#!/usr/bin/env python3

# 1) Run MadEvent using the options set in the input file to compute the widths and decays

from __future__ import print_function
import sys,os,glob
from configParserWrapper import ConfigParserExt
import logging,shutil
import subprocess
import tempfile
import time
import multiprocessing
import tqdm
import numpy as np

FORMAT = '%(levelname)s in %(module)s.%(funcName)s(): %(message)s at %(asctime)s'
logging.basicConfig(format=FORMAT,datefmt='%m/%d/%Y %I:%M:%S %p')
logger = logging.getLogger("MG5Scan")

def generateProcessFolder(parser):
    """
    Runs the madgraph process generation.
    This step just need to be performed once for a given
    model, since it is independent of the 
    numerical values of the model parameters.
    
    :param parser: Dictionary with parser sections.
    
    :return: True if successful. Otherwise False.
    """
    
    
    #Get run folder:    
    pars = parser["MadEventPars"]
    processCard = os.path.abspath(pars["proccard"])    
    if not os.path.isfile(processCard):
        logger.error("Process card %s not found" %processCard)
        raise ValueError()

    processFolder = os.path.abspath(pars["processFolder"])
    if os.path.isdir(processFolder):
        logger.warning("Process folder %s found. Skipping process generation." %processFolder)
        return False

    logger.info('Generating process using %s' %processCard)

    # Create copy of process card to replace output folder
    procCard = tempfile.mkstemp(suffix='.dat', prefix='procCard_')
    os.close(procCard[0])
    procCard = procCard[1]
    shutil.copy(processCard,procCard)
    with open(procCard,'r') as f:
        lines = f.readlines()
    lines = [l for l in lines[:] if l.strip()[:6] != 'output']
    lines.append('output %s\n' %processFolder)
    with open(procCard,'w') as f:
        for l in lines:
            f.write(l)
    
    #Generate process
    mg5Folder = os.path.abspath('./MG5')
    run = subprocess.Popen('./bin/mg5_aMC -f %s' %procCard,shell=True,
                                stdout=subprocess.PIPE,stderr=subprocess.PIPE,
                                cwd=mg5Folder)
         
    output,errorMsg = run.communicate()
    logger.debug('MG5 process error:\n %s \n' %errorMsg)
    logger.debug('MG5 process output:\n %s \n' %output)
    logger.info("Finished process generation")

    os.remove(procCard)
        
    return True

def runMadEvent(parser) -> str:
    
    """
    Runs the madgraph event generation and Pythia, if runPythia = True.
    
    :param parser: Dictionary with parser sections.
    :param runPythia: If True, will run the MG5-Pythia interface.
    :param runMadSpin: If True, will run the MadSpin.
    
    :return: Dictionary with run info. False if failed.
    """

    pars = parser["MadEventPars"]

    if not 'processFolder' in pars:
        logger.error('Process folder not defined.')
        return 'Error'        
    else:
        processFolder = pars['processFolder']
        if not os.path.isdir(processFolder):
            logger.error('Process folder %s not found.' %processFolder)
            return 'Error'
     
    if not 'runFolder' in pars:
        logger.error('Run folder not defined.')
        return 'Error'
    else:
        runFolder = pars['runFolder']

    if not 'outputFolder' in pars:
        logger.error('Output folder not defined.')
        return 'Error'
    else:
        outputFolder = pars['outputFolder']
        if not os.path.isdir(outputFolder):
            os.makedirs(outputFolder)

    # If run folder does not exist, create it using processFolder as a template:
    if not os.path.isdir(runFolder):
        runFolder = shutil.copytree(processFolder,runFolder,
                                    symlinks=True)
        logger.info("Created temporary folder %s" %runFolder) 

    if not os.path.isdir(runFolder):
        logger.error('Run folder %s not found.' %runFolder)
        return 'Error'

    # If run folder does not exist, create it using processFolder as a template:
    if 'paramcard' in pars:
        if os.path.isfile(pars['paramcard']):
            shutil.copyfile(pars['paramcard'],os.path.join(runFolder,'Cards/param_card.dat'))    
        else:
            raise ValueError("Param card %s not found" %pars['paramcard'])
    
    pdgs = str(pars['pdgs']).replace(',', ' ')
        
    #Generate commands file:       
    commandsFile = tempfile.mkstemp(suffix='.txt', prefix='MadEvent_commands_', dir=runFolder)
    os.close(commandsFile[0])
    commandsFileF = open(commandsFile[1],'w')
    commandsFileF.write(f'compute_widths {pdgs}\n')
    comms = parser["MadEventSet"]
    run_tag = comms.pop('run_tag')
    # Set the MadGraph parameters defined in the ini file
    for key,val in comms.items():
        commandsFileF.write('set %s %s\n' %(key,val))

    #Done setting up options
    commandsFileF.write('done\n')
    commandsFileF.close()
    commandsFile = commandsFile[1]      

    logger.debug("Computing widths events with command file %s" %commandsFile)
    run = subprocess.Popen('./bin/madevent < %s' %(commandsFile),
                           shell=True,stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE,cwd=runFolder)      
    output,errorMsg= run.communicate()    
    logger.debug('MG5 event error:\n %s \n' %errorMsg.decode())
    logger.debug('MG5 event output:\n %s \n' %output.decode())

    
    newFile = os.path.join(outputFolder,f'param_card_{run_tag}.dat')
    shutil.copyfile(os.path.join(runFolder,'Cards/param_card.dat'),newFile)
    logger.info(f'Copying {os.path.join(runFolder,'Cards/param_card.dat')} to {newFile}')
    
    if (runFolder != processFolder) and (runFolder != outputFolder):
        shutil.rmtree(runFolder)

    return newFile


def main(parfile,verbose):
   
    level = verbose
    levels = { "debug": logging.DEBUG, "info": logging.INFO,
               "warn": logging.WARNING,
               "warning": logging.WARNING, "error": logging.ERROR }
    if not level in levels:
        logger.error ( "Unknown log level ``%s'' supplied!" % level )
        sys.exit()
    logger.setLevel(level = levels[level])    

    parser = ConfigParserExt(inline_comment_prefixes="#")   
    ret = parser.read(parfile)
    if ret == []:
        logger.error( "No such file or directory: '%s'" % args.parfile)
        sys.exit()

    # Check if a parFile has been defined, if it has, create parsers from each row    
    if parser.has_option('AuxPars','parFile'):
        parFile = np.genfromtxt(parser.get('AuxPars','parFile'),
                                delimiter=',',names=True)
        parserList = []    
        for row in parFile:
            newParser = ConfigParserExt()
            newParser.read_dict(parser.toDict(raw=True))            
            for key,val in zip(parFile.dtype.names,row):
                newParser.set('AuxPars',key,str(val))
            parserList += newParser.expandLoops()
    else:
        #Get a list of parsers (in case loops have been defined)    
        parserList = parser.expandLoops()

    
    # Start multiprocessing pool
    ncpus = -1
    if parser.has_option("options","ncpu"):
        ncpus = int(parser.get("options","ncpu"))
    if ncpus  < 0:
        ncpus =  multiprocessing.cpu_count()
    ncpus = min(ncpus,len(parserList))
    pool = multiprocessing.Pool(processes=ncpus)
    if ncpus > 1:
        logger.info('Running %i jobs in parallel with %i processes' %(len(parserList),ncpus))
    else:
        logger.info('Running %i jobs in series with a single process' %(len(parserList)))
    
    children = []
    for irun,newParser in enumerate(parserList):
        processFolder = newParser.get('MadEventPars','processFolder')
        processFolder = os.path.abspath(processFolder)
        if processFolder[-1] == '/':
            processFolder = processFolder[:-1]
        if not os.path.isdir(processFolder):
            logger.info('Folder %s not found. Running MG5 to create folder.' %processFolder)
            generateProcessFolder(newParser)

        # Get largest existing events folder:
        run0 = 1
        eventsFolder = os.path.join(processFolder,'Events')
        if os.path.isdir(eventsFolder):
            for runF in glob.glob(os.path.join(eventsFolder,'run*')):
                run0 = max(run0,int(os.path.basename(runF).replace('run_',''))+1)

        # Create temporary folder names if running in parallel
        if ncpus > 1:
            # Create temporary folders
            runFolder = tempfile.mkdtemp(prefix='%s_'%(processFolder),suffix='_run_%02d' %(run0+irun))
            os.removedirs(runFolder)
        else:
            runFolder = processFolder

        newParser.set('MadEventPars','runFolder',runFolder)
        parserDict = newParser.toDict(raw=False)
        logger.debug('submitting with pars:\n %s \n' %parserDict)
        p = pool.apply_async(runMadEvent, args=(parserDict,))                       
        children.append(p)

        logger.debug('submitting with pars:\n %s \n' %parserDict)

    #     Wait for jobs to finish:
    output = []
    for p in tqdm.tqdm(children):
        output.append(p.get())
    
    return output
    


if __name__ == "__main__":
    
    import argparse    
    ap = argparse.ArgumentParser( description=
            "Run MadEvent in serial to compute widths for the parameters defined in the parameters file." )
    ap.add_argument('-p', '--parfile', default='parameters_lifetimeScan.ini',
            help='path to the parameters file [parameters_lifetimeScan.ini].')
    ap.add_argument('-v', '--verbose', default='info',
            help='verbose level (debug, info, warning or error). Default is info')


    
    t0 = time.time()

    args = ap.parse_args()
    output = main(args.parfile,args.verbose)
            
    print("\n\nDone in %3.2f min" %((time.time()-t0)/60.))
