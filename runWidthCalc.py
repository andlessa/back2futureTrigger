#!/usr/bin/env python3

# 1) Run MadEvent using the options set in the input file to compute the widths and decays

from __future__ import print_function
import sys,os,glob
from configParserWrapper import ConfigParserExt
import logging,shutil
import subprocess
import tempfile
import time,datetime
import tqdm

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
    pars = parser["MadGraphPars"]
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

    # If run folder does not exist, create it using processFolder as a template:
    if 'paramcard' in pars:
        if os.path.isfile(pars['paramcard']):
            shutil.copyfile(pars['paramcard'],os.path.join(processFolder,'Cards/param_card.dat'))    
        else:
            raise ValueError("Param card %s not found" %pars['paramcard'])
    
    pdgs = str(pars['pdgs']).replace(',', ' ')
        
    #Generate commands file:       
    commandsFile = tempfile.mkstemp(suffix='.txt', prefix='MadEvent_commands_', dir=processFolder)
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
                           stderr=subprocess.PIPE,cwd=processFolder)      
    output,errorMsg= run.communicate()    
    logger.debug('MG5 event error:\n %s \n' %errorMsg.decode())
    logger.debug('MG5 event output:\n %s \n' %output.decode())

    
    newFile = os.path.join(processFolder,f'width_results/param_card_{run_tag}.dat')
    if not os.path.isdir(os.path.join(processFolder,'width_results')):
        os.mkdir(os.path.join(processFolder,'width_results'))
    shutil.copyfile(os.path.join(processFolder,'Cards/param_card.dat'),newFile)
    
    os.remove(commandsFile)

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
            
    #Get a list of parsers (in case loops have been defined)    
    parserList = parser.expandLoops()

    
    for irun,newParser in enumerate(tqdm.tqdm(parserList)):
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

        
        parserDict = newParser.toDict(raw=False)
        logger.debug('submitting with pars:\n %s \n' %parserDict)
        _ = runMadEvent(parserDict)

    return
    


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
