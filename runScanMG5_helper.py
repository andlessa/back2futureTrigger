#!/usr/bin/env python3

# 1) Run MadGraph using the options set in the input file 
# (the proc_card.dat, parameter_card.dat and run_card.dat...).

from __future__ import print_function
import os,glob
import logging,shutil
import subprocess
import tempfile
import time,sys,datetime
from configParserWrapper import ConfigParserExt
from typing import Dict,Set,Union, NamedTuple
from collections import namedtuple
import pylhe
folderTuple = namedtuple('folderTuple', ['inputFolder', 'outputFolder'])

FORMAT = '%(levelname)s: %(message)s at %(asctime)s'
logging.basicConfig(format=FORMAT,datefmt='%m/%d/%Y %I:%M:%S %p')
logger = logging.getLogger("MG5Scan")

def generateInputFiles(parfile) -> Set[folderTuple]:
   
    parser = ConfigParserExt(inline_comment_prefixes="#")   
    ret = parser.read(parfile)
    if ret == []:
        logger.error( f"No such file or directory: {parfile}")
        return set([])
            
    #Get a list of parsers (in case loops have been defined)    
    parserList = parser.expandLoops()
    now = datetime.datetime.now()
    scanFolders = set([])
    for irun,newParser in enumerate(parserList):
        processFolder = newParser.get('MadGraphPars','processFolder')
        processFolder = os.path.abspath(processFolder)
        inputFolder = os.path.join(processFolder,'scan_inputFiles')
        outputFolder = os.path.join(processFolder,'scan_results')
        outputFolder = os.path.abspath(outputFolder)
        fTuple = folderTuple(inputFolder,outputFolder)
        if not os.path.isdir(processFolder):
            logger.info('Folder %s not found. Running MG5 to create folder.' %processFolder)
            generateProcess(newParser)
        # If it is the first time going through this process folder, create the input and output folders. 
        # Otherwise, skip this step to avoid overwriting existing input files.
        if fTuple not in scanFolders:
            if os.path.isdir(inputFolder):
                logger.warning(f'Input folder {inputFolder} already exists. It will be overwritten.')
                shutil.rmtree(inputFolder)
            if os.path.isdir(outputFolder):
                logger.warning(f'Output folder {outputFolder} already exists. It will be overwritten.')
                shutil.rmtree(outputFolder)

            logger.debug(f'Creating input folder {inputFolder} and output folder {outputFolder}')
            os.makedirs(inputFolder, exist_ok=True)
            os.makedirs(outputFolder, exist_ok=True)

        # Get largest existing events folder:
        run0 = 1
        eventsFolder = os.path.join(processFolder,'Events')
        if os.path.isdir(eventsFolder):
            for runF in glob.glob(os.path.join(eventsFolder,'run*')):
                run0 = max(run0,int(os.path.basename(runF).replace('run_',''))+1)

        # Create temporary folder
        runFolder = tempfile.mkdtemp(suffix='_run_%02d' %(run0+irun),
                                     dir=outputFolder)
        os.rmdir(runFolder)
        
        newParser.set('MadGraphPars','runFolder',os.path.abspath(runFolder))
        newParser.set('MadGraphPars','runNumber','%02d' %(run0+irun))
        newParser.set('MadGraphPars','resultsFolder',outputFolder)

        parserDict = newParser.toDict(raw=False,abspath_existing=True)
        # Create text config file to store the dictionary of parameters for this run.
        # This file can be read again by ConfigParserExt in the worker.
        outfile = f"{inputFolder}/job_{irun:05d}.ini"
        parserOut = ConfigParserExt()
        parserOut.read_dict(parserDict)
        with open(outfile, "w") as f:
            parserOut.write(f)
            
        scanFolders.add(fTuple)
    logger.info(f"Created {len(parserList)} input files at {now.strftime('%Y-%m-%d %H:%M')}")

    return scanFolders

def generateProcess(parser):
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
    logger.debug('MG5 process error:\n %s \n' %errorMsg.decode())
    logger.debug('MG5 process output:\n %s \n' %output.decode())
    logger.info("Finished process generation")

    os.remove(procCard)
        
    return True

def getInfoFromOutput(outputStr):
    """
    Try to fetch the run name, the number of events
    and the total cross-section from the MG5 output.

    :param outputStr: String containing the MG5 output

    :return: Dictionary with run summary
    """

    outputStr = str(outputStr)
    # Get summary block
    summaryBlock = outputStr.split("Results Summary for run")[-1]
    summaryBlock = summaryBlock.split("Done")[0]
    # Split block into lines
    summaryLines = [l for l in summaryBlock.split('\\n') if l.strip()]
    # Get info from lines:
    xsec, runNumb, runTag, nevts = None,None,None,-1
    for l in summaryLines:
        if 'tag' in l:
            runNumb,runTag = l.split('tag:')
            runNumb = runNumb.replace(':','').replace('=','').strip()
            runTag = runTag.replace(':','').replace('=','').strip()
        elif 'cross-section' in l.lower():
            xsec = eval(l.lower().split(':')[1].split('+-')[0])
        elif 'nb of events'in l.lower():
            nevts = eval(l.lower().split(':')[1])

    runInfo = {'run number' : runNumb, 'run tag' : runTag, 
               'cross-section (pb)' : xsec, 'Number of events' : int(nevts)}
    return runInfo

def runMG5(parser,runPythia=False, runMadSpin=False) -> Dict:
    
    """
    Runs the madgraph event generation and Pythia, if runPythia = True.
    
    :param parser: Dictionary with parser sections.
    :param runPythia: If True, will run the MG5-Pythia interface.
    :param runMadSpin: If True, will run the MadSpin.
    
    :return: Dictionary with run info. False if failed.
    """

    pars = parser["MadGraphPars"]
    if not 'runFolder' in pars:
        logger.error('Run folder not defined.')
        return {}        
    else:
        runFolder = pars['runFolder']

    if not 'processFolder' in pars:
        logger.error('Process folder not defined.')
        return {}        
    else:
        processFolder = pars['processFolder']
        if not os.path.isdir(processFolder):
            logger.error('Process folder %s not found.' %processFolder)
            return {}            

    # If run folder does not exist, create it using processFolder as a template:
    if not os.path.isdir(runFolder):
    # To avoid recursion when copying the process folder, we first copy the process folder 
    # to a temporary location and then move it to the results folder.
        run_tmp = shutil.copytree(processFolder,os.path.basename(runFolder),
                                    ignore=shutil.ignore_patterns('Events','*.lhe','scan_inputFiles',
                                                                  'scan_results','*.root','*.hepmc'),
                                    symlinks=True)
        os.makedirs(os.path.join(run_tmp,'Events'))
        shutil.move(run_tmp, runFolder)
       
        logger.info("Created temporary folder %s" %runFolder) 
            
    if not os.path.isdir(runFolder):
        logger.error('Run folder %s not found.' %runFolder)
        return {}

    if 'runcard' in pars:
        if os.path.isfile(pars['runcard']):    
            shutil.copyfile(pars['runcard'],os.path.join(runFolder,'Cards/run_card.dat'))
        else:
            raise ValueError("Run card %s not found" %pars['runcard'])
    if 'paramcard' in pars:
        if os.path.isfile(pars['paramcard']):
            shutil.copyfile(pars['paramcard'],os.path.join(runFolder,'Cards/param_card.dat'))    
        else:
            raise ValueError("Param card %s not found" %pars['paramcard'])
    

    # By default do not run Pythia or Delphes
    if runPythia and 'pythia8card' in pars:
        pythia8File = os.path.join(runFolder,'Cards/pythia8_card.dat')
        if os.path.isfile(pars['pythia8card']):
            shutil.copyfile(pars['pythia8card'],pythia8File) 

    if runMadSpin:
        if not 'madspincard' in pars or not os.path.isfile(pars['madspincard']):
            logger.error("MadSpin file not defined or not found.")
            return {}
        else:
            madspinFile = os.path.join(runFolder,'Cards/madspin_card.dat')
            shutil.copyfile(pars['madspincard'],madspinFile)   
        
    #Generate commands file:       
    commandsFile = tempfile.mkstemp(suffix='.txt', prefix='MG5_commands_', dir=runFolder)
    os.close(commandsFile[0])
    commandsFileF = open(commandsFile[1],'w')
    if runPythia:
        commandsFileF.write('shower=Pythia8\n')
    else:
        commandsFileF.write('shower=OFF\n')
    commandsFileF.write('detector=OFF\n')
    if runMadSpin:
        commandsFileF.write('madspin=ON\n')
    else:
        commandsFileF.write('madspin=OFF\n')
    commandsFileF.write('done\n')
    comms = parser["MadGraphSet"]
    # Set the MadGraph parameters defined in the ini file
    for key,val in comms.items():
        commandsFileF.write('set %s %s\n' %(key,val))

    #Done setting up options
    commandsFileF.write('done\n')
    commandsFileF.write('done\n')

    commandsFileF.close()
    commandsFile = commandsFile[1]      

    ncore = parser['options']['ncore']
    
    logger.debug("Generating MG5 events with command file %s" %commandsFile)
    # Fix the Python hash seed: MG5 stores the helicity filter (SubProcesses/P*/Hel/selection)
    # as str(set), whose ordering is random between Python sessions. This makes MG5 believe the
    # helicity selection changed and recompile all matrix*_optim.f files in every run. Only needed if hel_recycling = True
    #env = os.environ.copy()
    #env['PYTHONHASHSEED'] = '0'
    run = subprocess.Popen(f'./bin/generate_events --multicore --nb_core={ncore} < {commandsFile}',
                           shell=True,stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE,cwd=runFolder
                        #    ,env=env
                           )
      
    output,errorMsg= run.communicate()
    runInfo = {}
    runInfo.update(pars)

    logger.debug('MG5 event error:\n %s \n' %errorMsg.decode())
    logger.debug('MG5 event output:\n %s \n' %output.decode())
    # Get info from output
    runInfo.update(getInfoFromOutput(output))

    os.remove(commandsFile)

    return runInfo

def runDelphes(parser,runInfo,runDelphesPythia=True) -> Dict:
    """
    Runs Delphes on a hepmc file or directly (bypassing hepmc generation)
    on the LHE file.
    
    :param parser: Dictionary with parser sections.
    :param runInfo: Dictionary with MadGraph run information.
    :param runDelphesPythia: If True, run DelphesPythia8 on an input LHE file, otherwise run DelphesHepMC2 on an input HEPMC file.
    
    :return: Dictionary with run info. False if failed.
    """

    cleanOutput = parser['options']['cleanOutput']
    pars = parser["DelphesPars"]

    runFolder = parser['MadGraphPars']['runFolder']
    # Allow for running delphes multiple times with distinct cards
    delphescardList = pars['delphescard'].split(',')
    
    lheFile = None
    hepmcFile = None
    for delphescard in delphescardList:
        delphescard = os.path.abspath(delphescard)
        cardName = os.path.basename(delphescard)
        cardName = os.path.splitext(cardName)[0]
        rootFile = os.path.join(runFolder,'Events',runInfo['run number'], f"{runInfo['run tag']}_{cardName}.root")
        
        delphesFile = os.path.join(runFolder,'Events',runInfo['run number'],os.path.basename(delphescard))
        shutil.copyfile(delphescard,delphesFile)
        eventsFolder = os.path.join(runFolder,'Events',runInfo['run number'])

        if runDelphesPythia:
            # Get possible input files:
            lheFile = os.path.join(eventsFolder,'unweighted_events.lhe.gz')
            if not os.path.isfile(lheFile):
                logger.error(f'LHE file {lheFile} not found.')
                return runInfo
            if _getOpt(pars, 'matching', False):
                runInfo = runDelphesPythiaMatched(parser, runInfo, delphescard,
                                              lheFile, rootFile)
            else:
                runInfo = runDelphesPythiaDefault(parser, runInfo, delphescard, lheFile, rootFile)
        else:   
            hepmcFiles = list(glob.glob(os.path.join(eventsFolder,'*hepmc*')))
            if not hepmcFiles:
                logger.error(f'HEPMC file not found in {runFolder}')
                return runInfo
            if len(hepmcFiles) > 1:
                logger.warning(f'Found {len(hepmcFiles)} candidate files for the HEPMC file. Using {hepmcFiles[0]}.')
            hepmcFile = hepmcFiles[0]
            runInfo = runDelphesHepMC(parser, runInfo, delphescard, hepmcFile, rootFile)

        #Generate pythia log file
        if runDelphesPythia:
            pythia_log = rootFile.replace('.root','_pythia_delphes.log')
            with open(pythia_log,'w') as f:
                f.write(runInfo['DelphesOutput'])

    if cleanOutput:
        for f in [hepmcFile,lheFile]:
            if not f:
                continue
            if os.path.isfile(f):
                os.remove(f)
        

    return runInfo

def _getOpt(section, key, default=None):
    """
    Returns section[key] if present, otherwise default. Works for dictionaries and
    ConfigParserExt sections (whose .get(key,default) is not supported).
    """
    if key in section:
        return section[key]
    return default

def _readRunCardValue(runCardFile : str, name : str) -> Union[str,None]:
    """
    Reads the value of a parameter from a MadGraph run_card.dat (format: 'value = name ! comment').
    """
    if not os.path.isfile(runCardFile):
        return None
    with open(runCardFile,'r') as f:
        for line in f:
            line = line.split('!')[0].split('#')[0]
            if '=' not in line:
                continue
            val,key = line.split('=',1)
            if key.strip().lower() == name.lower():
                return val.strip()
    return None


def runDelphesPythiaMatched(parser, runInfo : Dict, delphescard : str,
                            lheFile : str, rootFile : str) -> Dict:
    """
    Runs Pythia8 with MLM matching (JetMatchingMadgraph, same settings as the MG5-Pythia8
    interface) and Delphes in a single step, without writing HepMC files.
    The LHE file is split into nchunks slices (via Beams:nSkipLHEFatInit), which are processed
    in parallel by DelphesPythia8Matched. The output ROOT files are merged with hadd and the
    matched cross section is written to <run_tag>_merged_xsecs.txt.

    Optional DelphesPars entries: delphesPythiaExe, nchunks (default = options:ncore),
    qcut (default = 1.5*xqcut), njetmax.

    :return: Updated runInfo dictionary.
    """

    pars = parser["DelphesPars"]
    mgSet = parser["MadGraphSet"]
    runFolder = parser['MadGraphPars']['runFolder']
    eventsFolder = os.path.dirname(rootFile)
    delphesDir = os.path.abspath(pars['delphesDir'])
    runCard = os.path.join(runFolder,'Cards','run_card.dat')

    defaultExe = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                              'delphesPythiaMatching','DelphesPythia8Matched')
    exe = os.path.abspath(_getOpt(pars, 'delphesPythiaExe', defaultExe))
    if not os.path.isfile(exe):
        logger.error(f'{exe} not found. Compile it with "make" in delphesPythiaMatching.')
        return runInfo

    # Pythia card (DelphesPars has priority over MadGraphPars)
    pythiaCard = _getOpt(pars, 'pythia8card', _getOpt(parser['MadGraphPars'], 'pythia8card'))
    pythiaCard = os.path.abspath(pythiaCard)
    with open(pythiaCard,'r') as f:
        baseCard = f.read()

    # Matching parameters, following MG5 (madevent_interface.py) defaults for MLM:
    xqcut = _getOpt(mgSet, 'xqcut', _readRunCardValue(runCard,'xqcut'))
    xqcut = float(xqcut)
    qcut = float(_getOpt(pars, 'qcut', 1.5*xqcut))
    if qcut < 1.5*xqcut:
        logger.warning(f'qCut = {qcut} is smaller than 1.5*xqcut = {1.5*xqcut}')
    nQmatch = int(_getOpt(mgSet, 'maxjetflavor', _readRunCardValue(runCard,'maxjetflavor') or 4))
    njetmax = int(pars['njetmax'])
    iseed = int(_getOpt(mgSet, 'iseed', 0))

    # Split the LHE file
    nLHE = pylhe.LHEFile.count_events(lheFile)
    nchunks = int(_getOpt(pars, 'nchunks', _getOpt(parser['options'], 'ncore', 1)))
    nchunks = max(1,min(nchunks,nLHE))
    chunkSize = -(-nLHE // nchunks)  # ceil
    logger.info(f'Running DelphesPythia8Matched on {nLHE} LHE events in {nchunks} chunks '
                f'(qCut = {qcut:.1f} GeV, nJetMax = {njetmax}, nQmatch = {nQmatch})')

    chunks = []
    for ichunk in range(nchunks):
        nskip = ichunk*chunkSize
        if nskip >= nLHE:
            break
        tag = os.path.splitext(rootFile)[0] + f'_chunk{ichunk}'
        chunk = {'card' : tag + '_pythia8.cmnd', 'root' : tag + '.root',
                 'summary' : tag + '_summary.txt', 'log' : tag + '.log'}
        with open(chunk['card'],'w') as f:
            f.write(baseCard)
            f.write('\n\n### Added by runDelphesPythiaMatched:\n')
            f.write('Beams:frameType = 4\n')
            f.write(f'Beams:LHEF = {lheFile}\n')
            f.write('Beams:setProductionScalesFromLHEF = on\n')
            f.write(f'Beams:nSkipLHEFatInit = {nskip}\n')
            f.write(f'Main:numberOfEvents = {min(chunkSize, nLHE-nskip)}\n') # number of LHE events
            f.write('Random:setSeed = on\n')
            f.write(f'Random:seed = {(iseed*1000 + ichunk + 1) % 900000000}\n')
            f.write('JetMatching:merge = on\n')
            f.write('JetMatching:scheme = 1\n')
            f.write('JetMatching:setMad = off\n')
            f.write(f'JetMatching:qCut = {qcut:.4f}\n')
            f.write(f'JetMatching:nJetMax = {njetmax}\n')
            f.write(f'JetMatching:nQmatch = {nQmatch}\n')
            f.write('JetMatching:coneRadius = 1.0\n')
            f.write('JetMatching:etaJetMax = 1000.0\n')
        for f in [chunk['root'],chunk['summary']]:
            if os.path.isfile(f):
                os.remove(f)
        chunk['logF'] = open(chunk['log'],'w')
        chunk['proc'] = subprocess.Popen([exe, delphescard, chunk['card'], chunk['root'], chunk['summary']],
                                         stdout=chunk['logF'], stderr=subprocess.STDOUT,
                                         cwd=delphesDir)
        chunks.append(chunk)

    # Wait for all chunks and collect summaries
    totals = {'nLHE' : 0, 'nAccepted' : 0, 'nVetoed' : 0, 'nErrors' : 0}
    sigmaLHE = None
    failed = False
    for chunk in chunks:
        ret = chunk['proc'].wait()
        chunk['logF'].close()
        if ret != 0 or not os.path.isfile(chunk['summary']):
            logger.error(f"DelphesPythia8Matched failed (exit code {ret}). See {chunk['log']}")
            failed = True
            continue
        summary = {}
        with open(chunk['summary'],'r') as f:
            for line in f:
                key,val = line.split('=')
                summary[key.strip()] = float(val)
        for key in totals:
            totals[key] += int(summary[key])
        sigmaLHE = summary['sigmaLHE_pb']

    # Concatenate logs
    pythia_log = rootFile.replace('.root','_pythia_delphes.log')
    with open(pythia_log,'w') as fout:
        for chunk in chunks:
            if os.path.isfile(chunk['log']):
                fout.write(f"\n######## {os.path.basename(chunk['card'])} ########\n")
                with open(chunk['log'],'r') as fin:
                    shutil.copyfileobj(fin,fout)

    if failed:
        runInfo.update({'DelphesError' : f'DelphesPythia8Matched failed. See {pythia_log}'})
        return runInfo

    # Merge ROOT files
    rootChunks = [chunk['root'] for chunk in chunks]
    if os.path.isfile(rootFile):
        os.remove(rootFile)
    if len(rootChunks) == 1:
        shutil.move(rootChunks[0],rootFile)
    else:
        run = subprocess.run(['hadd','-f',rootFile] + rootChunks,
                             capture_output=True, text=True)
        if run.returncode != 0:
            logger.error(f'hadd failed:\n{run.stderr}')
            runInfo.update({'DelphesError' : run.stderr})
            return runInfo

    # Matched cross section (unweighted LHE events): sigma_LHE * N_accepted/N_LHE
    eff = totals['nAccepted']/totals['nLHE']
    effErr = (eff*(1.0-eff)/totals['nLHE'])**0.5
    xsec, xsecErr = sigmaLHE*eff, sigmaLHE*effErr
    logger.info(f"Matching efficiency = {eff:.4f} +- {effErr:.4f} ({totals['nAccepted']}/{totals['nLHE']}, "
                f"{totals['nErrors']} errors). Matched cross section = {xsec:.4e} +- {xsecErr:.1e} pb")

    # Same file name and layout (2nd column = xsec in pb) as the MG5 merged_xsecs.txt
    xsecFile = os.path.join(eventsFolder, f"{runInfo['run tag']}_merged_xsecs.txt")
    with open(xsecFile,'w') as f:
        f.write('# Matched cross section from DelphesPythia8Matched (MLM, sigma_LHE*nAccepted/nLHE)\n')
        f.write('# qCut[GeV]  sigma[pb]  error[pb]  sigma_LHE[pb]  nLHE  nAccepted\n')
        f.write(f"{qcut:.2f}  {xsec:.6e}  {xsecErr:.6e}  {sigmaLHE:.6e}  {totals['nLHE']}  {totals['nAccepted']}\n")

    # Clean up chunk files
    for chunk in chunks:
        for key in ['card','root','summary','log']:
            if os.path.isfile(chunk[key]):
                os.remove(chunk[key])

    runInfo.update({'LHE cross-section (pb)' : sigmaLHE,
                    'cross-section (pb)' : xsec,
                    'cross-section error (pb)' : xsecErr,
                    'matching efficiency' : eff,
                    'Number of accepted events' : totals['nAccepted'],
                    'DelphesOutput' : pythia_log, 'DelphesError' : ''})

    return runInfo

def runDelphesPythiaDefault(parser,runInfo,delphescard,lheFile,rootFile) -> Dict:

    pars = parser["DelphesPars"]
    runFolder = parser['MadGraphPars']['runFolder']
    delphesDir = os.path.abspath(pars['delphesDir'])
    nevts = parser["MadGraphSet"]["nevents"]
    pythiaCard = os.path.abspath(pars['pythia8card'])
    pythiaFile = os.path.join(runFolder,'Events',runInfo['run number'],'pythia8_card.dat')
    shutil.copyfile(pythiaCard,pythiaFile)    
    #Generate pythia card with additional info
    with open(pythiaFile,'a') as f:
        f.write('\n\n')
        f.write('### Added lines:\n')
        f.write('Beams:frameType = 4\n')
        f.write('Beams:LHEF = %s\n' %lheFile)
        f.write('Beams:setProductionScalesFromLHEF=on\n')
        f.write('Main:numberOfEvents   = %i\n' %nevts)
        f.write('Main:timesAllowErrors   = %i\n' %int(1e-2*nevts)) # Allow at most 1% of failures
        # Matching specific parameters:
        if 'matching' in pars and pars['matching']:
            qcut = 1.5*float(parser["MadGraphSet"]['xqcut'])
            njet = pars['njetmax']
            f.write('JetMatching:merge=on \n')
            f.write('JetMatching:qCut         = %1.2f \n' %qcut)
            f.write('JetMatching:nJetMax      = %i \n' %njet)

    logger.debug("Running DelphesPythia8 with config files %s and %s" %(pythiaFile,delphescard))
    run = subprocess.Popen('./DelphesPythia8 %s %s %s' %(delphescard,pythiaFile,rootFile),shell=True,
                            stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True,
                            cwd=delphesDir)

    output,errorMsg = run.communicate()

    runInfo.update({'DelphesOutput' : output, 'DelphesError' : errorMsg})
    
    return runInfo

def runDelphesHepMC(parser,runInfo,delphescard,hepmcFile,rootFile) -> Dict:

    pars = parser["DelphesPars"]
    delphesDir = os.path.abspath(pars['delphesDir'])
    logger.debug("Running DelphesHepMC2 with files %s and %s" %(delphescard,hepmcFile))
    if os.path.splitext(hepmcFile)[-1] == '.gz':
        logger.info(f'HEPMC file {hepmcFile} is gzipped. Unzipping...')
        run = subprocess.Popen('gunzip -f %s' %hepmcFile,shell=True,
                                stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True)
        output,errorMsg = run.communicate()
        hepmcFile = os.path.splitext(hepmcFile)[0]
        logger.info(f'Unzipped HEPMC file to {hepmcFile}')
    run = subprocess.Popen('./DelphesHepMC2 %s %s %s' %(delphescard,rootFile,hepmcFile),shell=True,
                            stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True,
                            cwd=delphesDir)

    output,errorMsg = run.communicate()

    runInfo.update({'DelphesOutput' : output, 'DelphesError' : errorMsg})

    return runInfo


def generateEvents(parser):
    """
    Run MadGraph5 to generate LHE events and then DelphesPythia8.

    :param parser: Dictionary with parser sections.
    
    :return: Dictionary with run info. False if failed.
    """

    t0 = time.time()
    # Defaul settings:
    run_pythia = False
    run_delphespythia = False
    run_delphes = False
    run_madgraph = False
    run_madspin = False
    runInfo = {}

    if isinstance(parser,str):
        configFile = parser
        parser = ConfigParserExt(inline_comment_prefixes="#")
        ret = parser.read(configFile)
        if ret == []:
            logger.error(f"Could not read config file: {configFile}")
            return runInfo
        parserDict = parser.toDict(raw=False,abspath_existing=True)
    elif isinstance(parser,Dict):
        parserDict = parser
    else:
        logger.error(f"Parser should be either a file path or a dictionary and not  {type(parser).__name__}.")
        return runInfo

    assert isinstance(parserDict,Dict), f"Parser should be a dictionary and not {type(parserDict).__name__}."
    if isinstance(parserDict['options'], dict):        
        run_madgraph = parserDict["options"].get("runMadGraph", False)    
        run_pythia = parserDict["options"].get("runPythia", False)
        run_delphes = parserDict["options"].get("runDelphes", False)
        run_delphespythia = parserDict["options"].get("runDelphesPythia", False)
        run_madspin = parserDict["options"].get("runMadSpin", False)

    if run_pythia and run_delphespythia:
        logger.warning('Both runPythia and runDelphesPythia set to True. Setting runDelphesPythia to False')
        run_delphespythia = False
    

    if run_madgraph:
        runInfo = runMG5(parser,runPythia=run_pythia,
                         runMadSpin=run_madspin)
        if run_delphes or run_delphespythia:
            runInfo = runDelphes(parser,runInfo,
                             runDelphesPythia=run_delphespythia)
        
        runInfo.update({'time (s)' : time.time()-t0})
        return runInfo
    elif run_delphes or run_delphespythia or run_pythia:
        logger.error("Pythia and Delphes can only run if runMadGraph = True")
        return {}
    
def moveFolders(runInfo):
    """
    Move the run folders from the temporary running folder
    to the process folder.
    """

    logger.info('Finished event generation for run %i in %1.2f min' 
                %(int(runInfo['runNumber']),runInfo['time (s)']/60.))

    # Get run folder:
    runFolder = os.path.abspath(runInfo['runFolder'])
    runNumber = int(runInfo['runNumber'])
    processFolder = os.path.abspath(runInfo['processFolder'])
    # If run folder and process folder are the same, do nothing
    if runFolder == processFolder:
        return
    
    # Move run folder results to process folder
    eventFolder = os.path.join(processFolder,'Events')
    runDirs = list(glob.glob(os.path.join(runFolder,'Events','run_*')))
    if len(runDirs) != 1:
        logger.error('Something went wrong. Found %i run folders in %s' %(len(runDirs),runFolder))
        return False
    runDir = runDirs[0]
    finalRunDir = os.path.join(eventFolder,'run_%02d' %runNumber)
    logger.info('Moving %s to %s' %(runDir,finalRunDir))
    shutil.move(runDir,finalRunDir)
    logger.info('Deleting temporary folder %s' %(runFolder))
    shutil.rmtree(runFolder)

    
