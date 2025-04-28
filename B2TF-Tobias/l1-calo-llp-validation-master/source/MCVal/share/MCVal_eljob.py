#!/usr/bin/env python

# Read the submission directory as a command line argument. You can
# extend the list of arguments with your private ones later on.
import optparse
parser = optparse.OptionParser()
parser.add_option( '-s', '--submission-dir', dest = 'submission_dir',
                   action = 'store', type = 'string', default = 'submitDir',
                   help = 'Submission directory for EventLoop' )
parser.add_option( '-i', '--inputfile', dest = 'inputfile',
                   action = 'store', type = 'string', default = '',
                   help = 'Input filename' )
parser.add_option( '--scalarPDGID', dest = 'scalarPDGID',
                   action = 'store', type = 'int', default = 55,
                   help = 'PDGID of BSM scalar particle' )
parser.add_option( '--chi1PDGID', dest = 'chi1PDGID',
                   action = 'store', type = 'int', default = 4000023,
                   help = 'PDGID of BSM Chi1 particle' )
parser.add_option( '--chi0PDGID', dest = 'chi0PDGID',
                   action = 'store', type = 'int', default = 4000022,
                   help = 'PDGID of BSM Chi0 particle' )
( options, args ) = parser.parse_args()

# Set up (Py)ROOT.
import ROOT
ROOT.xAOD.Init().ignore()

# Set up the sample handler object. See comments from the C++ macro
# for the details about these lines.
import os
sh = ROOT.SH.SampleHandler()
sh.setMetaString( 'nc_tree', 'CollectionTree' )
inputFilePath = os.getcwd()+''

ROOT.SH.ScanDir().filePattern( options.inputfile ).scan( sh, inputFilePath )
sh.printContent()

# Create an EventLoop job.
job = ROOT.EL.Job()
job.sampleHandler( sh )
# This lets you run over a defined number of events. It is useful for debugging
#job.options().setDouble( ROOT.EL.Job.optMaxEvents, 100 )
job.outputAdd (ROOT.EL.OutputStream ('ANALYSIS'))

# Create the algorithm's configuration.
from AnaAlgorithm.DualUseConfig import createAlgorithm
alg = createAlgorithm ( 'MCValAlg', 'AnalysisAlg' )
alg.pdgIdBSM = options.scalarPDGID
alg.pdgIdHS = options.chi1PDGID
alg.pdgIdchi0 = options.chi0PDGID

print(options)

# later on we'll add some configuration options for our algorithm that go here

# Add our algorithm to the job
job.algsAdd( alg )

# Run the job using the direct driver.
driver = ROOT.EL.DirectDriver()
driver.submit( job, options.submission_dir )
