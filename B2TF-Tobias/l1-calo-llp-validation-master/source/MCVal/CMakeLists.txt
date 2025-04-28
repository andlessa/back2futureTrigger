# The name of the package:
atlas_subdir (MCVal)

# Add the shared library:
atlas_add_library (MCValLib
  MCVal/*.h Root/*.cxx
  PUBLIC_HEADERS MCVal
  LINK_LIBRARIES AnaAlgorithmLib xAODTruth xAODJet xAODMissingET TruthUtils)

if (XAOD_STANDALONE)
  # Add the dictionary (for AnalysisBase only):
  atlas_add_dictionary (MCValDict
    MCVal/MCValDict.h
    MCVal/selection.xml
    LINK_LIBRARIES MCValLib)
endif ()

if (NOT XAOD_STANDALONE)
  # Add a component library for AthAnalysis only:
  atlas_add_component (MCVal
    src/components/*.cxx
    LINK_LIBRARIES MCValLib)
endif ()

# Install files from the package:
#atlas_install_joboptions( share/*_jobOptions.py )
atlas_install_scripts( share/*_eljob.py )

