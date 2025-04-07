#!/usr/bin/env python3

from typing import Any,Union
import numpy as np
import os
import glob
import pyslha
from xml.etree import ElementTree as ET




def getModelDict(inputFile,model='minimalH',verbose=True,bannerFile=None):

    if model == 'minimalH':
        LLP = 4000023
        LSP = 4000022
        mother = 55        
    else:
        raise ValueError("Unreconized model %s" %model)

    modelInfoDict = {}
    f = inputFile
    if not os.path.isfile(f):
        print('File %s not found' %f)
        raise OSError()
    parsDict = {}    
    if bannerFile is None:
        bannerFile = list(glob.glob(os.path.join(os.path.dirname(f),'*banner*txt')))[0]
    with open(bannerFile,'r') as ff:
        slhaData = ff.read().split('<slha>')[1].split('</slha>')[0]
        slhaData = pyslha.readSLHA(slhaData)
            
    parsDict = {}
    parsDict['m1'] = slhaData.blocks['MASS'][LLP]
    parsDict['m0'] = slhaData.blocks['MASS'][LSP]
    parsDict['mS'] = slhaData.blocks['MASS'][mother]
    
    modelInfoDict.update(parsDict)
    if verbose:
        print('ms = ',parsDict['mS'])
        print('m1 = ',parsDict['m1'])
        print('m0 = ',parsDict['m0'])
    
    return modelInfoDict