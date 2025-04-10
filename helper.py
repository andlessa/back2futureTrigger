#!/usr/bin/env python3

from typing import Any,Union
import numpy as np
import os
import glob
import pyslha
from xml.etree import ElementTree as ET



texDict = {'$m_S$ = %1.0f GeV' : 'ms',
           '$m_1 = $ %1.0f GeV' : 'm1',
        #    '$m_0 = $ %1.2f GeV' : 'm0',
           '$\\Delta m_{10} = $ %1.0f GeV' : 'dm',
           '$c \\tau = $ %1.0f m' : 'ctau'}

axisDict = {'$m_S$ (GeV)' : 'ms',
           '$m_1$ (GeV)' : 'm1',
        #    '$m_0 = $ %1.2f GeV' : 'm0',
           '$\\Delta m_{10}$ (GeV)' : 'dm',
           '$c \\tau (m)$' : 'ctau',
           '$\\epsilon_{\\rm T}$' : 'eff',
           '$Z_{\\rm exp} = \\frac{n_s}{\sigma_b}$' : 'Z0',
           '$\sigma_{\\rm UL}^{\\rm exp} (fb)$' : 'sigmaUL',
           '$n_{s}$' : 'ns'}

varDict = {val : key for key,val in axisDict.items()}

def getTitleFromDF(df):

    title = []
    for texLabel,label in texDict.items():
        if len(df[label].unique()) > 1:
            continue
        val = df[label].unique()[0]
        title.append((texLabel % val))

    return r', '.join(title)

def getAxisLabel(df):
    
    axVars = []
    for texLabel,label in axisDict.items():
        if len(df[label].unique()) > 1:
            axVars.append(texLabel)
            
    if len(axVars) == 1:
        return (r'%s' %axVars[0])
    else:
        return 'x(?)'

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