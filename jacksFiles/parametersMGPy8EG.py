import os
import sys
import shutil
import math
from math import sqrt
from jinja2 import Template
from decimal import Decimal

def formatStr(w):
    if w == "'auto'":
        return 0
    else:
        return format(Decimal(w.strip("'")) , 'f')

def parametersMG(filename="../run/parameters.txt" , splitJobs=False, writeJOs=False, updateCouplings=False , initialDir="", outputDir=""):
    with open(filename, 'w') as file:

        directory_index = 100000
        nEvents = 5000
        nJobs = 1
        #parameters in the re-parametrised parameterspace

        doChi2 = True

        mSs = [450, 600, 800, 1000] #mass of scalar particle [GeV]
        beta1s = [0.1, 0.2, 0.25, 0.3, 0.4] #boost of the chi1 in the labframe
        Dm21s = [5] #mass difference between chi2 and chi1 [GeV]
        Dm10s = [60, 100] #mass difference between chi1 and chi0 [GeV] (come back and try 100)
        ctau0s = [20] #proper lifetime of Chi1 [m]

        Chxx10s = [0.1] #Chxx10 parameter ~lifetime
        lambdaUVs = [2.2e7] #LambdaUV ~lifetime
        sinAlphas = [0.2] #mixing angle scalar/ SM higgs

        lifetime_toys = [5., 10., 15., 20., 30., 40.] #m

        with open("../run/lifetimes.txt", "w") as f:
            for value in lifetime_toys:
                f.write(f"{value}\n")

        mass_s_chi2_chi1_chi0 = []
        for mS in mSs:
            for beta1 in beta1s:
                for Dm21 in Dm21s:
                    if not doChi2:
                        Dm21 = 0
                    #for the MG card, we need the mass of the chi2 and the chi1, hence we have to undo the reparimetrisation
                    #and invert the parameter for the squared mass difference between ms and mChi1
                    #m1 = (0.25 * (mS ** 2. - beta1 ** 2. * mS ** 2.)) ** 0.5 #that is the special case for Dm21 = 0 only
                    m1_solutions = [ (-beta1 ** 2 * Dm21 ** 3 + sgn * (-beta1 ** 2 * Dm21 ** 4 * mS ** 2 + 2 * beta1 ** 2 * Dm21 ** 2 * mS ** 4 - beta1 ** 2 * mS ** 6 + Dm21 ** 4 * mS ** 2 - 2 * Dm21 ** 2 * mS ** 4 + mS ** 6) ** .5 + beta1 ** 2 * Dm21 * mS ** 2 + Dm21 ** 3 - Dm21 * mS ** 2)/(2 * (beta1 ** 2 * Dm21 ** 2 - Dm21 ** 2 + mS ** 2)) for sgn in [+1 , -1]]
                    m1 = max(m1_solutions)
                    m2 = m1 + Dm21
                    for Dm10 in Dm10s:
                        m0 = m1 - Dm10
                        mass_s_chi2_chi1_chi0.append([mS,m2,m1,m0])

        for mass_pair in mass_s_chi2_chi1_chi0:
            for ctau0 in ctau0s:
                for Chxx10 in Chxx10s:
                    for lambdaUV in lambdaUVs:
                        for sinAlpha in sinAlphas:
                            if (mass_pair[0] < 0 or mass_pair[1] < 0 or mass_pair[2] < 0 or mass_pair[3] < 0):
                                print(f"ERROR       MC sample with (ms = {mass_pair[0]}, m2 = {mass_pair[1]}, m1 = {mass_pair[2]}, m0 = {mass_pair[3]}) not physical due to negative particle masses - Skipping it.")
                                directory_index += 1
                                continue
                            if not (mass_pair[0] > mass_pair[1] and mass_pair[1] >= mass_pair[2] and mass_pair[2] > mass_pair[3] and mass_pair[3] >= 125/2): #the latter condition means chi0 > mHiggs / 2
                                print(f"ERROR       MC sample with (ms = {mass_pair[0]}, m2 = {mass_pair[1]}, m1 = {mass_pair[2]}, m0 = {mass_pair[3]}) not physical due to wrong mass hierarchy - Skipping it.")
                                directory_index += 1
                                continue
                            if (mass_pair[1] + mass_pair[2] > mass_pair[0]):
                                print(f"ERROR       MC sample with (ms = {mass_pair[0]}, m2 = {mass_pair[1]}, m1 = {mass_pair[2]}, m0 = {mass_pair[3]}) is kinematically not allowed - Skipping it.")
                                directory_index += 1
                                continue

                            #if not stated otherwise, use auto-width calculation
                            wS , w2, w1 , w0 = "'auto'" , "'auto'" , "'auto'" , "'auto'"
                            
                            dir=f"{outputDir}/MCGeneration/{directory_index}"
                            os.makedirs(dir , exist_ok=True)
                            
                            if updateCouplings != "auto":
                                w1 = 6.58 * 3.e-17 / ctau0 #1/GeV
                                w1_previous = w1 #just initialise w1_previous as w1, we will update it below
                                if doChi2:
                                    paramCard = f"{dir}/0/PROC_DipoleDM2_minimalH_UFO_0/Cards/param_card.dat" 
                                else:
                                    paramCard = f"{dir}/0/PROC_backToTheFutureMinH_0/Cards/param_card.dat" 
                                with open(paramCard, "r") as card:
                                    for line in card:
                                        if line.strip().startswith("DECAY 4000023"):
                                            parts = line.strip().split()
                                            if len(parts) >= 3:
                                                w1_previous = float(parts[2])
                                        elif line.strip().startswith("DECAY 4000022"):
                                            parts = line.strip().split()
                                            if len(parts) >= 3:
                                                w0 = float(parts[2])
                                        elif line.strip().startswith("DECAY 55"):
                                            parts = line.strip().split()
                                            if len(parts) >= 3:
                                                wS = float(parts[2])
                                        elif line.strip().startswith("DECAY 4000024"):
                                            parts = line.strip().split()
                                            if len(parts) >= 3:
                                                w2 = float(parts[2])

                                #width scales with square of coupling. We can scale the current Chxx10 based on the ratio sqrt(newWidth / currentWidth)
                                #Chxx10 = Chxx10 * sqrt(w1 / w1_previous)
                                lambdaUV = lambdaUV * sqrt(w1_previous / w1)
                                wS = f"'{wS}'"
                                w1 = f"'{w1}'"
                                if doChi2:
                                    w2 = f"'{w2}'"
                                else:
                                    w2 = w1 #assign something arbitrary to w2
                                w0 = f"'{w0}'"
                               
                            if not splitJobs:
                                file.write(f"{directory_index} {mass_pair[0]} {mass_pair[1]} {mass_pair[2]} {mass_pair[3]} {ctau0} {Chxx10} {lambdaUV} {sinAlpha} {formatStr(wS)} {formatStr(w2)} {formatStr(w1)} {formatStr(w0)} {nEvents} {nJobs} {doChi2}\n")
                            if splitJobs:
                                for idx in range(1,nJobs+1):
                                    file.write(f"{directory_index} {mass_pair[0]} {mass_pair[1]} {mass_pair[2]} {mass_pair[3]} {ctau0} {Chxx10} {lambdaUV} {sinAlpha} {formatStr(wS)} {formatStr(w2)} {formatStr(w1)} {formatStr(w0)} {nEvents} {nJobs} {idx} {doChi2}\n")
                            if writeJOs:
                                #modivy JO Template accordingly
                                params = {
                                    "mS": mass_pair[0],
                                    "m2": mass_pair[1],
                                    "m1": mass_pair[2],
                                    "m0": mass_pair[3],
                                    "wS": "'auto'",
                                    "w2": "'auto'",
                                    "w1": "'auto'",
                                    "w0": "'auto'",
                                    "Chxx10": Chxx10,
                                    "lambdaUV": lambdaUV,
                                    "cTau": ctau0,
                                    "sinAlpha": sinAlpha
                                }
                                # Load template
                                if doChi2:
                                    template_path = f"{initialDir}/MGPy8EG_JOs_template.py.j2"
                                else:
                                    template_path = f"{initialDir}/MGPy8EG_JOs_template_chi1chi1.py.j2"
                                with open(template_path, "r") as f: #FIXME differences between chi1 chi1 and chi1 chi2 model
                                    template_content = f.read()
                                # Render with variables
                                template = Template(template_content)
                                rendered_script = template.render(params)
                                # Save rendered script
                                output_path=f"{dir}/mc.MGPy8EG_A14NNPDF23_slowLLP.py"
                                with open(output_path, "w") as f:
                                    f.write(rendered_script)

                            directory_index += 1

initialPath = sys.argv[1]
updateCouplings = sys.argv[2]

with open('../run/outputPath.txt', 'r') as file:
    outputDir = file.read().strip()

if not os.path.isdir(outputDir):
    os.makedirs(outputDir, exist_ok=True)

parametersMG(filename="../run/parameters.txt",splitJobs=False,writeJOs=False,updateCouplings=updateCouplings,initialDir=initialPath,outputDir=outputDir)
parametersMG(filename="../run/parameters_split.txt",splitJobs=True,writeJOs=False,updateCouplings=updateCouplings,initialDir=initialPath,outputDir=outputDir)

parametersMG(filename=f"{outputDir}/parameters.txt",splitJobs=False, writeJOs=True,updateCouplings=updateCouplings,initialDir=initialPath,outputDir=outputDir)
parametersMG(filename=f"{outputDir}/parameters_split.txt",splitJobs=True,writeJOs=False,updateCouplings=updateCouplings,initialDir=initialPath,outputDir=outputDir)
