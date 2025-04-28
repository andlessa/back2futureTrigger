# Slow LLP MC Generation Studies
This follows and adapts  the [MCGeneration](https://atlassoftwaredocs.web.cern.ch/analysis-software/AnalysisSWTutorial/mc_intro/) tutorials using a MadGraph/Pythia Card Provided by theorists Andre Lessa and Jose Zurita. 

## Get The Code

Set up ssh keys with gitlab if you haven't then clone recursively the repository.

```
git clone --recursive ssh://git@gitlab.cern.ch:7999/toheintz/edmgeneratorstudies.git
```

## Overview of workflow

The MC Production consists of three subsequent steps, namely MCGeneration, MCDerivation, and nTuple Dumping; thereafter, the trigger response is simulated and the trigger efficiencies are calculated. 
The procedure is parallelised using HTCondor and Condor DAGman with latter automatically processing the desired workflow (`batch_scripts/MCProduction<MG/PYTHIA>.dag`). 
After producing the MC samples, plots of the kinematics, and trigger efficiencies can be produced with a Jupyter Notebook (`batch_scripts/plot.py`)

Detailed comments on the different steps are provided below as well as suggestions which part of the code should/could be adjusted before running it.
**Explicit instructions how to actually run the code are given in the next section.**

#### Step 0: Preparation
The required file/directory structure is created with `parameters<MG/PYTHIA>.sh`, and the MG/Pythia cards are individually modified for each point in the parameter space.

Some "Hyper parameters" can be adjusted in `batch_scripts/parameters<MG/Pythia>.py`:
- `directory_index = 100000` should be a 6-digit number, and indicates the unique identifier of the first MC sample (if multiple samples with different parameters are produced, the ID is automatically counted upwards). If the MC Production is ran for a second time, the `directory_index` needs to be changed *manually* to a 6-digit number that is higher than the max ID of the previous production to avoid over-writing old samples.
- `nEvents = 10` indicates the number of events to be generated; should be at most 10000. 

Also, the desired parameter space can be specified in `parameters<MG/PYTHIA>.py`:
- `ms = [1000]` is a list of the scalar masses in GeV
- `betas = [0.01]` is a list of the $\beta^{\ast} = \sqrt{\frac{m_s^2 - 4m_1^2}{m_s^2}}$ values. (N.B. this variable is called `Dms1s = [0.01]` in parametersPYTHIA.py)
- `Dm10s =  [110]` is a list of the $\Delta m_{10} = m_{\chi_{1}} - m_{\chi_{0}}$ values in GeV.
- `lifetime_generated = [3.] #m` should contain only one value for the proper lifetime to be generated (one value shold be enough, samples with other lifetimes can be generated as toy MC samples according to `lifetime_toys = [...]`). Thoug, in Pythia multiple lifetimes can be explicitly generated in parallel by extending the `lifetime_generated` list.
- `lifetime_toys` is a list of lifetimes for which a toy MC approach should be used.

There are other parameters in the MadGraph+Pythia production:
- `Chxx10s = [1.5e-4]` should only contain one value that will be scaled automatically such that the width of the Chi1s sorresponds to `lifetime_generated`.
- `sinAlphas = [0.0001]` is a list of the $\sin \alpha$ values.
- The LambdaUV value is currently *fixed* to 5 TeV.


#### Step 1: MC Generation
##### a) MadGraph
The generation is made in `generationMG.sh` using the Athena Job Transform `Gen_tf.py`. 
Furthermore, we use the submodule `MadGraphModels` to link the Athena version of the backToTheFutureMinH model.

Instead of having the Chxx10 coupling as one of the parameters, it might be physically interesting to use the mean lifetime of the Chi1 (which is related to Chxx10) as parameter.
However, for off-shell Higgs there is no closed relation between lifetime, Chxx10, and the masses (three body decay).

Therefore, we first generate a small amount of events with the initial `Chxx10s` value as specified in Step 0; we use the result of the "auto width" calculation to scale the initial Chxx10 coupling to such a value that yields a Chi1 width corresponding to the `lifetime_generated` value.
We then use this scaled Chxx10 value in the second part of `generationMG.sh` to generate the desired number of events with the "aimed" proper lifetime. 

##### b) Pythia
The generation is made with `generationPYTHIA.sh`. The coupling Chxx10 is not a parameter, but Pythia does "force" the width according to the desired lifetime.

##### Toy MC approach
In general, it is sufficient to generate one lifetime value for each $m_s \times \beta \times \Delta m_{10}$ combination; we can use this sample and scale it in the trigger simulation (Step 4) to different lifetime values. 
(We decided to generate samples with proper lifetime of 3m as this is roughly where the eff. is at its maximum, and therefore scaling to the other lifetime points results in lowest uncertanties)

#### Step 2: MC Derivation
The TRUTH derivations are produced with `derivation.sh` using the Athena Job Transform `Derivation_tf.py`.

#### Step 3: nTuple dumping
Flat analysis-specific nTuples are produced using the submodule `l1-calo-llp-validation` (https://gitlab.cern.ch/toheintz/l1-calo-llp-validation).
The source code from `l1-calo-llp-validation` is compiled with `edmgeneratorstudies/batch_scripts/ntuplePre.sh`, and subsequently the eventloop `l1-calo-llp-validation/source/MCVal/share/MCVal_eljob.py` is run from `edmgeneratorstudies/batch_scripts/ntuple.sh`.

In this eventloop user-defined variables are written (see `l1-calo-llp-validation/source/MCVal/Root/MCValAlg.cxx`):
- Loop over truth particles: Scalar, BSM Chi1, BSM Chi0, leptons, quarks
- Loop over truth jets (all jet, and exclusively b jets)
- Higgs Boson: a) truth particles (on-shell higgs), b) reconstructed Higgs from opposite sign, same production vertex b quarks, and c) reclustered Higgs from b jets
- Delta Phi correction for out-of-Calorimeter decay (previously used, not sure if we should still use it later for the $\Delta \phi > \pi - 1$ cut)

#### Step 4: Trigger Simulation
The trigger simulation (`batch_scripts/trigger.py`) takes the nTuple with proper lifetime 3m as input. It simulatest the trigger decision, and returns a cutflow histogram.

- Both reclustered Higgs must be within $\eta < 3.2$ (acceptance of jFex)
- The transverse energy of the clustered Higgs Bosons (i.e. basically the sum of the ET of the jets) must meet the trigger threshold: for event N-1 $40 \text{GeV} < E_T < 100 \text{GeV}$, and for event N $40 \text{GeV} < E_T$.
- The clustered Higgs Bosons must be separated in $\phi$ by $\Delta \phi > \pi - 1$.  
- The LLP decay vertices must be within the acceptance of Tracker/ ECAL/ HCAL, i.e. $L_{xy} < 4 \text{m}$
- The "on-Time" decay should be within $0 < t < 10 \text{ns}$, whereas the "out-of-Time" decay takes place between $25 \text{ns} < t < 35 \text{ns}$.
 
This trigger simulation is not only performed for $c\tau = 3$ m, but is also repeated for different lifetimes $c\tau^\prime$. The cut values for $L_{xy}$ and $t$ are scaled by $1/R$ with $R = c\tau^\prime / 3 \text{m}$. 
Perhaps, a more sophisticated approach for the toy MCs would be as described here: https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/LifetimeReweighting 

#### Step 5: Calculation of Efficiencies
The trigger efficiencies are calculated in `batch_scripts/efficiency.py` 

1. `initialise_efficiencies`
First,  dictionaries for 1D ($m_s, \beta, \Delta_{10}, c\tau$) and 2D (6 combinations to combine 4 one-dim. variables) efficiencies are created, and the `root.TEfficiency` objects for all possible combinations are initialised.

2. `fill_efficiencies`
These efficiencies are filled binwise by looping over all samples, and calculating the ratio passed : total events based on the cutflow diagrams of the trigger simulation

3. `write_efficiencies`
The efficiencies are stored in .root files at `outputPath/efficiencies`.


#### Step 6: Post Processing (Plots)
A Jupyter Notebook `batch_scripts/plot.ipynb` is used to make plots of the efficiency curves and of some kinematics. 
(The notebook in the repo might me out-dated)

## Running Instructions

Before running the batch scripts for the first time, the local path of the repository needs to be added to each `<job>.sh` script. 
This can be done by typing
```
cd batch_scripts
source setupBATCH.sh
```
N.B. this script script only needs to be sourced before running the jobs for the very first time.

The data-heavy output files can be stored in a different directory which can be specified with
```
vim batch_scripts/outputPath.txt
...
```

The Condor jobs for the MG+Pythia production can be started with
```
cd batch_scripts
condor_submit_dag MCProductionMG.dag 
```
and the jobs for the pure Pythia generation with
```
cd batch_scripts
condor_submit_dag MCProductionPYTHIA.dag
```
MCProductionMG and MCProductionPYTHIA should not be ran in parallel, but subsequently after each other.

If the dag job has already been run sucessfully before and should be re-run (e.g. with other parameters), you must either delete the jobfiles which have been created automatically
```
rm MCProduction<MG/PYTHIA>.dag.*
```
or you can alternatively force DAGman to re-run the job with
```
condor_submit_dag -f MCProduction<MG/PYTHIA>.dag
```
In this case, the corresponding job files will be over-written.
The progress of the dag job can be checked with
```
condor_q -nobatch -dag
```

Plenty of job files will be created either directly after submitting the job or after sucess/failure of the job. In general, a file `MCProduction<MG/PYTHIA>.dag.metrics` will be created when the job has finished, and summarises whether the jobs were succesful or not. If `.rescue` files appear, something went wrong... after debugging, `condor_submit_dag MCProductionMG.dag`, will start the jobs automatically at the node which have failed (i.e. the successfull jobs are not ran again).

## Output Files
If the jobs run succesfully, several output files should be created (they can be found in the directory as specified in `batch_scripts/outputPath.txt`):
- EVNT files are stored in `outputPath/MCGeneration/evgen.xxxxxx.root` with xxxxxx being the 6-digit `directory_index`
- TRUTH Derivations in `outputPath/MCDerivation/xxxxxx/DAOD_TRUTH1.mcval.pool.xxxxxx.root`
- nTuples are stored in `outputPath/nTuples/xxxxxx/submitDirxxxxxx_mm_dd_hh_mm/data-ANALYSIS/xxxxxx.root` ("mm_dd_hh_mm" stands for 2-digits of month, day, hour, and minute): they contain an `analysis` tree with the respective branches as defined in `l1-calo-llp-validation/source/MCVal/Root/MCValAlg.cxx`
- The cutflow histograms of the trigger simulation can be found at `outputPath/nTuples/xxxxxx/submitDirxxxxxx_mm_dd_hh_mm/toys/cutflow_x.y.root` with `x.y` being the proper lifetime (in m) of the toy MC. Additionally, nTuples of the events surviving the trigger simulation are stored in `outputPath/nTuples/xxxxxx/submitDirxxxxxx_mm_dd_hh_mm/toys/trigger_x.y.root`.
- The 1D and 2D efficiency objects are stored in `outputPath/efficiencies/eff1D_var1.root` and `outputPath/efficiencies/eff2D_var1_var2.root` respectively. Within those root files, all possible combinations in the parameterspace are stored, the naming follows `var2_var3_var4` for 1D eff., and `var3_var4` for 2D eff (with `varx` being the value of the xth fixed variable). E.g. `2000.0_90.0_3.0` in `eff1D_beta.root` is the efficiency as function of beta for $m_s = 2000 \text{GeV}, $\Delta_{10} = 90 \text{GeV}, $c\tau = 3 \text{m}$.
- Plots are stored as specified in `batch_scripts/plot.ipynb`

If the jobs, however, failed, the log files might be useful:
- Log files for DAGman jobs are directly in `batch_scripts`
- Log/error/output files of the condor jobs are generally in the `edmgeneratorstudies` directory, e.g. `edmgeneratorstudies/nTuples`.
- Log files for the jobTransformers (e.g. `log.generate`, `log.Derivation` etc.) are in `outputPath/...`.

<!--
The workflow contains three main jobs (Node B-D) for the aforementioned three steps (Generation, Derivation, and Dumping), and a corresponding preprocessing script (Node A).
In the preprocessing script, the necessary files and directory structure will be created. Therefore, the parameters for the MCProduction should be defined here (before submitting to condorDag). This can currently be done inside the `batch_scripts/parameters.py` script. Also, the number of generated events can be modified in this script.
As mentioned earlier, three main jobs are carried out after the preprocessing:

1. MCGeneration.
The MCGeneration is currently carried out with AthGeneration 23.6.36 on Alma9. 
The `generation.sh` job basically generates the MCs for all parameter-sets with `Gen_tf.py`

2. MCDerivation.
The MCDerivation runs with Athena, 24.0.24 and produces Truth derivations by the `Derivation_tf.py` command.

3. NTuple Production.
The dumping of AODs to ntuples is implemented in `l1-calo-llp-validation`. After setting up Athena,25.2.2, the code is (if necessary) compiled and subsequently run. 

If necessary, the condor job options can be modified in the corresponding `.sub` files. This can be helpful if a certain node has an instable cvmfs connection. In this case, add the requirement 
```
Machine != "HOST_NAME"
```
Also, more CP memory or Disk space can be required in the job options if necessary.

Alternatively, when Condor Dugman is not working, the jobs can be submitted with `condor_submit <job>.sub` manually. However, successful completion of previous jobs is necessary for starting subsequent jobs.

## Post Processing Plots

The efficiency study and the plotting of some kinematic distributions is done with the jupyter notebook `/ultraSlowLLP/edmgeneratorstudies/eff.ipynb`. Further details are given in the notebook. 
Any modern pyROOT version should work. Personally I like to use miniconda to have my own local build but should also work on SWAN. 
Alternatively, a link to a Jupyter server can be created with  
```
source /cvmfs/sft.cern.ch/lcg/views/LCG_106/x86_64-el9-gcc13-opt/setup.sh
jupyter-notebook --no-browser --ip=* --port=8080
``` 
Thereafter, if working remotely create a local tunnel
```
ssh -L 8080:your-hep-desktop:8080 <hep-userid>@gw.hep.phy.cam.ac.uk
```
and open the link in VScode or a suitable Browser.


## ToDo

- [ ] Turn into a MC Request once simulation
- [x] Make $c\tau$ efficiency plots
- [ ] Add $\sigma$ baselines to MC estimates
--> 
