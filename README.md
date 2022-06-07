### Add to the existing ntuples, JES uncertainty source branches
Recommendation to split the JES uncertainty into all sources can be seen at:
```
https://twiki.cern.ch/twiki/bin/viewauth/CMS/JECUncertaintySources#Run_2_reduced_set_of_uncertainty
```
The procedure is done by reconstructing the jet observables by adding the individual JES source's uncertainty provided. 
For this we need Jet modules from latest `CMSSW`, hence, the macro is to be run in `src` directory of the latest release used for the ntuple production. To install the release and compile the packages, follow the instructions here:
```
https://github.com/qyguo/UFHZZAnalysisRun2/tree/UL_10_6_26
```

Once done, get the slimmed macro from `JES_split` branch. 

```
git clone -b JES_split git@github.com:qyguo/wow.git 
```

Update the ntuple's path in `slimNtuple_JEC.C` and then run it (as in `test_run.sh`):

```
root -l slimNtuple_JEC.C+\(2016,\"ttH_HToZZ_4LFilter_M124_TuneCP5_13TeV_powheg2_JHUGenV7011_pythia8\",true,true,false,true\)
```

arg[0]: default year is 2016

arg[1]: default sample is "ttH_HToZZ_4LFilter_M124_TuneCP5_13TeV_powheg2_JHUGenV7011_pythia8"

arg[2]: default isMC is true

arg[3]: default isSignal is true 

arg[4]: default Test is false (test true will have less branch)

arg[5]: default true to apply JEC 


if not target to slim the Ntuple, just add new muon SFs. Please run newSFmuon.C:

```
root -l newSFmuon.C\(2017,\"bbH_HToZZTo4L_M125_TuneCP2_13TeV-jhugenv7011-pythia8\"\)
```

arg[0]: default year is 2017
arg[1]: default sample is "bbH_HToZZTo4L_M125_TuneCP2_13TeV-jhugenv7011-pythia8"
