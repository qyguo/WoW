from ROOT import *
from array import array
import os
#from Utils import *

dirMC = {}

#dirMC['2018'] = '/eos/home-v/vmilosev/Skim_2018_HZZ/WoW/'
#dirMC['2017'] = '/eos/home-v/vmilosev/Skim_2018_HZZ/WoW/' # 2018 for now until we sort out the locations
#dirMC['2016'] = '/eos/home-v/vmilosev/Skim_2018_HZZ/WoW/' # 2018 for now until we sort out the locations
#dirMC['2018'] = '/eos/user/q/qguo/newNTuple_UL/2018/Slimmed_2p5/'
#dirMC['2018'] = '/publicfs/cms/data/hzz/guoqy/newNTuple_UL/2018/Slimmed_2p5_JES/'
dirMC['2018'] = '/publicfs/cms/data/hzz/guoqy/newNTuple_UL/2018/Slimmed_2p5_JES_3/'
dirMC['2018kgg'] = '/publicfs/cms/data/hzz/jtahir/UL2018/MC/Slimmed_2p5_JES_newkgg/'
dirMC['2018ggzz'] = '/publicfs/cms/data/hzz/guoqy/newNTuple_UL/2018/Slimmed_2p5/' # all ggZZ


#dirMC['2017'] = '/publicfs/cms/data/hzz/guoqy/newNTuple_UL/2017/Slimmed_2p5_JES/'
dirMC['2017'] = '/publicfs/cms/data/hzz/guoqy/newNTuple_UL/2017/Slimmed_2p5_JES_3/'
dirMC['2017kgg'] = '/publicfs/cms/data/hzz/jtahir/UL2017/MC/Slimmed_2p5_JES_newkgg/'
dirMC['2017ggzz'] = '/publicfs/cms/data/hzz/guoqy/newNTuple_UL/2017/Slimmed_2p5/' # all ggZZ
#dirMC['2016'] = '/eos/cms/store/group/phys_muon/TagAndProbe/HZZ4L/2016/UL/MC/postVFP/'
#dirMC['2016'] = '/eos/cms/store/group/phys_muon/TagAndProbe/HZZ4L/2016/UL/MC/Slimmed_2p5/'   # postVFP
#dirMC['2016'] = 'Slimmed_2p5/'   # postVFP
#dirMC['2016'] = 'Slimmed_2p5_JES/'   # postVFP
dirMC['2016'] = 'Slimmed_2p5_premj1j2fix/'   # postVFP
#dirMC['2016preVFP'] = 'Slimmed_2p5_premj1j2fix/'   # postVFP
dirMC['2016preVFP'] = '/publicfs/cms/data/hzz/jtahir/UL2016/MC/preVFP/Slimmed_2p5_JES/'   # preVFP
dirMC['2016preVFPkgg'] = '/publicfs/cms/data/hzz/jtahir/UL2016/MC/preVFP/Slimmed_2p5_JES_newkgg/'   # preVFP kgg
dirMC['2016preVFPggzz'] = '/publicfs/cms/data/hzz/jtahir/UL2016/MC/preVFP/Slimmed_2p5/'   # preVFP all ggZZ

dirMC['2016postVFP'] = '/publicfs/cms/data/hzz/jtahir/UL2016/MC/postVFP/Slimmed_2p5_JES/'   # postVFP
dirMC['2016postVFPkgg'] = '/publicfs/cms/data/hzz/jtahir/UL2016/MC/postVFP/Slimmed_2p5_JES_newkgg/'   # postVFP kgg
dirMC['2016postVFPggzz'] = '/publicfs/cms/data/hzz/jtahir/UL2016/MC/Slimmed_2p5/'   # postVFP all ggZZ

dirData_94 = '/eos/home-v/vmilosev/Skim_2018_HZZ/WoW/'  

#border_msg("samples directory: "+dirMC['2018'])

SamplesMC = {}

SamplesMC['2018'] = [
###
#'GluGluHToZZTo4L_M120_TuneCP5_13TeV_powheg2_JHUGenV7011_pythia8_slimmed_newMuSF_add2p5.root',
#'GluGluHToZZTo4L_M124_TuneCP5_13TeV_powheg2_JHUGenV7011_pythia8_slimmed_newMuSF_add2p5.root',
'GluGluHToZZTo4L_M125_TuneCP5_13TeV_powheg2_JHUGenV7011_pythia8_slimmed_newMuSF_add2p5.root',
#'GluGluHToZZTo4L_M126_TuneCP5_13TeV_powheg2_JHUGenV7011_pythia8_slimmed_newMuSF_add2p5.root',
#'GluGluHToZZTo4L_M130_TuneCP5_13TeV_powheg2_JHUGenV7011_pythia8_slimmed_newMuSF_add2p5.root',
#'VBF_HToZZTo4L_M120_TuneCP5_13TeV_powheg2_JHUGenV7011_pythia8_slimmed_newMuSF_add2p5.root',
#'VBF_HToZZTo4L_M124_TuneCP5_13TeV_powheg2_JHUGenV7011_pythia8_slimmed_newMuSF_add2p5.root',
'VBF_HToZZTo4L_M125_TuneCP5_13TeV_powheg2_JHUGenV7011_pythia8_slimmed_newMuSF_add2p5.root',
#'VBF_HToZZTo4L_M126_TuneCP5_13TeV_powheg2_JHUGenV7011_pythia8_slimmed_newMuSF_add2p5.root',
#'VBF_HToZZTo4L_M130_TuneCP5_13TeV_powheg2_JHUGenV7011_pythia8_slimmed_newMuSF_add2p5.root',
#'WH_HToZZTo4L_M120_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV7011_pythia8_slimmed_newMuSF_add2p5.root',
#'WH_HToZZTo4L_M124_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV7011_pythia8_slimmed_newMuSF_add2p5.root',
'WH_HToZZTo4L_M125_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV7011_pythia8_slimmed_newMuSF_add2p5.root',
#'WH_HToZZTo4L_M126_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV7011_pythia8_slimmed_newMuSF_add2p5.root',
#'WH_HToZZTo4L_M130_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV7011_pythia8_slimmed_newMuSF_add2p5.root',
#'ZH_HToZZ_4LFilter_M120_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV7011_pythia8_slimmed_newMuSF_add2p5.root',
#'ZH_HToZZ_4LFilter_M124_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV7011_pythia8_slimmed_newMuSF_add2p5.root',
'ZH_HToZZ_4LFilter_M125_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV7011_pythia8_slimmed_newMuSF_add2p5.root',
#'ZH_HToZZ_4LFilter_M126_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV7011_pythia8_slimmed_newMuSF_add2p5.root',
#'ZH_HToZZ_4LFilter_M130_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV7011_pythia8_slimmed_newMuSF_add2p5.root',
#'ttH_HToZZ_4LFilter_M120_TuneCP5_13TeV_powheg2_JHUGenV7011_pythia8_slimmed_newMuSF_add2p5.root',
#'ttH_HToZZ_4LFilter_M124_TuneCP5_13TeV_powheg2_JHUGenV7011_pythia8_slimmed_newMuSF_add2p5.root',
'ttH_HToZZ_4LFilter_M125_TuneCP5_13TeV_powheg2_JHUGenV7011_pythia8_slimmed_newMuSF_add2p5.root',
#'ttH_HToZZ_4LFilter_M126_TuneCP5_13TeV_powheg2_JHUGenV7011_pythia8_slimmed_newMuSF_add2p5.root',
#'ttH_HToZZ_4LFilter_M130_TuneCP5_13TeV_powheg2_JHUGenV7011_pythia8_slimmed_newMuSF_add2p5.root',
#'GluGluToContinToZZTo2e2mu_TuneCP5_13TeV-mcfm701-pythia8_slimmed_newMuSF_add2p5.root',
#'GluGluToContinToZZTo4e_TuneCP5_13TeV-mcfm701-pythia8_slimmed_newMuSF_add2p5.root',
#'GluGluToContinToZZTo4mu_TuneCP5_13TeV-mcfm701-pythia8_slimmed_newMuSF_add2p5.root',
#'GluGluToContinToZZTo2e2mu_TuneCP5_13TeV-mcfm701-pythia8_slimmed_newMuSF_add2p5.root',
#'GluGluToContinToZZTo4e_TuneCP5_13TeV-mcfm701-pythia8_slimmed_newMuSF_add2p5.root',
#'GluGluToContinToZZTo4mu_TuneCP5_13TeV-mcfm701-pythia8_slimmed_newMuSF_add2p5.root',
'GluGluToContinToZZTo4tau_TuneCP5_13TeV-mcfm701-pythia8_slimmed_newMuSF_add2p5.root',
'GluGluToContinToZZTo2e2tau_TuneCP5_13TeV-mcfm701-pythia8_slimmed_newMuSF_add2p5.root',
'GluGluToContinToZZTo2mu2tau_TuneCP5_13TeV-mcfm701-pythia8_slimmed_newMuSF_add2p5.root',
'GluGluToContinToZZTo2e2mu_TuneCP5_13TeV-mcfm701-pythia8_slimmed_newMuSF_add2p5.root',
'GluGluToContinToZZTo4e_TuneCP5_13TeV-mcfm701-pythia8_slimmed_newMuSF_add2p5.root',
'GluGluToContinToZZTo4mu_TuneCP5_13TeV-mcfm701-pythia8_slimmed_newMuSF_add2p5.root',
'ZZTo4L_TuneCP5_13TeV_powheg_pythia8_slimmed_newMuSF_add2p5.root'
]

SamplesMC['2017'] = [
#'GluGluHToZZTo4L_M120_TuneCP5_13TeV_powheg2_JHUGenV7011_pythia8_slimmed_newMuSF_add2p5.root',
#'GluGluHToZZTo4L_M124_TuneCP5_13TeV_powheg2_JHUGenV7011_pythia8_slimmed_newMuSF_add2p5.root',
'GluGluHToZZTo4L_M125_TuneCP5_13TeV_powheg2_JHUGenV7011_pythia8_slimmed_newMuSF_add2p5.root',
#'GluGluHToZZTo4L_M126_TuneCP5_13TeV_powheg2_JHUGenV7011_pythia8_slimmed_newMuSF_add2p5.root',
#'GluGluHToZZTo4L_M130_TuneCP5_13TeV_powheg2_JHUGenV7011_pythia8_slimmed_newMuSF_add2p5.root',
#'GluGluToZH_HToZZTo4L_M125_TuneCP5_13TeV-jhugenv723-pythia8_slimmed_newMuSF_add2p5.root',
#'VBF_HToZZTo4L_M120_TuneCP5_13TeV_powheg2_JHUGenV7011_pythia8_slimmed_newMuSF_add2p5.root',
#'VBF_HToZZTo4L_M124_TuneCP5_13TeV_powheg2_JHUGenV7011_pythia8_slimmed_newMuSF_add2p5.root',
'VBF_HToZZTo4L_M125_TuneCP5_13TeV_powheg2_JHUGenV7011_pythia8_slimmed_newMuSF_add2p5.root',
#'VBF_HToZZTo4L_M126_TuneCP5_13TeV_powheg2_JHUGenV7011_pythia8_slimmed_newMuSF_add2p5.root',
#'VBF_HToZZTo4L_M130_TuneCP5_13TeV_powheg2_JHUGenV7011_pythia8_slimmed_newMuSF_add2p5.root',
#'WH_HToZZTo4L_M120_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV7011_pythia8_slimmed_newMuSF_add2p5.root',
#'WH_HToZZTo4L_M124_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV7011_pythia8_slimmed_newMuSF_add2p5.root',
'WH_HToZZTo4L_M125_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV7011_pythia8_slimmed_newMuSF_add2p5.root',
#'WH_HToZZTo4L_M126_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV7011_pythia8_slimmed_newMuSF_add2p5.root',
#'WH_HToZZTo4L_M130_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV7011_pythia8_slimmed_newMuSF_add2p5.root',
##'WminusH_HToZZTo4L_M120_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV7011_pythia8_slimmed_newMuSF_add2p5.root',
##'WminusH_HToZZTo4L_M124_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV7011_pythia8_slimmed_newMuSF_add2p5.root',
##'WminusH_HToZZTo4L_M125_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV7011_pythia8_slimmed_newMuSF_add2p5.root',
##'WminusH_HToZZTo4L_M126_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV7011_pythia8_slimmed_newMuSF_add2p5.root',
##'WminusH_HToZZTo4L_M130_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV7011_pythia8_slimmed_newMuSF_add2p5.root',
##'WplusH_HToZZTo4L_M120_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV7011_pythia8_slimmed_newMuSF_add2p5.root',
##'WplusH_HToZZTo4L_M124_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV7011_pythia8_slimmed_newMuSF_add2p5.root',
##'WplusH_HToZZTo4L_M125_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV7011_pythia8_slimmed_newMuSF_add2p5.root',
##'WplusH_HToZZTo4L_M126_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV7011_pythia8_slimmed_newMuSF_add2p5.root',
##'WplusH_HToZZTo4L_M130_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV7011_pythia8_slimmed_newMuSF_add2p5.root',
#'ZH_HToZZ_4LFilter_M120_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV7011_pythia8_slimmed_newMuSF_add2p5.root',
#'ZH_HToZZ_4LFilter_M124_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV7011_pythia8_slimmed_newMuSF_add2p5.root',
'ZH_HToZZ_4LFilter_M125_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV7011_pythia8_slimmed_newMuSF_add2p5.root',
#'ZH_HToZZ_4LFilter_M126_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV7011_pythia8_slimmed_newMuSF_add2p5.root',
#'ZH_HToZZ_4LFilter_M130_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV7011_pythia8_slimmed_newMuSF_add2p5.root',
#'ttH_HToZZ_4LFilter_M120_TuneCP5_13TeV_powheg2_JHUGenV7011_pythia8_slimmed_newMuSF_add2p5.root',
#'ttH_HToZZ_4LFilter_M124_TuneCP5_13TeV_powheg2_JHUGenV7011_pythia8_slimmed_newMuSF_add2p5.root',
'ttH_HToZZ_4LFilter_M125_TuneCP5_13TeV_powheg2_JHUGenV7011_pythia8_slimmed_newMuSF_add2p5.root',
#'ttH_HToZZ_4LFilter_M126_TuneCP5_13TeV_powheg2_JHUGenV7011_pythia8_slimmed_newMuSF_add2p5.root',
#'ttH_HToZZ_4LFilter_M130_TuneCP5_13TeV_powheg2_JHUGenV7011_pythia8_slimmed_newMuSF_add2p5.root',
'GluGluToContinToZZTo4tau_TuneCP5_13TeV-mcfm701-pythia8_slimmed_newMuSF_add2p5.root',
'GluGluToContinToZZTo2e2tau_TuneCP5_13TeV-mcfm701-pythia8_slimmed_newMuSF_add2p5.root',
'GluGluToContinToZZTo2mu2tau_TuneCP5_13TeV-mcfm701-pythia8_slimmed_newMuSF_add2p5.root',
'GluGluToContinToZZTo2e2mu_TuneCP5_13TeV-mcfm701-pythia8_slimmed_newMuSF_add2p5.root',
'GluGluToContinToZZTo4e_TuneCP5_13TeV-mcfm701-pythia8_slimmed_newMuSF_add2p5.root',
'GluGluToContinToZZTo4mu_TuneCP5_13TeV-mcfm701-pythia8_slimmed_newMuSF_add2p5.root',
'ZZTo4L_TuneCP5_13TeV_powheg_pythia8_slimmed_newMuSF_add2p5.root'
]

# FIXME with actual file names for preVFP
SamplesMC['2016'] = [
'GluGluHToZZTo4L_M125_TuneCP5_13TeV_powheg2_JHUGenV7011_pythia8_slimmed_newMuSF_add2p5.root',
'VBF_HToZZTo4L_M125_TuneCP5_13TeV_powheg2_JHUGenV7011_pythia8_slimmed_newMuSF_add2p5.root',
'WH_HToZZTo4L_M125_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV7011_pythia8_slimmed_newMuSF_add2p5.root',
'ZH_HToZZ_4LFilter_M125_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV7011_pythia8_slimmed_newMuSF_add2p5.root',
'ttH_HToZZ_4LFilter_M125_TuneCP5_13TeV_powheg2_JHUGenV7011_pythia8_slimmed_newMuSF_add2p5.root',
'GluGluToContinToZZTo4tau_TuneCP5_13TeV-mcfm701-pythia8_slimmed_newMuSF_add2p5.root',
'GluGluToContinToZZTo2e2tau_TuneCP5_13TeV-mcfm701-pythia8_slimmed_newMuSF_add2p5.root',
'GluGluToContinToZZTo2mu2tau_TuneCP5_13TeV-mcfm701-pythia8_slimmed_newMuSF_add2p5.root',
'GluGluToContinToZZTo2e2mu_TuneCP5_13TeV-mcfm701-pythia8_slimmed_newMuSF_add2p5.root',
'GluGluToContinToZZTo4mu_TuneCP5_13TeV-mcfm701-pythia8_slimmed_newMuSF_add2p5.root',
'GluGluToContinToZZTo4e_TuneCP5_13TeV-mcfm701-pythia8_slimmed_newMuSF_add2p5.root',
'ZZTo4L_TuneCP5_13TeV_powheg_pythia8_slimmed_newMuSF_add2p5.root'
]

SamplesMC['2016preVFP'] = [
'GluGluHToZZTo4L_M125_TuneCP5_13TeV_powheg2_JHUGenV7011_pythia8_slimmed_newMuSF_add2p5.root',
'VBF_HToZZTo4L_M125_TuneCP5_13TeV_powheg2_JHUGenV7011_pythia8_slimmed_newMuSF_add2p5.root',
'WH_HToZZTo4L_M125_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV7011_pythia8_slimmed_newMuSF_add2p5.root',
'ZH_HToZZ_4LFilter_M125_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV7011_pythia8_slimmed_newMuSF_add2p5.root',
'ttH_HToZZ_4LFilter_M125_TuneCP5_13TeV_powheg2_JHUGenV7011_pythia8_slimmed_newMuSF_add2p5.root',
'GluGluToContinToZZTo4tau_TuneCP5_13TeV-mcfm701-pythia8_slimmed_newMuSF_add2p5.root',
'GluGluToContinToZZTo2e2tau_TuneCP5_13TeV-mcfm701-pythia8_slimmed_newMuSF_add2p5.root',
'GluGluToContinToZZTo2mu2tau_TuneCP5_13TeV-mcfm701-pythia8_slimmed_newMuSF_add2p5.root',
'GluGluToContinToZZTo2e2mu_TuneCP5_13TeV-mcfm701-pythia8_slimmed_newMuSF_add2p5.root',
'GluGluToContinToZZTo4mu_TuneCP5_13TeV-mcfm701-pythia8_slimmed_newMuSF_add2p5.root',
'GluGluToContinToZZTo4e_TuneCP5_13TeV-mcfm701-pythia8_slimmed_newMuSF_add2p5.root',
'ZZTo4L_TuneCP5_13TeV_powheg_pythia8_slimmed_newMuSF_add2p5.root'
]

SamplesMC['2016postVFP'] = [
'GluGluHToZZTo4L_M125_TuneCP5_13TeV_powheg2_JHUGenV7011_pythia8_slimmed_newMuSF_add2p5.root',
'VBF_HToZZTo4L_M125_TuneCP5_13TeV_powheg2_JHUGenV7011_pythia8_slimmed_newMuSF_add2p5.root',
'WH_HToZZTo4L_M125_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV7011_pythia8_slimmed_newMuSF_add2p5.root',
'ZH_HToZZ_4LFilter_M125_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV7011_pythia8_slimmed_newMuSF_add2p5.root',
'ttH_HToZZ_4LFilter_M125_TuneCP5_13TeV_powheg2_JHUGenV7011_pythia8_slimmed_newMuSF_add2p5.root',
'GluGluToContinToZZTo4tau_TuneCP5_13TeV-mcfm701-pythia8_slimmed_newMuSF_add2p5.root',
'GluGluToContinToZZTo2e2tau_TuneCP5_13TeV-mcfm701-pythia8_slimmed_newMuSF_add2p5.root',
'GluGluToContinToZZTo2mu2tau_TuneCP5_13TeV-mcfm701-pythia8_slimmed_newMuSF_add2p5.root',
'GluGluToContinToZZTo2e2mu_TuneCP5_13TeV-mcfm701-pythia8_slimmed_newMuSF_add2p5.root',
'GluGluToContinToZZTo4mu_TuneCP5_13TeV-mcfm701-pythia8_slimmed_newMuSF_add2p5.root',
'GluGluToContinToZZTo4e_TuneCP5_13TeV-mcfm701-pythia8_slimmed_newMuSF_add2p5.root',
'ZZTo4L_TuneCP5_13TeV_powheg_pythia8_slimmed_newMuSF_add2p5.root'
]



SamplesData = {}

#SamplesData['2016preVFP'] = ['data_UL2016_preVFP_noDuplicates_slimmed_add2p5_bkg.root']  # FIXME with actual file name
#SamplesData['2016postVFP'] = ['data_UL2016_postVFP_noDuplicates_slimmed_add2p5_bkg.root']
#SamplesData['2016'] = ['data_UL2016_noDuplicates_slimmed_add2p5_bkg.root']

#SamplesData['2016'] = ['data_UL2016_noDuplicates_slimmed_newMuSF_add2p5_bkg.root']
SamplesData['2016'] = ['data_UL2016_noDuplicates_slimmed_newMuSF_add2p5_sig.root']
SamplesData['2017'] = ['DataUL2017_all_noDuplicates_slimmed_newMuSF_add2p5.root']
SamplesData['2018'] = ['DataUL2018_all_noDuplicates_slimmed_newMuSF_add2p5.root']




def GrabMCTrees(era = '2018'):
    RootFile = {}
    Tree = {}
    nEvents = {}
    sumw = {}

    for i in range(0,len(SamplesMC[era])):

        sample = SamplesMC[era][i].rstrip('.root')


        if ("NNLOPS" in sample):
            RootFile[sample] = TFile(dirMC[era]+'/'+sample+'.root',"READ")
            Tree[sample]  = RootFile[sample].Get("Ana/passedEvents")

        elif ("ggH_amcatnloFXFX" in sample):
            RootFile[sample] = TFile(dirMC[era]+'/'+sample+'.root',"READ")
            Tree[sample]  = RootFile[sample].Get("Ana/passedEvents")

        else:
            RootFile[sample] = TFile.Open(dirMC[era]+'/'+sample+'.root',"READ")
            Tree[sample]  = RootFile[sample].Get("Ana/passedEvents")

        h_nevents = RootFile[sample].Get("Ana/nEvents")
        h_sumw = RootFile[sample].Get("Ana/sumWeights")

        if (h_nevents): nEvents[sample] = h_nevents.Integral()
        else: nEvents[sample] = 0.

        if (h_sumw): sumw[sample] = h_sumw.Integral()
        else: sumw[sample] = 0.

        if (not Tree[sample]): print sample+' has no passedEvents tree'
        else:
            print('{sample:37}\t nevents: {nevents:11}\t sumw: {sumw}'.format(
                sample = sample, nevents = nEvents[sample], sumw = sumw[sample]))

    return RootFile, Tree, nEvents, sumw
