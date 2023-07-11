#!/usr/bin/python
import optparse
import os, pwd, commands
import sys
#grootargs = []
#sys.argv = grootargs
import ROOT, re, string
from ROOT import *
from scipy import interpolate
import matplotlib.pyplot as plt

import math
from math import *

from tdrStyle import *
setTDRStyle()

#from LoadData import *
#from LoadData_yield import *
from LoadData_yield_tau import *
from zx_contribution import *
#from Utils import *


def parseOptions():

    global opt, args

    usage = ('usage: %prog [options]\n'
             + '%prog -h for help')
    parser = optparse.OptionParser(usage)

    parser.add_option('',   '--obsName',  dest='OBSNAME',  type='string',default='mass4l',   help='Name of the observable, supported: "inclusive", "pT4l", "eta4l", "massZ2", "nJets"')
    #parser.add_option('',   '--obsName',  dest='OBSNAME',  type='string',default='massZ1',   help='Name of the observable, supported: "inclusive", "pT4l", "eta4l", "massZ2", "nJets"')
    parser.add_option('',   '--obsBins',  dest='OBSBINS',  type='string',default='|105.0|160.0|',   help='Bin boundaries for the diff. measurement separated by "|", e.g. as "|0|50|100|", use the defalut if empty string')
    #parser.add_option('',   '--obsBins',  dest='OBSBINS',  type='string',default='|70.0|1000.0|',   help='Bin boundaries for the diff. measurement separated by "|", e.g. as "|0|50|100|", use the defalut if empty string')
    #parser.add_option('',   '--obsBins',  dest='OBSBINS',  type='string',default='|40.0|130.0|',   help='Bin boundaries for the diff. measurement separated by "|", e.g. as "|0|50|100|", use the defalut if empty string')
    #parser.add_option('',   '--era',  dest='ERA',  type='string',default='2018',   help='Era to analyze, e.g. 2016, 2017, 2018 or Full ')
    #parser.add_option('',   '--era',  dest='ERA',  type='string',default='2016',   help='Era to analyze, e.g. 2016, 2017, 2018 or Full ')
    #parser.add_option('',   '--era',  dest='ERA',  type='string',default='2017',   help='Era to analyze, e.g. 2016, 2017, 2018 or Full ')
    #parser.add_option('',   '--era',  dest='ERA',  type='string',default='Full',   help='Era to analyze, e.g. 2016, 2017, 2018 or Full ')
    parser.add_option('',   '--era',  dest='ERA',  type='string',default='',   help='Era to analyze, e.g. 2016, 2017, 2018 or Full ')
    #parser.add_option('', '--range',      dest='RNG',type='string',default='low', help='run with the type of mass window e.g. signal, full, low, high etc. ')
    #parser.add_option('', '--range',      dest='RNG',type='string',default='signal', help='run with the type of mass window e.g. signal, full, low, high etc. ')
    #parser.add_option('', '--range',      dest='RNG',type='string',default='high', help='run with the type of mass window e.g. signal, full, low, high etc. ')
    #parser.add_option('', '--range',      dest='RNG',type='string',default='full', help='run with the type of mass window e.g. signal, full, low, high etc. ')
    parser.add_option('', '--range',      dest='RNG',type='string',default='', help='run with the type of mass window e.g. signal, full, low, high etc. ')
    parser.add_option('',   '--debug',  dest='DEBUG',  type='int',default=0,   help='0 if debug false, else debug True')
#    parser.add_option("-l",action="callback",callback=callback_rootargs)
#    parser.add_option("-q",action="callback",callback=callback_rootargs)
#    parser.add_option("-b",action="callback",callback=callback_rootargs)

    global opt, args
    (opt, args) = parser.parse_args()

#grootargs = []
#def callback_rootargs(option, opt, value, parser):
#    grootargs.append(opt)

def extract_yields(nbins, obsName, obs_bins, channel, year, DEBUG = 0): 
#xlabel, xunits, prelim, setLogX, setLogY,EventCat, mh, DEBUG = 0):
#def plot_m4l(channel, var, mh, bin, low, high, m4llow, m4lhigh, xlabel, xunits, prelim, setLogX, setLogY,EventCat):


    x_points = [125]

    #procs=['trueH_XH','trueH','bkg_qqzz','obs_data','bkg_ggzz']
    #procs=['trueH_ggH','trueH_XH','trueH','bkg_qqzz','obs_data','bkg_ggzz']

    #procs=['bkg_ggzz','trueH_ggH','trueH_XH','trueH','bkg_qqzz','obs_data'] # OK
    procs=['bkg_ggzz','trueH_ggH','trueH_VBF','trueH_WH','trueH_ZH','trueH_ttH','trueH_XH','trueH','bkg_qqzz','obs_data'] # OK
    #procs=['trueH_VBF']
    #procs=['bkg_qqzz','bkg_ggzz']
    #procs=['bkg_qqzz','bkg_ggzz']
    #procs = ['obs_data']
#    procs = ['bkg_ggzz']
    #procs = ['trueH_ggH']
    #procs = ['trueH_XH']
    #procs = ['trueH_ggH','trueH','trueH_XH' ]


    #proc_selections = {"trueH":"passedFullSelection==1 ","out_trueH":"(passedFullSelection==1 && passedFiducialSelection!=1)","fakeH":"isH4l!=1","bkg_qqzz":"passedFullSelection==1","bkg_ggzz":"passedFullSelection==1","obs_data":"(passedZXCRSelection==1)"}
    #proc_selections = {"trueH":"passedFullSelection==1 ","out_trueH":"(passedFullSelection==1 && passedFiducialSelection!=1)","fakeH":"isH4l!=1","bkg_qqzz":"passedFullSelection==1","bkg_ggzz":"passedFullSelection==1","obs_data":"passedFullSelection==1"}
    #proc_selections = {"trueH":"passedFullSelection==1","trueH_XH":"passedFullSelection==1","out_trueH":"(passedFullSelection==1 && passedFiducialSelection!=1)","fakeH":"isH4l!=1","bkg_qqzz":"passedFullSelection==1","bkg_ggzz":"passedFullSelection==1","obs_data":"passedFullSelection==1"}
    #proc_selections = {"trueH":"passedFullSelection==1","trueH_XH":"passedFullSelection==1","trueH_ggH":"passedFullSelection==1","out_trueH":"(passedFullSelection==1 && passedFiducialSelection!=1)","fakeH":"isH4l!=1","bkg_qqzz":"passedFullSelection==1","bkg_ggzz":"passedFullSelection==1","obs_data":"passedFullSelection==1"}
    proc_selections = {"trueH":"passedFullSelection==1","trueH_XH":"passedFullSelection==1","trueH_ggH":"passedFullSelection==1","trueH_VBF":"passedFullSelection==1","trueH_WH":"passedFullSelection==1","trueH_ZH":"passedFullSelection==1","trueH_ttH":"passedFullSelection==1","out_trueH":"(passedFullSelection==1 && passedFiducialSelection!=1)","fakeH":"isH4l!=1","bkg_qqzz":"passedFullSelection==1","bkg_ggzz":"passedFullSelection==1","obs_data":"passedFullSelection==1"}
    #proc_selections = {"trueH":"passedFullSelection==1 ","out_trueH":"(passedFullSelection==1 && passedFiducialSelection!=1)","fakeH":"isH4l!=1","bkg_qqzz":"passedFullSelection==1","bkg_ggzz":"passedFullSelection==1","obs_data":"passedZXCRSelection==1"}


    cut_m4l_reco = {"2e2mu":"(mass2e2mu>0)","4mu":"(mass4mu>0)","4e":"(mass4e>0)","4l":"mass4l>0"}
    #cut_m4l_reco = {"2e2mu":"(mass2e2mu>105&&mass2e2mu<160)","4mu":"(mass4mu>105&&mass4mu<160)","4e":"(mass4e>105&&mass4e<160)","4l":"mass4l>105&&mass4l<160"}


    #lumi = {"2016":35867.0,"2017":41370.0,"2018":58800.0}
    #lumi = {"2016":35900.0,"2017":41500.0,"2018":59700.0}
    #lumi = {"2016":35900.0,"2017":41500.0,"2018":58970.48}  # old numbers
    lumi = {"2016":36330.0,"2017":41480.0,"2018":59830.0}  # https://twiki.cern.ch/twiki/bin/viewauth/CMS/PdmVRun2LegacyAnalysis
    lumi_preVFP=19520; lumi_postVFP=16810;

#    weights = {"sig":"genWeight*pileupWeight*prefiringWeight*dataMCWeight_new*crossSection","bkg_qqzz":"(genWeight*pileupWeight*prefiringWeight*dataMCWeight_new*crossSection)*(k_qqZZ_qcd_M*k_qqZZ_ewk)","bkg_ggzz":"(genWeight*pileupWeight*prefiringWeight*dataMCWeight_new*crossSection)*(k_ggZZ)","obs_data":"1.0"} #  orig
    weights = {"sig":"genWeight*pileupWeight*prefiringWeight*dataMCWeight_new*crossSection","bkg_qqzz":"(genWeight*pileupWeight*prefiringWeight*dataMCWeight_new*crossSection)*(k_qqZZ_qcd_M*k_qqZZ_ewk)","bkg_ggzz":"(genWeight*pileupWeight*prefiringWeight*dataMCWeight_new*crossSection)*(k_ggZZ_new)","obs_data":"1.0"} # new k_ggZZ from CJLST
    #weights = {"sig":"genWeight*pileupWeight*prefiringWeight*dataMCWeight_new*crossSection","bkg_qqzz":"(genWeight*pileupWeight*prefiringWeight*dataMCWeight_new*crossSection)*(k_qqZZ_qcd_M*k_qqZZ_ewk)","bkg_ggzz":"(genWeight*pileupWeight*prefiringWeight*dataMCWeight_new*crossSection)*(k_ggZZ_reco_new)","obs_data":"1.0"} # new k_ggZZ from CJLST, extracted from reco mass4l (not from GENmass4l)

    #weights = {"sig":"genWeight*pileupWeight*dataMCWeight_new*crossSection","bkg_qqzz":"(genWeight*pileupWeight*dataMCWeight_new*crossSection)*(k_qqZZ_qcd_M*k_qqZZ_ewk)","bkg_ggzz":"(genWeight*pileupWeight*dataMCWeight_new*crossSection)*(k_ggZZ)","obs_data":"1.0"}


    added_year = TH1D('a_year', 'a_year',m4l_bins,m4l_low,m4l_high)
    red_bkg_year = TH1D('red_bkg_year', 'red_bkg_year',m4l_bins,m4l_low,m4l_high)
    red_bkg_pre_year = TH1D('red_bkg_pre_year', 'red_bkg_pre_year',m4l_bins,m4l_low,m4l_high)
    red_bkg_post_year = TH1D('red_bkg_post_year', 'red_bkg_post_year',m4l_bins,m4l_low,m4l_high)
    irr_bkg_year = TH1D('irr_bkg_year', 'irr_bkg_year',m4l_bins,m4l_low,m4l_high)
    rare_bkg_year = TH1D('rare_bkg_year', 'rare_bkg_year',m4l_bins,m4l_low,m4l_high)
    qqzz_year = TH1D('qqzz_year', 'qqzz_year',m4l_bins,m4l_low,m4l_high)
    ggzz_year = TH1D('ggzz_year', 'ggzz_year',m4l_bins,m4l_low,m4l_high)
    sig_year = TH1D('sig_year', 'sig_year',m4l_bins,m4l_low,m4l_high)
    sig_XH_year = TH1D('sig_XH_year', 'sig_XH_year',m4l_bins,m4l_low,m4l_high)
    sig_ggH_year = TH1D('sig_ggH_year', 'sig_ggH_year',m4l_bins,m4l_low,m4l_high)
    data_year = TH1D('data_year', 'data_year',m4l_bins,m4l_low,m4l_high)

    sig_VBF_year = TH1D('sig_VBF_year', 'sig_VBF_year',m4l_bins,m4l_low,m4l_high)
    sig_WH_year = TH1D('sig_WH_year', 'sig_WH_year',m4l_bins,m4l_low,m4l_high)
    sig_ZH_year = TH1D('sig_ZH_year', 'sig_ZH_year',m4l_bins,m4l_low,m4l_high)
    sig_ttH_year = TH1D('sig_ttH_year', 'sig_ttH_year',m4l_bins,m4l_low,m4l_high)

    Histos = {}
    nom = {}

    RootFile = {}
    RootFile_pre = {}
    RootFile_post = {}
 
 
    if (year=='2018' or year=='2017'):
        RootFile, Tree, nEvents, sumw = GrabMCTrees(year)
        print "Tree", Tree
    else:
        print "getting trees from pre "
        RootFile_pre, Tree_pre, nEvents_pre, sumw_pre = GrabMCTrees('2016preVFP')
        print "getting trees from post "
        RootFile_post, Tree_post, nEvents_post, sumw_post = GrabMCTrees('2016postVFP')


    for x_point in x_points:
            for proc in procs:

                RootFile_list = []
                RootFile_list_pre = []
                RootFile_list_post = []

		#if ((proc=='trueH' or proc=='out_trueH' or proc=='fakeH')): 
		if ((proc=='trueH' or proc=='trueH_XH' or proc=='trueH_ggH' or proc=='out_trueH' or proc=='fakeH' or proc=='trueH_VBF' or proc=='trueH_WH' or proc=='trueH_ZH' or proc=='trueH_ttH')): 
		    recoweight = weights["sig"]
		else: 
		    recoweight = weights[proc]

		print "recoweight:    ",recoweight
		if (proc=='obs_data'): samples = SamplesData[year]
		else: samples =SamplesMC[year]


		sumw2 = []; sumw2_pre = []; sumw2_post = [];
                for i in range(0,len(samples)):
                    sample = samples[i].rstrip('.root')

                    if((proc!='obs_data' and proc!='bkg_ggzz') and (not sample.startswith('ZZTo4L') and (not str(x_point) in sample))): continue
                    if((proc=='bkg_qqzz') and (not sample.startswith('ZZTo4L'))): continue # and (not str(x_point) in sample))): continue
                    if((proc=='bkg_ggzz' and proc!='obs_data' and proc!='bkg_qqzz' and (proc!='trueH' or proc!='trueH_XH' or proc!='trueH_ggH' or proc!='out_trueH' or proc!='fakeH')) and (not sample.startswith('GluGluToContin'))): continue # and (not str(x_point) in sample))): continue
                    #if((proc=='bkg_ggzz' and proc!='obs_data' and proc!='bkg_qqzz' and (proc!='trueH' or proc!='out_trueH' or proc!='fakeH')) and (not channel in sample)): continue # and (not str(x_point) in sample))): continue
                    #if((proc=='bkg_ggzz' and proc!='obs_data' and proc!='bkg_qqzz' and (proc!='trueH' or proc!='trueH_XH' or proc!='trueH_ggH' or proc!='out_trueH' or proc!='fakeH')) and (channel!='4l' and not channel in sample)): continue # and (not str(x_point) in sample))): continue
                    if((proc=='bkg_ggzz' and proc!='obs_data' and proc!='bkg_qqzz' and (proc!='trueH' or proc!='trueH_XH' or proc!='trueH_ggH' or proc!='out_trueH' or proc!='fakeH') and (not 'tau' in sample)) and (channel!='4l' and not channel in sample)): continue # including the tau final states


                    if ((proc!='obs_data' and proc!='bkg_qqzz' and proc!='bkg_ggzz') and (proc=='trueH' or proc=='trueH_XH' or proc=='out_trueH' or proc=='fakeH') and ((sample.startswith('ZZTo4L') or sample.startswith('GluGluToContin')))): continue
                    if ((proc!='obs_data' and proc!='bkg_qqzz' and proc!='bkg_ggzz') and (proc=='trueH_XH') and ((sample.startswith('GluGluHToZZTo4L_M125') or (sample.startswith('ZZTo4L')) or sample.startswith('GluGluToContin')))): continue
                    if ((proc!='obs_data' and proc!='bkg_qqzz' and proc!='bkg_ggzz') and (proc=='trueH_ggH') and ((not sample.startswith('GluGluHToZZTo4L_M125') or (sample.startswith('ZZTo4L')) or sample.startswith('GluGluToContin')))): continue
                    if ((proc!='obs_data' and proc!='bkg_qqzz' and proc!='bkg_ggzz') and (proc=='trueH_VBF') and ((not sample.startswith('VBF') or (sample.startswith('ZZTo4L')) or sample.startswith('GluGluToContin')))): continue
                    if ((proc!='obs_data' and proc!='bkg_qqzz' and proc!='bkg_ggzz') and (proc=='trueH_WH') and ((not sample.startswith('WH') or (sample.startswith('ZZTo4L')) or sample.startswith('GluGluToContin')))): continue
                    if ((proc!='obs_data' and proc!='bkg_qqzz' and proc!='bkg_ggzz') and (proc=='trueH_ZH') and ((not sample.startswith('ZH') or (sample.startswith('ZZTo4L')) or sample.startswith('GluGluToContin')))): continue
                    if ((proc!='obs_data' and proc!='bkg_qqzz' and proc!='bkg_ggzz') and (proc=='trueH_ttH') and ((not sample.startswith('ttH') or (sample.startswith('ZZTo4L')) or sample.startswith('GluGluToContin')))): continue

		    if (year=='2018' or year=='2017'):
			if (proc=='bkg_ggzz'):
#                        if (proc=='bkg_ggzz_new'): # FIXME , temporary
                            RootFile[sample] = TFile(dirMC[year+'kgg']+'/'+sample+'.root',"READ")  # currently, all files in the same directory
            	            RootFile_list.append(dirMC[year+'kgg']+'/'+sample+'.root')
                            #RootFile[sample] = TFile(dirMC[year+'ggzz']+'/'+sample+'.root',"READ")  # currently, all files in the same directory
            	            #RootFile_list.append(dirMC[year+'ggzz']+'/'+sample+'.root')
		        else:
                            RootFile[sample] = TFile(dirMC[year]+'/'+sample+'.root',"READ")  # currently, all files in the same directory
                            RootFile_list.append(dirMC[year]+'/'+sample+'.root')
		    elif (year=='2016'):
                        if (proc=='bkg_ggzz'):
			    #print sample
		            #print "dirMC['2016preVFPggzz']: ", dirMC['2016preVFPggzz']	
			    #print "dirMC['2016postVFPggzz']: ", dirMC['2016postVFPggzz']
                        #if (proc=='bkg_ggzz_new'): # FIXME , temporary
    		            RootFile_pre[sample] = TFile(dirMC['2016preVFPkgg']+'/'+sample+'.root',"READ")
    			    RootFile_list_pre.append(dirMC['2016preVFPkgg']+'/'+sample+'.root')   
    		            RootFile_post[sample] = TFile(dirMC['2016postVFPkgg']+'/'+sample+'.root',"READ")
    			    RootFile_list_post.append(dirMC['2016postVFPkgg']+'/'+sample+'.root')   
#    		            RootFile_pre[sample] = TFile(dirMC['2016preVFPggzz']+'/'+sample+'.root',"READ")
#    			    RootFile_list_pre.append(dirMC['2016preVFPggzz']+'/'+sample+'.root')   
#    		            RootFile_post[sample] = TFile(dirMC['2016postVFPggzz']+'/'+sample+'.root',"READ")
#    			    RootFile_list_post.append(dirMC['2016postVFPggzz']+'/'+sample+'.root')   
                        else:
                            RootFile_pre[sample] = TFile(dirMC['2016preVFP']+'/'+sample+'.root',"READ")
                            RootFile_list_pre.append(dirMC['2016preVFP']+'/'+sample+'.root')
                            RootFile_post[sample] = TFile(dirMC['2016postVFP']+'/'+sample+'.root',"READ")
                            RootFile_list_post.append(dirMC['2016postVFP']+'/'+sample+'.root')

		    print "sample is :   ", sample
		    if (proc!='obs_data' and (year=='2018' or year=='2017')): 
			sumw2.append(sumw[sample])  # from keys to list, for onward use
		    elif (proc!='obs_data' and (year=='2016')): 
			sumw2_pre.append(sumw_pre[sample]); sumw2_post.append(sumw_post[sample])
		if (year=='2018' or year=='2017'):
                    files = [ROOT.TFile(i) for i in RootFile_list]
		    print "proc:  ", proc, "  channel:   ", channel
		    print "files are:  ", files
		elif (year=='2016'):
		    files_pre = [ROOT.TFile(i) for i in RootFile_list_pre]
		    files_post = [ROOT.TFile(i) for i in RootFile_list_post]
		    print "proc:  ", proc, "  channel:   ", channel
		    print "files_pre are:  ", files_pre
		    print "files_post are:  ", files_post

		if (year=='2018' or year=='2017'):
		    if(proc=='obs_data'): trees = [i.Get("passedEvents") for i in files]
                    else: trees = [i.Get("Ana/passedEvents") for i in files]
		    print "trees are:  ", trees
		elif (year=='2016'):
		    if(proc=='obs_data'): trees_pre = [i.Get("passedEvents") for i in files_pre]; trees_post = [j.Get("passedEvents") for j in files_post];
                    else: trees_pre = [i.Get("Ana/passedEvents") for i in files_pre]; trees_post = [j.Get("Ana/passedEvents") for j in files_post];

            
		for recobin in range(len(obs_bins)-1):
		    #obs_reco_low = obs_bins[recobin]
    		    #obs_reco_high = obs_bins[recobin+1]
		    obs_reco_low = m4l_low; #obs_bins[recobin]
    		    obs_reco_high = m4l_high; # obs_bins[recobin+1]

		    print "obs_reco_low", obs_reco_low, "     obs_reco_high", obs_reco_high
		    cutobs_reco = "("+obs_reco+">="+str(obs_reco_low)+" && "+obs_reco+"<"+str(obs_reco_high)+")"
		    cut_nom = "("+recoweight+")*("+cutobs_reco+" && "+proc_selections[proc]+" && "+cut_m4l_reco[channel]+")"

		    if (not (proc=='trueH' or proc=='trueH_XH' or proc=='out_trueH' or proc=='fakeH')): x_point=125 # temporary fix for bkg
		    #processBin_nom = proc+'_'+str(x_point)+'_'+obsName+'_'+channel+'_recobin'+str(recobin)
		    processBin_nom = proc+'_'+str(x_point)+'_'+obsName+'_'+channel+'_recobin'+str(recobin)+year  # to be individual to an year
		    nom[processBin_nom] = 0.0; 

					    
		    if (year=='2018' or year=='2017'):
    		        Histos[processBin_nom] = TH1D(processBin_nom,processBin_nom, m4l_bins, m4l_low, m4l_high)
    		        #Histos[processBin_nom+'kgg'] = TH1D(processBin_nom,processBin_nom, m4l_bins, 0, 2.6)
			#test_cut_kgg = "k_ggZZ!=0"
		        #Histos[processBin_nom].Sumw2()
		    elif (year=='2016'):
                        Histos[processBin_nom] = TH1D(processBin_nom,processBin_nom, m4l_bins, m4l_low, m4l_high)
                        Histos[processBin_nom].Sumw2()

			Histos[processBin_nom+'pre'] = TH1D(processBin_nom+'pre',processBin_nom+'pre', m4l_bins, m4l_low, m4l_high)
			Histos[processBin_nom+'pre'].Sumw2()
			Histos[processBin_nom+'post'] = TH1D(processBin_nom+'post',processBin_nom+'post', m4l_bins, m4l_low, m4l_high)
			Histos[processBin_nom+'post'].Sumw2()
                        Histos[processBin_nom+'comb'] = TH1D(processBin_nom+'comb',processBin_nom+'comb', m4l_bins, m4l_low, m4l_high)
                        Histos[processBin_nom+'comb'].Sumw2()


 		    print "cut_nom :   ", cut_nom 
#			print "cut_nuis_up :   ", cut_nuis_up
#			print "cut_nuis_dn :   ", cut_nuis_dn 
#
		    yield_nom = 0.;
    		    i=0;  # iterator to read sumw2 from corresponding root file, meant to work for signal proc. only
    
#		    red_bkg = TH1D('red_bkg', 'red_bkg',m4l_bins,m4l_low,m4l_high)
    		    if (year=='2018' or year=='2017'):
    		        for tree in trees:
    		            if(proc=='obs_data'):
                                #tree.Draw("mass4l >> "+processBin_nom,"("+cut_nom+")","goff")
                                tree.Draw("mass4l >> "+processBin_nom,"("+cut_nom+")","goff")
				data_year.Add(Histos[processBin_nom])
    		            else:
    		                print "proc is:   ", proc
    		                print "channel is:   ", channel
    		                print "i:  ",i
#    			        print "len(sumw2):    ", len(sumw2)
#    		                print "sumw2[",i,"] : ", sumw2[i] 
    		                #tree.Draw("mass4l >> "+processBin_nom,"("+cut_nom+"*("+str(lumi[year])+")/"+str(sumw2[i])+")","goff")
    		                tree.Draw("mass4l >> "+processBin_nom,"("+cut_nom+"*("+str(lumi[year])+")/"+str(sumw2[i])+")","goff")
## try
#				if (year=='2018' and proc=='bkg_ggzz'):
#				    with open("ggZZ_info_2018.txt","a") as file:
#					file.write("Run:    "+str(Run)+" \n")  # FIXME

				
#				tree.Draw("mass4l >> "+processBin_nom+'kgg',test_cut_kgg,"goff")
#				if(channel=='4mu'): Histos[processBin_nom+'kgg'].SaveAs('kgg_'+channel+'_'+year+'.root') 
				if(proc=='trueH'): sig_year.Add(Histos[processBin_nom]) 
				if(proc=='trueH_XH'): sig_XH_year.Add(Histos[processBin_nom]) 
				if(proc=='trueH_ggH'): sig_ggH_year.Add(Histos[processBin_nom]) 
				if(proc=='trueH_VBF'): sig_VBF_year.Add(Histos[processBin_nom]) 
				if(proc=='trueH_WH'): sig_WH_year.Add(Histos[processBin_nom]) 
				if(proc=='trueH_ZH'): sig_ZH_year.Add(Histos[processBin_nom]) 
				if(proc=='trueH_ttH'): sig_ttH_year.Add(Histos[processBin_nom]) 
				if(proc=='bkg_qqzz'): qqzz_year.Add(Histos[processBin_nom]); irr_bkg_year.Add(Histos[processBin_nom]) 
				if(proc=='bkg_ggzz'): ggzz_year.Add(Histos[processBin_nom]); irr_bkg_year.Add(Histos[processBin_nom]) 
				#if(proc=='bkg_qqzz'): irr_bkg.Add(Histos[processBin_nom])
				#if(proc=='bkg_qqzz'): irr_bkg.Add(Histos[processBin_nom])
				if(proc!='trueH'): added_year.Add(Histos[processBin_nom])
    		            i=i+1
    
    			    print "yield Integral", Histos[processBin_nom].Integral()
    		            yield_nom+=Histos[processBin_nom].Integral()
			#nom[processBin_nom] = yield_nom 
    
                    elif (year=='2016'):
    		        for pre_tree, post_tree in zip(trees_pre, trees_post):    
    		            print "pre_tree is:  ", pre_tree
    		            print "post_tree is:  ", post_tree
    
    		            if (proc=='obs_data'):
                                #pre_tree.Draw("mass4l >> "+processBin_nom+'pre',"("+cut_nom+")","goff")
                                pre_tree.Draw("mass4l >> "+processBin_nom+'pre',"("+cut_nom+")","goff")
                                #post_tree.Draw("mass4l >> "+processBin_nom+'post',"("+cut_nom+")","goff")
                                post_tree.Draw("mass4l >> "+processBin_nom+'post',"("+cut_nom+")","goff")
#                                data.Add(Histos[processBin_nom+'pre']); data.Add(Histos[processBin_nom+'post']); 		
				#data = Histos[processBin_nom+'pre'] + Histos[processBin_nom+'post']
				data_year.Add(Histos[processBin_nom+'pre'],Histos[processBin_nom+'post'])
				#data_year.Add(Histos[processBin_nom+'pre']); data_year.Add(Histos[processBin_nom+'post'])
				Histos[processBin_nom].Add(Histos[processBin_nom+'pre']); Histos[processBin_nom].Add(Histos[processBin_nom+'post']);
				Histos[processBin_nom+'comb'].Add(Histos[processBin_nom+'pre']); Histos[processBin_nom+'comb'].Add(Histos[processBin_nom+'post'])
				
                            else:
                                print "proc is:   ", proc
                                print "channel is:   ", channel
                                print "i:  ",i
#                                print "len(sumw2_pre):    ", len(sumw2_pre)
#                                print "sumw2_pre[",i,"] : ", sumw2_pre[i]
#                                print "len(sumw2_post):    ", len(sumw2_post)
#                                print "sumw2_post[",i,"] : ", sumw2_post[i]
				#if (proc=='trueH' or proc=='trueH_ggH' or proc=='trueH_XH'): # using only postVFP for signal reweighted to 36.33 
				if (proc=='trueH' or proc=='trueH_ggH' or proc=='trueH_XH' or proc=='trueH_VBF' or proc=='trueH_WH' or proc=='trueH_ZH' or proc=='trueH_ttH'): # using only postVFP for signal reweighted to 36.33 
				    post_tree.Draw("mass4l >> "+processBin_nom+'post',"("+cut_nom+"*("+str(lumi[year])+")/"+str(sumw2_post[i])+")","goff")
				elif (proc=='bkg_ggzz' or proc=='bkg_qqzz'):
			            pre_tree.Draw("mass4l >> "+processBin_nom+'pre',"("+cut_nom+"*("+str(lumi_preVFP)+")/"+str(sumw2_pre[i])+")","goff")
                                    post_tree.Draw("mass4l >> "+processBin_nom+'post',"("+cut_nom+"*("+str(lumi_postVFP)+")/"+str(sumw2_post[i])+")","goff")
				#Histos[processBin_nom+'comb'] = Histos[processBin_nom+'pre'].Add(Histos[processBin_nom+'post'])
				Histos[processBin_nom+'comb'].Add(Histos[processBin_nom+'pre'],Histos[processBin_nom+'post'])
				Histos[processBin_nom].Add(Histos[processBin_nom+'pre'],Histos[processBin_nom+'post']);

#				Histos[processBin_nom+'pre'].SaveAs(proc+'_'+channel+'_'+year+'_qqZZ_pre.root')
#				Histos[processBin_nom+'post'].SaveAs(proc+'_'+channel+'_'+year+'_qqZZ_post.root')
#				Histos[processBin_nom+'comb'].SaveAs(proc+'_'+channel+'_'+year+'_qqZZ_comb.root')
#				print "Histos[processBin_nom+'pre'].Integral()", Histos[processBin_nom+'pre'].Integral()
#				print "Histos[processBin_nom+'post'].Integral()", Histos[processBin_nom+'post'].Integral()
#				print "Histos[processBin_nom+'comb'].Integral()", Histos[processBin_nom+'comb'].Integral()

				if(proc!='trueH'): added_year.Add(Histos[processBin_nom])


				if(proc=='trueH'): sig_year.Add(Histos[processBin_nom]); #added.Add(sig);
				if(proc=='trueH_XH'): sig_XH_year.Add(Histos[processBin_nom]); #added.Add(sig);
				if(proc=='trueH_ggH'): sig_ggH_year.Add(Histos[processBin_nom]); #added.Add(sig);
                                if(proc=='trueH_VBF'): sig_VBF_year.Add(Histos[processBin_nom]) 
                                if(proc=='trueH_WH'): sig_WH_year.Add(Histos[processBin_nom]) 
                                if(proc=='trueH_ZH'): sig_ZH_year.Add(Histos[processBin_nom]) 
                                if(proc=='trueH_ttH'): sig_ttH_year.Add(Histos[processBin_nom])
				if(proc=="bkg_qqzz"):
				    qqzz_year.Add(Histos[processBin_nom]) 
				    irr_bkg_year.Add(Histos[processBin_nom]) 
                                if (proc=='bkg_ggzz'): 
				    ggzz_year.Add(Histos[processBin_nom]); 
				    irr_bkg_year.Add(Histos[processBin_nom]); 
				print "Histos[processBin_nom+'pre']:  ", Histos[processBin_nom+'pre'].Integral()
				print "Histos[processBin_nom+'post']:  ", Histos[processBin_nom+'post'].Integral()
				print "added_year integral:  ", added_year.Integral()
				#Histos[processBin_nom+'pre'].SaveAs(proc+'_'+channel+'_'+year+'pre.root')			
				#Histos[processBin_nom+'post'].SaveAs(proc+'_'+channel+'_'+year+'post.root')			
				#Histos[processBin_nom].SaveAs(proc+'_'+channel+'_'+year+'both.root')			
			#	sig.SaveAs(proc+'_'+channel+'_'+year+'both.root')			
	

                            i=i+1
    
    			    print "pre Integral", Histos[processBin_nom+'pre'].Integral()
    			    print "post Integral", Histos[processBin_nom+'post'].Integral()

		            yield_nom+=Histos[processBin_nom+'pre'].Integral()+Histos[processBin_nom+'post'].Integral() ## FIXME a check
			#nom[processBin_nom] = yield_nom
			#Histos[processBin_nom] = Histos[processBin_nom+'pre'].Add(Histos[processBin_nom+'post'])  # for later use to fit for Landau function
## commenting out
		#	Histos[processBin_nom].Add(Histos[processBin_nom+'pre'], Histos[processBin_nom+'post'])  # for later use to fit for Landau function

#		irr_bkg = Histos['bkg_qqzz_125_'+obsName+'_'+channel+'_recobin0'].Add(Histos['bkg_ggzz_125_'+obsName+'_'+channel+'_recobin0'])
#		print "irr_bkg intergral:    ", irr_bkg.Integral()

            
                    if(proc=='obs_data'):			    
			htemp = TH1F("htemp","htemp",1300,70,2000); # common to an year
			#htemp = TH1F("htemp","htemp",1300,70,770); # common to an year
                                ### red_bkg ###   # 2018 params..
			if (year=="2018"):
                            if (channel=="4e"):
                                fsum = TF1("fsum","landau(0)")
                                fsum.SetParameter(0,1.0)
                                #fsum.SetParameter(1,141.9)
                                fsum.SetParameter(1,134.395073)
                                #fsum.SetParameter(2,21.3)
                                fsum.SetParameter(2,26.093393)
                            if (channel=="4mu"):
                                fsum = TF1("fsum","landau(0)")
                                fsum.SetParameter(0,1.0)
                                #fsum.SetParameter(1,130.4)
                                fsum.SetParameter(1,122.830152)
                                #fsum.SetParameter(2,15.6)
                                fsum.SetParameter(2,21.622966)
                            if (channel=="2e2mu"):
                                fsum = TF1("fsum","landau(0)+landau(3)")
                                fsum.SetParameter(0,0.55)
                                #fsum.SetParameter(1,131.1)
                                fsum.SetParameter(1,127.821771)  # 2e2mu
                                #fsum.SetParameter(2,18.1)
                                fsum.SetParameter(2,20.097781)
                                fsum.SetParameter(3,0.45)
                                fsum.SetParameter(4,128.529522) # 2mu2e
                                fsum.SetParameter(5,21.337322)
                            if (channel=="4l"):
                                fsum = TF1("fsum","landau(0)+landau(3)+landau(6)+landau(9)")
                                fsum.SetParameter(0,9.8)
                                #fsum.SetParameter(1,141.9)
                                fsum.SetParameter(1,134.395073) # 4e
                                #fsum.SetParameter(2,21.3)
                                fsum.SetParameter(2,26.093393)
                                fsum.SetParameter(3,10.2) # ?
                                #fsum.SetParameter(4,130.4)
                                fsum.SetParameter(4,122.830152) # 4mu 
                                #fsum.SetParameter(5,15.6) # ?
                                fsum.SetParameter(5,21.622966) # ?
                                fsum.SetParameter(6,0.55*20.4) # ?
                                #fsum.SetParameter(7,131.1)
                                fsum.SetParameter(7,127.821771) # 2e2mu
                                #fsum.SetParameter(8,18.1)
                                fsum.SetParameter(8,20.097781) 
                                fsum.SetParameter(9,0.45*20.4) # ?
                                #fsum.SetParameter(10,133.8)
                                fsum.SetParameter(10,128.529522) # 2mu2e
                                #fsum.SetParameter(11,18.9)
                                fsum.SetParameter(11,21.337322)

			    nfill=100000
#                           htemp = TH1F("htemp","htemp",1300,70,2000);
                            htemp.FillRandom("fsum",nfill);
                            #htemp.Scale(21.1/htemp.Integral(htemp.FindBin(70),htemp.FindBin(2000)))
                            Histos[processBin_nom].FillRandom("fsum",nfill);
			    #red_bkg.FillRandom("fsum",nfill);
			    red_bkg_year.FillRandom("fsum",nfill);
			 
                        if (year=="2017"):
                            if (channel=="4e"):
                                fsum = TF1("fsum","landau(0)")
                                fsum.SetParameter(0,1.0)
                                #fsum.SetParameter(1,141.9)
                                fsum.SetParameter(1,138.203687) # 4e
                                #fsum.SetParameter(2,21.3)
                                fsum.SetParameter(2,26.366716)
                            if (channel=="4mu"):
                                fsum = TF1("fsum","landau(0)")
                                fsum.SetParameter(0,1.0)
                                #fsum.SetParameter(1,130.4)
                                fsum.SetParameter(1,121.343247) # 4mu
                                #fsum.SetParameter(2,15.6)
                                fsum.SetParameter(2,21.931045)
                            if (channel=="2e2mu"):
                                fsum = TF1("fsum","landau(0)+landau(3)")
                                fsum.SetParameter(0,0.55)
                                #fsum.SetParameter(1,131.1)
                                fsum.SetParameter(1,129.879999) # 2e2mu
                                #fsum.SetParameter(2,18.1)
                                fsum.SetParameter(2,20.149098)
                                fsum.SetParameter(3,0.45)
                                fsum.SetParameter(4,126.763438) # 2mu2e
                                fsum.SetParameter(5,22.033032)
                            if (channel=="4l"):
                                fsum = TF1("fsum","landau(0)+landau(3)+landau(6)+landau(9)")
                                fsum.SetParameter(0,9.8)
                                #fsum.SetParameter(1,141.9)
                                fsum.SetParameter(1,138.203687) # 4e
                                #fsum.SetParameter(2,21.3)
                                fsum.SetParameter(2,26.366716)
                                fsum.SetParameter(3,10.2) # ?
                                #fsum.SetParameter(4,130.4)
                                fsum.SetParameter(4,121.343247) # 4mu
                                fsum.SetParameter(5,21.931045) # 
                                fsum.SetParameter(6,0.55*20.4) # ?
                                #fsum.SetParameter(7,131.1)
                                fsum.SetParameter(7,129.879999) # 2e2mu
                                #fsum.SetParameter(8,18.1)
                                fsum.SetParameter(8,20.149098)
                                fsum.SetParameter(9,0.45*20.4) # ?
                                #fsum.SetParameter(10,133.8)
                                fsum.SetParameter(10,126.763438)  # 2mu2e
                                #fsum.SetParameter(11,18.9)
                                fsum.SetParameter(11,22.033032)

                            nfill=100000
#                           htemp = TH1F("htemp","htemp",1300,70,2000);
                            htemp.FillRandom("fsum",nfill);
                            #htemp.Scale(21.1/htemp.Integral(htemp.FindBin(70),htemp.FindBin(2000)))
                            Histos[processBin_nom].FillRandom("fsum",nfill);
			    #red_bkg.FillRandom("fsum",nfill);
			    red_bkg_year.FillRandom("fsum",nfill);

                        if (year=="2016"):
                        #############################
                             #####  preVFP  #####
                        ############################
			
                            if (channel=="4e"):  
                                fsum_pre = TF1("fsum_pre","landau(0)")
                                fsum_pre.SetParameter(0,1.0)
                                #fsum_pre.SetParameter(1,141.9)
                                fsum_pre.SetParameter(1,136.052353) # 4e
                                #fsum_pre.SetParameter(2,21.3)
                                fsum_pre.SetParameter(2,27.295657)
                            if (channel=="4mu"):
                                fsum_pre = TF1("fsum_pre","landau(0)")
                                fsum_pre.SetParameter(0,1.0)
                                #fsum_pre.SetParameter(1,130.4)
                                fsum_pre.SetParameter(1,131.436563) # 4mu
                                #fsum_pre.SetParameter(2,15.6)
                                fsum_pre.SetParameter(2,24.447732)
                            if (channel=="2e2mu"):
                                fsum_pre = TF1("fsum_pre","landau(0)+landau(3)")
                                fsum_pre.SetParameter(0,0.55)
                                #fsum_pre.SetParameter(1,131.1)
                                fsum_pre.SetParameter(1,121.165524) # 2e2mu
                                #fsum_pre.SetParameter(2,18.1)
                                fsum_pre.SetParameter(2,21.536749)
                                fsum_pre.SetParameter(3,0.45)
                                fsum_pre.SetParameter(4,132.162074) # 2mu2e
                                fsum_pre.SetParameter(5,24.088944)
                            if (channel=="4l"):
                                fsum_pre = TF1("fsum_pre","landau(0)+landau(3)+landau(6)+landau(9)")
                                fsum_pre.SetParameter(0,9.8)
                                #fsum_pre.SetParameter(1,141.9)
                                fsum_pre.SetParameter(1,136.052353) # 4e
                                #fsum_pre.SetParameter(2,21.3)
                                fsum_pre.SetParameter(2,27.295657)
                                fsum_pre.SetParameter(3,10.2) # ?
                                #fsum_pre.SetParameter(4,130.4)
                                fsum_pre.SetParameter(4,131.436563) # 4mu
                                fsum_pre.SetParameter(5,24.447732) # 
                                fsum_pre.SetParameter(6,0.55*20.4) # ?
                                #fsum_pre.SetParameter(7,131.1)
                                fsum_pre.SetParameter(7,121.165524) # 2e2mu
                                #fsum_pre.SetParameter(8,18.1)
                                fsum_pre.SetParameter(8,21.536749)
                                fsum_pre.SetParameter(9,0.45*20.4) # ?
                                #fsum_pre.SetParameter(10,133.8)
                                fsum_pre.SetParameter(10,132.162074)  # 2mu2e
                                #fsum_pre.SetParameter(11,18.9)
                                fsum_pre.SetParameter(11,24.088944)

                        #if (year=="2016post"):
                        #############################
                             #####  postVFP  #####
                        ############################
                            if (channel=="4e"):
                                fsum_post = TF1("fsum_post","landau(0)")
                                fsum_post.SetParameter(0,1.0)
                                #fsum_post.SetParameter(1,141.9)
                                fsum_post.SetParameter(1,131.729663) # 4e
                                #fsum_post.SetParameter(2,21.3)
                                fsum_post.SetParameter(2,25.804455)
                            if (channel=="4mu"):
                                fsum_post = TF1("fsum_post","landau(0)")
                                fsum_post.SetParameter(0,1.0)
                                #fsum_post.SetParameter(1,130.4)
                                fsum_post.SetParameter(1,120.000001) # 4mu
                                #fsum_post.SetParameter(2,15.6)
                                fsum_post.SetParameter(2,21.132331)
                            if (channel=="2e2mu"):
                                fsum_post = TF1("fsum_post","landau(0)+landau(3)")
                                fsum_post.SetParameter(0,0.55)
                                #fsum_post.SetParameter(1,131.1)
                                fsum_post.SetParameter(1,127.649715) # 2e2mu
                                #fsum_post.SetParameter(2,18.1)
                                fsum_post.SetParameter(2,18.023413)
                                fsum_post.SetParameter(3,0.45)
                                fsum_post.SetParameter(4,137.186696) # 2mu2e
                                fsum_post.SetParameter(5,24.621237)
                            if (channel=="4l"):
                                fsum_post = TF1("fsum_post","landau(0)+landau(3)+landau(6)+landau(9)")
                                fsum_post.SetParameter(0,9.8)
                                #fsum_post.SetParameter(1,141.9)
                                fsum_post.SetParameter(1,131.729663) # 4e
                                #fsum_post.SetParameter(2,21.3)
                                fsum_post.SetParameter(2,25.804455)
                                fsum_post.SetParameter(3,10.2) # ?
                                #fsum_post.SetParameter(4,130.4)
                                fsum_post.SetParameter(4,120.000001) # 4mu
                                fsum_post.SetParameter(5,21.132331) # 
                                fsum_post.SetParameter(6,0.55*20.4) # ?
                                #fsum_post.SetParameter(7,131.1)
                                fsum_post.SetParameter(7,127.649715) # 2e2mu
                                #fsum_post.SetParameter(8,18.1)
                                fsum_post.SetParameter(8,18.023413)
                                fsum_post.SetParameter(9,0.45*20.4) # ?
                                #fsum_post.SetParameter(10,133.8)
                                fsum_post.SetParameter(10,137.186696)  # 2mu2e
                                #fsum_post.SetParameter(11,18.9)
                                fsum_post.SetParameter(11,24.621237)
			    nfill=100000
			    htemp_pre = TH1F("htemp_pre","htemp_pre",1300,70,2000); htemp_pre.FillRandom("fsum_pre",nfill); #htemp.Add(htemp_pre)
			    #htemp_pre = TH1F("htemp_pre","htemp_pre",1300,70,770); htemp_pre.FillRandom("fsum_pre",nfill); #htemp.Add(htemp_pre)
			    htemp_post = TH1F("htemp_post","htemp_post",1300,70,2000); htemp_post.FillRandom("fsum_post",nfill); #htemp.Add(htemp_post)
			    #htemp_post = TH1F("htemp_post","htemp_post",1300,70,770); htemp_post.FillRandom("fsum_post",nfill); #htemp.Add(htemp_post)
			    htemp.Add(htemp_pre, htemp_pre);
			    #htemp.Add(htemp_pre); htemp.Add(htemp_pre);
			    #htemp.FillRandom("fsum",nfill);

			
                            Histos[processBin_nom+'pre'].FillRandom("fsum_pre",nfill); #Histos[processBin_nom].Add(Histos[processBin_nom+'pre'])
                            Histos[processBin_nom+'post'].FillRandom("fsum_post",nfill); #Histos[processBin_nom].Add(Histos[processBin_nom+'post'])
			    Histos[processBin_nom].Add(Histos[processBin_nom+'pre'], Histos[processBin_nom+'post'])
			    #Histos[processBin_nom].Add(Histos[processBin_nom+'pre']);
			    #Histos[processBin_nom].Add(Histos[processBin_nom+'post']);

			    #red_bkg_pre.FillRandom("fsum_pre",nfill); #red_bkg.Add(red_bkg_pre)
			    #red_bkg_post.FillRandom("fsum_post",nfill); #red_bkg.Add(red_bkg_post)
			    #red_bkg.Add(red_bkg_pre, red_bkg_post);

			    red_bkg_pre_year.FillRandom("fsum_pre",nfill); #red_bkg.Add(red_bkg_pre)
			    red_bkg_post_year.FillRandom("fsum_post",nfill); #red_bkg.Add(red_bkg_post)
			    red_bkg_year.Add(red_bkg_pre_year, red_bkg_post_year);
			    #red_bkg_year.Add(red_bkg_pre); red_bkg_year.Add(red_bkg_post);


#                        nfill=100000
#                        htemp = TH1F("htemp","htemp",1300,70,2000);
#                        htemp.FillRandom("fsum",nfill);
			#htemp.Scale(21.1/htemp.Integral(htemp.FindBin(70),htemp.FindBin(2000)))
			print "zx_contribution[year][channel]", zx_contribution[str(year)][channel]
			htemp.Scale(zx_contribution[str(year)][channel]/htemp.Integral(htemp.FindBin(70),htemp.FindBin(2000)))
			#htemp.Scale(zx_contribution[str(year)][channel]/htemp.Integral(htemp.FindBin(70),htemp.FindBin(770)))
# Filippo's selection 20,70,770

		
#                        if (channel=="4e"): htemp.Scale(21.1/htemp.Integral(htemp.FindBin(70),htemp.FindBin(2000)))
#                        if (channel=="4mu"): htemp.Scale(34.4/htemp.Integral(htemp.FindBin(70),htemp.FindBin(2000)))
#                        if (channel=="2e2mu"): htemp.Scale(59.9/htemp.Integral(htemp.FindBin(70),htemp.FindBin(2000)))
#                        if (channel=="4l"): htemp.Scale(115.4/htemp.Integral(htemp.FindBin(70),htemp.FindBin(2000)))
                        #Histos[processBin_nom].FillRandom("fsum",nfill);
                        #Histos[processBin_nom].SaveAs(proc+'_'+channel+'_'+year+'.root');


#			red_bkg.FillRandom("fsum",nfill);
			if (red_bkg_year.Integral(red_bkg_year.FindBin(m4l_low),red_bkg_year.FindBin(m4l_high)) >0 ):	
			    red_bkg_year.Scale(htemp.Integral(htemp.FindBin(m4l_low),htemp.FindBin(m4l_high))/red_bkg_year.Integral(red_bkg_year.FindBin(m4l_low),red_bkg_year.FindBin(m4l_high)))
			print 'Red Bkg year:  ', red_bkg_year.Integral(1,m4l_bins)
			yield_nom = red_bkg_year.Integral(1,m4l_bins)

			added_year.Add(red_bkg_year); # TJ should it be added ?

                        if (Histos[processBin_nom].Integral(Histos[processBin_nom].FindBin(m4l_low),Histos[processBin_nom].FindBin(m4l_high)) >0 ):
                            #Histos[processBin_nom].Scale(htemp.Integral(htemp.FindBin(105),htemp.FindBin(160))/Histos[processBin_nom].Integral(Histos[processBin_nom].FindBin(105),red_bkg.FindBin(160)))
                            Histos[processBin_nom].Scale(htemp.Integral(htemp.FindBin(m4l_low),htemp.FindBin(m4l_high))/Histos[processBin_nom].Integral(Histos[processBin_nom].FindBin(m4l_low),Histos[processBin_nom].FindBin(m4l_high)))
        	  	    #nom[processBin_nom] = Histos[processBin_nom].Integral(1,m4l_bins)	
        
#       		        print 'Red Bkg:  ',Histos[processBin_nom].Integral(1,m4l_bins)
			#print 'Red Bkg year:  ', red_bkg_year.Integral(1,m4l_bins)
                        
		    nom[processBin_nom] = yield_nom
		    #nom[processBin_nom] = Histos[processBin_nom].Integral(1,m4l_bins)  #  nominal rate for cards
                                    
            
    added.Add(added_year); red_bkg.Add(red_bkg_year); irr_bkg.Add(irr_bkg_year); qqzz.Add(qqzz_year); ggzz.Add(ggzz_year); sig.Add(sig_year); sig_XH.Add(sig_XH_year);sig_ggH.Add(sig_ggH_year); sig_VBF.Add(sig_VBF_year);sig_WH.Add(sig_WH_year);sig_ZH.Add(sig_ZH_year);sig_ttH.Add(sig_ttH_year);data.Add(data_year);
#    return added_year, red_bkg_year, irr_bkg_year, qqzz_year, ggzz_year, sig_year, data_year;
    if (not os.path.exists("datacardInputs/"+year+"/"+obs_reco)):
        os.system("mkdir -p datacardInputs/"+year+"/"+obs_reco+"")

    with open('datacardInputs/'+year+'/'+obs_reco+'/inputs_yields_'+opt.OBSNAME+'_'+channel+'_'+opt.RNG+'.py', 'w') as f:
        f.write('yield_nom = '+str(nom)+' \n')


            #Histos['bkg_qqzz_125_'+obsName+'_'+channel+'_recobin0'].Add(Histos['bkg_ggzz_125_'+obsName+'_'+channel+'_recobin0'])
#            print "irr_bkg intergral:    ", irr_bkg.Integral()



#    for x_point in x_points:
#        for channel in channels:
#            for recobin in range(len(obs_bins)-1):
#                processBin_nom_Hproc = 'Hproc_'+str(x_point)+'_'+obsName+'_'+channel+'_recobin'+str(recobin)
#                nom[processBin_nom_Hproc]=0;
#                yield_nom_Hproc=0;
#                for proc in procs:
#                    #if(not (proc=='trueH' or proc=='out_trueH' or proc=='fakeH')): continue
#                    if(not (proc=='trueH')): continue
#                    processBin_nom = proc+'_'+str(x_point)+'_'+obsName+'_'+channel+'_recobin'+str(recobin)
#                    yield_nom_Hproc+=nom[processBin_nom]
#                    #nom[processBin_nom_Hproc]=0;
##                   print "processBin_nom  ", processBin_nom
##                   print "nom[processBin_nom]  ", nom[processBin_nom]
##                   print "yield_nom_Hproc   ", yield_nom_Hproc
#                nom[processBin_nom_Hproc]=yield_nom_Hproc;
##               print "processBin_nom_Hproc  ", processBin_nom_Hproc
##               print "nom[processBin_nom_Hproc]  ", nom[processBin_nom_Hproc]
#
## interpolation part
#    for channel in channels:
#        for proc in procs:
#           #if(not (proc=='trueH' or proc=='out_trueH' or proc=='fakeH')): continue
#           if(not (proc=='trueH')): continue
#           for recobin in range(len(obs_bins)-1):
#               for nuis in JES_nuis:
#                   key_nom_125pt38 = proc+'_125.38_'+obsName+'_'+channel+'_recobin'+str(recobin)
#                   nom_points = []; 
#                   for x_point in x_points:
#                       key_nom = proc+'_'+str(x_point)+'_'+obsName+'_'+channel+'_recobin'+str(recobin)
#                       nom_points.append(nom[key_nom]); 
#    
#                   print "nom_points:   ", nom_points
#                   spl_nom_points  = interpolate.splrep(x_points, nom_points, k=2)
#                 
#                   nom[key_nom_125pt38] = float(interpolate.splev(125.38, spl_nom_points))


            #with open('datacardInputs/'+year+'/'+obs_reco+'/inputs_yields_'+opt.OBSNAME+'.py', 'w') as f:
#    with open('datacardInputs/'+year+'/'+obs_reco+'/inputs_yields_'+opt.OBSNAME+'.py', 'w') as f:
#            with open('datacardInputs/'+year+'/'+obs_reco+'/inputs_yields_'+opt.OBSNAME+'_'+channel+'.py', 'w') as f:
#        	f.write('yield_nom = '+str(nom)+' \n')

if __name__ == "__main__":
    global opt, args
    parseOptions()
    #sys.argv = grootargs
    observableBins = {0:(opt.OBSBINS.split("|")[1:(len(opt.OBSBINS.split("|"))-1)]),1:['0','inf']}[opt.OBSBINS=='inclusive']
    nbins = len(observableBins)

    obs_bins = opt.OBSBINS.split("|") 
    if (not (obs_bins[0] == '' and obs_bins[len(obs_bins)-1]=='')): 
        print 'BINS OPTION MUST START AND END WITH A |' 
    obs_bins.pop()
    obs_bins.pop(0)
    print "obs_bins:  ", obs_bins
    nbins = len(observableBins)
    
    obs_reco = opt.OBSNAME 
    obsName = opt.OBSNAME 

#    if (not os.path.exists("datacardInputs/"+year+"/"+obs_reco)):
#        os.system("mkdir -p datacardInputs/"+year+"/"+obs_reco+"")

    window = {"signal":{"lo":105.0,"hi":160.0,"bin":55},
   # window = {"signal":{"lo":105.0,"hi":160.0,"bin":500},
    #window = {"signal":{"lo":105.0,"hi":160.0,"bin":40},
    "full":{"lo":70.0,"hi":1000.0,"bin":90},
    #"full":{"lo":70.0,"hi":1000.0,"bin":930},
    #"low":{"lo":70.0,"hi":160.0,"bin":55},
    "low":{"lo":70.0,"hi":160.0,"bin":90},
    #"low":{"lo":70.0,"hi":160.0,"bin":95},
    #"high":{"lo":160,"hi":1000.0,"bin":75}
    "high":{"lo":160,"hi":1000.0,"bin":90}
    }

    m4l_low = window[opt.RNG]["lo"]
    m4l_high = window[opt.RNG]["hi"]
    m4l_bins = window[opt.RNG]["bin"]
 
 
    #m4l_low = 105.0
    #m4l_high = 160.0
    #m4l_bins = 55 # OK

    #m4l_bins = 35
    #m4l_bins = 40
    #m4l_bins = 55
    #m4l_bins = 70  ## temporary


    #m4l_low = 70.0
    #m4l_high = 1000.0
    #m4l_bins = 70  ## temporary

#    m4l_low = 40.0
#    m4l_high = 130.0
#    m4l_bins = 70  ## temporary


    print "era being processed: ", opt.ERA


    savevar = opt.OBSNAME

    #lumiplot = {"2018":"59.8 fb^{-1}","2017":"41.5 fb^{-1}","2016":"35.9 fb^{-1}", "Full":"138 fb^{-1}"}
    #lumiplot = {"2018":"59.8 fb^{-1}","2017":"41.5 fb^{-1}","2016":"36.3 fb^{-1}", "Full":"138 fb^{-1}"}
    lumiplot = {"2018":"59.8 fb^{-1}","2017":"41.5 fb^{-1}","2016":"36.31 fb^{-1}", "Full":"138 fb^{-1}"} # Reference : https://twiki.cern.ch/twiki/bin/viewauth/CMS/TWikiLUM
    lumiplot2018 = '58.8 fb^{-1}'

    fstates = {"4l":"m_{4l}","4e":"m_{4e}","4mu":"m_{4#mu}","2e2mu":"2e2#mu"}

    #save = savevar+'_'+str(m4l_low).replace('.','pt')+'to'+str(m4l_high).replace('.','pt')+'_'+channel+'_mh'+mh

    doratio = True
    #doratio = False



    print("Obs Name: {:15}  nBins: {:2}  bins: {}".format(opt.OBSNAME, nbins, observableBins))

    if (opt.ERA == '2016'): years = ['2016']
    if (opt.ERA == '2017'): years = ['2017']
    if (opt.ERA == '2018'): years = ['2018']
    if (opt.ERA == 'Full'): years = ['2016','2017','2018']

    channels=['4mu','2e2mu','4e', '4l']
    #channels=['4mu'] #
    #channels=['4l'] #
    #channels=['2e2mu'] #

    for channel in channels:

        stack = THStack('stack', 'stack')
        added = TH1D('a', 'a',m4l_bins,m4l_low,m4l_high)
        red_bkg = TH1D('red_bkg', 'red_bkg',m4l_bins,m4l_low,m4l_high)
        red_bkg_pre = TH1D('red_bkg_pre', 'red_bkg_pre',m4l_bins,m4l_low,m4l_high)
        red_bkg_post = TH1D('red_bkg_post', 'red_bkg_post',m4l_bins,m4l_low,m4l_high)
        irr_bkg = TH1D('irr_bkg', 'irr_bkg',m4l_bins,m4l_low,m4l_high)
        rare_bkg = TH1D('rare_bkg', 'rare_bkg',m4l_bins,m4l_low,m4l_high)
        qqzz = TH1D('qqzz', 'qqzz',m4l_bins,m4l_low,m4l_high)
        ggzz = TH1D('ggzz', 'ggzz',m4l_bins,m4l_low,m4l_high)
        sig = TH1D('sig', 'sig',m4l_bins,m4l_low,m4l_high)
        sig_XH = TH1D('sig_XH', 'sig_XH',m4l_bins,m4l_low,m4l_high)
        sig_VBF = TH1D('sig_VBF', 'sig_VBF',m4l_bins,m4l_low,m4l_high)
        sig_WH = TH1D('sig_WH', 'sig_WH',m4l_bins,m4l_low,m4l_high)
        sig_ZH = TH1D('sig_ZH', 'sig_ZH',m4l_bins,m4l_low,m4l_high)
        sig_ttH = TH1D('sig_ttH', 'sig_ttH',m4l_bins,m4l_low,m4l_high)
        sig_ggH = TH1D('sig_ggH', 'sig_ggH',m4l_bins,m4l_low,m4l_high)
        data = TH1D('data', 'data',m4l_bins,m4l_low,m4l_high)

        for year in years:
            #extract_yields(nbins, opt.OBSNAME, obs_bins, channel, years, 'm_{4l}', 'GeV', True, False, False,"-1","125", opt.DEBUG)
            #extract_yields(nbins, opt.OBSNAME, obs_bins, channel, year, 'm_{4l}', 'GeV', True, False, False,"-1","125", opt.DEBUG)
            extract_yields(nbins, opt.OBSNAME, obs_bins, channel, year, opt.DEBUG) #, 'm_{4l}', 'GeV', True, False, False,"-1","125", opt.DEBUG)
            #extract_yields(nbins, opt.OBSNAME, obs_bins, year, 'm_{4l}', 'GeV', True, False, False,"-1","125", opt.DEBUG)
    
    #        added.Add(added_year); red_bkg.Add(red_bkg_year); irr_bkg.Add(irr_bkg_year); qqzz.Add(qqzz_year); ggzz.Add(ggzz_year); sig.Add(sig_year); data.Add(data_year);


	#xlabel= 'm_{4l}'; 
	xlabel = fstates[channel];
	xunits= 'GeV'; prelim=True; setLogX=False; setLogY=False; EventCat='-1'; mh='125'
 
        save = savevar+'_'+str(m4l_low).replace('.','pt')+'to'+str(m4l_high).replace('.','pt')+'_'+channel+'_mh'+mh+opt.ERA
    
    
        if (obsName=='mass'+channel): c1 = TCanvas("c1","c1", 1200, 800)
        c1 = TCanvas("c1","c1", 800, 800)
        if (setLogX): c1.SetLogx()
        if (setLogY): c1.SetLogy()
        if (doratio): c1.SetBottomMargin(0.3)
        c1.SetRightMargin(0.03);
    
        ################################
        #### Reducible Background ###
        ################################
    
        #red_bkg.Smooth(1)
        red_bkg.SetLineColor(TColor.GetColor("#003300"))
        red_bkg.SetLineWidth(2)
        red_bkg.SetFillColor(TColor.GetColor("#669966"))
    ############################
     ###Rare Backgrounds ###
    ############################
        #rare_bkg.SetFillColor(kOrange+2)
        #rare_bkg.SetLineColor(kOrange+3)
        #rare_bkg.SetLineWidth(2)
        ################################
        #### Irreducible Backgrounds ###
        ################################
    
        irr_bkg.SetLineColor(TColor.GetColor("#000099"))
        irr_bkg.SetLineWidth(2)
        irr_bkg.SetFillColor(TColor.GetColor("#3366ff"))
    
        qqzz.SetLineColor(TColor.GetColor("#000099"))
        qqzz.SetLineWidth(2)
        qqzz.SetFillColor(TColor.GetColor("#99ccff"))
    
        ggzz.SetFillColor(TColor.GetColor("#3366ff"))
        ggzz.SetLineWidth(2)
        ggzz.SetLineColor(TColor.GetColor("#000099"))
    
        sig.SetFillColor(TColor.GetColor("#ffafaf"))
        sig.SetLineColor(TColor.GetColor("#cc0000"))
        sig.SetLineWidth(2)
 

        sig_ggH.SetFillColor(TColor.GetColor("#ffafaf"))
        sig_ggH.SetLineColor(TColor.GetColor("#cc0000"))
        sig_ggH.SetLineWidth(2)



   
        #sig_XH.SetFillColor(TColor.GetColor(ROOT.kGreen+3))
        #sig_XH.SetLineColor(TColor.GetColor(ROOT.kGreen+4))
        #sig_XH.SetLineWidth(2)

        #sig_XH.SetFillColor(TColor.GetColor("#886688"))
        #sig_XH.SetLineColor(TColor.GetColor("#5e2557"))
        #sig_XH.SetLineWidth(2)

#        sig_XH.SetLineColor(TColor.GetColor("#7ba487")) # light green, OK
#        sig_XH.SetFillColor(TColor.GetColor("#9acda9"))
#        sig_XH.SetLineWidth(2)


#        sig_XH.SetLineColor(TColor.GetColor("#e45c85")) # light pink
#        sig_XH.SetFillColor(TColor.GetColor("#ee99b3"))
#        sig_XH.SetLineWidth(2)
        #stack.Add(rare_bkg)

        #sig_XH.SetLineColor(TColor.GetColor("#60c3a2")) # light zink
        #sig_XH.SetFillColor(TColor.GetColor("#9bdac5"))
        #sig_XH.SetLineWidth(2)

#        sig_XH.SetLineColor(TColor.GetColor("#60c3a2")) # light zink
#        sig_XH.SetFillColor(TColor.GetColor("#cc7c7c"))
#        sig_XH.SetLineWidth(2)

        sig_XH.SetLineColor(TColor.GetColor("#7c3c5e")) # Dark pink 
        sig_XH.SetFillColor(TColor.GetColor("#e68cb0"))
        sig_XH.SetLineWidth(2)

        #stack.Add(rare_bkg)
        #if (obsName=='mass'+channel and EventCat=="-1"): 
        stack.Add(red_bkg)
        stack.Add(ggzz)
        stack.Add(qqzz)
        #stack.Add(irr_bkg)
        stack.Add(sig_XH)  # temporary commenting out.
        stack.Add(sig_ggH)  # temporary commenting out.
#        stack.Add(sig)
    
        stack.Draw("hist")
        stack.GetXaxis().SetMoreLogLabels(kTRUE)
        stack.GetXaxis().SetNoExponent(kTRUE)
        stack.SetMinimum(0.01)
        if (setLogY): stack.SetMaximum(50*max(stack.GetMaximum(),data.GetMaximum()+1))
        else: stack.SetMaximum(1.5*max(stack.GetMaximum(),data.GetMaximum()+1))
        if (xunits==''): stack.GetXaxis().SetTitle(xlabel)
        else: stack.GetXaxis().SetTitle(xlabel+' ['+xunits+']')
        if (doratio): stack.GetXaxis().SetTitleSize(0)
        if (doratio): stack.GetXaxis().SetLabelSize(0)
        #binsize = str(round(float((m4l_high-m4l_low)/nbins),2))
        binsize = str(round(float((m4l_high-m4l_low)/m4l_bins),2))
        if binsize.endswith('.0'): binsize=binsize.rstrip('.0')
        if (xunits==''): ylabel = 'Events / bin'
        else: ylabel = 'Events / '+binsize+' '+xunits
        stack.GetYaxis().SetTitle(ylabel)
    
    
        data.SetMarkerStyle(20)
        data.SetMarkerSize(1.2)
        data.SetBinErrorOption(TH1.kPoisson)
        data.Draw("ex0psame")
    
        #legend = TLegend(.62,.70,.90,.92)
        #legend = TLegend(.55,.70,.90,.92) # TJ
        #legend = TLegend(.55,.70,.90,.91) # TJ, OK
        #legend = TLegend(.52,.65,.92,.92) # TJ, looks OK
        #legend = TLegend(.55,.69,.92,.92) # TJ, looks OK 2
        #legend = TLegend(.54,.68,.92,.92) # TJ, looks OK 2
        legend = TLegend(.54,.68,.90,.92) # TJ, looks OK 2
        #legend = TLegend(.2,.30,.52,.60)
        legend.AddEntry(data, 'Data', "ep")
        #legend.AddEntry(sig, 'm_{H} = 125 GeV', "f")
        #legend.AddEntry(sig, 'm_{H} = 125 GeV (gg#rightarrowH (POWHEG) + XH)', "f")

        #legend.AddEntry(sig_ggH, 'gg#rightarrowH (POWHEG) + XH   (m_{H} = 125 GeV)', "f")

        #legend.AddEntry(sig_ggH, 'gg#rightarrowH (POWHEG)   (m_{H} = 125 GeV)', "f")
        legend.AddEntry(sig_ggH, 'gg #rightarrow H   (m_{H} = 125 GeV)', "f")

        #legend.AddEntry(sig, 'gg#rightarrowH (POWHEG) + XH   (m_{H} = 125 GeV)', "f")


        #legend.AddEntry(sig_XH, 'm_{H} = 125 GeV (sub-dominant signal only)', "f")
        #legend.AddEntry(sig_XH, 'XH = VBF + VH + ttH   (m_{H} = 125 GeV)', "f")
        legend.AddEntry(sig_XH, 'XH = VBF + VH + ttH ', "f")
#        legend.AddEntry(irr_bkg, 'Z#gamma*, ZZ', "f")
        legend.AddEntry(qqzz, 'qq #rightarrow ZZ', "f")
        legend.AddEntry(ggzz, 'gg #rightarrow ZZ', "f")
        #if (obsName=="mass4l"): legend.AddEntry(red_bkg, 'Z+X', "f")
        legend.AddEntry(red_bkg, 'Z+X', "f")
	#legend.SetTextFont(45);
        legend.SetShadowColor(0);
        legend.SetFillColor(0);
        legend.SetLineColor(0);
        legend.Draw("same")
    
        latex2 = TLatex()
        latex2.SetNDC()
        latex2.SetTextSize(0.6*c1.GetTopMargin())
        latex2.SetTextFont(42)
        #latex2.SetTextAlign(31) # align right
        latex2.SetTextAlign(31) # align right
    ##    latex2.DrawLatex(0.92, 0.94,lumiplot2016+" (13 TeV)")
    ##    latex2.DrawLatex(0.92, 0.94,lumiplot2017+" (13 TeV)")
        #latex2.DrawLatex(0.92, 0.94,lumiplot2018+" (13 TeV)")
        latex2.DrawLatex(0.92, 0.94,lumiplot[opt.ERA]+" (13 TeV)")
        #latex2.SetTextSize(1.0*c1.GetTopMargin())
        latex2.SetTextSize(0.7*c1.GetTopMargin())
        latex2.SetTextFont(62)
        latex2.SetTextAlign(11) # align right
        #latex2.DrawLatex(0.2, 0.84, "CMS")
        latex2.DrawLatex(0.156, 0.94, "CMS")
        #latex2.SetTextSize(0.8*c1.GetTopMargin())
        latex2.SetTextSize(0.6*c1.GetTopMargin())
        latex2.SetTextFont(52)
        latex2.SetTextAlign(11)
        #latex2.DrawLatex(0.28, 0.945, "Unpublished")    
        #latex2.DrawLatex(0.2, 0.78, "Preliminary")
        latex2.DrawLatex(0.27, 0.94, "Preliminary")
    
        lastbin=m4l_bins
        if (obsName=='mass4l' and m4l_high>700.0): lastbin=m4l_bins+1  # ??
    
        print 'Signal:  ',sig.Integral(1,lastbin)
        print 'Signal   XH:  ',sig_XH.Integral(1,lastbin)
        print 'Signal ggH:  ',sig_ggH.Integral(1,lastbin)
        #print 'Signal   XH:  ',sig_XH.Integral(1,lastbin)
#        print 'Signal - XH  :  ',(sig.Integral(1,lastbin)-sig_XH.Integral(1,lastbin));
        print 'Signal VBF:  ',sig_VBF.Integral(1,lastbin)
        print 'Signal WH:  ',sig_WH.Integral(1,lastbin)
        print 'Signal ZH:  ',sig_ZH.Integral(1,lastbin)
        print 'Signal ttH:  ',sig_ttH.Integral(1,lastbin)

        print 'ZZ Bkg:  ',irr_bkg.Integral(1,lastbin)
        print 'ggZZ Bkg:  ',ggzz.Integral(1,lastbin)
        print 'qqZZ Bkg:  ',qqzz.Integral(1,lastbin)
        print 'Red Bkg:  ',red_bkg.Integral(1,lastbin)
    ##    print 'Rare Bkg:  ',rare_bkg.Integral(1,lastbin)
    ##    print 'Signal:  ',sig.Integral(1,lastbin)
        print 'Data  :  ',data.Integral(1,lastbin)
    
        if (doratio):
    
            ratio = data.Clone('ratio')
            ratio.Divide(added)
            ratio.GetXaxis().SetMoreLogLabels(kTRUE)
            ratio.GetXaxis().SetNoExponent(kTRUE)
    
            pad = TPad("pad", "pad", 0.0, 0.0, 1.0, 1.0);
            if (setLogX): pad.SetLogx()
            pad.SetTopMargin(0.7);
            pad.SetRightMargin(0.03);
            pad.SetFillColor(0);
            pad.SetGridy(1);
            pad.SetFillStyle(0);
            pad.Draw();
            pad.cd(0);
    
            ratio.SetLineColor(1);
            ratio.SetMarkerColor(1);
            if (xunits==''):
                ratio.GetXaxis().SetTitle(xlabel)
            else:
                ratio.GetXaxis().SetTitle(xlabel+' ['+xunits+']')
            ratio.GetXaxis().SetTitle
            #ratio.GetYaxis().SetTitleSize(0.04);
            ratio.GetYaxis().SetTitleSize(0.03);
            ratio.GetYaxis().SetTitleOffset(1.8);
            #ratio.GetYaxis().SetTitle("Data/Bkg.");
            ratio.GetYaxis().SetTitle("Data/MC");
            ratio.GetYaxis().CenterTitle();
            ratio.GetYaxis().SetLabelSize(0.03);
            ratio.SetMarkerStyle(20);
            ratio.SetMarkerSize(1.2);
            ratio.SetLineWidth(2)
            ratio.SetLineStyle(1)
            ratio.SetLineColor(1)
            ratio.SetMarkerColor(1)
            #ratio.SetMinimum(0.51);
            ratio.SetMinimum(0.41);
            #ratio.SetMaximum(1.69);
            ratio.SetMaximum(1.79);
            ratio.Draw("ep");
    
        c1.SaveAs('distributions/Histo_' + save + '.pdf')
        c1.SaveAs('distributions/Histo_' + save + '.png')
    
    
        print '=========================================='
        print ''
    
        del c1
    



    print(" completed... :) ")
