import optparse
import os
import sys
import ROOT, re
from ROOT import *
from scipy import interpolate
import matplotlib.pyplot as plt
#from library import list_files as DATA
from LoadData import *
#from JES_branches import *
#from Input_Info import datacardInputs
#from Utils import *


if (not os.path.exists("datacardInputs")):
    os.system("mkdir -p datacardInputs")


def parseOptions():

    global opt, args, runAllSteps

    usage = ('usage: %prog [options]\n'
             + '%prog -h for help')
    parser = optparse.OptionParser(usage)
    #parser.add_option('',   '--obsName',  dest='OBSNAME',  type='string',default='mass4l',   help='Name of the observable, supported: "inclusive", "pT4l", "eta4l", "massZ2", "nJets"')
    #parser.add_option('',   '--obsName',  dest='OBSNAME',  type='string',default='pt_leadingjet_pt30_eta4p7',   help='Name of the observable, supported: "inclusive", "pT4l", "eta4l", "massZ2", "nJets"')
    #parser.add_option('',   '--obsName',  dest='OBSNAME',  type='string',default='mass4lj',   help='Name of the observable, supported: "inclusive", "pT4l", "eta4l", "massZ2", "nJets"')
    parser.add_option('',   '--obsName',  dest='OBSNAME',  type='string',default='mass4ljj',   help='Name of the observable, supported: "inclusive", "pT4l", "eta4l", "massZ2", "nJets"')
    #parser.add_option('',   '--obsName',  dest='OBSNAME',  type='string',default='TauB_Inc_0j_pTWgt',   help='Name of the observable, supported: "inclusive", "pT4l", "eta4l", "massZ2", "nJets"')
    #parser.add_option('',   '--obsBins',  dest='OBSBINS',  type='string',default='|105.0|160.0|',   help='Bin boundaries for the diff. measurement separated by "|", e.g. as "|0|50|100|", use the defalut if empty string')
    parser.add_option('',   '--obsBins',  dest='OBSBINS',  type='string',default='|105.0|1160.0|',   help='Bin boundaries for the diff. measurement separated by "|", e.g. as "|0|50|100|", use the defalut if empty string')
    #parser.add_option('',   '--obsBins',  dest='OBSBINS',  type='string',default='|105.0|160.0|',   help='Bin boundaries for the diff. measurement separated by "|", e.g. as "|0|50|100|", use the defalut if empty string')
    parser.add_option('',   '--year',  dest='YEAR',  type='string',default='2016',   help='Era to analyze, e.g. 2016, 2017, 2018 or Full ')
    parser.add_option('', '--bkg',      dest='BKG',type='string',default='', help='run with the type of zz background to float zz or qq_gg ')
    parser.add_option('',   '--debug',  dest='DEBUG',  type='int',default=0,   help='0 if debug false, else debug True')
    global opt, args
    (opt, args) = parser.parse_args()

def computeRatio(nominal, upVar, dnVar):
    if nominal == 0:
        return '-'
    else:
        dn_ratio = round(dnVar/nominal,3)
        up_ratio = round(upVar/nominal,3)
        if up_ratio==dn_ratio and up_ratio==1.000:
             return '-'
        elif up_ratio == dn_ratio:
            return str(dn_ratio)
        else:
            return str(dn_ratio)+'/'+str(up_ratio)
            #return str(up_ratio)+'/'+str(dn_ratio)



def extract_JES_nuis(nbins, obsName, obs_bins, DEBUG = 0):
#def extract_JES_nuis(x,nbins, obsName, obs_bins, DEBUG = 0):
    #x_points = [124, 126]
    x_points = [124,126,130]
    #x_points = [125]
#    x_points = [124]
    channels=['4mu','2e2mu','4e']
    #channels=['4e']
    #channels=['4mu']
#    modes=['ggH_powheg_JHUgen_','VBF_powheg_JHUgen_','ZH_powheg_JHUgen_','WH_powheg_JHUgen_','ttH_powheg_JHUgen_']
#    RootFile, Tree, nEvents, sumw = GrabMCTrees(opt.YEAR)
#    procs=['trueH','out_trueH','fakeH','bkg_qqzz','bkg_ggzz','bkg_zjets']
    #JES_nuis = ['Abs','Abs_year','BBEC1','BBEC1_year','EC2','EC2_year','FlavQCD','HF','HF_year','RelBal','RelSample_year', 'Total']
    #JES_nuis = ['Abs'] 
    #JES_nuis = ['FlavQCD'] 
#    procs = ['trueH','out_trueH','fakeH','bkg_qqzz','bkg_ggzz']
    #procs = ['trueH'] #,'out_trueH','fakeH','bkg_qqzz','bkg_ggzz']
   # procs = ['fakeH'] #,'out_trueH','fakeH','bkg_qqzz','bkg_ggzz']
    procs = ['trueH','out_trueH','fakeH','bkg_qqzz'] #,'bkg_ggzz']
    proc_selections = {"trueH":"passedFullSelection==1 ","out_trueH":"(passedFullSelection==1 && passedFiducialSelection!=1)","fakeH":"isH4l!=1","bkg_qqzz":"passedZ4lSelection==1","bkg_ggzz":"passedZ4lSelection==1","bkg_zjets":"(passedZXCRSelection==1 && nZXCRFailedLeptons==2)"}
    cut_m4l_reco = {"2e2mu":"(mass2e2mu>"+str(m4l_low)+" && mass2e2mu<"+str(m4l_high)+")","4mu":"(mass4mu>"+str(m4l_low)+" && mass4mu<"+str(m4l_high)+")","4e":"(mass4e>"+str(m4l_low)+" && mass4e<"+str(m4l_high)+")"}
#    recoweight = "genWeight*pileupWeight*dataMCWeight_new"

    weights = {"sig":"genWeight*pileupWeight*dataMCWeight_new","bkg_qqzz":"k_qqZZ_qcd_M*k_qqZZ_ewk","bkg_ggzz":"k_ggZZ","bkg_zjets":"1.0"}

    for x_point in x_points:
        for channel in channels:
            for proc in procs:
    	        RootFile = {}
    	        RootFile_list = []
		if ((proc=='trueH' or proc=='out_trueH' or proc=='fakeH')): 
		    recoweight = weights["sig"]
		else: 
		    recoweight = weights[proc]

		print "proc: ", proc, "    recoweight:    ",recoweight
                for i in range(0,len(SamplesMC[opt.YEAR])):
                    sample = SamplesMC[opt.YEAR][i].rstrip('.root')
		    if(sample.startswith("GluGluHToZZTo4")): continue # temporary, FIXME
		    print "sample is :   ", sample


                    if((not sample.startswith("ZZTo4L") and (not str(x_point) in sample))): continue
       	            if ((proc=='trueH' or proc=='out_trueH' or proc=='fakeH') and ((sample.startswith("ZZTo4L") or sample.startswith("GluGluToContin")))): continue
                    #if ((proc=='bkg_qqzz') and (("HToZZTo4L" in sample) or ("WH" in sample) or ("ZH" in sample) or ("ttH" in sample) or sample.startswith("GluGluToContin") )): continue
                    if ((proc=='bkg_qqzz') and (not sample.startswith("ZZTo4L") )): continue
                    if ((proc=='bkg_ggzz') and (not channel in sample)): continue
                    RootFile[sample] = TFile(dirMC[opt.YEAR]+'/'+sample+'.root',"READ")
        	    RootFile_list.append(dirMC[opt.YEAR]+'/'+sample+'.root')
                files = [ROOT.TFile(i) for i in RootFile_list]
		print "files are:  ", files

                trees = [i.Get("Ana/passedEvents") for i in files]
        	print "proc name:    ", proc 
        	print "RootFile_list:    ", RootFile_list
        	print "trees[0]		", trees[0]
		for recobin in range(len(obs_bins)-1):
		    obs_reco_low = obs_bins[recobin]
    		    obs_reco_high = obs_bins[recobin+1]
		    obs_gen_lowest = obs_bins[0]
    		    obs_gen_highest = obs_bins[len(obs_bins)-1]

		    cutobs_reco = "("+obs_reco+">="+str(obs_reco_low)+" && "+obs_reco+"<"+str(obs_reco_high)+")"
		    cut_nom = "("+recoweight+")*("+cutobs_reco+" && "+proc_selections[proc]+" && "+cut_m4l_reco[channel]+")"
		    #cut_nom = "("+weights["sig"]+")*("+cutobs_reco+" && "+proc_selections[proc]+" && "+cut_m4l_reco[channel]+")"

		    #processBin_nom = proc+'_'+obsName+'_'+channel+'_recobin'+str(recobin)
		    processBin_nom = proc+'_'+str(x_point)+'_'+obsName+'_'+channel+'_recobin'+str(recobin)
		    nom[processBin_nom] = 0.0; 
		
    		    Histos[processBin_nom] = TH1D(processBin_nom,processBin_nom, m4l_bins, m4l_low, m4l_high)
		    Histos[processBin_nom].Sumw2()

		    for nuis in JES_nuis:
		        processBin_jesup_nuis = proc+'_'+str(x_point)+'_'+obsName+'_'+channel+'_recobin'+str(recobin)+'_jesup_'+nuis
		        processBin_jesdn_nuis = proc+'_'+str(x_point)+'_'+obsName+'_'+channel+'_recobin'+str(recobin)+'_jesdn_'+nuis
		        processBin_ratio_nuis = proc+'_'+str(x_point)+'_'+obsName+'_'+channel+'_recobin'+str(recobin)+'_'+nuis

			Histos[processBin_jesup_nuis] = TH1D(processBin_jesup_nuis,processBin_jesup_nuis, m4l_bins, m4l_low, m4l_high)
			Histos[processBin_jesup_nuis].Sumw2()
			Histos[processBin_jesdn_nuis] = TH1D(processBin_jesdn_nuis,processBin_jesdn_nuis, m4l_bins, m4l_low, m4l_high)
			Histos[processBin_jesdn_nuis].Sumw2()

			#print "processBin_jesup_nuis:     ", processBin_jesup_nuis
			#print "processBin_jesdn_nuis:     ", processBin_jesdn_nuis

			cutobs_reco_nuis_up = "("+obs_reco+"_jesup_"+nuis+">="+str(obs_reco_low)+" && "+obs_reco+"_jesup_"+nuis+"<"+str(obs_reco_high)+")"
			cutobs_reco_nuis_dn = "("+obs_reco+"_jesdn_"+nuis+">="+str(obs_reco_low)+" && "+obs_reco+"_jesdn_"+nuis+"<"+str(obs_reco_high)+")"
			print "cutobs_reco_nuis_up:  " , cutobs_reco_nuis_up			
			print "cutobs_reco_nuis_dn:  " , cutobs_reco_nuis_dn			

			up[processBin_jesup_nuis] = 0.0; # Abs_year[processBin_jesup_nuis] = 0.0;
			dn[processBin_jesdn_nuis] = 0.0; # Abs_year[processBin_jesdn_nuis] = 0.0;
		        ratio[processBin_ratio_nuis] = 0.0;

			cut_nuis_up = "("+recoweight+")*("+cutobs_reco_nuis_up+" && "+proc_selections[proc]+" && "+cut_m4l_reco[channel]+")"
			cut_nuis_dn = "("+recoweight+")*("+cutobs_reco_nuis_dn+" && "+proc_selections[proc]+" && "+cut_m4l_reco[channel]+")"

			print "cut_nom :   ", cut_nom 
			print "cut_nuis_up :   ", cut_nuis_up
			print "cut_nuis_dn :   ", cut_nuis_dn 

			yield_nom = 0.;	yield_jesup_nuis = 0.; yield_jesdn_nuis = 0.;

			for tree in trees:
			    print "tree is:   ", tree

			    tree.Draw("mass4l >> "+processBin_nom,"("+cut_nom+")","goff")
			    tree.Draw("mass4l >> "+processBin_jesup_nuis,"("+cut_nuis_up+")","goff")
			    tree.Draw("mass4l >> "+processBin_jesdn_nuis,"("+cut_nuis_dn+")","goff")

			    yield_nom+=Histos[processBin_nom].Integral()
			    yield_jesup_nuis+=Histos[processBin_jesup_nuis].Integral()
			    yield_jesdn_nuis+=Histos[processBin_jesdn_nuis].Integral()

#			print "reading finished   ......."
#			print "yield_nom:   ", yield_nom
#			print "yield_jesup_nuis:   ", yield_jesup_nuis
#			print "yield_jesdn_nuis:   ", yield_jesdn_nuis
			
			nom[processBin_nom] = yield_nom 
			up[processBin_jesup_nuis] = yield_jesup_nuis 
			dn[processBin_jesdn_nuis] = yield_jesdn_nuis 
#		        print "ratio:    ", computeRatio(yield_nom,yield_jesup_nuis,yield_jesdn_nuis)	
			ratio[processBin_ratio_nuis] = computeRatio(yield_nom,yield_jesup_nuis,yield_jesdn_nuis)


    for channel in channels:
        for proc in procs:
	    if(not (proc=='trueH' or proc=='out_trueH' or proc=='fakeH')): continue
            for recobin in range(len(obs_bins)-1):
                for nuis in JES_nuis:
		    key_nom_125pt38 = proc+'_125.38_'+obsName+'_'+channel+'_recobin'+str(recobin)
                    key_jesup_125pt38 = proc+'_125.38_'+obsName+'_'+channel+'_recobin'+str(recobin)+'_jesup_'+nuis
                    key_jesdn_125pt38 = proc+'_125.38_'+obsName+'_'+channel+'_recobin'+str(recobin)+'_jesdn_'+nuis
		    key_ratio_nuis = proc+'_125.38_'+obsName+'_'+channel+'_recobin'+str(recobin)+'_'+nuis
    		    nom_points = []; jesUp_points = []; jesDn_points = [];
                    for x_point in x_points:
                        key_nom = proc+'_'+str(x_point)+'_'+obsName+'_'+channel+'_recobin'+str(recobin)
                        key_jesup = proc+'_'+str(x_point)+'_'+obsName+'_'+channel+'_recobin'+str(recobin)+'_jesup_'+nuis
                        key_jesdn = proc+'_'+str(x_point)+'_'+obsName+'_'+channel+'_recobin'+str(recobin)+'_jesdn_'+nuis
			nom_points.append(nom[key_nom]); jesUp_points.append(up[key_jesup]); jesDn_points.append(dn[key_jesdn]); 

		    print "nom_points:   ", nom_points
		    print "jesUp_points:   ", jesUp_points
		    print "jesDn_points:   ", jesDn_points			
		    spl_jesUp_points  = interpolate.splrep(x_points, jesUp_points, k=2)
		    spl_jesDn_points  = interpolate.splrep(x_points, jesDn_points, k=2)
		    spl_nom_points  = interpolate.splrep(x_points, nom_points, k=2)
		  
		    up[key_jesup_125pt38] = float(interpolate.splev(125.38, spl_jesUp_points))
		    dn[key_jesdn_125pt38] = float(interpolate.splev(125.38, spl_jesDn_points))
		    nom[key_nom_125pt38] = float(interpolate.splev(125.38, spl_nom_points))
		    ratio[key_ratio_nuis] = computeRatio(nom[key_nom_125pt38],up[key_jesup_125pt38],dn[key_jesdn_125pt38])

#    print "nom dictionary:  ", nom
#    print "up dictionary:  ",up
#    print "dn dictionary:  ",dn
    print "ratio dictionary:    ", ratio



if __name__ == "__main__":
    global opt, args
    parseOptions()
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

    m4l_low = 105.0
    m4l_high = 160.0
    m4l_bins = 35

    year = opt.YEAR
    print "year being processed: ",year

    Histos = {}
    nom = {}
    up = {}
    dn = {}
    ratio = {} 


    JES_nuis = ['Abs'] 
    print("Obs Name: {:15}  nBins: {:2}  bins: {}".format(opt.OBSNAME, nbins, observableBins))

    extract_JES_nuis(nbins, opt.OBSNAME, obs_bins, opt.DEBUG)
    #extract_JES_nuis(125.38, nbins, opt.OBSNAME, obs_bins, opt.DEBUG)

    ext=''
    with open('datacardInputs/inputs_JESnuis_'+opt.OBSNAME+'_'+year+'.py', 'w') as f:
        f.write('obsName = "'+str(opt.OBSNAME)+'" \n')
        f.write('obs_bins = '+str(obs_bins)+' \n')
        f.write('JES uncertainty sources = '+str(JES_nuis)+' \n')
        f.write('yield_nom = '+str(nom)+' \n')
        f.write('yield_up = '+str(up)+' \n')
        f.write('yield_dn = '+str(dn)+' \n')
        f.write('nuis_impact = '+str(ratio)+' \n')    


   


    print("Interpolation completed... :) ")
