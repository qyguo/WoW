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


#if (not os.path.exists("datacardInputs/"+opt.YEAR+"/"+obsName)):
#    #os.system("mkdir -p datacardInputs")
#    os.system("mkdir -p datacardInputs/"+opt.YEAR+"/"+obsname+"")


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
    parser.add_option('',   '--year',  dest='YEAR',  type='string',default='2017',   help='Era to analyze, e.g. 2016, 2017, 2018 or Full ')
    #parser.add_option('',   '--era',  dest='ERA',  type='string',default='2018',   help='Era to analyze, e.g. 2016, 2017, 2018 or Full ')
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
        elif dn_ratio>up_ratio and up_ratio>=1.000: #  
	    up_ratio= (1 - (dn_ratio - 1));
	    return str(dn_ratio)+'/'+str(up_ratio)
        elif up_ratio>dn_ratio and dn_ratio>=1.000: #  
	    dn_ratio= (1 - (up_ratio - 1));
	    return str(dn_ratio)+'/'+str(up_ratio)
	elif dn_ratio<1.000 and up_ratio<1.000:
	    dn_ratio=min(dn_ratio,up_ratio)
	    up_ratio=1+(1-dn_ratio)
	    return str(dn_ratio)+'/'+str(up_ratio)
        else:
            return str(dn_ratio)+'/'+str(up_ratio)
            #return str(up_ratio)+'/'+str(dn_ratio)



def extract_JES_nuis(nbins, obsName, obs_bins, DEBUG = 0):
    #x_points = [124,125,126]
    x_points = [125]
#    x_points = [124]
    channels=['4mu','2e2mu','4e']
    #channels=['2e2mu']
    #channels=['4mu']
    #procs=['trueH','out_trueH','fakeH','bkg_qqzz','bkg_ggzz','bkg_zjets']
    procs=['bkg_zjets','trueH','out_trueH','fakeH','bkg_qqzz','bkg_ggzz']
    #procs = ['bkg_qqzz']
    #procs = ['trueH']

    #proc_selections = {"trueH":"passedFullSelection==1 ","out_trueH":"(passedFullSelection==1 && passedFiducialSelection!=1)","fakeH":"isH4l!=1","bkg_qqzz":"passedZ4lSelection==1","bkg_ggzz":"passedZ4lSelection==1","bkg_zjets":"(passedZXCRSelection==1 && nZXCRFailedLeptons==2)"}

    #proc_selections = {"trueH":"passedFullSelection==1 ","out_trueH":"(passedFullSelection==1 && passedFiducialSelection!=1)","fakeH":"isH4l!=1","bkg_qqzz":"passedZ4lSelection==1","bkg_ggzz":"passedZ4lSelection==1","bkg_zjets":"(passedZXCRSelection==1)"}
    proc_selections = {"trueH":"passedFullSelection==1 ","out_trueH":"(passedFullSelection==1 && passedFiducialSelection!=1)","fakeH":"isH4l!=1","bkg_qqzz":"passedFullSelection==1","bkg_ggzz":"passedFullSelection==1","bkg_zjets":"(passedZXCRSelection==1)"}


    #proc_selections = {"trueH":"passedFullSelection==1 ","out_trueH":"(passedFullSelection==1 && passedFiducialSelection!=1)","fakeH":"isH4l!=1","bkg_qqzz":"(passedFullSelection!=1 && passedZ4lSelection==1)","bkg_ggzz":"(passedFullSelection!=1 && passedZ4lSelection==1)","bkg_zjets":"(passedFullSelection!=1 && passedZXCRSelection==1)"}  # test
    #proc_selections = {"trueH":"passedFullSelection==1 ","out_trueH":"(passedFullSelection==1 && passedFiducialSelection!=1)","fakeH":"isH4l!=1","bkg_qqzz":"(passedFullSelection!=1)","bkg_ggzz":"(passedFullSelection!=1)","bkg_zjets":"(passedFullSelection!=1)"}
    cut_m4l_reco = {"2e2mu":"(mass2e2mu>="+str(m4l_low)+" && mass2e2mu<="+str(m4l_high)+")","4mu":"(mass4mu>="+str(m4l_low)+" && mass4mu<="+str(m4l_high)+")","4e":"(mass4e>="+str(m4l_low)+" && mass4e<="+str(m4l_high)+")"}

    cut_zerobin = {"pt_leadingjet_pt30_eta4p7":"==0", "pTj2":"<=1", "mj1j2":"<=1", "pT4lj":"==0","mass4lj":"==0","dPhiHj1":"==0","dyHj1":"==0","TauB_Inc_0j_pTWgt":"==0","TauC_Inc_0j_EnergyWgt":"==0","pT4ljj":"<=1","mass4ljj":"<=1","dEtaj1j2":"<=1","dPhij1j2":"<=1","dPhiHj1j2":"<=1"}

    #lumi = {"2016":35867.0,"2017":41370.0,"2018":58800.0}
    #lumi = {"2016":35900.0,"2017":41500.0,"2018":59700.0}
    lumi = {"2016":35900.0,"2017":41480.0,"2018":59830.0}  # https://twiki.cern.ch/twiki/bin/viewauth/CMS/PdmVRun2LegacyAnalysis
    #lumi_part1=5.9+2.6 +2.3 +1.8 + 3.2;  #sum of luminosities of Run BCDEF 2016
    #lumi_preVFP=5.444+2.399 +4.256 +4.054 + 2.691;  #sum of luminosities of Run BCDEF 2016 (preVFP) # from processed lumi
    #lumi_postVFP=0.414 + 7.544 + 8.746;  # sum of luminosities of Run FGH (postVFP)
    lumi_preVFP=19520; lumi_postVFP=16810;  # https://twiki.cern.ch/twiki/bin/viewauth/CMS/PdmVRun2LegacyAnalysis

#    lumiWt1 = (lumi_part1)/(lumi_part1+lumi_part2)
#    lumiWt2 = (lumi_part2)/(lumi_part1+lumi_part2)
    #return((lumi_part1*sf1+ lumi_part2*sf2)/(lumi_part1+lumi_part2),sqrt( (lumiWt1 * sfe1)**2  + (lumiWt2 * sfe2)**2 ))



    #weights = {"sig":"genWeight*pileupWeight*prefiringWeight*dataMCWeight_new*crossSection*"+str(lumi[opt.YEAR]),"bkg_qqzz":"(genWeight*pileupWeight*prefiringWeight*dataMCWeight_new*crossSection)*(k_qqZZ_qcd_M*k_qqZZ_ewk)*"+str(lumi[opt.YEAR]),"bkg_ggzz":"(genWeight*pileupWeight*prefiringWeight*dataMCWeight_new*crossSection)*(k_ggZZ)*"+str(lumi[opt.YEAR]),"bkg_zjets":"1.0"}
    weights = {"sig":"genWeight*pileupWeight*prefiringWeight*dataMCWeight_new*crossSection","bkg_qqzz":"(genWeight*pileupWeight*prefiringWeight*dataMCWeight_new*crossSection)*(k_qqZZ_qcd_M*k_qqZZ_ewk)","bkg_ggzz":"(genWeight*pileupWeight*prefiringWeight*dataMCWeight_new*crossSection)*(k_ggZZ)","bkg_zjets":"1.0"}
#    int stop=0;
#    sumw = {};sumw_pre = {}; sumw_post = {} 
    for x_point in x_points:
        for channel in channels:
            for proc in procs:
    	        RootFile = {}
    	        RootFile_list = []
		RootFile_list_pre = []
		RootFile_list_post = []

		if ((proc=='trueH' or proc=='out_trueH' or proc=='fakeH')): 
		    recoweight = weights["sig"]
		else: 
		    recoweight = weights[proc]

		print "recoweight:    ",recoweight
		if (proc=='bkg_zjets'): samples = SamplesData[opt.YEAR]
		else: samples =SamplesMC[opt.YEAR]


		sumw2 = []; sumw2_pre = []; sumw2_post = [];
                for i in range(0,len(samples)):
                    sample = samples[i].rstrip('.root')

                    if((proc!='bkg_zjets' and proc!='bkg_ggzz') and (not sample.startswith('ZZTo4L') and (not str(x_point) in sample))): continue
                    if((proc=='bkg_qqzz') and (not sample.startswith('ZZTo4L'))): continue # and (not str(x_point) in sample))): continue
                    if((proc=='bkg_ggzz' and proc!='bkg_zjets' and proc!='bkg_qqzz' and (proc!='trueH' or proc!='out_trueH' or proc!='fakeH')) and (not sample.startswith('GluGluToContin'))): continue # and (not str(x_point) in sample))): continue
                    if((proc=='bkg_ggzz' and proc!='bkg_zjets' and proc!='bkg_qqzz' and (proc!='trueH' or proc!='out_trueH' or proc!='fakeH')) and (not channel in sample)): continue # and (not str(x_point) in sample))): continue
                    if ((proc!='bkg_zjets' and proc!='bkg_qqzz' and proc!='bkg_ggzz') and (proc=='trueH' or proc=='out_trueH' or proc=='fakeH') and ((sample.startswith('ZZTo4L') or sample.startswith('GluGluToContin')))): continue

#                    RootFile[sample] = TFile(dirMC[opt.YEAR]+'/'+sample+'.root',"READ")  # currently, all files in the same directory
#        	    RootFile_list.append(dirMC[opt.YEAR]+'/'+sample+'.root')
		    if (opt.YEAR=='2018' or opt.YEAR=='2017'):
                        RootFile[sample] = TFile(dirMC[opt.YEAR]+'/'+sample+'.root',"READ")  # currently, all files in the same directory
           	        RootFile_list.append(dirMC[opt.YEAR]+'/'+sample+'.root')
		    else:
		    #elif (opt.YEAR=='2016'):
		        RootFile_pre[sample] = TFile(dirMC['2016preVFP']+'/'+sample+'.root',"READ")
			RootFile_list_pre.append(dirMC['2016preVFP']+'/'+sample+'.root')   
		        RootFile_post[sample] = TFile(dirMC['2016postVFP']+'/'+sample+'.root',"READ")
			RootFile_list_post.append(dirMC['2016postVFP']+'/'+sample+'.root')   

		    print "sample is :   ", sample
		    #if (proc!='bkg_zjets'): sumw2.append(sumw[sample])
		    if (proc!='bkg_zjets' and (opt.YEAR=='2018' or opt.YEAR=='2017')): 
			sumw2.append(sumw[sample])  # from keys to list, for onward use
		    elif (proc!='bkg_zjets' and (opt.YEAR=='2016')): 
			sumw2_pre.append(sumw_pre[sample]); sumw2_post.append(sumw_post[sample])
#                files = [ROOT.TFile(i) for i in RootFile_list]
		if (opt.YEAR=='2018' or opt.YEAR=='2017'):
                    files = [ROOT.TFile(i) for i in RootFile_list]
		elif (opt.YEAR=='2016'):
		    files_pre = [ROOT.TFile(i) for i in RootFile_list_pre]
		    files_post = [ROOT.TFile(i) for i in RootFile_list_post]
		#    print "proc:  ", proc, "  channel:   ", channel
		    print "files_pre are:  ", files_pre
		    print "files_post are:  ", files_post
		print "proc:  ", proc, "  channel:   ", channel

		if (opt.YEAR=='2018' or opt.YEAR=='2017'):
		    if(proc=='bkg_zjets'): trees = [i.Get("passedEvents") for i in files]
                    else: trees = [i.Get("Ana/passedEvents") for i in files]
		elif (opt.YEAR=='2016'):
		    if(proc=='bkg_zjets'): trees_pre = [i.Get("passedEvents") for i in files_pre]; trees_post = [j.Get("passedEvents") for j in files_post];
                    else: trees_pre = [i.Get("Ana/passedEvents") for i in files_pre]; trees_post = [j.Get("Ana/passedEvents") for j in files_post];


		for recobin in range(len(obs_bins)-1):
		    obs_reco_low = obs_bins[recobin]
    		    obs_reco_high = obs_bins[recobin+1]
		    obs_gen_lowest = obs_bins[0]
    		    obs_gen_highest = obs_bins[len(obs_bins)-1]

		    # hack for bin0
		    if (recobin==0 and (not "Tau" in obs_reco and not "njets" in obs_reco)): 
		    #if (recobin==0 and (not "njets" in obs_reco)): 
			cutobs_reco = "njets_pt30_eta4p7"+cut_zerobin[obs_reco]
			#if ("Tau" in obs_reco):
			 #   cutobs_reco = "("+cutobs_reco+") && " + "("+obs_reco+">="+str(obs_reco_low)+" && "+obs_reco+"<"+str(obs_reco_high)+")" 
#			if ("Tau" in obs_reco):
#			    cutobs_reco = "("+cutobs_reco+") || " + "("+obs_reco+">="+str(obs_reco_low)+" && "+obs_reco+"<"+str(obs_reco_high)+")" 
			print "zero bin, obs name:  ", cutobs_reco
		    else:
		        cutobs_reco = "("+obs_reco+">="+str(obs_reco_low)+" && "+obs_reco+"<"+str(obs_reco_high)+")"
		    cut_nom = "("+recoweight+")*("+cutobs_reco+" && "+proc_selections[proc]+" && "+cut_m4l_reco[channel]+")"

		    #processBin_nom = proc+'_'+obsName+'_'+channel+'_recobin'+str(recobin)

		    if (not (proc=='trueH' or proc=='out_trueH' or proc=='fakeH')): x_point=125 # temporary fix for bkg

		    processBin_nom = proc+'_'+str(x_point)+'_'+obsName+'_'+channel+'_recobin'+str(recobin)
		    #processBin_nom_Hproc = 'Hproc_'+str(x_point)+'_'+obsName+'_'+channel+'_recobin'+str(recobin)
		    nom[processBin_nom] = 0.0; 
		    #nom[processBin_nom_Hproc] = 0.0;  # defined to sum-up the Higgs involved processes (trueH, fakeH, out_trueH) 
		
		    if (opt.YEAR=='2018' or opt.YEAR=='2017'):
    		        Histos[processBin_nom] = TH1D(processBin_nom,processBin_nom, m4l_bins, m4l_low, m4l_high)
		        Histos[processBin_nom].Sumw2()
		    elif (opt.YEAR=='2016'):
			Histos[processBin_nom+'pre'] = TH1D(processBin_nom+'pre',processBin_nom+'pre', m4l_bins, m4l_low, m4l_high)
			Histos[processBin_nom+'pre'].Sumw2()
			Histos[processBin_nom+'post'] = TH1D(processBin_nom+'post',processBin_nom+'post', m4l_bins, m4l_low, m4l_high)
			Histos[processBin_nom+'post'].Sumw2()
    		    #Histos[processBin_nom_Hproc] = TH1D(processBin_nom_Hproc,processBin_nom_Hproc, m4l_bins, m4l_low, m4l_high)
		    #Histos[processBin_nom_Hproc].Sumw2()

		    for nuis in JES_nuis:

		        processBin_jesup_nuis = proc+'_'+str(x_point)+'_'+obsName+'_'+channel+'_recobin'+str(recobin)+'_jesup_'+nuis
		        processBin_jesdn_nuis = proc+'_'+str(x_point)+'_'+obsName+'_'+channel+'_recobin'+str(recobin)+'_jesdn_'+nuis
		        processBin_ratio_nuis = proc+'_'+str(x_point)+'_'+obsName+'_'+channel+'_recobin'+str(recobin)+'_'+nuis

			if (opt.YEAR=='2018' or opt.YEAR=='2017'):
			    Histos[processBin_jesup_nuis] = TH1D(processBin_jesup_nuis,processBin_jesup_nuis, m4l_bins, m4l_low, m4l_high)
			    Histos[processBin_jesup_nuis].Sumw2()
			    Histos[processBin_jesdn_nuis] = TH1D(processBin_jesdn_nuis,processBin_jesdn_nuis, m4l_bins, m4l_low, m4l_high)
			    Histos[processBin_jesdn_nuis].Sumw2()
			elif (opt.YEAR=='2016'):
			    Histos[processBin_jesup_nuis+'pre'] = TH1D(processBin_jesup_nuis+'pre',processBin_jesup_nuis+'pre', m4l_bins, m4l_low, m4l_high)
                            Histos[processBin_jesup_nuis+'pre'].Sumw2()
                            Histos[processBin_jesdn_nuis+'pre'] = TH1D(processBin_jesdn_nuis+'pre',processBin_jesdn_nuis+'pre', m4l_bins, m4l_low, m4l_high)
                            Histos[processBin_jesdn_nuis+'pre'].Sumw2()

			    Histos[processBin_jesup_nuis+'post'] = TH1D(processBin_jesup_nuis+'post',processBin_jesup_nuis+'post', m4l_bins, m4l_low, m4l_high)
                            Histos[processBin_jesup_nuis+'post'].Sumw2()
                            Histos[processBin_jesdn_nuis+'post'] = TH1D(processBin_jesdn_nuis+'post',processBin_jesdn_nuis+'post', m4l_bins, m4l_low, m4l_high)
                            Histos[processBin_jesdn_nuis+'post'].Sumw2()


                    # hack for bin0
		        if (recobin==0 and (not "Tau" in obs_reco and not "njets" in obs_reco)): 
		        #if (recobin==0 and (not "njets" in obs_reco)): 
			    cutobs_reco_nuis_up = "njets_pt30_eta4p7_jesup_"+nuis+cut_zerobin[obs_reco]
			    cutobs_reco_nuis_dn = "njets_pt30_eta4p7_jesdn_"+nuis+cut_zerobin[obs_reco]
#                            if ("Tau" in obs_reco):
#                                #cutobs_reco_nuis_up = "("+cutobs_reco_nuis_up+") && " + "("+obs_reco+"_jesup_"+nuis+">="+str(obs_reco_low)+" && "+obs_reco+"_jesup_"+nuis+"<"+str(obs_reco_high)+")"
#                                #cutobs_reco_nuis_dn = "("+cutobs_reco_nuis_dn+") && " + "("+obs_reco+"_jesdn_"+nuis+">="+str(obs_reco_low)+" && "+obs_reco+"_jesdn_"+nuis+"<"+str(obs_reco_high)+")"
#                                cutobs_reco_nuis_up = "("+cutobs_reco_nuis_up+") || " + "("+obs_reco+"_jesup_"+nuis+">="+str(obs_reco_low)+" && "+obs_reco+"_jesup_"+nuis+"<"+str(obs_reco_high)+")"
#                                cutobs_reco_nuis_dn = "("+cutobs_reco_nuis_dn+") || " + "("+obs_reco+"_jesdn_"+nuis+">="+str(obs_reco_low)+" && "+obs_reco+"_jesdn_"+nuis+"<"+str(obs_reco_high)+")"
			else:
			    cutobs_reco_nuis_up = "("+obs_reco+"_jesup_"+nuis+">="+str(obs_reco_low)+" && "+obs_reco+"_jesup_"+nuis+"<"+str(obs_reco_high)+")"
			    cutobs_reco_nuis_dn = "("+obs_reco+"_jesdn_"+nuis+">="+str(obs_reco_low)+" && "+obs_reco+"_jesdn_"+nuis+"<"+str(obs_reco_high)+")"

			print "cutobs_reco_nuis_up:  " , cutobs_reco_nuis_up			
			print "cutobs_reco_nuis_dn:  " , cutobs_reco_nuis_dn			

		        ratio[processBin_ratio_nuis] = 0.0;
			if (nuis=="Abs"): up_Abs[processBin_jesup_nuis]=0.0; dn_Abs[processBin_jesdn_nuis]=0.0;ratio_Abs[processBin_ratio_nuis]=0.0;
			if (nuis=="Abs_year"): up_Abs_year[processBin_jesup_nuis]=0.0; dn_Abs_year[processBin_jesdn_nuis]=0.0;ratio_Abs_year[processBin_ratio_nuis]=0.0;
			if (nuis=="BBEC1"): up_BBEC1[processBin_jesup_nuis]=0.0; dn_BBEC1[processBin_jesdn_nuis]=0.0;ratio_BBEC1[processBin_ratio_nuis]=0.0;
			if (nuis=="BBEC1_year"): up_BBEC1_year[processBin_jesup_nuis]=0.0; dn_BBEC1_year[processBin_jesdn_nuis]=0.0;ratio_BBEC1_year[processBin_ratio_nuis]=0.0;
			if (nuis=="EC2"): up_EC2[processBin_jesup_nuis]=0.0; dn_EC2[processBin_jesdn_nuis]=0.0;ratio_EC2[processBin_ratio_nuis]=0.0;
			if (nuis=="EC2_year"): up_EC2_year[processBin_jesup_nuis]=0.0; dn_EC2_year[processBin_jesdn_nuis]=0.0;ratio_EC2_year[processBin_ratio_nuis]=0.0;
			if (nuis=="FlavQCD"): up_FlavQCD[processBin_jesup_nuis]=0.0; dn_FlavQCD[processBin_jesdn_nuis]=0.0;ratio_FlavQCD[processBin_ratio_nuis]=0.0;
			if (nuis=="HF"): up_HF[processBin_jesup_nuis]=0.0; dn_HF[processBin_jesdn_nuis]=0.0;ratio_HF[processBin_ratio_nuis]=0.0;
			if (nuis=="HF_year"): up_HF_year[processBin_jesup_nuis]=0.0; dn_HF_year[processBin_jesdn_nuis]=0.0;ratio_HF_year[processBin_ratio_nuis]=0.0;
			if (nuis=="RelBal"): up_RelBal[processBin_jesup_nuis]=0.0; dn_RelBal[processBin_jesdn_nuis]=0.0;ratio_RelBal[processBin_ratio_nuis]=0.0;
			if (nuis=="RelSample_year"): up_RelSample_year[processBin_jesup_nuis]=0.0; dn_RelSample_year[processBin_jesdn_nuis]=0.0;ratio_RelSample_year[processBin_ratio_nuis]=0.0;
			if (nuis=="Total"): up_Total[processBin_jesup_nuis]=0.0; dn_Total[processBin_jesdn_nuis]=0.0;ratio_Total[processBin_ratio_nuis]=0.0;


			cut_nuis_up = "("+recoweight+")*("+cutobs_reco_nuis_up+" && "+proc_selections[proc]+" && "+cut_m4l_reco[channel]+")"
			cut_nuis_dn = "("+recoweight+")*("+cutobs_reco_nuis_dn+" && "+proc_selections[proc]+" && "+cut_m4l_reco[channel]+")"

			print "cut_nom :   ", cut_nom 
			print "cut_nuis_up :   ", cut_nuis_up
			print "cut_nuis_dn :   ", cut_nuis_dn 

			yield_nom = 0.;	yield_jesup_nuis = 0.; yield_jesdn_nuis = 0.;
		#	yield_nom_Hproc=0;
			i=0;  # iterator to read sumw2 from corresponding root file, meant to work for signal proc. only
	
			if (opt.YEAR=='2018' or opt.YEAR=='2017'):
			    for tree in trees:
			        if(proc=='bkg_zjets'):
                                    tree.Draw("mass4l >> "+processBin_nom,"("+cut_nom+")","goff")
                                    tree.Draw("mass4l >> "+processBin_jesup_nuis,"("+cut_nuis_up+")","goff")
                                    tree.Draw("mass4l >> "+processBin_jesdn_nuis,"("+cut_nuis_dn+")","goff")				
			        else:
			            print "proc is:   ", proc
			            print "channel is:   ", channel
			            print "i:  ",i
				    print "len(sumw2):    ", len(sumw2)
			            print "sumw2[",i,"] : ", sumw2[i] 
			            tree.Draw("mass4l >> "+processBin_nom,"("+cut_nom+"*("+str(lumi[opt.YEAR])+")/"+str(sumw2[i])+")","goff")
			            tree.Draw("mass4l >> "+processBin_jesup_nuis,"("+cut_nuis_up+"*("+str(lumi[opt.YEAR])+")/"+str(sumw2[i])+")","goff")
			            tree.Draw("mass4l >> "+processBin_jesdn_nuis,"("+cut_nuis_dn+"*("+str(lumi[opt.YEAR])+")/"+str(sumw2[i])+")","goff")
			        i=i+1

			        yield_nom+=Histos[processBin_nom].Integral()
			        yield_jesup_nuis+=Histos[processBin_jesup_nuis].Integral()
			        yield_jesdn_nuis+=Histos[processBin_jesdn_nuis].Integral()

                        else:
			    for pre_tree, post_tree in zip(trees_pre, trees_post):    
				print "pre_tree is:  ", pre_tree
				print "post_tree is:  ", post_tree

				if(proc=='bkg_zjets'):
                                    pre_tree.Draw("mass4l >> "+processBin_nom+'pre',"("+cut_nom+")","goff")
                                    pre_tree.Draw("mass4l >> "+processBin_jesup_nuis+'pre',"("+cut_nuis_up+")","goff")
                                    pre_tree.Draw("mass4l >> "+processBin_jesdn_nuis+'pre',"("+cut_nuis_dn+")","goff")

                                    post_tree.Draw("mass4l >> "+processBin_nom+'post',"("+cut_nom+")","goff")
                                    post_tree.Draw("mass4l >> "+processBin_jesup_nuis+'post',"("+cut_nuis_up+")","goff")
                                    post_tree.Draw("mass4l >> "+processBin_jesdn_nuis+'post',"("+cut_nuis_dn+")","goff")
                                else:

                                    print "proc is:   ", proc
                                    print "channel is:   ", channel
                                    print "i:  ",i
                                    print "len(sumw2_pre):    ", len(sumw2_pre)
                                    print "sumw2_pre[",i,"] : ", sumw2_pre[i]
                                    print "len(sumw2_post):    ", len(sumw2_post)
                                    print "sumw2_post[",i,"] : ", sumw2_post[i]

                                    #pre_tree.Draw("mass4l >> "+processBin_nom+'pre',"("+cut_nom+"/"+str(sumw2_pre[i])+")","goff")
                                    pre_tree.Draw("mass4l >> "+processBin_nom+'pre',"("+cut_nom+"*("+str(lumi_preVFP)+")/"+str(sumw2_pre[i])+")","goff")
                                    #pre_tree.Draw("mass4l >> "+processBin_jesup_nuis+'pre',"("+cut_nuis_up+"/"+str(sumw2_pre[i])+")","goff")
                                    pre_tree.Draw("mass4l >> "+processBin_jesup_nuis+'pre',"("+cut_nuis_up+"*("+str(lumi_preVFP)+")/"+str(sumw2_pre[i])+")","goff")
                                    pre_tree.Draw("mass4l >> "+processBin_jesdn_nuis+'pre',"("+cut_nuis_dn+"*("+str(lumi_preVFP)+")/"+str(sumw2_pre[i])+")","goff")

                                    post_tree.Draw("mass4l >> "+processBin_nom+'post',"("+cut_nom+"*("+str(lumi_postVFP)+")/"+str(sumw2_post[i])+")","goff")
                                    post_tree.Draw("mass4l >> "+processBin_jesup_nuis+'post',"("+cut_nuis_up+"*("+str(lumi_postVFP)+")/"+str(sumw2_post[i])+")","goff")
                                    post_tree.Draw("mass4l >> "+processBin_jesdn_nuis+'post',"("+cut_nuis_dn+"*("+str(lumi_postVFP)+")/"+str(sumw2_post[i])+")","goff")

                                i=i+1

				print "pre Integral", Histos[processBin_nom+'pre'].Integral()
				print "post Integral", Histos[processBin_nom+'post'].Integral()

				
				yield_nom+= Histos[processBin_nom+'pre'].Integral() + Histos[processBin_nom+'post'].Integral()
				yield_jesup_nuis+= Histos[processBin_jesup_nuis+'pre'].Integral() + Histos[processBin_jesup_nuis+'post'].Integral()
				yield_jesdn_nuis+= Histos[processBin_jesdn_nuis+'pre'].Integral() + Histos[processBin_jesdn_nuis+'post'].Integral()


##       		
			nom[processBin_nom] = yield_nom 
			if (nuis=="Abs"): up_Abs[processBin_jesup_nuis] = yield_jesup_nuis; dn_Abs[processBin_jesdn_nuis] = yield_jesdn_nuis; ratio_Abs[processBin_ratio_nuis] = computeRatio(yield_nom,yield_jesup_nuis,yield_jesdn_nuis)
			if (nuis=="Abs_year"): up_Abs_year[processBin_jesup_nuis] = yield_jesup_nuis; dn_Abs_year[processBin_jesdn_nuis] = yield_jesdn_nuis; ratio_Abs_year[processBin_ratio_nuis] = computeRatio(yield_nom,yield_jesup_nuis,yield_jesdn_nuis)
			if (nuis=="BBEC1"): up_BBEC1[processBin_jesup_nuis] = yield_jesup_nuis; dn_BBEC1[processBin_jesdn_nuis] = yield_jesdn_nuis; ratio_BBEC1[processBin_ratio_nuis] = computeRatio(yield_nom,yield_jesup_nuis,yield_jesdn_nuis)
			if (nuis=="BBEC1_year"): up_BBEC1_year[processBin_jesup_nuis] = yield_jesup_nuis; dn_BBEC1_year[processBin_jesdn_nuis] = yield_jesdn_nuis; ratio_BBEC1_year[processBin_ratio_nuis] = computeRatio(yield_nom,yield_jesup_nuis,yield_jesdn_nuis)
			if (nuis=="EC2"): up_EC2[processBin_jesup_nuis] = yield_jesup_nuis; dn_EC2[processBin_jesdn_nuis] = yield_jesdn_nuis; ratio_EC2[processBin_ratio_nuis] = computeRatio(yield_nom,yield_jesup_nuis,yield_jesdn_nuis)
			if (nuis=="EC2_year"): up_EC2_year[processBin_jesup_nuis] = yield_jesup_nuis; dn_EC2_year[processBin_jesdn_nuis] = yield_jesdn_nuis; ratio_EC2_year[processBin_ratio_nuis] = computeRatio(yield_nom,yield_jesup_nuis,yield_jesdn_nuis)
			if (nuis=="FlavQCD"): up_FlavQCD[processBin_jesup_nuis] = yield_jesup_nuis; dn_FlavQCD[processBin_jesdn_nuis] = yield_jesdn_nuis; ratio_FlavQCD[processBin_ratio_nuis] = computeRatio(yield_nom,yield_jesup_nuis,yield_jesdn_nuis)
			if (nuis=="HF"): up_HF[processBin_jesup_nuis] = yield_jesup_nuis; dn_HF[processBin_jesdn_nuis] = yield_jesdn_nuis; ratio_HF[processBin_ratio_nuis] = computeRatio(yield_nom,yield_jesup_nuis,yield_jesdn_nuis)
			if (nuis=="HF_year"): up_HF_year[processBin_jesup_nuis] = yield_jesup_nuis; dn_HF_year[processBin_jesdn_nuis] = yield_jesdn_nuis; ratio_HF_year[processBin_ratio_nuis] = computeRatio(yield_nom,yield_jesup_nuis,yield_jesdn_nuis)
			if (nuis=="RelBal"): up_RelBal[processBin_jesup_nuis] = yield_jesup_nuis; dn_RelBal[processBin_jesdn_nuis] = yield_jesdn_nuis; ratio_RelBal[processBin_ratio_nuis] = computeRatio(yield_nom,yield_jesup_nuis,yield_jesdn_nuis)
			if (nuis=="RelSample_year"): up_RelSample_year[processBin_jesup_nuis] = yield_jesup_nuis; dn_RelSample_year[processBin_jesdn_nuis] = yield_jesdn_nuis; ratio_RelSample_year[processBin_ratio_nuis] = computeRatio(yield_nom,yield_jesup_nuis,yield_jesdn_nuis)
			if (nuis=="Total"): up_Total[processBin_jesup_nuis] = yield_jesup_nuis; dn_Total[processBin_jesdn_nuis] = yield_jesdn_nuis; ratio_Total[processBin_ratio_nuis] = computeRatio(yield_nom,yield_jesup_nuis,yield_jesdn_nuis)


## combing the H processes
    for x_point in x_points:
        for channel in channels:
            for recobin in range(len(obs_bins)-1):
		processBin_nom_Hproc = 'Hproc_'+str(x_point)+'_'+obsName+'_'+channel+'_recobin'+str(recobin)
		nom[processBin_nom_Hproc]=0;
		yield_nom_Hproc=0;
       	        for proc in procs:
                    if(not (proc=='trueH' or proc=='out_trueH' or proc=='fakeH')): continue
		    processBin_nom = proc+'_'+str(x_point)+'_'+obsName+'_'+channel+'_recobin'+str(recobin)
		    yield_nom_Hproc+=nom[processBin_nom]
		    #nom[processBin_nom_Hproc]=0;
#		    print "processBin_nom  ", processBin_nom
#		    print "nom[processBin_nom]  ", nom[processBin_nom]
#		    print "yield_nom_Hproc   ", yield_nom_Hproc
    	        nom[processBin_nom_Hproc]=yield_nom_Hproc;    
#		print "processBin_nom_Hproc  ", processBin_nom_Hproc
#		print "nom[processBin_nom_Hproc]  ", nom[processBin_nom_Hproc] 


## JES nuis
  	        #yield_jesup_nuis_Hproc=0; yield_jesdn_nuis_Hproc=0;
		for nuis in JES_nuis:
                    processBin_jesup_nuis_Hproc = 'Hproc_'+str(x_point)+'_'+obsName+'_'+channel+'_recobin'+str(recobin)+'_jesup_'+nuis
                    processBin_jesdn_nuis_Hproc = 'Hproc_'+str(x_point)+'_'+obsName+'_'+channel+'_recobin'+str(recobin)+'_jesdn_'+nuis
                    processBin_ratio_nuis_Hproc = 'Hproc_'+str(x_point)+'_'+obsName+'_'+channel+'_recobin'+str(recobin)+'_'+nuis

            	    ratio[processBin_ratio_nuis_Hproc] = 0.0;
                    if (nuis=="Abs"): up_Abs[processBin_jesup_nuis_Hproc]=0.0; dn_Abs[processBin_jesdn_nuis_Hproc]=0.0;ratio_Abs[processBin_ratio_nuis_Hproc]=0.0;
                    if (nuis=="Abs_year"): up_Abs_year[processBin_jesup_nuis_Hproc]=0.0; dn_Abs_year[processBin_jesdn_nuis_Hproc]=0.0;ratio_Abs_year[processBin_ratio_nuis_Hproc]=0.0;
                    if (nuis=="BBEC1"): up_BBEC1[processBin_jesup_nuis_Hproc]=0.0; dn_BBEC1[processBin_jesdn_nuis_Hproc]=0.0;ratio_BBEC1[processBin_ratio_nuis_Hproc]=0.0;
                    if (nuis=="BBEC1_year"): up_BBEC1_year[processBin_jesup_nuis_Hproc]=0.0; dn_BBEC1_year[processBin_jesdn_nuis_Hproc]=0.0;ratio_BBEC1_year[processBin_ratio_nuis_Hproc]=0.0;
                    if (nuis=="EC2"): up_EC2[processBin_jesup_nuis_Hproc]=0.0; dn_EC2[processBin_jesdn_nuis_Hproc]=0.0;ratio_EC2[processBin_ratio_nuis_Hproc]=0.0;
                    if (nuis=="EC2_year"): up_EC2_year[processBin_jesup_nuis_Hproc]=0.0; dn_EC2_year[processBin_jesdn_nuis_Hproc]=0.0;ratio_EC2_year[processBin_ratio_nuis_Hproc]=0.0;
                    if (nuis=="FlavQCD"): up_FlavQCD[processBin_jesup_nuis_Hproc]=0.0; dn_FlavQCD[processBin_jesdn_nuis_Hproc]=0.0;ratio_FlavQCD[processBin_ratio_nuis_Hproc]=0.0;
                    if (nuis=="HF"): up_HF[processBin_jesup_nuis_Hproc]=0.0; dn_HF[processBin_jesdn_nuis_Hproc]=0.0;ratio_HF[processBin_ratio_nuis_Hproc]=0.0;
                    if (nuis=="HF_year"): up_HF_year[processBin_jesup_nuis_Hproc]=0.0; dn_HF_year[processBin_jesdn_nuis_Hproc]=0.0;ratio_HF_year[processBin_ratio_nuis_Hproc]=0.0;
                    if (nuis=="RelBal"): up_RelBal[processBin_jesup_nuis_Hproc]=0.0; dn_RelBal[processBin_jesdn_nuis_Hproc]=0.0;ratio_RelBal[processBin_ratio_nuis_Hproc]=0.0;
                    if (nuis=="RelSample_year"): up_RelSample_year[processBin_jesup_nuis_Hproc]=0.0; dn_RelSample_year[processBin_jesdn_nuis_Hproc]=0.0;ratio_RelSample_year[processBin_ratio_nuis_Hproc]=0.0;
                    if (nuis=="Total"): up_Total[processBin_jesup_nuis_Hproc]=0.0; dn_Total[processBin_jesdn_nuis_Hproc]=0.0;ratio_Total[processBin_ratio_nuis_Hproc]=0.0;
 
		    yield_jesup_nuis_Hproc=0; yield_jesdn_nuis_Hproc=0;

                    for proc in procs:
                        if(not (proc=='trueH' or proc=='out_trueH' or proc=='fakeH')): continue
                        processBin_jesup_nuis = proc+'_'+str(x_point)+'_'+obsName+'_'+channel+'_recobin'+str(recobin)+'_jesup_'+nuis
                        processBin_jesdn_nuis = proc+'_'+str(x_point)+'_'+obsName+'_'+channel+'_recobin'+str(recobin)+'_jesdn_'+nuis
                        #processBin_ratio_nuis = proc+'_'+str(x_point)+'_'+obsName+'_'+channel+'_recobin'+str(recobin)+'_'+nuis
			#yield_jesup_nuis_Hproc+=	

			if (nuis=="Abs"): yield_jesup_nuis_Hproc+=up_Abs[processBin_jesup_nuis]; yield_jesdn_nuis_Hproc+=dn_Abs[processBin_jesdn_nuis];
			if (nuis=="Abs_year"): yield_jesup_nuis_Hproc+=up_Abs_year[processBin_jesup_nuis]; yield_jesdn_nuis_Hproc+=dn_Abs_year[processBin_jesdn_nuis];
			if (nuis=="BBEC1"): yield_jesup_nuis_Hproc+=up_BBEC1[processBin_jesup_nuis]; yield_jesdn_nuis_Hproc+=dn_BBEC1[processBin_jesdn_nuis];
			if (nuis=="BBEC1_year"): yield_jesup_nuis_Hproc+=up_BBEC1_year[processBin_jesup_nuis]; yield_jesdn_nuis_Hproc+=dn_BBEC1_year[processBin_jesdn_nuis];
			if (nuis=="EC2"): yield_jesup_nuis_Hproc+=up_EC2[processBin_jesup_nuis]; yield_jesdn_nuis_Hproc+=dn_EC2[processBin_jesdn_nuis];
			if (nuis=="EC2_year"): yield_jesup_nuis_Hproc+=up_EC2_year[processBin_jesup_nuis]; yield_jesdn_nuis_Hproc+=dn_EC2_year[processBin_jesdn_nuis];
			if (nuis=="FlavQCD"): yield_jesup_nuis_Hproc+=up_FlavQCD[processBin_jesup_nuis]; yield_jesdn_nuis_Hproc+=dn_FlavQCD[processBin_jesdn_nuis];
			if (nuis=="HF"): yield_jesup_nuis_Hproc+=up_HF[processBin_jesup_nuis]; yield_jesdn_nuis_Hproc+=dn_HF[processBin_jesdn_nuis];
			if (nuis=="HF_year"): yield_jesup_nuis_Hproc+=up_HF_year[processBin_jesup_nuis]; yield_jesdn_nuis_Hproc+=dn_HF_year[processBin_jesdn_nuis];
			if (nuis=="RelBal"): yield_jesup_nuis_Hproc+=up_RelBal[processBin_jesup_nuis]; yield_jesdn_nuis_Hproc+=dn_RelBal[processBin_jesdn_nuis];
			if (nuis=="RelSample_year"): yield_jesup_nuis_Hproc+=up_RelSample_year[processBin_jesup_nuis]; yield_jesdn_nuis_Hproc+=dn_RelSample_year[processBin_jesdn_nuis];
			if (nuis=="Total"): yield_jesup_nuis_Hproc+=up_Total[processBin_jesup_nuis]; yield_jesdn_nuis_Hproc+=dn_Total[processBin_jesdn_nuis];

   
    		    if (nuis=="Abs"): up_Abs[processBin_jesup_nuis_Hproc] = yield_jesup_nuis_Hproc; dn_Abs[processBin_jesdn_nuis_Hproc] = yield_jesdn_nuis_Hproc; ratio_Abs[processBin_ratio_nuis_Hproc] = computeRatio(yield_nom_Hproc,yield_jesup_nuis_Hproc,yield_jesdn_nuis_Hproc)
                    if (nuis=="Abs_year"): up_Abs_year[processBin_jesup_nuis_Hproc] = yield_jesup_nuis_Hproc; dn_Abs_year[processBin_jesdn_nuis_Hproc] = yield_jesdn_nuis_Hproc; ratio_Abs_year[processBin_ratio_nuis_Hproc] = computeRatio(yield_nom_Hproc,yield_jesup_nuis_Hproc,yield_jesdn_nuis_Hproc)
                    if (nuis=="BBEC1"): up_BBEC1[processBin_jesup_nuis_Hproc] = yield_jesup_nuis_Hproc; dn_BBEC1[processBin_jesdn_nuis_Hproc] = yield_jesdn_nuis_Hproc; ratio_BBEC1[processBin_ratio_nuis_Hproc] = computeRatio(yield_nom_Hproc,yield_jesup_nuis_Hproc,yield_jesdn_nuis_Hproc)
                    if (nuis=="BBEC1_year"): up_BBEC1_year[processBin_jesup_nuis_Hproc] = yield_jesup_nuis_Hproc; dn_BBEC1_year[processBin_jesdn_nuis_Hproc] = yield_jesdn_nuis_Hproc; ratio_BBEC1_year[processBin_ratio_nuis_Hproc] = computeRatio(yield_nom_Hproc,yield_jesup_nuis_Hproc,yield_jesdn_nuis_Hproc)
                    if (nuis=="EC2"): up_EC2[processBin_jesup_nuis_Hproc] = yield_jesup_nuis_Hproc; dn_EC2[processBin_jesdn_nuis_Hproc] = yield_jesdn_nuis_Hproc; ratio_EC2[processBin_ratio_nuis_Hproc] = computeRatio(yield_nom_Hproc,yield_jesup_nuis_Hproc,yield_jesdn_nuis_Hproc)
                    if (nuis=="EC2_year"): up_EC2_year[processBin_jesup_nuis_Hproc] = yield_jesup_nuis_Hproc; dn_EC2_year[processBin_jesdn_nuis_Hproc] = yield_jesdn_nuis_Hproc; ratio_EC2_year[processBin_ratio_nuis_Hproc] = computeRatio(yield_nom_Hproc,yield_jesup_nuis_Hproc,yield_jesdn_nuis_Hproc)
                    if (nuis=="FlavQCD"): up_FlavQCD[processBin_jesup_nuis_Hproc] = yield_jesup_nuis_Hproc; dn_FlavQCD[processBin_jesdn_nuis_Hproc] = yield_jesdn_nuis_Hproc; ratio_FlavQCD[processBin_ratio_nuis_Hproc] = computeRatio(yield_nom_Hproc,yield_jesup_nuis_Hproc,yield_jesdn_nuis_Hproc)
                    if (nuis=="HF"): up_HF[processBin_jesup_nuis_Hproc] = yield_jesup_nuis_Hproc; dn_HF[processBin_jesdn_nuis_Hproc] = yield_jesdn_nuis_Hproc; ratio_HF[processBin_ratio_nuis_Hproc] = computeRatio(yield_nom_Hproc,yield_jesup_nuis_Hproc,yield_jesdn_nuis_Hproc)
                    if (nuis=="HF_year"): up_HF_year[processBin_jesup_nuis_Hproc] = yield_jesup_nuis_Hproc; dn_HF_year[processBin_jesdn_nuis_Hproc] = yield_jesdn_nuis_Hproc; ratio_HF_year[processBin_ratio_nuis_Hproc] = computeRatio(yield_nom_Hproc,yield_jesup_nuis_Hproc,yield_jesdn_nuis_Hproc)
                    if (nuis=="RelBal"): up_RelBal[processBin_jesup_nuis_Hproc] = yield_jesup_nuis_Hproc; dn_RelBal[processBin_jesdn_nuis_Hproc] = yield_jesdn_nuis_Hproc; ratio_RelBal[processBin_ratio_nuis_Hproc] = computeRatio(yield_nom_Hproc,yield_jesup_nuis_Hproc,yield_jesdn_nuis_Hproc)
                    if (nuis=="RelSample_year"): up_RelSample_year[processBin_jesup_nuis_Hproc] = yield_jesup_nuis_Hproc; dn_RelSample_year[processBin_jesdn_nuis_Hproc] = yield_jesdn_nuis_Hproc; ratio_RelSample_year[processBin_ratio_nuis_Hproc] = computeRatio(yield_nom_Hproc,yield_jesup_nuis_Hproc,yield_jesdn_nuis_Hproc)
                    if (nuis=="Total"): up_Total[processBin_jesup_nuis_Hproc] = yield_jesup_nuis_Hproc; dn_Total[processBin_jesdn_nuis_Hproc] = yield_jesdn_nuis_Hproc; ratio_Total[processBin_ratio_nuis_Hproc] = computeRatio(yield_nom_Hproc,yield_jesup_nuis_Hproc,yield_jesdn_nuis_Hproc)
    



## interpolation part
#    for channel in channels:
#        for proc in procs:
#	    if(not (proc=='trueH' or proc=='out_trueH' or proc=='fakeH')): continue
#            for recobin in range(len(obs_bins)-1):
#                for nuis in JES_nuis:
#		    key_nom_125pt38 = proc+'_125.38_'+obsName+'_'+channel+'_recobin'+str(recobin)
#                    key_jesup_125pt38 = proc+'_125.38_'+obsName+'_'+channel+'_recobin'+str(recobin)+'_jesup_'+nuis
#                    key_jesdn_125pt38 = proc+'_125.38_'+obsName+'_'+channel+'_recobin'+str(recobin)+'_jesdn_'+nuis
#		    key_ratio_nuis = proc+'_125.38_'+obsName+'_'+channel+'_recobin'+str(recobin)+'_'+nuis
#    		    nom_points = []; jesUp_points = []; jesDn_points = [];
#                    for x_point in x_points:
#                        key_nom = proc+'_'+str(x_point)+'_'+obsName+'_'+channel+'_recobin'+str(recobin)
#                        key_jesup = proc+'_'+str(x_point)+'_'+obsName+'_'+channel+'_recobin'+str(recobin)+'_jesup_'+nuis
#                        key_jesdn = proc+'_'+str(x_point)+'_'+obsName+'_'+channel+'_recobin'+str(recobin)+'_jesdn_'+nuis
#			nom_points.append(nom[key_nom]); 
#			if (nuis=="Abs"): jesUp_points.append(up_Abs[key_jesup]); jesDn_points.append(dn_Abs[key_jesdn]); 
#			if (nuis=="Abs_year"): jesUp_points.append(up_Abs_year[key_jesup]); jesDn_points.append(dn_Abs_year[key_jesdn]); 
#			if (nuis=="BBEC1"): jesUp_points.append(up_BBEC1[key_jesup]); jesDn_points.append(dn_BBEC1[key_jesdn]); 
#			if (nuis=="BBEC1_year"): jesUp_points.append(up_BBEC1_year[key_jesup]); jesDn_points.append(dn_BBEC1_year[key_jesdn]); 
#			if (nuis=="EC2"): jesUp_points.append(up_EC2[key_jesup]); jesDn_points.append(dn_EC2[key_jesdn]); 
#			if (nuis=="EC2_year"): jesUp_points.append(up_EC2_year[key_jesup]); jesDn_points.append(dn_EC2_year[key_jesdn]); 
#			if (nuis=="FlavQCD"): jesUp_points.append(up_FlavQCD[key_jesup]); jesDn_points.append(dn_FlavQCD[key_jesdn]); 
#			if (nuis=="HF"): jesUp_points.append(up_HF[key_jesup]); jesDn_points.append(dn_HF[key_jesdn]); 
#			if (nuis=="HF_year"): jesUp_points.append(up_HF_year[key_jesup]); jesDn_points.append(dn_HF_year[key_jesdn]); 
#			if (nuis=="RelBal"): jesUp_points.append(up_RelBal[key_jesup]); jesDn_points.append(dn_RelBal[key_jesdn]); 
#			if (nuis=="RelSample_year"): jesUp_points.append(up_RelSample_year[key_jesup]); jesDn_points.append(dn_RelSample_year[key_jesdn]); 
#			if (nuis=="Total"): jesUp_points.append(up_Total[key_jesup]); jesDn_points.append(dn_Total[key_jesdn]); 
#
#		    print "nom_points:   ", nom_points
#		    print "jesUp_points:   ", jesUp_points
#		    print "jesDn_points:   ", jesDn_points			
#		    spl_jesUp_points  = interpolate.splrep(x_points, jesUp_points, k=2)
#		    spl_jesDn_points  = interpolate.splrep(x_points, jesDn_points, k=2)
#		    spl_nom_points  = interpolate.splrep(x_points, nom_points, k=2)
#		  
#		    #up[key_jesup_125pt38] = float(interpolate.splev(125.38, spl_jesUp_points))
#		    #dn[key_jesdn_125pt38] = float(interpolate.splev(125.38, spl_jesDn_points))
#		    nom[key_nom_125pt38] = float(interpolate.splev(125.38, spl_nom_points))
#		    #ratio[key_ratio_nuis] = computeRatio(nom[key_nom_125pt38],up[key_jesup_125pt38],dn[key_jesdn_125pt38])
#		    if (nuis=="Abs"): up_Abs[key_jesup_125pt38] = float(interpolate.splev(125.38, spl_jesUp_points)); dn_Abs[key_jesdn_125pt38] = float(interpolate.splev(125.38, spl_jesDn_points)); ratio_Abs[key_ratio_nuis] =  computeRatio(nom[key_nom_125pt38],up_Abs[key_jesup_125pt38],dn_Abs[key_jesdn_125pt38]);
#		    if (nuis=="Abs_year"): up_Abs_year[key_jesup_125pt38] = float(interpolate.splev(125.38, spl_jesUp_points)); dn_Abs_year[key_jesdn_125pt38] = float(interpolate.splev(125.38, spl_jesDn_points)); ratio_Abs_year[key_ratio_nuis] =  computeRatio(nom[key_nom_125pt38],up_Abs_year[key_jesup_125pt38],dn_Abs_year[key_jesdn_125pt38]);
#		    if (nuis=="BBEC1"): up_BBEC1[key_jesup_125pt38] = float(interpolate.splev(125.38, spl_jesUp_points)); dn_BBEC1[key_jesdn_125pt38] = float(interpolate.splev(125.38, spl_jesDn_points)); ratio_BBEC1[key_ratio_nuis] =  computeRatio(nom[key_nom_125pt38],up_BBEC1[key_jesup_125pt38],dn_BBEC1[key_jesdn_125pt38]);
#		    if (nuis=="BBEC1_year"): up_BBEC1_year[key_jesup_125pt38] = float(interpolate.splev(125.38, spl_jesUp_points)); dn_BBEC1_year[key_jesdn_125pt38] = float(interpolate.splev(125.38, spl_jesDn_points)); ratio_BBEC1_year[key_ratio_nuis] =  computeRatio(nom[key_nom_125pt38],up_BBEC1_year[key_jesup_125pt38],dn_BBEC1_year[key_jesdn_125pt38]);
#		    if (nuis=="EC2"): up_EC2[key_jesup_125pt38] = float(interpolate.splev(125.38, spl_jesUp_points)); dn_EC2[key_jesdn_125pt38] = float(interpolate.splev(125.38, spl_jesDn_points)); ratio_EC2[key_ratio_nuis] =  computeRatio(nom[key_nom_125pt38],up_EC2[key_jesup_125pt38],dn_EC2[key_jesdn_125pt38]);
#		    if (nuis=="EC2_year"): up_EC2_year[key_jesup_125pt38] = float(interpolate.splev(125.38, spl_jesUp_points)); dn_EC2_year[key_jesdn_125pt38] = float(interpolate.splev(125.38, spl_jesDn_points)); ratio_EC2_year[key_ratio_nuis] =  computeRatio(nom[key_nom_125pt38],up_EC2_year[key_jesup_125pt38],dn_EC2_year[key_jesdn_125pt38]);
#		    if (nuis=="FlavQCD"): up_FlavQCD[key_jesup_125pt38] = float(interpolate.splev(125.38, spl_jesUp_points)); dn_FlavQCD[key_jesdn_125pt38] = float(interpolate.splev(125.38, spl_jesDn_points)); ratio_FlavQCD[key_ratio_nuis] =  computeRatio(nom[key_nom_125pt38],up_FlavQCD[key_jesup_125pt38],dn_FlavQCD[key_jesdn_125pt38]);
#		    if (nuis=="HF"): up_HF[key_jesup_125pt38] = float(interpolate.splev(125.38, spl_jesUp_points)); dn_HF[key_jesdn_125pt38] = float(interpolate.splev(125.38, spl_jesDn_points)); ratio_HF[key_ratio_nuis] =  computeRatio(nom[key_nom_125pt38],up_HF[key_jesup_125pt38],dn_HF[key_jesdn_125pt38]);
#		    if (nuis=="HF_year"): up_HF_year[key_jesup_125pt38] = float(interpolate.splev(125.38, spl_jesUp_points)); dn_HF_year[key_jesdn_125pt38] = float(interpolate.splev(125.38, spl_jesDn_points)); ratio_HF_year[key_ratio_nuis] =  computeRatio(nom[key_nom_125pt38],up_HF_year[key_jesup_125pt38],dn_HF_year[key_jesdn_125pt38]);
#		    if (nuis=="RelBal"): up_RelBal[key_jesup_125pt38] = float(interpolate.splev(125.38, spl_jesUp_points)); dn_RelBal[key_jesdn_125pt38] = float(interpolate.splev(125.38, spl_jesDn_points)); ratio_RelBal[key_ratio_nuis] =  computeRatio(nom[key_nom_125pt38],up_RelBal[key_jesup_125pt38],dn_RelBal[key_jesdn_125pt38]);
#		    if (nuis=="RelSample_year"): up_RelSample_year[key_jesup_125pt38] = float(interpolate.splev(125.38, spl_jesUp_points)); dn_RelSample_year[key_jesdn_125pt38] = float(interpolate.splev(125.38, spl_jesDn_points)); ratio_RelSample_year[key_ratio_nuis] =  computeRatio(nom[key_nom_125pt38],up_RelSample_year[key_jesup_125pt38],dn_RelSample_year[key_jesdn_125pt38]);
#		    if (nuis=="Total"): up_Total[key_jesup_125pt38] = float(interpolate.splev(125.38, spl_jesUp_points)); dn_Total[key_jesdn_125pt38] = float(interpolate.splev(125.38, spl_jesDn_points)); ratio_Total[key_ratio_nuis] =  computeRatio(nom[key_nom_125pt38],up_Total[key_jesup_125pt38],dn_Total[key_jesdn_125pt38]);
#

        #with open('datacardInputs/inputs_JESnuis_'+opt.OBSNAME+'_'+channel+'_'+proc+'_'+opt.YEAR+'.py', 'w') as f:
    with open('datacardInputs/'+opt.YEAR+'/'+obs_reco+'/inputs_JESnuis_'+opt.OBSNAME+'.py', 'w') as f:
        #f.write('obsName = "'+str(opt.OBSNAME)+'" \n')
        #f.write('obs_bins = '+str(obs_bins)+' \n')
        #f.write('JES uncertainty sources = '+str(JES_nuis)+' \n')
        f.write('yield_nom = '+str(nom)+' \n')
        f.write('yield_up_Abs = '+str(up_Abs)+' \n')
        f.write('yield_dn_Abs = '+str(dn_Abs)+' \n')
        f.write('nuis_Abs = '+str(ratio_Abs)+' \n')
        f.write('yield_up_Abs_year = '+str(up_Abs_year)+' \n')
        f.write('yield_dn_Abs_year = '+str(dn_Abs_year)+' \n')
        f.write('nuis_Abs_year = '+str(ratio_Abs_year)+' \n')
        f.write('yield_up_BBEC1 = '+str(up_BBEC1)+' \n')
        f.write('yield_dn_BBEC1 = '+str(dn_BBEC1)+' \n')
        f.write('nuis_BBEC1 = '+str(ratio_BBEC1)+' \n')
        f.write('yield_up_BBEC1_year = '+str(up_BBEC1_year)+' \n')
        f.write('yield_dn_BBEC1_year = '+str(dn_BBEC1_year)+' \n')
        f.write('nuis_BBEC1_year = '+str(ratio_BBEC1_year)+' \n')
        f.write('yield_up_EC2 = '+str(up_EC2)+' \n')
        f.write('yield_dn_EC2 = '+str(dn_EC2)+' \n')
        f.write('nuis_EC2 = '+str(ratio_EC2)+' \n')
        f.write('yield_up_EC2_year = '+str(up_EC2_year)+' \n')
        f.write('yield_dn_EC2_year = '+str(dn_EC2_year)+' \n')
        f.write('nuis_EC2_year = '+str(ratio_EC2_year)+' \n')
        f.write('yield_up_FlavQCD = '+str(up_FlavQCD)+' \n')
        f.write('yield_dn_FlavQCD = '+str(dn_FlavQCD)+' \n')
        f.write('nuis_FlavQCD = '+str(ratio_FlavQCD)+' \n')
        f.write('yield_up_HF = '+str(up_HF)+' \n')
        f.write('yield_dn_HF = '+str(dn_HF)+' \n')
        f.write('nuis_HF = '+str(ratio_HF)+' \n')
        f.write('yield_up_HF_year = '+str(up_HF_year)+' \n')
        f.write('yield_dn_HF_year = '+str(dn_HF_year)+' \n')
        f.write('nuis_HF_year = '+str(ratio_HF_year)+' \n')
        f.write('yield_up_RelBal = '+str(up_RelBal)+' \n')
        f.write('yield_dn_RelBal = '+str(dn_RelBal)+' \n')
        f.write('nuis_RelBal = '+str(ratio_RelBal)+' \n')
        f.write('yield_up_RelSample_year = '+str(up_RelSample_year)+' \n')
        f.write('yield_dn_RelSample_year = '+str(dn_RelSample_year)+' \n')
        f.write('nuis_RelSample_year = '+str(ratio_RelSample_year)+' \n')
        f.write('yield_up_Total = '+str(up_Total)+' \n')
        f.write('yield_dn_Total = '+str(dn_Total)+' \n')
        f.write('nuis_Total = '+str(ratio_Total)+' \n')



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

    if (not os.path.exists("datacardInputs/"+opt.YEAR+"/"+obs_reco)):
    #os.system("mkdir -p datacardInputs")
        os.system("mkdir -p datacardInputs/"+opt.YEAR+"/"+obs_reco+"")


 
    m4l_low = 105.0
    m4l_high = 160.0
    m4l_bins = 35
    #m4l_bins = 70  ## temporary

    print "era being processed: ",opt.YEAR

    Histos = {}
    nom = {}
    #up = {}
    #dn = {}
    up_Abs = {}; up_Abs_year = {}; up_BBEC1 = {}; up_BBEC1_year = {}; up_EC2 = {}; up_EC2_year = {}; up_FlavQCD = {}; up_HF = {}; up_HF_year = {}; up_RelBal = {}; up_RelSample_year = {}; up_Total = {};
    dn_Abs = {}; dn_Abs_year = {}; dn_BBEC1 = {}; dn_BBEC1_year = {}; dn_EC2 = {}; dn_EC2_year = {}; dn_FlavQCD = {}; dn_HF = {}; dn_HF_year = {}; dn_RelBal = {}; dn_RelSample_year = {}; dn_Total = {};
    ratio_Abs = {}; ratio_Abs_year = {}; ratio_BBEC1 = {}; ratio_BBEC1_year = {}; ratio_EC2 = {}; ratio_EC2_year = {}; ratio_FlavQCD = {}; ratio_HF = {}; ratio_HF_year = {}; ratio_RelBal = {}; ratio_RelSample_year = {}; ratio_Total = {};

    ratio = {} 

#    RootFile, Tree, nEvents, sumw = GrabMCTrees(opt.YEAR)
#    print "sumw: ", sumw

    if (opt.YEAR=='2018' or opt.YEAR=='2017'):
        RootFile, Tree, nEvents, sumw = GrabMCTrees(opt.YEAR)
    else:
        print "getting trees from pre "
        RootFile_pre, Tree_pre, nEvents_pre, sumw_pre = GrabMCTrees('2016preVFP')
        print "getting trees from post "
        RootFile_post, Tree_post, nEvents_post, sumw_post = GrabMCTrees('2016postVFP')



    #JES_nuis = ['Abs','FlavQCD'] 
    #JES_nuis = ['Abs'] 
    JES_nuis = ['Abs','Abs_year','BBEC1','BBEC1_year','EC2','EC2_year','FlavQCD','HF','HF_year','RelBal','RelSample_year', 'Total']
    print("Obs Name: {:15}  nBins: {:2}  bins: {}".format(opt.OBSNAME, nbins, observableBins))

    extract_JES_nuis(nbins, opt.OBSNAME, obs_bins, opt.DEBUG)
   


    print("Interpolation completed... :) ")
