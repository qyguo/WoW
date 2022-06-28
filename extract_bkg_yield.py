import optparse
import os
import sys
import ROOT, re
from ROOT import *
from scipy import interpolate
import matplotlib.pyplot as plt
#from library import list_files as DATA
#from LoadData import *
from LoadData_yield import *
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
    parser.add_option('',   '--obsName',  dest='OBSNAME',  type='string',default='mass4l',   help='Name of the observable, supported: "inclusive", "pT4l", "eta4l", "massZ2", "nJets"')
    parser.add_option('',   '--obsBins',  dest='OBSBINS',  type='string',default='|105.0|160.0|',   help='Bin boundaries for the diff. measurement separated by "|", e.g. as "|0|50|100|", use the defalut if empty string')
    parser.add_option('',   '--year',  dest='YEAR',  type='string',default='2016',   help='Era to analyze, e.g. 2016, 2017, 2018 or Full ')
    parser.add_option('', '--bkg',      dest='BKG',type='string',default='', help='run with the type of zz background to float zz or qq_gg ')
    parser.add_option('',   '--debug',  dest='DEBUG',  type='int',default=0,   help='0 if debug false, else debug True')
    global opt, args
    (opt, args) = parser.parse_args()

def extract_JES_nuis(nbins, obsName, obs_bins, DEBUG = 0):
    x_points = [125]
    channels=['4mu','2e2mu','4e']
    #procs=['bkg_zjets','trueH','out_trueH','fakeH','bkg_qqzz','bkg_ggzz']
    #procs=['bkg_zjets','bkg_qqzz','bkg_ggzz']
    procs=['bkg_qqzz','bkg_ggzz']
    #procs = ['bkg_zjets']
    #procs = ['bkg_ggzz']


    #proc_selections = {"trueH":"passedFullSelection==1 ","out_trueH":"(passedFullSelection==1 && passedFiducialSelection!=1)","fakeH":"isH4l!=1","bkg_qqzz":"passedFullSelection==1","bkg_ggzz":"passedFullSelection==1","bkg_zjets":"(passedZXCRSelection==1)"}
    proc_selections = {"trueH":"passedFullSelection==1 ","out_trueH":"(passedFullSelection==1 && passedFiducialSelection!=1)","fakeH":"isH4l!=1","bkg_qqzz":"passedFullSelection==1","bkg_ggzz":"passedFullSelection==1","bkg_zjets":"passedFullSelection==1"}


    cut_m4l_reco = {"2e2mu":"(mass2e2mu>0)","4mu":"(mass4mu>0)","4e":"(mass4e>0)"}


    #lumi = {"2016":35867.0,"2017":41370.0,"2018":58800.0}
    #lumi = {"2016":35900.0,"2017":41500.0,"2018":59700.0}
    lumi = {"2016":35900.0,"2017":41500.0,"2018":58970.48}
    lumi_preVFP=19520; lumi_postVFP=16810;

    weights = {"sig":"genWeight*pileupWeight*prefiringWeight*dataMCWeight_new*crossSection","bkg_qqzz":"(genWeight*pileupWeight*prefiringWeight*dataMCWeight_new*crossSection)*(k_qqZZ_qcd_M*k_qqZZ_ewk)","bkg_ggzz":"(genWeight*pileupWeight*prefiringWeight*dataMCWeight_new*crossSection)*(k_ggZZ)","bkg_zjets":"1.0"}
    #weights = {"sig":"genWeight*pileupWeight*dataMCWeight_new*crossSection","bkg_qqzz":"(genWeight*pileupWeight*dataMCWeight_new*crossSection)*(k_qqZZ_qcd_M*k_qqZZ_ewk)","bkg_ggzz":"(genWeight*pileupWeight*dataMCWeight_new*crossSection)*(k_ggZZ)","bkg_zjets":"1.0"}
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

		    if (opt.YEAR=='2018' or opt.YEAR=='2017'):
                        RootFile[sample] = TFile(dirMC[opt.YEAR]+'/'+sample+'.root',"READ")  # currently, all files in the same directory
           	        RootFile_list.append(dirMC[opt.YEAR]+'/'+sample+'.root')
		    #else:
		    elif (opt.YEAR=='2016'):
		        RootFile_pre[sample] = TFile(dirMC['2016preVFP']+'/'+sample+'.root',"READ")
			RootFile_list_pre.append(dirMC['2016preVFP']+'/'+sample+'.root')   
		        RootFile_post[sample] = TFile(dirMC['2016postVFP']+'/'+sample+'.root',"READ")
			RootFile_list_post.append(dirMC['2016postVFP']+'/'+sample+'.root')   

		    print "sample is :   ", sample
		    if (proc!='bkg_zjets' and (opt.YEAR=='2018' or opt.YEAR=='2017')): 
			sumw2.append(sumw[sample])  # from keys to list, for onward use
		    elif (proc!='bkg_zjets' and (opt.YEAR=='2016')): 
			sumw2_pre.append(sumw_pre[sample]); sumw2_post.append(sumw_post[sample])
		if (opt.YEAR=='2018' or opt.YEAR=='2017'):
                    files = [ROOT.TFile(i) for i in RootFile_list]
		    print "proc:  ", proc, "  channel:   ", channel
		    print "files are:  ", files
		elif (opt.YEAR=='2016'):
		    files_pre = [ROOT.TFile(i) for i in RootFile_list_pre]
		    files_post = [ROOT.TFile(i) for i in RootFile_list_post]
		    print "proc:  ", proc, "  channel:   ", channel
		    print "files_pre are:  ", files_pre
		    print "files_post are:  ", files_post

		if (opt.YEAR=='2018' or opt.YEAR=='2017'):
		    if(proc=='bkg_zjets'): trees = [i.Get("passedEvents") for i in files]
                    else: trees = [i.Get("Ana/passedEvents") for i in files]
		    print "trees are:  ", trees
		elif (opt.YEAR=='2016'):
		    if(proc=='bkg_zjets'): trees_pre = [i.Get("passedEvents") for i in files_pre]; trees_post = [j.Get("passedEvents") for j in files_post];
                    else: trees_pre = [i.Get("Ana/passedEvents") for i in files_pre]; trees_post = [j.Get("Ana/passedEvents") for j in files_post];


		for recobin in range(len(obs_bins)-1):
		    obs_reco_low = obs_bins[recobin]
    		    obs_reco_high = obs_bins[recobin+1]

		    print "obs_reco_low", obs_reco_low, "     obs_reco_high", obs_reco_high
		    cutobs_reco = "("+obs_reco+">="+str(obs_reco_low)+" && "+obs_reco+"<"+str(obs_reco_high)+")"
		    cut_nom = "("+recoweight+")*("+cutobs_reco+" && "+proc_selections[proc]+" && "+cut_m4l_reco[channel]+")"

		    processBin_nom = proc+'_'+str(x_point)+'_'+obsName+'_'+channel+'_recobin'+str(recobin)
		    nom[processBin_nom] = 0.0; 
		
		    if (opt.YEAR=='2018' or opt.YEAR=='2017'):
    		        Histos[processBin_nom] = TH1D(processBin_nom,processBin_nom, m4l_bins, m4l_low, m4l_high)
		        Histos[processBin_nom].Sumw2()
		    elif (opt.YEAR=='2016'):
			Histos[processBin_nom+'pre'] = TH1D(processBin_nom+'pre',processBin_nom+'pre', m4l_bins, m4l_low, m4l_high)
			Histos[processBin_nom+'pre'].Sumw2()
			Histos[processBin_nom+'post'] = TH1D(processBin_nom+'post',processBin_nom+'post', m4l_bins, m4l_low, m4l_high)
			Histos[processBin_nom+'post'].Sumw2()

 		    print "cut_nom :   ", cut_nom 
#			print "cut_nuis_up :   ", cut_nuis_up
#			print "cut_nuis_dn :   ", cut_nuis_dn 
#
		    yield_nom = 0.;
    		    i=0;  # iterator to read sumw2 from corresponding root file, meant to work for signal proc. only
    
		    red_bkg = TH1D('red_bkg', 'red_bkg',m4l_bins,m4l_low,m4l_high)
    		    if (opt.YEAR=='2018' or opt.YEAR=='2017'):
    		        for tree in trees:
    		            if(proc=='bkg_zjets'):
                                #tree.Draw("mass4l >> "+processBin_nom,"("+cut_nom+")","goff")
                                tree.Draw("mass4l >> "+processBin_nom,"("+cut_nom+")","goff")
    		            else:
    		                print "proc is:   ", proc
    		                print "channel is:   ", channel
    		                print "i:  ",i
    			        print "len(sumw2):    ", len(sumw2)
    		                print "sumw2[",i,"] : ", sumw2[i] 
    		                #tree.Draw("mass4l >> "+processBin_nom,"("+cut_nom+"*("+str(lumi[opt.YEAR])+")/"+str(sumw2[i])+")","goff")
    		                tree.Draw("mass4l >> "+processBin_nom,"("+cut_nom+"*("+str(lumi[opt.YEAR])+")/"+str(sumw2[i])+")","goff")
    		            i=i+1
    
    			    print "yield Integral", Histos[processBin_nom].Integral()
    		            yield_nom+=Histos[processBin_nom].Integral()
			#nom[processBin_nom] = yield_nom 
    
                    elif (opt.YEAR=='2016'):
    		        for pre_tree, post_tree in zip(trees_pre, trees_post):    
    		            print "pre_tree is:  ", pre_tree
    		            print "post_tree is:  ", post_tree
    
    		            if(proc=='bkg_zjets'):
                                #pre_tree.Draw("mass4l >> "+processBin_nom+'pre',"("+cut_nom+")","goff")
                                pre_tree.Draw("mass4l >> "+processBin_nom+'pre',"("+cut_nom+")","goff")
    
                                #post_tree.Draw("mass4l >> "+processBin_nom+'post',"("+cut_nom+")","goff")
                                post_tree.Draw("mass4l >> "+processBin_nom+'post',"("+cut_nom+")","goff")
				
                            else:
    
                                print "proc is:   ", proc
                                print "channel is:   ", channel
                                print "i:  ",i
                                print "len(sumw2_pre):    ", len(sumw2_pre)
                                print "sumw2_pre[",i,"] : ", sumw2_pre[i]
                                print "len(sumw2_post):    ", len(sumw2_post)
                                print "sumw2_post[",i,"] : ", sumw2_post[i]
    
                                pre_tree.Draw("mass4l >> "+processBin_nom+'pre',"("+cut_nom+"*("+str(lumi_preVFP)+")/"+str(sumw2_pre[i])+")","goff")
                                post_tree.Draw("mass4l >> "+processBin_nom+'post',"("+cut_nom+"*("+str(lumi_postVFP)+")/"+str(sumw2_post[i])+")","goff")
    
                            i=i+1
    
    			    print "pre Integral", Histos[processBin_nom+'pre'].Integral()
    			    print "post Integral", Histos[processBin_nom+'post'].Integral()

		            yield_nom+=Histos[processBin_nom+'pre'].Integral()+Histos[processBin_nom+'post'].Integral() ## FIXME a check
			#nom[processBin_nom] = yield_nom
			Histos[processBin_nom] = Histos[processBin_nom+'pre'] + Histos[processBin_nom+'post']  # for later use to fit for Landau function
                
                    if(proc=='bkg_zjets'):			    
                                ### red_bkg ###
                        if (channel=="4e"):
                            fsum = TF1("fsum","landau(0)")
                            fsum.SetParameter(0,1.0)
                            fsum.SetParameter(1,141.9)
                            fsum.SetParameter(2,21.3)
                        if (channel=="4mu"):
                            fsum = TF1("fsum","landau(0)")
                            fsum.SetParameter(0,1.0)
                            fsum.SetParameter(1,130.4)
                            fsum.SetParameter(2,15.6)
                        if (channel=="2e2mu"):
                            fsum = TF1("fsum","landau(0)+landau(3)")
                            fsum.SetParameter(0,0.55)
                            fsum.SetParameter(1,131.1)
                            fsum.SetParameter(2,18.1)
                            fsum.SetParameter(3,0.45)
                            fsum.SetParameter(4,133.8)
                            fsum.SetParameter(5,18.9)
                        if (channel=="4l"):
                            fsum = TF1("fsum","landau(0)+landau(3)+landau(6)+landau(9)")
                            fsum.SetParameter(0,9.8)
                            fsum.SetParameter(1,141.9)
                            fsum.SetParameter(2,21.3)
                            fsum.SetParameter(3,10.2)
                            fsum.SetParameter(4,130.4)
                            fsum.SetParameter(5,15.6)
                            fsum.SetParameter(6,0.55*20.4)
                            fsum.SetParameter(7,131.1)
                            fsum.SetParameter(8,18.1)
                            fsum.SetParameter(9,0.45*20.4)
                            fsum.SetParameter(10,133.8)
                            fsum.SetParameter(11,18.9)
                        
                        nfill=100000
                        htemp = TH1F("htemp","htemp",1300,70,2000);
                        htemp.FillRandom("fsum",nfill);
                        if (channel=="4e"): htemp.Scale(21.1/htemp.Integral(htemp.FindBin(70),htemp.FindBin(2000)))
                        if (channel=="4mu"): htemp.Scale(34.4/htemp.Integral(htemp.FindBin(70),htemp.FindBin(2000)))
                        if (channel=="2e2mu"): htemp.Scale(59.9/htemp.Integral(htemp.FindBin(70),htemp.FindBin(2000)))
                        if (channel=="4l"): htemp.Scale(115.4/htemp.Integral(htemp.FindBin(70),htemp.FindBin(2000)))
                        Histos[processBin_nom].FillRandom("fsum",nfill);
                        if (Histos[processBin_nom].Integral(Histos[processBin_nom].FindBin(105),Histos[processBin_nom].FindBin(160)) >0 ):
                            #Histos[processBin_nom].Scale(htemp.Integral(htemp.FindBin(105),htemp.FindBin(160))/Histos[processBin_nom].Integral(Histos[processBin_nom].FindBin(105),red_bkg.FindBin(160)))
                            Histos[processBin_nom].Scale(htemp.Integral(htemp.FindBin(m4l_low),htemp.FindBin(m4l_high))/Histos[processBin_nom].Integral(Histos[processBin_nom].FindBin(m4l_low),red_bkg.FindBin(m4l_high)))
        	  	    nom[processBin_nom] = Histos[processBin_nom].Integral(1,m4l_bins)	
        
        		    print 'Red Bkg:  ',Histos[processBin_nom].Integral(1,m4l_bins)
                        
		    #nom[processBin_nom] = yield_nom
		    nom[processBin_nom] = Histos[processBin_nom].Integral(1,m4l_bins)



    			#(lumi_part1*sf1+ lumi_part2*sf2)/(lumi_part1+lumi_part2)
                            #yield_nom+=Histos[processBin_nom].Integral()
#                            yield_nom+=(((Histos[processBin_nom+'pre'].Integral())*lumi_preVFP)+((Histos[processBin_nom+'post'].Integral())*lumi_postVFP))/(lumi_preVFP+lumi_postVFP)
	                #yield_nom+=Histos[processBin_nom+'pre'].Integral()+Histos[processBin_nom+'post'].Integral() ## FIXME a check
		    #nom[processBin_nom] = yield_nom 
#
#
    with open('datacardInputs/'+opt.YEAR+'/'+obs_reco+'/inputs_yields_'+opt.OBSNAME+'.py', 'w') as f:
        f.write('yield_nom = '+str(nom)+' \n')

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
        os.system("mkdir -p datacardInputs/"+opt.YEAR+"/"+obs_reco+"")

    m4l_low = 105.0
    m4l_high = 160.0
    #m4l_bins = 35
    m4l_bins = 40
    #m4l_bins = 70  ## temporary

    print "era being processed: ",opt.YEAR

    Histos = {}
    nom = {}

    if (opt.YEAR=='2018' or opt.YEAR=='2017'):
        RootFile, Tree, nEvents, sumw = GrabMCTrees(opt.YEAR)
    else:
        print "getting trees from pre "
        RootFile_pre, Tree_pre, nEvents_pre, sumw_pre = GrabMCTrees('2016preVFP')
        print "getting trees from post "
        RootFile_post, Tree_post, nEvents_post, sumw_post = GrabMCTrees('2016postVFP')


    print("Obs Name: {:15}  nBins: {:2}  bins: {}".format(opt.OBSNAME, nbins, observableBins))

    extract_JES_nuis(nbins, opt.OBSNAME, obs_bins, opt.DEBUG)
   


    print("Interpolation completed... :) ")
