#include <iostream>
#include <set>
#include <TString.h>
#include <TFile.h>
#include <TTree.h>
#include <TMath.h>
#include <memory>
#include <string>
#include <map>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <algorithm>
#include <stdlib.h>
#include <stdio.h>
#include <cmath>
#include <iomanip>
#include <set>

#define PI 3.14159

// user include files 
#include "TROOT.h"
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TSpline.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TString.h"
#include "TLorentzVector.h"
#include "Math/VectorUtil.h"
#include "TClonesArray.h"
#include "TCanvas.h"

float dataMC_(float pt_, float eta_, TH2F* hMuScaleFac)
{
    float pt = TMath::Min(pt_,(float)199.0);
    float eta = eta_;
    //float eta = TMath::Abs(eta_);
    return hMuScaleFac->GetBinContent(hMuScaleFac->FindBin(eta,pt));
}

float dataMCErr_(float pt_, float eta_, TH2F* hMuScaleFacUnc)
{
    float pt = TMath::Min(pt_,(float)199.0);
    float eta = eta_;
    //float eta = TMath::Abs(eta_);
    return hMuScaleFacUnc->GetBinContent(hMuScaleFacUnc->FindBin(eta,pt));
}

void slimNtuple(const int & _year_=2017, const string & _name_DS_="bbH_HToZZTo4L_M125_TuneCP2_13TeV-jhugenv7011-pythia8", const bool & isMC = true,  const bool & isSignal = true, const bool & _Test=false) {
    std::cout<<"year: "<<_year_<<"; name_DS: "<<_name_DS_<<"; isMC: "<<isMC<<"; isSignal: "<<isSignal<<"; _Test: "<<_Test<<std::endl;

    int year = _year_;
    int Year = TMath::Abs(year);
    //string pre_name = "root://cmsio5.rc.ufl.edu:1094//store/user/t2/users/ferrico/Qianying/";
    //string pre_name = "/publicfs/cms/data/hzz/guoqy/newNTuple_UL/" + whichYear[year-2016] + "/MC/";
    string pre_name = "/publicfs/cms/data/hzz/guoqy/newNTuple_UL/";
    if (year==2016)    pre_name="";
    else if (year==-2016) //UL16APV
        pre_name="";

    string whichYear[3] = {"2016/", "2017/", "2018/"};
    string MC__ = "MC/";
    if (!isMC) MC__ = "Data/";
    //string name_DS = "GluGluHToZZTo4L_M125_2017";
    string name_DS = _name_DS_.c_str();

    //TString filename = prefix+".root";
    string filename = (pre_name + whichYear[year-2016] + MC__ + name_DS + ".root").c_str();
    string outputFile_Slimmed = (pre_name + whichYear[year-2016] + "Slimmed/").c_str(); 

    std::cout<<"Year: "<<year<<std::endl;
    std::cout<<filename<<std::endl;

    // UL new muon SFs
    string mu_scalefac_name_161718[3] = {"final_HZZ_SF_2016UL_mupogsysts_newLoose.root", "final_HZZ_SF_2017UL_mupogsysts_newLoose.root", "final_HZZ_SF_2018UL_mupogsysts_newLoose.root"};
    //string mu_scalefacFileInPath = ("/publicfs/cms/data/hzz/guoqy/newNTuple_UL/newMuonSF/" + mu_scalefac_name_161718[Year-2016]).c_str();
    string mu_scalefacFileInPath = ("newMuonSF/" + mu_scalefac_name_161718[Year-2016]).c_str();
    TFile *fMuScalFac = TFile::Open(mu_scalefacFileInPath.c_str());
    TH2F *hMuScaleFac = (TH2F*)fMuScalFac->Get("FINAL");
    TH2F *hMuScaleFacUnc = (TH2F*)fMuScalFac->Get("ERROR");

    TFile *oldfile = TFile::Open(filename.c_str());
    if (isMC) oldfile->cd("Ana");
    //TTree *oldtree = (TTree*)oldfile->Get("Ana/passedEvents");
    TTree *oldtree = (TTree*)gDirectory->Get("passedEvents");
    TH1F *th[7];

    if (isMC)
    {
        th[0]=(TH1F*)oldfile->Get("Ana/nEvents");
        th[1]=(TH1F*)oldfile->Get("Ana/sumWeights");
        th[2]=(TH1F*)oldfile->Get("Ana/sumWeightsPU");
        th[3]=(TH1F*)oldfile->Get("Ana/nVtx");
        th[4]=(TH1F*)oldfile->Get("Ana/nVtx_ReWeighted");
        th[5]=(TH1F*)oldfile->Get("Ana/nInteractions");
        th[6]=(TH1F*)oldfile->Get("Ana/nInteraction_ReWeighted");
    }
    Long64_t nentries = oldtree->GetEntries();
    std::cout<<nentries<<" total entries."<<std::endl;
    //ULong64_t Run, LumiSect, Event;
    //Bool_t          passedFullSelection;
    Float_t         dataMCWeight;
    Float_t         eventWeight;
    vector<int>     *lep_id;
    vector<float>   *lep_pt;
    vector<float>   *lep_eta;
    Int_t           lep_Hindex[4];
    vector<float>   *lep_dataMC;
    vector<float>   *lep_dataMCErr;
    vector<float>   lep_dataMC_new;
    vector<float>   lep_dataMCErr_new;
    Float_t         dataMCWeight_new;
    Float_t         eventWeight_new;
    lep_id = 0; 
    lep_pt = 0;    
    lep_eta = 0;  
    lep_dataMC = 0;
    lep_dataMCErr = 0;
    if (isMC)
    {
        //oldtree->SetBranchAddress("Run",&Run);
        //oldtree->SetBranchAddress("LumiSect",&LumiSect);
        //oldtree->SetBranchAddress("Event",&Event);
        //oldtree->SetBranchAddress("passedFullSelection", &passedFullSelection);
        oldtree->SetBranchAddress("dataMCWeight", &dataMCWeight);
        oldtree->SetBranchAddress("eventWeight", &eventWeight);
        oldtree->SetBranchAddress("lep_id", &lep_id);
        oldtree->SetBranchAddress("lep_pt", &lep_pt);
        oldtree->SetBranchAddress("lep_Hindex", lep_Hindex);
        oldtree->SetBranchAddress("lep_eta", &lep_eta);
        oldtree->SetBranchAddress("lep_dataMC", &lep_dataMC);
        oldtree->SetBranchAddress("lep_dataMCErr", &lep_dataMCErr);
    }
    
    Bool_t passedZ4lSelection;
    Bool_t passedZXCRSelection;
    oldtree->SetBranchAddress("passedZ4lSelection", &passedZ4lSelection);
    oldtree->SetBranchAddress("passedZXCRSelection", &passedZXCRSelection);

    oldtree->SetBranchStatus("*",0); //Disables All Branches
    //Then enables only select branches
    //
    if (isMC){
        oldtree->SetBranchStatus("lep_id",1);
        oldtree->SetBranchStatus("lep_pt",1);
        oldtree->SetBranchStatus("lep_eta",1);
        oldtree->SetBranchStatus("lep_dataMC",1);
        oldtree->SetBranchStatus("lep_dataMCErr",1);
    }

    oldtree->SetBranchStatus("Run",1);
    oldtree->SetBranchStatus("LumiSect",1);
    oldtree->SetBranchStatus("Event",1);
    oldtree->SetBranchStatus("passedFullSelection",1);
    oldtree->SetBranchStatus("passedFiducialSelection",1);

    oldtree->SetBranchStatus("passedZ4lSelection",1);
    oldtree->SetBranchStatus("passedZ1LSelection",1);
    oldtree->SetBranchStatus("passedZ4lZ1LSelection",1);
    oldtree->SetBranchStatus("passedZ4lZXCRSelection",1);
    oldtree->SetBranchStatus("passedZXCRSelection",1);
    oldtree->SetBranchStatus("nZXCRFailedLeptons",1);

    oldtree->SetBranchStatus("genWeight",1);
    oldtree->SetBranchStatus("pileupWeight",1);
    oldtree->SetBranchStatus("prefiringWeight",1);
    oldtree->SetBranchStatus("dataMCWeight",1);
    oldtree->SetBranchStatus("eventWeight",1);
    oldtree->SetBranchStatus("crossSection",1);
    oldtree->SetBranchStatus("k_qqZZ_qcd_M",1);
    oldtree->SetBranchStatus("k_qqZZ_ewk",1);
    oldtree->SetBranchStatus("k_ggZZ",1);
    oldtree->SetBranchStatus("mass4lREFIT",1);

    oldtree->SetBranchStatus("finalState",1);
    oldtree->SetBranchStatus("idL1",1);
    oldtree->SetBranchStatus("idL2",1);
    oldtree->SetBranchStatus("idL3",1);
    oldtree->SetBranchStatus("idL4",1);
    oldtree->SetBranchStatus("pTL1",1);
    oldtree->SetBranchStatus("pTL2",1);
    oldtree->SetBranchStatus("pTL3",1);
    oldtree->SetBranchStatus("pTL4",1);
    oldtree->SetBranchStatus("etaL1",1);
    oldtree->SetBranchStatus("etaL2",1);
    oldtree->SetBranchStatus("etaL3",1);
    oldtree->SetBranchStatus("etaL4",1);
    oldtree->SetBranchStatus("pTZ1",1);
    oldtree->SetBranchStatus("pTZ2",1);

    oldtree->SetBranchStatus("mass4l",1);
    oldtree->SetBranchStatus("mass2e2mu",1);
    oldtree->SetBranchStatus("mass4mu",1);
    oldtree->SetBranchStatus("mass4e",1);
    oldtree->SetBranchStatus("lep_genindex",1);
    oldtree->SetBranchStatus("lep_Hindex",1);
    oldtree->SetBranchStatus("massZ1",1);
    oldtree->SetBranchStatus("massZ2",1);
    oldtree->SetBranchStatus("pT4l",1);
    oldtree->SetBranchStatus("eta4l",1);
    oldtree->SetBranchStatus("phi4l",1);
    oldtree->SetBranchStatus("rapidity4l",1);
    oldtree->SetBranchStatus("Phi",1);
    oldtree->SetBranchStatus("Phi1",1);
    oldtree->SetBranchStatus("cosTheta1",1);
    oldtree->SetBranchStatus("cosTheta2",1);
    oldtree->SetBranchStatus("cosThetaStar",1);
    oldtree->SetBranchStatus("njets_pt30_eta4p7",1);
    oldtree->SetBranchStatus("njets_pt30_eta2p5",1);
    oldtree->SetBranchStatus("pt_leadingjet_pt30_eta4p7",1);
    oldtree->SetBranchStatus("pt_leadingjet_pt30_eta2p5",1);
    oldtree->SetBranchStatus("njets_pt30_eta4p7_jesup",1);
    oldtree->SetBranchStatus("njets_pt30_eta4p7_jesdn",1);
    oldtree->SetBranchStatus("pt_leadingjet_pt30_eta4p7_jesup",1);
    oldtree->SetBranchStatus("pt_leadingjet_pt30_eta4p7_jesdn",1);
    oldtree->SetBranchStatus("njets_pt30_eta2p5_jesup",1);
    oldtree->SetBranchStatus("njets_pt30_eta2p5_jesdn",1);
    oldtree->SetBranchStatus("pt_leadingjet_pt30_eta2p5_jesup",1);
    oldtree->SetBranchStatus("pt_leadingjet_pt30_eta2p5_jesdn",1);

    if (!_Test) {
        oldtree->SetBranchStatus("D_0m",1);
        oldtree->SetBranchStatus("D_CP",1);
        oldtree->SetBranchStatus("D_0hp",1);
        oldtree->SetBranchStatus("D_int",1);
        oldtree->SetBranchStatus("D_L1",1);
        oldtree->SetBranchStatus("D_L1Zg",1);
        oldtree->SetBranchStatus("TauC_Inc_0j_EnergyWgt",1);
        oldtree->SetBranchStatus("TauB_Inc_0j_pTWgt",1);
        oldtree->SetBranchStatus("pTj2",1);
        oldtree->SetBranchStatus("pTj2_jesup",1);
        oldtree->SetBranchStatus("pTj2_jesdn",1);
        oldtree->SetBranchStatus("mj1j2",1);
        oldtree->SetBranchStatus("mj1j2_jesup",1);
        oldtree->SetBranchStatus("mj1j2_jesdn",1);
        oldtree->SetBranchStatus("dEtaj1j2",1);
        oldtree->SetBranchStatus("dEtaj1j2_jesup",1);
        oldtree->SetBranchStatus("dEtaj1j2_jesdn",1);
        oldtree->SetBranchStatus("dPhij1j2",1);
        oldtree->SetBranchStatus("dPhij1j2_jesup",1);
        oldtree->SetBranchStatus("dPhij1j2_jesdn",1);
        oldtree->SetBranchStatus("pT4lj",1);
        oldtree->SetBranchStatus("pT4lj_jesup",1);
        oldtree->SetBranchStatus("pT4lj_jesdn",1);
        oldtree->SetBranchStatus("mass4lj",1);
        oldtree->SetBranchStatus("mass4lj_jesup",1);
        oldtree->SetBranchStatus("mass4lj_jesdn",1);
        oldtree->SetBranchStatus("pT4ljj",1);
        oldtree->SetBranchStatus("pT4ljj_jesup",1);
        oldtree->SetBranchStatus("pT4ljj_jesdn",1);
        oldtree->SetBranchStatus("mass4ljj",1);
        oldtree->SetBranchStatus("mass4ljj_jesup",1);
        oldtree->SetBranchStatus("mass4ljj_jesdn",1);
    }

    oldtree->SetBranchStatus("GENZ_DaughtersId",1);
    oldtree->SetBranchStatus("GENZ_MomId",1);
    oldtree->SetBranchStatus("GENlep_MomMomId",1);
    oldtree->SetBranchStatus("GENlep_MomId",1);
    oldtree->SetBranchStatus("GENlep_Hindex",1);
    oldtree->SetBranchStatus("GENlep_id",1);

    oldtree->SetBranchStatus("GENmass4l",1);
    oldtree->SetBranchStatus("GENmassZ1",1);
    oldtree->SetBranchStatus("GENmassZ2",1);
    oldtree->SetBranchStatus("GENpT4l",1);
    oldtree->SetBranchStatus("GENeta4l",1);
    oldtree->SetBranchStatus("GENphi4l",1);
    oldtree->SetBranchStatus("GENrapidity4l",1);
    oldtree->SetBranchStatus("GENPhi",1);
    oldtree->SetBranchStatus("GENPhi1",1);
    oldtree->SetBranchStatus("GENcosTheta1",1);
    oldtree->SetBranchStatus("GENcosTheta2",1);
    oldtree->SetBranchStatus("GENcosThetaStar",1);
    oldtree->SetBranchStatus("GENnjets_pt30_eta4p7",1);
    oldtree->SetBranchStatus("GENnjets_pt30_eta2p5",1);
    oldtree->SetBranchStatus("GENpt_leadingjet_pt30_eta4p7",1);
    oldtree->SetBranchStatus("GENpt_leadingjet_pt30_eta2p5",1);

// Weights for the uncertainty module
    oldtree->SetBranchStatus("nnloWeights",1);
    oldtree->SetBranchStatus("qcdWeights",1);
    oldtree->SetBranchStatus("pdfENVup",1);


    if (!_Test) {
        oldtree->SetBranchStatus("GEND_0m",1);
        oldtree->SetBranchStatus("GEND_CP",1);
        oldtree->SetBranchStatus("GEND_0hp",1);
        oldtree->SetBranchStatus("GEND_int",1);
        oldtree->SetBranchStatus("GEND_L1",1);
        oldtree->SetBranchStatus("GEND_L1Zg",1);
        oldtree->SetBranchStatus("GEN_TauC_Inc_0j_EnergyWgt",1);
        oldtree->SetBranchStatus("GEN_TauB_Inc_0j_pTWgt",1);

        oldtree->SetBranchStatus("GENpTj2",1);
        oldtree->SetBranchStatus("GENmj1j2",1);
        oldtree->SetBranchStatus("GENdEtaj1j2",1);
        oldtree->SetBranchStatus("GENdPhij1j2",1);
        oldtree->SetBranchStatus("GENpT4lj",1);
        oldtree->SetBranchStatus("GENmass4lj",1);
        oldtree->SetBranchStatus("GENpT4ljj",1);
        oldtree->SetBranchStatus("GENmass4ljj",1);
    }

    //Create a new file + a clone of old tree in new file
    TFile *newfile = new TFile(
            (outputFile_Slimmed+name_DS+"_slimmed_newMuSF.root").c_str()
            ,"recreate");
    cout<<"Output file: "<<newfile->GetName()<<endl;
    if (isMC)
    {
        newfile->mkdir("Ana");
        newfile->cd("Ana");
    }
    if (isMC){
        for (int ii=0; ii<7; ii++)
        {
            th[ii]->Write();
        }
    }
    TTree *newtree = oldtree->CloneTree(0);
    newtree->Branch("lep_dataMC_new", &lep_dataMC_new);
    newtree->Branch("lep_dataMCErr_new", &lep_dataMCErr_new);
    newtree->Branch("dataMCWeight_new", &dataMCWeight_new);
    newtree->Branch("eventWeight_new", &eventWeight_new);

    std::set<TString> runlumieventSet;
    for (Long64_t i=0;i<nentries; i++) {
        lep_dataMC_new.clear();
        lep_dataMCErr_new.clear();
        //if (i>=2000000) continue;
        //if (i>=2000000&&_Test) break;
        if (i%100000==0) std::cout<<i<<"/"<<nentries<<std::endl;
        oldtree->GetEntry(i);

        if (isMC)
        {
            if (isSignal || (!isSignal && (passedZ4lSelection ||passedZXCRSelection)))
            {
                for(Long64_t j=0; j<lep_pt->size(); j++)
                {
                    if(TMath::Abs((*lep_id)[j])==13)
                    {
                        lep_dataMC_new.push_back(dataMC_((*lep_pt)[j], (*lep_eta)[j], hMuScaleFac));
                        lep_dataMCErr_new.push_back(dataMCErr_((*lep_pt)[j], (*lep_eta)[j], hMuScaleFacUnc));
                    }
                    else
                    {
                        lep_dataMC_new.push_back((*lep_dataMC)[j]);
                        lep_dataMCErr_new.push_back((*lep_dataMCErr)[j]);
                    }
                }

                int sum_lep_Hindex=lep_Hindex[0]+lep_Hindex[1]+lep_Hindex[2]+lep_Hindex[3];//lep_Hindex default value is -1
                if (sum_lep_Hindex!=-4)
                {
                    dataMCWeight_new = lep_dataMC_new[lep_Hindex[0]]*lep_dataMC_new[lep_Hindex[1]]*lep_dataMC_new[lep_Hindex[2]]*lep_dataMC_new[lep_Hindex[3]]; 
                    eventWeight_new = eventWeight/dataMCWeight*dataMCWeight_new;
                }
                else
                {
                    dataMCWeight_new = 1.0;
                    eventWeight_new = eventWeight/dataMCWeight*dataMCWeight_new;
                }
            }
        }

        if (isSignal)
        {
            newtree->Fill();
        }
        if (!isSignal && (passedZ4lSelection ||passedZXCRSelection))
        {
            newtree->Fill();////if only store the event passedFullSelection
        }
        
    }

    newtree->AutoSave();
    newfile->Close();
    oldfile->Close();
    //delete oldfile;
    delete newfile;
    exit(EXIT_FAILURE);
}

