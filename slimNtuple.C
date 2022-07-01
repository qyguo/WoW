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
    //string filename = (pre_name + whichYear[year-2016] + MC__ + name_DS + ".root").c_str();
    //string outputFile_Slimmed = (pre_name + whichYear[year-2016] + "Slimmed_2p5/").c_str(); 
    string filename = "/publicfs/cms/data/hzz/guoqy/newNTuple_UL/ggh_amcatnloFXFX_ntuple/ggh_amcatnloFXFX_gen.root";
    string outputFile_Slimmed = "/publicfs/cms/data/hzz/guoqy/newNTuple_UL/ggh_amcatnloFXFX_ntuple/";

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

    Float_t         GENpT4lj_2p5;
    Float_t         GENmass4lj_2p5;
    Float_t         GENpT4ljj_2p5;
    Float_t         GENmass4ljj_2p5;
    Float_t         pT4lj_2p5;
    Float_t         mass4lj_2p5;
    Float_t         pT4ljj_2p5;
    Float_t         mass4ljj_2p5;
    Float_t         pT4lj_2p5_jesup;
    Float_t         mass4lj_2p5_jesup;
    Float_t         pT4ljj_2p5_jesup;
    Float_t         mass4ljj_2p5_jesup;
    Float_t         pT4lj_2p5_jesdn;
    Float_t         mass4lj_2p5_jesdn;
    Float_t         pT4ljj_2p5_jesdn;
    Float_t         mass4ljj_2p5_jesdn;

    Float_t         GENpT4l;
    Float_t         GENeta4l;
    Float_t         GENphi4l;
    Float_t         GENmass4l;
    Float_t         GENpTj1_2p5;
    Float_t         GENetaj1_2p5;
    Float_t         GENphij1_2p5;
    Float_t         GENmj1_2p5;
    Float_t         GENpTj2_2p5;
    Float_t         GENetaj2_2p5;
    Float_t         GENphij2_2p5;
    Float_t         GENmj2_2p5;
    Int_t           GENnjets_pt30_eta2p5;
    Float_t         pT4l;
    Float_t         eta4l;
    Float_t         phi4l;
    Float_t         mass4l;
    Float_t         pTj1_2p5;
    Float_t         etaj1_2p5;
    Float_t         phij1_2p5;
    Float_t         mj1_2p5;
    Float_t         pTj2_2p5;
    Float_t         etaj2_2p5;
    Float_t         phij2_2p5;
    Float_t         mj2_2p5;
    Float_t         pTj1_2p5_jesup;
    Float_t         etaj1_2p5_jesup;
    Float_t         phij1_2p5_jesup;
    Float_t         mj1_2p5_jesup;
    Float_t         pTj2_2p5_jesup;
    Float_t         etaj2_2p5_jesup;
    Float_t         phij2_2p5_jesup;
    Float_t         mj2_2p5_jesup;
    Float_t         pTj1_2p5_jesdn;
    Float_t         etaj1_2p5_jesdn;
    Float_t         phij1_2p5_jesdn;
    Float_t         mj1_2p5_jesdn;
    Float_t         pTj2_2p5_jesdn;
    Float_t         etaj2_2p5_jesdn;
    Float_t         phij2_2p5_jesdn;
    Float_t         mj2_2p5_jesdn;
    Int_t           njets_pt30_eta2p5;
    Int_t           njets_pt30_eta2p5_jesup;
    Int_t           njets_pt30_eta2p5_jesdn;

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
        oldtree->SetBranchAddress("GENpT4l", &GENpT4l);
        oldtree->SetBranchAddress("GENeta4l", &GENeta4l);
        oldtree->SetBranchAddress("GENphi4l", &GENphi4l);
        oldtree->SetBranchAddress("GENmass4l", &GENmass4l);
        oldtree->SetBranchAddress("GENpTj1_2p5", &GENpTj1_2p5);
        oldtree->SetBranchAddress("GENetaj1_2p5", &GENetaj1_2p5);
        oldtree->SetBranchAddress("GENphij1_2p5", &GENphij1_2p5);
        oldtree->SetBranchAddress("GENmj1_2p5", &GENmj1_2p5);
        oldtree->SetBranchAddress("GENpTj2_2p5", &GENpTj2_2p5);
        oldtree->SetBranchAddress("GENetaj2_2p5", &GENetaj2_2p5);
        oldtree->SetBranchAddress("GENphij2_2p5", &GENphij2_2p5);
        oldtree->SetBranchAddress("GENmj2_2p5", &GENmj2_2p5);
        oldtree->SetBranchAddress("GENnjets_pt30_eta2p5", &GENnjets_pt30_eta2p5);
    }
    oldtree->SetBranchAddress("pTj1_2p5", &pTj1_2p5);
    oldtree->SetBranchAddress("etaj1_2p5", &etaj1_2p5);
    oldtree->SetBranchAddress("phij1_2p5", &phij1_2p5);
    oldtree->SetBranchAddress("mj1_2p5", &mj1_2p5);
    oldtree->SetBranchAddress("pTj2_2p5", &pTj2_2p5);
    oldtree->SetBranchAddress("etaj2_2p5", &etaj2_2p5);
    oldtree->SetBranchAddress("phij2_2p5", &phij2_2p5);
    oldtree->SetBranchAddress("mj2_2p5", &mj2_2p5);
    oldtree->SetBranchAddress("pTj1_2p5_jesup", &pTj1_2p5_jesup);
    oldtree->SetBranchAddress("etaj1_2p5_jesup", &etaj1_2p5_jesup);
    oldtree->SetBranchAddress("phij1_2p5_jesup", &phij1_2p5_jesup);
    oldtree->SetBranchAddress("mj1_2p5_jesup", &mj1_2p5_jesup);
    oldtree->SetBranchAddress("pTj2_2p5_jesup", &pTj2_2p5_jesup);
    oldtree->SetBranchAddress("etaj2_2p5_jesup", &etaj2_2p5_jesup);
    oldtree->SetBranchAddress("phij2_2p5_jesup", &phij2_2p5_jesup);
    oldtree->SetBranchAddress("mj2_2p5_jesup", &mj2_2p5_jesup);
    oldtree->SetBranchAddress("pTj1_2p5_jesdn", &pTj1_2p5_jesdn);
    oldtree->SetBranchAddress("etaj1_2p5_jesdn", &etaj1_2p5_jesdn);
    oldtree->SetBranchAddress("phij1_2p5_jesdn", &phij1_2p5_jesdn);
    oldtree->SetBranchAddress("mj1_2p5_jesdn", &mj1_2p5_jesdn);
    oldtree->SetBranchAddress("pTj2_2p5_jesdn", &pTj2_2p5_jesdn);
    oldtree->SetBranchAddress("etaj2_2p5_jesdn", &etaj2_2p5_jesdn);
    oldtree->SetBranchAddress("phij2_2p5_jesdn", &phij2_2p5_jesdn);
    oldtree->SetBranchAddress("mj2_2p5_jesdn", &mj2_2p5_jesdn);
    oldtree->SetBranchAddress("pT4l", &pT4l);
    oldtree->SetBranchAddress("eta4l", &eta4l);
    oldtree->SetBranchAddress("phi4l", &phi4l);
    oldtree->SetBranchAddress("mass4l", &mass4l);
    oldtree->SetBranchAddress("njets_pt30_eta2p5", &njets_pt30_eta2p5);
    oldtree->SetBranchAddress("njets_pt30_eta2p5_jesup", &njets_pt30_eta2p5_jesup);
    oldtree->SetBranchAddress("njets_pt30_eta2p5_jesdn", &njets_pt30_eta2p5_jesdn);
    
    Bool_t passedZ4lSelection;
    Bool_t passedZXCRSelection;
    oldtree->SetBranchAddress("passedZ4lSelection", &passedZ4lSelection);
    oldtree->SetBranchAddress("passedZXCRSelection", &passedZXCRSelection);

    oldtree->SetBranchStatus("*",0); //Disables All Branches
    //Then enables only select branches
    //
    //
    oldtree->SetBranchStatus("lep_id",1);
    oldtree->SetBranchStatus("lep_pt",1);
    oldtree->SetBranchStatus("lep_eta",1);
    oldtree->SetBranchStatus("lep_dataMC",1);
    oldtree->SetBranchStatus("lep_dataMCErr",1);

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

    oldtree->SetBranchStatus("met",1);
    oldtree->SetBranchStatus("lep_phi",1);
    oldtree->SetBranchStatus("lep_mass",1);
    oldtree->SetBranchStatus("lep_tightId",1);
    oldtree->SetBranchStatus("lep_RelIsoNoFSR",1);
    oldtree->SetBranchStatus("lep_matchedR03_PdgId",1);
    oldtree->SetBranchStatus("lep_matchedR03_MomId",1);
    oldtree->SetBranchStatus("lep_matchedR03_MomMomId",1);

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
    oldtree->SetBranchStatus("mass4lErr",1);
    oldtree->SetBranchStatus("mass4lErrREFIT",1);

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
//
        oldtree->SetBranchStatus("pTj2_2p5",1);
        oldtree->SetBranchStatus("pTj2_2p5_jesup",1);
        oldtree->SetBranchStatus("pTj2_2p5_jesdn",1);
        oldtree->SetBranchStatus("mj1j2_2p5",1);
        oldtree->SetBranchStatus("mj1j2_2p5_jesup",1);
        oldtree->SetBranchStatus("mj1j2_2p5_jesdn",1);
        oldtree->SetBranchStatus("dEtaj1j2_2p5",1);
        oldtree->SetBranchStatus("dEtaj1j2_2p5_jesup",1);
        oldtree->SetBranchStatus("dEtaj1j2_2p5_jesdn",1);
        oldtree->SetBranchStatus("dPhij1j2_2p5",1);
        oldtree->SetBranchStatus("dPhij1j2_2p5_jesup",1);
        oldtree->SetBranchStatus("dPhij1j2_2p5_jesdn",1);
        oldtree->SetBranchStatus("pTj1_2p5",1);
        oldtree->SetBranchStatus("etaj1_2p5",1);
        oldtree->SetBranchStatus("phij1_2p5",1);
        oldtree->SetBranchStatus("mj1_2p5",1);
        //oldtree->SetBranchStatus("pTj2_2p5",1);
        oldtree->SetBranchStatus("etaj2_2p5",1);
        oldtree->SetBranchStatus("phij2_2p5",1);
        oldtree->SetBranchStatus("mj2_2p5",1);
        oldtree->SetBranchStatus("pTj1_2p5_jesup",1);
        oldtree->SetBranchStatus("etaj1_2p5_jesup",1);
        oldtree->SetBranchStatus("phij1_2p5_jesup",1);
        oldtree->SetBranchStatus("mj1_2p5_jesup",1);
        //oldtree->SetBranchStatus("pTj2_2p5_jesup",1);
        oldtree->SetBranchStatus("etaj2_2p5_jesup",1);
        oldtree->SetBranchStatus("phij2_2p5_jesup",1);
        oldtree->SetBranchStatus("mj2_2p5_jesup",1);
        oldtree->SetBranchStatus("pTj1_2p5_jesdn",1);
        oldtree->SetBranchStatus("etaj1_2p5_jesdn",1);
        oldtree->SetBranchStatus("phij1_2p5_jesdn",1);
        oldtree->SetBranchStatus("mj1_2p5_jesdn",1);
        //oldtree->SetBranchStatus("pTj2_2p5_jesdn",1);
        oldtree->SetBranchStatus("etaj2_2p5_jesdn",1);
        oldtree->SetBranchStatus("phij2_2p5_jesdn",1);
        oldtree->SetBranchStatus("mj2_2p5_jesdn",1);
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
        oldtree->SetBranchStatus("GENpTj2_2p5",1);
        oldtree->SetBranchStatus("GENmj1j2_2p5",1);
        oldtree->SetBranchStatus("GENdEtaj1j2_2p5",1);
        oldtree->SetBranchStatus("GENdPhij1j2_2p5",1);

        oldtree->SetBranchStatus("GENpTj1_2p5",1);
        oldtree->SetBranchStatus("GENetaj1_2p5",1);
        oldtree->SetBranchStatus("GENphij1_2p5",1);
        oldtree->SetBranchStatus("GENmj1_2p5",1);
        //oldtree->SetBranchStatus("GENpTj2_2p5",1);
        oldtree->SetBranchStatus("GENetaj2_2p5",1);
        oldtree->SetBranchStatus("GENphij2_2p5",1);
        oldtree->SetBranchStatus("GENmj2_2p5",1);
    }

    //Create a new file + a clone of old tree in new file
    TFile *newfile = new TFile(
            (outputFile_Slimmed+name_DS+"_slimmed_newMuSF_add2p5.root").c_str()
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
    newtree->Branch("GENpT4lj_2p5", &GENpT4lj_2p5);
    newtree->Branch("GENmass4lj_2p5", &GENmass4lj_2p5);
    newtree->Branch("GENpT4ljj_2p5", &GENpT4ljj_2p5);
    newtree->Branch("GENmass4ljj_2p5", &GENmass4ljj_2p5);
    newtree->Branch("pT4lj_2p5", &pT4lj_2p5);
    newtree->Branch("mass4lj_2p5", &mass4lj_2p5);
    newtree->Branch("pT4ljj_2p5", &pT4ljj_2p5);
    newtree->Branch("mass4ljj_2p5", &mass4ljj_2p5);
    newtree->Branch("pT4lj_2p5_jesup", &pT4lj_2p5_jesup);
    newtree->Branch("mass4lj_2p5_jesup", &mass4lj_2p5_jesup);
    newtree->Branch("pT4ljj_2p5_jesup", &pT4ljj_2p5_jesup);
    newtree->Branch("mass4ljj_2p5_jesup", &mass4ljj_2p5_jesup);
    newtree->Branch("pT4lj_2p5_jesdn", &pT4lj_2p5_jesdn);
    newtree->Branch("mass4lj_2p5_jesdn", &mass4lj_2p5_jesdn);
    newtree->Branch("pT4ljj_2p5_jesdn", &pT4ljj_2p5_jesdn);
    newtree->Branch("mass4ljj_2p5_jesdn", &mass4ljj_2p5_jesdn);

    std::set<TString> runlumieventSet;
    for (Long64_t i=0;i<nentries; i++) {
        lep_dataMC_new.clear();
        lep_dataMCErr_new.clear();
        eventWeight_new=-1.0;
        dataMCWeight_new=-1.0;
        GENpT4lj_2p5=-1.0;
        GENmass4lj_2p5=-1.0;
        GENpT4ljj_2p5=-1.0;
        GENmass4ljj_2p5=-1.0;
        pT4lj_2p5=-1.0;
        mass4lj_2p5=-1.0;
        pT4ljj_2p5=-1.0;
        mass4ljj_2p5=-1.0;
        pT4lj_2p5_jesup=-1.0;
        mass4lj_2p5_jesup=-1.0;
        pT4ljj_2p5_jesup=-1.0;
        mass4ljj_2p5_jesup=-1.0;
        pT4lj_2p5_jesdn=-1.0;
        mass4lj_2p5_jesdn=-1.0;
        pT4ljj_2p5_jesdn=-1.0;
        mass4ljj_2p5_jesdn=-1.0;
        //if (i>=2000000) continue;
        //if (i>=2000000&&_Test) break;
        if (i%100000==0) std::cout<<i<<"/"<<nentries<<std::endl;
        oldtree->GetEntry(i);

        if (!isMC || (isMC && (isSignal || (!isSignal && ( passedZXCRSelection||passedZ1LSelection )))) )
        {
            if (pT4l>0)
            {
                TLorentzVector j1, j2, Higgs, hj, jj, hjj, j1_jesup, j1_jesdn, j2_jesup, j2_jesdn, hj_jesup, hj_jesdn, jj_jesup, jj_jesdn, hjj_jesup, hjj_jesdn;
                Higgs.SetPtEtaPhiM(pT4l, eta4l, phi4l, mass4l);
                if(njets_pt30_eta2p5>0)
                {
                    j1.SetPtEtaPhiM(pTj1_2p5,etaj1_2p5,phij1_2p5,mj1_2p5);
                    hj = Higgs + j1;
                    pT4lj_2p5 = hj.Pt();
                    mass4lj_2p5 = hj.M();
                }
                if(njets_pt30_eta2p5>1)
                {
                    j2.SetPtEtaPhiM(pTj2_2p5,etaj2_2p5,phij2_2p5,mj2_2p5);
                    jj = j1 + j2;
                    hjj = Higgs + jj;
                    pT4ljj_2p5 = hjj.Pt();
                    mass4ljj_2p5 = hjj.M();
                }
                if(njets_pt30_eta2p5_jesup>0)
                {
                    j1_jesup.SetPtEtaPhiM(pTj1_2p5_jesup,etaj1_2p5_jesup,phij1_2p5_jesup,mj1_2p5_jesup);
                    hj_jesup = Higgs + j1_jesup;
                    pT4lj_2p5_jesup = hj_jesup.Pt();
                    mass4lj_2p5_jesup = hj_jesup.M();
                }
                if(njets_pt30_eta2p5_jesup>1)
                {
                    j2_jesup.SetPtEtaPhiM(pTj2_2p5_jesup,etaj2_2p5_jesup,phij2_2p5_jesup,mj2_2p5_jesup);
                    jj_jesup = j1_jesup + j2_jesup;
                    hjj_jesup = Higgs + jj_jesup;
                    pT4ljj_2p5_jesup = hjj_jesup.Pt();
                    mass4ljj_2p5_jesup = hjj_jesup.M();
                }
                if(njets_pt30_eta2p5_jesdn>0)
                {
                    j1_jesdn.SetPtEtaPhiM(pTj1_2p5_jesdn,etaj1_2p5_jesdn,phij1_2p5_jesdn,mj1_2p5_jesdn);
                    hj_jesdn = Higgs + j1_jesdn;
                    pT4lj_2p5_jesdn = hj_jesdn.Pt();
                    mass4lj_2p5_jesdn = hj_jesdn.M();
                }
                if(njets_pt30_eta2p5_jesdn>1)
                {
                    j2_jesdn.SetPtEtaPhiM(pTj2_2p5_jesdn,etaj2_2p5_jesdn,phij2_2p5_jesdn,mj2_2p5_jesdn);
                    jj_jesdn = j1_jesdn + j2_jesdn;
                    hjj_jesdn = Higgs + jj_jesdn;
                    pT4ljj_2p5_jesdn = hjj_jesdn.Pt();
                    mass4ljj_2p5_jesdn = hjj_jesdn.M();
                }
            }
        }

        if (isMC)
        {
            if (isSignal || (!isSignal && (passedZXCRSelection || passedZ1LSelection)))
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

                if(GENpT4l>0)
                {
                    TLorentzVector GENj1, GENj2, GENHiggs, GENhj, GENjj, GENhjj;
                    GENHiggs.SetPtEtaPhiM(GENpT4l, GENeta4l, GENphi4l, GENmass4l);

                    if(GENnjets_pt30_eta2p5>0)
                    {
                        GENj1.SetPtEtaPhiM(GENpTj1_2p5, GENetaj1_2p5, GENphij1_2p5, GENmj1_2p5);
                        GENhj = GENj1 + GENHiggs;
                        GENpT4lj_2p5=GENhj.Pt();
                        GENmass4lj_2p5=GENhj.M();
                    }
                    if(GENnjets_pt30_eta2p5>1)
                    {
                        GENj2.SetPtEtaPhiM(GENpTj2_2p5, GENetaj2_2p5, GENphij2_2p5, GENmj2_2p5);
                        GENjj = GENj1 + GENj2;
                        GENhjj = GENjj + GENHiggs;
                        GENpT4ljj_2p5 = GENhjj.Pt();
                        GENmass4ljj_2p5 = GENhjj.M();
                    }
                }
            }
        }

        if (isSignal)
        {
            newtree->Fill();
        }
        if (!isSignal && (passedZ4lSelection || passedZXCRSelection || passedZ1LSelection))
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

