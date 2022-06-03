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
#include "TSpline.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TLorentzVector.h"
#include "Math/VectorUtil.h"
#include "TClonesArray.h"
#include "TCanvas.h"
#include "slimNtuple_JES.h"
#include "variables.h"
//#include "JetMETCorrections/Modules/interface/JetResolution.h"
//#include "DataFormats/PatCandidates/interface/Jet.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"
#include "DataFormats/Math/interface/deltaPhi.h"
// Jet energy correction
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"
#include "JetMETCorrections/Modules/interface/JetResolution.h"
// NJettiness module
//#include "UFHZZAnalysisRun2/UFHZZ4LAna/interface/NJettiness.h"

using namespace edm;
using namespace std;



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

float TauC(vector<float> pt, vector<float> eta, vector<float> phi, vector<float> mass, TLorentzVector H)
{
	float TauC_j = 0; float TauC_jmax = 0;
	for( unsigned int k = 0; k< pt.size(); k++) {
	    TLorentzVector theJet;
	    theJet.SetPtEtaPhiM(pt[k],eta[k],phi[k],mass[k]);
	    TauC_j = sqrt(theJet.Pt()*theJet.Pt() + theJet.M()*theJet.M())/(2*cosh(theJet.Rapidity() - H.Rapidity()));
	    if (TauC_j > TauC_jmax) {
		TauC_jmax = TauC_j;
		}
	}
    return TauC_jmax;
}

float TauB(vector<float> pt, vector<float> eta, vector<float> phi, vector<float> mass, TLorentzVector H)
{
        float TauB_j = 0; float TauB_jmax = 0;
        for( unsigned int k = 0; k< pt.size(); k++) {
            TLorentzVector theJet;
            theJet.SetPtEtaPhiM(pt[k],eta[k],phi[k],mass[k]);
	    TauB_j = sqrt(theJet.Pt()*theJet.Pt() + theJet.M()*theJet.M())*exp(-1*(abs(theJet.Rapidity() - H.Rapidity())));   // equivalent ?
	    //TauB_j = sqrt(theJet.Pt()*theJet.Pt() + theJet.M()*theJet.M())*exp(-1*(theJet.Rapidity() - H.Rapidity()));   // equivalent ?
	    //TauB_j = sqrt(theJet.energy()*theJet.energy() - theJet.pz()*theJet.pz())*exp(-1*(theJet.Rapidity() - H.Rapidity()));
            if (TauB_j > TauB_jmax) {
                TauB_jmax = TauB_j;
                }
        }
    return TauB_jmax;
}



//void slimNtuple_JES(const int & _year_=2016, const string & _name_DS_="GluGluHToZZTo4L_M125_TuneCP5_13TeV_powheg2_JHUGenV7011_pythia8", const bool & isMC = true,  const bool & isSignal = true, const bool & _Test=false) {
void slimNtuple_JES(const int & _year_=2016, const string & _name_DS_="ttH_HToZZ_4LFilter_M124_TuneCP5_13TeV_powheg2_JHUGenV7011_pythia8", const bool & isMC = true,  const bool & isSignal = true, const bool & _Test=false, const bool & applyJEC_ = true) {
    std::cout<<"year: "<<_year_<<"; name_DS: "<<_name_DS_<<"; isMC: "<<isMC<<"; isSignal: "<<isSignal<<"; _Test: "<<_Test<<std::endl;

    //gROOT->ProcessLine("gSystem->Load(\"AddIncludePath(\"-I$CMSSW_BASE/src/WoW/JECUncertaintySources/\")");
    gROOT->ProcessLine("gSystem->Load(\"libRecoJetsJetProducers.so\")");
    gROOT->ProcessLine("gSystem->Load(\"libSimDataFormatsHTXS.so\")");
    gROOT->ProcessLine("gSystem->Load(\"libPhysicsToolsPatAlgos.so\")");
    gROOT->ProcessLine("gSystem->Load(\"libCondFormatsJetMETObjects.so\")");
    gROOT->ProcessLine("gSystem->Load(\"libFWCoreParameterSet.so\")");
 
    int year = _year_;
    int Year = TMath::Abs(year);
    //string pre_name = "/publicfs/cms/data/hzz/guoqy/newNTuple_UL/" + whichYear[year-2016] + "/MC/";
    //string pre_name = "/publicfs/cms/data/hzz/guoqy/newNTuple_UL/";
    //string pre_name = "/publicfs/cms/data/hzz/jtahir/UL2016/MC/preVFP/GluGluHToZZTo4L_M125_TuneCP5_13TeV_powheg2_JHUGenV7011_pythia8/";
    string pre_name1 = "/publicfs/cms/data/hzz/jtahir/UL2016/MC/";
    string pre_name = (pre_name1 + _name_DS_).c_str();
    //string pre_name = "/afs/cern.ch/user/t/tjavaid/workspace/BUF_ntuplizer/UL/fresh/CMSSW_10_6_26/src/";
//    if (year==2016)    pre_name="";
//    else if (year==-2016) //UL16APV
//        pre_name="";

    string whichYear[3] = {"2016/", "2017/", "2018/"};
    string MC__ = "MC/";
    if (!isMC) MC__ = "Data/";
    //string name_DS = "GluGluHToZZTo4L_M125_2017";
    string name_DS = _name_DS_.c_str();

    //TString filename = prefix+".root";
    //string filename = (pre_name + whichYear[year-2016] + MC__ + name_DS + ".root").c_str();
//    string filename = (pre_name + name_DS + ".root").c_str();
    string filename = (pre_name + "/" + name_DS + ".root").c_str();
    //string outputFile_Slimmed = (pre_name + whichYear[year-2016] + "Slimmed_2p5/").c_str(); 
    string outputFile_Slimmed = ("Slimmed_2p5/"); //.c_str(); 

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
    ULong64_t Run, LumiSect, Event;
//    Bool_t          passedFullSelection;
    Bool_t          isH4l;
    vector<int> *lep_genindex;
    vector<int> *GENlep_MomId;
    vector<int> *GENlep_MomMomId;
    //int GENlep_Hindex[4];
//    int lep_Hindex[4];

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
    vector<float> pt_nom {};
// JES split vars
/*    vector<float> jes_unc_split {};
    vector<float> pt_jesup_split {};
    vector<float> pt_jesdn_split {};
    float singleContr_jes_unc = 0; 
*/
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
    Float_t         GENmassZZ;
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
    Float_t         rapidity4l;
    Float_t         phi4l;
    Float_t         mass4l;
    Float_t         pt_leadingjet_pt30_eta4p7;
    Float_t         TauC_Inc_0j_EnergyWgt;
    Float_t         TauB_Inc_0j_pTWgt;
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
    float pTj1;
    Float_t         TauC_Inc_0j_EnergyWgt_nom, TauB_Inc_0j_pTWgt_nom;


// try
/*    std::vector<> *jet_mass;
    std::vector<> *jet_pt; 
    std::vector<float> *jet_eta; 
    std::vector<float> *jet_phi;
    std::vector<float> *jet_iscleanH4l; 
*/
    vector<float> *jet_mass;
    vector<float> *jet_pt; 
    vector<float> *jet_eta; 
    vector<float> *jet_phi;
    vector<float> *jet_iscleanH4l; 

    vector<float> jes_unc_split_Total {}; 
    vector<float> jes_unc_split_Abs {}; 
    vector<float> jes_unc_split_Abs_year {}; 
    vector<float> jes_unc_split_BBEC1 {}; 
    vector<float> jes_unc_split_BBEC1_year {};
    vector<float> jes_unc_split_EC2 {}; 
    vector<float> jes_unc_split_EC2_year {}; 
    vector<float> jes_unc_split_FlavQCD {}; 
    vector<float> jes_unc_split_HF {}; 
    vector<float> jes_unc_split_HF_year {};
    vector<float> jes_unc_split_RelBal {}; 
    vector<float> jes_unc_split_RelSample_year {}; 
    // up
    vector<float> pt_jesup_split_Total {}; 
    vector<float> pt_jesup_split_Abs {}; 
    vector<float> pt_jesup_split_Abs_year {}; 
    vector<float> pt_jesup_split_BBEC1 {}; 
    vector<float> pt_jesup_split_BBEC1_year {};
    vector<float> pt_jesup_split_EC2 {}; 
    vector<float> pt_jesup_split_EC2_year {}; 
    vector<float> pt_jesup_split_FlavQCD {}; 
    vector<float> pt_jesup_split_HF {}; 
    vector<float> pt_jesup_split_HF_year {};
    vector<float> pt_jesup_split_RelBal {}; 
    vector<float> pt_jesup_split_RelSample_year {}; 
    //dn
    vector<float> pt_jesdn_split_Total {}; 
    vector<float> pt_jesdn_split_Abs {}; 
    vector<float> pt_jesdn_split_Abs_year {}; 
    vector<float> pt_jesdn_split_BBEC1 {}; 
    vector<float> pt_jesdn_split_BBEC1_year {};
    vector<float> pt_jesdn_split_EC2 {}; 
    vector<float> pt_jesdn_split_EC2_year {}; 
    vector<float> pt_jesdn_split_FlavQCD {}; 
    vector<float> pt_jesdn_split_HF {}; 
    vector<float> pt_jesdn_split_HF_year {};
    vector<float> pt_jesdn_split_RelBal {}; 
    vector<float> pt_jesdn_split_RelSample_year {};

    vector<float> eta_jes_split_Total {};
    vector<float> eta_jes_split_Abs {};
    vector<float> eta_jes_split_Abs_year {};
    vector<float> eta_jes_split_BBEC1 {};
    vector<float> eta_jes_split_BBEC1_year {};
    vector<float> eta_jes_split_EC2 {};
    vector<float> eta_jes_split_EC2_year {};
    vector<float> eta_jes_split_FlavQCD {};
    vector<float> eta_jes_split_HF {};
    vector<float> eta_jes_split_HF_year {};
    vector<float> eta_jes_split_RelBal {};
    vector<float> eta_jes_split_RelSample_year {};

    vector<float> phi_jes_split_Total {};
    vector<float> phi_jes_split_Abs {};
    vector<float> phi_jes_split_Abs_year {};
    vector<float> phi_jes_split_BBEC1 {};
    vector<float> phi_jes_split_BBEC1_year {};
    vector<float> phi_jes_split_EC2 {};
    vector<float> phi_jes_split_EC2_year {};
    vector<float> phi_jes_split_FlavQCD {};
    vector<float> phi_jes_split_HF {};
    vector<float> phi_jes_split_HF_year {};
    vector<float> phi_jes_split_RelBal {};
    vector<float> phi_jes_split_RelSample_year {};

    vector<float> mass_jes_split_Total {};
    vector<float> mass_jes_split_Abs {};
    vector<float> mass_jes_split_Abs_year {};
    vector<float> mass_jes_split_BBEC1 {};
    vector<float> mass_jes_split_BBEC1_year {};
    vector<float> mass_jes_split_EC2 {};
    vector<float> mass_jes_split_EC2_year {};
    vector<float> mass_jes_split_FlavQCD {};
    vector<float> mass_jes_split_HF {};
    vector<float> mass_jes_split_HF_year {};
    vector<float> mass_jes_split_RelBal {};
    vector<float> mass_jes_split_RelSample_year {};


    lep_id = 0; 
    lep_pt = 0;    
    lep_eta = 0;  
    lep_dataMC = 0;
    lep_dataMCErr = 0;
    jet_pt = 0; 
    jet_eta = 0; 
    jet_phi = 0; 
    jet_mass = 0; 
    jet_iscleanH4l = 0;
    lep_genindex = 0;
    GENlep_MomId = 0;
    GENlep_MomMomId = 0;

 
    std::string jecUncFile_;
    std::vector<string> uncSources {};
    std::vector<JetCorrectionUncertainty*> splittedUncerts_;

    bool doM4lSkim = true;
    float GENmass4l_lo = 105.0;
    float mass4l_lo = 105.0;
    float GENmass4l_hi = 160.0; 
    float mass4l_hi = 160.0; 


// JEC split
   bool redoJets=true;
   if (redoJets) {
   if (year == 2016) {
         edm::FileInPath jecUncFile("JECUncertaintySources/RegroupedV2_Summer16_07Aug2017_V11_MC_UncertaintySources_AK4PFchs.txt");
         jecUncFile_ = jecUncFile.fullPath();
         uncSources.push_back("Total");
         uncSources.push_back("Absolute"); uncSources.push_back("Absolute_2016");
         uncSources.push_back("BBEC1"); uncSources.push_back("BBEC1_2016");
         uncSources.push_back("EC2"); uncSources.push_back("EC2_2016");
         uncSources.push_back("FlavorQCD");
         uncSources.push_back("HF"); uncSources.push_back("HF_2016");
         uncSources.push_back("RelativeBal");
         uncSources.push_back("RelativeSample_2016");
    }
   if (year == 2017)   {
         edm::FileInPath jecUncFile("JECUncertaintySources/RegroupedV2_Fall17_17Nov2017_V32_MC_UncertaintySources_AK4PFchs.txt");
         jecUncFile_ = jecUncFile.fullPath();
         uncSources.push_back("Total");
         uncSources.push_back("Absolute"); uncSources.push_back("Absolute_2017");
         uncSources.push_back("BBEC1"); uncSources.push_back("BBEC1_2017");
         uncSources.push_back("EC2"); uncSources.push_back("EC2_2017");
         uncSources.push_back("FlavorQCD");
         uncSources.push_back("HF"); uncSources.push_back("HF_2017");
         uncSources.push_back("RelativeBal");
         uncSources.push_back("RelativeSample_2017");
         }
   if (year == 2018)   {
         edm::FileInPath jecUncFile("JECUncertaintySources/RegroupedV2_Autumn18_V19_MC_UncertaintySources_AK4PFchs.txt");
         jecUncFile_ = jecUncFile.fullPath();
         uncSources.push_back("Total");
         uncSources.push_back("Absolute"); uncSources.push_back("Absolute_2018");
         uncSources.push_back("BBEC1"); uncSources.push_back("BBEC1_2018");
         uncSources.push_back("EC2"); uncSources.push_back("EC2_2018");
         uncSources.push_back("FlavorQCD");
         uncSources.push_back("HF"); uncSources.push_back("HF_2018");
         uncSources.push_back("RelativeBal");
         uncSources.push_back("RelativeSample_2018");
         }
    }
    if(applyJEC_ && isMC)
    {
      for (unsigned s_unc = 0; s_unc < uncSources.size(); s_unc++)
	{
	  JetCorrectorParameters corrParams = JetCorrectorParameters(jecUncFile_, uncSources[s_unc]);
	  splittedUncerts_.push_back(new JetCorrectionUncertainty(corrParams));
          //cout<<"splittedUncerts_[0]      "<<splittedUncerts_[0]<<endl;
	}
    }

   

        Bool_t passedFullSelection;
        oldtree->SetBranchAddress("passedFullSelection", &passedFullSelection);
        oldtree->SetBranchAddress("lep_genindex", &lep_genindex);
        oldtree->SetBranchAddress("GENlep_MomId", &GENlep_MomId);
        oldtree->SetBranchAddress("GENlep_MomMomId", &GENlep_MomMomId);
	oldtree->SetBranchAddress("lep_Hindex", lep_Hindex);

    if (isMC)
    {
        //oldtree->SetBranchAddress("Run",&Run);
        //oldtree->SetBranchAddress("LumiSect",&LumiSect);
        oldtree->SetBranchAddress("Event",&Event);
        oldtree->SetBranchAddress("dataMCWeight", &dataMCWeight);
        oldtree->SetBranchAddress("eventWeight", &eventWeight);
        oldtree->SetBranchAddress("lep_id", &lep_id);
        oldtree->SetBranchAddress("lep_pt", &lep_pt);
//        oldtree->SetBranchAddress("lep_Hindex", lep_Hindex);
        oldtree->SetBranchAddress("lep_eta", &lep_eta);
        oldtree->SetBranchAddress("lep_dataMC", &lep_dataMC);
        oldtree->SetBranchAddress("lep_dataMCErr", &lep_dataMCErr);
        oldtree->SetBranchAddress("GENpT4l", &GENpT4l);
        oldtree->SetBranchAddress("GENmassZZ", &GENmassZZ);
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
/*        oldtree->SetBranchAddress("lep_genindex", &lep_genindex);
        oldtree->SetBranchAddress("GENlep_MomId", &GENlep_MomId);
        oldtree->SetBranchAddress("GENlep_MomMomId", &GENlep_MomMomId);
        oldtree->SetBranchAddress("jet_pt", &jet_pt);
        oldtree->SetBranchAddress("jet_eta", &jet_eta);
        oldtree->SetBranchAddress("jet_phi", &jet_phi);
        oldtree->SetBranchAddress("jet_mass", &jet_mass);
        oldtree->SetBranchAddress("jet_iscleanH4l", &jet_iscleanH4l); */


    }
    oldtree->SetBranchAddress("TauC_Inc_0j_EnergyWgt", &TauC_Inc_0j_EnergyWgt);
    oldtree->SetBranchAddress("TauB_Inc_0j_pTWgt", &TauB_Inc_0j_pTWgt);
    oldtree->SetBranchAddress("pt_leadingjet_pt30_eta4p7", &pt_leadingjet_pt30_eta4p7);
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
    oldtree->SetBranchAddress("rapidity4l", &rapidity4l);
    oldtree->SetBranchAddress("phi4l", &phi4l);
    oldtree->SetBranchAddress("mass4l", &mass4l);
    oldtree->SetBranchAddress("njets_pt30_eta2p5", &njets_pt30_eta2p5);
    oldtree->SetBranchAddress("njets_pt30_eta2p5_jesup", &njets_pt30_eta2p5_jesup);
    oldtree->SetBranchAddress("njets_pt30_eta2p5_jesdn", &njets_pt30_eta2p5_jesdn);
    oldtree->SetBranchAddress("jet_pt", &jet_pt);
    oldtree->SetBranchAddress("jet_eta", &jet_eta);
    oldtree->SetBranchAddress("jet_phi", &jet_phi);
    oldtree->SetBranchAddress("jet_mass", &jet_mass);
    oldtree->SetBranchAddress("jet_iscleanH4l", &jet_iscleanH4l); 

    
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


    oldtree->SetBranchStatus("jet_pt",1);
    oldtree->SetBranchStatus("jet_eta",1);
    oldtree->SetBranchStatus("jet_phi",1);
    oldtree->SetBranchStatus("jet_mass",1);
    oldtree->SetBranchStatus("jet_iscleanH4l",1);

    oldtree->SetBranchStatus("Run",1);
    oldtree->SetBranchStatus("LumiSect",1);
    oldtree->SetBranchStatus("Event",1);
    oldtree->SetBranchStatus("passedFiducialSelection",1);

    oldtree->SetBranchStatus("passedZ4lSelection",1);
    oldtree->SetBranchStatus("passedZ1LSelection",1);
    oldtree->SetBranchStatus("passedZ4lZ1LSelection",1);
    oldtree->SetBranchStatus("passedZ4lZXCRSelection",1);
    oldtree->SetBranchStatus("passedZXCRSelection",1);
    oldtree->SetBranchStatus("nZXCRFailedLeptons",1);
    oldtree->SetBranchStatus("passedFullSelection",1);

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
    oldtree->SetBranchStatus("GENlep_MomId",1);
    oldtree->SetBranchStatus("GENlep_MomMomId",1);
    oldtree->SetBranchStatus("massZ1",1);
    oldtree->SetBranchStatus("massZ2",1);
    oldtree->SetBranchStatus("pT4l",1);
    oldtree->SetBranchStatus("eta4l",1);
    oldtree->SetBranchStatus("rapidity4l",1);
    oldtree->SetBranchStatus("phi4l",1);
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
        oldtree->SetBranchStatus("dPhiHj1j2",1);
        oldtree->SetBranchStatus("dPhiHj1j2_jesup",1);
        oldtree->SetBranchStatus("dPhiHj1j2_jesdn",1);
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
        oldtree->SetBranchStatus("dPhiHj1j2_2p5",1);
        oldtree->SetBranchStatus("dPhiHj1j2_2p5_jesup",1);
        oldtree->SetBranchStatus("dPhiHj1j2_2p5_jesdn",1);
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
    oldtree->SetBranchStatus("lep_Hindex",1);
    oldtree->SetBranchStatus("GENlep_id",1);

    oldtree->SetBranchStatus("GENmass4l",1);
    oldtree->SetBranchStatus("GENmassZ1",1);
    oldtree->SetBranchStatus("GENmassZ2",1);
    oldtree->SetBranchStatus("GENmassZZ",1);
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

    newtree->Branch("isH4l",&isH4l,"isH4l/O");
// validation branches
    newtree->Branch("pTj1",&pTj1,"pTj1/F");
    newtree->Branch("TauC_Inc_0j_EnergyWgt_nom", &TauC_Inc_0j_EnergyWgt_nom, "TauC_Inc_0j_EnergyWgt_nom/F");
    newtree->Branch("TauB_Inc_0j_pTWgt_nom", &TauB_Inc_0j_pTWgt_nom, "TauB_Inc_0j_pTWgt_nom/F");

// Abs JES 
// jesup_Abs
    newtree->Branch("njets_pt30_eta4p7_jesup_Abs", &njets_pt30_eta4p7_jesup_Abs, "njets_pt30_eta4p7_jesup_Abs/F");
    newtree->Branch("TauC_Inc_0j_EnergyWgt_jesup_Abs", &TauC_Inc_0j_EnergyWgt_jesup_Abs, "TauC_Inc_0j_EnergyWgt_jesup_Abs/F");
    newtree->Branch("TauB_Inc_0j_pTWgt_jesup_Abs", &TauB_Inc_0j_pTWgt_jesup_Abs, "TauB_Inc_0j_pTWgt_jesup_Abs/F");
    newtree->Branch("njets_pt30_eta2p5_jesup_Abs", &njets_pt30_eta2p5_jesup_Abs, "njets_pt30_eta2p5_jesup_Abs/F");
    newtree->Branch("pt_leadingjet_pt30_eta4p7_jesup_Abs",&pt_leadingjet_pt30_eta4p7_jesup_Abs,"pt_leadingjet_pt30_eta4p7_jesup_Abs/F");
    newtree->Branch("pTj1_jesup_Abs",&pTj1_jesup_Abs,"pTj1_jesup_Abs/F");
    newtree->Branch("etaj1_jesup_Abs",&etaj1_jesup_Abs,"etaj1_jesup_Abs/F");
    newtree->Branch("pTj2_jesup_Abs",&pTj2_jesup_Abs,"pTj2_jesup_Abs/F");
    newtree->Branch("etaj2_jesup_Abs",&etaj2_jesup_Abs,"etaj2_jesup_Abs/F");
    newtree->Branch("yj1_jesup_Abs",&yj1_jesup_Abs,"yj1_jesup_Abs/F");
    newtree->Branch("yj2_jesup_Abs",&yj2_jesup_Abs,"yj2_jesup_Abs/F");
    newtree->Branch("dPhiHj1_jesup_Abs",&dPhiHj1_jesup_Abs,"dPhiHj1_jesup_Abs/F"); 
    newtree->Branch("mass4lj_jesup_Abs",&mass4lj_jesup_Abs,"mass4lj_jesup_Abs/F"); 
    newtree->Branch("mass4ljj_jesup_Abs",&mass4ljj_jesup_Abs,"mass4ljj_jesup_Abs/F"); 
    newtree->Branch("pT4lj_jesup_Abs",&pT4lj_jesup_Abs,"pT4lj_jesup_Abs/F"); 
    newtree->Branch("pT4ljj_jesup_Abs",&pT4ljj_jesup_Abs,"pT4ljj_jesup_Abs/F"); 
    newtree->Branch("dyHj1_jesup_Abs",&dyHj1_jesup_Abs,"dyHj1_jesup_Abs/F");
    newtree->Branch("mj1j2_jesup_Abs",&mj1j2_jesup_Abs,"mj1j2_jesup_Abs/F"); 
    newtree->Branch("dEtaj1j2_jesup_Abs",&dEtaj1j2_jesup_Abs,"dEtaj1j2_jesup_Abs/F");
    newtree->Branch("dPhij1j2_jesup_Abs",&dPhij1j2_jesup_Abs,"dPhij1j2_jesup_Abs/F"); 
    newtree->Branch("dPhiHj1j2_jesup_Abs",&dPhiHj1j2_jesup_Abs,"dPhiHj1j2_jesup_Abs/F");
    newtree->Branch("pTj1_2p5_jesup_Abs",&pTj1_2p5_jesup_Abs,"pTj1_2p5_jesup_Abs/F"); 
    newtree->Branch("pt_leadingjet_pt30_eta2p5_jesup_Abs",&pt_leadingjet_pt30_eta2p5_jesup_Abs,"pt_leadingjet_pt30_eta2p5_jesup_Abs/F");
    newtree->Branch("yj1_2p5_jesup_Abs",&yj1_2p5_jesup_Abs,"yj1_2p5_jesup_Abs/F");
    newtree->Branch("pTj2_2p5_jesup_Abs",&pTj2_2p5_jesup_Abs,"pTj2_2p5_jesup_Abs/F"); 
    newtree->Branch("yj2_2p5_jesup_Abs",&yj2_2p5_jesup_Abs,"yj2_2p5_jesup_Abs/F");
    newtree->Branch("dPhiHj1_2p5_jesup_Abs",&dPhiHj1_2p5_jesup_Abs,"dPhiHj1_2p5_jesup_Abs/F"); 
    newtree->Branch("mass4lj_2p5_jesup_Abs",&mass4lj_2p5_jesup_Abs,"mass4lj_2p5_jesup_Abs/F");
    newtree->Branch("mass4ljj_2p5_jesup_Abs",&mass4ljj_2p5_jesup_Abs,"mass4ljj_2p5_jesup_Abs/F");
    newtree->Branch("pT4lj_2p5_jesup_Abs",&pT4lj_2p5_jesup_Abs,"pT4lj_2p5_jesup_Abs/F");
    newtree->Branch("pT4ljj_2p5_jesup_Abs",&pT4ljj_2p5_jesup_Abs,"pT4ljj_2p5_jesup_Abs/F");
    newtree->Branch("dyHj1_2p5_jesup_Abs",&dyHj1_2p5_jesup_Abs,"dyHj1_2p5_jesup_Abs/F");
    newtree->Branch("mj1j2_2p5_jesup_Abs",&mj1j2_2p5_jesup_Abs,"mj1j2_2p5_jesup_Abs/F"); 
    newtree->Branch("dEtaj1j2_2p5_jesup_Abs",&dEtaj1j2_2p5_jesup_Abs,"dEtaj1j2_2p5_jesup_Abs/F");
    newtree->Branch("dPhij1j2_2p5_jesup_Abs",&dPhij1j2_2p5_jesup_Abs,"dPhij1j2_2p5_jesup_Abs/F"); 
    newtree->Branch("dPhiHj1j2_2p5_jesup_Abs",&dPhiHj1j2_2p5_jesup_Abs,"dPhiHj1j2_2p5_jesup_Abs/F");
// jesdn_Abs
    newtree->Branch("njets_pt30_eta4p7_jesdn_Abs", &njets_pt30_eta4p7_jesdn_Abs, "njets_pt30_eta4p7_jesdn_Abs/F");
    newtree->Branch("TauC_Inc_0j_EnergyWgt_jesdn_Abs", &TauC_Inc_0j_EnergyWgt_jesdn_Abs, "TauC_Inc_0j_EnergyWgt_jesdn_Abs/F");
    newtree->Branch("TauB_Inc_0j_pTWgt_jesdn_Abs", &TauB_Inc_0j_pTWgt_jesdn_Abs, "TauB_Inc_0j_pTWgt_jesdn_Abs/F");
    newtree->Branch("njets_pt30_eta2p5_jesdn_Abs", &njets_pt30_eta2p5_jesdn_Abs, "njets_pt30_eta2p5_jesdn_Abs/F");
    newtree->Branch("pt_leadingjet_pt30_eta4p7_jesdn_Abs",&pt_leadingjet_pt30_eta4p7_jesdn_Abs,"pt_leadingjet_pt30_eta4p7_jesdn_Abs/F");
    newtree->Branch("pTj1_jesdn_Abs",&pTj1_jesdn_Abs,"pTj1_jesdn_Abs/F");
    newtree->Branch("etaj1_jesdn_Abs",&etaj1_jesdn_Abs,"etaj1_jesdn_Abs/F");
    newtree->Branch("pTj2_jesdn_Abs",&pTj2_jesdn_Abs,"pTj2_jesdn_Abs/F");
    newtree->Branch("etaj2_jesdn_Abs",&etaj2_jesdn_Abs,"etaj2_jesdn_Abs/F");
    newtree->Branch("yj1_jesdn_Abs",&yj1_jesdn_Abs,"yj1_jesdn_Abs/F");
    newtree->Branch("yj2_jesdn_Abs",&yj2_jesdn_Abs,"yj2_jesdn_Abs/F");
    newtree->Branch("dPhiHj1_jesdn_Abs",&dPhiHj1_jesdn_Abs,"dPhiHj1_jesdn_Abs/F"); 
    newtree->Branch("dyHj1_jesdn_Abs",&dyHj1_jesdn_Abs,"dyHj1_jesdn_Abs/F");
    newtree->Branch("mass4lj_jesdn_Abs",&mass4lj_jesdn_Abs,"mass4lj_jesdn_Abs/F");
    newtree->Branch("mass4ljj_jesdn_Abs",&mass4ljj_jesdn_Abs,"mass4ljj_jesdn_Abs/F");
    newtree->Branch("pT4lj_jesdn_Abs",&pT4lj_jesdn_Abs,"pT4lj_jesdn_Abs/F");
    newtree->Branch("pT4ljj_jesdn_Abs",&pT4ljj_jesdn_Abs,"pT4ljj_jesdn_Abs/F");
    newtree->Branch("mj1j2_jesdn_Abs",&mj1j2_jesdn_Abs,"mj1j2_jesdn_Abs/F"); 
    newtree->Branch("dEtaj1j2_jesdn_Abs",&dEtaj1j2_jesdn_Abs,"dEtaj1j2_jesdn_Abs/F");
    newtree->Branch("dPhij1j2_jesdn_Abs",&dPhij1j2_jesdn_Abs,"dPhij1j2_jesdn_Abs/F"); 
    newtree->Branch("dPhiHj1j2_jesdn_Abs",&dPhiHj1j2_jesdn_Abs,"dPhiHj1j2_jesdn_Abs/F");
    newtree->Branch("pTj1_2p5_jesdn_Abs",&pTj1_2p5_jesdn_Abs,"pTj1_2p5_jesdn_Abs/F"); 
    newtree->Branch("pt_leadingjet_pt30_eta2p5_jesdn_Abs",&pt_leadingjet_pt30_eta2p5_jesdn_Abs,"pt_leadingjet_pt30_eta2p5_2p5_jesdn_Abs/F"); 
    newtree->Branch("yj1_2p5_jesdn_Abs",&yj1_2p5_jesdn_Abs,"yj1_2p5_jesdn_Abs/F");
    newtree->Branch("pTj2_2p5_jesdn_Abs",&pTj2_2p5_jesdn_Abs,"pTj2_2p5_jesdn_Abs/F"); 
    newtree->Branch("yj2_2p5_jesdn_Abs",&yj2_2p5_jesdn_Abs,"yj2_2p5_jesdn_Abs/F");
    newtree->Branch("mass4lj_2p5_jesdn_Abs",&mass4lj_2p5_jesdn_Abs,"mass4lj_2p5_jesdn_Abs/F");
    newtree->Branch("mass4ljj_2p5_jesdn_Abs",&mass4ljj_2p5_jesdn_Abs,"mass4ljj_2p5_jesdn_Abs/F");
    newtree->Branch("pT4lj_2p5_jesdn_Abs",&pT4lj_2p5_jesdn_Abs,"pT4lj_2p5_jesdn_Abs/F");
    newtree->Branch("pT4ljj_2p5_jesdn_Abs",&pT4ljj_2p5_jesdn_Abs,"pT4ljj_2p5_jesdn_Abs/F");
    newtree->Branch("dPhiHj1_2p5_jesdn_Abs",&dPhiHj1_2p5_jesdn_Abs,"dPhiHj1_2p5_jesdn_Abs/F"); 
    newtree->Branch("dyHj1_2p5_jesdn_Abs",&dyHj1_2p5_jesdn_Abs,"dyHj1_2p5_jesdn_Abs/F");
    newtree->Branch("mj1j2_2p5_jesdn_Abs",&mj1j2_2p5_jesdn_Abs,"mj1j2_2p5_jesdn_Abs/F"); 
    newtree->Branch("dEtaj1j2_2p5_jesdn_Abs",&dEtaj1j2_2p5_jesdn_Abs,"dEtaj1j2_2p5_jesdn_Abs/F");
    newtree->Branch("dPhij1j2_2p5_jesdn_Abs",&dPhij1j2_2p5_jesdn_Abs,"dPhij1j2_2p5_jesdn_Abs/F"); 
    newtree->Branch("dPhiHj1j2_2p5_jesdn_Abs",&dPhiHj1j2_2p5_jesdn_Abs,"dPhiHj1j2_2p5_jesdn_Abs/F");
//////////////////////////////////////////////////////////////////////////////////////////////////////

// Abs_year JES 
// jesup_Abs_year
    newtree->Branch("njets_pt30_eta4p7_jesup_Abs_year", &njets_pt30_eta4p7_jesup_Abs_year, "njets_pt30_eta4p7_jesup_Abs_year/F");
    newtree->Branch("TauC_Inc_0j_EnergyWgt_jesup_Abs_year", &TauC_Inc_0j_EnergyWgt_jesup_Abs_year, "TauC_Inc_0j_EnergyWgt_jesup_Abs_year/F");
    newtree->Branch("TauB_Inc_0j_pTWgt_jesup_Abs_year", &TauB_Inc_0j_pTWgt_jesup_Abs_year, "TauB_Inc_0j_pTWgt_jesup_Abs_year/F");
    newtree->Branch("njets_pt30_eta2p5_jesup_Abs_year", &njets_pt30_eta2p5_jesup_Abs_year, "njets_pt30_eta2p5_jesup_Abs_year/F");
    newtree->Branch("pt_leadingjet_pt30_eta4p7_jesup_Abs_year",&pt_leadingjet_pt30_eta4p7_jesup_Abs_year,"pt_leadingjet_pt30_eta4p7_jesup_Abs_year/F");
    newtree->Branch("pTj1_jesup_Abs_year",&pTj1_jesup_Abs_year,"pTj1_jesup_Abs_year/F");
    newtree->Branch("etaj1_jesup_Abs_year",&etaj1_jesup_Abs_year,"etaj1_jesup_Abs_year/F");
    newtree->Branch("pTj2_jesup_Abs_year",&pTj2_jesup_Abs_year,"pTj2_jesup_Abs_year/F");
    newtree->Branch("etaj2_jesup_Abs_year",&etaj2_jesup_Abs_year,"etaj2_jesup_Abs_year/F");
    newtree->Branch("yj1_jesup_Abs_year",&yj1_jesup_Abs_year,"yj1_jesup_Abs_year/F");
    newtree->Branch("yj2_jesup_Abs_year",&yj2_jesup_Abs_year,"yj2_jesup_Abs_year/F");
    newtree->Branch("dPhiHj1_jesup_Abs_year",&dPhiHj1_jesup_Abs_year,"dPhiHj1_jesup_Abs_year/F"); 
    newtree->Branch("mass4lj_jesup_Abs_year",&mass4lj_jesup_Abs_year,"mass4lj_jesup_Abs_year/F"); 
    newtree->Branch("mass4ljj_jesup_Abs_year",&mass4ljj_jesup_Abs_year,"mass4ljj_jesup_Abs_year/F"); 
    newtree->Branch("pT4lj_jesup_Abs_year",&pT4lj_jesup_Abs_year,"pT4lj_jesup_Abs_year/F"); 
    newtree->Branch("pT4ljj_jesup_Abs_year",&pT4ljj_jesup_Abs_year,"pT4ljj_jesup_Abs_year/F"); 
    newtree->Branch("dyHj1_jesup_Abs_year",&dyHj1_jesup_Abs_year,"dyHj1_jesup_Abs_year/F");
    newtree->Branch("mj1j2_jesup_Abs_year",&mj1j2_jesup_Abs_year,"mj1j2_jesup_Abs_year/F"); 
    newtree->Branch("dEtaj1j2_jesup_Abs_year",&dEtaj1j2_jesup_Abs_year,"dEtaj1j2_jesup_Abs_year/F");
    newtree->Branch("dPhij1j2_jesup_Abs_year",&dPhij1j2_jesup_Abs_year,"dPhij1j2_jesup_Abs_year/F"); 
    newtree->Branch("dPhiHj1j2_jesup_Abs_year",&dPhiHj1j2_jesup_Abs_year,"dPhiHj1j2_jesup_Abs_year/F");
    newtree->Branch("pTj1_2p5_jesup_Abs_year",&pTj1_2p5_jesup_Abs_year,"pTj1_2p5_jesup_Abs_year/F"); 
    newtree->Branch("pt_leadingjet_pt30_eta2p5_jesup_Abs_year",&pt_leadingjet_pt30_eta2p5_jesup_Abs_year,"pt_leadingjet_pt30_eta2p5_jesup_Abs_year/F");
    newtree->Branch("yj1_2p5_jesup_Abs_year",&yj1_2p5_jesup_Abs_year,"yj1_2p5_jesup_Abs_year/F");
    newtree->Branch("pTj2_2p5_jesup_Abs_year",&pTj2_2p5_jesup_Abs_year,"pTj2_2p5_jesup_Abs_year/F"); 
    newtree->Branch("yj2_2p5_jesup_Abs_year",&yj2_2p5_jesup_Abs_year,"yj2_2p5_jesup_Abs_year/F");
    newtree->Branch("dPhiHj1_2p5_jesup_Abs_year",&dPhiHj1_2p5_jesup_Abs_year,"dPhiHj1_2p5_jesup_Abs_year/F"); 
    newtree->Branch("mass4lj_2p5_jesup_Abs_year",&mass4lj_2p5_jesup_Abs_year,"mass4lj_2p5_jesup_Abs_year/F");
    newtree->Branch("mass4ljj_2p5_jesup_Abs_year",&mass4ljj_2p5_jesup_Abs_year,"mass4ljj_2p5_jesup_Abs_year/F");
    newtree->Branch("pT4lj_2p5_jesup_Abs_year",&pT4lj_2p5_jesup_Abs_year,"pT4lj_2p5_jesup_Abs_year/F");
    newtree->Branch("pT4ljj_2p5_jesup_Abs_year",&pT4ljj_2p5_jesup_Abs_year,"pT4ljj_2p5_jesup_Abs_year/F");
    newtree->Branch("dyHj1_2p5_jesup_Abs_year",&dyHj1_2p5_jesup_Abs_year,"dyHj1_2p5_jesup_Abs_year/F");
    newtree->Branch("mj1j2_2p5_jesup_Abs_year",&mj1j2_2p5_jesup_Abs_year,"mj1j2_2p5_jesup_Abs_year/F"); 
    newtree->Branch("dEtaj1j2_2p5_jesup_Abs_year",&dEtaj1j2_2p5_jesup_Abs_year,"dEtaj1j2_2p5_jesup_Abs_year/F");
    newtree->Branch("dPhij1j2_2p5_jesup_Abs_year",&dPhij1j2_2p5_jesup_Abs_year,"dPhij1j2_2p5_jesup_Abs_year/F"); 
    newtree->Branch("dPhiHj1j2_2p5_jesup_Abs_year",&dPhiHj1j2_2p5_jesup_Abs_year,"dPhiHj1j2_2p5_jesup_Abs_year/F");
// jesdn_Abs_year
    newtree->Branch("njets_pt30_eta4p7_jesdn_Abs_year", &njets_pt30_eta4p7_jesdn_Abs_year, "njets_pt30_eta4p7_jesdn_Abs_year/F");
    newtree->Branch("TauC_Inc_0j_EnergyWgt_jesdn_Abs_year", &TauC_Inc_0j_EnergyWgt_jesdn_Abs_year, "TauC_Inc_0j_EnergyWgt_jesdn_Abs_year/F");
    newtree->Branch("TauB_Inc_0j_pTWgt_jesdn_Abs_year", &TauB_Inc_0j_pTWgt_jesdn_Abs_year, "TauB_Inc_0j_pTWgt_jesdn_Abs_year/F");
    newtree->Branch("njets_pt30_eta2p5_jesdn_Abs_year", &njets_pt30_eta2p5_jesdn_Abs_year, "njets_pt30_eta2p5_jesdn_Abs_year/F");
    newtree->Branch("pt_leadingjet_pt30_eta4p7_jesdn_Abs_year",&pt_leadingjet_pt30_eta4p7_jesdn_Abs_year,"pt_leadingjet_pt30_eta4p7_jesdn_Abs_year/F");
    newtree->Branch("pTj1_jesdn_Abs_year",&pTj1_jesdn_Abs_year,"pTj1_jesdn_Abs_year/F");
    newtree->Branch("etaj1_jesdn_Abs_year",&etaj1_jesdn_Abs_year,"etaj1_jesdn_Abs_year/F");
    newtree->Branch("pTj2_jesdn_Abs_year",&pTj2_jesdn_Abs_year,"pTj2_jesdn_Abs_year/F");
    newtree->Branch("etaj2_jesdn_Abs_year",&etaj2_jesdn_Abs_year,"etaj2_jesdn_Abs_year/F");
    newtree->Branch("yj1_jesdn_Abs_year",&yj1_jesdn_Abs_year,"yj1_jesdn_Abs_year/F");
    newtree->Branch("yj2_jesdn_Abs_year",&yj2_jesdn_Abs_year,"yj2_jesdn_Abs_year/F");
    newtree->Branch("dPhiHj1_jesdn_Abs_year",&dPhiHj1_jesdn_Abs_year,"dPhiHj1_jesdn_Abs_year/F"); 
    newtree->Branch("dyHj1_jesdn_Abs_year",&dyHj1_jesdn_Abs_year,"dyHj1_jesdn_Abs_year/F");
    newtree->Branch("mass4lj_jesdn_Abs_year",&mass4lj_jesdn_Abs_year,"mass4lj_jesdn_Abs_year/F");
    newtree->Branch("mass4ljj_jesdn_Abs_year",&mass4ljj_jesdn_Abs_year,"mass4ljj_jesdn_Abs_year/F");
    newtree->Branch("pT4lj_jesdn_Abs_year",&pT4lj_jesdn_Abs_year,"pT4lj_jesdn_Abs_year/F");
    newtree->Branch("pT4ljj_jesdn_Abs_year",&pT4ljj_jesdn_Abs_year,"pT4ljj_jesdn_Abs_year/F");
    newtree->Branch("mj1j2_jesdn_Abs_year",&mj1j2_jesdn_Abs_year,"mj1j2_jesdn_Abs_year/F"); 
    newtree->Branch("dEtaj1j2_jesdn_Abs_year",&dEtaj1j2_jesdn_Abs_year,"dEtaj1j2_jesdn_Abs_year/F");
    newtree->Branch("dPhij1j2_jesdn_Abs_year",&dPhij1j2_jesdn_Abs_year,"dPhij1j2_jesdn_Abs_year/F"); 
    newtree->Branch("dPhiHj1j2_jesdn_Abs_year",&dPhiHj1j2_jesdn_Abs_year,"dPhiHj1j2_jesdn_Abs_year/F");
    newtree->Branch("pTj1_2p5_jesdn_Abs_year",&pTj1_2p5_jesdn_Abs_year,"pTj1_2p5_jesdn_Abs_year/F"); 
    newtree->Branch("pt_leadingjet_pt30_eta2p5_jesdn_Abs_year",&pt_leadingjet_pt30_eta2p5_jesdn_Abs_year,"pt_leadingjet_pt30_eta2p5_2p5_jesdn_Abs_year/F"); 
    newtree->Branch("yj1_2p5_jesdn_Abs_year",&yj1_2p5_jesdn_Abs_year,"yj1_2p5_jesdn_Abs_year/F");
    newtree->Branch("pTj2_2p5_jesdn_Abs_year",&pTj2_2p5_jesdn_Abs_year,"pTj2_2p5_jesdn_Abs_year/F"); 
    newtree->Branch("yj2_2p5_jesdn_Abs_year",&yj2_2p5_jesdn_Abs_year,"yj2_2p5_jesdn_Abs_year/F");
    newtree->Branch("mass4lj_2p5_jesdn_Abs_year",&mass4lj_2p5_jesdn_Abs_year,"mass4lj_2p5_jesdn_Abs_year/F");
    newtree->Branch("mass4ljj_2p5_jesdn_Abs_year",&mass4ljj_2p5_jesdn_Abs_year,"mass4ljj_2p5_jesdn_Abs_year/F");
    newtree->Branch("pT4lj_2p5_jesdn_Abs_year",&pT4lj_2p5_jesdn_Abs_year,"pT4lj_2p5_jesdn_Abs_year/F");
    newtree->Branch("pT4ljj_2p5_jesdn_Abs_year",&pT4ljj_2p5_jesdn_Abs_year,"pT4ljj_2p5_jesdn_Abs_year/F");
    newtree->Branch("dPhiHj1_2p5_jesdn_Abs_year",&dPhiHj1_2p5_jesdn_Abs_year,"dPhiHj1_2p5_jesdn_Abs_year/F"); 
    newtree->Branch("dyHj1_2p5_jesdn_Abs_year",&dyHj1_2p5_jesdn_Abs_year,"dyHj1_2p5_jesdn_Abs_year/F");
    newtree->Branch("mj1j2_2p5_jesdn_Abs_year",&mj1j2_2p5_jesdn_Abs_year,"mj1j2_2p5_jesdn_Abs_year/F"); 
    newtree->Branch("dEtaj1j2_2p5_jesdn_Abs_year",&dEtaj1j2_2p5_jesdn_Abs_year,"dEtaj1j2_2p5_jesdn_Abs_year/F");
    newtree->Branch("dPhij1j2_2p5_jesdn_Abs_year",&dPhij1j2_2p5_jesdn_Abs_year,"dPhij1j2_2p5_jesdn_Abs_year/F"); 
    newtree->Branch("dPhiHj1j2_2p5_jesdn_Abs_year",&dPhiHj1j2_2p5_jesdn_Abs_year,"dPhiHj1j2_2p5_jesdn_Abs_year/F");
/////////////////////////////////////////////////////////////////////////////////////////////////////

// BBEC1 JES 
// jesup_BBEC1
    newtree->Branch("njets_pt30_eta4p7_jesup_BBEC1", &njets_pt30_eta4p7_jesup_BBEC1, "njets_pt30_eta4p7_jesup_BBEC1/F");
    newtree->Branch("TauC_Inc_0j_EnergyWgt_jesup_BBEC1", &TauC_Inc_0j_EnergyWgt_jesup_BBEC1, "TauC_Inc_0j_EnergyWgt_jesup_BBEC1/F");
    newtree->Branch("TauB_Inc_0j_pTWgt_jesup_BBEC1", &TauB_Inc_0j_pTWgt_jesup_BBEC1, "TauB_Inc_0j_pTWgt_jesup_BBEC1/F");
    newtree->Branch("njets_pt30_eta2p5_jesup_BBEC1", &njets_pt30_eta2p5_jesup_BBEC1, "njets_pt30_eta2p5_jesup_BBEC1/F");
    newtree->Branch("pt_leadingjet_pt30_eta4p7_jesup_BBEC1",&pt_leadingjet_pt30_eta4p7_jesup_BBEC1,"pt_leadingjet_pt30_eta4p7_jesup_BBEC1/F");
    newtree->Branch("pTj1_jesup_BBEC1",&pTj1_jesup_BBEC1,"pTj1_jesup_BBEC1/F");
    newtree->Branch("etaj1_jesup_BBEC1",&etaj1_jesup_BBEC1,"etaj1_jesup_BBEC1/F");
    newtree->Branch("pTj2_jesup_BBEC1",&pTj2_jesup_BBEC1,"pTj2_jesup_BBEC1/F");
    newtree->Branch("etaj2_jesup_BBEC1",&etaj2_jesup_BBEC1,"etaj2_jesup_BBEC1/F");
    newtree->Branch("yj1_jesup_BBEC1",&yj1_jesup_BBEC1,"yj1_jesup_BBEC1/F");
    newtree->Branch("yj2_jesup_BBEC1",&yj2_jesup_BBEC1,"yj2_jesup_BBEC1/F");
    newtree->Branch("dPhiHj1_jesup_BBEC1",&dPhiHj1_jesup_BBEC1,"dPhiHj1_jesup_BBEC1/F"); 
    newtree->Branch("mass4lj_jesup_BBEC1",&mass4lj_jesup_BBEC1,"mass4lj_jesup_BBEC1/F"); 
    newtree->Branch("mass4ljj_jesup_BBEC1",&mass4ljj_jesup_BBEC1,"mass4ljj_jesup_BBEC1/F"); 
    newtree->Branch("pT4lj_jesup_BBEC1",&pT4lj_jesup_BBEC1,"pT4lj_jesup_BBEC1/F"); 
    newtree->Branch("pT4ljj_jesup_BBEC1",&pT4ljj_jesup_BBEC1,"pT4ljj_jesup_BBEC1/F"); 
    newtree->Branch("dyHj1_jesup_BBEC1",&dyHj1_jesup_BBEC1,"dyHj1_jesup_BBEC1/F");
    newtree->Branch("mj1j2_jesup_BBEC1",&mj1j2_jesup_BBEC1,"mj1j2_jesup_BBEC1/F"); 
    newtree->Branch("dEtaj1j2_jesup_BBEC1",&dEtaj1j2_jesup_BBEC1,"dEtaj1j2_jesup_BBEC1/F");
    newtree->Branch("dPhij1j2_jesup_BBEC1",&dPhij1j2_jesup_BBEC1,"dPhij1j2_jesup_BBEC1/F"); 
    newtree->Branch("dPhiHj1j2_jesup_BBEC1",&dPhiHj1j2_jesup_BBEC1,"dPhiHj1j2_jesup_BBEC1/F");
    newtree->Branch("pTj1_2p5_jesup_BBEC1",&pTj1_2p5_jesup_BBEC1,"pTj1_2p5_jesup_BBEC1/F"); 
    newtree->Branch("pt_leadingjet_pt30_eta2p5_jesup_BBEC1",&pt_leadingjet_pt30_eta2p5_jesup_BBEC1,"pt_leadingjet_pt30_eta2p5_jesup_BBEC1/F");
    newtree->Branch("yj1_2p5_jesup_BBEC1",&yj1_2p5_jesup_BBEC1,"yj1_2p5_jesup_BBEC1/F");
    newtree->Branch("pTj2_2p5_jesup_BBEC1",&pTj2_2p5_jesup_BBEC1,"pTj2_2p5_jesup_BBEC1/F"); 
    newtree->Branch("yj2_2p5_jesup_BBEC1",&yj2_2p5_jesup_BBEC1,"yj2_2p5_jesup_BBEC1/F");
    newtree->Branch("dPhiHj1_2p5_jesup_BBEC1",&dPhiHj1_2p5_jesup_BBEC1,"dPhiHj1_2p5_jesup_BBEC1/F"); 
    newtree->Branch("mass4lj_2p5_jesup_BBEC1",&mass4lj_2p5_jesup_BBEC1,"mass4lj_2p5_jesup_BBEC1/F");
    newtree->Branch("mass4ljj_2p5_jesup_BBEC1",&mass4ljj_2p5_jesup_BBEC1,"mass4ljj_2p5_jesup_BBEC1/F");
    newtree->Branch("pT4lj_2p5_jesup_BBEC1",&pT4lj_2p5_jesup_BBEC1,"pT4lj_2p5_jesup_BBEC1/F");
    newtree->Branch("pT4ljj_2p5_jesup_BBEC1",&pT4ljj_2p5_jesup_BBEC1,"pT4ljj_2p5_jesup_BBEC1/F");
    newtree->Branch("dyHj1_2p5_jesup_BBEC1",&dyHj1_2p5_jesup_BBEC1,"dyHj1_2p5_jesup_BBEC1/F");
    newtree->Branch("mj1j2_2p5_jesup_BBEC1",&mj1j2_2p5_jesup_BBEC1,"mj1j2_2p5_jesup_BBEC1/F"); 
    newtree->Branch("dEtaj1j2_2p5_jesup_BBEC1",&dEtaj1j2_2p5_jesup_BBEC1,"dEtaj1j2_2p5_jesup_BBEC1/F");
    newtree->Branch("dPhij1j2_2p5_jesup_BBEC1",&dPhij1j2_2p5_jesup_BBEC1,"dPhij1j2_2p5_jesup_BBEC1/F"); 
    newtree->Branch("dPhiHj1j2_2p5_jesup_BBEC1",&dPhiHj1j2_2p5_jesup_BBEC1,"dPhiHj1j2_2p5_jesup_BBEC1/F");
// jesdn_BBEC1
    newtree->Branch("njets_pt30_eta4p7_jesdn_BBEC1", &njets_pt30_eta4p7_jesdn_BBEC1, "njets_pt30_eta4p7_jesdn_BBEC1/F");
    newtree->Branch("TauC_Inc_0j_EnergyWgt_jesdn_BBEC1", &TauC_Inc_0j_EnergyWgt_jesdn_BBEC1, "TauC_Inc_0j_EnergyWgt_jesdn_BBEC1/F");
    newtree->Branch("TauB_Inc_0j_pTWgt_jesdn_BBEC1", &TauB_Inc_0j_pTWgt_jesdn_BBEC1, "TauB_Inc_0j_pTWgt_jesdn_BBEC1/F");
    newtree->Branch("njets_pt30_eta2p5_jesdn_BBEC1", &njets_pt30_eta2p5_jesdn_BBEC1, "njets_pt30_eta2p5_jesdn_BBEC1/F");
    newtree->Branch("pt_leadingjet_pt30_eta4p7_jesdn_BBEC1",&pt_leadingjet_pt30_eta4p7_jesdn_BBEC1,"pt_leadingjet_pt30_eta4p7_jesdn_BBEC1/F");
    newtree->Branch("pTj1_jesdn_BBEC1",&pTj1_jesdn_BBEC1,"pTj1_jesdn_BBEC1/F");
    newtree->Branch("etaj1_jesdn_BBEC1",&etaj1_jesdn_BBEC1,"etaj1_jesdn_BBEC1/F");
    newtree->Branch("pTj2_jesdn_BBEC1",&pTj2_jesdn_BBEC1,"pTj2_jesdn_BBEC1/F");
    newtree->Branch("etaj2_jesdn_BBEC1",&etaj2_jesdn_BBEC1,"etaj2_jesdn_BBEC1/F");
    newtree->Branch("yj1_jesdn_BBEC1",&yj1_jesdn_BBEC1,"yj1_jesdn_BBEC1/F");
    newtree->Branch("yj2_jesdn_BBEC1",&yj2_jesdn_BBEC1,"yj2_jesdn_BBEC1/F");
    newtree->Branch("dPhiHj1_jesdn_BBEC1",&dPhiHj1_jesdn_BBEC1,"dPhiHj1_jesdn_BBEC1/F"); 
    newtree->Branch("dyHj1_jesdn_BBEC1",&dyHj1_jesdn_BBEC1,"dyHj1_jesdn_BBEC1/F");
    newtree->Branch("mass4lj_jesdn_BBEC1",&mass4lj_jesdn_BBEC1,"mass4lj_jesdn_BBEC1/F");
    newtree->Branch("mass4ljj_jesdn_BBEC1",&mass4ljj_jesdn_BBEC1,"mass4ljj_jesdn_BBEC1/F");
    newtree->Branch("pT4lj_jesdn_BBEC1",&pT4lj_jesdn_BBEC1,"pT4lj_jesdn_BBEC1/F");
    newtree->Branch("pT4ljj_jesdn_BBEC1",&pT4ljj_jesdn_BBEC1,"pT4ljj_jesdn_BBEC1/F");
    newtree->Branch("mj1j2_jesdn_BBEC1",&mj1j2_jesdn_BBEC1,"mj1j2_jesdn_BBEC1/F"); 
    newtree->Branch("dEtaj1j2_jesdn_BBEC1",&dEtaj1j2_jesdn_BBEC1,"dEtaj1j2_jesdn_BBEC1/F");
    newtree->Branch("dPhij1j2_jesdn_BBEC1",&dPhij1j2_jesdn_BBEC1,"dPhij1j2_jesdn_BBEC1/F"); 
    newtree->Branch("dPhiHj1j2_jesdn_BBEC1",&dPhiHj1j2_jesdn_BBEC1,"dPhiHj1j2_jesdn_BBEC1/F");
    newtree->Branch("pTj1_2p5_jesdn_BBEC1",&pTj1_2p5_jesdn_BBEC1,"pTj1_2p5_jesdn_BBEC1/F"); 
    newtree->Branch("pt_leadingjet_pt30_eta2p5_jesdn_BBEC1",&pt_leadingjet_pt30_eta2p5_jesdn_BBEC1,"pt_leadingjet_pt30_eta2p5_2p5_jesdn_BBEC1/F"); 
    newtree->Branch("yj1_2p5_jesdn_BBEC1",&yj1_2p5_jesdn_BBEC1,"yj1_2p5_jesdn_BBEC1/F");
    newtree->Branch("pTj2_2p5_jesdn_BBEC1",&pTj2_2p5_jesdn_BBEC1,"pTj2_2p5_jesdn_BBEC1/F"); 
    newtree->Branch("yj2_2p5_jesdn_BBEC1",&yj2_2p5_jesdn_BBEC1,"yj2_2p5_jesdn_BBEC1/F");
    newtree->Branch("mass4lj_2p5_jesdn_BBEC1",&mass4lj_2p5_jesdn_BBEC1,"mass4lj_2p5_jesdn_BBEC1/F");
    newtree->Branch("mass4ljj_2p5_jesdn_BBEC1",&mass4ljj_2p5_jesdn_BBEC1,"mass4ljj_2p5_jesdn_BBEC1/F");
    newtree->Branch("pT4lj_2p5_jesdn_BBEC1",&pT4lj_2p5_jesdn_BBEC1,"pT4lj_2p5_jesdn_BBEC1/F");
    newtree->Branch("pT4ljj_2p5_jesdn_BBEC1",&pT4ljj_2p5_jesdn_BBEC1,"pT4ljj_2p5_jesdn_BBEC1/F");
    newtree->Branch("dPhiHj1_2p5_jesdn_BBEC1",&dPhiHj1_2p5_jesdn_BBEC1,"dPhiHj1_2p5_jesdn_BBEC1/F"); 
    newtree->Branch("dyHj1_2p5_jesdn_BBEC1",&dyHj1_2p5_jesdn_BBEC1,"dyHj1_2p5_jesdn_BBEC1/F");
    newtree->Branch("mj1j2_2p5_jesdn_BBEC1",&mj1j2_2p5_jesdn_BBEC1,"mj1j2_2p5_jesdn_BBEC1/F"); 
    newtree->Branch("dEtaj1j2_2p5_jesdn_BBEC1",&dEtaj1j2_2p5_jesdn_BBEC1,"dEtaj1j2_2p5_jesdn_BBEC1/F");
    newtree->Branch("dPhij1j2_2p5_jesdn_BBEC1",&dPhij1j2_2p5_jesdn_BBEC1,"dPhij1j2_2p5_jesdn_BBEC1/F"); 
    newtree->Branch("dPhiHj1j2_2p5_jesdn_BBEC1",&dPhiHj1j2_2p5_jesdn_BBEC1,"dPhiHj1j2_2p5_jesdn_BBEC1/F");
//////////////////////////////////////////////////////////////////////////////////////////////////////////

// BBEC1_year JES 
// jesup_BBEC1_year
    newtree->Branch("njets_pt30_eta4p7_jesup_BBEC1_year", &njets_pt30_eta4p7_jesup_BBEC1_year, "njets_pt30_eta4p7_jesup_BBEC1_year/F");
    newtree->Branch("TauC_Inc_0j_EnergyWgt_jesup_BBEC1_year", &TauC_Inc_0j_EnergyWgt_jesup_BBEC1_year, "TauC_Inc_0j_EnergyWgt_jesup_BBEC1_year/F");
    newtree->Branch("TauB_Inc_0j_pTWgt_jesup_BBEC1_year", &TauB_Inc_0j_pTWgt_jesup_BBEC1_year, "TauB_Inc_0j_pTWgt_jesup_BBEC1_year/F");
    newtree->Branch("njets_pt30_eta2p5_jesup_BBEC1_year", &njets_pt30_eta2p5_jesup_BBEC1_year, "njets_pt30_eta2p5_jesup_BBEC1_year/F");
    newtree->Branch("pt_leadingjet_pt30_eta4p7_jesup_BBEC1_year",&pt_leadingjet_pt30_eta4p7_jesup_BBEC1_year,"pt_leadingjet_pt30_eta4p7_jesup_BBEC1_year/F");
    newtree->Branch("pTj1_jesup_BBEC1_year",&pTj1_jesup_BBEC1_year,"pTj1_jesup_BBEC1_year/F");
    newtree->Branch("etaj1_jesup_BBEC1_year",&etaj1_jesup_BBEC1_year,"etaj1_jesup_BBEC1_year/F");
    newtree->Branch("pTj2_jesup_BBEC1_year",&pTj2_jesup_BBEC1_year,"pTj2_jesup_BBEC1_year/F");
    newtree->Branch("etaj2_jesup_BBEC1_year",&etaj2_jesup_BBEC1_year,"etaj2_jesup_BBEC1_year/F");
    newtree->Branch("yj1_jesup_BBEC1_year",&yj1_jesup_BBEC1_year,"yj1_jesup_BBEC1_year/F");
    newtree->Branch("yj2_jesup_BBEC1_year",&yj2_jesup_BBEC1_year,"yj2_jesup_BBEC1_year/F");
    newtree->Branch("dPhiHj1_jesup_BBEC1_year",&dPhiHj1_jesup_BBEC1_year,"dPhiHj1_jesup_BBEC1_year/F"); 
    newtree->Branch("mass4lj_jesup_BBEC1_year",&mass4lj_jesup_BBEC1_year,"mass4lj_jesup_BBEC1_year/F"); 
    newtree->Branch("mass4ljj_jesup_BBEC1_year",&mass4ljj_jesup_BBEC1_year,"mass4ljj_jesup_BBEC1_year/F"); 
    newtree->Branch("pT4lj_jesup_BBEC1_year",&pT4lj_jesup_BBEC1_year,"pT4lj_jesup_BBEC1_year/F"); 
    newtree->Branch("pT4ljj_jesup_BBEC1_year",&pT4ljj_jesup_BBEC1_year,"pT4ljj_jesup_BBEC1_year/F"); 
    newtree->Branch("dyHj1_jesup_BBEC1_year",&dyHj1_jesup_BBEC1_year,"dyHj1_jesup_BBEC1_year/F");
    newtree->Branch("mj1j2_jesup_BBEC1_year",&mj1j2_jesup_BBEC1_year,"mj1j2_jesup_BBEC1_year/F"); 
    newtree->Branch("dEtaj1j2_jesup_BBEC1_year",&dEtaj1j2_jesup_BBEC1_year,"dEtaj1j2_jesup_BBEC1_year/F");
    newtree->Branch("dPhij1j2_jesup_BBEC1_year",&dPhij1j2_jesup_BBEC1_year,"dPhij1j2_jesup_BBEC1_year/F"); 
    newtree->Branch("dPhiHj1j2_jesup_BBEC1_year",&dPhiHj1j2_jesup_BBEC1_year,"dPhiHj1j2_jesup_BBEC1_year/F");
    newtree->Branch("pTj1_2p5_jesup_BBEC1_year",&pTj1_2p5_jesup_BBEC1_year,"pTj1_2p5_jesup_BBEC1_year/F"); 
    newtree->Branch("pt_leadingjet_pt30_eta2p5_jesup_BBEC1_year",&pt_leadingjet_pt30_eta2p5_jesup_BBEC1_year,"pt_leadingjet_pt30_eta2p5_jesup_BBEC1_year/F");
    newtree->Branch("yj1_2p5_jesup_BBEC1_year",&yj1_2p5_jesup_BBEC1_year,"yj1_2p5_jesup_BBEC1_year/F");
    newtree->Branch("pTj2_2p5_jesup_BBEC1_year",&pTj2_2p5_jesup_BBEC1_year,"pTj2_2p5_jesup_BBEC1_year/F"); 
    newtree->Branch("yj2_2p5_jesup_BBEC1_year",&yj2_2p5_jesup_BBEC1_year,"yj2_2p5_jesup_BBEC1_year/F");
    newtree->Branch("dPhiHj1_2p5_jesup_BBEC1_year",&dPhiHj1_2p5_jesup_BBEC1_year,"dPhiHj1_2p5_jesup_BBEC1_year/F"); 
    newtree->Branch("mass4lj_2p5_jesup_BBEC1_year",&mass4lj_2p5_jesup_BBEC1_year,"mass4lj_2p5_jesup_BBEC1_year/F");
    newtree->Branch("mass4ljj_2p5_jesup_BBEC1_year",&mass4ljj_2p5_jesup_BBEC1_year,"mass4ljj_2p5_jesup_BBEC1_year/F");
    newtree->Branch("pT4lj_2p5_jesup_BBEC1_year",&pT4lj_2p5_jesup_BBEC1_year,"pT4lj_2p5_jesup_BBEC1_year/F");
    newtree->Branch("pT4ljj_2p5_jesup_BBEC1_year",&pT4ljj_2p5_jesup_BBEC1_year,"pT4ljj_2p5_jesup_BBEC1_year/F");
    newtree->Branch("dyHj1_2p5_jesup_BBEC1_year",&dyHj1_2p5_jesup_BBEC1_year,"dyHj1_2p5_jesup_BBEC1_year/F");
    newtree->Branch("mj1j2_2p5_jesup_BBEC1_year",&mj1j2_2p5_jesup_BBEC1_year,"mj1j2_2p5_jesup_BBEC1_year/F"); 
    newtree->Branch("dEtaj1j2_2p5_jesup_BBEC1_year",&dEtaj1j2_2p5_jesup_BBEC1_year,"dEtaj1j2_2p5_jesup_BBEC1_year/F");
    newtree->Branch("dPhij1j2_2p5_jesup_BBEC1_year",&dPhij1j2_2p5_jesup_BBEC1_year,"dPhij1j2_2p5_jesup_BBEC1_year/F"); 
    newtree->Branch("dPhiHj1j2_2p5_jesup_BBEC1_year",&dPhiHj1j2_2p5_jesup_BBEC1_year,"dPhiHj1j2_2p5_jesup_BBEC1_year/F");
// jesdn_BBEC1_year
    newtree->Branch("njets_pt30_eta4p7_jesdn_BBEC1_year", &njets_pt30_eta4p7_jesdn_BBEC1_year, "njets_pt30_eta4p7_jesdn_BBEC1_year/F");
    newtree->Branch("TauC_Inc_0j_EnergyWgt_jesdn_BBEC1_year", &TauC_Inc_0j_EnergyWgt_jesdn_BBEC1_year, "TauC_Inc_0j_EnergyWgt_jesdn_BBEC1_year/F");
    newtree->Branch("TauB_Inc_0j_pTWgt_jesdn_BBEC1_year", &TauB_Inc_0j_pTWgt_jesdn_BBEC1_year, "TauB_Inc_0j_pTWgt_jesdn_BBEC1_year/F");
    newtree->Branch("njets_pt30_eta2p5_jesdn_BBEC1_year", &njets_pt30_eta2p5_jesdn_BBEC1_year, "njets_pt30_eta2p5_jesdn_BBEC1_year/F");
    newtree->Branch("pt_leadingjet_pt30_eta4p7_jesdn_BBEC1_year",&pt_leadingjet_pt30_eta4p7_jesdn_BBEC1_year,"pt_leadingjet_pt30_eta4p7_jesdn_BBEC1_year/F");
    newtree->Branch("pTj1_jesdn_BBEC1_year",&pTj1_jesdn_BBEC1_year,"pTj1_jesdn_BBEC1_year/F");
    newtree->Branch("etaj1_jesdn_BBEC1_year",&etaj1_jesdn_BBEC1_year,"etaj1_jesdn_BBEC1_year/F");
    newtree->Branch("pTj2_jesdn_BBEC1_year",&pTj2_jesdn_BBEC1_year,"pTj2_jesdn_BBEC1_year/F");
    newtree->Branch("etaj2_jesdn_BBEC1_year",&etaj2_jesdn_BBEC1_year,"etaj2_jesdn_BBEC1_year/F");
    newtree->Branch("yj1_jesdn_BBEC1_year",&yj1_jesdn_BBEC1_year,"yj1_jesdn_BBEC1_year/F");
    newtree->Branch("yj2_jesdn_BBEC1_year",&yj2_jesdn_BBEC1_year,"yj2_jesdn_BBEC1_year/F");
    newtree->Branch("dPhiHj1_jesdn_BBEC1_year",&dPhiHj1_jesdn_BBEC1_year,"dPhiHj1_jesdn_BBEC1_year/F"); 
    newtree->Branch("dyHj1_jesdn_BBEC1_year",&dyHj1_jesdn_BBEC1_year,"dyHj1_jesdn_BBEC1_year/F");
    newtree->Branch("mass4lj_jesdn_BBEC1_year",&mass4lj_jesdn_BBEC1_year,"mass4lj_jesdn_BBEC1_year/F");
    newtree->Branch("mass4ljj_jesdn_BBEC1_year",&mass4ljj_jesdn_BBEC1_year,"mass4ljj_jesdn_BBEC1_year/F");
    newtree->Branch("pT4lj_jesdn_BBEC1_year",&pT4lj_jesdn_BBEC1_year,"pT4lj_jesdn_BBEC1_year/F");
    newtree->Branch("pT4ljj_jesdn_BBEC1_year",&pT4ljj_jesdn_BBEC1_year,"pT4ljj_jesdn_BBEC1_year/F");
    newtree->Branch("mj1j2_jesdn_BBEC1_year",&mj1j2_jesdn_BBEC1_year,"mj1j2_jesdn_BBEC1_year/F"); 
    newtree->Branch("dEtaj1j2_jesdn_BBEC1_year",&dEtaj1j2_jesdn_BBEC1_year,"dEtaj1j2_jesdn_BBEC1_year/F");
    newtree->Branch("dPhij1j2_jesdn_BBEC1_year",&dPhij1j2_jesdn_BBEC1_year,"dPhij1j2_jesdn_BBEC1_year/F"); 
    newtree->Branch("dPhiHj1j2_jesdn_BBEC1_year",&dPhiHj1j2_jesdn_BBEC1_year,"dPhiHj1j2_jesdn_BBEC1_year/F");
    newtree->Branch("pTj1_2p5_jesdn_BBEC1_year",&pTj1_2p5_jesdn_BBEC1_year,"pTj1_2p5_jesdn_BBEC1_year/F"); 
    newtree->Branch("pt_leadingjet_pt30_eta2p5_jesdn_BBEC1_year",&pt_leadingjet_pt30_eta2p5_jesdn_BBEC1_year,"pt_leadingjet_pt30_eta2p5_2p5_jesdn_BBEC1_year/F"); 
    newtree->Branch("yj1_2p5_jesdn_BBEC1_year",&yj1_2p5_jesdn_BBEC1_year,"yj1_2p5_jesdn_BBEC1_year/F");
    newtree->Branch("pTj2_2p5_jesdn_BBEC1_year",&pTj2_2p5_jesdn_BBEC1_year,"pTj2_2p5_jesdn_BBEC1_year/F"); 
    newtree->Branch("yj2_2p5_jesdn_BBEC1_year",&yj2_2p5_jesdn_BBEC1_year,"yj2_2p5_jesdn_BBEC1_year/F");
    newtree->Branch("mass4lj_2p5_jesdn_BBEC1_year",&mass4lj_2p5_jesdn_BBEC1_year,"mass4lj_2p5_jesdn_BBEC1_year/F");
    newtree->Branch("mass4ljj_2p5_jesdn_BBEC1_year",&mass4ljj_2p5_jesdn_BBEC1_year,"mass4ljj_2p5_jesdn_BBEC1_year/F");
    newtree->Branch("pT4lj_2p5_jesdn_BBEC1_year",&pT4lj_2p5_jesdn_BBEC1_year,"pT4lj_2p5_jesdn_BBEC1_year/F");
    newtree->Branch("pT4ljj_2p5_jesdn_BBEC1_year",&pT4ljj_2p5_jesdn_BBEC1_year,"pT4ljj_2p5_jesdn_BBEC1_year/F");
    newtree->Branch("dPhiHj1_2p5_jesdn_BBEC1_year",&dPhiHj1_2p5_jesdn_BBEC1_year,"dPhiHj1_2p5_jesdn_BBEC1_year/F"); 
    newtree->Branch("dyHj1_2p5_jesdn_BBEC1_year",&dyHj1_2p5_jesdn_BBEC1_year,"dyHj1_2p5_jesdn_BBEC1_year/F");
    newtree->Branch("mj1j2_2p5_jesdn_BBEC1_year",&mj1j2_2p5_jesdn_BBEC1_year,"mj1j2_2p5_jesdn_BBEC1_year/F"); 
    newtree->Branch("dEtaj1j2_2p5_jesdn_BBEC1_year",&dEtaj1j2_2p5_jesdn_BBEC1_year,"dEtaj1j2_2p5_jesdn_BBEC1_year/F");
    newtree->Branch("dPhij1j2_2p5_jesdn_BBEC1_year",&dPhij1j2_2p5_jesdn_BBEC1_year,"dPhij1j2_2p5_jesdn_BBEC1_year/F"); 
    newtree->Branch("dPhiHj1j2_2p5_jesdn_BBEC1_year",&dPhiHj1j2_2p5_jesdn_BBEC1_year,"dPhiHj1j2_2p5_jesdn_BBEC1_year/F");
////////////////////////////////////////////////////////////////////////////////////////////////////////////

// EC2 JES 
// jesup_EC2
    newtree->Branch("njets_pt30_eta4p7_jesup_EC2", &njets_pt30_eta4p7_jesup_EC2, "njets_pt30_eta4p7_jesup_EC2/F");
    newtree->Branch("TauC_Inc_0j_EnergyWgt_jesup_EC2", &TauC_Inc_0j_EnergyWgt_jesup_EC2, "TauC_Inc_0j_EnergyWgt_jesup_EC2/F");
    newtree->Branch("TauB_Inc_0j_pTWgt_jesup_EC2", &TauB_Inc_0j_pTWgt_jesup_EC2, "TauB_Inc_0j_pTWgt_jesup_EC2/F");
    newtree->Branch("njets_pt30_eta2p5_jesup_EC2", &njets_pt30_eta2p5_jesup_EC2, "njets_pt30_eta2p5_jesup_EC2/F");
    newtree->Branch("pt_leadingjet_pt30_eta4p7_jesup_EC2",&pt_leadingjet_pt30_eta4p7_jesup_EC2,"pt_leadingjet_pt30_eta4p7_jesup_EC2/F");
    newtree->Branch("pTj1_jesup_EC2",&pTj1_jesup_EC2,"pTj1_jesup_EC2/F");
    newtree->Branch("etaj1_jesup_EC2",&etaj1_jesup_EC2,"etaj1_jesup_EC2/F");
    newtree->Branch("pTj2_jesup_EC2",&pTj2_jesup_EC2,"pTj2_jesup_EC2/F");
    newtree->Branch("etaj2_jesup_EC2",&etaj2_jesup_EC2,"etaj2_jesup_EC2/F");
    newtree->Branch("yj1_jesup_EC2",&yj1_jesup_EC2,"yj1_jesup_EC2/F");
    newtree->Branch("yj2_jesup_EC2",&yj2_jesup_EC2,"yj2_jesup_EC2/F");
    newtree->Branch("dPhiHj1_jesup_EC2",&dPhiHj1_jesup_EC2,"dPhiHj1_jesup_EC2/F"); 
    newtree->Branch("mass4lj_jesup_EC2",&mass4lj_jesup_EC2,"mass4lj_jesup_EC2/F"); 
    newtree->Branch("mass4ljj_jesup_EC2",&mass4ljj_jesup_EC2,"mass4ljj_jesup_EC2/F"); 
    newtree->Branch("pT4lj_jesup_EC2",&pT4lj_jesup_EC2,"pT4lj_jesup_EC2/F"); 
    newtree->Branch("pT4ljj_jesup_EC2",&pT4ljj_jesup_EC2,"pT4ljj_jesup_EC2/F"); 
    newtree->Branch("dyHj1_jesup_EC2",&dyHj1_jesup_EC2,"dyHj1_jesup_EC2/F");
    newtree->Branch("mj1j2_jesup_EC2",&mj1j2_jesup_EC2,"mj1j2_jesup_EC2/F"); 
    newtree->Branch("dEtaj1j2_jesup_EC2",&dEtaj1j2_jesup_EC2,"dEtaj1j2_jesup_EC2/F");
    newtree->Branch("dPhij1j2_jesup_EC2",&dPhij1j2_jesup_EC2,"dPhij1j2_jesup_EC2/F"); 
    newtree->Branch("dPhiHj1j2_jesup_EC2",&dPhiHj1j2_jesup_EC2,"dPhiHj1j2_jesup_EC2/F");
    newtree->Branch("pTj1_2p5_jesup_EC2",&pTj1_2p5_jesup_EC2,"pTj1_2p5_jesup_EC2/F"); 
    newtree->Branch("pt_leadingjet_pt30_eta2p5_jesup_EC2",&pt_leadingjet_pt30_eta2p5_jesup_EC2,"pt_leadingjet_pt30_eta2p5_jesup_EC2/F");
    newtree->Branch("yj1_2p5_jesup_EC2",&yj1_2p5_jesup_EC2,"yj1_2p5_jesup_EC2/F");
    newtree->Branch("pTj2_2p5_jesup_EC2",&pTj2_2p5_jesup_EC2,"pTj2_2p5_jesup_EC2/F"); 
    newtree->Branch("yj2_2p5_jesup_EC2",&yj2_2p5_jesup_EC2,"yj2_2p5_jesup_EC2/F");
    newtree->Branch("dPhiHj1_2p5_jesup_EC2",&dPhiHj1_2p5_jesup_EC2,"dPhiHj1_2p5_jesup_EC2/F"); 
    newtree->Branch("mass4lj_2p5_jesup_EC2",&mass4lj_2p5_jesup_EC2,"mass4lj_2p5_jesup_EC2/F");
    newtree->Branch("mass4ljj_2p5_jesup_EC2",&mass4ljj_2p5_jesup_EC2,"mass4ljj_2p5_jesup_EC2/F");
    newtree->Branch("pT4lj_2p5_jesup_EC2",&pT4lj_2p5_jesup_EC2,"pT4lj_2p5_jesup_EC2/F");
    newtree->Branch("pT4ljj_2p5_jesup_EC2",&pT4ljj_2p5_jesup_EC2,"pT4ljj_2p5_jesup_EC2/F");
    newtree->Branch("dyHj1_2p5_jesup_EC2",&dyHj1_2p5_jesup_EC2,"dyHj1_2p5_jesup_EC2/F");
    newtree->Branch("mj1j2_2p5_jesup_EC2",&mj1j2_2p5_jesup_EC2,"mj1j2_2p5_jesup_EC2/F"); 
    newtree->Branch("dEtaj1j2_2p5_jesup_EC2",&dEtaj1j2_2p5_jesup_EC2,"dEtaj1j2_2p5_jesup_EC2/F");
    newtree->Branch("dPhij1j2_2p5_jesup_EC2",&dPhij1j2_2p5_jesup_EC2,"dPhij1j2_2p5_jesup_EC2/F"); 
    newtree->Branch("dPhiHj1j2_2p5_jesup_EC2",&dPhiHj1j2_2p5_jesup_EC2,"dPhiHj1j2_2p5_jesup_EC2/F");
// jesdn_EC2
    newtree->Branch("njets_pt30_eta4p7_jesdn_EC2", &njets_pt30_eta4p7_jesdn_EC2, "njets_pt30_eta4p7_jesdn_EC2/F");
    newtree->Branch("TauC_Inc_0j_EnergyWgt_jesdn_EC2", &TauC_Inc_0j_EnergyWgt_jesdn_EC2, "TauC_Inc_0j_EnergyWgt_jesdn_EC2/F");
    newtree->Branch("TauB_Inc_0j_pTWgt_jesdn_EC2", &TauB_Inc_0j_pTWgt_jesdn_EC2, "TauB_Inc_0j_pTWgt_jesdn_EC2/F");
    newtree->Branch("njets_pt30_eta2p5_jesdn_EC2", &njets_pt30_eta2p5_jesdn_EC2, "njets_pt30_eta2p5_jesdn_EC2/F");
    newtree->Branch("pt_leadingjet_pt30_eta4p7_jesdn_EC2",&pt_leadingjet_pt30_eta4p7_jesdn_EC2,"pt_leadingjet_pt30_eta4p7_jesdn_EC2/F");
    newtree->Branch("pTj1_jesdn_EC2",&pTj1_jesdn_EC2,"pTj1_jesdn_EC2/F");
    newtree->Branch("etaj1_jesdn_EC2",&etaj1_jesdn_EC2,"etaj1_jesdn_EC2/F");
    newtree->Branch("pTj2_jesdn_EC2",&pTj2_jesdn_EC2,"pTj2_jesdn_EC2/F");
    newtree->Branch("etaj2_jesdn_EC2",&etaj2_jesdn_EC2,"etaj2_jesdn_EC2/F");
    newtree->Branch("yj1_jesdn_EC2",&yj1_jesdn_EC2,"yj1_jesdn_EC2/F");
    newtree->Branch("yj2_jesdn_EC2",&yj2_jesdn_EC2,"yj2_jesdn_EC2/F");
    newtree->Branch("dPhiHj1_jesdn_EC2",&dPhiHj1_jesdn_EC2,"dPhiHj1_jesdn_EC2/F"); 
    newtree->Branch("dyHj1_jesdn_EC2",&dyHj1_jesdn_EC2,"dyHj1_jesdn_EC2/F");
    newtree->Branch("mass4lj_jesdn_EC2",&mass4lj_jesdn_EC2,"mass4lj_jesdn_EC2/F");
    newtree->Branch("mass4ljj_jesdn_EC2",&mass4ljj_jesdn_EC2,"mass4ljj_jesdn_EC2/F");
    newtree->Branch("pT4lj_jesdn_EC2",&pT4lj_jesdn_EC2,"pT4lj_jesdn_EC2/F");
    newtree->Branch("pT4ljj_jesdn_EC2",&pT4ljj_jesdn_EC2,"pT4ljj_jesdn_EC2/F");
    newtree->Branch("mj1j2_jesdn_EC2",&mj1j2_jesdn_EC2,"mj1j2_jesdn_EC2/F"); 
    newtree->Branch("dEtaj1j2_jesdn_EC2",&dEtaj1j2_jesdn_EC2,"dEtaj1j2_jesdn_EC2/F");
    newtree->Branch("dPhij1j2_jesdn_EC2",&dPhij1j2_jesdn_EC2,"dPhij1j2_jesdn_EC2/F"); 
    newtree->Branch("dPhiHj1j2_jesdn_EC2",&dPhiHj1j2_jesdn_EC2,"dPhiHj1j2_jesdn_EC2/F");
    newtree->Branch("pTj1_2p5_jesdn_EC2",&pTj1_2p5_jesdn_EC2,"pTj1_2p5_jesdn_EC2/F"); 
    newtree->Branch("pt_leadingjet_pt30_eta2p5_jesdn_EC2",&pt_leadingjet_pt30_eta2p5_jesdn_EC2,"pt_leadingjet_pt30_eta2p5_2p5_jesdn_EC2/F"); 
    newtree->Branch("yj1_2p5_jesdn_EC2",&yj1_2p5_jesdn_EC2,"yj1_2p5_jesdn_EC2/F");
    newtree->Branch("pTj2_2p5_jesdn_EC2",&pTj2_2p5_jesdn_EC2,"pTj2_2p5_jesdn_EC2/F"); 
    newtree->Branch("yj2_2p5_jesdn_EC2",&yj2_2p5_jesdn_EC2,"yj2_2p5_jesdn_EC2/F");
    newtree->Branch("mass4lj_2p5_jesdn_EC2",&mass4lj_2p5_jesdn_EC2,"mass4lj_2p5_jesdn_EC2/F");
    newtree->Branch("mass4ljj_2p5_jesdn_EC2",&mass4ljj_2p5_jesdn_EC2,"mass4ljj_2p5_jesdn_EC2/F");
    newtree->Branch("pT4lj_2p5_jesdn_EC2",&pT4lj_2p5_jesdn_EC2,"pT4lj_2p5_jesdn_EC2/F");
    newtree->Branch("pT4ljj_2p5_jesdn_EC2",&pT4ljj_2p5_jesdn_EC2,"pT4ljj_2p5_jesdn_EC2/F");
    newtree->Branch("dPhiHj1_2p5_jesdn_EC2",&dPhiHj1_2p5_jesdn_EC2,"dPhiHj1_2p5_jesdn_EC2/F"); 
    newtree->Branch("dyHj1_2p5_jesdn_EC2",&dyHj1_2p5_jesdn_EC2,"dyHj1_2p5_jesdn_EC2/F");
    newtree->Branch("mj1j2_2p5_jesdn_EC2",&mj1j2_2p5_jesdn_EC2,"mj1j2_2p5_jesdn_EC2/F"); 
    newtree->Branch("dEtaj1j2_2p5_jesdn_EC2",&dEtaj1j2_2p5_jesdn_EC2,"dEtaj1j2_2p5_jesdn_EC2/F");
    newtree->Branch("dPhij1j2_2p5_jesdn_EC2",&dPhij1j2_2p5_jesdn_EC2,"dPhij1j2_2p5_jesdn_EC2/F"); 
    newtree->Branch("dPhiHj1j2_2p5_jesdn_EC2",&dPhiHj1j2_2p5_jesdn_EC2,"dPhiHj1j2_2p5_jesdn_EC2/F");
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// EC2_year JES 
// jesup_EC2_year
    newtree->Branch("njets_pt30_eta4p7_jesup_EC2_year", &njets_pt30_eta4p7_jesup_EC2_year, "njets_pt30_eta4p7_jesup_EC2_year/F");
    newtree->Branch("TauC_Inc_0j_EnergyWgt_jesup_EC2_year", &TauC_Inc_0j_EnergyWgt_jesup_EC2_year, "TauC_Inc_0j_EnergyWgt_jesup_EC2_year/F");
    newtree->Branch("TauB_Inc_0j_pTWgt_jesup_EC2_year", &TauB_Inc_0j_pTWgt_jesup_EC2_year, "TauB_Inc_0j_pTWgt_jesup_EC2_year/F");
    newtree->Branch("njets_pt30_eta2p5_jesup_EC2_year", &njets_pt30_eta2p5_jesup_EC2_year, "njets_pt30_eta2p5_jesup_EC2_year/F");
    newtree->Branch("pt_leadingjet_pt30_eta4p7_jesup_EC2_year",&pt_leadingjet_pt30_eta4p7_jesup_EC2_year,"pt_leadingjet_pt30_eta4p7_jesup_EC2_year/F");
    newtree->Branch("pTj1_jesup_EC2_year",&pTj1_jesup_EC2_year,"pTj1_jesup_EC2_year/F");
    newtree->Branch("etaj1_jesup_EC2_year",&etaj1_jesup_EC2_year,"etaj1_jesup_EC2_year/F");
    newtree->Branch("pTj2_jesup_EC2_year",&pTj2_jesup_EC2_year,"pTj2_jesup_EC2_year/F");
    newtree->Branch("etaj2_jesup_EC2_year",&etaj2_jesup_EC2_year,"etaj2_jesup_EC2_year/F");
    newtree->Branch("yj1_jesup_EC2_year",&yj1_jesup_EC2_year,"yj1_jesup_EC2_year/F");
    newtree->Branch("yj2_jesup_EC2_year",&yj2_jesup_EC2_year,"yj2_jesup_EC2_year/F");
    newtree->Branch("dPhiHj1_jesup_EC2_year",&dPhiHj1_jesup_EC2_year,"dPhiHj1_jesup_EC2_year/F"); 
    newtree->Branch("mass4lj_jesup_EC2_year",&mass4lj_jesup_EC2_year,"mass4lj_jesup_EC2_year/F"); 
    newtree->Branch("mass4ljj_jesup_EC2_year",&mass4ljj_jesup_EC2_year,"mass4ljj_jesup_EC2_year/F"); 
    newtree->Branch("pT4lj_jesup_EC2_year",&pT4lj_jesup_EC2_year,"pT4lj_jesup_EC2_year/F"); 
    newtree->Branch("pT4ljj_jesup_EC2_year",&pT4ljj_jesup_EC2_year,"pT4ljj_jesup_EC2_year/F"); 
    newtree->Branch("dyHj1_jesup_EC2_year",&dyHj1_jesup_EC2_year,"dyHj1_jesup_EC2_year/F");
    newtree->Branch("mj1j2_jesup_EC2_year",&mj1j2_jesup_EC2_year,"mj1j2_jesup_EC2_year/F"); 
    newtree->Branch("dEtaj1j2_jesup_EC2_year",&dEtaj1j2_jesup_EC2_year,"dEtaj1j2_jesup_EC2_year/F");
    newtree->Branch("dPhij1j2_jesup_EC2_year",&dPhij1j2_jesup_EC2_year,"dPhij1j2_jesup_EC2_year/F"); 
    newtree->Branch("dPhiHj1j2_jesup_EC2_year",&dPhiHj1j2_jesup_EC2_year,"dPhiHj1j2_jesup_EC2_year/F");
    newtree->Branch("pTj1_2p5_jesup_EC2_year",&pTj1_2p5_jesup_EC2_year,"pTj1_2p5_jesup_EC2_year/F"); 
    newtree->Branch("pt_leadingjet_pt30_eta2p5_jesup_EC2_year",&pt_leadingjet_pt30_eta2p5_jesup_EC2_year,"pt_leadingjet_pt30_eta2p5_jesup_EC2_year/F");
    newtree->Branch("yj1_2p5_jesup_EC2_year",&yj1_2p5_jesup_EC2_year,"yj1_2p5_jesup_EC2_year/F");
    newtree->Branch("pTj2_2p5_jesup_EC2_year",&pTj2_2p5_jesup_EC2_year,"pTj2_2p5_jesup_EC2_year/F"); 
    newtree->Branch("yj2_2p5_jesup_EC2_year",&yj2_2p5_jesup_EC2_year,"yj2_2p5_jesup_EC2_year/F");
    newtree->Branch("dPhiHj1_2p5_jesup_EC2_year",&dPhiHj1_2p5_jesup_EC2_year,"dPhiHj1_2p5_jesup_EC2_year/F"); 
    newtree->Branch("mass4lj_2p5_jesup_EC2_year",&mass4lj_2p5_jesup_EC2_year,"mass4lj_2p5_jesup_EC2_year/F");
    newtree->Branch("mass4ljj_2p5_jesup_EC2_year",&mass4ljj_2p5_jesup_EC2_year,"mass4ljj_2p5_jesup_EC2_year/F");
    newtree->Branch("pT4lj_2p5_jesup_EC2_year",&pT4lj_2p5_jesup_EC2_year,"pT4lj_2p5_jesup_EC2_year/F");
    newtree->Branch("pT4ljj_2p5_jesup_EC2_year",&pT4ljj_2p5_jesup_EC2_year,"pT4ljj_2p5_jesup_EC2_year/F");
    newtree->Branch("dyHj1_2p5_jesup_EC2_year",&dyHj1_2p5_jesup_EC2_year,"dyHj1_2p5_jesup_EC2_year/F");
    newtree->Branch("mj1j2_2p5_jesup_EC2_year",&mj1j2_2p5_jesup_EC2_year,"mj1j2_2p5_jesup_EC2_year/F"); 
    newtree->Branch("dEtaj1j2_2p5_jesup_EC2_year",&dEtaj1j2_2p5_jesup_EC2_year,"dEtaj1j2_2p5_jesup_EC2_year/F");
    newtree->Branch("dPhij1j2_2p5_jesup_EC2_year",&dPhij1j2_2p5_jesup_EC2_year,"dPhij1j2_2p5_jesup_EC2_year/F"); 
    newtree->Branch("dPhiHj1j2_2p5_jesup_EC2_year",&dPhiHj1j2_2p5_jesup_EC2_year,"dPhiHj1j2_2p5_jesup_EC2_year/F");
// jesdn_EC2_year
    newtree->Branch("njets_pt30_eta4p7_jesdn_EC2_year", &njets_pt30_eta4p7_jesdn_EC2_year, "njets_pt30_eta4p7_jesdn_EC2_year/F");
    newtree->Branch("TauC_Inc_0j_EnergyWgt_jesdn_EC2_year", &TauC_Inc_0j_EnergyWgt_jesdn_EC2_year, "TauC_Inc_0j_EnergyWgt_jesdn_EC2_year/F");
    newtree->Branch("TauB_Inc_0j_pTWgt_jesdn_EC2_year", &TauB_Inc_0j_pTWgt_jesdn_EC2_year, "TauB_Inc_0j_pTWgt_jesdn_EC2_year/F");
    newtree->Branch("njets_pt30_eta2p5_jesdn_EC2_year", &njets_pt30_eta2p5_jesdn_EC2_year, "njets_pt30_eta2p5_jesdn_EC2_year/F");
    newtree->Branch("pt_leadingjet_pt30_eta4p7_jesdn_EC2_year",&pt_leadingjet_pt30_eta4p7_jesdn_EC2_year,"pt_leadingjet_pt30_eta4p7_jesdn_EC2_year/F");
    newtree->Branch("pTj1_jesdn_EC2_year",&pTj1_jesdn_EC2_year,"pTj1_jesdn_EC2_year/F");
    newtree->Branch("etaj1_jesdn_EC2_year",&etaj1_jesdn_EC2_year,"etaj1_jesdn_EC2_year/F");
    newtree->Branch("pTj2_jesdn_EC2_year",&pTj2_jesdn_EC2_year,"pTj2_jesdn_EC2_year/F");
    newtree->Branch("etaj2_jesdn_EC2_year",&etaj2_jesdn_EC2_year,"etaj2_jesdn_EC2_year/F");
    newtree->Branch("yj1_jesdn_EC2_year",&yj1_jesdn_EC2_year,"yj1_jesdn_EC2_year/F");
    newtree->Branch("yj2_jesdn_EC2_year",&yj2_jesdn_EC2_year,"yj2_jesdn_EC2_year/F");
    newtree->Branch("dPhiHj1_jesdn_EC2_year",&dPhiHj1_jesdn_EC2_year,"dPhiHj1_jesdn_EC2_year/F"); 
    newtree->Branch("dyHj1_jesdn_EC2_year",&dyHj1_jesdn_EC2_year,"dyHj1_jesdn_EC2_year/F");
    newtree->Branch("mass4lj_jesdn_EC2_year",&mass4lj_jesdn_EC2_year,"mass4lj_jesdn_EC2_year/F");
    newtree->Branch("mass4ljj_jesdn_EC2_year",&mass4ljj_jesdn_EC2_year,"mass4ljj_jesdn_EC2_year/F");
    newtree->Branch("pT4lj_jesdn_EC2_year",&pT4lj_jesdn_EC2_year,"pT4lj_jesdn_EC2_year/F");
    newtree->Branch("pT4ljj_jesdn_EC2_year",&pT4ljj_jesdn_EC2_year,"pT4ljj_jesdn_EC2_year/F");
    newtree->Branch("mj1j2_jesdn_EC2_year",&mj1j2_jesdn_EC2_year,"mj1j2_jesdn_EC2_year/F"); 
    newtree->Branch("dEtaj1j2_jesdn_EC2_year",&dEtaj1j2_jesdn_EC2_year,"dEtaj1j2_jesdn_EC2_year/F");
    newtree->Branch("dPhij1j2_jesdn_EC2_year",&dPhij1j2_jesdn_EC2_year,"dPhij1j2_jesdn_EC2_year/F"); 
    newtree->Branch("dPhiHj1j2_jesdn_EC2_year",&dPhiHj1j2_jesdn_EC2_year,"dPhiHj1j2_jesdn_EC2_year/F");
    newtree->Branch("pTj1_2p5_jesdn_EC2_year",&pTj1_2p5_jesdn_EC2_year,"pTj1_2p5_jesdn_EC2_year/F"); 
    newtree->Branch("pt_leadingjet_pt30_eta2p5_jesdn_EC2_year",&pt_leadingjet_pt30_eta2p5_jesdn_EC2_year,"pt_leadingjet_pt30_eta2p5_2p5_jesdn_EC2_year/F"); 
    newtree->Branch("yj1_2p5_jesdn_EC2_year",&yj1_2p5_jesdn_EC2_year,"yj1_2p5_jesdn_EC2_year/F");
    newtree->Branch("pTj2_2p5_jesdn_EC2_year",&pTj2_2p5_jesdn_EC2_year,"pTj2_2p5_jesdn_EC2_year/F"); 
    newtree->Branch("yj2_2p5_jesdn_EC2_year",&yj2_2p5_jesdn_EC2_year,"yj2_2p5_jesdn_EC2_year/F");
    newtree->Branch("mass4lj_2p5_jesdn_EC2_year",&mass4lj_2p5_jesdn_EC2_year,"mass4lj_2p5_jesdn_EC2_year/F");
    newtree->Branch("mass4ljj_2p5_jesdn_EC2_year",&mass4ljj_2p5_jesdn_EC2_year,"mass4ljj_2p5_jesdn_EC2_year/F");
    newtree->Branch("pT4lj_2p5_jesdn_EC2_year",&pT4lj_2p5_jesdn_EC2_year,"pT4lj_2p5_jesdn_EC2_year/F");
    newtree->Branch("pT4ljj_2p5_jesdn_EC2_year",&pT4ljj_2p5_jesdn_EC2_year,"pT4ljj_2p5_jesdn_EC2_year/F");
    newtree->Branch("dPhiHj1_2p5_jesdn_EC2_year",&dPhiHj1_2p5_jesdn_EC2_year,"dPhiHj1_2p5_jesdn_EC2_year/F"); 
    newtree->Branch("dyHj1_2p5_jesdn_EC2_year",&dyHj1_2p5_jesdn_EC2_year,"dyHj1_2p5_jesdn_EC2_year/F");
    newtree->Branch("mj1j2_2p5_jesdn_EC2_year",&mj1j2_2p5_jesdn_EC2_year,"mj1j2_2p5_jesdn_EC2_year/F"); 
    newtree->Branch("dEtaj1j2_2p5_jesdn_EC2_year",&dEtaj1j2_2p5_jesdn_EC2_year,"dEtaj1j2_2p5_jesdn_EC2_year/F");
    newtree->Branch("dPhij1j2_2p5_jesdn_EC2_year",&dPhij1j2_2p5_jesdn_EC2_year,"dPhij1j2_2p5_jesdn_EC2_year/F"); 
    newtree->Branch("dPhiHj1j2_2p5_jesdn_EC2_year",&dPhiHj1j2_2p5_jesdn_EC2_year,"dPhiHj1j2_2p5_jesdn_EC2_year/F");
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// FlavQCD JES 
// jesup_FlavQCD
    newtree->Branch("njets_pt30_eta4p7_jesup_FlavQCD", &njets_pt30_eta4p7_jesup_FlavQCD, "njets_pt30_eta4p7_jesup_FlavQCD/F");
    newtree->Branch("TauC_Inc_0j_EnergyWgt_jesup_FlavQCD", &TauC_Inc_0j_EnergyWgt_jesup_FlavQCD, "TauC_Inc_0j_EnergyWgt_jesup_FlavQCD/F");
    newtree->Branch("TauB_Inc_0j_pTWgt_jesup_FlavQCD", &TauB_Inc_0j_pTWgt_jesup_FlavQCD, "TauB_Inc_0j_pTWgt_jesup_FlavQCD/F");
    newtree->Branch("njets_pt30_eta2p5_jesup_FlavQCD", &njets_pt30_eta2p5_jesup_FlavQCD, "njets_pt30_eta2p5_jesup_FlavQCD/F");
    newtree->Branch("pt_leadingjet_pt30_eta4p7_jesup_FlavQCD",&pt_leadingjet_pt30_eta4p7_jesup_FlavQCD,"pt_leadingjet_pt30_eta4p7_jesup_FlavQCD/F");
    newtree->Branch("pTj1_jesup_FlavQCD",&pTj1_jesup_FlavQCD,"pTj1_jesup_FlavQCD/F");
    newtree->Branch("etaj1_jesup_FlavQCD",&etaj1_jesup_FlavQCD,"etaj1_jesup_FlavQCD/F");
    newtree->Branch("pTj2_jesup_FlavQCD",&pTj2_jesup_FlavQCD,"pTj2_jesup_FlavQCD/F");
    newtree->Branch("etaj2_jesup_FlavQCD",&etaj2_jesup_FlavQCD,"etaj2_jesup_FlavQCD/F");
    newtree->Branch("yj1_jesup_FlavQCD",&yj1_jesup_FlavQCD,"yj1_jesup_FlavQCD/F");
    newtree->Branch("yj2_jesup_FlavQCD",&yj2_jesup_FlavQCD,"yj2_jesup_FlavQCD/F");
    newtree->Branch("dPhiHj1_jesup_FlavQCD",&dPhiHj1_jesup_FlavQCD,"dPhiHj1_jesup_FlavQCD/F"); 
    newtree->Branch("mass4lj_jesup_FlavQCD",&mass4lj_jesup_FlavQCD,"mass4lj_jesup_FlavQCD/F"); 
    newtree->Branch("mass4ljj_jesup_FlavQCD",&mass4ljj_jesup_FlavQCD,"mass4ljj_jesup_FlavQCD/F"); 
    newtree->Branch("pT4lj_jesup_FlavQCD",&pT4lj_jesup_FlavQCD,"pT4lj_jesup_FlavQCD/F"); 
    newtree->Branch("pT4ljj_jesup_FlavQCD",&pT4ljj_jesup_FlavQCD,"pT4ljj_jesup_FlavQCD/F"); 
    newtree->Branch("dyHj1_jesup_FlavQCD",&dyHj1_jesup_FlavQCD,"dyHj1_jesup_FlavQCD/F");
    newtree->Branch("mj1j2_jesup_FlavQCD",&mj1j2_jesup_FlavQCD,"mj1j2_jesup_FlavQCD/F"); 
    newtree->Branch("dEtaj1j2_jesup_FlavQCD",&dEtaj1j2_jesup_FlavQCD,"dEtaj1j2_jesup_FlavQCD/F");
    newtree->Branch("dPhij1j2_jesup_FlavQCD",&dPhij1j2_jesup_FlavQCD,"dPhij1j2_jesup_FlavQCD/F"); 
    newtree->Branch("dPhiHj1j2_jesup_FlavQCD",&dPhiHj1j2_jesup_FlavQCD,"dPhiHj1j2_jesup_FlavQCD/F");
    newtree->Branch("pTj1_2p5_jesup_FlavQCD",&pTj1_2p5_jesup_FlavQCD,"pTj1_2p5_jesup_FlavQCD/F"); 
    newtree->Branch("pt_leadingjet_pt30_eta2p5_jesup_FlavQCD",&pt_leadingjet_pt30_eta2p5_jesup_FlavQCD,"pt_leadingjet_pt30_eta2p5_jesup_FlavQCD/F");
    newtree->Branch("yj1_2p5_jesup_FlavQCD",&yj1_2p5_jesup_FlavQCD,"yj1_2p5_jesup_FlavQCD/F");
    newtree->Branch("pTj2_2p5_jesup_FlavQCD",&pTj2_2p5_jesup_FlavQCD,"pTj2_2p5_jesup_FlavQCD/F"); 
    newtree->Branch("yj2_2p5_jesup_FlavQCD",&yj2_2p5_jesup_FlavQCD,"yj2_2p5_jesup_FlavQCD/F");
    newtree->Branch("dPhiHj1_2p5_jesup_FlavQCD",&dPhiHj1_2p5_jesup_FlavQCD,"dPhiHj1_2p5_jesup_FlavQCD/F"); 
    newtree->Branch("mass4lj_2p5_jesup_FlavQCD",&mass4lj_2p5_jesup_FlavQCD,"mass4lj_2p5_jesup_FlavQCD/F");
    newtree->Branch("mass4ljj_2p5_jesup_FlavQCD",&mass4ljj_2p5_jesup_FlavQCD,"mass4ljj_2p5_jesup_FlavQCD/F");
    newtree->Branch("pT4lj_2p5_jesup_FlavQCD",&pT4lj_2p5_jesup_FlavQCD,"pT4lj_2p5_jesup_FlavQCD/F");
    newtree->Branch("pT4ljj_2p5_jesup_FlavQCD",&pT4ljj_2p5_jesup_FlavQCD,"pT4ljj_2p5_jesup_FlavQCD/F");
    newtree->Branch("dyHj1_2p5_jesup_FlavQCD",&dyHj1_2p5_jesup_FlavQCD,"dyHj1_2p5_jesup_FlavQCD/F");
    newtree->Branch("mj1j2_2p5_jesup_FlavQCD",&mj1j2_2p5_jesup_FlavQCD,"mj1j2_2p5_jesup_FlavQCD/F"); 
    newtree->Branch("dEtaj1j2_2p5_jesup_FlavQCD",&dEtaj1j2_2p5_jesup_FlavQCD,"dEtaj1j2_2p5_jesup_FlavQCD/F");
    newtree->Branch("dPhij1j2_2p5_jesup_FlavQCD",&dPhij1j2_2p5_jesup_FlavQCD,"dPhij1j2_2p5_jesup_FlavQCD/F"); 
    newtree->Branch("dPhiHj1j2_2p5_jesup_FlavQCD",&dPhiHj1j2_2p5_jesup_FlavQCD,"dPhiHj1j2_2p5_jesup_FlavQCD/F");
// jesdn_FlavQCD
    newtree->Branch("njets_pt30_eta4p7_jesdn_FlavQCD", &njets_pt30_eta4p7_jesdn_FlavQCD, "njets_pt30_eta4p7_jesdn_FlavQCD/F");
    newtree->Branch("TauC_Inc_0j_EnergyWgt_jesdn_FlavQCD", &TauC_Inc_0j_EnergyWgt_jesdn_FlavQCD, "TauC_Inc_0j_EnergyWgt_jesdn_FlavQCD/F");
    newtree->Branch("TauB_Inc_0j_pTWgt_jesdn_FlavQCD", &TauB_Inc_0j_pTWgt_jesdn_FlavQCD, "TauB_Inc_0j_pTWgt_jesdn_FlavQCD/F");
    newtree->Branch("njets_pt30_eta2p5_jesdn_FlavQCD", &njets_pt30_eta2p5_jesdn_FlavQCD, "njets_pt30_eta2p5_jesdn_FlavQCD/F");
    newtree->Branch("pt_leadingjet_pt30_eta4p7_jesdn_FlavQCD",&pt_leadingjet_pt30_eta4p7_jesdn_FlavQCD,"pt_leadingjet_pt30_eta4p7_jesdn_FlavQCD/F");
    newtree->Branch("pTj1_jesdn_FlavQCD",&pTj1_jesdn_FlavQCD,"pTj1_jesdn_FlavQCD/F");
    newtree->Branch("etaj1_jesdn_FlavQCD",&etaj1_jesdn_FlavQCD,"etaj1_jesdn_FlavQCD/F");
    newtree->Branch("pTj2_jesdn_FlavQCD",&pTj2_jesdn_FlavQCD,"pTj2_jesdn_FlavQCD/F");
    newtree->Branch("etaj2_jesdn_FlavQCD",&etaj2_jesdn_FlavQCD,"etaj2_jesdn_FlavQCD/F");
    newtree->Branch("yj1_jesdn_FlavQCD",&yj1_jesdn_FlavQCD,"yj1_jesdn_FlavQCD/F");
    newtree->Branch("yj2_jesdn_FlavQCD",&yj2_jesdn_FlavQCD,"yj2_jesdn_FlavQCD/F");
    newtree->Branch("dPhiHj1_jesdn_FlavQCD",&dPhiHj1_jesdn_FlavQCD,"dPhiHj1_jesdn_FlavQCD/F"); 
    newtree->Branch("dyHj1_jesdn_FlavQCD",&dyHj1_jesdn_FlavQCD,"dyHj1_jesdn_FlavQCD/F");
    newtree->Branch("mass4lj_jesdn_FlavQCD",&mass4lj_jesdn_FlavQCD,"mass4lj_jesdn_FlavQCD/F");
    newtree->Branch("mass4ljj_jesdn_FlavQCD",&mass4ljj_jesdn_FlavQCD,"mass4ljj_jesdn_FlavQCD/F");
    newtree->Branch("pT4lj_jesdn_FlavQCD",&pT4lj_jesdn_FlavQCD,"pT4lj_jesdn_FlavQCD/F");
    newtree->Branch("pT4ljj_jesdn_FlavQCD",&pT4ljj_jesdn_FlavQCD,"pT4ljj_jesdn_FlavQCD/F");
    newtree->Branch("mj1j2_jesdn_FlavQCD",&mj1j2_jesdn_FlavQCD,"mj1j2_jesdn_FlavQCD/F"); 
    newtree->Branch("dEtaj1j2_jesdn_FlavQCD",&dEtaj1j2_jesdn_FlavQCD,"dEtaj1j2_jesdn_FlavQCD/F");
    newtree->Branch("dPhij1j2_jesdn_FlavQCD",&dPhij1j2_jesdn_FlavQCD,"dPhij1j2_jesdn_FlavQCD/F"); 
    newtree->Branch("dPhiHj1j2_jesdn_FlavQCD",&dPhiHj1j2_jesdn_FlavQCD,"dPhiHj1j2_jesdn_FlavQCD/F");
    newtree->Branch("pTj1_2p5_jesdn_FlavQCD",&pTj1_2p5_jesdn_FlavQCD,"pTj1_2p5_jesdn_FlavQCD/F"); 
    newtree->Branch("pt_leadingjet_pt30_eta2p5_jesdn_FlavQCD",&pt_leadingjet_pt30_eta2p5_jesdn_FlavQCD,"pt_leadingjet_pt30_eta2p5_2p5_jesdn_FlavQCD/F"); 
    newtree->Branch("yj1_2p5_jesdn_FlavQCD",&yj1_2p5_jesdn_FlavQCD,"yj1_2p5_jesdn_FlavQCD/F");
    newtree->Branch("pTj2_2p5_jesdn_FlavQCD",&pTj2_2p5_jesdn_FlavQCD,"pTj2_2p5_jesdn_FlavQCD/F"); 
    newtree->Branch("yj2_2p5_jesdn_FlavQCD",&yj2_2p5_jesdn_FlavQCD,"yj2_2p5_jesdn_FlavQCD/F");
    newtree->Branch("mass4lj_2p5_jesdn_FlavQCD",&mass4lj_2p5_jesdn_FlavQCD,"mass4lj_2p5_jesdn_FlavQCD/F");
    newtree->Branch("mass4ljj_2p5_jesdn_FlavQCD",&mass4ljj_2p5_jesdn_FlavQCD,"mass4ljj_2p5_jesdn_FlavQCD/F");
    newtree->Branch("pT4lj_2p5_jesdn_FlavQCD",&pT4lj_2p5_jesdn_FlavQCD,"pT4lj_2p5_jesdn_FlavQCD/F");
    newtree->Branch("pT4ljj_2p5_jesdn_FlavQCD",&pT4ljj_2p5_jesdn_FlavQCD,"pT4ljj_2p5_jesdn_FlavQCD/F");
    newtree->Branch("dPhiHj1_2p5_jesdn_FlavQCD",&dPhiHj1_2p5_jesdn_FlavQCD,"dPhiHj1_2p5_jesdn_FlavQCD/F"); 
    newtree->Branch("dyHj1_2p5_jesdn_FlavQCD",&dyHj1_2p5_jesdn_FlavQCD,"dyHj1_2p5_jesdn_FlavQCD/F");
    newtree->Branch("mj1j2_2p5_jesdn_FlavQCD",&mj1j2_2p5_jesdn_FlavQCD,"mj1j2_2p5_jesdn_FlavQCD/F"); 
    newtree->Branch("dEtaj1j2_2p5_jesdn_FlavQCD",&dEtaj1j2_2p5_jesdn_FlavQCD,"dEtaj1j2_2p5_jesdn_FlavQCD/F");
    newtree->Branch("dPhij1j2_2p5_jesdn_FlavQCD",&dPhij1j2_2p5_jesdn_FlavQCD,"dPhij1j2_2p5_jesdn_FlavQCD/F"); 
    newtree->Branch("dPhiHj1j2_2p5_jesdn_FlavQCD",&dPhiHj1j2_2p5_jesdn_FlavQCD,"dPhiHj1j2_2p5_jesdn_FlavQCD/F");
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// HF JES 
// jesup_HF
    newtree->Branch("njets_pt30_eta4p7_jesup_HF", &njets_pt30_eta4p7_jesup_HF, "njets_pt30_eta4p7_jesup_HF/F");
    newtree->Branch("TauC_Inc_0j_EnergyWgt_jesup_HF", &TauC_Inc_0j_EnergyWgt_jesup_HF, "TauC_Inc_0j_EnergyWgt_jesup_HF/F");
    newtree->Branch("TauB_Inc_0j_pTWgt_jesup_HF", &TauB_Inc_0j_pTWgt_jesup_HF, "TauB_Inc_0j_pTWgt_jesup_HF/F");
    newtree->Branch("njets_pt30_eta2p5_jesup_HF", &njets_pt30_eta2p5_jesup_HF, "njets_pt30_eta2p5_jesup_HF/F");
    newtree->Branch("pt_leadingjet_pt30_eta4p7_jesup_HF",&pt_leadingjet_pt30_eta4p7_jesup_HF,"pt_leadingjet_pt30_eta4p7_jesup_HF/F");
    newtree->Branch("pTj1_jesup_HF",&pTj1_jesup_HF,"pTj1_jesup_HF/F");
    newtree->Branch("etaj1_jesup_HF",&etaj1_jesup_HF,"etaj1_jesup_HF/F");
    newtree->Branch("pTj2_jesup_HF",&pTj2_jesup_HF,"pTj2_jesup_HF/F");
    newtree->Branch("etaj2_jesup_HF",&etaj2_jesup_HF,"etaj2_jesup_HF/F");
    newtree->Branch("yj1_jesup_HF",&yj1_jesup_HF,"yj1_jesup_HF/F");
    newtree->Branch("yj2_jesup_HF",&yj2_jesup_HF,"yj2_jesup_HF/F");
    newtree->Branch("dPhiHj1_jesup_HF",&dPhiHj1_jesup_HF,"dPhiHj1_jesup_HF/F"); 
    newtree->Branch("mass4lj_jesup_HF",&mass4lj_jesup_HF,"mass4lj_jesup_HF/F"); 
    newtree->Branch("mass4ljj_jesup_HF",&mass4ljj_jesup_HF,"mass4ljj_jesup_HF/F"); 
    newtree->Branch("pT4lj_jesup_HF",&pT4lj_jesup_HF,"pT4lj_jesup_HF/F"); 
    newtree->Branch("pT4ljj_jesup_HF",&pT4ljj_jesup_HF,"pT4ljj_jesup_HF/F"); 
    newtree->Branch("dyHj1_jesup_HF",&dyHj1_jesup_HF,"dyHj1_jesup_HF/F");
    newtree->Branch("mj1j2_jesup_HF",&mj1j2_jesup_HF,"mj1j2_jesup_HF/F"); 
    newtree->Branch("dEtaj1j2_jesup_HF",&dEtaj1j2_jesup_HF,"dEtaj1j2_jesup_HF/F");
    newtree->Branch("dPhij1j2_jesup_HF",&dPhij1j2_jesup_HF,"dPhij1j2_jesup_HF/F"); 
    newtree->Branch("dPhiHj1j2_jesup_HF",&dPhiHj1j2_jesup_HF,"dPhiHj1j2_jesup_HF/F");
    newtree->Branch("pTj1_2p5_jesup_HF",&pTj1_2p5_jesup_HF,"pTj1_2p5_jesup_HF/F"); 
    newtree->Branch("pt_leadingjet_pt30_eta2p5_jesup_HF",&pt_leadingjet_pt30_eta2p5_jesup_HF,"pt_leadingjet_pt30_eta2p5_jesup_HF/F");
    newtree->Branch("yj1_2p5_jesup_HF",&yj1_2p5_jesup_HF,"yj1_2p5_jesup_HF/F");
    newtree->Branch("pTj2_2p5_jesup_HF",&pTj2_2p5_jesup_HF,"pTj2_2p5_jesup_HF/F"); 
    newtree->Branch("yj2_2p5_jesup_HF",&yj2_2p5_jesup_HF,"yj2_2p5_jesup_HF/F");
    newtree->Branch("dPhiHj1_2p5_jesup_HF",&dPhiHj1_2p5_jesup_HF,"dPhiHj1_2p5_jesup_HF/F"); 
    newtree->Branch("mass4lj_2p5_jesup_HF",&mass4lj_2p5_jesup_HF,"mass4lj_2p5_jesup_HF/F");
    newtree->Branch("mass4ljj_2p5_jesup_HF",&mass4ljj_2p5_jesup_HF,"mass4ljj_2p5_jesup_HF/F");
    newtree->Branch("pT4lj_2p5_jesup_HF",&pT4lj_2p5_jesup_HF,"pT4lj_2p5_jesup_HF/F");
    newtree->Branch("pT4ljj_2p5_jesup_HF",&pT4ljj_2p5_jesup_HF,"pT4ljj_2p5_jesup_HF/F");
    newtree->Branch("dyHj1_2p5_jesup_HF",&dyHj1_2p5_jesup_HF,"dyHj1_2p5_jesup_HF/F");
    newtree->Branch("mj1j2_2p5_jesup_HF",&mj1j2_2p5_jesup_HF,"mj1j2_2p5_jesup_HF/F"); 
    newtree->Branch("dEtaj1j2_2p5_jesup_HF",&dEtaj1j2_2p5_jesup_HF,"dEtaj1j2_2p5_jesup_HF/F");
    newtree->Branch("dPhij1j2_2p5_jesup_HF",&dPhij1j2_2p5_jesup_HF,"dPhij1j2_2p5_jesup_HF/F"); 
    newtree->Branch("dPhiHj1j2_2p5_jesup_HF",&dPhiHj1j2_2p5_jesup_HF,"dPhiHj1j2_2p5_jesup_HF/F");
// jesdn_HF
    newtree->Branch("njets_pt30_eta4p7_jesdn_HF", &njets_pt30_eta4p7_jesdn_HF, "njets_pt30_eta4p7_jesdn_HF/F");
    newtree->Branch("TauC_Inc_0j_EnergyWgt_jesdn_HF", &TauC_Inc_0j_EnergyWgt_jesdn_HF, "TauC_Inc_0j_EnergyWgt_jesdn_HF/F");
    newtree->Branch("TauB_Inc_0j_pTWgt_jesdn_HF", &TauB_Inc_0j_pTWgt_jesdn_HF, "TauB_Inc_0j_pTWgt_jesdn_HF/F");
    newtree->Branch("njets_pt30_eta2p5_jesdn_HF", &njets_pt30_eta2p5_jesdn_HF, "njets_pt30_eta2p5_jesdn_HF/F");
    newtree->Branch("pt_leadingjet_pt30_eta4p7_jesdn_HF",&pt_leadingjet_pt30_eta4p7_jesdn_HF,"pt_leadingjet_pt30_eta4p7_jesdn_HF/F");
    newtree->Branch("pTj1_jesdn_HF",&pTj1_jesdn_HF,"pTj1_jesdn_HF/F");
    newtree->Branch("etaj1_jesdn_HF",&etaj1_jesdn_HF,"etaj1_jesdn_HF/F");
    newtree->Branch("pTj2_jesdn_HF",&pTj2_jesdn_HF,"pTj2_jesdn_HF/F");
    newtree->Branch("etaj2_jesdn_HF",&etaj2_jesdn_HF,"etaj2_jesdn_HF/F");
    newtree->Branch("yj1_jesdn_HF",&yj1_jesdn_HF,"yj1_jesdn_HF/F");
    newtree->Branch("yj2_jesdn_HF",&yj2_jesdn_HF,"yj2_jesdn_HF/F");
    newtree->Branch("dPhiHj1_jesdn_HF",&dPhiHj1_jesdn_HF,"dPhiHj1_jesdn_HF/F"); 
    newtree->Branch("dyHj1_jesdn_HF",&dyHj1_jesdn_HF,"dyHj1_jesdn_HF/F");
    newtree->Branch("mass4lj_jesdn_HF",&mass4lj_jesdn_HF,"mass4lj_jesdn_HF/F");
    newtree->Branch("mass4ljj_jesdn_HF",&mass4ljj_jesdn_HF,"mass4ljj_jesdn_HF/F");
    newtree->Branch("pT4lj_jesdn_HF",&pT4lj_jesdn_HF,"pT4lj_jesdn_HF/F");
    newtree->Branch("pT4ljj_jesdn_HF",&pT4ljj_jesdn_HF,"pT4ljj_jesdn_HF/F");
    newtree->Branch("mj1j2_jesdn_HF",&mj1j2_jesdn_HF,"mj1j2_jesdn_HF/F"); 
    newtree->Branch("dEtaj1j2_jesdn_HF",&dEtaj1j2_jesdn_HF,"dEtaj1j2_jesdn_HF/F");
    newtree->Branch("dPhij1j2_jesdn_HF",&dPhij1j2_jesdn_HF,"dPhij1j2_jesdn_HF/F"); 
    newtree->Branch("dPhiHj1j2_jesdn_HF",&dPhiHj1j2_jesdn_HF,"dPhiHj1j2_jesdn_HF/F");
    newtree->Branch("pTj1_2p5_jesdn_HF",&pTj1_2p5_jesdn_HF,"pTj1_2p5_jesdn_HF/F"); 
    newtree->Branch("pt_leadingjet_pt30_eta2p5_jesdn_HF",&pt_leadingjet_pt30_eta2p5_jesdn_HF,"pt_leadingjet_pt30_eta2p5_2p5_jesdn_HF/F"); 
    newtree->Branch("yj1_2p5_jesdn_HF",&yj1_2p5_jesdn_HF,"yj1_2p5_jesdn_HF/F");
    newtree->Branch("pTj2_2p5_jesdn_HF",&pTj2_2p5_jesdn_HF,"pTj2_2p5_jesdn_HF/F"); 
    newtree->Branch("yj2_2p5_jesdn_HF",&yj2_2p5_jesdn_HF,"yj2_2p5_jesdn_HF/F");
    newtree->Branch("mass4lj_2p5_jesdn_HF",&mass4lj_2p5_jesdn_HF,"mass4lj_2p5_jesdn_HF/F");
    newtree->Branch("mass4ljj_2p5_jesdn_HF",&mass4ljj_2p5_jesdn_HF,"mass4ljj_2p5_jesdn_HF/F");
    newtree->Branch("pT4lj_2p5_jesdn_HF",&pT4lj_2p5_jesdn_HF,"pT4lj_2p5_jesdn_HF/F");
    newtree->Branch("pT4ljj_2p5_jesdn_HF",&pT4ljj_2p5_jesdn_HF,"pT4ljj_2p5_jesdn_HF/F");
    newtree->Branch("dPhiHj1_2p5_jesdn_HF",&dPhiHj1_2p5_jesdn_HF,"dPhiHj1_2p5_jesdn_HF/F"); 
    newtree->Branch("dyHj1_2p5_jesdn_HF",&dyHj1_2p5_jesdn_HF,"dyHj1_2p5_jesdn_HF/F");
    newtree->Branch("mj1j2_2p5_jesdn_HF",&mj1j2_2p5_jesdn_HF,"mj1j2_2p5_jesdn_HF/F"); 
    newtree->Branch("dEtaj1j2_2p5_jesdn_HF",&dEtaj1j2_2p5_jesdn_HF,"dEtaj1j2_2p5_jesdn_HF/F");
    newtree->Branch("dPhij1j2_2p5_jesdn_HF",&dPhij1j2_2p5_jesdn_HF,"dPhij1j2_2p5_jesdn_HF/F"); 
    newtree->Branch("dPhiHj1j2_2p5_jesdn_HF",&dPhiHj1j2_2p5_jesdn_HF,"dPhiHj1j2_2p5_jesdn_HF/F");
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// HF_year JES 
// jesup_HF_year
    newtree->Branch("njets_pt30_eta4p7_jesup_HF_year", &njets_pt30_eta4p7_jesup_HF_year, "njets_pt30_eta4p7_jesup_HF_year/F");
    newtree->Branch("TauC_Inc_0j_EnergyWgt_jesup_HF_year", &TauC_Inc_0j_EnergyWgt_jesup_HF_year, "TauC_Inc_0j_EnergyWgt_jesup_HF_year/F");
    newtree->Branch("TauB_Inc_0j_pTWgt_jesup_HF_year", &TauB_Inc_0j_pTWgt_jesup_HF_year, "TauB_Inc_0j_pTWgt_jesup_HF_year/F");
    newtree->Branch("njets_pt30_eta2p5_jesup_HF_year", &njets_pt30_eta2p5_jesup_HF_year, "njets_pt30_eta2p5_jesup_HF_year/F");
    newtree->Branch("pt_leadingjet_pt30_eta4p7_jesup_HF_year",&pt_leadingjet_pt30_eta4p7_jesup_HF_year,"pt_leadingjet_pt30_eta4p7_jesup_HF_year/F");
    newtree->Branch("pTj1_jesup_HF_year",&pTj1_jesup_HF_year,"pTj1_jesup_HF_year/F");
    newtree->Branch("etaj1_jesup_HF_year",&etaj1_jesup_HF_year,"etaj1_jesup_HF_year/F");
    newtree->Branch("pTj2_jesup_HF_year",&pTj2_jesup_HF_year,"pTj2_jesup_HF_year/F");
    newtree->Branch("etaj2_jesup_HF_year",&etaj2_jesup_HF_year,"etaj2_jesup_HF_year/F");
    newtree->Branch("yj1_jesup_HF_year",&yj1_jesup_HF_year,"yj1_jesup_HF_year/F");
    newtree->Branch("yj2_jesup_HF_year",&yj2_jesup_HF_year,"yj2_jesup_HF_year/F");
    newtree->Branch("dPhiHj1_jesup_HF_year",&dPhiHj1_jesup_HF_year,"dPhiHj1_jesup_HF_year/F"); 
    newtree->Branch("mass4lj_jesup_HF_year",&mass4lj_jesup_HF_year,"mass4lj_jesup_HF_year/F"); 
    newtree->Branch("mass4ljj_jesup_HF_year",&mass4ljj_jesup_HF_year,"mass4ljj_jesup_HF_year/F"); 
    newtree->Branch("pT4lj_jesup_HF_year",&pT4lj_jesup_HF_year,"pT4lj_jesup_HF_year/F"); 
    newtree->Branch("pT4ljj_jesup_HF_year",&pT4ljj_jesup_HF_year,"pT4ljj_jesup_HF_year/F"); 
    newtree->Branch("dyHj1_jesup_HF_year",&dyHj1_jesup_HF_year,"dyHj1_jesup_HF_year/F");
    newtree->Branch("mj1j2_jesup_HF_year",&mj1j2_jesup_HF_year,"mj1j2_jesup_HF_year/F"); 
    newtree->Branch("dEtaj1j2_jesup_HF_year",&dEtaj1j2_jesup_HF_year,"dEtaj1j2_jesup_HF_year/F");
    newtree->Branch("dPhij1j2_jesup_HF_year",&dPhij1j2_jesup_HF_year,"dPhij1j2_jesup_HF_year/F"); 
    newtree->Branch("dPhiHj1j2_jesup_HF_year",&dPhiHj1j2_jesup_HF_year,"dPhiHj1j2_jesup_HF_year/F");
    newtree->Branch("pTj1_2p5_jesup_HF_year",&pTj1_2p5_jesup_HF_year,"pTj1_2p5_jesup_HF_year/F"); 
    newtree->Branch("pt_leadingjet_pt30_eta2p5_jesup_HF_year",&pt_leadingjet_pt30_eta2p5_jesup_HF_year,"pt_leadingjet_pt30_eta2p5_jesup_HF_year/F");
    newtree->Branch("yj1_2p5_jesup_HF_year",&yj1_2p5_jesup_HF_year,"yj1_2p5_jesup_HF_year/F");
    newtree->Branch("pTj2_2p5_jesup_HF_year",&pTj2_2p5_jesup_HF_year,"pTj2_2p5_jesup_HF_year/F"); 
    newtree->Branch("yj2_2p5_jesup_HF_year",&yj2_2p5_jesup_HF_year,"yj2_2p5_jesup_HF_year/F");
    newtree->Branch("dPhiHj1_2p5_jesup_HF_year",&dPhiHj1_2p5_jesup_HF_year,"dPhiHj1_2p5_jesup_HF_year/F"); 
    newtree->Branch("mass4lj_2p5_jesup_HF_year",&mass4lj_2p5_jesup_HF_year,"mass4lj_2p5_jesup_HF_year/F");
    newtree->Branch("mass4ljj_2p5_jesup_HF_year",&mass4ljj_2p5_jesup_HF_year,"mass4ljj_2p5_jesup_HF_year/F");
    newtree->Branch("pT4lj_2p5_jesup_HF_year",&pT4lj_2p5_jesup_HF_year,"pT4lj_2p5_jesup_HF_year/F");
    newtree->Branch("pT4ljj_2p5_jesup_HF_year",&pT4ljj_2p5_jesup_HF_year,"pT4ljj_2p5_jesup_HF_year/F");
    newtree->Branch("dyHj1_2p5_jesup_HF_year",&dyHj1_2p5_jesup_HF_year,"dyHj1_2p5_jesup_HF_year/F");
    newtree->Branch("mj1j2_2p5_jesup_HF_year",&mj1j2_2p5_jesup_HF_year,"mj1j2_2p5_jesup_HF_year/F"); 
    newtree->Branch("dEtaj1j2_2p5_jesup_HF_year",&dEtaj1j2_2p5_jesup_HF_year,"dEtaj1j2_2p5_jesup_HF_year/F");
    newtree->Branch("dPhij1j2_2p5_jesup_HF_year",&dPhij1j2_2p5_jesup_HF_year,"dPhij1j2_2p5_jesup_HF_year/F"); 
    newtree->Branch("dPhiHj1j2_2p5_jesup_HF_year",&dPhiHj1j2_2p5_jesup_HF_year,"dPhiHj1j2_2p5_jesup_HF_year/F");
// jesdn_HF_year
    newtree->Branch("njets_pt30_eta4p7_jesdn_HF_year", &njets_pt30_eta4p7_jesdn_HF_year, "njets_pt30_eta4p7_jesdn_HF_year/F");
    newtree->Branch("TauC_Inc_0j_EnergyWgt_jesdn_HF_year", &TauC_Inc_0j_EnergyWgt_jesdn_HF_year, "TauC_Inc_0j_EnergyWgt_jesdn_HF_year/F");
    newtree->Branch("TauB_Inc_0j_pTWgt_jesdn_HF_year", &TauB_Inc_0j_pTWgt_jesdn_HF_year, "TauB_Inc_0j_pTWgt_jesdn_HF_year/F");
    newtree->Branch("njets_pt30_eta2p5_jesdn_HF_year", &njets_pt30_eta2p5_jesdn_HF_year, "njets_pt30_eta2p5_jesdn_HF_year/F");
    newtree->Branch("pt_leadingjet_pt30_eta4p7_jesdn_HF_year",&pt_leadingjet_pt30_eta4p7_jesdn_HF_year,"pt_leadingjet_pt30_eta4p7_jesdn_HF_year/F");
    newtree->Branch("pTj1_jesdn_HF_year",&pTj1_jesdn_HF_year,"pTj1_jesdn_HF_year/F");
    newtree->Branch("etaj1_jesdn_HF_year",&etaj1_jesdn_HF_year,"etaj1_jesdn_HF_year/F");
    newtree->Branch("pTj2_jesdn_HF_year",&pTj2_jesdn_HF_year,"pTj2_jesdn_HF_year/F");
    newtree->Branch("etaj2_jesdn_HF_year",&etaj2_jesdn_HF_year,"etaj2_jesdn_HF_year/F");
    newtree->Branch("yj1_jesdn_HF_year",&yj1_jesdn_HF_year,"yj1_jesdn_HF_year/F");
    newtree->Branch("yj2_jesdn_HF_year",&yj2_jesdn_HF_year,"yj2_jesdn_HF_year/F");
    newtree->Branch("dPhiHj1_jesdn_HF_year",&dPhiHj1_jesdn_HF_year,"dPhiHj1_jesdn_HF_year/F"); 
    newtree->Branch("dyHj1_jesdn_HF_year",&dyHj1_jesdn_HF_year,"dyHj1_jesdn_HF_year/F");
    newtree->Branch("mass4lj_jesdn_HF_year",&mass4lj_jesdn_HF_year,"mass4lj_jesdn_HF_year/F");
    newtree->Branch("mass4ljj_jesdn_HF_year",&mass4ljj_jesdn_HF_year,"mass4ljj_jesdn_HF_year/F");
    newtree->Branch("pT4lj_jesdn_HF_year",&pT4lj_jesdn_HF_year,"pT4lj_jesdn_HF_year/F");
    newtree->Branch("pT4ljj_jesdn_HF_year",&pT4ljj_jesdn_HF_year,"pT4ljj_jesdn_HF_year/F");
    newtree->Branch("mj1j2_jesdn_HF_year",&mj1j2_jesdn_HF_year,"mj1j2_jesdn_HF_year/F"); 
    newtree->Branch("dEtaj1j2_jesdn_HF_year",&dEtaj1j2_jesdn_HF_year,"dEtaj1j2_jesdn_HF_year/F");
    newtree->Branch("dPhij1j2_jesdn_HF_year",&dPhij1j2_jesdn_HF_year,"dPhij1j2_jesdn_HF_year/F"); 
    newtree->Branch("dPhiHj1j2_jesdn_HF_year",&dPhiHj1j2_jesdn_HF_year,"dPhiHj1j2_jesdn_HF_year/F");
    newtree->Branch("pTj1_2p5_jesdn_HF_year",&pTj1_2p5_jesdn_HF_year,"pTj1_2p5_jesdn_HF_year/F"); 
    newtree->Branch("pt_leadingjet_pt30_eta2p5_jesdn_HF_year",&pt_leadingjet_pt30_eta2p5_jesdn_HF_year,"pt_leadingjet_pt30_eta2p5_2p5_jesdn_HF_year/F"); 
    newtree->Branch("yj1_2p5_jesdn_HF_year",&yj1_2p5_jesdn_HF_year,"yj1_2p5_jesdn_HF_year/F");
    newtree->Branch("pTj2_2p5_jesdn_HF_year",&pTj2_2p5_jesdn_HF_year,"pTj2_2p5_jesdn_HF_year/F"); 
    newtree->Branch("yj2_2p5_jesdn_HF_year",&yj2_2p5_jesdn_HF_year,"yj2_2p5_jesdn_HF_year/F");
    newtree->Branch("mass4lj_2p5_jesdn_HF_year",&mass4lj_2p5_jesdn_HF_year,"mass4lj_2p5_jesdn_HF_year/F");
    newtree->Branch("mass4ljj_2p5_jesdn_HF_year",&mass4ljj_2p5_jesdn_HF_year,"mass4ljj_2p5_jesdn_HF_year/F");
    newtree->Branch("pT4lj_2p5_jesdn_HF_year",&pT4lj_2p5_jesdn_HF_year,"pT4lj_2p5_jesdn_HF_year/F");
    newtree->Branch("pT4ljj_2p5_jesdn_HF_year",&pT4ljj_2p5_jesdn_HF_year,"pT4ljj_2p5_jesdn_HF_year/F");
    newtree->Branch("dPhiHj1_2p5_jesdn_HF_year",&dPhiHj1_2p5_jesdn_HF_year,"dPhiHj1_2p5_jesdn_HF_year/F"); 
    newtree->Branch("dyHj1_2p5_jesdn_HF_year",&dyHj1_2p5_jesdn_HF_year,"dyHj1_2p5_jesdn_HF_year/F");
    newtree->Branch("mj1j2_2p5_jesdn_HF_year",&mj1j2_2p5_jesdn_HF_year,"mj1j2_2p5_jesdn_HF_year/F"); 
    newtree->Branch("dEtaj1j2_2p5_jesdn_HF_year",&dEtaj1j2_2p5_jesdn_HF_year,"dEtaj1j2_2p5_jesdn_HF_year/F");
    newtree->Branch("dPhij1j2_2p5_jesdn_HF_year",&dPhij1j2_2p5_jesdn_HF_year,"dPhij1j2_2p5_jesdn_HF_year/F"); 
    newtree->Branch("dPhiHj1j2_2p5_jesdn_HF_year",&dPhiHj1j2_2p5_jesdn_HF_year,"dPhiHj1j2_2p5_jesdn_HF_year/F");
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// RelBal JES 
// jesup_RelBal
    newtree->Branch("njets_pt30_eta4p7_jesup_RelBal", &njets_pt30_eta4p7_jesup_RelBal, "njets_pt30_eta4p7_jesup_RelBal/F");
    newtree->Branch("TauC_Inc_0j_EnergyWgt_jesup_RelBal", &TauC_Inc_0j_EnergyWgt_jesup_RelBal, "TauC_Inc_0j_EnergyWgt_jesup_RelBal/F");
    newtree->Branch("TauB_Inc_0j_pTWgt_jesup_RelBal", &TauB_Inc_0j_pTWgt_jesup_RelBal, "TauB_Inc_0j_pTWgt_jesup_RelBal/F");
    newtree->Branch("njets_pt30_eta2p5_jesup_RelBal", &njets_pt30_eta2p5_jesup_RelBal, "njets_pt30_eta2p5_jesup_RelBal/F");
    newtree->Branch("pt_leadingjet_pt30_eta4p7_jesup_RelBal",&pt_leadingjet_pt30_eta4p7_jesup_RelBal,"pt_leadingjet_pt30_eta4p7_jesup_RelBal/F");
    newtree->Branch("pTj1_jesup_RelBal",&pTj1_jesup_RelBal,"pTj1_jesup_RelBal/F");
    newtree->Branch("etaj1_jesup_RelBal",&etaj1_jesup_RelBal,"etaj1_jesup_RelBal/F");
    newtree->Branch("pTj2_jesup_RelBal",&pTj2_jesup_RelBal,"pTj2_jesup_RelBal/F");
    newtree->Branch("etaj2_jesup_RelBal",&etaj2_jesup_RelBal,"etaj2_jesup_RelBal/F");
    newtree->Branch("yj1_jesup_RelBal",&yj1_jesup_RelBal,"yj1_jesup_RelBal/F");
    newtree->Branch("yj2_jesup_RelBal",&yj2_jesup_RelBal,"yj2_jesup_RelBal/F");
    newtree->Branch("dPhiHj1_jesup_RelBal",&dPhiHj1_jesup_RelBal,"dPhiHj1_jesup_RelBal/F"); 
    newtree->Branch("mass4lj_jesup_RelBal",&mass4lj_jesup_RelBal,"mass4lj_jesup_RelBal/F"); 
    newtree->Branch("mass4ljj_jesup_RelBal",&mass4ljj_jesup_RelBal,"mass4ljj_jesup_RelBal/F"); 
    newtree->Branch("pT4lj_jesup_RelBal",&pT4lj_jesup_RelBal,"pT4lj_jesup_RelBal/F"); 
    newtree->Branch("pT4ljj_jesup_RelBal",&pT4ljj_jesup_RelBal,"pT4ljj_jesup_RelBal/F"); 
    newtree->Branch("dyHj1_jesup_RelBal",&dyHj1_jesup_RelBal,"dyHj1_jesup_RelBal/F");
    newtree->Branch("mj1j2_jesup_RelBal",&mj1j2_jesup_RelBal,"mj1j2_jesup_RelBal/F"); 
    newtree->Branch("dEtaj1j2_jesup_RelBal",&dEtaj1j2_jesup_RelBal,"dEtaj1j2_jesup_RelBal/F");
    newtree->Branch("dPhij1j2_jesup_RelBal",&dPhij1j2_jesup_RelBal,"dPhij1j2_jesup_RelBal/F"); 
    newtree->Branch("dPhiHj1j2_jesup_RelBal",&dPhiHj1j2_jesup_RelBal,"dPhiHj1j2_jesup_RelBal/F");
    newtree->Branch("pTj1_2p5_jesup_RelBal",&pTj1_2p5_jesup_RelBal,"pTj1_2p5_jesup_RelBal/F"); 
    newtree->Branch("pt_leadingjet_pt30_eta2p5_jesup_RelBal",&pt_leadingjet_pt30_eta2p5_jesup_RelBal,"pt_leadingjet_pt30_eta2p5_jesup_RelBal/F");
    newtree->Branch("yj1_2p5_jesup_RelBal",&yj1_2p5_jesup_RelBal,"yj1_2p5_jesup_RelBal/F");
    newtree->Branch("pTj2_2p5_jesup_RelBal",&pTj2_2p5_jesup_RelBal,"pTj2_2p5_jesup_RelBal/F"); 
    newtree->Branch("yj2_2p5_jesup_RelBal",&yj2_2p5_jesup_RelBal,"yj2_2p5_jesup_RelBal/F");
    newtree->Branch("dPhiHj1_2p5_jesup_RelBal",&dPhiHj1_2p5_jesup_RelBal,"dPhiHj1_2p5_jesup_RelBal/F"); 
    newtree->Branch("mass4lj_2p5_jesup_RelBal",&mass4lj_2p5_jesup_RelBal,"mass4lj_2p5_jesup_RelBal/F");
    newtree->Branch("mass4ljj_2p5_jesup_RelBal",&mass4ljj_2p5_jesup_RelBal,"mass4ljj_2p5_jesup_RelBal/F");
    newtree->Branch("pT4lj_2p5_jesup_RelBal",&pT4lj_2p5_jesup_RelBal,"pT4lj_2p5_jesup_RelBal/F");
    newtree->Branch("pT4ljj_2p5_jesup_RelBal",&pT4ljj_2p5_jesup_RelBal,"pT4ljj_2p5_jesup_RelBal/F");
    newtree->Branch("dyHj1_2p5_jesup_RelBal",&dyHj1_2p5_jesup_RelBal,"dyHj1_2p5_jesup_RelBal/F");
    newtree->Branch("mj1j2_2p5_jesup_RelBal",&mj1j2_2p5_jesup_RelBal,"mj1j2_2p5_jesup_RelBal/F"); 
    newtree->Branch("dEtaj1j2_2p5_jesup_RelBal",&dEtaj1j2_2p5_jesup_RelBal,"dEtaj1j2_2p5_jesup_RelBal/F");
    newtree->Branch("dPhij1j2_2p5_jesup_RelBal",&dPhij1j2_2p5_jesup_RelBal,"dPhij1j2_2p5_jesup_RelBal/F"); 
    newtree->Branch("dPhiHj1j2_2p5_jesup_RelBal",&dPhiHj1j2_2p5_jesup_RelBal,"dPhiHj1j2_2p5_jesup_RelBal/F");
// jesdn_RelBal
    newtree->Branch("njets_pt30_eta4p7_jesdn_RelBal", &njets_pt30_eta4p7_jesdn_RelBal, "njets_pt30_eta4p7_jesdn_RelBal/F");
    newtree->Branch("TauC_Inc_0j_EnergyWgt_jesdn_RelBal", &TauC_Inc_0j_EnergyWgt_jesdn_RelBal, "TauC_Inc_0j_EnergyWgt_jesdn_RelBal/F");
    newtree->Branch("TauB_Inc_0j_pTWgt_jesdn_RelBal", &TauB_Inc_0j_pTWgt_jesdn_RelBal, "TauB_Inc_0j_pTWgt_jesdn_RelBal/F");
    newtree->Branch("njets_pt30_eta2p5_jesdn_RelBal", &njets_pt30_eta2p5_jesdn_RelBal, "njets_pt30_eta2p5_jesdn_RelBal/F");
    newtree->Branch("pt_leadingjet_pt30_eta4p7_jesdn_RelBal",&pt_leadingjet_pt30_eta4p7_jesdn_RelBal,"pt_leadingjet_pt30_eta4p7_jesdn_RelBal/F");
    newtree->Branch("pTj1_jesdn_RelBal",&pTj1_jesdn_RelBal,"pTj1_jesdn_RelBal/F");
    newtree->Branch("etaj1_jesdn_RelBal",&etaj1_jesdn_RelBal,"etaj1_jesdn_RelBal/F");
    newtree->Branch("pTj2_jesdn_RelBal",&pTj2_jesdn_RelBal,"pTj2_jesdn_RelBal/F");
    newtree->Branch("etaj2_jesdn_RelBal",&etaj2_jesdn_RelBal,"etaj2_jesdn_RelBal/F");
    newtree->Branch("yj1_jesdn_RelBal",&yj1_jesdn_RelBal,"yj1_jesdn_RelBal/F");
    newtree->Branch("yj2_jesdn_RelBal",&yj2_jesdn_RelBal,"yj2_jesdn_RelBal/F");
    newtree->Branch("dPhiHj1_jesdn_RelBal",&dPhiHj1_jesdn_RelBal,"dPhiHj1_jesdn_RelBal/F"); 
    newtree->Branch("dyHj1_jesdn_RelBal",&dyHj1_jesdn_RelBal,"dyHj1_jesdn_RelBal/F");
    newtree->Branch("mass4lj_jesdn_RelBal",&mass4lj_jesdn_RelBal,"mass4lj_jesdn_RelBal/F");
    newtree->Branch("mass4ljj_jesdn_RelBal",&mass4ljj_jesdn_RelBal,"mass4ljj_jesdn_RelBal/F");
    newtree->Branch("pT4lj_jesdn_RelBal",&pT4lj_jesdn_RelBal,"pT4lj_jesdn_RelBal/F");
    newtree->Branch("pT4ljj_jesdn_RelBal",&pT4ljj_jesdn_RelBal,"pT4ljj_jesdn_RelBal/F");
    newtree->Branch("mj1j2_jesdn_RelBal",&mj1j2_jesdn_RelBal,"mj1j2_jesdn_RelBal/F"); 
    newtree->Branch("dEtaj1j2_jesdn_RelBal",&dEtaj1j2_jesdn_RelBal,"dEtaj1j2_jesdn_RelBal/F");
    newtree->Branch("dPhij1j2_jesdn_RelBal",&dPhij1j2_jesdn_RelBal,"dPhij1j2_jesdn_RelBal/F"); 
    newtree->Branch("dPhiHj1j2_jesdn_RelBal",&dPhiHj1j2_jesdn_RelBal,"dPhiHj1j2_jesdn_RelBal/F");
    newtree->Branch("pTj1_2p5_jesdn_RelBal",&pTj1_2p5_jesdn_RelBal,"pTj1_2p5_jesdn_RelBal/F"); 
    newtree->Branch("pt_leadingjet_pt30_eta2p5_jesdn_RelBal",&pt_leadingjet_pt30_eta2p5_jesdn_RelBal,"pt_leadingjet_pt30_eta2p5_2p5_jesdn_RelBal/F"); 
    newtree->Branch("yj1_2p5_jesdn_RelBal",&yj1_2p5_jesdn_RelBal,"yj1_2p5_jesdn_RelBal/F");
    newtree->Branch("pTj2_2p5_jesdn_RelBal",&pTj2_2p5_jesdn_RelBal,"pTj2_2p5_jesdn_RelBal/F"); 
    newtree->Branch("yj2_2p5_jesdn_RelBal",&yj2_2p5_jesdn_RelBal,"yj2_2p5_jesdn_RelBal/F");
    newtree->Branch("mass4lj_2p5_jesdn_RelBal",&mass4lj_2p5_jesdn_RelBal,"mass4lj_2p5_jesdn_RelBal/F");
    newtree->Branch("mass4ljj_2p5_jesdn_RelBal",&mass4ljj_2p5_jesdn_RelBal,"mass4ljj_2p5_jesdn_RelBal/F");
    newtree->Branch("pT4lj_2p5_jesdn_RelBal",&pT4lj_2p5_jesdn_RelBal,"pT4lj_2p5_jesdn_RelBal/F");
    newtree->Branch("pT4ljj_2p5_jesdn_RelBal",&pT4ljj_2p5_jesdn_RelBal,"pT4ljj_2p5_jesdn_RelBal/F");
    newtree->Branch("dPhiHj1_2p5_jesdn_RelBal",&dPhiHj1_2p5_jesdn_RelBal,"dPhiHj1_2p5_jesdn_RelBal/F"); 
    newtree->Branch("dyHj1_2p5_jesdn_RelBal",&dyHj1_2p5_jesdn_RelBal,"dyHj1_2p5_jesdn_RelBal/F");
    newtree->Branch("mj1j2_2p5_jesdn_RelBal",&mj1j2_2p5_jesdn_RelBal,"mj1j2_2p5_jesdn_RelBal/F"); 
    newtree->Branch("dEtaj1j2_2p5_jesdn_RelBal",&dEtaj1j2_2p5_jesdn_RelBal,"dEtaj1j2_2p5_jesdn_RelBal/F");
    newtree->Branch("dPhij1j2_2p5_jesdn_RelBal",&dPhij1j2_2p5_jesdn_RelBal,"dPhij1j2_2p5_jesdn_RelBal/F"); 
    newtree->Branch("dPhiHj1j2_2p5_jesdn_RelBal",&dPhiHj1j2_2p5_jesdn_RelBal,"dPhiHj1j2_2p5_jesdn_RelBal/F");
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// RelSample_year JES 
// jesup_RelSample_year
    newtree->Branch("njets_pt30_eta4p7_jesup_RelSample_year", &njets_pt30_eta4p7_jesup_RelSample_year, "njets_pt30_eta4p7_jesup_RelSample_year/F");
    newtree->Branch("TauC_Inc_0j_EnergyWgt_jesup_RelSample_year", &TauC_Inc_0j_EnergyWgt_jesup_RelSample_year, "TauC_Inc_0j_EnergyWgt_jesup_RelSample_year/F");
    newtree->Branch("TauB_Inc_0j_pTWgt_jesup_RelSample_year", &TauB_Inc_0j_pTWgt_jesup_RelSample_year, "TauB_Inc_0j_pTWgt_jesup_RelSample_year/F");
    newtree->Branch("njets_pt30_eta2p5_jesup_RelSample_year", &njets_pt30_eta2p5_jesup_RelSample_year, "njets_pt30_eta2p5_jesup_RelSample_year/F");
    newtree->Branch("pt_leadingjet_pt30_eta4p7_jesup_RelSample_year",&pt_leadingjet_pt30_eta4p7_jesup_RelSample_year,"pt_leadingjet_pt30_eta4p7_jesup_RelSample_year/F");
    newtree->Branch("pTj1_jesup_RelSample_year",&pTj1_jesup_RelSample_year,"pTj1_jesup_RelSample_year/F");
    newtree->Branch("etaj1_jesup_RelSample_year",&etaj1_jesup_RelSample_year,"etaj1_jesup_RelSample_year/F");
    newtree->Branch("pTj2_jesup_RelSample_year",&pTj2_jesup_RelSample_year,"pTj2_jesup_RelSample_year/F");
    newtree->Branch("etaj2_jesup_RelSample_year",&etaj2_jesup_RelSample_year,"etaj2_jesup_RelSample_year/F");
    newtree->Branch("yj1_jesup_RelSample_year",&yj1_jesup_RelSample_year,"yj1_jesup_RelSample_year/F");
    newtree->Branch("yj2_jesup_RelSample_year",&yj2_jesup_RelSample_year,"yj2_jesup_RelSample_year/F");
    newtree->Branch("dPhiHj1_jesup_RelSample_year",&dPhiHj1_jesup_RelSample_year,"dPhiHj1_jesup_RelSample_year/F"); 
    newtree->Branch("mass4lj_jesup_RelSample_year",&mass4lj_jesup_RelSample_year,"mass4lj_jesup_RelSample_year/F"); 
    newtree->Branch("mass4ljj_jesup_RelSample_year",&mass4ljj_jesup_RelSample_year,"mass4ljj_jesup_RelSample_year/F"); 
    newtree->Branch("pT4lj_jesup_RelSample_year",&pT4lj_jesup_RelSample_year,"pT4lj_jesup_RelSample_year/F"); 
    newtree->Branch("pT4ljj_jesup_RelSample_year",&pT4ljj_jesup_RelSample_year,"pT4ljj_jesup_RelSample_year/F"); 
    newtree->Branch("dyHj1_jesup_RelSample_year",&dyHj1_jesup_RelSample_year,"dyHj1_jesup_RelSample_year/F");
    newtree->Branch("mj1j2_jesup_RelSample_year",&mj1j2_jesup_RelSample_year,"mj1j2_jesup_RelSample_year/F"); 
    newtree->Branch("dEtaj1j2_jesup_RelSample_year",&dEtaj1j2_jesup_RelSample_year,"dEtaj1j2_jesup_RelSample_year/F");
    newtree->Branch("dPhij1j2_jesup_RelSample_year",&dPhij1j2_jesup_RelSample_year,"dPhij1j2_jesup_RelSample_year/F"); 
    newtree->Branch("dPhiHj1j2_jesup_RelSample_year",&dPhiHj1j2_jesup_RelSample_year,"dPhiHj1j2_jesup_RelSample_year/F");
    newtree->Branch("pTj1_2p5_jesup_RelSample_year",&pTj1_2p5_jesup_RelSample_year,"pTj1_2p5_jesup_RelSample_year/F"); 
    newtree->Branch("pt_leadingjet_pt30_eta2p5_jesup_RelSample_year",&pt_leadingjet_pt30_eta2p5_jesup_RelSample_year,"pt_leadingjet_pt30_eta2p5_jesup_RelSample_year/F");
    newtree->Branch("yj1_2p5_jesup_RelSample_year",&yj1_2p5_jesup_RelSample_year,"yj1_2p5_jesup_RelSample_year/F");
    newtree->Branch("pTj2_2p5_jesup_RelSample_year",&pTj2_2p5_jesup_RelSample_year,"pTj2_2p5_jesup_RelSample_year/F"); 
    newtree->Branch("yj2_2p5_jesup_RelSample_year",&yj2_2p5_jesup_RelSample_year,"yj2_2p5_jesup_RelSample_year/F");
    newtree->Branch("dPhiHj1_2p5_jesup_RelSample_year",&dPhiHj1_2p5_jesup_RelSample_year,"dPhiHj1_2p5_jesup_RelSample_year/F"); 
    newtree->Branch("mass4lj_2p5_jesup_RelSample_year",&mass4lj_2p5_jesup_RelSample_year,"mass4lj_2p5_jesup_RelSample_year/F");
    newtree->Branch("mass4ljj_2p5_jesup_RelSample_year",&mass4ljj_2p5_jesup_RelSample_year,"mass4ljj_2p5_jesup_RelSample_year/F");
    newtree->Branch("pT4lj_2p5_jesup_RelSample_year",&pT4lj_2p5_jesup_RelSample_year,"pT4lj_2p5_jesup_RelSample_year/F");
    newtree->Branch("pT4ljj_2p5_jesup_RelSample_year",&pT4ljj_2p5_jesup_RelSample_year,"pT4ljj_2p5_jesup_RelSample_year/F");
    newtree->Branch("dyHj1_2p5_jesup_RelSample_year",&dyHj1_2p5_jesup_RelSample_year,"dyHj1_2p5_jesup_RelSample_year/F");
    newtree->Branch("mj1j2_2p5_jesup_RelSample_year",&mj1j2_2p5_jesup_RelSample_year,"mj1j2_2p5_jesup_RelSample_year/F"); 
    newtree->Branch("dEtaj1j2_2p5_jesup_RelSample_year",&dEtaj1j2_2p5_jesup_RelSample_year,"dEtaj1j2_2p5_jesup_RelSample_year/F");
    newtree->Branch("dPhij1j2_2p5_jesup_RelSample_year",&dPhij1j2_2p5_jesup_RelSample_year,"dPhij1j2_2p5_jesup_RelSample_year/F"); 
    newtree->Branch("dPhiHj1j2_2p5_jesup_RelSample_year",&dPhiHj1j2_2p5_jesup_RelSample_year,"dPhiHj1j2_2p5_jesup_RelSample_year/F");
// jesdn_RelSample_year
    newtree->Branch("njets_pt30_eta4p7_jesdn_RelSample_year", &njets_pt30_eta4p7_jesdn_RelSample_year, "njets_pt30_eta4p7_jesdn_RelSample_year/F");
    newtree->Branch("TauC_Inc_0j_EnergyWgt_jesdn_RelSample_year", &TauC_Inc_0j_EnergyWgt_jesdn_RelSample_year, "TauC_Inc_0j_EnergyWgt_jesdn_RelSample_year/F");
    newtree->Branch("TauB_Inc_0j_pTWgt_jesdn_RelSample_year", &TauB_Inc_0j_pTWgt_jesdn_RelSample_year, "TauB_Inc_0j_pTWgt_jesdn_RelSample_year/F");
    newtree->Branch("njets_pt30_eta2p5_jesdn_RelSample_year", &njets_pt30_eta2p5_jesdn_RelSample_year, "njets_pt30_eta2p5_jesdn_RelSample_year/F");
    newtree->Branch("pt_leadingjet_pt30_eta4p7_jesdn_RelSample_year",&pt_leadingjet_pt30_eta4p7_jesdn_RelSample_year,"pt_leadingjet_pt30_eta4p7_jesdn_RelSample_year/F");
    newtree->Branch("pTj1_jesdn_RelSample_year",&pTj1_jesdn_RelSample_year,"pTj1_jesdn_RelSample_year/F");
    newtree->Branch("etaj1_jesdn_RelSample_year",&etaj1_jesdn_RelSample_year,"etaj1_jesdn_RelSample_year/F");
    newtree->Branch("pTj2_jesdn_RelSample_year",&pTj2_jesdn_RelSample_year,"pTj2_jesdn_RelSample_year/F");
    newtree->Branch("etaj2_jesdn_RelSample_year",&etaj2_jesdn_RelSample_year,"etaj2_jesdn_RelSample_year/F");
    newtree->Branch("yj1_jesdn_RelSample_year",&yj1_jesdn_RelSample_year,"yj1_jesdn_RelSample_year/F");
    newtree->Branch("yj2_jesdn_RelSample_year",&yj2_jesdn_RelSample_year,"yj2_jesdn_RelSample_year/F");
    newtree->Branch("dPhiHj1_jesdn_RelSample_year",&dPhiHj1_jesdn_RelSample_year,"dPhiHj1_jesdn_RelSample_year/F"); 
    newtree->Branch("dyHj1_jesdn_RelSample_year",&dyHj1_jesdn_RelSample_year,"dyHj1_jesdn_RelSample_year/F");
    newtree->Branch("mass4lj_jesdn_RelSample_year",&mass4lj_jesdn_RelSample_year,"mass4lj_jesdn_RelSample_year/F");
    newtree->Branch("mass4ljj_jesdn_RelSample_year",&mass4ljj_jesdn_RelSample_year,"mass4ljj_jesdn_RelSample_year/F");
    newtree->Branch("pT4lj_jesdn_RelSample_year",&pT4lj_jesdn_RelSample_year,"pT4lj_jesdn_RelSample_year/F");
    newtree->Branch("pT4ljj_jesdn_RelSample_year",&pT4ljj_jesdn_RelSample_year,"pT4ljj_jesdn_RelSample_year/F");
    newtree->Branch("mj1j2_jesdn_RelSample_year",&mj1j2_jesdn_RelSample_year,"mj1j2_jesdn_RelSample_year/F"); 
    newtree->Branch("dEtaj1j2_jesdn_RelSample_year",&dEtaj1j2_jesdn_RelSample_year,"dEtaj1j2_jesdn_RelSample_year/F");
    newtree->Branch("dPhij1j2_jesdn_RelSample_year",&dPhij1j2_jesdn_RelSample_year,"dPhij1j2_jesdn_RelSample_year/F"); 
    newtree->Branch("dPhiHj1j2_jesdn_RelSample_year",&dPhiHj1j2_jesdn_RelSample_year,"dPhiHj1j2_jesdn_RelSample_year/F");
    newtree->Branch("pTj1_2p5_jesdn_RelSample_year",&pTj1_2p5_jesdn_RelSample_year,"pTj1_2p5_jesdn_RelSample_year/F"); 
    newtree->Branch("pt_leadingjet_pt30_eta2p5_jesdn_RelSample_year",&pt_leadingjet_pt30_eta2p5_jesdn_RelSample_year,"pt_leadingjet_pt30_eta2p5_2p5_jesdn_RelSample_year/F"); 
    newtree->Branch("yj1_2p5_jesdn_RelSample_year",&yj1_2p5_jesdn_RelSample_year,"yj1_2p5_jesdn_RelSample_year/F");
    newtree->Branch("pTj2_2p5_jesdn_RelSample_year",&pTj2_2p5_jesdn_RelSample_year,"pTj2_2p5_jesdn_RelSample_year/F"); 
    newtree->Branch("yj2_2p5_jesdn_RelSample_year",&yj2_2p5_jesdn_RelSample_year,"yj2_2p5_jesdn_RelSample_year/F");
    newtree->Branch("mass4lj_2p5_jesdn_RelSample_year",&mass4lj_2p5_jesdn_RelSample_year,"mass4lj_2p5_jesdn_RelSample_year/F");
    newtree->Branch("mass4ljj_2p5_jesdn_RelSample_year",&mass4ljj_2p5_jesdn_RelSample_year,"mass4ljj_2p5_jesdn_RelSample_year/F");
    newtree->Branch("pT4lj_2p5_jesdn_RelSample_year",&pT4lj_2p5_jesdn_RelSample_year,"pT4lj_2p5_jesdn_RelSample_year/F");
    newtree->Branch("pT4ljj_2p5_jesdn_RelSample_year",&pT4ljj_2p5_jesdn_RelSample_year,"pT4ljj_2p5_jesdn_RelSample_year/F");
    newtree->Branch("dPhiHj1_2p5_jesdn_RelSample_year",&dPhiHj1_2p5_jesdn_RelSample_year,"dPhiHj1_2p5_jesdn_RelSample_year/F"); 
    newtree->Branch("dyHj1_2p5_jesdn_RelSample_year",&dyHj1_2p5_jesdn_RelSample_year,"dyHj1_2p5_jesdn_RelSample_year/F");
    newtree->Branch("mj1j2_2p5_jesdn_RelSample_year",&mj1j2_2p5_jesdn_RelSample_year,"mj1j2_2p5_jesdn_RelSample_year/F"); 
    newtree->Branch("dEtaj1j2_2p5_jesdn_RelSample_year",&dEtaj1j2_2p5_jesdn_RelSample_year,"dEtaj1j2_2p5_jesdn_RelSample_year/F");
    newtree->Branch("dPhij1j2_2p5_jesdn_RelSample_year",&dPhij1j2_2p5_jesdn_RelSample_year,"dPhij1j2_2p5_jesdn_RelSample_year/F"); 
    newtree->Branch("dPhiHj1j2_2p5_jesdn_RelSample_year",&dPhiHj1j2_2p5_jesdn_RelSample_year,"dPhiHj1j2_2p5_jesdn_RelSample_year/F");
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Total JES 
// jesup_Total
    newtree->Branch("njets_pt30_eta4p7_jesup_Total", &njets_pt30_eta4p7_jesup_Total, "njets_pt30_eta4p7_jesup_Total/F");
    newtree->Branch("TauC_Inc_0j_EnergyWgt_jesup_Total", &TauC_Inc_0j_EnergyWgt_jesup_Total, "TauC_Inc_0j_EnergyWgt_jesup_Total/F");
    newtree->Branch("TauB_Inc_0j_pTWgt_jesup_Total", &TauB_Inc_0j_pTWgt_jesup_Total, "TauB_Inc_0j_pTWgt_jesup_Total/F");
    newtree->Branch("njets_pt30_eta2p5_jesup_Total", &njets_pt30_eta2p5_jesup_Total, "njets_pt30_eta2p5_jesup_Total/F");
    newtree->Branch("pt_leadingjet_pt30_eta4p7_jesup_Total",&pt_leadingjet_pt30_eta4p7_jesup_Total,"pt_leadingjet_pt30_eta4p7_jesup_Total/F");
    newtree->Branch("pTj1_jesup_Total",&pTj1_jesup_Total,"pTj1_jesup_Total/F");
    newtree->Branch("etaj1_jesup_Total",&etaj1_jesup_Total,"etaj1_jesup_Total/F");
    newtree->Branch("pTj2_jesup_Total",&pTj2_jesup_Total,"pTj2_jesup_Total/F");
    newtree->Branch("etaj2_jesup_Total",&etaj2_jesup_Total,"etaj2_jesup_Total/F");
    newtree->Branch("yj1_jesup_Total",&yj1_jesup_Total,"yj1_jesup_Total/F");
    newtree->Branch("yj2_jesup_Total",&yj2_jesup_Total,"yj2_jesup_Total/F");
    newtree->Branch("dPhiHj1_jesup_Total",&dPhiHj1_jesup_Total,"dPhiHj1_jesup_Total/F"); 
    newtree->Branch("mass4lj_jesup_Total",&mass4lj_jesup_Total,"mass4lj_jesup_Total/F"); 
    newtree->Branch("mass4ljj_jesup_Total",&mass4ljj_jesup_Total,"mass4ljj_jesup_Total/F"); 
    newtree->Branch("pT4lj_jesup_Total",&pT4lj_jesup_Total,"pT4lj_jesup_Total/F"); 
    newtree->Branch("pT4ljj_jesup_Total",&pT4ljj_jesup_Total,"pT4ljj_jesup_Total/F"); 
    newtree->Branch("dyHj1_jesup_Total",&dyHj1_jesup_Total,"dyHj1_jesup_Total/F");
    newtree->Branch("mj1j2_jesup_Total",&mj1j2_jesup_Total,"mj1j2_jesup_Total/F"); 
    newtree->Branch("dEtaj1j2_jesup_Total",&dEtaj1j2_jesup_Total,"dEtaj1j2_jesup_Total/F");
    newtree->Branch("dPhij1j2_jesup_Total",&dPhij1j2_jesup_Total,"dPhij1j2_jesup_Total/F"); 
    newtree->Branch("dPhiHj1j2_jesup_Total",&dPhiHj1j2_jesup_Total,"dPhiHj1j2_jesup_Total/F");
    newtree->Branch("pTj1_2p5_jesup_Total",&pTj1_2p5_jesup_Total,"pTj1_2p5_jesup_Total/F"); 
    newtree->Branch("pt_leadingjet_pt30_eta2p5_jesup_Total",&pt_leadingjet_pt30_eta2p5_jesup_Total,"pt_leadingjet_pt30_eta2p5_jesup_Total/F");
    newtree->Branch("yj1_2p5_jesup_Total",&yj1_2p5_jesup_Total,"yj1_2p5_jesup_Total/F");
    newtree->Branch("pTj2_2p5_jesup_Total",&pTj2_2p5_jesup_Total,"pTj2_2p5_jesup_Total/F"); 
    newtree->Branch("yj2_2p5_jesup_Total",&yj2_2p5_jesup_Total,"yj2_2p5_jesup_Total/F");
    newtree->Branch("dPhiHj1_2p5_jesup_Total",&dPhiHj1_2p5_jesup_Total,"dPhiHj1_2p5_jesup_Total/F"); 
    newtree->Branch("mass4lj_2p5_jesup_Total",&mass4lj_2p5_jesup_Total,"mass4lj_2p5_jesup_Total/F");
    newtree->Branch("mass4ljj_2p5_jesup_Total",&mass4ljj_2p5_jesup_Total,"mass4ljj_2p5_jesup_Total/F");
    newtree->Branch("pT4lj_2p5_jesup_Total",&pT4lj_2p5_jesup_Total,"pT4lj_2p5_jesup_Total/F");
    newtree->Branch("pT4ljj_2p5_jesup_Total",&pT4ljj_2p5_jesup_Total,"pT4ljj_2p5_jesup_Total/F");
    newtree->Branch("dyHj1_2p5_jesup_Total",&dyHj1_2p5_jesup_Total,"dyHj1_2p5_jesup_Total/F");
    newtree->Branch("mj1j2_2p5_jesup_Total",&mj1j2_2p5_jesup_Total,"mj1j2_2p5_jesup_Total/F"); 
    newtree->Branch("dEtaj1j2_2p5_jesup_Total",&dEtaj1j2_2p5_jesup_Total,"dEtaj1j2_2p5_jesup_Total/F");
    newtree->Branch("dPhij1j2_2p5_jesup_Total",&dPhij1j2_2p5_jesup_Total,"dPhij1j2_2p5_jesup_Total/F"); 
    newtree->Branch("dPhiHj1j2_2p5_jesup_Total",&dPhiHj1j2_2p5_jesup_Total,"dPhiHj1j2_2p5_jesup_Total/F");
// jesdn_Total
    newtree->Branch("njets_pt30_eta4p7_jesdn_Total", &njets_pt30_eta4p7_jesdn_Total, "njets_pt30_eta4p7_jesdn_Total/F");
    newtree->Branch("TauC_Inc_0j_EnergyWgt_jesdn_Total", &TauC_Inc_0j_EnergyWgt_jesdn_Total, "TauC_Inc_0j_EnergyWgt_jesdn_Total/F");
    newtree->Branch("TauB_Inc_0j_pTWgt_jesdn_Total", &TauB_Inc_0j_pTWgt_jesdn_Total, "TauB_Inc_0j_pTWgt_jesdn_Total/F");
    newtree->Branch("njets_pt30_eta2p5_jesdn_Total", &njets_pt30_eta2p5_jesdn_Total, "njets_pt30_eta2p5_jesdn_Total/F");
    newtree->Branch("pt_leadingjet_pt30_eta4p7_jesdn_Total",&pt_leadingjet_pt30_eta4p7_jesdn_Total,"pt_leadingjet_pt30_eta4p7_jesdn_Total/F");
    newtree->Branch("pTj1_jesdn_Total",&pTj1_jesdn_Total,"pTj1_jesdn_Total/F");
    newtree->Branch("etaj1_jesdn_Total",&etaj1_jesdn_Total,"etaj1_jesdn_Total/F");
    newtree->Branch("pTj2_jesdn_Total",&pTj2_jesdn_Total,"pTj2_jesdn_Total/F");
    newtree->Branch("etaj2_jesdn_Total",&etaj2_jesdn_Total,"etaj2_jesdn_Total/F");
    newtree->Branch("yj1_jesdn_Total",&yj1_jesdn_Total,"yj1_jesdn_Total/F");
    newtree->Branch("yj2_jesdn_Total",&yj2_jesdn_Total,"yj2_jesdn_Total/F");
    newtree->Branch("dPhiHj1_jesdn_Total",&dPhiHj1_jesdn_Total,"dPhiHj1_jesdn_Total/F"); 
    newtree->Branch("dyHj1_jesdn_Total",&dyHj1_jesdn_Total,"dyHj1_jesdn_Total/F");
    newtree->Branch("mass4lj_jesdn_Total",&mass4lj_jesdn_Total,"mass4lj_jesdn_Total/F");
    newtree->Branch("mass4ljj_jesdn_Total",&mass4ljj_jesdn_Total,"mass4ljj_jesdn_Total/F");
    newtree->Branch("pT4lj_jesdn_Total",&pT4lj_jesdn_Total,"pT4lj_jesdn_Total/F");
    newtree->Branch("pT4ljj_jesdn_Total",&pT4ljj_jesdn_Total,"pT4ljj_jesdn_Total/F");
    newtree->Branch("mj1j2_jesdn_Total",&mj1j2_jesdn_Total,"mj1j2_jesdn_Total/F"); 
    newtree->Branch("dEtaj1j2_jesdn_Total",&dEtaj1j2_jesdn_Total,"dEtaj1j2_jesdn_Total/F");
    newtree->Branch("dPhij1j2_jesdn_Total",&dPhij1j2_jesdn_Total,"dPhij1j2_jesdn_Total/F"); 
    newtree->Branch("dPhiHj1j2_jesdn_Total",&dPhiHj1j2_jesdn_Total,"dPhiHj1j2_jesdn_Total/F");
    newtree->Branch("pTj1_2p5_jesdn_Total",&pTj1_2p5_jesdn_Total,"pTj1_2p5_jesdn_Total/F"); 
    newtree->Branch("pt_leadingjet_pt30_eta2p5_jesdn_Total",&pt_leadingjet_pt30_eta2p5_jesdn_Total,"pt_leadingjet_pt30_eta2p5_2p5_jesdn_Total/F"); 
    newtree->Branch("yj1_2p5_jesdn_Total",&yj1_2p5_jesdn_Total,"yj1_2p5_jesdn_Total/F");
    newtree->Branch("pTj2_2p5_jesdn_Total",&pTj2_2p5_jesdn_Total,"pTj2_2p5_jesdn_Total/F"); 
    newtree->Branch("yj2_2p5_jesdn_Total",&yj2_2p5_jesdn_Total,"yj2_2p5_jesdn_Total/F");
    newtree->Branch("mass4lj_2p5_jesdn_Total",&mass4lj_2p5_jesdn_Total,"mass4lj_2p5_jesdn_Total/F");
    newtree->Branch("mass4ljj_2p5_jesdn_Total",&mass4ljj_2p5_jesdn_Total,"mass4ljj_2p5_jesdn_Total/F");
    newtree->Branch("pT4lj_2p5_jesdn_Total",&pT4lj_2p5_jesdn_Total,"pT4lj_2p5_jesdn_Total/F");
    newtree->Branch("pT4ljj_2p5_jesdn_Total",&pT4ljj_2p5_jesdn_Total,"pT4ljj_2p5_jesdn_Total/F");
    newtree->Branch("dPhiHj1_2p5_jesdn_Total",&dPhiHj1_2p5_jesdn_Total,"dPhiHj1_2p5_jesdn_Total/F"); 
    newtree->Branch("dyHj1_2p5_jesdn_Total",&dyHj1_2p5_jesdn_Total,"dyHj1_2p5_jesdn_Total/F");
    newtree->Branch("mj1j2_2p5_jesdn_Total",&mj1j2_2p5_jesdn_Total,"mj1j2_2p5_jesdn_Total/F"); 
    newtree->Branch("dEtaj1j2_2p5_jesdn_Total",&dEtaj1j2_2p5_jesdn_Total,"dEtaj1j2_2p5_jesdn_Total/F");
    newtree->Branch("dPhij1j2_2p5_jesdn_Total",&dPhij1j2_2p5_jesdn_Total,"dPhij1j2_2p5_jesdn_Total/F"); 
    newtree->Branch("dPhiHj1j2_2p5_jesdn_Total",&dPhiHj1j2_2p5_jesdn_Total,"dPhiHj1j2_2p5_jesdn_Total/F");



    std::set<TString> runlumieventSet;
    for (Long64_t i=0;i<nentries; i++) {
    //for (Long64_t i=0;i<10000; i++) {  // temp

	// do a Skim
	if ( doM4lSkim && (GENmass4l<GENmass4l_lo || GENmass4l>GENmass4l_hi) && !(passedFullSelection==1 && mass4l>mass4l_lo && mass4l<mass4l_hi) ) continue;

        lep_dataMC_new.clear();
        lep_dataMCErr_new.clear();
/*	jes_unc_split.clear();
	pt_jesup_split.clear();
	pt_jesdn_split.clear();
*/
	isH4l = false; 
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

// validation branches
        pTj1=-9999.0; 
        TauC_Inc_0j_EnergyWgt_nom=-9999.0, TauB_Inc_0j_pTWgt_nom=-9999.0;	
// initialize Abs variables
        jet1index_jesup_Abs=-1, jet2index_jesup_Abs=-1, jet1index2p5_jesup_Abs=-1, jet2index2p5_jesup_Abs=-1, jet1pt_jesup_Abs=0.0, jet2pt_jesup_Abs=0.0, jet1pt2p5_jesup_Abs=0.0, jet2pt2p5_jesup_Abs=0.0,jet1index_jesdn_Abs=-1, jet2index_jesdn_Abs=-1,jet1index2p5_jesdn_Abs=-1, jet2index2p5_jesdn_Abs=-1, jet1pt_jesdn_Abs=0.0, jet2pt_jesdn_Abs=0.0,jet1pt2p5_jesdn_Abs=0.0, jet2pt2p5_jesdn_Abs=0.0;
	TauC_Inc_0j_EnergyWgt_jesup_Abs=-9999.0, TauB_Inc_0j_pTWgt_jesup_Abs=-9999.0, TauC_Inc_0j_EnergyWgt_jesdn_Abs=-9999.0, TauB_Inc_0j_pTWgt_jesdn_Abs=-9999.0, njets_pt30_eta4p7_jesup_Abs=0, njets_pt30_eta2p5_jesup_Abs=0, njets_pt30_eta4p7_jesdn_Abs=0, njets_pt30_eta2p5_jesdn_Abs=0, pTj1_jesup_Abs=-9999.0, pt_leadingjet_pt30_eta4p7_jesup_Abs=-9999.0, pTj1_2p5_jesup_Abs=-9999.0, pt_leadingjet_pt30_eta2p5_jesup_Abs=-9999.0, pTj2_jesup_Abs=-9999.0, dEtaj1j2_jesup_Abs=-9999.0, yj1_jesup_Abs=-9999.0, yj2_jesup_Abs=-9999.0, dPhiHj1_jesup_Abs=-9999.0, dyHj1_jesup_Abs=-9999.0, mass4lj_jesup_Abs=-9999.0, pT4lj_jesup_Abs=-9999.0, mass4ljj_jesup_Abs=-9999.0, pT4ljj_jesup_Abs=-9999.0, mass4lj_2p5_jesup_Abs=-9999.0, pT4lj_2p5_jesup_Abs=-9999.0, mass4ljj_2p5_jesup_Abs=-9999.0, pT4ljj_2p5_jesup_Abs=-9999.0, mj1j2_jesup_Abs=-9999.0, dEtaj1j2_jesup_Abs=-9999.0, dPhij1j2_jesup_Abs=-9999.0, dPhiHj1j2_jesup_Abs=-9999.0, yj1_2p5_jesup_Abs=-9999.0, yj2_2p5_jesup_Abs=-9999.0, dPhiHj1_2p5_jesup_Abs=-9999.0, dyHj1_2p5_jesup_Abs=-9999.0, mj1j2_2p5_jesup_Abs=-9999.0, dEtaj1j2_2p5_jesup_Abs=-9999.0, dPhij1j2_2p5_jesup_Abs=-9999.0, dPhiHj1j2_2p5_jesup_Abs=-9999.0, pTj1_jesdn_Abs=-9999.0, pTj2_jesdn_Abs=-9999.0, pTj1_2p5_jesdn_Abs=-9999.0, pt_leadingjet_pt30_eta2p5_jesdn_Abs=-9999.0, pTj2_jesdn_Abs=-9999.0, mj1j2_jesdn_Abs=-9999.0, dEtaj1j2_jesdn_Abs=-9999.0, yj1_jesdn_Abs=-9999.0, yj2_jesdn_Abs=-9999.0, dPhiHj1_jesdn_Abs=-9999.0, dyHj1_jesdn_Abs=-9999.0, mj1j2_jesdn_Abs=-9999.0, dEtaj1j2_jesdn_Abs=-9999.0, dPhij1j2_jesdn_Abs=-9999.0, dPhiHj1j2_jesdn_Abs=-9999.0, yj1_2p5_jesdn_Abs=-9999.0, yj2_2p5_jesdn_Abs=-9999.0, dPhiHj1_2p5_jesdn_Abs=-9999.0, dyHj1_2p5_jesdn_Abs=-9999.0, mj1j2_2p5_jesdn_Abs=-9999.0, dEtaj1j2_2p5_jesdn_Abs=-9999.0, dPhij1j2_2p5_jesdn_Abs=-9999.0, dPhiHj1j2_2p5_jesdn_Abs=-9999.0, mass4lj_jesdn_Abs=-9999.0, pT4lj_jesdn_Abs=-9999.0, mass4ljj_jesdn_Abs=-9999.0, pT4ljj_jesdn_Abs=-9999.0, mass4lj_2p5_jesdn_Abs=-9999.0, pT4lj_2p5_jesdn_Abs=-9999.0, mass4ljj_2p5_jesdn_Abs=-9999.0, pT4ljj_2p5_jesdn_Abs=-9999.0;  // FIXME 


// initialize Abs_year variables
        jet1index_jesup_Abs_year=-1, jet2index_jesup_Abs_year=-1, jet1index2p5_jesup_Abs_year=-1, jet2index2p5_jesup_Abs_year=-1, jet1pt_jesup_Abs_year=0.0, jet2pt_jesup_Abs_year=0.0, jet1pt2p5_jesup_Abs_year=0.0, jet2pt2p5_jesup_Abs_year=0.0,jet1index_jesdn_Abs_year=-1, jet2index_jesdn_Abs_year=-1,jet1index2p5_jesdn_Abs_year=-1, jet2index2p5_jesdn_Abs_year=-1, jet1pt_jesdn_Abs_year=0.0, jet2pt_jesdn_Abs_year=0.0,jet1pt2p5_jesdn_Abs_year=0.0, jet2pt2p5_jesdn_Abs_year=0.0;
	TauC_Inc_0j_EnergyWgt_jesup_Abs_year=-9999.0, TauB_Inc_0j_pTWgt_jesup_Abs_year=-9999.0, TauC_Inc_0j_EnergyWgt_jesdn_Abs_year=-9999.0, TauB_Inc_0j_pTWgt_jesdn_Abs_year=-9999.0, njets_pt30_eta4p7_jesup_Abs_year=0, njets_pt30_eta2p5_jesup_Abs_year=0, njets_pt30_eta4p7_jesdn_Abs_year=0, njets_pt30_eta2p5_jesdn_Abs_year=0, pTj1_jesup_Abs_year=-9999.0, pt_leadingjet_pt30_eta4p7_jesup_Abs_year=-9999.0, pTj1_2p5_jesup_Abs_year=-9999.0, pt_leadingjet_pt30_eta2p5_jesup_Abs_year=-9999.0, pTj2_jesup_Abs_year=-9999.0, dEtaj1j2_jesup_Abs_year=-9999.0, yj1_jesup_Abs_year=-9999.0, yj2_jesup_Abs_year=-9999.0, dPhiHj1_jesup_Abs_year=-9999.0, dyHj1_jesup_Abs_year=-9999.0, mass4lj_jesup_Abs_year=-9999.0, pT4lj_jesup_Abs_year=-9999.0, mass4ljj_jesup_Abs_year=-9999.0, pT4ljj_jesup_Abs_year=-9999.0, mass4lj_2p5_jesup_Abs_year=-9999.0, pT4lj_2p5_jesup_Abs_year=-9999.0, mass4ljj_2p5_jesup_Abs_year=-9999.0, pT4ljj_2p5_jesup_Abs_year=-9999.0, mj1j2_jesup_Abs_year=-9999.0, dEtaj1j2_jesup_Abs_year=-9999.0, dPhij1j2_jesup_Abs_year=-9999.0, dPhiHj1j2_jesup_Abs_year=-9999.0, yj1_2p5_jesup_Abs_year=-9999.0, yj2_2p5_jesup_Abs_year=-9999.0, dPhiHj1_2p5_jesup_Abs_year=-9999.0, dyHj1_2p5_jesup_Abs_year=-9999.0, mj1j2_2p5_jesup_Abs_year=-9999.0, dEtaj1j2_2p5_jesup_Abs_year=-9999.0, dPhij1j2_2p5_jesup_Abs_year=-9999.0, dPhiHj1j2_2p5_jesup_Abs_year=-9999.0, pTj1_jesdn_Abs_year=-9999.0, pTj2_jesdn_Abs_year=-9999.0, pTj1_2p5_jesdn_Abs_year=-9999.0, pt_leadingjet_pt30_eta2p5_jesdn_Abs_year=-9999.0, pTj2_jesdn_Abs_year=-9999.0, mj1j2_jesdn_Abs_year=-9999.0, dEtaj1j2_jesdn_Abs_year=-9999.0, yj1_jesdn_Abs_year=-9999.0, yj2_jesdn_Abs_year=-9999.0, dPhiHj1_jesdn_Abs_year=-9999.0, dyHj1_jesdn_Abs_year=-9999.0, mj1j2_jesdn_Abs_year=-9999.0, dEtaj1j2_jesdn_Abs_year=-9999.0, dPhij1j2_jesdn_Abs_year=-9999.0, dPhiHj1j2_jesdn_Abs_year=-9999.0, yj1_2p5_jesdn_Abs_year=-9999.0, yj2_2p5_jesdn_Abs_year=-9999.0, dPhiHj1_2p5_jesdn_Abs_year=-9999.0, dyHj1_2p5_jesdn_Abs_year=-9999.0, mj1j2_2p5_jesdn_Abs_year=-9999.0, dEtaj1j2_2p5_jesdn_Abs_year=-9999.0, dPhij1j2_2p5_jesdn_Abs_year=-9999.0, dPhiHj1j2_2p5_jesdn_Abs_year=-9999.0, mass4lj_jesdn_Abs_year=-9999.0, pT4lj_jesdn_Abs_year=-9999.0, mass4ljj_jesdn_Abs_year=-9999.0, pT4ljj_jesdn_Abs_year=-9999.0, mass4lj_2p5_jesdn_Abs_year=-9999.0, pT4lj_2p5_jesdn_Abs_year=-9999.0, mass4ljj_2p5_jesdn_Abs_year=-9999.0, pT4ljj_2p5_jesdn_Abs_year=-9999.0;  // FIXME 


// initialize Total variables
        jet1index_jesup_Total=-1, jet2index_jesup_Total=-1, jet1index2p5_jesup_Total=-1, jet2index2p5_jesup_Total=-1, jet1pt_jesup_Total=0.0, jet2pt_jesup_Total=0.0, jet1pt2p5_jesup_Total=0.0, jet2pt2p5_jesup_Total=0.0,jet1index_jesdn_Total=-1, jet2index_jesdn_Total=-1,jet1index2p5_jesdn_Total=-1, jet2index2p5_jesdn_Total=-1, jet1pt_jesdn_Total=0.0, jet2pt_jesdn_Total=0.0,jet1pt2p5_jesdn_Total=0.0, jet2pt2p5_jesdn_Total=0.0;
	TauC_Inc_0j_EnergyWgt_jesup_Total=-9999.0, TauB_Inc_0j_pTWgt_jesup_Total=-9999.0, TauC_Inc_0j_EnergyWgt_jesdn_Total=-9999.0, TauB_Inc_0j_pTWgt_jesdn_Total=-9999.0, njets_pt30_eta4p7_jesup_Total=0, njets_pt30_eta2p5_jesup_Total=0, njets_pt30_eta4p7_jesdn_Total=0, njets_pt30_eta2p5_jesdn_Total=0, pTj1_jesup_Total=-9999.0, pt_leadingjet_pt30_eta4p7_jesup_Total=-9999.0, pTj1_2p5_jesup_Total=-9999.0, pt_leadingjet_pt30_eta2p5_jesup_Total=-9999.0, pTj2_jesup_Total=-9999.0, dEtaj1j2_jesup_Total=-9999.0, yj1_jesup_Total=-9999.0, yj2_jesup_Total=-9999.0, dPhiHj1_jesup_Total=-9999.0, dyHj1_jesup_Total=-9999.0, mass4lj_jesup_Total=-9999.0, pT4lj_jesup_Total=-9999.0, mass4ljj_jesup_Total=-9999.0, pT4ljj_jesup_Total=-9999.0, mass4lj_2p5_jesup_Total=-9999.0, pT4lj_2p5_jesup_Total=-9999.0, mass4ljj_2p5_jesup_Total=-9999.0, pT4ljj_2p5_jesup_Total=-9999.0, mj1j2_jesup_Total=-9999.0, dEtaj1j2_jesup_Total=-9999.0, dPhij1j2_jesup_Total=-9999.0, dPhiHj1j2_jesup_Total=-9999.0, yj1_2p5_jesup_Total=-9999.0, yj2_2p5_jesup_Total=-9999.0, dPhiHj1_2p5_jesup_Total=-9999.0, dyHj1_2p5_jesup_Total=-9999.0, mj1j2_2p5_jesup_Total=-9999.0, dEtaj1j2_2p5_jesup_Total=-9999.0, dPhij1j2_2p5_jesup_Total=-9999.0, dPhiHj1j2_2p5_jesup_Total=-9999.0, pTj1_jesdn_Total=-9999.0, pTj2_jesdn_Total=-9999.0, pTj1_2p5_jesdn_Total=-9999.0, pt_leadingjet_pt30_eta2p5_jesdn_Total=-9999.0, pTj2_jesdn_Total=-9999.0, mj1j2_jesdn_Total=-9999.0, dEtaj1j2_jesdn_Total=-9999.0, yj1_jesdn_Total=-9999.0, yj2_jesdn_Total=-9999.0, dPhiHj1_jesdn_Total=-9999.0, dyHj1_jesdn_Total=-9999.0, mj1j2_jesdn_Total=-9999.0, dEtaj1j2_jesdn_Total=-9999.0, dPhij1j2_jesdn_Total=-9999.0, dPhiHj1j2_jesdn_Total=-9999.0, yj1_2p5_jesdn_Total=-9999.0, yj2_2p5_jesdn_Total=-9999.0, dPhiHj1_2p5_jesdn_Total=-9999.0, dyHj1_2p5_jesdn_Total=-9999.0, mj1j2_2p5_jesdn_Total=-9999.0, dEtaj1j2_2p5_jesdn_Total=-9999.0, dPhij1j2_2p5_jesdn_Total=-9999.0, dPhiHj1j2_2p5_jesdn_Total=-9999.0, mass4lj_jesdn_Total=-9999.0, pT4lj_jesdn_Total=-9999.0, mass4ljj_jesdn_Total=-9999.0, pT4ljj_jesdn_Total=-9999.0, mass4lj_2p5_jesdn_Total=-9999.0, pT4lj_2p5_jesdn_Total=-9999.0, mass4ljj_2p5_jesdn_Total=-9999.0, pT4ljj_2p5_jesdn_Total=-9999.0;  // FIXME 


// initialize BBEC1 variables
        jet1index_jesup_BBEC1=-1, jet2index_jesup_BBEC1=-1, jet1index2p5_jesup_BBEC1=-1, jet2index2p5_jesup_BBEC1=-1, jet1pt_jesup_BBEC1=0.0, jet2pt_jesup_BBEC1=0.0, jet1pt2p5_jesup_BBEC1=0.0, jet2pt2p5_jesup_BBEC1=0.0,jet1index_jesdn_BBEC1=-1, jet2index_jesdn_BBEC1=-1,jet1index2p5_jesdn_BBEC1=-1, jet2index2p5_jesdn_BBEC1=-1, jet1pt_jesdn_BBEC1=0.0, jet2pt_jesdn_BBEC1=0.0,jet1pt2p5_jesdn_BBEC1=0.0, jet2pt2p5_jesdn_BBEC1=0.0;
	TauC_Inc_0j_EnergyWgt_jesup_BBEC1=-9999.0, TauB_Inc_0j_pTWgt_jesup_BBEC1=-9999.0, TauC_Inc_0j_EnergyWgt_jesdn_BBEC1=-9999.0, TauB_Inc_0j_pTWgt_jesdn_BBEC1=-9999.0, njets_pt30_eta4p7_jesup_BBEC1=0, njets_pt30_eta2p5_jesup_BBEC1=0, njets_pt30_eta4p7_jesdn_BBEC1=0, njets_pt30_eta2p5_jesdn_BBEC1=0, pTj1_jesup_BBEC1=-9999.0, pt_leadingjet_pt30_eta4p7_jesup_BBEC1=-9999.0, pTj1_2p5_jesup_BBEC1=-9999.0, pt_leadingjet_pt30_eta2p5_jesup_BBEC1=-9999.0, pTj2_jesup_BBEC1=-9999.0, dEtaj1j2_jesup_BBEC1=-9999.0, yj1_jesup_BBEC1=-9999.0, yj2_jesup_BBEC1=-9999.0, dPhiHj1_jesup_BBEC1=-9999.0, dyHj1_jesup_BBEC1=-9999.0, mass4lj_jesup_BBEC1=-9999.0, pT4lj_jesup_BBEC1=-9999.0, mass4ljj_jesup_BBEC1=-9999.0, pT4ljj_jesup_BBEC1=-9999.0, mass4lj_2p5_jesup_BBEC1=-9999.0, pT4lj_2p5_jesup_BBEC1=-9999.0, mass4ljj_2p5_jesup_BBEC1=-9999.0, pT4ljj_2p5_jesup_BBEC1=-9999.0, mj1j2_jesup_BBEC1=-9999.0, dEtaj1j2_jesup_BBEC1=-9999.0, dPhij1j2_jesup_BBEC1=-9999.0, dPhiHj1j2_jesup_BBEC1=-9999.0, yj1_2p5_jesup_BBEC1=-9999.0, yj2_2p5_jesup_BBEC1=-9999.0, dPhiHj1_2p5_jesup_BBEC1=-9999.0, dyHj1_2p5_jesup_BBEC1=-9999.0, mj1j2_2p5_jesup_BBEC1=-9999.0, dEtaj1j2_2p5_jesup_BBEC1=-9999.0, dPhij1j2_2p5_jesup_BBEC1=-9999.0, dPhiHj1j2_2p5_jesup_BBEC1=-9999.0, pTj1_jesdn_BBEC1=-9999.0, pTj2_jesdn_BBEC1=-9999.0, pTj1_2p5_jesdn_BBEC1=-9999.0, pt_leadingjet_pt30_eta2p5_jesdn_BBEC1=-9999.0, pTj2_jesdn_BBEC1=-9999.0, mj1j2_jesdn_BBEC1=-9999.0, dEtaj1j2_jesdn_BBEC1=-9999.0, yj1_jesdn_BBEC1=-9999.0, yj2_jesdn_BBEC1=-9999.0, dPhiHj1_jesdn_BBEC1=-9999.0, dyHj1_jesdn_BBEC1=-9999.0, mj1j2_jesdn_BBEC1=-9999.0, dEtaj1j2_jesdn_BBEC1=-9999.0, dPhij1j2_jesdn_BBEC1=-9999.0, dPhiHj1j2_jesdn_BBEC1=-9999.0, yj1_2p5_jesdn_BBEC1=-9999.0, yj2_2p5_jesdn_BBEC1=-9999.0, dPhiHj1_2p5_jesdn_BBEC1=-9999.0, dyHj1_2p5_jesdn_BBEC1=-9999.0, mj1j2_2p5_jesdn_BBEC1=-9999.0, dEtaj1j2_2p5_jesdn_BBEC1=-9999.0, dPhij1j2_2p5_jesdn_BBEC1=-9999.0, dPhiHj1j2_2p5_jesdn_BBEC1=-9999.0, mass4lj_jesdn_BBEC1=-9999.0, pT4lj_jesdn_BBEC1=-9999.0, mass4ljj_jesdn_BBEC1=-9999.0, pT4ljj_jesdn_BBEC1=-9999.0, mass4lj_2p5_jesdn_BBEC1=-9999.0, pT4lj_2p5_jesdn_BBEC1=-9999.0, mass4ljj_2p5_jesdn_BBEC1=-9999.0, pT4ljj_2p5_jesdn_BBEC1=-9999.0;  // FIXME 


// initialize BBEC1_year variables
        jet1index_jesup_BBEC1_year=-1, jet2index_jesup_BBEC1_year=-1, jet1index2p5_jesup_BBEC1_year=-1, jet2index2p5_jesup_BBEC1_year=-1, jet1pt_jesup_BBEC1_year=0.0, jet2pt_jesup_BBEC1_year=0.0, jet1pt2p5_jesup_BBEC1_year=0.0, jet2pt2p5_jesup_BBEC1_year=0.0,jet1index_jesdn_BBEC1_year=-1, jet2index_jesdn_BBEC1_year=-1,jet1index2p5_jesdn_BBEC1_year=-1, jet2index2p5_jesdn_BBEC1_year=-1, jet1pt_jesdn_BBEC1_year=0.0, jet2pt_jesdn_BBEC1_year=0.0,jet1pt2p5_jesdn_BBEC1_year=0.0, jet2pt2p5_jesdn_BBEC1_year=0.0;
	TauC_Inc_0j_EnergyWgt_jesup_BBEC1_year=-9999.0, TauB_Inc_0j_pTWgt_jesup_BBEC1_year=-9999.0, TauC_Inc_0j_EnergyWgt_jesdn_BBEC1_year=-9999.0, TauB_Inc_0j_pTWgt_jesdn_BBEC1_year=-9999.0, njets_pt30_eta4p7_jesup_BBEC1_year=0, njets_pt30_eta2p5_jesup_BBEC1_year=0, njets_pt30_eta4p7_jesdn_BBEC1_year=0, njets_pt30_eta2p5_jesdn_BBEC1_year=0, pTj1_jesup_BBEC1_year=-9999.0, pt_leadingjet_pt30_eta4p7_jesup_BBEC1_year=-9999.0, pTj1_2p5_jesup_BBEC1_year=-9999.0, pt_leadingjet_pt30_eta2p5_jesup_BBEC1_year=-9999.0, pTj2_jesup_BBEC1_year=-9999.0, dEtaj1j2_jesup_BBEC1_year=-9999.0, yj1_jesup_BBEC1_year=-9999.0, yj2_jesup_BBEC1_year=-9999.0, dPhiHj1_jesup_BBEC1_year=-9999.0, dyHj1_jesup_BBEC1_year=-9999.0, mass4lj_jesup_BBEC1_year=-9999.0, pT4lj_jesup_BBEC1_year=-9999.0, mass4ljj_jesup_BBEC1_year=-9999.0, pT4ljj_jesup_BBEC1_year=-9999.0, mass4lj_2p5_jesup_BBEC1_year=-9999.0, pT4lj_2p5_jesup_BBEC1_year=-9999.0, mass4ljj_2p5_jesup_BBEC1_year=-9999.0, pT4ljj_2p5_jesup_BBEC1_year=-9999.0, mj1j2_jesup_BBEC1_year=-9999.0, dEtaj1j2_jesup_BBEC1_year=-9999.0, dPhij1j2_jesup_BBEC1_year=-9999.0, dPhiHj1j2_jesup_BBEC1_year=-9999.0, yj1_2p5_jesup_BBEC1_year=-9999.0, yj2_2p5_jesup_BBEC1_year=-9999.0, dPhiHj1_2p5_jesup_BBEC1_year=-9999.0, dyHj1_2p5_jesup_BBEC1_year=-9999.0, mj1j2_2p5_jesup_BBEC1_year=-9999.0, dEtaj1j2_2p5_jesup_BBEC1_year=-9999.0, dPhij1j2_2p5_jesup_BBEC1_year=-9999.0, dPhiHj1j2_2p5_jesup_BBEC1_year=-9999.0, pTj1_jesdn_BBEC1_year=-9999.0, pTj2_jesdn_BBEC1_year=-9999.0, pTj1_2p5_jesdn_BBEC1_year=-9999.0, pt_leadingjet_pt30_eta2p5_jesdn_BBEC1_year=-9999.0, pTj2_jesdn_BBEC1_year=-9999.0, mj1j2_jesdn_BBEC1_year=-9999.0, dEtaj1j2_jesdn_BBEC1_year=-9999.0, yj1_jesdn_BBEC1_year=-9999.0, yj2_jesdn_BBEC1_year=-9999.0, dPhiHj1_jesdn_BBEC1_year=-9999.0, dyHj1_jesdn_BBEC1_year=-9999.0, mj1j2_jesdn_BBEC1_year=-9999.0, dEtaj1j2_jesdn_BBEC1_year=-9999.0, dPhij1j2_jesdn_BBEC1_year=-9999.0, dPhiHj1j2_jesdn_BBEC1_year=-9999.0, yj1_2p5_jesdn_BBEC1_year=-9999.0, yj2_2p5_jesdn_BBEC1_year=-9999.0, dPhiHj1_2p5_jesdn_BBEC1_year=-9999.0, dyHj1_2p5_jesdn_BBEC1_year=-9999.0, mj1j2_2p5_jesdn_BBEC1_year=-9999.0, dEtaj1j2_2p5_jesdn_BBEC1_year=-9999.0, dPhij1j2_2p5_jesdn_BBEC1_year=-9999.0, dPhiHj1j2_2p5_jesdn_BBEC1_year=-9999.0, mass4lj_jesdn_BBEC1_year=-9999.0, pT4lj_jesdn_BBEC1_year=-9999.0, mass4ljj_jesdn_BBEC1_year=-9999.0, pT4ljj_jesdn_BBEC1_year=-9999.0, mass4lj_2p5_jesdn_BBEC1_year=-9999.0, pT4lj_2p5_jesdn_BBEC1_year=-9999.0, mass4ljj_2p5_jesdn_BBEC1_year=-9999.0, pT4ljj_2p5_jesdn_BBEC1_year=-9999.0;  // FIXME 


// initialize EC2 variables
        jet1index_jesup_EC2=-1, jet2index_jesup_EC2=-1, jet1index2p5_jesup_EC2=-1, jet2index2p5_jesup_EC2=-1, jet1pt_jesup_EC2=0.0, jet2pt_jesup_EC2=0.0, jet1pt2p5_jesup_EC2=0.0, jet2pt2p5_jesup_EC2=0.0,jet1index_jesdn_EC2=-1, jet2index_jesdn_EC2=-1,jet1index2p5_jesdn_EC2=-1, jet2index2p5_jesdn_EC2=-1, jet1pt_jesdn_EC2=0.0, jet2pt_jesdn_EC2=0.0,jet1pt2p5_jesdn_EC2=0.0, jet2pt2p5_jesdn_EC2=0.0;
	TauC_Inc_0j_EnergyWgt_jesup_EC2=-9999.0, TauB_Inc_0j_pTWgt_jesup_EC2=-9999.0, TauC_Inc_0j_EnergyWgt_jesdn_EC2=-9999.0, TauB_Inc_0j_pTWgt_jesdn_EC2=-9999.0, njets_pt30_eta4p7_jesup_EC2=0, njets_pt30_eta2p5_jesup_EC2=0, njets_pt30_eta4p7_jesdn_EC2=0, njets_pt30_eta2p5_jesdn_EC2=0, pTj1_jesup_EC2=-9999.0, pt_leadingjet_pt30_eta4p7_jesup_EC2=-9999.0, pTj1_2p5_jesup_EC2=-9999.0, pt_leadingjet_pt30_eta2p5_jesup_EC2=-9999.0, pTj2_jesup_EC2=-9999.0, dEtaj1j2_jesup_EC2=-9999.0, yj1_jesup_EC2=-9999.0, yj2_jesup_EC2=-9999.0, dPhiHj1_jesup_EC2=-9999.0, dyHj1_jesup_EC2=-9999.0, mass4lj_jesup_EC2=-9999.0, pT4lj_jesup_EC2=-9999.0, mass4ljj_jesup_EC2=-9999.0, pT4ljj_jesup_EC2=-9999.0, mass4lj_2p5_jesup_EC2=-9999.0, pT4lj_2p5_jesup_EC2=-9999.0, mass4ljj_2p5_jesup_EC2=-9999.0, pT4ljj_2p5_jesup_EC2=-9999.0, mj1j2_jesup_EC2=-9999.0, dEtaj1j2_jesup_EC2=-9999.0, dPhij1j2_jesup_EC2=-9999.0, dPhiHj1j2_jesup_EC2=-9999.0, yj1_2p5_jesup_EC2=-9999.0, yj2_2p5_jesup_EC2=-9999.0, dPhiHj1_2p5_jesup_EC2=-9999.0, dyHj1_2p5_jesup_EC2=-9999.0, mj1j2_2p5_jesup_EC2=-9999.0, dEtaj1j2_2p5_jesup_EC2=-9999.0, dPhij1j2_2p5_jesup_EC2=-9999.0, dPhiHj1j2_2p5_jesup_EC2=-9999.0, pTj1_jesdn_EC2=-9999.0, pTj2_jesdn_EC2=-9999.0, pTj1_2p5_jesdn_EC2=-9999.0, pt_leadingjet_pt30_eta2p5_jesdn_EC2=-9999.0, pTj2_jesdn_EC2=-9999.0, mj1j2_jesdn_EC2=-9999.0, dEtaj1j2_jesdn_EC2=-9999.0, yj1_jesdn_EC2=-9999.0, yj2_jesdn_EC2=-9999.0, dPhiHj1_jesdn_EC2=-9999.0, dyHj1_jesdn_EC2=-9999.0, mj1j2_jesdn_EC2=-9999.0, dEtaj1j2_jesdn_EC2=-9999.0, dPhij1j2_jesdn_EC2=-9999.0, dPhiHj1j2_jesdn_EC2=-9999.0, yj1_2p5_jesdn_EC2=-9999.0, yj2_2p5_jesdn_EC2=-9999.0, dPhiHj1_2p5_jesdn_EC2=-9999.0, dyHj1_2p5_jesdn_EC2=-9999.0, mj1j2_2p5_jesdn_EC2=-9999.0, dEtaj1j2_2p5_jesdn_EC2=-9999.0, dPhij1j2_2p5_jesdn_EC2=-9999.0, dPhiHj1j2_2p5_jesdn_EC2=-9999.0, mass4lj_jesdn_EC2=-9999.0, pT4lj_jesdn_EC2=-9999.0, mass4ljj_jesdn_EC2=-9999.0, pT4ljj_jesdn_EC2=-9999.0, mass4lj_2p5_jesdn_EC2=-9999.0, pT4lj_2p5_jesdn_EC2=-9999.0, mass4ljj_2p5_jesdn_EC2=-9999.0, pT4ljj_2p5_jesdn_EC2=-9999.0;  // FIXME 


// initialize EC2_year variables
        jet1index_jesup_EC2_year=-1, jet2index_jesup_EC2_year=-1, jet1index2p5_jesup_EC2_year=-1, jet2index2p5_jesup_EC2_year=-1, jet1pt_jesup_EC2_year=0.0, jet2pt_jesup_EC2_year=0.0, jet1pt2p5_jesup_EC2_year=0.0, jet2pt2p5_jesup_EC2_year=0.0,jet1index_jesdn_EC2_year=-1, jet2index_jesdn_EC2_year=-1,jet1index2p5_jesdn_EC2_year=-1, jet2index2p5_jesdn_EC2_year=-1, jet1pt_jesdn_EC2_year=0.0, jet2pt_jesdn_EC2_year=0.0,jet1pt2p5_jesdn_EC2_year=0.0, jet2pt2p5_jesdn_EC2_year=0.0;
	TauC_Inc_0j_EnergyWgt_jesup_EC2_year=-9999.0, TauB_Inc_0j_pTWgt_jesup_EC2_year=-9999.0, TauC_Inc_0j_EnergyWgt_jesdn_EC2_year=-9999.0, TauB_Inc_0j_pTWgt_jesdn_EC2_year=-9999.0, njets_pt30_eta4p7_jesup_EC2_year=0, njets_pt30_eta2p5_jesup_EC2_year=0, njets_pt30_eta4p7_jesdn_EC2_year=0, njets_pt30_eta2p5_jesdn_EC2_year=0, pTj1_jesup_EC2_year=-9999.0, pt_leadingjet_pt30_eta4p7_jesup_EC2_year=-9999.0, pTj1_2p5_jesup_EC2_year=-9999.0, pt_leadingjet_pt30_eta2p5_jesup_EC2_year=-9999.0, pTj2_jesup_EC2_year=-9999.0, dEtaj1j2_jesup_EC2_year=-9999.0, yj1_jesup_EC2_year=-9999.0, yj2_jesup_EC2_year=-9999.0, dPhiHj1_jesup_EC2_year=-9999.0, dyHj1_jesup_EC2_year=-9999.0, mass4lj_jesup_EC2_year=-9999.0, pT4lj_jesup_EC2_year=-9999.0, mass4ljj_jesup_EC2_year=-9999.0, pT4ljj_jesup_EC2_year=-9999.0, mass4lj_2p5_jesup_EC2_year=-9999.0, pT4lj_2p5_jesup_EC2_year=-9999.0, mass4ljj_2p5_jesup_EC2_year=-9999.0, pT4ljj_2p5_jesup_EC2_year=-9999.0, mj1j2_jesup_EC2_year=-9999.0, dEtaj1j2_jesup_EC2_year=-9999.0, dPhij1j2_jesup_EC2_year=-9999.0, dPhiHj1j2_jesup_EC2_year=-9999.0, yj1_2p5_jesup_EC2_year=-9999.0, yj2_2p5_jesup_EC2_year=-9999.0, dPhiHj1_2p5_jesup_EC2_year=-9999.0, dyHj1_2p5_jesup_EC2_year=-9999.0, mj1j2_2p5_jesup_EC2_year=-9999.0, dEtaj1j2_2p5_jesup_EC2_year=-9999.0, dPhij1j2_2p5_jesup_EC2_year=-9999.0, dPhiHj1j2_2p5_jesup_EC2_year=-9999.0, pTj1_jesdn_EC2_year=-9999.0, pTj2_jesdn_EC2_year=-9999.0, pTj1_2p5_jesdn_EC2_year=-9999.0, pt_leadingjet_pt30_eta2p5_jesdn_EC2_year=-9999.0, pTj2_jesdn_EC2_year=-9999.0, mj1j2_jesdn_EC2_year=-9999.0, dEtaj1j2_jesdn_EC2_year=-9999.0, yj1_jesdn_EC2_year=-9999.0, yj2_jesdn_EC2_year=-9999.0, dPhiHj1_jesdn_EC2_year=-9999.0, dyHj1_jesdn_EC2_year=-9999.0, mj1j2_jesdn_EC2_year=-9999.0, dEtaj1j2_jesdn_EC2_year=-9999.0, dPhij1j2_jesdn_EC2_year=-9999.0, dPhiHj1j2_jesdn_EC2_year=-9999.0, yj1_2p5_jesdn_EC2_year=-9999.0, yj2_2p5_jesdn_EC2_year=-9999.0, dPhiHj1_2p5_jesdn_EC2_year=-9999.0, dyHj1_2p5_jesdn_EC2_year=-9999.0, mj1j2_2p5_jesdn_EC2_year=-9999.0, dEtaj1j2_2p5_jesdn_EC2_year=-9999.0, dPhij1j2_2p5_jesdn_EC2_year=-9999.0, dPhiHj1j2_2p5_jesdn_EC2_year=-9999.0, mass4lj_jesdn_EC2_year=-9999.0, pT4lj_jesdn_EC2_year=-9999.0, mass4ljj_jesdn_EC2_year=-9999.0, pT4ljj_jesdn_EC2_year=-9999.0, mass4lj_2p5_jesdn_EC2_year=-9999.0, pT4lj_2p5_jesdn_EC2_year=-9999.0, mass4ljj_2p5_jesdn_EC2_year=-9999.0, pT4ljj_2p5_jesdn_EC2_year=-9999.0;  // FIXME 


// initialize FlavQCD variables
        jet1index_jesup_FlavQCD=-1, jet2index_jesup_FlavQCD=-1, jet1index2p5_jesup_FlavQCD=-1, jet2index2p5_jesup_FlavQCD=-1, jet1pt_jesup_FlavQCD=0.0, jet2pt_jesup_FlavQCD=0.0, jet1pt2p5_jesup_FlavQCD=0.0, jet2pt2p5_jesup_FlavQCD=0.0,jet1index_jesdn_FlavQCD=-1, jet2index_jesdn_FlavQCD=-1,jet1index2p5_jesdn_FlavQCD=-1, jet2index2p5_jesdn_FlavQCD=-1, jet1pt_jesdn_FlavQCD=0.0, jet2pt_jesdn_FlavQCD=0.0,jet1pt2p5_jesdn_FlavQCD=0.0, jet2pt2p5_jesdn_FlavQCD=0.0;
	TauC_Inc_0j_EnergyWgt_jesup_FlavQCD=-9999.0, TauB_Inc_0j_pTWgt_jesup_FlavQCD=-9999.0, TauC_Inc_0j_EnergyWgt_jesdn_FlavQCD=-9999.0, TauB_Inc_0j_pTWgt_jesdn_FlavQCD=-9999.0, njets_pt30_eta4p7_jesup_FlavQCD=0, njets_pt30_eta2p5_jesup_FlavQCD=0, njets_pt30_eta4p7_jesdn_FlavQCD=0, njets_pt30_eta2p5_jesdn_FlavQCD=0, pTj1_jesup_FlavQCD=-9999.0, pt_leadingjet_pt30_eta4p7_jesup_FlavQCD=-9999.0, pTj1_2p5_jesup_FlavQCD=-9999.0, pt_leadingjet_pt30_eta2p5_jesup_FlavQCD=-9999.0, pTj2_jesup_FlavQCD=-9999.0, dEtaj1j2_jesup_FlavQCD=-9999.0, yj1_jesup_FlavQCD=-9999.0, yj2_jesup_FlavQCD=-9999.0, dPhiHj1_jesup_FlavQCD=-9999.0, dyHj1_jesup_FlavQCD=-9999.0, mass4lj_jesup_FlavQCD=-9999.0, pT4lj_jesup_FlavQCD=-9999.0, mass4ljj_jesup_FlavQCD=-9999.0, pT4ljj_jesup_FlavQCD=-9999.0, mass4lj_2p5_jesup_FlavQCD=-9999.0, pT4lj_2p5_jesup_FlavQCD=-9999.0, mass4ljj_2p5_jesup_FlavQCD=-9999.0, pT4ljj_2p5_jesup_FlavQCD=-9999.0, mj1j2_jesup_FlavQCD=-9999.0, dEtaj1j2_jesup_FlavQCD=-9999.0, dPhij1j2_jesup_FlavQCD=-9999.0, dPhiHj1j2_jesup_FlavQCD=-9999.0, yj1_2p5_jesup_FlavQCD=-9999.0, yj2_2p5_jesup_FlavQCD=-9999.0, dPhiHj1_2p5_jesup_FlavQCD=-9999.0, dyHj1_2p5_jesup_FlavQCD=-9999.0, mj1j2_2p5_jesup_FlavQCD=-9999.0, dEtaj1j2_2p5_jesup_FlavQCD=-9999.0, dPhij1j2_2p5_jesup_FlavQCD=-9999.0, dPhiHj1j2_2p5_jesup_FlavQCD=-9999.0, pTj1_jesdn_FlavQCD=-9999.0, pTj2_jesdn_FlavQCD=-9999.0, pTj1_2p5_jesdn_FlavQCD=-9999.0, pt_leadingjet_pt30_eta2p5_jesdn_FlavQCD=-9999.0, pTj2_jesdn_FlavQCD=-9999.0, mj1j2_jesdn_FlavQCD=-9999.0, dEtaj1j2_jesdn_FlavQCD=-9999.0, yj1_jesdn_FlavQCD=-9999.0, yj2_jesdn_FlavQCD=-9999.0, dPhiHj1_jesdn_FlavQCD=-9999.0, dyHj1_jesdn_FlavQCD=-9999.0, mj1j2_jesdn_FlavQCD=-9999.0, dEtaj1j2_jesdn_FlavQCD=-9999.0, dPhij1j2_jesdn_FlavQCD=-9999.0, dPhiHj1j2_jesdn_FlavQCD=-9999.0, yj1_2p5_jesdn_FlavQCD=-9999.0, yj2_2p5_jesdn_FlavQCD=-9999.0, dPhiHj1_2p5_jesdn_FlavQCD=-9999.0, dyHj1_2p5_jesdn_FlavQCD=-9999.0, mj1j2_2p5_jesdn_FlavQCD=-9999.0, dEtaj1j2_2p5_jesdn_FlavQCD=-9999.0, dPhij1j2_2p5_jesdn_FlavQCD=-9999.0, dPhiHj1j2_2p5_jesdn_FlavQCD=-9999.0, mass4lj_jesdn_FlavQCD=-9999.0, pT4lj_jesdn_FlavQCD=-9999.0, mass4ljj_jesdn_FlavQCD=-9999.0, pT4ljj_jesdn_FlavQCD=-9999.0, mass4lj_2p5_jesdn_FlavQCD=-9999.0, pT4lj_2p5_jesdn_FlavQCD=-9999.0, mass4ljj_2p5_jesdn_FlavQCD=-9999.0, pT4ljj_2p5_jesdn_FlavQCD=-9999.0;  // FIXME 


// initialize HF variables
        jet1index_jesup_HF=-1, jet2index_jesup_HF=-1, jet1index2p5_jesup_HF=-1, jet2index2p5_jesup_HF=-1, jet1pt_jesup_HF=0.0, jet2pt_jesup_HF=0.0, jet1pt2p5_jesup_HF=0.0, jet2pt2p5_jesup_HF=0.0,jet1index_jesdn_HF=-1, jet2index_jesdn_HF=-1,jet1index2p5_jesdn_HF=-1, jet2index2p5_jesdn_HF=-1, jet1pt_jesdn_HF=0.0, jet2pt_jesdn_HF=0.0,jet1pt2p5_jesdn_HF=0.0, jet2pt2p5_jesdn_HF=0.0;
	TauC_Inc_0j_EnergyWgt_jesup_HF=-9999.0, TauB_Inc_0j_pTWgt_jesup_HF=-9999.0, TauC_Inc_0j_EnergyWgt_jesdn_HF=-9999.0, TauB_Inc_0j_pTWgt_jesdn_HF=-9999.0, njets_pt30_eta4p7_jesup_HF=0, njets_pt30_eta2p5_jesup_HF=0, njets_pt30_eta4p7_jesdn_HF=0, njets_pt30_eta2p5_jesdn_HF=0, pTj1_jesup_HF=-9999.0, pt_leadingjet_pt30_eta4p7_jesup_HF=-9999.0, pTj1_2p5_jesup_HF=-9999.0, pt_leadingjet_pt30_eta2p5_jesup_HF=-9999.0, pTj2_jesup_HF=-9999.0, dEtaj1j2_jesup_HF=-9999.0, yj1_jesup_HF=-9999.0, yj2_jesup_HF=-9999.0, dPhiHj1_jesup_HF=-9999.0, dyHj1_jesup_HF=-9999.0, mass4lj_jesup_HF=-9999.0, pT4lj_jesup_HF=-9999.0, mass4ljj_jesup_HF=-9999.0, pT4ljj_jesup_HF=-9999.0, mass4lj_2p5_jesup_HF=-9999.0, pT4lj_2p5_jesup_HF=-9999.0, mass4ljj_2p5_jesup_HF=-9999.0, pT4ljj_2p5_jesup_HF=-9999.0, mj1j2_jesup_HF=-9999.0, dEtaj1j2_jesup_HF=-9999.0, dPhij1j2_jesup_HF=-9999.0, dPhiHj1j2_jesup_HF=-9999.0, yj1_2p5_jesup_HF=-9999.0, yj2_2p5_jesup_HF=-9999.0, dPhiHj1_2p5_jesup_HF=-9999.0, dyHj1_2p5_jesup_HF=-9999.0, mj1j2_2p5_jesup_HF=-9999.0, dEtaj1j2_2p5_jesup_HF=-9999.0, dPhij1j2_2p5_jesup_HF=-9999.0, dPhiHj1j2_2p5_jesup_HF=-9999.0, pTj1_jesdn_HF=-9999.0, pTj2_jesdn_HF=-9999.0, pTj1_2p5_jesdn_HF=-9999.0, pt_leadingjet_pt30_eta2p5_jesdn_HF=-9999.0, pTj2_jesdn_HF=-9999.0, mj1j2_jesdn_HF=-9999.0, dEtaj1j2_jesdn_HF=-9999.0, yj1_jesdn_HF=-9999.0, yj2_jesdn_HF=-9999.0, dPhiHj1_jesdn_HF=-9999.0, dyHj1_jesdn_HF=-9999.0, mj1j2_jesdn_HF=-9999.0, dEtaj1j2_jesdn_HF=-9999.0, dPhij1j2_jesdn_HF=-9999.0, dPhiHj1j2_jesdn_HF=-9999.0, yj1_2p5_jesdn_HF=-9999.0, yj2_2p5_jesdn_HF=-9999.0, dPhiHj1_2p5_jesdn_HF=-9999.0, dyHj1_2p5_jesdn_HF=-9999.0, mj1j2_2p5_jesdn_HF=-9999.0, dEtaj1j2_2p5_jesdn_HF=-9999.0, dPhij1j2_2p5_jesdn_HF=-9999.0, dPhiHj1j2_2p5_jesdn_HF=-9999.0, mass4lj_jesdn_HF=-9999.0, pT4lj_jesdn_HF=-9999.0, mass4ljj_jesdn_HF=-9999.0, pT4ljj_jesdn_HF=-9999.0, mass4lj_2p5_jesdn_HF=-9999.0, pT4lj_2p5_jesdn_HF=-9999.0, mass4ljj_2p5_jesdn_HF=-9999.0, pT4ljj_2p5_jesdn_HF=-9999.0;  // FIXME 


// initialize HF_year variables
        jet1index_jesup_HF_year=-1, jet2index_jesup_HF_year=-1, jet1index2p5_jesup_HF_year=-1, jet2index2p5_jesup_HF_year=-1, jet1pt_jesup_HF_year=0.0, jet2pt_jesup_HF_year=0.0, jet1pt2p5_jesup_HF_year=0.0, jet2pt2p5_jesup_HF_year=0.0,jet1index_jesdn_HF_year=-1, jet2index_jesdn_HF_year=-1,jet1index2p5_jesdn_HF_year=-1, jet2index2p5_jesdn_HF_year=-1, jet1pt_jesdn_HF_year=0.0, jet2pt_jesdn_HF_year=0.0,jet1pt2p5_jesdn_HF_year=0.0, jet2pt2p5_jesdn_HF_year=0.0;
	TauC_Inc_0j_EnergyWgt_jesup_HF_year=-9999.0, TauB_Inc_0j_pTWgt_jesup_HF_year=-9999.0, TauC_Inc_0j_EnergyWgt_jesdn_HF_year=-9999.0, TauB_Inc_0j_pTWgt_jesdn_HF_year=-9999.0, njets_pt30_eta4p7_jesup_HF_year=0, njets_pt30_eta2p5_jesup_HF_year=0, njets_pt30_eta4p7_jesdn_HF_year=0, njets_pt30_eta2p5_jesdn_HF_year=0, pTj1_jesup_HF_year=-9999.0, pt_leadingjet_pt30_eta4p7_jesup_HF_year=-9999.0, pTj1_2p5_jesup_HF_year=-9999.0, pt_leadingjet_pt30_eta2p5_jesup_HF_year=-9999.0, pTj2_jesup_HF_year=-9999.0, dEtaj1j2_jesup_HF_year=-9999.0, yj1_jesup_HF_year=-9999.0, yj2_jesup_HF_year=-9999.0, dPhiHj1_jesup_HF_year=-9999.0, dyHj1_jesup_HF_year=-9999.0, mass4lj_jesup_HF_year=-9999.0, pT4lj_jesup_HF_year=-9999.0, mass4ljj_jesup_HF_year=-9999.0, pT4ljj_jesup_HF_year=-9999.0, mass4lj_2p5_jesup_HF_year=-9999.0, pT4lj_2p5_jesup_HF_year=-9999.0, mass4ljj_2p5_jesup_HF_year=-9999.0, pT4ljj_2p5_jesup_HF_year=-9999.0, mj1j2_jesup_HF_year=-9999.0, dEtaj1j2_jesup_HF_year=-9999.0, dPhij1j2_jesup_HF_year=-9999.0, dPhiHj1j2_jesup_HF_year=-9999.0, yj1_2p5_jesup_HF_year=-9999.0, yj2_2p5_jesup_HF_year=-9999.0, dPhiHj1_2p5_jesup_HF_year=-9999.0, dyHj1_2p5_jesup_HF_year=-9999.0, mj1j2_2p5_jesup_HF_year=-9999.0, dEtaj1j2_2p5_jesup_HF_year=-9999.0, dPhij1j2_2p5_jesup_HF_year=-9999.0, dPhiHj1j2_2p5_jesup_HF_year=-9999.0, pTj1_jesdn_HF_year=-9999.0, pTj2_jesdn_HF_year=-9999.0, pTj1_2p5_jesdn_HF_year=-9999.0, pt_leadingjet_pt30_eta2p5_jesdn_HF_year=-9999.0, pTj2_jesdn_HF_year=-9999.0, mj1j2_jesdn_HF_year=-9999.0, dEtaj1j2_jesdn_HF_year=-9999.0, yj1_jesdn_HF_year=-9999.0, yj2_jesdn_HF_year=-9999.0, dPhiHj1_jesdn_HF_year=-9999.0, dyHj1_jesdn_HF_year=-9999.0, mj1j2_jesdn_HF_year=-9999.0, dEtaj1j2_jesdn_HF_year=-9999.0, dPhij1j2_jesdn_HF_year=-9999.0, dPhiHj1j2_jesdn_HF_year=-9999.0, yj1_2p5_jesdn_HF_year=-9999.0, yj2_2p5_jesdn_HF_year=-9999.0, dPhiHj1_2p5_jesdn_HF_year=-9999.0, dyHj1_2p5_jesdn_HF_year=-9999.0, mj1j2_2p5_jesdn_HF_year=-9999.0, dEtaj1j2_2p5_jesdn_HF_year=-9999.0, dPhij1j2_2p5_jesdn_HF_year=-9999.0, dPhiHj1j2_2p5_jesdn_HF_year=-9999.0, mass4lj_jesdn_HF_year=-9999.0, pT4lj_jesdn_HF_year=-9999.0, mass4ljj_jesdn_HF_year=-9999.0, pT4ljj_jesdn_HF_year=-9999.0, mass4lj_2p5_jesdn_HF_year=-9999.0, pT4lj_2p5_jesdn_HF_year=-9999.0, mass4ljj_2p5_jesdn_HF_year=-9999.0, pT4ljj_2p5_jesdn_HF_year=-9999.0;  // FIXME 


// initialize RelBal variables
        jet1index_jesup_RelBal=-1, jet2index_jesup_RelBal=-1, jet1index2p5_jesup_RelBal=-1, jet2index2p5_jesup_RelBal=-1, jet1pt_jesup_RelBal=0.0, jet2pt_jesup_RelBal=0.0, jet1pt2p5_jesup_RelBal=0.0, jet2pt2p5_jesup_RelBal=0.0,jet1index_jesdn_RelBal=-1, jet2index_jesdn_RelBal=-1,jet1index2p5_jesdn_RelBal=-1, jet2index2p5_jesdn_RelBal=-1, jet1pt_jesdn_RelBal=0.0, jet2pt_jesdn_RelBal=0.0,jet1pt2p5_jesdn_RelBal=0.0, jet2pt2p5_jesdn_RelBal=0.0;
	TauC_Inc_0j_EnergyWgt_jesup_RelBal=-9999.0, TauB_Inc_0j_pTWgt_jesup_RelBal=-9999.0, TauC_Inc_0j_EnergyWgt_jesdn_RelBal=-9999.0, TauB_Inc_0j_pTWgt_jesdn_RelBal=-9999.0, njets_pt30_eta4p7_jesup_RelBal=0, njets_pt30_eta2p5_jesup_RelBal=0, njets_pt30_eta4p7_jesdn_RelBal=0, njets_pt30_eta2p5_jesdn_RelBal=0, pTj1_jesup_RelBal=-9999.0, pt_leadingjet_pt30_eta4p7_jesup_RelBal=-9999.0, pTj1_2p5_jesup_RelBal=-9999.0, pt_leadingjet_pt30_eta2p5_jesup_RelBal=-9999.0, pTj2_jesup_RelBal=-9999.0, dEtaj1j2_jesup_RelBal=-9999.0, yj1_jesup_RelBal=-9999.0, yj2_jesup_RelBal=-9999.0, dPhiHj1_jesup_RelBal=-9999.0, dyHj1_jesup_RelBal=-9999.0, mass4lj_jesup_RelBal=-9999.0, pT4lj_jesup_RelBal=-9999.0, mass4ljj_jesup_RelBal=-9999.0, pT4ljj_jesup_RelBal=-9999.0, mass4lj_2p5_jesup_RelBal=-9999.0, pT4lj_2p5_jesup_RelBal=-9999.0, mass4ljj_2p5_jesup_RelBal=-9999.0, pT4ljj_2p5_jesup_RelBal=-9999.0, mj1j2_jesup_RelBal=-9999.0, dEtaj1j2_jesup_RelBal=-9999.0, dPhij1j2_jesup_RelBal=-9999.0, dPhiHj1j2_jesup_RelBal=-9999.0, yj1_2p5_jesup_RelBal=-9999.0, yj2_2p5_jesup_RelBal=-9999.0, dPhiHj1_2p5_jesup_RelBal=-9999.0, dyHj1_2p5_jesup_RelBal=-9999.0, mj1j2_2p5_jesup_RelBal=-9999.0, dEtaj1j2_2p5_jesup_RelBal=-9999.0, dPhij1j2_2p5_jesup_RelBal=-9999.0, dPhiHj1j2_2p5_jesup_RelBal=-9999.0, pTj1_jesdn_RelBal=-9999.0, pTj2_jesdn_RelBal=-9999.0, pTj1_2p5_jesdn_RelBal=-9999.0, pt_leadingjet_pt30_eta2p5_jesdn_RelBal=-9999.0, pTj2_jesdn_RelBal=-9999.0, mj1j2_jesdn_RelBal=-9999.0, dEtaj1j2_jesdn_RelBal=-9999.0, yj1_jesdn_RelBal=-9999.0, yj2_jesdn_RelBal=-9999.0, dPhiHj1_jesdn_RelBal=-9999.0, dyHj1_jesdn_RelBal=-9999.0, mj1j2_jesdn_RelBal=-9999.0, dEtaj1j2_jesdn_RelBal=-9999.0, dPhij1j2_jesdn_RelBal=-9999.0, dPhiHj1j2_jesdn_RelBal=-9999.0, yj1_2p5_jesdn_RelBal=-9999.0, yj2_2p5_jesdn_RelBal=-9999.0, dPhiHj1_2p5_jesdn_RelBal=-9999.0, dyHj1_2p5_jesdn_RelBal=-9999.0, mj1j2_2p5_jesdn_RelBal=-9999.0, dEtaj1j2_2p5_jesdn_RelBal=-9999.0, dPhij1j2_2p5_jesdn_RelBal=-9999.0, dPhiHj1j2_2p5_jesdn_RelBal=-9999.0, mass4lj_jesdn_RelBal=-9999.0, pT4lj_jesdn_RelBal=-9999.0, mass4ljj_jesdn_RelBal=-9999.0, pT4ljj_jesdn_RelBal=-9999.0, mass4lj_2p5_jesdn_RelBal=-9999.0, pT4lj_2p5_jesdn_RelBal=-9999.0, mass4ljj_2p5_jesdn_RelBal=-9999.0, pT4ljj_2p5_jesdn_RelBal=-9999.0;  // FIXME 


// initialize RelSample_year variables
        jet1index_jesup_RelSample_year=-1, jet2index_jesup_RelSample_year=-1, jet1index2p5_jesup_RelSample_year=-1, jet2index2p5_jesup_RelSample_year=-1, jet1pt_jesup_RelSample_year=0.0, jet2pt_jesup_RelSample_year=0.0, jet1pt2p5_jesup_RelSample_year=0.0, jet2pt2p5_jesup_RelSample_year=0.0,jet1index_jesdn_RelSample_year=-1, jet2index_jesdn_RelSample_year=-1,jet1index2p5_jesdn_RelSample_year=-1, jet2index2p5_jesdn_RelSample_year=-1, jet1pt_jesdn_RelSample_year=0.0, jet2pt_jesdn_RelSample_year=0.0,jet1pt2p5_jesdn_RelSample_year=0.0, jet2pt2p5_jesdn_RelSample_year=0.0;
	TauC_Inc_0j_EnergyWgt_jesup_RelSample_year=-9999.0, TauB_Inc_0j_pTWgt_jesup_RelSample_year=-9999.0, TauC_Inc_0j_EnergyWgt_jesdn_RelSample_year=-9999.0, TauB_Inc_0j_pTWgt_jesdn_RelSample_year=-9999.0, njets_pt30_eta4p7_jesup_RelSample_year=0, njets_pt30_eta2p5_jesup_RelSample_year=0, njets_pt30_eta4p7_jesdn_RelSample_year=0, njets_pt30_eta2p5_jesdn_RelSample_year=0, pTj1_jesup_RelSample_year=-9999.0, pt_leadingjet_pt30_eta4p7_jesup_RelSample_year=-9999.0, pTj1_2p5_jesup_RelSample_year=-9999.0, pt_leadingjet_pt30_eta2p5_jesup_RelSample_year=-9999.0, pTj2_jesup_RelSample_year=-9999.0, dEtaj1j2_jesup_RelSample_year=-9999.0, yj1_jesup_RelSample_year=-9999.0, yj2_jesup_RelSample_year=-9999.0, dPhiHj1_jesup_RelSample_year=-9999.0, dyHj1_jesup_RelSample_year=-9999.0, mass4lj_jesup_RelSample_year=-9999.0, pT4lj_jesup_RelSample_year=-9999.0, mass4ljj_jesup_RelSample_year=-9999.0, pT4ljj_jesup_RelSample_year=-9999.0, mass4lj_2p5_jesup_RelSample_year=-9999.0, pT4lj_2p5_jesup_RelSample_year=-9999.0, mass4ljj_2p5_jesup_RelSample_year=-9999.0, pT4ljj_2p5_jesup_RelSample_year=-9999.0, mj1j2_jesup_RelSample_year=-9999.0, dEtaj1j2_jesup_RelSample_year=-9999.0, dPhij1j2_jesup_RelSample_year=-9999.0, dPhiHj1j2_jesup_RelSample_year=-9999.0, yj1_2p5_jesup_RelSample_year=-9999.0, yj2_2p5_jesup_RelSample_year=-9999.0, dPhiHj1_2p5_jesup_RelSample_year=-9999.0, dyHj1_2p5_jesup_RelSample_year=-9999.0, mj1j2_2p5_jesup_RelSample_year=-9999.0, dEtaj1j2_2p5_jesup_RelSample_year=-9999.0, dPhij1j2_2p5_jesup_RelSample_year=-9999.0, dPhiHj1j2_2p5_jesup_RelSample_year=-9999.0, pTj1_jesdn_RelSample_year=-9999.0, pTj2_jesdn_RelSample_year=-9999.0, pTj1_2p5_jesdn_RelSample_year=-9999.0, pt_leadingjet_pt30_eta2p5_jesdn_RelSample_year=-9999.0, pTj2_jesdn_RelSample_year=-9999.0, mj1j2_jesdn_RelSample_year=-9999.0, dEtaj1j2_jesdn_RelSample_year=-9999.0, yj1_jesdn_RelSample_year=-9999.0, yj2_jesdn_RelSample_year=-9999.0, dPhiHj1_jesdn_RelSample_year=-9999.0, dyHj1_jesdn_RelSample_year=-9999.0, mj1j2_jesdn_RelSample_year=-9999.0, dEtaj1j2_jesdn_RelSample_year=-9999.0, dPhij1j2_jesdn_RelSample_year=-9999.0, dPhiHj1j2_jesdn_RelSample_year=-9999.0, yj1_2p5_jesdn_RelSample_year=-9999.0, yj2_2p5_jesdn_RelSample_year=-9999.0, dPhiHj1_2p5_jesdn_RelSample_year=-9999.0, dyHj1_2p5_jesdn_RelSample_year=-9999.0, mj1j2_2p5_jesdn_RelSample_year=-9999.0, dEtaj1j2_2p5_jesdn_RelSample_year=-9999.0, dPhij1j2_2p5_jesdn_RelSample_year=-9999.0, dPhiHj1j2_2p5_jesdn_RelSample_year=-9999.0, mass4lj_jesdn_RelSample_year=-9999.0, pT4lj_jesdn_RelSample_year=-9999.0, mass4ljj_jesdn_RelSample_year=-9999.0, pT4ljj_jesdn_RelSample_year=-9999.0, mass4lj_2p5_jesdn_RelSample_year=-9999.0, pT4lj_2p5_jesdn_RelSample_year=-9999.0, mass4ljj_2p5_jesdn_RelSample_year=-9999.0, pT4ljj_2p5_jesdn_RelSample_year=-9999.0;  // FIXME 





        //if (i>=2000000) continue;
        //if (i>=2000000&&_Test) break;
        if (i%100000==0) std::cout<<i<<"/"<<nentries<<std::endl;

        oldtree->GetEntry(i);
//	cout <<"passedFullSelection:   "<<passedFullSelection<<endl;
//	cout <<"lep_Hindex[0]:   "<<lep_Hindex[0]<<endl;
//	cout <<"(*GENlep_MomId)[0]:    "<<(*GENlep_MomId)[0]<<endl;
//	cout <<"(*GENlep_MomMomId)[0]:    "<<(*GENlep_MomMomId)[0]<<endl;

	if (passedFullSelection) {
	isH4l = ( (((*lep_genindex)[passedFullSelection*lep_Hindex[0]]>-0.5)*(*GENlep_MomMomId)[max(0,(*lep_genindex)[passedFullSelection*lep_Hindex[0]])]==25 && ((*lep_genindex)[passedFullSelection*lep_Hindex[0]]>-0.5)*(*GENlep_MomId)[max(0,(*lep_genindex)[passedFullSelection*lep_Hindex[0]])]==23 && ((*lep_genindex)[passedFullSelection*lep_Hindex[1]]>-0.5)*(*GENlep_MomMomId)[max(0,(*lep_genindex)[passedFullSelection*lep_Hindex[1]])]==25 && ((*lep_genindex)[passedFullSelection*lep_Hindex[1]]>-0.5)*(*GENlep_MomId)[max(0,(*lep_genindex)[passedFullSelection*lep_Hindex[1]])]==23 && ((*lep_genindex)[passedFullSelection*lep_Hindex[2]]>-0.5)*(*GENlep_MomMomId)[max(0,(*lep_genindex)[passedFullSelection*lep_Hindex[2]])]==25 && ((*lep_genindex)[passedFullSelection*lep_Hindex[2]]>-0.5)*(*GENlep_MomId)[max(0,(*lep_genindex)[passedFullSelection*lep_Hindex[2]])]==23 && ((*lep_genindex)[passedFullSelection*lep_Hindex[3]]>-0.5)*(*GENlep_MomMomId)[max(0,(*lep_genindex)[passedFullSelection*lep_Hindex[3]])]==25 && ((*lep_genindex)[passedFullSelection*lep_Hindex[3]]>-0.5)*(*GENlep_MomId)[max(0,(*lep_genindex)[passedFullSelection*lep_Hindex[3]])]==23) );
	}


	if (redoJets) {

           //     vector<float> goodJets_JECJER_pt30_eta4p7 {}; 
                vector<float> jes_unc_split {};
//up    
                vector<float> pt_jesup_split {};
//dn
                vector<float> pt_jesdn_split {};


                vector<float> eta_jes_split {};
                vector<float> phi_jes_split {};
                vector<float> mass_jes_split {};

                float singleContr_jes_unc ;

                TLorentzVector thisJet;
// tempo variables for TauC TauB nominal values

                for( unsigned int k = 0; k<(*jet_iscleanH4l).size(); k++) {
                    if ((*jet_pt)[k]<30.0 || abs((*jet_eta)[k])>4.7) continue;
                    thisJet.SetPtEtaPhiM((*jet_pt)[((*jet_iscleanH4l)[k])],(*jet_eta)[((*jet_iscleanH4l)[k])],(*jet_phi)[((*jet_iscleanH4l)[k])],(*jet_mass)[((*jet_iscleanH4l)[k])]);
		    // validation variable pTj1
		    if ((*jet_pt)[((*jet_iscleanH4l)[k])] > pTj1) {
			pTj1 = (*jet_pt)[((*jet_iscleanH4l)[k])];
			}
			
			pt_nom.push_back(thisJet.Pt());
                    //if(applyJEC_ && isMC)
                    if(applyJEC_)
                      {
                        for (unsigned s_unc = 0; s_unc < uncSources.size(); s_unc++)
                          {
                            singleContr_jes_unc = 0;
                            splittedUncerts_[s_unc]->setJetEta(thisJet.Eta());
                            splittedUncerts_[s_unc]->setJetPt(thisJet.Pt());
                            singleContr_jes_unc = splittedUncerts_[s_unc]->getUncertainty(true);
                            jes_unc_split.push_back(singleContr_jes_unc);
                            pt_jesup_split.push_back(thisJet.Pt() * (1.0 + singleContr_jes_unc));
                            pt_jesdn_split.push_back(thisJet.Pt() * (1.0 - singleContr_jes_unc));
                            eta_jes_split.push_back(thisJet.Eta());
                            phi_jes_split.push_back(thisJet.Phi());
                            mass_jes_split.push_back(thisJet.M());

                            cout<<"Uncertainty source is:  uncSources["<<s_unc<<"]  "<<uncSources[s_unc]<<endl;
                            cout<<"thisJet.Pt():   "<<thisJet.Pt()<<endl;
                            cout<<"jes_unc_split["<<s_unc<<"]  "<<jes_unc_split[s_unc]<<endl;
                            cout<<"pt_jesup_split["<<s_unc<<"]  "<<pt_jesup_split[s_unc]<<endl;
                            cout<<"pt_jesdn_split["<<s_unc<<"]  "<<pt_jesdn_split[s_unc]<<endl;
                          }
                      }
                    else
                      {
                        for(unsigned s_unc = 0; s_unc < uncSources.size(); s_unc++)
                          {
                            jes_unc_split.push_back(-999.);
                            pt_jesup_split.push_back(-999.);
                            pt_jesdn_split.push_back(-999.);

                            eta_jes_split.push_back(-999.); 
                            phi_jes_split.push_back(-999.);
                            mass_jes_split.push_back(-999.);			
                          }
                      }
//filling variables
        //unc, pt_up, pt_dn,eta, phi, mass
                      jes_unc_split_Total.push_back(jes_unc_split[0]);jes_unc_split_Abs.push_back(jes_unc_split[1]); jes_unc_split_Abs_year.push_back(jes_unc_split[2]); jes_unc_split_BBEC1.push_back(jes_unc_split[3]); jes_unc_split_BBEC1_year.push_back(jes_unc_split[4]);jes_unc_split_EC2.push_back(jes_unc_split[5]); jes_unc_split_EC2_year.push_back(jes_unc_split[6]); jes_unc_split_FlavQCD.push_back(jes_unc_split[7]); jes_unc_split_HF.push_back(jes_unc_split[8]); jes_unc_split_HF_year.push_back(jes_unc_split[9]);jes_unc_split_RelBal.push_back(jes_unc_split[10]); jes_unc_split_RelSample_year.push_back(jes_unc_split[11]);
                      pt_jesup_split_Total.push_back(pt_jesup_split[0]);pt_jesup_split_Abs.push_back(pt_jesup_split[1]); pt_jesup_split_Abs_year.push_back(pt_jesup_split[2]); pt_jesup_split_BBEC1.push_back(pt_jesup_split[3]); pt_jesup_split_BBEC1_year.push_back(pt_jesup_split[4]);pt_jesup_split_EC2.push_back(pt_jesup_split[5]); pt_jesup_split_EC2_year.push_back(pt_jesup_split[6]); pt_jesup_split_FlavQCD.push_back(pt_jesup_split[7]); pt_jesup_split_HF.push_back(pt_jesup_split[8]); pt_jesup_split_HF_year.push_back(pt_jesup_split[9]);pt_jesup_split_RelBal.push_back(pt_jesup_split[10]); pt_jesup_split_RelSample_year.push_back(pt_jesup_split[11]);
                      pt_jesdn_split_Total.push_back(pt_jesdn_split[0]);pt_jesdn_split_Abs.push_back(pt_jesdn_split[1]); pt_jesdn_split_Abs_year.push_back(pt_jesdn_split[2]); pt_jesdn_split_BBEC1.push_back(pt_jesdn_split[3]); pt_jesdn_split_BBEC1_year.push_back(pt_jesdn_split[4]);pt_jesdn_split_EC2.push_back(pt_jesdn_split[5]); pt_jesdn_split_EC2_year.push_back(pt_jesdn_split[6]); pt_jesdn_split_FlavQCD.push_back(pt_jesdn_split[7]); pt_jesdn_split_HF.push_back(pt_jesdn_split[8]); pt_jesdn_split_HF_year.push_back(pt_jesdn_split[9]);pt_jesdn_split_RelBal.push_back(pt_jesdn_split[10]); pt_jesdn_split_RelSample_year.push_back(pt_jesdn_split[11]);
                      eta_jes_split_Total.push_back(eta_jes_split[0]);eta_jes_split_Abs.push_back(eta_jes_split[1]); eta_jes_split_Abs_year.push_back(eta_jes_split[2]); eta_jes_split_BBEC1.push_back(eta_jes_split[3]); eta_jes_split_BBEC1_year.push_back(eta_jes_split[4]);eta_jes_split_EC2.push_back(eta_jes_split[5]); eta_jes_split_EC2_year.push_back(eta_jes_split[6]); eta_jes_split_FlavQCD.push_back(eta_jes_split[7]); eta_jes_split_HF.push_back(eta_jes_split[8]); eta_jes_split_HF_year.push_back(eta_jes_split[9]);eta_jes_split_RelBal.push_back(eta_jes_split[10]); eta_jes_split_RelSample_year.push_back(eta_jes_split[11]);		
                      phi_jes_split_Total.push_back(phi_jes_split[0]);phi_jes_split_Abs.push_back(phi_jes_split[1]); phi_jes_split_Abs_year.push_back(phi_jes_split[2]); phi_jes_split_BBEC1.push_back(phi_jes_split[3]); phi_jes_split_BBEC1_year.push_back(phi_jes_split[4]);phi_jes_split_EC2.push_back(phi_jes_split[5]); phi_jes_split_EC2_year.push_back(phi_jes_split[6]); phi_jes_split_FlavQCD.push_back(phi_jes_split[7]); phi_jes_split_HF.push_back(phi_jes_split[8]); phi_jes_split_HF_year.push_back(phi_jes_split[9]);phi_jes_split_RelBal.push_back(phi_jes_split[10]); phi_jes_split_RelSample_year.push_back(phi_jes_split[11]);		
                      mass_jes_split_Total.push_back(mass_jes_split[0]);mass_jes_split_Abs.push_back(mass_jes_split[1]); mass_jes_split_Abs_year.push_back(mass_jes_split[2]); mass_jes_split_BBEC1.push_back(mass_jes_split[3]); mass_jes_split_BBEC1_year.push_back(mass_jes_split[4]);mass_jes_split_EC2.push_back(mass_jes_split[5]); mass_jes_split_EC2_year.push_back(mass_jes_split[6]); mass_jes_split_FlavQCD.push_back(mass_jes_split[7]); mass_jes_split_HF.push_back(mass_jes_split[8]); mass_jes_split_HF_year.push_back(mass_jes_split[9]);mass_jes_split_RelBal.push_back(mass_jes_split[10]); mass_jes_split_RelSample_year.push_back(mass_jes_split[11]);		

	              jes_unc_split.clear(); pt_jesup_split.clear();pt_jesdn_split.clear();eta_jes_split.clear();phi_jes_split.clear();mass_jes_split.clear();

		      } // end jet loop

       	       }  // redoJets
////////////////////////////////////////////////////////////////////////////
                TLorentzVector Higgs;
                Higgs.SetPtEtaPhiM(pT4l, eta4l, phi4l, mass4l);

                TauC_Inc_0j_EnergyWgt_nom = TauC(pt_nom, eta_jes_split_Abs,phi_jes_split_Abs,mass_jes_split_Abs, Higgs);
                TauB_Inc_0j_pTWgt_nom = TauB(pt_nom, eta_jes_split_Abs,phi_jes_split_Abs,mass_jes_split_Abs, Higgs);
		pt_nom.clear();  
////////////////////////////////////////////////////////////////////////////////
////    Starting "Abs" JES uncertainty source
////////////////////////////////////////////////////////////////////////////////

                TauC_Inc_0j_EnergyWgt_jesup_Abs = TauC(pt_jesup_split_Abs, eta_jes_split_Abs,phi_jes_split_Abs,mass_jes_split_Abs, Higgs);
                TauB_Inc_0j_pTWgt_jesup_Abs = TauB(pt_jesup_split_Abs, eta_jes_split_Abs,phi_jes_split_Abs,mass_jes_split_Abs, Higgs);


                for( unsigned int k = 0; k< pt_jesup_split_Abs.size(); k++) {
                    if (pt_jesup_split_Abs[k]<30.0 || abs(eta_jes_split_Abs[k])>4.7) continue;
		    TLorentzVector thisJet_jesup_Abs;
		    thisJet_jesup_Abs.SetPtEtaPhiM(pt_jesup_split_Abs[k],eta_jes_split_Abs[k],phi_jes_split_Abs[k],mass_jes_split_Abs[k]);

                        njets_pt30_eta4p7_jesup_Abs+=1;  

                        if (thisJet_jesup_Abs.Pt()>jet1pt_jesup_Abs) {
                            jet2pt_jesup_Abs=jet1pt_jesup_Abs; jet2index_jesup_Abs=jet1index_jesup_Abs;
                            jet1pt_jesup_Abs=thisJet_jesup_Abs.Pt(); jet1index_jesup_Abs=k;
                        } else if (thisJet_jesup_Abs.Pt()>jet2pt_jesup_Abs) {
                            jet2pt_jesup_Abs=thisJet_jesup_Abs.Pt(); jet2index_jesup_Abs=k;
                        }
                        if (abs(thisJet_jesup_Abs.Eta())<2.5) {
                            njets_pt30_eta2p5_jesup_Abs+=1;
                            if (thisJet_jesup_Abs.Pt()>jet1pt2p5_jesup_Abs) {
                                jet2pt2p5_jesup_Abs=jet1pt2p5_jesup_Abs; jet2index2p5_jesup_Abs=jet1index2p5_jesup_Abs;
                                jet1pt2p5_jesup_Abs=thisJet_jesup_Abs.Pt(); jet1index2p5_jesup_Abs=k;
                            } else if (thisJet_jesup_Abs.Pt()>jet2pt2p5_jesup_Abs) {
                                jet2pt2p5_jesup_Abs=thisJet_jesup_Abs.Pt(); jet2index2p5_jesup_Abs=k;
                            }
                        }
                }
                cout<<njets_pt30_eta4p7_jesup_Abs<<" jets (jesup_Abs)"<<endl;


		TLorentzVector Jet1_jesup_Abs, Jet1_2p5_jesup_Abs, Jet2_jesup_Abs, Jet2_2p5_jesup_Abs;
                if (njets_pt30_eta4p7_jesup_Abs > 0) {
		    Jet1_jesup_Abs.SetPtEtaPhiM(pt_jesup_split_Abs[jet1index_jesup_Abs],eta_jes_split_Abs[jet1index_jesup_Abs],phi_jes_split_Abs[jet1index_jesup_Abs], mass_jes_split_Abs[jet1index_jesup_Abs]);

                    pt_leadingjet_pt30_eta4p7_jesup_Abs=Jet1_jesup_Abs.Pt(); 
                    pTj1_jesup_Abs=Jet1_jesup_Abs.Pt(); 
                    etaj1_jesup_Abs=Jet1_jesup_Abs.Eta();
                    yj1_jesup_Abs=Jet1_jesup_Abs.Rapidity();
		    pT4lj_jesup_Abs=(Higgs+Jet1_jesup_Abs).Pt();
		    mass4lj_jesup_Abs=(Higgs+Jet1_jesup_Abs).M();
		    dPhiHj1_jesup_Abs=deltaPhi(Higgs.Phi(),Jet1_jesup_Abs.Phi());
                    dyHj1_jesup_Abs=TMath::Abs(rapidity4l-yj1_jesup_Abs);
                }
                if (njets_pt30_eta4p7_jesup_Abs > 1) {
                    Jet2_jesup_Abs.SetPtEtaPhiM(pt_jesup_split_Abs[jet2index_jesup_Abs],eta_jes_split_Abs[jet2index_jesup_Abs],phi_jes_split_Abs[jet2index_jesup_Abs], mass_jes_split_Abs[jet2index_jesup_Abs]);

                    pTj2_jesup_Abs=Jet2_jesup_Abs.Pt();
                    etaj2_jesup_Abs=Jet2_jesup_Abs.Eta();
                    yj2_jesup_Abs=Jet2_jesup_Abs.Rapidity();
                    pT4ljj_jesup_Abs=(Higgs+Jet1_jesup_Abs+Jet2_jesup_Abs).Pt();
                    mass4ljj_jesup_Abs=(Higgs+Jet1_jesup_Abs+Jet2_jesup_Abs).M();
		    mj1j2_jesup_Abs=(Jet1_jesup_Abs+Jet2_jesup_Abs).M();
                    dEtaj1j2_jesup_Abs=TMath::Abs(Jet1_jesup_Abs.Eta()-Jet2_jesup_Abs.Eta());
                    dPhij1j2_jesup_Abs=deltaPhi(Jet1_jesup_Abs.Phi(),Jet2_jesup_Abs.Phi());
                    dPhiHj1j2_jesup_Abs=deltaPhi(Higgs.Phi(),(Jet1_jesup_Abs+Jet2_jesup_Abs).Phi());
                }
                if (njets_pt30_eta2p5_jesup_Abs > 0) {
                    Jet1_2p5_jesup_Abs.SetPtEtaPhiM(pt_jesup_split_Abs[jet1index2p5_jesup_Abs],eta_jes_split_Abs[jet1index2p5_jesup_Abs],phi_jes_split_Abs[jet1index2p5_jesup_Abs], mass_jes_split_Abs[jet1index2p5_jesup_Abs]);
                    pt_leadingjet_pt30_eta2p5_jesup_Abs=Jet1_2p5_jesup_Abs.Pt();
                    pTj1_2p5_jesup_Abs=Jet1_2p5_jesup_Abs.Pt();
                    yj1_2p5_jesup_Abs=Jet1_2p5_jesup_Abs.Rapidity();
                    pT4lj_2p5_jesup_Abs=(Higgs+Jet1_2p5_jesup_Abs).Pt();
                    mass4lj_2p5_jesup_Abs=(Higgs+Jet1_2p5_jesup_Abs).M();
		    dPhiHj1_2p5_jesup_Abs=deltaPhi(Higgs.Phi(),Jet1_2p5_jesup_Abs.Phi());
                    dyHj1_2p5_jesup_Abs=TMath::Abs(rapidity4l-yj1_2p5_jesup_Abs);
                }
                if (njets_pt30_eta2p5_jesup_Abs > 1) {
                    Jet2_2p5_jesup_Abs.SetPtEtaPhiM(pt_jesup_split_Abs[jet2index2p5_jesup_Abs],eta_jes_split_Abs[jet2index2p5_jesup_Abs],phi_jes_split_Abs[jet2index2p5_jesup_Abs], mass_jes_split_Abs[jet2index2p5_jesup_Abs]);
                    pTj2_2p5_jesup_Abs=Jet2_2p5_jesup_Abs.Pt();
                    yj2_2p5_jesup_Abs=Jet2_2p5_jesup_Abs.Rapidity();
                    pT4ljj_2p5_jesup_Abs=(Higgs+Jet1_2p5_jesup_Abs+Jet2_2p5_jesup_Abs).Pt();
                    mass4ljj_2p5_jesup_Abs=(Higgs+Jet1_2p5_jesup_Abs+Jet2_2p5_jesup_Abs).M();
                    mj1j2_2p5_jesup_Abs=(Jet1_2p5_jesup_Abs+Jet2_2p5_jesup_Abs).M();                    
                    dEtaj1j2_2p5_jesup_Abs=TMath::Abs(Jet1_2p5_jesup_Abs.Eta()-Jet2_2p5_jesup_Abs.Eta());
                    dPhij1j2_2p5_jesup_Abs=deltaPhi(Jet1_2p5_jesup_Abs.Phi(),Jet2_2p5_jesup_Abs.Phi());
                    dPhiHj1j2_2p5_jesup_Abs=deltaPhi(Higgs.Phi(),(Jet1_2p5_jesup_Abs+Jet2_2p5_jesup_Abs).Phi());

                }

		pt_jesup_split_Abs.clear();  

/////////////////
// Abs dn start
                TauC_Inc_0j_EnergyWgt_jesdn_Abs = TauC(pt_jesdn_split_Abs, eta_jes_split_Abs,phi_jes_split_Abs,mass_jes_split_Abs, Higgs);
                TauB_Inc_0j_pTWgt_jesdn_Abs = TauB(pt_jesdn_split_Abs, eta_jes_split_Abs,phi_jes_split_Abs,mass_jes_split_Abs, Higgs);


                for( unsigned int k = 0; k< pt_jesdn_split_Abs.size(); k++) {
                    if (pt_jesdn_split_Abs[k]<30.0 || abs(eta_jes_split_Abs[k])>4.7) continue;
                    TLorentzVector thisJet_jesdn_Abs;
		    thisJet_jesdn_Abs.SetPtEtaPhiM(pt_jesdn_split_Abs[k],eta_jes_split_Abs[k],phi_jes_split_Abs[k],mass_jes_split_Abs[k]);

                    bool isclean_H4l=true;
                        njets_pt30_eta4p7_jesdn_Abs+=1;
			if (thisJet_jesdn_Abs.Pt()>jet1pt_jesdn_Abs) {
                            jet2pt_jesdn_Abs=jet1pt_jesdn_Abs; jet2index_jesdn_Abs=jet1index_jesdn_Abs;
                            jet1pt_jesdn_Abs=thisJet_jesdn_Abs.Pt(); jet1index_jesdn_Abs=k;
                        } else if (thisJet_jesdn_Abs.Pt()>jet2pt_jesdn_Abs) {
                            jet2pt_jesdn_Abs=thisJet_jesdn_Abs.Pt(); jet2index_jesdn_Abs=k;
                        }
                        if (abs(thisJet_jesdn_Abs.Eta())<2.5) {
                            njets_pt30_eta2p5_jesdn_Abs+=1;
                            if (thisJet_jesdn_Abs.Pt()>jet1pt2p5_jesdn_Abs) {
                                jet2pt2p5_jesdn_Abs=jet1pt2p5_jesdn_Abs; jet2index2p5_jesdn_Abs=jet1index2p5_jesdn_Abs;
                                jet1pt2p5_jesdn_Abs=thisJet_jesdn_Abs.Pt(); jet1index2p5_jesdn_Abs=k;
                            } else if (thisJet_jesdn_Abs.Pt()>jet2pt2p5_jesdn_Abs) {
                                jet2pt2p5_jesdn_Abs=thisJet_jesdn_Abs.Pt(); jet2index2p5_jesdn_Abs=k;
                            }
                        }
                }

// Filling Abs dn variables

                TLorentzVector Jet1_jesdn_Abs, Jet1_2p5_jesdn_Abs, Jet2_jesdn_Abs, Jet2_2p5_jesdn_Abs;

                if (njets_pt30_eta4p7_jesdn_Abs > 0) { 
                    Jet1_jesdn_Abs.SetPtEtaPhiM(pt_jesdn_split_Abs[jet1index_jesdn_Abs],eta_jes_split_Abs[jet1index_jesdn_Abs],phi_jes_split_Abs[jet1index_jesdn_Abs], mass_jes_split_Abs[jet1index_jesdn_Abs]);
                    pt_leadingjet_pt30_eta4p7_jesdn_Abs=Jet1_jesdn_Abs.Pt(); 
                    pTj1_jesdn_Abs=Jet1_jesdn_Abs.Pt(); 
                    etaj1_jesdn_Abs=Jet1_jesdn_Abs.Eta();
                    yj1_jesdn_Abs=Jet1_jesdn_Abs.Rapidity();
                    pT4lj_jesdn_Abs=(Higgs+Jet1_jesdn_Abs).Pt();
                    mass4lj_jesdn_Abs=(Higgs+Jet1_jesdn_Abs).M();
                    dPhiHj1_jesdn_Abs=deltaPhi(Higgs.Phi(),Jet1_jesdn_Abs.Phi());
                    dyHj1_jesdn_Abs=TMath::Abs(rapidity4l-yj1_jesdn_Abs);
                }    
                if (njets_pt30_eta4p7_jesdn_Abs > 1) { 
                    Jet2_jesdn_Abs.SetPtEtaPhiM(pt_jesdn_split_Abs[jet2index_jesdn_Abs],eta_jes_split_Abs[jet2index_jesdn_Abs],phi_jes_split_Abs[jet2index_jesdn_Abs], mass_jes_split_Abs[jet2index_jesdn_Abs]);

                    pTj2_jesdn_Abs=Jet2_jesdn_Abs.Pt();
                    etaj2_jesdn_Abs=Jet2_jesdn_Abs.Eta();
                    yj2_jesdn_Abs=Jet2_jesdn_Abs.Rapidity();
                    pT4ljj_jesdn_Abs=(Higgs+Jet1_jesdn_Abs+Jet2_jesdn_Abs).Pt();
                    mass4ljj_jesdn_Abs=(Higgs+Jet1_jesdn_Abs+Jet2_jesdn_Abs).M();
                    dEtaj1j2_jesdn_Abs=TMath::Abs(Jet1_jesdn_Abs.Eta()-Jet2_jesdn_Abs.Eta());
                    dPhij1j2_jesdn_Abs=deltaPhi(Jet1_jesdn_Abs.Phi(),Jet2_jesdn_Abs.Phi());
                    dPhiHj1j2_jesdn_Abs=deltaPhi(Higgs.Phi(),(Jet1_jesdn_Abs+Jet2_jesdn_Abs).Phi());
                }    
                if (njets_pt30_eta2p5_jesdn_Abs > 0) { 
                    Jet1_2p5_jesdn_Abs.SetPtEtaPhiM(pt_jesdn_split_Abs[jet1index2p5_jesdn_Abs],eta_jes_split_Abs[jet1index2p5_jesdn_Abs],phi_jes_split_Abs[jet1index2p5_jesdn_Abs], mass_jes_split_Abs[jet1index2p5_jesdn_Abs]);
                    pt_leadingjet_pt30_eta2p5_jesdn_Abs=Jet1_2p5_jesdn_Abs.Pt();
                    pTj1_2p5_jesdn_Abs=Jet1_2p5_jesdn_Abs.Pt();
                    yj1_2p5_jesdn_Abs=Jet1_2p5_jesdn_Abs.Rapidity();
                    pT4lj_2p5_jesdn_Abs=(Higgs+Jet1_2p5_jesdn_Abs).Pt();
                    mass4lj_2p5_jesdn_Abs=(Higgs+Jet1_2p5_jesdn_Abs).M();
		    dPhiHj1_2p5_jesdn_Abs=deltaPhi(Higgs.Phi(),Jet1_2p5_jesdn_Abs.Phi());
                    dyHj1_2p5_jesdn_Abs=TMath::Abs(rapidity4l-yj1_2p5_jesdn_Abs);
                }    
                if (njets_pt30_eta2p5_jesdn_Abs > 1) { 
                    Jet2_2p5_jesdn_Abs.SetPtEtaPhiM(pt_jesdn_split_Abs[jet2index2p5_jesdn_Abs],eta_jes_split_Abs[jet2index2p5_jesdn_Abs],phi_jes_split_Abs[jet2index2p5_jesdn_Abs], mass_jes_split_Abs[jet2index2p5_jesdn_Abs]);
                    pTj2_2p5_jesdn_Abs=Jet2_2p5_jesdn_Abs.Pt();
                    yj2_2p5_jesdn_Abs=Jet2_2p5_jesdn_Abs.Rapidity();
                    pT4ljj_2p5_jesdn_Abs=(Higgs+Jet1_2p5_jesdn_Abs+Jet2_2p5_jesdn_Abs).Pt();
                    mass4ljj_2p5_jesdn_Abs=(Higgs+Jet1_2p5_jesdn_Abs+Jet2_2p5_jesdn_Abs).M();
                    mj1j2_2p5_jesdn_Abs=(Jet1_2p5_jesdn_Abs+Jet2_2p5_jesdn_Abs).M();     
                    dEtaj1j2_2p5_jesdn_Abs=TMath::Abs(Jet1_2p5_jesdn_Abs.Eta()-Jet2_2p5_jesdn_Abs.Eta());
	            dPhij1j2_2p5_jesdn_Abs=deltaPhi(Jet1_2p5_jesdn_Abs.Phi(),Jet2_2p5_jesdn_Abs.Phi());
                    dPhiHj1j2_2p5_jesdn_Abs=deltaPhi(Higgs.Phi(),(Jet1_2p5_jesdn_Abs+Jet2_2p5_jesdn_Abs).Phi());

                }    

                pt_jesdn_split_Abs.clear();

		eta_jes_split_Abs.clear();
		phi_jes_split_Abs.clear();
		mass_jes_split_Abs.clear();


////////////////////////////////////////////////////////////////////////////////
////    Starting "Abs_year" JES uncertainty source
////////////////////////////////////////////////////////////////////////////////

                TauC_Inc_0j_EnergyWgt_jesup_Abs_year = TauC(pt_jesup_split_Abs_year, eta_jes_split_Abs_year,phi_jes_split_Abs_year,mass_jes_split_Abs_year, Higgs);
                TauB_Inc_0j_pTWgt_jesup_Abs_year = TauB(pt_jesup_split_Abs_year, eta_jes_split_Abs_year,phi_jes_split_Abs_year,mass_jes_split_Abs_year, Higgs);

                for( unsigned int k = 0; k< pt_jesup_split_Abs_year.size(); k++) {
                    if (pt_jesup_split_Abs_year[k]<30.0 || abs(eta_jes_split_Abs_year[k])>4.7) continue;
		    TLorentzVector thisJet_jesup_Abs_year;
		    thisJet_jesup_Abs_year.SetPtEtaPhiM(pt_jesup_split_Abs_year[k],eta_jes_split_Abs_year[k],phi_jes_split_Abs_year[k],mass_jes_split_Abs_year[k]);

                        njets_pt30_eta4p7_jesup_Abs_year+=1;  

                        if (thisJet_jesup_Abs_year.Pt()>jet1pt_jesup_Abs_year) {
                            jet2pt_jesup_Abs_year=jet1pt_jesup_Abs_year; jet2index_jesup_Abs_year=jet1index_jesup_Abs_year;
                            jet1pt_jesup_Abs_year=thisJet_jesup_Abs_year.Pt(); jet1index_jesup_Abs_year=k;
                        } else if (thisJet_jesup_Abs_year.Pt()>jet2pt_jesup_Abs_year) {
                            jet2pt_jesup_Abs_year=thisJet_jesup_Abs_year.Pt(); jet2index_jesup_Abs_year=k;
                        }
                        if (abs(thisJet_jesup_Abs_year.Eta())<2.5) {
                            njets_pt30_eta2p5_jesup_Abs_year+=1;
                            if (thisJet_jesup_Abs_year.Pt()>jet1pt2p5_jesup_Abs_year) {
                                jet2pt2p5_jesup_Abs_year=jet1pt2p5_jesup_Abs_year; jet2index2p5_jesup_Abs_year=jet1index2p5_jesup_Abs_year;
                                jet1pt2p5_jesup_Abs_year=thisJet_jesup_Abs_year.Pt(); jet1index2p5_jesup_Abs_year=k;
                            } else if (thisJet_jesup_Abs_year.Pt()>jet2pt2p5_jesup_Abs_year) {
                                jet2pt2p5_jesup_Abs_year=thisJet_jesup_Abs_year.Pt(); jet2index2p5_jesup_Abs_year=k;
                            }
                        }
                }
                cout<<njets_pt30_eta4p7_jesup_Abs_year<<" jets (jesup_Abs_year)"<<endl;


		TLorentzVector Jet1_jesup_Abs_year, Jet1_2p5_jesup_Abs_year, Jet2_jesup_Abs_year, Jet2_2p5_jesup_Abs_year;
                if (njets_pt30_eta4p7_jesup_Abs_year > 0) {
		    Jet1_jesup_Abs_year.SetPtEtaPhiM(pt_jesup_split_Abs_year[jet1index_jesup_Abs_year],eta_jes_split_Abs_year[jet1index_jesup_Abs_year],phi_jes_split_Abs_year[jet1index_jesup_Abs_year], mass_jes_split_Abs_year[jet1index_jesup_Abs_year]);

                    pt_leadingjet_pt30_eta4p7_jesup_Abs_year=Jet1_jesup_Abs_year.Pt(); 
                    pTj1_jesup_Abs_year=Jet1_jesup_Abs_year.Pt(); 
                    etaj1_jesup_Abs_year=Jet1_jesup_Abs_year.Eta();
                    yj1_jesup_Abs_year=Jet1_jesup_Abs_year.Rapidity();
		    pT4lj_jesup_Abs_year=(Higgs+Jet1_jesup_Abs_year).Pt();
		    mass4lj_jesup_Abs_year=(Higgs+Jet1_jesup_Abs_year).M();
		    dPhiHj1_jesup_Abs_year=deltaPhi(Higgs.Phi(),Jet1_jesup_Abs_year.Phi());
                    dyHj1_jesup_Abs_year=TMath::Abs(rapidity4l-yj1_jesup_Abs_year);
                }
                if (njets_pt30_eta4p7_jesup_Abs_year > 1) {
                    Jet2_jesup_Abs_year.SetPtEtaPhiM(pt_jesup_split_Abs_year[jet2index_jesup_Abs_year],eta_jes_split_Abs_year[jet2index_jesup_Abs_year],phi_jes_split_Abs_year[jet2index_jesup_Abs_year], mass_jes_split_Abs_year[jet2index_jesup_Abs_year]);

                    pTj2_jesup_Abs_year=Jet2_jesup_Abs_year.Pt();
                    etaj2_jesup_Abs_year=Jet2_jesup_Abs_year.Eta();
                    yj2_jesup_Abs_year=Jet2_jesup_Abs_year.Rapidity();
                    pT4ljj_jesup_Abs_year=(Higgs+Jet1_jesup_Abs_year+Jet2_jesup_Abs_year).Pt();
                    mass4ljj_jesup_Abs_year=(Higgs+Jet1_jesup_Abs_year+Jet2_jesup_Abs_year).M();
		    mj1j2_jesup_Abs_year=(Jet1_jesup_Abs_year+Jet2_jesup_Abs_year).M();
                    dEtaj1j2_jesup_Abs_year=TMath::Abs(Jet1_jesup_Abs_year.Eta()-Jet2_jesup_Abs_year.Eta());
                    dPhij1j2_jesup_Abs_year=deltaPhi(Jet1_jesup_Abs_year.Phi(),Jet2_jesup_Abs_year.Phi());
                    dPhiHj1j2_jesup_Abs_year=deltaPhi(Higgs.Phi(),(Jet1_jesup_Abs_year+Jet2_jesup_Abs_year).Phi());
                }
                if (njets_pt30_eta2p5_jesup_Abs_year > 0) {
                    Jet1_2p5_jesup_Abs_year.SetPtEtaPhiM(pt_jesup_split_Abs_year[jet1index2p5_jesup_Abs_year],eta_jes_split_Abs_year[jet1index2p5_jesup_Abs_year],phi_jes_split_Abs_year[jet1index2p5_jesup_Abs_year], mass_jes_split_Abs_year[jet1index2p5_jesup_Abs_year]);
                    pt_leadingjet_pt30_eta2p5_jesup_Abs_year=Jet1_2p5_jesup_Abs_year.Pt();
                    pTj1_2p5_jesup_Abs_year=Jet1_2p5_jesup_Abs_year.Pt();
                    yj1_2p5_jesup_Abs_year=Jet1_2p5_jesup_Abs_year.Rapidity();
                    pT4lj_2p5_jesup_Abs_year=(Higgs+Jet1_2p5_jesup_Abs_year).Pt();
                    mass4lj_2p5_jesup_Abs_year=(Higgs+Jet1_2p5_jesup_Abs_year).M();
		    dPhiHj1_2p5_jesup_Abs_year=deltaPhi(Higgs.Phi(),Jet1_2p5_jesup_Abs_year.Phi());
                    dyHj1_2p5_jesup_Abs_year=TMath::Abs(rapidity4l-yj1_2p5_jesup_Abs_year);
                }
                if (njets_pt30_eta2p5_jesup_Abs_year > 1) {
                    Jet2_2p5_jesup_Abs_year.SetPtEtaPhiM(pt_jesup_split_Abs_year[jet2index2p5_jesup_Abs_year],eta_jes_split_Abs_year[jet2index2p5_jesup_Abs_year],phi_jes_split_Abs_year[jet2index2p5_jesup_Abs_year], mass_jes_split_Abs_year[jet2index2p5_jesup_Abs_year]);
                    pTj2_2p5_jesup_Abs_year=Jet2_2p5_jesup_Abs_year.Pt();
                    yj2_2p5_jesup_Abs_year=Jet2_2p5_jesup_Abs_year.Rapidity();
                    pT4ljj_2p5_jesup_Abs_year=(Higgs+Jet1_2p5_jesup_Abs_year+Jet2_2p5_jesup_Abs_year).Pt();
                    mass4ljj_2p5_jesup_Abs_year=(Higgs+Jet1_2p5_jesup_Abs_year+Jet2_2p5_jesup_Abs_year).M();
                    mj1j2_2p5_jesup_Abs_year=(Jet1_2p5_jesup_Abs_year+Jet2_2p5_jesup_Abs_year).M();                    
                    dEtaj1j2_2p5_jesup_Abs_year=TMath::Abs(Jet1_2p5_jesup_Abs_year.Eta()-Jet2_2p5_jesup_Abs_year.Eta());
                    dPhij1j2_2p5_jesup_Abs_year=deltaPhi(Jet1_2p5_jesup_Abs_year.Phi(),Jet2_2p5_jesup_Abs_year.Phi());
                    dPhiHj1j2_2p5_jesup_Abs_year=deltaPhi(Higgs.Phi(),(Jet1_2p5_jesup_Abs_year+Jet2_2p5_jesup_Abs_year).Phi());

                }

		pt_jesup_split_Abs_year.clear();  

/////////////////
// Abs_year dn start
                TauC_Inc_0j_EnergyWgt_jesdn_Abs_year = TauC(pt_jesdn_split_Abs_year, eta_jes_split_Abs_year,phi_jes_split_Abs_year,mass_jes_split_Abs_year, Higgs);
                TauB_Inc_0j_pTWgt_jesdn_Abs_year = TauB(pt_jesdn_split_Abs_year, eta_jes_split_Abs_year,phi_jes_split_Abs_year,mass_jes_split_Abs_year, Higgs);


                for( unsigned int k = 0; k< pt_jesdn_split_Abs_year.size(); k++) {
                    if (pt_jesdn_split_Abs_year[k]<30.0 || abs(eta_jes_split_Abs_year[k])>4.7) continue;
                    TLorentzVector thisJet_jesdn_Abs_year;
		    thisJet_jesdn_Abs_year.SetPtEtaPhiM(pt_jesdn_split_Abs_year[k],eta_jes_split_Abs_year[k],phi_jes_split_Abs_year[k],mass_jes_split_Abs_year[k]);

                    bool isclean_H4l=true;
                        njets_pt30_eta4p7_jesdn_Abs_year+=1;
			if (thisJet_jesdn_Abs_year.Pt()>jet1pt_jesdn_Abs_year) {
                            jet2pt_jesdn_Abs_year=jet1pt_jesdn_Abs_year; jet2index_jesdn_Abs_year=jet1index_jesdn_Abs_year;
                            jet1pt_jesdn_Abs_year=thisJet_jesdn_Abs_year.Pt(); jet1index_jesdn_Abs_year=k;
                        } else if (thisJet_jesdn_Abs_year.Pt()>jet2pt_jesdn_Abs_year) {
                            jet2pt_jesdn_Abs_year=thisJet_jesdn_Abs_year.Pt(); jet2index_jesdn_Abs_year=k;
                        }
                        if (abs(thisJet_jesdn_Abs_year.Eta())<2.5) {
                            njets_pt30_eta2p5_jesdn_Abs_year+=1;
                            if (thisJet_jesdn_Abs_year.Pt()>jet1pt2p5_jesdn_Abs_year) {
                                jet2pt2p5_jesdn_Abs_year=jet1pt2p5_jesdn_Abs_year; jet2index2p5_jesdn_Abs_year=jet1index2p5_jesdn_Abs_year;
                                jet1pt2p5_jesdn_Abs_year=thisJet_jesdn_Abs_year.Pt(); jet1index2p5_jesdn_Abs_year=k;
                            } else if (thisJet_jesdn_Abs_year.Pt()>jet2pt2p5_jesdn_Abs_year) {
                                jet2pt2p5_jesdn_Abs_year=thisJet_jesdn_Abs_year.Pt(); jet2index2p5_jesdn_Abs_year=k;
                            }
                        }
                }

// Filling Abs_year dn variables

                TLorentzVector Jet1_jesdn_Abs_year, Jet1_2p5_jesdn_Abs_year, Jet2_jesdn_Abs_year, Jet2_2p5_jesdn_Abs_year;

                if (njets_pt30_eta4p7_jesdn_Abs_year > 0) { 
                    Jet1_jesdn_Abs_year.SetPtEtaPhiM(pt_jesdn_split_Abs_year[jet1index_jesdn_Abs_year],eta_jes_split_Abs_year[jet1index_jesdn_Abs_year],phi_jes_split_Abs_year[jet1index_jesdn_Abs_year], mass_jes_split_Abs_year[jet1index_jesdn_Abs_year]);
                    pt_leadingjet_pt30_eta4p7_jesdn_Abs_year=Jet1_jesdn_Abs_year.Pt(); 
                    pTj1_jesdn_Abs_year=Jet1_jesdn_Abs_year.Pt(); 
                    etaj1_jesdn_Abs_year=Jet1_jesdn_Abs_year.Eta();
                    yj1_jesdn_Abs_year=Jet1_jesdn_Abs_year.Rapidity();
                    pT4lj_jesdn_Abs_year=(Higgs+Jet1_jesdn_Abs_year).Pt();
                    mass4lj_jesdn_Abs_year=(Higgs+Jet1_jesdn_Abs_year).M();
                    dPhiHj1_jesdn_Abs_year=deltaPhi(Higgs.Phi(),Jet1_jesdn_Abs_year.Phi());
                    dyHj1_jesdn_Abs_year=TMath::Abs(rapidity4l-yj1_jesdn_Abs_year);
                }    
                if (njets_pt30_eta4p7_jesdn_Abs_year > 1) { 
                    Jet2_jesdn_Abs_year.SetPtEtaPhiM(pt_jesdn_split_Abs_year[jet2index_jesdn_Abs_year],eta_jes_split_Abs_year[jet2index_jesdn_Abs_year],phi_jes_split_Abs_year[jet2index_jesdn_Abs_year], mass_jes_split_Abs_year[jet2index_jesdn_Abs_year]);

                    pTj2_jesdn_Abs_year=Jet2_jesdn_Abs_year.Pt();
                    etaj2_jesdn_Abs_year=Jet2_jesdn_Abs_year.Eta();
                    yj2_jesdn_Abs_year=Jet2_jesdn_Abs_year.Rapidity();
                    pT4ljj_jesdn_Abs_year=(Higgs+Jet1_jesdn_Abs_year+Jet2_jesdn_Abs_year).Pt();
                    mass4ljj_jesdn_Abs_year=(Higgs+Jet1_jesdn_Abs_year+Jet2_jesdn_Abs_year).M();
                    dEtaj1j2_jesdn_Abs_year=TMath::Abs(Jet1_jesdn_Abs_year.Eta()-Jet2_jesdn_Abs_year.Eta());
                    dPhij1j2_jesdn_Abs_year=deltaPhi(Jet1_jesdn_Abs_year.Phi(),Jet2_jesdn_Abs_year.Phi());
                    dPhiHj1j2_jesdn_Abs_year=deltaPhi(Higgs.Phi(),(Jet1_jesdn_Abs_year+Jet2_jesdn_Abs_year).Phi());
                }    
                if (njets_pt30_eta2p5_jesdn_Abs_year > 0) { 
                    Jet1_2p5_jesdn_Abs_year.SetPtEtaPhiM(pt_jesdn_split_Abs_year[jet1index2p5_jesdn_Abs_year],eta_jes_split_Abs_year[jet1index2p5_jesdn_Abs_year],phi_jes_split_Abs_year[jet1index2p5_jesdn_Abs_year], mass_jes_split_Abs_year[jet1index2p5_jesdn_Abs_year]);
                    pt_leadingjet_pt30_eta2p5_jesdn_Abs_year=Jet1_2p5_jesdn_Abs_year.Pt();
                    pTj1_2p5_jesdn_Abs_year=Jet1_2p5_jesdn_Abs_year.Pt();
                    yj1_2p5_jesdn_Abs_year=Jet1_2p5_jesdn_Abs_year.Rapidity();
                    pT4lj_2p5_jesdn_Abs_year=(Higgs+Jet1_2p5_jesdn_Abs_year).Pt();
                    mass4lj_2p5_jesdn_Abs_year=(Higgs+Jet1_2p5_jesdn_Abs_year).M();
		    dPhiHj1_2p5_jesdn_Abs_year=deltaPhi(Higgs.Phi(),Jet1_2p5_jesdn_Abs_year.Phi());
                    dyHj1_2p5_jesdn_Abs_year=TMath::Abs(rapidity4l-yj1_2p5_jesdn_Abs_year);
                }    
                if (njets_pt30_eta2p5_jesdn_Abs_year > 1) { 
                    Jet2_2p5_jesdn_Abs_year.SetPtEtaPhiM(pt_jesdn_split_Abs_year[jet2index2p5_jesdn_Abs_year],eta_jes_split_Abs_year[jet2index2p5_jesdn_Abs_year],phi_jes_split_Abs_year[jet2index2p5_jesdn_Abs_year], mass_jes_split_Abs_year[jet2index2p5_jesdn_Abs_year]);
                    pTj2_2p5_jesdn_Abs_year=Jet2_2p5_jesdn_Abs_year.Pt();
                    yj2_2p5_jesdn_Abs_year=Jet2_2p5_jesdn_Abs_year.Rapidity();
                    pT4ljj_2p5_jesdn_Abs_year=(Higgs+Jet1_2p5_jesdn_Abs_year+Jet2_2p5_jesdn_Abs_year).Pt();
                    mass4ljj_2p5_jesdn_Abs_year=(Higgs+Jet1_2p5_jesdn_Abs_year+Jet2_2p5_jesdn_Abs_year).M();
                    mj1j2_2p5_jesdn_Abs_year=(Jet1_2p5_jesdn_Abs_year+Jet2_2p5_jesdn_Abs_year).M();     
                    dEtaj1j2_2p5_jesdn_Abs_year=TMath::Abs(Jet1_2p5_jesdn_Abs_year.Eta()-Jet2_2p5_jesdn_Abs_year.Eta());
	            dPhij1j2_2p5_jesdn_Abs_year=deltaPhi(Jet1_2p5_jesdn_Abs_year.Phi(),Jet2_2p5_jesdn_Abs_year.Phi());
                    dPhiHj1j2_2p5_jesdn_Abs_year=deltaPhi(Higgs.Phi(),(Jet1_2p5_jesdn_Abs_year+Jet2_2p5_jesdn_Abs_year).Phi());

                }    

                pt_jesdn_split_Abs_year.clear();

		eta_jes_split_Abs_year.clear();
		phi_jes_split_Abs_year.clear();
		mass_jes_split_Abs_year.clear();


////////////////////////////////////////////////////////////////////////////////
////    Starting "BBEC1" JES uncertainty source
////////////////////////////////////////////////////////////////////////////////

                TauC_Inc_0j_EnergyWgt_jesup_BBEC1 = TauC(pt_jesup_split_BBEC1, eta_jes_split_BBEC1,phi_jes_split_BBEC1,mass_jes_split_BBEC1, Higgs);
                TauB_Inc_0j_pTWgt_jesup_BBEC1 = TauB(pt_jesup_split_BBEC1, eta_jes_split_BBEC1,phi_jes_split_BBEC1,mass_jes_split_BBEC1, Higgs);

                for( unsigned int k = 0; k< pt_jesup_split_BBEC1.size(); k++) {
                    if (pt_jesup_split_BBEC1[k]<30.0 || abs(eta_jes_split_BBEC1[k])>4.7) continue;
		    TLorentzVector thisJet_jesup_BBEC1;
		    thisJet_jesup_BBEC1.SetPtEtaPhiM(pt_jesup_split_BBEC1[k],eta_jes_split_BBEC1[k],phi_jes_split_BBEC1[k],mass_jes_split_BBEC1[k]);

                        njets_pt30_eta4p7_jesup_BBEC1+=1;  

                        if (thisJet_jesup_BBEC1.Pt()>jet1pt_jesup_BBEC1) {
                            jet2pt_jesup_BBEC1=jet1pt_jesup_BBEC1; jet2index_jesup_BBEC1=jet1index_jesup_BBEC1;
                            jet1pt_jesup_BBEC1=thisJet_jesup_BBEC1.Pt(); jet1index_jesup_BBEC1=k;
                        } else if (thisJet_jesup_BBEC1.Pt()>jet2pt_jesup_BBEC1) {
                            jet2pt_jesup_BBEC1=thisJet_jesup_BBEC1.Pt(); jet2index_jesup_BBEC1=k;
                        }
                        if (abs(thisJet_jesup_BBEC1.Eta())<2.5) {
                            njets_pt30_eta2p5_jesup_BBEC1+=1;
                            if (thisJet_jesup_BBEC1.Pt()>jet1pt2p5_jesup_BBEC1) {
                                jet2pt2p5_jesup_BBEC1=jet1pt2p5_jesup_BBEC1; jet2index2p5_jesup_BBEC1=jet1index2p5_jesup_BBEC1;
                                jet1pt2p5_jesup_BBEC1=thisJet_jesup_BBEC1.Pt(); jet1index2p5_jesup_BBEC1=k;
                            } else if (thisJet_jesup_BBEC1.Pt()>jet2pt2p5_jesup_BBEC1) {
                                jet2pt2p5_jesup_BBEC1=thisJet_jesup_BBEC1.Pt(); jet2index2p5_jesup_BBEC1=k;
                            }
                        }
                }
                cout<<njets_pt30_eta4p7_jesup_BBEC1<<" jets (jesup_BBEC1)"<<endl;


		TLorentzVector Jet1_jesup_BBEC1, Jet1_2p5_jesup_BBEC1, Jet2_jesup_BBEC1, Jet2_2p5_jesup_BBEC1;
                if (njets_pt30_eta4p7_jesup_BBEC1 > 0) {
		    Jet1_jesup_BBEC1.SetPtEtaPhiM(pt_jesup_split_BBEC1[jet1index_jesup_BBEC1],eta_jes_split_BBEC1[jet1index_jesup_BBEC1],phi_jes_split_BBEC1[jet1index_jesup_BBEC1], mass_jes_split_BBEC1[jet1index_jesup_BBEC1]);

                    pt_leadingjet_pt30_eta4p7_jesup_BBEC1=Jet1_jesup_BBEC1.Pt(); 
                    pTj1_jesup_BBEC1=Jet1_jesup_BBEC1.Pt(); 
                    etaj1_jesup_BBEC1=Jet1_jesup_BBEC1.Eta();
                    yj1_jesup_BBEC1=Jet1_jesup_BBEC1.Rapidity();
		    pT4lj_jesup_BBEC1=(Higgs+Jet1_jesup_BBEC1).Pt();
		    mass4lj_jesup_BBEC1=(Higgs+Jet1_jesup_BBEC1).M();
		    dPhiHj1_jesup_BBEC1=deltaPhi(Higgs.Phi(),Jet1_jesup_BBEC1.Phi());
                    dyHj1_jesup_BBEC1=TMath::Abs(rapidity4l-yj1_jesup_BBEC1);
                }
                if (njets_pt30_eta4p7_jesup_BBEC1 > 1) {
                    Jet2_jesup_BBEC1.SetPtEtaPhiM(pt_jesup_split_BBEC1[jet2index_jesup_BBEC1],eta_jes_split_BBEC1[jet2index_jesup_BBEC1],phi_jes_split_BBEC1[jet2index_jesup_BBEC1], mass_jes_split_BBEC1[jet2index_jesup_BBEC1]);

                    pTj2_jesup_BBEC1=Jet2_jesup_BBEC1.Pt();
                    etaj2_jesup_BBEC1=Jet2_jesup_BBEC1.Eta();
                    yj2_jesup_BBEC1=Jet2_jesup_BBEC1.Rapidity();
                    pT4ljj_jesup_BBEC1=(Higgs+Jet1_jesup_BBEC1+Jet2_jesup_BBEC1).Pt();
                    mass4ljj_jesup_BBEC1=(Higgs+Jet1_jesup_BBEC1+Jet2_jesup_BBEC1).M();
		    mj1j2_jesup_BBEC1=(Jet1_jesup_BBEC1+Jet2_jesup_BBEC1).M();
                    dEtaj1j2_jesup_BBEC1=TMath::Abs(Jet1_jesup_BBEC1.Eta()-Jet2_jesup_BBEC1.Eta());
                    dPhij1j2_jesup_BBEC1=deltaPhi(Jet1_jesup_BBEC1.Phi(),Jet2_jesup_BBEC1.Phi());
                    dPhiHj1j2_jesup_BBEC1=deltaPhi(Higgs.Phi(),(Jet1_jesup_BBEC1+Jet2_jesup_BBEC1).Phi());
                }
                if (njets_pt30_eta2p5_jesup_BBEC1 > 0) {
                    Jet1_2p5_jesup_BBEC1.SetPtEtaPhiM(pt_jesup_split_BBEC1[jet1index2p5_jesup_BBEC1],eta_jes_split_BBEC1[jet1index2p5_jesup_BBEC1],phi_jes_split_BBEC1[jet1index2p5_jesup_BBEC1], mass_jes_split_BBEC1[jet1index2p5_jesup_BBEC1]);
                    pt_leadingjet_pt30_eta2p5_jesup_BBEC1=Jet1_2p5_jesup_BBEC1.Pt();
                    pTj1_2p5_jesup_BBEC1=Jet1_2p5_jesup_BBEC1.Pt();
                    yj1_2p5_jesup_BBEC1=Jet1_2p5_jesup_BBEC1.Rapidity();
                    pT4lj_2p5_jesup_BBEC1=(Higgs+Jet1_2p5_jesup_BBEC1).Pt();
                    mass4lj_2p5_jesup_BBEC1=(Higgs+Jet1_2p5_jesup_BBEC1).M();
		    dPhiHj1_2p5_jesup_BBEC1=deltaPhi(Higgs.Phi(),Jet1_2p5_jesup_BBEC1.Phi());
                    dyHj1_2p5_jesup_BBEC1=TMath::Abs(rapidity4l-yj1_2p5_jesup_BBEC1);
                }
                if (njets_pt30_eta2p5_jesup_BBEC1 > 1) {
                    Jet2_2p5_jesup_BBEC1.SetPtEtaPhiM(pt_jesup_split_BBEC1[jet2index2p5_jesup_BBEC1],eta_jes_split_BBEC1[jet2index2p5_jesup_BBEC1],phi_jes_split_BBEC1[jet2index2p5_jesup_BBEC1], mass_jes_split_BBEC1[jet2index2p5_jesup_BBEC1]);
                    pTj2_2p5_jesup_BBEC1=Jet2_2p5_jesup_BBEC1.Pt();
                    yj2_2p5_jesup_BBEC1=Jet2_2p5_jesup_BBEC1.Rapidity();
                    pT4ljj_2p5_jesup_BBEC1=(Higgs+Jet1_2p5_jesup_BBEC1+Jet2_2p5_jesup_BBEC1).Pt();
                    mass4ljj_2p5_jesup_BBEC1=(Higgs+Jet1_2p5_jesup_BBEC1+Jet2_2p5_jesup_BBEC1).M();
                    mj1j2_2p5_jesup_BBEC1=(Jet1_2p5_jesup_BBEC1+Jet2_2p5_jesup_BBEC1).M();                    
                    dEtaj1j2_2p5_jesup_BBEC1=TMath::Abs(Jet1_2p5_jesup_BBEC1.Eta()-Jet2_2p5_jesup_BBEC1.Eta());
                    dPhij1j2_2p5_jesup_BBEC1=deltaPhi(Jet1_2p5_jesup_BBEC1.Phi(),Jet2_2p5_jesup_BBEC1.Phi());
                    dPhiHj1j2_2p5_jesup_BBEC1=deltaPhi(Higgs.Phi(),(Jet1_2p5_jesup_BBEC1+Jet2_2p5_jesup_BBEC1).Phi());

                }

		pt_jesup_split_BBEC1.clear();  

/////////////////
// BBEC1 dn start
                TauC_Inc_0j_EnergyWgt_jesdn_BBEC1 = TauC(pt_jesdn_split_BBEC1, eta_jes_split_BBEC1,phi_jes_split_BBEC1,mass_jes_split_BBEC1, Higgs);
                TauB_Inc_0j_pTWgt_jesdn_BBEC1 = TauB(pt_jesdn_split_BBEC1, eta_jes_split_BBEC1,phi_jes_split_BBEC1,mass_jes_split_BBEC1, Higgs);


                for( unsigned int k = 0; k< pt_jesdn_split_BBEC1.size(); k++) {
                    if (pt_jesdn_split_BBEC1[k]<30.0 || abs(eta_jes_split_BBEC1[k])>4.7) continue;
                    TLorentzVector thisJet_jesdn_BBEC1;
		    thisJet_jesdn_BBEC1.SetPtEtaPhiM(pt_jesdn_split_BBEC1[k],eta_jes_split_BBEC1[k],phi_jes_split_BBEC1[k],mass_jes_split_BBEC1[k]);

                    bool isclean_H4l=true;
                        njets_pt30_eta4p7_jesdn_BBEC1+=1;
			if (thisJet_jesdn_BBEC1.Pt()>jet1pt_jesdn_BBEC1) {
                            jet2pt_jesdn_BBEC1=jet1pt_jesdn_BBEC1; jet2index_jesdn_BBEC1=jet1index_jesdn_BBEC1;
                            jet1pt_jesdn_BBEC1=thisJet_jesdn_BBEC1.Pt(); jet1index_jesdn_BBEC1=k;
                        } else if (thisJet_jesdn_BBEC1.Pt()>jet2pt_jesdn_BBEC1) {
                            jet2pt_jesdn_BBEC1=thisJet_jesdn_BBEC1.Pt(); jet2index_jesdn_BBEC1=k;
                        }
                        if (abs(thisJet_jesdn_BBEC1.Eta())<2.5) {
                            njets_pt30_eta2p5_jesdn_BBEC1+=1;
                            if (thisJet_jesdn_BBEC1.Pt()>jet1pt2p5_jesdn_BBEC1) {
                                jet2pt2p5_jesdn_BBEC1=jet1pt2p5_jesdn_BBEC1; jet2index2p5_jesdn_BBEC1=jet1index2p5_jesdn_BBEC1;
                                jet1pt2p5_jesdn_BBEC1=thisJet_jesdn_BBEC1.Pt(); jet1index2p5_jesdn_BBEC1=k;
                            } else if (thisJet_jesdn_BBEC1.Pt()>jet2pt2p5_jesdn_BBEC1) {
                                jet2pt2p5_jesdn_BBEC1=thisJet_jesdn_BBEC1.Pt(); jet2index2p5_jesdn_BBEC1=k;
                            }
                        }
                }

// Filling BBEC1 dn variables

                TLorentzVector Jet1_jesdn_BBEC1, Jet1_2p5_jesdn_BBEC1, Jet2_jesdn_BBEC1, Jet2_2p5_jesdn_BBEC1;

                if (njets_pt30_eta4p7_jesdn_BBEC1 > 0) { 
                    Jet1_jesdn_BBEC1.SetPtEtaPhiM(pt_jesdn_split_BBEC1[jet1index_jesdn_BBEC1],eta_jes_split_BBEC1[jet1index_jesdn_BBEC1],phi_jes_split_BBEC1[jet1index_jesdn_BBEC1], mass_jes_split_BBEC1[jet1index_jesdn_BBEC1]);
                    pt_leadingjet_pt30_eta4p7_jesdn_BBEC1=Jet1_jesdn_BBEC1.Pt(); 
                    pTj1_jesdn_BBEC1=Jet1_jesdn_BBEC1.Pt(); 
                    etaj1_jesdn_BBEC1=Jet1_jesdn_BBEC1.Eta();
                    yj1_jesdn_BBEC1=Jet1_jesdn_BBEC1.Rapidity();
                    pT4lj_jesdn_BBEC1=(Higgs+Jet1_jesdn_BBEC1).Pt();
                    mass4lj_jesdn_BBEC1=(Higgs+Jet1_jesdn_BBEC1).M();
                    dPhiHj1_jesdn_BBEC1=deltaPhi(Higgs.Phi(),Jet1_jesdn_BBEC1.Phi());
                    dyHj1_jesdn_BBEC1=TMath::Abs(rapidity4l-yj1_jesdn_BBEC1);
                }    
                if (njets_pt30_eta4p7_jesdn_BBEC1 > 1) { 
                    Jet2_jesdn_BBEC1.SetPtEtaPhiM(pt_jesdn_split_BBEC1[jet2index_jesdn_BBEC1],eta_jes_split_BBEC1[jet2index_jesdn_BBEC1],phi_jes_split_BBEC1[jet2index_jesdn_BBEC1], mass_jes_split_BBEC1[jet2index_jesdn_BBEC1]);

                    pTj2_jesdn_BBEC1=Jet2_jesdn_BBEC1.Pt();
                    etaj2_jesdn_BBEC1=Jet2_jesdn_BBEC1.Eta();
                    yj2_jesdn_BBEC1=Jet2_jesdn_BBEC1.Rapidity();
                    pT4ljj_jesdn_BBEC1=(Higgs+Jet1_jesdn_BBEC1+Jet2_jesdn_BBEC1).Pt();
                    mass4ljj_jesdn_BBEC1=(Higgs+Jet1_jesdn_BBEC1+Jet2_jesdn_BBEC1).M();
                    dEtaj1j2_jesdn_BBEC1=TMath::Abs(Jet1_jesdn_BBEC1.Eta()-Jet2_jesdn_BBEC1.Eta());
                    dPhij1j2_jesdn_BBEC1=deltaPhi(Jet1_jesdn_BBEC1.Phi(),Jet2_jesdn_BBEC1.Phi());
                    dPhiHj1j2_jesdn_BBEC1=deltaPhi(Higgs.Phi(),(Jet1_jesdn_BBEC1+Jet2_jesdn_BBEC1).Phi());
                }    
                if (njets_pt30_eta2p5_jesdn_BBEC1 > 0) { 
                    Jet1_2p5_jesdn_BBEC1.SetPtEtaPhiM(pt_jesdn_split_BBEC1[jet1index2p5_jesdn_BBEC1],eta_jes_split_BBEC1[jet1index2p5_jesdn_BBEC1],phi_jes_split_BBEC1[jet1index2p5_jesdn_BBEC1], mass_jes_split_BBEC1[jet1index2p5_jesdn_BBEC1]);
                    pt_leadingjet_pt30_eta2p5_jesdn_BBEC1=Jet1_2p5_jesdn_BBEC1.Pt();
                    pTj1_2p5_jesdn_BBEC1=Jet1_2p5_jesdn_BBEC1.Pt();
                    yj1_2p5_jesdn_BBEC1=Jet1_2p5_jesdn_BBEC1.Rapidity();
                    pT4lj_2p5_jesdn_BBEC1=(Higgs+Jet1_2p5_jesdn_BBEC1).Pt();
                    mass4lj_2p5_jesdn_BBEC1=(Higgs+Jet1_2p5_jesdn_BBEC1).M();
		    dPhiHj1_2p5_jesdn_BBEC1=deltaPhi(Higgs.Phi(),Jet1_2p5_jesdn_BBEC1.Phi());
                    dyHj1_2p5_jesdn_BBEC1=TMath::Abs(rapidity4l-yj1_2p5_jesdn_BBEC1);
                }    
                if (njets_pt30_eta2p5_jesdn_BBEC1 > 1) { 
                    Jet2_2p5_jesdn_BBEC1.SetPtEtaPhiM(pt_jesdn_split_BBEC1[jet2index2p5_jesdn_BBEC1],eta_jes_split_BBEC1[jet2index2p5_jesdn_BBEC1],phi_jes_split_BBEC1[jet2index2p5_jesdn_BBEC1], mass_jes_split_BBEC1[jet2index2p5_jesdn_BBEC1]);
                    pTj2_2p5_jesdn_BBEC1=Jet2_2p5_jesdn_BBEC1.Pt();
                    yj2_2p5_jesdn_BBEC1=Jet2_2p5_jesdn_BBEC1.Rapidity();
                    pT4ljj_2p5_jesdn_BBEC1=(Higgs+Jet1_2p5_jesdn_BBEC1+Jet2_2p5_jesdn_BBEC1).Pt();
                    mass4ljj_2p5_jesdn_BBEC1=(Higgs+Jet1_2p5_jesdn_BBEC1+Jet2_2p5_jesdn_BBEC1).M();
                    mj1j2_2p5_jesdn_BBEC1=(Jet1_2p5_jesdn_BBEC1+Jet2_2p5_jesdn_BBEC1).M();     
                    dEtaj1j2_2p5_jesdn_BBEC1=TMath::Abs(Jet1_2p5_jesdn_BBEC1.Eta()-Jet2_2p5_jesdn_BBEC1.Eta());
	            dPhij1j2_2p5_jesdn_BBEC1=deltaPhi(Jet1_2p5_jesdn_BBEC1.Phi(),Jet2_2p5_jesdn_BBEC1.Phi());
                    dPhiHj1j2_2p5_jesdn_BBEC1=deltaPhi(Higgs.Phi(),(Jet1_2p5_jesdn_BBEC1+Jet2_2p5_jesdn_BBEC1).Phi());

                }    

                pt_jesdn_split_BBEC1.clear();

		eta_jes_split_BBEC1.clear();
		phi_jes_split_BBEC1.clear();
		mass_jes_split_BBEC1.clear();


////////////////////////////////////////////////////////////////////////////////
////    Starting "BBEC1_year" JES uncertainty source
////////////////////////////////////////////////////////////////////////////////

                TauC_Inc_0j_EnergyWgt_jesup_BBEC1_year = TauC(pt_jesup_split_BBEC1_year, eta_jes_split_BBEC1_year,phi_jes_split_BBEC1_year,mass_jes_split_BBEC1_year, Higgs);
                TauB_Inc_0j_pTWgt_jesup_BBEC1_year = TauB(pt_jesup_split_BBEC1_year, eta_jes_split_BBEC1_year,phi_jes_split_BBEC1_year,mass_jes_split_BBEC1_year, Higgs);


                for( unsigned int k = 0; k< pt_jesup_split_BBEC1_year.size(); k++) {
                    if (pt_jesup_split_BBEC1_year[k]<30.0 || abs(eta_jes_split_BBEC1_year[k])>4.7) continue;
		    TLorentzVector thisJet_jesup_BBEC1_year;
		    thisJet_jesup_BBEC1_year.SetPtEtaPhiM(pt_jesup_split_BBEC1_year[k],eta_jes_split_BBEC1_year[k],phi_jes_split_BBEC1_year[k],mass_jes_split_BBEC1_year[k]);

                        njets_pt30_eta4p7_jesup_BBEC1_year+=1;  

                        if (thisJet_jesup_BBEC1_year.Pt()>jet1pt_jesup_BBEC1_year) {
                            jet2pt_jesup_BBEC1_year=jet1pt_jesup_BBEC1_year; jet2index_jesup_BBEC1_year=jet1index_jesup_BBEC1_year;
                            jet1pt_jesup_BBEC1_year=thisJet_jesup_BBEC1_year.Pt(); jet1index_jesup_BBEC1_year=k;
                        } else if (thisJet_jesup_BBEC1_year.Pt()>jet2pt_jesup_BBEC1_year) {
                            jet2pt_jesup_BBEC1_year=thisJet_jesup_BBEC1_year.Pt(); jet2index_jesup_BBEC1_year=k;
                        }
                        if (abs(thisJet_jesup_BBEC1_year.Eta())<2.5) {
                            njets_pt30_eta2p5_jesup_BBEC1_year+=1;
                            if (thisJet_jesup_BBEC1_year.Pt()>jet1pt2p5_jesup_BBEC1_year) {
                                jet2pt2p5_jesup_BBEC1_year=jet1pt2p5_jesup_BBEC1_year; jet2index2p5_jesup_BBEC1_year=jet1index2p5_jesup_BBEC1_year;
                                jet1pt2p5_jesup_BBEC1_year=thisJet_jesup_BBEC1_year.Pt(); jet1index2p5_jesup_BBEC1_year=k;
                            } else if (thisJet_jesup_BBEC1_year.Pt()>jet2pt2p5_jesup_BBEC1_year) {
                                jet2pt2p5_jesup_BBEC1_year=thisJet_jesup_BBEC1_year.Pt(); jet2index2p5_jesup_BBEC1_year=k;
                            }
                        }
                }
                cout<<njets_pt30_eta4p7_jesup_BBEC1_year<<" jets (jesup_BBEC1_year)"<<endl;


		TLorentzVector Jet1_jesup_BBEC1_year, Jet1_2p5_jesup_BBEC1_year, Jet2_jesup_BBEC1_year, Jet2_2p5_jesup_BBEC1_year;
                if (njets_pt30_eta4p7_jesup_BBEC1_year > 0) {
		    Jet1_jesup_BBEC1_year.SetPtEtaPhiM(pt_jesup_split_BBEC1_year[jet1index_jesup_BBEC1_year],eta_jes_split_BBEC1_year[jet1index_jesup_BBEC1_year],phi_jes_split_BBEC1_year[jet1index_jesup_BBEC1_year], mass_jes_split_BBEC1_year[jet1index_jesup_BBEC1_year]);

                    pt_leadingjet_pt30_eta4p7_jesup_BBEC1_year=Jet1_jesup_BBEC1_year.Pt(); 
                    pTj1_jesup_BBEC1_year=Jet1_jesup_BBEC1_year.Pt(); 
                    etaj1_jesup_BBEC1_year=Jet1_jesup_BBEC1_year.Eta();
                    yj1_jesup_BBEC1_year=Jet1_jesup_BBEC1_year.Rapidity();
		    pT4lj_jesup_BBEC1_year=(Higgs+Jet1_jesup_BBEC1_year).Pt();
		    mass4lj_jesup_BBEC1_year=(Higgs+Jet1_jesup_BBEC1_year).M();
		    dPhiHj1_jesup_BBEC1_year=deltaPhi(Higgs.Phi(),Jet1_jesup_BBEC1_year.Phi());
                    dyHj1_jesup_BBEC1_year=TMath::Abs(rapidity4l-yj1_jesup_BBEC1_year);
                }
                if (njets_pt30_eta4p7_jesup_BBEC1_year > 1) {
                    Jet2_jesup_BBEC1_year.SetPtEtaPhiM(pt_jesup_split_BBEC1_year[jet2index_jesup_BBEC1_year],eta_jes_split_BBEC1_year[jet2index_jesup_BBEC1_year],phi_jes_split_BBEC1_year[jet2index_jesup_BBEC1_year], mass_jes_split_BBEC1_year[jet2index_jesup_BBEC1_year]);

                    pTj2_jesup_BBEC1_year=Jet2_jesup_BBEC1_year.Pt();
                    etaj2_jesup_BBEC1_year=Jet2_jesup_BBEC1_year.Eta();
                    yj2_jesup_BBEC1_year=Jet2_jesup_BBEC1_year.Rapidity();
                    pT4ljj_jesup_BBEC1_year=(Higgs+Jet1_jesup_BBEC1_year+Jet2_jesup_BBEC1_year).Pt();
                    mass4ljj_jesup_BBEC1_year=(Higgs+Jet1_jesup_BBEC1_year+Jet2_jesup_BBEC1_year).M();
		    mj1j2_jesup_BBEC1_year=(Jet1_jesup_BBEC1_year+Jet2_jesup_BBEC1_year).M();
                    dEtaj1j2_jesup_BBEC1_year=TMath::Abs(Jet1_jesup_BBEC1_year.Eta()-Jet2_jesup_BBEC1_year.Eta());
                    dPhij1j2_jesup_BBEC1_year=deltaPhi(Jet1_jesup_BBEC1_year.Phi(),Jet2_jesup_BBEC1_year.Phi());
                    dPhiHj1j2_jesup_BBEC1_year=deltaPhi(Higgs.Phi(),(Jet1_jesup_BBEC1_year+Jet2_jesup_BBEC1_year).Phi());
                }
                if (njets_pt30_eta2p5_jesup_BBEC1_year > 0) {
                    Jet1_2p5_jesup_BBEC1_year.SetPtEtaPhiM(pt_jesup_split_BBEC1_year[jet1index2p5_jesup_BBEC1_year],eta_jes_split_BBEC1_year[jet1index2p5_jesup_BBEC1_year],phi_jes_split_BBEC1_year[jet1index2p5_jesup_BBEC1_year], mass_jes_split_BBEC1_year[jet1index2p5_jesup_BBEC1_year]);
                    pt_leadingjet_pt30_eta2p5_jesup_BBEC1_year=Jet1_2p5_jesup_BBEC1_year.Pt();
                    pTj1_2p5_jesup_BBEC1_year=Jet1_2p5_jesup_BBEC1_year.Pt();
                    yj1_2p5_jesup_BBEC1_year=Jet1_2p5_jesup_BBEC1_year.Rapidity();
                    pT4lj_2p5_jesup_BBEC1_year=(Higgs+Jet1_2p5_jesup_BBEC1_year).Pt();
                    mass4lj_2p5_jesup_BBEC1_year=(Higgs+Jet1_2p5_jesup_BBEC1_year).M();
		    dPhiHj1_2p5_jesup_BBEC1_year=deltaPhi(Higgs.Phi(),Jet1_2p5_jesup_BBEC1_year.Phi());
                    dyHj1_2p5_jesup_BBEC1_year=TMath::Abs(rapidity4l-yj1_2p5_jesup_BBEC1_year);
                }
                if (njets_pt30_eta2p5_jesup_BBEC1_year > 1) {
                    Jet2_2p5_jesup_BBEC1_year.SetPtEtaPhiM(pt_jesup_split_BBEC1_year[jet2index2p5_jesup_BBEC1_year],eta_jes_split_BBEC1_year[jet2index2p5_jesup_BBEC1_year],phi_jes_split_BBEC1_year[jet2index2p5_jesup_BBEC1_year], mass_jes_split_BBEC1_year[jet2index2p5_jesup_BBEC1_year]);
                    pTj2_2p5_jesup_BBEC1_year=Jet2_2p5_jesup_BBEC1_year.Pt();
                    yj2_2p5_jesup_BBEC1_year=Jet2_2p5_jesup_BBEC1_year.Rapidity();
                    pT4ljj_2p5_jesup_BBEC1_year=(Higgs+Jet1_2p5_jesup_BBEC1_year+Jet2_2p5_jesup_BBEC1_year).Pt();
                    mass4ljj_2p5_jesup_BBEC1_year=(Higgs+Jet1_2p5_jesup_BBEC1_year+Jet2_2p5_jesup_BBEC1_year).M();
                    mj1j2_2p5_jesup_BBEC1_year=(Jet1_2p5_jesup_BBEC1_year+Jet2_2p5_jesup_BBEC1_year).M();                    
                    dEtaj1j2_2p5_jesup_BBEC1_year=TMath::Abs(Jet1_2p5_jesup_BBEC1_year.Eta()-Jet2_2p5_jesup_BBEC1_year.Eta());
                    dPhij1j2_2p5_jesup_BBEC1_year=deltaPhi(Jet1_2p5_jesup_BBEC1_year.Phi(),Jet2_2p5_jesup_BBEC1_year.Phi());
                    dPhiHj1j2_2p5_jesup_BBEC1_year=deltaPhi(Higgs.Phi(),(Jet1_2p5_jesup_BBEC1_year+Jet2_2p5_jesup_BBEC1_year).Phi());

                }

		pt_jesup_split_BBEC1_year.clear();  

/////////////////
// BBEC1_year dn start
                TauC_Inc_0j_EnergyWgt_jesdn_BBEC1_year = TauC(pt_jesdn_split_BBEC1_year, eta_jes_split_BBEC1_year,phi_jes_split_BBEC1_year,mass_jes_split_BBEC1_year, Higgs);
                TauB_Inc_0j_pTWgt_jesdn_BBEC1_year = TauB(pt_jesdn_split_BBEC1_year, eta_jes_split_BBEC1_year,phi_jes_split_BBEC1_year,mass_jes_split_BBEC1_year, Higgs);


                for( unsigned int k = 0; k< pt_jesdn_split_BBEC1_year.size(); k++) {
                    if (pt_jesdn_split_BBEC1_year[k]<30.0 || abs(eta_jes_split_BBEC1_year[k])>4.7) continue;
                    TLorentzVector thisJet_jesdn_BBEC1_year;
		    thisJet_jesdn_BBEC1_year.SetPtEtaPhiM(pt_jesdn_split_BBEC1_year[k],eta_jes_split_BBEC1_year[k],phi_jes_split_BBEC1_year[k],mass_jes_split_BBEC1_year[k]);

                    bool isclean_H4l=true;
                        njets_pt30_eta4p7_jesdn_BBEC1_year+=1;
			if (thisJet_jesdn_BBEC1_year.Pt()>jet1pt_jesdn_BBEC1_year) {
                            jet2pt_jesdn_BBEC1_year=jet1pt_jesdn_BBEC1_year; jet2index_jesdn_BBEC1_year=jet1index_jesdn_BBEC1_year;
                            jet1pt_jesdn_BBEC1_year=thisJet_jesdn_BBEC1_year.Pt(); jet1index_jesdn_BBEC1_year=k;
                        } else if (thisJet_jesdn_BBEC1_year.Pt()>jet2pt_jesdn_BBEC1_year) {
                            jet2pt_jesdn_BBEC1_year=thisJet_jesdn_BBEC1_year.Pt(); jet2index_jesdn_BBEC1_year=k;
                        }
                        if (abs(thisJet_jesdn_BBEC1_year.Eta())<2.5) {
                            njets_pt30_eta2p5_jesdn_BBEC1_year+=1;
                            if (thisJet_jesdn_BBEC1_year.Pt()>jet1pt2p5_jesdn_BBEC1_year) {
                                jet2pt2p5_jesdn_BBEC1_year=jet1pt2p5_jesdn_BBEC1_year; jet2index2p5_jesdn_BBEC1_year=jet1index2p5_jesdn_BBEC1_year;
                                jet1pt2p5_jesdn_BBEC1_year=thisJet_jesdn_BBEC1_year.Pt(); jet1index2p5_jesdn_BBEC1_year=k;
                            } else if (thisJet_jesdn_BBEC1_year.Pt()>jet2pt2p5_jesdn_BBEC1_year) {
                                jet2pt2p5_jesdn_BBEC1_year=thisJet_jesdn_BBEC1_year.Pt(); jet2index2p5_jesdn_BBEC1_year=k;
                            }
                        }
                }

// Filling BBEC1_year dn variables

                TLorentzVector Jet1_jesdn_BBEC1_year, Jet1_2p5_jesdn_BBEC1_year, Jet2_jesdn_BBEC1_year, Jet2_2p5_jesdn_BBEC1_year;

                if (njets_pt30_eta4p7_jesdn_BBEC1_year > 0) { 
                    Jet1_jesdn_BBEC1_year.SetPtEtaPhiM(pt_jesdn_split_BBEC1_year[jet1index_jesdn_BBEC1_year],eta_jes_split_BBEC1_year[jet1index_jesdn_BBEC1_year],phi_jes_split_BBEC1_year[jet1index_jesdn_BBEC1_year], mass_jes_split_BBEC1_year[jet1index_jesdn_BBEC1_year]);
                    pt_leadingjet_pt30_eta4p7_jesdn_BBEC1_year=Jet1_jesdn_BBEC1_year.Pt(); 
                    pTj1_jesdn_BBEC1_year=Jet1_jesdn_BBEC1_year.Pt(); 
                    etaj1_jesdn_BBEC1_year=Jet1_jesdn_BBEC1_year.Eta();
                    yj1_jesdn_BBEC1_year=Jet1_jesdn_BBEC1_year.Rapidity();
                    pT4lj_jesdn_BBEC1_year=(Higgs+Jet1_jesdn_BBEC1_year).Pt();
                    mass4lj_jesdn_BBEC1_year=(Higgs+Jet1_jesdn_BBEC1_year).M();
                    dPhiHj1_jesdn_BBEC1_year=deltaPhi(Higgs.Phi(),Jet1_jesdn_BBEC1_year.Phi());
                    dyHj1_jesdn_BBEC1_year=TMath::Abs(rapidity4l-yj1_jesdn_BBEC1_year);
                }    
                if (njets_pt30_eta4p7_jesdn_BBEC1_year > 1) { 
                    Jet2_jesdn_BBEC1_year.SetPtEtaPhiM(pt_jesdn_split_BBEC1_year[jet2index_jesdn_BBEC1_year],eta_jes_split_BBEC1_year[jet2index_jesdn_BBEC1_year],phi_jes_split_BBEC1_year[jet2index_jesdn_BBEC1_year], mass_jes_split_BBEC1_year[jet2index_jesdn_BBEC1_year]);

                    pTj2_jesdn_BBEC1_year=Jet2_jesdn_BBEC1_year.Pt();
                    etaj2_jesdn_BBEC1_year=Jet2_jesdn_BBEC1_year.Eta();
                    yj2_jesdn_BBEC1_year=Jet2_jesdn_BBEC1_year.Rapidity();
                    pT4ljj_jesdn_BBEC1_year=(Higgs+Jet1_jesdn_BBEC1_year+Jet2_jesdn_BBEC1_year).Pt();
                    mass4ljj_jesdn_BBEC1_year=(Higgs+Jet1_jesdn_BBEC1_year+Jet2_jesdn_BBEC1_year).M();
                    dEtaj1j2_jesdn_BBEC1_year=TMath::Abs(Jet1_jesdn_BBEC1_year.Eta()-Jet2_jesdn_BBEC1_year.Eta());
                    dPhij1j2_jesdn_BBEC1_year=deltaPhi(Jet1_jesdn_BBEC1_year.Phi(),Jet2_jesdn_BBEC1_year.Phi());
                    dPhiHj1j2_jesdn_BBEC1_year=deltaPhi(Higgs.Phi(),(Jet1_jesdn_BBEC1_year+Jet2_jesdn_BBEC1_year).Phi());
                }    
                if (njets_pt30_eta2p5_jesdn_BBEC1_year > 0) { 
                    Jet1_2p5_jesdn_BBEC1_year.SetPtEtaPhiM(pt_jesdn_split_BBEC1_year[jet1index2p5_jesdn_BBEC1_year],eta_jes_split_BBEC1_year[jet1index2p5_jesdn_BBEC1_year],phi_jes_split_BBEC1_year[jet1index2p5_jesdn_BBEC1_year], mass_jes_split_BBEC1_year[jet1index2p5_jesdn_BBEC1_year]);
                    pt_leadingjet_pt30_eta2p5_jesdn_BBEC1_year=Jet1_2p5_jesdn_BBEC1_year.Pt();
                    pTj1_2p5_jesdn_BBEC1_year=Jet1_2p5_jesdn_BBEC1_year.Pt();
                    yj1_2p5_jesdn_BBEC1_year=Jet1_2p5_jesdn_BBEC1_year.Rapidity();
                    pT4lj_2p5_jesdn_BBEC1_year=(Higgs+Jet1_2p5_jesdn_BBEC1_year).Pt();
                    mass4lj_2p5_jesdn_BBEC1_year=(Higgs+Jet1_2p5_jesdn_BBEC1_year).M();
		    dPhiHj1_2p5_jesdn_BBEC1_year=deltaPhi(Higgs.Phi(),Jet1_2p5_jesdn_BBEC1_year.Phi());
                    dyHj1_2p5_jesdn_BBEC1_year=TMath::Abs(rapidity4l-yj1_2p5_jesdn_BBEC1_year);
                }    
                if (njets_pt30_eta2p5_jesdn_BBEC1_year > 1) { 
                    Jet2_2p5_jesdn_BBEC1_year.SetPtEtaPhiM(pt_jesdn_split_BBEC1_year[jet2index2p5_jesdn_BBEC1_year],eta_jes_split_BBEC1_year[jet2index2p5_jesdn_BBEC1_year],phi_jes_split_BBEC1_year[jet2index2p5_jesdn_BBEC1_year], mass_jes_split_BBEC1_year[jet2index2p5_jesdn_BBEC1_year]);
                    pTj2_2p5_jesdn_BBEC1_year=Jet2_2p5_jesdn_BBEC1_year.Pt();
                    yj2_2p5_jesdn_BBEC1_year=Jet2_2p5_jesdn_BBEC1_year.Rapidity();
                    pT4ljj_2p5_jesdn_BBEC1_year=(Higgs+Jet1_2p5_jesdn_BBEC1_year+Jet2_2p5_jesdn_BBEC1_year).Pt();
                    mass4ljj_2p5_jesdn_BBEC1_year=(Higgs+Jet1_2p5_jesdn_BBEC1_year+Jet2_2p5_jesdn_BBEC1_year).M();
                    mj1j2_2p5_jesdn_BBEC1_year=(Jet1_2p5_jesdn_BBEC1_year+Jet2_2p5_jesdn_BBEC1_year).M();     
                    dEtaj1j2_2p5_jesdn_BBEC1_year=TMath::Abs(Jet1_2p5_jesdn_BBEC1_year.Eta()-Jet2_2p5_jesdn_BBEC1_year.Eta());
	            dPhij1j2_2p5_jesdn_BBEC1_year=deltaPhi(Jet1_2p5_jesdn_BBEC1_year.Phi(),Jet2_2p5_jesdn_BBEC1_year.Phi());
                    dPhiHj1j2_2p5_jesdn_BBEC1_year=deltaPhi(Higgs.Phi(),(Jet1_2p5_jesdn_BBEC1_year+Jet2_2p5_jesdn_BBEC1_year).Phi());

                }    

                pt_jesdn_split_BBEC1_year.clear();

		eta_jes_split_BBEC1_year.clear();
		phi_jes_split_BBEC1_year.clear();
		mass_jes_split_BBEC1_year.clear();


////////////////////////////////////////////////////////////////////////////////
////    Starting "EC2" JES uncertainty source
////////////////////////////////////////////////////////////////////////////////

                TauC_Inc_0j_EnergyWgt_jesup_EC2 = TauC(pt_jesup_split_EC2, eta_jes_split_EC2,phi_jes_split_EC2,mass_jes_split_EC2, Higgs);
                TauB_Inc_0j_pTWgt_jesup_EC2 = TauB(pt_jesup_split_EC2, eta_jes_split_EC2,phi_jes_split_EC2,mass_jes_split_EC2, Higgs);


                for( unsigned int k = 0; k< pt_jesup_split_EC2.size(); k++) {
                    if (pt_jesup_split_EC2[k]<30.0 || abs(eta_jes_split_EC2[k])>4.7) continue;
		    TLorentzVector thisJet_jesup_EC2;
		    thisJet_jesup_EC2.SetPtEtaPhiM(pt_jesup_split_EC2[k],eta_jes_split_EC2[k],phi_jes_split_EC2[k],mass_jes_split_EC2[k]);

                        njets_pt30_eta4p7_jesup_EC2+=1;  

                        if (thisJet_jesup_EC2.Pt()>jet1pt_jesup_EC2) {
                            jet2pt_jesup_EC2=jet1pt_jesup_EC2; jet2index_jesup_EC2=jet1index_jesup_EC2;
                            jet1pt_jesup_EC2=thisJet_jesup_EC2.Pt(); jet1index_jesup_EC2=k;
                        } else if (thisJet_jesup_EC2.Pt()>jet2pt_jesup_EC2) {
                            jet2pt_jesup_EC2=thisJet_jesup_EC2.Pt(); jet2index_jesup_EC2=k;
                        }
                        if (abs(thisJet_jesup_EC2.Eta())<2.5) {
                            njets_pt30_eta2p5_jesup_EC2+=1;
                            if (thisJet_jesup_EC2.Pt()>jet1pt2p5_jesup_EC2) {
                                jet2pt2p5_jesup_EC2=jet1pt2p5_jesup_EC2; jet2index2p5_jesup_EC2=jet1index2p5_jesup_EC2;
                                jet1pt2p5_jesup_EC2=thisJet_jesup_EC2.Pt(); jet1index2p5_jesup_EC2=k;
                            } else if (thisJet_jesup_EC2.Pt()>jet2pt2p5_jesup_EC2) {
                                jet2pt2p5_jesup_EC2=thisJet_jesup_EC2.Pt(); jet2index2p5_jesup_EC2=k;
                            }
                        }
                }
                cout<<njets_pt30_eta4p7_jesup_EC2<<" jets (jesup_EC2)"<<endl;


		TLorentzVector Jet1_jesup_EC2, Jet1_2p5_jesup_EC2, Jet2_jesup_EC2, Jet2_2p5_jesup_EC2;
                if (njets_pt30_eta4p7_jesup_EC2 > 0) {
		    Jet1_jesup_EC2.SetPtEtaPhiM(pt_jesup_split_EC2[jet1index_jesup_EC2],eta_jes_split_EC2[jet1index_jesup_EC2],phi_jes_split_EC2[jet1index_jesup_EC2], mass_jes_split_EC2[jet1index_jesup_EC2]);

                    pt_leadingjet_pt30_eta4p7_jesup_EC2=Jet1_jesup_EC2.Pt(); 
                    pTj1_jesup_EC2=Jet1_jesup_EC2.Pt(); 
                    etaj1_jesup_EC2=Jet1_jesup_EC2.Eta();
                    yj1_jesup_EC2=Jet1_jesup_EC2.Rapidity();
		    pT4lj_jesup_EC2=(Higgs+Jet1_jesup_EC2).Pt();
		    mass4lj_jesup_EC2=(Higgs+Jet1_jesup_EC2).M();
		    dPhiHj1_jesup_EC2=deltaPhi(Higgs.Phi(),Jet1_jesup_EC2.Phi());
                    dyHj1_jesup_EC2=TMath::Abs(rapidity4l-yj1_jesup_EC2);
                }
                if (njets_pt30_eta4p7_jesup_EC2 > 1) {
                    Jet2_jesup_EC2.SetPtEtaPhiM(pt_jesup_split_EC2[jet2index_jesup_EC2],eta_jes_split_EC2[jet2index_jesup_EC2],phi_jes_split_EC2[jet2index_jesup_EC2], mass_jes_split_EC2[jet2index_jesup_EC2]);

                    pTj2_jesup_EC2=Jet2_jesup_EC2.Pt();
                    etaj2_jesup_EC2=Jet2_jesup_EC2.Eta();
                    yj2_jesup_EC2=Jet2_jesup_EC2.Rapidity();
                    pT4ljj_jesup_EC2=(Higgs+Jet1_jesup_EC2+Jet2_jesup_EC2).Pt();
                    mass4ljj_jesup_EC2=(Higgs+Jet1_jesup_EC2+Jet2_jesup_EC2).M();
		    mj1j2_jesup_EC2=(Jet1_jesup_EC2+Jet2_jesup_EC2).M();
                    dEtaj1j2_jesup_EC2=TMath::Abs(Jet1_jesup_EC2.Eta()-Jet2_jesup_EC2.Eta());
                    dPhij1j2_jesup_EC2=deltaPhi(Jet1_jesup_EC2.Phi(),Jet2_jesup_EC2.Phi());
                    dPhiHj1j2_jesup_EC2=deltaPhi(Higgs.Phi(),(Jet1_jesup_EC2+Jet2_jesup_EC2).Phi());
                }
                if (njets_pt30_eta2p5_jesup_EC2 > 0) {
                    Jet1_2p5_jesup_EC2.SetPtEtaPhiM(pt_jesup_split_EC2[jet1index2p5_jesup_EC2],eta_jes_split_EC2[jet1index2p5_jesup_EC2],phi_jes_split_EC2[jet1index2p5_jesup_EC2], mass_jes_split_EC2[jet1index2p5_jesup_EC2]);
                    pt_leadingjet_pt30_eta2p5_jesup_EC2=Jet1_2p5_jesup_EC2.Pt();
                    pTj1_2p5_jesup_EC2=Jet1_2p5_jesup_EC2.Pt();
                    yj1_2p5_jesup_EC2=Jet1_2p5_jesup_EC2.Rapidity();
                    pT4lj_2p5_jesup_EC2=(Higgs+Jet1_2p5_jesup_EC2).Pt();
                    mass4lj_2p5_jesup_EC2=(Higgs+Jet1_2p5_jesup_EC2).M();
		    dPhiHj1_2p5_jesup_EC2=deltaPhi(Higgs.Phi(),Jet1_2p5_jesup_EC2.Phi());
                    dyHj1_2p5_jesup_EC2=TMath::Abs(rapidity4l-yj1_2p5_jesup_EC2);
                }
                if (njets_pt30_eta2p5_jesup_EC2 > 1) {
                    Jet2_2p5_jesup_EC2.SetPtEtaPhiM(pt_jesup_split_EC2[jet2index2p5_jesup_EC2],eta_jes_split_EC2[jet2index2p5_jesup_EC2],phi_jes_split_EC2[jet2index2p5_jesup_EC2], mass_jes_split_EC2[jet2index2p5_jesup_EC2]);
                    pTj2_2p5_jesup_EC2=Jet2_2p5_jesup_EC2.Pt();
                    yj2_2p5_jesup_EC2=Jet2_2p5_jesup_EC2.Rapidity();
                    pT4ljj_2p5_jesup_EC2=(Higgs+Jet1_2p5_jesup_EC2+Jet2_2p5_jesup_EC2).Pt();
                    mass4ljj_2p5_jesup_EC2=(Higgs+Jet1_2p5_jesup_EC2+Jet2_2p5_jesup_EC2).M();
                    mj1j2_2p5_jesup_EC2=(Jet1_2p5_jesup_EC2+Jet2_2p5_jesup_EC2).M();                    
                    dEtaj1j2_2p5_jesup_EC2=TMath::Abs(Jet1_2p5_jesup_EC2.Eta()-Jet2_2p5_jesup_EC2.Eta());
                    dPhij1j2_2p5_jesup_EC2=deltaPhi(Jet1_2p5_jesup_EC2.Phi(),Jet2_2p5_jesup_EC2.Phi());
                    dPhiHj1j2_2p5_jesup_EC2=deltaPhi(Higgs.Phi(),(Jet1_2p5_jesup_EC2+Jet2_2p5_jesup_EC2).Phi());

                }

		pt_jesup_split_EC2.clear();  

/////////////////
// EC2 dn start
                TauC_Inc_0j_EnergyWgt_jesdn_EC2 = TauC(pt_jesdn_split_EC2, eta_jes_split_EC2,phi_jes_split_EC2,mass_jes_split_EC2, Higgs);
                TauB_Inc_0j_pTWgt_jesdn_EC2 = TauB(pt_jesdn_split_EC2, eta_jes_split_EC2,phi_jes_split_EC2,mass_jes_split_EC2, Higgs);


                for( unsigned int k = 0; k< pt_jesdn_split_EC2.size(); k++) {
                    if (pt_jesdn_split_EC2[k]<30.0 || abs(eta_jes_split_EC2[k])>4.7) continue;
                    TLorentzVector thisJet_jesdn_EC2;
		    thisJet_jesdn_EC2.SetPtEtaPhiM(pt_jesdn_split_EC2[k],eta_jes_split_EC2[k],phi_jes_split_EC2[k],mass_jes_split_EC2[k]);

                    bool isclean_H4l=true;
                        njets_pt30_eta4p7_jesdn_EC2+=1;
			if (thisJet_jesdn_EC2.Pt()>jet1pt_jesdn_EC2) {
                            jet2pt_jesdn_EC2=jet1pt_jesdn_EC2; jet2index_jesdn_EC2=jet1index_jesdn_EC2;
                            jet1pt_jesdn_EC2=thisJet_jesdn_EC2.Pt(); jet1index_jesdn_EC2=k;
                        } else if (thisJet_jesdn_EC2.Pt()>jet2pt_jesdn_EC2) {
                            jet2pt_jesdn_EC2=thisJet_jesdn_EC2.Pt(); jet2index_jesdn_EC2=k;
                        }
                        if (abs(thisJet_jesdn_EC2.Eta())<2.5) {
                            njets_pt30_eta2p5_jesdn_EC2+=1;
                            if (thisJet_jesdn_EC2.Pt()>jet1pt2p5_jesdn_EC2) {
                                jet2pt2p5_jesdn_EC2=jet1pt2p5_jesdn_EC2; jet2index2p5_jesdn_EC2=jet1index2p5_jesdn_EC2;
                                jet1pt2p5_jesdn_EC2=thisJet_jesdn_EC2.Pt(); jet1index2p5_jesdn_EC2=k;
                            } else if (thisJet_jesdn_EC2.Pt()>jet2pt2p5_jesdn_EC2) {
                                jet2pt2p5_jesdn_EC2=thisJet_jesdn_EC2.Pt(); jet2index2p5_jesdn_EC2=k;
                            }
                        }
                }

// Filling EC2 dn variables

                TLorentzVector Jet1_jesdn_EC2, Jet1_2p5_jesdn_EC2, Jet2_jesdn_EC2, Jet2_2p5_jesdn_EC2;

                if (njets_pt30_eta4p7_jesdn_EC2 > 0) { 
                    Jet1_jesdn_EC2.SetPtEtaPhiM(pt_jesdn_split_EC2[jet1index_jesdn_EC2],eta_jes_split_EC2[jet1index_jesdn_EC2],phi_jes_split_EC2[jet1index_jesdn_EC2], mass_jes_split_EC2[jet1index_jesdn_EC2]);
                    pt_leadingjet_pt30_eta4p7_jesdn_EC2=Jet1_jesdn_EC2.Pt(); 
                    pTj1_jesdn_EC2=Jet1_jesdn_EC2.Pt(); 
                    etaj1_jesdn_EC2=Jet1_jesdn_EC2.Eta();
                    yj1_jesdn_EC2=Jet1_jesdn_EC2.Rapidity();
                    pT4lj_jesdn_EC2=(Higgs+Jet1_jesdn_EC2).Pt();
                    mass4lj_jesdn_EC2=(Higgs+Jet1_jesdn_EC2).M();
                    dPhiHj1_jesdn_EC2=deltaPhi(Higgs.Phi(),Jet1_jesdn_EC2.Phi());
                    dyHj1_jesdn_EC2=TMath::Abs(rapidity4l-yj1_jesdn_EC2);
                }    
                if (njets_pt30_eta4p7_jesdn_EC2 > 1) { 
                    Jet2_jesdn_EC2.SetPtEtaPhiM(pt_jesdn_split_EC2[jet2index_jesdn_EC2],eta_jes_split_EC2[jet2index_jesdn_EC2],phi_jes_split_EC2[jet2index_jesdn_EC2], mass_jes_split_EC2[jet2index_jesdn_EC2]);

                    pTj2_jesdn_EC2=Jet2_jesdn_EC2.Pt();
                    etaj2_jesdn_EC2=Jet2_jesdn_EC2.Eta();
                    yj2_jesdn_EC2=Jet2_jesdn_EC2.Rapidity();
                    pT4ljj_jesdn_EC2=(Higgs+Jet1_jesdn_EC2+Jet2_jesdn_EC2).Pt();
                    mass4ljj_jesdn_EC2=(Higgs+Jet1_jesdn_EC2+Jet2_jesdn_EC2).M();
                    dEtaj1j2_jesdn_EC2=TMath::Abs(Jet1_jesdn_EC2.Eta()-Jet2_jesdn_EC2.Eta());
                    dPhij1j2_jesdn_EC2=deltaPhi(Jet1_jesdn_EC2.Phi(),Jet2_jesdn_EC2.Phi());
                    dPhiHj1j2_jesdn_EC2=deltaPhi(Higgs.Phi(),(Jet1_jesdn_EC2+Jet2_jesdn_EC2).Phi());
                }    
                if (njets_pt30_eta2p5_jesdn_EC2 > 0) { 
                    Jet1_2p5_jesdn_EC2.SetPtEtaPhiM(pt_jesdn_split_EC2[jet1index2p5_jesdn_EC2],eta_jes_split_EC2[jet1index2p5_jesdn_EC2],phi_jes_split_EC2[jet1index2p5_jesdn_EC2], mass_jes_split_EC2[jet1index2p5_jesdn_EC2]);
                    pt_leadingjet_pt30_eta2p5_jesdn_EC2=Jet1_2p5_jesdn_EC2.Pt();
                    pTj1_2p5_jesdn_EC2=Jet1_2p5_jesdn_EC2.Pt();
                    yj1_2p5_jesdn_EC2=Jet1_2p5_jesdn_EC2.Rapidity();
                    pT4lj_2p5_jesdn_EC2=(Higgs+Jet1_2p5_jesdn_EC2).Pt();
                    mass4lj_2p5_jesdn_EC2=(Higgs+Jet1_2p5_jesdn_EC2).M();
		    dPhiHj1_2p5_jesdn_EC2=deltaPhi(Higgs.Phi(),Jet1_2p5_jesdn_EC2.Phi());
                    dyHj1_2p5_jesdn_EC2=TMath::Abs(rapidity4l-yj1_2p5_jesdn_EC2);
                }    
                if (njets_pt30_eta2p5_jesdn_EC2 > 1) { 
                    Jet2_2p5_jesdn_EC2.SetPtEtaPhiM(pt_jesdn_split_EC2[jet2index2p5_jesdn_EC2],eta_jes_split_EC2[jet2index2p5_jesdn_EC2],phi_jes_split_EC2[jet2index2p5_jesdn_EC2], mass_jes_split_EC2[jet2index2p5_jesdn_EC2]);
                    pTj2_2p5_jesdn_EC2=Jet2_2p5_jesdn_EC2.Pt();
                    yj2_2p5_jesdn_EC2=Jet2_2p5_jesdn_EC2.Rapidity();
                    pT4ljj_2p5_jesdn_EC2=(Higgs+Jet1_2p5_jesdn_EC2+Jet2_2p5_jesdn_EC2).Pt();
                    mass4ljj_2p5_jesdn_EC2=(Higgs+Jet1_2p5_jesdn_EC2+Jet2_2p5_jesdn_EC2).M();
                    mj1j2_2p5_jesdn_EC2=(Jet1_2p5_jesdn_EC2+Jet2_2p5_jesdn_EC2).M();     
                    dEtaj1j2_2p5_jesdn_EC2=TMath::Abs(Jet1_2p5_jesdn_EC2.Eta()-Jet2_2p5_jesdn_EC2.Eta());
	            dPhij1j2_2p5_jesdn_EC2=deltaPhi(Jet1_2p5_jesdn_EC2.Phi(),Jet2_2p5_jesdn_EC2.Phi());
                    dPhiHj1j2_2p5_jesdn_EC2=deltaPhi(Higgs.Phi(),(Jet1_2p5_jesdn_EC2+Jet2_2p5_jesdn_EC2).Phi());

                }    

                pt_jesdn_split_EC2.clear();

		eta_jes_split_EC2.clear();
		phi_jes_split_EC2.clear();
		mass_jes_split_EC2.clear();


////////////////////////////////////////////////////////////////////////////////
////    Starting "EC2_year" JES uncertainty source
////////////////////////////////////////////////////////////////////////////////

                TauC_Inc_0j_EnergyWgt_jesup_EC2_year = TauC(pt_jesup_split_EC2_year, eta_jes_split_EC2_year,phi_jes_split_EC2_year,mass_jes_split_EC2_year, Higgs);
                TauB_Inc_0j_pTWgt_jesup_EC2_year = TauB(pt_jesup_split_EC2_year, eta_jes_split_EC2_year,phi_jes_split_EC2_year,mass_jes_split_EC2_year, Higgs);


                for( unsigned int k = 0; k< pt_jesup_split_EC2_year.size(); k++) {
                    if (pt_jesup_split_EC2_year[k]<30.0 || abs(eta_jes_split_EC2_year[k])>4.7) continue;
		    TLorentzVector thisJet_jesup_EC2_year;
		    thisJet_jesup_EC2_year.SetPtEtaPhiM(pt_jesup_split_EC2_year[k],eta_jes_split_EC2_year[k],phi_jes_split_EC2_year[k],mass_jes_split_EC2_year[k]);

                        njets_pt30_eta4p7_jesup_EC2_year+=1;  

                        if (thisJet_jesup_EC2_year.Pt()>jet1pt_jesup_EC2_year) {
                            jet2pt_jesup_EC2_year=jet1pt_jesup_EC2_year; jet2index_jesup_EC2_year=jet1index_jesup_EC2_year;
                            jet1pt_jesup_EC2_year=thisJet_jesup_EC2_year.Pt(); jet1index_jesup_EC2_year=k;
                        } else if (thisJet_jesup_EC2_year.Pt()>jet2pt_jesup_EC2_year) {
                            jet2pt_jesup_EC2_year=thisJet_jesup_EC2_year.Pt(); jet2index_jesup_EC2_year=k;
                        }
                        if (abs(thisJet_jesup_EC2_year.Eta())<2.5) {
                            njets_pt30_eta2p5_jesup_EC2_year+=1;
                            if (thisJet_jesup_EC2_year.Pt()>jet1pt2p5_jesup_EC2_year) {
                                jet2pt2p5_jesup_EC2_year=jet1pt2p5_jesup_EC2_year; jet2index2p5_jesup_EC2_year=jet1index2p5_jesup_EC2_year;
                                jet1pt2p5_jesup_EC2_year=thisJet_jesup_EC2_year.Pt(); jet1index2p5_jesup_EC2_year=k;
                            } else if (thisJet_jesup_EC2_year.Pt()>jet2pt2p5_jesup_EC2_year) {
                                jet2pt2p5_jesup_EC2_year=thisJet_jesup_EC2_year.Pt(); jet2index2p5_jesup_EC2_year=k;
                            }
                        }
                }
                cout<<njets_pt30_eta4p7_jesup_EC2_year<<" jets (jesup_EC2_year)"<<endl;


		TLorentzVector Jet1_jesup_EC2_year, Jet1_2p5_jesup_EC2_year, Jet2_jesup_EC2_year, Jet2_2p5_jesup_EC2_year;
                if (njets_pt30_eta4p7_jesup_EC2_year > 0) {
		    Jet1_jesup_EC2_year.SetPtEtaPhiM(pt_jesup_split_EC2_year[jet1index_jesup_EC2_year],eta_jes_split_EC2_year[jet1index_jesup_EC2_year],phi_jes_split_EC2_year[jet1index_jesup_EC2_year], mass_jes_split_EC2_year[jet1index_jesup_EC2_year]);

                    pt_leadingjet_pt30_eta4p7_jesup_EC2_year=Jet1_jesup_EC2_year.Pt(); 
                    pTj1_jesup_EC2_year=Jet1_jesup_EC2_year.Pt(); 
                    etaj1_jesup_EC2_year=Jet1_jesup_EC2_year.Eta();
                    yj1_jesup_EC2_year=Jet1_jesup_EC2_year.Rapidity();
		    pT4lj_jesup_EC2_year=(Higgs+Jet1_jesup_EC2_year).Pt();
		    mass4lj_jesup_EC2_year=(Higgs+Jet1_jesup_EC2_year).M();
		    dPhiHj1_jesup_EC2_year=deltaPhi(Higgs.Phi(),Jet1_jesup_EC2_year.Phi());
                    dyHj1_jesup_EC2_year=TMath::Abs(rapidity4l-yj1_jesup_EC2_year);
                }
                if (njets_pt30_eta4p7_jesup_EC2_year > 1) {
                    Jet2_jesup_EC2_year.SetPtEtaPhiM(pt_jesup_split_EC2_year[jet2index_jesup_EC2_year],eta_jes_split_EC2_year[jet2index_jesup_EC2_year],phi_jes_split_EC2_year[jet2index_jesup_EC2_year], mass_jes_split_EC2_year[jet2index_jesup_EC2_year]);

                    pTj2_jesup_EC2_year=Jet2_jesup_EC2_year.Pt();
                    etaj2_jesup_EC2_year=Jet2_jesup_EC2_year.Eta();
                    yj2_jesup_EC2_year=Jet2_jesup_EC2_year.Rapidity();
                    pT4ljj_jesup_EC2_year=(Higgs+Jet1_jesup_EC2_year+Jet2_jesup_EC2_year).Pt();
                    mass4ljj_jesup_EC2_year=(Higgs+Jet1_jesup_EC2_year+Jet2_jesup_EC2_year).M();
		    mj1j2_jesup_EC2_year=(Jet1_jesup_EC2_year+Jet2_jesup_EC2_year).M();
                    dEtaj1j2_jesup_EC2_year=TMath::Abs(Jet1_jesup_EC2_year.Eta()-Jet2_jesup_EC2_year.Eta());
                    dPhij1j2_jesup_EC2_year=deltaPhi(Jet1_jesup_EC2_year.Phi(),Jet2_jesup_EC2_year.Phi());
                    dPhiHj1j2_jesup_EC2_year=deltaPhi(Higgs.Phi(),(Jet1_jesup_EC2_year+Jet2_jesup_EC2_year).Phi());
                }
                if (njets_pt30_eta2p5_jesup_EC2_year > 0) {
                    Jet1_2p5_jesup_EC2_year.SetPtEtaPhiM(pt_jesup_split_EC2_year[jet1index2p5_jesup_EC2_year],eta_jes_split_EC2_year[jet1index2p5_jesup_EC2_year],phi_jes_split_EC2_year[jet1index2p5_jesup_EC2_year], mass_jes_split_EC2_year[jet1index2p5_jesup_EC2_year]);
                    pt_leadingjet_pt30_eta2p5_jesup_EC2_year=Jet1_2p5_jesup_EC2_year.Pt();
                    pTj1_2p5_jesup_EC2_year=Jet1_2p5_jesup_EC2_year.Pt();
                    yj1_2p5_jesup_EC2_year=Jet1_2p5_jesup_EC2_year.Rapidity();
                    pT4lj_2p5_jesup_EC2_year=(Higgs+Jet1_2p5_jesup_EC2_year).Pt();
                    mass4lj_2p5_jesup_EC2_year=(Higgs+Jet1_2p5_jesup_EC2_year).M();
		    dPhiHj1_2p5_jesup_EC2_year=deltaPhi(Higgs.Phi(),Jet1_2p5_jesup_EC2_year.Phi());
                    dyHj1_2p5_jesup_EC2_year=TMath::Abs(rapidity4l-yj1_2p5_jesup_EC2_year);
                }
                if (njets_pt30_eta2p5_jesup_EC2_year > 1) {
                    Jet2_2p5_jesup_EC2_year.SetPtEtaPhiM(pt_jesup_split_EC2_year[jet2index2p5_jesup_EC2_year],eta_jes_split_EC2_year[jet2index2p5_jesup_EC2_year],phi_jes_split_EC2_year[jet2index2p5_jesup_EC2_year], mass_jes_split_EC2_year[jet2index2p5_jesup_EC2_year]);
                    pTj2_2p5_jesup_EC2_year=Jet2_2p5_jesup_EC2_year.Pt();
                    yj2_2p5_jesup_EC2_year=Jet2_2p5_jesup_EC2_year.Rapidity();
                    pT4ljj_2p5_jesup_EC2_year=(Higgs+Jet1_2p5_jesup_EC2_year+Jet2_2p5_jesup_EC2_year).Pt();
                    mass4ljj_2p5_jesup_EC2_year=(Higgs+Jet1_2p5_jesup_EC2_year+Jet2_2p5_jesup_EC2_year).M();
                    mj1j2_2p5_jesup_EC2_year=(Jet1_2p5_jesup_EC2_year+Jet2_2p5_jesup_EC2_year).M();                    
                    dEtaj1j2_2p5_jesup_EC2_year=TMath::Abs(Jet1_2p5_jesup_EC2_year.Eta()-Jet2_2p5_jesup_EC2_year.Eta());
                    dPhij1j2_2p5_jesup_EC2_year=deltaPhi(Jet1_2p5_jesup_EC2_year.Phi(),Jet2_2p5_jesup_EC2_year.Phi());
                    dPhiHj1j2_2p5_jesup_EC2_year=deltaPhi(Higgs.Phi(),(Jet1_2p5_jesup_EC2_year+Jet2_2p5_jesup_EC2_year).Phi());

                }

		pt_jesup_split_EC2_year.clear();  

/////////////////
// EC2_year dn start
                TauC_Inc_0j_EnergyWgt_jesdn_EC2_year = TauC(pt_jesdn_split_EC2_year, eta_jes_split_EC2_year,phi_jes_split_EC2_year,mass_jes_split_EC2_year, Higgs);
                TauB_Inc_0j_pTWgt_jesdn_EC2_year = TauB(pt_jesdn_split_EC2_year, eta_jes_split_EC2_year,phi_jes_split_EC2_year,mass_jes_split_EC2_year, Higgs);


                for( unsigned int k = 0; k< pt_jesdn_split_EC2_year.size(); k++) {
                    if (pt_jesdn_split_EC2_year[k]<30.0 || abs(eta_jes_split_EC2_year[k])>4.7) continue;
                    TLorentzVector thisJet_jesdn_EC2_year;
		    thisJet_jesdn_EC2_year.SetPtEtaPhiM(pt_jesdn_split_EC2_year[k],eta_jes_split_EC2_year[k],phi_jes_split_EC2_year[k],mass_jes_split_EC2_year[k]);

                    bool isclean_H4l=true;
                        njets_pt30_eta4p7_jesdn_EC2_year+=1;
			if (thisJet_jesdn_EC2_year.Pt()>jet1pt_jesdn_EC2_year) {
                            jet2pt_jesdn_EC2_year=jet1pt_jesdn_EC2_year; jet2index_jesdn_EC2_year=jet1index_jesdn_EC2_year;
                            jet1pt_jesdn_EC2_year=thisJet_jesdn_EC2_year.Pt(); jet1index_jesdn_EC2_year=k;
                        } else if (thisJet_jesdn_EC2_year.Pt()>jet2pt_jesdn_EC2_year) {
                            jet2pt_jesdn_EC2_year=thisJet_jesdn_EC2_year.Pt(); jet2index_jesdn_EC2_year=k;
                        }
                        if (abs(thisJet_jesdn_EC2_year.Eta())<2.5) {
                            njets_pt30_eta2p5_jesdn_EC2_year+=1;
                            if (thisJet_jesdn_EC2_year.Pt()>jet1pt2p5_jesdn_EC2_year) {
                                jet2pt2p5_jesdn_EC2_year=jet1pt2p5_jesdn_EC2_year; jet2index2p5_jesdn_EC2_year=jet1index2p5_jesdn_EC2_year;
                                jet1pt2p5_jesdn_EC2_year=thisJet_jesdn_EC2_year.Pt(); jet1index2p5_jesdn_EC2_year=k;
                            } else if (thisJet_jesdn_EC2_year.Pt()>jet2pt2p5_jesdn_EC2_year) {
                                jet2pt2p5_jesdn_EC2_year=thisJet_jesdn_EC2_year.Pt(); jet2index2p5_jesdn_EC2_year=k;
                            }
                        }
                }

// Filling EC2_year dn variables

                TLorentzVector Jet1_jesdn_EC2_year, Jet1_2p5_jesdn_EC2_year, Jet2_jesdn_EC2_year, Jet2_2p5_jesdn_EC2_year;

                if (njets_pt30_eta4p7_jesdn_EC2_year > 0) { 
                    Jet1_jesdn_EC2_year.SetPtEtaPhiM(pt_jesdn_split_EC2_year[jet1index_jesdn_EC2_year],eta_jes_split_EC2_year[jet1index_jesdn_EC2_year],phi_jes_split_EC2_year[jet1index_jesdn_EC2_year], mass_jes_split_EC2_year[jet1index_jesdn_EC2_year]);
                    pt_leadingjet_pt30_eta4p7_jesdn_EC2_year=Jet1_jesdn_EC2_year.Pt(); 
                    pTj1_jesdn_EC2_year=Jet1_jesdn_EC2_year.Pt(); 
                    etaj1_jesdn_EC2_year=Jet1_jesdn_EC2_year.Eta();
                    yj1_jesdn_EC2_year=Jet1_jesdn_EC2_year.Rapidity();
                    pT4lj_jesdn_EC2_year=(Higgs+Jet1_jesdn_EC2_year).Pt();
                    mass4lj_jesdn_EC2_year=(Higgs+Jet1_jesdn_EC2_year).M();
                    dPhiHj1_jesdn_EC2_year=deltaPhi(Higgs.Phi(),Jet1_jesdn_EC2_year.Phi());
                    dyHj1_jesdn_EC2_year=TMath::Abs(rapidity4l-yj1_jesdn_EC2_year);
                }    
                if (njets_pt30_eta4p7_jesdn_EC2_year > 1) { 
                    Jet2_jesdn_EC2_year.SetPtEtaPhiM(pt_jesdn_split_EC2_year[jet2index_jesdn_EC2_year],eta_jes_split_EC2_year[jet2index_jesdn_EC2_year],phi_jes_split_EC2_year[jet2index_jesdn_EC2_year], mass_jes_split_EC2_year[jet2index_jesdn_EC2_year]);

                    pTj2_jesdn_EC2_year=Jet2_jesdn_EC2_year.Pt();
                    etaj2_jesdn_EC2_year=Jet2_jesdn_EC2_year.Eta();
                    yj2_jesdn_EC2_year=Jet2_jesdn_EC2_year.Rapidity();
                    pT4ljj_jesdn_EC2_year=(Higgs+Jet1_jesdn_EC2_year+Jet2_jesdn_EC2_year).Pt();
                    mass4ljj_jesdn_EC2_year=(Higgs+Jet1_jesdn_EC2_year+Jet2_jesdn_EC2_year).M();
                    dEtaj1j2_jesdn_EC2_year=TMath::Abs(Jet1_jesdn_EC2_year.Eta()-Jet2_jesdn_EC2_year.Eta());
                    dPhij1j2_jesdn_EC2_year=deltaPhi(Jet1_jesdn_EC2_year.Phi(),Jet2_jesdn_EC2_year.Phi());
                    dPhiHj1j2_jesdn_EC2_year=deltaPhi(Higgs.Phi(),(Jet1_jesdn_EC2_year+Jet2_jesdn_EC2_year).Phi());
                }    
                if (njets_pt30_eta2p5_jesdn_EC2_year > 0) { 
                    Jet1_2p5_jesdn_EC2_year.SetPtEtaPhiM(pt_jesdn_split_EC2_year[jet1index2p5_jesdn_EC2_year],eta_jes_split_EC2_year[jet1index2p5_jesdn_EC2_year],phi_jes_split_EC2_year[jet1index2p5_jesdn_EC2_year], mass_jes_split_EC2_year[jet1index2p5_jesdn_EC2_year]);
                    pt_leadingjet_pt30_eta2p5_jesdn_EC2_year=Jet1_2p5_jesdn_EC2_year.Pt();
                    pTj1_2p5_jesdn_EC2_year=Jet1_2p5_jesdn_EC2_year.Pt();
                    yj1_2p5_jesdn_EC2_year=Jet1_2p5_jesdn_EC2_year.Rapidity();
                    pT4lj_2p5_jesdn_EC2_year=(Higgs+Jet1_2p5_jesdn_EC2_year).Pt();
                    mass4lj_2p5_jesdn_EC2_year=(Higgs+Jet1_2p5_jesdn_EC2_year).M();
		    dPhiHj1_2p5_jesdn_EC2_year=deltaPhi(Higgs.Phi(),Jet1_2p5_jesdn_EC2_year.Phi());
                    dyHj1_2p5_jesdn_EC2_year=TMath::Abs(rapidity4l-yj1_2p5_jesdn_EC2_year);
                }    
                if (njets_pt30_eta2p5_jesdn_EC2_year > 1) { 
                    Jet2_2p5_jesdn_EC2_year.SetPtEtaPhiM(pt_jesdn_split_EC2_year[jet2index2p5_jesdn_EC2_year],eta_jes_split_EC2_year[jet2index2p5_jesdn_EC2_year],phi_jes_split_EC2_year[jet2index2p5_jesdn_EC2_year], mass_jes_split_EC2_year[jet2index2p5_jesdn_EC2_year]);
                    pTj2_2p5_jesdn_EC2_year=Jet2_2p5_jesdn_EC2_year.Pt();
                    yj2_2p5_jesdn_EC2_year=Jet2_2p5_jesdn_EC2_year.Rapidity();
                    pT4ljj_2p5_jesdn_EC2_year=(Higgs+Jet1_2p5_jesdn_EC2_year+Jet2_2p5_jesdn_EC2_year).Pt();
                    mass4ljj_2p5_jesdn_EC2_year=(Higgs+Jet1_2p5_jesdn_EC2_year+Jet2_2p5_jesdn_EC2_year).M();
                    mj1j2_2p5_jesdn_EC2_year=(Jet1_2p5_jesdn_EC2_year+Jet2_2p5_jesdn_EC2_year).M();     
                    dEtaj1j2_2p5_jesdn_EC2_year=TMath::Abs(Jet1_2p5_jesdn_EC2_year.Eta()-Jet2_2p5_jesdn_EC2_year.Eta());
	            dPhij1j2_2p5_jesdn_EC2_year=deltaPhi(Jet1_2p5_jesdn_EC2_year.Phi(),Jet2_2p5_jesdn_EC2_year.Phi());
                    dPhiHj1j2_2p5_jesdn_EC2_year=deltaPhi(Higgs.Phi(),(Jet1_2p5_jesdn_EC2_year+Jet2_2p5_jesdn_EC2_year).Phi());

                }    

                pt_jesdn_split_EC2_year.clear();

		eta_jes_split_EC2_year.clear();
		phi_jes_split_EC2_year.clear();
		mass_jes_split_EC2_year.clear();



////////////////////////////////////////////////////////////////////////////////
////    Starting "FlavQCD" JES uncertainty source
////////////////////////////////////////////////////////////////////////////////

                TauC_Inc_0j_EnergyWgt_jesup_FlavQCD = TauC(pt_jesup_split_FlavQCD, eta_jes_split_FlavQCD,phi_jes_split_FlavQCD,mass_jes_split_FlavQCD, Higgs);
                TauB_Inc_0j_pTWgt_jesup_FlavQCD = TauB(pt_jesup_split_FlavQCD, eta_jes_split_FlavQCD,phi_jes_split_FlavQCD,mass_jes_split_FlavQCD, Higgs);


                for( unsigned int k = 0; k< pt_jesup_split_FlavQCD.size(); k++) {
                    if (pt_jesup_split_FlavQCD[k]<30.0 || abs(eta_jes_split_FlavQCD[k])>4.7) continue;
		    TLorentzVector thisJet_jesup_FlavQCD;
		    thisJet_jesup_FlavQCD.SetPtEtaPhiM(pt_jesup_split_FlavQCD[k],eta_jes_split_FlavQCD[k],phi_jes_split_FlavQCD[k],mass_jes_split_FlavQCD[k]);

                        njets_pt30_eta4p7_jesup_FlavQCD+=1;  

                        if (thisJet_jesup_FlavQCD.Pt()>jet1pt_jesup_FlavQCD) {
                            jet2pt_jesup_FlavQCD=jet1pt_jesup_FlavQCD; jet2index_jesup_FlavQCD=jet1index_jesup_FlavQCD;
                            jet1pt_jesup_FlavQCD=thisJet_jesup_FlavQCD.Pt(); jet1index_jesup_FlavQCD=k;
                        } else if (thisJet_jesup_FlavQCD.Pt()>jet2pt_jesup_FlavQCD) {
                            jet2pt_jesup_FlavQCD=thisJet_jesup_FlavQCD.Pt(); jet2index_jesup_FlavQCD=k;
                        }
                        if (abs(thisJet_jesup_FlavQCD.Eta())<2.5) {
                            njets_pt30_eta2p5_jesup_FlavQCD+=1;
                            if (thisJet_jesup_FlavQCD.Pt()>jet1pt2p5_jesup_FlavQCD) {
                                jet2pt2p5_jesup_FlavQCD=jet1pt2p5_jesup_FlavQCD; jet2index2p5_jesup_FlavQCD=jet1index2p5_jesup_FlavQCD;
                                jet1pt2p5_jesup_FlavQCD=thisJet_jesup_FlavQCD.Pt(); jet1index2p5_jesup_FlavQCD=k;
                            } else if (thisJet_jesup_FlavQCD.Pt()>jet2pt2p5_jesup_FlavQCD) {
                                jet2pt2p5_jesup_FlavQCD=thisJet_jesup_FlavQCD.Pt(); jet2index2p5_jesup_FlavQCD=k;
                            }
                        }
                }
                cout<<njets_pt30_eta4p7_jesup_FlavQCD<<" jets (jesup_FlavQCD)"<<endl;


		TLorentzVector Jet1_jesup_FlavQCD, Jet1_2p5_jesup_FlavQCD, Jet2_jesup_FlavQCD, Jet2_2p5_jesup_FlavQCD;
                if (njets_pt30_eta4p7_jesup_FlavQCD > 0) {
		    Jet1_jesup_FlavQCD.SetPtEtaPhiM(pt_jesup_split_FlavQCD[jet1index_jesup_FlavQCD],eta_jes_split_FlavQCD[jet1index_jesup_FlavQCD],phi_jes_split_FlavQCD[jet1index_jesup_FlavQCD], mass_jes_split_FlavQCD[jet1index_jesup_FlavQCD]);

                    pt_leadingjet_pt30_eta4p7_jesup_FlavQCD=Jet1_jesup_FlavQCD.Pt(); 
                    pTj1_jesup_FlavQCD=Jet1_jesup_FlavQCD.Pt(); 
                    etaj1_jesup_FlavQCD=Jet1_jesup_FlavQCD.Eta();
                    yj1_jesup_FlavQCD=Jet1_jesup_FlavQCD.Rapidity();
		    pT4lj_jesup_FlavQCD=(Higgs+Jet1_jesup_FlavQCD).Pt();
		    mass4lj_jesup_FlavQCD=(Higgs+Jet1_jesup_FlavQCD).M();
		    dPhiHj1_jesup_FlavQCD=deltaPhi(Higgs.Phi(),Jet1_jesup_FlavQCD.Phi());
                    dyHj1_jesup_FlavQCD=TMath::Abs(rapidity4l-yj1_jesup_FlavQCD);
                }
                if (njets_pt30_eta4p7_jesup_FlavQCD > 1) {
                    Jet2_jesup_FlavQCD.SetPtEtaPhiM(pt_jesup_split_FlavQCD[jet2index_jesup_FlavQCD],eta_jes_split_FlavQCD[jet2index_jesup_FlavQCD],phi_jes_split_FlavQCD[jet2index_jesup_FlavQCD], mass_jes_split_FlavQCD[jet2index_jesup_FlavQCD]);

                    pTj2_jesup_FlavQCD=Jet2_jesup_FlavQCD.Pt();
                    etaj2_jesup_FlavQCD=Jet2_jesup_FlavQCD.Eta();
                    yj2_jesup_FlavQCD=Jet2_jesup_FlavQCD.Rapidity();
                    pT4ljj_jesup_FlavQCD=(Higgs+Jet1_jesup_FlavQCD+Jet2_jesup_FlavQCD).Pt();
                    mass4ljj_jesup_FlavQCD=(Higgs+Jet1_jesup_FlavQCD+Jet2_jesup_FlavQCD).M();
		    mj1j2_jesup_FlavQCD=(Jet1_jesup_FlavQCD+Jet2_jesup_FlavQCD).M();
                    dEtaj1j2_jesup_FlavQCD=TMath::Abs(Jet1_jesup_FlavQCD.Eta()-Jet2_jesup_FlavQCD.Eta());
                    dPhij1j2_jesup_FlavQCD=deltaPhi(Jet1_jesup_FlavQCD.Phi(),Jet2_jesup_FlavQCD.Phi());
                    dPhiHj1j2_jesup_FlavQCD=deltaPhi(Higgs.Phi(),(Jet1_jesup_FlavQCD+Jet2_jesup_FlavQCD).Phi());
                }
                if (njets_pt30_eta2p5_jesup_FlavQCD > 0) {
                    Jet1_2p5_jesup_FlavQCD.SetPtEtaPhiM(pt_jesup_split_FlavQCD[jet1index2p5_jesup_FlavQCD],eta_jes_split_FlavQCD[jet1index2p5_jesup_FlavQCD],phi_jes_split_FlavQCD[jet1index2p5_jesup_FlavQCD], mass_jes_split_FlavQCD[jet1index2p5_jesup_FlavQCD]);
                    pt_leadingjet_pt30_eta2p5_jesup_FlavQCD=Jet1_2p5_jesup_FlavQCD.Pt();
                    pTj1_2p5_jesup_FlavQCD=Jet1_2p5_jesup_FlavQCD.Pt();
                    yj1_2p5_jesup_FlavQCD=Jet1_2p5_jesup_FlavQCD.Rapidity();
                    pT4lj_2p5_jesup_FlavQCD=(Higgs+Jet1_2p5_jesup_FlavQCD).Pt();
                    mass4lj_2p5_jesup_FlavQCD=(Higgs+Jet1_2p5_jesup_FlavQCD).M();
		    dPhiHj1_2p5_jesup_FlavQCD=deltaPhi(Higgs.Phi(),Jet1_2p5_jesup_FlavQCD.Phi());
                    dyHj1_2p5_jesup_FlavQCD=TMath::Abs(rapidity4l-yj1_2p5_jesup_FlavQCD);
                }
                if (njets_pt30_eta2p5_jesup_FlavQCD > 1) {
                    Jet2_2p5_jesup_FlavQCD.SetPtEtaPhiM(pt_jesup_split_FlavQCD[jet2index2p5_jesup_FlavQCD],eta_jes_split_FlavQCD[jet2index2p5_jesup_FlavQCD],phi_jes_split_FlavQCD[jet2index2p5_jesup_FlavQCD], mass_jes_split_FlavQCD[jet2index2p5_jesup_FlavQCD]);
                    pTj2_2p5_jesup_FlavQCD=Jet2_2p5_jesup_FlavQCD.Pt();
                    yj2_2p5_jesup_FlavQCD=Jet2_2p5_jesup_FlavQCD.Rapidity();
                    pT4ljj_2p5_jesup_FlavQCD=(Higgs+Jet1_2p5_jesup_FlavQCD+Jet2_2p5_jesup_FlavQCD).Pt();
                    mass4ljj_2p5_jesup_FlavQCD=(Higgs+Jet1_2p5_jesup_FlavQCD+Jet2_2p5_jesup_FlavQCD).M();
                    mj1j2_2p5_jesup_FlavQCD=(Jet1_2p5_jesup_FlavQCD+Jet2_2p5_jesup_FlavQCD).M();                    
                    dEtaj1j2_2p5_jesup_FlavQCD=TMath::Abs(Jet1_2p5_jesup_FlavQCD.Eta()-Jet2_2p5_jesup_FlavQCD.Eta());
                    dPhij1j2_2p5_jesup_FlavQCD=deltaPhi(Jet1_2p5_jesup_FlavQCD.Phi(),Jet2_2p5_jesup_FlavQCD.Phi());
                    dPhiHj1j2_2p5_jesup_FlavQCD=deltaPhi(Higgs.Phi(),(Jet1_2p5_jesup_FlavQCD+Jet2_2p5_jesup_FlavQCD).Phi());

                }

		pt_jesup_split_FlavQCD.clear();  

/////////////////
// FlavQCD dn start
                TauC_Inc_0j_EnergyWgt_jesdn_FlavQCD = TauC(pt_jesdn_split_FlavQCD, eta_jes_split_FlavQCD,phi_jes_split_FlavQCD,mass_jes_split_FlavQCD, Higgs);
                TauB_Inc_0j_pTWgt_jesdn_FlavQCD = TauB(pt_jesdn_split_FlavQCD, eta_jes_split_FlavQCD,phi_jes_split_FlavQCD,mass_jes_split_FlavQCD, Higgs);


                for( unsigned int k = 0; k< pt_jesdn_split_FlavQCD.size(); k++) {
                    if (pt_jesdn_split_FlavQCD[k]<30.0 || abs(eta_jes_split_FlavQCD[k])>4.7) continue;
                    TLorentzVector thisJet_jesdn_FlavQCD;
		    thisJet_jesdn_FlavQCD.SetPtEtaPhiM(pt_jesdn_split_FlavQCD[k],eta_jes_split_FlavQCD[k],phi_jes_split_FlavQCD[k],mass_jes_split_FlavQCD[k]);

                    bool isclean_H4l=true;
                        njets_pt30_eta4p7_jesdn_FlavQCD+=1;
			if (thisJet_jesdn_FlavQCD.Pt()>jet1pt_jesdn_FlavQCD) {
                            jet2pt_jesdn_FlavQCD=jet1pt_jesdn_FlavQCD; jet2index_jesdn_FlavQCD=jet1index_jesdn_FlavQCD;
                            jet1pt_jesdn_FlavQCD=thisJet_jesdn_FlavQCD.Pt(); jet1index_jesdn_FlavQCD=k;
                        } else if (thisJet_jesdn_FlavQCD.Pt()>jet2pt_jesdn_FlavQCD) {
                            jet2pt_jesdn_FlavQCD=thisJet_jesdn_FlavQCD.Pt(); jet2index_jesdn_FlavQCD=k;
                        }
                        if (abs(thisJet_jesdn_FlavQCD.Eta())<2.5) {
                            njets_pt30_eta2p5_jesdn_FlavQCD+=1;
                            if (thisJet_jesdn_FlavQCD.Pt()>jet1pt2p5_jesdn_FlavQCD) {
                                jet2pt2p5_jesdn_FlavQCD=jet1pt2p5_jesdn_FlavQCD; jet2index2p5_jesdn_FlavQCD=jet1index2p5_jesdn_FlavQCD;
                                jet1pt2p5_jesdn_FlavQCD=thisJet_jesdn_FlavQCD.Pt(); jet1index2p5_jesdn_FlavQCD=k;
                            } else if (thisJet_jesdn_FlavQCD.Pt()>jet2pt2p5_jesdn_FlavQCD) {
                                jet2pt2p5_jesdn_FlavQCD=thisJet_jesdn_FlavQCD.Pt(); jet2index2p5_jesdn_FlavQCD=k;
                            }
                        }
                }

// Filling FlavQCD dn variables

                TLorentzVector Jet1_jesdn_FlavQCD, Jet1_2p5_jesdn_FlavQCD, Jet2_jesdn_FlavQCD, Jet2_2p5_jesdn_FlavQCD;

                if (njets_pt30_eta4p7_jesdn_FlavQCD > 0) { 
                    Jet1_jesdn_FlavQCD.SetPtEtaPhiM(pt_jesdn_split_FlavQCD[jet1index_jesdn_FlavQCD],eta_jes_split_FlavQCD[jet1index_jesdn_FlavQCD],phi_jes_split_FlavQCD[jet1index_jesdn_FlavQCD], mass_jes_split_FlavQCD[jet1index_jesdn_FlavQCD]);
                    pt_leadingjet_pt30_eta4p7_jesdn_FlavQCD=Jet1_jesdn_FlavQCD.Pt(); 
                    pTj1_jesdn_FlavQCD=Jet1_jesdn_FlavQCD.Pt(); 
                    etaj1_jesdn_FlavQCD=Jet1_jesdn_FlavQCD.Eta();
                    yj1_jesdn_FlavQCD=Jet1_jesdn_FlavQCD.Rapidity();
                    pT4lj_jesdn_FlavQCD=(Higgs+Jet1_jesdn_FlavQCD).Pt();
                    mass4lj_jesdn_FlavQCD=(Higgs+Jet1_jesdn_FlavQCD).M();
                    dPhiHj1_jesdn_FlavQCD=deltaPhi(Higgs.Phi(),Jet1_jesdn_FlavQCD.Phi());
                    dyHj1_jesdn_FlavQCD=TMath::Abs(rapidity4l-yj1_jesdn_FlavQCD);
                }    
                if (njets_pt30_eta4p7_jesdn_FlavQCD > 1) { 
                    Jet2_jesdn_FlavQCD.SetPtEtaPhiM(pt_jesdn_split_FlavQCD[jet2index_jesdn_FlavQCD],eta_jes_split_FlavQCD[jet2index_jesdn_FlavQCD],phi_jes_split_FlavQCD[jet2index_jesdn_FlavQCD], mass_jes_split_FlavQCD[jet2index_jesdn_FlavQCD]);

                    pTj2_jesdn_FlavQCD=Jet2_jesdn_FlavQCD.Pt();
                    etaj2_jesdn_FlavQCD=Jet2_jesdn_FlavQCD.Eta();
                    yj2_jesdn_FlavQCD=Jet2_jesdn_FlavQCD.Rapidity();
                    pT4ljj_jesdn_FlavQCD=(Higgs+Jet1_jesdn_FlavQCD+Jet2_jesdn_FlavQCD).Pt();
                    mass4ljj_jesdn_FlavQCD=(Higgs+Jet1_jesdn_FlavQCD+Jet2_jesdn_FlavQCD).M();
                    dEtaj1j2_jesdn_FlavQCD=TMath::Abs(Jet1_jesdn_FlavQCD.Eta()-Jet2_jesdn_FlavQCD.Eta());
                    dPhij1j2_jesdn_FlavQCD=deltaPhi(Jet1_jesdn_FlavQCD.Phi(),Jet2_jesdn_FlavQCD.Phi());
                    dPhiHj1j2_jesdn_FlavQCD=deltaPhi(Higgs.Phi(),(Jet1_jesdn_FlavQCD+Jet2_jesdn_FlavQCD).Phi());
                }    
                if (njets_pt30_eta2p5_jesdn_FlavQCD > 0) { 
                    Jet1_2p5_jesdn_FlavQCD.SetPtEtaPhiM(pt_jesdn_split_FlavQCD[jet1index2p5_jesdn_FlavQCD],eta_jes_split_FlavQCD[jet1index2p5_jesdn_FlavQCD],phi_jes_split_FlavQCD[jet1index2p5_jesdn_FlavQCD], mass_jes_split_FlavQCD[jet1index2p5_jesdn_FlavQCD]);
                    pt_leadingjet_pt30_eta2p5_jesdn_FlavQCD=Jet1_2p5_jesdn_FlavQCD.Pt();
                    pTj1_2p5_jesdn_FlavQCD=Jet1_2p5_jesdn_FlavQCD.Pt();
                    yj1_2p5_jesdn_FlavQCD=Jet1_2p5_jesdn_FlavQCD.Rapidity();
                    pT4lj_2p5_jesdn_FlavQCD=(Higgs+Jet1_2p5_jesdn_FlavQCD).Pt();
                    mass4lj_2p5_jesdn_FlavQCD=(Higgs+Jet1_2p5_jesdn_FlavQCD).M();
		    dPhiHj1_2p5_jesdn_FlavQCD=deltaPhi(Higgs.Phi(),Jet1_2p5_jesdn_FlavQCD.Phi());
                    dyHj1_2p5_jesdn_FlavQCD=TMath::Abs(rapidity4l-yj1_2p5_jesdn_FlavQCD);
                }    
                if (njets_pt30_eta2p5_jesdn_FlavQCD > 1) { 
                    Jet2_2p5_jesdn_FlavQCD.SetPtEtaPhiM(pt_jesdn_split_FlavQCD[jet2index2p5_jesdn_FlavQCD],eta_jes_split_FlavQCD[jet2index2p5_jesdn_FlavQCD],phi_jes_split_FlavQCD[jet2index2p5_jesdn_FlavQCD], mass_jes_split_FlavQCD[jet2index2p5_jesdn_FlavQCD]);
                    pTj2_2p5_jesdn_FlavQCD=Jet2_2p5_jesdn_FlavQCD.Pt();
                    yj2_2p5_jesdn_FlavQCD=Jet2_2p5_jesdn_FlavQCD.Rapidity();
                    pT4ljj_2p5_jesdn_FlavQCD=(Higgs+Jet1_2p5_jesdn_FlavQCD+Jet2_2p5_jesdn_FlavQCD).Pt();
                    mass4ljj_2p5_jesdn_FlavQCD=(Higgs+Jet1_2p5_jesdn_FlavQCD+Jet2_2p5_jesdn_FlavQCD).M();
                    mj1j2_2p5_jesdn_FlavQCD=(Jet1_2p5_jesdn_FlavQCD+Jet2_2p5_jesdn_FlavQCD).M();     
                    dEtaj1j2_2p5_jesdn_FlavQCD=TMath::Abs(Jet1_2p5_jesdn_FlavQCD.Eta()-Jet2_2p5_jesdn_FlavQCD.Eta());
	            dPhij1j2_2p5_jesdn_FlavQCD=deltaPhi(Jet1_2p5_jesdn_FlavQCD.Phi(),Jet2_2p5_jesdn_FlavQCD.Phi());
                    dPhiHj1j2_2p5_jesdn_FlavQCD=deltaPhi(Higgs.Phi(),(Jet1_2p5_jesdn_FlavQCD+Jet2_2p5_jesdn_FlavQCD).Phi());

                }    

                pt_jesdn_split_FlavQCD.clear();

		eta_jes_split_FlavQCD.clear();
		phi_jes_split_FlavQCD.clear();
		mass_jes_split_FlavQCD.clear();


////////////////////////////////////////////////////////////////////////////////
////    Starting "HF" JES uncertainty source
////////////////////////////////////////////////////////////////////////////////

                TauC_Inc_0j_EnergyWgt_jesup_HF = TauC(pt_jesup_split_HF, eta_jes_split_HF,phi_jes_split_HF,mass_jes_split_HF, Higgs);
                TauB_Inc_0j_pTWgt_jesup_HF = TauB(pt_jesup_split_HF, eta_jes_split_HF,phi_jes_split_HF,mass_jes_split_HF, Higgs);


                for( unsigned int k = 0; k< pt_jesup_split_HF.size(); k++) {
                    if (pt_jesup_split_HF[k]<30.0 || abs(eta_jes_split_HF[k])>4.7) continue;
		    TLorentzVector thisJet_jesup_HF;
		    thisJet_jesup_HF.SetPtEtaPhiM(pt_jesup_split_HF[k],eta_jes_split_HF[k],phi_jes_split_HF[k],mass_jes_split_HF[k]);

                        njets_pt30_eta4p7_jesup_HF+=1;  

                        if (thisJet_jesup_HF.Pt()>jet1pt_jesup_HF) {
                            jet2pt_jesup_HF=jet1pt_jesup_HF; jet2index_jesup_HF=jet1index_jesup_HF;
                            jet1pt_jesup_HF=thisJet_jesup_HF.Pt(); jet1index_jesup_HF=k;
                        } else if (thisJet_jesup_HF.Pt()>jet2pt_jesup_HF) {
                            jet2pt_jesup_HF=thisJet_jesup_HF.Pt(); jet2index_jesup_HF=k;
                        }
                        if (abs(thisJet_jesup_HF.Eta())<2.5) {
                            njets_pt30_eta2p5_jesup_HF+=1;
                            if (thisJet_jesup_HF.Pt()>jet1pt2p5_jesup_HF) {
                                jet2pt2p5_jesup_HF=jet1pt2p5_jesup_HF; jet2index2p5_jesup_HF=jet1index2p5_jesup_HF;
                                jet1pt2p5_jesup_HF=thisJet_jesup_HF.Pt(); jet1index2p5_jesup_HF=k;
                            } else if (thisJet_jesup_HF.Pt()>jet2pt2p5_jesup_HF) {
                                jet2pt2p5_jesup_HF=thisJet_jesup_HF.Pt(); jet2index2p5_jesup_HF=k;
                            }
                        }
                }
                cout<<njets_pt30_eta4p7_jesup_HF<<" jets (jesup_HF)"<<endl;


		TLorentzVector Jet1_jesup_HF, Jet1_2p5_jesup_HF, Jet2_jesup_HF, Jet2_2p5_jesup_HF;
                if (njets_pt30_eta4p7_jesup_HF > 0) {
		    Jet1_jesup_HF.SetPtEtaPhiM(pt_jesup_split_HF[jet1index_jesup_HF],eta_jes_split_HF[jet1index_jesup_HF],phi_jes_split_HF[jet1index_jesup_HF], mass_jes_split_HF[jet1index_jesup_HF]);

                    pt_leadingjet_pt30_eta4p7_jesup_HF=Jet1_jesup_HF.Pt(); 
                    pTj1_jesup_HF=Jet1_jesup_HF.Pt(); 
                    etaj1_jesup_HF=Jet1_jesup_HF.Eta();
                    yj1_jesup_HF=Jet1_jesup_HF.Rapidity();
		    pT4lj_jesup_HF=(Higgs+Jet1_jesup_HF).Pt();
		    mass4lj_jesup_HF=(Higgs+Jet1_jesup_HF).M();
		    dPhiHj1_jesup_HF=deltaPhi(Higgs.Phi(),Jet1_jesup_HF.Phi());
                    dyHj1_jesup_HF=TMath::Abs(rapidity4l-yj1_jesup_HF);
                }
                if (njets_pt30_eta4p7_jesup_HF > 1) {
                    Jet2_jesup_HF.SetPtEtaPhiM(pt_jesup_split_HF[jet2index_jesup_HF],eta_jes_split_HF[jet2index_jesup_HF],phi_jes_split_HF[jet2index_jesup_HF], mass_jes_split_HF[jet2index_jesup_HF]);

                    pTj2_jesup_HF=Jet2_jesup_HF.Pt();
                    etaj2_jesup_HF=Jet2_jesup_HF.Eta();
                    yj2_jesup_HF=Jet2_jesup_HF.Rapidity();
                    pT4ljj_jesup_HF=(Higgs+Jet1_jesup_HF+Jet2_jesup_HF).Pt();
                    mass4ljj_jesup_HF=(Higgs+Jet1_jesup_HF+Jet2_jesup_HF).M();
		    mj1j2_jesup_HF=(Jet1_jesup_HF+Jet2_jesup_HF).M();
                    dEtaj1j2_jesup_HF=TMath::Abs(Jet1_jesup_HF.Eta()-Jet2_jesup_HF.Eta());
                    dPhij1j2_jesup_HF=deltaPhi(Jet1_jesup_HF.Phi(),Jet2_jesup_HF.Phi());
                    dPhiHj1j2_jesup_HF=deltaPhi(Higgs.Phi(),(Jet1_jesup_HF+Jet2_jesup_HF).Phi());
                }
                if (njets_pt30_eta2p5_jesup_HF > 0) {
                    Jet1_2p5_jesup_HF.SetPtEtaPhiM(pt_jesup_split_HF[jet1index2p5_jesup_HF],eta_jes_split_HF[jet1index2p5_jesup_HF],phi_jes_split_HF[jet1index2p5_jesup_HF], mass_jes_split_HF[jet1index2p5_jesup_HF]);
                    pt_leadingjet_pt30_eta2p5_jesup_HF=Jet1_2p5_jesup_HF.Pt();
                    pTj1_2p5_jesup_HF=Jet1_2p5_jesup_HF.Pt();
                    yj1_2p5_jesup_HF=Jet1_2p5_jesup_HF.Rapidity();
                    pT4lj_2p5_jesup_HF=(Higgs+Jet1_2p5_jesup_HF).Pt();
                    mass4lj_2p5_jesup_HF=(Higgs+Jet1_2p5_jesup_HF).M();
		    dPhiHj1_2p5_jesup_HF=deltaPhi(Higgs.Phi(),Jet1_2p5_jesup_HF.Phi());
                    dyHj1_2p5_jesup_HF=TMath::Abs(rapidity4l-yj1_2p5_jesup_HF);
                }
                if (njets_pt30_eta2p5_jesup_HF > 1) {
                    Jet2_2p5_jesup_HF.SetPtEtaPhiM(pt_jesup_split_HF[jet2index2p5_jesup_HF],eta_jes_split_HF[jet2index2p5_jesup_HF],phi_jes_split_HF[jet2index2p5_jesup_HF], mass_jes_split_HF[jet2index2p5_jesup_HF]);
                    pTj2_2p5_jesup_HF=Jet2_2p5_jesup_HF.Pt();
                    yj2_2p5_jesup_HF=Jet2_2p5_jesup_HF.Rapidity();
                    pT4ljj_2p5_jesup_HF=(Higgs+Jet1_2p5_jesup_HF+Jet2_2p5_jesup_HF).Pt();
                    mass4ljj_2p5_jesup_HF=(Higgs+Jet1_2p5_jesup_HF+Jet2_2p5_jesup_HF).M();
                    mj1j2_2p5_jesup_HF=(Jet1_2p5_jesup_HF+Jet2_2p5_jesup_HF).M();                    
                    dEtaj1j2_2p5_jesup_HF=TMath::Abs(Jet1_2p5_jesup_HF.Eta()-Jet2_2p5_jesup_HF.Eta());
                    dPhij1j2_2p5_jesup_HF=deltaPhi(Jet1_2p5_jesup_HF.Phi(),Jet2_2p5_jesup_HF.Phi());
                    dPhiHj1j2_2p5_jesup_HF=deltaPhi(Higgs.Phi(),(Jet1_2p5_jesup_HF+Jet2_2p5_jesup_HF).Phi());

                }

		pt_jesup_split_HF.clear();  

/////////////////
// HF dn start
                TauC_Inc_0j_EnergyWgt_jesdn_HF = TauC(pt_jesdn_split_HF, eta_jes_split_HF,phi_jes_split_HF,mass_jes_split_HF, Higgs);
                TauB_Inc_0j_pTWgt_jesdn_HF = TauB(pt_jesdn_split_HF, eta_jes_split_HF,phi_jes_split_HF,mass_jes_split_HF, Higgs);


                for( unsigned int k = 0; k< pt_jesdn_split_HF.size(); k++) {
                    if (pt_jesdn_split_HF[k]<30.0 || abs(eta_jes_split_HF[k])>4.7) continue;
                    TLorentzVector thisJet_jesdn_HF;
		    thisJet_jesdn_HF.SetPtEtaPhiM(pt_jesdn_split_HF[k],eta_jes_split_HF[k],phi_jes_split_HF[k],mass_jes_split_HF[k]);

                    bool isclean_H4l=true;
                        njets_pt30_eta4p7_jesdn_HF+=1;
			if (thisJet_jesdn_HF.Pt()>jet1pt_jesdn_HF) {
                            jet2pt_jesdn_HF=jet1pt_jesdn_HF; jet2index_jesdn_HF=jet1index_jesdn_HF;
                            jet1pt_jesdn_HF=thisJet_jesdn_HF.Pt(); jet1index_jesdn_HF=k;
                        } else if (thisJet_jesdn_HF.Pt()>jet2pt_jesdn_HF) {
                            jet2pt_jesdn_HF=thisJet_jesdn_HF.Pt(); jet2index_jesdn_HF=k;
                        }
                        if (abs(thisJet_jesdn_HF.Eta())<2.5) {
                            njets_pt30_eta2p5_jesdn_HF+=1;
                            if (thisJet_jesdn_HF.Pt()>jet1pt2p5_jesdn_HF) {
                                jet2pt2p5_jesdn_HF=jet1pt2p5_jesdn_HF; jet2index2p5_jesdn_HF=jet1index2p5_jesdn_HF;
                                jet1pt2p5_jesdn_HF=thisJet_jesdn_HF.Pt(); jet1index2p5_jesdn_HF=k;
                            } else if (thisJet_jesdn_HF.Pt()>jet2pt2p5_jesdn_HF) {
                                jet2pt2p5_jesdn_HF=thisJet_jesdn_HF.Pt(); jet2index2p5_jesdn_HF=k;
                            }
                        }
                }

// Filling HF dn variables

                TLorentzVector Jet1_jesdn_HF, Jet1_2p5_jesdn_HF, Jet2_jesdn_HF, Jet2_2p5_jesdn_HF;

                if (njets_pt30_eta4p7_jesdn_HF > 0) { 
                    Jet1_jesdn_HF.SetPtEtaPhiM(pt_jesdn_split_HF[jet1index_jesdn_HF],eta_jes_split_HF[jet1index_jesdn_HF],phi_jes_split_HF[jet1index_jesdn_HF], mass_jes_split_HF[jet1index_jesdn_HF]);
                    pt_leadingjet_pt30_eta4p7_jesdn_HF=Jet1_jesdn_HF.Pt(); 
                    pTj1_jesdn_HF=Jet1_jesdn_HF.Pt(); 
                    etaj1_jesdn_HF=Jet1_jesdn_HF.Eta();
                    yj1_jesdn_HF=Jet1_jesdn_HF.Rapidity();
                    pT4lj_jesdn_HF=(Higgs+Jet1_jesdn_HF).Pt();
                    mass4lj_jesdn_HF=(Higgs+Jet1_jesdn_HF).M();
                    dPhiHj1_jesdn_HF=deltaPhi(Higgs.Phi(),Jet1_jesdn_HF.Phi());
                    dyHj1_jesdn_HF=TMath::Abs(rapidity4l-yj1_jesdn_HF);
                }    
                if (njets_pt30_eta4p7_jesdn_HF > 1) { 
                    Jet2_jesdn_HF.SetPtEtaPhiM(pt_jesdn_split_HF[jet2index_jesdn_HF],eta_jes_split_HF[jet2index_jesdn_HF],phi_jes_split_HF[jet2index_jesdn_HF], mass_jes_split_HF[jet2index_jesdn_HF]);

                    pTj2_jesdn_HF=Jet2_jesdn_HF.Pt();
                    etaj2_jesdn_HF=Jet2_jesdn_HF.Eta();
                    yj2_jesdn_HF=Jet2_jesdn_HF.Rapidity();
                    pT4ljj_jesdn_HF=(Higgs+Jet1_jesdn_HF+Jet2_jesdn_HF).Pt();
                    mass4ljj_jesdn_HF=(Higgs+Jet1_jesdn_HF+Jet2_jesdn_HF).M();
                    dEtaj1j2_jesdn_HF=TMath::Abs(Jet1_jesdn_HF.Eta()-Jet2_jesdn_HF.Eta());
                    dPhij1j2_jesdn_HF=deltaPhi(Jet1_jesdn_HF.Phi(),Jet2_jesdn_HF.Phi());
                    dPhiHj1j2_jesdn_HF=deltaPhi(Higgs.Phi(),(Jet1_jesdn_HF+Jet2_jesdn_HF).Phi());
                }    
                if (njets_pt30_eta2p5_jesdn_HF > 0) { 
                    Jet1_2p5_jesdn_HF.SetPtEtaPhiM(pt_jesdn_split_HF[jet1index2p5_jesdn_HF],eta_jes_split_HF[jet1index2p5_jesdn_HF],phi_jes_split_HF[jet1index2p5_jesdn_HF], mass_jes_split_HF[jet1index2p5_jesdn_HF]);
                    pt_leadingjet_pt30_eta2p5_jesdn_HF=Jet1_2p5_jesdn_HF.Pt();
                    pTj1_2p5_jesdn_HF=Jet1_2p5_jesdn_HF.Pt();
                    yj1_2p5_jesdn_HF=Jet1_2p5_jesdn_HF.Rapidity();
                    pT4lj_2p5_jesdn_HF=(Higgs+Jet1_2p5_jesdn_HF).Pt();
                    mass4lj_2p5_jesdn_HF=(Higgs+Jet1_2p5_jesdn_HF).M();
		    dPhiHj1_2p5_jesdn_HF=deltaPhi(Higgs.Phi(),Jet1_2p5_jesdn_HF.Phi());
                    dyHj1_2p5_jesdn_HF=TMath::Abs(rapidity4l-yj1_2p5_jesdn_HF);
                }    
                if (njets_pt30_eta2p5_jesdn_HF > 1) { 
                    Jet2_2p5_jesdn_HF.SetPtEtaPhiM(pt_jesdn_split_HF[jet2index2p5_jesdn_HF],eta_jes_split_HF[jet2index2p5_jesdn_HF],phi_jes_split_HF[jet2index2p5_jesdn_HF], mass_jes_split_HF[jet2index2p5_jesdn_HF]);
                    pTj2_2p5_jesdn_HF=Jet2_2p5_jesdn_HF.Pt();
                    yj2_2p5_jesdn_HF=Jet2_2p5_jesdn_HF.Rapidity();
                    pT4ljj_2p5_jesdn_HF=(Higgs+Jet1_2p5_jesdn_HF+Jet2_2p5_jesdn_HF).Pt();
                    mass4ljj_2p5_jesdn_HF=(Higgs+Jet1_2p5_jesdn_HF+Jet2_2p5_jesdn_HF).M();
                    mj1j2_2p5_jesdn_HF=(Jet1_2p5_jesdn_HF+Jet2_2p5_jesdn_HF).M();     
                    dEtaj1j2_2p5_jesdn_HF=TMath::Abs(Jet1_2p5_jesdn_HF.Eta()-Jet2_2p5_jesdn_HF.Eta());
	            dPhij1j2_2p5_jesdn_HF=deltaPhi(Jet1_2p5_jesdn_HF.Phi(),Jet2_2p5_jesdn_HF.Phi());
                    dPhiHj1j2_2p5_jesdn_HF=deltaPhi(Higgs.Phi(),(Jet1_2p5_jesdn_HF+Jet2_2p5_jesdn_HF).Phi());

                }    

                pt_jesdn_split_HF.clear();

		eta_jes_split_HF.clear();
		phi_jes_split_HF.clear();
		mass_jes_split_HF.clear();


////////////////////////////////////////////////////////////////////////////////
////    Starting "HF_year" JES uncertainty source
////////////////////////////////////////////////////////////////////////////////

                TauC_Inc_0j_EnergyWgt_jesup_HF_year = TauC(pt_jesup_split_HF_year, eta_jes_split_HF_year,phi_jes_split_HF_year,mass_jes_split_HF_year, Higgs);
                TauB_Inc_0j_pTWgt_jesup_HF_year = TauB(pt_jesup_split_HF_year, eta_jes_split_HF_year,phi_jes_split_HF_year,mass_jes_split_HF_year, Higgs);


                for( unsigned int k = 0; k< pt_jesup_split_HF_year.size(); k++) {
                    if (pt_jesup_split_HF_year[k]<30.0 || abs(eta_jes_split_HF_year[k])>4.7) continue;
		    TLorentzVector thisJet_jesup_HF_year;
		    thisJet_jesup_HF_year.SetPtEtaPhiM(pt_jesup_split_HF_year[k],eta_jes_split_HF_year[k],phi_jes_split_HF_year[k],mass_jes_split_HF_year[k]);

                        njets_pt30_eta4p7_jesup_HF_year+=1;  

                        if (thisJet_jesup_HF_year.Pt()>jet1pt_jesup_HF_year) {
                            jet2pt_jesup_HF_year=jet1pt_jesup_HF_year; jet2index_jesup_HF_year=jet1index_jesup_HF_year;
                            jet1pt_jesup_HF_year=thisJet_jesup_HF_year.Pt(); jet1index_jesup_HF_year=k;
                        } else if (thisJet_jesup_HF_year.Pt()>jet2pt_jesup_HF_year) {
                            jet2pt_jesup_HF_year=thisJet_jesup_HF_year.Pt(); jet2index_jesup_HF_year=k;
                        }
                        if (abs(thisJet_jesup_HF_year.Eta())<2.5) {
                            njets_pt30_eta2p5_jesup_HF_year+=1;
                            if (thisJet_jesup_HF_year.Pt()>jet1pt2p5_jesup_HF_year) {
                                jet2pt2p5_jesup_HF_year=jet1pt2p5_jesup_HF_year; jet2index2p5_jesup_HF_year=jet1index2p5_jesup_HF_year;
                                jet1pt2p5_jesup_HF_year=thisJet_jesup_HF_year.Pt(); jet1index2p5_jesup_HF_year=k;
                            } else if (thisJet_jesup_HF_year.Pt()>jet2pt2p5_jesup_HF_year) {
                                jet2pt2p5_jesup_HF_year=thisJet_jesup_HF_year.Pt(); jet2index2p5_jesup_HF_year=k;
                            }
                        }
                }
                cout<<njets_pt30_eta4p7_jesup_HF_year<<" jets (jesup_HF_year)"<<endl;


		TLorentzVector Jet1_jesup_HF_year, Jet1_2p5_jesup_HF_year, Jet2_jesup_HF_year, Jet2_2p5_jesup_HF_year;
                if (njets_pt30_eta4p7_jesup_HF_year > 0) {
		    Jet1_jesup_HF_year.SetPtEtaPhiM(pt_jesup_split_HF_year[jet1index_jesup_HF_year],eta_jes_split_HF_year[jet1index_jesup_HF_year],phi_jes_split_HF_year[jet1index_jesup_HF_year], mass_jes_split_HF_year[jet1index_jesup_HF_year]);

                    pt_leadingjet_pt30_eta4p7_jesup_HF_year=Jet1_jesup_HF_year.Pt(); 
                    pTj1_jesup_HF_year=Jet1_jesup_HF_year.Pt(); 
                    etaj1_jesup_HF_year=Jet1_jesup_HF_year.Eta();
                    yj1_jesup_HF_year=Jet1_jesup_HF_year.Rapidity();
		    pT4lj_jesup_HF_year=(Higgs+Jet1_jesup_HF_year).Pt();
		    mass4lj_jesup_HF_year=(Higgs+Jet1_jesup_HF_year).M();
		    dPhiHj1_jesup_HF_year=deltaPhi(Higgs.Phi(),Jet1_jesup_HF_year.Phi());
                    dyHj1_jesup_HF_year=TMath::Abs(rapidity4l-yj1_jesup_HF_year);
                }
                if (njets_pt30_eta4p7_jesup_HF_year > 1) {
                    Jet2_jesup_HF_year.SetPtEtaPhiM(pt_jesup_split_HF_year[jet2index_jesup_HF_year],eta_jes_split_HF_year[jet2index_jesup_HF_year],phi_jes_split_HF_year[jet2index_jesup_HF_year], mass_jes_split_HF_year[jet2index_jesup_HF_year]);

                    pTj2_jesup_HF_year=Jet2_jesup_HF_year.Pt();
                    etaj2_jesup_HF_year=Jet2_jesup_HF_year.Eta();
                    yj2_jesup_HF_year=Jet2_jesup_HF_year.Rapidity();
                    pT4ljj_jesup_HF_year=(Higgs+Jet1_jesup_HF_year+Jet2_jesup_HF_year).Pt();
                    mass4ljj_jesup_HF_year=(Higgs+Jet1_jesup_HF_year+Jet2_jesup_HF_year).M();
		    mj1j2_jesup_HF_year=(Jet1_jesup_HF_year+Jet2_jesup_HF_year).M();
                    dEtaj1j2_jesup_HF_year=TMath::Abs(Jet1_jesup_HF_year.Eta()-Jet2_jesup_HF_year.Eta());
                    dPhij1j2_jesup_HF_year=deltaPhi(Jet1_jesup_HF_year.Phi(),Jet2_jesup_HF_year.Phi());
                    dPhiHj1j2_jesup_HF_year=deltaPhi(Higgs.Phi(),(Jet1_jesup_HF_year+Jet2_jesup_HF_year).Phi());
                }
                if (njets_pt30_eta2p5_jesup_HF_year > 0) {
                    Jet1_2p5_jesup_HF_year.SetPtEtaPhiM(pt_jesup_split_HF_year[jet1index2p5_jesup_HF_year],eta_jes_split_HF_year[jet1index2p5_jesup_HF_year],phi_jes_split_HF_year[jet1index2p5_jesup_HF_year], mass_jes_split_HF_year[jet1index2p5_jesup_HF_year]);
                    pt_leadingjet_pt30_eta2p5_jesup_HF_year=Jet1_2p5_jesup_HF_year.Pt();
                    pTj1_2p5_jesup_HF_year=Jet1_2p5_jesup_HF_year.Pt();
                    yj1_2p5_jesup_HF_year=Jet1_2p5_jesup_HF_year.Rapidity();
                    pT4lj_2p5_jesup_HF_year=(Higgs+Jet1_2p5_jesup_HF_year).Pt();
                    mass4lj_2p5_jesup_HF_year=(Higgs+Jet1_2p5_jesup_HF_year).M();
		    dPhiHj1_2p5_jesup_HF_year=deltaPhi(Higgs.Phi(),Jet1_2p5_jesup_HF_year.Phi());
                    dyHj1_2p5_jesup_HF_year=TMath::Abs(rapidity4l-yj1_2p5_jesup_HF_year);
                }
                if (njets_pt30_eta2p5_jesup_HF_year > 1) {
                    Jet2_2p5_jesup_HF_year.SetPtEtaPhiM(pt_jesup_split_HF_year[jet2index2p5_jesup_HF_year],eta_jes_split_HF_year[jet2index2p5_jesup_HF_year],phi_jes_split_HF_year[jet2index2p5_jesup_HF_year], mass_jes_split_HF_year[jet2index2p5_jesup_HF_year]);
                    pTj2_2p5_jesup_HF_year=Jet2_2p5_jesup_HF_year.Pt();
                    yj2_2p5_jesup_HF_year=Jet2_2p5_jesup_HF_year.Rapidity();
                    pT4ljj_2p5_jesup_HF_year=(Higgs+Jet1_2p5_jesup_HF_year+Jet2_2p5_jesup_HF_year).Pt();
                    mass4ljj_2p5_jesup_HF_year=(Higgs+Jet1_2p5_jesup_HF_year+Jet2_2p5_jesup_HF_year).M();
                    mj1j2_2p5_jesup_HF_year=(Jet1_2p5_jesup_HF_year+Jet2_2p5_jesup_HF_year).M();                    
                    dEtaj1j2_2p5_jesup_HF_year=TMath::Abs(Jet1_2p5_jesup_HF_year.Eta()-Jet2_2p5_jesup_HF_year.Eta());
                    dPhij1j2_2p5_jesup_HF_year=deltaPhi(Jet1_2p5_jesup_HF_year.Phi(),Jet2_2p5_jesup_HF_year.Phi());
                    dPhiHj1j2_2p5_jesup_HF_year=deltaPhi(Higgs.Phi(),(Jet1_2p5_jesup_HF_year+Jet2_2p5_jesup_HF_year).Phi());

                }

		pt_jesup_split_HF_year.clear();  

/////////////////
// HF_year dn start
                TauC_Inc_0j_EnergyWgt_jesdn_HF_year = TauC(pt_jesdn_split_HF_year, eta_jes_split_HF_year,phi_jes_split_HF_year,mass_jes_split_HF_year, Higgs);
                TauB_Inc_0j_pTWgt_jesdn_HF_year = TauB(pt_jesdn_split_HF_year, eta_jes_split_HF_year,phi_jes_split_HF_year,mass_jes_split_HF_year, Higgs);


                for( unsigned int k = 0; k< pt_jesdn_split_HF_year.size(); k++) {
                    if (pt_jesdn_split_HF_year[k]<30.0 || abs(eta_jes_split_HF_year[k])>4.7) continue;
                    TLorentzVector thisJet_jesdn_HF_year;
		    thisJet_jesdn_HF_year.SetPtEtaPhiM(pt_jesdn_split_HF_year[k],eta_jes_split_HF_year[k],phi_jes_split_HF_year[k],mass_jes_split_HF_year[k]);

                    bool isclean_H4l=true;
                        njets_pt30_eta4p7_jesdn_HF_year+=1;
			if (thisJet_jesdn_HF_year.Pt()>jet1pt_jesdn_HF_year) {
                            jet2pt_jesdn_HF_year=jet1pt_jesdn_HF_year; jet2index_jesdn_HF_year=jet1index_jesdn_HF_year;
                            jet1pt_jesdn_HF_year=thisJet_jesdn_HF_year.Pt(); jet1index_jesdn_HF_year=k;
                        } else if (thisJet_jesdn_HF_year.Pt()>jet2pt_jesdn_HF_year) {
                            jet2pt_jesdn_HF_year=thisJet_jesdn_HF_year.Pt(); jet2index_jesdn_HF_year=k;
                        }
                        if (abs(thisJet_jesdn_HF_year.Eta())<2.5) {
                            njets_pt30_eta2p5_jesdn_HF_year+=1;
                            if (thisJet_jesdn_HF_year.Pt()>jet1pt2p5_jesdn_HF_year) {
                                jet2pt2p5_jesdn_HF_year=jet1pt2p5_jesdn_HF_year; jet2index2p5_jesdn_HF_year=jet1index2p5_jesdn_HF_year;
                                jet1pt2p5_jesdn_HF_year=thisJet_jesdn_HF_year.Pt(); jet1index2p5_jesdn_HF_year=k;
                            } else if (thisJet_jesdn_HF_year.Pt()>jet2pt2p5_jesdn_HF_year) {
                                jet2pt2p5_jesdn_HF_year=thisJet_jesdn_HF_year.Pt(); jet2index2p5_jesdn_HF_year=k;
                            }
                        }
                }

// Filling HF_year dn variables

                TLorentzVector Jet1_jesdn_HF_year, Jet1_2p5_jesdn_HF_year, Jet2_jesdn_HF_year, Jet2_2p5_jesdn_HF_year;

                if (njets_pt30_eta4p7_jesdn_HF_year > 0) { 
                    Jet1_jesdn_HF_year.SetPtEtaPhiM(pt_jesdn_split_HF_year[jet1index_jesdn_HF_year],eta_jes_split_HF_year[jet1index_jesdn_HF_year],phi_jes_split_HF_year[jet1index_jesdn_HF_year], mass_jes_split_HF_year[jet1index_jesdn_HF_year]);
                    pt_leadingjet_pt30_eta4p7_jesdn_HF_year=Jet1_jesdn_HF_year.Pt(); 
                    pTj1_jesdn_HF_year=Jet1_jesdn_HF_year.Pt(); 
                    etaj1_jesdn_HF_year=Jet1_jesdn_HF_year.Eta();
                    yj1_jesdn_HF_year=Jet1_jesdn_HF_year.Rapidity();
                    pT4lj_jesdn_HF_year=(Higgs+Jet1_jesdn_HF_year).Pt();
                    mass4lj_jesdn_HF_year=(Higgs+Jet1_jesdn_HF_year).M();
                    dPhiHj1_jesdn_HF_year=deltaPhi(Higgs.Phi(),Jet1_jesdn_HF_year.Phi());
                    dyHj1_jesdn_HF_year=TMath::Abs(rapidity4l-yj1_jesdn_HF_year);
                }    
                if (njets_pt30_eta4p7_jesdn_HF_year > 1) { 
                    Jet2_jesdn_HF_year.SetPtEtaPhiM(pt_jesdn_split_HF_year[jet2index_jesdn_HF_year],eta_jes_split_HF_year[jet2index_jesdn_HF_year],phi_jes_split_HF_year[jet2index_jesdn_HF_year], mass_jes_split_HF_year[jet2index_jesdn_HF_year]);

                    pTj2_jesdn_HF_year=Jet2_jesdn_HF_year.Pt();
                    etaj2_jesdn_HF_year=Jet2_jesdn_HF_year.Eta();
                    yj2_jesdn_HF_year=Jet2_jesdn_HF_year.Rapidity();
                    pT4ljj_jesdn_HF_year=(Higgs+Jet1_jesdn_HF_year+Jet2_jesdn_HF_year).Pt();
                    mass4ljj_jesdn_HF_year=(Higgs+Jet1_jesdn_HF_year+Jet2_jesdn_HF_year).M();
                    dEtaj1j2_jesdn_HF_year=TMath::Abs(Jet1_jesdn_HF_year.Eta()-Jet2_jesdn_HF_year.Eta());
                    dPhij1j2_jesdn_HF_year=deltaPhi(Jet1_jesdn_HF_year.Phi(),Jet2_jesdn_HF_year.Phi());
                    dPhiHj1j2_jesdn_HF_year=deltaPhi(Higgs.Phi(),(Jet1_jesdn_HF_year+Jet2_jesdn_HF_year).Phi());
                }    
                if (njets_pt30_eta2p5_jesdn_HF_year > 0) { 
                    Jet1_2p5_jesdn_HF_year.SetPtEtaPhiM(pt_jesdn_split_HF_year[jet1index2p5_jesdn_HF_year],eta_jes_split_HF_year[jet1index2p5_jesdn_HF_year],phi_jes_split_HF_year[jet1index2p5_jesdn_HF_year], mass_jes_split_HF_year[jet1index2p5_jesdn_HF_year]);
                    pt_leadingjet_pt30_eta2p5_jesdn_HF_year=Jet1_2p5_jesdn_HF_year.Pt();
                    pTj1_2p5_jesdn_HF_year=Jet1_2p5_jesdn_HF_year.Pt();
                    yj1_2p5_jesdn_HF_year=Jet1_2p5_jesdn_HF_year.Rapidity();
                    pT4lj_2p5_jesdn_HF_year=(Higgs+Jet1_2p5_jesdn_HF_year).Pt();
                    mass4lj_2p5_jesdn_HF_year=(Higgs+Jet1_2p5_jesdn_HF_year).M();
		    dPhiHj1_2p5_jesdn_HF_year=deltaPhi(Higgs.Phi(),Jet1_2p5_jesdn_HF_year.Phi());
                    dyHj1_2p5_jesdn_HF_year=TMath::Abs(rapidity4l-yj1_2p5_jesdn_HF_year);
                }    
                if (njets_pt30_eta2p5_jesdn_HF_year > 1) { 
                    Jet2_2p5_jesdn_HF_year.SetPtEtaPhiM(pt_jesdn_split_HF_year[jet2index2p5_jesdn_HF_year],eta_jes_split_HF_year[jet2index2p5_jesdn_HF_year],phi_jes_split_HF_year[jet2index2p5_jesdn_HF_year], mass_jes_split_HF_year[jet2index2p5_jesdn_HF_year]);
                    pTj2_2p5_jesdn_HF_year=Jet2_2p5_jesdn_HF_year.Pt();
                    yj2_2p5_jesdn_HF_year=Jet2_2p5_jesdn_HF_year.Rapidity();
                    pT4ljj_2p5_jesdn_HF_year=(Higgs+Jet1_2p5_jesdn_HF_year+Jet2_2p5_jesdn_HF_year).Pt();
                    mass4ljj_2p5_jesdn_HF_year=(Higgs+Jet1_2p5_jesdn_HF_year+Jet2_2p5_jesdn_HF_year).M();
                    mj1j2_2p5_jesdn_HF_year=(Jet1_2p5_jesdn_HF_year+Jet2_2p5_jesdn_HF_year).M();     
                    dEtaj1j2_2p5_jesdn_HF_year=TMath::Abs(Jet1_2p5_jesdn_HF_year.Eta()-Jet2_2p5_jesdn_HF_year.Eta());
	            dPhij1j2_2p5_jesdn_HF_year=deltaPhi(Jet1_2p5_jesdn_HF_year.Phi(),Jet2_2p5_jesdn_HF_year.Phi());
                    dPhiHj1j2_2p5_jesdn_HF_year=deltaPhi(Higgs.Phi(),(Jet1_2p5_jesdn_HF_year+Jet2_2p5_jesdn_HF_year).Phi());

                }    

                pt_jesdn_split_HF_year.clear();

		eta_jes_split_HF_year.clear();
		phi_jes_split_HF_year.clear();
		mass_jes_split_HF_year.clear();


////////////////////////////////////////////////////////////////////////////////
////    Starting "RelBal" JES uncertainty source
////////////////////////////////////////////////////////////////////////////////

                TauC_Inc_0j_EnergyWgt_jesup_RelBal = TauC(pt_jesup_split_RelBal, eta_jes_split_RelBal,phi_jes_split_RelBal,mass_jes_split_RelBal, Higgs);
                TauB_Inc_0j_pTWgt_jesup_RelBal = TauB(pt_jesup_split_RelBal, eta_jes_split_RelBal,phi_jes_split_RelBal,mass_jes_split_RelBal, Higgs);


                for( unsigned int k = 0; k< pt_jesup_split_RelBal.size(); k++) {
                    if (pt_jesup_split_RelBal[k]<30.0 || abs(eta_jes_split_RelBal[k])>4.7) continue;
		    TLorentzVector thisJet_jesup_RelBal;
		    thisJet_jesup_RelBal.SetPtEtaPhiM(pt_jesup_split_RelBal[k],eta_jes_split_RelBal[k],phi_jes_split_RelBal[k],mass_jes_split_RelBal[k]);

                        njets_pt30_eta4p7_jesup_RelBal+=1;  

                        if (thisJet_jesup_RelBal.Pt()>jet1pt_jesup_RelBal) {
                            jet2pt_jesup_RelBal=jet1pt_jesup_RelBal; jet2index_jesup_RelBal=jet1index_jesup_RelBal;
                            jet1pt_jesup_RelBal=thisJet_jesup_RelBal.Pt(); jet1index_jesup_RelBal=k;
                        } else if (thisJet_jesup_RelBal.Pt()>jet2pt_jesup_RelBal) {
                            jet2pt_jesup_RelBal=thisJet_jesup_RelBal.Pt(); jet2index_jesup_RelBal=k;
                        }
                        if (abs(thisJet_jesup_RelBal.Eta())<2.5) {
                            njets_pt30_eta2p5_jesup_RelBal+=1;
                            if (thisJet_jesup_RelBal.Pt()>jet1pt2p5_jesup_RelBal) {
                                jet2pt2p5_jesup_RelBal=jet1pt2p5_jesup_RelBal; jet2index2p5_jesup_RelBal=jet1index2p5_jesup_RelBal;
                                jet1pt2p5_jesup_RelBal=thisJet_jesup_RelBal.Pt(); jet1index2p5_jesup_RelBal=k;
                            } else if (thisJet_jesup_RelBal.Pt()>jet2pt2p5_jesup_RelBal) {
                                jet2pt2p5_jesup_RelBal=thisJet_jesup_RelBal.Pt(); jet2index2p5_jesup_RelBal=k;
                            }
                        }
                }
                cout<<njets_pt30_eta4p7_jesup_RelBal<<" jets (jesup_RelBal)"<<endl;


		TLorentzVector Jet1_jesup_RelBal, Jet1_2p5_jesup_RelBal, Jet2_jesup_RelBal, Jet2_2p5_jesup_RelBal;
                if (njets_pt30_eta4p7_jesup_RelBal > 0) {
		    Jet1_jesup_RelBal.SetPtEtaPhiM(pt_jesup_split_RelBal[jet1index_jesup_RelBal],eta_jes_split_RelBal[jet1index_jesup_RelBal],phi_jes_split_RelBal[jet1index_jesup_RelBal], mass_jes_split_RelBal[jet1index_jesup_RelBal]);

                    pt_leadingjet_pt30_eta4p7_jesup_RelBal=Jet1_jesup_RelBal.Pt(); 
                    pTj1_jesup_RelBal=Jet1_jesup_RelBal.Pt(); 
                    etaj1_jesup_RelBal=Jet1_jesup_RelBal.Eta();
                    yj1_jesup_RelBal=Jet1_jesup_RelBal.Rapidity();
		    pT4lj_jesup_RelBal=(Higgs+Jet1_jesup_RelBal).Pt();
		    mass4lj_jesup_RelBal=(Higgs+Jet1_jesup_RelBal).M();
		    dPhiHj1_jesup_RelBal=deltaPhi(Higgs.Phi(),Jet1_jesup_RelBal.Phi());
                    dyHj1_jesup_RelBal=TMath::Abs(rapidity4l-yj1_jesup_RelBal);
                }
                if (njets_pt30_eta4p7_jesup_RelBal > 1) {
                    Jet2_jesup_RelBal.SetPtEtaPhiM(pt_jesup_split_RelBal[jet2index_jesup_RelBal],eta_jes_split_RelBal[jet2index_jesup_RelBal],phi_jes_split_RelBal[jet2index_jesup_RelBal], mass_jes_split_RelBal[jet2index_jesup_RelBal]);

                    pTj2_jesup_RelBal=Jet2_jesup_RelBal.Pt();
                    etaj2_jesup_RelBal=Jet2_jesup_RelBal.Eta();
                    yj2_jesup_RelBal=Jet2_jesup_RelBal.Rapidity();
                    pT4ljj_jesup_RelBal=(Higgs+Jet1_jesup_RelBal+Jet2_jesup_RelBal).Pt();
                    mass4ljj_jesup_RelBal=(Higgs+Jet1_jesup_RelBal+Jet2_jesup_RelBal).M();
		    mj1j2_jesup_RelBal=(Jet1_jesup_RelBal+Jet2_jesup_RelBal).M();
                    dEtaj1j2_jesup_RelBal=TMath::Abs(Jet1_jesup_RelBal.Eta()-Jet2_jesup_RelBal.Eta());
                    dPhij1j2_jesup_RelBal=deltaPhi(Jet1_jesup_RelBal.Phi(),Jet2_jesup_RelBal.Phi());
                    dPhiHj1j2_jesup_RelBal=deltaPhi(Higgs.Phi(),(Jet1_jesup_RelBal+Jet2_jesup_RelBal).Phi());
                }
                if (njets_pt30_eta2p5_jesup_RelBal > 0) {
                    Jet1_2p5_jesup_RelBal.SetPtEtaPhiM(pt_jesup_split_RelBal[jet1index2p5_jesup_RelBal],eta_jes_split_RelBal[jet1index2p5_jesup_RelBal],phi_jes_split_RelBal[jet1index2p5_jesup_RelBal], mass_jes_split_RelBal[jet1index2p5_jesup_RelBal]);
                    pt_leadingjet_pt30_eta2p5_jesup_RelBal=Jet1_2p5_jesup_RelBal.Pt();
                    pTj1_2p5_jesup_RelBal=Jet1_2p5_jesup_RelBal.Pt();
                    yj1_2p5_jesup_RelBal=Jet1_2p5_jesup_RelBal.Rapidity();
                    pT4lj_2p5_jesup_RelBal=(Higgs+Jet1_2p5_jesup_RelBal).Pt();
                    mass4lj_2p5_jesup_RelBal=(Higgs+Jet1_2p5_jesup_RelBal).M();
		    dPhiHj1_2p5_jesup_RelBal=deltaPhi(Higgs.Phi(),Jet1_2p5_jesup_RelBal.Phi());
                    dyHj1_2p5_jesup_RelBal=TMath::Abs(rapidity4l-yj1_2p5_jesup_RelBal);
                }
                if (njets_pt30_eta2p5_jesup_RelBal > 1) {
                    Jet2_2p5_jesup_RelBal.SetPtEtaPhiM(pt_jesup_split_RelBal[jet2index2p5_jesup_RelBal],eta_jes_split_RelBal[jet2index2p5_jesup_RelBal],phi_jes_split_RelBal[jet2index2p5_jesup_RelBal], mass_jes_split_RelBal[jet2index2p5_jesup_RelBal]);
                    pTj2_2p5_jesup_RelBal=Jet2_2p5_jesup_RelBal.Pt();
                    yj2_2p5_jesup_RelBal=Jet2_2p5_jesup_RelBal.Rapidity();
                    pT4ljj_2p5_jesup_RelBal=(Higgs+Jet1_2p5_jesup_RelBal+Jet2_2p5_jesup_RelBal).Pt();
                    mass4ljj_2p5_jesup_RelBal=(Higgs+Jet1_2p5_jesup_RelBal+Jet2_2p5_jesup_RelBal).M();
                    mj1j2_2p5_jesup_RelBal=(Jet1_2p5_jesup_RelBal+Jet2_2p5_jesup_RelBal).M();                    
                    dEtaj1j2_2p5_jesup_RelBal=TMath::Abs(Jet1_2p5_jesup_RelBal.Eta()-Jet2_2p5_jesup_RelBal.Eta());
                    dPhij1j2_2p5_jesup_RelBal=deltaPhi(Jet1_2p5_jesup_RelBal.Phi(),Jet2_2p5_jesup_RelBal.Phi());
                    dPhiHj1j2_2p5_jesup_RelBal=deltaPhi(Higgs.Phi(),(Jet1_2p5_jesup_RelBal+Jet2_2p5_jesup_RelBal).Phi());

                }

		pt_jesup_split_RelBal.clear();  

/////////////////
// RelBal dn start
                TauC_Inc_0j_EnergyWgt_jesdn_RelBal = TauC(pt_jesdn_split_RelBal, eta_jes_split_RelBal,phi_jes_split_RelBal,mass_jes_split_RelBal, Higgs);
                TauB_Inc_0j_pTWgt_jesdn_RelBal = TauB(pt_jesdn_split_RelBal, eta_jes_split_RelBal,phi_jes_split_RelBal,mass_jes_split_RelBal, Higgs);


                for( unsigned int k = 0; k< pt_jesdn_split_RelBal.size(); k++) {
                    if (pt_jesdn_split_RelBal[k]<30.0 || abs(eta_jes_split_RelBal[k])>4.7) continue;
                    TLorentzVector thisJet_jesdn_RelBal;
		    thisJet_jesdn_RelBal.SetPtEtaPhiM(pt_jesdn_split_RelBal[k],eta_jes_split_RelBal[k],phi_jes_split_RelBal[k],mass_jes_split_RelBal[k]);

                    bool isclean_H4l=true;
                        njets_pt30_eta4p7_jesdn_RelBal+=1;
			if (thisJet_jesdn_RelBal.Pt()>jet1pt_jesdn_RelBal) {
                            jet2pt_jesdn_RelBal=jet1pt_jesdn_RelBal; jet2index_jesdn_RelBal=jet1index_jesdn_RelBal;
                            jet1pt_jesdn_RelBal=thisJet_jesdn_RelBal.Pt(); jet1index_jesdn_RelBal=k;
                        } else if (thisJet_jesdn_RelBal.Pt()>jet2pt_jesdn_RelBal) {
                            jet2pt_jesdn_RelBal=thisJet_jesdn_RelBal.Pt(); jet2index_jesdn_RelBal=k;
                        }
                        if (abs(thisJet_jesdn_RelBal.Eta())<2.5) {
                            njets_pt30_eta2p5_jesdn_RelBal+=1;
                            if (thisJet_jesdn_RelBal.Pt()>jet1pt2p5_jesdn_RelBal) {
                                jet2pt2p5_jesdn_RelBal=jet1pt2p5_jesdn_RelBal; jet2index2p5_jesdn_RelBal=jet1index2p5_jesdn_RelBal;
                                jet1pt2p5_jesdn_RelBal=thisJet_jesdn_RelBal.Pt(); jet1index2p5_jesdn_RelBal=k;
                            } else if (thisJet_jesdn_RelBal.Pt()>jet2pt2p5_jesdn_RelBal) {
                                jet2pt2p5_jesdn_RelBal=thisJet_jesdn_RelBal.Pt(); jet2index2p5_jesdn_RelBal=k;
                            }
                        }
                }

// Filling RelBal dn variables

                TLorentzVector Jet1_jesdn_RelBal, Jet1_2p5_jesdn_RelBal, Jet2_jesdn_RelBal, Jet2_2p5_jesdn_RelBal;

                if (njets_pt30_eta4p7_jesdn_RelBal > 0) { 
                    Jet1_jesdn_RelBal.SetPtEtaPhiM(pt_jesdn_split_RelBal[jet1index_jesdn_RelBal],eta_jes_split_RelBal[jet1index_jesdn_RelBal],phi_jes_split_RelBal[jet1index_jesdn_RelBal], mass_jes_split_RelBal[jet1index_jesdn_RelBal]);
                    pt_leadingjet_pt30_eta4p7_jesdn_RelBal=Jet1_jesdn_RelBal.Pt(); 
                    pTj1_jesdn_RelBal=Jet1_jesdn_RelBal.Pt(); 
                    etaj1_jesdn_RelBal=Jet1_jesdn_RelBal.Eta();
                    yj1_jesdn_RelBal=Jet1_jesdn_RelBal.Rapidity();
                    pT4lj_jesdn_RelBal=(Higgs+Jet1_jesdn_RelBal).Pt();
                    mass4lj_jesdn_RelBal=(Higgs+Jet1_jesdn_RelBal).M();
                    dPhiHj1_jesdn_RelBal=deltaPhi(Higgs.Phi(),Jet1_jesdn_RelBal.Phi());
                    dyHj1_jesdn_RelBal=TMath::Abs(rapidity4l-yj1_jesdn_RelBal);
                }    
                if (njets_pt30_eta4p7_jesdn_RelBal > 1) { 
                    Jet2_jesdn_RelBal.SetPtEtaPhiM(pt_jesdn_split_RelBal[jet2index_jesdn_RelBal],eta_jes_split_RelBal[jet2index_jesdn_RelBal],phi_jes_split_RelBal[jet2index_jesdn_RelBal], mass_jes_split_RelBal[jet2index_jesdn_RelBal]);

                    pTj2_jesdn_RelBal=Jet2_jesdn_RelBal.Pt();
                    etaj2_jesdn_RelBal=Jet2_jesdn_RelBal.Eta();
                    yj2_jesdn_RelBal=Jet2_jesdn_RelBal.Rapidity();
                    pT4ljj_jesdn_RelBal=(Higgs+Jet1_jesdn_RelBal+Jet2_jesdn_RelBal).Pt();
                    mass4ljj_jesdn_RelBal=(Higgs+Jet1_jesdn_RelBal+Jet2_jesdn_RelBal).M();
                    dEtaj1j2_jesdn_RelBal=TMath::Abs(Jet1_jesdn_RelBal.Eta()-Jet2_jesdn_RelBal.Eta());
                    dPhij1j2_jesdn_RelBal=deltaPhi(Jet1_jesdn_RelBal.Phi(),Jet2_jesdn_RelBal.Phi());
                    dPhiHj1j2_jesdn_RelBal=deltaPhi(Higgs.Phi(),(Jet1_jesdn_RelBal+Jet2_jesdn_RelBal).Phi());
                }    
                if (njets_pt30_eta2p5_jesdn_RelBal > 0) { 
                    Jet1_2p5_jesdn_RelBal.SetPtEtaPhiM(pt_jesdn_split_RelBal[jet1index2p5_jesdn_RelBal],eta_jes_split_RelBal[jet1index2p5_jesdn_RelBal],phi_jes_split_RelBal[jet1index2p5_jesdn_RelBal], mass_jes_split_RelBal[jet1index2p5_jesdn_RelBal]);
                    pt_leadingjet_pt30_eta2p5_jesdn_RelBal=Jet1_2p5_jesdn_RelBal.Pt();
                    pTj1_2p5_jesdn_RelBal=Jet1_2p5_jesdn_RelBal.Pt();
                    yj1_2p5_jesdn_RelBal=Jet1_2p5_jesdn_RelBal.Rapidity();
                    pT4lj_2p5_jesdn_RelBal=(Higgs+Jet1_2p5_jesdn_RelBal).Pt();
                    mass4lj_2p5_jesdn_RelBal=(Higgs+Jet1_2p5_jesdn_RelBal).M();
		    dPhiHj1_2p5_jesdn_RelBal=deltaPhi(Higgs.Phi(),Jet1_2p5_jesdn_RelBal.Phi());
                    dyHj1_2p5_jesdn_RelBal=TMath::Abs(rapidity4l-yj1_2p5_jesdn_RelBal);
                }    
                if (njets_pt30_eta2p5_jesdn_RelBal > 1) { 
                    Jet2_2p5_jesdn_RelBal.SetPtEtaPhiM(pt_jesdn_split_RelBal[jet2index2p5_jesdn_RelBal],eta_jes_split_RelBal[jet2index2p5_jesdn_RelBal],phi_jes_split_RelBal[jet2index2p5_jesdn_RelBal], mass_jes_split_RelBal[jet2index2p5_jesdn_RelBal]);
                    pTj2_2p5_jesdn_RelBal=Jet2_2p5_jesdn_RelBal.Pt();
                    yj2_2p5_jesdn_RelBal=Jet2_2p5_jesdn_RelBal.Rapidity();
                    pT4ljj_2p5_jesdn_RelBal=(Higgs+Jet1_2p5_jesdn_RelBal+Jet2_2p5_jesdn_RelBal).Pt();
                    mass4ljj_2p5_jesdn_RelBal=(Higgs+Jet1_2p5_jesdn_RelBal+Jet2_2p5_jesdn_RelBal).M();
                    mj1j2_2p5_jesdn_RelBal=(Jet1_2p5_jesdn_RelBal+Jet2_2p5_jesdn_RelBal).M();     
                    dEtaj1j2_2p5_jesdn_RelBal=TMath::Abs(Jet1_2p5_jesdn_RelBal.Eta()-Jet2_2p5_jesdn_RelBal.Eta());
	            dPhij1j2_2p5_jesdn_RelBal=deltaPhi(Jet1_2p5_jesdn_RelBal.Phi(),Jet2_2p5_jesdn_RelBal.Phi());
                    dPhiHj1j2_2p5_jesdn_RelBal=deltaPhi(Higgs.Phi(),(Jet1_2p5_jesdn_RelBal+Jet2_2p5_jesdn_RelBal).Phi());

                }    

                pt_jesdn_split_RelBal.clear();

		eta_jes_split_RelBal.clear();
		phi_jes_split_RelBal.clear();
		mass_jes_split_RelBal.clear();


////////////////////////////////////////////////////////////////////////////////
////    Starting "RelSample_year" JES uncertainty source
////////////////////////////////////////////////////////////////////////////////

                TauC_Inc_0j_EnergyWgt_jesup_RelSample_year = TauC(pt_jesup_split_RelSample_year, eta_jes_split_RelSample_year,phi_jes_split_RelSample_year,mass_jes_split_RelSample_year, Higgs);
                TauB_Inc_0j_pTWgt_jesup_RelSample_year = TauB(pt_jesup_split_RelSample_year, eta_jes_split_RelSample_year,phi_jes_split_RelSample_year,mass_jes_split_RelSample_year, Higgs);


                for( unsigned int k = 0; k< pt_jesup_split_RelSample_year.size(); k++) {
                    if (pt_jesup_split_RelSample_year[k]<30.0 || abs(eta_jes_split_RelSample_year[k])>4.7) continue;
		    TLorentzVector thisJet_jesup_RelSample_year;
		    thisJet_jesup_RelSample_year.SetPtEtaPhiM(pt_jesup_split_RelSample_year[k],eta_jes_split_RelSample_year[k],phi_jes_split_RelSample_year[k],mass_jes_split_RelSample_year[k]);

                        njets_pt30_eta4p7_jesup_RelSample_year+=1;  

                        if (thisJet_jesup_RelSample_year.Pt()>jet1pt_jesup_RelSample_year) {
                            jet2pt_jesup_RelSample_year=jet1pt_jesup_RelSample_year; jet2index_jesup_RelSample_year=jet1index_jesup_RelSample_year;
                            jet1pt_jesup_RelSample_year=thisJet_jesup_RelSample_year.Pt(); jet1index_jesup_RelSample_year=k;
                        } else if (thisJet_jesup_RelSample_year.Pt()>jet2pt_jesup_RelSample_year) {
                            jet2pt_jesup_RelSample_year=thisJet_jesup_RelSample_year.Pt(); jet2index_jesup_RelSample_year=k;
                        }
                        if (abs(thisJet_jesup_RelSample_year.Eta())<2.5) {
                            njets_pt30_eta2p5_jesup_RelSample_year+=1;
                            if (thisJet_jesup_RelSample_year.Pt()>jet1pt2p5_jesup_RelSample_year) {
                                jet2pt2p5_jesup_RelSample_year=jet1pt2p5_jesup_RelSample_year; jet2index2p5_jesup_RelSample_year=jet1index2p5_jesup_RelSample_year;
                                jet1pt2p5_jesup_RelSample_year=thisJet_jesup_RelSample_year.Pt(); jet1index2p5_jesup_RelSample_year=k;
                            } else if (thisJet_jesup_RelSample_year.Pt()>jet2pt2p5_jesup_RelSample_year) {
                                jet2pt2p5_jesup_RelSample_year=thisJet_jesup_RelSample_year.Pt(); jet2index2p5_jesup_RelSample_year=k;
                            }
                        }
                }
                cout<<njets_pt30_eta4p7_jesup_RelSample_year<<" jets (jesup_RelSample_year)"<<endl;


		TLorentzVector Jet1_jesup_RelSample_year, Jet1_2p5_jesup_RelSample_year, Jet2_jesup_RelSample_year, Jet2_2p5_jesup_RelSample_year;
                if (njets_pt30_eta4p7_jesup_RelSample_year > 0) {
		    Jet1_jesup_RelSample_year.SetPtEtaPhiM(pt_jesup_split_RelSample_year[jet1index_jesup_RelSample_year],eta_jes_split_RelSample_year[jet1index_jesup_RelSample_year],phi_jes_split_RelSample_year[jet1index_jesup_RelSample_year], mass_jes_split_RelSample_year[jet1index_jesup_RelSample_year]);

                    pt_leadingjet_pt30_eta4p7_jesup_RelSample_year=Jet1_jesup_RelSample_year.Pt(); 
                    pTj1_jesup_RelSample_year=Jet1_jesup_RelSample_year.Pt(); 
                    etaj1_jesup_RelSample_year=Jet1_jesup_RelSample_year.Eta();
                    yj1_jesup_RelSample_year=Jet1_jesup_RelSample_year.Rapidity();
		    pT4lj_jesup_RelSample_year=(Higgs+Jet1_jesup_RelSample_year).Pt();
		    mass4lj_jesup_RelSample_year=(Higgs+Jet1_jesup_RelSample_year).M();
		    dPhiHj1_jesup_RelSample_year=deltaPhi(Higgs.Phi(),Jet1_jesup_RelSample_year.Phi());
                    dyHj1_jesup_RelSample_year=TMath::Abs(rapidity4l-yj1_jesup_RelSample_year);
                }
                if (njets_pt30_eta4p7_jesup_RelSample_year > 1) {
                    Jet2_jesup_RelSample_year.SetPtEtaPhiM(pt_jesup_split_RelSample_year[jet2index_jesup_RelSample_year],eta_jes_split_RelSample_year[jet2index_jesup_RelSample_year],phi_jes_split_RelSample_year[jet2index_jesup_RelSample_year], mass_jes_split_RelSample_year[jet2index_jesup_RelSample_year]);

                    pTj2_jesup_RelSample_year=Jet2_jesup_RelSample_year.Pt();
                    etaj2_jesup_RelSample_year=Jet2_jesup_RelSample_year.Eta();
                    yj2_jesup_RelSample_year=Jet2_jesup_RelSample_year.Rapidity();
                    pT4ljj_jesup_RelSample_year=(Higgs+Jet1_jesup_RelSample_year+Jet2_jesup_RelSample_year).Pt();
                    mass4ljj_jesup_RelSample_year=(Higgs+Jet1_jesup_RelSample_year+Jet2_jesup_RelSample_year).M();
		    mj1j2_jesup_RelSample_year=(Jet1_jesup_RelSample_year+Jet2_jesup_RelSample_year).M();
                    dEtaj1j2_jesup_RelSample_year=TMath::Abs(Jet1_jesup_RelSample_year.Eta()-Jet2_jesup_RelSample_year.Eta());
                    dPhij1j2_jesup_RelSample_year=deltaPhi(Jet1_jesup_RelSample_year.Phi(),Jet2_jesup_RelSample_year.Phi());
                    dPhiHj1j2_jesup_RelSample_year=deltaPhi(Higgs.Phi(),(Jet1_jesup_RelSample_year+Jet2_jesup_RelSample_year).Phi());
                }
                if (njets_pt30_eta2p5_jesup_RelSample_year > 0) {
                    Jet1_2p5_jesup_RelSample_year.SetPtEtaPhiM(pt_jesup_split_RelSample_year[jet1index2p5_jesup_RelSample_year],eta_jes_split_RelSample_year[jet1index2p5_jesup_RelSample_year],phi_jes_split_RelSample_year[jet1index2p5_jesup_RelSample_year], mass_jes_split_RelSample_year[jet1index2p5_jesup_RelSample_year]);
                    pt_leadingjet_pt30_eta2p5_jesup_RelSample_year=Jet1_2p5_jesup_RelSample_year.Pt();
                    pTj1_2p5_jesup_RelSample_year=Jet1_2p5_jesup_RelSample_year.Pt();
                    yj1_2p5_jesup_RelSample_year=Jet1_2p5_jesup_RelSample_year.Rapidity();
                    pT4lj_2p5_jesup_RelSample_year=(Higgs+Jet1_2p5_jesup_RelSample_year).Pt();
                    mass4lj_2p5_jesup_RelSample_year=(Higgs+Jet1_2p5_jesup_RelSample_year).M();
		    dPhiHj1_2p5_jesup_RelSample_year=deltaPhi(Higgs.Phi(),Jet1_2p5_jesup_RelSample_year.Phi());
                    dyHj1_2p5_jesup_RelSample_year=TMath::Abs(rapidity4l-yj1_2p5_jesup_RelSample_year);
                }
                if (njets_pt30_eta2p5_jesup_RelSample_year > 1) {
                    Jet2_2p5_jesup_RelSample_year.SetPtEtaPhiM(pt_jesup_split_RelSample_year[jet2index2p5_jesup_RelSample_year],eta_jes_split_RelSample_year[jet2index2p5_jesup_RelSample_year],phi_jes_split_RelSample_year[jet2index2p5_jesup_RelSample_year], mass_jes_split_RelSample_year[jet2index2p5_jesup_RelSample_year]);
                    pTj2_2p5_jesup_RelSample_year=Jet2_2p5_jesup_RelSample_year.Pt();
                    yj2_2p5_jesup_RelSample_year=Jet2_2p5_jesup_RelSample_year.Rapidity();
                    pT4ljj_2p5_jesup_RelSample_year=(Higgs+Jet1_2p5_jesup_RelSample_year+Jet2_2p5_jesup_RelSample_year).Pt();
                    mass4ljj_2p5_jesup_RelSample_year=(Higgs+Jet1_2p5_jesup_RelSample_year+Jet2_2p5_jesup_RelSample_year).M();
                    mj1j2_2p5_jesup_RelSample_year=(Jet1_2p5_jesup_RelSample_year+Jet2_2p5_jesup_RelSample_year).M();                    
                    dEtaj1j2_2p5_jesup_RelSample_year=TMath::Abs(Jet1_2p5_jesup_RelSample_year.Eta()-Jet2_2p5_jesup_RelSample_year.Eta());
                    dPhij1j2_2p5_jesup_RelSample_year=deltaPhi(Jet1_2p5_jesup_RelSample_year.Phi(),Jet2_2p5_jesup_RelSample_year.Phi());
                    dPhiHj1j2_2p5_jesup_RelSample_year=deltaPhi(Higgs.Phi(),(Jet1_2p5_jesup_RelSample_year+Jet2_2p5_jesup_RelSample_year).Phi());

                }

		pt_jesup_split_RelSample_year.clear();  

/////////////////
// RelSample_year dn start
                TauC_Inc_0j_EnergyWgt_jesdn_RelSample_year = TauC(pt_jesdn_split_RelSample_year, eta_jes_split_RelSample_year,phi_jes_split_RelSample_year,mass_jes_split_RelSample_year, Higgs);
                TauB_Inc_0j_pTWgt_jesdn_RelSample_year = TauB(pt_jesdn_split_RelSample_year, eta_jes_split_RelSample_year,phi_jes_split_RelSample_year,mass_jes_split_RelSample_year, Higgs);


                for( unsigned int k = 0; k< pt_jesdn_split_RelSample_year.size(); k++) {
                    if (pt_jesdn_split_RelSample_year[k]<30.0 || abs(eta_jes_split_RelSample_year[k])>4.7) continue;
                    TLorentzVector thisJet_jesdn_RelSample_year;
		    thisJet_jesdn_RelSample_year.SetPtEtaPhiM(pt_jesdn_split_RelSample_year[k],eta_jes_split_RelSample_year[k],phi_jes_split_RelSample_year[k],mass_jes_split_RelSample_year[k]);

                    bool isclean_H4l=true;
                        njets_pt30_eta4p7_jesdn_RelSample_year+=1;
			if (thisJet_jesdn_RelSample_year.Pt()>jet1pt_jesdn_RelSample_year) {
                            jet2pt_jesdn_RelSample_year=jet1pt_jesdn_RelSample_year; jet2index_jesdn_RelSample_year=jet1index_jesdn_RelSample_year;
                            jet1pt_jesdn_RelSample_year=thisJet_jesdn_RelSample_year.Pt(); jet1index_jesdn_RelSample_year=k;
                        } else if (thisJet_jesdn_RelSample_year.Pt()>jet2pt_jesdn_RelSample_year) {
                            jet2pt_jesdn_RelSample_year=thisJet_jesdn_RelSample_year.Pt(); jet2index_jesdn_RelSample_year=k;
                        }
                        if (abs(thisJet_jesdn_RelSample_year.Eta())<2.5) {
                            njets_pt30_eta2p5_jesdn_RelSample_year+=1;
                            if (thisJet_jesdn_RelSample_year.Pt()>jet1pt2p5_jesdn_RelSample_year) {
                                jet2pt2p5_jesdn_RelSample_year=jet1pt2p5_jesdn_RelSample_year; jet2index2p5_jesdn_RelSample_year=jet1index2p5_jesdn_RelSample_year;
                                jet1pt2p5_jesdn_RelSample_year=thisJet_jesdn_RelSample_year.Pt(); jet1index2p5_jesdn_RelSample_year=k;
                            } else if (thisJet_jesdn_RelSample_year.Pt()>jet2pt2p5_jesdn_RelSample_year) {
                                jet2pt2p5_jesdn_RelSample_year=thisJet_jesdn_RelSample_year.Pt(); jet2index2p5_jesdn_RelSample_year=k;
                            }
                        }
                }

// Filling RelSample_year dn variables

                TLorentzVector Jet1_jesdn_RelSample_year, Jet1_2p5_jesdn_RelSample_year, Jet2_jesdn_RelSample_year, Jet2_2p5_jesdn_RelSample_year;

                if (njets_pt30_eta4p7_jesdn_RelSample_year > 0) { 
                    Jet1_jesdn_RelSample_year.SetPtEtaPhiM(pt_jesdn_split_RelSample_year[jet1index_jesdn_RelSample_year],eta_jes_split_RelSample_year[jet1index_jesdn_RelSample_year],phi_jes_split_RelSample_year[jet1index_jesdn_RelSample_year], mass_jes_split_RelSample_year[jet1index_jesdn_RelSample_year]);
                    pt_leadingjet_pt30_eta4p7_jesdn_RelSample_year=Jet1_jesdn_RelSample_year.Pt(); 
                    pTj1_jesdn_RelSample_year=Jet1_jesdn_RelSample_year.Pt(); 
                    etaj1_jesdn_RelSample_year=Jet1_jesdn_RelSample_year.Eta();
                    yj1_jesdn_RelSample_year=Jet1_jesdn_RelSample_year.Rapidity();
                    pT4lj_jesdn_RelSample_year=(Higgs+Jet1_jesdn_RelSample_year).Pt();
                    mass4lj_jesdn_RelSample_year=(Higgs+Jet1_jesdn_RelSample_year).M();
                    dPhiHj1_jesdn_RelSample_year=deltaPhi(Higgs.Phi(),Jet1_jesdn_RelSample_year.Phi());
                    dyHj1_jesdn_RelSample_year=TMath::Abs(rapidity4l-yj1_jesdn_RelSample_year);
                }    
                if (njets_pt30_eta4p7_jesdn_RelSample_year > 1) { 
                    Jet2_jesdn_RelSample_year.SetPtEtaPhiM(pt_jesdn_split_RelSample_year[jet2index_jesdn_RelSample_year],eta_jes_split_RelSample_year[jet2index_jesdn_RelSample_year],phi_jes_split_RelSample_year[jet2index_jesdn_RelSample_year], mass_jes_split_RelSample_year[jet2index_jesdn_RelSample_year]);

                    pTj2_jesdn_RelSample_year=Jet2_jesdn_RelSample_year.Pt();
                    etaj2_jesdn_RelSample_year=Jet2_jesdn_RelSample_year.Eta();
                    yj2_jesdn_RelSample_year=Jet2_jesdn_RelSample_year.Rapidity();
                    pT4ljj_jesdn_RelSample_year=(Higgs+Jet1_jesdn_RelSample_year+Jet2_jesdn_RelSample_year).Pt();
                    mass4ljj_jesdn_RelSample_year=(Higgs+Jet1_jesdn_RelSample_year+Jet2_jesdn_RelSample_year).M();
                    dEtaj1j2_jesdn_RelSample_year=TMath::Abs(Jet1_jesdn_RelSample_year.Eta()-Jet2_jesdn_RelSample_year.Eta());
                    dPhij1j2_jesdn_RelSample_year=deltaPhi(Jet1_jesdn_RelSample_year.Phi(),Jet2_jesdn_RelSample_year.Phi());
                    dPhiHj1j2_jesdn_RelSample_year=deltaPhi(Higgs.Phi(),(Jet1_jesdn_RelSample_year+Jet2_jesdn_RelSample_year).Phi());
                }    
                if (njets_pt30_eta2p5_jesdn_RelSample_year > 0) { 
                    Jet1_2p5_jesdn_RelSample_year.SetPtEtaPhiM(pt_jesdn_split_RelSample_year[jet1index2p5_jesdn_RelSample_year],eta_jes_split_RelSample_year[jet1index2p5_jesdn_RelSample_year],phi_jes_split_RelSample_year[jet1index2p5_jesdn_RelSample_year], mass_jes_split_RelSample_year[jet1index2p5_jesdn_RelSample_year]);
                    pt_leadingjet_pt30_eta2p5_jesdn_RelSample_year=Jet1_2p5_jesdn_RelSample_year.Pt();
                    pTj1_2p5_jesdn_RelSample_year=Jet1_2p5_jesdn_RelSample_year.Pt();
                    yj1_2p5_jesdn_RelSample_year=Jet1_2p5_jesdn_RelSample_year.Rapidity();
                    pT4lj_2p5_jesdn_RelSample_year=(Higgs+Jet1_2p5_jesdn_RelSample_year).Pt();
                    mass4lj_2p5_jesdn_RelSample_year=(Higgs+Jet1_2p5_jesdn_RelSample_year).M();
		    dPhiHj1_2p5_jesdn_RelSample_year=deltaPhi(Higgs.Phi(),Jet1_2p5_jesdn_RelSample_year.Phi());
                    dyHj1_2p5_jesdn_RelSample_year=TMath::Abs(rapidity4l-yj1_2p5_jesdn_RelSample_year);
                }    
                if (njets_pt30_eta2p5_jesdn_RelSample_year > 1) { 
                    Jet2_2p5_jesdn_RelSample_year.SetPtEtaPhiM(pt_jesdn_split_RelSample_year[jet2index2p5_jesdn_RelSample_year],eta_jes_split_RelSample_year[jet2index2p5_jesdn_RelSample_year],phi_jes_split_RelSample_year[jet2index2p5_jesdn_RelSample_year], mass_jes_split_RelSample_year[jet2index2p5_jesdn_RelSample_year]);
                    pTj2_2p5_jesdn_RelSample_year=Jet2_2p5_jesdn_RelSample_year.Pt();
                    yj2_2p5_jesdn_RelSample_year=Jet2_2p5_jesdn_RelSample_year.Rapidity();
                    pT4ljj_2p5_jesdn_RelSample_year=(Higgs+Jet1_2p5_jesdn_RelSample_year+Jet2_2p5_jesdn_RelSample_year).Pt();
                    mass4ljj_2p5_jesdn_RelSample_year=(Higgs+Jet1_2p5_jesdn_RelSample_year+Jet2_2p5_jesdn_RelSample_year).M();
                    mj1j2_2p5_jesdn_RelSample_year=(Jet1_2p5_jesdn_RelSample_year+Jet2_2p5_jesdn_RelSample_year).M();     
                    dEtaj1j2_2p5_jesdn_RelSample_year=TMath::Abs(Jet1_2p5_jesdn_RelSample_year.Eta()-Jet2_2p5_jesdn_RelSample_year.Eta());
	            dPhij1j2_2p5_jesdn_RelSample_year=deltaPhi(Jet1_2p5_jesdn_RelSample_year.Phi(),Jet2_2p5_jesdn_RelSample_year.Phi());
                    dPhiHj1j2_2p5_jesdn_RelSample_year=deltaPhi(Higgs.Phi(),(Jet1_2p5_jesdn_RelSample_year+Jet2_2p5_jesdn_RelSample_year).Phi());

                }    

                pt_jesdn_split_RelSample_year.clear();

		eta_jes_split_RelSample_year.clear();
		phi_jes_split_RelSample_year.clear();
		mass_jes_split_RelSample_year.clear();


////////////////////////////////////////////////////////////////////////////////
////    Starting "Total" JES uncertainty source
////////////////////////////////////////////////////////////////////////////////

                TauC_Inc_0j_EnergyWgt_jesup_Total = TauC(pt_jesup_split_Total, eta_jes_split_Total,phi_jes_split_Total,mass_jes_split_Total, Higgs);
                TauB_Inc_0j_pTWgt_jesup_Total = TauB(pt_jesup_split_Total, eta_jes_split_Total,phi_jes_split_Total,mass_jes_split_Total, Higgs);


                for( unsigned int k = 0; k< pt_jesup_split_Total.size(); k++) {
                    if (pt_jesup_split_Total[k]<30.0 || abs(eta_jes_split_Total[k])>4.7) continue;
		    TLorentzVector thisJet_jesup_Total;
		    thisJet_jesup_Total.SetPtEtaPhiM(pt_jesup_split_Total[k],eta_jes_split_Total[k],phi_jes_split_Total[k],mass_jes_split_Total[k]);

                        njets_pt30_eta4p7_jesup_Total+=1;  

                        if (thisJet_jesup_Total.Pt()>jet1pt_jesup_Total) {
                            jet2pt_jesup_Total=jet1pt_jesup_Total; jet2index_jesup_Total=jet1index_jesup_Total;
                            jet1pt_jesup_Total=thisJet_jesup_Total.Pt(); jet1index_jesup_Total=k;
                        } else if (thisJet_jesup_Total.Pt()>jet2pt_jesup_Total) {
                            jet2pt_jesup_Total=thisJet_jesup_Total.Pt(); jet2index_jesup_Total=k;
                        }
                        if (abs(thisJet_jesup_Total.Eta())<2.5) {
                            njets_pt30_eta2p5_jesup_Total+=1;
                            if (thisJet_jesup_Total.Pt()>jet1pt2p5_jesup_Total) {
                                jet2pt2p5_jesup_Total=jet1pt2p5_jesup_Total; jet2index2p5_jesup_Total=jet1index2p5_jesup_Total;
                                jet1pt2p5_jesup_Total=thisJet_jesup_Total.Pt(); jet1index2p5_jesup_Total=k;
                            } else if (thisJet_jesup_Total.Pt()>jet2pt2p5_jesup_Total) {
                                jet2pt2p5_jesup_Total=thisJet_jesup_Total.Pt(); jet2index2p5_jesup_Total=k;
                            }
                        }
                }
                cout<<njets_pt30_eta4p7_jesup_Total<<" jets (jesup_Total)"<<endl;


		TLorentzVector Jet1_jesup_Total, Jet1_2p5_jesup_Total, Jet2_jesup_Total, Jet2_2p5_jesup_Total;
                if (njets_pt30_eta4p7_jesup_Total > 0) {
		    Jet1_jesup_Total.SetPtEtaPhiM(pt_jesup_split_Total[jet1index_jesup_Total],eta_jes_split_Total[jet1index_jesup_Total],phi_jes_split_Total[jet1index_jesup_Total], mass_jes_split_Total[jet1index_jesup_Total]);

                    pt_leadingjet_pt30_eta4p7_jesup_Total=Jet1_jesup_Total.Pt(); 
                    pTj1_jesup_Total=Jet1_jesup_Total.Pt(); 
                    etaj1_jesup_Total=Jet1_jesup_Total.Eta();
                    yj1_jesup_Total=Jet1_jesup_Total.Rapidity();
		    pT4lj_jesup_Total=(Higgs+Jet1_jesup_Total).Pt();
		    mass4lj_jesup_Total=(Higgs+Jet1_jesup_Total).M();
		    dPhiHj1_jesup_Total=deltaPhi(Higgs.Phi(),Jet1_jesup_Total.Phi());
                    dyHj1_jesup_Total=TMath::Abs(rapidity4l-yj1_jesup_Total);
                }
                if (njets_pt30_eta4p7_jesup_Total > 1) {
                    Jet2_jesup_Total.SetPtEtaPhiM(pt_jesup_split_Total[jet2index_jesup_Total],eta_jes_split_Total[jet2index_jesup_Total],phi_jes_split_Total[jet2index_jesup_Total], mass_jes_split_Total[jet2index_jesup_Total]);

                    pTj2_jesup_Total=Jet2_jesup_Total.Pt();
                    etaj2_jesup_Total=Jet2_jesup_Total.Eta();
                    yj2_jesup_Total=Jet2_jesup_Total.Rapidity();
                    pT4ljj_jesup_Total=(Higgs+Jet1_jesup_Total+Jet2_jesup_Total).Pt();
                    mass4ljj_jesup_Total=(Higgs+Jet1_jesup_Total+Jet2_jesup_Total).M();
		    mj1j2_jesup_Total=(Jet1_jesup_Total+Jet2_jesup_Total).M();
                    dEtaj1j2_jesup_Total=TMath::Abs(Jet1_jesup_Total.Eta()-Jet2_jesup_Total.Eta());
                    dPhij1j2_jesup_Total=deltaPhi(Jet1_jesup_Total.Phi(),Jet2_jesup_Total.Phi());
                    dPhiHj1j2_jesup_Total=deltaPhi(Higgs.Phi(),(Jet1_jesup_Total+Jet2_jesup_Total).Phi());
                }
                if (njets_pt30_eta2p5_jesup_Total > 0) {
                    Jet1_2p5_jesup_Total.SetPtEtaPhiM(pt_jesup_split_Total[jet1index2p5_jesup_Total],eta_jes_split_Total[jet1index2p5_jesup_Total],phi_jes_split_Total[jet1index2p5_jesup_Total], mass_jes_split_Total[jet1index2p5_jesup_Total]);
                    pt_leadingjet_pt30_eta2p5_jesup_Total=Jet1_2p5_jesup_Total.Pt();
                    pTj1_2p5_jesup_Total=Jet1_2p5_jesup_Total.Pt();
                    yj1_2p5_jesup_Total=Jet1_2p5_jesup_Total.Rapidity();
                    pT4lj_2p5_jesup_Total=(Higgs+Jet1_2p5_jesup_Total).Pt();
                    mass4lj_2p5_jesup_Total=(Higgs+Jet1_2p5_jesup_Total).M();
		    dPhiHj1_2p5_jesup_Total=deltaPhi(Higgs.Phi(),Jet1_2p5_jesup_Total.Phi());
                    dyHj1_2p5_jesup_Total=TMath::Abs(rapidity4l-yj1_2p5_jesup_Total);
                }
                if (njets_pt30_eta2p5_jesup_Total > 1) {
                    Jet2_2p5_jesup_Total.SetPtEtaPhiM(pt_jesup_split_Total[jet2index2p5_jesup_Total],eta_jes_split_Total[jet2index2p5_jesup_Total],phi_jes_split_Total[jet2index2p5_jesup_Total], mass_jes_split_Total[jet2index2p5_jesup_Total]);
                    pTj2_2p5_jesup_Total=Jet2_2p5_jesup_Total.Pt();
                    yj2_2p5_jesup_Total=Jet2_2p5_jesup_Total.Rapidity();
                    pT4ljj_2p5_jesup_Total=(Higgs+Jet1_2p5_jesup_Total+Jet2_2p5_jesup_Total).Pt();
                    mass4ljj_2p5_jesup_Total=(Higgs+Jet1_2p5_jesup_Total+Jet2_2p5_jesup_Total).M();
                    mj1j2_2p5_jesup_Total=(Jet1_2p5_jesup_Total+Jet2_2p5_jesup_Total).M();                    
                    dEtaj1j2_2p5_jesup_Total=TMath::Abs(Jet1_2p5_jesup_Total.Eta()-Jet2_2p5_jesup_Total.Eta());
                    dPhij1j2_2p5_jesup_Total=deltaPhi(Jet1_2p5_jesup_Total.Phi(),Jet2_2p5_jesup_Total.Phi());
                    dPhiHj1j2_2p5_jesup_Total=deltaPhi(Higgs.Phi(),(Jet1_2p5_jesup_Total+Jet2_2p5_jesup_Total).Phi());

                }

		pt_jesup_split_Total.clear();  

/////////////////
// Total dn start
                TauC_Inc_0j_EnergyWgt_jesdn_Total = TauC(pt_jesdn_split_Total, eta_jes_split_Total,phi_jes_split_Total,mass_jes_split_Total, Higgs);
                TauB_Inc_0j_pTWgt_jesdn_Total = TauB(pt_jesdn_split_Total, eta_jes_split_Total,phi_jes_split_Total,mass_jes_split_Total, Higgs);


                for( unsigned int k = 0; k< pt_jesdn_split_Total.size(); k++) {
                    if (pt_jesdn_split_Total[k]<30.0 || abs(eta_jes_split_Total[k])>4.7) continue;
                    TLorentzVector thisJet_jesdn_Total;
		    thisJet_jesdn_Total.SetPtEtaPhiM(pt_jesdn_split_Total[k],eta_jes_split_Total[k],phi_jes_split_Total[k],mass_jes_split_Total[k]);

                    bool isclean_H4l=true;
                        njets_pt30_eta4p7_jesdn_Total+=1;
			if (thisJet_jesdn_Total.Pt()>jet1pt_jesdn_Total) {
                            jet2pt_jesdn_Total=jet1pt_jesdn_Total; jet2index_jesdn_Total=jet1index_jesdn_Total;
                            jet1pt_jesdn_Total=thisJet_jesdn_Total.Pt(); jet1index_jesdn_Total=k;
                        } else if (thisJet_jesdn_Total.Pt()>jet2pt_jesdn_Total) {
                            jet2pt_jesdn_Total=thisJet_jesdn_Total.Pt(); jet2index_jesdn_Total=k;
                        }
                        if (abs(thisJet_jesdn_Total.Eta())<2.5) {
                            njets_pt30_eta2p5_jesdn_Total+=1;
                            if (thisJet_jesdn_Total.Pt()>jet1pt2p5_jesdn_Total) {
                                jet2pt2p5_jesdn_Total=jet1pt2p5_jesdn_Total; jet2index2p5_jesdn_Total=jet1index2p5_jesdn_Total;
                                jet1pt2p5_jesdn_Total=thisJet_jesdn_Total.Pt(); jet1index2p5_jesdn_Total=k;
                            } else if (thisJet_jesdn_Total.Pt()>jet2pt2p5_jesdn_Total) {
                                jet2pt2p5_jesdn_Total=thisJet_jesdn_Total.Pt(); jet2index2p5_jesdn_Total=k;
                            }
                        }
                }

// Filling Total dn variables

                TLorentzVector Jet1_jesdn_Total, Jet1_2p5_jesdn_Total, Jet2_jesdn_Total, Jet2_2p5_jesdn_Total;

                if (njets_pt30_eta4p7_jesdn_Total > 0) { 
                    Jet1_jesdn_Total.SetPtEtaPhiM(pt_jesdn_split_Total[jet1index_jesdn_Total],eta_jes_split_Total[jet1index_jesdn_Total],phi_jes_split_Total[jet1index_jesdn_Total], mass_jes_split_Total[jet1index_jesdn_Total]);
                    pt_leadingjet_pt30_eta4p7_jesdn_Total=Jet1_jesdn_Total.Pt(); 
                    pTj1_jesdn_Total=Jet1_jesdn_Total.Pt(); 
                    etaj1_jesdn_Total=Jet1_jesdn_Total.Eta();
                    yj1_jesdn_Total=Jet1_jesdn_Total.Rapidity();
                    pT4lj_jesdn_Total=(Higgs+Jet1_jesdn_Total).Pt();
                    mass4lj_jesdn_Total=(Higgs+Jet1_jesdn_Total).M();
                    dPhiHj1_jesdn_Total=deltaPhi(Higgs.Phi(),Jet1_jesdn_Total.Phi());
                    dyHj1_jesdn_Total=TMath::Abs(rapidity4l-yj1_jesdn_Total);
                }    
                if (njets_pt30_eta4p7_jesdn_Total > 1) { 
                    Jet2_jesdn_Total.SetPtEtaPhiM(pt_jesdn_split_Total[jet2index_jesdn_Total],eta_jes_split_Total[jet2index_jesdn_Total],phi_jes_split_Total[jet2index_jesdn_Total], mass_jes_split_Total[jet2index_jesdn_Total]);

                    pTj2_jesdn_Total=Jet2_jesdn_Total.Pt();
                    etaj2_jesdn_Total=Jet2_jesdn_Total.Eta();
                    yj2_jesdn_Total=Jet2_jesdn_Total.Rapidity();
                    pT4ljj_jesdn_Total=(Higgs+Jet1_jesdn_Total+Jet2_jesdn_Total).Pt();
                    mass4ljj_jesdn_Total=(Higgs+Jet1_jesdn_Total+Jet2_jesdn_Total).M();
                    dEtaj1j2_jesdn_Total=TMath::Abs(Jet1_jesdn_Total.Eta()-Jet2_jesdn_Total.Eta());
                    dPhij1j2_jesdn_Total=deltaPhi(Jet1_jesdn_Total.Phi(),Jet2_jesdn_Total.Phi());
                    dPhiHj1j2_jesdn_Total=deltaPhi(Higgs.Phi(),(Jet1_jesdn_Total+Jet2_jesdn_Total).Phi());
                }    
                if (njets_pt30_eta2p5_jesdn_Total > 0) { 
                    Jet1_2p5_jesdn_Total.SetPtEtaPhiM(pt_jesdn_split_Total[jet1index2p5_jesdn_Total],eta_jes_split_Total[jet1index2p5_jesdn_Total],phi_jes_split_Total[jet1index2p5_jesdn_Total], mass_jes_split_Total[jet1index2p5_jesdn_Total]);
                    pt_leadingjet_pt30_eta2p5_jesdn_Total=Jet1_2p5_jesdn_Total.Pt();
                    pTj1_2p5_jesdn_Total=Jet1_2p5_jesdn_Total.Pt();
                    yj1_2p5_jesdn_Total=Jet1_2p5_jesdn_Total.Rapidity();
                    pT4lj_2p5_jesdn_Total=(Higgs+Jet1_2p5_jesdn_Total).Pt();
                    mass4lj_2p5_jesdn_Total=(Higgs+Jet1_2p5_jesdn_Total).M();
		    dPhiHj1_2p5_jesdn_Total=deltaPhi(Higgs.Phi(),Jet1_2p5_jesdn_Total.Phi());
                    dyHj1_2p5_jesdn_Total=TMath::Abs(rapidity4l-yj1_2p5_jesdn_Total);
                }    
                if (njets_pt30_eta2p5_jesdn_Total > 1) { 
                    Jet2_2p5_jesdn_Total.SetPtEtaPhiM(pt_jesdn_split_Total[jet2index2p5_jesdn_Total],eta_jes_split_Total[jet2index2p5_jesdn_Total],phi_jes_split_Total[jet2index2p5_jesdn_Total], mass_jes_split_Total[jet2index2p5_jesdn_Total]);
                    pTj2_2p5_jesdn_Total=Jet2_2p5_jesdn_Total.Pt();
                    yj2_2p5_jesdn_Total=Jet2_2p5_jesdn_Total.Rapidity();
                    pT4ljj_2p5_jesdn_Total=(Higgs+Jet1_2p5_jesdn_Total+Jet2_2p5_jesdn_Total).Pt();
                    mass4ljj_2p5_jesdn_Total=(Higgs+Jet1_2p5_jesdn_Total+Jet2_2p5_jesdn_Total).M();
                    mj1j2_2p5_jesdn_Total=(Jet1_2p5_jesdn_Total+Jet2_2p5_jesdn_Total).M();     
                    dEtaj1j2_2p5_jesdn_Total=TMath::Abs(Jet1_2p5_jesdn_Total.Eta()-Jet2_2p5_jesdn_Total.Eta());
	            dPhij1j2_2p5_jesdn_Total=deltaPhi(Jet1_2p5_jesdn_Total.Phi(),Jet2_2p5_jesdn_Total.Phi());
                    dPhiHj1j2_2p5_jesdn_Total=deltaPhi(Higgs.Phi(),(Jet1_2p5_jesdn_Total+Jet2_2p5_jesdn_Total).Phi());

                }    

                pt_jesdn_split_Total.clear();

		eta_jes_split_Total.clear();
		phi_jes_split_Total.clear();
		mass_jes_split_Total.clear();








//////////////////////////////////////////////////////////////////////////////
//    END filling all JES uncertainty sources
//////////////////////////////////////////////////////////////////////////////





/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        if (!isMC || (isMC && (isSignal || (!isSignal && (passedZ4lSelection ||passedZXCRSelection)))) )
        {
            if (pT4l>0)
            {
                //TLorentzVector j1, j2, Higgs, hj, jj, hjj, j1_jesup, j1_jesdn, j2_jesup, j2_jesdn, hj_jesup, hj_jesdn, jj_jesup, jj_jesdn, hjj_jesup, hjj_jesdn;
                TLorentzVector j1, j2, hj, jj, hjj, j1_jesup, j1_jesdn, j2_jesup, j2_jesdn, hj_jesup, hj_jesdn, jj_jesup, jj_jesdn, hjj_jesup, hjj_jesdn;
                //Higgs.SetPtEtaPhiM(pT4l, eta4l, phi4l, mass4l);
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
    cout<<"skimming done !!  "<<endl;
    exit(EXIT_FAILURE);
}

