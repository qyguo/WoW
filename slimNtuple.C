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

void order_pt(vector<float> *GENlep_pt, vector<int> &order_index)
{
    vector<float> GENlep_pt_ordered;
    GENlep_pt_ordered.clear();
    order_index.clear();
    for(unsigned int i=0; i<GENlep_pt->size(); i++)
    {
        if(GENlep_pt_ordered.size()==0 || (*GENlep_pt)[i]<GENlep_pt_ordered[GENlep_pt_ordered.size()-1])
        {
            GENlep_pt_ordered.push_back((*GENlep_pt)[i]);
            order_index.push_back(i);
            continue;
        }
        for(unsigned int j=0; j<GENlep_pt_ordered.size(); j++)
        {
            if((*GENlep_pt)[i]>=GENlep_pt_ordered[j])
            {
                GENlep_pt_ordered.insert(GENlep_pt_ordered.begin()+j,(*GENlep_pt)[i]);
                order_index.insert(order_index.begin()+j, i);
                break;
            }
        }
    }
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
            TauB_j = sqrt(theJet.Pt()*theJet.Pt() + theJet.M()*theJet.M())*exp(-1*(abs(theJet.Rapidity() - H.Rapidity())));
            if (TauB_j > TauB_jmax) {
                TauB_jmax = TauB_j;
                }
        }
    return TauB_jmax;
}



//void slimNtuple(const int & _year_=2017, const string & _name_DS_="GluGluHToZZTo4L_M125_2017", const bool & isMC = true,  const bool & isSignal = true, const bool & _OldNtuple=false, const bool & _addNewVar=false, const bool & _Test = false) {
void slimNtuple(const int & _year_=2016, const string & _name_DS_="testGGH_nnlops_GENonly", const bool & isMC = true,  const bool & isSignal = true, const bool & _OldNtuple=true, const bool & _addNewVar=true, const bool & _Test = false) {

    int year = _year_;
    //string pre_name = "root://cmsio5.rc.ufl.edu:1094//store/user/t2/users/ferrico/Qianying/";
    string pre_name = "/publicfs/cms/data/hzz/jtahir/UL2016/CMSSW_10_6_26/src/WoW/";
    string whichYear[3] = {"2016/", "2017/", "2018/"};
    string MC__ = "MC/";
    if (!isMC) MC__ = "DATA/";
    //string name_DS = "GluGluHToZZTo4L_M125_2017";
    string name_DS = _name_DS_.c_str();

    //TString filename = prefix+".root";
    //string filename = (pre_name + whichYear[year-2016] + MC__ + name_DS + ".root").c_str();
    string filename = (pre_name + name_DS + ".root").c_str();
    //filename = "Sync_106X_2018UL_30112021_-1ggH_V2.root";

    std::cout<<"Year: "<<year<<std::endl;
    std::cout<<filename<<std::endl;

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

    vector<int> order_index;
    vector<float> pt_nom {};
    vector<float> eta_nom {};
    vector<float> phi_nom {};
    vector<float> mass_nom {};

    TLorentzVector j1, j2, Higgs, hj, jj, hjj;
    TLorentzVector GENj1, GENj2, GENl1, GENl2, GENl3, GENl4, GENHiggs, GENhj, GENjj, GENhjj;
    TLorentzVector GENj1_2p5, GENj2_2p5, GENhj_2p5, GENjj_2p5, GENhjj_2p5;
    float pTj2(-1.0), mj1j2(-1.0), dEtaj1j2(-99.0), dPhij1j2(-99.0), pT4lj(-1.0), mass4lj(-1.0), pT4ljj(-1.0), mass4ljj(-1.0);
    float GENpTj2(-1.0), GENmj1j2(-1.0), GENdEtaj1j2(-99.0), GENdPhij1j2(-99.0), GENpT4lj(-1.0), GENmass4lj(-1.0), GENpT4ljj(-1.0), GENmass4ljj(-1.0), GEN_TauB_Inc_0j_pTWgt(-99), GEN_TauC_Inc_0j_EnergyWgt(-99);
    float GENpTj2_2p5(-1.0), GENmj1j2_2p5(-1.0), GENdEtaj1j2_2p5(-99.0), GENdPhij1j2_2p5(-99.0), GENpT4lj_2p5(-1.0), GENmass4lj_2p5(-1.0), GENpT4ljj_2p5(-1.0), GENmass4ljj_2p5(-1.0);


    //ULong64_t Run, LumiSect, Event;
    //Bool_t          passedFullSelection;
    //Float_t         dataMCWeight;
    //Float_t         eventWeight;
    //vector<int>     *lep_id;
    //vector<float>   *lep_pt;
    //vector<float>   *lep_eta;
    //Int_t           lep_Hindex[4];
    //vector<float>   *lep_dataMC;
    //vector<float>   lep_dataMC_new;
    //vector<float>   *lep_dataMCErr;
    //lep_id = 0; 
    //lep_pt = 0;    
    //lep_eta = 0;  
    //lep_dataMC = 0;
    //lep_dataMCErr = 0;
    //oldtree->SetBranchAddress("Run",&Run);
    //oldtree->SetBranchAddress("LumiSect",&LumiSect);
    //oldtree->SetBranchAddress("Event",&Event);
    //oldtree->SetBranchAddress("passedFullSelection", &passedFullSelection);
    //oldtree->SetBranchAddress("dataMCWeight", &dataMCWeight);
    //oldtree->SetBranchAddress("eventWeight", &eventWeight);
    //oldtree->SetBranchAddress("lep_id", &lep_id);
    //oldtree->SetBranchAddress("lep_pt", &lep_pt);
    //oldtree->SetBranchAddress("lep_Hindex", lep_Hindex);
    //oldtree->SetBranchAddress("lep_eta", &lep_eta);
    //oldtree->SetBranchAddress("lep_dataMC", &lep_dataMC);
    //oldtree->SetBranchAddress("lep_dataMCErr", &lep_dataMCErr);
    
    Bool_t          passedFullSelection;
    Bool_t          passedFiducialSelection;
    Int_t           njets_pt30_eta4p7;
    vector<float>   *jet_pt;
    vector<float>   *jet_eta;
    vector<float>   *jet_phi;
    vector<float>   *jet_mass;
    Int_t           jet1index;
    Int_t           jet2index;
    Float_t         pT4l;
    Float_t         eta4l;
    Float_t         phi4l;
    Float_t         mass4l;

    jet_pt = 0;
    jet_eta = 0;
    jet_phi = 0;
    jet_mass = 0;

    Int_t           GENnjets_pt30_eta4p7;
    Int_t           GENnjets_pt30_eta2p5;
    vector<float>   *GENjet_pt;
    vector<float>   *GENjet_eta;
    vector<float>   *GENjet_phi;
    vector<float>   *GENjet_mass;
    //Int_t           GENjet1index;
    //Int_t           GENjet2index;
    vector<float>   *GENlep_pt;
    vector<float>   *GENlep_eta;
    vector<float>   *GENlep_phi;
    vector<float>   *GENlep_mass;
    Int_t           GENlep_Hindex[4];
    Float_t         GENmass4l;
    GENjet_pt = 0;
    GENjet_eta = 0;
    GENjet_phi = 0;
    GENjet_mass = 0;
    GENlep_pt = 0;
    GENlep_eta = 0;
    GENlep_phi = 0;
    GENlep_mass = 0;

    if (_addNewVar && _OldNtuple)
    {
        oldtree->SetBranchAddress("passedFullSelection", &passedFullSelection);
        oldtree->SetBranchAddress("passedFiducialSelection", &passedFiducialSelection);
        oldtree->SetBranchAddress("njets_pt30_eta4p7", &njets_pt30_eta4p7);
        oldtree->SetBranchAddress("jet_pt", &jet_pt);
        oldtree->SetBranchAddress("jet_eta", &jet_eta);
        oldtree->SetBranchAddress("jet_phi", &jet_phi);
        oldtree->SetBranchAddress("jet_mass", &jet_mass);
        oldtree->SetBranchAddress("jet1index", &jet1index);
        oldtree->SetBranchAddress("jet2index", &jet2index);
        oldtree->SetBranchAddress("pT4l", &pT4l);
        oldtree->SetBranchAddress("eta4l", &eta4l);
        oldtree->SetBranchAddress("phi4l", &phi4l);
        oldtree->SetBranchAddress("mass4l", &mass4l);

        oldtree->SetBranchAddress("GENnjets_pt30_eta4p7", &GENnjets_pt30_eta4p7);
        oldtree->SetBranchAddress("GENnjets_pt30_eta2p5", &GENnjets_pt30_eta2p5);
        oldtree->SetBranchAddress("GENjet_pt", &GENjet_pt);
        oldtree->SetBranchAddress("GENjet_eta", &GENjet_eta);
        oldtree->SetBranchAddress("GENjet_phi", &GENjet_phi);
        oldtree->SetBranchAddress("GENjet_mass", &GENjet_mass);
        //oldtree->SetBranchAddress("GENjet1index", &GENjet1index);
        //oldtree->SetBranchAddress("GENjet2index", &GENjet2index);
        oldtree->SetBranchAddress("GENlep_pt", &GENlep_pt);
        oldtree->SetBranchAddress("GENlep_eta", &GENlep_eta);
        oldtree->SetBranchAddress("GENlep_phi", &GENlep_phi);
        oldtree->SetBranchAddress("GENlep_mass", &GENlep_mass);
        oldtree->SetBranchAddress("GENlep_Hindex", &GENlep_Hindex);
        oldtree->SetBranchAddress("GENmass4l", &GENmass4l);
    }

    Bool_t passedZ4lSelection;
    Bool_t passedZXCRSelection;
    oldtree->SetBranchAddress("passedZ4lSelection", &passedZ4lSelection);
    oldtree->SetBranchAddress("passedZXCRSelection", &passedZXCRSelection);

    oldtree->SetBranchStatus("*",0); //Disables All Branches
    //Then enables only select branches

    if (_addNewVar && _OldNtuple)
    {
        oldtree->SetBranchStatus("jet_pt", 1);
        oldtree->SetBranchStatus("jet_eta", 1);
        oldtree->SetBranchStatus("jet_phi", 1);
        oldtree->SetBranchStatus("jet_mass", 1);
        oldtree->SetBranchStatus("jet1index", 1);
        oldtree->SetBranchStatus("jet2index", 1);

        oldtree->SetBranchStatus("GENjet_pt", 1);
        oldtree->SetBranchStatus("GENjet_eta", 1);
        oldtree->SetBranchStatus("GENjet_phi", 1);
        oldtree->SetBranchStatus("GENjet_mass", 1);
        //oldtree->SetBranchStatus("GENjet1index", 1);
        //oldtree->SetBranchStatus("GENjet2index", 1);
        oldtree->SetBranchStatus("GENlep_pt", 1);
        oldtree->SetBranchStatus("GENlep_eta", 1);
        oldtree->SetBranchStatus("GENlep_phi", 1);
        oldtree->SetBranchStatus("GENlep_mass", 1);
        oldtree->SetBranchStatus("GENmass4l", 1);
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
//    oldtree->SetBranchStatus("prefiringWeight",1);
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

    if (!_OldNtuple) {
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


    if (!_OldNtuple) {
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
            (name_DS+"_slimmed.root").c_str()
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
    //newtree->Branch("dataMCWeight_new", &dataMCWeight_new);
    //newtree->Branch("eventWeight_new", &eventWeight_new);
    if (_addNewVar && _OldNtuple)
    {
        newtree->Branch("pTj2", &pTj2);
        newtree->Branch("mj1j2", &mj1j2);
        newtree->Branch("dEtaj1j2", &dEtaj1j2);
        newtree->Branch("dPhij1j2", &dPhij1j2);
        newtree->Branch("pT4lj", &pT4lj);
        newtree->Branch("mass4lj", &mass4lj);
        newtree->Branch("pT4ljj", &pT4ljj);
        newtree->Branch("mass4ljj", &mass4ljj);

        newtree->Branch("GENpTj2", &GENpTj2);
        newtree->Branch("GENmj1j2", &GENmj1j2);
        newtree->Branch("GENdEtaj1j2", &GENdEtaj1j2);
        newtree->Branch("GENdPhij1j2", &GENdPhij1j2);
        newtree->Branch("GENpT4lj", &GENpT4lj);
        newtree->Branch("GENmass4lj", &GENmass4lj);
        newtree->Branch("GENpT4ljj", &GENpT4ljj);
        newtree->Branch("GENmass4ljj", &GENmass4ljj);
        newtree->Branch("GENpTj2_2p5", &GENpTj2_2p5);
        newtree->Branch("GENmj1j2_2p5", &GENmj1j2_2p5);
        newtree->Branch("GENdEtaj1j2_2p5", &GENdEtaj1j2_2p5);
        newtree->Branch("GENdPhij1j2_2p5", &GENdPhij1j2_2p5);
        newtree->Branch("GENpT4lj_2p5", &GENpT4lj_2p5);
        newtree->Branch("GENmass4lj_2p5", &GENmass4lj_2p5);
        newtree->Branch("GENpT4ljj_2p5", &GENpT4ljj_2p5);
        newtree->Branch("GENmass4ljj_2p5", &GENmass4ljj_2p5);
        newtree->Branch("GEN_TauC_Inc_0j_EnergyWgt", &GEN_TauB_Inc_0j_pTWgt);
        newtree->Branch("GEN_TauB_Inc_0j_pTWgt", &GEN_TauB_Inc_0j_pTWgt);

    }

    std::set<TString> runlumieventSet;
    for (Long64_t i=0;i<nentries; i++) {
        //lep_dataMC_new.clear();
        if (_addNewVar && _OldNtuple)
        {
            order_index.clear();
/*
	    pt_nom.clear();
	    eta_nom.clear();
	    phi_nom.clear();
	    mass_nom.clear();
*/

            pTj2=-1.0; mj1j2=-1.0; dEtaj1j2=-99.0; dPhij1j2=-99.0;
            pT4lj=-1.0; mass4lj=-1.0; pT4ljj=-1.0; mass4ljj=-1.0;
            GENpTj2=-1.0; GENmj1j2=-1.0; GENdEtaj1j2=-99.0; GENdPhij1j2=-99.0;
            GENpT4lj=-1.0; GENmass4lj=-1.0; GENpT4ljj=-1.0; GENmass4ljj=-1.0;
            GENpTj2_2p5=-1.0; GENmj1j2_2p5=-1.0; GENdEtaj1j2_2p5=-99.0; GENdPhij1j2_2p5=-99.0;
            GENpT4lj_2p5=-1.0; GENmass4lj_2p5=-1.0; GENpT4ljj_2p5=-1.0; GENmass4ljj_2p5=-1.0;
	    GEN_TauB_Inc_0j_pTWgt=-99, GEN_TauC_Inc_0j_EnergyWgt=-99;

        }
        //if (i>=2000000) continue;
        if (i>=20000&&_Test) break;
        if (i%100000==0) std::cout<<i<<"/"<<nentries<<std::endl;
        oldtree->GetEntry(i);

        if (_addNewVar && _OldNtuple)
        {

           if(njets_pt30_eta4p7>0)
            {
                j1.SetPtEtaPhiM((*jet_pt)[jet1index],(*jet_eta)[jet1index],(*jet_phi)[jet1index],(*jet_mass)[jet1index]);
            }
            if(njets_pt30_eta4p7>1)
            {
                j2.SetPtEtaPhiM((*jet_pt)[jet2index],(*jet_eta)[jet2index],(*jet_phi)[jet2index],(*jet_mass)[jet2index]);
                jj = j1 + j2;
                pTj2 = j2.Pt();
                mj1j2 = jj.M();
                dEtaj1j2 = j1.Eta()-j2.Eta();
                dPhij1j2 = j1.DeltaPhi(j2);
            }
            if(GENnjets_pt30_eta4p7>0)
            {
                order_pt(GENjet_pt, order_index);
                int GENjet1index = order_index[0];
                GENj1.SetPtEtaPhiM((*GENjet_pt)[GENjet1index],(*GENjet_eta)[GENjet1index],(*GENjet_phi)[GENjet1index],(*GENjet_mass)[GENjet1index]);
            }
            if(GENnjets_pt30_eta4p7>1)
            {
                int GENjet2index = order_index[1];
                GENj2.SetPtEtaPhiM((*GENjet_pt)[GENjet2index],(*GENjet_eta)[GENjet2index],(*GENjet_phi)[GENjet2index],(*GENjet_mass)[GENjet2index]);
                GENjj = GENj1 + GENj2;

                GENpTj2 = GENj2.Pt();
                GENmj1j2 = GENjj.M();
                GENdEtaj1j2 = GENj1.Eta()-GENj2.Eta();
                GENdPhij1j2 = GENj1.DeltaPhi(GENj2);
            }
// 2p5
            if(GENnjets_pt30_eta2p5>0)
            {
                order_pt(GENjet_pt, order_index);
                int GENjet1index_2p5 = order_index[0];
                GENj1_2p5.SetPtEtaPhiM((*GENjet_pt)[GENjet1index_2p5],(*GENjet_eta)[GENjet1index_2p5],(*GENjet_phi)[GENjet1index_2p5],(*GENjet_mass)[GENjet1index_2p5]);
            }
            if(GENnjets_pt30_eta2p5>1)
            {
                int GENjet2index_2p5 = order_index[1];
                GENj2_2p5.SetPtEtaPhiM((*GENjet_pt)[GENjet2index_2p5],(*GENjet_eta)[GENjet2index_2p5],(*GENjet_phi)[GENjet2index_2p5],(*GENjet_mass)[GENjet2index_2p5]);
                GENjj_2p5 = GENj1_2p5 + GENj2_2p5;

                GENpTj2_2p5 = GENj2_2p5.Pt();
                GENmj1j2_2p5 = GENjj_2p5.M();
                GENdEtaj1j2_2p5 = GENj1_2p5.Eta()-GENj2_2p5.Eta();
                GENdPhij1j2_2p5 = GENj1_2p5.DeltaPhi(GENj2_2p5);
            }





            if(passedFiducialSelection)
            {
                for (unsigned int k=0; k<(*GENjet_pt).size(); k++) {

                    TLorentzVector thisGENJet;
                    thisGENJet.SetPtEtaPhiM((*GENjet_pt)[k],(*GENjet_eta)[k],(*GENjet_phi)[k],(*GENjet_mass)[k]);

                    if ((*GENjet_pt)[k]<30.0 || abs((*GENjet_eta)[k])>4.7) continue;
                    
                    pt_nom.push_back(thisGENJet.Pt());
                    eta_nom.push_back(thisGENJet.Eta());
                    phi_nom.push_back(thisGENJet.Phi());
                    mass_nom.push_back(thisGENJet.M());
                } 


                GENl1.SetPtEtaPhiM((*GENlep_pt)[GENlep_Hindex[0]],(*GENlep_eta)[GENlep_Hindex[0]],(*GENlep_phi)[GENlep_Hindex[0]],(*GENlep_mass)[GENlep_Hindex[0]]);
                GENl2.SetPtEtaPhiM((*GENlep_pt)[GENlep_Hindex[1]],(*GENlep_eta)[GENlep_Hindex[1]],(*GENlep_phi)[GENlep_Hindex[1]],(*GENlep_mass)[GENlep_Hindex[1]]);
                GENl3.SetPtEtaPhiM((*GENlep_pt)[GENlep_Hindex[2]],(*GENlep_eta)[GENlep_Hindex[2]],(*GENlep_phi)[GENlep_Hindex[2]],(*GENlep_mass)[GENlep_Hindex[2]]);
                GENl4.SetPtEtaPhiM((*GENlep_pt)[GENlep_Hindex[3]],(*GENlep_eta)[GENlep_Hindex[3]],(*GENlep_phi)[GENlep_Hindex[3]],(*GENlep_mass)[GENlep_Hindex[3]]);

                GENHiggs = GENl1 + GENl2 + GENl3+ GENl4;

                GEN_TauC_Inc_0j_EnergyWgt = TauC(pt_nom, eta_nom, phi_nom, mass_nom, GENHiggs);
                GEN_TauB_Inc_0j_pTWgt = TauB(pt_nom, eta_nom, phi_nom, mass_nom, GENHiggs);

                pt_nom.clear();
                eta_nom.clear();
                phi_nom.clear();
                mass_nom.clear();


                Float_t GENmass4l_ = GENHiggs.M();
                //if((GENmass4l-GENmass4l_)>0.01)
                //    cout<<"GENmass4l: "<<GENmass4l<<"; GENmass4l_: "<<GENmass4l_<<endl;
                if(GENnjets_pt30_eta4p7>0)
                {
                    GENhj = GENj1 + GENHiggs;
                    GENpT4lj = GENhj.Pt();
                    GENmass4lj = GENhj.M();
                }
                if(GENnjets_pt30_eta4p7>1)
                {
                    GENhjj = GENjj + GENHiggs;
                    GENpT4ljj = GENhjj.Pt();
                    GENmass4ljj = GENhjj.M();
                }
//2p5
                if(GENnjets_pt30_eta2p5>0)
                {
                    GENhj_2p5 = GENj1_2p5 + GENHiggs;
                    GENpT4lj_2p5 = GENhj_2p5.Pt();
                    GENmass4lj_2p5 = GENhj_2p5.M();
                }
                if(GENnjets_pt30_eta2p5>1)
                {
                    GENhjj_2p5 = GENjj_2p5 + GENHiggs;
                    GENpT4ljj_2p5 = GENhjj_2p5.Pt();
                    GENmass4ljj_2p5 = GENhjj_2p5.M();
                }


            }

            if(passedFullSelection)
            {
                Higgs.SetPtEtaPhiM(pT4l, eta4l, phi4l, mass4l);
                if(njets_pt30_eta4p7>0)
                {
                    hj = j1 + Higgs;

                    pT4lj = hj.Pt();
                    mass4lj = hj.M();
                }
                if(njets_pt30_eta4p7>1)
                {
                    //jj = j1 + j2;
                    hjj = jj + Higgs;

                    pT4ljj = hjj.Pt();
                    mass4ljj = hjj.M();
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

