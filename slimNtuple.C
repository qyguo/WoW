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

void slimNtuple(const int & _year_=2017, const string & _name_DS_="GluGluHToZZTo4L_M125_2017") {
    int year = _year_;
    string pre_name = "root://cmsio5.rc.ufl.edu:1094//store/user/t2/users/ferrico/Qianying/";
    string whichYear[3] = {"2016/", "2017/", "2018/"};
    string MC__ = "MC/";
    //string name_DS = "GluGluHToZZTo4L_M125_2017";
    string name_DS = _name_DS_.c_str();

    //TString filename = prefix+".root";
    string filename = (pre_name + whichYear[year-2016] + MC__ + name_DS + ".root").c_str();
    //filename = "Sync_106X_2018UL_30112021_-1ggH_V2.root";

    std::cout<<"Year: "<<year<<std::endl;
    std::cout<<filename<<std::endl;

    TFile *oldfile = TFile::Open(filename.c_str());
    oldfile->cd("Ana");
    //TTree *oldtree = (TTree*)oldfile->Get("Ana/passedEvents");
    TTree *oldtree = (TTree*)gDirectory->Get("passedEvents");

    TH1F *th[7];
    th[0]=(TH1F*)oldfile->Get("Ana/nEvents");
    th[1]=(TH1F*)oldfile->Get("Ana/sumWeights");
    th[2]=(TH1F*)oldfile->Get("Ana/sumWeightsPU");
    th[3]=(TH1F*)oldfile->Get("Ana/nVtx");
    th[4]=(TH1F*)oldfile->Get("Ana/nVtx_ReWeighted");
    th[5]=(TH1F*)oldfile->Get("Ana/nInteractions");
    th[6]=(TH1F*)oldfile->Get("Ana/nInteraction_ReWeighted");

    Long64_t nentries = oldtree->GetEntries();
    std::cout<<nentries<<" total entries."<<std::endl;
    ULong64_t Run, LumiSect, Event;
    Bool_t          passedFullSelection;
    Float_t         dataMCWeight;
    Float_t         eventWeight;
    vector<int>     *lep_id;
    vector<float>   *lep_pt;
    vector<float>   *lep_eta;
    Int_t           lep_Hindex[4];
    vector<float>   *lep_dataMC;
    vector<float>   lep_dataMC_new;
    vector<float>   *lep_dataMCErr;
    lep_id = 0; 
    lep_pt = 0;    
    lep_eta = 0;  
    lep_dataMC = 0;
    lep_dataMCErr = 0;
    oldtree->SetBranchAddress("Run",&Run);
    oldtree->SetBranchAddress("LumiSect",&LumiSect);
    oldtree->SetBranchAddress("Event",&Event);
    oldtree->SetBranchAddress("passedFullSelection", &passedFullSelection);
    oldtree->SetBranchAddress("dataMCWeight", &dataMCWeight);
    oldtree->SetBranchAddress("eventWeight", &eventWeight);
    oldtree->SetBranchAddress("lep_id", &lep_id);
    oldtree->SetBranchAddress("lep_pt", &lep_pt);
    oldtree->SetBranchAddress("lep_Hindex", lep_Hindex);
    oldtree->SetBranchAddress("lep_eta", &lep_eta);
    oldtree->SetBranchAddress("lep_dataMC", &lep_dataMC);
    oldtree->SetBranchAddress("lep_dataMCErr", &lep_dataMCErr);

    oldtree->SetBranchStatus("*",0); //Disables All Branches
    //Then enables only select branches
    oldtree->SetBranchStatus("Run",1);
    oldtree->SetBranchStatus("LumiSect",1);
    oldtree->SetBranchStatus("Event",1);
    oldtree->SetBranchStatus("passedFullSelection",1);
    oldtree->SetBranchStatus("eventWeight",1);
    oldtree->SetBranchStatus("passedFiducialSelection",1);
    oldtree->SetBranchStatus("mass4l",1);
    oldtree->SetBranchStatus("pT4l",1);
    oldtree->SetBranchStatus("eta4l",1);
    oldtree->SetBranchStatus("phi4l",1);
    oldtree->SetBranchStatus("Phi",1);
    oldtree->SetBranchStatus("Phi1",1);
    oldtree->SetBranchStatus("cosTheta1",1);
    oldtree->SetBranchStatus("cosTheta2",1);
    oldtree->SetBranchStatus("cosThetaStar",1);
    oldtree->SetBranchStatus("lep_pt",1);
    oldtree->SetBranchStatus("lep_eta",1);

    //Create a new file + a clone of old tree in new file
    TFile *newfile = new TFile(
            (name_DS+"_slimmed.root").c_str()
            ,"recreate");
    cout<<"Output file: "<<newfile->GetName()<<endl;
    newfile->mkdir("Ana");
    newfile->cd("Ana");
    for (int ii=0; ii<7; ii++)
    {
      th[ii]->Write();
    }
    TTree *newtree = oldtree->CloneTree(0);
    //newtree->Branch("dataMCWeight_new", &dataMCWeight_new);
    //newtree->Branch("eventWeight_new", &eventWeight_new);

    std::set<TString> runlumieventSet;
    for (Long64_t i=0;i<nentries; i++) {
        lep_dataMC_new.clear();
        //if (i>=2000000) continue;
        if (i%100000==0) std::cout<<i<<"/"<<nentries<<std::endl;
        oldtree->GetEntry(i);
        if (passedFullSelection)
        {
            newtree->Fill();////if only store the event passedFullSelection
        }
        
        //newtree->Fill();
    }

    newtree->AutoSave();
    newfile->Close();
    oldfile->Close();
    //delete oldfile;
    delete newfile;
    exit(EXIT_FAILURE);
}

