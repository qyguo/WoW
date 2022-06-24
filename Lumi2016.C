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

string Path_(string file_name, string str="_")
{
    string path = file_name.c_str();
    //auto pos = path.find_last_of ('_');
    auto pos = path.rfind(str); ///reverse find
    string dir;
    dir.clear();
    if (pos != string::npos) {
        auto directory = path.substr (0,pos);
        auto file = path.substr (pos);
        dir = directory;
    }
    return dir;
}

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

void Lumi2016(const int & _year_=2017, const string & _name_DS_="bbH_HToZZTo4L_M125_TuneCP2_13TeV-jhugenv7011-pythia8", const bool & isMC = true,  const bool & isSignal = true, const bool & _Test=false) {
    std::cout<<"year: "<<_year_<<"; name_DS: "<<_name_DS_<<"; isMC: "<<isMC<<"; isSignal: "<<isSignal<<"; _Test: "<<_Test<<std::endl;

    int year = _year_;
    int Year = TMath::Abs(year);
    if (Year!=2016)    continue;

    string pre_name = "/publicfs/cms/data/hzz/guoqy/newNTuple_UL/";
    if (year==2016)    pre_name="/eos/cms/store/group/phys_muon/TagAndProbe/HZZ4L/2016/UL/MC/postVFP/";
    else if (year==-2016) //UL16APV
        pre_name="/eos/cms/store/group/phys_muon/TagAndProbe/HZZ4L/2016/UL/MC/preVFP/";

    string whichYear[3] = {"2016/", "2017/", "2018/"};
    string MC__ = "MC/";
    //string files=Path_(_name_DS_);
    //string MC__ = (Path_(Path_(_name_DS_, "_slimmed_newMuSF_"))+"/").c_str();
    if (!isMC) MC__ = "Data/";
    //string name_DS = "GluGluHToZZTo4L_M125_2017";
    string name_DS = _name_DS_.c_str();
    MC__="";

    //string filename = (pre_name + whichYear[year-2016] + MC__ + "Slimmed_2p5_JES/" + name_DS + ".root").c_str();
    //string outputFile_Slimmed = (pre_name + whichYear[year-2016] + MC__  + "Slimmed_2p5_JES_4More/").c_str(); 
    string filename = (pre_name + name_DS + ".root").c_str();
    string outputFile_Slimmed = (pre_name).c_str();

    std::cout<<"Year: "<<year<<std::endl;
    std::cout<<filename<<std::endl;

    TFile *oldfile = TFile::Open(filename.c_str());
    if (isMC) oldfile->cd("Ana");
    TTree *oldtree = (TTree*)gDirectory->Get("passedEvents");
    TH1F *th[7];
    TH1F *th_LimiWeight2016;
    th_LimiWeight2016=new TH1F("Lumi_Weight","lumi weight for 2016",2,0,2);

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

    Float_t Lumi_Weight;
    Float_t Lumi_Curr;
    Float_t TotalLumi_2016;
    if (year==2016)    Lumi_Curr=16.81;
    else if (year=-2016)    Lumi_Curr=19.52;
    TotalLumi_2016 = 16.81+19.52;
    Lumi_Weight=Lumi_Curr/TotalLumi_2016;
    oldtree->SetBranchStatus("*",1);

    th_LimiWeight2016->SetBinContent(1,Lumi_Weight);

    //Create a new file + a clone of old tree in new file
    TFile *newfile = new TFile(
            (outputFile_Slimmed+name_DS+"_slimmed_newMuSF_add2p5_LumiWeight.root").c_str()
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
            th_LimiWeight2016->Write();
        }
    }
    TTree *newtree = oldtree->CloneTree(0);
    newtree->Branch("Lumi_Weight", &Lumi_Weight);

    std::set<TString> runlumieventSet;
    for (Long64_t i=0;i<nentries; i++) {
        if (i%100000==0) std::cout<<i<<"/"<<nentries<<std::endl;
        oldtree->GetEntry(i);
        newtree->Fill();
    }
    newtree->AutoSave();
    newfile->Close();
    oldfile->Close();
    //delete oldfile;
    delete newfile;
    exit(EXIT_FAILURE);
}
