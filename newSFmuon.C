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
    return hMuScaleFac->GetBinContent(hMuScaleFac->FindBin(eta,pt));
}

float dataMCErr_(float pt_, float eta_, TH2F* hMuScaleFacUnc)
{
    float pt = TMath::Min(pt_,(float)199.0);
    float eta = eta_;
    return hMuScaleFacUnc->GetBinContent(hMuScaleFacUnc->FindBin(eta,pt));
}

void newSFmuon(const int & _year_=2017, const string & _name_DS_="bbH_HToZZTo4L_M125_TuneCP2_13TeV-jhugenv7011-pythia8") {
    int year = _year_;
    int Year = TMath::Abs(year);
    //string pre_name = "root://cmsio5.rc.ufl.edu:1094//store/user/t2/users/ferrico/Qianying/";
    string whichYear[3] = {"2016/", "2017/", "2018/"};
    string MC__ = "MC/";
    //string name_DS = "GluGluHToZZTo4L_M125_2017";
    string name_DS = _name_DS_.c_str();

    //TString filename = prefix+".root";
    //string filename = (pre_name + whichYear[year-2016] + MC__ + name_DS + ".root").c_str();
    //pre_name = "/raid/raid9/qguo/Run2/after/Run2_2/new/CMSSW_10_2_18/src/resultsAna_2018_ggHMad/crab_ggH_amcatnloFXFX/results/";
    string pre_name = "/publicfs/cms/data/hzz/guoqy/newNTuple_UL/" + whichYear[year-2016] + "/MC/";
    if (year==2016)    pre_name="";
    else if (year==-2016) //UL16APV
        pre_name="";
    string filename = (pre_name + name_DS + ".root").c_str();

    std::cout<<"year: "<<year<<std::endl;
    std::cout<<filename<<std::endl;

    //old
    //string mu_scalefac_name_161718[3] = {"final_HZZ_SF_2016_legacy_mupogsysts.root", "final_HZZ_SF_2017_rereco_mupogsysts_3010.root", "final_HZZ_SF_2018_rereco_mupogsysts_3010.root"};
    //new
    //string mu_scalefac_name_161718[3] = {"final_HZZ_muon_SF_2016RunB2H_legacy_newLoose_newIso_paper.root", "final_HZZ_muon_SF_2017_newLooseIso_mupogSysts_paper.root", "final_HZZ_muon_SF_2018RunA2D_ER_newLoose_newIso_paper.root"};
    //new UL
    string mu_scalefac_name_161718[3] = {"final_HZZ_SF_2016UL_mupogsysts_newLoose.root", "final_HZZ_SF_2017UL_mupogsysts_newLoose.root", "final_HZZ_SF_2018UL_mupogsysts_newLoose.root"};
    //edm::FileInPath mu_scalefacFileInPath(("UFHZZAnalysisRun2/UFHZZ4LAna/data/"+mu_scalefac_name_161718[year-2016]).c_str());

    //string mu_scalefacFileInPath = ("/raid/raid9/qguo/Run2/after/Run2_2/test3/CMSSW_10_2_18/src/UFHZZAnalysisRun2/UFHZZ4LAna/data/"+mu_scalefac_name_161718[year-2016]).c_str();
    string mu_scalefacFileInPath = ("/publicfs/cms/data/hzz/guoqy/newNTuple_UL/newMuonSF/" + mu_scalefac_name_161718[Year-2016]).c_str();
    TFile *fMuScalFac = TFile::Open(mu_scalefacFileInPath.c_str());
    TH2F *hMuScaleFac = (TH2F*)fMuScalFac->Get("FINAL");
    TH2F *hMuScaleFacUnc = (TH2F*)fMuScalFac->Get("ERROR");

    TFile *oldfile = TFile::Open(filename.c_str());
    oldfile->cd("Ana");
    //TTree *oldtree = (TTree*)oldfile->Get("Ana/passedEvents");
    
    TTree *oldtree = (TTree*)gDirectory->Get("passedEvents");
    //TTree *oldtree = (TTree*)oldfile->Get("passedEvents");

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
    vector<float>   *lep_dataMCErr;
    vector<float>   lep_dataMC_new;
    vector<float>   lep_dataMCErr_new;
    Float_t         crossSection;
    Float_t         dataMCWeight_new;
    Float_t         eventWeight_new;
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
    oldtree->SetBranchAddress("crossSection", &crossSection);

    //Create a new file + a clone of old tree in new file
    TFile *newfile = new TFile(
            (pre_name+"/WithNewMuonSF/"+name_DS+"_newMuonSF_UL.root").c_str()
            ,"recreate");
    cout<<"Output file: "<<newfile->GetName()<<endl;
    newfile->mkdir("Ana");
    newfile->cd("Ana");
    TTree *newtree = oldtree->CloneTree(0);
    newtree->Branch("lep_dataMC_new", &lep_dataMC_new);
    newtree->Branch("lep_dataMCErr_new", &lep_dataMCErr_new);
    newtree->Branch("dataMCWeight_new", &dataMCWeight_new);
    newtree->Branch("eventWeight_new", &eventWeight_new);
    //newtree->Branch("crossSection", &crossSection);

    std::set<TString> runlumieventSet;
    for (Long64_t i=0;i<nentries; i++) {
        //lep_id = 0; 
        //lep_pt = 0;    
        //lep_eta = 0;  
        //lep_dataMC = 0;
        //lep_dataMCErr = 0;
        lep_dataMC_new.clear();
        lep_dataMCErr_new.clear();
        //if (i>70000000) continue;
        if (i%100000==0) std::cout<<i<<"/"<<nentries<<std::endl;
        //if (i>30000001) continue;
        oldtree->GetEntry(i);
        //Float_t crossSection_new = 49.27*0.0002502;//124 ggH
        //Float_t crossSection_new = 47.89*0.0003001;//126ggH
        //Float_t crossSection_old = crossSection;
        //crossSection = crossSection_new;

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
            //dataMCWeight = lep_dataMC[lep_Hindex[0]]*lep_dataMC[lep_Hindex[1]]*lep_dataMC[lep_Hindex[2]]*lep_dataMC[lep_Hindex[3]];
            //eventWeight = genWeight*crossSection*pileupWeight*dataMCWeight*prefiringWeight;
            int sum_lep_Hindex=lep_Hindex[0]+lep_Hindex[1]+lep_Hindex[2]+lep_Hindex[3];//lep_Hindex default value is -1
            if (sum_lep_Hindex!=-4)
            //if (passedFullSelection)
            {
                //eventWeight = eventWeight/crossSection_old*crossSection;
                dataMCWeight_new = lep_dataMC_new[lep_Hindex[0]]*lep_dataMC_new[lep_Hindex[1]]*lep_dataMC_new[lep_Hindex[2]]*lep_dataMC_new[lep_Hindex[3]]; 
                eventWeight_new = eventWeight/dataMCWeight*dataMCWeight_new;
            }
            else
            {

                dataMCWeight_new = 1.0;
                eventWeight_new = eventWeight/dataMCWeight*dataMCWeight_new;
            }

        
        newtree->Fill();
    }

    //newtree->Print();
    newtree->AutoSave();
    //delete oldfile;
    delete newfile;
    fMuScalFac->Close();
    exit(EXIT_FAILURE);
}

