#include <TFile.h>
#include <TTree.h>
#include <iostream>

void createdata(bool useRefit=false) {

   //TFile *oldfile = new TFile("Data_Mar14_bestCandLegacy_NoDuplicates.root","READ");
   TFile *oldfile = new TFile("data_UL2016_noDuplicates.root","READ");
   //TFile *oldfile = new TFile("/publicfs/cms/data/hzz/guoqy/newNTuple_UL/2017/Data/DataUL2017_all_noDuplicates.root","READ");
   //TFile *oldfile = new TFile("/publicfs/cms/data/hzz/guoqy/newNTuple_UL/2018/Data/DataUL2018_all_noDuplicates.root","READ");

   TTree *oldtree = (TTree*)oldfile->Get("passedEvents");

   ULong64_t Run;
   oldtree->SetBranchAddress("Run",&Run);
   ULong64_t LumiSect;
   oldtree->SetBranchAddress("LumiSect",&LumiSect);
   ULong64_t Event;
   oldtree->SetBranchAddress("Event",&Event);

   Bool_t passedFullSelection;
   oldtree->SetBranchAddress("passedFullSelection",&passedFullSelection);

   Float_t mass4l;
   oldtree->SetBranchAddress("mass4l",&mass4l);
   Float_t mass4lErr;
   oldtree->SetBranchAddress("mass4lErr",&mass4lErr);
   Float_t mass4lREFIT;
   oldtree->SetBranchAddress("mass4lREFIT",&mass4lREFIT);
   Float_t mass4lErrREFIT;
   oldtree->SetBranchAddress("mass4lErrREFIT",&mass4lErrREFIT);
   Float_t mass4e;
   oldtree->SetBranchAddress("mass4e",&mass4e);
   Float_t mass4mu;
   oldtree->SetBranchAddress("mass4mu",&mass4mu);
   Float_t mass2e2mu;
   oldtree->SetBranchAddress("mass2e2mu",&mass2e2mu);

   Float_t pT4l;
   oldtree->SetBranchAddress("pT4l",&pT4l);
   Float_t rapidity4l;
   oldtree->SetBranchAddress("rapidity4l",&rapidity4l);
   Int_t njets_pt30_eta4p7;
   oldtree->SetBranchAddress("njets_pt30_eta4p7",&njets_pt30_eta4p7);
   Int_t njets_pt30_eta2p5;
   oldtree->SetBranchAddress("njets_pt30_eta2p5",&njets_pt30_eta2p5);
   Float_t pt_leadingjet_pt30_eta4p7;
   oldtree->SetBranchAddress("pt_leadingjet_pt30_eta4p7",&pt_leadingjet_pt30_eta4p7);
   Float_t pt_leadingjet_pt30_eta2p5;
   oldtree->SetBranchAddress("pt_leadingjet_pt30_eta2p5",&pt_leadingjet_pt30_eta2p5);
// additional
   Float_t massZ1;
   oldtree->SetBranchAddress("massZ1",&massZ1);
   Float_t massZ2;
   oldtree->SetBranchAddress("massZ2",&massZ2);
   Float_t cosThetaStar;
   oldtree->SetBranchAddress("cosThetaStar",&cosThetaStar);
   Float_t cosTheta1;
   oldtree->SetBranchAddress("cosTheta1",&cosTheta1);
   Float_t cosTheta2;
   oldtree->SetBranchAddress("cosTheta2",&cosTheta2);
   Float_t Phi;
   oldtree->SetBranchAddress("Phi",&Phi);
   Float_t Phi1;
   oldtree->SetBranchAddress("Phi1",&Phi1);
//4p7
   Float_t pTj2;
   oldtree->SetBranchAddress("pTj2",&pTj2);
   Float_t mj1j2;
   oldtree->SetBranchAddress("mj1j2",&mj1j2);
   Float_t dEtaj1j2;
   oldtree->SetBranchAddress("dEtaj1j2",&dEtaj1j2);
   Float_t dPhij1j2;
   oldtree->SetBranchAddress("dPhij1j2",&dPhij1j2);
   Float_t pT4lj;
   oldtree->SetBranchAddress("pT4lj",&pT4lj);
   Float_t mass4lj;
   oldtree->SetBranchAddress("mass4lj",&mass4lj);
   Float_t pT4ljj;
   oldtree->SetBranchAddress("pT4ljj",&pT4ljj);
   Float_t mass4ljj;
   oldtree->SetBranchAddress("mass4ljj",&mass4ljj);
//2p5
   Float_t pTj2_2p5;
   oldtree->SetBranchAddress("pTj2_2p5",&pTj2_2p5);
   Float_t mj1j2_2p5;
   oldtree->SetBranchAddress("mj1j2_2p5",&mj1j2_2p5);
   Float_t dEtaj1j2_2p5;
   oldtree->SetBranchAddress("dEtaj1j2_2p5",&dEtaj1j2_2p5);
   Float_t dPhij1j2_2p5;
   oldtree->SetBranchAddress("dPhij1j2_2p5",&dPhij1j2_2p5);
   Float_t pT4lj_2p5;
   oldtree->SetBranchAddress("pT4lj_2p5",&pT4lj_2p5);
   Float_t mass4lj_2p5;
   oldtree->SetBranchAddress("mass4lj_2p5",&mass4lj_2p5);
   Float_t pT4ljj_2p5;
   oldtree->SetBranchAddress("pT4ljj_2p5",&pT4ljj_2p5);
   Float_t mass4ljj_2p5;
   oldtree->SetBranchAddress("mass4ljj_2p5",&mass4ljj_2p5);
//KD
   Float_t D_0m;
   oldtree->SetBranchAddress("D_0m",&D_0m);
   Float_t D_CP;
   oldtree->SetBranchAddress("D_CP",&D_CP);
   Float_t D_0hp;
   oldtree->SetBranchAddress("D_0hp",&D_0hp);
   Float_t D_int;
   oldtree->SetBranchAddress("D_int",&D_int);
   Float_t D_L1;
   oldtree->SetBranchAddress("D_L1",&D_L1);
   Float_t D_L1Zg;
   oldtree->SetBranchAddress("D_L1Zg",&D_L1Zg);
//Tau
   Float_t TauC_Inc_0j_EnergyWgt;
   oldtree->SetBranchAddress("TauC_Inc_0j_EnergyWgt",&TauC_Inc_0j_EnergyWgt);
   Float_t TauB_Inc_0j_pTWgt;
   oldtree->SetBranchAddress("TauB_Inc_0j_pTWgt",&TauB_Inc_0j_pTWgt);

   Int_t finalState;
   oldtree->SetBranchAddress("finalState",&finalState);

   Float_t D_bkg_kin;
   oldtree->SetBranchAddress("D_bkg_kin",&D_bkg_kin);
  
   std::vector<ULong64_t> runVec, lumiVec, eventVec;

   TFile *newfile = new TFile("data_UL2016_noDuplicates_created.root","recreate");
   //TFile *newfile = new TFile("data_UL2017_noDuplicates_created.root","recreate");
   //TFile *newfile = new TFile("data_UL2018_noDuplicates_created.root","recreate");
   TTree *newtree = new TTree("passedEvents","passedEvents");
   
   Float_t CMS_zz4l_mass;
   Float_t CMS_zz4l_massErr;
   Float_t melaLD;
   newtree->Branch("CMS_zz4l_mass",&CMS_zz4l_mass);
   newtree->Branch("CMS_zz4l_massErr",&CMS_zz4l_massErr);
   newtree->Branch("melaLD",&melaLD);
   newtree->Branch("mass4e",&mass4e);
   newtree->Branch("mass4mu",&mass4mu);
   newtree->Branch("mass2e2mu",&mass2e2mu);
   newtree->Branch("pT4l",&pT4l);
   newtree->Branch("rapidity4l",&rapidity4l);
   newtree->Branch("njets_pt30_eta4p7",&njets_pt30_eta4p7);
   newtree->Branch("njets_pt30_eta2p5",&njets_pt30_eta2p5);
   newtree->Branch("pt_leadingjet_pt30_eta4p7",&pt_leadingjet_pt30_eta4p7);
   newtree->Branch("pt_leadingjet_pt30_eta2p5",&pt_leadingjet_pt30_eta2p5);
//additional
   newtree->Branch("massZ1",&massZ1);
   newtree->Branch("massZ2",&massZ2);
   newtree->Branch("rapidity4l",&rapidity4l);
   newtree->Branch("cosThetaStar",&cosThetaStar);
   newtree->Branch("cosTheta1",&cosTheta1);
   newtree->Branch("cosTheta2",&cosTheta2);
   newtree->Branch("Phi",&Phi);
   newtree->Branch("Phi1",&Phi1);
//4p7
   newtree->Branch("pTj2",&pTj2);
   newtree->Branch("mj1j2",&mj1j2);
   newtree->Branch("dEtaj1j2",&dEtaj1j2);
   newtree->Branch("dPhij1j2",&dPhij1j2);
   newtree->Branch("pT4lj",&pT4lj);
   newtree->Branch("mass4lj",&mass4lj);
   newtree->Branch("pT4ljj",&pT4ljj);
   newtree->Branch("mass4ljj",&mass4ljj);
//2p5
   newtree->Branch("pTj2_2p5",&pTj2_2p5);
   newtree->Branch("mj1j2_2p5",&mj1j2_2p5);
   newtree->Branch("dEtaj1j2_2p5",&dEtaj1j2_2p5);
   newtree->Branch("dPhij1j2_2p5",&dPhij1j2_2p5);
   newtree->Branch("pT4lj_2p5",&pT4lj_2p5);
   newtree->Branch("mass4lj_2p5",&mass4lj_2p5);
   newtree->Branch("pT4ljj_2p5",&pT4ljj_2p5);
   newtree->Branch("mass4ljj_2p5",&mass4ljj_2p5);
//KD
   newtree->Branch("D_0m",&D_0m);
   newtree->Branch("D_CP",&D_CP);
   newtree->Branch("D_0hp",&D_0hp);
   newtree->Branch("D_int",&D_int);
   newtree->Branch("D_L1",&D_L1);
   newtree->Branch("D_CP",&D_CP);
   newtree->Branch("D_L1Zg",&D_L1Zg);
// Tau
   newtree->Branch("TauC_Inc_0j_EnergyWgt",&TauC_Inc_0j_EnergyWgt);
   newtree->Branch("TauB_Inc_0j_pTWgt",&TauB_Inc_0j_pTWgt);

   newtree->Branch("finalState",&finalState);

   int n4e = 0;
   int n2e2mu = 0;
   int n4mu = 0;

   int nEntries = oldtree->GetEntries();

   std::cout<<"using refit? "<<useRefit<<std::endl;

   for(int i = 0; i < oldtree->GetEntries(); i++) {

       oldtree->GetEntry(i);

       if (i%1000==0) std::cout<<i<<"/"<<nEntries<<std::endl;

       if (!passedFullSelection) continue;

       if (useRefit) {
           CMS_zz4l_mass = mass4lREFIT;
           CMS_zz4l_massErr = mass4lErrREFIT;           
           if (mass4e>0.0) mass4e = mass4lREFIT;
           if (mass4mu>0.0) mass4mu = mass4lREFIT;
           if (mass2e2mu>0.0) mass2e2mu = mass4lREFIT;
       }
       else {
           CMS_zz4l_mass = mass4l;
           CMS_zz4l_massErr = mass4lErr;
       }

       melaLD = D_bkg_kin;

       if (mass4e>105.0 && mass4e<160.0) n4e++;
       if (mass4mu>105.0 && mass4mu<160.0) n4mu++;
       if (mass2e2mu>105.0 && mass2e2mu<160.0) n2e2mu++;
      
       newtree->Fill();

   }

   std::cout<<"105.0-160.0: "<<n4e<<" 4e, "<<n4mu<<" 4mu, "<<n2e2mu<<" 2e2mu."<<std::endl;

   newfile->Write();

}
