#include "hades.h"
#include "htool.h"
#include "hloop.h"
#include "hphysicsconstants.h"
#include "hrootsource.h"
#include "hiterator.h"

#include "TMacro.h"

#include "haddef.h"
#include "hgeantdef.h"

#include "hparticledef.h"
#include "hparticlestructs.h"
#include "hparticlecand.h"
#include "hparticlecandsim.h"
#include "hcategorymanager.h"
#include "hparticletracksorter.h"
#include "hgeantkine.h"
#include "hparticleevtinfo.h"
#include "hparticletool.h"
#include "henergylosscorrpar.h"
#include "hparticlebtring.h"

#include "hmdcdef.h"
#include "hmdctrackddef.h"
#include "hmdctrackgdef.h"
#include "horadef.h"
#include "horasimdef.h"
#include "hstartdef.h"
#include "richdef.h"
#include "rpcdef.h"
#include "showerdef.h"
#include "simulationdef.h"
#include "tofdef.h"
#include "walldef.h"

#include "TROOT.h"
#include "TSystem.h"
#include "TString.h"
#include "TTree.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH3F.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TFile.h"
#include "TStopwatch.h"
#include "TCutG.h"
#include "TF1.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TNtuple.h"

#include <iostream>
#include <vector>
#include <algorithm>
#include <map>
#include <stdio.h>

using namespace std;

Bool_t isGoodNonFitted(HParticleCand* cand, Float_t nonFittedCut4){
    // return kTRUE if candidate has no non fitted
    // neighbour closer than cut
    if(     cand->getAngleToNearbyUnfittedInner()<-90) return kTRUE; // no non fitted in within 15 deg
    if(fabs(cand->getAngleToNearbyUnfittedInner()) > nonFittedCut4) return kTRUE;
    return kFALSE;
}
Bool_t isGoodFitted(HParticleCand* cand, Float_t fittedCut4){
    // return kTRUE if candidate has no  fitted
    // neighbour closer than cut
    if(     cand->getAngleToNearbyFittedInner()<-90) return kTRUE; // no fitted in within 15 deg
    if(fabs(cand->getAngleToNearbyFittedInner()) > fittedCut4) return kTRUE;
    return kFALSE;
}




Int_t dieleAna(TString, TString, Int_t );
int main(int argc,char **argv)
{
  
  cout << "argc = " << argc <<endl;
  switch (argc) {
           
  case 4:
      return  dieleAna(TString(argv[1]),TString(argv[2]),atoi(argv[3]));
      break;
   
  default:
      {
      cout<<"default"<<endl;
      return dieleAna("","",0);
      }
    break;
       
  }
  return 1;
}


Int_t dieleAna(TString inputlist, TString outfile, Int_t nev=-1)

{
    HLoop* loop = new HLoop(kTRUE);  // kTRUE : create Hades  (needed to work with standard eventstructure)


    TString readCategories = "-*,+HParticleCand,+HParticleEvtInfo,+HStart2Hit,+HStart2Cal,+HRichHit,+HParticleBtRing";
    loop->addMultFiles(inputlist);

    if(!loop->setInput(readCategories)) { exit(1); }
    loop->printCategories();

    HParticleEvtInfo* evtinfo;

    HParticleCand* cand1;
    HParticleCand* cand2;


    HParticleBtRing* btRing;
    HParticleBtRing* btRing2;



    HCategory* candCat = (HCategory*)HCategoryManager::getCategory(catParticleCand);
    if(!candCat) { exit(1); }

    HCategory* btCat = (HCategory*)HCategoryManager::getCategory(catParticleBtRing);
    if(!btCat)   { exit(1); }




    Int_t binsx    = 200;
    Int_t binsxmom = 300;
    Float_t momMax = 2000;
    Float_t momMin =-2000;
    //TH1::SetDefaultSumw2(kTRUE);
    //-------------------------------------------------------------------------------------------------
    // Define histograms
    //-------------------------------------------------------------------------------------------------
     
    TH1F* hCounter = new TH1F("hCounter","hCounter",5,0,5);
    hCounter->GetXaxis()->SetBinLabel(1,"all");
    hCounter->GetXaxis()->SetBinLabel(4,"events with lepton candidates");
   
    TH2F* htheta_phi_raw_pair= new TH2F("htheta_phi_raw_pair","",90,0,360,40,10,90);
    TH2F* htheta_phi_raw_pair_ele= new TH2F("htheta_phi_raw_pair_ele","",90,0,360,40,10,90);
    TH2F* htheta_phi_raw_pair_pos= new TH2F("htheta_phi_raw_pair_pos","",90,0,360,40,10,90);
    
    TH2F* htheta_phi_raw_single_ele= new TH2F("htheta_phi_raw_single_ele","",45,0,360,20,10,90);
    TH2F* htheta_phi_raw_single_pos= new TH2F("htheta_phi_raw_single_pos","",45,0,360,20,10,90);
    TH2F* htheta_phi_raw_single_ele_pid= new TH2F("htheta_phi_raw_single_ele_pid","",45,0,360,20,10,90);
    TH2F* htheta_phi_raw_single_pos_pid= new TH2F("htheta_phi_raw_single_pos_pid","",45,0,360,20,10,90);

    TH2F* htheta_phi_raw_single_ele_rebin= new TH2F("htheta_phi_raw_single_ele_rebin","",6,0,360,4,10,90);
    TH2F* htheta_phi_raw_single_pos_rebin= new TH2F("htheta_phi_raw_single_pos_rebin","",6,0,360,4,10,90);

    TH2F* htheta_phi_raw_pair_ele_rebin= new TH2F("htheta_phi_raw_pair_ele_rebin","",6,0,360,4,10,90);
    TH2F* htheta_phi_raw_pair_pos_rebin= new TH2F("htheta_phi_raw_pair_pos_rebin","",6,0,360,4,10,90);



    TH2F *hmassNP_missmass_vs_hmissmassNP_oa9deg = new TH2F ("hmassNP_missmass_vs_hmissmassNP_oa9deg","",1500,0,700,1500,600,1500);
    TH2F *hmassPP_missmass_vs_hmissmassPP_oa9deg = new TH2F ("hmassPP_missmass_vs_hmissmassPP_oa9deg","",1500,0,700,1500,600,1500);
    TH2F *hmassNN_missmass_vs_hmissmassNN_oa9deg = new TH2F ("hmassNN_missmass_vs_hmissmassNN_oa9deg","",1500,0,700,1500,600,1500);




    TH2F *hmassNP_missmass_vs_hmissmassNP_oa9deg_invm_new = new TH2F ("hmassNP_missmass_vs_hmissmassNP_oa9deg_invm_new","",1500,0,700,1500,600,1500);
    TH2F *hmassPP_missmass_vs_hmissmassPP_oa9deg_invm_new = new TH2F ("hmassPP_missmass_vs_hmissmassPP_oa9deg_invm_new","",1500,0,700,1500,600,1500);
    TH2F *hmassNN_missmass_vs_hmissmassNN_oa9deg_invm_new = new TH2F ("hmassNN_missmass_vs_hmissmassNN_oa9deg_invm_new","",1500,0,700,1500,600,1500);



    TH2F *hmassNP_missmass_vs_hmissmassNP_oa9deg_cut_neutron_peak = new TH2F ("hmassNP_missmass_vs_hmissmassNP_oa9deg_cut_neutron_peak","",1500,0,700,1500,600,1500);
    TH2F *hmassPP_missmass_vs_hmissmassPP_oa9deg_cut_neutron_peak = new TH2F ("hmassPP_missmass_vs_hmissmassPP_oa9deg_cut_neutron_peak","",1500,0,700,1500,600,1500);
    TH2F *hmassNN_missmass_vs_hmissmassNN_oa9deg_cut_neutron_peak = new TH2F ("hmassNN_missmass_vs_hmissmassNN_oa9deg_cut_neutron_peak","",1500,0,700,1500,600,1500);



  
    TH1F* hsec_single_pos = new TH1F("hsec_single_pos","",6,0,6);
    TH1F* hsec_single_ele = new TH1F("hsec_single_ele","",6,0,6);
    TH1F* hsec_pair_pos = new TH1F("hsec_pair_pos","",6,0,6);
    TH1F* hsec_pair_ele = new TH1F("hsec_pair_ele","",6,0,6);
    TH1F* hsec_pair_pos_op9 = new TH1F("hsec_pair_pos_op9","",6,0,6);
    TH1F* hsec_pair_ele_op9 = new TH1F("hsec_pair_ele_op9","",6,0,6);

    TH1F* hringPadNr[6];
    TH1F* hringPadNr_15_30[6];
    TH1F* hringPadNr_30_50[6];
    TH1F* hringPadNr_50_70[6];
    TH1F* hringPadNr_70_90[6];
    TH1F* hringClusterNr[6];
    TH1F* hringClusterNr_15_30[6];
    TH1F* hringClusterNr_30_50[6];
    TH1F* hringClusterNr_50_70[6];
    TH1F* hringClusterNr_70_90[6];
    TH1F* hringAC[6];
    TH1F* hringAC_15_30[6];
    TH1F* hringAC_30_50[6];
    TH1F* hringAC_50_70[6];
    TH1F* hringAC_70_90[6];

    for(Int_t j = 0; j < 6; j++) {
      hringPadNr[j]= new TH1F(Form("hringPadNr_sec%i",j),"",50,0,50);
      hringPadNr_15_30[j]= new TH1F(Form("hringPadNr_15_30_sec%i",j),"",50,0,50);
      hringPadNr_30_50[j]= new TH1F(Form("hringPadNr_30_50_sec%i",j),"",50,0,50);
      hringPadNr_50_70[j]= new TH1F(Form("hringPadNr_50_70_sec%i",j),"",50,0,50);
      hringPadNr_70_90[j]= new TH1F(Form("hringPadNr_70_90_sec%i",j),"",50,0,50);
      hringClusterNr[j]= new TH1F(Form("hringClusterNr_sec%i",j),"",30,0,30);
      hringClusterNr_15_30[j]= new TH1F(Form("hringClusterNr_15_30_sec%i",j),"",30,0,30);
      hringClusterNr_30_50[j]= new TH1F(Form("hringClusterNr_30_50_sec%i",j),"",30,0,30);
      hringClusterNr_50_70[j]= new TH1F(Form("hringClusterNr_50_70_sec%i",j),"",30,0,30);
      hringClusterNr_70_90[j]= new TH1F(Form("hringClusterNr_70_90_sec%i",j),"",30,0,30);
      hringAC[j]= new TH1F(Form("hringAC_sec%i",j),"",400,0,400);
      hringAC_15_30[j]= new TH1F(Form("hringAC_15_30_sec%i",j),"",400,0,400);
      hringAC_30_50[j]= new TH1F(Form("hringAC_30_50_sec%i",j),"",400,0,400);
      hringAC_50_70[j]= new TH1F(Form("hringAC_50_70_sec%i",j),"",400,0,400);
      hringAC_70_90[j]= new TH1F(Form("hringAC_70_90_sec%i",j),"",400,0,400);
    }

    TH1F* hringPadNr_single[6];
    TH1F* hringPadNr_15_30_single[6];
    TH1F* hringPadNr_30_50_single[6];
    TH1F* hringPadNr_50_70_single[6];
    TH1F* hringPadNr_70_90_single[6];
    TH1F* hringClusterNr_single[6];
    TH1F* hringClusterNr_15_30_single[6];
    TH1F* hringClusterNr_30_50_single[6];
    TH1F* hringClusterNr_50_70_single[6];
    TH1F* hringClusterNr_70_90_single[6];
    TH1F* hringAC_single[6];
    TH1F* hringAC_15_30_single[6];
    TH1F* hringAC_30_50_single[6];
    TH1F* hringAC_50_70_single[6];
    TH1F* hringAC_70_90_single[6];
    
    TH1F* hrich_sector_single = new TH1F("hrich_sector_single","",6,0,6);
    TH1F* hrich_sector_single_pion = new TH1F("hrich_sector_single_pion","",6,0,6);
    
    for(Int_t j = 0; j < 6; j++) {
      hringPadNr_single[j]= new TH1F(Form("hringPadNr_single_sec%i",j),"",50,0,50);
      hringPadNr_15_30_single[j]= new TH1F(Form("hringPadNr_15_30_single_sec%i",j),"",50,0,50);
      hringPadNr_30_50_single[j]= new TH1F(Form("hringPadNr_30_50_single_sec%i",j),"",50,0,50);
      hringPadNr_50_70_single[j]= new TH1F(Form("hringPadNr_50_70_single_sec%i",j),"",50,0,50);
      hringPadNr_70_90_single[j]= new TH1F(Form("hringPadNr_70_90_single_sec%i",j),"",50,0,50);
      hringClusterNr_single[j]= new TH1F(Form("hringClusterNr_single_sec%i",j),"",30,0,30);
      hringClusterNr_15_30_single[j]= new TH1F(Form("hringClusterNr_15_30_single_sec%i",j),"",30,0,30);
      hringClusterNr_30_50_single[j]= new TH1F(Form("hringClusterNr_30_50_single_sec%i",j),"",30,0,30);
      hringClusterNr_50_70_single[j]= new TH1F(Form("hringClusterNr_50_70_single_sec%i",j),"",30,0,30);
      hringClusterNr_70_90_single[j]= new TH1F(Form("hringClusterNr_70_90_single_sec%i",j),"",30,0,30);
      hringAC_single[j]= new TH1F(Form("hringAC_single_sec%i",j),"",400,0,400);
      hringAC_15_30_single[j]= new TH1F(Form("hringAC_15_30_single_sec%i",j),"",400,0,400);
      hringAC_30_50_single[j]= new TH1F(Form("hringAC_30_50_single_sec%i",j),"",400,0,400);
      hringAC_50_70_single[j]= new TH1F(Form("hringAC_50_70_single_sec%i",j),"",400,0,400);
      hringAC_70_90_single[j]= new TH1F(Form("hringAC_70_90_single_sec%i",j),"",400,0,400);
    }

    TH1F* hmassNP_oa9deg_no_eff = new TH1F("hmassNP_oa9deg_no_eff","Invariant mass",60,0,600);
    TH1F* hmassPP_oa9deg_no_eff = new TH1F("hmassPP_oa9deg_no_eff","Invariant mass",60,0,600);
    TH1F* hmassNN_oa9deg_no_eff = new TH1F("hmassNN_oa9deg_no_eff","Invariant mass",60,0,600);
    // hmassReco->SetXTitle("M_{inv} [MeV/c^{2}]");
    // hmassReco->SetYTitle("dN/dM_{inv} 1/[MeV/c^{2}]");
    // hmassReco->SetMarkerSize(0.5);
    // hmassReco->SetMarkerStyle(1);


   
    TH1F* hmissmassNP_oa9deg = new TH1F("hmissmassNP_oa9deg" ,"",160,600,1400);
    TH1F* hmissmassPP_oa9deg = new TH1F("hmissmassPP_oa9deg" ,"",160,600,1400);
    TH1F* hmissmassNN_oa9deg = new TH1F("hmissmassNN_oa9deg" ,"",160,600,1400);
    
    TH1F* hmassNP_missmass = new TH1F("hmassNP_missmass" ,"",60,0,600);
    TH1F* hmassPP_missmass = new TH1F("hmassPP_missmass" ,"",60,0,600);
    TH1F* hmassNN_missmass = new TH1F("hmassNN_missmass" ,"",60,0,600);

    
    TH1F* hmissmassNP_Reco_pi0 = new TH1F("hmissmassNP_Reco_pi0" ,"",160,600,1400);
    TH1F* hcos_theta = new TH1F("hcos_theta","",100,-1,1);

    TH1F* hmom_pos_reco = new TH1F("hmom_pos_reco","",70,0,700);
    TH1F* hmom_ele_reco = new TH1F("hmom_ele_reco","",70,0,700);

    TH1F* hmom_pos_reco_PP = new TH1F("hmom_pos_reco_PP","",70,0,700);
    TH1F* hmom_ele_reco_NN = new TH1F("hmom_ele_reco_NN","",70,0,700);

    TH1F *htheta_pos_reco = new TH1F("htheta_pos_reco","",40,10,90);
    TH1F *htheta_ele_reco = new TH1F("htheta_ele_reco","",40,10,90);

    TH1F *htheta_pos_reco_PP = new TH1F("htheta_pos_reco_PP","",40,10,90);
    TH1F *htheta_ele_reco_NN = new TH1F("htheta_ele_reco_NN","",40,10,90);
    
    TH1F* hoangle_reco = new TH1F("hoangle_reco","",180,0,180);
    TH1F* hoangle_reco_PP = new TH1F("hoangle_reco_PP","",180,0,180);
    TH1F* hoangle_reco_NN = new TH1F("hoangle_reco_NN","",180,0,180);
    TH1F* hpt = new TH1F("hpt","",350,0,700);
    TH1F* hpt_PP = new TH1F("hpt_PP","",350,0,700);
    TH1F* hpt_NN = new TH1F("hpt_NN","",350,0,700);

    TH1F* hmom1mom2 = new TH1F("mom1mom2","",1000,0,1000);
    TH1F* hmom1mom2_PP = new TH1F("mom1mom2_PP","",1000,0,1000);
    TH1F* hmom1mom2_NN = new TH1F("mom1mom2_NN","",1000,0,1000);

    TH1F* hpt_0_140 = new TH1F("hpt_0_140","",350,0,700);
    TH1F* hpt_PP_0_140 = new TH1F("hpt_PP_0_140","",350,0,700);
    TH1F* hpt_NN_0_140 = new TH1F("hpt_NN_0_140","",350,0,700);

    TH1F* hpt_140_300 = new TH1F("hpt_140_300","",350,0,700);
    TH1F* hpt_PP_140_300 = new TH1F("hpt_PP_140_300","",350,0,700);
    TH1F* hpt_NN_140_300 = new TH1F("hpt_NN_140_300","",350,0,700);

    TH1F* hpt_300_600 = new TH1F("hpt_300_600","",350,0,700);
    TH1F* hpt_PP_300_600 = new TH1F("hpt_PP_300_600","",350,0,700);
    TH1F* hpt_NN_300_600 = new TH1F("hpt_NN_300_600","",350,0,700);


     // Vertex
    TH1F* hvertexZ   = new TH1F("hvertexZ", "", 1000,-200,200);
    TH1F* hvertexX   = new TH1F("hvertexX", "", 1000,-200,200);
    TH1F* hvertexY   = new TH1F("hvertexY", "", 1000,-200,200);
    TH2F* hvertexXZ  = new TH2F("hvertexXZ","",360,-100,20,500,-10,10);
    TH2F* hvertexXY  = new TH2F("hvertexXY","",500,-10,10,500,-10,10);
    TH2F* hvertexYZ  = new TH2F("hvertexYZ","",360,-100,20,500,-10,10);
     // Vertex reco
    TH1F* hvertexZ_reco   = new TH1F("hvertexZ_reco", "", 1000,-200,200);
    TH1F* hvertexX_reco   = new TH1F("hvertexX_reco", "", 1000,-200,200);
    TH1F* hvertexY_reco   = new TH1F("hvertexY_reco", "", 1000,-200,200);
    TH2F* hvertexXZ_reco  = new TH2F("hvertexXZ_reco","",360,-100,20,500,-10,10);
    TH2F* hvertexXY_reco  = new TH2F("hvertexXY_reco","",500,-10,10,500,-10,10);
    TH2F* hvertexYZ_reco  = new TH2F("hvertexYZ_reco","",360,-1000,200,500,-10,10);




    TH1F* hnCand     = new TH1F("hnCand","",20,0,20);
    TH1F* hnCand_best     = new TH1F("hnCand_best","",20,0,20);
    TH1F* hnCand_best_excl     = new TH1F("hnCand_best_excl","",20,0,20);
    TH2F* hr_vs_z = new TH2F("hr_vs_z", "",200,-500,300, 100,-100,100);
    TH2F* hBetaMom_dilep1_op9 = new TH2F("hBetaMom_dilep1_op9","", 200,-1000,1000,500,0.7,1.4);
    TH2F* hBetaMom = new TH2F("hBetaMom","", 200,-1000,1000,500,0,1.4);



    TH2F* hBetaMom_ring = new TH2F("hBetaMom_ring","", 200,-1000,1000,500,0.7,1.4);
    TH2F* hBetaMom_ring_no_backtracking = new TH2F("hBetaMom_ring_no_backtracking","", 200,-1000,1000,500,0.7,1.4);
    TH2F* hBetaMom_backtracking = new TH2F("hBetaMom_backtracking","", 200,-1000,1000,500,0.7,1.4);
    TH2F* hBetaMom_backtracking_no_ring = new TH2F("hBetaMom_backtracking_no_ring","", 200,-1000,1000,500,0.7,1.4);
    
    TH2F* hrichQa_mom_ring = new TH2F("hrichQa_mom_ring","",400,-800,800,200,-2,18);
    TH2F* hrichQa_mom_backtracking = new TH2F("hrichQa_mom_backtracking","",400,-800,800,200,-2,18);
    TH2F* hrichQa_mom_ring_no_backtracking = new TH2F("hrichQa_mom_ring_no_backtracking","",400,-800,800,200,-2,18);
    TH2F* hrichQa_mom_ring_no_backtracking_high_beta = new TH2F("hrichQa_mom_ring_no_backtracking_high_beta","",400,-800,800,200,-2,18);
    TH2F* hrichQa_mom_backtracking_no_ring = new TH2F("hrichQa_mom_backtracking_no_ring","",400,-800,800,200,-2,18);
    
    //-----------------------------------------------------------------------------------------//**********************************************    
    Int_t nbytes = 0;
    Int_t entries = loop->getEntries();  // all entries of all files in chain
    Int_t nSecondary=0;
    Bool_t printHistoryFlag=kTRUE;
      
    TStopwatch timer;
    timer.Reset();
    timer.Start();

    TLorentzVector vpi(0,0,685,699);
    TLorentzVector vnu(0,0,0,938);
    TLorentzVector in = vpi+ vnu;

    HEnergyLossCorrPar energylosspar;
    //energylosspar.setDefaultPar("aug14_PE");
    energylosspar.setDefaultPar("aug14_pe");


    
    HParticleTrackSorter sorter;
    sorter.setUseYMatching(kTRUE);
    //sorter.setBetaLeptonCut(0.5);
    //sorter.setRICHMatching(HParticleTrackSorter::kUseRKRICHWindow,8.); // select matching RICH-MDC for selectLeptons() function
    sorter.setMETAQACut(4.);
    //sorter.setDebug();                                            // for debug
    //sorter.setPrintLevel(3);                                      // max prints
    //sorter.setRICHMatching(HParticleTrackSorter::kUseRKRICHWindow,8.); // select matching RICH-MDC for selectLeptons() function
    //sorter.setIgnoreInnerMDC();                                   // do not reject Double_t inner MDC hits
    //sorter.setIgnoreOuterMDC();                                   // do not reject Double_t outer MDC hits
    //sorter.setIgnoreMETA();                                       // do not reject Double_t META hits
    //sorter.setIgnorePreviousIndex();                              // do not reject indices from previous selctions
    sorter.init();                                                  // get catgegory pointers etc...

    Float_t oAngleCut      = 9.;
    //Float_t momCut         = 80.;
    Float_t momCut         = 100.;


    Float_t nonFittedCut4  = 6;
    Float_t fittedCut4     = 6;

    Float_t Chi2RkCut      = 500.;
    Float_t Chi2InCut      = 500.;
    Float_t Chi2OutCut     = 500.;
    
    Int_t counter_lep_event=0;
    
    
    //-------------------------------------------------------------------------------------------------
    // Event LOOP begin
    for(Int_t i = 1; i <= nev; i++)
    {
	//----------break if last event is reached-------------
	if(loop->nextEvent(i) <= 0) { cout<<" end recieved "<<endl; break; } // last event reached
	HTool::printProgress(i,nev,1,"Analysing evt# :");
	hCounter->Fill(0);

        evtinfo = HCategoryManager::getObject(evtinfo,catParticleEvtInfo,0);

        if(!evtinfo->isGoodEvent(1)) continue;
	
	sorter.cleanUp();
	//------------------------------------------------------------------------

	sorter.resetFlags(kTRUE,kTRUE,kTRUE,kTRUE);            // reset all flags for flags (0-28) ,reject,used,lepton
	// Int_t nCandLep     = sorter.fill(selectLeptonsBeta);   // fill only good leptons
	//Int_t nCandLep       = sorter.fill(HParticleTrackSorter::selectLeptons);
	//Int_t nCandLepBest = sorter.selectBest(HParticleTrackSorter::kIsBestRKRKMETA,HParticleTrackSorter::kIsLepton);
	Int_t nCandHad     = sorter.fill(HParticleTrackSorter::selectHadrons);   // fill only good hadrons (already marked good leptons will be skipped)
	Int_t nCandHadBest = sorter.selectBest(HParticleTrackSorter::kIsBestRKRKMETA,HParticleTrackSorter::kIsHadron);


	// Get vertex from combined fit of fitted inner segments



	Double_t vertexX = gHades->getCurrentEvent()->getHeader()->getVertex().getX();
        Double_t vertexY = gHades->getCurrentEvent()->getHeader()->getVertex().getY();
        Double_t vertexZ = gHades->getCurrentEvent()->getHeader()->getVertex().getZ();

	Double_t vertexX_reco = gHades->getCurrentEvent()->getHeader()->getVertexReco().getX();
        Double_t vertexY_reco = gHades->getCurrentEvent()->getHeader()->getVertexReco().getY();
        Double_t vertexZ_reco = gHades->getCurrentEvent()->getHeader()->getVertexReco().getZ();
	

	//if(gHades->getCurrentEvent()->getHeader()->getVertexReco().getChi2()==-1) continue;
	
	hvertexXZ->Fill(vertexZ,vertexX);
	hvertexXY->Fill(vertexY,vertexX);
	hvertexYZ->Fill(vertexZ,vertexY);
	hvertexX->Fill(vertexX);
	hvertexY->Fill(vertexY);
	hvertexZ->Fill(vertexZ);



	hvertexXZ_reco->Fill(vertexZ_reco,vertexX_reco);
	hvertexXY_reco->Fill(vertexY_reco,vertexX_reco);
	hvertexYZ_reco->Fill(vertexZ_reco,vertexY_reco);
	hvertexX_reco->Fill(vertexX_reco);
	hvertexY_reco->Fill(vertexY_reco);
	hvertexZ_reco->Fill(vertexZ_reco);







	
      	//----------------looping data-------------------------
	Int_t index = 0;
	Int_t index1 = 0;
        Int_t index2 = 0;
	Int_t size = candCat->getEntries();

	Bool_t leptonFlag[size];

	for(Int_t j = 0; j < size; j ++){
	  leptonFlag[j]  = kFALSE;	  
	  cand1 = HCategoryManager::getObject(cand1,candCat,j);
	  
	  cand1->calc4vectorProperties();

	  Int_t sys1          =cand1->getSystem();
	  Float_t richQa1     = cand1->getRichMatchingQuality();
	  Float_t metaQa1     = cand1->getMetaMatchQuality();
	  Int_t charge1       = cand1->getCharge();
	  //Int_t pid          = cand1->getGeantPID();
	  Int_t sec1          =  cand1->getSector();
	  
	  Float_t chi2RK1      = cand1->getChi2();
	  Float_t chi2In1      = cand1->getInnerSegmentChi2();
	  Float_t chi2Out1     = cand1->getOuterSegmentChi2();
	  
	  Float_t beta1       = cand1->getBeta();
	  Float_t mom1       = cand1->getMomentum();
	  Float_t theta1      = cand1->getTheta();
	  Float_t phi1        = cand1->getPhi();
	  
	  Float_t z1           = cand1->getZ();
	  Float_t rad1         = cand1->getR();

	  Int_t richind_1       = cand1->getRichInd();

	  Float_t momCor = energylosspar.getCorrMom(3,mom1,theta1);

	  mom1=momCor;


	  hBetaMom->Fill(mom1*charge1,beta1);

          //hVertexZ->Fill(z1);


	  //if(beta1<0.80) continue;
	  if(z1 <-100 || z1> 50) continue;
	  if(TMath::Abs(rad1)>40) continue;
	  if(cand1->isAtAnyMdcEdge()==kTRUE) continue;
	  if(mom1<momCut) continue;

	  //if(richQa1>4 && charge1==-1) continue;
          //RICH
	  //if(richQa1>4 ) continue;
	  //if(mom1>300) continue;
	  if(chi2RK1>Chi2RkCut || chi2In1>Chi2InCut || chi2Out1>Chi2OutCut) continue;
          //RICH
	  //if(!cand1->isFlagBit(Particle::kIsAcceptedHitRICH) && charge1==-1)continue;
	  //if(!cand1->isFlagBit(Particle::kIsAcceptedHitRICH) )continue;
          //
	  if(!cand1->isFlagBit(Particle::kIsAcceptedHitInnerMDC)) continue;
	  if(!cand1->isFlagBit(Particle::kIsAcceptedHitOuterMDC)) continue;
	  if(!cand1->isFlagBit(Particle::kIsAcceptedHitMETA))  continue;
	  if(!cand1->isFlagBit(Particle::kIsAcceptedRK)) continue;

          //BACK

	  if(cand1->getRichBTInd()==-1) continue;
	  btRing = HCategoryManager::getObject(btRing,btCat,cand1->getRichBTInd());
	  if(btRing->getMaxima()<1)  continue;

	  hr_vs_z->Fill(z1,rad1);
	  leptonFlag[j]=kTRUE;
	  //if(cand1->getGeantParentTrackNum()!=0) continue;
	  
	  if(charge1 ==1) {
	    htheta_phi_raw_single_pos->Fill(phi1,theta1);
	    htheta_phi_raw_single_pos_rebin->Fill(phi1,theta1);
	  }
	  else if(charge1 ==-1){
	    htheta_phi_raw_single_ele->Fill(phi1,theta1);
	    htheta_phi_raw_single_ele_rebin->Fill(phi1,theta1);
	  }

	  //if(pid ==9 || pid ==8 ) hrich_sector_single_pion->Fill(sec1);
	  
	  if(richind_1 !=-1 ){
	    HRichHit* richhit_1 = HParticleTool::getRichHit(richind_1);
	    
	    Int_t   ringSector_1   = richhit_1->getSector();
	    Float_t ringPadNr_1    = (float)richhit_1->getRingPadNr();
	   
	    Int_t   cluster_num_1  = richhit_1->getRingClusterNr();
	    
	    Float_t ringAmpl_1     = (float)richhit_1->getRingAmplitude();
	    
	    Float_t ringAC_1       = (float)ringAmpl_1/ringPadNr_1;
	    
	    hrich_sector_single->Fill(ringSector_1);
	    hringPadNr_single[ringSector_1]->Fill(ringPadNr_1);
	    hringClusterNr_single[ringSector_1]->Fill(cluster_num_1);
	    hringAC_single[ringSector_1]->Fill(ringAC_1);
	    if(theta1> 15 && theta1<30) {
	      hringPadNr_15_30_single[ringSector_1]->Fill(ringPadNr_1);
	      hringClusterNr_15_30_single[ringSector_1]->Fill(cluster_num_1);
	      hringAC_15_30_single[ringSector_1]->Fill(ringAC_1);
	    }
	    else if(theta1> 30 && theta1<50) {
	      hringPadNr_30_50_single[ringSector_1]->Fill(ringPadNr_1);
	      hringClusterNr_30_50_single[ringSector_1]->Fill(cluster_num_1);
	      hringAC_30_50_single[ringSector_1]->Fill(ringAC_1);
	    }
	    else if(theta1> 50 && theta1<70) {
	      hringPadNr_50_70_single[ringSector_1]->Fill(ringPadNr_1);
	      hringClusterNr_50_70_single[ringSector_1]->Fill(cluster_num_1);
	      hringAC_50_70_single[ringSector_1]->Fill(ringAC_1);
	    }
	    else if(theta1> 70 && theta1<90) {
	      hringPadNr_70_90_single[ringSector_1]->Fill(ringPadNr_1);
	      hringClusterNr_70_90_single[ringSector_1]->Fill(cluster_num_1);
	      hringAC_70_90_single[ringSector_1]->Fill(ringAC_1);
	    }
	  }
	}  
	  
	//------------------------------------------------------------------
        // find good pairs surviving oA cut
	vector<Int_t> indexOAfailed;
	for(Int_t j = 0; j < size-1; j ++){
	    cand1 = HCategoryManager::getObject(cand1,candCat,j);
	    cand1->calc4vectorProperties();
	    //if(!cand1->isFlagBit(Particle::kIsLepton)) continue;
	    index1 = cand1 ->getIndex();
	    if(leptonFlag[j]==kFALSE) continue;

	    for(Int_t k = j+1; k < size; k ++){
		cand2 = HCategoryManager::getObject(cand2,candCat,k);
		cand2->calc4vectorProperties();
		//if(!cand2->isFlagBit(Particle::kIsLepton)) continue;
		index2 = cand2 ->getIndex();
		if(leptonFlag[k]==kFALSE) continue;

		if(isGoodNonFitted(cand1,nonFittedCut4) &&
		   isGoodNonFitted(cand2,nonFittedCut4) &&
		   isGoodFitted   (cand1,fittedCut4   ) &&
		   isGoodFitted   (cand2,fittedCut4   )
		  ) {
		    if((1==cand1->getCharge() && -1==cand2->getCharge()) ||
		       (1==cand2->getCharge() && -1==cand1->getCharge())) {

			Float_t oAngle   = cand1->Angle(cand2->Vect())*TMath::RadToDeg();
			if(oAngle<=oAngleCut){
			    if(find(indexOAfailed.begin(),indexOAfailed.end(),index1)==indexOAfailed.end()){
				indexOAfailed.push_back(index1);
			    }
			    if(find(indexOAfailed.begin(),indexOAfailed.end(),index2)==indexOAfailed.end()){
				indexOAfailed.push_back(index2);
			    }
			}
		    }
		}
	    }
	}


	//================================================================
	// Pairs loop begin
	for(Int_t m = 0; m < size-1; m ++)
	  {
	    cand1 = HCategoryManager::getObject(cand1,candCat,m);
	    cand1->calc4vectorProperties();
	    
	    if(!cand1->isFlagBit(Particle::kIsUsed)) continue;
	    //if(!cand1->isFlagBit(Particle::kIsLepton)) continue;
	  
	    index1 = cand1 ->getIndex();
	    	      
	    
	    for(Int_t k = m+1; k < size; k ++)
	      {
		cand2 = HCategoryManager::getObject(cand2,candCat,k);
		cand2->calc4vectorProperties();
		
		if(!cand2->isFlagBit(Particle::kIsUsed)) continue;
		//if(!cand2->isFlagBit(Particle::kIsLepton)) continue;
						
		//if(cand1->isGhostTrack()==1 || cand2->isGhostTrack()==1)continue;

		Int_t sec1           = cand1->getSector();
		Int_t sec2           = cand2->getSector();
		Int_t charge1        = cand1->getCharge();
		Int_t charge2        = cand2->getCharge();
		Float_t beta1        = cand1->getBeta();
		Float_t beta2        = cand2->getBeta();
		Float_t mom1         = cand1->getMomentum();
		Float_t mom2         = cand2->getMomentum();
		Float_t richQa1      = cand1->getRichMatchingQuality();
		Float_t richQa2      = cand2->getRichMatchingQuality();
		Float_t z1           = cand1->getZ();
		Float_t rad1         = cand1->getR();
		Float_t z2           = cand2->getZ();
		Float_t rad2         = cand2->getR();
		Float_t theta1       = cand1->getTheta();
		Float_t phi1         = cand1->getPhi();
		Float_t theta2       = cand2->getTheta();
		Float_t phi2         = cand2->getPhi();
		Int_t richind_1       = cand1->getRichInd();
		Int_t richind_2       = cand2->getRichInd();
		Float_t chi2RK1      = cand1->getChi2();
		Float_t chi2In1      = cand1->getInnerSegmentChi2();
		Float_t chi2Out1     = cand1->getOuterSegmentChi2();
		Float_t chi2RK2      = cand2->getChi2();
		Float_t chi2In2      = cand2->getInnerSegmentChi2();
		Float_t chi2Out2     = cand2->getOuterSegmentChi2();
		Float_t metaQa1      = cand1->getMetaMatchQuality();
		Float_t metaQa2      = cand2->getMetaMatchQuality();

		Float_t momCor = energylosspar.getCorrMom(3,mom1,theta1);
		Float_t momCor2 = energylosspar.getCorrMom(3,mom2,theta2);

		mom1=momCor;
		mom2=momCor2;

		cand1->setMomentum(mom1);
		cand2->setMomentum(mom2);
			









		//if(beta1<0.80) continue;
		//if(beta2<0.80) continue;
		if(z1 <-100 || z1> 50) continue;
		if(z2 <-100 || z2> 50) continue;
		if(TMath::Abs(rad1)>40) continue;
		if(TMath::Abs(rad2)>40) continue;
		if(cand1->isAtAnyMdcEdge()==kTRUE) continue;
		if(cand2->isAtAnyMdcEdge()==kTRUE) continue;
                if(mom1<momCut) continue;
		if(mom2<momCut) continue;

		//if(richQa1>8 ) continue;
	        //if(richQa2>8 ) continue;
	        // if(richQa1>8 && charge1==-1) continue;
		// if(richQa2>8 && charge2==-1) continue;

	        //if(richQa1>4 ) continue;
		//if(richQa2>4 ) continue;

	        //if(richQa1>4 && charge1==-1) continue;
	        //if(richQa2>4 && charge2==-1) continue;
		//if(mom1>300) continue;
		//if(charge1==-1 && mom1>250) continue;
		//if(mom2>300) continue;
		//if(charge2==-1 && mom2>250) continue;
		//if(cand1->getRichBTInd()==-1 && charge1==-1) continue;
		//if(cand2->getRichBTInd()==-1 && charge2==-1) continue;
	        
		// if(metaQa1 > 4) continue;
		// if(metaQa2 > 4) continue;
		if(chi2RK1>Chi2RkCut || chi2In1>Chi2InCut || chi2Out1>Chi2OutCut) continue;
		if(chi2RK2>Chi2RkCut || chi2In2>Chi2InCut || chi2Out2>Chi2OutCut) continue;
                //RICH
		//if(!cand1->isFlagBit(Particle::kIsAcceptedHitRICH) && charge1==-1)continue;
	       // if(!cand1->isFlagBit(Particle::kIsAcceptedHitRICH) )continue;
                //
		if(!cand1->isFlagBit(Particle::kIsAcceptedHitInnerMDC)) continue;
		if(!cand1->isFlagBit(Particle::kIsAcceptedHitOuterMDC)) continue;
		if(!cand1->isFlagBit(Particle::kIsAcceptedHitMETA))  continue;
		if(!cand1->isFlagBit(Particle::kIsAcceptedRK)) continue;
		//RICH
		//if(!cand2->isFlagBit(Particle::kIsAcceptedHitRICH) && charge2==-1 )continue;
		//if(!cand2->isFlagBit(Particle::kIsAcceptedHitRICH) )continue;
                //
		if(!cand2->isFlagBit(Particle::kIsAcceptedHitInnerMDC)) continue;
		if(!cand2->isFlagBit(Particle::kIsAcceptedHitOuterMDC)) continue;
		if(!cand2->isFlagBit(Particle::kIsAcceptedHitMETA))  continue;
		if(!cand2->isFlagBit(Particle::kIsAcceptedRK)) continue;

		if(cand1->getRichBTInd()==-1 && cand1->isFlagBit(Particle::kIsAcceptedHitRICH))  {
		  hrichQa_mom_ring_no_backtracking->Fill(mom1*charge1,richQa1);
		  if(beta1>0.93)hrichQa_mom_ring_no_backtracking_high_beta->Fill(mom1*charge1,richQa1);
		  hBetaMom_ring_no_backtracking->Fill(mom1*charge1,beta1);
		  
		}
		if(cand2->getRichBTInd()==-1 && cand2->isFlagBit(Particle::kIsAcceptedHitRICH)) {
		  hrichQa_mom_ring_no_backtracking->Fill(mom2*charge2,richQa2);
		  if(beta2>0.93)hrichQa_mom_ring_no_backtracking_high_beta->Fill(mom2*charge2,richQa2);
		  hBetaMom_ring_no_backtracking->Fill(mom2*charge2,beta2);
		}
		if(cand1->isFlagBit(Particle::kIsAcceptedHitRICH))  {
		  hrichQa_mom_ring->Fill(mom1*charge1,richQa1);
		  hBetaMom_ring->Fill(mom1*charge1,beta1);
		}
		if(cand2->isFlagBit(Particle::kIsAcceptedHitRICH))  {
		  hrichQa_mom_ring->Fill(mom2*charge2,richQa2);
		  hBetaMom_ring->Fill(mom1*charge1,beta1);
		}


		//BACK

		if(cand1->getRichBTInd()==-1 ) continue;
	        if(cand2->getRichBTInd()==-1 ) continue;
		btRing = HCategoryManager::getObject(btRing,btCat,cand1->getRichBTInd());
		if(btRing->getMaxima()<1) continue;
		//}

                  //if(cand2->getRichBTInd()!=-1 && charge2==-1){
		btRing2 = HCategoryManager::getObject(btRing2,btCat,cand2->getRichBTInd());
		if(btRing2->getMaxima()<1) continue;

               //

		hrichQa_mom_backtracking->Fill(mom1*charge1,richQa1);
		hrichQa_mom_backtracking->Fill(mom2*charge2,richQa2);
		hBetaMom_backtracking->Fill(mom1*charge1,beta1);
		hBetaMom_backtracking->Fill(mom2*charge2,beta2);
		if( !cand1->isFlagBit(Particle::kIsAcceptedHitRICH))  {
		  hrichQa_mom_backtracking_no_ring->Fill(mom1*charge1,richQa1);
		  hBetaMom_backtracking_no_ring->Fill(mom1*charge1,beta1);
		}
		if(!cand2->isFlagBit(Particle::kIsAcceptedHitRICH)) {
		  hrichQa_mom_backtracking_no_ring->Fill(mom2*charge2,richQa2);
		  hBetaMom_backtracking_no_ring->Fill(mom2*charge2,beta2);
		}
		
		// if(((cand1->getGeantPID()==2 && cand2->getGeantPID()==3 ) ||
		// 	(cand1->getGeantPID()==3 && cand2->getGeantPID()==2 )) &&
		//        (cand1->getGeantParentTrackNum()==0 && cand2->getGeantParentTrackNum()==0 ) )
		//    {
		//      if(cand1->isFlagBit(Particle::kIsAcceptedHitRICH) && charge1==1 )hsec_pair_pos->Fill(sec1);
		//      else if(cand1->isFlagBit(Particle::kIsAcceptedHitRICH) && charge1==-1 )hsec_pair_ele->Fill(sec1);
		     
		//      if(cand2->isFlagBit(Particle::kIsAcceptedHitRICH) && charge2==1 )hsec_pair_pos->Fill(sec2);
		//      else if(cand2->isFlagBit(Particle::kIsAcceptedHitRICH) && charge2==-1 )hsec_pair_ele->Fill(sec2);
		//    }
				
		//if(!cand1->isFlagBit(Particle::kIsAcceptedHitRICH)||!cand2->isFlagBit(Particle::kIsAcceptedHitRICH)) continue;
		index2 = cand2 ->getIndex();
		
		// Pair observables
		TLorentzVector pair = (*cand1) + (*cand2);
					
		Float_t invM     = pair.M();
		Float_t oAngle   = cand1->Angle(cand2->Vect())*TMath::RadToDeg();
                Float_t angle_pair = pair.Theta();
		Float_t pt       = pair.Pt();
				
		Bool_t goodNonFitted4 = kFALSE;
		Bool_t goodFitted4    = kFALSE;
		
		if(isGoodNonFitted(cand1,nonFittedCut4) &&
		   isGoodNonFitted(cand2,nonFittedCut4)) goodNonFitted4 = kTRUE;
		if(isGoodFitted   (cand1,fittedCut4   ) &&
		   isGoodFitted   (cand2,fittedCut4   )) goodFitted4    = kTRUE;
		
		if(oAngle>9 )
		{
		    //Select true e+e- pairs
		    if(goodNonFitted4 && goodFitted4)
		    {
			//----------------------------------------------------------------------
			// NP  recursive remove cand from pairs failing oA cut
		        if(find(indexOAfailed.begin(),indexOAfailed.end(),index1)==indexOAfailed.end() &&
		           find(indexOAfailed.begin(),indexOAfailed.end(),index2)==indexOAfailed.end()
			  ){
			    if(((1==cand1->getCharge() && -1==cand2->getCharge()) ||
				(1==cand2->getCharge() && -1==cand1->getCharge())) )
			    {
				  hmassNP_oa9deg_no_eff->Fill(invM);
				  hBetaMom_dilep1_op9->Fill(mom1*charge1,beta1);
				  hBetaMom_dilep1_op9->Fill(mom2*charge2,beta2);
				  hpt->Fill(pt);

				  if(invM<140){
				      hoangle_reco->Fill(oAngle);
				      hpt_0_140->Fill(pt);
				      hmom1mom2->Fill(TMath::Sqrt(mom1*mom2));
				      hcos_theta->Fill(TMath::Cos(angle_pair));
				      TLorentzVector sub = in - pair;
				      Float_t missmass = sub.M();
				      hmissmassNP_Reco_pi0->Fill(missmass);

				      if(1==cand1->getCharge())  {
					  hmom_pos_reco->Fill(mom1);
					  htheta_phi_raw_pair_pos->Fill(phi1,theta1);
					  htheta_phi_raw_pair_pos_rebin->Fill(phi1,theta1);
				      }
				      else if(-1==cand1->getCharge())  	{
					  hmom_ele_reco->Fill(mom1);
					  htheta_phi_raw_pair_ele->Fill(phi1,theta1);
					  htheta_phi_raw_pair_ele_rebin->Fill(phi1,theta1);
				      }
				      if(1==cand2->getCharge())  {
					  hmom_pos_reco->Fill(mom2);
					  htheta_phi_raw_pair_pos->Fill(phi2,theta2);
					  htheta_phi_raw_pair_pos_rebin->Fill(phi2,theta2);
				      }
				      else if(-1==cand2->getCharge())  {
					  hmom_ele_reco->Fill(mom2);
					  htheta_phi_raw_pair_ele->Fill(phi2,theta2);
					  htheta_phi_raw_pair_ele_rebin->Fill(phi2,theta2);
				      }

				  }
                                  TLorentzVector sub = in - pair;
				Float_t missmass = sub.M();
				hmassNP_missmass_vs_hmissmassNP_oa9deg_invm_new->Fill(invM,missmass);



				  if(invM>140 && invM <300) hpt_140_300->Fill(pt);
				  else if(invM>300) hpt_300_600->Fill(pt);

				  if(invM>140){
				      TLorentzVector sub = in - pair;
				      Float_t missmass = sub.M();
				      hmissmassNP_oa9deg->Fill(missmass);
				    
				  }
				  if((missmass >900 &&  missmass < 1030)){
				      hmassNP_missmass->Fill(invM);
				      if(invM>300) hnCand_best_excl->Fill(nCandHadBest);
				  }

				  if(richind_1 !=-1 && richind_2 !=-1){
				      HRichHit* richhit_1 = HParticleTool::getRichHit(richind_1);
				      HRichHit* richhit_2 = HParticleTool::getRichHit(richind_2);
				      Int_t   ringSector_1   = richhit_1->getSector();
				      Float_t ringPadNr_1    = (float)richhit_1->getRingPadNr();
				      Int_t   ringSector_2   = richhit_2->getSector();
				      Float_t ringPadNr_2    = (float)richhit_2->getRingPadNr();
				      Int_t   cluster_num_1  = richhit_1->getRingClusterNr();
				      Int_t   cluster_num_2  = richhit_2->getRingClusterNr();
				      Float_t ringAmpl_1     = (float)richhit_1->getRingAmplitude();
				      Float_t ringAmpl_2     = (float)richhit_2->getRingAmplitude();
				      Float_t ringAC_1       = (float)ringAmpl_1/ringPadNr_1;
				      Float_t ringAC_2       = (float)ringAmpl_2/ringPadNr_2;

				      if(invM<150){
					  hringPadNr[ringSector_1]->Fill(ringPadNr_1);
					  hringPadNr[ringSector_2]->Fill(ringPadNr_2);
					  hringClusterNr[ringSector_1]->Fill(cluster_num_1);
					  hringClusterNr[ringSector_2]->Fill(cluster_num_2);
					  hringAC[ringSector_1]->Fill(ringAC_1);
					  hringAC[ringSector_2]->Fill(ringAC_2);
					  if(theta1> 15 && theta1<30) {
					      hringPadNr_15_30[ringSector_1]->Fill(ringPadNr_1);
					      hringClusterNr_15_30[ringSector_1]->Fill(cluster_num_1);
					      hringAC_15_30[ringSector_1]->Fill(ringAC_1);
					  }
					  else if(theta1> 30 && theta1<50) {
					      hringPadNr_30_50[ringSector_1]->Fill(ringPadNr_1);
					      hringClusterNr_30_50[ringSector_1]->Fill(cluster_num_1);
					      hringAC_30_50[ringSector_1]->Fill(ringAC_1);
					  }
					  else if(theta1> 50 && theta1<70) {
					      hringPadNr_50_70[ringSector_1]->Fill(ringPadNr_1);
					      hringClusterNr_50_70[ringSector_1]->Fill(cluster_num_1);
					      hringAC_50_70[ringSector_1]->Fill(ringAC_1);
					  }
					  else if(theta1> 70 && theta1<90) {
					      hringPadNr_70_90[ringSector_1]->Fill(ringPadNr_1);
					      hringClusterNr_70_90[ringSector_1]->Fill(cluster_num_1);
					      hringAC_70_90[ringSector_1]->Fill(ringAC_1);
					  }
					  if(theta2> 15 && theta2<30) {
					      hringPadNr_15_30[ringSector_2]->Fill(ringPadNr_2);
					      hringClusterNr_15_30[ringSector_2]->Fill(cluster_num_2);
					      hringAC_15_30[ringSector_2]->Fill(ringAC_2);
					  }
					  else if(theta2> 30 && theta2<50) {
					      hringPadNr_30_50[ringSector_2]->Fill(ringPadNr_2);
					      hringClusterNr_30_50[ringSector_2]->Fill(cluster_num_2);
					      hringAC_30_50[ringSector_2]->Fill(ringAC_2);
					  }
					  else if(theta2> 50 && theta2<70) {
					      hringPadNr_50_70[ringSector_2]->Fill(ringPadNr_2);
					      hringClusterNr_50_70[ringSector_2]->Fill(cluster_num_2);
					      hringAC_50_70[ringSector_2]->Fill(ringAC_2);
					  }
					  else if(theta2> 70 && theta2<90) {
					      hringPadNr_70_90[ringSector_2]->Fill(ringPadNr_2);
					      hringClusterNr_50_70[ringSector_2]->Fill(cluster_num_2);
					      hringAC_70_90[ringSector_2]->Fill(ringAC_2);
					  }
				      }
				  }
			      }
			      if(1==cand1->getCharge() && 1==cand2->getCharge()){
				  hmassPP_oa9deg_no_eff->Fill(invM);
				  hoangle_reco_PP->Fill(oAngle);
				  hpt_PP->Fill(pt);
				  hmom1mom2_PP->Fill(TMath::Sqrt(mom1*mom2));
				  hmom_pos_reco_PP->Fill(mom1);
				  hmom_pos_reco_PP->Fill(mom2);
				  htheta_pos_reco_PP->Fill(theta1);
				  htheta_pos_reco_PP->Fill(theta2);
				  if(invM<140)hpt_PP_0_140->Fill(pt);
				  else if(invM>140 && invM <300) hpt_PP_140_300->Fill(pt);
				  else if(invM>300) hpt_PP_300_600->Fill(pt);
				  TLorentzVector sub = in - pair;
				  Float_t missmass = sub.M();
				  if(invM>140){
				      hmissmassPP_oa9deg->Fill(missmass);
				  }
				  if((missmass >900 &&  missmass < 1030)){
				      hmassPP_missmass->Fill(invM);
				  }
				  hmassPP_missmass_vs_hmissmassPP_oa9deg_invm_new->Fill(invM,missmass);


			      }

			      if(-1==cand1->getCharge() && -1==cand2->getCharge() ) {
				  hmassNN_oa9deg_no_eff->Fill(invM);
				  hoangle_reco_NN->Fill(oAngle);
				  hpt_NN->Fill(pt);
				  hmom1mom2_NN->Fill(TMath::Sqrt(mom1*mom2));
				  hmom_ele_reco_NN->Fill(mom1);
				  hmom_ele_reco_NN->Fill(mom2);
				  htheta_ele_reco_NN->Fill(theta1);
				  htheta_ele_reco_NN->Fill(theta2);
				  if(invM<140)hpt_NN_0_140->Fill(pt);
				  else if(invM>140 && invM <300) hpt_NN_140_300->Fill(pt);
				  else if(invM>300) hpt_NN_300_600->Fill(pt);
				  TLorentzVector sub = in - pair;
				  Float_t missmass = sub.M();
				  if(invM>140){
				      hmissmassNN_oa9deg->Fill(missmass);
				  }
				  if((missmass >900 &&  missmass < 1030)){
				      hmassNN_missmass->Fill(invM);
				  }

				  hmassNN_missmass_vs_hmissmassNN_oa9deg_invm_new->Fill(invM,missmass);

			      }




			      //}

			      if(counter_lep_event!=i) hCounter->Fill(3);
			      counter_lep_event==i;
					}                                      ////opening angle
			      }
			  }        //                                //fit/unfit cut
		      } // inner loop
		  } // outer loop

		hnCand->Fill(nCandHad);
		hnCand_best->Fill(nCandHadBest);
    } // end event loop
   
    sorter.finalize();
    timer.Stop();
   
    
    //-------------------------------------------------
   
    TFile* out = new TFile(outfile.Data(),"RECREATE");

    out->cd();




    hmassNP_missmass_vs_hmissmassNP_oa9deg_invm_new->Write();
    hmassNN_missmass_vs_hmissmassNN_oa9deg_invm_new->Write();
    hmassPP_missmass_vs_hmissmassPP_oa9deg_invm_new->Write();


    hBetaMom_ring->Write(); 
    hBetaMom_ring_no_backtracking->Write(); 
    hBetaMom_backtracking->Write(); 
    hBetaMom_backtracking_no_ring->Write(); 
    
    hrichQa_mom_ring->Write();
    hrichQa_mom_ring_no_backtracking->Write();
    hrichQa_mom_ring_no_backtracking_high_beta->Write();
    hrichQa_mom_backtracking->Write(); 
    hrichQa_mom_backtracking_no_ring->Write(); 
        
    hBetaMom_dilep1_op9->Write();
    hBetaMom->Write();



    hvertexXZ->Write();
    hvertexXY->Write();
    hvertexYZ->Write();
    hvertexX->Write();
    hvertexY->Write();
    hvertexZ->Write();

    hvertexXZ_reco->Write();
    hvertexXY_reco->Write();
    hvertexYZ_reco->Write();
    hvertexX_reco->Write();
    hvertexY_reco->Write();
    hvertexZ_reco->Write();
















    hr_vs_z->Write();
    hnCand->Write();
    hnCand_best->Write();
    hnCand_best_excl->Write();
    hpt->Write();
    hpt_PP->Write();
    hpt_NN->Write();
    
    hpt_0_140->Write();
    hpt_PP_0_140->Write();
    hpt_NN_0_140->Write();
    
    hpt_140_300->Write();
    hpt_PP_140_300->Write();
    hpt_NN_140_300->Write();

    hpt_300_600->Write();
    hpt_PP_300_600->Write();
    hpt_NN_300_600->Write();

    hmom1mom2->Write();
    hoangle_reco->Write();
    hmom_pos_reco->Write();
    hmom_ele_reco->Write();
    hmom1mom2_PP->Write();
    hoangle_reco_PP->Write();
    hmom1mom2_NN->Write();
    hoangle_reco_NN->Write();
    hmom_pos_reco_PP->Write();
    hmom_ele_reco_NN->Write();
    htheta_pos_reco_PP->Write();
    htheta_ele_reco_NN->Write();
    hmissmassNP_Reco_pi0->Write();
    hcos_theta->Write(); 
    hmassNP_oa9deg_no_eff->Write();
    hmassPP_oa9deg_no_eff->Write();
    hmassNN_oa9deg_no_eff->Write();
    hmissmassNP_oa9deg->Write();
    hmassNP_missmass->Write();
    hmissmassPP_oa9deg->Write();
    hmassPP_missmass->Write();
    hmissmassNN_oa9deg->Write();
    hmassNN_missmass->Write();
    hrich_sector_single_pion->Write();
    htheta_phi_raw_pair_ele_rebin->Write();
    htheta_phi_raw_pair_pos_rebin->Write();
    
    htheta_phi_raw_single_pos_rebin->Write();
    htheta_phi_raw_single_ele_rebin->Write();
    
    htheta_phi_raw_single_pos->Write();
    htheta_phi_raw_single_ele->Write();
    
    for(Int_t j = 0; j < 6; j++) {
      hringPadNr_single[j]->Write();
      hringPadNr_15_30_single[j]->Write();
      hringPadNr_30_50_single[j]->Write();
      hringPadNr_50_70_single[j]->Write();
      hringPadNr_70_90_single[j]->Write();
      hringClusterNr_single[j]->Write();
      hringClusterNr_15_30_single[j]->Write();
      hringClusterNr_30_50_single[j]->Write();
      hringClusterNr_50_70_single[j]->Write();
      hringClusterNr_70_90_single[j]->Write();
      hringAC_single[j]->Write();;
      hringAC_15_30_single[j]->Write();;
      hringAC_30_50_single[j]->Write();;
      hringAC_50_70_single[j]->Write();;
      hringAC_70_90_single[j]->Write();;
    }
    
    for(Int_t j = 0; j < 6; j++) {
      hringPadNr[j]->Write();
      hringPadNr_15_30[j]->Write();
      hringPadNr_30_50[j]->Write();
      hringPadNr_50_70[j]->Write();
      hringPadNr_70_90[j]->Write();
      hringClusterNr[j]->Write();
      hringClusterNr_15_30[j]->Write();
      hringClusterNr_30_50[j]->Write();
      hringClusterNr_50_70[j]->Write();
      hringClusterNr_70_90[j]->Write();
      hringAC[j]->Write();;
      hringAC_15_30[j]->Write();;
      hringAC_30_50[j]->Write();;
      hringAC_50_70[j]->Write();;
      hringAC_70_90[j]->Write();;
      
    }
    
    hsec_single_pos->Write(); 
    hsec_single_ele->Write(); 
    hsec_pair_pos->Write();
    hsec_pair_ele->Write();
    hsec_pair_pos_op9->Write();
    hsec_pair_ele_op9->Write();
    
    hCounter->Write();
    
    htheta_phi_raw_pair->Write();
    htheta_phi_raw_pair_ele->Write();
    htheta_phi_raw_pair_pos->Write();


    TMacro m("dieleAna.cc");
    m.Write();




    out->Save();
    out->Close();

    cout<<"####################################################"<<endl;

    return 0;
}
