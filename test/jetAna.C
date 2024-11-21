#include "TTree.h"
#include "TChain.h"
#include "TFile.h"
#include <iostream>
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TMath.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "/phenix/mpcex/lajoie/ISU/new_unfolding/RooUnfold/src/RooUnfoldResponse.h"
#include "TF1.h"
#include <TGraphAsymmErrors.h>
#include "TSystem.h"
#include "TVector3.h"

#include <iostream>
#include <algorithm>
#include <string>
#include <vector>

#include "bins.h"

#include "SpinDBOutput.hh"
#include "SpinDBContent.hh"

#ifdef __CINT__ 
#pragma link C++ nestedclasses;
#pragma link C++ nestedtypedefs;
#pragma link C++ class std::vector<std::vector<float> >+;
#pragma link C++ class std::vector<std::vector<int> >+;
#endif

// ----------------------------------------------------------------------------
// EEC calculations
// ----------------------------------------------------------------------------
/*  Make sure the EEC path is pointed to wherever
 *  the repo is located. We'll need <utility> for
 *  std::pair.
 */

#include <utility>
#include "./include/PHEnergyCorrelator.h"

// ----------------------------------------------------------------------------

const int maxRecoJets = 100;
int nRecoJets;
float r_centrality;
float r_vertex; 
int r_spinPat;
int r_ip12_clock_cross;
int r_level1_clock_cross;
int r_runNumber;
float r_pT[maxRecoJets];
float r_ml_pT[maxRecoJets];
float r_phi[maxRecoJets];
float r_eta[maxRecoJets];
float r_zg[maxRecoJets];
float r_oang[maxRecoJets];
float r_cf[maxRecoJets];
int r_nc[maxRecoJets]; 

bool arm0FiredTrigger; 
bool arm1FiredTrigger; 

std::vector< std::vector<float> > *re_cs_charge; 
std::vector< std::vector<float> > *re_cs_z;
std::vector< std::vector<float> > *re_cs_jT; 
std::vector< std::vector<float> > *re_cs_dR; 
std::vector< std::vector<float> > *re_cs_eta; 
std::vector< std::vector<float> > *re_cs_phi; 
  
void BranchAddressRecoJets(TTree *tReco, int doFakeJets){

  if(doFakeJets){

    tReco->SetBranchAddress("nFakeJets", &nRecoJets);
    tReco->SetBranchAddress("centrality", &r_centrality);
    tReco->SetBranchAddress("vertex", &r_vertex);
    tReco->SetBranchAddress("ip12_clock_cross", &r_ip12_clock_cross);
    tReco->SetBranchAddress("level1_clock_cross", &r_level1_clock_cross);
    tReco->SetBranchAddress("spinPattern", &r_spinPat);
    tReco->SetBranchAddress("runNumber", &r_runNumber);
    tReco->SetBranchAddress("pT", r_pT);
    tReco->SetBranchAddress("ml_pT", r_ml_pT);
    tReco->SetBranchAddress("phi", r_phi);
    tReco->SetBranchAddress("eta", r_eta);
    tReco->SetBranchAddress("cf", r_cf);
    tReco->SetBranchAddress("nc", r_nc);
    tReco->SetBranchAddress("arm0FiredTrigger", &arm0FiredTrigger);
    tReco->SetBranchAddress("arm1FiredTrigger", &arm1FiredTrigger);
    
    for(int i=0; i<maxRecoJets; i++){
      r_zg[i] = 0.0;
    }

    // no constituents array
    re_cs_charge = new std::vector< std::vector<float> >; 
    re_cs_charge->clear();
    re_cs_z =  new std::vector< std::vector<float> >;
    re_cs_z->clear();
    re_cs_jT =  new std::vector< std::vector<float> >;
    re_cs_jT->clear();
    re_cs_dR =  new std::vector< std::vector<float> >;
    re_cs_dR->clear();
    re_cs_eta =  new std::vector< std::vector<float> >;
    re_cs_eta->clear();
    re_cs_phi =  new std::vector< std::vector<float> >;
    re_cs_phi->clear();

  }
  else{

    tReco->SetBranchAddress("nRecoJets", &nRecoJets);
    tReco->SetBranchAddress("centrality", &r_centrality);
    tReco->SetBranchAddress("vertex", &r_vertex);
    tReco->SetBranchAddress("ip12_clock_cross", &r_ip12_clock_cross);
    tReco->SetBranchAddress("level1_clock_cross", &r_level1_clock_cross);
    tReco->SetBranchAddress("spinPattern", &r_spinPat);
    tReco->SetBranchAddress("runNumber", &r_runNumber);
    tReco->SetBranchAddress("pT", r_pT);
    tReco->SetBranchAddress("ml_pT", r_ml_pT);
    tReco->SetBranchAddress("phi", r_phi);
    tReco->SetBranchAddress("eta", r_eta);
    tReco->SetBranchAddress("z_g", r_zg);
    tReco->SetBranchAddress("oang", r_oang);
    tReco->SetBranchAddress("cf", r_cf);
    tReco->SetBranchAddress("nc", r_nc);
    tReco->SetBranchAddress("arm0FiredTrigger", &arm0FiredTrigger);
    tReco->SetBranchAddress("arm1FiredTrigger", &arm1FiredTrigger);

    tReco -> SetBranchAddress("cs_charge",&re_cs_charge);
    tReco -> SetBranchAddress("cs_z",&re_cs_z);
    tReco -> SetBranchAddress("cs_jT",&re_cs_jT);
    tReco -> SetBranchAddress("cs_dR",&re_cs_dR);
    tReco -> SetBranchAddress("cs_eta",&re_cs_eta);
    tReco -> SetBranchAddress("cs_phi",&re_cs_phi);

  }

  return;

}//BranchAddressRecoJets

double JetCharge( int jetIdx ){

  // Calculate jet charge - arXiv:1209.2421v2

  const double kappa = 0.5; 
  double jc = 0; 

  for(unsigned int i=0; i<re_cs_z->at(jetIdx).size(); i++){
    if(re_cs_charge->at(jetIdx).at(i)==0.0) continue;   

    double ph = sqrt(pow(re_cs_z->at(jetIdx).at(i)*r_pT[jetIdx],2) + 
		     pow(re_cs_jT->at(jetIdx).at(i),2) ); 
    double htheta = 2.0*atan(exp(-re_cs_eta->at(jetIdx).at(i)));
    double cpt = ph*sin(htheta); 

    jc += (re_cs_charge->at(jetIdx).at(i))*pow(cpt, kappa); 

  }

  jc /= pow(r_pT[jetIdx],kappa); 

  return jc; 

}

int getArm(float phi) {

  if (phi >= -TMath::PiOver2() && phi < TMath::PiOver2()) return 0; 
  if (phi >=  TMath::PiOver2() && phi < 3.0*TMath::PiOver2()) return 1; 
  return -1;

}

int GetBinReco(float pT, int NPTBINS, double *PTBINS){

  if( (pT <  PTBINS[0]) || (pT >=  PTBINS[NPTBINS]) ){
    //cout << "out of range pT =  " << pT << " " << PTBINS[0] << " " << PTBINS[NPTBINS] << endl;
    return -999;
  }

  int retVal = -1;
  for(int i = 0; i < NPTBINS; i++){
    if(pT >= PTBINS[i] && pT < PTBINS[i+1]){ 
      return i;
    }
  }//i

}//GetBinReco

int GetMaxPtIndex(float *pT, float *cf, float *phi, int nRecoJets, int runno, int NPTBINS, double *PTBINS, int AcceptFlag){

  float maxPt = 0.0;
  int indexMax = -1;
  for(int i = 0; i < nRecoJets; i++){

    if((runno==8) && (getArm(phi[i])==0)) continue; // East arm only
    //if((runno==8) && (cf[i]>=0.75)) continue;       // PPG184 cut
    if((runno==8) && ((cf[i]<0.3)||(cf[i]>=0.7))) continue;       // std cut
    
    // If AcceptFlag >=0 then use it to select the arm
    // That will be used for jet finding. 

    if(AcceptFlag>=0){
      if(getArm(phi[i])!=AcceptFlag) continue; 
    }

    if(pT[i] > maxPt){
      maxPt = pT[i];
      indexMax = i;
    }

  }

  return indexMax;

}//GetMaxPt

int GetSpinPattern(int bspin, int yspin){

  int ret_spinPat = -9999; 
  if(bspin == 1 && yspin == 1) {

    ret_spinPat = 0;
  } 
  else if(bspin == -1 && yspin == 1) { 
 
   ret_spinPat = 1;

  } 
  else if( bspin == 1 && yspin == -1) { 
   
    ret_spinPat = 2;
  }
  else if( bspin == -1 && yspin == -1) { 
  
    ret_spinPat = 3;
  }
  else{                     
    
    ret_spinPat = 4;
  }

  return ret_spinPat;
}

void jetAna(int RUNNUM = 12, int isHI = 0, float R = 0.3, float centLow = 0.0, float centHigh = 20.0, int MB = 0, int useML = 0,  
	    int AcceptFlag = -1, std::string inSuffix = "", int doFakeJets = 0){

  // --------------------------------------------------------------------------
  // EEC calculations
  // --------------------------------------------------------------------------
  /*  Illustration of how to declare and initialize
   *  the EEC calculator. The exact choice of ptjet,
   *  CF, and spin bins are TBD.
   */ 

// define flags to turn off certain calculations
#define doRecoEEC 1

  // pt jet bins
  std::vector< std::pair<float, float> > ptJetBins;
  ptJetBins.push_back( std::make_pair(5., 10.) );
  ptJetBins.push_back( std::make_pair(10., 15.) );
  ptJetBins.push_back( std::make_pair(15., 20.) );

  // cf jet bins
  std::vector< std::pair<float, float> > cfJetBins;
  cfJetBins.push_back( std::make_pair(0., 0.5) );
  cfJetBins.push_back( std::make_pair(0.5, 1.) );

  // now declare calculators
  //   - the 1st argument is what's going to be used
  //     to calculate the weights (Pt, Et, or E)
  //   - and the 2nd argument is what power we're
  //     going to raise the weights to (by default
  //     it's set to 1.0, the usual EEC definition)
  PHEC::Calculator recoEEC( PHEC::Type::Pt, 1.0 );

  // set histogram tags
  recoEEC.SetHistTag( "DataJet" );

  // set bins (spin bins are set in a similar way)
  recoEEC.SetPtJetBins( ptJetBins );
  recoEEC.SetCFJetBins( cfJetBins );

  // run initialization routine to generate 
  // desired histograms
  //   - the 1st argument turns on/off 2-point 
  //     histograms (by default set to ON)
  //   - the 2nd argument turns on/off 3-point
  //     histograms (TODO)
  //   - and the 3rd argument turns on/off
  //     lambda-tagged EEC histograms (TODO)
  recoEEC.Init(true, false, false);

  // --------------------------------------------------------------------------
 

  float vertexCut = 0.0;
  if((RUNNUM == 12)||(RUNNUM==15)) {
    vertexCut = 10.0;
    if(isHI!=0) vertexCut = 10.0; 
  }
  else if(RUNNUM == 13) vertexCut = 30.0;
  else if(RUNNUM == 8) vertexCut = 20.0;
  else{ 
    cout << "ERROR in runnumber" << endl;
    return; 
  }

  TString Suffix = inSuffix; 

  if(isHI && useML){
    cout << "ERROR: Machine Learning not currently supported for HI" << endl;  
    return; 
  }

  if(useML && !((RUNNUM==12)||(RUNNUM==13))){
    cout << "ERROR: Machine Learning not currently supported for run = " << RUNNUM << endl;  
    return; 
  }

  TString HIString = ""; 
  if(isHI){
    if(RUNNUM==8) HIString = "dAu"; 
    if(RUNNUM==12) HIString = "CuAu"; 
    if(RUNNUM==15) HIString = "pAu"; 
  }
  
  int ptBins = 0;
  if((RUNNUM ==  12)||(RUNNUM==8)||(RUNNUM==15)){
    ptBins = NPTBINS_RECO_12;
  }
  else if(RUNNUM == 13){
    ptBins = NPTBINS_RECO_13;
  }
  else{
    cout << "Wrong Run" << endl;
  }
  
  const int NPTBINS = ptBins;
  double PTBINS[NPTBINS + 1];
  for(unsigned int i = 0; i < (unsigned int)(NPTBINS+1); i++){
    PTBINS[i] = 0.0;
  }

  if((RUNNUM == 12)||(RUNNUM == 8)||(RUNNUM==15)){
    memcpy(PTBINS, PTBINS_RECO_12, sizeof(PTBINS));
  }
  else if(RUNNUM == 13){
    memcpy(PTBINS, PTBINS_RECO_13, sizeof(PTBINS));
  }
  else{
    cout << "WRONG RUNNUMBER FOR BINNING" << endl;
  }
  for(unsigned int i = 0; i < (unsigned int)(NPTBINS+1); i++){
    cout << PTBINS[i] << endl;
  }

  TString MLString = ""; 
  if(useML) MLString = "_ML"; 

  // Update the output name based on arm selection
  TString oSuffix = Suffix; 
  if(AcceptFlag==0) oSuffix = "_Arm0" + Suffix; 
  if(AcceptFlag==1) oSuffix = "_Arm1" + Suffix; 

  // BBC efficiency weights
  TGraphAsymmErrors *bbcEff = NULL;  
  TFile *fEff = NULL; 
  if(((RUNNUM==8)||(RUNNUM==12)||(RUNNUM==13)||(RUNNUM==15))&&(isHI==0)){
    fEff = new TFile(Form("Run%iBBCEff.root",RUNNUM), "read");
    bbcEff = (TGraphAsymmErrors *)fEff->Get("BBCEff");
    if(!bbcEff){
      cout << " ERROR - missing bbc efficiency graph!" << endl; 
      return;
    }
    else{
      cout << " BBC efficiency loaded." << endl; 
    }
  }

  TString fakeString = ""; 
  if(doFakeJets) fakeString = "Fake"; 

  TFile *fOut = NULL; 
  if(isHI==0){
    if(MB==1)
      fOut = new TFile(Form("jetAna%sRun%i_ppMB_R%4.2f%s%s.root", fakeString.Data(),RUNNUM,R,MLString.Data(),oSuffix.Data()), "Recreate");
    else
      fOut = new TFile(Form("jetAna%sRun%i_pp_R%4.2f%s%s.root", fakeString.Data(),RUNNUM,R,MLString.Data(),oSuffix.Data()), "Recreate");
  }
  else{
    if(MB==1)
      fOut = new TFile(Form("jetAna%sRun%i_%sMB_R%4.2f_%i_%i%s%s.root", fakeString.Data(), RUNNUM, HIString.Data(), R, (int)centLow, (int)centHigh, MLString.Data(),oSuffix.Data()), "Recreate");
    else
      fOut = new TFile(Form("jetAna%sRun%i_%s_R%4.2f_%i_%i%s%s.root", fakeString.Data(), RUNNUM, HIString.Data(), R, (int)centLow, (int)centHigh, MLString.Data(),oSuffix.Data()), "Recreate");
  }

  TH1D *hReco = new TH1D("hReco", "Highest pT Reco Jets", NPTBINS, PTBINS);
  hReco->Sumw2();
  TH1D *hReco_ERTfired = new TH1D("hReco_ERTfired", "Highest pT Reco Jets", NPTBINS, PTBINS);
  hReco_ERTfired->Sumw2();
  TH1D *hReco_ERTfiredZg0 = new TH1D("hReco_ERTfiredZg0", "Highest pT Reco Jets (no Zg)", NPTBINS, PTBINS);
  hReco_ERTfiredZg0->Sumw2();
  TH1D *hReco_ERTfiredZg1 = new TH1D("hReco_ERTfiredZg1", "Highest pT Reco Jets (valid Zg)", NPTBINS, PTBINS);
  hReco_ERTfiredZg1->Sumw2();

  TH2D *hRecoJetCharge_ERTfired = new TH2D("hRecoJetCharge_ERTfired", "Jet Charge vs. Reco Pt", NPTBINS, PTBINS, NJCBINS, JCBINS );
  hRecoJetCharge_ERTfired->Sumw2();

  TH1D *hRecoSamePat = new TH1D("hRecoSamePat_ERTfired", "Same spin: ++/--", NPTBINS, PTBINS);
  hRecoSamePat->Sumw2();
  TH1D *hRecoOppPat = new TH1D("hRecoOppPat_ERTfired", "Opp spin: +-/-+", NPTBINS, PTBINS);
  hRecoOppPat->Sumw2();

  TH1F *hRecoSamePatP = new TH1F("hRecoSamePat_ERTfiredP", "Same spin: ++/--", NPTBINS, PTBINS);
  hRecoSamePatP->Sumw2();
  TH1F *hRecoOppPatP = new TH1F("hRecoOppPat_ERTfiredP", "Opp spin: +-/-+", NPTBINS, PTBINS);
  hRecoOppPatP->Sumw2(); 

  TH1F *hRecoSamePatP2 = new TH1F("hRecoSamePat_ERTfiredP2", "Same spin: ++/--", NPTBINS, PTBINS);
  hRecoSamePatP2->Sumw2();
  TH1F *hRecoOppPatP2 = new TH1F("hRecoOppPat_ERTfiredP2", "Opp spin: +-/-+", NPTBINS, PTBINS);
  hRecoOppPatP2->Sumw2();

  TH1D *hRecoSamePat_even = new TH1D("hRecoSamePat_even_ERTfired", "Same spin: ++/--, even crossing", NPTBINS, PTBINS);
  hRecoSamePat_even->Sumw2();
  TH1D *hRecoSamePat_odd = new TH1D("hRecoSamePat_odd_ERTfired", "Same spin: ++/--, odd crossing", NPTBINS, PTBINS);
  hRecoSamePat_odd->Sumw2();
  TH1D *hRecoOppPat_even = new TH1D("hRecoOppPat_even_ERTfired", "Opp spin: +-/-+, even crossing", NPTBINS, PTBINS);
  hRecoOppPat_even->Sumw2();
  TH1D *hRecoOppPat_odd = new TH1D("hRecoOppPat_odd_ERTfired", "Opp spin: +-/-+, odd crossing", NPTBINS, PTBINS);
  hRecoOppPat_odd->Sumw2();

  TH1F *hRecoSamePatP_even = new TH1F("hRecoSamePatP_even_ERTfired", "Same spin: ++/--, even crossing", NPTBINS, PTBINS);
  hRecoSamePatP_even->Sumw2();
  TH1F *hRecoSamePatP_odd = new TH1F("hRecoSamePatP_odd_ERTfired", "Same spin: ++/--, odd crossing", NPTBINS, PTBINS);
  hRecoSamePatP_odd->Sumw2();
  TH1F *hRecoOppPatP_even = new TH1F("hRecoOppPatP_even_ERTfired", "Opp spin: +-/-+, even crossing", NPTBINS, PTBINS);
  hRecoOppPatP_even->Sumw2();
  TH1F *hRecoOppPatP_odd = new TH1F("hRecoOppPatP_odd_ERTfired", "Opp spin: +-/-+, odd crossing", NPTBINS, PTBINS);
  hRecoOppPatP_odd->Sumw2();

  TH1F *hRecoSamePatP2_even = new TH1F("hRecoSamePatP2_even_ERTfired", "Same spin: ++/--, even crossing", NPTBINS, PTBINS);
  hRecoSamePatP2_even->Sumw2();
  TH1F *hRecoSamePatP2_odd = new TH1F("hRecoSamePatP2_odd_ERTfired", "Same spin: ++/--, odd crossing", NPTBINS, PTBINS);
  hRecoSamePatP2_odd->Sumw2();
  TH1F *hRecoOppPatP2_even = new TH1F("hRecoOppPatP2_even_ERTfired", "Opp spin: +-/-+, even crossing", NPTBINS, PTBINS);
  hRecoOppPatP2_even->Sumw2();
  TH1F *hRecoOppPatP2_odd = new TH1F("hRecoOppPatP2_odd_ERTfired", "Opp spin: +-/-+, odd crossing", NPTBINS, PTBINS);
  hRecoOppPatP2_odd->Sumw2();

  TH1D *hRecoSpinPat[2][4];
  TH1F *hRecoSpinPatP[2][4];
  TH1F *hRecoSpinPatP2[2][4];
  for(unsigned int j = 0; j < 2; j++){
    for(unsigned int i = 0; i < 4; i++){
      hRecoSpinPat[j][i] = new TH1D(Form("hRecoSpinPat_%i_%i", j, i), Form("Reco Jets, spinPat: %i %i", j, i), NPTBINS, PTBINS);
      hRecoSpinPat[j][i]->Sumw2();
      hRecoSpinPatP[j][i] = new TH1F(Form("hRecoSpinPatP_%i_%i", j, i), Form("Reco Jets, spinPat: %i %i", j, i), NPTBINS, PTBINS);
      hRecoSpinPatP[j][i]->Sumw2();
      hRecoSpinPatP2[j][i] = new TH1F(Form("hRecoSpinPatP2_%i_%i", j, i), Form("Reco Jets, spinPat: %i %i", j, i), NPTBINS, PTBINS);
      hRecoSpinPatP2[j][i]->Sumw2();
    }
  }

  TH1D *hReco_bin[NPTBINS];
  TH1D *hReco_bin_ERTfired[NPTBINS];
  for(unsigned int i = 0; i < (unsigned int) NPTBINS; i++){
    hReco_bin[i] = new TH1D(Form("hRecoBin_%i", i), Form("Bin %i", i+1), 100, PTBINS[i], PTBINS[i+1]);
    hReco_bin_ERTfired[i] = new TH1D(Form("hRecoBinERTfired_%i", i), Form("Bin %i", i+1), 100, PTBINS[i], PTBINS[i+1]);
  }//i

  TH1D *hEventCount = new TH1D("hEventCount", "Count of all events", 1, -0.5, 0.5);
  TH1D *hEventCountValidVertex = new TH1D("hEventCountValidVertex", "Count of all events with valid vertex", 1, -0.5, 0.5);
  TH1D *hVertex = new TH1D("hVertex", "Vertex Distribution; Vertex (cm)", 150, -150.0, 150.0);
  TH1D *hVertexCut = new TH1D("hVertexCut", "Vertex Distribution after cut; Vertex (cm)", 30, -30.0, 30.0);
  TH1D *hVertexJets = new TH1D("hVertexJets", "Vertex Distribution; Vertex (cm)", 30, -30.0, 30.0);

  TH2D *hZg = new TH2D("hZg", "z_{g} vs. reco p_{T}", NPTBINS, PTBINS, NZGBINS, ZGBINS);
  hZg->Sumw2(); 
  TH2D *hOang = new TH2D("hOang", "oang vs. reco p_{T}", NPTBINS, PTBINS, NOANGBINS, OANGBINS);
  hOang->Sumw2(); 

  TH2D *hZgSpinPat[2][4];
  for(unsigned int j = 0; j < 2; j++){
    for(unsigned int i = 0; i < 4; i++){
      hZgSpinPat[j][i] = new TH2D(Form("hZgSpinPat_%i_%i", j, i), Form("hZgSpinPat_%i_%i", j, i), NPTBINS, PTBINS, NZGBINS, ZGBINS);
      hZgSpinPat[j][i]->Sumw2(); 
    }
  }

  TH2D *hZgSamePat_even = new TH2D("hZgSamePat_even", "Same spin: ++/--, even crossing", NPTBINS, PTBINS, NZGBINS, ZGBINS);
  hZgSamePat_even->Sumw2();   
  TH2D *hZgOppPat_even = new TH2D("hZgOppPat_even", "Opposite spin: +-/-+, even crossing", NPTBINS, PTBINS, NZGBINS, ZGBINS);
  hZgOppPat_even->Sumw2();   
  TH2D *hZgSamePat_odd = new TH2D("hZgSamePat_odd", "Same spin: ++/--, odd crossing", NPTBINS, PTBINS, NZGBINS, ZGBINS);
  hZgSamePat_odd->Sumw2();   
  TH2D *hZgOppPat_odd = new TH2D("hZgOppPat_odd", "Opposite spin: +-/-+, odd crossing", NPTBINS, PTBINS, NZGBINS, ZGBINS);
  hZgOppPat_odd->Sumw2();   

  TH2D *hFFZ = new TH2D("hFFZ", "z vs. reco p_{T}", NPTBINS, PTBINS, NFFZBINS, FFZBINS);
  hFFZ->Sumw2(); 

  TH2D *hFFZSpinPat[2][4];
  for(unsigned int j = 0; j < 2; j++){
    for(unsigned int i = 0; i < 4; i++){
      hFFZSpinPat[j][i] = new TH2D(Form("hFFZSpinPat_%i_%i", j, i), Form("hFFZSpinPat_%i_%i", j, i), NPTBINS, PTBINS, NFFZBINS, FFZBINS);
      hFFZSpinPat[j][i]->Sumw2(); 
    }
  }

  TH2D *hFFZSamePat_even = new TH2D("hFFZSamePat_even", "Same spin: ++/--, even crossing", NPTBINS, PTBINS, NFFZBINS, FFZBINS);
  hFFZSamePat_even->Sumw2();   
  TH2D *hFFZOppPat_even = new TH2D("hFFZOppPat_even", "Opposite spin: +-/-+, even crossing", NPTBINS, PTBINS, NFFZBINS, FFZBINS);
  hFFZOppPat_even->Sumw2();   
  TH2D *hFFZSamePat_odd = new TH2D("hFFZSamePat_odd", "Same spin: ++/--, odd crossing", NPTBINS, PTBINS, NFFZBINS, FFZBINS);
  hFFZSamePat_odd->Sumw2();   
  TH2D *hFFZOppPat_odd = new TH2D("hFFZOppPat_odd", "Opposite spin: +-/-+, odd crossing", NPTBINS, PTBINS, NFFZBINS, FFZBINS);
  hFFZOppPat_odd->Sumw2();   

  TH2D *hFFXI = new TH2D("hFFXI", "z vs. reco p_{T}", NPTBINS, PTBINS, NFFXIBINS, FFXIBINS);
  hFFXI->Sumw2(); 

  TH2D *hFFXISpinPat[2][4];
  for(unsigned int j = 0; j < 2; j++){
    for(unsigned int i = 0; i < 4; i++){
      hFFXISpinPat[j][i] = new TH2D(Form("hFFXISpinPat_%i_%i", j, i), Form("hFFXISpinPat_%i_%i", j, i), NPTBINS, PTBINS, NFFXIBINS, FFXIBINS);
      hFFXISpinPat[j][i]->Sumw2(); 
    }
  }

  TH2D *hFFXISamePat_even = new TH2D("hFFXISamePat_even", "Same spin: ++/--, even crossing", NPTBINS, PTBINS, NFFXIBINS, FFXIBINS);
  hFFXISamePat_even->Sumw2();   
  TH2D *hFFXIOppPat_even = new TH2D("hFFXIOppPat_even", "Opposite spin: +-/-+, even crossing", NPTBINS, PTBINS, NFFXIBINS, FFXIBINS);
  hFFXIOppPat_even->Sumw2();   
  TH2D *hFFXISamePat_odd = new TH2D("hFFXISamePat_odd", "Same spin: ++/--, odd crossing", NPTBINS, PTBINS, NFFXIBINS, FFXIBINS);
  hFFXISamePat_odd->Sumw2();   
  TH2D *hFFXIOppPat_odd = new TH2D("hFFXIOppPat_odd", "Opposite spin: +-/-+, odd crossing", NPTBINS, PTBINS, NFFXIBINS, FFXIBINS);
  hFFXIOppPat_odd->Sumw2();   


  //TH2D *hFFJT = new TH2D("hFFJT", "jT vs. reco p_{T}", NPTBINS, PTBINS, NFFJTBINS, FFJTBINS);
  TH2D *hFFJT = new TH2D("hFFJT", "jT/pT vs. reco p_{T}", NPTBINS, PTBINS, NFFJTBINS, FFJTBINS);
  hFFJT->Sumw2(); 

  TH2D *hFFJTSpinPat[2][4];
  for(unsigned int j = 0; j < 2; j++){
    for(unsigned int i = 0; i < 4; i++){
      hFFJTSpinPat[j][i] = new TH2D(Form("hFFJTSpinPat_%i_%i", j, i), Form("hFFJTSpinPat_%i_%i", j, i), NPTBINS, PTBINS, NFFJTBINS, FFJTBINS);
      hFFJTSpinPat[j][i]->Sumw2(); 
    }
  }

  TH2D *hFFJTSamePat_even = new TH2D("hFFJTSamePat_even", "Same spin: ++/--, even crossing", NPTBINS, PTBINS, NFFJTBINS, FFJTBINS);
  hFFJTSamePat_even->Sumw2();   
  TH2D *hFFJTOppPat_even = new TH2D("hFFJTOppPat_even", "Opposite spin: +-/-+, even crossing", NPTBINS, PTBINS, NFFJTBINS, FFJTBINS);
  hFFJTOppPat_even->Sumw2();   
  TH2D *hFFJTSamePat_odd = new TH2D("hFFJTSamePat_odd", "Same spin: ++/--, odd crossing", NPTBINS, PTBINS, NFFJTBINS, FFJTBINS);
  hFFJTSamePat_odd->Sumw2();   
  TH2D *hFFJTOppPat_odd = new TH2D("hFFJTOppPat_odd", "Opposite spin: +-/-+, odd crossing", NPTBINS, PTBINS, NFFJTBINS, FFJTBINS);
  hFFJTOppPat_odd->Sumw2();   

  TH2D *hFFdR = new TH2D("hFFdR", "dR vs. reco p_{T}", NPTBINS, PTBINS, NFFDRBINS, FFDRBINS);
  hFFdR->Sumw2(); 

  TH2D *hFFdRSpinPat[2][4];
  for(unsigned int j = 0; j < 2; j++){
    for(unsigned int i = 0; i < 4; i++){
      hFFdRSpinPat[j][i] = new TH2D(Form("hFFdRSpinPat_%i_%i", j, i), Form("hFFdRSpinPat_%i_%i", j, i), NPTBINS, PTBINS, NFFDRBINS, FFDRBINS);
      hFFdRSpinPat[j][i]->Sumw2(); 
    }
  }

  TH2D *hFFdRSamePat_even = new TH2D("hFFdRSamePat_even", "Same spin: ++/--, even crossing", NPTBINS, PTBINS, NFFDRBINS, FFDRBINS);
  hFFdRSamePat_even->Sumw2();   
  TH2D *hFFdROppPat_even = new TH2D("hFFdROppPat_even", "Opposite spin: +-/-+, even crossing", NPTBINS, PTBINS, NFFDRBINS, FFDRBINS);
  hFFdROppPat_even->Sumw2();   
  TH2D *hFFdRSamePat_odd = new TH2D("hFFdRSamePat_odd", "Same spin: ++/--, odd crossing", NPTBINS, PTBINS, NFFDRBINS, FFDRBINS);
  hFFdRSamePat_odd->Sumw2();   
  TH2D *hFFdROppPat_odd = new TH2D("hFFdROppPat_odd", "Opposite spin: +-/-+, odd crossing", NPTBINS, PTBINS, NFFDRBINS, FFDRBINS);
  hFFdROppPat_odd->Sumw2();   

  // Phi distribution w.r.t the jet axis as a function of jet pT, z
  // Only valid for spin 

  TH2D *hJetPhiBluePol[2][2];
  TH2D *hJetPhiYellowPol[2][2];
  TH2D *hJetPhiBluePolJCPos[2][2];
  TH2D *hJetPhiYellowPolJCPos[2][2];
  TH2D *hJetPhiBluePolJCNeg[2][2];
  TH2D *hJetPhiYellowPolJCNeg[2][2];
  for(unsigned int i = 0; i < 2; i++){
    for(unsigned int j = 0; j < 2; j++){
      hJetPhiBluePol[i][j] = new TH2D(Form("hJetPhiBluePol_%i_%i", i, j), Form("hFFPhiBluePol_%i_%i", i, j), NPTBINS, PTBINS, NFFPHIBINS, FFPHIBINS);
      hJetPhiBluePol[i][j]->Sumw2(); 
      hJetPhiYellowPol[i][j] = new TH2D(Form("hJetPhiYellowPol_%i_%i", i, j), Form("hFFPhiYellowPol_%i_%i", i, j), NPTBINS, PTBINS, NFFPHIBINS, FFPHIBINS);
      hJetPhiYellowPol[i][j]->Sumw2(); 
      hJetPhiBluePolJCPos[i][j] = new TH2D(Form("hJetPhiBluePolJCPos_%i_%i", i, j), Form("hFFPhiBluePolJCPos_%i_%i", i, j), NPTBINS, PTBINS, NFFPHIBINS, FFPHIBINS);
      hJetPhiBluePolJCPos[i][j]->Sumw2(); 
      hJetPhiYellowPolJCPos[i][j] = new TH2D(Form("hJetPhiYellowPolJCPos_%i_%i", i, j), Form("hFFPhiYellowPolJCPos_%i_%i", i, j), NPTBINS, PTBINS, NFFPHIBINS, FFPHIBINS);
      hJetPhiYellowPolJCPos[i][j]->Sumw2(); 
      hJetPhiBluePolJCNeg[i][j] = new TH2D(Form("hJetPhiBluePolJCNeg_%i_%i", i, j), Form("hFFPhiBluePolJCNeg_%i_%i", i, j), NPTBINS, PTBINS, NFFPHIBINS, FFPHIBINS);
      hJetPhiBluePolJCNeg[i][j]->Sumw2(); 
      hJetPhiYellowPolJCNeg[i][j] = new TH2D(Form("hJetPhiYellowPolJCNeg_%i_%i", i, j), Form("hFFPhiYellowPolJCNeg_%i_%i", i, j), NPTBINS, PTBINS, NFFPHIBINS, FFPHIBINS);
      hJetPhiYellowPolJCNeg[i][j]->Sumw2(); 
    }
  }

  TH3D *hFFPhiBluePolPos[2][2];
  TH3D *hFFPhiYellowPolPos[2][2];
  TH3D *hFFPhiBluePolNeg[2][2];
  TH3D *hFFPhiYellowPolNeg[2][2];
  for(unsigned int i = 0; i < 2; i++){
    for(unsigned int j = 0; j < 2; j++){
      hFFPhiBluePolPos[i][j] = new TH3D(Form("hFFPhiBluePolPos_%i_%i", i, j), Form("hFFPhiBluePolPos_%i_%i", i, j), NPTBINS, PTBINS, NFFZBINS, FFZBINS, NFFPHIBINS, FFPHIBINS);
      hFFPhiBluePolPos[i][j]->Sumw2(); 
      hFFPhiYellowPolPos[i][j] = new TH3D(Form("hFFPhiYellowPolPos_%i_%i", i, j), Form("hFFPhiYellowPolPos_%i_%i", i, j), NPTBINS, PTBINS, NFFZBINS, FFZBINS, NFFPHIBINS, FFPHIBINS);
      hFFPhiYellowPolPos[i][j]->Sumw2(); 
      hFFPhiBluePolNeg[i][j] = new TH3D(Form("hFFPhiBluePolNeg_%i_%i", i, j), Form("hFFPhiBluePolPosNeg_%i_%i", i, j), NPTBINS, PTBINS, NFFZBINS, FFZBINS, NFFPHIBINS, FFPHIBINS);
      hFFPhiBluePolNeg[i][j]->Sumw2(); 
      hFFPhiYellowPolNeg[i][j] = new TH3D(Form("hFFPhiYellowPolNeg_%i_%i", i, j), Form("hFFPhiYellowPolNeg_%i_%i", i, j), NPTBINS, PTBINS, NFFZBINS, FFZBINS, NFFPHIBINS, FFPHIBINS);
      hFFPhiYellowPolNeg[i][j]->Sumw2(); 
    }
  }

  TH3D *hFFTwoPhiBluePolPos[2][2];
  TH3D *hFFTwoPhiYellowPolPos[2][2];
  TH3D *hFFTwoPhiBluePolNeg[2][2];
  TH3D *hFFTwoPhiYellowPolNeg[2][2];
  for(unsigned int i = 0; i < 2; i++){
    for(unsigned int j = 0; j < 2; j++){
      hFFTwoPhiBluePolPos[i][j] = new TH3D(Form("hFFTwoPhiBluePolPos_%i_%i", i, j), Form("hFFTwoPhiBluePolPos_%i_%i", i, j), NPTBINS, PTBINS, NFFZBINS, FFZBINS, NFFPHIBINS, FFPHIBINS);
      hFFTwoPhiBluePolPos[i][j]->Sumw2(); 
      hFFTwoPhiYellowPolPos[i][j] = new TH3D(Form("hFFTwoPhiYellowPolPos_%i_%i", i, j), Form("hFFTwoPhiYellowPolPos_%i_%i", i, j), NPTBINS, PTBINS, NFFZBINS, FFZBINS, NFFPHIBINS, FFPHIBINS);
      hFFTwoPhiYellowPolPos[i][j]->Sumw2(); 
      hFFTwoPhiBluePolNeg[i][j] = new TH3D(Form("hFFTwoPhiBluePolNeg_%i_%i", i, j), Form("hFFTwoPhiBluePolPosNeg_%i_%i", i, j), NPTBINS, PTBINS, NFFZBINS, FFZBINS, NFFPHIBINS, FFPHIBINS);
      hFFTwoPhiBluePolNeg[i][j]->Sumw2(); 
      hFFTwoPhiYellowPolNeg[i][j] = new TH3D(Form("hFFTwoPhiYellowPolNeg_%i_%i", i, j), Form("hFFTwoPhiYellowPolNeg_%i_%i", i, j), NPTBINS, PTBINS, NFFZBINS, FFZBINS, NFFPHIBINS, FFPHIBINS);
      hFFTwoPhiYellowPolNeg[i][j]->Sumw2(); 
    }
  }

  // jet yield in run number vs. IP12 crossing

  TH2D *JetYield = new TH2D("JetYield","yield in runnumber vs. ip12 crossing",121,-0.5,120.5,2208,387538.5,396367.5); 
  TH2D *JetYieldNar = new TH2D("JetYieldNar","yield in runnumber vs. ip12 crossing",121,-0.5,120.5,20,391859.5,391879.5); 
  
  TH2D *BluePol = new TH2D("BluePol","Blue Pol vs. ip12 crossing",120,0,120,100,0.0,1.0); 
  TH2D *YellowPol = new TH2D("YellowPol","Yellow Pol vs. ip12 crossing",120,0,120,100,0.0,1.0); 

  // charged fraction, etc.

  int ncfbins = 100; 
  double cfbins[101]; 
  for(unsigned int i=0; i<101; i++) cfbins[i] = (1.0/100.0)*i;

  TH2D *hCF = new TH2D("hCF", "Charged Fraction",NPTBINS, PTBINS, ncfbins, cfbins);
  hCF->Sumw2();

  int ncbins = 21; 
  double cbins[22]; 
  for(unsigned int i=0; i<22; i++) cbins[i] = -0.5 + 1.0*i;

  TH2D *hNC = new TH2D("hNC", "Jet NC",NPTBINS,PTBINS,ncbins,cbins);
  hNC->Sumw2();

  TH1D *hPhi = new TH1D("hPhi", "Jet Phi",100,-TMath::PiOver2(),3.0*TMath::PiOver2());
  hPhi->Sumw2();

  // Variables for polarization tracking 

  int currRunNumber = -1;
  float Pblue = 1.0; 
  float Pyellow = 1.0; 
  float bpolIndexed[120] = {0.0}; 
  float ypolIndexed[120] = {0.0}; 
  float Rweight = 1.0; 
  float RweightEO[2] = {1.0,1.0}; 
  SpinDBOutput spin_out("phnxrc");

  int tMu_runnum = 0;
  double meannum = 0;
  double bbcnovtxraw = 0;
  double bbcnovtxlive = 0;
  long long ss_bbcnovtxlive[120]; 
  int ss_ok = 0; 

  TFile *fMu = NULL; 
  TTree* tMu = NULL; 
  if(RUNNUM==13){
    fMu = new TFile("lumitree_run13Fin.root", "read");
    tMu = (TTree*)fMu->Get("T");
    if(tMu){
      tMu->SetBranchAddress("runnum", &tMu_runnum);
      tMu->SetBranchAddress("mu", &meannum);
      tMu->SetBranchAddress("bbcnovtxraw", &bbcnovtxraw);
      tMu->SetBranchAddress("bbcnovtxlive", &bbcnovtxlive);
      tMu->SetBranchAddress("ss_ok", &ss_ok);
      tMu->SetBranchAddress("ss_bbcnovtxlive", &ss_bbcnovtxlive);
    }
    else{
      cout << "Unable to access luminosity tree!" << endl; 
    }
  }

  // Input chains
  
  TFile *fIn = NULL; 
  if(isHI==0){
    if(MB==1)
      fIn = new TFile(Form("./chainedRun%i_ppMB_R%4.2f%s.root", RUNNUM, R, Suffix.Data()), "read");
    else
      fIn = new TFile(Form("./chainedRun%i_pp_R%4.2f%s.root", RUNNUM, R, Suffix.Data()), "read");
  }
  else {
    if(MB==1)
      fIn = new TFile(Form("./chainedRun%i_%sMB_R%4.2f%s.root", RUNNUM, HIString.Data(), R, Suffix.Data()), "read");
    else
      fIn = new TFile(Form("./chainedRun%i_%s_R%4.2f%s.root", RUNNUM, HIString.Data(), R, Suffix.Data()), "read");
  }

  //TTree *tReco = (TTree*)fIn->Get("rJets");
  TTree *tReco = NULL; 
  if(doFakeJets)
    tReco = (TTree*)fIn->Get("fJets");
  else
    tReco = (TTree*)fIn->Get("rJets");

  BranchAddressRecoJets(tReco, doFakeJets);

  unsigned long long nEntriesReco = tReco->GetEntries(); 
  cout << " nEntriesReco = " << nEntriesReco << endl;

  // Random number generator - for tests and debug
  TRandom3 rand; 

  //   LOOP OVER RECO JETS TREE ///
  for(unsigned long long iRecoEntry = 0; iRecoEntry < nEntriesReco; iRecoEntry++){

    tReco->GetEntry(iRecoEntry);

    if(nRecoJets>maxRecoJets){
      cout << "ERROR! nRecoJets = " << nRecoJets << " (> maxRecoJets = " << maxRecoJets << ")" << endl; 
      return;
    }
    
    if((RUNNUM==12) && (isHI==1) && ((r_centrality<centLow) || (r_centrality>=centHigh))) continue; 
    if((RUNNUM==8) && (isHI==1) && ((r_centrality<centLow) || (r_centrality>=centHigh))) continue; 
    if((RUNNUM==15) && (isHI==1) && ((r_centrality<centLow) || (r_centrality>=centHigh))) continue; 

    // TEST - eliminate two crossings after empties
    //if((RUNNUM==13) && 
    //   (((r_runNumber<392200) && ((r_ip12_clock_cross==31)||(r_ip12_clock_cross==32)||(r_ip12_clock_cross==71)||(r_ip12_clock_cross==72))) ||
    //    (r_ip12_clock_cross<=3)) ) continue;
    
    // Check for run number change, update polarization
    if((RUNNUM==13) && (r_runNumber!=currRunNumber)){

     // Update polarization values 

      SpinDBContent spin_cont;
 
      int storedbcont = spin_out.StoreDBContent(r_runNumber, r_runNumber);
      if(storedbcont != 1) {
	cout << "ERROR: StoreDBContent failed for run = " << r_runNumber << endl;      
	return; 
      }

      int getdbcontstore = spin_out.GetDBContentStore(spin_cont,r_runNumber);
      if( getdbcontstore != 1) {
	cout << "ERROR: GetDBContentStore failed for run =  " << r_runNumber << endl;      
	continue; 
      }

      for(int ip12_clock_cross = 0; ip12_clock_cross<120; ip12_clock_cross++){

	double bpol, bpol_err; 
	int getbluepol = spin_cont.GetPolarizationBlue(ip12_clock_cross, bpol, bpol_err);
	if((getbluepol==1)&&(bpol<0.7)){ // eliminate unusually high blue polarization values
	  bpolIndexed[ip12_clock_cross] = bpol; 
	  BluePol->Fill(ip12_clock_cross,bpol); 
	}
	else{
	  cerr << "ERROR getting bpol for " << r_runNumber << " " << ip12_clock_cross << endl; 
	  bpolIndexed[ip12_clock_cross] = 0.0; 
	}

	double ypol, ypol_err; 
	int getpolyellow = spin_cont.GetPolarizationYellow(ip12_clock_cross, ypol, ypol_err);
	if(getpolyellow==1){
	  ypolIndexed[ip12_clock_cross] = ypol; 
	  YellowPol->Fill(ip12_clock_cross,ypol); 
	}
        else{
	  cerr << "ERROR getting ypol for " << r_runNumber << " " << ip12_clock_cross << endl; 
	  ypolIndexed[ip12_clock_cross] = 0.0; 
	}
      
      }

      // Update relative luminosity for this run

      int entries = tMu->GetEntries();
      bool chk = false;
      for(int i = 0; i < entries; i++){//loop over tree   
	tMu->GetEntry(i);
	if(tMu_runnum == r_runNumber){
	  chk = true;
	  break;
	}
      }

      double R_lumi[2][4] = {{0.0, 0.0, 0.0, 0.0},{0.0, 0.0, 0.0, 0.0}}; 

      for(int ip12_clock_cross = 0; ip12_clock_cross<120; ip12_clock_cross++){

	int even_odd = ip12_clock_cross%2; 

	int bspin = spin_cont.GetSpinPatternBlue(ip12_clock_cross);
	int yspin = spin_cont.GetSpinPatternYellow(ip12_clock_cross);

	int spinPattern = GetSpinPattern(bspin, yspin);

	if(chk && (ss_ok!=0) && (spinPattern>=0) && (spinPattern<4)){
	  R_lumi[even_odd][spinPattern] += (double)(ss_bbcnovtxlive[ip12_clock_cross]);
	}
	else{
	  R_lumi[even_odd][spinPattern] += (double)spin_cont.GetScalerBbcNoCut(ip12_clock_cross);
	}
      }

      RweightEO[0] = (R_lumi[0][1]+R_lumi[0][2])/(R_lumi[0][0]+R_lumi[0][3]); 
      RweightEO[1] = (R_lumi[1][1]+R_lumi[1][2])/(R_lumi[1][0]+R_lumi[1][3]);   
      Rweight = (R_lumi[0][1]+R_lumi[0][2]+R_lumi[1][1]+R_lumi[1][2])/(R_lumi[0][0]+R_lumi[0][3]+R_lumi[1][0]+R_lumi[1][3]); 

      // Notify user

      cout << " Run = " << r_runNumber << " RweightEO[0] = " << RweightEO[0] << " RweightEO[1] = " << RweightEO[1] << endl; 

      // update current run number
      
      currRunNumber = r_runNumber; 

    }
    
    // Set polarization 

    if(RUNNUM==13){
      Pblue = bpolIndexed[r_ip12_clock_cross]; 
      Pyellow = ypolIndexed[r_ip12_clock_cross]; 
    }
    else{
      Pblue = 1.0; 
      Pyellow = 1.0; 
    }

    hEventCount->Fill(0.0);
    if(r_vertex<-90000) continue; // vertex not valid
    hEventCountValidVertex->Fill(0.0); 

    hVertex->Fill(r_vertex);
    //vertex cut
    if(fabs(r_vertex) > vertexCut) continue;
    hVertexCut->Fill(r_vertex);

    if(nRecoJets > 0){

      hVertexJets->Fill(r_vertex);

      int indexMax = -1;
      if(useML)
        indexMax = GetMaxPtIndex(r_ml_pT, r_cf, r_phi, nRecoJets, RUNNUM, NPTBINS, PTBINS, AcceptFlag);
      else
        indexMax = GetMaxPtIndex(r_pT, r_cf, r_phi, nRecoJets, RUNNUM, NPTBINS, PTBINS, AcceptFlag);

      // do not consider jets with unphysical reco pT
      if(indexMax>=0){
	if(useML){
	  if(r_ml_pT[indexMax]>=PTBINS[NPTBINS]) indexMax = -1; 
	  if(r_ml_pT[indexMax]<PTBINS[0]) indexMax = -1; 
	}
	else{
	  if(r_pT[indexMax]>=PTBINS[NPTBINS]) indexMax = -1; 
	  if(r_pT[indexMax]<PTBINS[0]) indexMax = -1; 
	}
      }

      if(indexMax<0) continue; 
      float r_maxPt = -9999.0; 
      if(useML)
	r_maxPt = r_ml_pT[indexMax];
      else
	r_maxPt = r_pT[indexMax];
      
      if(fabs(r_eta[indexMax])>0.15) continue; 

      int binReco = GetBinReco(r_maxPt, NPTBINS, PTBINS);
      if(binReco<0) continue; 

      // BBC efficiency weighting
      double weight = 1.0; 
      if(bbcEff) weight = (1.0/bbcEff->Eval(r_maxPt)); 
      
      hReco->Fill(r_maxPt,weight);
      hCF->Fill(r_maxPt, r_cf[indexMax],weight); 
      hPhi->Fill(r_phi[indexMax],weight); 
      hNC->Fill(r_maxPt, r_nc[indexMax],weight); 

      bool trigFired = false;
      if((arm0FiredTrigger && (getArm(r_phi[indexMax])==0)) || (arm1FiredTrigger && (getArm(r_phi[indexMax])==1))) trigFired = true;  
      if((RUNNUM==12) && (isHI==1)) trigFired = true; 

      hReco_bin[binReco]->Fill(r_maxPt,weight);

      if(trigFired){

	hReco_ERTfired->Fill(r_maxPt,weight);

	hReco_bin_ERTfired[binReco]->Fill(r_maxPt,weight);
 	
	double jetCharge = 0.0; 
	if(doFakeJets==0) jetCharge = JetCharge(indexMax); 
	hRecoJetCharge_ERTfired->Fill(r_maxPt,jetCharge,weight);
	
	// yields unweighted
	JetYield->Fill(r_ip12_clock_cross, r_runNumber); 
	JetYieldNar->Fill(r_ip12_clock_cross, r_runNumber);

	if(r_zg[indexMax]>=0.1 && r_zg[indexMax]<=0.5 ) { // valid Zg
	  hZg->Fill(r_maxPt, r_zg[indexMax],weight);
	  hOang->Fill(r_maxPt, r_oang[indexMax],weight);
	  if(r_zg[indexMax]>=0.1 && r_zg[indexMax]<0.3 )
	    hReco_ERTfiredZg0->Fill(r_maxPt,weight);
	  else
	    hReco_ERTfiredZg1->Fill(r_maxPt,weight);
	}

	if(doFakeJets==0){
	  for(unsigned int i=0; i<re_cs_z->at(indexMax).size(); i++){
	    if(re_cs_charge->at(indexMax).at(i)==0.0) continue; 
	    hFFZ->Fill(r_maxPt, re_cs_z->at(indexMax).at(i),weight); 
	    hFFXI->Fill(r_maxPt, -log(re_cs_z->at(indexMax).at(i)),weight); 
	  }

	  for(unsigned int i=0; i<re_cs_jT->at(indexMax).size(); i++){
	    if(re_cs_charge->at(indexMax).at(i)==0.0) continue;
	    //hFFJT->Fill(r_maxPt, re_cs_jT->at(indexMax).at(i),weight); 
	    hFFJT->Fill(r_maxPt, re_cs_jT->at(indexMax).at(i)/r_maxPt,weight); 
	  }

	  for(unsigned int i=0; i<re_cs_dR->at(indexMax).size(); i++){
	    if(re_cs_charge->at(indexMax).at(i)==0.0) continue;
	    hFFdR->Fill(r_maxPt, re_cs_dR->at(indexMax).at(i),weight); 
	  }
	}

	//SPIN pattern histos
	// 0 = even, 1 = odd
	int even_odd = r_ip12_clock_cross%2;

	// DEBUG - Randomize spin pattern
	//r_spinPat = (int)rand.Uniform(0.0,4.0); 

	// DEBUG - swap spin pattern
	/*
	if((r_spinPat==0)||(r_spinPat==3)){

	  if(r_spinPat==0) 
	    r_spinPat = 1; 
	  else
	    r_spinPat = 2; 

	}
	else{

	  if(r_spinPat==1) 
	    r_spinPat = 0; 
	  else
	    r_spinPat = 3; 

	}
	*/

        // ---------------------------------------------------------------------
        // EEC calculations over max pt jet
        // ---------------------------------------------------------------------
        /* Here we actually run the relevant calculations
         * on the jets and their constituents. Note that
         * if you haven't turned on histograms for a
         * particular calculation (e.g. the 3-point), then
         * the code won't try to fill the corresponding
         * histograms.
         */
        if (doMaxJetEEC) {

          // collect jet information into a handy struct
          //   - NOTE: the spin for the bunch x-ing is
          //     bundled w/ the jets (the last argument)
          //   - For now, it's just a dummy value
          PHEC::Type::Jet jet(
            r_cf[indexMax],
            r_pT[indexMax],
            r_eta[indexMax],
            r_phi[indexMax],
            1.
          );

          // loop through pairs of constituents
          for (std::size_t iCstA = 0; iCstA < re_cs_z->at(indexMax).size(); ++iCstA) {
            for (std::size_t iCstB = 0; iCstB < re_cs_z->at(indexMax).size(); ++iCstB) {

              // skip diagonal
              if (iCstA == iCstB) continue;

              // collect cst information into a handy struct
              PHEC::Type::Cst cstA(
                re_cs_z->at(indexMax).at(iCstA),
                re_cs_jT->at(indexMax).at(iCstA),
                re_cs_eta->at(indexMax).at(iCstA),
                re_cs_phi->at(indexMax).at(iCstA),
                re_cs_charge->at(indexMax).at(iCstA)
              );
              PHEC::Type::Cst cstB(
                re_cs_z->at(indexMax).at(iCstB),
                re_cs_jT->at(indexMax).at(iCstB),
                re_cs_eta->at(indexMax).at(iCstB),
                re_cs_phi->at(indexMax).at(iCstB),
                re_cs_charge->at(indexMax).at(iCstB)
              );

              // run 2-point calculation for pair
              recoEEC.CalcEEC( jet, std::make_pair(cstA, cstB) );

            }  // end 2nd cst loop
          }  // end 1st cst loop
        }  // end max jet eec calculation

        // ---------------------------------------------------------------------

	if(r_spinPat < 4){
	  hRecoSpinPat[even_odd][r_spinPat]->Fill(r_maxPt,weight);
	  if((r_spinPat==0)||(r_spinPat==3)){ //ss
	    hRecoSpinPatP[even_odd][r_spinPat]->Fill(r_maxPt,weight*Pblue*Pyellow);
	    hRecoSpinPatP2[even_odd][r_spinPat]->Fill(r_maxPt,weight*Pblue*Pblue*Pyellow*Pyellow);
	  }
	  else{ //os
	    hRecoSpinPatP[even_odd][r_spinPat]->Fill(r_maxPt,weight*Pblue*Pyellow*RweightEO[even_odd]);
	    hRecoSpinPatP2[even_odd][r_spinPat]->Fill(r_maxPt,weight*Pblue*Pblue*Pyellow*Pyellow*RweightEO[even_odd]);
	  }
	  if(r_zg[indexMax]>0.0) hZgSpinPat[even_odd][r_spinPat]->Fill(r_maxPt, r_zg[indexMax],weight); 

	  if(doFakeJets==0){
	    for(unsigned int i=0; i<re_cs_z->at(indexMax).size(); i++){
	      if(re_cs_charge->at(indexMax).at(i)==0.0) continue; 
	      hFFZSpinPat[even_odd][r_spinPat]->Fill(r_maxPt, re_cs_z->at(indexMax).at(i),weight); 
	      hFFXISpinPat[even_odd][r_spinPat]->Fill(r_maxPt, -log(re_cs_z->at(indexMax).at(i)),weight); 
	    }

	    for(unsigned int i=0; i<re_cs_jT->at(indexMax).size(); i++){
	      if(re_cs_charge->at(indexMax).at(i)==0.0) continue;  
	      //hFFJTSpinPat[even_odd][r_spinPat]->Fill(r_maxPt, re_cs_jT->at(indexMax).at(i),weight); 
	      hFFJTSpinPat[even_odd][r_spinPat]->Fill(r_maxPt, re_cs_jT->at(indexMax).at(i)/r_maxPt,weight); 
	    }

	    for(unsigned int i=0; i<re_cs_dR->at(indexMax).size(); i++){
	      if(re_cs_charge->at(indexMax).at(i)==0.0) continue;   
	      hFFdRSpinPat[even_odd][r_spinPat]->Fill(r_maxPt, re_cs_dR->at(indexMax).at(i),weight); 
	    }
	  }

	  if(((RUNNUM==12)||(RUNNUM==15))&&(isHI==0)){

	    // Collins histograms
	    // define unit vectors for angles

	    TVector3 blue_beam(0.0,0.0,1.0); 
	    TVector3 yellow_beam(0.0,0.0,-1.0); 

	    TVector3 blue_spin(0.0,1.0,0.0);
	    if((r_spinPat==2)||(r_spinPat==3)) blue_spin.SetY(-1.0); 
	    TVector3 yellow_spin(0.0,1.0,0.0);
	    if((r_spinPat==1)||(r_spinPat==3)) yellow_spin.SetY(-1.0); 

	    double theta = 2.0*atan(exp(-r_eta[indexMax]));
	    double pz = r_maxPt/tan(theta); 
	    double px = (pz/cos(theta))*cos(r_phi[indexMax]); 
	    double py = (pz/cos(theta))*sin(r_phi[indexMax]); 

	    TVector3 jet(px,py,pz); 
	    double jet_norm = jet.Mag(); 
	    jet.SetX(jet.X()/jet_norm); 
	    jet.SetY(jet.Y()/jet_norm); 
	    jet.SetZ(jet.Z()/jet_norm); 

	    // PhiSpin is the angle between the plane defined by the
	    // polarized beam and the jet

	    TVector3 jet_blue_beam_perp = (blue_beam.Cross(jet)).Unit();  
	    TVector3 jet_yellow_beam_perp = (yellow_beam.Cross(jet)).Unit();  
	 
	    double bluePhiSpin = TMath::PiOver2() - acos(jet_blue_beam_perp.Dot(blue_spin)); 
	    double yellowPhiSpin = TMath::PiOver2() - acos(jet_yellow_beam_perp.Dot(yellow_spin)); 

	    if((r_spinPat==2)||(r_spinPat==3)){
	      hJetPhiBluePol[1][even_odd]->Fill(r_maxPt,  bluePhiSpin, weight); // spin down 
	      if(jetCharge>0.0)
		hJetPhiBluePolJCPos[1][even_odd]->Fill(r_maxPt,  bluePhiSpin, weight); // spin down 
	      else
		hJetPhiBluePolJCNeg[1][even_odd]->Fill(r_maxPt,  bluePhiSpin, weight); // spin down 
	    }
	    else{
	      hJetPhiBluePol[0][even_odd]->Fill(r_maxPt,  bluePhiSpin, weight); // spin up 
	      if(jetCharge>0.0)
	        hJetPhiBluePolJCPos[0][even_odd]->Fill(r_maxPt,  bluePhiSpin, weight); // spin up 
	      else
	        hJetPhiBluePolJCNeg[0][even_odd]->Fill(r_maxPt,  bluePhiSpin, weight); // spin up 
	    }
	
	    if((r_spinPat==1)||(r_spinPat==3)){
	      hJetPhiYellowPol[1][even_odd]->Fill(r_maxPt,  yellowPhiSpin, weight); // spin down
	      if(jetCharge>0.0)
		hJetPhiYellowPolJCPos[1][even_odd]->Fill(r_maxPt,  yellowPhiSpin, weight); // spin down
	      else
	        hJetPhiYellowPolJCNeg[1][even_odd]->Fill(r_maxPt,  yellowPhiSpin, weight); // spin down
	    }
	    else{
	      hJetPhiYellowPol[0][even_odd]->Fill(r_maxPt,  yellowPhiSpin, weight); // spin up 
	      if(jetCharge>0.0)
	        hJetPhiYellowPolJCPos[0][even_odd]->Fill(r_maxPt,  yellowPhiSpin, weight); // spin up 
	      else
	        hJetPhiYellowPolJCNeg[0][even_odd]->Fill(r_maxPt,  yellowPhiSpin, weight); // spin up 
	    }
	  
	    if(doFakeJets==0){
	      for(unsigned int i=0; i<re_cs_z->at(indexMax).size(); i++){
		if(re_cs_charge->at(indexMax).at(i)==0.0) continue;   

		// PhiHadron is the angle between the plane defined by the 
		// polarized beam and the jet (above) and the plane containing the 
		// jet and the hadron

		double ph = sqrt(pow(re_cs_z->at(indexMax).at(i)*r_maxPt,2) + 
				 pow(re_cs_jT->at(indexMax).at(i),2) ); 

		double htheta = 2.0*atan(exp(-re_cs_eta->at(indexMax).at(i)));
		double px = ph*sin(htheta)*cos(re_cs_phi->at(indexMax).at(i));
		double py = ph*sin(htheta)*sin(re_cs_phi->at(indexMax).at(i));
		double pz = ph*cos(htheta);

		TVector3 hadron(px,py,pz);
		double hadron_norm = hadron.Mag(); 
		hadron.SetX(hadron.X()/hadron_norm); 
		hadron.SetY(hadron.Y()/hadron_norm); 
		hadron.SetZ(hadron.Z()/hadron_norm); 

		TVector3 hadron_jet_perp = (jet.Cross(hadron)).Unit(); 
	    
		double bluePhiHadron = acos(jet_blue_beam_perp.Dot(hadron_jet_perp)); 
		double yellowPhiHadron = acos(jet_yellow_beam_perp.Dot(hadron_jet_perp)); 
		if(bluePhiHadron>TMath::PiOver2()) bluePhiHadron -= TMath::Pi(); 
		if(yellowPhiHadron>TMath::PiOver2()) yellowPhiHadron -= TMath::Pi(); 

		double blue_dPhiSpin = bluePhiSpin - bluePhiHadron; 
		if (blue_dPhiSpin > TMath::Pi())
		  blue_dPhiSpin -= TMath::Pi(); 
		else if (blue_dPhiSpin < 0 )
		  blue_dPhiSpin += TMath::Pi(); 

		double blue2PhiHadron = 2.0*bluePhiHadron; 
		if( (blue2PhiHadron>TMath::PiOver2()) && (blue2PhiHadron<=3.0*TMath::PiOver2()) ) 
		  blue2PhiHadron -= TMath::Pi(); 
		else if( blue2PhiHadron>3.0*TMath::PiOver2() ) 
		  blue2PhiHadron -= 2.0*TMath::Pi(); 
		else if( (blue2PhiHadron<-TMath::PiOver2()) && (blue2PhiHadron>=-3.0*TMath::PiOver2()) )
		  blue2PhiHadron += TMath::Pi(); 
		else if( blue2PhiHadron<-3.0*TMath::PiOver2() ) 
		  blue2PhiHadron += 2.0*TMath::Pi(); 

		double blue_d2PhiSpin = bluePhiSpin - blue2PhiHadron; 
		if (blue_d2PhiSpin > TMath::Pi())
		  blue_d2PhiSpin -= TMath::Pi(); 
		else if (blue_d2PhiSpin < 0 )
		  blue_d2PhiSpin += TMath::Pi(); 

		double yellow_dPhiSpin = yellowPhiSpin - yellowPhiHadron; 
		if (yellow_dPhiSpin > TMath::Pi())
		  yellow_dPhiSpin -= TMath::Pi(); 
		else if (yellow_dPhiSpin < 0 )
		  yellow_dPhiSpin += TMath::Pi(); 

		double yellow2PhiHadron = 2.0*yellowPhiHadron; 
		if( (yellow2PhiHadron>TMath::PiOver2()) && (yellow2PhiHadron<=3.0*TMath::PiOver2()) ) 
		  yellow2PhiHadron -= TMath::Pi(); 
		else if( yellow2PhiHadron>3.0*TMath::PiOver2() ) 
		  yellow2PhiHadron -= 2.0*TMath::Pi(); 
		else if( (yellow2PhiHadron<-TMath::PiOver2()) && (yellow2PhiHadron>=-3.0*TMath::PiOver2()) )
		  yellow2PhiHadron += TMath::Pi(); 
		else if( yellow2PhiHadron<-3.0*TMath::PiOver2() ) 
		  yellow2PhiHadron += 2.0*TMath::Pi(); 

		double yellow_d2PhiSpin = yellowPhiSpin - yellow2PhiHadron; 
		if (yellow_d2PhiSpin > TMath::Pi())
		  yellow_d2PhiSpin -= TMath::Pi(); 
		else if (yellow_d2PhiSpin < 0 )
		  yellow_d2PhiSpin += TMath::Pi(); 

		if(re_cs_charge->at(indexMax).at(i)>0.0){

		  if((r_spinPat==2)||(r_spinPat==3)){
		    hFFPhiBluePolPos[1][even_odd]->Fill(r_maxPt, re_cs_z->at(indexMax).at(i), blue_dPhiSpin, weight); // spin down
		    hFFTwoPhiBluePolPos[1][even_odd]->Fill(r_maxPt, re_cs_z->at(indexMax).at(i), blue_d2PhiSpin, weight); // spin down
		  }
		  else{
		    hFFPhiBluePolPos[0][even_odd]->Fill(r_maxPt, re_cs_z->at(indexMax).at(i), blue_dPhiSpin, weight); // spin up
		    hFFTwoPhiBluePolPos[0][even_odd]->Fill(r_maxPt, re_cs_z->at(indexMax).at(i), blue_d2PhiSpin, weight); // spin up
		  }
		
		  if((r_spinPat==1)||(r_spinPat==3)){
		    hFFPhiYellowPolPos[1][even_odd]->Fill(r_maxPt, re_cs_z->at(indexMax).at(i), yellow_dPhiSpin, weight); // spin down 
		    hFFTwoPhiYellowPolPos[1][even_odd]->Fill(r_maxPt, re_cs_z->at(indexMax).at(i), yellow_d2PhiSpin, weight); // spin down 
		  }
		  else{ 
		    hFFPhiYellowPolPos[0][even_odd]->Fill(r_maxPt, re_cs_z->at(indexMax).at(i), yellow_dPhiSpin, weight); // spin up
		    hFFTwoPhiYellowPolPos[0][even_odd]->Fill(r_maxPt, re_cs_z->at(indexMax).at(i), yellow_d2PhiSpin, weight); // spin up
		  }
	    
		}
		else if(re_cs_charge->at(indexMax).at(i)<0.0){

		  if((r_spinPat==2)||(r_spinPat==3)){
		    hFFPhiBluePolNeg[1][even_odd]->Fill(r_maxPt, re_cs_z->at(indexMax).at(i), blue_dPhiSpin, weight); // spin down
		    hFFTwoPhiBluePolNeg[1][even_odd]->Fill(r_maxPt, re_cs_z->at(indexMax).at(i), blue_d2PhiSpin, weight); // spin down
		  }
		  else{
		    hFFPhiBluePolNeg[0][even_odd]->Fill(r_maxPt, re_cs_z->at(indexMax).at(i), blue_dPhiSpin, weight); // spin up
		    hFFTwoPhiBluePolNeg[0][even_odd]->Fill(r_maxPt, re_cs_z->at(indexMax).at(i), blue_d2PhiSpin, weight); // spin up
		  }

		  if((r_spinPat==1)||(r_spinPat==3)){
		    hFFPhiYellowPolNeg[1][even_odd]->Fill(r_maxPt, re_cs_z->at(indexMax).at(i), yellow_dPhiSpin, weight); // spin down
		    hFFTwoPhiYellowPolNeg[1][even_odd]->Fill(r_maxPt, re_cs_z->at(indexMax).at(i), yellow_d2PhiSpin, weight); // spin down
		  }
		  else{
		    hFFPhiYellowPolNeg[0][even_odd]->Fill(r_maxPt, re_cs_z->at(indexMax).at(i), yellow_dPhiSpin, weight); // spin up
		    hFFTwoPhiYellowPolNeg[0][even_odd]->Fill(r_maxPt, re_cs_z->at(indexMax).at(i), yellow_d2PhiSpin, weight); // spin up
		  }

		}

	      }

	    }

	  }

	}
 
      }
     
    }

  }//iRecoEntry

  //add same and opp spin histos
  //unweighted
  hRecoSamePat_even->Add(hRecoSpinPat[0][0]);
  hRecoSamePat_even->Add(hRecoSpinPat[0][3]);
  hRecoSamePat_odd->Add(hRecoSpinPat[1][0]);
  hRecoSamePat_odd->Add(hRecoSpinPat[1][3]);

  hRecoSamePat->Add(hRecoSamePat_even);
  hRecoSamePat->Add(hRecoSamePat_odd);

  hRecoOppPat_even->Add(hRecoSpinPat[0][1]);
  hRecoOppPat_even->Add(hRecoSpinPat[0][2]);
  hRecoOppPat_odd->Add(hRecoSpinPat[1][1]);
  hRecoOppPat_odd->Add(hRecoSpinPat[1][2]);

  hRecoOppPat->Add(hRecoOppPat_even);
  hRecoOppPat->Add(hRecoOppPat_odd);

  // pol weighted
  hRecoSamePatP_even->Add(hRecoSpinPatP[0][0]);
  hRecoSamePatP_even->Add(hRecoSpinPatP[0][3]);
  hRecoSamePatP_odd->Add(hRecoSpinPatP[1][0]);
  hRecoSamePatP_odd->Add(hRecoSpinPatP[1][3]);

  hRecoSamePatP->Add(hRecoSamePatP_even);
  hRecoSamePatP->Add(hRecoSamePatP_odd);

  hRecoOppPatP_even->Add(hRecoSpinPatP[0][1]);
  hRecoOppPatP_even->Add(hRecoSpinPatP[0][2]);
  hRecoOppPatP_odd->Add(hRecoSpinPatP[1][1]);
  hRecoOppPatP_odd->Add(hRecoSpinPatP[1][2]);

  hRecoOppPatP->Add(hRecoOppPatP_even);
  hRecoOppPatP->Add(hRecoOppPatP_odd);

  // pol^2 weighted
  hRecoSamePatP2_even->Add(hRecoSpinPatP2[0][0]);
  hRecoSamePatP2_even->Add(hRecoSpinPatP2[0][3]);
  hRecoSamePatP2_odd->Add(hRecoSpinPatP2[1][0]);
  hRecoSamePatP2_odd->Add(hRecoSpinPatP2[1][3]);

  hRecoSamePatP2->Add(hRecoSamePatP2_even);
  hRecoSamePatP2->Add(hRecoSamePatP2_odd);

  hRecoOppPatP2_even->Add(hRecoSpinPatP2[0][1]);
  hRecoOppPatP2_even->Add(hRecoSpinPatP2[0][2]);
  hRecoOppPatP2_odd->Add(hRecoSpinPatP2[1][1]);
  hRecoOppPatP2_odd->Add(hRecoSpinPatP2[1][2]);

  hRecoOppPatP2->Add(hRecoOppPatP2_even);
  hRecoOppPatP2->Add(hRecoOppPatP2_odd);


  hZgSamePat_even->Add(hZgSpinPat[0][0]);
  hZgSamePat_even->Add(hZgSpinPat[0][3]);
  hZgSamePat_odd->Add(hZgSpinPat[1][0]);
  hZgSamePat_odd->Add(hZgSpinPat[1][3]);
  
  hZgOppPat_even->Add(hZgSpinPat[0][1]);
  hZgOppPat_even->Add(hZgSpinPat[0][2]);
  hZgOppPat_odd->Add(hZgSpinPat[1][1]);
  hZgOppPat_odd->Add(hZgSpinPat[1][2]);

  hFFZSamePat_even->Add(hFFZSpinPat[0][0]);
  hFFZSamePat_even->Add(hFFZSpinPat[0][3]);
  hFFZSamePat_odd->Add(hFFZSpinPat[1][0]);
  hFFZSamePat_odd->Add(hFFZSpinPat[1][3]);

  hFFZOppPat_even->Add(hFFZSpinPat[0][1]);
  hFFZOppPat_even->Add(hFFZSpinPat[0][2]);
  hFFZOppPat_odd->Add(hFFZSpinPat[1][1]);
  hFFZOppPat_odd->Add(hFFZSpinPat[1][2]);

  hFFXISamePat_even->Add(hFFXISpinPat[0][0]);
  hFFXISamePat_even->Add(hFFXISpinPat[0][3]);
  hFFXISamePat_odd->Add(hFFXISpinPat[1][0]);
  hFFXISamePat_odd->Add(hFFXISpinPat[1][3]);

  hFFXIOppPat_even->Add(hFFXISpinPat[0][1]);
  hFFXIOppPat_even->Add(hFFXISpinPat[0][2]);
  hFFXIOppPat_odd->Add(hFFXISpinPat[1][1]);
  hFFXIOppPat_odd->Add(hFFXISpinPat[1][2]);

  hFFJTSamePat_even->Add(hFFJTSpinPat[0][0]);
  hFFJTSamePat_even->Add(hFFJTSpinPat[0][3]);
  hFFJTSamePat_odd->Add(hFFJTSpinPat[1][0]);
  hFFJTSamePat_odd->Add(hFFJTSpinPat[1][3]);

  hFFJTOppPat_even->Add(hFFJTSpinPat[0][1]);
  hFFJTOppPat_even->Add(hFFJTSpinPat[0][2]);
  hFFJTOppPat_odd->Add(hFFJTSpinPat[1][1]);
  hFFJTOppPat_odd->Add(hFFJTSpinPat[1][2]);

  hFFdRSamePat_even->Add(hFFdRSpinPat[0][0]);
  hFFdRSamePat_even->Add(hFFdRSpinPat[0][3]);
  hFFdRSamePat_odd->Add(hFFdRSpinPat[1][0]);
  hFFdRSamePat_odd->Add(hFFdRSpinPat[1][3]);

  hFFdROppPat_even->Add(hFFdRSpinPat[0][1]);
  hFFdROppPat_even->Add(hFFdRSpinPat[0][2]);
  hFFdROppPat_odd->Add(hFFdRSpinPat[1][1]);
  hFFdROppPat_odd->Add(hFFdRSpinPat[1][2]);

  // --------------------------------------------------------------------------
  // EEC calculations
  // --------------------------------------------------------------------------
  /*  And lastly, here we save all histograms
   *  to whatever file is passed to the End()
   *  method.
   */ 

  if (doRecoEEC) recoEEC.End( fOut );

  // --------------------------------------------------------------------------

  fOut->Write(); 
  fOut->Close();

  return;

}//jetAna
