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

int crossingShift = 0; 

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
  else if( yspin == 0 ) { // Run-15 p+Au (only p is polarized)               
    
    if(bspin == 1)
      ret_spinPat = 4;
    else if(bspin = -1) 
      ret_spinPat = 5; 

  }
  else
    ret_spinPat = 6; 

  return ret_spinPat;
}

void jetAna(int RUNNUM = 12, int isHI = 0, float R = 0.3, float centLow = 0.0, float centHigh = 20.0, int MB = 0, int useML = 0,  
	    int AcceptFlag = -1, std::string inSuffix = "", int doFakeJets = 0){

  // --------------------------------------------------------------------------
  // EEC calculations
  // --------------------------------------------------------------------------

// define flags to turn off eec calculations
#define doDataEEC 1
#define doDataEECChargedOnly 0

// define flags to turn on/off certain binnings
#define doJetCFBins 0
#define doJetChargeBins 0

  // pt jet bins
  std::vector< std::pair<float, float> > ptJetBins;
  ptJetBins.push_back( std::make_pair(5., 10.) );
  ptJetBins.push_back( std::make_pair(10., 15.) );
  ptJetBins.push_back( std::make_pair(15., 20.) );

  // cf jet bins
  std::vector< std::pair<float, float> > cfJetBins;
  cfJetBins.push_back( std::make_pair(0., 1.0) );

  // jet charge bins
  std::vector< std::pair<float, float> > chJetBins;
  chJetBins.push_back( std::make_pair(-100., 0.0) );
  chJetBins.push_back( std::make_pair(0.0, 100.) );

  // now declare calculators
  //   - 1st argument: quantity used for weights (pt, et, or e)
  //   - 2nd argument: power to raise weights to
  PHEC::Calculator dataEEC( PHEC::Type::Pt, 1.0 );

  // set histogram tags
  dataEEC.SetHistTag( "DataJet" );

  // set pt, cf, and charge jet bins
  dataEEC.SetPtJetBins( ptJetBins );
  if (doJetCFBins) {
    dataEEC.SetCFJetBins( cfJetBins );
  }
  if (doJetChargeBins) {
    dataEEC.SetChargeBins( chJetBins );
  }

  // turn on spin sorting
  dataEEC.SetDoSpinBins( true );

  // run initialization routine to generate 
  // desired histograms
  //   - 1st argument: turn on/off 2-point histograms 
  //   - 2nd argument: turn on/off 3-point histograms (TODO)
  //   - 3rd argument: turn on/off lambda EEC histograms (TODO)
  dataEEC.Init(true, false, false);

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

  TH2D *hRecoJetCharge_ERTfired = new TH2D("hRecoJetCharge_ERTfired", "Jet Charge vs. Reco Pt", NPTBINS, PTBINS, NJCBINS, JCBINS );
  hRecoJetCharge_ERTfired->Sumw2();

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

  TH2D *hFFZ = new TH2D("hFFZ", "z vs. reco p_{T}", NPTBINS, PTBINS, NFFZBINS, FFZBINS);
  hFFZ->Sumw2(); 

  TH2D *hFFXI = new TH2D("hFFXI", "z vs. reco p_{T}", NPTBINS, PTBINS, NFFXIBINS, FFXIBINS);
  hFFXI->Sumw2(); 

  TH2D *hFFJT = new TH2D("hFFJT", "jT/pT vs. reco p_{T}", NPTBINS, PTBINS, NFFJTBINS, FFJTBINS);
  hFFJT->Sumw2(); 

  TH2D *hFFdR = new TH2D("hFFdR", "dR vs. reco p_{T}", NPTBINS, PTBINS, NFFDRBINS, FFDRBINS);
  hFFdR->Sumw2(); 

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

    // Bad pAu runs - not in spin database
    if((r_runNumber==434147)||(r_runNumber==434148)||(r_runNumber==434150)||(r_runNumber==434151)) continue; 

    // Check for run number change, update polarization
    if(((RUNNUM==13)||(RUNNUM==15)) && (r_runNumber!=currRunNumber)){

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
	currRunNumber = r_runNumber; // don't keep hammering the DB
	continue; 
      }

      crossingShift = spin_cont.GetCrossingShift();

      for(int ip12_clock_cross = 0; ip12_clock_cross<120; ip12_clock_cross++){

	double bpol, bpol_err; 
	int getbluepol = spin_cont.GetPolarizationBlue(ip12_clock_cross, bpol, bpol_err);
	//if((getbluepol==1)&&(bpol<0.7)){ // eliminate unusually high blue polarization values
	if(getbluepol==1){ 
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

      if(RUNNUM==13){

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

	double R_lumi[2][7] = {{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}}; 

	for(int ip12_clock_cross = 0; ip12_clock_cross<120; ip12_clock_cross++){

	  int even_odd = ip12_clock_cross%2; 

	  int bspin = spin_cont.GetSpinPatternBlue(ip12_clock_cross);
	  int yspin = spin_cont.GetSpinPatternYellow(ip12_clock_cross);

	  int spinPattern = GetSpinPattern(bspin, yspin);

	  if(chk && (ss_ok!=0) && (spinPattern>=0) && (spinPattern<6)){
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

      }

      // update current run number
      
      currRunNumber = r_runNumber; 

    }

    // Check for valid spin information
    // For Run-15 pAu this needs to be recalculated

    if((RUNNUM==15) && (isHI==1)){
      r_ip12_clock_cross = (r_level1_clock_cross + crossingShift)%120;
    }

    // sanity check clock crossing number

    if( (r_ip12_clock_cross<0) || (r_ip12_clock_cross>120) ){
      cout << "Bad clock cross, event skipped, run = " << r_runNumber << " ip12_clock_cross = " << r_ip12_clock_cross << endl; 
      continue; 
    }
        
    // Set polarization 

    if((RUNNUM==13)||(RUNNUM==15)){
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

	if(r_spinPat < 6){

            // ---------------------------------------------------------------------
            // EEC calculations over max pt jet
            // ---------------------------------------------------------------------
            /* N.B. cuts on jet pt, eta, and CNF are baked into the requirement
             *   that indexMax >= 0. A max reco jet has to satisfy these cuts.
             */
            if (doDataEEC) {

              // collect jet and spin information into a handy struct
              PHEC::Type::Jet jet_data(
                r_cf[indexMax],
                r_pT[indexMax],
                r_eta[indexMax],
                r_phi[indexMax],
                jetCharge,
                r_spinPat
              );

              // loop through pairs of constituents
              for (
                std::size_t iCstA = 0;
                iCstA < re_cs_z->at(indexMax).size();
                ++iCstA
              ) {

                // keep only charged csts.s
                if (doDataEECChargedOnly && (re_cs_charge->at(indexMax).at(iCstA) == 0.0)) {
                  continue;
                }

                for (
                  std::size_t iCstB = 0;
                  iCstB <= iCstA;
                  ++iCstB
                ) {

                  // keep only charged csts.s
                  if (doDataEECChargedOnly && (re_cs_charge->at(indexMax).at(iCstB) == 0.0)) {
                    continue;
                  }

                  // collect cst information into a handy struct
                  PHEC::Type::Cst cstA_data(
                    re_cs_z->at(indexMax).at(iCstA),
                    re_cs_jT->at(indexMax).at(iCstA),
                    re_cs_eta->at(indexMax).at(iCstA),
                    re_cs_phi->at(indexMax).at(iCstA),
                    re_cs_charge->at(indexMax).at(iCstA)
                  );
                  PHEC::Type::Cst cstB_data(
                    re_cs_z->at(indexMax).at(iCstB),
                    re_cs_jT->at(indexMax).at(iCstB),
                    re_cs_eta->at(indexMax).at(iCstB),
                    re_cs_phi->at(indexMax).at(iCstB),
                    re_cs_charge->at(indexMax).at(iCstB)
                  );

                  // run 2-point calculation for pair
                  dataEEC.CalcEEC( jet_data, std::make_pair(cstA_data, cstB_data) );

                }  // end 2nd cst loop
              }  // end 1st cst loop
            }  // end max jet eec calculation

            // ---------------------------------------------------------------------
	  
	}

      }
 
    }
     
  }

  // --------------------------------------------------------------------------
  // EEC calculations
  // --------------------------------------------------------------------------
  /*  And lastly, here we save all histograms
   *  to whatever file is passed to the End()
   *  method.
   */ 

  if (doDataEEC) dataEEC.End( fOut );

  // --------------------------------------------------------------------------

  fOut->Write(); 
  fOut->Close();

  return;

}//jetAna
