#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TMath.h>


void QuickModelMod() {

  TFile* f = new TFile("checkingQuickModel_step3_randomizeJet.nIter10K.d28m2y2025.root", "recreate");

  TRandom *rand = new TRandom3(); 

  // beam axis
  TVector3 beam(0,0,1);

  // Histogram

  TH1D *testAngle = new TH1D("testAngle","",100,-TMath::TwoPi(),TMath::TwoPi()); 
  TH2D *testAngle_v_Phi = new TH2D("testAngle_v_Phi","",100,-TMath::TwoPi(),TMath::TwoPi(),100,-TMath::TwoPi(),TMath::TwoPi()); 
  TH2D *testAngle_v_Sense = new TH2D("testAngle_v_sense", "", 100, -5., 5., 100, -TMath::TwoPi(), TMath::TwoPi());

  for(int i=0; i<10000; i++){

    // -------------------------------

#if 1
    // completely random jet axis
    double xj = rand->Uniform(-1.,1.);
    double yj = rand->Uniform(-1.,1.);
    double zj = rand->Uniform(-1.,1.);
    TVector3 jet(xj,yj,zj); 

    // normal to jet-beam plane
    TVector3 njetbeam = beam.Cross(jet);

#else

    // jet axis along the beam axis

    TVector3 jet(0,0,1); 

    // normal to jet-beam plane (in x-y plane by definition)
    double xjb = rand->Uniform(-1.,1.);
    double yjb = rand->Uniform(-1.,1.);
    TVector3 njetbeam(xjb,yjb,0); 
#endif

    // --------------------------------

#if 1
    // completely random hadron
    double xh = rand->Uniform(-1.,1.);
    double yh = rand->Uniform(-1.,1.);
    double zh = rand->Uniform(-1.,1.);
    TVector3 hadron_unnorm(xh, yh, zh);
    TVector3 hadron = hadron_unnorm.Unit();
    double phi = hadron.Phi();
#else
    // random hadron around the jet axis
    double phi   = rand->Uniform()*TMath::TwoPi(); 
    TVector3 hadron_unnorm(cos(phi),sin(phi),rand->Uniform()); // +z only
    TVector3 hadron = hadron_unnorm.Unit(); 
#endif

    // ------------------------------

    // Normal to hadron-jet plane

    TVector3 hjetnorm = jet.Cross(hadron).Unit(); 

    // Calculate the signed angle between the two planes 

    double normAngle = atan2( njetbeam.Cross(hjetnorm).Mag(), njetbeam.Dot(hjetnorm) );

    // If "sense" is negative then the direction of rotation is opposite
    // use this to get the full 2pi range and define a zero and sense
    // of direction. 

    double sense = njetbeam.Dot(hadron);
    if(sense<0.0) normAngle = TMath::TwoPi() - normAngle; 

    // Histogram

    testAngle->Fill(normAngle); 
    testAngle_v_Phi->Fill(phi,normAngle); 
    testAngle_v_Sense->Fill(sense,normAngle); 
  }

  f -> cd();
  testAngle->Write();
  testAngle_v_Phi->Write();
  testAngle_v_sense->Write();
  f -> Close();

}
