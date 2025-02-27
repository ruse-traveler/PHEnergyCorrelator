#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TMath.h>


void QuickModel() {

  TFile* f = new TFile("checkingQuickModel.nIter10K.d27m2y2025.root", "recreate");

  TRandom *rand = new TRandom3(); 

  // jet axis along the beam axis

  TVector3 jet(0,0,1); 

  // normal to jet-beam plane (artificial y-axis)

  TVector3 njetbeam(0,1,0); 

  // Histogram

  TH1D *testAngle = new TH1D("testAngle","",100,-TMath::TwoPi(),TMath::TwoPi()); 
  TH2D *testAngle_v_Phi = new TH2D("testAngle_v_Phi","",100,0,TMath::TwoPi(),100,-TMath::TwoPi(),TMath::TwoPi()); 
  TH2D *testAngle_v_Sense = new TH2D("testAngle_v_sense", "", 100, -5., 5., 100, -TMath::TwoPi(), TMath::TwoPi());

  for(int i=0; i<10000; i++){

    // random hadron around the jet axis

    double phi   = rand->Uniform()*TMath::TwoPi(); 
    TVector3 hadron_unnorm(cos(phi),sin(phi),rand->Uniform()); // +z only
    TVector3 hadron = hadron_unnorm.Unit(); 

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
