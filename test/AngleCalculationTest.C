/// ============================================================================
/*! \file    AngleCalculationTest.C
 *  \authors Derek Anderson
 *  \date    02.13.2025
 *
 *  Macro to test collins/boer-mulders angle calculation. Randomly
 *  generates "spin", "jet", and "hadron" vectors and calculates
 *  the above angles between them.
 */
/// ============================================================================

#define ANGLECALCULATIONTEST_C

// c++ utilities
#include <iostream>
#include <utility>
// root libraries
#include <TDatime.h>
#include <TFile.h>
#include <TH1.h>
#include <TMath.h>
#include <TRandom.h>
#include <TRandom3.h>
#include <TVector3.h>

// convenience types
typedef std::tuple<std::size_t, float, float> Bins;



// ============================================================================
//! Test macro for anglue calculations
// ============================================================================
void AngleCalculationTest(const std::size_t nIter = 100) {

  // initialize rng
  TDatime*  time  = new TDatime();
  TRandom3* rando = new TRandom3();
  rando -> SetSeed(time -> Get());

  // histogram binning
  Bins bAngles = std::make_tuple(180, -12.60, 12.60);

  // turn on errors
  TH1::SetDefaultSumw2(true);

  // initialize histograms
  TH1D* hPhiSpinB = new TH1D("hPhiSpinB", "", get<0>(bAngles), get<1>(bAngles), get<2>(bAngles));
  TH1D* hPhiSpinY = new TH1D("hPhiSpinY", "", get<0>(bAngles), get<1>(bAngles), get<2>(bAngles));
  TH1D* hPhiHadB  = new TH1D("hPhiHadB", "", get<0>(bAngles), get<1>(bAngles), get<2>(bAngles));
  TH1D* hPhiHadY  = new TH1D("hPhiHadY", "", get<0>(bAngles), get<1>(bAngles), get<2>(bAngles));
  TH1D* hPhiHad2B = new TH1D("hPhiHad2B", "", get<0>(bAngles), get<1>(bAngles), get<2>(bAngles));
  TH1D* hPhiHad2Y = new TH1D("hPhiHad2Y", "", get<0>(bAngles), get<1>(bAngles), get<2>(bAngles));
  TH1D* hPhiCollB = new TH1D("hPhiCollB", "", get<0>(bAngles), get<1>(bAngles), get<2>(bAngles));
  TH1D* hPhiCollY = new TH1D("hPhiCollY", "", get<0>(bAngles), get<1>(bAngles), get<2>(bAngles));
  TH1D* hPhiBoerB = new TH1D("hPhiBoerB", "", get<0>(bAngles), get<1>(bAngles), get<2>(bAngles));
  TH1D* hPhiBoerY = new TH1D("hPhiBoerY", "", get<0>(bAngles), get<1>(bAngles), get<2>(bAngles));

  // --------------------------------------------------------------------------
  // run calculation nIter times
  // --------------------------------------------------------------------------
  for (std::size_t iIter = 0; iIter < nIter; ++iIter) {

    // set up vectors ---------------------------------------------------------

    // set beam directions
    TVector3 vecBeamB3(0.0, 0.0, 1.0);
    TVector3 vecBeamY3(0.0, 0.0, -1.0);

    // get random spin directions
    TVector3 vecSpinB3(rando -> Uniform(), rando -> Uniform(), 0.0);
    TVector3 vecSpinY3(rando -> Uniform(), rando -> Uniform(), 0.0);

    // get random jet/hadron directions
    TVector3 vecJet3(rando -> Uniform(), rando -> Uniform(), rando -> Uniform());
    TVector3 vecHad3(rando -> Uniform(), rando -> Uniform(), rando -> Uniform());

    // now normalize spin/jet/hadron vectors
    TVector3 unitSpinB3 = vecSpinB3.Unit();
    TVector3 unitSpinY3 = vecSpinY3.Unit();
    TVector3 unitJet3   = vecJet3.Unit();
    TVector3 unitHad3   = vecHad3.Unit();

    // collins & boer-mulders angle calculations --------------------------

    // (1) get vectors normal to the jet-beam plane
    std::pair<TVector3, TVector3> normJetBeam3 = std::make_pair(
      ( vecBeamB3.Cross(unitJet3) ).Unit(),
      ( vecBeamY3.Cross(unitJet3) ).Unit()
    );

    // (2) get phiSpin: angles between the jet-beam plane and spin
    //   - n.b. for spin pattern >= 4, the yellow spin is randomized
    //   - angle between jet plane and spin 
    //   - note that we get the full [0,2pi) range for these angles
    double phiSpinBlue = TMath::PiOver2() - atan2( normJetBeam3.first.Cross(unitSpinB3).Mag(), normJetBeam3.first.Dot(unitSpinB3) );
    double phiSpinYell = TMath::PiOver2() - atan2( normJetBeam3.second.Cross(unitSpinY3).Mag(), normJetBeam3.second.Dot(unitSpinY3) );
    if (phiSpinBlue < 0)               phiSpinBlue += TMath::TwoPi();
    if (phiSpinBlue >= TMath::TwoPi()) phiSpinBlue -= TMath::TwoPi();
    if (phiSpinYell < 0)               phiSpinYell += TMath::TwoPi();
    if (phiSpinYell >= TMath::TwoPi()) phiSpinYell -= TMath::TwoPi();

    // fill spin histograms
    hPhiSpinB -> Fill(phiSpinBlue);
    hPhiSpinY -> Fill(phiSpinYell);

    // (3) get vector normal to hadron average-jet plane
    TVector3 normHadJet3 = ( unitJet3.Cross(unitHad3) ).Unit();

    // (4) get phiHadron: angle between the jet-beam plane and the
    //   - angle between jet-hadron plane
    //   - constrain to range [0,2pi)
    double phiHadBlue = atan2( normJetBeam3.first.Cross(normHadJet3).Mag(), normJetBeam3.first.Dot(normHadJet3) );
    double phiHadYell = atan2( normJetBeam3.second.Cross(normHadJet3).Mag(), normJetBeam3.second.Dot(normHadJet3) );
    if (phiHadBlue < 0)               phiHadBlue += TMath::TwoPi();
    if (phiHadBlue >= TMath::TwoPi()) phiHadBlue -= TMath::TwoPi();
    if (phiHadYell < 0)               phiHadYell += TMath::TwoPi();
    if (phiHadYell >= TMath::TwoPi()) phiHadYell -= TMath::TwoPi();

    // fill histograms
    hPhiHadBlue -> Fill(phiHadBlue);
    hPhiHadYell -> Fill(phiHadYell);

    // (5) double phiHadron for boer-mulders,
    //   - constrain to [0, 2pi)
    double phiHadBlue2 = 2.0 * phiHadBlue;
    double phiHadYell2 = 2.0 * phiHadYell;
    if (phiHadBlue2 < 0)               phiHadBlue2 += TMath::TwoPi();
    if (phiHadBlue2 >= TMath::TwoPi()) phiHadBlue2 -= TMath::TwoPi();
    if (phiHadYell2 < 0)               phiHadYell2 += TMath::TwoPi();
    if (phiHadYell2 >= TMath::TwoPi()) phiHadYell2 -= TMath::TwoPi();

    // fill histograms
    hPhiHadBlue2 -> Fill(phiHadBlue2);
    hPhiHadYell2 -> Fill(phiHadYell2);

    // (6) now calculate phiColl: phiSpin - phiHadron,
    //   - constrain to [0, 2pi)
    double phiCollBlue = phiSpinBlue - phiHadBlue;
    double phiCollYell = phiSpinYell - phiHadYell;
    if (phiCollBlue < 0)               phiCollBlue += TMath::TwoPi();
    if (phiCollBlue >= TMath::TwoPi()) phiCollBlue -= TMath::TwoPi();
    if (phiCollYell < 0)               phiCollYell += TMath::TwoPi();
    if (phiCollYell >= TMath::TwoPi()) phiCollYell -= TMath::TwoPi();

    // fill histograms
    hPhiCollBlue -> Fill(phiCollBlue);
    hPhiCollYell -> Fill(phiCollYell);

    // (7) now calculate phiBoer: phiSpin - (2 * phiHadron),
    //   - constrain phiBoerBlue to [0, 2pi)
    double phiBoerBlue = phiSpinBlue - phiHadBlue2;
    double phiBoerYell = phiSpinYell - phiHadYell2;
    if (phiBoerBlue < 0)               phiBoerBlue += TMath::TwoPi();
    if (phiBoerBlue >= TMath::TwoPi()) phiBoerBlue -= TMath::TwoPi();
    if (phiBoerYell < 0)               phiBoerYell += TMath::TwoPi();
    if (phiBoerYell >= TMath::TwoPi()) phiBoerYell -= TMath::TwoPi();

    // fill histograms
    hPhiBoerBlue -> Fill(phiBoerBlue);
    hPhiBoerYell -> Fill(phiBoerYell);

  }  // end iter loop

  /* TODO
   *   - normalize hists
   *   - make plot
   *   - write output
   */
  return;

}  // end AngleCalculationTest

// end ========================================================================
