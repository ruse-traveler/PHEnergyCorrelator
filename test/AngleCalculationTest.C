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
#include <cmath>
#include <iostream>
#include <string>
#include <utility>
// root libraries
#include <TCanvas.h>
#include <TDatime.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TLegend.h>
#include <TMath.h>
#include <TPad.h>
#include <TRandom.h>
#include <TRandom3.h>
#include <TVector3.h>



// ============================================================================
//! Test macro for anglue calculations
// ============================================================================
/*! Note: wrapMode sets what wrapping to do on the
 *  collins angle
 *    0 = wrap to [0, 2pi)
 *    1 = wrap to [0, pi)
 */
void AngleCalculationTest(
  const std::string oFile = "angleCalcTest.nIter10K_wrapCollByPi_noWrap.d23m2y2025.root",
  const std::size_t nIter = 1000,
  const std::size_t wrapMode = 1,
  const bool doWrap = true,
  const bool doBatch = false
) {

  // announce start
  std::cout << "\n  Starting angle calculation test..." << std::endl;

  // create output file
  TFile* fOutput = new TFile(oFile.data(), "recreate");
  if (!fOutput) {
    std::cerr << "PANIC: couldn't create output file!\n" << std::endl;
    return;
  }
  std::cout << "    Created output file." << std::endl;

  // initialize rng
  TDatime*  time  = new TDatime();
  TRandom3* rando = new TRandom3();
  rando -> SetSeed(time -> Get());
  std::cout << "    Initialize RNG." << std::endl;

  // histogram binning
  const int   nAngBins  = 180;
  const int   nCosBins  = 40;
  const int   nXYBins   = 400;
  const int   nZBins    = 400;
  const int   nMagBins  = 100;
  const float xAngStart = -12.60;
  const float xCosStart = -1.0;
  const float xXYStart  = -2.0;
  const float xZStart   = -2.0;
  const float xMagStart = -5.;
  const float xAngStop  = 12.60; 
  const float xCosStop  = 1.0;
  const float xXYStop   = 2.0;
  const float xZStop    = 2.0;
  const float xMagStop  = 5.;

  // turn on errors
  TH1::SetDefaultSumw2(true);
  TH2::SetDefaultSumw2(true);

  // initialize histograms
  TH1D* hInputPhiSpinB     = new TH1D("hInputPhiSpinB", "#phi_{spin}^{B} input", nAngBins, xAngStart, xAngStop);
  TH1D* hInputPhiSpinY     = new TH1D("hInputPhiSpinY", "#phi_{spin}^{Y} input", nAngBins, xAngStart, xAngStop);
  TH2D* hInputXYSpinB      = new TH2D("hInputXYSpinB", "(x,y) sampled (spin B)", nXYBins, xXYStart, xXYStop, nXYBins, xXYStart, xXYStop);
  TH2D* hInputXYSpinY      = new TH2D("hInputXYSpinY", "(x,y) sampled (spin Y)", nXYBins, xXYStart, xXYStop, nXYBins, xXYStart, xXYStop);
  TH1D* hInputPhiJet       = new TH1D("hInputPhiJet", "#phi_{jet} input", nAngBins, xAngStart, xAngStop);
  TH1D* hInputThetaJet     = new TH1D("hInputThetaJet", "#theta_{jet} input", nAngBins, xAngStart, xAngStop);
  TH1D* hInputCosThJet     = new TH1D("hInputCosThJet", "cos#theta_{jet} input", nCosBins, xCosStart, xCosStop);
  TH2D* hInputXYJet        = new TH2D("hInputXYJet", "(x,y) sampled (jet)", nXYBins, xXYStart, xXYStop, nXYBins, xXYStart, xXYStop);
  TH2D* hInputZThJet       = new TH2D("hInputZThJet", "(z,#theta) sampled (jet)", nZBins, xZStart, xZStop, nAngBins, xAngStart, xAngStop);
  TH1D* hInputPhiHad       = new TH1D("hInputPhiHad", "#phi_{h} input", nAngBins, xAngStart, xAngStop);
  TH1D* hInputThetaHad     = new TH1D("hInputThetaHad", "#theta_{h} input", nAngBins, xAngStart, xAngStop);
  TH1D* hInputCosThHad     = new TH1D("hInputCosThHad", "cos#theta_{h} input", nCosBins, xCosStart, xCosStop);
  TH2D* hInputXYHad        = new TH2D("hInputXYHad", "(x,y) sampled (had)", nXYBins, xXYStart, xXYStop, nXYBins, xXYStart, xXYStop);
  TH2D* hInputZThHad       = new TH2D("hInputZThHad", "(z,#theta) sampled (had)", nZBins, xZStart, xZStop, nAngBins, xAngStart, xAngStop);
  TH1D* hCalcPhiJetBeamB   = new TH1D("hCalcPhiJetBeamB", "#phi_{jet-beam}^{B}", nAngBins, xAngStart, xAngStop);
  TH1D* hCalcThetaJetBeamB = new TH1D("hCalcThetaJetBeamB", "#theta_{jet-beam}^{B}", nAngBins, xAngStart, xAngStop);
  TH1D* hCalcCosThJetBeamB = new TH1D("hCalcCosThJetBeamB", "cos#theta_{jet-beam}^{B}", nCosBins, xCosStart, xCosStop);
  TH2D* hCalcZThJetBeamB   = new TH2D("hCalcZThJetBeamB", "(z,#theta) sampled (jet-blue beam)", nZBins, xZStart, xZStop, nAngBins, xAngStart, xAngStop);
  TH1D* hCalcPhiJetBeamY   = new TH1D("hCalcPhiJetBeamY", "#phi_{jet-beam}^{Y}", nAngBins, xAngStart, xAngStop);
  TH1D* hCalcThetaJetBeamY = new TH1D("hCalcThetaJetBeamY", "#theta_{jet-beam}^{Y}", nAngBins, xAngStart, xAngStop);
  TH1D* hCalcCosThJetBeamY = new TH1D("hCalcCosThJetBeamY", "cos#theta_{jet-beam}^{Y}", nCosBins, xCosStart, xCosStop);
  TH2D* hCalcZThJetBeamY   = new TH2D("hCalcZThJetBeamY", "(z,#theta) sampled (jet-yellow beam)", nZBins, xZStart, xZStop, nAngBins, xAngStart, xAngStop);
  TH1D* hCalcPhiHadJet     = new TH1D("hCalcPhiJetHad", "#phi_{jet-h}", nAngBins, xAngStart, xAngStop);
  TH1D* hCalcThetaHadJet   = new TH1D("hCalcThetaJetHad", "#theta_{jet-h}", nAngBins, xAngStart, xAngStop);
  TH1D* hCalcCosThHadJet   = new TH1D("hCalcCosThJetHad", "cos#theta_{jet-h}", nCosBins, xCosStart, xCosStop);
  TH2D* hCalcZThHadJet     = new TH2D("hCalcZThHadJet", "(z,#theta) sampled (jet-hadron)", nZBins, xZStart, xZStop, nAngBins, xAngStart, xAngStop);
  TH1D* hCalcJHBCrossMagB  = new TH1D("hCalcJHBCrossMagB", "|#bf{v}_{jet-beam}^{B} #times #bf{v}_{jet-had}|", nMagBins, xMagStart, xMagStop);
  TH1D* hCalcJHBDotMagB    = new TH1D("hCalcJHBDotMagB", "#bf{v}_{jet-beam}^{B} #upoint #bf{v}_{jet-had}", nMagBins, xMagStart, xMagStop);
  TH1D* hCalcJHBCrossMagY  = new TH1D("hCalcJHBCrossMagY", "|#bf{v}_{jet-beam}^{Y} #times #bf{v}_{jet-had}|", nMagBins, xMagStart, xMagStop);
  TH1D* hCalcJHBDotMagY    = new TH1D("hCalcJHBDotMagY", "#bf{v}_{jet-beam}^{Y} #upoint #bf{v}_{jet-had}", nMagBins, xMagStart, xMagStop);
  TH1D* hPhiSpinB          = new TH1D("hPhiSpinB", "#phi_{spin}^{B}", nAngBins, xAngStart, xAngStop);
  TH1D* hPhiSpinY          = new TH1D("hPhiSpinY", "#phi_{spin}^{Y}", nAngBins, xAngStart, xAngStop);
  TH1D* hPhiHadB           = new TH1D("hPhiHadB", "#phi_{h}^{B}", nAngBins, xAngStart, xAngStop);
  TH1D* hPhiHadY           = new TH1D("hPhiHadY", "#phi_{h}^{Y}", nAngBins, xAngStart, xAngStop);
  TH1D* hPhiHad2B          = new TH1D("hPhiHad2B", "2#phi_{h}^{B}", nAngBins, xAngStart, xAngStop);
  TH1D* hPhiHad2Y          = new TH1D("hPhiHad2Y", "2#phi_{h}^{Y}", nAngBins, xAngStart, xAngStop);
  TH1D* hPhiCollB          = new TH1D("hPhiCollB", "#phi_{collins}^{B}", nAngBins, xAngStart, xAngStop);
  TH1D* hPhiCollY          = new TH1D("hPhiCollY", "#phi_{collins}^{Y}", nAngBins, xAngStart, xAngStop);
  TH1D* hPhiBoerB          = new TH1D("hPhiBoerB", "#phi_{boer}^{B}", nAngBins, xAngStart, xAngStop);
  TH1D* hPhiBoerY          = new TH1D("hPhiBoerY", "#phi_{boer}^{Y}", nAngBins, xAngStart, xAngStop);
  TH1D* hAltPhiSpinB       = new TH1D("hAltPhiSpinB", "#phi_{spin}^{B}", nAngBins, xAngStart, xAngStop);
  TH1D* hAltPhiSpinY       = new TH1D("hAltPhiSpinY", "#phi_{spin}^{Y}", nAngBins, xAngStart, xAngStop);
  TH1D* hAltPhiHadB        = new TH1D("hAltPhiHadB", "#phi_{h}^{B}", nAngBins, xAngStart, xAngStop);
  TH1D* hAltPhiHadY        = new TH1D("hAltPhiHadY", "#phi_{h}^{Y}", nAngBins, xAngStart, xAngStop);
  TH1D* hAltPhiHad2B       = new TH1D("hAltPhiHad2B", "2#phi_{h}^{B}", nAngBins, xAngStart, xAngStop);
  TH1D* hAltPhiHad2Y       = new TH1D("hAltPhiHad2Y", "2#phi_{h}^{Y}", nAngBins, xAngStart, xAngStop);
  TH1D* hAltPhiCollB       = new TH1D("hAltPhiCollB", "#phi_{collins}^{B}", nAngBins, xAngStart, xAngStop);
  TH1D* hAltPhiCollY       = new TH1D("hAltPhiCollY", "#phi_{collins}^{Y}", nAngBins, xAngStart, xAngStop);
  TH1D* hAltPhiBoerB       = new TH1D("hAltPhiBoerB", "#phi_{boer}^{B}", nAngBins, xAngStart, xAngStop);
  TH1D* hAltPhiBoerY       = new TH1D("hAltPhiBoerY", "#phi_{boer}^{Y}", nAngBins, xAngStart, xAngStop);
  std::cout << "    Created histograms.\n"
            << "    MC loop: running " << nIter << " iterations:"
            << std::endl;

  // --------------------------------------------------------------------------
  // run calculation nIter times
  // --------------------------------------------------------------------------
  for (std::size_t iIter = 0; iIter < nIter; ++iIter) {

    // announce progress
    std::size_t iProg = iIter + 1;
    if (iProg == nIter) {
      std::cout << "      iter " << iProg << "/" << nIter << "..." << std::endl;
    } else if (doBatch) {
      std::cout << "      iter " << iProg << "/" << nIter << "..." << std::endl;
    } else {
      std::cout << "      iter " << iProg << "/" << nIter << "...\r" << std::flush;
    }

    // set up vectors ---------------------------------------------------------

    // set beam directions
    TVector3 vecBeamB3(0.0, 0.0, 1.0);
    TVector3 vecBeamY3(0.0, 0.0, -1.0);

    // get random phi = (0, 2pi) for spins
    const double phiRandSpinB = rando -> Uniform(0.0, TMath::TwoPi());
    const double phiRandSpinY = rando -> Uniform(0.0, TMath::TwoPi());

    // now translate into (x, y) values (assuming
    // r = 1.0)
    const double xRandSpinB = cos(phiRandSpinB);
    const double yRandSpinB = sin(phiRandSpinB);
    const double xRandSpinY = cos(phiRandSpinY);
    const double yRandSpinY = sin(phiRandSpinY);

    // set spin directions, fill input histograms
    TVector3 vecSpinB3(xRandSpinB, yRandSpinB, 0.0);
    TVector3 vecSpinY3(xRandSpinY, yRandSpinY, 0.0);
    hInputPhiSpinB -> Fill( vecSpinB3.Phi() );
    hInputPhiSpinY -> Fill( vecSpinY3.Phi() );
    hInputXYSpinB  -> Fill( vecSpinB3.X(), vecSpinB3.Y() );
    hInputXYSpinY  -> Fill( vecSpinY3.X(), vecSpinY3.Y() );

    // get random phi = (0, 2pi), cos-theta = (-1, 1), and u
    // for jet/hadron
    const double phiRandJet = rando -> Uniform(0.0, TMath::TwoPi());
    const double phiRandHad = rando -> Uniform(0.0, TMath::TwoPi());
    const double cosRandJet = rando -> Uniform(-1.0, 1.0);
    const double cosRandHad = rando -> Uniform(-1.0, 1.0);
    const double uRandJet   = rando -> Uniform(0.0, 1.0);
    const double uRandHad   = rando -> Uniform(0.0, 1.0);

    // translate cos-theta, u into theta, r
    const double thetaRandJet = acos(cosRandJet);
    const double thetaRandHad = acos(cosRandHad);
    const double rRandJet     = pow(uRandJet, 1.0/3.0);
    const double rRandHad     = pow(uRandHad, 1.0/3.0);

    // now translate into (x, y, z) values (assuming
    // r = 1.0)
    const double xRandJet = rRandJet * sin(thetaRandJet) * cos(phiRandJet);
    const double yRandJet = rRandJet * sin(thetaRandJet) * sin(phiRandJet);
    const double zRandJet = rRandJet * cos(thetaRandJet);
    const double xRandHad = rRandHad * sin(thetaRandHad) * cos(phiRandHad);
    const double yRandHad = rRandHad * sin(thetaRandHad) * sin(phiRandHad);
    const double zRandHad = rRandHad * cos(thetaRandHad);

    // get random jet/hadron directions, fill input histograms
    TVector3 vecJet3(xRandJet, yRandJet, zRandJet);
    TVector3 vecHad3(xRandHad, yRandHad, zRandHad);
    hInputPhiJet   -> Fill( vecJet3.Phi() );
    hInputThetaJet -> Fill( vecJet3.Theta() );
    hInputCosThJet -> Fill( cos(vecJet3.Theta()) );
    hInputXYJet    -> Fill( vecJet3.X(), vecJet3.Y() );
    hInputZThJet   -> Fill( vecJet3.Z(), vecJet3.Theta() );
    hInputPhiHad   -> Fill( vecHad3.Phi() );
    hInputThetaHad -> Fill( vecHad3.Theta() );
    hInputCosThHad -> Fill( cos(vecHad3.Theta()) );
    hInputXYHad    -> Fill( vecHad3.X(), vecHad3.Y() );
    hInputZThHad   -> Fill( vecHad3.Z(), vecHad3.Theta() );

    // now normalize spin/jet/hadron vectors
    TVector3 unitSpinB3 = vecSpinB3.Unit();
    TVector3 unitSpinY3 = vecSpinY3.Unit();
    TVector3 unitJet3   = vecJet3.Unit();
    TVector3 unitHad3   = vecHad3.Unit();

    // collins & boer-mulders angle calculations --------------------------


    // (1) get vectors normal to the jet-beam plane,
    //     fill intermediate histograms
    std::pair<TVector3, TVector3> normJetBeam3 = std::make_pair(
      ( vecBeamB3.Cross(unitJet3) ).Unit(),
      ( vecBeamY3.Cross(unitJet3) ).Unit()
    );
    hCalcPhiJetBeamB   -> Fill( normJetBeam3.first.Phi() );
    hCalcThetaJetBeamB -> Fill( normJetBeam3.first.Theta() );
    hCalcCosThJetBeamB -> Fill( cos(normJetBeam3.first.Theta()) );
    hCalcZThJetBeamB   -> Fill( normJetBeam3.first.Z(), normJetBeam3.first.Theta() );
    hCalcPhiJetBeamY   -> Fill( normJetBeam3.second.Phi() );
    hCalcThetaJetBeamY -> Fill( normJetBeam3.second.Theta() );
    hCalcCosThJetBeamY -> Fill( cos(normJetBeam3.second.Theta()) );
    hCalcZThJetBeamY   -> Fill( normJetBeam3.second.Z(), normJetBeam3.second.Theta() );

    // (2) get phiSpin: angles between the jet-beam plane and spin
    //   - n.b. for spin pattern >= 4, the yellow spin is randomized
    //   - angle between jet plane and spin 
    //   - note that we get the full [0,2pi) range for these angles
    double phiSpinBlue = TMath::PiOver2() - atan2( normJetBeam3.first.Cross(unitSpinB3).Mag(), normJetBeam3.first.Dot(unitSpinB3) );
    double phiSpinYell = TMath::PiOver2() - atan2( normJetBeam3.second.Cross(unitSpinY3).Mag(), normJetBeam3.second.Dot(unitSpinY3) );
    if (doWrap) {
      if (phiSpinBlue < 0)               phiSpinBlue += TMath::TwoPi();
      if (phiSpinBlue >= TMath::TwoPi()) phiSpinBlue -= TMath::TwoPi();
      if (phiSpinYell < 0)               phiSpinYell += TMath::TwoPi();
      if (phiSpinYell >= TMath::TwoPi()) phiSpinYell -= TMath::TwoPi();
    }

    // fill spin histograms
    hPhiSpinB -> Fill(phiSpinBlue);
    hPhiSpinY -> Fill(phiSpinYell);

    // (3) get vector normal to hadron average-jet plane,
    //     fill intermediate histograms
    TVector3 normHadJet3 = ( unitJet3.Cross(unitHad3) ).Unit();
    hCalcPhiHadJet   -> Fill( normHadJet3.Phi() );
    hCalcThetaHadJet -> Fill( normHadJet3.Theta() );
    hCalcCosThHadJet -> Fill( cos(normHadJet3.Theta()) );
    hCalcZThHadJet   -> Fill( normHadJet3.Z(), normHadJet3.Theta() );

    // cross/dot jet-beam and jet-hadron normals
    double jhbCrossB = normJetBeam3.first.Cross(normHadJet3).Mag();
    double jhbDotB   = normJetBeam3.first.Dot(normHadJet3);
    double jhbCrossY = normJetBeam3.second.Cross(normHadJet3).Mag();
    double jhbDotY   = normJetBeam3.second.Dot(normHadJet3);

    // fill intermediate histograms
    hCalcJHBCrossMagB -> Fill(jhbCrossB); 
    hCalcJHBDotMagB   -> Fill(jhbDotB);
    hCalcJHBCrossMagY -> Fill(jhbCrossY);
    hCalcJHBDotMagY   -> Fill(jhbDotY);

    // (4) get phiHadron: angle between the jet-beam plane and the
    //   - angle between jet-hadron plane
    //   - constrain to range [0,2pi)
    double phiHadBlue = atan2( jhbCrossB, jhbDotB );
    double phiHadYell = atan2( jhbCrossY, jhbDotY );
    if (doWrap) {
      if (phiHadBlue < 0)               phiHadBlue += TMath::TwoPi();
      if (phiHadBlue >= TMath::TwoPi()) phiHadBlue -= TMath::TwoPi();
      if (phiHadYell < 0)               phiHadYell += TMath::TwoPi();
      if (phiHadYell >= TMath::TwoPi()) phiHadYell -= TMath::TwoPi();
    }

    // fill histograms
    hPhiHadB -> Fill(phiHadBlue);
    hPhiHadY -> Fill(phiHadYell);

    // (5) double phiHadron for boer-mulders,
    //   - constrain to [0, 2pi)
    double phiHadBlue2 = 2.0 * phiHadBlue;
    double phiHadYell2 = 2.0 * phiHadYell;
    if (doWrap) {
      if (phiHadBlue2 < 0)               phiHadBlue2 += TMath::TwoPi();
      if (phiHadBlue2 >= TMath::TwoPi()) phiHadBlue2 -= TMath::TwoPi();
      if (phiHadYell2 < 0)               phiHadYell2 += TMath::TwoPi();
      if (phiHadYell2 >= TMath::TwoPi()) phiHadYell2 -= TMath::TwoPi();
    }

    // fill histograms
    hPhiHad2B -> Fill(phiHadBlue2);
    hPhiHad2Y -> Fill(phiHadYell2);

    // (6) now calculate phiColl: phiSpin - phiHadron,
    //   - constrain to [0, 2pi) OR [0, pi)
    double phiCollBlue = phiSpinBlue - phiHadBlue;
    double phiCollYell = phiSpinYell - phiHadYell;
    if (doWrap) {
      if (wrapMode == 1) {
        if (phiCollBlue < 0)            phiCollBlue += TMath::Pi();
        if (phiCollBlue >= TMath::Pi()) phiCollBlue -= TMath::Pi();
        if (phiCollYell < 0)            phiCollYell += TMath::Pi();
        if (phiCollYell >= TMath::Pi()) phiCollYell -= TMath::Pi();
      } else {
        if (phiCollBlue < 0)               phiCollBlue += TMath::TwoPi();
        if (phiCollBlue >= TMath::TwoPi()) phiCollBlue -= TMath::TwoPi();
        if (phiCollYell < 0)               phiCollYell += TMath::TwoPi();
        if (phiCollYell >= TMath::TwoPi()) phiCollYell -= TMath::TwoPi();
      }
    }

    // fill histograms
    hPhiCollB -> Fill(phiCollBlue);
    hPhiCollY -> Fill(phiCollYell);

    // (7) now calculate phiBoer: phiSpin - (2 * phiHadron),
    //   - constrain phiBoerBlue to [0, 2pi)
    double phiBoerBlue = phiSpinBlue - phiHadBlue2;
    double phiBoerYell = phiSpinYell - phiHadYell2;
    if (doWrap) {
      if (phiBoerBlue < 0)               phiBoerBlue += TMath::TwoPi();
      if (phiBoerBlue >= TMath::TwoPi()) phiBoerBlue -= TMath::TwoPi();
      if (phiBoerYell < 0)               phiBoerYell += TMath::TwoPi();
      if (phiBoerYell >= TMath::TwoPi()) phiBoerYell -= TMath::TwoPi();
    }

    // fill histograms
    hPhiBoerB -> Fill(phiBoerBlue);
    hPhiBoerY -> Fill(phiBoerYell);

    // alternate calc: what if we used acos to get the angles? ----------------

    // (2) starting at get phiSpin: angles between the jet-beam plane and spin
    //   - n.b. for spin pattern >= 4, the yellow spin is randomized
    //   - angle between jet plane and spin
    //   - note that we get the full [0,2pi) range for these angles
    double phiSpinBlueAlt = TMath::PiOver2() - acos( normJetBeam3.first.Dot(unitSpinB3) / (normJetBeam3.first.Mag() * unitSpinB3.Mag()) );
    double phiSpinYellAlt = TMath::PiOver2() - acos( normJetBeam3.second.Dot(unitSpinY3) / (normJetBeam3.second.Mag() * unitSpinY3.Mag()) );
    if (doWrap) {
      if (phiSpinBlueAlt < 0)               phiSpinBlueAlt += TMath::TwoPi();
      if (phiSpinBlueAlt >= TMath::TwoPi()) phiSpinBlueAlt -= TMath::TwoPi();
      if (phiSpinYellAlt < 0)               phiSpinYellAlt += TMath::TwoPi();
      if (phiSpinYellAlt >= TMath::TwoPi()) phiSpinYellAlt -= TMath::TwoPi();
    }

    // fill spin histograms
    hAltPhiSpinB -> Fill(phiSpinBlueAlt);
    hAltPhiSpinY -> Fill(phiSpinYellAlt);

    // (4) now jump to phiHadron: angle between the jet-beam plane and the
    //     jet-hadron plane
    //   - angle between jet-hadron plane
    //   - constrain to range [0,2pi)
    double phiHadBlueAlt = acos( normJetBeam3.first.Dot(normHadJet3) / (normJetBeam3.first.Mag() * normHadJet3.Mag()) );
    double phiHadYellAlt = acos( normJetBeam3.second.Dot(normHadJet3) /  (normJetBeam3.second.Mag() * normHadJet3.Mag()) );
    if (doWrap) {
      if (phiHadBlueAlt < 0)               phiHadBlueAlt += TMath::TwoPi();
      if (phiHadBlueAlt >= TMath::TwoPi()) phiHadBlueAlt -= TMath::TwoPi();
      if (phiHadYellAlt < 0)               phiHadYellAlt += TMath::TwoPi();
      if (phiHadYellAlt >= TMath::TwoPi()) phiHadYellAlt -= TMath::TwoPi();
    }

    // fill histograms
    hAltPhiHadB -> Fill(phiHadBlueAlt);
    hAltPhiHadY -> Fill(phiHadYellAlt);

    // (5) double phiHadron for boer-mulders,
    //   - constrain to [0, 2pi)
    double phiHadBlueAlt2 = 2.0 * phiHadBlueAlt;
    double phiHadYellAlt2 = 2.0 * phiHadYellAlt;
    if (doWrap) {
      if (phiHadBlueAlt2 < 0)               phiHadBlueAlt2 += TMath::TwoPi();
      if (phiHadBlueAlt2 >= TMath::TwoPi()) phiHadBlueAlt2 -= TMath::TwoPi();
      if (phiHadYellAlt2 < 0)               phiHadYellAlt2 += TMath::TwoPi();
      if (phiHadYellAlt2 >= TMath::TwoPi()) phiHadYellAlt2 -= TMath::TwoPi();
    }

    // fill histograms
    hAltPhiHad2B -> Fill(phiHadBlueAlt2);
    hAltPhiHad2Y -> Fill(phiHadYellAlt2);

    // (6) now calculate phiColl: phiSpin - phiHadron,
    //   - constrain to [0, 2pi) OR [0, pi)
    double phiCollBlueAlt = phiSpinBlueAlt - phiHadBlueAlt;
    double phiCollYellAlt = phiSpinYellAlt - phiHadYellAlt;
    if (doWrap) {
      if (wrapMode == 1) {
        if (phiCollBlueAlt < 0)            phiCollBlueAlt += TMath::Pi();
        if (phiCollBlueAlt >= TMath::Pi()) phiCollBlueAlt -= TMath::Pi();
        if (phiCollYellAlt < 0)            phiCollYellAlt += TMath::Pi();
        if (phiCollYellAlt >= TMath::Pi()) phiCollYellAlt -= TMath::Pi();
      } else {
        if (phiCollBlueAlt < 0)               phiCollBlueAlt += TMath::TwoPi();
        if (phiCollBlueAlt >= TMath::TwoPi()) phiCollBlueAlt -= TMath::TwoPi();
        if (phiCollYellAlt < 0)               phiCollYellAlt += TMath::TwoPi();
        if (phiCollYellAlt >= TMath::TwoPi()) phiCollYellAlt -= TMath::TwoPi();
      }
    }

    // fill histograms
    hAltPhiCollB -> Fill(phiCollBlueAlt);
    hAltPhiCollY -> Fill(phiCollYellAlt);

    // (7) now calculate phiBoer: phiSpin - (2 * phiHadron),
    //   - constrain phiBoerBlue to [0, 2pi)
    double phiBoerBlueAlt = phiSpinBlueAlt - phiHadBlueAlt2;
    double phiBoerYellAlt = phiSpinYellAlt - phiHadYellAlt2;
    if (doWrap) {
      if (phiBoerBlueAlt < 0)               phiBoerBlueAlt += TMath::TwoPi();
      if (phiBoerBlueAlt >= TMath::TwoPi()) phiBoerBlueAlt -= TMath::TwoPi();
      if (phiBoerYellAlt < 0)               phiBoerYellAlt += TMath::TwoPi();
      if (phiBoerYellAlt >= TMath::TwoPi()) phiBoerYellAlt -= TMath::TwoPi();
    }

    // fill histograms
    hAltPhiBoerB -> Fill(phiBoerBlueAlt);
    hAltPhiBoerY -> Fill(phiBoerYellAlt);

  }  // end iter loop
  std::cout << "    MC loop finished!" << std::endl;

  // normalize histograms
  //   - n.b. (x,y), (z,theta), and magnitude
  //     hists NOT normalized
  hInputPhiSpinB     -> Scale(1. / (double) nIter);
  hInputPhiSpinY     -> Scale(1. / (double) nIter);
  hInputPhiJet       -> Scale(1. / (double) nIter);
  hInputThetaJet     -> Scale(1. / (double) nIter);
  hInputCosThJet     -> Scale(1. / (double) nIter);
  hInputPhiHad       -> Scale(1. / (double) nIter);
  hInputThetaHad     -> Scale(1. / (double) nIter);
  hInputCosThHad     -> Scale(1. / (double) nIter);
  hCalcPhiJetBeamB   -> Scale(1. / (double) nIter);
  hCalcThetaJetBeamB -> Scale(1. / (double) nIter);
  hCalcCosThJetBeamB -> Scale(1. / (double) nIter);
  hCalcPhiJetBeamY   -> Scale(1. / (double) nIter);
  hCalcThetaJetBeamY -> Scale(1. / (double) nIter);
  hCalcCosThJetBeamY -> Scale(1. / (double) nIter);
  hCalcPhiHadJet     -> Scale(1. / (double) nIter);
  hCalcThetaHadJet   -> Scale(1. / (double) nIter);
  hCalcCosThHadJet   -> Scale(1. / (double) nIter);
  hPhiSpinB          -> Scale(1. / (double) nIter);
  hPhiSpinY          -> Scale(1. / (double) nIter);
  hPhiHadB           -> Scale(1. / (double) nIter);
  hPhiHadY           -> Scale(1. / (double) nIter);
  hPhiHad2B          -> Scale(1. / (double) nIter);
  hPhiHad2Y          -> Scale(1. / (double) nIter);
  hPhiCollB          -> Scale(1. / (double) nIter);
  hPhiCollY          -> Scale(1. / (double) nIter);
  hPhiBoerB          -> Scale(1. / (double) nIter);
  hPhiBoerY          -> Scale(1. / (double) nIter);
  hAltPhiSpinB       -> Scale(1. / (double) nIter);
  hAltPhiSpinY       -> Scale(1. / (double) nIter);
  hAltPhiHadB        -> Scale(1. / (double) nIter);
  hAltPhiHadY        -> Scale(1. / (double) nIter);
  hAltPhiHad2B       -> Scale(1. / (double) nIter);
  hAltPhiHad2Y       -> Scale(1. / (double) nIter);
  hAltPhiCollB       -> Scale(1. / (double) nIter);
  hAltPhiCollY       -> Scale(1. / (double) nIter);
  hAltPhiBoerB       -> Scale(1. / (double) nIter);
  hAltPhiBoerY       -> Scale(1. / (double) nIter);
  std::cout << "    Normalized histograms." << std::endl;

  // create frame histograms
  TH1D* hPhiFrame    = hInputPhiSpinB -> Clone();
  TH1D* hAltPhiFrame = hInputPhiSpinB -> Clone();
  TH1D* hThetaFrame  = hInputThetaJet -> Clone();
  TH1D* hCosThFrame  = hInputCosThJet -> Clone();
  hPhiFrame    -> Reset("ICES");
  hPhiFrame    -> SetName("hPhiFrame");
  hPhiFrame    -> SetTitle(";#phi [rad]");
  hAltPhiFrame -> Reset("ICES");
  hAltPhiFrame -> SetName("hAltPhiFrame");
  hAltPhiFrame -> SetTitle("Calculated using cos#phi = A#dotB/|A||B|;#phi [rad]");
  hThetaFrame  -> Reset("ICES");
  hThetaFrame  -> SetName("hThetaFrame");
  hThetaFrame  -> SetTitle(";#theta [rad]");
  hCosThFrame  -> Reset("ICES");
  hCosThFrame  -> SetName("hCosThFrame");
  hCosThFrame  -> SetTitle(";cos#theta");

  // create header
  TString sHeader("#bf{");
  sHeader += nIter;
  sHeader += " iterations}";

  // styles for histograms
  const int col[10] = {
    923,
    921,
    799,
    899,
    617,
    879,
    859,
    839,
    819,
    402
  };
  const int mar[10] = {
    29,
    24,
    25,
    26,
    32,
    27,
    28,
    46,
    30,
    24
  };

  // set styles
  hInputPhiSpinB     -> SetLineColor(col[0]);
  hInputPhiSpinB     -> SetMarkerColor(col[0]);
  hInputPhiSpinB     -> SetMarkerStyle(mar[0]);
  hInputPhiSpinY     -> SetLineColor(col[0]);
  hInputPhiSpinY     -> SetMarkerColor(col[0]);
  hInputPhiSpinY     -> SetMarkerStyle(mar[0]);
  hInputPhiJet       -> SetLineColor(col[1]);
  hInputPhiJet       -> SetMarkerColor(col[1]);
  hInputPhiJet       -> SetMarkerStyle(mar[1]);
  hInputThetaJet     -> SetLineColor(col[1]);
  hInputThetaJet     -> SetMarkerColor(col[1]);
  hInputThetaJet     -> SetMarkerStyle(mar[1]);
  hInputCosThJet     -> SetLineColor(col[1]);
  hInputCosThJet     -> SetMarkerColor(col[1]);
  hInputCosThJet     -> SetMarkerStyle(mar[1]);
  hInputPhiHad       -> SetLineColor(col[2]);
  hInputPhiHad       -> SetMarkerColor(col[2]);
  hInputPhiHad       -> SetMarkerStyle(mar[2]);
  hInputThetaHad     -> SetLineColor(col[2]);
  hInputThetaHad     -> SetMarkerColor(col[2]);
  hInputThetaHad     -> SetMarkerStyle(mar[2]);
  hInputCosThHad     -> SetLineColor(col[2]);
  hInputCosThHad     -> SetMarkerColor(col[2]);
  hInputCosThHad     -> SetMarkerStyle(mar[2]);
  hCalcPhiJetBeamB   -> SetLineColor(col[3]);
  hCalcPhiJetBeamB   -> SetMarkerColor(col[3]);
  hCalcPhiJetBeamB   -> SetMarkerStyle(mar[3]);
  hCalcThetaJetBeamB -> SetLineColor(col[3]);
  hCalcThetaJetBeamB -> SetMarkerColor(col[3]);
  hCalcThetaJetBeamB -> SetMarkerStyle(mar[3]);
  hCalcCosThJetBeamB -> SetLineColor(col[3]);
  hCalcCosThJetBeamB -> SetMarkerColor(col[3]);
  hCalcCosThJetBeamB -> SetMarkerStyle(mar[3]);
  hCalcPhiJetBeamY   -> SetLineColor(col[3]);
  hCalcPhiJetBeamY   -> SetMarkerColor(col[3]);
  hCalcPhiJetBeamY   -> SetMarkerStyle(mar[3]);
  hCalcThetaJetBeamY -> SetLineColor(col[3]);
  hCalcThetaJetBeamY -> SetMarkerColor(col[3]);
  hCalcThetaJetBeamY -> SetMarkerStyle(mar[3]);
  hCalcCosThJetBeamY -> SetLineColor(col[3]);
  hCalcCosThJetBeamY -> SetMarkerColor(col[3]);
  hCalcCosThJetBeamY -> SetMarkerStyle(mar[3]);
  hCalcPhiHadJet     -> SetLineColor(col[4]);
  hCalcPhiHadJet     -> SetMarkerColor(col[4]);
  hCalcPhiHadJet     -> SetMarkerStyle(mar[4]);
  hCalcThetaHadJet   -> SetLineColor(col[4]);
  hCalcThetaHadJet   -> SetMarkerColor(col[4]);
  hCalcThetaHadJet   -> SetMarkerStyle(mar[4]);
  hCalcCosThHadJet   -> SetLineColor(col[4]);
  hCalcCosThHadJet   -> SetMarkerColor(col[4]);
  hCalcCosThHadJet   -> SetMarkerStyle(mar[4]);
  hPhiSpinB          -> SetLineColor(col[5]);
  hPhiSpinB          -> SetMarkerColor(col[5]);
  hPhiSpinB          -> SetMarkerStyle(mar[5]);
  hPhiSpinY          -> SetLineColor(col[5]);
  hPhiSpinY          -> SetMarkerColor(col[5]);
  hPhiSpinY          -> SetMarkerStyle(mar[5]);
  hAltPhiSpinB       -> SetLineColor(col[5]);
  hAltPhiSpinB       -> SetMarkerColor(col[5]);
  hAltPhiSpinB       -> SetMarkerStyle(mar[5]);
  hAltPhiSpinY       -> SetLineColor(col[5]);
  hAltPhiSpinY       -> SetMarkerColor(col[5]);
  hAltPhiSpinY       -> SetMarkerStyle(mar[5]);
  hPhiHadB           -> SetLineColor(col[6]);
  hPhiHadB           -> SetMarkerColor(col[6]);
  hPhiHadB           -> SetMarkerStyle(mar[6]);
  hPhiHadY           -> SetLineColor(col[6]);
  hPhiHadY           -> SetMarkerColor(col[6]);
  hPhiHadY           -> SetMarkerStyle(mar[6]);
  hAltPhiHadB        -> SetLineColor(col[6]);
  hAltPhiHadB        -> SetMarkerColor(col[6]);
  hAltPhiHadB        -> SetMarkerStyle(mar[6]);
  hAltPhiHadY        -> SetLineColor(col[6]);
  hAltPhiHadY        -> SetMarkerColor(col[6]);
  hAltPhiHadY        -> SetMarkerStyle(mar[6]);
  hPhiHad2B          -> SetLineColor(col[7]);
  hPhiHad2B          -> SetMarkerColor(col[7]);
  hPhiHad2B          -> SetMarkerStyle(mar[7]);
  hPhiHad2Y          -> SetLineColor(col[7]);
  hPhiHad2Y          -> SetMarkerColor(col[7]);
  hPhiHad2Y          -> SetMarkerStyle(mar[7]);
  hAltPhiHad2B       -> SetLineColor(col[7]);
  hAltPhiHad2B       -> SetMarkerColor(col[7]);
  hAltPhiHad2B       -> SetMarkerStyle(mar[7]);
  hAltPhiHad2Y       -> SetLineColor(col[7]);
  hAltPhiHad2Y       -> SetMarkerColor(col[7]);
  hAltPhiHad2Y       -> SetMarkerStyle(mar[7]);
  hPhiCollB          -> SetLineColor(col[8]);
  hPhiCollB          -> SetMarkerColor(col[8]);
  hPhiCollB          -> SetMarkerStyle(mar[8]);
  hPhiCollY          -> SetLineColor(col[8]);
  hPhiCollY          -> SetMarkerColor(col[8]);
  hPhiCollY          -> SetMarkerStyle(mar[8]);
  hAltPhiCollB       -> SetLineColor(col[8]);
  hAltPhiCollB       -> SetMarkerColor(col[8]);
  hAltPhiCollB       -> SetMarkerStyle(mar[8]);
  hAltPhiCollY       -> SetLineColor(col[8]);
  hAltPhiCollY       -> SetMarkerColor(col[8]);
  hAltPhiCollY       -> SetMarkerStyle(mar[8]);
  hPhiBoerB          -> SetLineColor(col[9]);
  hPhiBoerB          -> SetMarkerColor(col[9]);
  hPhiBoerB          -> SetMarkerStyle(mar[9]);
  hPhiBoerY          -> SetLineColor(col[9]);
  hPhiBoerY          -> SetMarkerColor(col[9]);
  hPhiBoerY          -> SetMarkerStyle(mar[9]);
  hAltPhiBoerB       -> SetLineColor(col[9]);
  hAltPhiBoerB       -> SetMarkerColor(col[9]);
  hAltPhiBoerB       -> SetMarkerStyle(mar[9]);
  hAltPhiBoerY       -> SetLineColor(col[9]);
  hAltPhiBoerY       -> SetMarkerColor(col[9]);
  hAltPhiBoerY       -> SetMarkerStyle(mar[9]);
  std::cout << "    Set styles." << std::endl;

  // create phi input plots
  {

    // create legend for blue phi
    TLegend* lInputPhiB = new TLegend(0.1, 0.1, 0.3, 0.4, sHeader.Data());
    lInputPhiB -> SetFillColor(0);
    lInputPhiB -> SetLineColor(0);
    lInputPhiB -> SetTextFont(42);
    lInputPhiB -> SetTextAlign(12);
    lInputPhiB -> AddEntry(hInputPhiSpinB, hInputPhiSpinB -> GetTitle(), "P");
    lInputPhiB -> AddEntry(hInputPhiJet, hInputPhiJet -> GetTitle(), "P");
    lInputPhiB -> AddEntry(hInputPhiHad, hInputPhiHad -> GetTitle(), "P");
    lInputPhiB -> AddEntry(hCalcPhiJetBeamB, hCalcPhiJetBeamB -> GetTitle(), "P");
    lInputPhiB -> AddEntry(hCalcPhiHadJet, hCalcPhiHadJet -> GetTitle(), "P");

    // create plot for blue phi
    TCanvas* cInputPhiB = new TCanvas("cInputPhiBlue", "", 750, 750);
    cInputPhiB       -> SetGrid(0, 0);
    cInputPhiB       -> cd();
    hPhiFrame        -> Draw();
    hInputPhiSpinB   -> Draw("same");
    hInputPhiJet     -> Draw("same");
    hInputPhiHad     -> Draw("same");
    hCalcPhiJetBeamB -> Draw("same");
    hCalcPhiHadJet   -> Draw("same");
    lInputPhiB       -> Draw();
    fOutput          -> cd();
    cInputPhiB       -> Write();
    cInputPhiB       -> Close();

    // create legend for yellow phi
    TLegend* lInputPhiY = new TLegend(0.1, 0.1, 0.3, 0.4, sHeader.Data());
    lInputPhiY -> SetFillColor(0);
    lInputPhiY -> SetLineColor(0);
    lInputPhiY -> SetTextFont(42);
    lInputPhiY -> SetTextAlign(12);
    lInputPhiY -> AddEntry(hInputPhiSpinY, hInputPhiSpinY -> GetTitle(), "P");
    lInputPhiY -> AddEntry(hInputPhiJet, hInputPhiJet -> GetTitle(), "P");
    lInputPhiY -> AddEntry(hInputPhiHad, hInputPhiHad -> GetTitle(), "P");
    lInputPhiY -> AddEntry(hCalcPhiJetBeamY, hCalcPhiJetBeamY -> GetTitle(), "P");
    lInputPhiY -> AddEntry(hCalcPhiHadJet, hCalcPhiHadJet -> GetTitle(), "P");

    // create plot for yellow phi
    TCanvas* cInputPhiY = new TCanvas("cInputPhiYellow", "", 750, 750);
    cInputPhiY       -> SetGrid(0, 0);
    cInputPhiY       -> cd();
    hPhiFrame        -> Draw();
    hInputPhiSpinY   -> Draw("same");
    hInputPhiJet     -> Draw("same");
    hInputPhiHad     -> Draw("same");
    hCalcPhiJetBeamY -> Draw("same");
    hCalcPhiHadJet   -> Draw("same");
    lInputPhiY       -> Draw();
    fOutput          -> cd();
    cInputPhiY       -> Write();
    cInputPhiY       -> Close();

  }  // end input phi plot making
  std::cout << "    Created phi input plots." << std::endl;

  // create theta input plots
  {

    // create legend for blue theta
    TLegend* lInputThetaB = new TLegend(0.1, 0.1, 0.3, 0.35, sHeader.Data());
    lInputThetaB -> SetFillColor(0);
    lInputThetaB -> SetLineColor(0);
    lInputThetaB -> SetTextFont(42);
    lInputThetaB -> SetTextAlign(12);
    lInputThetaB -> AddEntry(hInputThetaJet, hInputThetaJet -> GetTitle(), "P");
    lInputThetaB -> AddEntry(hInputThetaHad, hInputThetaHad -> GetTitle(), "P");
    lInputThetaB -> AddEntry(hCalcThetaJetBeamB, hCalcThetaJetBeamB -> GetTitle(), "P");
    lInputThetaB -> AddEntry(hCalcThetaHadJet, hCalcThetaHadJet -> GetTitle(), "P");

    // create plot for blue theta
    TCanvas* cInputThetaB = new TCanvas("cInputThetaBlue", "", 750, 750);
    cInputThetaB       -> SetGrid(0, 0);
    cInputThetaB       -> cd();
    hThetaFrame        -> Draw();
    hInputThetaJet     -> Draw("same");
    hInputThetaHad     -> Draw("same");
    hCalcThetaJetBeamB -> Draw("same");
    hCalcThetaHadJet   -> Draw("same");
    lInputThetaB       -> Draw();
    fOutput            -> cd();
    cInputThetaB       -> Write();
    cInputThetaB       -> Close();

    // create legend for yellow theta
    TLegend* lInputThetaY = new TLegend(0.1, 0.1, 0.3, 0.35, sHeader.Data());
    lInputThetaY -> SetFillColor(0);
    lInputThetaY -> SetLineColor(0);
    lInputThetaY -> SetTextFont(42);
    lInputThetaY -> SetTextAlign(12);
    lInputThetaY -> AddEntry(hInputThetaJet, hInputThetaJet -> GetTitle(), "P");
    lInputThetaY -> AddEntry(hInputThetaHad, hInputThetaHad -> GetTitle(), "P");
    lInputThetaY -> AddEntry(hCalcThetaJetBeamY, hCalcThetaJetBeamY -> GetTitle(), "P");
    lInputThetaY -> AddEntry(hCalcThetaHadJet, hCalcThetaHadJet -> GetTitle(), "P");

    // create plot for blue theta
    TCanvas* cInputThetaY = new TCanvas("cInputThetaYellow", "", 750, 750);
    cInputThetaY       -> SetGrid(0, 0);
    cInputThetaY       -> cd();
    hThetaFrame        -> Draw();
    hInputThetaJet     -> Draw("same");
    hInputThetaHad     -> Draw("same");
    hCalcThetaJetBeamY -> Draw("same");
    hCalcThetaHadJet   -> Draw("same");
    lInputThetaY       -> Draw();
    fOutput            -> cd();
    cInputThetaY       -> Write();
    cInputThetaY       -> Close();

    // create legend for blue cos-th
    TLegend* lInputCosThB = new TLegend(0.1, 0.1, 0.3, 0.35, sHeader.Data());
    lInputCosThB -> SetFillColor(0);
    lInputCosThB -> SetLineColor(0);
    lInputCosThB -> SetTextFont(42);
    lInputCosThB -> SetTextAlign(12);
    lInputCosThB -> AddEntry(hInputCosThJet, hInputCosThJet -> GetTitle(), "P");
    lInputCosThB -> AddEntry(hInputCosThHad, hInputCosThHad -> GetTitle(), "P");
    lInputCosThB -> AddEntry(hCalcCosThJetBeamB, hCalcCosThJetBeamB -> GetTitle(), "P");
    lInputCosThB -> AddEntry(hCalcCosThHadJet, hCalcCosThHadJet -> GetTitle(), "P");

    // create plot for blue cos-th
    TCanvas* cInputCosThB = new TCanvas("cInputCosThBlue", "", 750, 750);
    cInputCosThB       -> SetGrid(0, 0);
    cInputCosThB       -> cd();
    hCosThFrame        -> Draw();
    hInputCosThJet     -> Draw("same");
    hInputCosThHad     -> Draw("same");
    hCalcCosThJetBeamB -> Draw("same");
    hCalcCosThHadJet   -> Draw("same");
    lInputCosThB       -> Draw();
    fOutput            -> cd();
    cInputCosThB       -> Write();

    // create legend for blue cos-th
    TLegend* lInputCosThY = new TLegend(0.1, 0.1, 0.3, 0.35, sHeader.Data());
    lInputCosThY -> SetFillColor(0);
    lInputCosThY -> SetLineColor(0);
    lInputCosThY -> SetTextFont(42);
    lInputCosThY -> SetTextAlign(12);
    lInputCosThY -> AddEntry(hInputCosThJet, hInputCosThJet -> GetTitle(), "P");
    lInputCosThY -> AddEntry(hInputCosThHad, hInputCosThHad -> GetTitle(), "P");
    lInputCosThY -> AddEntry(hCalcCosThJetBeamY, hCalcCosThJetBeamY -> GetTitle(), "P");
    lInputCosThY -> AddEntry(hCalcCosThHadJet, hCalcCosThHadJet -> GetTitle(), "P");

    // create plot for blue cos-th
    TCanvas* cInputCosThY = new TCanvas("cInputCosThYellow", "", 750, 750);
    cInputCosThY       -> SetGrid(0, 0);
    cInputCosThY       -> cd();
    hCosThFrame        -> Draw();
    hInputCosThJet     -> Draw("same");
    hInputCosThHad     -> Draw("same");
    hCalcCosThJetBeamY -> Draw("same");
    hCalcCosThHadJet   -> Draw("same");
    lInputCosThY       -> Draw();
    fOutput            -> cd();
    cInputCosThY       -> Write();
    cInputCosThY       -> Close();
    cInputCosThY       -> Close();

    // create (z,theta) input plot
    TCanvas* cInputZTh = new TCanvas("cInputZTh", "", 1500, 1250);
    TPad*    pJet      = new TPad("pJet", "", 0.00, 0.50, 0.33, 1.00);
    TPad*    pHad      = new TPad("pHad", "", 0.33, 0.50, 0.66, 1.00);
    TPad*    pJetHad   = new TPad("pJetHad", "", 0.66, 0.50, 1.00, 1.00);
    TPad*    pBeamB    = new TPad("pBeamB", "",  0.00, 0.00, 0.33, 0.50);
    TPad*    pBeamY    = new TPad("pBeamY", "", 0.33, 0.00, 0.66, 0.50);
    cInputZTh        -> SetGrid(0, 0);
    pJet             -> SetGrid(0, 0);
    pHad             -> SetGrid(0, 0);
    pJetHad          -> SetGrid(0, 0);
    pBeamB           -> SetGrid(0, 0);
    pBeamY           -> SetGrid(0, 0);
    cInputZTh        -> cd();
    pJet             -> Draw();
    pHad             -> Draw();
    pJetHad          -> Draw();
    pBeamB           -> Draw();
    pBeamY           -> Draw();
    pJet             -> cd();
    hInputZThJet     -> Draw();
    pHad             -> cd();
    hInputZThHad     -> Draw();
    pJetHad          -> cd();
    hCalcZThHadJet   -> Draw();
    pBeamB           -> cd();
    hCalcZThJetBeamB -> Draw();
    pBeamY           -> cd();
    hCalcZThJetBeamY -> Draw();
    fOutput          -> cd();
    cInputZTh        -> Write();
    cInputZTh        -> Close();

  }  // end input theta plot making
  std::cout << "    Created theta input plots." << std::endl;

  // create (x,y) inuput plot
  {

    // create (x,y) plot
    TCanvas* cInputXY = new TCanvas("cInputXY", "", 1250, 1250);
    TPad*    pSpinB   = new TPad("pSpinB", "", 0.00, 0.50, 0.50, 1.00);
    TPad*    pSpinY   = new TPad("pSpinY", "", 0.50, 0.50, 1.00, 1.00);
    TPad*    pJet     = new TPad("pJet", "", 0.00, 0.00, 0.50, 0.50);
    TPad*    pHad     = new TPad("pHad", "", 0.50, 0.00, 1.00, 0.50);
    cInputXY      -> SetGrid(0, 0);
    pSpinB        -> SetGrid(0, 0);
    pSpinY        -> SetGrid(0, 0);
    pJet          -> SetGrid(0, 0);
    pHad          -> SetGrid(0, 0);
    cInputXY      -> cd();
    pSpinB        -> Draw();
    pSpinY        -> Draw();
    pJet          -> Draw();
    pHad          -> Draw();
    pSpinB        -> cd();
    hInputXYSpinB -> Draw();
    pSpinY        -> cd();
    hInputXYSpinY -> Draw();
    pJet          -> cd();
    hInputXYJet   -> Draw();
    pHad          -> cd();
    hInputXYHad   -> Draw();
    fOutput       -> cd();
    cInputXY      -> Write();
    cInputXY      -> Close();

  }
  std::cout << "    Created (x,y) input plot." << std::endl;

  // create phi output plots
  {

    // create legend for blue phi
    TLegend* lOutputPhiB = new TLegend(0.1, 0.1, 0.3, 0.4, sHeader.Data());
    lOutputPhiB -> SetFillColor(0);
    lOutputPhiB -> SetLineColor(0);
    lOutputPhiB -> SetTextFont(42);
    lOutputPhiB -> SetTextAlign(12);
    lOutputPhiB -> AddEntry(hPhiSpinB, hPhiSpinB -> GetTitle(), "P");
    lOutputPhiB -> AddEntry(hPhiHadB, hPhiHadB -> GetTitle(), "P");
    lOutputPhiB -> AddEntry(hPhiHad2B, hPhiHad2B -> GetTitle(), "P");
    lOutputPhiB -> AddEntry(hPhiCollB, hPhiCollB -> GetTitle(), "P");
    lOutputPhiB -> AddEntry(hPhiBoerB, hPhiBoerB -> GetTitle(), "P");

    // create plot for blue phi
    TCanvas* cOutputPhiB = new TCanvas("cOutputPhiBlue", "", 750, 750);
    cOutputPhiB -> SetGrid(0, 0);
    cOutputPhiB -> cd();
    hPhiFrame   -> Draw();
    hPhiSpinB   -> Draw("same");
    hPhiHadB    -> Draw("same");
    hPhiHad2B   -> Draw("same");
    hPhiCollB   -> Draw("same");
    hPhiBoerB   -> Draw("same");
    lOutputPhiB -> Draw();
    fOutput     -> cd();
    cOutputPhiB -> Write();
    cOutputPhiB -> Close();

    // create legend for yellow phi
    TLegend* lOutputPhiY = new TLegend(0.1, 0.1, 0.3, 0.4, sHeader.Data());
    lOutputPhiY -> SetFillColor(0);
    lOutputPhiY -> SetLineColor(0);
    lOutputPhiY -> SetTextFont(42);
    lOutputPhiY -> SetTextAlign(12);
    lOutputPhiY -> AddEntry(hPhiSpinY, hPhiSpinY -> GetTitle(), "P");
    lOutputPhiY -> AddEntry(hPhiHadY, hPhiHadY -> GetTitle(), "P");
    lOutputPhiY -> AddEntry(hPhiHad2Y, hPhiHad2Y -> GetTitle(), "P");
    lOutputPhiY -> AddEntry(hPhiCollY, hPhiCollY -> GetTitle(), "P");
    lOutputPhiY -> AddEntry(hPhiBoerY, hPhiBoerY -> GetTitle(), "P");

    // create plot for yellow phi
    TCanvas* cOutputPhiY = new TCanvas("cOutputPhiYellow", "", 750, 750);
    cOutputPhiY -> SetGrid(0, 0);
    cOutputPhiY -> cd();
    hPhiFrame   -> Draw();
    hPhiSpinY   -> Draw("same");
    hPhiHadY    -> Draw("same");
    hPhiHad2Y   -> Draw("same");
    hPhiCollY   -> Draw("same");
    hPhiBoerY   -> Draw("same");
    lOutputPhiY -> Draw();
    fOutput     -> cd();
    cOutputPhiY -> Write();
    cOutputPhiY -> Close();

  }  // end output phi plot making
  std::cout << "    Created phi output plots." << std::endl;

  // create alt phi output plots
  {

    // create legend for blue alt phi
    TLegend* lOutputAltPhiB = new TLegend(0.1, 0.1, 0.3, 0.4, sHeader.Data());
    lOutputAltPhiB -> SetFillColor(0);
    lOutputAltPhiB -> SetLineColor(0);
    lOutputAltPhiB -> SetTextFont(42);
    lOutputAltPhiB -> SetTextAlign(12);
    lOutputAltPhiB -> AddEntry(hAltPhiSpinB, hAltPhiSpinB -> GetTitle(), "P");
    lOutputAltPhiB -> AddEntry(hAltPhiHadB, hAltPhiHadB -> GetTitle(), "P");
    lOutputAltPhiB -> AddEntry(hAltPhiHad2B, hAltPhiHad2B -> GetTitle(), "P");
    lOutputAltPhiB -> AddEntry(hAltPhiCollB, hAltPhiCollB -> GetTitle(), "P");
    lOutputAltPhiB -> AddEntry(hAltPhiBoerB, hAltPhiBoerB -> GetTitle(), "P");

    // create plot for blue alt phi
    TCanvas* cOutputAltPhiB = new TCanvas("cOutputAltPhiBlue", "", 750, 750);
    cOutputAltPhiB -> SetGrid(0, 0);
    cOutputAltPhiB -> cd();
    hAltPhiFrame   -> Draw();
    hAltPhiSpinB   -> Draw("same");
    hAltPhiHadB    -> Draw("same");
    hAltPhiHad2B   -> Draw("same");
    hAltPhiCollB   -> Draw("same");
    hAltPhiBoerB   -> Draw("same");
    lOutputAltPhiB -> Draw();
    fOutput        -> cd();
    cOutputAltPhiB -> Write();
    cOutputAltPhiB -> Close();

    // create legend for yellow alt phi
    TLegend* lOutputAltPhiY = new TLegend(0.1, 0.1, 0.3, 0.4, sHeader.Data());
    lOutputAltPhiY -> SetFillColor(0);
    lOutputAltPhiY -> SetLineColor(0);
    lOutputAltPhiY -> SetTextFont(42);
    lOutputAltPhiY -> SetTextAlign(12);
    lOutputAltPhiY -> AddEntry(hAltPhiSpinY, hAltPhiSpinY -> GetTitle(), "P");
    lOutputAltPhiY -> AddEntry(hAltPhiHadY, hAltPhiHadY -> GetTitle(), "P");
    lOutputAltPhiY -> AddEntry(hAltPhiHad2Y, hAltPhiHad2Y -> GetTitle(), "P");
    lOutputAltPhiY -> AddEntry(hAltPhiCollY, hAltPhiCollY -> GetTitle(), "P");
    lOutputAltPhiY -> AddEntry(hAltPhiBoerY, hAltPhiBoerY -> GetTitle(), "P");

    // create plot for yellow alt phi
    TCanvas* cOutputAltPhiY = new TCanvas("cOutputAltPhiYellow", "", 750, 750);
    cOutputAltPhiY -> SetGrid(0, 0);
    cOutputAltPhiY -> cd();
    hAltPhiFrame   -> Draw();
    hAltPhiSpinY   -> Draw("same");
    hAltPhiHadY    -> Draw("same");
    hAltPhiHad2Y   -> Draw("same");
    hAltPhiCollY   -> Draw("same");
    hAltPhiBoerY   -> Draw("same");
    lOutputAltPhiY -> Draw();
    fOutput        -> cd();
    cOutputAltPhiY -> Write();
    cOutputAltPhiY -> Close();

  }  // end output alt phi plot making
  std::cout << "    Created alt phi output plots." << std::endl;

  // create everything plots
  {

    // create legend for blue phi
    TLegend* lEverythingPhiB = new TLegend(0.1, 0.1, 0.3, 0.6, sHeader.Data());
    lEverythingPhiB -> SetFillColor(0);
    lEverythingPhiB -> SetLineColor(0);
    lEverythingPhiB -> SetTextFont(42);
    lEverythingPhiB -> SetTextAlign(12);
    lEverythingPhiB -> AddEntry(hInputPhiSpinB, hInputPhiSpinB -> GetTitle(), "P");
    lEverythingPhiB -> AddEntry(hInputPhiJet, hInputPhiJet -> GetTitle(), "P");
    lEverythingPhiB -> AddEntry(hInputPhiHad, hInputPhiHad -> GetTitle(), "P");
    lEverythingPhiB -> AddEntry(hCalcPhiJetBeamB, hCalcPhiJetBeamB -> GetTitle(), "P");
    lEverythingPhiB -> AddEntry(hCalcPhiHadJet, hCalcPhiHadJet -> GetTitle(), "P");
    lEverythingPhiB -> AddEntry(hPhiSpinB, hPhiSpinB -> GetTitle(), "P");
    lEverythingPhiB -> AddEntry(hPhiHadB, hPhiHadB -> GetTitle(), "P");
    lEverythingPhiB -> AddEntry(hPhiHad2B, hPhiHad2B -> GetTitle(), "P");
    lEverythingPhiB -> AddEntry(hPhiCollB, hPhiCollB -> GetTitle(), "P");
    lEverythingPhiB -> AddEntry(hPhiBoerB, hPhiBoerB -> GetTitle(), "P");

    // create plot for blue phi
    TCanvas* cEverythingPhiB = new TCanvas("cEverythingPhiBlue", "", 750, 750);
    cEverythingPhiB  -> SetGrid(0, 0);
    cEverythingPhiB  -> cd();
    hPhiFrame        -> Draw();
    hInputPhiSpinB   -> Draw("same");
    hInputPhiJet     -> Draw("same");
    hInputPhiHad     -> Draw("same");
    hCalcPhiJetBeamB -> Draw("same");
    hCalcPhiHadJet   -> Draw("same");
    hPhiSpinB        -> Draw("same");
    hPhiHadB         -> Draw("same");
    hPhiHad2B        -> Draw("same");
    hPhiCollB        -> Draw("same");
    hPhiBoerB        -> Draw("same");
    lEverythingPhiB  -> Draw();
    fOutput          -> cd();
    cEverythingPhiB  -> Write();
    cEverythingPhiB  -> Close();

    // create legend for yellow phi
    TLegend* lEverythingPhiY = new TLegend(0.1, 0.1, 0.3, 0.6, sHeader.Data());
    lEverythingPhiY -> SetFillColor(0);
    lEverythingPhiY -> SetLineColor(0);
    lEverythingPhiY -> SetTextFont(42);
    lEverythingPhiY -> SetTextAlign(12);
    lEverythingPhiY -> AddEntry(hInputPhiSpinY, hInputPhiSpinY -> GetTitle(), "P");
    lEverythingPhiY -> AddEntry(hInputPhiJet, hInputPhiJet -> GetTitle(), "P");
    lEverythingPhiY -> AddEntry(hInputPhiHad, hInputPhiHad -> GetTitle(), "P");
    lEverythingPhiY -> AddEntry(hCalcPhiJetBeamY, hCalcPhiJetBeamY -> GetTitle(), "P");
    lEverythingPhiY -> AddEntry(hCalcPhiHadJet, hCalcPhiHadJet -> GetTitle(), "P");
    lEverythingPhiY -> AddEntry(hPhiSpinY, hPhiSpinY -> GetTitle(), "P");
    lEverythingPhiY -> AddEntry(hPhiHadY, hPhiHadY -> GetTitle(), "P");
    lEverythingPhiY -> AddEntry(hPhiHad2Y, hPhiHad2Y -> GetTitle(), "P");
    lEverythingPhiY -> AddEntry(hPhiCollY, hPhiCollY -> GetTitle(), "P");
    lEverythingPhiY -> AddEntry(hPhiBoerY, hPhiBoerY -> GetTitle(), "P");

    // create plot for yellow phi
    TCanvas* cEverythingPhiY = new TCanvas("cEverythingPhiYellow", "", 750, 750);
    cEverythingPhiY  -> SetGrid(0, 0);
    cEverythingPhiY  -> cd();
    hPhiFrame        -> Draw();
    hInputPhiSpinY   -> Draw("same");
    hInputPhiJet     -> Draw("same");
    hInputPhiHad     -> Draw("same");
    hCalcPhiJetBeamY -> Draw("same");
    hCalcPhiHadJet   -> Draw("same");
    hPhiSpinB        -> Draw("same");
    hPhiHadB         -> Draw("same");
    hPhiHad2B        -> Draw("same");
    hPhiCollB        -> Draw("same");
    hPhiBoerB        -> Draw("same");
    lEverythingPhiY  -> Draw();
    fOutput          -> cd();
    cEverythingPhiY  -> Write();
    cEverythingPhiY  -> Close();

  }  // end everything plot making
  std::cout << "    Created phi everything plots." << std::endl;

  // save histograms
  fOutput            -> cd();
  hInputPhiSpinB     -> Write();
  hInputPhiSpinY     -> Write();
  hInputXYSpinB      -> Write();
  hInputXYSpinY      -> Write();
  hInputPhiJet       -> Write();
  hInputThetaJet     -> Write();
  hInputCosThJet     -> Write();
  hInputZThJet       -> Write();
  hInputXYJet        -> Write();
  hInputPhiHad       -> Write();
  hInputThetaHad     -> Write();
  hInputCosThHad     -> Write();
  hInputXYHad        -> Write();
  hInputZThHad       -> Write();
  hCalcPhiJetBeamB   -> Write();
  hCalcThetaJetBeamB -> Write();
  hCalcCosThJetBeamB -> Write();
  hCalcZThJetBeamB   -> Write();
  hCalcPhiJetBeamY   -> Write();
  hCalcThetaJetBeamY -> Write();
  hCalcCosThJetBeamY -> Write();
  hCalcZThJetBeamY   -> Write();
  hCalcPhiHadJet     -> Write();
  hCalcThetaHadJet   -> Write();
  hCalcCosThHadJet   -> Write();
  hCalcZThHadJet     -> Write();
  hCalcJHBCrossMagB  -> Write(); 
  hCalcJHBDotMagB    -> Write();
  hCalcJHBCrossMagY  -> Write();
  hCalcJHBDotMagY    -> Write();
  hPhiSpinB          -> Write();
  hPhiSpinY          -> Write();
  hPhiHadB           -> Write();
  hPhiHadY           -> Write();
  hPhiHad2B          -> Write();
  hPhiHad2Y          -> Write();
  hPhiCollB          -> Write();
  hPhiCollY          -> Write();
  hPhiBoerB          -> Write();
  hPhiBoerY          -> Write();
  hAltPhiSpinB       -> Write();
  hAltPhiSpinY       -> Write();
  hAltPhiHadB        -> Write();
  hAltPhiHadY        -> Write();
  hAltPhiHad2B       -> Write();
  hAltPhiHad2Y       -> Write();
  hAltPhiCollB       -> Write();
  hAltPhiCollY       -> Write();
  hAltPhiBoerB       -> Write();
  hAltPhiBoerY       -> Write();
  std::cout << "    Saved histograms." << std::endl;

  // close output file
  fOutput -> cd();
  fOutput -> Close();

  // announce end and exit
  std::cout << "  Finished angle calculation test!\n" << std::endl;
  return;

}  // end AngleCalculationTest

// end ========================================================================
