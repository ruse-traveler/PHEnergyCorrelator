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
#include <limits>
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



// ===========================================================================
//! Namespace for acceptance-related constants
// ===========================================================================
namespace Accept {

  ///! PHENX central arm eta-acceptacne
  float EtaMin() {return -0.35;}
  float EtaMax() {return 0.35;}

  ///! PHENIX central arm west phi-acceptance
  float PhiWestMin() {return (30.0 / 180.0) * TMath::Pi();}
  float PhiWestMax() {return (120.0 / 180.0) * TMath::Pi();}

  ///! PHENIX central arm east phi-acceptance
  float PhiEastMin() {return (30.0 / 180.0) * TMath::Pi();}
  float PhiEastMax() {return (120.0 / 180.0) * TMath::Pi();}

  // -------------------------------------------------------------------------
  //! Helper method to check if (theta, phi) are in acceptance
  // -------------------------------------------------------------------------
  bool IsIn(const double theta, const double phi) {

    // calculate eta and check acceptance
    const double eta     = -1.0 * std::log(std::tan(theta / 2.0));
    const bool   isInEta = ((eta > EtaMin()) && (eta < EtaMax()));

    // wrap phi to be in (0, 2pi)
    double phiUse = phi;
    if (phi <= 0.0)           phiUse += TMath::TwoPi();
    if (phi > TMath::TwoPi()) phiUse -= TMath::TwoPi();  

    // check if phi is in west acceptance
    const bool isInWestPhi = ((phiUse > PhiWestMin()) && (phiUse < PhiWestMax()));

    // check if phi is in east acceptance
    const double phiEast     = phiUse - TMath::Pi();
    const double isInEastPhi = ((phiEast > PhiEastMin()) && (phiEast < PhiEastMax()));

    // return if in eta + phi accept
    const bool isInAccept = (isInEta && (isInWestPhi || isInEastPhi));
    return isInAccept;

  }  // end 'IsIn(double x 2)'

}  // end Const namespace



// ============================================================================
//! Test macro for anglue calculations
// ============================================================================
/*! Note: wrapMode sets what wrapping to do on the
 *  collins angle
 *    0 = wrap to [0, 2pi)
 *    1 = wrap to [0, pi)
 *  
 *  Note: doDot turns on adjusting phiHad based on
 *  the sign of the dot-product between the jet-
 *  beam and jet-hadron normals
 */
void AngleCalculationTest(
  const std::string oFile = "angleCalcTest.nIter10K_addDiFFAngles.d7m7y2025.root",
  const std::size_t nIter = 10000,
  const std::size_t wrapMode = 2,
  const bool doWrap = false,
  const bool doDot = true,
  const bool doAccept = false,
  const bool doBatch = true
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

  // initialize collins & BM histograms
  TH1D* hInputPhiSpinB        = new TH1D("hInputPhiSpinB", "#phi_{spin}^{B} input", nAngBins, xAngStart, xAngStop);
  TH1D* hInputPhiSpinY        = new TH1D("hInputPhiSpinY", "#phi_{spin}^{Y} input", nAngBins, xAngStart, xAngStop);
  TH2D* hInputXYSpinB         = new TH2D("hInputXYSpinB", "(x,y) sampled (spin B)", nXYBins, xXYStart, xXYStop, nXYBins, xXYStart, xXYStop);
  TH2D* hInputXYSpinY         = new TH2D("hInputXYSpinY", "(x,y) sampled (spin Y)", nXYBins, xXYStart, xXYStop, nXYBins, xXYStart, xXYStop);
  TH1D* hInputPhiJet          = new TH1D("hInputPhiJet", "#phi_{jet} input", nAngBins, xAngStart, xAngStop);
  TH1D* hInputThetaJet        = new TH1D("hInputThetaJet", "#theta_{jet} input", nAngBins, xAngStart, xAngStop);
  TH1D* hInputCosThJet        = new TH1D("hInputCosThJet", "cos#theta_{jet} input", nCosBins, xCosStart, xCosStop);
  TH2D* hInputXYJet           = new TH2D("hInputXYJet", "(x,y) sampled (jet)", nXYBins, xXYStart, xXYStop, nXYBins, xXYStart, xXYStop);
  TH2D* hInputZThJet          = new TH2D("hInputZThJet", "(z,#theta) sampled (jet)", nZBins, xZStart, xZStop, nAngBins, xAngStart, xAngStop);
  TH1D* hInputPhiHad          = new TH1D("hInputPhiHad", "#phi_{h} input", nAngBins, xAngStart, xAngStop);
  TH1D* hInputThetaHad        = new TH1D("hInputThetaHad", "#theta_{h} input", nAngBins, xAngStart, xAngStop);
  TH1D* hInputCosThHad        = new TH1D("hInputCosThHad", "cos#theta_{h} input", nCosBins, xCosStart, xCosStop);
  TH2D* hInputXYHad           = new TH2D("hInputXYHad", "(x,y) sampled (had)", nXYBins, xXYStart, xXYStop, nXYBins, xXYStart, xXYStop);
  TH2D* hInputZThHad          = new TH2D("hInputZThHad", "(z,#theta) sampled (had)", nZBins, xZStart, xZStop, nAngBins, xAngStart, xAngStop);
  TH1D* hCalcPhiJetBeamB      = new TH1D("hCalcPhiJetBeamB", "#phi_{jet-beam}^{B}", nAngBins, xAngStart, xAngStop);
  TH1D* hCalcThetaJetBeamB    = new TH1D("hCalcThetaJetBeamB", "#theta_{jet-beam}^{B}", nAngBins, xAngStart, xAngStop);
  TH1D* hCalcCosThJetBeamB    = new TH1D("hCalcCosThJetBeamB", "cos#theta_{jet-beam}^{B}", nCosBins, xCosStart, xCosStop);
  TH2D* hCalcZThJetBeamB      = new TH2D("hCalcZThJetBeamB", "(z,#theta) sampled (jet-blue beam)", nZBins, xZStart, xZStop, nAngBins, xAngStart, xAngStop);
  TH1D* hCalcPhiJetBeamY      = new TH1D("hCalcPhiJetBeamY", "#phi_{jet-beam}^{Y}", nAngBins, xAngStart, xAngStop);
  TH1D* hCalcThetaJetBeamY    = new TH1D("hCalcThetaJetBeamY", "#theta_{jet-beam}^{Y}", nAngBins, xAngStart, xAngStop);
  TH1D* hCalcCosThJetBeamY    = new TH1D("hCalcCosThJetBeamY", "cos#theta_{jet-beam}^{Y}", nCosBins, xCosStart, xCosStop);
  TH2D* hCalcZThJetBeamY      = new TH2D("hCalcZThJetBeamY", "(z,#theta) sampled (jet-yellow beam)", nZBins, xZStart, xZStop, nAngBins, xAngStart, xAngStop);
  TH1D* hCalcPhiHadJet        = new TH1D("hCalcPhiJetHad", "#phi_{jet-h}", nAngBins, xAngStart, xAngStop);
  TH1D* hCalcThetaHadJet      = new TH1D("hCalcThetaJetHad", "#theta_{jet-h}", nAngBins, xAngStart, xAngStop);
  TH1D* hCalcCosThHadJet      = new TH1D("hCalcCosThJetHad", "cos#theta_{jet-h}", nCosBins, xCosStart, xCosStop);
  TH2D* hCalcZThHadJet        = new TH2D("hCalcZThHadJet", "(z,#theta) sampled (jet-hadron)", nZBins, xZStart, xZStop, nAngBins, xAngStart, xAngStop);
  TH1D* hCalcJHBCrossMagB     = new TH1D("hCalcJHBCrossMagB", "|#bf{v}_{jet-beam}^{B}#times#bf{v}_{jet-had}|", nMagBins, xMagStart, xMagStop);
  TH1D* hCalcJHBDotB          = new TH1D("hCalcJHBDotB", "#bf{v}_{jet-beam}^{B} #upoint#bf{v}_{jet-had}", nMagBins, xMagStart, xMagStop);
  TH1D* hCalcJHBCrossMagY     = new TH1D("hCalcJHBCrossMagY", "|#bf{v}_{jet-beam}^{Y}#times#bf{v}_{jet-had}|", nMagBins, xMagStart, xMagStop);
  TH1D* hCalcJHBDotY          = new TH1D("hCalcJHBDotY", "#bf{v}_{jet-beam}^{Y}#upoint#bf{v}_{jet-had}", nMagBins, xMagStart, xMagStop);
  TH1D* hPhiSpinB             = new TH1D("hPhiSpinB", "#phi_{spin}^{B}", nAngBins, xAngStart, xAngStop);
  TH1D* hPhiSpinY             = new TH1D("hPhiSpinY", "#phi_{spin}^{Y}", nAngBins, xAngStart, xAngStop);
  TH1D* hPhiHadB              = new TH1D("hPhiHadB", "#phi_{h}^{B}", nAngBins, xAngStart, xAngStop);
  TH1D* hPhiHadY              = new TH1D("hPhiHadY", "#phi_{h}^{Y}", nAngBins, xAngStart, xAngStop);
  TH1D* hPhiHad2B             = new TH1D("hPhiHad2B", "2#phi_{h}^{B}", nAngBins, xAngStart, xAngStop);
  TH1D* hPhiHad2Y             = new TH1D("hPhiHad2Y", "2#phi_{h}^{Y}", nAngBins, xAngStart, xAngStop);
  TH1D* hPhiCollB             = new TH1D("hPhiCollB", "#phi_{collins}^{B}", nAngBins, xAngStart, xAngStop);
  TH1D* hPhiCollY             = new TH1D("hPhiCollY", "#phi_{collins}^{Y}", nAngBins, xAngStart, xAngStop);
  TH1D* hPhiBoerB             = new TH1D("hPhiBoerB", "#phi_{boer}^{B}", nAngBins, xAngStart, xAngStop);
  TH1D* hPhiBoerY             = new TH1D("hPhiBoerY", "#phi_{boer}^{Y}", nAngBins, xAngStart, xAngStop);
  TH1D* hAltPhiSpinB          = new TH1D("hAltPhiSpinB", "#phi_{spin}^{B}", nAngBins, xAngStart, xAngStop);
  TH1D* hAltPhiSpinY          = new TH1D("hAltPhiSpinY", "#phi_{spin}^{Y}", nAngBins, xAngStart, xAngStop);
  TH1D* hAltPhiHadB           = new TH1D("hAltPhiHadB", "#phi_{h}^{B}", nAngBins, xAngStart, xAngStop);
  TH1D* hAltPhiHadY           = new TH1D("hAltPhiHadY", "#phi_{h}^{Y}", nAngBins, xAngStart, xAngStop);
  TH1D* hAltPhiHad2B          = new TH1D("hAltPhiHad2B", "2#phi_{h}^{B}", nAngBins, xAngStart, xAngStop);
  TH1D* hAltPhiHad2Y          = new TH1D("hAltPhiHad2Y", "2#phi_{h}^{Y}", nAngBins, xAngStart, xAngStop);
  TH1D* hAltPhiCollB          = new TH1D("hAltPhiCollB", "#phi_{collins}^{B}", nAngBins, xAngStart, xAngStop);
  TH1D* hAltPhiCollY          = new TH1D("hAltPhiCollY", "#phi_{collins}^{Y}", nAngBins, xAngStart, xAngStop);
  TH1D* hAltPhiBoerB          = new TH1D("hAltPhiBoerB", "#phi_{boer}^{B}", nAngBins, xAngStart, xAngStop);
  TH1D* hAltPhiBoerY          = new TH1D("hAltPhiBoerY", "#phi_{boer}^{Y}", nAngBins, xAngStart, xAngStop);
  TH1D* hCheckPhiHadVsAltB    = new TH1D("hCheckPhiHadVsAltB", "|#phi_{h}^{B}-#acute{#phi}_{h}^{B}|", nAngBins, xAngStart, xAngStop);
  TH1D* hCheckPhiHadVsAltY    = new TH1D("hCheckPhiHadVsAltY", "|#phi_{h}^{Y}-#acute{#phi}_{h}^{Y}|", nAngBins, xAngStart, xAngStop);
  TH2D* hCheckPhiHadVsDotB    = new TH2D("hCheckPhiHadVsDotB", ";#bf{v}_{jet-beam}^{B}#upoint#bf{v}_{jet-hadron};#phi_{had}^{B}", nMagBins, xMagStart, xMagStop, nAngBins, xAngStart, xAngStop);
  TH2D* hCheckPhiHadVsDotY    = new TH2D("hCheckPhiHadVsDotY", ";#bf{v}_{jet-beam}^{Y}#upoint#bf{v}_{jet-hadron};#phi_{had}^{Y}", nMagBins, xMagStart, xMagStop, nAngBins, xAngStart, xAngStop);
  TH2D* hCheckAltPhiHadVsDotB = new TH2D("hCheckAltPhiHadVsDotB", ";#bf{v}_{jet-beam}^{B}#upoint#bf{v}_{jet-hadron};#phi_{had}^{B}", nMagBins, xMagStart, xMagStop, nAngBins, xAngStart, xAngStop);
  TH2D* hCheckAltPhiHadVsDotY = new TH2D("hCheckAltPhiHadVsDotY", ";#bf{v}_{jet-beam}^{Y}#upoint#bf{v}_{jet-hadron};#phi_{had}^{Y}", nMagBins, xMagStart, xMagStop, nAngBins, xAngStart, xAngStop);

  // initialize DiFF hists
  TH1D* hInputPhiSec   = new TH1D("hInputPhiSec", "#phi_{sh} input", nAngBins, xAngStart, xAngStop);
  TH1D* hInputThetaSec = new TH1D("hInputThetaSec", "#theta_{sh} input", nAngBins, xAngStart, xAngStop);
  TH1D* hInputCosThSec = new TH1D("hInputCosThSec", "cos#theta_{sh} input", nCosBins, xCosStart, xCosStop);
  TH2D* hInputXYSec    = new TH2D("hInputXYSec", "(x,y) sampled (2nd had)", nXYBins, xXYStart, xXYStop, nXYBins, xXYStart, xXYStop);
  TH2D* hInputZThSec   = new TH2D("hInputZThSec", "(z,#theta) sampled (2nd had)", nZBins, xZStart, xZStop, nAngBins, xAngStart, xAngStop);
  TH1D* hThSpinBeamB   = new TH1D("hThSpinBeamB", "#theta(spin, beam)_{B}", nAngBins, xAngStart, xAngStop);
  TH1D* hThSpinPCB     = new TH1D("hThSpinPCB", "#theta(spin, PC)_{B}", nAngBins, xAngStart, xAngStop);
  TH1D* hThSpinBeamY   = new TH1D("hThSpinBeamY", "#theta(spin, beam)_{Y}", nAngBins, xAngStart, xAngStop);
  TH1D* hThSpinPCY     = new TH1D("hThSpinPCY", "#theta(spin, PC)_{Y}", nAngBins, xAngStart, xAngStop);
  TH1D* hThPCRCY       = new TH1D("hThPCRCY", "#theta(PC, RC)_{Y}", nAngBins, xAngStart, xAngStop);
  TH1D* hThBeamRCY     = new TH1D("hThBeamRCY", "#theta(beam, RC)_{Y}", nAngBins, xAngStart, xAngStop);
  TH1D* hThetaSB       = new TH1D("hThetaSB", "#theta_{SB}", nAngBins, xAngStart, xAngStop);
  TH1D* hThetaSY       = new TH1D("hThetaSY", "#theta_{SY}", nAngBins, xAngStart, xAngStop);
  TH1D* hThetaRC       = new TH1D("hThetaRC", "#theta_{RC}", nAngBins, xAngStart, xAngStop);
  TH1D* hThetaSBRC     = new TH1D("hThetaSBRC", "#theta_{SB} - #theta_{RC}", nAngBins, xAngStart, xAngStop);
  TH1D* hThetaSYRC     = new TH1D("hThetaSYRC", "#theta_{SY} - #theta_{RC}", nAngBins, xAngStart, xAngStop);
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

    // get random u for jet/hadron vectors, translate to r
    const double uRandJet = rando -> Uniform(0.0, 1.0);
    const double uRandHad = rando -> Uniform(0.0, 1.0);
    const double uRandSec = rando -> Uniform(0.0, 1.0);
    const double rRandJet = pow(uRandJet, 1.0/3.0);
    const double rRandHad = pow(uRandHad, 1.0/3.0);
    const double rRandSec = pow(uRandSec, 1.0/3.0);

    // get jet, hadron vectors and constrain to
    // PHENIX acceptance if need be
    bool   goodJetVals = false;
    bool   goodHadVals = false;
    bool   goodSecVals = false;
    double phiUseJet   = 0.0;
    double phiUseHad   = 0.0;
    double phiUseSec   = 0.0;
    double thetaUseJet = 0.0;
    double thetaUseHad = 0.0;
    double thetaUseSec = 0.0;
    do {

      // get random phi = (0, 2pi), cos-theta = (-1, 1)
      const double phiRandJet = rando -> Uniform(0.0, TMath::TwoPi());
      const double phiRandHad = rando -> Uniform(0.0, TMath::TwoPi());
      const double phiRandSec = rando -> Uniform(0.0, TMath::TwoPi());
      const double cosRandJet = rando -> Uniform(-1.0, 1.0);
      const double cosRandHad = rando -> Uniform(-1.0, 1.0);
      const double cosRandSec = rando -> Uniform(-1.0, 1.0);

      // translate cos-theta
      const double thetaRandJet = acos(cosRandJet);
      const double thetaRandHad = acos(cosRandHad);
      const double thetaRandSec = acos(cosRandSec);

      // now check if in acceptance if needed
      goodJetVals = doAccept ? Accept::IsIn(thetaRandJet, phiRandJet) : true;
      goodHadVals = doAccept ? Accept::IsIn(thetaRandHad, phiRandHad) : true;
      goodSecVals = doAccept ? Accept::IsIn(thetaRandSec, phiRandSec) : true;

      // and set use values if good
      if (goodJetVals && goodHadVals && goodSecVals) {
        phiUseJet   = phiRandJet;
        phiUseHad   = phiRandHad;
        phiUseSec   = phiRandSec;
        thetaUseJet = thetaRandJet;
        thetaUseHad = thetaRandHad;
        thetaUseSec = thetaRandSec;
      }

    } while (!goodJetVals || !goodHadVals);

    // now translate into (x, y, z) values (assuming
    // r = 1.0)
    const double xUseJet = rRandJet * sin(thetaUseJet) * cos(phiUseJet);
    const double yUseJet = rRandJet * sin(thetaUseJet) * sin(phiUseJet);
    const double zUseJet = rRandJet * cos(thetaUseJet);
    const double xUseHad = rRandHad * sin(thetaUseHad) * cos(phiUseHad);
    const double yUseHad = rRandHad * sin(thetaUseHad) * sin(phiUseHad);
    const double zUseHad = rRandHad * cos(thetaUseHad);
    const double xUseSec = rRandSec * sin(thetaUseSec) * cos(phiUseSec);
    const double yUseSec = rRandSec * sin(thetaUseSec) * sin(phiUseSec);
    const double zUseSec = rRandSec * cos(thetaUseSec);

    // get random jet/hadron directions, fill input histograms
    TVector3 vecJet3(xUseJet, yUseJet, zUseJet);
    TVector3 vecHad3(xUseHad, yUseHad, zUseHad);
    TVector3 vecSec3(xUseSec, yUseSec, zUseSec);
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
    hInputPhiSec   -> Fill( vecSec3.Phi() );
    hInputThetaSec -> Fill( vecSec3.Theta() );
    hInputCosThSec -> Fill( cos(vecSec3.Theta()) );
    hInputXYSec    -> Fill( vecSec3.X(), vecSec3.Y() );
    hInputZThSec   -> Fill( vecSec3.Z(), vecSec3.Theta() );

    // now normalize spin/jet/hadron vectors
    TVector3 unitSpinB3 = vecSpinB3.Unit();
    TVector3 unitSpinY3 = vecSpinY3.Unit();
    TVector3 unitJet3   = vecJet3.Unit();
    TVector3 unitHad3   = vecHad3.Unit();
    TVector3 unitSec3   = vecSec3.Unit();

    // collins & boer-mulders angle calculations --------------------------

    // (0) get vectors normal to the spin-beam & jet-beam plane
    std::pair<TVector3, TVector3> normSpinBeam3 = std::make_pair(
      ( vecBeamB3.Cross(unitSpinB3) ).Unit(),
      ( vecBeamY3.Cross(unitSpinY3) ).Unit()
    );
    std::pair<TVector3, TVector3> normJetBeam3 = std::make_pair(
      ( vecBeamB3.Cross(unitJet3) ).Unit(),
      ( vecBeamY3.Cross(unitJet3) ).Unit()
    );

    // fill intermediate histograms
    hCalcPhiJetBeamB   -> Fill( normJetBeam3.first.Phi() );
    hCalcThetaJetBeamB -> Fill( normJetBeam3.first.Theta() );
    hCalcCosThJetBeamB -> Fill( cos(normJetBeam3.first.Theta()) );
    hCalcZThJetBeamB   -> Fill( normJetBeam3.first.Z(), normJetBeam3.first.Theta() );
    hCalcPhiJetBeamY   -> Fill( normJetBeam3.second.Phi() );
    hCalcThetaJetBeamY -> Fill( normJetBeam3.second.Theta() );
    hCalcCosThJetBeamY -> Fill( cos(normJetBeam3.second.Theta()) );
    hCalcZThJetBeamY   -> Fill( normJetBeam3.second.Z(), normJetBeam3.second.Theta() );

    // ------------------------------------------------------------------------

    // (1) cross spin into jet-beam plane
    std::pair<TVector3, TVector3> crossJetBeamSpin3 = std::make_pair(
      normJetBeam3.first.Cross( normSpinBeam3.first ),
      normJetBeam3.second.Cross( normSpinBeam3.second )
    );

    // check dot product of spin into jet-beam
    const double spinJetBeamDotB = normJetBeam3.first.Dot( unitSpinB3 );
    const double spinJetBeamDotY = normJetBeam3.second.Dot( unitSpinY3 );

    // (2) get phiSpin: angles between the jet-beam plane and spin
    double phiSpinBlue = atan2( crossJetBeamSpin3.first.Mag(), normJetBeam3.first.Dot(normSpinBeam3.first) );
    double phiSpinYell = atan2( crossJetBeamSpin3.second.Mag(), normJetBeam3.second.Dot(normSpinBeam3.second) );

    // if using the dot product, take the outer angle if the
    // spin-jet-beam dot product is negative
    if (doDot) {
      if (spinJetBeamDotB < 0.0) phiSpinBlue = TMath::TwoPi() - phiSpinBlue;
      if (spinJetBeamDotY < 0.0) phiSpinYell = TMath::TwoPi() - phiSpinYell;
    }

    // if doing wrapping, constrain to range [0, 2pi)
    if (doWrap) {
      if (phiSpinBlue < 0)               phiSpinBlue += TMath::TwoPi();
      if (phiSpinBlue >= TMath::TwoPi()) phiSpinBlue -= TMath::TwoPi();
      if (phiSpinYell < 0)               phiSpinYell += TMath::TwoPi();
      if (phiSpinYell >= TMath::TwoPi()) phiSpinYell -= TMath::TwoPi();
    }

    // fill spin histograms
    hPhiSpinB -> Fill(phiSpinBlue);
    hPhiSpinY -> Fill(phiSpinYell);

    // ------------------------------------------------------------------------

    // (3) get vector normal to hadron average-jet plane,
    TVector3 normHadJet3 = ( unitJet3.Cross(unitHad3) ).Unit();

    // fill intermediate histograms
    hCalcPhiHadJet   -> Fill( normHadJet3.Phi() );
    hCalcThetaHadJet -> Fill( normHadJet3.Theta() );
    hCalcCosThHadJet -> Fill( cos(normHadJet3.Theta()) );
    hCalcZThHadJet   -> Fill( normHadJet3.Z(), normHadJet3.Theta() );

    // (4) cross jet-hadron normal into jet-beam normal
    std::pair<TVector3, TVector3> crossJetBeamHadron3 = std::make_pair(
      normJetBeam3.first.Cross(normHadJet3),
      normJetBeam3.second.Cross(normHadJet3)
    );

    // check dot product of hadron into jet-beam
    const double hadJetBeamDotB = normJetBeam3.first.Dot( unitHad3 );
    const double hadJetBeamDotY = normJetBeam3.second.Dot( unitHad3 );

    // fill intermediate histograms
    hCalcJHBCrossMagB -> Fill( crossJetBeamHadron3.first.Mag() ); 
    hCalcJHBDotB      -> Fill( normJetBeam3.first.Dot(normHadJet3) );
    hCalcJHBCrossMagY -> Fill( crossJetBeamHadron3.second.Mag() );
    hCalcJHBDotY      -> Fill( normJetBeam3.second.Dot(normHadJet3) );

    // (5) get phiHadron: angle between the jet-beam plane and
    //     jet-hadron plane
    double phiHadBlue = atan2( crossJetBeamHadron3.first.Mag(), normJetBeam3.first.Dot(normHadJet3) );
    double phiHadYell = atan2( crossJetBeamHadron3.second.Mag(), normJetBeam3.second.Dot(normHadJet3) );

    // if using the dot product, take the outer angle if the dot
    // product is negative
    if (doDot) {
      if (hadJetBeamDotB < 0.0) phiHadBlue = TMath::TwoPi() - phiHadBlue;
      if (hadJetBeamDotY < 0.0) phiHadYell = TMath::TwoPi() - phiHadYell;
    }

    // if doing wrapping, constrain to range [0,2pi)
    if (doWrap) {
      if (phiHadBlue < 0)               phiHadBlue += TMath::TwoPi();
      if (phiHadBlue >= TMath::TwoPi()) phiHadBlue -= TMath::TwoPi();
      if (phiHadYell < 0)               phiHadYell += TMath::TwoPi();
      if (phiHadYell >= TMath::TwoPi()) phiHadYell -= TMath::TwoPi();
    }

    // fill histograms
    hPhiHadB           -> Fill(phiHadBlue);
    hPhiHadY           -> Fill(phiHadYell);
    hCheckPhiHadVsDotB -> Fill(hadJetBeamDotB, phiHadBlue);
    hCheckPhiHadVsDotY -> Fill(hadJetBeamDotY, phiHadYell);

    // (6) double phiHadron for boer-mulders,
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

    // ------------------------------------------------------------------------

    // (7) now calculate phiColl: phiSpin - phiHadron,
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

    // (8) now calculate phiBoer: phiSpin - (2 * phiHadron),
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

    // (2) starting at get phiSpin: angles between the jet-beam plane and
    //     spin-beam plane
    //   - constrain to [0,2pi)
    double phiSpinBlueAlt = acos( normJetBeam3.first.Dot(normSpinBeam3.first) / (normJetBeam3.first.Mag() * normSpinBeam3.first.Mag()) );
    double phiSpinYellAlt = acos( normJetBeam3.second.Dot(normSpinBeam3.second) / (normJetBeam3.second.Mag() * normSpinBeam3.second.Mag()) );

    // if using the dot product, take the outer angle if the
    // spin-jet-beam dot product is negative
    if (doDot) {
      if (spinJetBeamDotB < 0.0) phiSpinBlueAlt = TMath::TwoPi() - phiSpinBlueAlt;
      if (spinJetBeamDotY < 0.0) phiSpinYellAlt = TMath::TwoPi() - phiSpinYellAlt;
    }

    // if doing wrapping, constrain to range [0, 2pi)
    if (doWrap) {
      if (phiSpinBlueAlt < 0)               phiSpinBlueAlt += TMath::TwoPi();
      if (phiSpinBlueAlt >= TMath::TwoPi()) phiSpinBlueAlt -= TMath::TwoPi();
      if (phiSpinYellAlt < 0)               phiSpinYellAlt += TMath::TwoPi();
      if (phiSpinYellAlt >= TMath::TwoPi()) phiSpinYellAlt -= TMath::TwoPi();
    }

    // fill spin histograms
    hAltPhiSpinB -> Fill(phiSpinBlueAlt);
    hAltPhiSpinY -> Fill(phiSpinYellAlt);

    // ------------------------------------------------------------------------

    // (4) get dot products and magnitudes
    double jhAltMag   = normHadJet3.Mag();
    double jbAltMagB  = normJetBeam3.first.Mag();
    double jbAltMagY  = normJetBeam3.second.Mag();
    double jhbAltDotB = normJetBeam3.first.Dot(normHadJet3);
    double jhbAltDotY = normJetBeam3.second.Dot(normHadJet3);


    // (5) now jump to phiHadron: angle between the jet-beam plane and the
    //     jet-hadron plane
    double phiHadBlueAlt = acos( normJetBeam3.first.Dot(normHadJet3) / (normJetBeam3.first.Mag() * normHadJet3.Mag()) );
    double phiHadYellAlt = acos( normJetBeam3.second.Dot(normHadJet3) / (normJetBeam3.second.Mag() * normHadJet3.Mag()) );

    // if using the dot product, take the outer angle if the dot
    // product is negative
    if (doDot) {
      if (hadJetBeamDotB < 0) phiHadBlueAlt = TMath::TwoPi() - phiHadBlueAlt;
      if (hadJetBeamDotY < 0) phiHadYellAlt = TMath::TwoPi() - phiHadYellAlt;
    }

    // if doing wrapping, constrain to range [0,2pi)
    if (doWrap) {
      if (phiHadBlueAlt < 0)               phiHadBlueAlt += TMath::TwoPi();
      if (phiHadBlueAlt >= TMath::TwoPi()) phiHadBlueAlt -= TMath::TwoPi();
      if (phiHadYellAlt < 0)               phiHadYellAlt += TMath::TwoPi();
      if (phiHadYellAlt >= TMath::TwoPi()) phiHadYellAlt -= TMath::TwoPi();
    }

    // cross-check: take difference between normal and alternate
    // phi hadron calculations
    double diffPhiHadBlue = abs(phiHadBlue - phiHadBlueAlt);
    double diffPhiHadYell = abs(phiHadYell - phiHadYellAlt);
    if (diffPhiHadBlue < std::numeric_limits<float>::epsilon()) diffPhiHadBlue = 0.0;
    if (diffPhiHadYell < std::numeric_limits<float>::epsilon()) diffPhiHadYell = 0.0;

    // fill intermediate histograms
    hCheckPhiHadVsAltB -> Fill(diffPhiHadBlue); 
    hCheckPhiHadVsAltY -> Fill(diffPhiHadYell);

    // fill histograms
    hAltPhiHadB           -> Fill(phiHadBlueAlt);
    hAltPhiHadY           -> Fill(phiHadYellAlt);
    hCheckAltPhiHadVsDotB -> Fill(hadJetBeamDotB, phiHadBlueAlt);
    hCheckAltPhiHadVsDotY -> Fill(hadJetBeamDotY, phiHadYellAlt);

    // (6) double phiHadron for boer-mulders,
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

    // ------------------------------------------------------------------------

    // (7) now calculate phiColl: phiSpin - phiHadron,
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

    // DiFF angle calculation -------------------------------------------------

    // (0) calculate PC & RC vectors 
    //     - note difference in notation
    //       * PC --> vecHadPC3
    //       * RC --> vecHadRC3
    //       * PB --> vecBeamB3
    //       * PA --> vecBeamY3
    //       * SB --> vecSpinB3
    //       * SA --> vecSpinY3
    TVector3 vecHadPC3  = vecHad3 + vecSec3;
    TVector3 vecHadRC3  = 0.5 * (vecHad3 - vecSec3);
    TVector3 unitHadPC3 = vecHadPC3.Unit();

    // get unit beam vectors (just to be safe)
    TVector3 unitBeamB3 = vecBeamB3.Unit();
    TVector3 unitBeamY3 = vecBeamY3.Unit();

    // (1) calculate blue spin-beam & spin-PC angles
    //     - cThetaSB --> thetaSpinBeamB
    //     - sThetaSB --> thetaSpinPCB
    double thetaSpinBeamB = (
      unitBeamB3.Cross(vecHadPC3) * (
        1.0 / (
          unitBeamB3.Cross(vecHadPC3).Mag()
        )
      )
    ).Dot(
      unitBeamB3.Cross(vecSpinB3) * (
        1.0 / (
          unitBeamB3.Cross(vecSpinB3).Mag()
        )
      )
    );
    double thetaSpinPCB = (
      vecHadPC3.Cross(vecSpinB3)
    ).Dot(unitBeamB3) * (
      1.0 / (
        (unitBeamB3.Cross(vecHadPC3).Mag()) *
        (unitBeamB3.Cross(vecSpinB3).Mag())
      )
    );

    // and calculate yellow spin-beam & spin-PC angles
    // - cThetaSA --> thetaSpinBeamY
    // - sThetaSA --> thetaSpinPCY
    double thetaSpinBeamY = (
      unitBeamY3.Cross(vecHadPC3) * (
        1.0 / (
          unitBeamY3.Cross(vecHadPC3).Mag()
        )
      )
    ).Dot(
      unitBeamY3.Cross(vecSpinY3) * (
        1.0 / (
          unitBeamY3.Cross(vecSpinY3).Mag()
        )
      )
    );
    double thetaSpinPCY = (
      vecHadPC3.Cross(vecSpinY3)
    ).Dot(unitBeamY3) * (
      1.0 / (
        (unitBeamY3.Cross(vecHadPC3).Mag()) * 
        (unitBeamY3.Cross(vecSpinY3).Mag())
      )
    );

    // (2) calculate RC angles
    //     - cThetaRC --> thetaPCRCY
    //     - sThetaRC --> thetaBeamRCY
    double thetaPCRCY = (
      unitHadPC3.Cross(vecBeamY3) * (
        1.0 / (
          unitHadPC3.Cross(vecBeamY3).Mag()
        )
      )
    ).Dot(
      unitHadPC3.Cross(vecHadRC3) * (
        1.0 / (
          unitHadPC3.Cross(vecHadRC3).Mag()
        )
      )
    );
    double thetaBeamRCY = (
      vecBeamY3.Cross(vecHadRC3)
    ).Dot(unitHadPC3) * (
      1.0 / (
        (unitHadPC3.Cross(vecBeamY3).Mag()) *
        (unitHadPC3.Cross(vecHadRC3).Mag())
      )
    );

    // fill histograms
    hThSpinBeamB -> Fill( thetaSpinBeamB );
    hThSpinPCB   -> Fill( thetaSpinPCB );
    hThSpinBeamY -> Fill( thetaSpinBeamY );
    hThSpinPCY   -> Fill( thetaSpinPCY );
    hThPCRCY     -> Fill( thetaPCRCY );
    hThBeamRCY   -> Fill( thetaBeamRCY );

    // ------------------------------------------------------------------------

    // (3) make sure angles span full [-pi, pi]
    //     - cThetaSB --> thetaSpinBeamB
    //     - sThetaSB --> thetaSpinPCB
    //     - cThetaSA --> thetaSpinBeamY
    //     - sThetaSA --> thetaSpinPCY
    //     - cThetaRC --> thetaPCRCY
    //     - sThetaRC --> thetaBeamRCY
    //     - ThetaSB  --> thetaSB
    //     - ThetaSA  --> thetaSY
    //     - ThetaRC  --> thetaRC
    double thetaSB = (thetaSpinPCB > 0.0) ? acos(thetaSpinBeamB) : -acos(thetaSpinBeamB);
    double thetaSY = (thetaSpinPCY > 0.0) ? acos(thetaSpinBeamY) : -acos(thetaSpinBeamY);
    double thetaRC = (thetaBeamRCY > 0.0) ? acos(thetaPCRCY) : -cos(thetaPCRCY);

    // (4) wrap angles
    if (doWrap) {
      if (wrapMode == 1) {
        if (thetaSB < 0)            thetaSB += TMath::Pi();
        if (thetaSB >= TMath::Pi()) thetaSB -= TMath::Pi();
        if (thetaSY < 0)            thetaSY += TMath::Pi();
        if (thetaSY >= TMath::Pi()) thetaSY -= TMath::Pi();
        if (thetaRC < 0)            thetaRC += TMath::Pi();
        if (thetaRC >= TMath::Pi()) thetaRC -= TMath::Pi();
      } else {
        if (thetaSB < 0)               thetaSB += TMath::TwoPi();
        if (thetaSB >= TMath::TwoPi()) thetaSB -= TMath::TwoPi();
        if (thetaSY < 0)               thetaSY += TMath::TwoPi();
        if (thetaSY >= TMath::TwoPi()) thetaSY -= TMath::TwoPi();
        if (thetaRC < 0)               thetaRC += TMath::TwoPi();
        if (thetaRC >= TMath::TwoPi()) thetaRC -= TMath::TwoPi();
      }
    }

    // fill histograms
    hThetaSB -> Fill( thetaSB );
    hThetaSY -> Fill( thetaSY );
    hThetaRC -> Fill( thetaRC );

    // ------------------------------------------------------------------------

    // (5) take difference
    double thetaSBRC = thetaSB - thetaRC;
    double thetaSYRC = thetaSY - thetaRC;

    // (6) and wrap differences
    if (doWrap) {
      if (wrapMode == 1) {
        if (thetaSBRC < 0)            thetaSBRC += TMath::Pi();
        if (thetaSBRC >= TMath::Pi()) thetaSBRC -= TMath::Pi();
        if (thetaSYRC < 0)            thetaSBRC += TMath::Pi();
        if (thetaSYRC >= TMath::Pi()) thetaSBRC -= TMath::Pi();
      } else {
        if (thetaSBRC < 0)               thetaSBRC += TMath::TwoPi();
        if (thetaSBRC >= TMath::TwoPi()) thetaSBRC -= TMath::TwoPi();
        if (thetaSYRC < 0)               thetaSBRC += TMath::TwoPi();
        if (thetaSYRC >= TMath::TwoPi()) thetaSBRC -= TMath::TwoPi();
      }
    }

    // fill histograms
    hThetaSBRC -> Fill( thetaSBRC );
    hThetaSYRC -> Fill( thetaSYRC );

  }  // end iter loop
  std::cout << "    MC loop finished!" << std::endl;

  // --------------------------------------------------------------------------

  // normalize collins & BM histograms
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

  // normalize DiFF histograms
  hInputPhiSec   -> Scale( 1. / (double) nIter);
  hInputThetaSec -> Scale( 1. / (double) nIter);
  hInputCosThSec -> Scale( 1. / (double) nIter);
  hThSpinBeamB   -> Scale( 1. / (double) nIter);
  hThSpinPCB     -> Scale( 1. / (double) nIter);
  hThSpinBeamY   -> Scale( 1. / (double) nIter);
  hThSpinPCY     -> Scale( 1. / (double) nIter);
  hThPCRCY       -> Scale( 1. / (double) nIter);
  hThBeamRCY     -> Scale( 1. / (double) nIter);
  hThetaSB       -> Scale( 1. / (double) nIter);
  hThetaSY       -> Scale( 1. / (double) nIter);
  hThetaRC       -> Scale( 1. / (double) nIter);
  hThetaSBRC     -> Scale( 1. / (double) nIter);
  hThetaSYRC     -> Scale( 1. / (double) nIter);
  std::cout << "    Normalized histograms." << std::endl;

  // --------------------------------------------------------------------------

  // create frame histograms
  TH1D* hPhiFrame    = (TH1D*) hInputPhiSpinB -> Clone();
  TH1D* hAltPhiFrame = (TH1D*) hInputPhiSpinB -> Clone();
  TH1D* hThetaFrame  = (TH1D*) hInputThetaJet -> Clone();
  TH1D* hCosThFrame  = (TH1D*) hInputCosThJet -> Clone();
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

  // set collins & BM styles
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

  // set DiFF styles
  hInputPhiSec   -> SetLineColor(col[3]);
  hInputPhiSec   -> SetMarkerColor(col[3]);
  hInputPhiSec   -> SetMarkerStyle(mar[3]);
  hInputThetaSec -> SetLineColor(col[3]);
  hInputThetaSec -> SetMarkerColor(col[3]);
  hInputThetaSec -> SetMarkerStyle(mar[3]);
  hInputCosThSec -> SetLineColor(col[3]);
  hInputCosThSec -> SetMarkerColor(col[3]);
  hInputCosThSec -> SetMarkerStyle(mar[3]);
  hThSpinBeamB   -> SetLineColor(col[0]);
  hThSpinBeamB   -> SetMarkerColor(col[0]);
  hThSpinBeamB   -> SetMarkerStyle(mar[0]);
  hThSpinPCB     -> SetLineColor(col[1]);
  hThSpinPCB     -> SetMarkerColor(col[1]);
  hThSpinPCB     -> SetMarkerStyle(mar[1]);
  hThSpinBeamY   -> SetLineColor(col[0]);
  hThSpinBeamY   -> SetMarkerColor(col[0]);
  hThSpinBeamY   -> SetMarkerStyle(mar[0]);
  hThSpinPCY     -> SetLineColor(col[1]);
  hThSpinPCY     -> SetMarkerColor(col[1]);
  hThSpinPCY     -> SetMarkerStyle(mar[1]);
  hThPCRCY       -> SetLineColor(col[2]);
  hThPCRCY       -> SetMarkerColor(col[2]);
  hThPCRCY       -> SetMarkerStyle(mar[2]);
  hThBeamRCY     -> SetLineColor(col[3]);
  hThBeamRCY     -> SetMarkerColor(col[3]);
  hThBeamRCY     -> SetMarkerStyle(mar[3]);
  hThetaSB       -> SetLineColor(col[4]);
  hThetaSB       -> SetMarkerColor(col[4]);
  hThetaSB       -> SetMarkerStyle(mar[4]);
  hThetaSY       -> SetLineColor(col[5]);
  hThetaSY       -> SetMarkerColor(col[5]);
  hThetaSY       -> SetMarkerStyle(mar[5]);
  hThetaRC       -> SetLineColor(col[6]);
  hThetaRC       -> SetMarkerColor(col[6]);
  hThetaRC       -> SetMarkerStyle(mar[6]);
  hThetaSBRC     -> SetLineColor(col[7]);
  hThetaSBRC     -> SetMarkerColor(col[7]);
  hThetaSBRC     -> SetMarkerStyle(mar[7]);
  hThetaSYRC     -> SetLineColor(col[8]);
  hThetaSYRC     -> SetMarkerColor(col[8]);
  hThetaSYRC     -> SetMarkerStyle(mar[8]);
  std::cout << "    Set styles." << std::endl;

  // --------------------------------------------------------------------------

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

  // --------------------------------------------------------------------------

  // create DiFF input plots
  {

    // create legend for blue DiFF phi input
    TLegend* lDiffPhiInputB = new TLegend(0.1, 0.1, 0.3, 0.4, sHeader.Data());
    lDiffPhiInputB -> SetFillColor(0);
    lDiffPhiInputB -> SetLineColor(0);
    lDiffPhiInputB -> SetTextFont(42);
    lDiffPhiInputB -> SetTextAlign(12);
    lDiffPhiInputB -> AddEntry(hInputPhiSpinB, hInputPhiSpinB -> GetTitle(), "P");
    lDiffPhiInputB -> AddEntry(hInputPhiHad, hInputPhiHad -> GetTitle(), "P");
    lDiffPhiInputB -> AddEntry(hInputPhiSec, hInputPhiSec-> GetTitle(), "P");

    // create plot for blue DiFF phi input
    TCanvas* cDiffPhiInputB = new TCanvas("cDiffPhiInputB", "", 750, 750);
    cDiffPhiInputB -> SetGrid(0, 0);
    cDiffPhiInputB -> cd();
    hPhiFrame      -> Draw();
    hInputPhiSpinB -> Draw("same");
    hInputPhiHad   -> Draw("same");
    hInputPhiSec   -> Draw("same");
    lDiffPhiInputB -> Draw();
    fOutput        -> cd();
    cDiffPhiInputB -> Write();
    cDiffPhiInputB -> Close();

    // create legend for yellow DiFF phi input
    TLegend* lDiffPhiInputY = new TLegend(0.1, 0.1, 0.3, 0.4, sHeader.Data());
    lDiffPhiInputY -> SetFillColor(0);
    lDiffPhiInputY -> SetLineColor(0);
    lDiffPhiInputY -> SetTextFont(42);
    lDiffPhiInputY -> SetTextAlign(12);
    lDiffPhiInputY -> AddEntry(hInputPhiSpinY, hInputPhiSpinY -> GetTitle(), "P");
    lDiffPhiInputY -> AddEntry(hInputPhiHad, hInputPhiHad -> GetTitle(), "P");
    lDiffPhiInputY -> AddEntry(hInputPhiSec, hInputPhiSec-> GetTitle(), "P");

    // create plot for yellow DiFF phi input
    TCanvas* cDiffPhiInputY = new TCanvas("cDiffPhiInputY", "", 750, 750);
    cDiffPhiInputY -> SetGrid(0, 0);
    cDiffPhiInputY -> cd();
    hPhiFrame      -> Draw();
    hInputPhiSpinY -> Draw("same");
    hInputPhiHad   -> Draw("same");
    hInputPhiSec   -> Draw("same");
    lDiffPhiInputY -> Draw();
    fOutput        -> cd();
    cDiffPhiInputY -> Write();
    cDiffPhiInputY -> Close();

  }  // end DiFF phi input plot making
  std::cout << "    Created DiFF phi input plots." << std::endl;

  // create (x,y) inuput plot
  {

    // create (x,y) plot
    TCanvas* cDiffInputXY = new TCanvas("cDiffInputXY", "", 1250, 1250);
    TPad*    pDiffSpinB   = new TPad("pDiffSpinB", "", 0.00, 0.50, 0.50, 1.00);
    TPad*    pDiffSpinY   = new TPad("pDiffSpinY", "", 0.50, 0.50, 1.00, 1.00);
    TPad*    pDiffHad     = new TPad("pDiffHad", "", 0.00, 0.00, 0.50, 0.50);
    TPad*    pDiffSec     = new TPad("pDiffSec", "", 0.50, 0.00, 1.00, 0.50);
    cDiffInputXY      -> SetGrid(0, 0);
    pDiffSpinB        -> SetGrid(0, 0);
    pDiffSpinY        -> SetGrid(0, 0);
    pDiffHad          -> SetGrid(0, 0);
    pDiffSec          -> SetGrid(0, 0);
    cDiffInputXY      -> cd();
    pDiffSpinB        -> Draw();
    pDiffSpinY        -> Draw();
    pDiffHad          -> Draw();
    pDiffSec          -> Draw();
    pDiffSpinB        -> cd();
    hInputXYSpinB     -> Draw();
    pDiffSpinY        -> cd();
    hInputXYSpinY     -> Draw();
    pDiffHad          -> cd();
    hInputXYHad       -> Draw();
    pDiffSec          -> cd();
    hInputXYSec       -> Draw();
    fOutput           -> cd();
    cDiffInputXY      -> Write();
    cDiffInputXY      -> Close();

  }  // end making DiFF (x,y) input plots
  std::cout << "    Created DiFF (x,y) input plot." << std::endl;

  // create Diff input theta plots
  {

    TCanvas* cDiffThetaInput = new TCanvas("cDiffThetaInput", "", 750, 750);
    TPad*    pDiffTheta      = new TPad("pDiffTheta", "", 0.00, 0.50, 0.50, 1.00);
    TPad*    pDiffCosTh      = new TPad("pDiffCosTh", "", 0.50, 0.50, 1.00, 1.00);
    TPad*    pDiffZTh        = new TPad("pDiffZTh",   "", 0.00, 0.00, 0.50, 0.50);
    cDiffThetaInput -> SetGrid(0, 0);
    pDiffTheta      -> SetGrid(0, 0);
    pDiffCosTh      -> SetGrid(0, 0);
    pDiffZTh        -> SetGrid(0, 0);
    cDiffThetaInput -> cd();
    pDiffTheta      -> Draw();
    pDiffCosTh      -> Draw();
    pDiffZTh        -> Draw();
    pDiffTheta      -> cd();
    hInputThetaSec  -> Draw();
    pDiffCosTh      -> cd();
    hInputCosThSec  -> Draw();
    pDiffZTh        -> cd();
    hInputZThSec    -> Draw();
    fOutput         -> cd();
    cDiffThetaInput -> Write();
    cDiffThetaInput -> Close();

  }  // end making DiFF input theta plots
  std::cout << "    Create DiFF theta input plots." << std::endl;

  // create DiFF calc plots
  {

    // create legend for blue DiFF calculation
    TLegend* lDiffThetaCalcB = new TLegend(0.1, 0.1, 0.3, 0.4, sHeader.Data());
    lDiffThetaCalcB -> SetFillColor(0);
    lDiffThetaCalcB -> SetLineColor(0);
    lDiffThetaCalcB -> SetTextFont(42);
    lDiffThetaCalcB -> SetTextAlign(12);
    lDiffThetaCalcB -> AddEntry(hThSpinBeamB, hThSpinBeamB -> GetTitle(), "p");
    lDiffThetaCalcB -> AddEntry(hThSpinPCB, hThSpinPCB -> GetTitle(), "p");
    lDiffThetaCalcB -> AddEntry(hThPCRCY, hThPCRCY -> GetTitle(), "p");
    lDiffThetaCalcB -> AddEntry(hThBeamRCY, hThBeamRCY -> GetTitle(), "p");
    lDiffThetaCalcB -> AddEntry(hThetaSB, hThetaSB -> GetTitle(), "p");
    lDiffThetaCalcB -> AddEntry(hThetaRC, hThetaRC -> GetTitle(), "p");
    lDiffThetaCalcB -> AddEntry(hThetaSBRC, hThetaSBRC -> GetTitle(), "p");

    // create plot for blue DiFF calculation
    TCanvas *cDiffThetaCalcB = new TCanvas("cDiffThetaCalcB", "", 750, 750);
    cDiffThetaCalcB -> SetGrid(0, 0);
    cDiffThetaCalcB -> cd();
    hThetaFrame     -> Draw();
    hThSpinBeamB    -> Draw("same");
    hThSpinPCB      -> Draw("same");
    hThPCRCY        -> Draw("same");
    hThBeamRCY      -> Draw("same");
    hThetaSB        -> Draw("same");
    hThetaRC        -> Draw("same");
    hThetaSBRC      -> Draw("same");
    lDiffThetaCalcB -> Draw();
    fOutput         -> cd();
    cDiffThetaCalcB -> Write();
    cDiffThetaCalcB -> Close();

    // create legend for blue DiFF calculation
    TLegend* lDiffThetaCalcY = new TLegend(0.1, 0.1, 0.3, 0.4, sHeader.Data());
    lDiffThetaCalcY -> SetFillColor(0);
    lDiffThetaCalcY -> SetLineColor(0);
    lDiffThetaCalcY -> SetTextFont(42);
    lDiffThetaCalcY -> SetTextAlign(12);
    lDiffThetaCalcY -> AddEntry(hThSpinBeamY, hThSpinBeamY -> GetTitle(), "p");
    lDiffThetaCalcY -> AddEntry(hThSpinPCY, hThSpinPCY -> GetTitle(), "p");
    lDiffThetaCalcY -> AddEntry(hThPCRCY, hThPCRCY -> GetTitle(), "p");
    lDiffThetaCalcY -> AddEntry(hThBeamRCY, hThBeamRCY -> GetTitle(), "p");
    lDiffThetaCalcY -> AddEntry(hThetaSY, hThetaSY -> GetTitle(), "p");
    lDiffThetaCalcB -> AddEntry(hThetaRC, hThetaRC -> GetTitle(), "p");
    lDiffThetaCalcY -> AddEntry(hThetaSYRC, hThetaSYRC -> GetTitle(), "p");

    // create plot for blue DiFF calculation
    TCanvas *cDiffThetaCalcY = new TCanvas("cDiffThetaCalcY", "", 750, 750);
    cDiffThetaCalcY -> SetGrid(0, 0);
    cDiffThetaCalcY -> cd();
    hThetaFrame     -> Draw();
    hThSpinBeamY    -> Draw("same");
    hThSpinPCY      -> Draw("same");
    hThPCRCY        -> Draw("same");
    hThBeamRCY      -> Draw("same");
    hThetaSY        -> Draw("same");
    hThetaRC        -> Draw("same");
    hThetaSYRC      -> Draw("same");
    lDiffThetaCalcY -> Draw();
    fOutput         -> cd();
    cDiffThetaCalcY -> Write();
    cDiffThetaCalcY -> Close();

  }  // end DiFF calc plot making
  std::cout << "    Created DiFF calculation plots." << std::endl;

  // save collins & BM histograms
  fOutput               -> cd();
  hInputPhiSpinB        -> Write();
  hInputPhiSpinY        -> Write();
  hInputXYSpinB         -> Write();
  hInputXYSpinY         -> Write();
  hInputPhiJet          -> Write();
  hInputThetaJet        -> Write();
  hInputCosThJet        -> Write();
  hInputZThJet          -> Write();
  hInputXYJet           -> Write();
  hInputPhiHad          -> Write();
  hInputThetaHad        -> Write();
  hInputCosThHad        -> Write();
  hInputXYHad           -> Write();
  hInputZThHad          -> Write();
  hCalcPhiJetBeamB      -> Write();
  hCalcThetaJetBeamB    -> Write();
  hCalcCosThJetBeamB    -> Write();
  hCalcZThJetBeamB      -> Write();
  hCalcPhiJetBeamY      -> Write();
  hCalcThetaJetBeamY    -> Write();
  hCalcCosThJetBeamY    -> Write();
  hCalcZThJetBeamY      -> Write();
  hCalcPhiHadJet        -> Write();
  hCalcThetaHadJet      -> Write();
  hCalcCosThHadJet      -> Write();
  hCalcZThHadJet        -> Write();
  hCalcJHBCrossMagB     -> Write(); 
  hCalcJHBDotB          -> Write();
  hCalcJHBCrossMagY     -> Write();
  hCalcJHBDotY          -> Write();
  hPhiSpinB             -> Write();
  hPhiSpinY             -> Write();
  hPhiHadB              -> Write();
  hPhiHadY              -> Write();
  hPhiHad2B             -> Write();
  hPhiHad2Y             -> Write();
  hPhiCollB             -> Write();
  hPhiCollY             -> Write();
  hPhiBoerB             -> Write();
  hPhiBoerY             -> Write();
  hAltPhiSpinB          -> Write();
  hAltPhiSpinY          -> Write();
  hAltPhiHadB           -> Write();
  hAltPhiHadY           -> Write();
  hAltPhiHad2B          -> Write();
  hAltPhiHad2Y          -> Write();
  hAltPhiCollB          -> Write();
  hAltPhiCollY          -> Write();
  hAltPhiBoerB          -> Write();
  hAltPhiBoerY          -> Write();
  hCheckPhiHadVsAltB    -> Write(); 
  hCheckPhiHadVsAltY    -> Write();
  hCheckPhiHadVsDotB    -> Write();
  hCheckPhiHadVsDotY    -> Write();
  hCheckAltPhiHadVsDotB -> Write();
  hCheckAltPhiHadVsDotY -> Write();

  // save DiFF histograms
  hInputPhiSec   -> Write();
  hInputThetaSec -> Write();
  hInputCosThSec -> Write();
  hInputXYSec    -> Write();
  hInputZThSec   -> Write();
  hThSpinBeamB   -> Write();
  hThSpinPCB     -> Write();
  hThSpinBeamY   -> Write();
  hThSpinPCY     -> Write();
  hThPCRCY       -> Write();
  hThBeamRCY     -> Write();
  hThetaSB       -> Write();
  hThetaSY       -> Write();
  hThetaRC       -> Write();
  hThetaSBRC     -> Write();
  hThetaSYRC     -> Write();
  std::cout << "    Saved histograms." << std::endl;

  // close output file
  fOutput -> cd();
  fOutput -> Close();

  // announce end and exit
  std::cout << "  Finished angle calculation test!\n" << std::endl;
  return;

}  // end AngleCalculationTest

// end ========================================================================
