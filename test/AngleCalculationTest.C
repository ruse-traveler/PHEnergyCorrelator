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
void AngleCalculationTest(
  const std::string oFile = "angleCalcTest.nIter10K_uniformSphereWithXYPlots_doWrap.d19m2y2025.root",
  const std::size_t nIter = 10000,
  const bool doWrap = true,
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
  const int   nXYBins   = 400;
  const float xAngStart = -12.60;
  const float xXYStart  = -2.0;
  const float xAngStop  = 12.60; 
  const float xXYStop   = 2.0;

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
  TH2D* hInputXYJet        = new TH2D("hInputXYJet", "(x,y) sampled (jet)", nXYBins, xXYStart, xXYStop, nXYBins, xXYStart, xXYStop);
  TH1D* hInputPhiHad       = new TH1D("hInputPhiHad", "#phi_{h} input", nAngBins, xAngStart, xAngStop);
  TH1D* hInputThetaHad     = new TH1D("hInputThetaHad", "#theta_{h} input", nAngBins, xAngStart, xAngStop);
  TH2D* hInputXYHad        = new TH2D("hInputXYHad", "(x,y) sampled (had)", nXYBins, xXYStart, xXYStop, nXYBins, xXYStart, xXYStop);
  TH1D* hCalcPhiJetBeamB   = new TH1D("hCalcPhiJetBeamB", "#phi_{jet-beam}^{B}", nAngBins, xAngStart, xAngStop);
  TH1D* hCalcThetaJetBeamB = new TH1D("hCalcThetaJetBeamB", "#theta_{jet-beam}^{B}", nAngBins, xAngStart, xAngStop);
  TH1D* hCalcPhiJetBeamY   = new TH1D("hCalcPhiJetBeamY", "#phi_{jet-beam}^{Y}", nAngBins, xAngStart, xAngStop);
  TH1D* hCalcThetaJetBeamY = new TH1D("hCalcThetaJetBeamY", "#theta_{jet-beam}^{Y}", nAngBins, xAngStart, xAngStop);
  TH1D* hCalcPhiHadJet     = new TH1D("hCalcPhiJetHad", "#phi_{jet-h}", nAngBins, xAngStart, xAngStop);
  TH1D* hCalcThetaHadJet   = new TH1D("hCalcThetaJetHad", "#theta_{jet-h}", nAngBins, xAngStart, xAngStop);
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
    TVector3 vecSpinY3(xRandSpinY, yRandSpinB, 0.0);
    hInputPhiSpinB -> Fill( vecSpinB3.Phi() );
    hInputPhiSpinY -> Fill( vecSpinY3.Phi() );
    hInputXYSpinB  -> Fill( vecSpinB3.X(), vecSpinB3.Y() );
    hInputXYSpinY  -> Fill( vecSpinB3.X(), vecSpinB3.Y() );

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
    hInputXYJet    -> Fill( vecJet3.X(), vecJet3.Y() );
    hInputPhiHad   -> Fill( vecHad3.Phi() );
    hInputThetaHad -> Fill( vecHad3.Theta() );
    hInputXYHad    -> Fill( vecHad3.X(), vecHad3.Y() );

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
    hCalcPhiJetBeamY   -> Fill( normJetBeam3.second.Phi() );
    hCalcThetaJetBeamY -> Fill( normJetBeam3.second.Theta() );

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

    // (4) get phiHadron: angle between the jet-beam plane and the
    //   - angle between jet-hadron plane
    //   - constrain to range [0,2pi)
    double phiHadBlue = atan2( normJetBeam3.first.Cross(normHadJet3).Mag(), normJetBeam3.first.Dot(normHadJet3) );
    double phiHadYell = atan2( normJetBeam3.second.Cross(normHadJet3).Mag(), normJetBeam3.second.Dot(normHadJet3) );
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
    //   - constrain to [0, 2pi)
    double phiCollBlue = phiSpinBlue - phiHadBlue;
    double phiCollYell = phiSpinYell - phiHadYell;
      if (doWrap) {
      if (phiCollBlue < 0)               phiCollBlue += TMath::TwoPi();
      if (phiCollBlue >= TMath::TwoPi()) phiCollBlue -= TMath::TwoPi();
      if (phiCollYell < 0)               phiCollYell += TMath::TwoPi();
      if (phiCollYell >= TMath::TwoPi()) phiCollYell -= TMath::TwoPi();
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

  }  // end iter loop
  std::cout << "    MC loop finished!" << std::endl;

  // normalize histograms
  //   - n.b. (x,y) hists NOT normalized
  hInputPhiSpinB     -> Scale(1. / (double) nIter);
  hInputPhiSpinY     -> Scale(1. / (double) nIter);
  hInputPhiJet       -> Scale(1. / (double) nIter);
  hInputThetaJet     -> Scale(1. / (double) nIter);
  hInputPhiHad       -> Scale(1. / (double) nIter);
  hInputThetaHad     -> Scale(1. / (double) nIter);
  hCalcPhiJetBeamB   -> Scale(1. / (double) nIter);
  hCalcThetaJetBeamB -> Scale(1. / (double) nIter);
  hCalcPhiJetBeamY   -> Scale(1. / (double) nIter);
  hCalcThetaJetBeamY -> Scale(1. / (double) nIter);
  hCalcPhiHadJet     -> Scale(1. / (double) nIter);
  hCalcThetaHadJet   -> Scale(1. / (double) nIter);
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
  std::cout << "    Normalized histograms." << std::endl;

  // create frame histograms
  TH1D* hPhiFrame   = hInputPhiSpinB -> Clone();
  TH1D* hThetaFrame = hInputThetaJet -> Clone();
  hPhiFrame   -> Reset("ICES");
  hPhiFrame   -> SetName("hPhiFrame");
  hPhiFrame   -> SetTitle(";#phi [rad]");
  hThetaFrame -> Reset("ICES");
  hThetaFrame -> SetName("hThetaFrame");
  hThetaFrame -> SetTitle(";#theta [rad]");

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
  hInputPhiHad       -> SetLineColor(col[2]);
  hInputPhiHad       -> SetMarkerColor(col[2]);
  hInputPhiHad       -> SetMarkerStyle(mar[2]);
  hInputThetaHad     -> SetLineColor(col[2]);
  hInputThetaHad     -> SetMarkerColor(col[2]);
  hInputThetaHad     -> SetMarkerStyle(mar[2]);
  hCalcPhiJetBeamB   -> SetLineColor(col[3]);
  hCalcPhiJetBeamB   -> SetMarkerColor(col[3]);
  hCalcPhiJetBeamB   -> SetMarkerStyle(mar[3]);
  hCalcThetaJetBeamB -> SetLineColor(col[3]);
  hCalcThetaJetBeamB -> SetMarkerColor(col[3]);
  hCalcThetaJetBeamB -> SetMarkerStyle(mar[3]);
  hCalcPhiJetBeamY   -> SetLineColor(col[3]);
  hCalcPhiJetBeamY   -> SetMarkerColor(col[3]);
  hCalcPhiJetBeamY   -> SetMarkerStyle(mar[3]);
  hCalcThetaJetBeamY -> SetLineColor(col[3]);
  hCalcThetaJetBeamY -> SetMarkerColor(col[3]);
  hCalcThetaJetBeamY -> SetMarkerStyle(mar[3]);
  hCalcPhiHadJet     -> SetLineColor(col[4]);
  hCalcPhiHadJet     -> SetMarkerColor(col[4]);
  hCalcPhiHadJet     -> SetMarkerStyle(mar[4]);
  hCalcThetaHadJet   -> SetLineColor(col[4]);
  hCalcThetaHadJet   -> SetMarkerColor(col[4]);
  hCalcThetaHadJet   -> SetMarkerStyle(mar[4]);
  hPhiSpinB          -> SetLineColor(col[5]);
  hPhiSpinB          -> SetMarkerColor(col[5]);
  hPhiSpinB          -> SetMarkerStyle(mar[5]);
  hPhiSpinY          -> SetLineColor(col[5]);
  hPhiSpinY          -> SetMarkerColor(col[5]);
  hPhiSpinY          -> SetMarkerStyle(mar[5]);
  hPhiHadB           -> SetLineColor(col[6]);
  hPhiHadB           -> SetMarkerColor(col[6]);
  hPhiHadB           -> SetMarkerStyle(mar[6]);
  hPhiHadY           -> SetLineColor(col[6]);
  hPhiHadY           -> SetMarkerColor(col[6]);
  hPhiHadY           -> SetMarkerStyle(mar[6]);
  hPhiHad2B          -> SetLineColor(col[7]);
  hPhiHad2B          -> SetMarkerColor(col[7]);
  hPhiHad2B          -> SetMarkerStyle(mar[7]);
  hPhiHad2Y          -> SetLineColor(col[7]);
  hPhiHad2Y          -> SetMarkerColor(col[7]);
  hPhiHad2Y          -> SetMarkerStyle(mar[7]);
  hPhiCollB          -> SetLineColor(col[8]);
  hPhiCollB          -> SetMarkerColor(col[8]);
  hPhiCollB          -> SetMarkerStyle(mar[8]);
  hPhiCollY          -> SetLineColor(col[8]);
  hPhiCollY          -> SetMarkerColor(col[8]);
  hPhiCollY          -> SetMarkerStyle(mar[8]);
  hPhiBoerB          -> SetLineColor(col[9]);
  hPhiBoerB          -> SetMarkerColor(col[9]);
  hPhiBoerB          -> SetMarkerStyle(mar[9]);
  hPhiBoerY          -> SetLineColor(col[9]);
  hPhiBoerY          -> SetMarkerColor(col[9]);
  hPhiBoerY          -> SetMarkerStyle(mar[9]);
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
  hInputXYJet        -> Write();
  hInputPhiHad       -> Write();
  hInputThetaHad     -> Write();
  hInputXYHad        -> Write();
  hCalcPhiJetBeamB   -> Write();
  hCalcThetaJetBeamB -> Write();
  hCalcPhiJetBeamY   -> Write();
  hCalcThetaJetBeamY -> Write();
  hCalcPhiHadJet     -> Write();
  hCalcThetaHadJet   -> Write();
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
  std::cout << "    Saved histograms." << std::endl;

  // close output file
  fOutput -> cd();
  fOutput -> Close();

  // announce end and exit
  std::cout << "  Finished angle calculation test!\n" << std::endl;
  return;

}  // end AngleCalculationTest

// end ========================================================================
