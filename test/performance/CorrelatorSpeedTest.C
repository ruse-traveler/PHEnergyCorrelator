/// ============================================================================
/*! \file    CorrelatorSpeedTest.C
 *  \authors Derek Anderson
 *  \date    02.24.2024
 *
 *  Macro to check speed of ENC calculations
 */
/// ============================================================================

#define CORRELATORSPEEDTEST_C

// c++ utilities
#include <iostream>
#include <utility>
#include <vector>
// root libraries
#include <TDatime.h>
#include <TMath.h>
#include <TRandom.h>
#include <TRandom3.h>
#include <TStopwatch.h>
// analysis header
#include "../../include/PHEnergyCorrelator.h"



// ============================================================================
//! Test speed of PHEnergyCorrelator library.
// ============================================================================
void CorrelatorSpeedTest(
  const std::string outfile = "speedTest.change3_fixHasTypo.nIter10KnJet1nCst3.d6m3y2025.root",
  const std::size_t nIter = 10000,
  const std::size_t nJet = 1,
  const std::size_t nCst = 3,
  const bool doBatch = true
) {

  // announce start
  std::cout << "\n  Beginning correlator speed test." << std::endl;

  // jet ranges for rng
  const double ptJetRNG[2]   = {3.0, 40.0};
  const double cfJetRNG[2]   = {0.3, 0.9};
  const double chJetRNG[2]   = {-20., 20.};
  const double etaJetRNG[2]  = {-0.5, 0.5};
  const double phiJetRNG[2]  = {-TMath::Pi(), TMath::Pi()};
  const int    spinJetRNG[2] = {0, 6};

  // cst ranges for rng
  const double zCstRNG[2]   = {0.0, 1.0};
  const double jtCstRNG[2]  = {0.1, 20.0};
  const double chCstRNG[2]  = {-1.0, 1.0};
  const double etaCstRNG[2] = {-0.5, 0.5};
  const double phiCstRNG[2] = {-TMath::Pi(), TMath::Pi()};

  // --------------------------------------------------------------------------
  // Set up calculation
  // --------------------------------------------------------------------------

  // pt jet bins
  std::vector< std::pair<float, float> > ptjetbins;
  ptjetbins.push_back( std::make_pair(5., 10.) );
  ptjetbins.push_back( std::make_pair(10., 15.) );
  ptjetbins.push_back( std::make_pair(15., 20.) );

  // charge jet bins
  std::vector< std::pair<float, float> > chjetbins;
  chjetbins.push_back( std::make_pair(-100., 0.0) );
  chjetbins.push_back( std::make_pair(0.0, 100.) );

  // instantiate calculator
  PHEC::Calculator calc(PHEC::Type::Pt);
  calc.SetPtJetBins(ptjetbins);
  calc.SetChargeBins(chjetbins);
  calc.SetDoSpinBins(true);
  calc.Init(true);
  std::cout << "    Initialized calculator." << std::endl;

  // initialize rng
  TDatime*  time  = new TDatime();
  TRandom3* rando = new TRandom3();
  rando -> SetSeed(time -> Get());
  std::cout << "    Initialized RNG." << std::endl;

  // --------------------------------------------------------------------------
  // Create nIter fake events
  // --------------------------------------------------------------------------

  // start timer
  TStopwatch* watch = new TStopwatch();
  watch -> Start();
  std::cout << "    MC loop: running " << nIter << " iterations:"
            << std::endl;

  // mc loop
  std::vector<PHEC::Type::Cst> csts;
  for (std::size_t iIter = 0; iIter < nIter; ++iIter) {

    // print progress every 100 iterations
    std::size_t iProg = iIter + 1;
    if (iIter % 100 == 0) {
      if (doBatch || (iProg + 100 >= nIter)) {
        std::cout << "      iter " << iProg << "..." << std::endl;
      } else {
        std::cout << "      iter " << iProg << "...\r" << std::flush;
      }
    }

    // generate nJet jets each with nCst cst.s
    for (std::size_t iJet = 0; iJet < nJet; ++iJet) {

      // generate random jet
      PHEC::Type::Jet jet(
        rando -> Uniform(cfJetRNG[0], cfJetRNG[1]),
        rando -> Uniform(ptJetRNG[0], ptJetRNG[1]),
        rando -> Uniform(etaJetRNG[0], etaJetRNG[1]),
        rando -> Uniform(phiJetRNG[0], phiJetRNG[1]),
        rando -> Uniform(chJetRNG[0], chJetRNG[1]),
        rando -> Uniform(spinJetRNG[0], spinJetRNG[1])
      );

      // make sure cst vector is clear
      csts.clear();

      // generate random cst.s
      for (std::size_t iCst = 0; iCst < nCst; ++iCst) {
        csts.push_back(
          PHEC::Type::Cst(
            rando -> Uniform(zCstRNG[0], zCstRNG[1]),
            rando -> Uniform(jtCstRNG[0], jtCstRNG[1]),
            rando -> Uniform(etaCstRNG[0], etaCstRNG[1]),
            rando -> Uniform(phiCstRNG[0], phiCstRNG[1]),
            rando -> Uniform(chCstRNG[0], chCstRNG[1])
          )
        );
      }

      // now run calculation
      for (std::size_t iCstA = 0; iCstA < nCst; ++iCstA) {
        for (std::size_t iCstB = 0; iCstB <= iCstA; ++iCstB) {
          calc.CalcEEC(
            jet,
            std::make_pair(csts[iCstA], csts[iCstB])
          );
        }  // end cst B loop 
      }  // end cst A loop
    }  // end jet loop
  }  // end mc loop

  // stop timer
  watch -> Stop();
  std::cout << "    MC loop finished! Time elapsed:\n"
            << "    -------------------------------\n"
            << std::endl;

  // print time
  watch -> Print();
  std::cout << "\n    -------------------------------" << std::endl;

  // --------------------------------------------------------------------------
  // Save output and end
  // --------------------------------------------------------------------------

  // create output file
  TFile* output = new TFile(outfile.data(), "recreate");

  // save histograms to output
  calc.End(output);
  std::cout << "    Saved output." << std::endl;

  // announce end & exit
  std::cout << "  Correlator speed test complete!\n" << std::endl;
  return;

}

// end ========================================================================
