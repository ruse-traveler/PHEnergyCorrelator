/// ============================================================================
/*! \file    PHAngleCalculatorTest.C
 *  \authors Derek Anderson
 *  \date    08.01.2025
 *
 *  Macro for unit test of angle calculation class. Intended
 *  to confirm that calculations are identical.
 */
/// ============================================================================

#define PHANGLECALCULATORTEST_C

// c++ utilities
#include <iostream>
#include <string>
#include <utility>
#include <vector>
// root libraries
#include <TMath.h>
#include <TRandom.h>
#include <TRandom3.h>
// analysis header (adjust as needed)
#include "../PHEnergyCorrelator/include/PHEnergyCorrelator.h"
//#include "../PHECJohnFork/include/PHEnergyCorrelator.h"



// ============================================================================
//! Macro for angle calculator unit test
// ============================================================================
void PHAngleCalculatorTest(
  const std::string outfile = "angleUnitTest.factorAngles_run0_addPhiSpinCalc.d1m8y2025.root"
) {

  // announce start
  std::cout << "\n  Beginning angle calculation unit test." << std::endl;

  // --------------------------------------------------------------------------
  // Set up calculation
  // --------------------------------------------------------------------------

  // calculation parameters
  //   - fixed for reproducibility
  const std::size_t nIter   = 1000;
  const std::size_t nJet    = 1;
  const std::size_t nCst    = 3;
  const bool        doBatch = true;

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
  //   - use fixed seed for reproducibility
  TRandom3* rando = new TRandom3();
  rando -> SetSeed(1);
  std::cout << "    Initialized RNG." << std::endl;

  // --------------------------------------------------------------------------
  // Create nIter fake events
  // --------------------------------------------------------------------------

  // announce start
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

  // announce end 
  std::cout << "    MC loop finished!" << std::endl;

  // --------------------------------------------------------------------------
  // Save output and end
  // --------------------------------------------------------------------------

  // create output file
  TFile* output = new TFile(outfile.data(), "recreate");

  // save histograms to output
  calc.End(output);
  std::cout << "    Saved output." << std::endl;

  // announce end & exit
  std::cout << "  Angle calculation unit test complete!\n" << std::endl;
  return;

}

// end ========================================================================
