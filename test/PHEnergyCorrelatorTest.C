/// ============================================================================
/*! \file    PHEnergyCorrelatorTest.C
 *  \authors Derek Anderson
 *  \date    09.21.2024
 *
 *  Macro to test compiling/running PHEnergyCorrelator library.
 */
/// ============================================================================

#define PHENERGYCORRELATORTEST_C

// c++ utilities
#include <iostream>
#include <utility>
#include <vector>
// root libraries
#include <TMath.h>
// analysis header
#include "../include/PHEnergyCorrelator.h"



// ============================================================================
//! Test macro for PHEnergyCorrelator library.
// ============================================================================
void PHEnergyCorrelatorTest() {

  // announce start
  std::cout << "\n  Beginning PHEnergyCorrelator test macro." << std::endl;

  // --------------------------------------------------------------------------
  // Test binning
  // --------------------------------------------------------------------------
  std::cout << "    Case [0]: test binning" << std::endl;

  // instantiate bin database
  PHEC::Bins bins;

  // loop over R_{L/S/M} bin edges and print them
  const std::vector<double> edges = bins.Get("side").GetBins();
  for (std::size_t iedge = 0; iedge < edges.size(); ++iedge) {
    std::cout << "      --- R_{L}: edge[" << iedge << "] = " << edges[iedge] << std::endl;
  }

  // --------------------------------------------------------------------------
  // Test histogram manager
  // --------------------------------------------------------------------------
  std::cout << "    Case [1]: test histogram manager" << std::endl;

  // pt jet bins
  std::vector< std::pair<float, float> > ptjetbins;
  ptjetbins.push_back( std::make_pair(5., 10.) );
  ptjetbins.push_back( std::make_pair(10., 15.) );
  ptjetbins.push_back( std::make_pair(15., 20.) );

  // cf jet bins
  std::vector< std::pair<float, float> > cfjetbins;
  cfjetbins.push_back( std::make_pair(0., 0.5) );
  cfjetbins.push_back( std::make_pair(0.5, 1.) );

  // instantiate calculator
  PHEC::Calculator calc(PHEC::Type::Pt);
  calc.SetPtJetBins(ptjetbins);
  calc.SetCFJetBins(cfjetbins);
  calc.SetHistTag("Test");

  // check no. of bins
  std::cout << "      --- N pt bins = " << calc.GetManager().GetNPtJetBins() << "\n"
            << "      --- N CF bins = " << calc.GetManager().GetNCFJetBins()
            << std::endl;

  // create histograms
  calc.Init(true);
  std::cout << "      --- N tags    = " << calc.GetManager().GetNIndexTags() << "\n"
            << "      --- N hists   = " << calc.GetManager().GetNHists()
            << std::endl;

  // --------------------------------------------------------------------------
  // Test calculator
  // --------------------------------------------------------------------------
  /* TODO changes values to something actually physical,
   *   these are just random...
   */
  std::cout << "    Case [2]: test calculator." << std::endl;

  // pi/4, 5pi/4
  const double piDiv4  = TMath::PiOver4();
  const double pi5Div4 = 5. * TMath::PiOver4();

  // jet values
  std::vector<PHEC::Type::Jet> jets;
  jets.push_back( PHEC::Type::Jet(0.25, 8.0, 0.2, piDiv4, 1.)   );
  jets.push_back( PHEC::Type::Jet(0.75, 13., -0.2, pi5Div4, 1.) );

  // cst values
  std::vector< std::vector< PHEC::Type::Cst > > csts( jets.size() );
  csts[0].push_back( PHEC::Type::Cst(0.25, 1.0, 0.25, piDiv4 + 0.02, 1.)   );
  csts[0].push_back( PHEC::Type::Cst(0.10, 0.3, 0.18, piDiv4 - 0.02, -1.)  );
  csts[1].push_back( PHEC::Type::Cst(0.57, 2.1, -0.16, pi5Div4 + 0.1, 1.)  );
  csts[1].push_back( PHEC::Type::Cst(0.20, 0.9, -0.21, pi5Div4 - 0.1, -1.) );
  csts[1].push_back( PHEC::Type::Cst(0.09, 0.1, -0.19, pi5Div4, 1.)        );

  // run calculations
  for (std::size_t ijet = 0; ijet < jets.size(); ++ijet) {
    for (std::size_t icst_a = 0; icst_a < csts[ijet].size(); ++icst_a) {
      for (std::size_t icst_b = 0; icst_b < csts[ijet].size(); ++icst_b) {

        // skip diagonal
        if (icst_a == icst_b) continue;

        // do calculation
        calc.CalcEEC(
          jets[ijet],
          std::make_pair(csts[ijet][icst_a], csts[ijet][icst_b])
        );
      }
    }
  }

  // --------------------------------------------------------------------------
  // Save histograms
  // --------------------------------------------------------------------------
  std::cout << "    Case [3]: test saving histograms" << std::endl;

  // create output file
  TFile* output = new TFile("test.root", "recreate");

  // save histograms to output
  calc.End(output);

  // --------------------------------------------------------------------------
  // Tests complete
  // --------------------------------------------------------------------------

  // close output file
  output -> cd();
  output -> Close();
  std::cout << "    Output file closed." << std::endl;

  // announce end & exit
  std::cout << "  PHEnergyCorrelator test complete!\n" << std::endl;
  return;

}

// end ========================================================================
