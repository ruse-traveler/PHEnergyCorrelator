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

  // check no. of bins
  std::cout << "      --- N pt bins = " << calc.GetManager().GetNPtJetBins() << "\n"
            << "      --- N CF bins = " << calc.GetManager().GetNCFJetBins()
            << std::endl;

  // create histograms
  calc.Init(true);
  std::cout << "      --- N tags    = " << calc.GetManager().GetNTags() << "\n"
            << "      --- N hists   = " << calc.GetManager().GetNHists()
            << std::endl;

  // --------------------------------------------------------------------------
  // Save histograms
  // --------------------------------------------------------------------------
  std::cout << "    Case [2]: test saving histograms" << std::endl;

  // create output file
  TFile* output = new TFile("test.root", "recreate");

  // save histograms to output
  calc.GetManager().SaveHists(output);

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
