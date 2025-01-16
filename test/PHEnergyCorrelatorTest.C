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
  std::cout << "      --- [PASS] tested R_{L} binning" << std::endl;

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
  cfjetbins.push_back( std::make_pair(0.3, 0.6) );

  // charge jet bins
  std::vector< std::pair<float, float> > chjetbins;
  chjetbins.push_back( std::make_pair(-100., 0.0) );
  chjetbins.push_back( std::make_pair(0.0, 100.) );

  // instantiate calculator
  PHEC::Calculator calc_a(PHEC::Type::Pt);
  calc_a.SetPtJetBins(ptjetbins);
  calc_a.SetCFJetBins(cfjetbins);
  calc_a.SetChargeBins(chjetbins);
  calc_a.SetHistTag("FirstCalculation");

  // check no. of bins
  std::cout << "      --- N pt bins = " << calc_a.GetManager().GetNPtJetBins() << "\n"
            << "      --- N CF bins = " << calc_a.GetManager().GetNCFJetBins() << "\n"
            << "      --- N ch bins = " << calc_a.GetManager().GetNChargeBins()
            << std::endl;

  // create histograms
  calc_a.Init(true);
  std::cout << "      --- N tags    = " << calc_a.GetManager().GetNIndexTags() << "\n"
            << "      --- N hists   = " << calc_a.GetManager().GetNHists()     << "\n"
            << "      --- [PASS] tested histogram manager"
            << std::endl;

  // --------------------------------------------------------------------------
  // Prepare dummy jet/cst values
  // --------------------------------------------------------------------------
  /*! N.B. the jet/cst values are random: they're here just to
   *  provide a unit test and make sure everything works.
   */
  std::cout << "    Preparing dummy jet/cst values for tests" << std::endl;

  // pi/4, pi/3, 5pi/4, 4pi/3
  const double piDiv4  = TMath::PiOver4();
  const double piDiv3  = TMath::Pi() / 3.;
  const double pi5Div4 = 5. * TMath::PiOver4();
  const double pi4Div3 = 4. * piDiv3;

  // jet values (enough to cover all spin patterns)
  // jet ctor arguments:
  //   0 --> cf
  //   1 --> pt
  //   2 --> eta
  //   3 --> phi
  //   4 --> charge
  //   5 --> spin pattern
  std::vector<PHEC::Type::Jet> jets;
  jets.push_back(
    PHEC::Type::Jet(0.40, 9.0, 0.75, piDiv3, -0.5, 0)
  );
  jets.push_back(
    PHEC::Type::Jet(0.25, 8.0, 0.2, piDiv4, 0.5, 1)
  );
  jets.push_back(
    PHEC::Type::Jet(0.75, 13., -0.2, pi5Div4, -2.3, 2)
  );
  jets.push_back(
    PHEC::Type::Jet(0.90, 5.0, 0.1, piDiv3, 3.2, 3)
  );
  jets.push_back(
    PHEC::Type::Jet(1.0, 7.0, -0.05, pi4Div3, -1.7, 4)
  );
  jets.push_back(
    PHEC::Type::Jet(0.50, 3.0, -0.75, pi4Div3, 1.3, 5)
  );
  jets.push_back(
    PHEC::Type::Jet(0.66, 4.0, -0.05, pi4Div3, 0.4, 6)
  );

  // cst values (enough to make sure calculations run)
  // cst ctor arguments:
  //   0 --> z
  //   1 --> jt
  //   2 --> eta
  //   3 --> phi
  //   4 --> chrg
  std::vector< std::vector< PHEC::Type::Cst > > csts( jets.size() );
  for (std::size_t ijet = 0; ijet < jets.size(); ++ijet) {
    if ((ijet % 2) == 0) {
      csts[ijet].push_back( PHEC::Type::Cst(0.25, 1.0, 0.25, piDiv4 + 0.02, 1.)   );
      csts[ijet].push_back( PHEC::Type::Cst(0.10, 0.3, 0.18, piDiv4 - 0.02, -1.)  );
    } else {
      csts[ijet].push_back( PHEC::Type::Cst(0.57, 2.1, -0.16, pi5Div4 + 0.1, 1.)  );
      csts[ijet].push_back( PHEC::Type::Cst(0.20, 0.9, -0.21, pi5Div4 - 0.1, -1.) );
      csts[ijet].push_back( PHEC::Type::Cst(0.09, 0.1, -0.19, pi5Div4, 1.)        );
    }
  }
  std::cout << "    Dummy jet/cst values prepared" << std::endl;

  // --------------------------------------------------------------------------
  // Test calculator
  // --------------------------------------------------------------------------
  std::cout << "    Case [2]: test calculator" << std::endl;

  // run calculations
  for (std::size_t ijet = 0; ijet < jets.size(); ++ijet) {
    for (std::size_t icst_a = 0; icst_a < csts[ijet].size(); ++icst_a) {
      for (std::size_t icst_b = 0; icst_b < csts[ijet].size(); ++icst_b) {
        calc_a.CalcEEC(
          jets[ijet],
          std::make_pair(csts[ijet][icst_a], csts[ijet][icst_b])
        );
      }
    }
  }
  std::cout << "      --- [PASS] ran 1st calculation" << std::endl;

  // --------------------------------------------------------------------------
  // Test adding 2nd calculator
  // --------------------------------------------------------------------------
  std::cout << "    Case [3]: test adding a 2nd calculator" << std::endl;

  // instantiate calculator
  PHEC::Calculator calc_b(PHEC::Type::Pt);
  calc_b.SetPtJetBins(ptjetbins);
  calc_b.SetCFJetBins(cfjetbins);
  calc_b.SetChargeBins(chjetbins);
  calc_b.SetHistTag("SecondCalculation");
  calc_b.Init(true);

  // run calculations
  for (std::size_t ijet = 0; ijet < jets.size(); ++ijet) {
    for (std::size_t icst_a = 0; icst_a < csts[ijet].size(); ++icst_a) {
      for (std::size_t icst_b = 0; icst_b < csts[ijet].size(); ++icst_b) {

        // do calculation
        calc_b.CalcEEC(
          jets[ijet],
          std::make_pair(csts[ijet][icst_a], csts[ijet][icst_b])
        );
      }
    }
  }
  std::cout << "      --- [PASS] ran second calculation" << std::endl;

  // --------------------------------------------------------------------------
  // Test spin sorting
  // --------------------------------------------------------------------------
  std::cout << "    Case [4]: test spin sorting" << std::endl;

  // instantiate calculator
  PHEC::Calculator calc_c(PHEC::Type::Pt);
  calc_c.SetPtJetBins(ptjetbins);
  calc_c.SetCFJetBins(cfjetbins);
  calc_c.SetChargeBins(chjetbins);
  calc_c.SetDoSpinBins(true);
  calc_c.SetHistTag("ThirdCalculation");
  calc_c.Init(true);
  std::cout << "      --- initialized third calculation" << std::endl;

  // run calculations
  for (std::size_t ijet = 0; ijet < jets.size(); ++ijet) {
    for (std::size_t icst_a = 0; icst_a < csts[ijet].size(); ++icst_a) {
      for (std::size_t icst_b = 0; icst_b < csts[ijet].size(); ++icst_b) {

        // do calculation
        calc_c.CalcEEC(
          jets[ijet],
          std::make_pair(csts[ijet][icst_a], csts[ijet][icst_b])
        );
      }
    }
    std::cout << "      --- ran calculation for spin pattern " << jets[ijet].pattern << std::endl;
  }
  std::cout << "      --- [PASS] ran third calculation" << std::endl;

  // --------------------------------------------------------------------------
  // Save histograms
  // --------------------------------------------------------------------------
  std::cout << "    Case [5]: test saving histograms" << std::endl;

  // create output file
  TFile* output = new TFile("test.root", "recreate");

  // save histograms to output
  calc_a.End(output);
  calc_b.End(output);
  calc_c.End(output);
  std::cout << "      --- [PASS] histograms saved" << std::endl;

  // --------------------------------------------------------------------------
  // Tests complete
  // --------------------------------------------------------------------------

  // close output file
  output -> cd();
  output -> Close();
  std::cout << "    Output file closed" << std::endl;

  // announce end & exit
  std::cout << "  PHEnergyCorrelator test complete!\n" << std::endl;
  return;

}

// end ========================================================================
