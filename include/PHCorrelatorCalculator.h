/// ============================================================================
/*! \file    PHCorrelatorCalculator.h
 *  \authors Derek Anderson, Alex Clarke
 *  \date    09.21.2024
 *
 *  Cclass to run n-point energy-energy correlator
 *  calculations on inputs.
 */
/// ============================================================================

#ifndef PHCORRELATORCALCULATOR_H
#define PHCORRELATORCALCULATOR_H

// c++ utilities
#include <utility>
#include <vector>
// root libraries
#include <TLorentzVector.h>
#include <TVector3.h>



namespace PHEnergyCorrelator {

  // ==========================================================================
  //! ENC Calculator
  // ==========================================================================
  class Calculator {

    private:

      // data members
      std::vector< std::pair<float, float> > m_ptjet_bins;
      std::vector< std::pair<float, float> > m_cfjet_bins;
      std::vector< std::pair<float, float> > m_spin_bins;

      /* TODO
       *   - add flags
       *   - GetPtJetIndex()
       *   - GetCFJetIndex()
       *   - GetSpinIndex()
       *   - GetEnergyWeight()
       */

    public:

      // default ctor/dtor
      Calculator()  {};
      ~Calculator() {};

      /* TODO
       *   - DoEECCalc(jet, {cst, cst})
       *   - DoE3CCalc(jet, {cst, cst, cst})
       *   - DoLECCalc(jet, {lambda, cst})
       *   - FindLambda({cst...})
       */

  };  // end PHEnergyCorrelator::Calculator

}  // end PHEnergyCorrelator namespace

#endif

// end ========================================================================
