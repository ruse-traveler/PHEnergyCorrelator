/// ============================================================================
/*! \file    PHCorrelatorCalculator.h
 *  \authors Derek Anderson, Alex Clarke
 *  \date    09.21.2024
 *
 *  Driver class to run n-point energy-energy correlator
 *  calculations on inputs.
 */
/// ============================================================================

#ifndef PHCORRELATORCALCULATOR_H
#define PHCORRELATORCALCULATOR_H

// c++ utilities
#include <algorithm>
#include <utility>
#include <vector>
// root libraries
#include <TLorentzVector.h>
#include <TVector3.h>
// analysis componenets
#include "PHCorrelatorManager.h"
#include "PHCorrelatorTypes.h"
#include "PHCorrelatorTools.h"



namespace PHEnergyCorrelator {

  // ==========================================================================
  //! ENC Calculator
  // ==========================================================================
  class Calculator {

    private:

      /* TODO
       *   - GetPtJetIndex()
       *   - GetCFJetIndex()
       *   - GetSpinIndex()
       *   - GetEnergyWeight()
       */

      // data members (options)
      bool m_do_eec_hist;
      bool m_do_e3c_hist;
      bool m_do_lec_hist;
      bool m_do_pt_bins;
      bool m_do_cf_bins;
      bool m_do_sp_bins;

      // data members (bins)
      std::vector< std::pair<float, float> > m_ptjet_bins;
      std::vector< std::pair<float, float> > m_cfjet_bins;
      std::vector< std::pair<float, float> > m_spin_bins;

      // data member (hist manager)
      Manager m_manager;

    public:

      /* TODO
       *   - DoEECCalc(jet, {cst, cst})
       *   - DoE3CCalc(jet, {cst, cst, cst})
       *   - DoLECCalc(jet, {lambda, cst})
       *   - FindLambda({cst...})
       */

      // ----------------------------------------------------------------------
      //! Getters
      // ----------------------------------------------------------------------
      Manager& GetManager() {return m_manager;}

      // ----------------------------------------------------------------------
      //! Set jet pt bins
      // ----------------------------------------------------------------------
      void SetPtJetBins(const std::vector< std::pair<float, float> >& bins) {

        // turn on pt binning
        m_do_pt_bins = true;

        // copy bins to member
        m_ptjet_bins.resize( bins.size() );
        std::copy(bins.begin(), bins.end(), m_ptjet_bins.begin());

        // update hist manager and exit
        m_manager.DoPtJetBins( m_ptjet_bins.size() );
        return;

      }  // end 'SetPtJetBins(std::vector<std::pair<float, float>>&)'

      // ----------------------------------------------------------------------
      //! Set jet CF bins
      // ----------------------------------------------------------------------
      void SetCFJetBins(const std::vector< std::pair<float, float> >& bins) {

        // turn on cf binning
        m_do_cf_bins = true;

        // copy bins to member
        m_cfjet_bins.resize( bins.size() );
        std::copy(bins.begin(), bins.end(), m_cfjet_bins.begin());

        // update hist manager and exit
        m_manager.DoCFJetBins( m_cfjet_bins.size() );
        return;

      }  // end 'SetCFJetBins(std::vector<std::pair<float, float>>&)'

      // ----------------------------------------------------------------------
      //! Set spin bins
      // ----------------------------------------------------------------------
      void SetSpinBins(const std::vector< std::pair<float, float> >& bins) {

        // turn on spin binning
        m_do_sp_bins = true;

        // copy bins to member
        m_spin_bins.resize( bins.size() );
        std::copy(bins.begin(), bins.end(), m_spin_bins.begin());

        // update hist manager and exit
        m_manager.DoSpinBins( m_spin_bins.size() );
        return;

      }  // end 'SetCFJetBins(std::vector<std::pair<float, float>>&)'

      // ----------------------------------------------------------------------
      //! Initialize calculator
      // ----------------------------------------------------------------------
      void Init(const bool do_eec, const bool do_e3c = false, const bool do_lec = false) {

        // turn on/off relevent hists
        m_do_eec_hist = do_eec;
        m_do_e3c_hist = do_e3c;
        m_do_lec_hist = do_lec;

        // generate histograms
        m_manager.DoEECHists(m_do_eec_hist);
        m_manager.DoE3CHists(m_do_e3c_hist);
        m_manager.DoLECHists(m_do_lec_hist);
        m_manager.GenerateHists();
        return;

      } // end 'Init(bool, bool, bool)'

      // ----------------------------------------------------------------------
      //! End calculations
      // ----------------------------------------------------------------------
      void End(TFile* file) {

        // save histograms to file
        m_manager.SaveHists(file);
        return;

      }  // end 'End(TFile*)'

      // ----------------------------------------------------------------------
      //! default ctor
      // ----------------------------------------------------------------------
      Calculator()  {

        m_do_eec_hist = false;
        m_do_e3c_hist = false;
        m_do_lec_hist = false;
        m_do_pt_bins  = false;
        m_do_cf_bins  = false;
        m_do_sp_bins  = false;

      }  // end default ctor

      // ----------------------------------------------------------------------
      //! default dtor
      // ----------------------------------------------------------------------
      ~Calculator() {};

  };  // end PHEnergyCorrelator::Calculator

}  // end PHEnergyCorrelator namespace

#endif

// end ========================================================================
