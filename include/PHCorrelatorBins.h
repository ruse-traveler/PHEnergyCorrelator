/// ============================================================================
/*! \file    PHCorrelatorBins.h
 *  \authors Derek Anderson
 *  \date    09.24.2024
 *
 *  Class to centralize binning of histograms filed
 *  during ENC calculations.
 */
/// ============================================================================

#ifndef PHCORRELATORBINS_H
#define PHCORRELATORBINS_H

// c++ utilities
#include <map>
#include <string>
// analysis components
#include "PHCorrelatorBinning.h"
#include "PHCorrelatorTools.h"
#include "PHCorrelatorTypes.h"



namespace PHEnergyCorrelator {

  // ==========================================================================
  //! Bin Database
  // ==========================================================================
  /*! A class to centralize binning for various quantities like R_{L}, etc.
   *  Methods are provided to update exisiting/add new bin definitions
   *  on the fly.
   */
  class Bins {

    private:

      // data members
      std::map<std::string, Binning> m_bins;

    public:

      // ----------------------------------------------------------------------
      //! Add a binning
      // ----------------------------------------------------------------------
      void Add(const std::string& name, const Binning& binning) {

        // throw error if binning already exists
        if (m_bins.count(name) >= 1) {
          assert(m_bins.count(name) == 0);
        }

        // otherwise insert new binning
        m_bins[name] = binning;
        return;

      }  // end 'Add(std::string&, Binning&)'

      // ----------------------------------------------------------------------
      //! Change a binning
      // ----------------------------------------------------------------------
      void Set(const std::string& variable, const Binning& binning) {

        // throw error if binning doesn't exist
        if (m_bins.count(variable) == 0) {
          assert(m_bins.count(variable) >= 1);
        }

        // otherwise update binning
        m_bins[variable] = binning;
        return;

      }  // end 'Set(std::string&, Binning&)'

      // ----------------------------------------------------------------------
      //! Get a binning
      // ----------------------------------------------------------------------
      Binning Get(const std::string& variable) {

        // throw error if binning doesn't exist
        if (m_bins.count(variable) == 0) {
          assert(m_bins.count(variable) >= 1);
        }

        // otherwise return binning
        return m_bins[variable];

      }  // end 'Get(std::string&)'

      // ----------------------------------------------------------------------
      //! default ctor
      // ----------------------------------------------------------------------
      Bins() {

        m_bins["energy"]   = Binning(202, -1., 100.);
        m_bins["side"]     = Binning(75, 1e-5, 1., Type::Log);
        m_bins["logside"]  = Binning(75, -5., 0.);
        m_bins["angle"]    = Binning(180, -3.15, 3.15);
        m_bins["cosangle"] = Binning(20, -1., 1.);
        m_bins["xi"]       = Binning(100, 0., 1.);
        m_bins["pattern"]  = Binning(10, -0.5, 9.5);

      }  // end 'ctor()'

      // ----------------------------------------------------------------------
      //! default dtor
      // ----------------------------------------------------------------------
      ~Bins() {};

  };  // end Bins

}  // end PHEnergyCorrelator namespace

#endif

// end ========================================================================
