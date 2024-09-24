/// ============================================================================
/*! \file    PHCorrelatorBins.h
 *  \authors Derek Anderson
 *  \date    09.24.2024
 *
 *  Classes to handle binning of histograms filed
 *  during ENC calculations.
 */
/// ============================================================================

#ifndef PHCORRELATORBINS_H
#define PHCORRELATORBINS_H

// c++ utilities
#include <map>
#include <string>
#include <vector>
// root libraries
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
// analysis components
#include "PHCorrelatorTools.h"
#include "PHCorrelatorTypes.h"



namespace PHEnergyCorrelator {

  // ==========================================================================
  //! Binning definition
  // ==========================================================================
  /*! A small class to consolidate data
   *  for defining histogram bins.
   */ 
  class Binning {

    private:

      // data members
      uint32_t            m_num;
      double              m_start;
      double              m_stop;
      std::vector<double> m_bins;

    public:

      // ----------------------------------------------------------------------
      //! Uniform bin getters
      //-----------------------------------------------------------------------
      uint32_t GetNum()   const {return m_num;}
      double   GetStart() const {return m_start;}
      double   GetStop()  const {return m_stop;}

      // ----------------------------------------------------------------------
      //! Variable bin getter
      // ----------------------------------------------------------------------
      std::vector<double> GetBins() const {return m_bins;}

      // ----------------------------------------------------------------------
      //! default ctor/dtor
      // ----------------------------------------------------------------------
      Binning()  {};
      ~Binning() {};

      // ----------------------------------------------------------------------
      //! ctor accepting uniform parameters
      // ----------------------------------------------------------------------
      /* TODO add flag to set log bins vs. not */ 
      Binning(const uint32_t num, const double start, const double stop) {

        m_num   = num;
        m_start = start;
        m_stop  = stop;

      }  // end ctor(uint32_t, double, double)

      // ----------------------------------------------------------------------
      //! ctor accepting non-uniform parameters
      // ----------------------------------------------------------------------
      Binning(const std::vector<double> edges) {

        m_bins  = edges;
        m_num   = edges.size() - 1;
        m_start = edges.front();
        m_stop  = edges.back();

      }  // end ctor(std::vector<double>)

  };  // end Binning



  // ==========================================================================
  //! Bin Database
  // ==========================================================================
  /*! A class to centralize binning for various quantities like R_{L}, etc.
   *  Methods are provided to update exisiting/add new bin definitions
   *  on the fly.
   */
  class Bins {

    private:

      /* TODO
       *  - Add R_{L/S/M}
       *  - Add cos(angle)
       *  - Add angle
       *  - Add xi
       */
      std::map<std::string, Binning> m_bins = {
        { "energy", {202, -1., 100.} }
      };

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
        m_bins.insert( {name, binning} );
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
          assert(m_bins.count(variable) == 0);
        }

        // otherwise return binning
        return m_bins[variable];

      }  // end 'Get(std::string&)'

      // ----------------------------------------------------------------------
      //! default ctor/dtor
      // ----------------------------------------------------------------------
      Bins()  {};
      ~Bins() {};

  };  // end Bins

}  // end PHEnergyCorrelator namespace

#endif

// end ========================================================================
