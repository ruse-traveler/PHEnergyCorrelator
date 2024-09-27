/// ============================================================================
/*! \file    PHCorrelatorConstants.h
 *  \authors Derek Anderson
 *  \date    09.27.2024
 *
 *  Namespace to collect constants in PHENIX
 *  ENC analyses.
 */
/// ============================================================================

#ifndef PHCORRELATORCONSTANTS_H
#define PHCORRELATORCONSTANTS_H

// c++ utilities
#include <string>



namespace PHEnergyCorrelator {

  // ==========================================================================
  //! PHEnergyCorrelator Constants
  // ==========================================================================
  namespace Const {

    // ------------------------------------------------------------------------
    //! Base for log axes
    // ------------------------------------------------------------------------
    inline double Base() {
      const double base = 10.;
      return base;
    }

    // ------------------------------------------------------------------------
    //! Tag for PtJet bins
    // ------------------------------------------------------------------------
    inline std::string PtTag() {
      const std::string pttag = "pt";
      return pttag;
    }

    // ------------------------------------------------------------------------
    //! Tag for CFJet bins
    // ------------------------------------------------------------------------
    inline std::string CFTag() {
      const std::string cftag = "cf";
      return cftag;
    }

    // ------------------------------------------------------------------------
    //! Tag for spin bins
    // ------------------------------------------------------------------------
    inline std::string SpinTag() {
      const std::string sptag = "sp";
      return sptag;
    }

  }   // end Const namespace
}  // end PHEnergyCorrelator namespace

#endif

// end =========================================================================
