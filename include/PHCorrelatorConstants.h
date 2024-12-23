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
// root libraries
#include <TVector3.h>



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

    // ------------------------------------------------------------------------
    //! Spin up in lab coordinates
    // ------------------------------------------------------------------------
    inline TVector3 SpinUp() {
      const TVector3 up(0.0, 1.0, 0.0);
      return up;
    }

    // ------------------------------------------------------------------------
    //! Spin down in lab coordinates
    // ------------------------------------------------------------------------
    inline TVector3 SpinDown() {
      const TVector3 down(0.0, -1.0, 0.0);
      return down;
    }

    // ------------------------------------------------------------------------
    //! Null spin in lab coordinates
    // ------------------------------------------------------------------------
    inline TVector3 SpinNull() {
      const TVector3 null(0.0, 0.0, 0.0);
      return null;
    }

  }   // end Const namespace
}  // end PHEnergyCorrelator namespace

#endif

// end =========================================================================
