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
#include <limits>
#include <string>
// root libraries
#include <TVector3.h>



namespace PHEnergyCorrelator {

  // ==========================================================================
  //! PHEnergyCorrelator Constants
  // ==========================================================================
  namespace Const {

    // ------------------------------------------------------------------------
    // Default value for index arguments
    // ------------------------------------------------------------------------
    inline std::size_t IndexDefault() {
      const std::size_t def = 0;
      return def;
    }

    // ------------------------------------------------------------------------
    // Default value for int arguments
    // ------------------------------------------------------------------------
    inline int IntDefault() {
      const int def = std::numeric_limits<int>::max();
      return def;
    }

    // ------------------------------------------------------------------------
    // Default value for double arguments
    // ------------------------------------------------------------------------
    inline double DoubleDefault() {
      const int def = std::numeric_limits<double>::max();
      return def;
    }

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
    //! Blue beam direction in lab coordinates
    // ------------------------------------------------------------------------
    inline TVector3 BlueBeam() {
      const TVector3 blue(0.0, 0.0, 1.0);
      return blue;
    }

    // ------------------------------------------------------------------------
    //! Yellow beam direction in lab coordinates
    // ------------------------------------------------------------------------
    inline TVector3 YellowBeam() {
      const TVector3 yellow(0.0, 0.0, -1.0);
      return yellow;
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
