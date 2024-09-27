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



namespace PHEnergyCorrelator {

  // ==========================================================================
  //! PHEnergyCorrelator Constants
  // ==========================================================================
  namespace Const {

    // ------------------------------------------------------------------------
    //! Base for log axes
    // ------------------------------------------------------------------------
    inline double const& Base() {
       static double base = 10.;
       return base;
    }

  }   // end Const namespace
}  // end PHEnergyCorrelator namespace

#endif

// end =========================================================================
