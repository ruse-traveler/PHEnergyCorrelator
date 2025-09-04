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
#include <TRandom3.h>

namespace PHEnergyCorrelator {

  // Random generator for null spin vector
  TRandom3 *rand = new TRandom3();

  // ==========================================================================
  //! PHEnergyCorrelator Constants
  // ==========================================================================
  namespace Const {

    // ------------------------------------------------------------------------
    //! Total number of bins per spin case
    // -----------------------------------------------------------------------
    /*! This is the number of "bins" you'll get for each possible spin:
     *    [0] integrated pt and charge;
     *    [1] binned pt, integrated charge;
     *    [2] integrated pt, binned charge;
     *    [3] binned pt and charge
     */
    inline std::size_t NBinsPerSpin() {
      const std::size_t num = 4;
      return num;
    }

    // ------------------------------------------------------------------------
    //! Blue spin start point
    // ------------------------------------------------------------------------
    /*! Starting point in index vector for blue spins
     */
    inline std::size_t BlueSpinStart() {
      const std::size_t start = NBinsPerSpin();
      return start;
    }

    // ------------------------------------------------------------------------
    //! Yellow spin start point
    // ------------------------------------------------------------------------
    /*! Starting point in index vector for blue spins
     */
    inline std::size_t YellSpinStart() {
      const std::size_t start = 2 * NBinsPerSpin();
      return start;
    }

    // ------------------------------------------------------------------------
    //! Default value for index arguments
    // ------------------------------------------------------------------------
    inline std::size_t IndexDefault() {
      const std::size_t def = 0;
      return def;
    }

    // ------------------------------------------------------------------------
    //! Default value for int arguments
    // ------------------------------------------------------------------------
    inline int IntDefault() {
      const int def = std::numeric_limits<int>::max();
      return def;
    }

    // ------------------------------------------------------------------------
    //! Default value for double arguments
    // ------------------------------------------------------------------------
    inline double DoubleDefault() {
      const double def = std::numeric_limits<double>::max();
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
    //! Tag for integrated "bins"
    // ------------------------------------------------------------------------
    inline std::string IntTag() {
      const std::string inttag = "INT";
      return inttag;
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
    //! Tag for charge bins
    // ------------------------------------------------------------------------
    inline std::string ChrgTag() {
      const std::string chtag = "ch";
      return chtag;
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
    /*! We actually return a random normalized spin in the x-y plane to avoid
     *  problems with the vector calculations.  Since there is in principle no
     *  spin dependence this should be a valid thing to do.
     */
    inline TVector3 SpinNull() {
      const TVector3 null(rand->Uniform(), rand->Uniform(), 0.0);
      return null.Unit();
    }

  }   // end Const namespace
}  // end PHEnergyCorrelator namespace

#endif

// end =========================================================================
