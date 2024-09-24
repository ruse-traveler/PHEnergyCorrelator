/// ============================================================================
/*! \file    PHCorrelatorTools.h
 *  \authors Derek Anderson
 *  \date    09.23.2024
 *
 *  Namespace to collect useful methods for ENC calculations.
 */
/// ============================================================================

#ifndef PHCORRELATORTOOLS_H
#define PHCORRELATORTOOLS_H

// c++ utilities
#include <algorithm>
#include <cmath>
#include <vector>
// root libraries
#include <TLorentzVector.h>
#include <TVector3.h>



namespace PHEnergyCorrelator {

  // ==========================================================================
  //! PHEnergyCorrelator Tools
  // ==========================================================================
  namespace Tools {

    // ------------------------------------------------------------------------
    //! Divide a range into a certain number of bins
    // ------------------------------------------------------------------------
    /* TODO add flag for log vs. not */
    std::vector<double> GetBinEdges(
      const uint32_t num,
      const double start,
      const double stop
    ) {

      std::vector<double> bins;

      /* TODO calculations go here */
      return bins;

    }  // end 'GetBinEdges(uint32_t, double, double)'



    // ------------------------------------------------------------------------
    //! Get constituent 4-vector in lab frame from 3-momenta in jet frame
    // ------------------------------------------------------------------------
    /* TODO might be good to collect (z, jt, eta, phi) into a struct */
    TLorentzVector GetCstPxPyPz(
      const TVector3 pjet,
      const float z,
      const float jt,
      const float eta,
      const float phi
    ) {

      // get fractional momenta wrt jet
      TVector3 pcst = z * pjet;

      // calculate momentum components in lab frame
      const float ptot = hypot(jt, pcst.Mag());
      const float px   = ptot * std::cosh(eta) * std::cos(phi);
      const float py   = ptot * std::cosh(eta) * std::sin(phi);
      const float pz   = ptot * std::sinh(eta);

      // return 4-vector
      return TLorentzVector(px, py, pz, ptot);

    }  // end 'GetCstpxPyPz(TVector3, float, float, float, float)'

  }  // end Tools namespace
}  // end PHEnergyCorrelator namespace

#endif

// end ========================================================================
