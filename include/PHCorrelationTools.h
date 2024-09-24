/// ============================================================================
/*! \file    PHCorrelationTools.h
 *  \authors Derek Anderson
 *  \date    09.23.2024
 *
 *  Namespace to collect useful methods.
 */
/// ============================================================================

#ifndef PHCORRELATIONTOOLS_H
#define PHCORRELATIONTOOLS_H

// c++ utilities
#include <cmath>
// root libraries
#include <TLorentzVector.h>
#include <TVector3.h>



namespace PHEnergyCorrelator {

  // ==========================================================================
  //! PHEnergyCorrelator Tools
  // ==========================================================================
  namespace Tools {

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
      const float ptot = std::hypot(jt, pcst.Mag());
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
