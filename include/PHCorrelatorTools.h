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
#include <cassert>
#include <cmath>
#include <vector>
// root libraries
#include <TLorentzVector.h>
#include <TString.h>
#include <TVector3.h>
// analysis components
#include "PHCorrelatorConstants.h"
#include "PHCorrelatorTypes.h"



namespace PHEnergyCorrelator {

  // ==========================================================================
  //! PHEnergyCorrelator Tools
  // ==========================================================================
  namespace Tools {

    // ------------------------------------------------------------------------
    //! Wrapper function for std::pow
    // ------------------------------------------------------------------------
    double Exponentiate(const double arg) {

      return std::pow(Const::Base(), arg);

    }  // end 'Exponentiate(double&)'



    // ------------------------------------------------------------------------
    //! Wrapper function for log
    // ------------------------------------------------------------------------
    double Log(const double arg) {

      return std::log10(arg) / std::log10(Const::Base());

    }  // end 'Log(double&)'



    // ------------------------------------------------------------------------
    //! Convert an index to a string
    // ------------------------------------------------------------------------
    std::string StringifyIndex(const std::size_t index) {

      TString str((Ssiz_t) index);
      return str.Data();

    }  // end 'StringifyIndex(std::size_t)'



    // ------------------------------------------------------------------------
    //! Divide a range into a certain number of bins
    // ------------------------------------------------------------------------
    std::vector<double> GetBinEdges(
      const std::size_t num,
      const double start,
      const double stop,
      const Type::Axis axis = Type::Log
    ) {

      // throw error if start/stop are out of order
      // or if num is zero
      if (num <= 0)     assert(num > 0);
      if (start > stop) assert(start <= stop);

      // set start/stop, calculate bin steps
      const double start_use = (axis == Type::Log) ? Log(start) : start;
      const double stop_use  = (axis == Type::Log) ? Log(stop)  : stop;
      const double step      = (stop_use - start_use) / num;
 
      // instantiate vector to hold bins
      std::vector<double> bins;

      // and fill vector
      double edge = start_use;
      for (std::size_t inum = 0; inum < num; ++inum) {
        bins.push_back( edge );
        edge += step;
      }
      bins.push_back( edge );

      // if need be, transform back from log
      if (axis == Type::Log) {
        std::transform(
          bins.begin(),
          bins.end(),
          bins.begin(),
          Exponentiate
        ); 
      }
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
