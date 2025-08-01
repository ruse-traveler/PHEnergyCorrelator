/// ============================================================================
/*! \file    PHCorrelatorAngler.h
 *  \authors Derek Anderson
 *  \date    08.01.2024
 *
 *  Class to calculate angles during
 *  ENC calculations.
 */
/// ============================================================================

#ifndef PHCORRELATORANGLER_H
#define PHCORRELATORANGLER_H

// c++ utilities
#include <cmath>
// root libraries
#include <TLorentzVector.h>
#include <TMath.h>
#include <TVector3.h>



namespace PHEnergyCorrelator {

  // ==========================================================================
  //! ENC Angle Calculator 
  // ==========================================================================
  class Angler {

    private:

      bool   m_do_wrap;  /*!< turn angle wrapping on/off */
      double m_wrap;     /*!< value to wrap angles by */

    public:

      // ----------------------------------------------------------------------
      //! Getter
      // ----------------------------------------------------------------------
      bool   GetDoWrap() const {return m_do_wrap;}
      double GetWrap()   const {return m_wrap;}

      // ----------------------------------------------------------------------
      //! Setter
      // ----------------------------------------------------------------------
      void SetDoWrap(const bool do_wrap) {m_do_wrap = do_wrap;}
      void SetWrap(const double wrap)    {m_wrap    = wrap;}

      // ----------------------------------------------------------------------
      //! Get angle between two planes
      // ----------------------------------------------------------------------
      /*! Calculate angle between two planes (e.g. jet-beam and spin-beam
       *  planes)
       *
       *    \param unitPlaneA3 normalized normal vector for 1st plane
       *    \param unitPlaneB3 normalized normal vector for 2nd plane
       *    \param vecSensor3  unnormalized vector to determine sign of angle
       */
      double GetTwoPlaneAngle(
        const TVector3& unitPlaneA3,
        const TVector3& unitPlaneB3,
        const TVector3& vecSensor3
      ) {

        // get angle
        double phi = atan2(
          unitPlaneA3.Cross(unitPlaneB3).Mag(),
          unitPlaneA3.Dot(unitPlaneB3)
        );

        // determine correct sign of angle and return
        //   - by definition, atan2 only returns values
        //     between [0, pi]
	//   - we'll then define the zero of phi as plane A so
	//     that if vecSensor3 is in the opposite hemisphere
	//     of plane A's normal, we take the supplementary
	//     angle
	//   - that way the result will then be in [0, 2pi] rather
	//     than [0, pi]
        if (unitPlaneA3.Dot(vecSensor3) < 0.0) {
          phi = TMath::TwoPi() - phi;
        }
        return phi;

      }  // end 'GetTwoPlaneAngle(TVector3& x 3)'

      // ----------------------------------------------------------------------
      //! default ctor
      // ----------------------------------------------------------------------
      Angler()  {

        m_do_wrap = true;
        m_wrap    = TMath::Pi();

      }  // end default ctor

      // ----------------------------------------------------------------------
      //! default dtor
      // ----------------------------------------------------------------------
      ~Angler() {};

      // ----------------------------------------------------------------------
      //! ctor accepting arguments
      // ----------------------------------------------------------------------
      Angler(const bool do_wrap, const double wrap = TMath::Pi()) {

        m_do_wrap = do_wrap;
        m_wrap    = wrap;

      }  // end ctor(bool, double)

  };  // end PHEnergyCorrelator::Angler

}  // end PHEnergyCorrelator namespace

#endif

// end ========================================================================
