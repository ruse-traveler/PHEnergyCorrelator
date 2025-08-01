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
      //! Get angle between jet-beam and spin-beam planes
      // ----------------------------------------------------------------------
      /*! Calculate angle between jet-beam and the spin-beam planes.
       *
       *    \param unitJetBeam3  normalized jet-beam plane normal vector
       *    \param unitSpinBeam3 normalized spin-beam plane normal vector
       *    \param vecSpin3      unnormalized spin vector
       */
      double GetPhiSpin(
        const TVector3& unitJetBeam3,
        const TVector3& unitSpinBeam3,
        const TVector3& vecSpin3
      ) {

        // get angle
        double phiSpin = atan2(
          unitJetBeam3.Cross(unitSpinBeam3).Mag(),
          unitJetBeam3.Dot(unitSpinBeam3)
        );

        // determine correct sign of angle and return
	//   - we'll define the zero of phiSpin as the jet-beam plane
	//     and the sense of rotation
	//   - the result will then be in [0, 2pi] rather than
	//     [0, pi]
	if (unitJetBeam3.Dot(vecSpin3) < 0.0) {
          phiSpin = TMath::TwoPi() - phiSpin; 
        }
        return phiSpin;

      }  // end 'GetPhiSpin(TVector3& x 3)'

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
