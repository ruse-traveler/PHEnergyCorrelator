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
      double m_min;      /*!< min value to wrap to */
      double m_max;      /*!< max value to wrap to */

    public:

      // ----------------------------------------------------------------------
      //! Getter
      // ----------------------------------------------------------------------
      bool   GetDoWrap() const {return m_do_wrap;}
      double GetWrap()   const {return m_wrap;}
      double GetMin()    const {return m_min;}
      double GetMax()    const {return m_max;}

      // ----------------------------------------------------------------------
      //! Setter
      // ----------------------------------------------------------------------
      void SetDoWrap(const bool do_wrap) {m_do_wrap = do_wrap;}
      void SetWrap(const double wrap)    {m_wrap    = wrap;}
      void SetMin(const double min)      {m_min     = min;}
      void SetMax(const double max)      {m_max     = max;}

      // ----------------------------------------------------------------------
      //! Wrap an angle
      // ----------------------------------------------------------------------
      void Wrap(double& angle) const {

        if (angle <  m_min) angle += m_wrap;
        if (angle >= m_max) angle -= m_wrap;
        return;

      }  // end 'Wrap(double&)'

      // ----------------------------------------------------------------------
      //! Get angle between two planes
      // ----------------------------------------------------------------------
      /*! Calculate angle between two planes (e.g. jet-beam and spin-beam
       *  planes)
       *
       *    \param unit_plane_a3 normalized normal vector for 1st plane
       *    \param unit_plane_b3 normalized normal vector for 2nd plane
       *    \param vec_sensor3  unnormalized vector to determine sign of angle
       */
      double GetTwoPlaneAngle(
        const TVector3& unit_plane_a3,
        const TVector3& unit_plane_b3,
        const TVector3& vec_sensor3
      ) const {

        // get angle
        double phi = atan2(
          unit_plane_a3.Cross(unit_plane_b3).Mag(),
          unit_plane_a3.Dot(unit_plane_b3)
        );

        // determine correct sign of angle and return
        //   - by definition, atan2 only returns values
        //     between [0, pi]
	//   - we'll then define the zero of phi as plane A so
	//     that if vec_sensor3 is in the opposite hemisphere
	//     of plane A's normal, we take the supplementary
	//     angle
	//   - that way the result will then be in [0, 2pi] rather
	//     than [0, pi]
        if (unit_plane_a3.Dot(vec_sensor3) < 0.0) {
          phi = TMath::TwoPi() - phi;
        }
        return phi;

      }  // end 'GetTwoPlaneAngle(TVector3& x 3)'

      // ----------------------------------------------------------------------
      //! Get Collins angle
      // ----------------------------------------------------------------------
      /*! Calculates collins angle given a phi_spin and a phi_hadron.
       *  Angle is wrapped accordingly based on m_do_wrap and m_wrap.
       *
       *    \params phi_spin angle between spin-beam plane and jet-beam plane
       *    \params phi_had  angle between hadron-jet plane and jet-beam plane
       */
      double GetCollinsAngle(const double phi_spin, const double phi_had) const {

        double coll = phi_spin - phi_had;
        if (m_do_wrap) {
          Wrap(coll);
        }
        return coll;

      }  // end 'GetCollinsAngle(double, double)'

      // ----------------------------------------------------------------------
      //! Get Boer-Mulders angle
      // ----------------------------------------------------------------------
      /*! Calculates boer-mulders angle given a phi_spin and a phi_hadron.
       *  Angle is wrapped accordingly based on m_do_wrap and m_wrap.
       *
       *    \params phi_spin angle between spin-beam plane and jet-beam plane
       *    \params phi_had  angle between hadron-jet plane and jet-beam plane
       */
      double GetBoerMuldersAngle(const double phi_spin, const double phi_had) const {

        // double phi_hadron
        double had2 = 2.0 * phi_had;
        if (m_do_wrap) {
          Wrap(had2);
        }

        // get bm angle
        double boer = phi_spin - had2;
        if (m_do_wrap) {
          Wrap(boer);
        }
        return boer;

      }  // end 'GetBoerMuldersAngle(double, double)'

      // ----------------------------------------------------------------------
      //! Get DiFF spin-beam angle
      // ----------------------------------------------------------------------
      /*! Calculates theta_s, angle between spin and the scattering plane, the
       *  plane defined by the beam and PC (dihadron center-of-mass).  Angle is
       *  wrapped accordingly based on m_do_wrap and m_wrap.  For more details,
       *  see arXiv:hep-ph/0409174.
       *
       *  \param unit_beam     normalized beam direction (i.e. PB or PA)
       *  \param vec_spin      spin direction of beam (i.e. SB or SA)
       *  \param vec_dihad_com dihadron center-of-mass (i.e. PC)
       */
      double GetDiFFSpinBeamAngle(
        const TVector3& unit_beam,
        const TVector3& vec_spin,
        const TVector3& vec_dihad_com
      ) const {

        // calculate cosine of theta_s
	double costheta = (
          unit_beam.Cross(vec_dihad_com) * (
            1.0 / (
              unit_beam.Cross(vec_dihad_com).Mag()
            )
          )
        ).Dot(
          unit_beam.Cross(vec_spin) * (
            1.0 / unit_beam.Cross(vec_spin).Mag()
          )
        );

        // calculate sine of theta_s
	double sintheta = (
          vec_dihad_com.Cross(vec_spin)
        ).Dot(unit_beam) * (
          1.0 / (
            (unit_beam.Cross(vec_dihad_com).Mag()) *
            (unit_beam.Cross(vec_spin).Mag())
          )
        );

        // get thetas (thetaSB, thetaSA)
        //   - by definition, acos only returns values
        //     between [0, pi]
        //   - so the sign of sintheta will determine if
        //     angle is in [-pi, 0] or [0, pi]
        double thetas = (sintheta > 0.0) ? acos(costheta) : -acos(costheta);
        if (m_do_wrap) {
          Wrap(thetas);
        }
        return thetas;

      }  // end 'GetDiFFSpinBeamAngle(TVector3& x 3)'

      // ----------------------------------------------------------------------
      //! default ctor
      // ----------------------------------------------------------------------
      Angler()  {

        m_do_wrap = true;
        m_wrap    = TMath::TwoPi();
        m_min     = 0.0;
        m_max     = TMath::TwoPi();

      }  // end default ctor

      // ----------------------------------------------------------------------
      //! default dtor
      // ----------------------------------------------------------------------
      ~Angler() {};

      // ----------------------------------------------------------------------
      //! ctor accepting arguments
      // ----------------------------------------------------------------------
      Angler(
        const bool do_wrap,
        const double wrap = TMath::TwoPi(),
        const double min = 0.0,
        const double max = TMath::TwoPi()
      ) {

        m_do_wrap = do_wrap;
        m_wrap    = wrap;
        m_min     = min;
        m_max     = max;

      }  // end ctor(bool, double x 3)

  };  // end PHEnergyCorrelator::Angler

}  // end PHEnergyCorrelator namespace

#endif

// end ========================================================================
