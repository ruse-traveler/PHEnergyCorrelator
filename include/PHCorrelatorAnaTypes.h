/// ============================================================================
/*! \file    PHCorrelatorAnaTypes.h
 *  \authors Derek Anderson
 *  \date    09.21.2024
 *
 *  Useful types related to analysis/histogramming.
 */
/// ============================================================================

#ifndef PHCORRELATORANATYPES_H
#define PHCORRELATORANATYPES_H

// analysis components
#include "PHCorrelatorConstants.h"



namespace PHEnergyCorrelator {
  namespace Type {

    // ------------------------------------------------------------------------
    //! Types of binning
    // ------------------------------------------------------------------------
    enum Axis {Log, Norm};

    // ------------------------------------------------------------------------
    //! Weight types
    // ------------------------------------------------------------------------
    enum Weight {
      E,   /*!< weight by energy */
      Et,  /*!< weight by transverse energy */
      Pt   /*!< weight by transverse momentum (lab frame) */
    };

    // ----------------------------------------------------------------------
    //! Possible spin patterns
    // ----------------------------------------------------------------------
    enum Pattern {
      PPBUYU = 0,  /*!< blue up, yellow up (pp) */
      PPBDYU = 1,  /*!< blue down, yellow up (pp) */
      PPBUYD = 2,  /*!< blue up, yellow down (pp) */
      PPBDYD = 3,  /*!< blue down, yellow down (pp) */
      PABU   = 4,  /*!< blue up (pAu) */
      PABD   = 5   /*!< blue down (pAu) */
    };



    // ------------------------------------------------------------------------
    //! Jet information
    // ------------------------------------------------------------------------
    struct Jet {

      // data members
      double cf;
      double pt;
      double eta;
      double phi;
      int    pattern;

      //! default ctor/dtor
      Jet()  {};
      ~Jet() {};

      //! ctor accepting arguments
      Jet(
        const double carg,
        const double parg,
        const double harg,
        const double farg,
        const int    aarg = Const::IntDefault()
      ) {
        cf      = carg;
        pt      = parg;
        eta     = harg;
        phi     = farg;
        pattern = aarg; 
      };

    };  // end Jet



    // ------------------------------------------------------------------------
    //! Constituent information
    // ------------------------------------------------------------------------
    struct Cst {

      // data members
      double z;
      double jt;
      double eta;
      double phi;
      double chrg;

      //! default ctor/dtor
      Cst()  {};
      ~Cst() {};

      //! ctor accpeting arguments
      Cst(
        const double zarg,
        const double jarg,
        const double harg,
        const double farg,
        const double carg
      ) {
        z    = zarg;
        jt   = jarg;
        eta  = harg;
        phi  = farg;
        chrg = carg;
      }  // end ctor(double x 5)

    };  // end Cst



    // ------------------------------------------------------------------------
    //! Histogram index
    // ------------------------------------------------------------------------
    struct HistIndex {

      // data members
      std::size_t pt;
      std::size_t cf;
      std::size_t spin;

      //! default ctor/dtor
      HistIndex()  {};
      ~HistIndex() {};

      //! ctor accepting arguments
      HistIndex(
        const std::size_t ipt,
        const std::size_t icf = Const::IndexDefault(),
        const std::size_t isp = Const::IndexDefault()
      ) {
        pt   = ipt;
        cf   = icf;
        spin = isp;
      }  // end ctor(std::size_t, std::size_t, std::size_t)

    };  // end HistIndex



    // ------------------------------------------------------------------------
    //! Histogram content
    // ------------------------------------------------------------------------
    struct HistContent {

      // data members
      double weight;   //!< energy weight
      double rl;       //!< longest side length
      double rm;       //!< medium side length (for E3C)
      double rs;       //!< shortest side length (for E3C and greater)
      double xi;       //!< \f$R_{s}/R_{m}\f$ (for E3C)
      double theta;    //!< \f$asin\sqrt(1-\frac{(R_{l}-R_{m})^{2}}{R_{s}^{2}})\f$ (for E3C)
      double phiHAvgB; //!< hadron-hadron angle wrt blue spin
      double phiHAvgY; //!< hadron-hadron angle wrt yellow spin
      double phiCollB; //!< collins (hadron-hadron-jet) angle wrt blue spin
      double phiCollY; //!< collins (hadron-hadron-jet) angle wrt yellow spin
      double spinB;    //!< blue spin y-component
      double spinY;    //!< yellow spin y-component
      int    pattern;  //!< spin pattern

      //! default ctor/dtor
      HistContent()  {};
      ~HistContent() {};

      //! ctor accepting only 2-point arguments
      HistContent(
        const double w,
        const double l,
        const double ab = Const::DoubleDefault(),
        const double ay = Const::DoubleDefault(),
        const double cb = Const::DoubleDefault(),
        const double cy = Const::DoubleDefault(),
        const double sb = Const::DoubleDefault(),
        const double sy = Const::DoubleDefault(),
        const int    p = Const::IntDefault()
      ) {
        weight   = w;
        rl       = l;
        phiHAvgB = ab;
        phiHAvgY = ay;
        phiCollB = cb;
        phiCollY = cy;
        spinB    = sb;
        spinY    = sy;
        pattern  = p;
      }  // end ctor(double x 8, int)

      //! ctor accepting only 3-point arguments
      HistContent(
        const double w,
        const double x,
        const double t,
        const double l,
        const double m = Const::DoubleDefault(),
        const double s = Const::DoubleDefault()
      ) {
        weight = w;
        xi     = x;
        theta  = t;
        rl     = l;
        rm     = m;
        rs     = s;
      }  // end ctor(double x 6)

      //! ctor accepting all arguments
      HistContent(
        const double w,
        const double l,
        const double m,
        const double s,
        const double x,
        const double t,
        const double ab = Const::DoubleDefault(),
        const double ay = Const::DoubleDefault(),
        const double cb = Const::DoubleDefault(),
        const double cy = Const::DoubleDefault(),
        const double sb = Const::DoubleDefault(),
        const double sy = Const::DoubleDefault(),
        const int    p = Const::IntDefault()
      ) {
        weight   = w;
        rl       = l;
        rm       = m;
        rs       = s;
        xi       = x;
        theta    = t;
        phiHAvgB = ab;
        phiHAvgY = ay;
        phiCollB = cb;
        phiCollY = cy;
        spinB    = sb;
        spinY    = sy;
        pattern  = p;
      }  // end ctor(double x 11, int x 1)

    };  // end HistContent

  }  // end Type namespace
}  // end PHEnergyCorrelator namespace

#endif

// end =========================================================================
