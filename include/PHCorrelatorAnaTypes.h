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
      double charge;
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
        const double qarg = Const::DoubleDefault(),
        const int    aarg = Const::IntDefault()
      ) {
        cf      = carg;
        pt      = parg;
        eta     = harg;
        phi     = farg;
        charge  = qarg;
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
      std::size_t chrg;
      std::size_t spin;

      //! default ctor/dtor
      HistIndex()  {};
      ~HistIndex() {};

      //! ctor accepting arguments
      HistIndex(
        const std::size_t ipt,
        const std::size_t icf = Const::IndexDefault(),
        const std::size_t ich = Const::IndexDefault(),
        const std::size_t isp = Const::IndexDefault()
      ) {
        pt   = ipt;
        cf   = icf;
        chrg = ich;
        spin = isp;
      }  // end ctor(std::size_t x 4)

    };  // end HistIndex



    // ------------------------------------------------------------------------
    //! Histogram content
    // ------------------------------------------------------------------------
    struct HistContent {

      // data members
      double weight;     //!< energy weight
      double rl;         //!< longest side length
      double rm;         //!< medium side length (for E3C)
      double rs;         //!< shortest side length (for E3C and greater)
      double xi;         //!< \f$R_{s}/R_{m}\f$ (for E3C)
      double theta;      //!< \f$asin\sqrt(1-\frac{(R_{l}-R_{m})^{2}}{R_{s}^{2}})\f$ (for E3C)
      double phiCollB;   //!< blue collins angle: phiSpin - phiHadron
      double phiCollY;   //!< yellow collins angle: phiSpin - phiHadron
      double phiBoerB;   //!< blue boer-mulders angle: phiSpin - (2*phiHadron)
      double phiBoerY;   //!< yellow boer-mulders angle: phiSpin - (2*phiHadron)
      double thSpinRCB;  //!< blue DiFF angle: phiSpin - RC
      double thSpinRCY;  //!< yellow DiFF angle: phiSpin - RC
      double spinB;      //!< blue spin y-component
      double spinY;      //!< yellow spin y-component
      int    pattern;    //!< spin pattern

      //! default ctor/dtor
      HistContent()  {};
      ~HistContent() {};

      //! ctor accepting only 2-point arguments
      HistContent(
        const double w,
        const double l,
        const double db = Const::DoubleDefault(),
        const double dy = Const::DoubleDefault(),
        const double sb = Const::DoubleDefault(),
        const double sy = Const::DoubleDefault(),
        const int    p = Const::IntDefault()
      ) {
        weight    = w;
        rl        = l;
        thSpinRCB = db;
        thSpinRCY = dy;
        spinB     = sb;
        spinY     = sy;
        pattern   = p;
      }  // end ctor(double x 7, int)

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
        const double cb = Const::DoubleDefault(),
        const double cy = Const::DoubleDefault(),
        const double bb = Const::DoubleDefault(),
        const double by = Const::DoubleDefault(),
        const double db = Const::DoubleDefault(),
        const double dy = Const::DoubleDefault(),
        const double sb = Const::DoubleDefault(),
        const double sy = Const::DoubleDefault(),
        const int    p = Const::IntDefault()
      ) {
        weight    = w;
        rl        = l;
        rm        = m;
        rs        = s;
        xi        = x;
        theta     = t;
        phiCollB  = cb;
        phiCollY  = cy;
        phiBoerB  = bb;
        phiBoerY  = by;
        thSpinRCB = db;
        thSpinRCY = dy;
        spinB     = sb;
        spinY     = sy;
        pattern   = p;
      }  // end ctor(double x 13, int x 1)

    };  // end HistContent

  }  // end Type namespace
}  // end PHEnergyCorrelator namespace

#endif

// end =========================================================================
