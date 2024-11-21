/// ============================================================================
/*! \file    PHCorrelatorTypes.h
 *  \authors Derek Anderson
 *  \date    09.21.2024
 *
 *  Namespace to collect types used in PHENIX
 *  ENC analyses.
 */
/// ============================================================================

#ifndef PHCORRELATORTYPES_H
#define PHCORRELATORTYPES_H



namespace PHEnergyCorrelator {

  // ==========================================================================
  //! PHEnergyCorrelator Types
  // ==========================================================================
  namespace Type {

    // ------------------------------------------------------------------------
    //! Types of binning
    // ------------------------------------------------------------------------
    enum Axis {Log, Norm};

    // ------------------------------------------------------------------------
    //! Weight types
    // ------------------------------------------------------------------------
    enum Weight {E, Et, Pt};



    // ------------------------------------------------------------------------
    //! Jet information
    // ------------------------------------------------------------------------
    struct Jet {

      // data members
      double cf;
      double pt;
      double eta;
      double phi;
      double spin;

      //! default ctor/dtor
      Jet()  {};
      ~Jet() {};

      //! ctor accepting arguments
      Jet(
        const double carg,
        const double parg,
        const double harg,
        const double farg,
        const double sarg
      ) {
        cf   = carg;
        pt   = parg;
        eta  = harg;
        phi  = farg;
        spin = sarg;
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
        const std::size_t icf = 0,
        const std::size_t isp = 0
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
      double weight;
      double rl;
      double rm;
      double rs;
      double xi;
      double theta;
      double phi;

      //! default ctor/dtor
      HistContent()  {};
      ~HistContent() {};

      //! ctor accepting only 2-point arguments
      HistContent(
        const double w,
        const double l,
        const double f = 0.
      ) {
        weight = w;
        rl     = l;
        phi    = f;
      }  // end ctor(double x 3)

      //! ctor accepting only 3-point arguments
      HistContent(
        const double w,
        const double x,
        const double t,
        const double l,
        const double m = 0.,
        const double s = 0.
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
        const double f = 0.
      ) {
        weight = w;
        rl     = l;
        rm     = m;
        rs     = s;
        xi     = x;
        theta  = t;
        phi    = f;
      }  // end ctor(double x 7)

    };  // end HistContent

  }  // end Type namespace
}  // end PHEnergyCorrelator namespace

#endif

// end =========================================================================
