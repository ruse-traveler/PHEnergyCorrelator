/// ============================================================================
/*! \file    PHCorrelatorBinning.h
 *  \authors Derek Anderson
 *  \date    09.27.2024
 *
 *  Class to define histogram binnings.
 */
/// ============================================================================

#ifndef PHCORRELATORBINNING_H
#define PHCORRELATORBINNING_H

// c++ utilities
#include <string>
#include <vector>
// analysis components
#include "PHCorrelatorTools.h"
#include "PHCorrelatorTypes.h"



namespace PHEnergyCorrelator {

  // ==========================================================================
  //! Binning definition
  // ==========================================================================
  /*! A small class to consolidate data
   *  for defining histogram bins.
   */ 
  class Binning {

    private:

      // data members
      double              m_start;
      double              m_stop;
      std::size_t         m_num;
      std::vector<double> m_bins;

    public:

      // ----------------------------------------------------------------------
      //! Uniform bin getters
      //-----------------------------------------------------------------------
      double      GetStart() const {return m_start;}
      double      GetStop()  const {return m_stop;}
      std::size_t GetNum()   const {return m_num;}

      // ----------------------------------------------------------------------
      //! Variable bin getter
      // ----------------------------------------------------------------------
      std::vector<double> GetBins() const {return m_bins;}

      // ----------------------------------------------------------------------
      //! default ctor/dtor
      // ----------------------------------------------------------------------
      Binning()  {};
      ~Binning() {};

      // ----------------------------------------------------------------------
      //! ctor accepting uniform parameters
      // ----------------------------------------------------------------------
      Binning(
        const std::size_t num,
        const double start,
        const double stop,
        const Type::Axis axis = Type::Norm
      ) {

        m_num   = num;
        m_start = start;
        m_stop  = stop;
        m_bins  = Tools::GetBinEdges(m_num, m_start, m_stop, axis);

      }  // end ctor(uint32_t, double, double, Type::Axis)

      // ----------------------------------------------------------------------
      //! ctor accepting non-uniform parameters
      // ----------------------------------------------------------------------
      Binning(const std::vector<double> edges) {

        m_bins  = edges;
        m_num   = edges.size() - 1;
        m_start = edges.front();
        m_stop  = edges.back();

      }  // end ctor(std::vector<double>)

  };  // end Binning

}  // end PHEnergyCorrelator namespace

#endif

// end ========================================================================
