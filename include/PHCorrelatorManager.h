/// ============================================================================
/*! \file    PHCorrelatorManager.h
 *  \authors Derek Anderson
 *  \date    09.24.2024
 *
 *  Class to manage histograms to be filled during
 *  ENC calculations.
 */
/// ============================================================================

#ifndef PHCORRELATORMANAGER_H
#define PHCORRELATORMANAGER_H

// c++ utilities
#include <map>
#include <string>
#include <vector>
// root libraries
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
// analysis components
#include "PHCorrelatorBins.h"



namespace PHEnergyCorrelator {

  // ==========================================================================
  //! ENC Histogram Manager
  // ==========================================================================
  class Manager {

    private:

      // data members
      std::map<std::string, TH1D*> m_map_1d;
      std::map<std::string, TH2D*> m_map_2d;
      std::map<std::string, TH3D*> m_map_3d;

      /* TODO
       *   - add flags
       *   - CreateTags()
       *   - GenerateEECHists()
       *   - GenerateE3CHists()
       *   - GenerateLECHists()
       */

    public:

      // default ctor/dtor
      Manager()  {};
      ~Manager() {};

      /* TODO
       *   - SetOptions()
       *   - GetHistTag()
       *   - GetHist1D(string) 
       *   - GetHist2D(string) 
       *   - GetHist3D(string) 
       *   - GenerateHists()
       *   - SaveHists()
       *   - FillHists(Tag, Weight, RL, RM, RS)
       */

  };  // end PHEnergyCorrelator::Manager

}  // end PHEnergyCorrelator namespace

#endif

// end ========================================================================
