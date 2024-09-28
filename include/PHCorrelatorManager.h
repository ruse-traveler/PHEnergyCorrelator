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
#include <cassert>
#include <map>
#include <string>
#include <vector>
// root libraries
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
// analysis components
#include "PHCorrelatorBins.h"
#include "PHCorrelatorConstants.h"
#include "PHCorrelatorHistogram.h"
#include "PHCorrelatorTools.h"
#include "PHCorrelatorTypes.h"



namespace PHEnergyCorrelator {

  // ==========================================================================
  //! ENC Histogram Manager
  // ==========================================================================
  class Manager {

    private:

      /* TODO
       *   - GenerateGenHists() (?) 
       *   - GenerateE3CHists()
       *   - GenerateLECHists()
       */

      // data members (options)
      bool m_do_eec_hist;
      bool m_do_e3c_hist;
      bool m_do_lec_hist;
      bool m_do_pt_bins;
      bool m_do_cf_bins;
      bool m_do_sp_bins;

      // data members (no. of bins)
      std::size_t m_nbins_pt;
      std::size_t m_nbins_cf;
      std::size_t m_nbins_sp;

      // data members (index tags)
      std::vector<std::string> m_index_tags;

      // data members (histograms)
      std::map<std::string, TH1D*> m_hist_1d;
      std::map<std::string, TH2D*> m_hist_2d;
      std::map<std::string, TH3D*> m_hist_3d;

      // data members (bins)
      Bins m_bins;

      // ----------------------------------------------------------------------
      //! Make a tag from a histogram index 
      // ----------------------------------------------------------------------
      std::string MakeTag(const Type::HistIndex& index) {

        // by default, return empty string
        std::string tag = "";

        if (m_do_pt_bins) tag += Const::PtTag() + Tools::StringifyIndex(index.pt);
        if (m_do_cf_bins) tag += Const::CFTag() + Tools::StringifyIndex(index.cf);
        if (m_do_sp_bins) tag += Const::SpinTag() + Tools::StringifyIndex(index.spin);
        return tag;

      }  // end 'GetTag(Type::HistIndex&)'

      // ----------------------------------------------------------------------
      //! Create tags for bins
      // ----------------------------------------------------------------------
      void CreateBinTags() {

        // build list of indices
        std::vector<Type::HistIndex> indices;
        for (std::size_t ipt = 0; ipt < m_nbins_pt; ++ipt) {
          for (std::size_t icf = 0; icf < m_nbins_cf; ++icf) {
            for (std::size_t isp = 0; isp < m_nbins_sp; ++isp) {
              indices.push_back( Type::HistIndex(ipt, icf, isp) );
            }
          }
        }

        // create tags and return
        m_index_tags.clear();
        for (std::size_t index = 0; index < indices.size(); ++index) {
          m_index_tags.push_back( MakeTag(indices[index]) );
        }
        return;

      }  // end 'CreateBinTags()'

      // ----------------------------------------------------------------------
      //! Generate 2-point histograms
      // ----------------------------------------------------------------------
      void GenerateEECHists() {

        // 1d histogram definitions
        std::vector<Histogram> def_1d;
        def_1d.push_back( Histogram("EECStat",     "", {"R_{L}"}, {m_bins.Get("side")})    );
        def_1d.push_back( Histogram("EECWidth",    "", {"R_{L}"}, {m_bins.Get("side")})    );
        def_1d.push_back( Histogram("LogEECStat",  "", {"R_{L}"}, {m_bins.Get("logside")}) );
        def_1d.push_back( Histogram("LogEECWidth", "", {"R_{L}"}, {m_bins.Get("logside")}) );

        // create histograms
        for (std::size_t ihist = 0; ihist < def_1d.size(); ++ihist) {
          for (std::size_t index = 0; index < m_index_tags.size(); ++index) {

            // grab definition, adjust name
            Histogram hist = def_1d[ihist];
            hist.PrependToName("h");
            hist.AppendToName("_" + m_index_tags[index]);

            // create histogram, set errors
            m_hist_1d[hist.GetName()] = hist.MakeTH1();
            m_hist_1d[hist.GetName()] -> Sumw2();

          }
        }
        return;

      }  // end 'GenerateEECHists()'

    public:

      /* TODO
       *   - SetOptions()
       *   - GenerateHists()
       *   - SaveHists()
       *   - FillHists(Tag, Weight, RL, RM, RS)
       */

      // ----------------------------------------------------------------------
      //! Bin on jet pt
      // ----------------------------------------------------------------------
      void DoPtJetBins(const std::size_t nbins) {

        m_nbins_pt   = nbins;
        m_do_pt_bins = true;
        return;

      }  // end 'DoPtJetBins(std::size_t)'

      // ----------------------------------------------------------------------
      //! Bin on jet CF
      // ----------------------------------------------------------------------
      void DoCFJetBins(const std::size_t nbins) {

        m_nbins_pt   = nbins;
        m_do_cf_bins = true;
        return;

      }  // end 'DoCFJetBins(std::size_t)'

      // ----------------------------------------------------------------------
      //! Bin on spin
      // ----------------------------------------------------------------------
      void DoSpinBins(const std::size_t nbins) {

        m_nbins_sp   = nbins;
        m_do_sp_bins = true;
        return;

      }  // end 'DoSpinBins(std::size_t)'

      // ----------------------------------------------------------------------
      //! Get a 1D histogram
      // ----------------------------------------------------------------------
      TH1D* GetHist1D(const std::string& tag) {

        // throw error if binning doesn't exist
        if (m_hist_1d.count(tag) == 0) {
          assert(m_hist_1d.count(tag) >= 1);
        }

        // otherwise return hist
        return m_hist_1d[tag];

      }  // end 'GetHist1D(std::string&)'

      // ----------------------------------------------------------------------
      //! Get a 2D histogram
      // ----------------------------------------------------------------------
      TH2D* GetHist2D(const std::string& tag) {

        // throw error if binning doesn't exist
        if (m_hist_2d.count(tag) == 0) {
          assert(m_hist_2d.count(tag) >= 1);
        }

        // otherwise return hist
        return m_hist_2d[tag];

      }  // end 'GetHist2D(std::string&)'

      // ----------------------------------------------------------------------
      //! Get a 3D histogram
      // ----------------------------------------------------------------------
      TH3D* GetHist3D(const std::string& tag) {

        // throw error if binning doesn't exist
        if (m_hist_3d.count(tag) == 0) {
          assert(m_hist_3d.count(tag) >= 1);
        }

        // otherwise return hist
        return m_hist_3d[tag];

      }  // end 'GetHist3D(std::string&)'

      // ----------------------------------------------------------------------
      //! Get a histogram tag from a histogram index
      // ----------------------------------------------------------------------
      std::string GetTag(const Type::HistIndex& index) {

        return MakeTag(index);

      }  // end 'GetTag(Type::HistIndex&)'

      // ----------------------------------------------------------------------
      //! default ctor/dtor
      // ----------------------------------------------------------------------
      Manager()  {};
      ~Manager() {};

  };  // end PHEnergyCorrelator::Manager

}  // end PHEnergyCorrelator namespace

#endif

// end ========================================================================
