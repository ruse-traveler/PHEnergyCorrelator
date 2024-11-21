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
#include <algorithm>
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

  // -------------------------------------------------------------------------
  //! Iterators for iterating over histogram maps
  // -------------------------------------------------------------------------
  typedef std::map<std::string, TH1D*>::iterator it_th1;
  typedef std::map<std::string, TH2D*>::iterator it_th2;
  typedef std::map<std::string, TH3D*>::iterator it_th3;



  // ==========================================================================
  //! ENC Histogram Manager
  // ==========================================================================
  class Manager {

    private:

      /* TODO
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

      // data members (tags)
      std::string              m_hist_tag;
      std::string              m_hist_pref;
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
      std::string MakeIndexTag(const Type::HistIndex& index) const {

        // by default, return empty string
        std::string tag = "";

        // add appropriate indices
        if (m_do_pt_bins) tag += Const::PtTag() + Tools::StringifyIndex(index.pt);
        if (m_do_cf_bins) tag += Const::CFTag() + Tools::StringifyIndex(index.cf);
        if (m_do_sp_bins) tag += Const::SpinTag() + Tools::StringifyIndex(index.spin);
        return tag;

      }  // end 'MakeIndexTag(Type::HistIndex&)'

      // ----------------------------------------------------------------------
      //! Make a histogram suffix from a tag
      // ----------------------------------------------------------------------
      std::string MakeHistSuffix(const std::string& tag) const {

        return "_" + tag;

      }  // end 'MakeHistSuffix(std::string&)'

      // ----------------------------------------------------------------------
      //! Make a histogram name from a base and a tag
      // ----------------------------------------------------------------------
      std::string MakeHistName(const std::string& base, const std::string& tag) const {

        return m_hist_pref + base + MakeHistSuffix(tag);

      }  // end 'MakeHistName(std::string&, std::string&)'

      // ----------------------------------------------------------------------
      //! Make histograms out of a list of definitions
      // ----------------------------------------------------------------------
      void MakeHistograms(const std::vector<Histogram>& defs, const int dim) {

        for (std::size_t ihist = 0; ihist < defs.size(); ++ihist) {
          for (std::size_t index = 0; index < m_index_tags.size(); ++index) {

            // grab definition, adjust name
            Histogram hist = defs[ihist];
            hist.PrependToName( m_hist_pref );
            hist.AppendToName( MakeHistSuffix(m_index_tags[index]) );

            // create histogram, set errors
            switch (dim) {
              case 1:
                m_hist_1d[hist.GetName()] = hist.MakeTH1();
                m_hist_1d[hist.GetName()] -> Sumw2();
                break;
              case 2:
                m_hist_2d[hist.GetName()] = hist.MakeTH2();
                m_hist_2d[hist.GetName()] -> Sumw2();
                break;
              case 3:
                m_hist_3d[hist.GetName()] = hist.MakeTH3();
                m_hist_3d[hist.GetName()] -> Sumw2();
                break;
              default:
                assert((dim >= 1) && (dim <= 3));
                break;
            }
          }
        }
        return;

      }  // end 'MakeHistograms(std::vector<Histogram>&, uint32_t)'

      // ----------------------------------------------------------------------
      //! Generate tags for bins
      // ----------------------------------------------------------------------
      void GenerateIndexTags() {

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
          m_index_tags.push_back( MakeIndexTag(indices[index]) );
        }
        return;

      }  // end 'GenerateIndexTags()'

      // ----------------------------------------------------------------------
      //! Generate prefix for histograms
      // ----------------------------------------------------------------------
      void GenerateHistPrefix() {

        m_hist_pref = "h" + m_hist_tag;
        return;

      }  // end 'GenerateHistPrefix()'

      // ----------------------------------------------------------------------
      //! Generate 2-point histograms
      // ----------------------------------------------------------------------
      void GenerateEECHists() {

        // 1d histogram definitions
        std::vector<Histogram> def_1d;
        def_1d.push_back(Histogram("EECStat", "", "R_{L}", m_bins.Get("side")));
        def_1d.push_back(Histogram("EECWidth", "", "R_{L}", m_bins.Get("side")));
        def_1d.push_back(Histogram("LogEECStat", "", "log R_{L}", m_bins.Get("logside")));
        def_1d.push_back(Histogram("LogEECWidth", "", "log R_{L}", m_bins.Get("logside")));
        def_1d.push_back(Histogram("SpinPhiStat", "", "#varphi", m_bins.Get("angle")));

        // vectors of binnings for 2d histograms
        std::vector<Binning> spinside_bins;
        std::vector<Binning> spinlogside_bins;
        spinside_bins.push_back(m_bins.Get("side"));
        spinside_bins.push_back(m_bins.Get("angle"));
        spinlogside_bins.push_back(m_bins.Get("logside"));
        spinlogside_bins.push_back(m_bins.Get("angle"));

        // vectors of axis titles for 2d histograms
        std::vector<std::string> spinside_titles;
        std::vector<std::string> spinlogside_titles;
        spinside_titles.push_back("R_{L}");
        spinside_titles.push_back("#varphi");
        spinlogside_titles.push_back("log R_{L}");
        spinlogside_titles.push_back("#varphi");

        // 2D histogram definitions
        std::vector<Histogram> def_2d;
        def_2d.push_back(
          Histogram("EECPhiVsRStat", "", spinside_titles, spinside_bins)
        );
        def_2d.push_back(
          Histogram("EECPhiVsLogRStat", "", spinlogside_titles, spinlogside_bins)
        );

        // create histograms
        MakeHistograms(def_1d, 1);
        MakeHistograms(def_2d, 2);
        return;

      }  // end 'GenerateEECHists()'

      // ----------------------------------------------------------------------
      //! Set variances of relevant 2-point histograms
      // ----------------------------------------------------------------------
      void SetEECVariances() {

        // names of EEC hists to set variances of
        std::vector<std::string> to_set;
        to_set.push_back( "EECWidth" );
        to_set.push_back( "LogEECWidth" );

        for (std::size_t index = 0; index < m_index_tags.size(); ++index) {
          for (std::size_t iset = 0; iset < to_set.size(); ++iset) {
            Histogram::SetHist1DErrToVar(
              m_hist_1d[ MakeHistName(to_set[iset], m_index_tags[index]) ]
            );
          }
        }
        return;

      }  // end 'SetEECVariances()'

    public:

      // ----------------------------------------------------------------------
      //! Getters
      // ----------------------------------------------------------------------
      std::string GetHistTag()    const {return m_hist_tag;}
      std::size_t GetNPtJetBins() const {return m_nbins_pt;}
      std::size_t GetNCFJetBins() const {return m_nbins_cf;}
      std::size_t GetNSpinBins()  const {return m_nbins_sp;}
      std::size_t GetNIndexTags() const {return m_index_tags.size();}
      std::size_t GetNHist1D()    const {return m_hist_1d.size();}
      std::size_t GetNHist2D()    const {return m_hist_2d.size();}
      std::size_t GetNHist3D()    const {return m_hist_3d.size();}
      std::size_t GetNHists()     const {return GetNHist1D() + GetNHist2D() + GetNHist3D();}

      // ----------------------------------------------------------------------
      //! Set hist tag
      // ----------------------------------------------------------------------
      void SetHistTag(const std::string& tag) {m_hist_tag = tag;}

      // ----------------------------------------------------------------------
      //! Set histogram options
      // ----------------------------------------------------------------------
      void DoEECHists(const bool dohists) {m_do_eec_hist = dohists;}
      void DoE3CHists(const bool dohists) {m_do_e3c_hist = dohists;}
      void DoLECHists(const bool dohists) {m_do_lec_hist = dohists;}

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

        m_nbins_cf   = nbins;
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
      //! Generate histograms
      // ----------------------------------------------------------------------
      void GenerateHists() {

        // 1st make sure there'll be at least 1 index
        m_nbins_pt = m_do_pt_bins ? m_nbins_pt : (std::size_t) 1;
        m_nbins_cf = m_do_cf_bins ? m_nbins_cf : (std::size_t) 1;
        m_nbins_sp = m_do_sp_bins ? m_nbins_sp : (std::size_t) 1;

        // then create tags for each bin and histogrma prefixes
        GenerateIndexTags();
        GenerateHistPrefix();

        // finally generate appropriate histograms
        //   - TODO add others when ready
        if (m_do_eec_hist) GenerateEECHists();
        return;

      }  // end 'GenerateHists()'

      // ----------------------------------------------------------------------
      //! Fill EEC histograms
      // ----------------------------------------------------------------------
      void FillEECHists(const Type::HistIndex& index, const Type::HistContent& content) {

        // grab hist tag from index
        const std::string tag = MakeIndexTag(index);

        // fill 1d histograms
        m_hist_1d[MakeHistName("EECStat", tag)] -> Fill(content.rl, content.weight);
        m_hist_1d[MakeHistName("EECWidth", tag)] -> Fill(content.rl, content.weight);
        m_hist_1d[MakeHistName("LogEECStat", tag)] -> Fill(Tools::Log(content.rl), content.weight);
        m_hist_1d[MakeHistName("LogEECWidth", tag)] -> Fill(Tools::Log(content.rl), content.weight);
        m_hist_1d[MakeHistName("SpinPhiStat", tag)] -> Fill(content.phi);

        // fill 2d histograms
        m_hist_2d[MakeHistName("EECPhiVsRStat", tag)] -> Fill(content.rl, content.phi, content.weight);
        m_hist_2d[MakeHistName("EECPhiVsLogRStat", tag)] -> Fill(Tools::Log(content.rl), content.phi, content.weight);
        return;

      }  // end 'FillEECHists(Type::HistIndex&, Type::HistContent&)'

      // ----------------------------------------------------------------------
      //! Save histograms to a file
      // ----------------------------------------------------------------------
      void SaveHists(TFile* file) {

        // set variances on relevant histograms
        //   - TODO add others when ready
        if (m_do_eec_hist) SetEECVariances();

        // throw error if cd failed
        const bool good_cd = file -> cd();
        if (!good_cd) {
          assert(good_cd);
        }

        // then save histograms
        for (it_th1 it1d = m_hist_1d.begin(); it1d != m_hist_1d.end(); ++it1d) {
          it1d -> second -> Write();
        }
        for (it_th2 it2d = m_hist_2d.begin(); it2d != m_hist_2d.end(); ++it2d) {
          it2d -> second -> Write();
        }
        for (it_th3 it3d = m_hist_3d.begin(); it3d != m_hist_3d.end(); ++it3d) {
          it3d -> second -> Write();
        }
        return;

      }  // end 'SaveHists(TFile*)'

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
      std::string GetIndexTag(const Type::HistIndex& index) {

        return MakeIndexTag(index);

      }  // end 'GetIndexTag(Type::HistIndex&)'

      // ----------------------------------------------------------------------
      //! default ctor
      // ----------------------------------------------------------------------
      Manager()  {

        m_do_eec_hist = false;
        m_do_e3c_hist = false;
        m_do_lec_hist = false;
        m_do_pt_bins  = false;
        m_do_cf_bins  = false;
        m_do_sp_bins  = false;
        m_nbins_pt    = 1;
        m_nbins_cf    = 1;
        m_nbins_sp    = 1;
        m_hist_tag    = "";
        m_hist_pref   = "";

      }  // end default ctor

      // ----------------------------------------------------------------------
      //! default dtor
      // ----------------------------------------------------------------------
      ~Manager() {};

      // ----------------------------------------------------------------------
      //! ctor accepting arguments
      // ----------------------------------------------------------------------
      Manager(const bool do_eec, const bool do_e3c = false, const bool do_lec = false) {

        m_do_eec_hist = do_eec;
        m_do_e3c_hist = do_e3c;
        m_do_lec_hist = do_lec;
        m_do_pt_bins  = false;
        m_do_cf_bins  = false;
        m_do_sp_bins  = false;
        m_nbins_pt    = 1;
        m_nbins_cf    = 1;
        m_nbins_sp    = 1;
        m_hist_tag    = "";
        m_hist_pref   = "";

      }  // end 'Manager(bool, bool, bool)'

  };  // end PHEnergyCorrelator::Manager

}  // end PHEnergyCorrelator namespace

#endif

// end ========================================================================
