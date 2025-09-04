/// ============================================================================
/*! \file    PHCorrelatorHistManager.h
 *  \authors Derek Anderson
 *  \date    09.24.2024
 *
 *  Class to manage histograms to be filled during
 *  ENC calculations.
 */
/// ============================================================================

#ifndef PHCORRELATORHISTMANAGER_H
#define PHCORRELATORHISTMANAGER_H

// c++ utilities
#include <algorithm>
#include <cassert>
#include <map>
#include <stdint.h>
#include <string>
#include <vector>
// root libraries
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TString.h>
// analysis components
#include "PHCorrelatorAnaTools.h"
#include "PHCorrelatorAnaTypes.h"
#include "PHCorrelatorBins.h"
#include "PHCorrelatorConstants.h"
#include "PHCorrelatorHistogram.h"

//! Turns on/off width calculation.
#define DO_WIDTH_CALC 0
 

namespace PHEnergyCorrelator {

  // -------------------------------------------------------------------------
  //! Iterators for iterating over histogram maps
  // -------------------------------------------------------------------------
  typedef std::map<unsigned int, TH1D*>::iterator it_th1;
  typedef std::map<unsigned int, TH2D*>::iterator it_th2;
  typedef std::map<unsigned int, TH3D*>::iterator it_th3;



  // ==========================================================================
  //! ENC Histogram Manager
  // ==========================================================================
  class HistManager {

    public:

      // ----------------------------------------------------------------------
      //! Indices of spin configurations
      // ----------------------------------------------------------------------
      enum Spin {
        Int  = 0,  /*!< integrated over spin */
        BU   = 1,  /*!< blue up (int. over yellow) */
        BD   = 2,  /*!< blue down (int. over yellow) */
        YU   = 3,  /*!< yellow up (int. over blue) */
        YD   = 4,  /*!< yellow down (int. over blue) */
        BUYU = 5,  /*!< blue up, yellow up */
        BUYD = 6,  /*!< blue up, yellow up */
        BDYU = 7,  /*!< blue down, yellow up */
        BDYD = 8   /*!< blue down, yellow down */
      };

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
      bool m_do_ch_bins;
      bool m_do_sp_bins;

      // data members (no. of bins)
      std::size_t m_nbins_pt;
      std::size_t m_nbins_cf;
      std::size_t m_nbins_ch;
      std::size_t m_nbins_sp;

      // data members (tags)
      std::string              m_hist_tag;
      std::string              m_hist_pref;
      std::vector<std::string> m_index_tags;

      // data members (histograms)
      std::map<unsigned int, TH1D*> m_hist_1d;
      std::map<unsigned int, TH2D*> m_hist_2d;
      std::map<unsigned int, TH3D*> m_hist_3d;

      // data members (bins)
      Bins m_bins;

      // ----------------------------------------------------------------------
      //! Convert an index to a string
      // ----------------------------------------------------------------------
      std::string StringifyIndex(const std::size_t index) const {

        // create TString, add index
        TString tstr;
        tstr += index;

        // create std::string and return
        const std::string sstr(tstr.Data());
        return sstr;

      }  // end 'StringifyIndex(std::size_t)'

      // ----------------------------------------------------------------------
      //! Translate pt index into a label
      // ----------------------------------------------------------------------
      std::string GetPtLabel(const std::size_t ipt) const {

        if (ipt < m_nbins_pt) {
          return StringifyIndex(ipt);
        } else {
          return Const::IntTag();
        }

      }  // end 'GetPtLabel(std::size_t, std::size_t)'

      // ----------------------------------------------------------------------
      //! Translate charge index into a label
      // ----------------------------------------------------------------------
      std::string GetChrgLabel(const std::size_t ich) const {

        if (ich < m_nbins_ch) {
          return StringifyIndex(ich);
        } else {
          return Const::IntTag();
        }

      }  // end 'GetChrgLabel(std::size_t, std::size_t)'

      // ----------------------------------------------------------------------
      //! Translate spin index into a label
      // ----------------------------------------------------------------------
      /*! Spin "bins" correspond to different combinations of
       *  blue/yellow polarizations, indexed by the enum `Spin`.
       *  They are
       */ 
      std::string GetSpinLabel(const std::size_t isp) const {

        std::string label = "";
        switch (isp) {
          case Int:
            label = Const::IntTag();
            break;
          case BU:
            label = "BU";
            break;
          case BD:
            label = "BD";
            break;
          case YU:
            label = "YU";
            break;
          case YD:
            label = "YD";
            break;
          case BUYU:
            label = "BUYU";
            break;
          case BUYD:
            label = "BUYD";
            break;
          case BDYU:
            label = "BDYU";
            break;
          case BDYD:
            label = "BDYD";
            break;
          default:
            label = "";
            break;
        }
        return label;

      }  // end 'GetSpinLabel(std::size_t)'

      // ----------------------------------------------------------------------
      //! Make a tag from a histogram index 
      // ----------------------------------------------------------------------
      std::string MakeIndexTag(const Type::HistIndex& index) const {

        // by default, return empty string
        std::string tag = "";

        // add appropriate indices
        if (m_do_pt_bins) tag += Const::PtTag() + GetPtLabel(index.pt);
        if (m_do_cf_bins) tag += Const::CFTag() + StringifyIndex(index.cf);
        if (m_do_ch_bins) tag += Const::ChrgTag() + GetChrgLabel(index.chrg);
        if (m_do_sp_bins) tag += Const::SpinTag() + GetSpinLabel(index.spin);
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
      //! Jenkins' one-at-a-time hash for strings
      // ----------------------------------------------------------------------
      /*! From this stack overflow post:
       *    - https://stackoverflow.com/questions/114085/fast-string-hashing-algorithm-with-low-collision-rates-with-32-bit-integer
       *  For the Jenkins hashes, see here:
       *    - https://en.wikipedia.org/wiki/Jenkins_hash_function
       *    - https://www.burtleburtle.net/bob/hash/doobs.html
       */
      unsigned int HashString(const char* string) const {

        unsigned int hash = 0;
        for (; *string; ++string) {
          hash += *string;
          hash += (hash << 10);
          hash ^= (hash >> 6);
        }
        hash += (hash << 3);
        hash ^= (hash >> 11);
        hash += (hash << 15);
        return hash;

      }  // end 'HashString(char*)'

      // ----------------------------------------------------------------------
      //! Make a histogram name and hash it
      // ----------------------------------------------------------------------
      unsigned int MakeHashedName(const std::string& base, const std::string& tag) const {

        const std::string name = MakeHistName(base, tag);
        return HashString( name.data() );

      }

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
                m_hist_1d[ HashString(hist.GetName().data()) ] = hist.MakeTH1();
                break;
              case 2:
                m_hist_2d[ HashString(hist.GetName().data()) ] = hist.MakeTH2();
                break;
              case 3:
                m_hist_3d[ HashString(hist.GetName().data()) ] = hist.MakeTH3();
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

        // set n bins for pt & charge
        //   - n.b. once more, these will have 1 additional "bin",
        //     for integration
        const std::size_t nbins_pt_use = m_nbins_pt + 1;
        const std::size_t nbins_ch_use = m_nbins_ch + 1;

        // build list of indices
        std::vector<Type::HistIndex> indices;
        for (std::size_t ipt = 0; ipt < nbins_pt_use; ++ipt) {
          for (std::size_t icf = 0; icf < m_nbins_cf; ++icf) {
            for (std::size_t ich = 0; ich < nbins_ch_use; ++ich) {
              for (std::size_t isp = 0; isp < m_nbins_sp; ++isp) {
                indices.push_back( Type::HistIndex(ipt, icf, ich, isp) );
              }
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

        // 1d spin axis titles
        const std::string diffB_title("#theta_{s_{B}} - #theta_{R_{c}}");
        const std::string diffY_title("#theta_{s_{Y}} - #theta_{R_{c}}");

        // 1d histogram definitions
        std::vector<Histogram> def_1d;
        def_1d.push_back(
          Histogram("EECStat", "", "R_{L}", m_bins.Get("side"))
        );
        def_1d.push_back(
          Histogram("DiFFBlueSpinRC", "", diffB_title, m_bins.Get("angle"))
        );
        def_1d.push_back(
          Histogram("DiFFYellSpinRC", "", diffY_title, m_bins.Get("angle"))
        );

        // vectors of binnings for 2d histograms
        std::vector<Binning> angleXside_bins;
        angleXside_bins.push_back(m_bins.Get("side"));
        angleXside_bins.push_back(m_bins.Get("angle"));

        // vectors of axis titles for 2d histograms
        std::vector<std::string> diffXsideB_titles;
        std::vector<std::string> diffXsideY_titles;
        diffXsideB_titles.push_back("R_{L}");
        diffXsideB_titles.push_back(diffB_title);
        diffXsideY_titles.push_back("R_{L}");
        diffXsideY_titles.push_back(diffY_title);

        // 2D histogram definitions
        std::vector<Histogram> def_2d;
        def_2d.push_back(
          Histogram("DiFFBlueSpinRCVsRStat", "", diffXsideB_titles, angleXside_bins)
        );
        def_2d.push_back(
          Histogram("DiFFYellSpinRCVsRStat", "", diffXsideY_titles, angleXside_bins)
        );

        // create histograms
        MakeHistograms(def_1d, 1);
        MakeHistograms(def_2d, 2);
        return;

      }  // end 'GenerateEECHists()'

#if DO_WIDTH_CALC
      // ----------------------------------------------------------------------
      //! Set variances of relevant 2-point histograms
      // ----------------------------------------------------------------------
      /*! N.B. this function is currently unused. If/when we decide to include
       *  histograms with widths reported rather than statisical errors, this
       *  will come into play.
       */
      void SetEECVariances() {

        // names of EEC hists to set variances of
        std::vector<std::string> to_set;

        for (std::size_t index = 0; index < m_index_tags.size(); ++index) {
          for (std::size_t iset = 0; iset < to_set.size(); ++iset) {
            Histogram::SetHist1DErrToVar(
              m_hist_1d.at( MakeHashedName(to_set[iset], m_index_tags[index]) )
            );
          }
        }
        return;

      }  // end 'SetEECVariances()'
#endif

    public:

      // ----------------------------------------------------------------------
      //! Getters
      // ----------------------------------------------------------------------
      std::string GetHistTag()      const {return m_hist_tag;}
      std::size_t GetNPtJetBins()   const {return m_nbins_pt;}
      std::size_t GetNCFJetBins()   const {return m_nbins_cf;}
      std::size_t GetNChargeBins()  const {return m_nbins_ch;}
      std::size_t GetNSpinBins()    const {return m_nbins_sp;}
      std::size_t GetNIndexTags()   const {return m_index_tags.size();}
      std::size_t GetNHist1D()      const {return m_hist_1d.size();}
      std::size_t GetNHist2D()      const {return m_hist_2d.size();}
      std::size_t GetNHist3D()      const {return m_hist_3d.size();}
      std::size_t GetNHists()       const {return GetNHist1D() + GetNHist2D() + GetNHist3D();}
      bool        GetDoPtJetBins()  const {return m_do_pt_bins;}
      bool        GetDoCFJetBins()  const {return m_do_cf_bins;}
      bool        GetDoChargeBins() const {return m_do_ch_bins;}
      bool        GetDoSpinBins()   const {return m_do_sp_bins;}
      bool        GetDoEECHists()   const {return m_do_eec_hist;}
      bool        GetDoE3CHists()   const {return m_do_e3c_hist;}
      bool        GetDoLECHists()   const {return m_do_lec_hist;}

      // ----------------------------------------------------------------------
      //! Setters
      // ----------------------------------------------------------------------
      void SetHistTag(const std::string& tag) {m_hist_tag    = tag;}
      void SetDoEECHists(const bool dohists)  {m_do_eec_hist = dohists;}
      void SetDoE3CHists(const bool dohists)  {m_do_e3c_hist = dohists;}
      void SetDoLECHists(const bool dohists)  {m_do_lec_hist = dohists;}

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
      //! Bin on jet charge
      // ----------------------------------------------------------------------
      void DoChargeBins(const std::size_t nbins) {

        m_nbins_ch   = nbins;
        m_do_ch_bins = true;
        return;

      }  // end 'DoChargeBins(std::size_t)'

      // ----------------------------------------------------------------------
      //! Bin on spin
      // ----------------------------------------------------------------------
      void DoSpinBins(const bool spin) {

        m_do_sp_bins = spin;
        return;

      }  // end 'DoSpinBins()'

      // ----------------------------------------------------------------------
      //! Generate histograms
      // ----------------------------------------------------------------------
      void GenerateHists() {

        // 1st make sure there'll be at least 1 index
        //   - n.b. as always, not that pt and chrg will add one
        //     "bin" for integration
        m_nbins_pt = m_do_pt_bins ? m_nbins_pt : (std::size_t) 0;
        m_nbins_cf = m_do_cf_bins ? m_nbins_cf : (std::size_t) 1;
        m_nbins_ch = m_do_ch_bins ? m_nbins_ch : (std::size_t) 0;
        m_nbins_sp = m_do_sp_bins ? m_nbins_sp : (std::size_t) 1;

        // then create tags for each bin and histogrma prefixes
        GenerateIndexTags();
        GenerateHistPrefix();

        // turn on errors ahead of hist generation
        TH1::SetDefaultSumw2(true);
        TH2::SetDefaultSumw2(true);
        TH3::SetDefaultSumw2(true);

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
        m_hist_1d[ MakeHashedName("EECStat", tag) ] -> Fill(content.rl, content.weight);
        m_hist_1d[ MakeHashedName("DiFFBlueSpinRC", tag) ] -> Fill(content.thSpinRCB);
        m_hist_1d[ MakeHashedName("DiFFYellSpinRC", tag) ] -> Fill(content.thSpinRCY);

        // fill 2d histograms
        m_hist_2d[ MakeHashedName("DiFFBlueSpinRCVsRStat", tag) ] -> Fill(
          content.rl, content.thSpinRCB, content.weight
        );
        m_hist_2d[ MakeHashedName("DiFFYellSpinRCVsRStat", tag) ] -> Fill(
          content.rl, content.thSpinRCY, content.weight
        );
        return;

      }  // end 'FillEECHists(Type::HistIndex&, Type::HistContent&)'

      // ----------------------------------------------------------------------
      //! Save histograms to a file
      // ----------------------------------------------------------------------
      void SaveHists(TFile* file) {

#if DO_WIDTH_CALC
        // set variances on relevant histograms
        if (m_do_eec_hist) SetEECVariances();
#endif

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
        const unsigned int key = HashString(tag.data());
        if (m_hist_1d.count(key) == 0) {
          assert(m_hist_1d.count(key) >= 1);
        }

        // otherwise return hist
        return m_hist_1d[key];

      }  // end 'GetHist1D(std::string&)'

      // ----------------------------------------------------------------------
      //! Get a 2D histogram
      // ----------------------------------------------------------------------
      TH2D* GetHist2D(const std::string& tag) {

        // throw error if binning doesn't exist
        const unsigned int key = HashString(tag.data());
        if (m_hist_2d.count(key) == 0) {
          assert(m_hist_2d.count(key) >= 1);
        }

        // otherwise return hist
        return m_hist_2d[key];

      }  // end 'GetHist2D(std::string&)'

      // ----------------------------------------------------------------------
      //! Get a 3D histogram
      // ----------------------------------------------------------------------
      TH3D* GetHist3D(const std::string& tag) {

        // throw error if binning doesn't exist
        const unsigned int key = HashString(tag.data());
        if (m_hist_3d.count(key) == 0) {
          assert(m_hist_3d.count(key) >= 1);
        }

        // otherwise return hist
        return m_hist_3d[key];

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
      HistManager()  {

        m_do_eec_hist = false;
        m_do_e3c_hist = false;
        m_do_lec_hist = false;
        m_do_pt_bins  = false;
        m_do_cf_bins  = false;
        m_do_ch_bins  = false;
        m_do_sp_bins  = false;
        m_nbins_pt    = 0;  // n.b. there will always be 1 additional integrated "bin"
        m_nbins_cf    = 1;
        m_nbins_ch    = 0;  // n.b. there will always be 1 additional integrated "bin"
        m_nbins_sp    = 9;
        m_hist_tag    = "";
        m_hist_pref   = "";

      }  // end default ctor

      // ----------------------------------------------------------------------
      //! default dtor
      // ----------------------------------------------------------------------
      ~HistManager() {};

      // ----------------------------------------------------------------------
      //! ctor accepting arguments
      // ----------------------------------------------------------------------
      HistManager(const bool do_eec, const bool do_e3c = false, const bool do_lec = false) {

        m_do_eec_hist = do_eec;
        m_do_e3c_hist = do_e3c;
        m_do_lec_hist = do_lec;
        m_do_pt_bins  = false;
        m_do_cf_bins  = false;
        m_do_ch_bins  = false;
        m_do_sp_bins  = false;
        m_nbins_pt    = 0;  // n.b. there will always be 1 additional integrated "bin"
        m_nbins_cf    = 1;
        m_nbins_ch    = 0;  // n.b. there will always be 1 additional integrated "bin"
        m_nbins_sp    = 9;
        m_hist_tag    = "";
        m_hist_pref   = "";

      }  // end 'HistManager(bool, bool, bool)'

  };  // end PHEnergyCorrelator::HistManager

}  // end PHEnergyCorrelator namespace

#endif

// end ========================================================================
