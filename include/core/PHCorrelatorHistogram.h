/// ============================================================================
/*! \file    PHCorrelatorHistogram.h
 *  \authors Derek Anderson
 *  \date    09.27.2024
 *
 *  Class to define histograms.
 */
/// ============================================================================

#ifndef PHCORRELATORHISTOGRAM_H
#define PHCORRELATORHISTOGRAM_H

// c++ utilities
#include <cassert>
#include <string>
#include <vector>
// root libraries
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
// analysis components
#include "PHCorrelatorAnaTools.h"
#include "PHCorrelatorAnaTypes.h"
#include "PHCorrelatorBinning.h"



namespace PHEnergyCorrelator {

  // ==========================================================================
  //! Histogram definition
  // ==========================================================================
  /*! A small class to consolidate necessary
   *  data to define a TH1, TH2, or TH3. Used
   *  to produce a THND. 
   */ 
  class Histogram {

    private:

      // data members
      std::string m_name;
      std::string m_title;
      std::string m_title_x;
      std::string m_title_y;
      std::string m_title_z;
      Binning     m_bins_x;
      Binning     m_bins_y;
      Binning     m_bins_z;

      // ----------------------------------------------------------------------
      //! Make histogram title
      // ----------------------------------------------------------------------
      std::string MakeTitle() const {

        std::string title = m_title;
        title.append( ";" + m_title_x );
        title.append( ";" + m_title_y );
        title.append( ";" + m_title_z );
        return title;

      }  // end 'MakeTitle()'

    public:

      // ----------------------------------------------------------------------
      //! Getters
      // ----------------------------------------------------------------------
      std::string GetName()   const {return m_name;}
      std::string GetTitle()  const {return m_title;}
      std::string GetTitleX() const {return m_title_x;}
      std::string GetTitleY() const {return m_title_y;}
      std::string GetTitleZ() const {return m_title_z;}
      Binning     GetBinsX()  const {return m_bins_x;}
      Binning     GetBinsY()  const {return m_bins_y;}
      Binning     GetBinsZ()  const {return m_bins_z;}

      // ----------------------------------------------------------------------
      //! Set and modify histogram title/name
      // ----------------------------------------------------------------------
      void SetHistTitle(const std::string& title)   {m_title = title;}
      void SetHistName(const std::string& name)     {m_name = name;}
      void PrependToName(const std::string& prefix) {m_name = prefix + m_name;}
      void AppendToName(const std::string& suffix)  {m_name = m_name + suffix;}

      // ----------------------------------------------------------------------
      //! Set specific axis title
      // ----------------------------------------------------------------------
      void SetAxisTitle(const std::string& title, const std::size_t axis) {

        switch (axis) {
          case 0:
            m_title_x = title;
            break;
          case 1:
            m_title_y = title;
            break;
          case 2:
            m_title_z = title;
            break;
          default:
            assert(axis < 3);
            break;
        }
        return;

      }  // end 'SetAxisTitle(std::string&, std::size_t)'

      // ----------------------------------------------------------------------
      //! Set specific binning
      // ----------------------------------------------------------------------
      void SetAxisBins(const Binning& bins, const std::size_t axis) {

        switch (axis) {
          case 0:
            m_bins_x = bins;
            break;
          case 1:
            m_bins_y = bins;
            break;
          case 2:
            m_bins_z = bins;
            break;
          default:
            assert(axis < 3);
            break;
        }
        return;

      }  // end 'SetAxisBins(Binning&, std::size_t)'

      // ----------------------------------------------------------------------
      //! Set axis titles via list
      // ----------------------------------------------------------------------
      void SetAxisTitles(const std::vector<std::string>& titles) {

        if (titles.size() >= 1) m_title_x = titles[0];
        if (titles.size() >= 2) m_title_y = titles[1];
        if (titles.size() >= 3) m_title_z = titles[2];
        return;

      }  // end 'SetAxisTitles(std::vector<std::string>&)'

      // ----------------------------------------------------------------------
      //! Set binning via list
      // ----------------------------------------------------------------------
      void SetAxisBins(const std::vector<Binning>& bins) {

        if (bins.size() >= 1) m_bins_x = bins[0];
        if (bins.size() >= 2) m_bins_y = bins[1];
        if (bins.size() >= 3) m_bins_z = bins[2];
        return;

      }  // end 'SetAxisBins(std::vector<Binning>&)'

      // ----------------------------------------------------------------------
      //! Generate a TH1D
      // ----------------------------------------------------------------------
      TH1D* MakeTH1() const {

        // make hist + axis titles
        const std::string title = MakeTitle();

        TH1D* hist = new TH1D(
          m_name.data(),
          title.data(),
          m_bins_x.GetNum(),
          m_bins_x.GetBins().data()
        );
        return hist;

      }  // end 'MakeTH1()'

      // ----------------------------------------------------------------------
      //! Generate a TH2D
      // ----------------------------------------------------------------------
      TH2D* MakeTH2() const {

        // make hist + axis titles
        const std::string title = MakeTitle();

        TH2D* hist = new TH2D(
          m_name.data(),
          title.data(),
          m_bins_x.GetNum(),
          m_bins_x.GetBins().data(),
          m_bins_y.GetNum(),
          m_bins_y.GetBins().data()
        );
        return hist;

      }  // end 'MakeTH2()'

      // ----------------------------------------------------------------------
      //! Generate a TH3D
      // ----------------------------------------------------------------------
      TH3D* MakeTH3() const {

        // make hist + axis titles
        const std::string title = MakeTitle();

        TH3D* hist = new TH3D(
          m_name.data(),
          title.data(),
          m_bins_x.GetNum(),
          m_bins_x.GetBins().data(),
          m_bins_y.GetNum(),
          m_bins_y.GetBins().data(),
          m_bins_z.GetNum(),
          m_bins_z.GetBins().data()
        );
        return hist;

      }  // end 'MakeTH3()'

      // -----------------------------------------------------------------------
      //! Helper method to set 1d histogram errors to be variances
      // -----------------------------------------------------------------------
      static void SetHist1DErrToVar(TH1D* hist) {

        for (int ibin = 1; ibin <= hist -> GetNbinsX(); ++ibin) {
          const double var = Tools::GetVariance(
            hist -> GetBinError(ibin),
            hist -> GetBinContent(ibin)
          );
          hist -> SetBinError(ibin, var);
        }  // end bin loop
        return;

      }  // end 'SetHist1DErrTOVar(TH1D*)'

      // ----------------------------------------------------------------------
      //! default ctor/dtor
      // ----------------------------------------------------------------------
      Histogram()   {};
       ~Histogram() {};

      // ----------------------------------------------------------------------
      //! ctor accepting arguments (1 binning, 1 axis title)
      // ----------------------------------------------------------------------
      Histogram(
        const std::string& hist_name,
        const std::string& hist_title,
        const std::string& axis_title,
        const Binning& axis_bins
      ) {

        SetHistName(hist_name);
        SetHistTitle(hist_title);
        SetAxisTitle(axis_title, 0);
        SetAxisBins(axis_bins, 0);

      }  // end ctor(std::string&, std::string&, std::string&, Binning&)'

      // ----------------------------------------------------------------------
      //! ctor accepting arguments (vector of axis titles/bins)
      // ----------------------------------------------------------------------
      Histogram(
        const std::string& hist_name,
        const std::string& hist_title,
        const std::vector<std::string>& axis_titles,
        const std::vector<Binning>& axis_bins
      ) {

        SetHistName(hist_name);
        SetHistTitle(hist_title);
        SetAxisTitles(axis_titles);
        SetAxisBins(axis_bins);

      }  // end ctor(std::string&, std::string&, std::vector<std::string>&, std::vector<Binning>&)

  };  // end Histogram
}  // end PHEnergyCorrelator namespace

#endif

// end ========================================================================
