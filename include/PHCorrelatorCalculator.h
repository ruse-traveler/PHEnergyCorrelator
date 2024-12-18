/// ============================================================================
/*! \file    PHCorrelatorCalculator.h
 *  \authors Derek Anderson, Alex Clarke
 *  \date    09.21.2024
 *
 *  Driver class to run n-point energy-energy correlator
 *  calculations on inputs.
 */
/// ============================================================================

#ifndef PHCORRELATORCALCULATOR_H
#define PHCORRELATORCALCULATOR_H

// c++ utilities
#include <algorithm>
#include <cmath>
#include <string>
#include <utility>
#include <vector>
// root libraries
#include <TLorentzVector.h>
#include <TMath.h>
#include <TVector3.h>
// analysis componenets
#include "PHCorrelatorManager.h"
#include "PHCorrelatorTypes.h"
#include "PHCorrelatorTools.h"

// TEST
#include <iostream>


namespace PHEnergyCorrelator {

  // ==========================================================================
  //! ENC Calculator
  // ==========================================================================
  class Calculator {

    private:

      // data members (calc options)
      double       m_weight_power;
      Type::Weight m_weight_type;

      // data members (bins)
      std::vector< std::pair<float, float> > m_ptjet_bins;
      std::vector< std::pair<float, float> > m_cfjet_bins;

      // data member (hist manager)
      Manager m_manager;

      // ---------------------------------------------------------------------=
      //! Get weight of a constituent
      // ----------------------------------------------------------------------
      double GetCstWeight(TLorentzVector& cst, TLorentzVector& jet) {

        // grab relevant cst & jet values
        double numer = 1.0;
        double denom = 1.0;
        switch (m_weight_type) {

          case Type::E:
            numer = cst.E();
            denom = jet.E();
            break;

          case Type::Et:
            numer = cst.Et();
            denom = jet.Et();
            break;

          case Type::Pt:
            numer = cst.Pt();
            denom = jet.Pt();
            break;

          default:
            numer = cst.Pt();
            denom = jet.Pt();
            break;

        }  // end switch

        // raise cst, jet values to specified value (defualt is 1.0)
        numer = pow(numer, m_weight_power);
        denom = pow(denom, m_weight_power);

        // calculate weight and exit
        const double weight = numer / denom;
        return weight;

      }  // end 'GetCstWeight(TLorentzVector&, TLorentzVector&)'

      // ----------------------------------------------------------------------
      //! Get hist index/indices
      // ----------------------------------------------------------------------
      /*! N.B. when doing spin binning, if an unexpected spin pattern is
       *  provided then only ONE index will be returned, the spin integrated
       *  index.
       */
      std::vector<Type::HistIndex> GetHistIndices(const Type::Jet& jet) {

        // for pt and cf, index will correspond to what bin
        // the jet falls in
        Type::HistIndex iptcf(0, 0, 0);

        // determine pt bin
        if (m_manager.GetDoPtJetBins()) {
          for (std::size_t ipt = 0; ipt < m_ptjet_bins.size(); ++ipt) {
            if ((jet.pt >= m_ptjet_bins[ipt].first) && (jet.pt < m_ptjet_bins[ipt].second)) {
              iptcf.pt = ipt;
            }
          }  // end pt bin loop
        }

        // determine cf bin
        if (m_manager.GetDoCFJetBins()) {
          for (std::size_t icf = 0; icf < m_cfjet_bins.size(); ++icf) {
            if ((jet.cf >= m_cfjet_bins[icf].first) && (jet.cf < m_cfjet_bins[icf].second)) {
              iptcf.cf = icf;
            }
          }  // end cf bin loop
        }

        // but for spin, we'll have 2 indices
        // for each spin pattern
        //   - pattern = 1 --> spin indices = {0, 3}
        //   - pattern = 2 --> spin indices = {1, 2}
        //   - pattern = 3 --> spin indices = {1, 3}
        //   - pattern = 4 --> spin indices = {0, 2}
        // plus the index for integrating over
        // spins (sp = 4)
        std::vector<Type::HistIndex> indices;

        // determine spin bins
        if (m_manager.GetDoSpinBins()) {
            switch (jet.pattern) {

              // blue up, yellow down
              case 1:
                indices.push_back( Type::HistIndex(iptcf.pt, iptcf.cf, 0) );
                indices.push_back( Type::HistIndex(iptcf.pt, iptcf.cf, 3) );
                break;

              // blue down, yellow up
              case 2:
                indices.push_back( Type::HistIndex(iptcf.pt, iptcf.cf, 1) );
                indices.push_back( Type::HistIndex(iptcf.pt, iptcf.cf, 2) );
                break;

              // blue down, yellow down
              case 3:
                indices.push_back( Type::HistIndex(iptcf.pt, iptcf.cf, 1) );
                indices.push_back( Type::HistIndex(iptcf.pt, iptcf.cf, 3) );
                break;

              // blue up, yellow up
              case 4:
                indices.push_back( Type::HistIndex(iptcf.pt, iptcf.cf, 0) );
                indices.push_back( Type::HistIndex(iptcf.pt, iptcf.cf, 2) );
                break;

              // by default, only add integrated
              default:
                break;
          }
        }
        indices.push_back( Type::HistIndex(iptcf.pt, iptcf.cf, 4) );
        return indices;

      }  // end 'GetHistIndices(Type::Jet&)'

    public:

      /* TODO
       *   - DoE3CCalc(jet, {cst, cst, cst})
       *   - DoLECCalc(jet, {lambda, cst})
       *   - FindLambda({cst...})
       */

      // ----------------------------------------------------------------------
      //! Getters
      // ----------------------------------------------------------------------
      Manager& GetManager() {return m_manager;}

      // ----------------------------------------------------------------------
      //! Setters
      // ----------------------------------------------------------------------
      void SetWeightPower(const double power)       {m_weight_power = power;}
      void SetWeightType(const Type::Weight weight) {m_weight_type  = weight;}
      void SetHistTag(const std::string& tag)       {m_manager.SetHistTag(tag);}

      // ----------------------------------------------------------------------
      //! Set jet pt bins
      // ----------------------------------------------------------------------
      void SetPtJetBins(const std::vector< std::pair<float, float> >& bins) {

        // copy bins to member
        m_ptjet_bins.resize( bins.size() );
        std::copy(bins.begin(), bins.end(), m_ptjet_bins.begin());

        // update hist manager and exit
        m_manager.DoPtJetBins( m_ptjet_bins.size() );
        return;

      }  // end 'SetPtJetBins(std::vector<std::pair<float, float>>&)'

      // ----------------------------------------------------------------------
      //! Set jet CF bins
      // ----------------------------------------------------------------------
      void SetCFJetBins(const std::vector< std::pair<float, float> >& bins) {

        // copy bins to member
        m_cfjet_bins.resize( bins.size() );
        std::copy(bins.begin(), bins.end(), m_cfjet_bins.begin());

        // update hist manager and exit
        m_manager.DoCFJetBins( m_cfjet_bins.size() );
        return;

      }  // end 'SetCFJetBins(std::vector<std::pair<float, float>>&)'

      // ----------------------------------------------------------------------
      //! Turn on/off spin binning
      // ----------------------------------------------------------------------
      void SetDoSpinBins(const bool spin) {

        // and set corresponding flag in manager
        m_manager.DoSpinBins(spin);
        return;

      }  // end 'SetDoSpinBins(bool)'

      // ----------------------------------------------------------------------
      //! Initialize calculator
      // ----------------------------------------------------------------------
      void Init(const bool do_eec, const bool do_e3c = false, const bool do_lec = false) {

        // turn on/off relevant histograms
        m_manager.SetDoEECHists(do_eec);
        m_manager.SetDoE3CHists(do_e3c);
        m_manager.SetDoLECHists(do_lec);

        // then generate necessary histograms
        m_manager.GenerateHists();
        return;

      } // end 'Init(bool, bool, bool)'

      // ----------------------------------------------------------------------
      //! Do EEC calculation
      // ----------------------------------------------------------------------
      /*! Note that the optional third argument, `evtweight`, is there
       *  to allow for weighting by ckin, spin, etc. By default, it's
       *  set to 1.
       */ 
      void CalcEEC(
        const Type::Jet& jet,
        const std::pair<Type::Cst, Type::Cst>& csts,
        const double evt_weight = 1.0
      ) {

        // get jet 4-momenta
        TLorentzVector jet_vec = Tools::GetJetLorentz(jet);

        // get cst 4-momenta
        std::pair<TLorentzVector, TLorentzVector> cst_vecs = std::make_pair(
          Tools::GetCstLorentz(jet_vec.Vect(), csts.first),
          Tools::GetCstLorentz(jet_vec.Vect(), csts.second)
        );

        // get vector distance b/n average of cst.s and spin direction
        TLorentzVector cst_avg = Tools::GetWeightedAvgLorentz(
          cst_vecs.first,
          cst_vecs.second
        );

        // now get weights
        std::pair<double, double> cst_weights = std::make_pair(
          GetCstWeight(cst_vecs.first, jet_vec),
          GetCstWeight(cst_vecs.second, jet_vec)
        );

        // calculate RL (dist b/n cst.s for EEC), EEC, and
        // angle b/n the cst average and spin
        const double dist    = Tools::GetCstDist(csts);
        const double weight  = cst_weights.first * cst_weights.second * evt_weight;
        const double dphiblu = remainder(cst_avg.Phi() - jet.phiblu, TMath::TwoPi());
        const double dphiyel = remainder(cst_avg.Phi() - jet.phiyel, TMath::TwoPi());

        // bundle results for histogram filling
        Type::HistContent content_int(weight, dist);
        Type::HistContent content_blu(weight, dist, dphiblu, jet.pattern);
        Type::HistContent content_yel(weight, dist, dphiyel, jet.pattern);

        // fill histograms and exit
        if (m_manager.GetDoEECHists()) {

          // grab hist indices; if doing spin binning, the order
          // of the vector will always be
          //   [0] = blue beam index
          //   [1] = yellow beam index
          //   [2] = integrated
          // otherwise vector will have 1 entry corresponding
          // to JUST the integrated case
          std::vector<Type::HistIndex> indices = GetHistIndices(jet);

          // fill relevant histograms
          if (m_manager.GetDoSpinBins() && (indices.size() > 1)) {
            m_manager.FillEECHists(indices[0], content_blu);
            m_manager.FillEECHists(indices[1], content_yel);
            m_manager.FillEECHists(indices[2], content_int);
          } else {
            m_manager.FillEECHists(indices[0], content_int);
          }
        }
        return;

      }  // end 'CalcEEC(Type::Jet&, std::pair<Type::Cst, Type::Cst>&, double)'

      // ----------------------------------------------------------------------
      //! End calculations
      // ----------------------------------------------------------------------
      void End(TFile* file) {

        // save histograms to file
        m_manager.SaveHists(file);
        return;

      }  // end 'End(TFile*)'

      // ----------------------------------------------------------------------
      //! default ctor
      // ----------------------------------------------------------------------
      Calculator()  {

        m_weight_power = 1.0;
        m_weight_type  = Type::Pt;

      }  // end default ctor

      // ----------------------------------------------------------------------
      //! default dtor
      // ----------------------------------------------------------------------
      ~Calculator() {};

      // ----------------------------------------------------------------------
      //! ctor accepting arguments
      // ----------------------------------------------------------------------
      Calculator(const Type::Weight weight, const double power = 1.0) {

        m_weight_power = power;
        m_weight_type  = weight;

      }  // end ctor(Type::Weight, double)

  };  // end PHEnergyCorrelator::Calculator

}  // end PHEnergyCorrelator namespace

#endif

// end ========================================================================
