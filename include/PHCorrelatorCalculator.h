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
#include <limits>
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
      /*! Returns a vector of indices of histograms to be filled. If not doing
       *  spin sorting, vector will only have ONE entry, which is the pt (and
       *  maybe CF) index.
       *
       *  If doing spin sorting, the vector will have either 1, 2, or 4
       *  entries.
       *    - `size() == 4`: pp case; entries correspond to spin-integrated,
       *       blue-only, yellow-only, and blue-and-yellow indices.
       *    - `size() == 2`: pAu case; entries correspond to spin-integrated
       *       and blue-only indices.
       *    - `size() == 1`: an unexpected spin pattern was provided, only
       *       spin-integrated case returned.
       *
       *  The order of the vector will always be
       *    [0] = integrated
       *    [1] = blue beam index
       *    [2] = yellow beam index
       *    [3] = blue and yellow index
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

        // By default, add spin-integrated bin
        std::vector<Type::HistIndex> indices;
        indices.push_back( Type::HistIndex(iptcf.pt, iptcf.cf, Manager::Int) );

        // if needed, determine spin bins
        if (m_manager.GetDoSpinBins()) {
            switch (jet.pattern) {

              // blue up, yellow up (pp)
              case Type::PPBUYU:
                indices.push_back( Type::HistIndex(iptcf.pt, iptcf.cf, Manager::BU) );
                indices.push_back( Type::HistIndex(iptcf.pt, iptcf.cf, Manager::YU) );
                indices.push_back( Type::HistIndex(iptcf.pt, iptcf.cf, Manager::BUYU) );
                break;

              // blue down, yellow up (pp)
              case Type::PPBDYU:
                indices.push_back( Type::HistIndex(iptcf.pt, iptcf.cf, Manager::BD) );
                indices.push_back( Type::HistIndex(iptcf.pt, iptcf.cf, Manager::YU) );
                indices.push_back( Type::HistIndex(iptcf.pt, iptcf.cf, Manager::BDYU) );
                break;

              // blue up, yellow down (pp)
              case Type::PPBUYD:
                indices.push_back( Type::HistIndex(iptcf.pt, iptcf.cf, Manager::BU) );
                indices.push_back( Type::HistIndex(iptcf.pt, iptcf.cf, Manager::YD) );
                indices.push_back( Type::HistIndex(iptcf.pt, iptcf.cf, Manager::BUYD) );
                break;

              // blue down, yellow down (pp)
              case Type::PPBDYD:
                indices.push_back( Type::HistIndex(iptcf.pt, iptcf.cf, Manager::BD) );
                indices.push_back( Type::HistIndex(iptcf.pt, iptcf.cf, Manager::YD) );
                indices.push_back( Type::HistIndex(iptcf.pt, iptcf.cf, Manager::BDYD) );
                break;

              // blue up (pAu)
              case Type::PABU:
                indices.push_back( Type::HistIndex(iptcf.pt, iptcf.cf, Manager::BU) );
                break;

              // blue down (pAu)
              case Type::PABD:
                indices.push_back( Type::HistIndex(iptcf.pt, iptcf.cf, Manager::BD) );
                break;

              // by default, only add integrated
              default:
                break;

          }
        }
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
        TLorentzVector vecJet4  = Tools::GetJetLorentz(jet, false);
        TLorentzVector unitJet4 = Tools::GetJetLorentz(jet, true);

        // get cst 4-momenta
        std::pair<TLorentzVector, TLorentzVector> vecCst4 = std::make_pair(
          Tools::GetCstLorentz(csts.first, jet.pt, false),
          Tools::GetCstLorentz(csts.second, jet.pt, false)
        );
        std::pair<TLorentzVector, TLorentzVector> unitCst4 = std::make_pair(
          Tools::GetCstLorentz(csts.first, jet.pt, true),
          Tools::GetCstLorentz(csts.second, jet.pt, true)
        );

        // get average of cst 3-vectors
        TVector3 avgCst3 = Tools::GetWeightedAvgVector(
          vecCst4.first.Vect(),
          vecCst4.second.Vect()
        );

        // get spins
        std::pair<TVector3, TVector3> vecSpin3 = Tools::GetSpins( jet.pattern );

        // now calculate vectors normal to hadron-spin and jet-spin planes
        std::pair<TVector3, TVector3> normHadSpin3 = std::make_pair(
          ( vecSpin3.first.Cross(avgCst3) ).Unit(),
          ( vecSpin3.second.Cross(avgCst3) ).Unit()
        );
        std::pair<TVector3, TVector3> normJetSpin3 = std::make_pair(
          ( vecSpin3.first.Cross(unitJet4.Vect()) ).Unit(),
          ( vecSpin3.second.Cross(unitJet4.Vect()) ).Unit()
        );

        // next calculate vectors normal to hadron-hadron and hadron-jet planes
        TVector3 normHadHad3 = ( unitCst4.first.Vect().Cross(unitCst4.second.Vect()) ).Unit();
        TVector3 normJetHad3 = ( unitJet4.Vect().Cross(avgCst3) ).Unit();

        // and finally, compute angles wrt to spins
        double phiHadBlue = acos( normHadSpin3.first.Dot(normHadHad3) );
        double phiHadYell = acos( normHadSpin3.second.Dot(normHadHad3) );
        double phiJetBlue = acos( normJetSpin3.first.Dot(normJetHad3) );
        double phiJetYell = acos( normJetSpin3.first.Dot(normJetHad3) );
        if (phiHadBlue > TMath::PiOver2()) phiHadBlue -= TMath::Pi();
        if (phiHadYell > TMath::PiOver2()) phiHadYell -= TMath::Pi();
        if (phiJetBlue > TMath::PiOver2()) phiJetBlue -= TMath::Pi();
        if (phiJetYell > TMath::PiOver2()) phiJetYell -= TMath::Pi();

        // now get EEC weights
        std::pair<double, double> cst_weights = std::make_pair(
          GetCstWeight(vecCst4.first, vecJet4),
          GetCstWeight(vecCst4.second, vecJet4)
        );

        // then calculate RL (dist b/n cst.s for EEC) and overall EEC weight
        const double dist    = Tools::GetCstDist(csts);
        const double weight  = cst_weights.first * cst_weights.second * evt_weight;

        // fill histograms if needed
        //   - FIXME this can be simplified...
        if (m_manager.GetDoEECHists()) {

          // grab hist indices
          std::vector<Type::HistIndex> indices = GetHistIndices(jet);

          // fill spin integrated case
          Type::HistContent content_int(weight, dist);
          m_manager.FillEECHists(indices[0], content_int);

          // if needed, fill spin sorted histograms
          if (m_manager.GetDoSpinBins() && (indices.size() > 1)) {

            // fill blue-only case
            Type::HistContent content_blu(
              weight,
              dist,
              phiHadBlue,
              std::numeric_limits<double>::max(),
              phiJetBlue,
              std::numeric_limits<double>::max(),
              vecSpin3.first.Y(),
              vecSpin3.second.Y(),
              jet.pattern
            );
            m_manager.FillEECHists(indices[1], content_blu);

            // fill yellow-only, both cases
            if (indices.size() > 2) {
              Type::HistContent content_yel(
                weight,
                dist,
                std::numeric_limits<double>::max(),
                phiHadYell,
                std::numeric_limits<double>::max(),
                phiJetYell,
                vecSpin3.first.Y(),
                vecSpin3.second.Y(),
                jet.pattern
              );
              Type::HistContent content_both(
                weight,
                dist,
                phiHadBlue,
                phiHadYell,
                phiJetBlue,
                phiJetYell,
                vecSpin3.first.Y(),
                vecSpin3.second.Y(),
                jet.pattern
              );
              m_manager.FillEECHists(indices[2], content_yel);
              m_manager.FillEECHists(indices[3], content_both);
            }
          }  // end spin hist filling
        }  // end hist filing
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
