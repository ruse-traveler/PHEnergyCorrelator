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
#include "PHCorrelatorAnaTools.h"
#include "PHCorrelatorAnaTypes.h"
#include "PHCorrelatorAngler.h"
#include "PHCorrelatorHistManager.h"



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
      std::vector< std::pair<float, float> > m_chrg_bins;

      // data members (phec components)
      HistManager m_manager;
      Angler      m_angler;

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
       *  spin sorting, vector will only have FOUR entries, which are
       *    [0] = integrated pt, charge (+cf bin)
       *    [1] = pt bin, integrated charge (+cf bin)
       *    [2] = charge bin, integrated pt (+cf bin)
       *    [3] = pt bin, charge bin (+cf bin)
       *
       *  If doing spin sorting, the vector will have either 4, 8, or 16
       *  entries.
       *    - `size() == 16`: pp case; entries correspond to spin-integrated,
       *       blue-only, yellow-only, and blue-and-yellow indices.
       *    - `size() == 8`: pAu case; entries correspond to spin-integrated
       *       and blue-only indices.
       *    - `size() == 4`: an unexpected spin pattern was provided, only
       *       spin-integrated case returned.
       *
       *  The order of the vector of indices will always be
       *    [0 - 3]   = spin integrated
       *    [4 - 7]   = blue beam index
       *    [8 - 11]  = yellow beam index
       *    [12 - 15] = blue and yellow index
       */
      std::vector<Type::HistIndex> GetHistIndices(const Type::Jet& jet) {

        // for pt and cf, index will correspond to what bin
        // the jet falls in PLUS integrated pt/charge bin
        Type::HistIndex base_index(0, 0, 0, 0);

        // determine pt bin
        if (m_manager.GetDoPtJetBins()) {
          for (std::size_t ipt = 0; ipt < m_ptjet_bins.size(); ++ipt) {
            if ((jet.pt >= m_ptjet_bins[ipt].first) && (jet.pt < m_ptjet_bins[ipt].second)) {
              base_index.pt = ipt;
            }
          }  // end pt bin loop
        }

        // determine cf bin
        //   - n.b. there is NO integrated bin for cf
        if (m_manager.GetDoCFJetBins()) {
          for (std::size_t icf = 0; icf < m_cfjet_bins.size(); ++icf) {
            if ((jet.cf >= m_cfjet_bins[icf].first) && (jet.cf < m_cfjet_bins[icf].second)) {
              base_index.cf = icf;
            }
          }  // end charge bin loop
        }

        // determine charge bin
        if (m_manager.GetDoChargeBins()) {
          for (std::size_t ich = 0; ich < m_chrg_bins.size(); ++ich) {
            if ((jet.charge >= m_chrg_bins[ich].first) && (jet.charge < m_chrg_bins[ich].second)) {
              base_index.chrg = ich;
            }
          }  // end charge bin loop
        }

        // By default, only add spin-integrated bin
        std::vector<std::size_t> spin_indices;
        spin_indices.push_back( HistManager::Int );

        // if needed, determine spin bins
        if (m_manager.GetDoSpinBins()) {
            switch (jet.pattern) {

              // blue up, yellow up (pp)
              case Type::PPBUYU:
                spin_indices.push_back( HistManager::BU );
                spin_indices.push_back( HistManager::YU );
                spin_indices.push_back( HistManager::BUYU );
                break;

              // blue down, yellow up (pp)
              case Type::PPBDYU:
                spin_indices.push_back( HistManager::BD );
                spin_indices.push_back( HistManager::YU );
                spin_indices.push_back( HistManager::BDYU );
                break;

              // blue up, yellow down (pp)
              case Type::PPBUYD:
                spin_indices.push_back( HistManager::BU );
                spin_indices.push_back( HistManager::YD );
                spin_indices.push_back( HistManager::BUYD );
                break;

              // blue down, yellow down (pp)
              case Type::PPBDYD:
                spin_indices.push_back( HistManager::BD );
                spin_indices.push_back( HistManager::YD );
                spin_indices.push_back( HistManager::BDYD );
                break;

              // blue up (pAu)
              case Type::PABU:
                spin_indices.push_back( HistManager::BU );
                break;

              // blue down (pAu)
              case Type::PABD:
                spin_indices.push_back( HistManager::BD );
                break;

              // by default, only add integrated
              default:
                break;

          }
        }

        // now assemble list of indices to fill
        std::vector<Type::HistIndex> indices;
        for (std::size_t isp = 0; isp < spin_indices.size(); ++isp) {

          // integrated everything (except cf)
          indices.push_back(
            Type::HistIndex(
              m_ptjet_bins.size(),
              base_index.cf,
              m_chrg_bins.size(),
              spin_indices[isp]
            )
          );

          // binned pt, integrated charge
          indices.push_back(
            Type::HistIndex(
              base_index.pt,
              base_index.cf,
              m_chrg_bins.size(),
              spin_indices[isp]
            )
          );

          // binned charge, integrated pt
          indices.push_back(
            Type::HistIndex(
              m_ptjet_bins.size(),
              base_index.cf,
              base_index.chrg,
              spin_indices[isp]
            )
          );

          // binned everything
          indices.push_back(
            Type::HistIndex(
              base_index.pt,
              base_index.cf,
              base_index.chrg,
              spin_indices[isp]
            )
          );
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
      HistManager& GetManager() {return m_manager;}

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
      //! Set jet charge bins
      // ----------------------------------------------------------------------
      void SetChargeBins(const std::vector< std::pair<float, float> >& bins) {

        // copy bins to member
        m_chrg_bins.resize( bins.size() );
        std::copy(bins.begin(), bins.end(), m_chrg_bins.begin());

        // update hist manager and exit
        m_manager.DoChargeBins( m_chrg_bins.size() );
        return;

      }  // end 'SetChargeBins(std::vector<std::pair<float, float>>&)'

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

        // calculate jet, cst kinematics --------------------------------------

        // get jet 4-momenta
        TLorentzVector vecJet4  = Tools::GetJetLorentz(jet, false);
        TLorentzVector unitJet4 = Tools::GetJetLorentz(jet, true);

        // get cst 4-momenta
        std::pair<TLorentzVector, TLorentzVector> vecCst4 = std::make_pair(
          Tools::GetCstLorentz(csts.first, jet.pt, false),
          Tools::GetCstLorentz(csts.second, jet.pt, false)
        );

        // get average of cst 3-vectors
        TVector3 vecAvgCst3 = Tools::GetWeightedAvgVector(
          vecCst4.first.Vect(),
          vecCst4.second.Vect(),
          false
        );
        TVector3 unitAvgCst3 = Tools::GetWeightedAvgVector(
          vecCst4.first.Vect(),
          vecCst4.second.Vect(),
          true
        );

        // collins & boer-mulders angle calculations --------------------------

        // (0) get beam and spin directions
        //     - components:
        //         first  = blue beam/spin
        //         second = yellow beam/spin
        //     - n.b. for spin pattern >= 4, the
        //       yellow spin is randomized
        std::pair<TVector3, TVector3> vecBeam3 = Tools::GetBeams();
        std::pair<TVector3, TVector3> vecSpin3 = Tools::GetSpins( jet.pattern );

        // (1) get vectors normal to the jet-beam plane
        std::pair<TVector3, TVector3> normJetBeam3 = std::make_pair(
          ( vecBeam3.first.Cross(unitJet4.Vect()) ).Unit(),
          ( vecBeam3.second.Cross(unitJet4.Vect()) ).Unit()
        );

        // (2) get vectors normal to the spin-beam plane
        std::pair<TVector3, TVector3> normJetSpin = std::make_pair(
          ( vecBeam3.first.Cross(vecSpin3.first) ).Unit(),
          ( vecBeam3.second.Cross(vecSpin3.second) ).Unit()
        );

        // (3) get phiSpin: angles between jet-beam and spin-beam planes
	//     - between [0, 2pi] by definition 
        double phiSpinBlue = m_angler.GetTwoPlaneAngle(normJetBeam3.first, normJetSpin.first, vecSpin3.first);
        double phiSpinYell = m_angler.GetTwoPlaneAngle(normJetBeam3.second, normJetSpin.second, vecSpin3.second);

        // (4) get vector normal to hadron average-jet plane
        TVector3 normHadJet3 = ( unitJet4.Vect().Cross(unitAvgCst3) ).Unit();

        // (5) get phiHadron: angle between jet-beam and hadron-jet planes
        //     - between [0, 2pi] by definition
        double phiHadBlue = m_angler.GetTwoPlaneAngle(normJetBeam3.first, normHadJet3, unitAvgCst3);
        double phiHadYell = m_angler.GetTwoPlaneAngle(normJetBeam3.second, normHadJet3, unitAvgCst3);

        // (6) now calculate phiColl: phiSpin - phiHadron,
        //     - constrain to [0, 2pi)
        double phiCollBlue = m_angler.GetCollinsAngle(phiSpinBlue, phiHadBlue);
        double phiCollYell = m_angler.GetCollinsAngle(phiSpinYell, phiHadYell);

        // (7) now calculate phiBoer: phiSpin - (2 * phiHadron),
        double phiBoerBlue = m_angler.GetBoerMuldersAngle(phiSpinBlue, phiHadBlue);
        double phiBoerYell = m_angler.GetBoerMuldersAngle(phiSpinYell, phiHadYell);
 
        // calculate eec quantities -------------------------------------------

        // get EEC weights
        std::pair<double, double> cst_weights = std::make_pair(
          GetCstWeight(vecCst4.first, vecJet4),
          GetCstWeight(vecCst4.second, vecJet4)
        );

        // and then calculate RL (dist b/n cst.s for EEC) and overall EEC weight
        const double dist    = Tools::GetCstDist(csts);
        const double weight  = cst_weights.first * cst_weights.second * evt_weight;

        // fill histograms ---------------------------------------------------=

        // fill histograms if needed
        if (m_manager.GetDoEECHists()) {

          // grab hist indices
          std::vector<Type::HistIndex> indices = GetHistIndices(jet);

          // collect quantities to be histogrammed
          Type::HistContent content(weight, dist);
          if (m_manager.GetDoSpinBins()) {
            content.phiCollB = phiCollBlue;
            content.phiCollY = phiCollYell;
            content.phiBoerB = phiBoerBlue;
            content.phiBoerY = phiBoerYell;
            content.spinB    = vecSpin3.first.Y();
            content.spinY    = vecSpin3.second.Y();
            content.pattern  = jet.pattern;
          }

          // fill spin-integrated histograms
          for (std::size_t idx = 0; idx < Const::NBinsPerSpin(); ++idx) {
            m_manager.FillEECHists(indices[idx], content);
          }

          // if needed, fill spin sorted histograms
          if (m_manager.GetDoSpinBins() && (indices.size() > Const::BlueSpinStart())) {

            // fill blue spins
            for (
              std::size_t idx = Const::BlueSpinStart();
              idx < Const::YellSpinStart();
              ++idx
            ) {
              m_manager.FillEECHists(indices[idx], content);
            }

            // fill yellow and both spins
            if (indices.size() > Const::YellSpinStart()) {
              for (
                std::size_t idx = Const::YellSpinStart();
                idx < indices.size();
                ++idx
              ) {
                m_manager.FillEECHists(indices[idx], content);
              }
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
