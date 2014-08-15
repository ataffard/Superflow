#pragma once

#include <vector>
#include <cstddef> // size_t

#include "TLorentzVector.h"
#include "TVector2.h"

#include "SusyNtuple/SusyNt.h" // Lepton, Jet, Met, and all that
#include "SusyNtuple/SusyDefs.h" // LeptonVector, JetVector and all that
#include "SusyNtuple/SusyNtTools.h"

// namespace Susy
// {
//   class Lepton;
//   class Jet;
//   class Met;
// }
// typedef std::vector<Susy::Lepton*> LeptonVector;
// typedef std::vector<Susy::Jet*> JetVector;

namespace sflow
{
  /// Kinematic variables for the sflow study
  /**
     The highest-pt lepton is stored as l0, and the second one as l1.
     Only works with electrons and muons (not taus).

     davide.gerbaudo@gmail.com
     Jun 2014
   */
  struct DileptonVariables {
    DileptonVariables() { reset(); }
    bool isMu0, isMu1;
    bool hasFiredTrig, hasTrigMatch;
    float q0, q1;
    float pt0, pt1;
    float phi0, phi1;
    float eta0, eta1;
    float mll, detall;
    float mcoll01, mcoll10, mcoll;
    size_t numCentralLightJets;
    float j0pt,  j1pt,  j2pt;
    float j0eta, j1eta, j2eta;
    float j0phi, j1phi, j2phi;
    float mt0, mt1, mtllmet;
    float met, metPhi;
    float metrel;
    bool isSs() const { return (q0*q1)>0.0; }
    bool isOs() const { return !isSs(); }
    bool isEe() const { return (!isMu0 && !isMu1); }
    bool isMumu() const { return (isMu0 && isMu1); }
    bool isEmu() const { return (!isMu0 && isMu1); }
    bool isMue() const { return (isMu0 && !isMu1); }
    float mtmin() const { return mt0<mt1 ? mt0 : mt1; }
    float mtmax() const { return mt0>mt1 ? mt0 : mt1; }
    /// reset all variables to their default value
    void reset() {
      isMu0 = isMu1 = false;
      hasFiredTrig = hasTrigMatch = false;
      q0 = q1 = 0.0;
      pt0 = pt1 =  0.0;
      eta0 = eta1 = 0.0;
      phi0 = phi1 = 0.0;
      mll = detall = 0.0;
      mcoll01 = mcoll10 = mcoll = 0.0;
      numCentralLightJets = 0;
      j0pt  = j1pt  = j2pt  = 0.0;
      j0eta = j1eta = j2eta = 0.0;
      j0phi = j1phi = j2phi = 0.0;
      mt0 = mt1 = mtllmet = 0.0;
      met = metPhi = 0.0;
      metrel = 0.0;
    }
    /// pass all selection criteria below
    bool passAllSelectionCriteria() const;
    bool passL0Pt() const;
    bool passL1Pt() const;
    bool passLeptonPt() const;
    bool passJetVeto() const;
    bool passDeltaPhiLl() const;
    bool passDeltaPhiL1Met() const;
  };
  /// compute and assign all DilepVars attributes
  DileptonVariables computeDileptonVariables(const LeptonVector &leptons, const Susy::Met *met,
                                             const JetVector &jets);
  /// transverse W mass variable
  float transverseMass(const TLorentzVector &lep, const TLorentzVector &met);
  /// invariant mass under the assumption that one lepton is from a leptonic tau
  /**
     The second lepton (l1) is assumed to be the one from the tau.
     See definition in  hep-ph/1405.4545.
   */
  float computeCollinearMzLepTau(const TLorentzVector &l0,
                                 const TLorentzVector &l1,
                                 const TLorentzVector &met);
  /// invariant mass under the assumption that both leptons are from leptonic taus
  /**
     Apply the collinear approximation to both decays.
     Re-written based on HWWlvlvCode::calculate_METBasedVariables
   */
  float computeCollinearMzTauTau(const TLorentzVector &l0,
                                 const TLorentzVector &l1,
                                 const TLorentzVector &met);
}
