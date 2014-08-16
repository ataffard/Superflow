#include "Superflow/DileptonVariables.h"

#include <cmath> // sqrt
#include <math.h> // cos, fabs
#include <cassert>

using sflow::DileptonVariables;

//-----------------------------------------
DileptonVariables sflow::computeDileptonVariables(const LeptonVector &leptons, const Susy::Met *met,
    const JetVector &jets)
{
    DileptonVariables v;
    if (leptons.size() > 1) {
        bool sorted = leptons[0]->Pt() > leptons[1]->Pt();
        Susy::Lepton &l0 = (sorted ? *leptons[0] : *leptons[1]);
        Susy::Lepton &l1 = (sorted ? *leptons[1] : *leptons[0]);
        bool isEl0(l0.isEle()), isEl1(l1.isEle()), isMu0(l0.isMu()), isMu1(l1.isMu());
        assert(isEl0 != isMu0 && isEl1 != isMu1); // assuming we're only dealing with electrons or muons
        v.isMu0 = isMu0;
        v.isMu1 = isMu1;
        v.q0 = l0.q;
        v.q1 = l1.q;
        v.pt0 = l0.Pt();
        v.pt1 = l1.Pt();
        v.eta0 = l0.Eta();
        v.eta1 = l1.Eta();
        v.phi0 = l0.Phi();
        v.phi1 = l1.Phi();
        TLorentzVector ll(l0 + l1);
        v.mll = ll.M();
        v.mcoll01 = sflow::computeCollinearMzLepTau(l0, l1, met->lv());
        v.mcoll10 = sflow::computeCollinearMzLepTau(l1, l0, met->lv());
        v.mcoll = sflow::computeCollinearMzTauTau(l0, l1, met->lv());
        v.detall = fabs(l0.Eta() - l1.Eta());
        LeptonVector lepts;
        lepts.push_back(&l0);
        lepts.push_back(&l1);
        v.metrel = SusyNtTools::getMetRel(met, lepts, jets);
        v.mt0 = sflow::transverseMass(l0, met->lv());
        v.mt1 = sflow::transverseMass(l1, met->lv());
        v.mtllmet = transverseMass(ll, met->lv());
        v.met = met->Et;
        v.metPhi = met->phi;
        v.numCentralLightJets = jets.size();
        if (jets.size() > 0) { const Susy::Jet& j = (*jets[0]); v.j0pt = j.Pt(); v.j0eta = j.Eta(); v.j0phi = j.Phi(); }
        if (jets.size() > 1) { const Susy::Jet& j = (*jets[1]); v.j1pt = j.Pt(); v.j1eta = j.Eta(); v.j1phi = j.Phi(); }
        if (jets.size() > 2) { const Susy::Jet& j = (*jets[2]); v.j2pt = j.Pt(); v.j2eta = j.Eta(); v.j2phi = j.Phi(); }
    }
    return v;
}
//-----------------------------------------
float sflow::transverseMass(const TLorentzVector &lep, const TLorentzVector &met)
{
    return std::sqrt(2.0 * lep.Pt() * met.Et() *(1 - cos(lep.DeltaPhi(met))));
}
//-----------------------------------------
bool DileptonVariables::passAllSelectionCriteria() const
{
    return passLeptonPt() && passJetVeto() && passDeltaPhiLl() && passDeltaPhiL1Met();
}
//-----------------------------------------
bool DileptonVariables::passL0Pt() const { return pt0 > 30.0; }
bool DileptonVariables::passL1Pt() const { return pt1 > 20.0; }
bool DileptonVariables::passLeptonPt() const { return passL0Pt() && passL1Pt(); }
bool DileptonVariables::passJetVeto() const { return numCentralLightJets == 0; }
//-----------------------------------------
bool DileptonVariables::passDeltaPhiLl() const
{
    float dphill = fabs(TVector2::Phi_mpi_pi(TVector2(1.0, 0.0).Rotate(phi0)
        .DeltaPhi(TVector2(1.0, 0.0).Rotate(phi1))));
    return dphill > 2.5;
}
//-----------------------------------------
bool DileptonVariables::passDeltaPhiL1Met() const
{
    float dphil1met = fabs(TVector2::Phi_mpi_pi(TVector2(1.0, 0.0).Rotate(phi1)
        .DeltaPhi(TVector2(1.0, 0.0).Rotate(metPhi))));
    return dphil1met < 0.5;
}
//-----------------------------------------
float sflow::computeCollinearMzLepTau(const TLorentzVector &l0,
    const TLorentzVector &l1,
    const TLorentzVector &met)
{
    return std::sqrt(2.0 * l0.Pt() *
        (l1.Pt() + met.Et()) *
        (cosh(l0.Eta() - l1.Eta()) - cos(l0.DeltaPhi(l1))));
}
//-----------------------------------------
float sflow::computeCollinearMzTauTau(const TLorentzVector &l0,
    const TLorentzVector &l1,
    const TLorentzVector &met)
{
    float px0(l0.Px()), py0(l0.Py());
    float px1(l1.Px()), py1(l1.Py());
    float pxm(met.Px()), pym(met.Py());
    float num(px0*py1 - py0*px1);
    float den1(py1*pxm - px1*pym + px0*py1 - py0*px1);
    float den2(px0*pym - py0*pxm + px0*py1 - py0*px1);
    float x1 = (den1 != 0.0 ? (num / den1) : 0.0);
    float x2 = (den2 != 0.0 ? (num / den2) : 0.0);
    bool kinematicallyPossible(x1*x2 > 0.0);
    return (kinematicallyPossible ? (l0 + l1).M() / std::sqrt(x1*x2) : -1.0);
}
//-----------------------------------------
