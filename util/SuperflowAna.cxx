// SuperflowAna.cxx
//

#include <cstdlib>
#include <cmath>
#include <fstream> 
#include <iostream>
#include <string>
#include <getopt.h>

#include "TChain.h"
#include "TVectorD.h"
//#include "Cintex/Cintex.h"

#include "SusyNtuple/ChainHelper.h"
#include "SusyNtuple/string_utils.h"

#include "Superflow/Superflow.h"
#include "Superflow/Superlink.h"
#include "Superflow/Cut.h"
#include "Superflow/StringTools.h"
#include "Superflow/PhysicsTools.h"
#include "Superflow/LeptonTruthDefinitions.h"

#include "Mt2/mt2_bisect.h"

using namespace std;
using namespace sflow;

// constants
const double GeV_to_MeV = 1000.0;

// function prototypes
void print_usage(const char *exeName);
void read_options(int argc, char* argv[], TChain* chain, int& n_skip_, int& num_events_, string& sample_, SuperflowRunMode& run_mode_, SusyNtSys& nt_sysnt_);

int main(int argc, char* argv[])
{
    // START read-in
    int n_skip_ = 0;
    int num_events_ = -1;
    string sample_;
    SuperflowRunMode run_mode = SuperflowRunMode::nominal;
    SusyNtSys nt_sys_ = NtSys_NOM;

    TChain* chain = new TChain("susyNt");
    chain->SetDirectory(0);

    read_options(argc, argv, chain, n_skip_, num_events_, sample_, run_mode, nt_sys_); // defined below
    // END read-in

    Superflow* cutflow = new Superflow(); // initialize the cutflow
    cutflow->setSampleName(sample_);
    cutflow->setRunMode(run_mode);
    cout << "Analysis    Total Entries: " << chain->GetEntries() << endl;

    if (run_mode == SuperflowRunMode::single_event_syst) cutflow->setSingleEventSyst(nt_sys_);

    // START Setup cuts
    // START Setup cuts
    // START Setup cuts

    *cutflow << CutName("exactly two signal leptons") << [](Superlink* sl) -> bool {
        please return sl->leptons->size() == 2;
    }; // First cut to help speed.

    *cutflow << CutName("opposite sign") << [](Superlink* sl) -> bool {
        return (sl->leptons->at(0)->q * sl->leptons->at(1)->q < 0);
    };

    *cutflow << CutName("m_ll > 20 GeV") << [](Superlink* sl) -> bool {
        return (*sl->leptons->at(0) + *sl->leptons->at(1)).M() > 20.0;
    };

    *cutflow << CutName("leading lepton Pt > 25 GeV") << [](Superlink* sl) -> bool {
        return sl->leptons->at(0)->Pt() > 25.0;
    };

    *cutflow << CutName("subleading lepton Pt > 10 GeV") << [](Superlink* sl) -> bool {
        please return sl->leptons->at(1)->Pt() > 10.0;
    };

    *cutflow << CutName("muon eta < 2.4") << [](Superlink* sl) -> bool {
        return (!sl->leptons->at(0)->isMu() || abs(sl->leptons->at(0)->Eta()) < 2.4)
            && (!sl->leptons->at(1)->isMu() || abs(sl->leptons->at(1)->Eta()) < 2.4);
    };

    *cutflow << CutName("tau veto") << [](Superlink* sl) -> bool {
        return sl->taus->size() == 0;
    };

    *cutflow << CutName("pass dilepton trigger") << [](Superlink* sl) -> bool {
        return sl->dileptonTrigger->passDilTrig(*sl->leptons, sl->met->lv().Pt(), sl->nt->evt());
    };

    *cutflow << CutName("prompt leptons") << [](Superlink* sl) -> bool {
        bool pass_ = true;

        if (sl->isMC && !sl->doFake) {
            for (int l_ = 0; l_ < sl->leptons->size(); l_++) {
                bool isReal = sl->leptons->at(l_)->truthType == LeptonTruthType::PROMPT;

                bool isChargeFlip = sl->leptons->at(l_)->isEle()
                    && static_cast<Electron*>(sl->leptons->at(l_))->isChargeFlip;

                if (!isReal || isChargeFlip) pass_ = false;
            }
        }
        return pass_;
    };

    *cutflow << CutName("forward jets veto") << [](Superlink* sl) -> bool {
        return SusyNtTools::numberOfFJets(*sl->jets) == 0;
    };

    // END Setup cuts
    // END Setup cuts
    // END Setup cuts

    // GAP //
    // GAP //
    // GAP //

    // START Setup output trees
    // START Setup output trees
    // START Setup output trees

    *cutflow << NewVar("event weight"); {
        *cutflow << HFTname("eventweight");
        *cutflow << [](Superlink* sl, var_double*) -> double { please return sl->weights->product(); };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("run number"); {
        *cutflow << HFTname("runNumber");
        *cutflow << [](Superlink* sl, var_int*) -> int { return sl->nt->evt()->run; };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("event number"); {
        *cutflow << HFTname("eventNumber");
        *cutflow << [](Superlink* sl, var_int*) -> int { return sl->nt->evt()->event; };
        *cutflow << SaveVar();
    }

    // LEPTONS
    // LEPTONS
    // LEPTONS

    *cutflow << NewVar("is e + e"); {
        *cutflow << HFTname("isEE");
        *cutflow << [](Superlink* sl, var_bool*) -> bool { return sl->leptons->at(0)->isEle() && sl->leptons->at(1)->isEle(); };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("is e + mu"); {
        *cutflow << HFTname("isEMU");
        *cutflow << [](Superlink* sl, var_bool*) -> bool { return sl->leptons->at(0)->isEle() ^ sl->leptons->at(1)->isEle(); };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("is mu + mu"); {
        *cutflow << HFTname("isMUMU");
        *cutflow << [](Superlink* sl, var_bool*) -> bool { return sl->leptons->at(0)->isMu() && sl->leptons->at(1)->isMu(); };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("is opposite-sign"); {
        *cutflow << HFTname("isOS");
        *cutflow << [](Superlink* sl, var_bool*) -> bool { return sl->leptons->at(0)->q * sl->leptons->at(1)->q < 0; };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("lepton-1 Pt"); {
        *cutflow << HFTname("lept1Pt");
        *cutflow << [](Superlink* sl, var_float*) -> double { return sl->leptons->at(0)->Pt() * GeV_to_MeV; };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("lepton-2 Pt"); {
        *cutflow << HFTname("lept2Pt");
        *cutflow << [](Superlink* sl, var_float*) -> double { return sl->leptons->at(1)->Pt() * GeV_to_MeV; };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("lepton-1 Eta"); {
        *cutflow << HFTname("lept1Eta");
        *cutflow << [](Superlink* sl, var_float*) -> double { return sl->leptons->at(0)->Eta(); };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("lepton-2 Eta"); {
        *cutflow << HFTname("lept2Eta");
        *cutflow << [](Superlink* sl, var_float*) -> double { return sl->leptons->at(1)->Eta(); };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("lepton-1 Phi"); {
        *cutflow << HFTname("lept1Phi");
        *cutflow << [](Superlink* sl, var_float*) -> double { return sl->leptons->at(0)->Phi(); };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("lepton-2 Phi"); {
        *cutflow << HFTname("lept2Phi");
        *cutflow << [](Superlink* sl, var_float*) -> double { return sl->leptons->at(1)->Phi(); };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("lepton-1 E"); {
        *cutflow << HFTname("lept1E");
        *cutflow << [](Superlink* sl, var_float*) -> double { return sl->leptons->at(0)->E() * GeV_to_MeV; };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("lepton-2 E"); {
        *cutflow << HFTname("lept2E");
        *cutflow << [](Superlink* sl, var_float*) -> double { return sl->leptons->at(1)->E() * GeV_to_MeV; };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("lepton-1 charge"); {
        *cutflow << HFTname("lept1q");
        *cutflow << [](Superlink* sl, var_int*) -> int { return sl->leptons->at(0)->q; };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("lepton-2 charge"); {
        *cutflow << HFTname("lept2q");
        *cutflow << [](Superlink* sl, var_int*) -> int { return sl->leptons->at(1)->q; };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("lepton-1 flavor"); { // 0=el, 1=mu, see HistFitterTree.h
        *cutflow << HFTname("lept1Flav");
        *cutflow << [](Superlink* sl, var_int*) -> int { return sl->leptons->at(0)->isEle() ? 0 : 1; };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("lepton-2 flavor"); {
        *cutflow << HFTname("lept2Flav");
        *cutflow << [](Superlink* sl, var_int*) -> int { return sl->leptons->at(1)->isEle() ? 0 : 1; };
        *cutflow << SaveVar();
    }

    // JETS
    // JETS
    // JETS

    *cutflow << NewVar("number of central light jets"); {
        *cutflow << HFTname("L2nCentralLightJets");
        *cutflow << [](Superlink* sl, var_int*) -> int {
            return SusyNtTools::numberOfCLJets(*sl->jets, sl->jvfTool, sl->nt_sys, sl->anaType);
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("number of central b jets"); {
        *cutflow << HFTname("L2nCentralBJets");
        *cutflow << [](Superlink* sl, var_int*) -> int { return SusyNtTools::numberOfCBJets(*sl->jets); };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("number of forward jets"); {
        *cutflow << HFTname("L2nForwardJets");
        *cutflow << [](Superlink* sl, var_int*) -> int { return SusyNtTools::numberOfFJets(*sl->jets); };
        *cutflow << SaveVar();
    }

    // MET
    // MET
    // MET

    *cutflow << NewVar("transverse missing energy (Et)"); {
        *cutflow << HFTname("met");
        *cutflow << [](Superlink* sl, var_float*) -> double { return sl->met->Et * GeV_to_MeV; };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("transverse missing energy (Phi)"); {
        *cutflow << HFTname("metPhi");
        *cutflow << [](Superlink* sl, var_float*) -> double { return sl->met->phi; };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("Etmiss Rel"); {
        *cutflow << HFTname("L2METrel");
        *cutflow << [](Superlink* sl, var_float*) -> double { return SusyNtTools::getMetRel(sl->met, *sl->leptons, *sl->jets); };
        *cutflow << SaveVar();
    }


    // VARS
    // VARS
    // VARS

    TLorentzVector local_ll;

    *cutflow << NewVar("mass of di-lepton system, M_ll"); {
        *cutflow << HFTname("L2Mll");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            local_ll = (*sl->leptons->at(0) + *sl->leptons->at(1));
            return local_ll.M() * GeV_to_MeV;
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("delta Eta(ll)"); {
        *cutflow << HFTname("deltaEtaLl");
        *cutflow << [](Superlink* sl, var_float*) -> double { return abs(sl->leptons->at(0)->Eta() - sl->leptons->at(1)->Eta()); };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("mTWW or mtllmet"); {
        *cutflow << HFTname("mtllmet");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            TLorentzVector met_ = sl->met->lv();
            return sqrt(2.0 * local_ll.Pt() * met_.Et() *(1 - cos(local_ll.DeltaPhi(met_)))) * GeV_to_MeV;
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("Ht (m_Eff)"); {
        *cutflow << HFTname("Ht");
        *cutflow << [](Superlink* sl, var_float*) -> double {
            double ht = 0.0;

            ht += sl->leptons->at(0)->Pt() + sl->leptons->at(1)->Pt();
            ht += sl->met->Et;
            for (int i = 0; i < sl->jets->size(); i++) {
                if (sl->jets->at(i)->Pt() > 20.0) {
                    ht += sl->jets->at(i)->Pt();
                }
            }
            return ht * GeV_to_MeV;
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("stransverse mass"); {
        *cutflow << HFTname("MT2");
        *cutflow << [](Superlink* sl, var_float*) -> double {

            mt2_bisect::mt2 mt2_event;
            double *pa, *pb, *pmiss;
            pa = new double[3]; pa[0] = sl->leptons->at(0)->M(); pa[1] = sl->leptons->at(0)->Px(), pa[2] = sl->leptons->at(0)->Py();
            pb = new double[3]; pb[0] = sl->leptons->at(1)->M(); pb[1] = sl->leptons->at(1)->Px(), pb[2] = sl->leptons->at(1)->Py();
            pmiss = new double[3]; pmiss[0] = 0.0; pmiss[1] = sl->met->Et * cos(sl->met->phi); pmiss[2] = sl->met->Et * sin(sl->met->phi);

            mt2_event.set_momenta(pa, pb, pmiss);
            mt2_event.set_mn(0.0); // LSP mass = 0 is Generic

            double mt2_ = mt2_event.get_mt2();
            // SUPRESSED messages "Deltasq_high not found at event 0" (~1 per 10000 events)

            delete[] pa;
            delete[] pb;
            delete[] pmiss;

            return mt2_ * GeV_to_MeV;
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("collinear mass"); {
        *cutflow << HFTname("mcoll");
        *cutflow << [](Superlink* sl, var_float*) -> double { 
            return PhysicsTools::mColl(sl->leptons->at(0), sl->leptons->at(1), sl->met->lv()) * GeV_to_MeV; };
        *cutflow << SaveVar();
    }

    TVector2 metCorr;

    *cutflow << NewVar("ETmiss Corr (magnitude)"); {
        *cutflow << HFTname("metCorr");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            double met_x = sl->met->lv().Px() - sl->met->softTerm_etx;
            double met_y = sl->met->lv().Py() - sl->met->softTerm_ety;
            metCorr = TVector2(met_x, met_y);

            return metCorr.Mod() * GeV_to_MeV;
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("ETmiss Corr (Phi)"); {
        *cutflow << HFTname("metCorrPhi");
        *cutflow << [&](Superlink* sl, var_float*) -> double { return metCorr.Phi(); };
        *cutflow << SaveVar();
    }

    Susy::Met corrMet;

    *cutflow << NewVar("Etmiss Rel Corr"); {
        *cutflow << HFTname("METrelCorr");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            corrMet = *sl->met;
            corrMet.Et = metCorr.Mod();
            corrMet.phi = metCorr.Phi();

            please return SusyNtTools::getMetRel(&corrMet, *sl->leptons, *sl->jets) * GeV_to_MeV;
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("collinear mass Corr"); {
        *cutflow << HFTname("mcollCorr");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            TLorentzVector metLV = TLorentzVector(metCorr.Mod(), 0, metCorr.Phi(), metCorr.Mod());
            return PhysicsTools::mColl(sl->leptons->at(0), sl->leptons->at(1), metLV) * GeV_to_MeV;
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("mT l1 Corr"); {
        *cutflow << HFTname("mTl1Corr");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            please return PhysicsTools::mT(sl->leptons->at(0), corrMet.lv()) * GeV_to_MeV;
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("mT l2 Corr"); {
        *cutflow << HFTname("mTl2Corr");
        *cutflow << [&](Superlink* sl, var_float*) -> double { return PhysicsTools::mT(sl->leptons->at(1), corrMet.lv()) * GeV_to_MeV; };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("mCT"); {
        *cutflow << HFTname("mct");
        *cutflow << [&](Superlink* sl, var_float*) -> double { please return PhysicsTools::mCT(*sl->leptons->at(0), *sl->leptons->at(1)) * GeV_to_MeV; };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("mCT perpendicular"); {
        *cutflow << HFTname("mctPerp");
        *cutflow << [&](Superlink* sl, var_float*) -> double { return PhysicsTools::mCTperp(*sl->leptons->at(0), *sl->leptons->at(1), sl->met->lv()) * GeV_to_MeV; };
        *cutflow << SaveVar();
    }

    // END Setup output trees
    // END Setup output trees
    // END Setup output trees

    // GAP //
    // GAP //
    // GAP //

    // START Setup systematics
    // START Setup systematics
    // START Setup systematics

    // weight variation systematics
    // *cutflow << NewSystematic("shift in electron trigger weights"); {
    //     *cutflow << WeightSystematic(SupersysWeight::ETRIGREWUP, SupersysWeight::ETRIGREWDOWN);
    //     *cutflow << TreeName("ETRIGREW");
    //     *cutflow << SaveSystematic();
    // } // This is now included as an event variation systematic
    // 
    // *cutflow << NewSystematic("shift in muon trigger weights"); {
    //     *cutflow << WeightSystematic(SupersysWeight::MTRIGREWUP, SupersysWeight::MTRIGREWDOWN);
    //     *cutflow << TreeName("MTRIGREW");
    //     *cutflow << SaveSystematic();
    // } // This is now included as an event variation systematic

    *cutflow << NewSystematic("shift in b-tag scale factor"); {
        *cutflow << WeightSystematic(SupersysWeight::BJETUP, SupersysWeight::BJETDOWN);
        *cutflow << TreeName("BJET");
        *cutflow << SaveSystematic();
    }

    *cutflow << NewSystematic("shift in c-tag scale factor"); {
        *cutflow << WeightSystematic(SupersysWeight::CJETUP, SupersysWeight::CJETDOWN);
        *cutflow << TreeName("CJET");
        *cutflow << SaveSystematic();
    }

    *cutflow << NewSystematic("shift in light-tag scale factor"); {
        *cutflow << WeightSystematic(SupersysWeight::BMISTAGUP, SupersysWeight::BMISTAGDOWN);
        *cutflow << TreeName("BMISTAG");
        *cutflow << SaveSystematic();
    }

    *cutflow << NewSystematic("shift in electron efficiency"); {
        *cutflow << WeightSystematic(SupersysWeight::ESFUP, SupersysWeight::ESFDOWN);
        *cutflow << TreeName("ESF");
        *cutflow << SaveSystematic();
    }

    *cutflow << NewSystematic("shift in muon efficiency"); {
        *cutflow << WeightSystematic(SupersysWeight::MEFFUP, SupersysWeight::MEFFDOWN);
        *cutflow << TreeName("MEFF");
        *cutflow << SaveSystematic();
    }

    *cutflow << NewSystematic("shift in pileup"); {
        *cutflow << WeightSystematic(SupersysWeight::PILEUPUP, SupersysWeight::PILEUPDOWN);
        *cutflow << TreeName("PILEUP");
        *cutflow << SaveSystematic();
    }

    *cutflow << NewSystematic("shift in cross section"); {
        *cutflow << WeightSystematic(SupersysWeight::XSUP, SupersysWeight::XSDOWN);
        *cutflow << TreeName("XS");
        *cutflow << SaveSystematic();
    }

    // event variation systematics
    *cutflow << NewSystematic("Electron Scale Zsys + sigma"); {
        *cutflow << EventSystematic(NtSys_EES_Z_UP);
        *cutflow << TreeName("EESZUP");
        *cutflow << SaveSystematic();
    }

    *cutflow << NewSystematic("Electron Scale Zsys - sigma"); {
        *cutflow << EventSystematic(NtSys_EES_Z_DN);
        *cutflow << TreeName("EESZDOWN");
        *cutflow << SaveSystematic();
    }

    *cutflow << NewSystematic("Electron Scale Material + sigma"); {
        *cutflow << EventSystematic(NtSys_EES_MAT_UP);
        *cutflow << TreeName("EESMATUP");
        *cutflow << SaveSystematic();
    }

    *cutflow << NewSystematic("Electron Scale Material - sigma"); {
        *cutflow << EventSystematic(NtSys_EES_MAT_DN);
        *cutflow << TreeName("EESMATDOWN");
        *cutflow << SaveSystematic();
    }

    *cutflow << NewSystematic("Electron Scale Presampler + sigma"); {
        *cutflow << EventSystematic(NtSys_EES_PS_UP);
        *cutflow << TreeName("EESPSUP");
        *cutflow << SaveSystematic();
    }

    *cutflow << NewSystematic("Electron Scale Presampler - sigma"); {
        *cutflow << EventSystematic(NtSys_EES_PS_DN);
        *cutflow << TreeName("EESPSDOWN");
        *cutflow << SaveSystematic();
    }

    *cutflow << NewSystematic("Electron Scale Low Pt + sigma"); {
        *cutflow << EventSystematic(NtSys_EES_LOW_UP);
        *cutflow << TreeName("EESLOWUP");
        *cutflow << SaveSystematic();
    }

    *cutflow << NewSystematic("Electron Scale Low Pt - sigma"); {
        *cutflow << EventSystematic(NtSys_EES_LOW_DN);
        *cutflow << TreeName("EESLOWDOWN");
        *cutflow << SaveSystematic();
    }

    *cutflow << NewSystematic("Electron Resolution + sigma"); {
        *cutflow << EventSystematic(NtSys_EER_UP);
        *cutflow << TreeName("EERUP");
        *cutflow << SaveSystematic();
    }

    *cutflow << NewSystematic("Electron Resolution - sigma"); {
        *cutflow << EventSystematic(NtSys_EER_DN);
        *cutflow << TreeName("EERDOWN");
        *cutflow << SaveSystematic();
    }

    *cutflow << NewSystematic("Muon MS track + sigma"); {
        *cutflow << EventSystematic(NtSys_MS_UP);
        *cutflow << TreeName("MSUP");
        *cutflow << SaveSystematic();
    }

    *cutflow << NewSystematic("Muon MS track - sigma"); {
        *cutflow << EventSystematic(NtSys_MS_DN);
        *cutflow << TreeName("MSDOWN");
        *cutflow << SaveSystematic();
    }

    *cutflow << NewSystematic("Muon ID track + sigma"); {
        *cutflow << EventSystematic(NtSys_ID_UP);
        *cutflow << TreeName("IDUP");
        *cutflow << SaveSystematic();
    }

    *cutflow << NewSystematic("Muon ID track - sigma"); {
        *cutflow << EventSystematic(NtSys_ID_DN);
        *cutflow << TreeName("IDDOWN");
        *cutflow << SaveSystematic();
    }

    *cutflow << NewSystematic("Jet Energy Scale + sigma"); {
        *cutflow << EventSystematic(NtSys_JES_UP);
        *cutflow << TreeName("JESUP");
        *cutflow << SaveSystematic();
    }

    *cutflow << NewSystematic("Jet Energy Scale - sigma"); {
        *cutflow << EventSystematic(NtSys_JES_DN);
        *cutflow << TreeName("JESDOWN");
        *cutflow << SaveSystematic();
    }

    *cutflow << NewSystematic("Jet Energy Resolution (gaussian)"); {
        *cutflow << EventSystematic(NtSys_JER);
        *cutflow << TreeName("JER");
        *cutflow << SaveSystematic();
    }

    *cutflow << NewSystematic("Met scale soft term + sigma"); {
        *cutflow << EventSystematic(NtSys_SCALEST_UP);
        *cutflow << TreeName("SCALESTUP");
        *cutflow << SaveSystematic();
    }

    *cutflow << NewSystematic("Met scale soft term - sigma"); {
        *cutflow << EventSystematic(NtSys_SCALEST_DN);
        *cutflow << TreeName("SCALESTDOWN");
        *cutflow << SaveSystematic();
    }

    *cutflow << NewSystematic("Met resolution soft term + sigma"); {
        *cutflow << EventSystematic(NtSys_RESOST);
        *cutflow << TreeName("RESOST");
        *cutflow << SaveSystematic();
    }

    *cutflow << NewSystematic("Trigger Scale factor + error for el"); {
        *cutflow << EventSystematic(NtSys_TRIGSF_EL_UP);
        *cutflow << TreeName("TRIGSFELUP");
        *cutflow << SaveSystematic();
    }

    *cutflow << NewSystematic("Trigger Scale factor - error for el"); {
        *cutflow << EventSystematic(NtSys_TRIGSF_EL_DN);
        *cutflow << TreeName("TRIGSFELDN");
        *cutflow << SaveSystematic();
    }

    *cutflow << NewSystematic("Trigger Scale factor + error for mu"); {
        *cutflow << EventSystematic(NtSys_TRIGSF_MU_UP);
        *cutflow << TreeName("TRIGSFMUUP");
        *cutflow << SaveSystematic();
    }

    *cutflow << NewSystematic("Trigger Scale factor - error for mu"); {
        *cutflow << EventSystematic(NtSys_TRIGSF_MU_DN);
        *cutflow << TreeName("TRIGSFMUDN");
        *cutflow << SaveSystematic();
    }

    *cutflow << NewSystematic("Tau energy scale + sigma"); {
        *cutflow << EventSystematic(NtSys_TES_UP);
        *cutflow << TreeName("TESUP");
        *cutflow << SaveSystematic();
    }

    *cutflow << NewSystematic("Tau energy scale - sigma"); {
        *cutflow << EventSystematic(NtSys_TES_DN);
        *cutflow << TreeName("TESDOWN");
        *cutflow << SaveSystematic();
    }

    *cutflow << NewSystematic("Jet JVF cut + sigma"); {
        *cutflow << EventSystematic(NtSys_JVF_UP);
        *cutflow << TreeName("JVFUP");
        *cutflow << SaveSystematic();
    }

    *cutflow << NewSystematic("Jet JVF cut - sigma"); {
        *cutflow << EventSystematic(NtSys_JVF_DN);
        *cutflow << TreeName("JVFDOWN");
        *cutflow << SaveSystematic();
    }

    // END Setup systematics
    // END Setup systematics
    // END Setup systematics

    // Initialize the cutflow and start the event loop.
    chain->Process(cutflow, sample_.c_str(), num_events_, n_skip_);

    delete cutflow;
    delete chain;

    cout << "Done." << endl;
    exit(0);
}

void print_usage(const char *exeName) // Courtesy of Davide
{
    cout << "Usage: " << endl
        << exeName << " options" << endl
        << "\t -h [--help]                   : print this message \n"
        << "\t -v [--verbose]                : toggle verbose \n"
        << "\t -d [--debug]                  : toggle debug \n"
        << "\t -i [--input]   <file.txt>     : input root file (or filelist or dir) \n"
        << "\t -n [--num-events] <int>       : number of events to be processed \n"
        << "\t -s [--sample_]  <samplename>   : sample_ name \n"
        << "\t -e [--event-list] <file.root> : file where the eventlist is cached \n"
        << "Example: \n"
        << exeName << "\n"
        << " -i filelist/Sherpa_CT10_lllnu_WZ_MassiveCB.txt \n"
        << " -s Sherpa_CT10_lllnu_WZ_MassiveCB \n"
        << " -e out/selection/eventlist_files \n"
        << endl;
}

void read_options(int argc, char* argv[], TChain* chain, int& n_skip_, int& num_events_, string& sample_,
    SuperflowRunMode& run_mode_, SusyNtSys& nt_sys) // Courtesy of Davide
{
    bool debug = false;
    bool nominal_ = false;
    bool nominal_and_weight_syst_ = false;
    bool all_syst_ = false;
    bool single_event_syst_ = false;

    string systematic_ = "undefined";

    string input;
    string output;
    string eventlist;
    string tuple_out;

    /** Read inputs to program */
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "/n") == 0)
            num_events_ = atoi(argv[++i]);
        else if (strcmp(argv[i], "/c") == 0)
            nominal_ = true;
        else if (strcmp(argv[i], "/w") == 0)
            nominal_and_weight_syst_ = true;
        else if (strcmp(argv[i], "/e") == 0) {
            single_event_syst_ = true;
        }
        else if (strcmp(argv[i], "/a") == 0) {
            all_syst_ = true;
        }
        else if (strcmp(argv[i], "/i") == 0) {
            input = argv[++i];
        }
        else if (strcmp(argv[i], "/s") == 0) {
            systematic_ = argv[++i];
        }
        else {
            cout << "Analysis    Error (fatal): Bad arguments." << endl;
            exit(1);
        }
    }

    // TChain* chain = new TChain("susyNt");
    bool inputIsFile = susy::utils::endswith(input, ".root");
    bool inputIsList = susy::utils::endswith(input, ".txt");
    bool inputIsDir = susy::utils::endswith(input, "/");
    bool validInput(inputIsFile || inputIsList || inputIsDir);
    if (!validInput) {
        cout << "Analysis    invalid input '" << input << "'" << endl;
        print_usage(argv[0]);
        exit(1);
    }
    if (inputIsFile) {
        ChainHelper::addFile(chain, input);
        cout << "Analysis    file: " << input << endl;
        cout << "Analysis    file: " << input << endl;
        cout << "Analysis    file: " << input << endl;
        sample_ = input;
    }
    if (inputIsList) {
        ChainHelper::addFileList(chain, input);
        cout << "Analysis    list: " << input << endl;
        cout << "Analysis    list: " << input << endl;
        cout << "Analysis    list: " << input << endl;
        ifstream infile(input.c_str());
        if (infile.good()) {
            string sLine;
            getline(infile, sLine);
            sample_ = sLine;
        }
        else {
            sample_ = input;
        }
        infile.close();
    }
    if (inputIsDir) {
        ChainHelper::addFileDir(chain, input);
        cout << "Analysis    dir: " << input << endl;
        cout << "Analysis    dir: " << input << endl;
        cout << "Analysis    dir: " << input << endl;
        sample_ = input;
    }
    Long64_t tot_num_events = chain->GetEntries();
    num_events_ = (num_events_ < 0 ? tot_num_events : num_events_);
    if (debug) chain->ls();

    if (nominal_) {
        run_mode_ = SuperflowRunMode::nominal;
        cout << "Analysis    run mode: SuperflowRunMode::nominal" << endl;
    }
    if (nominal_and_weight_syst_) {
        run_mode_ = SuperflowRunMode::nominal_and_weight_syst;
        cout << "Analysis    run mode: SuperflowRunMode::nominal_and_weight_syst" << endl;
    }
    if (single_event_syst_) {
        run_mode_ = SuperflowRunMode::single_event_syst;
        cout << "Analysis    run mode: SuperflowRunMode::single_event_syst" << endl;
    }

    if (all_syst_) {
        run_mode_ = SuperflowRunMode::all_syst;
        cout << "Analysis    run mode: SuperflowRunMode::all_syst" << endl;
    }

    map <string, SusyNtSys> event_syst_map;
    event_syst_map["EESZUP"] = NtSys_EES_Z_UP;
    event_syst_map["EESZDOWN"] = NtSys_EES_Z_DN;
    event_syst_map["EESMATUP"] = NtSys_EES_MAT_UP;
    event_syst_map["EESMATDOWN"] = NtSys_EES_MAT_DN;
    event_syst_map["EESPSUP"] = NtSys_EES_PS_UP;
    event_syst_map["EESPSDOWN"] = NtSys_EES_PS_DN;
    event_syst_map["EESLOWUP"] = NtSys_EES_LOW_UP;
    event_syst_map["EESLOWDOWN"] = NtSys_EES_LOW_DN;
    event_syst_map["EERUP"] = NtSys_EER_UP;
    event_syst_map["EERDOWN"] = NtSys_EER_DN;
    event_syst_map["MSUP"] = NtSys_MS_UP;
    event_syst_map["MSDOWN"] = NtSys_MS_DN;
    event_syst_map["IDUP"] = NtSys_ID_UP;
    event_syst_map["IDDOWN"] = NtSys_ID_DN;
    event_syst_map["JESUP"] = NtSys_JES_UP;
    event_syst_map["JESDOWN"] = NtSys_JES_DN;
    event_syst_map["JER"] = NtSys_JER;
    event_syst_map["SCALESTUP"] = NtSys_SCALEST_UP;
    event_syst_map["SCALESTDOWN"] = NtSys_SCALEST_DN;
    event_syst_map["RESOST"] = NtSys_RESOST;
    event_syst_map["TRIGSFELUP"] = NtSys_TRIGSF_EL_UP;
    event_syst_map["TRIGSFELDN"] = NtSys_TRIGSF_EL_DN;
    event_syst_map["TRIGSFMUUP"] = NtSys_TRIGSF_MU_UP;
    event_syst_map["TRIGSFMUDN"] = NtSys_TRIGSF_MU_DN;
    event_syst_map["TESUP"] = NtSys_TES_UP;
    event_syst_map["TESDOWN"] = NtSys_TES_DN;
    event_syst_map["JVFUP"] = NtSys_JVF_UP;
    event_syst_map["JVFDOWN"] = NtSys_JVF_DN;

    if (single_event_syst_) {
        if (event_syst_map.count(systematic_) == 1) {
            nt_sys = event_syst_map[systematic_];
        }
        else {
            cout << "Analysis" << "    ERROR (fatal): Event systematic option -t " << systematic_ << " / not found." << endl;
            exit(1);
        }
    }
}


// List of event systematics (for scripting)
// 
// "EESZUP",
// "EESZDOWN",
// "EESMATUP",
// "EESMATDOWN",
// "EESPSUP",
// "EESPSDOWN",
// "EESLOWUP",
// "EESLOWDOWN",
// "EERUP",
// "EERDOWN",
// "MSUP",
// "MSDOWN",
// "IDUP",
// "IDDOWN",
// "JESUP",
// "JESDOWN",
// "JER",
// "SCALESTUP",
// "SCALESTDOWN",
// "RESOST",
// "TRIGSFELUP",
// "TRIGSFELDN",
// "TRIGSFMUUP",
// "TRIGSFMUDN",
// "TESUP",
// "TESDOWN",
// "JVFUP",
// "JVFDOWN",
// 


// A SELECTION OF CUTS
// A SELECTION OF CUTS
// A SELECTION OF CUTS
//
// *cutflow << CutName("exactly two base leptons") << [](Superlink* sl) -> bool {
//     return sl->baseLeptons->size() == 2;
// };
//
// *cutflow << CutName("m_ll > 20 GeV") << [](Superlink* sl) -> bool {
//     double m_ll = (*sl->baseLeptons->at(0) + *sl->baseLeptons->at(1)).M();
//     return m_ll > 20.0;
// };
//
// *cutflow << CutName("is MM") << [](Superlink* sl) -> bool {
//     return sl->baseLeptons->at(0)->isMu() && sl->baseLeptons->at(1)->isMu();
// }; // debug only !!!
//
// *cutflow << CutName("exactly two signal leptons") << [](Superlink* sl) -> bool {
//     return sl->leptons->size() == 2;
// };
//
// A SELECTION OF CUTS
// A SELECTION OF CUTS
// A SELECTION OF CUTS


//rem// Cut* is_e_mu = new IsEMu(); // initialize cut
//rem// Cut* is_two_lepton = new Is2Lepton();
//rem// Cut* is_e_e_and_os = new IsEEOS();
//rem// *cutflow << (new IsEEOS()); // push cut to cutflow
//rem// class Is2Lepton : public Cut {
//rem// public:
//rem//     Is2Lepton()
//rem//     {
//rem//         name = "base is 2-lepton";
//rem//     }
//rem//     bool operator() (Superlink* sl) // return true to pass the cut
//rem//     {
//rem//         return sl->baseLeptons->size() == 2; // exactly two base leptons
//rem//     }
//rem// };
//rem// 
//rem// class IsEMu : public Cut {
//rem// public:
//rem//     IsEMu()
//rem//     {
//rem//         name = "base is e + mu";
//rem//     }
//rem//     bool operator() (Superlink* sl) // return true to pass the cut
//rem//     {
//rem//         if (sl->baseLeptons->size() == 2) { // exactly two base leptons
//rem//             return sl->baseLeptons->at(0)->isEle() ^ sl->baseLeptons->at(1)->isEle(); // e + mu
//rem//         }
//rem//         else {
//rem//             return false;
//rem//         }
//rem//     }
//rem// };
//rem// 
//rem// class IsEEOS : public Cut {
//rem// public:
//rem//     IsEEOS()
//rem//     {
//rem//         name = "base is e + e + OS";
//rem//     }
//rem//     bool operator() (Superlink* sl) // return true to pass the cut
//rem//     {
//rem//         if (sl->baseLeptons->size() == 2) { // exactly two base leptons
//rem//             return (sl->baseLeptons->at(0)->q * sl->baseLeptons->at(1)->q < 0)
//rem//                 && (sl->baseLeptons->at(0)->isEle() && sl->baseLeptons->at(1)->isEle()); // e + mu
//rem//         }
//rem//         else {
//rem//             return false;
//rem//         }
//rem//     }
//rem// };
