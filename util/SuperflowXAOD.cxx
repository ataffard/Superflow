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

#include "TMath.h"
#include "TVector2.h"

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
    SusyNtSys nt_sys_ = Susy::NtSys::NOM;

    TChain* chain = new TChain("susyNt");
    chain->SetDirectory(0);

    read_options(argc, argv, chain, n_skip_, num_events_, sample_, run_mode, nt_sys_); // defined below
    // END read-in

    Superflow* cutflow = new Superflow(); // initialize the cutflow
    cutflow->setAnaType(Ana_2Lep, true); // Ana_2Lep Ana_2LepWH 
    cutflow->setSampleName(sample_);
    cutflow->setRunMode(run_mode);
    cutflow->setChain(chain);
    cutflow->setCountWeights(true);

    cout << "Analysis    Total Entries: " << chain->GetEntries() << endl;

    if (run_mode == SuperflowRunMode::single_event_syst) cutflow->setSingleEventSyst(nt_sys_);

    // START Setup cuts
    // START Setup cuts
    // START Setup cuts

    *cutflow << CutName("read in") << [](Superlink* sl) -> bool { return true; };

    *cutflow << CutName("exactly two base leptons") << [](Superlink* sl) -> bool {
        return sl->baseLeptons->size() == 2;
    };

    *cutflow << CutName("m_ll > 20 GeV") << [](Superlink* sl) -> bool {
        return (*sl->baseLeptons->at(0) + *sl->baseLeptons->at(1)).M() > 20.0;
    };

    *cutflow << CutName("leading lepton eta < 2.4") << [](Superlink* sl) -> bool {
        return abs(sl->baseLeptons->at(0)->Eta()) < 2.4;
    };

    *cutflow << CutName("sub-leading lepton eta < 2.4") << [](Superlink* sl) -> bool {
        return abs(sl->baseLeptons->at(1)->Eta()) < 2.4;
    };

    *cutflow << CutName("tau veto") << [](Superlink* sl) -> bool {
        return sl->taus->size() == 0;
    };

    *cutflow << CutName("leading lepton Pt > 30 GeV") << [](Superlink* sl) -> bool {
        return sl->baseLeptons->at(0)->Pt() > 30.0; // THIS WAS RUN AT 30 GeV IN PREVIOUS PRODUCTIONS
    };

    *cutflow << CutName("subleading lepton Pt > 20 GeV") << [](Superlink* sl) -> bool {
        return sl->baseLeptons->at(1)->Pt() > 20.0;
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
        *cutflow << [](Superlink* sl, var_double*) -> double { return sl->weights->product(); };
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

    *cutflow << NewVar("is Monte Carlo"); {
        *cutflow << HFTname("isMC");
        *cutflow << [](Superlink* sl, var_bool*) -> bool { return sl->nt->evt()->isMC ? true : false; };
        *cutflow << SaveVar();
    }

    // LEPTONS
    // LEPTONS
    // LEPTONS

    *cutflow << NewVar("is e + e"); {
        *cutflow << HFTname("isElEl");
        *cutflow << [](Superlink* sl, var_bool*) -> bool { return sl->baseLeptons->at(0)->isEle() && sl->baseLeptons->at(1)->isEle(); };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("is e + mu"); {
        *cutflow << HFTname("isElMu");
        *cutflow << [](Superlink* sl, var_bool*) -> bool { return sl->baseLeptons->at(0)->isEle() ^ sl->baseLeptons->at(1)->isEle(); };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("is mu + mu"); {
        *cutflow << HFTname("isMuMu");
        *cutflow << [](Superlink* sl, var_bool*) -> bool { return sl->baseLeptons->at(0)->isMu() && sl->baseLeptons->at(1)->isMu(); };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("is opposite-sign"); {
        *cutflow << HFTname("isOS");
        *cutflow << [](Superlink* sl, var_bool*) -> bool { return sl->baseLeptons->at(0)->q * sl->baseLeptons->at(1)->q < 0; };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("is Mu (lead) + E (sub)"); {
        *cutflow << HFTname("isME");
        *cutflow << [](Superlink* sl, var_bool*) -> bool { return sl->baseLeptons->at(0)->isMu() && sl->baseLeptons->at(1)->isEle(); };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("is E (lead) + Mu (sub)"); {
        *cutflow << HFTname("isEM");
        *cutflow << [](Superlink* sl, var_bool*) -> bool { return sl->baseLeptons->at(0)->isEle() && sl->baseLeptons->at(1)->isMu(); };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("lepton-1 Pt"); {
        *cutflow << HFTname("lept1Pt");
        *cutflow << [](Superlink* sl, var_float*) -> double { return sl->baseLeptons->at(0)->Pt(); };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("lepton-1 Eta"); {
        *cutflow << HFTname("lept1Eta");
        *cutflow << [](Superlink* sl, var_float*) -> double { return sl->baseLeptons->at(0)->Eta(); };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("lepton-1 Phi"); {
        *cutflow << HFTname("lept1Phi");
        *cutflow << [](Superlink* sl, var_float*) -> double { return sl->baseLeptons->at(0)->Phi(); };
        *cutflow << SaveVar();
    }

    // *cutflow << NewVar("lepton-1 Energy"); {
    //     *cutflow << HFTname("lept1E");
    //     *cutflow << [](Superlink* sl, var_float*) -> double { return sl->baseLeptons->at(0)->E(); };
    //     *cutflow << SaveVar();
    // }

    *cutflow << NewVar("lepton-1 charge"); {
        *cutflow << HFTname("lept1q");
        *cutflow << [](Superlink* sl, var_int*) -> int { return sl->baseLeptons->at(0)->q; };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("lepton-1 flavor"); {
        *cutflow << HFTname("lept1Flav");
        *cutflow << [](Superlink* sl, var_int*) -> int { return sl->baseLeptons->at(0)->isEle() ? 0 : 1; };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("lepton-1 isEle"); {
        *cutflow << HFTname("lept1isEle");
        *cutflow << [](Superlink* sl, var_bool*) -> bool { return sl->baseLeptons->at(0)->isEle(); };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("lepton-1 isMu"); {
        *cutflow << HFTname("lept1isMu");
        *cutflow << [](Superlink* sl, var_bool*) -> bool { return sl->baseLeptons->at(0)->isMu(); };
        *cutflow << SaveVar();
    }

    // START isolation variables

    *cutflow << NewVar("lepton-1 ptcone20"); {
        *cutflow << HFTname("lept1ptcone20");
        *cutflow << [](Superlink* sl, var_float*) -> double { return sl->baseLeptons->at(0)->ptcone20; };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("lepton-1 ptcone20 / Pt"); {
        *cutflow << HFTname("lept1ptcone20_Pt");
        *cutflow << [](Superlink* sl, var_float*) -> double { return abs(sl->baseLeptons->at(0)->ptcone20 / sl->baseLeptons->at(0)->Pt()); };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("lepton-1 ptcone30"); {
        *cutflow << HFTname("lept1ptcone30");
        *cutflow << [](Superlink* sl, var_float*) -> double { return sl->baseLeptons->at(0)->ptcone30; };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("lepton-1 ptcone30 / Pt"); {
        *cutflow << HFTname("lept1ptcone30_Pt");
        *cutflow << [](Superlink* sl, var_float*) -> double { return abs(sl->baseLeptons->at(0)->ptcone30 / sl->baseLeptons->at(0)->Pt()); };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("lepton-1 etcone20"); {
        *cutflow << HFTname("lept1etcone20");
        *cutflow << [](Superlink* sl, var_float*) -> double { return sl->baseLeptons->at(0)->etcone20; };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("lepton-1 etcone20 / Pt"); {
        *cutflow << HFTname("lept1etcone20_Pt");
        *cutflow << [](Superlink* sl, var_float*) -> double { return abs(sl->baseLeptons->at(0)->etcone20 / sl->baseLeptons->at(0)->Pt()); };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("lepton-1 z0"); {
        *cutflow << HFTname("lept1z0");
        *cutflow << [](Superlink* sl, var_float*) -> double { return sl->baseLeptons->at(0)->z0; };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("lepton-1 z0 error"); {
        *cutflow << HFTname("lept1z0_err");
        *cutflow << [](Superlink* sl, var_float*) -> double { return sl->baseLeptons->at(0)->errZ0; };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("lepton-1 z0 * sin(theta)"); {
        *cutflow << HFTname("lept1z0sin_theta");
        *cutflow << [](Superlink* sl, var_float*) -> double { return sl->baseLeptons->at(0)->z0SinTheta(); };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("lepton-1 d0"); {
        *cutflow << HFTname("lept1d0");
        *cutflow << [](Superlink* sl, var_float*) -> double { return sl->baseLeptons->at(0)->d0; };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("lepton-1 d0 error"); {
        *cutflow << HFTname("lept1d0_err");
        *cutflow << [](Superlink* sl, var_float*) -> double { return sl->baseLeptons->at(0)->errD0; };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("lepton-1 d0 significance"); {
        *cutflow << HFTname("lept1d0_sig");
        *cutflow << [](Superlink* sl, var_float*) -> double { return abs(sl->baseLeptons->at(0)->d0 / sl->baseLeptons->at(0)->errD0); };
        *cutflow << SaveVar();
    }

    // END isolation variables

    *cutflow << NewVar("lepton-2 Pt"); {
        *cutflow << HFTname("lept2Pt");
        *cutflow << [](Superlink* sl, var_float*) -> double { return sl->baseLeptons->at(1)->Pt(); };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("lepton-2 Eta"); {
        *cutflow << HFTname("lept2Eta");
        *cutflow << [](Superlink* sl, var_float*) -> double { return sl->baseLeptons->at(1)->Eta(); };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("lepton-2 Phi"); {
        *cutflow << HFTname("lept2Phi");
        *cutflow << [](Superlink* sl, var_float*) -> double { return sl->baseLeptons->at(1)->Phi(); };
        *cutflow << SaveVar();
    }

    // *cutflow << NewVar("lepton-2 Energy"); {
    //     *cutflow << HFTname("lept2E");
    //     *cutflow << [](Superlink* sl, var_float*) -> double { return sl->baseLeptons->at(1)->E(); };
    //     *cutflow << SaveVar();
    // }

    *cutflow << NewVar("lepton-2 charge"); {
        *cutflow << HFTname("lept2q");
        *cutflow << [](Superlink* sl, var_int*) -> int { return sl->baseLeptons->at(1)->q; };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("lepton-2 isEle"); {
        *cutflow << HFTname("lept2isEle");
        *cutflow << [](Superlink* sl, var_bool*) -> bool { return sl->baseLeptons->at(1)->isEle(); };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("lepton-2 isMu"); {
        *cutflow << HFTname("lept2isMu");
        *cutflow << [](Superlink* sl, var_bool*) -> bool { return sl->baseLeptons->at(1)->isMu(); };
        *cutflow << SaveVar();
    }

    // START isolation variables

    *cutflow << NewVar("lepton-2 ptcone20"); {
        *cutflow << HFTname("lept2ptcone20");
        *cutflow << [](Superlink* sl, var_float*) -> double { return sl->baseLeptons->at(1)->ptcone20; };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("lepton-2 ptcone20 / Pt"); {
        *cutflow << HFTname("lept2ptcone20_Pt");
        *cutflow << [](Superlink* sl, var_float*) -> double { return abs(sl->baseLeptons->at(1)->ptcone20 / sl->baseLeptons->at(1)->Pt()); };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("lepton-2 ptcone30"); {
        *cutflow << HFTname("lept2ptcone30");
        *cutflow << [](Superlink* sl, var_float*) -> double { return sl->baseLeptons->at(1)->ptcone30; };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("lepton-2 ptcone30 / Pt"); {
        *cutflow << HFTname("lept2ptcone30_Pt");
        *cutflow << [](Superlink* sl, var_float*) -> double { return abs(sl->baseLeptons->at(1)->ptcone30 / sl->baseLeptons->at(1)->Pt()); };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("lepton-2 etcone20"); {
        *cutflow << HFTname("lept2etcone20");
        *cutflow << [](Superlink* sl, var_float*) -> double { return sl->baseLeptons->at(1)->etcone20; };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("lepton-2 etcone20 / Pt"); {
        *cutflow << HFTname("lept2etcone20_Pt");
        *cutflow << [](Superlink* sl, var_float*) -> double { return abs(sl->baseLeptons->at(1)->etcone20 / sl->baseLeptons->at(1)->Pt()); };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("lepton-2 z0"); {
        *cutflow << HFTname("lept2z0");
        *cutflow << [](Superlink* sl, var_float*) -> double { return sl->baseLeptons->at(1)->z0; };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("lepton-2 z0 error"); {
        *cutflow << HFTname("lept2z0_err");
        *cutflow << [](Superlink* sl, var_float*) -> double { return sl->baseLeptons->at(1)->errZ0; };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("lepton-2 z0 * sin(theta)"); {
        *cutflow << HFTname("lept2z0sin_theta");
        *cutflow << [](Superlink* sl, var_float*) -> double { return sl->baseLeptons->at(1)->z0SinTheta(); };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("lepton-2 d0"); {
        *cutflow << HFTname("lept2d0");
        *cutflow << [](Superlink* sl, var_float*) -> double { return sl->baseLeptons->at(1)->d0; };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("lepton-2 d0 error"); {
        *cutflow << HFTname("lept2d0_err");
        *cutflow << [](Superlink* sl, var_float*) -> double { return sl->baseLeptons->at(1)->errD0; };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("lepton-2 d0 significance"); {
        *cutflow << HFTname("lept2d0_sig");
        *cutflow << [](Superlink* sl, var_float*) -> double { return abs(sl->baseLeptons->at(1)->d0 / sl->baseLeptons->at(1)->errD0); };
        *cutflow << SaveVar();
    }

    // END isolation variables

    // JETS
    // JETS
    // JETS

    JetVector central_light_jets; // local variable!

    *cutflow << NewVar("number of central light jets"); {
        *cutflow << HFTname("nCentralLightJets");
        *cutflow << [&](Superlink* sl, var_int*) -> int {
            for (int i = 0; i < sl->jets->size(); i++) {
                if (sl->tools->isCentralLightJet(sl->jets->at(i), sl->jvfTool, sl->nt_sys, sl->anaType)) {
                    central_light_jets.push_back(sl->jets->at(i));
                }
            }
            return central_light_jets.size();
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("number of central b jets"); {
        *cutflow << HFTname("nCentralBJets");
        *cutflow << [](Superlink* sl, var_int*) -> int { return sl->tools->numberOfCBJets(*sl->jets); };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("number of forward jets"); {
        *cutflow << HFTname("nForwardJets");
        *cutflow << [](Superlink* sl, var_int*) -> int { return sl->tools->numberOfFJets(*sl->jets); };
        *cutflow << SaveVar();
    }


    *cutflow << NewVar("jet-1 Pt"); {
        *cutflow << HFTname("jet1Pt");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return central_light_jets.size() >= 1 ? central_light_jets[0]->Pt() : 0.0;
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("jet-1 Eta"); {
        *cutflow << HFTname("jet1Eta");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return central_light_jets.size() >= 1 ? central_light_jets[0]->Eta() : 0.0;
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("jet-1 Phi"); {
        *cutflow << HFTname("jet1Phi");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return central_light_jets.size() >= 1 ? central_light_jets[0]->Phi() : 0.0;
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("jet-2 Pt"); {
        *cutflow << HFTname("jet2Pt");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return central_light_jets.size() >= 2 ? central_light_jets[1]->Pt() : 0.0;
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("jet-2 Eta"); {
        *cutflow << HFTname("jet2Eta");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return central_light_jets.size() >= 2 ? central_light_jets[1]->Eta() : 0.0;
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("jet-2 Phi"); {
        *cutflow << HFTname("jet2Phi");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return central_light_jets.size() >= 2 ? central_light_jets[1]->Phi() : 0.0;
        };
        *cutflow << SaveVar();
    }

    *cutflow << [&](Superlink* sl, var_void*) { central_light_jets.clear(); };

    // MET
    // MET
    // MET

    *cutflow << NewVar("transverse missing energy (Et)"); {
        *cutflow << HFTname("met");
        *cutflow << [](Superlink* sl, var_float*) -> double { return sl->met->Et; };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("transverse missing energy (Phi)"); {
        *cutflow << HFTname("metPhi");
        *cutflow << [](Superlink* sl, var_float*) -> double { return sl->met->phi; };
        *cutflow << SaveVar();
    }

    // VARS
    // VARS
    // VARS

    TLorentzVector local_ll; // local variable!

    *cutflow << NewVar("mass of di-lepton system, M_ll"); {
        *cutflow << HFTname("ll_M");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            local_ll = (*sl->baseLeptons->at(0) + *sl->baseLeptons->at(1));
            return local_ll.M();
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("Pt of di-lepton system, Pt_ll"); {
        *cutflow << HFTname("ll_Pt");
        *cutflow << [&](Superlink* sl, var_float*) -> double { return local_ll.Pt(); };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("delta Eta of di-lepton system"); {
        *cutflow << HFTname("ll_deltaEta");
        *cutflow << [](Superlink* sl, var_float*) -> double { return abs(sl->baseLeptons->at(0)->Eta() - sl->baseLeptons->at(1)->Eta()); };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("delta Phi of di-lepton system"); {
        *cutflow << HFTname("ll_deltaPhi");
        *cutflow << [](Superlink* sl, var_float*) -> double { return acos(cos(sl->baseLeptons->at(0)->Phi() - sl->baseLeptons->at(1)->Phi())); };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("stransverse mass"); {
        *cutflow << HFTname("MT2");
        *cutflow << [](Superlink* sl, var_float*) -> double {

            mt2_bisect::mt2 mt2_event;
            double *pa, *pb, *pmiss;
            pa = new double[3]; pa[0] = sl->baseLeptons->at(0)->M(); pa[1] = sl->baseLeptons->at(0)->Px(), pa[2] = sl->baseLeptons->at(0)->Py();
            pb = new double[3]; pb[0] = sl->baseLeptons->at(1)->M(); pb[1] = sl->baseLeptons->at(1)->Px(), pb[2] = sl->baseLeptons->at(1)->Py();
            pmiss = new double[3]; pmiss[0] = 0.0; pmiss[1] = sl->met->Et * cos(sl->met->phi); pmiss[2] = sl->met->Et * sin(sl->met->phi);

            mt2_event.set_momenta(pa, pb, pmiss);
            mt2_event.set_mn(0.0); // LSP mass = 0 is Generic

            double mt2_ = mt2_event.get_mt2();
            // SUPRESSED messages "Deltasq_high not found at event 0" (~1 per 10000 events)

            delete[] pa;
            delete[] pb;
            delete[] pmiss;

            return mt2_;
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("Ht (m_Eff)"); {
        *cutflow << HFTname("HTcorr");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            double ht = 0.0;

            ht += sl->baseLeptons->at(0)->Pt() + sl->baseLeptons->at(1)->Pt();
            ht += sl->met->Et;
            for (int i = 0; i < sl->jets->size(); i++) {
                if (sl->jets->at(i)->Pt() > 20.0) {
                    ht += sl->jets->at(i)->Pt();
                }
            }
            return ht;
        };
        *cutflow << SaveVar();
    }

    // Initialize the cutflow and start the event loop.
    chain->Process(cutflow, sample_.c_str(), num_events_, n_skip_);

    delete cutflow;
    delete chain;

    cout << "Done." << endl;
    exit(0);
}

void read_options(int argc, char* argv[], TChain* chain, int& n_skip_, int& num_events_, string& sample_,
                  SuperflowRunMode& run_mode_, SusyNtSys& nt_sys)
{
    bool nominal_ = false;
    bool nominal_and_weight_syst_ = false;
    bool all_syst_ = false;
    bool single_event_syst_ = false;

    string systematic_ = "undefined";

    string input;

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

    bool inputIsFile = susy::utils::endswith(input, ".root");
    bool inputIsList = susy::utils::endswith(input, ".txt");
    bool inputIsDir = susy::utils::endswith(input, "/");
    bool validInput(inputIsFile || inputIsList || inputIsDir);
    if (!validInput) {
        cout << "Analysis    invalid input '" << input << "'" << endl;
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
    // if (debug) chain->ls();

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
    event_syst_map["EESZUP"] = Susy::NtSys::EES_Z_UP;
    event_syst_map["EESZDOWN"] = Susy::NtSys::EES_Z_DN;
    event_syst_map["EESMATUP"] = Susy::NtSys::EES_MAT_UP;
    event_syst_map["EESMATDOWN"] = Susy::NtSys::EES_MAT_DN;
    event_syst_map["EESPSUP"] = Susy::NtSys::EES_PS_UP;
    event_syst_map["EESPSDOWN"] = Susy::NtSys::EES_PS_DN;
    event_syst_map["EESLOWUP"] = Susy::NtSys::EES_LOW_UP;
    event_syst_map["EESLOWDOWN"] = Susy::NtSys::EES_LOW_DN;
    event_syst_map["EERUP"] = Susy::NtSys::EER_UP;
    event_syst_map["EERDOWN"] = Susy::NtSys::EER_DN;
    event_syst_map["MSUP"] = Susy::NtSys::MS_UP;
    event_syst_map["MSDOWN"] = Susy::NtSys::MS_DN;
    event_syst_map["IDUP"] = Susy::NtSys::ID_UP;
    event_syst_map["IDDOWN"] = Susy::NtSys::ID_DN;
    event_syst_map["JESUP"] = Susy::NtSys::JES_UP;
    event_syst_map["JESDOWN"] = Susy::NtSys::JES_DN;
    event_syst_map["JER"] = Susy::NtSys::JER;
    event_syst_map["SCALESTUP"] = Susy::NtSys::SCALEST_UP;
    event_syst_map["SCALESTDOWN"] = Susy::NtSys::SCALEST_DN;
    event_syst_map["RESOST"] = Susy::NtSys::RESOST;
    // event_syst_map["TRIGSFELUP"] = Susy::NtSys::TRIGSF_EL_UP; // doesn't exist
    // event_syst_map["TRIGSFELDN"] = Susy::NtSys::TRIGSF_EL_DN;
    // event_syst_map["TRIGSFMUUP"] = Susy::NtSys::TRIGSF_MU_UP;
    // event_syst_map["TRIGSFMUDN"] = Susy::NtSys::TRIGSF_MU_DN;
    event_syst_map["TESUP"] = Susy::NtSys::TES_UP;
    event_syst_map["TESDOWN"] = Susy::NtSys::TES_DN;
    event_syst_map["JVFUP"] = Susy::NtSys::JVF_UP;
    event_syst_map["JVFDOWN"] = Susy::NtSys::JVF_DN;

    if (single_event_syst_) {
        if (event_syst_map.count(systematic_) == 1) {
            nt_sys = event_syst_map[systematic_];
        }
        else {
            cout << "Analysis" << "    ERROR (fatal): Event systematic option /s " << systematic_ << " -> not found." << endl;
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


// A SELECTION OF CUTS and Variables
// A SELECTION OF CUTS and Variables
// A SELECTION OF CUTS and Variables
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
// *cutflow << NewVar("mCT"); {
//     *cutflow << HFTname("mct");
//     *cutflow << [&](Superlink* sl, var_float*) -> double {  return PhysicsTools::mCT(*sl->leptons->at(0), *sl->leptons->at(1)); };
//     *cutflow << SaveVar();
// }
// 
// *cutflow << NewVar("mCT perpendicular"); {
//     *cutflow << HFTname("mctPerp");
//     *cutflow << [&](Superlink* sl, var_float*) -> double { return PhysicsTools::mCTperp(*sl->leptons->at(0), *sl->leptons->at(1), sl->met->lv()); };
//     *cutflow << SaveVar();
// }
//
// *cutflow << NewSystematic("Trigger Scale factor + error for el"); {
//     *cutflow << EventSystematic(TRIGSF_EL_UP);
//     *cutflow << TreeName("TRIGSFELUP");
//     *cutflow << SaveSystematic();
// } // These are correctly included as event variation systematics
// 
// *cutflow << NewSystematic("Trigger Scale factor - error for el"); {
//     *cutflow << EventSystematic(TRIGSF_EL_DN);
//     *cutflow << TreeName("TRIGSFELDN");
//     *cutflow << SaveSystematic();
// } // These are correctly included as event variation systematics
// 
// *cutflow << NewSystematic("Trigger Scale factor + error for mu"); {
//     *cutflow << EventSystematic(TRIGSF_MU_UP);
//     *cutflow << TreeName("TRIGSFMUUP");
//     *cutflow << SaveSystematic();
// } // These are correctly included as event variation systematics
// 
// *cutflow << NewSystematic("Trigger Scale factor - error for mu"); {
//     *cutflow << EventSystematic(TRIGSF_MU_DN);
//     *cutflow << TreeName("TRIGSFMUDN");
//     *cutflow << SaveSystematic();
// } // These are correctly included as event variation systematics
//
// A SELECTION OF CUTS and Variables
// A SELECTION OF CUTS and Variables
// A SELECTION OF CUTS and Variables


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
