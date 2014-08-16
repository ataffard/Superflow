#pragma once

#include <vector>

#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TEntryList.h"

#include "TVector2.h"

#include "SusyNtuple/SusyDefs.h"
#include "SusyNtuple/SusyNtAna.h"
#include "SusyNtuple/SusyNtTools.h"

#include "SusyNtuple/MCWeighter.h"
#include "SusyNtuple/DilTrigLogic.h"

#include "Superflow/DataDefinitions.h"

#include "Superflow/Cut.h"
#include "Superflow/Superlink.h"
#include "Superflow/Supervar.h"
#include "Superflow/Supersys.h"
#include "Superflow/EventFlags.h"

using namespace DataDefinitions;

namespace sflow {

    enum class SuperflowRunMode {
        nominal,
        nominal_and_weight_syst,
        single_event_syst,
        all_syst,
        data,
        null
    };

    enum class ATLAS_period { A, B, C, D, E, F, G, H, I, J, K, L, M, N, O, P, Q, R, S, T, U, V, W, X, Y, Z, null };
    enum class ATLAS_stream { Egamma, Muons, null };

    class Superflow : public SusyNtAna {

    public:
        Superflow();
        ~Superflow();

        // Cut Operators
        Superflow& operator<<(CutName cut_);
        void operator<<(std::function<bool(Superlink*)> cut_);

        // Var Operators
        void operator<<(std::function<double(Superlink*, var_float*)> var_);
        void operator<<(std::function<double(Superlink*, var_double*)> var_);
        void operator<<(std::function<int(Superlink*, var_int*)> var_);
        void operator<<(std::function<bool(Superlink*, var_bool*)> var_);
        void operator<<(std::function<void(Superlink*, var_void*)> var_);

        void operator<<(NewSystematic new_sys);
        void operator<<(TreeName tree_name);
        void operator<<(EventSystematic obj_);
        void operator<<(WeightSystematic obj_);
        void operator<<(SaveSystematic save_var);

        void operator<<(NewVar new_var_name);
        void operator<<(HFTname hft_name);
        void operator<<(SaveVar save_var);

        void Begin(TTree *tree); ///< called before looping on entries
        void Init(TTree *tree); ///< called when the TChain is attached
        void Terminate(); ///< called after looping is finished
        Bool_t Notify(); ///< called at each event
        Bool_t Process(Long64_t entry); ///< called at each event

        void setCountWeights(bool value); ///< Toggle the display of the weighted cuts. (default off)
        void setRunMode(SuperflowRunMode run_mode_);
        void setSingleEventSyst(SusyNtSys nt_syst_);
        void setChain(TChain* input_chain_);

    protected:
        void attach_superlink(Superlink* sl_);

        DilTrigLogic* m_trigObj; ///< trigger logic class
        MCWeighter* m_mcWeighter; ///< tool to determine the normalization

        bool assignNonStaticWeightComponents(
            Susy::SusyNtObject &ntobj,
            MCWeighter &weighter,
            const LeptonVector& leptons,
            const JetVector& jets,
            Supersys* super_sys,
            Superweight* weightComponents);

        double computeDileptonTriggerWeight(const LeptonVector &leptons, const SusyNtSys sys);

        double computeBtagWeight(const JetVector& jets, const Susy::Event* evt, SupersysWeight sys);

        double computeLeptonEfficiencySf(const Susy::Lepton &lep, const SupersysWeight sys);

        EventFlags computeEventFlags();

        // vector<Cut* > m_CutStore;
        vector<std::function<bool(Superlink*)>> m_LambdaCutStore;
        vector<string> m_LambdaCutStoreNames;

        vector<double> m_RawCounter;
        vector<double> m_WeightCounter;

        double m_passed;
        double m_weighted;
        bool m_countWeights;

        Superweight* m_weights;

        bool m_super_isData;

        string m_outputFileName;
        string m_entry_list_FileName;
        string m_tree_name_auto;
        TFile* m_outputFile;
        TFile* m_entryListFile;
        TTree* m_HFT;

        TFile** m_output_array;
        TTree** m_HFT_array;

        TEntryList* m_entry_list_total;
        TEntryList* m_entry_list_single_tree;

        ATLAS_period m_period;
        ATLAS_stream m_stream;

        std::function<double(Superlink*, var_float*)> m_nullExprFloat;
        std::function<double(Superlink*, var_double*)> m_nullExprDouble;
        std::function<int(Superlink*, var_int*)> m_nullExprInt;
        std::function<bool(Superlink*, var_bool*)> m_nullExprBool;
        std::function<void(Superlink*, var_void*)> m_nullExprVoid;

        vector<std::function<double(Superlink*, var_float*)>> m_varExprFloat;
        vector<std::function<double(Superlink*, var_double*)>> m_varExprDouble;
        vector<std::function<int(Superlink*, var_int*)>> m_varExprInt;
        vector<std::function<bool(Superlink*, var_bool*)>> m_varExprBool;
        vector<std::function<void(Superlink*, var_void*)>> m_varExprVoid;
        //vector<Float_t> m_varFloat;
        //vector<Double_t> m_varDouble;
        //vector<Int_t> m_varInt;

        Float_t* m_varFloat;
        Double_t* m_varDouble;
        Int_t* m_varInt;
        Bool_t* m_varBool;

        Float_t** m_varFloat_array;
        Double_t** m_varDouble_array;
        Int_t** m_varInt_array;
        Bool_t** m_varBool_array;

        SupervarState m_varState;

        vector<SupervarType> m_varType;
        vector<string> m_varNiceName;
        vector<string> m_varHFTName;
        bool m_superVar_hasFunction;
        bool m_superVar_hasNiceName;
        bool m_superVar_hasHFTName;
        int m_superVar_Untitled;

        SuperflowRunMode m_runMode;
        SupersysState m_sysState;

        Supersys m_sysTemplate;
        vector<Supersys> m_sysStore;
        bool m_sys_hasNiceName;
        bool m_sys_hasTreeName;
        bool m_sys_hasType;
        bool m_sys_hasSystematic;

        SusyNtSys m_singleEventSyst;
        Supersys* m_RunSyst;

        vector<int> index_weight_sys;
        vector<int> index_event_sys;

        int m_tree_leafs_size;
        int m_weight_leaf_offset;

    private:
        /// initialize weighter used for normalization
        bool initMcWeighter(TTree *tree);
        string app_name = "Superflow    ";
        int m_LambdaCutStoreUntitled;
        bool m_LambdaCutStore_Name_Exists;

        map<ATLAS_stream, map<ATLAS_period, string>> m_data_periods;
        map<SusyNtSys, string> m_NtSys_to_string;

        TChain* m_input_chain;

        const double m_epsilon = 1e-12;
    };

};