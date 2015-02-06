#include <cstdlib>
#include <string>
#include <cassert>
#include <cmath> // isnan
#include <cfloat> // FLT_MAX, FLT_MIN
#include <iomanip> // setw, setprecision
#include <iostream>
#include <fstream> 
#include <sstream>  // std::ostringstream
#include <dirent.h> // UNIX
#include <vector>
#include <map>

#include "TROOT.h"
//#include "TDirectory.h"
#include "TMath.h"

#include "TTree.h"
#include "TChain.h"
#include "TFile.h"
#include "TCut.h"
#include "TEventList.h"

#include "THStack.h"
#include "TStyle.h"
#include "TPad.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TH1.h"
#include "TH1F.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "TAxis.h"
#include "TLine.h"
#include "TLatex.h"
#include "TStyle.h"
using namespace std;

void listFiles(vector<string>& list_, string dir_);

#define MC_TEXT_DIR "/gdata/atlas/suneetu/Documents/LFV_Higgs2014/S9_Superflow/Superflow/data/"
#define RAW_SAMPLES_DIR "/data7/atlas/suneetu/Documents/LFV_Higgs2014/output/R12_September_15/"
//#define NEW_SAMPLES_DIR "/data7/atlas/suneetu/Documents/LFV_Higgs2014/output/R4_August_30/Processed/"
#define NEW_SAMPLES_DIR "/local/scratch/suneetu/R12_Output/"

#define OUT_FILENAME "LFVHSIG8TeV.root"

#define OUT_PREFIX "file_list_"

// const int n_files = 6;
// string files[] = { "Higgs", "top_MCNLO", "WW_Sherpa", "WZ_ZZ_Sherpa", "Zjets_AlpgenPythia", "Data" };
// string ext = ".txt";
// 
// string new_tree_names[] = { "Higgs", "Top", "WW", "ZV", "Zjets", "Data" };

const int n_files = 4;
string files[] = { "LFV_169670", "LFV_169671", "LFV_169672", "LFV_169673"};
string ext = ".txt";

string new_tree_names[] = { "169670", "169671", "169672", "169673" };

enum hft_sys_type {
    k_sys_object,
    k_sys_one_sided_object,
    k_sys_weight,
    k_sys_adhoc,
    k_N_sys_type
};

struct hft_systematic {
    hft_sys_type sys_type;
    string basename;
    string up_name;
    string down_name;
    string preface;
    string adhoc_tcut;
    string adhoc_tcut2;
    int binding_group;// inside a binding group, *trees* add linearly
};

vector<hft_systematic> defineSystematics()
{
    vector<hft_systematic> hft_syst;

    hft_systematic syst;

    syst.adhoc_tcut2 = "";

    // CENTRAL
    // these apply to all
    syst.sys_type = k_sys_one_sided_object;
    syst.up_name = "";
    syst.down_name = "";
    syst.preface = "";
    syst.adhoc_tcut = "";
    syst.binding_group = 0;
    syst.basename = "CENTRAL";
    hft_syst.push_back(syst);

    // OBJECT SYSTEMATICS
    // these apply to all
    syst.sys_type = k_sys_object;
    syst.up_name = "UP";
    syst.down_name = "DOWN";
    syst.preface = "";
    syst.adhoc_tcut = "";
    syst.binding_group = 0;

    // EESZ
    syst.basename = "EESZ";
    hft_syst.push_back(syst);

    // EER
    syst.basename = "EER";
    hft_syst.push_back(syst);

    // EESLOW
    syst.basename = "EESLOW";
    hft_syst.push_back(syst);

    // EESMAT
    syst.basename = "EESMAT";
    hft_syst.push_back(syst);

    // EESPS
    syst.basename = "EESPS";
    hft_syst.push_back(syst);

    // ID
    syst.basename = "ID";
    hft_syst.push_back(syst);

    // JES
    syst.basename = "JES";
    hft_syst.push_back(syst);

    // MS
    syst.basename = "MS";
    hft_syst.push_back(syst);

    // SCALEST
    syst.basename = "SCALEST";
    hft_syst.push_back(syst);

    // TES
    syst.basename = "TES";
    hft_syst.push_back(syst);

    // ONE-SIDED SYSTEMATICS
    // these apply to all
    syst.sys_type = k_sys_one_sided_object;
    syst.up_name = "";
    syst.down_name = "";
    syst.preface = "";
    syst.adhoc_tcut = "";
    syst.binding_group = 0;

    // JER
    syst.basename = "JER";
    hft_syst.push_back(syst);

    // RESOST
    syst.basename = "RESOST";
    hft_syst.push_back(syst);

    // DONE //
    return hft_syst;
}

string trim(string s)
{
    string val = s.erase(s.find_last_not_of(" \n\r\t") + 1);
    val.erase(0, s.find_first_not_of(" \n\r\t"));
    return val;
}
void listFiles(vector<string>& list_, string dir_);

int main(int argc, char** argv)
{
    // ROOT::Cintex::Cintex::Enable();

    vector<string> sample_dir_files;
    listFiles(sample_dir_files, RAW_SAMPLES_DIR);

    vector<hft_systematic> hft_syst = defineSystematics();

    vector<string> syst_strings;

    string outputFileName = string(NEW_SAMPLES_DIR) + OUT_FILENAME;
    TFile* outputFile = new TFile(outputFileName.data(), "RECREATE");
    outputFile->Close();
    delete outputFile;

    for (int i = 0; i < hft_syst.size(); i++) {
        if (hft_syst[i].sys_type == k_sys_object) {
            syst_strings.push_back(hft_syst[i].basename + hft_syst[i].up_name);
            syst_strings.push_back(hft_syst[i].basename + hft_syst[i].down_name);
        }
        else if (hft_syst[i].sys_type == k_sys_one_sided_object) {
            syst_strings.push_back(hft_syst[i].basename);
        }
    }

    for (int s_ = 0; s_ < syst_strings.size(); s_++) { // systematics
        for (int g_ = 0; g_ < n_files; g_++) { // known groups
            if (syst_strings[s_].compare("CENTRAL") != 0 && files[g_].compare("Data") == 0) continue;

            vector<string> sample_list;

            stringstream in_filename;
            in_filename << MC_TEXT_DIR << files[g_] << ext;

            ifstream mc_in_(in_filename.str().data());

            if (!mc_in_.is_open()) {
                cout << "File not open. Exiting." << endl;
                return 0;
            }
            else {
                cout << endl << "Opened " << in_filename.str() << endl;
            }

            while (!mc_in_.eof()) {
                string line_ = "";
                getline(mc_in_, line_);
                line_ = trim(line_);

                size_t found_ = line_.find_first_of('.');

                if (found_ != string::npos && line_ != "") {
                    string sample_ = line_.substr(found_ + 1);
                    sample_list.push_back(sample_);
                }
            }

            string tree_name = new_tree_names[g_] + "_" + syst_strings[s_];
            cout << "tree: " << tree_name << endl;

            outputFile = new TFile(outputFileName.data(), "UPDATE");
            outputFile->cd();

            TChain* merge_chain = new TChain(tree_name.data());
            TTree::SetMaxTreeSize(137438953472LL);

            int sum_trees = 0;

            for (int i = 0; i < sample_list.size(); i++) {
                string sample_treename = string("id_") + sample_list[i];
                string sample_filename = syst_strings[s_] + "_" + sample_list[i] + ".root";
                bool flag_found_root = false;

                auto file_nm = find(sample_dir_files.begin(), sample_dir_files.end(), sample_filename);
                if (file_nm != sample_dir_files.end()) {
                    flag_found_root = true;
                    // cout << "found: " << *file_nm << endl;
                }
                else {
                    cout << sample_filename << " not found           #%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#" << endl;
                }

                if (!flag_found_root) continue;

                string root_file_name = string(RAW_SAMPLES_DIR) + sample_filename;

                TFile* in_file = new TFile(root_file_name.data());
                TTree* in_tree = static_cast<TTree*>(in_file->Get(sample_treename.data()));

                int n_entries = 0;

                if (in_tree != nullptr) {
                    n_entries = in_tree->GetEntries();
                    cout << sample_filename << " " << to_string(n_entries) << endl;
                    sum_trees += n_entries;
                }

                if (n_entries == 0) continue;

                merge_chain->AddFile(root_file_name.data(), 0, sample_treename.data());

                delete in_tree;
                delete in_file;
                // cout << sample_treename << endl;
                // cout << sample_filename << endl;
                // cout << endl;
            }

            outputFile->cd();
            merge_chain->Merge(outputFile, 0, "fast");

            delete merge_chain;
            cout << "sum: " << sum_trees << endl;
            mc_in_.close();
        }
        // no code allowed here
    }



    // outputFile->Write();
    // outputFile->Close();
    // delete outputFile;

    return 0;
}

void listFiles(vector<string>& list_, string dir_)
{
    DIR *pDIR;
    struct dirent *entry;

    if ((pDIR = opendir(dir_.data()))) {
        while ((entry = readdir(pDIR))) {
            if (strcmp(entry->d_name, ".") != 0 && strcmp(entry->d_name, "..") != 0)
                list_.push_back(entry->d_name);
        }
        closedir(pDIR);
    }
}