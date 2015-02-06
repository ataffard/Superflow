// Copyright 2014, The ATLAS Collaboration
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <cmath>
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
#include "TColor.h"
#include "TStyle.h"

#include "Superflow/StringTools.h"

//#include "Cintex/Cintex.h"

using namespace std;
using namespace sflow;

#define SUFFIX_NAME "_3"

#define DEBUG_NTC true

#define DIR_0 "/scratch/suneetu/R35_Output_XAOD/"
#define DIR_1 "/data7/atlas/suneetu/Documents/UCI_XAOD_2014/C3_Superflow/susynt-write/prod/Superflow/plots/x7_testing/"

#define BG_FILE "LFVHBKG8TeV.root"
#define DATA_FILE "LFVHBKG8TeV.root"
//#define DATA_FILE "LFVHDAT8TeV.root"
#define FAKE_FILE "LFVHFAK8TeV.root"
#define TOP_FILE "LFVHBKG8TeV.root"
//#define TOP_FILE "LFVHBKG8TeV.root"

//#define BG_FILE "LFVHBKG8TeV.root"
//#define DATA_FILE "LFVHBKG8TeV.root"
//#define FAKE_FILE "LFVHFAK8TeV.root"
//#define TOP_FILE "LFVHBKG8TeV.root"

const string m_central = "CENTRAL";

const int n_hists = 1;
string hists_entries[] = {};

const double Pi = 3.141592653589793238462;

struct plot {
    string root_get_member;
    string root_get_member_name;
    string root_get_factor;
    string signal_region;
    string x_label;
    string y_label;
    string request_cuts;
    string reweight;
    string blind_data;
    bool is_Log;
    double x_bin_width;
    double x_range_max;
    double x_range_min;
    double y_range_max;
    double y_range_min;
    string long_name;
    int run_group;
    double x_red_arrow;
};

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

enum processed_tree_type {
    k_tree_mc,
    k_tree_signal,
    k_tree_fakes,
    k_tree_data,
    k_N_tree_type
};

struct processed_tree {
    processed_tree_type tree_type;
    string treename;
    string displayname;
    string blinding_cuts;
    string special_cuts;
    string file;
    int color;
    vector<int> systematic_groups;
    TEventList central_list;
    map<string, TEventList> all_lists;
    bool all_lists_exist;
};

void convertErrorsToPoisson(TH1* histo, TGraphAsymmErrors* graph)
{
    // Needed variables
    double value = 0;
    double error_poisson_up = 0;
    double error_poisson_down = 0;
    double alpha = 0.158655, beta = 0.158655; // 68%

    // loop over bins and overwrite values
    for (int i = 1; i < histo->GetNbinsX() + 1; i++) {
        value = histo->GetBinContent(i);
        if (value != 0) {
            error_poisson_up = 0.5*TMath::ChisquareQuantile(1 - beta, 2 * (value + 1)) - value;
            error_poisson_down = value - 0.5*TMath::ChisquareQuantile(alpha, 2 * value);
            graph->SetPoint(i - 1, histo->GetBinCenter(i), value);
            graph->SetPointError(i - 1, 0., 0., error_poisson_down, error_poisson_up);
        }
        else {
            graph->SetPoint(i - 1, histo->GetBinCenter(i), -10);
            graph->SetPointError(i - 1, 0., 0., 0., 0.);
        }
    }
}

void myText(Double_t x, Double_t y, Color_t color, const char *text)
{
    Double_t tsize = 0.055;
    TLatex l; //l.SetTextAlign(12); 
    l.SetTextSize(tsize);
    l.SetNDC();
    l.SetTextColor(color);
    l.DrawLatex(x, y, text);
}

vector<processed_tree> defineTrees()
{
    vector<processed_tree> hft_trees;

    processed_tree tree_;

    vector<int> grouplist;
    grouplist.clear();

    TEventList blank;
    tree_.central_list = blank;
    tree_.all_lists_exist = false;

    // monte-carlo trees
    // these apply to the group
    tree_.tree_type = k_tree_mc;
    tree_.blinding_cuts = "";
    tree_.special_cuts = "";

    tree_.file = BG_FILE;

    // Higgs
    tree_.treename = "Higgs";
    tree_.displayname = "Higgs";
    tree_.color = TColor::GetColor("#ffff99");
    grouplist.clear();
    grouplist.push_back(0);
    grouplist.push_back(7);// XS
    tree_.systematic_groups = grouplist;
    hft_trees.push_back(tree_);

    // ZV
    tree_.treename = "ZV";
    tree_.displayname = "ZV";
    tree_.color = TColor::GetColor("#ccff66");
    grouplist.clear();
    grouplist.push_back(0);
    grouplist.push_back(3);// XS
    tree_.systematic_groups = grouplist;
    hft_trees.push_back(tree_);

    // Zjets
    tree_.treename = "Zjets";
    tree_.displayname = "Z+jets";
    tree_.color = TColor::GetColor("#ffcc00");
    grouplist.clear();
    grouplist.push_back(0);
    grouplist.push_back(6);// XS
    tree_.systematic_groups = grouplist;
    hft_trees.push_back(tree_);

    // Top
    tree_.treename = "Top";
    tree_.displayname = "tt+Wt";
    tree_.color = TColor::GetColor("#a50021");
    grouplist.clear();
    grouplist.push_back(0);
    grouplist.push_back(5);// XS
    tree_.systematic_groups = grouplist;
    hft_trees.push_back(tree_);

    // WW
    tree_.treename = "WW";
    tree_.displayname = "WW";
    tree_.color = TColor::GetColor("#336699");
    grouplist.clear();
    grouplist.push_back(0);
    grouplist.push_back(4);// XS
    tree_.systematic_groups = grouplist;
    hft_trees.push_back(tree_);

    tree_.file = DATA_FILE;

    // DATA trees
    tree_.treename = "Data";
    tree_.displayname = "Data";
    tree_.tree_type = k_tree_data;
    tree_.blinding_cuts = "";// 
    tree_.special_cuts = "";
    tree_.color = kBlack;
    grouplist.clear();// no systematics on data
    tree_.systematic_groups = grouplist;
    hft_trees.push_back(tree_);

    return hft_trees;
}

vector<hft_systematic> defineSystematics()
{
    vector<hft_systematic> hft_syst;

    hft_systematic syst;

    syst.adhoc_tcut2 = "";

    // DONE //
    return hft_syst;
}

vector<plot> definePlots()
{
    vector<plot> vect_plot;
    plot plot_;

    plot_.is_Log = false;
    plot_.root_get_member_name = "";
    plot_.reweight = "";
    plot_.x_red_arrow = 0.0;
    plot_.run_group = 0;

    // plot_.is_Log = true;


#include "Superflow/define_test_SR.cxx"
// #include "Superflow/define_test_VR_5.cxx"
// #include "Superflow/define_test_ptconeplots.cxx"

    return vect_plot;
}

int main(int argc, char* argv[])
{
    string request_plot = "";
    bool plot_requested = false;

    string request_var = "";
    bool var_requested = false;

    int request_group = 0;
    bool group_requested = false;

    bool no_systematics = false;

    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "/p") == 0) {
            request_plot = argv[++i];
            if (request_plot != "") plot_requested = true;
        }
        else if (strcmp(argv[i], "/g") == 0) {
            request_group = atoi(argv[++i]);
            if (request_group > 0) group_requested = true;
        }
        else if (strcmp(argv[i], "/v") == 0) {
            request_var = argv[++i];
            if (request_var != "") var_requested = true;
        }
        else if (strcmp(argv[i], "/nosys") == 0) {
            no_systematics = true;
        }
        else {
            cout << "LFVH_Plots    Error (fatal): Bad arguments." << endl;
            exit(1);
        }
    }

    vector<plot> request_plots = definePlots();
    vector<processed_tree> hft_trees = defineTrees();
    vector<hft_systematic> hft_syst = defineSystematics();

    // define the all-important TCut strings
    map <string, string> tcuts;
    stringstream tcuts_join;

    //SR_OPT_no_jets
    tcuts_join << "isOS && lept1Pt>30. && lept2Pt>20. && nCentralBJets>0 && nForwardJets==0 && MT2 > 40.0";
    tcuts_join << " && ((lept1isEle && lept1d0_sig<5.0 && lept1z0sin_theta<0.4) || (lept1isMu && lept1d0_sig<3.0 && lept1z0sin_theta<1.0)) && ((lept2isEle && lept2d0_sig<5.0 && lept2z0sin_theta<0.4) || (lept2isMu && lept2d0_sig<3.0 && lept2z0sin_theta<1.0))";
    tcuts_join << " && ((lept1isEle && lept1ptcone30_Pt<0.16 && lept1etcone20_Pt<0.18) || (lept1isMu && lept1ptcone30_Pt<0.12)) && ((lept2isEle && lept2ptcone30_Pt<0.16 && lept2etcone20_Pt<0.18) || (lept2isMu && lept2ptcone30_Pt<0.12))";
    tcuts["VR_Top"] = tcuts_join.str(); tcuts_join.str("");// clear stringstream

    tcuts_join << "isOS && lept1Pt>30. && lept2Pt>20. && nCentralBJets==0 && nForwardJets==0 && nCentralLightJets<2 && MT2>40.0";
    tcuts_join << " && ((lept1isEle && lept1d0_sig<5.0 && lept1z0sin_theta<0.4) || (lept1isMu && lept1d0_sig<3.0 && lept1z0sin_theta<1.0)) && ((lept2isEle && lept2d0_sig<5.0 && lept2z0sin_theta<0.4) || (lept2isMu && lept2d0_sig<3.0 && lept2z0sin_theta<1.0))";
    tcuts_join << " && ((lept1isEle && lept1ptcone30_Pt<0.16 && lept1etcone20_Pt<0.18) || (lept1isMu && lept1ptcone30_Pt<0.12)) && ((lept2isEle && lept2ptcone30_Pt<0.16 && lept2etcone20_Pt<0.18) || (lept2isMu && lept2ptcone30_Pt<0.12))";
    tcuts["VR_WW"] = tcuts_join.str(); tcuts_join.str("");// clear stringstream

    tcuts_join << "isMuMu && isOS && lept1Pt>35. && lept2Pt>20. && nCentralBJets==0 && nForwardJets==0 && nCentralLightJets<2 && ll_M>(91.2-10.) && ll_M<(91.2+10.)";
    tcuts_join << " && ((lept1isEle && lept1d0_sig<5.0 && lept1z0sin_theta<0.4) || (lept1isMu && lept1d0_sig<3.0 && lept1z0sin_theta<1.0)) && ((lept2isEle && lept2d0_sig<5.0 && lept2z0sin_theta<0.4) || (lept2isMu && lept2d0_sig<3.0 && lept2z0sin_theta<1.0))";
    tcuts_join << " && ((lept1isEle && lept1ptcone30_Pt<0.16 && lept1etcone20_Pt<0.18) || (lept1isMu && lept1ptcone30_Pt<0.12)) && ((lept2isEle && lept2ptcone30_Pt<0.16 && lept2etcone20_Pt<0.18) || (lept2isMu && lept2ptcone30_Pt<0.12))";
    tcuts["VR_Zmm"] = tcuts_join.str(); tcuts_join.str("");// clear stringstream

    tcuts_join << "isElEl && isOS && lept1Pt>35. && lept2Pt>20. && nCentralBJets==0 && nForwardJets==0 && nCentralLightJets<2 && ll_M>(91.2-10.) && ll_M<(91.2+10.)";
    tcuts_join << " && ((lept1isEle && lept1d0_sig<5.0 && lept1z0sin_theta<0.4) || (lept1isMu && lept1d0_sig<3.0 && lept1z0sin_theta<1.0)) && ((lept2isEle && lept2d0_sig<5.0 && lept2z0sin_theta<0.4) || (lept2isMu && lept2d0_sig<3.0 && lept2z0sin_theta<1.0))";
    tcuts_join << " && ((lept1isEle && lept1ptcone30_Pt<0.16 && lept1etcone20_Pt<0.18) || (lept1isMu && lept1ptcone30_Pt<0.12)) && ((lept2isEle && lept2ptcone30_Pt<0.16 && lept2etcone20_Pt<0.18) || (lept2isMu && lept2ptcone30_Pt<0.12))";
    tcuts["VR_Zee"] = tcuts_join.str(); tcuts_join.str("");// clear stringstream


    tcuts_join << "isOS && (isEM|| isME) && !(ll_M>(91.2-10.) && ll_M<(91.2+10.)) && met>50. && nCentralBJets==0";
    tcuts["VR_X5"] = tcuts_join.str(); tcuts_join.str("");// clear stringstream

    tcuts_join << "isOS && isMuMu && nCentralBJets==0";
    tcuts["VR_MuMu"] = tcuts_join.str(); tcuts_join.str("");// clear stringstream

    tcuts["none"] = "1";

    TCanvas* ttemp = new TCanvas("canvas_1", "", 768, 768);

    vector<vector<vector<double>>> table_columns;
    vector<string> sr;

    string last_sr = "";

    for (int r_ = 0; r_ < request_plots.size(); r_++) {

        plot plot_ = request_plots[r_];
        sr.push_back(plot_.signal_region);

        if (plot_requested && plot_.long_name.compare(request_plot) != 0) { continue; };
        if (group_requested && plot_.run_group != request_group) { continue; };

        if (plot_.signal_region.compare(last_sr) != 0) {
            for (int t_ = 0; t_ < hft_trees.size(); t_++) {
                hft_trees[t_].all_lists_exist = false;
            }
            last_sr = plot_.signal_region;
        }

        vector<vector<vector<double>>> table_col_;// column is per process

        cout << setprecision(4);
        cout << std::fixed;
        cout << "+++++++++++++ " << plot_.signal_region << " +++++++++++++" << endl;

        // Complete/finalize TCut strings.
        stringstream signal_region_stream;
        signal_region_stream << "( ";
        if (plot_.request_cuts.compare("") != 0) {
            cout << "cut requested: " << plot_.request_cuts << endl;
            signal_region_stream << plot_.request_cuts << " && ";
        }
        signal_region_stream << tcuts[plot_.signal_region];
        signal_region_stream << " )";// Don't forget the close parentheses

        string signal_tcut = signal_region_stream.str();

        cout << signal_region_stream.str() << endl;

        TH1F** hist_central = new TH1F*[hft_trees.size()];
        for (int t_ = 0; t_ < hft_trees.size(); t_++) {
            stringstream hist_;
            if (plot_.root_get_member_name.compare("") == 0) {
                hist_ << "h_" << plot_.signal_region << "_" << plot_.root_get_member << "_" << hft_trees[t_].treename;
            }
            else {
                hist_ << "h_" << plot_.signal_region << "_" << plot_.root_get_member_name << "_" << hft_trees[t_].treename;
            }
            // calculate the number of bins
            int n_bins = floor((plot_.x_range_max - plot_.x_range_min) / (plot_.x_bin_width) + 0.5);
            hist_central[t_] = new TH1F(hist_.str().data(), hist_.str().data(), n_bins, plot_.x_range_min, plot_.x_range_max);
            hist_central[t_]->Sumw2();
        }

        int num_bins = hist_central[0]->GetNbinsX();

        double*** systematic_up = new double**[num_bins];
        for (int b_ = 0; b_ < num_bins; b_++) {
            systematic_up[b_] = new double*[hft_trees.size()];
            for (int t_ = 0; t_ < hft_trees.size(); t_++) {
                systematic_up[b_][t_] = new double[hft_syst.size()];
                for (int s_ = 0; s_ < hft_syst.size(); s_++) {
                    systematic_up[b_][t_][s_] = 0;
                }
            }
        }

        double*** systematic_down = new double**[num_bins];
        for (int b_ = 0; b_ < num_bins; b_++) {
            systematic_down[b_] = new double*[hft_trees.size()];
            for (int t_ = 0; t_ < hft_trees.size(); t_++) {
                systematic_down[b_][t_] = new double[hft_syst.size()];
                for (int s_ = 0; s_ < hft_syst.size(); s_++) {
                    systematic_down[b_][t_][s_] = 0;
                }
            }
        }

        vector<TH1F> total_ststematics;
        vector<vector<double>> tree_total;
        vector<vector<double>> tree_stat_error;

        vector<double> single_tree_total;

        for (int b_ = 0; b_ < num_bins; b_++) {
            vector<double> entry;
            tree_total.push_back(entry);
            tree_stat_error.push_back(entry);
        }

        // define central yields table
        for (int t_ = 0; t_ < hft_trees.size(); t_++) {
            stringstream filename_stream;
            filename_stream << DIR_0 << hft_trees[t_].file;
            string filename_ = filename_stream.str();

            stringstream chain_name;// select CENTRAL
            chain_name << hft_trees[t_].treename << "_" << m_central;

            TChain* chain_;
            //TChain* un_chain;
            chain_ = new TChain(chain_name.str().data());
            chain_->Add(filename_.data());

            //un_chain = new TChain(chain_name.str().data());
            //un_chain->Add(filename_.data());

            if (chain_->IsZombie()) {
                cout << chain_name.str() << "not found." << endl;
                continue;
            }

            TCut sel = TCut(signal_tcut.data());// key TCut

            TCut weight = TCut("eventweight");// selects eventweight leaf
            TCut blind = TCut(hft_trees[t_].blinding_cuts.data());

            // START reweighting segment
            TCut tc_reweight = TCut("1.0");
            if (plot_.reweight.compare("") != 0 && hft_trees[t_].tree_type == k_tree_mc) {
                tc_reweight = TCut(plot_.reweight.data());
            }
            // END reweighting segment

            // START selection and event-list
            if (!hft_trees[t_].all_lists_exist || plot_.request_cuts.compare("") != 0) {
                string list_name = string("list_") + hist_central[t_]->GetName();

                stringstream draw_list;
                draw_list << ">> " << list_name; // hist_central[t_]->GetName();
                stringstream draw_hist;
                draw_hist << plot_.root_get_member << plot_.root_get_factor << " >>+ " << hist_central[t_]->GetName();

                chain_->Draw(draw_list.str().data(), (sel + blind)); // drawing list


                TEventList* eventList = (TEventList*)gDirectory->Get(list_name.data());
                chain_->SetEventList(eventList);

                hft_trees[t_].central_list = *eventList;

                chain_->Draw(draw_hist.str().data(), weight * tc_reweight); // drawing histogram
                //un_chain->Draw(draw_hist.str().data(), (sel + blind)); // drawing list
            }
            else {
                stringstream draw_hist;
                draw_hist << plot_.root_get_member << plot_.root_get_factor << " >>+ " << hist_central[t_]->GetName();

                chain_->SetEventList(&hft_trees[t_].central_list);
                chain_->Draw(draw_hist.str().data(), weight * tc_reweight); // drawing histogram
            }
            // END selection and event-list

            double error_ = 0;
            double integral_ = hist_central[t_]->IntegralAndError(0, hist_central[t_]->GetNbinsX() + 1, error_);

            cout << hist_central[t_]->GetName() << ":" << integral_ << " (" << hist_central[t_]->GetEntries() << ") " << endl;

            for (int b_ = 1; b_ <= num_bins; b_++) {
                tree_total[b_ - 1].push_back(hist_central[t_]->GetBinContent(b_));
                tree_stat_error[b_ - 1].push_back(hist_central[t_]->GetBinError(b_));
            }



            single_tree_total.push_back(integral_);

            //re tree_total.push_back(integral_);
            //re tree_stat_error.push_back(error_);

            // /if (hft_trees[t_].tree_type == k_tree_data) {
            // /	//Int_t runNumber = 0;
            // /	//Int_t eventNumber = 0;
            // /	//chain_->SetBranchAddress("runNumber", &runNumber);
            // /	//chain_->SetBranchAddress("eventNumber", &eventNumber);
            // /
            // /	stringstream select_;
            // /	select_ << ">> data_events";
            // /	chain_->Draw(select_.str().data(), sel * weight);
            // /	TEventList* eventList = (TEventList*)gDirectory->Get("data_events");
            // /	chain_->SetEventList(eventList);
            // /
            // /	chain_->Scan("runNumber:eventNumber");
            // /
            // /	//int n_entries = chain_->GetEntries();
            // /	//cout << chain_name.str() << ": " << n_entries << endl;
            // /	//for (int i = 0; i < n_entries; i++) {
            // /	//	chain_->GetEntry(i);
            // /	//	cout << "DATA: runNumber-" << runNumber << " eventNumber-" << eventNumber << endl;
            // /	//}
            // /
            // /	delete eventList;
            // /} // Direct reading

            delete chain_;
            //delete un_chain;
            // delete eventList; in gDirectory, not needed.draw_hist
        }

        cout << ">";

        for (int t_ = 0; t_ < hft_trees.size(); t_++) {
            stringstream filename_stream;
            filename_stream << DIR_0 << hft_trees[t_].file;
            string filename_ = filename_stream.str();

            if (!DEBUG_NTC) {
                cout << "<< " << hft_trees[t_].treename;
            }
            int lcount = 0;

            for (int s_ = 0; s_ < hft_syst.size(); s_++) {

                if (no_systematics) continue;

                if (hft_trees[t_].systematic_groups.size() == 0) continue;

                bool add_this_systematic = false;
                for (int i = 0; i < hft_trees[t_].systematic_groups.size(); i++) {
                    if (hft_trees[t_].systematic_groups[i] == hft_syst[s_].binding_group) add_this_systematic = true;
                }

                if (!add_this_systematic) continue;

                stringstream chain_name_up;
                stringstream chain_name_down;

                // classify systematics
                if (hft_syst[s_].sys_type == k_sys_object) {
                    chain_name_up << hft_trees[t_].treename << "_" << hft_syst[s_].basename << hft_syst[s_].up_name;
                    chain_name_down << hft_trees[t_].treename << "_" << hft_syst[s_].basename << hft_syst[s_].down_name;
                }
                else if (hft_syst[s_].sys_type == k_sys_one_sided_object) {
                    chain_name_up << hft_trees[t_].treename << "_" << hft_syst[s_].basename;
                    chain_name_down << hft_trees[t_].treename << "_" << hft_syst[s_].basename;
                }
                else if (hft_syst[s_].sys_type == k_sys_weight) {
                    chain_name_up << hft_trees[t_].treename << "_" << m_central;
                    chain_name_down << hft_trees[t_].treename << "_" << m_central;
                    //reweight_flag = true;
                }
                else if (hft_syst[s_].sys_type == k_sys_adhoc) {
                    chain_name_up << hft_trees[t_].treename << "_" << m_central;
                    chain_name_down << hft_trees[t_].treename << "_" << m_central;
                    //reweight_flag = true;
                }
                else {
                    continue;
                }

                TChain* chain_up;
                chain_up = new TChain(chain_name_up.str().data());
                chain_up->Add(filename_.data());

                TChain* chain_down;
                chain_down = new TChain(chain_name_down.str().data());
                chain_down->Add(filename_.data());

                if (chain_up->IsZombie()) {
                    cout << chain_name_up.str() << "not found." << endl;
                    continue;
                }
                if (chain_down->IsZombie()) {
                    cout << chain_name_down.str() << "not found." << endl;
                    continue;
                }

                TCut sel = TCut(signal_tcut.data());// key TCut
                TCut weight = TCut("eventweight");// selects eventweight leaf
                //TCut blind = TCut(hft_trees[t_].blinding_cuts.data());

                // apply direct weighting
                TCut syst_var_up;
                TCut syst_var_down;
                if (hft_syst[s_].sys_type == k_sys_weight) {
                    stringstream weight_name_up;
                    weight_name_up << hft_syst[s_].preface << hft_syst[s_].basename << hft_syst[s_].up_name;
                    stringstream weight_name_down;
                    weight_name_down << hft_syst[s_].preface << hft_syst[s_].basename << hft_syst[s_].down_name;

                    syst_var_up = TCut(weight_name_up.str().data());
                    syst_var_down = TCut(weight_name_down.str().data());
                }
                else if (hft_syst[s_].sys_type == k_sys_adhoc) {
                    syst_var_up = TCut(hft_syst[s_].adhoc_tcut.data());
                    syst_var_down = TCut(hft_syst[s_].adhoc_tcut2.data());
                }
                else {
                    syst_var_up = TCut("1.0");
                    syst_var_down = TCut("1.0");
                }

                // START reweighting segment
                TCut tc_reweight = TCut("1.0");
                if (plot_.reweight.compare("") != 0 && hft_trees[t_].tree_type == k_tree_mc) {
                    tc_reweight = TCut(plot_.reweight.data());
                }
                // END reweighting segment

                int n_bins = floor((plot_.x_range_max - plot_.x_range_min) / (plot_.x_bin_width) + 0.5);

                // temp histogram to get the variation
                TH1F* temp_hist_up = new TH1F("temp_hist_up", "temp_hist_up", n_bins, plot_.x_range_min, plot_.x_range_max);
                TH1F* temp_hist_down = new TH1F("temp_hist_down", "temp_hist_down", n_bins, plot_.x_range_min, plot_.x_range_max);

                stringstream draw_up;
                draw_up << plot_.root_get_member << plot_.root_get_factor << " >> +" << "temp_hist_up";
                stringstream draw_down;
                draw_down << plot_.root_get_member << plot_.root_get_factor << " >> +" << "temp_hist_down";

                if (hft_syst[s_].sys_type == k_sys_weight || hft_syst[s_].sys_type == k_sys_adhoc) {
                    chain_up->SetEventList(&hft_trees[t_].central_list);
                    chain_down->SetEventList(&hft_trees[t_].central_list);
                    chain_up->Draw(draw_up.str().data(), weight * syst_var_up * tc_reweight);
                    chain_down->Draw(draw_down.str().data(), weight * syst_var_down * tc_reweight);
                }
                // else if (hft_syst[s_].sys_type == k_sys_adhoc) {
                //     chain_up->Draw(draw_up.str().data(), sel * weight * syst_var_up * tc_reweight);
                //     chain_down->Draw(draw_down.str().data(), sel * weight * syst_var_down * tc_reweight);
                // } // needless duplication
                else if (hft_trees[t_].all_lists_exist && plot_.request_cuts.compare("") == 0) {
                    string list_name_up = string("eventlist_") + hft_syst[s_].basename + hft_syst[s_].up_name;
                    string list_name_down = string("eventlist_") + hft_syst[s_].basename + hft_syst[s_].down_name;

                    chain_up->SetEventList(&hft_trees[t_].all_lists[list_name_up]);
                    chain_down->SetEventList(&hft_trees[t_].all_lists[list_name_down]);

                    chain_up->Draw(draw_up.str().data(), weight * tc_reweight);
                    chain_down->Draw(draw_down.str().data(), weight * tc_reweight);
                }
                else {
                    string list_name_up = string("eventlist_") + hft_syst[s_].basename + hft_syst[s_].up_name;
                    string list_name_down = string("eventlist_") + hft_syst[s_].basename + hft_syst[s_].down_name;

                    string draw_list_up = string(">> ") + list_name_up;
                    string draw_list_down = string(">> ") + list_name_down;

                    chain_up->Draw(draw_list_up.data(), sel); // drawing list
                    chain_down->Draw(draw_list_down.data(), sel); // drawing list

                    TEventList* eventList_up = (TEventList*)gDirectory->Get(list_name_up.data());
                    TEventList* eventList_down = (TEventList*)gDirectory->Get(list_name_down.data());

                    chain_up->SetEventList(eventList_up);
                    chain_down->SetEventList(eventList_down);

                    hft_trees[t_].all_lists[list_name_up] = *eventList_up;
                    hft_trees[t_].all_lists[list_name_down] = *eventList_down;

                    chain_up->Draw(draw_up.str().data(), weight * tc_reweight);
                    chain_down->Draw(draw_down.str().data(), weight * tc_reweight);
                }

                double int_up = temp_hist_up->Integral(0, -1) - single_tree_total[t_];
                double int_down = temp_hist_down->Integral(0, -1) - single_tree_total[t_];

                bool use_up = false;
                bool use_down = false;
                // sign and sidedness adjustments
                if (hft_syst[s_].sys_type == k_sys_one_sided_object) {
                    if (int_up > 0) {
                        int_down = 0.0;
                        use_up = true;
                    }
                    else {
                        int_up = 0.0;
                        use_down = true;
                    }
                }
                else if (hft_syst[s_].sys_type == k_sys_adhoc) {
                    int_up = abs(int_up);
                    int_down = -abs(int_down);
                }

                // systematic_up[t_][s_] = int_up;
                // systematic_down[t_][s_] = int_down;

                for (int b_ = 1; b_ <= num_bins; b_++) {
                    systematic_up[b_ - 1][t_][s_] = temp_hist_up->GetBinContent(b_) - tree_total[b_ - 1][t_];
                    systematic_down[b_ - 1][t_][s_] = temp_hist_down->GetBinContent(b_) - tree_total[b_ - 1][t_];

                    if (use_up)  systematic_down[b_ - 1][t_][s_] = 0;
                    if (use_down) systematic_up[b_ - 1][t_][s_] = 0;
                }

                if (DEBUG_NTC) {
                    // cout << hft_trees[t_].treename << " > " << hft_syst[s_].basename << " ";
                    // cout << tree_total[t_] << ",+" << int_up / tree_total[t_] << ",-" << int_down / tree_total[t_] << endl;
                    //cout << hft_trees[t_].treename << " > " << hft_syst[s_].basename << " ";
                    //cout << single_tree_total[t_] << ",+" << int_up << ",-" << int_down << endl;

                    cout << pad_width(hft_trees[t_].treename + " > " + hft_syst[s_].basename, 20);
                    cout << pad_width(to_string(single_tree_total[t_]) + ",", 16);
                    cout << pad_width("+" + to_string(int_up), 12) << pad_width("-" + to_string(int_down), 12) << endl;

                }

                delete temp_hist_up;
                temp_hist_up = 0;
                delete temp_hist_down;
                temp_hist_down = 0;
                delete chain_up;
                delete chain_down;

                lcount++;

                if (!(lcount % 10 == 0 || lcount == 0)) cout << "+";
                else cout << "/";
                cout << flush;
            }

            hft_trees[t_].all_lists_exist = true; // set when you've run through the tree once

            if (!DEBUG_NTC) {
                cout << endl;
            }
        }

        for (int b_ = 1; b_ <= num_bins; b_++) {
            vector<vector<double> > entry;
            // for (int t_ = 0; t_ < hft_trees.size(); t_++) {
            // 
            //     vector<double> single_entry(4, 0.0);
            //     entry.push_back(single_entry);
            // 
            // }
            table_col_.push_back(entry);
        }

        // per grouping table
        for (int b_ = 0; b_ < num_bins; b_++) {
            for (int t_ = 0; t_ < hft_trees.size(); t_++) {
                vector<double> single_entry(4, 0.0);

                double syst_up = 0;
                double syst_down = 0;

                for (int s_ = 0; s_ < hft_syst.size(); s_++) {
                    double up_var = systematic_up[b_][t_][s_];
                    double down_var = systematic_down[b_][t_][s_];

                    if (up_var > 0) syst_up += up_var * up_var;
                    else syst_down += up_var * up_var;

                    if (down_var < 0) syst_down += down_var * down_var;
                    else syst_up += down_var * down_var;
                }

                syst_up = sqrt(syst_up);
                syst_down = sqrt(syst_down);

                single_entry[0] = tree_total[b_][t_];
                single_entry[1] = tree_stat_error[b_][t_];
                single_entry[2] = syst_up;
                single_entry[3] = syst_down;

                table_col_[b_].push_back(single_entry);
            }
        }

        // gather systematic groups
        vector<int> syst_groups;
        for (int s_ = 0; s_ < hft_syst.size(); s_++) {
            bool add_group = true;
            for (int g_ = 0; g_ < syst_groups.size(); g_++) {
                if (syst_groups[g_] == hft_syst[s_].binding_group) add_group = false;
            }

            if (add_group) syst_groups.push_back(hft_syst[s_].binding_group);
        }

        vector<double> total_stat_up;
        vector<double> total_stat_down;

        for (int b_ = 0; b_ < num_bins; b_++) {
            total_stat_up.push_back(0);
            total_stat_down.push_back(0);
        }

        // totals
        for (int b_ = 0; b_ < num_bins; b_++) {
            for (int g_ = 0; g_ < syst_groups.size(); g_++) {
                for (int s_ = 0; s_ < hft_syst.size(); s_++) {
                    if (hft_syst[s_].binding_group != syst_groups[g_]) continue;

                    double syst_up = 0;
                    double syst_down = 0;

                    for (int t_ = 0; t_ < hft_trees.size(); t_++) {

                        bool add_this_systematic = false;
                        for (int i = 0; i < hft_trees[t_].systematic_groups.size(); i++) {
                            if (hft_trees[t_].systematic_groups[i] == hft_syst[s_].binding_group) add_this_systematic = true;
                        }

                        if (!add_this_systematic) continue;

                        double up_var = systematic_up[b_][t_][s_];
                        double down_var = systematic_down[b_][t_][s_];

                        syst_up += up_var;
                        syst_down += down_var;
                    }

                    //total_stat_up += syst_up * syst_up;
                    //total_stat_down += syst_down * syst_down;

                    if (syst_up > 0) {
                        total_stat_up[b_] += syst_up * syst_up;
                    }
                    else {
                        total_stat_down[b_] += syst_up * syst_up;
                    }

                    if (syst_down < 0) {
                        total_stat_down[b_] += syst_down * syst_down;
                    }
                    else {
                        total_stat_up[b_] += syst_down * syst_down;
                    }
                }
            }
        }

        for (int b_ = 0; b_ < num_bins; b_++) {
            total_stat_up[b_] = sqrt(total_stat_up[b_]);
            total_stat_down[b_] = sqrt(total_stat_down[b_]);
        }

        vector<double> realtotal;
        vector<double> staterror;
        for (int b_ = 0; b_ < num_bins; b_++) {
            realtotal.push_back(0);
            staterror.push_back(0);
            for (int t_ = 0; t_ < hft_trees.size(); t_++) {
                if (hft_trees[t_].tree_type == k_tree_data) continue;
                realtotal[b_] += tree_total[b_][t_];
                staterror[b_] += tree_stat_error[b_][t_] * tree_stat_error[b_][t_];
            }
            staterror[b_] = sqrt(staterror[b_]);
        }

        //cout << "Total" << " > " << realtotal << " \xB1" << staterror;
        //cout << " + " << total_stat_up << " - " << total_stat_down << endl;

        vector<double> total_entry(4, 0.0);

        // realtotal[b_];        // these are the ultimate results.
        // staterror[b_];        // these are the ultimate results.
        // total_stat_up[b_];    // these are the ultimate results.
        // total_stat_down[b_];  // these are the ultimate results.

        /// /// /// /// ///
        /// /// /// /// ///
        /// /// /// /// ///
        /// /// /// /// ///
        /// /// /// /// ///

        double* nominal_yield = new double[hist_central[0]->GetNbinsX()];
        for (int j = 0; j < hist_central[0]->GetNbinsX(); j++) nominal_yield[j] = 0.0;
        // calculate and store, per bin nominal_yield;
        TH1F* hist_prediction = new TH1F();
        delete hist_prediction;
        for (int t_ = 0; t_ < hft_trees.size(); t_++) {
            if (hft_trees[t_].tree_type == k_tree_mc) {
                hist_prediction = (TH1F*)hist_central[t_]->Clone();
                break;
            }
        }
        bool skip_init = true;
        for (int t_ = 0; t_ < hft_trees.size(); t_++) {
            if (hft_trees[t_].tree_type == k_tree_mc || hft_trees[t_].tree_type == k_tree_fakes) {
                if (skip_init) {
                    skip_init = false;
                }
                else {
                    hist_prediction->Add(hist_central[t_], 1.0);
                }
            }
        }
        for (int j = 1; j <= hist_central[0]->GetNbinsX(); j++) {// note bin numbering convention
            nominal_yield[j - 1] = hist_prediction->GetBinContent(j);
        }
        delete hist_prediction;

        // X BIN CENTER
        double* x_bin_center = new double[hist_central[0]->GetNbinsX()];
        for (int j = 1; j <= hist_central[0]->GetNbinsX(); j++) x_bin_center[j - 1] = hist_central[0]->GetBinCenter(j);
        // X BIN DOWN
        double* x_bin_down = new double[hist_central[0]->GetNbinsX()];
        for (int j = 0; j < hist_central[0]->GetNbinsX(); j++) x_bin_down[j] = hist_central[0]->GetBinWidth(j) / 2;
        // X BIN UP
        double* x_bin_up = new double[hist_central[0]->GetNbinsX()];
        for (int j = 0; j < hist_central[0]->GetNbinsX(); j++) x_bin_up[j] = hist_central[0]->GetBinWidth(j) / 2;
        // Y BIN DOWN
        double* y_bin_down = new double[hist_central[0]->GetNbinsX()];
        for (int j = 0; j < hist_central[0]->GetNbinsX(); j++) y_bin_down[j] = sqrt(total_stat_down[j] * total_stat_down[j] + staterror[j] * staterror[j]);
        // Y BIN UP
        double* y_bin_up = new double[hist_central[0]->GetNbinsX()];
        for (int j = 0; j < hist_central[0]->GetNbinsX(); j++) y_bin_up[j] = sqrt(total_stat_up[j] * total_stat_up[j] + staterror[j] * staterror[j]);

        // order histograms smallest to largest
        vector<int> tree_ordering;
        {
            map<double, int> tree_map;
            for (int t_ = 0; t_ < hft_trees.size(); t_++) {
                tree_map[single_tree_total[t_]] = t_;
            }
            for (auto t_ = tree_map.begin(); t_ != tree_map.end(); t_++) {
                tree_ordering.push_back(t_->second);
            }
        }

        cout << "-------------" << endl;

        // Output to pdf or eps

        gROOT->SetStyle("Plain");
        gStyle->SetOptStat(false);
        gStyle->SetTitleBorderSize(0);
        gStyle->SetLineWidth(1);
        gStyle->SetTextFont(42);

        gStyle->SetPadTickX(1);
        gStyle->SetPadTickY(1);

        gStyle->SetErrorX(0.0);

        stringstream out_name;
        out_name << plot_.signal_region << "_" << plot_.root_get_member;
        stringstream out_file;
        out_file << DIR_1 << plot_.long_name << SUFFIX_NAME << ".eps";

        stringstream out_file_image;
        out_file_image << DIR_1 << plot_.long_name << SUFFIX_NAME << ".png";

        TCanvas* tc = new TCanvas("canvas", "", 768, 768);
        TLegend* tl = new TLegend(0.68, 0.49, 0.88, 0.87);

        tc->Divide(1, 2);

        TPad* canvas_up = (TPad*)tc->GetListOfPrimitives()->FindObject("canvas_1");
        TPad* canvas_dw = (TPad*)tc->GetListOfPrimitives()->FindObject("canvas_2");

        // define the size
        double up_height = 0.75;
        double font_size_dw = 0.73;
        double dw_height = 0.31;

        // set pad size
        canvas_up->SetPad(0., 1 - up_height, 1., 1.);
        canvas_dw->SetPad(0., 0., 1., dw_height);
        canvas_up->SetFrameFillColor(0);
        canvas_up->SetFillColor(0);
        canvas_dw->SetFrameFillColor(0);
        canvas_dw->SetFillColor(0);

        canvas_up->SetLeftMargin(0.1);
        canvas_dw->SetLeftMargin(0.1);
        canvas_up->SetRightMargin(0.075);
        canvas_dw->SetRightMargin(0.075);
        canvas_dw->SetBottomMargin(0.4);
        canvas_dw->SetTopMargin(0.05);

        // draw top figure
        canvas_up->cd();

        if (plot_.is_Log) canvas_up->SetLogy(true);

        tl->SetBorderSize(0);
        tl->SetLineColor(0);
        tl->SetLineWidth(0);
        tl->SetTextFont(62);
        tl->SetTextSize(0.035);
        tl->SetFillStyle(0);

        // declare early for legend
        TGraphAsymmErrors* Data = new TGraphAsymmErrors();
        // data legend FIRST
        tl->AddEntry(Data, "Data 2012", "p");

        THStack* hist_stack = new THStack(out_name.str().data(), "");

        //for (int t_ = 0; t_ < hft_trees.size(); t_++) {
        for (int i = 0; i < tree_ordering.size(); i++) {
            int t_ = tree_ordering[i];

            if (hft_trees[t_].tree_type == k_tree_mc) {
                hist_central[t_]->SetMarkerStyle(1);
                hist_central[t_]->SetLineColor(kBlack);
                hist_central[t_]->SetLineWidth(2);
                hist_central[t_]->SetFillColor(hft_trees[t_].color);
                hist_stack->Add(hist_central[t_]);
            }
        }

        for (int i = 0; i < tree_ordering.size(); i++) {
            int t_ = tree_ordering[i];

            if (hft_trees[t_].tree_type == k_tree_fakes) {
                hist_central[t_]->SetMarkerStyle(1);
                hist_central[t_]->SetLineColor(kBlack);
                hist_central[t_]->SetLineWidth(2);
                hist_central[t_]->SetFillColor(hft_trees[t_].color);
                hist_stack->Add(hist_central[t_]);
            }
        }
        //}

        // Signal
        for (int t_ = 0; t_ < hft_trees.size(); t_++) {
            if (hft_trees[t_].tree_type == k_tree_signal) {
                hist_central[t_]->SetMarkerStyle(1);
                hist_central[t_]->SetLineColor(hft_trees[t_].color);
                hist_central[t_]->SetLineWidth(3);
                hist_central[t_]->SetLineStyle(2);
            }
        }

        // MC Legend
        for (int t_ = 0; t_ < hft_trees.size(); t_++) {
            if (hft_trees[t_].tree_type == k_tree_mc) tl->AddEntry(hist_central[t_], hft_trees[t_].displayname.data(), "f");
        }

        for (int t_ = 0; t_ < hft_trees.size(); t_++) {
            if (hft_trees[t_].tree_type == k_tree_fakes) tl->AddEntry(hist_central[t_], hft_trees[t_].displayname.data(), "f");
        }

        // Signal legend
        for (int t_ = 0; t_ < hft_trees.size(); t_++) {
            if (hft_trees[t_].tree_type == k_tree_signal) {
                tl->AddEntry(hist_central[t_], hft_trees[t_].displayname.data(), "l");
            }

        }

        hist_stack->Draw("HIST");
        hist_stack->GetXaxis()->SetTitle(plot_.x_label.data());
        hist_stack->GetYaxis()->SetTitle(plot_.y_label.data());
        hist_stack->GetXaxis()->SetTitleFont(42);
        hist_stack->GetYaxis()->SetTitleFont(42);
        hist_stack->GetXaxis()->SetLabelFont(42);
        hist_stack->GetYaxis()->SetLabelFont(42);
        hist_stack->GetXaxis()->SetLabelOffset(999);
        hist_stack->GetXaxis()->SetTitleOffset(999);
        hist_stack->SetMinimum(plot_.y_range_min);
        hist_stack->SetMaximum(plot_.y_range_max);
        hist_stack->GetXaxis()->SetLabelSize(0.046);
        hist_stack->GetYaxis()->SetLabelSize(0.046);
        hist_stack->GetXaxis()->SetTitleSize(0.13);
        hist_stack->GetYaxis()->SetTitleSize(0.048);

        TGraphAsymmErrors* asym_errors = new TGraphAsymmErrors(hist_central[0]->GetNbinsX(), x_bin_center, nominal_yield, x_bin_down, x_bin_up, y_bin_down, y_bin_up);
        gStyle->SetHatchesSpacing(0.9);
        gStyle->SetHatchesLineWidth(1);
        asym_errors->SetFillStyle(3354); // 3004
        asym_errors->SetFillColor(kGray + 3);
        tl->AddEntry(asym_errors, "Bkg. Uncert.", "f");
        asym_errors->Draw("option same 02");

        // Decoration
        char annoyingLabel1[100] = "#bf{#it{ATLAS}} Internal", annoyingLabel2[100] = "#scale[0.6]{#int} L dt = 20.3 fb^{-1} #sqrt{s} = 8 TeV";
        myText(0.14, 0.82, kBlack, annoyingLabel1);
        myText(0.14, 0.74, kBlack, annoyingLabel2);
        string var_name = plot_.x_label;
        size_t find_g = var_name.find("[GeV]");
        if (find_g != string::npos) {
            var_name = var_name.substr(0, find_g);
        }
        string title_ = plot_.signal_region;
        if (plot_.request_cuts.compare("") != 0) {
            title_ += string("_") + plot_.request_cuts;
        }
        title_ += string(": ") + var_name;

        myText(0.14, 0.940, kBlack, title_.c_str());

        tl->Draw();


        for (int t_ = 0; t_ < hft_trees.size(); t_++) {
            if (hft_trees[t_].tree_type == k_tree_data) {
                convertErrorsToPoisson(hist_central[t_], Data);
                break;
            }
        }

        // Signal
        for (int t_ = 0; t_ < hft_trees.size(); t_++) {
            if (hft_trees[t_].tree_type == k_tree_signal) {
                hist_central[t_]->Draw("same HIST");
            }
        }

        Data->SetLineWidth(2);
        Data->SetMarkerStyle(8);
        Data->SetMarkerSize(1.5);
        Data->SetLineColor(1);
        Data->Draw("option same pz");

        // draw bottom figure
        // draw bottom figure
        // draw bottom figure

        canvas_dw->cd();

        TLine line_ = TLine(plot_.x_range_min, 1.0, plot_.x_range_max, 1.0);
        TLine line_up = TLine(plot_.x_range_min, 1.5, plot_.x_range_max, 1.5);
        TLine line_down = TLine(plot_.x_range_min, 0.5, plot_.x_range_max, 0.5);

        line_up.SetLineStyle(3);
        line_down.SetLineStyle(3);

        line_.SetLineColor(kRed);
        line_.SetLineStyle(2);

        //
        double* nominal_one = new double[hist_central[0]->GetNbinsX()];
        for (int j = 0; j < hist_central[0]->GetNbinsX(); j++) nominal_one[j] = 1.0;

        double* value_zero = new double[hist_central[0]->GetNbinsX()];
        for (int j = 0; j < hist_central[0]->GetNbinsX(); j++) value_zero[j] = 0.0;

        double* y_bin_down_ratio = new double[hist_central[0]->GetNbinsX()];
        for (int j = 0; j < hist_central[0]->GetNbinsX(); j++) {
            if (nominal_yield[j] != 0) {
                y_bin_down_ratio[j] = y_bin_down[j] / nominal_yield[j];
            }
            else {
                y_bin_down_ratio[j] = 0;
            }
        }

        double* y_bin_up_ratio = new double[hist_central[0]->GetNbinsX()];
        for (int j = 0; j < hist_central[0]->GetNbinsX(); j++) {
            if (nominal_yield[j] != 0) {
                y_bin_up_ratio[j] = y_bin_up[j] / nominal_yield[j];
            }
            else {
                y_bin_up_ratio[j] = 0;
            }
        }

        double* data_ratio = new double[hist_central[0]->GetNbinsX()];
        for (int j = 0; j < hist_central[0]->GetNbinsX(); j++) {
            if (nominal_yield[j] > 0 && hist_central[hft_trees.size() - 1]->GetBinContent(j + 1) > 0) {
                data_ratio[j] = hist_central[hft_trees.size() - 1]->GetBinContent(j + 1) / nominal_yield[j];
            }
            else {
                data_ratio[j] = 100.0;
            }
        }

        double* data_err_down_ratio = new double[hist_central[0]->GetNbinsX()];
        double* data_err_up_ratio = new double[hist_central[0]->GetNbinsX()];
        for (int j = 0; j < hist_central[0]->GetNbinsX(); j++) {
            if (nominal_yield[j] > 0) {

                data_err_down_ratio[j] = 1.0 - ((nominal_yield[j] - Data->GetErrorYlow(j)) / nominal_yield[j]);

                data_err_up_ratio[j] = ((nominal_yield[j] + Data->GetErrorYhigh(j)) / nominal_yield[j]) - 1.0;
            }
            else {
                data_err_down_ratio[j] = 0.0;
                data_err_up_ratio[j] = 0.0;
            }
        }

        TGraphAsymmErrors* asym_errors2 = new TGraphAsymmErrors(hist_central[0]->GetNbinsX(), x_bin_center, nominal_one, x_bin_down, x_bin_up, y_bin_down_ratio, y_bin_up_ratio);

        asym_errors2->SetFillStyle(3354); // 3004
        asym_errors2->SetFillColor(kGray + 3);

        TH1F* hist_ratio = (TH1F*)hist_central[0]->Clone();
        TGraphAsymmErrors* Data2 = new TGraphAsymmErrors(hist_central[0]->GetNbinsX(), x_bin_center, data_ratio, value_zero, value_zero, data_err_down_ratio, data_err_up_ratio);

        hist_ratio->SetName("ratio_");
        hist_ratio->SetMaximum(2.0);
        hist_ratio->SetMinimum(0.0);
        hist_ratio->GetXaxis()->SetTitleFont(42);
        hist_ratio->GetXaxis()->SetTitleFont(42);
        hist_ratio->GetYaxis()->SetTitleFont(42);
        hist_ratio->GetXaxis()->SetLabelFont(42);
        hist_ratio->GetYaxis()->SetLabelFont(42);
        hist_ratio->GetXaxis()->SetTitle(plot_.x_label.data());
        hist_ratio->GetYaxis()->SetTitle("Data / SM");
        hist_ratio->GetXaxis()->SetTickLength(0.06);
        hist_ratio->GetXaxis()->SetLabelSize(0.15 * font_size_dw);
        hist_ratio->GetYaxis()->SetLabelSize(0.15 * font_size_dw);
        hist_ratio->GetXaxis()->SetTitleSize(0.19 * font_size_dw);
        hist_ratio->GetYaxis()->SetTitleSize(0.16 * font_size_dw);
        hist_ratio->GetYaxis()->SetTitleOffset(0.40);
        hist_ratio->GetYaxis()->SetNdivisions(4, false);
        hist_ratio->Draw("AXIS");

        asym_errors2->Draw("option same 02");

        Data2->SetTitle("");
        Data2->SetLineWidth(2);
        Data2->SetMarkerStyle(8);
        Data2->SetMarkerSize(1.5);
        Data2->SetLineColor(1);
        Data2->Draw("option same 0pz");

        line_.Draw("same");
        line_up.Draw("same");
        line_down.Draw("same");












        stringstream out_file_root;
        out_file_root << DIR_1 << plot_.long_name << ".root";

        tc->Print(out_file.str().data());
        tc->Print(out_file_image.str().data());
        // tc->SaveAs(out_file_root.str().data());
        tc->Close();

        /// /// /// /// ///
        /// /// /// /// ///
        /// /// /// /// ///
        /// /// /// /// ///
        /// /// /// /// ///

        /// /// /// /// ///
        /// /// /// /// ///
        /// /// /// /// ///
        /// /// /// /// ///
        /// /// /// /// ///


        delete asym_errors;
        delete asym_errors2;
        delete Data;
        delete Data2;
        delete hist_stack;

        delete[] data_ratio;
        delete[] data_err_up_ratio;
        delete[] data_err_down_ratio;

        delete[] nominal_yield;
        delete[] x_bin_center;
        delete[] x_bin_down;
        delete[] x_bin_up;
        delete[] y_bin_down;
        delete[] y_bin_up;

        delete[] nominal_one;
        delete[] value_zero;
        delete[] y_bin_down_ratio;
        delete[] y_bin_up_ratio;

        for (int t_ = 0; t_ < hft_trees.size(); t_++) {
            delete hist_central[t_];
        }
        delete[] hist_central;

        for (int b_ = 0; b_ < num_bins; b_++) {
            for (int t_ = 0; t_ < hft_trees.size(); t_++) {
                delete[] systematic_up[b_][t_];
            }
            delete[] systematic_up[b_];
        }
        delete[] systematic_up;

        for (int b_ = 0; b_ < num_bins; b_++) {
            for (int t_ = 0; t_ < hft_trees.size(); t_++) {
                delete[] systematic_down[b_][t_];
            }
            delete[] systematic_down[b_];
        }
        delete[] systematic_down;


        delete tl;
        delete tc;

        ttemp->Clear();
    }

    delete ttemp;

    cout << "Done." << endl;
    return 0;
}






// // // VR_Etmiss_Zmm
// // // Plot: MT2
// // plot_.root_get_member = "MT2";
// // plot_.root_get_factor = "";
// // plot_.signal_region = "VR_Etmiss_Zmm";
// // plot_.request_cuts = "";
// // plot_.blind_data = "";
// // plot_.x_label = "m_{T2} [GeV]";
// // plot_.y_label = "Events / 10 GeV";
// // plot_.x_bin_width = 10.0;
// // plot_.x_range_min = 0.0;
// // plot_.x_range_max = 300.0;
// // plot_.y_range_min = 1.0e-1;
// // plot_.y_range_max = 1.0e8;
// // plot_.long_name = "VR_Etmiss_Zmm_MT2";
// // vect_plot.push_back(plot_);