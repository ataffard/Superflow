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
#include "TStyle.h"

//#include "Cintex/Cintex.h"

using namespace std;

#define DEBUG_NTC true

#define DIR_0 "/gdata/atlas/suneetu/Documents/LFV_Higgs2014/output/R3_August_22_patch/Processed/"
#define DIR_1 "/gdata/atlas/suneetu/Documents/LFV_Higgs2014/S9_Superflow/Superflow/plots/"

#define BG_FILE "LFVHBKG8TeV.root"
#define DATA_FILE "LFVHBKG8TeV.root"

const string m_central = "CENTRAL";

const int n_hists = 1;
string hists_entries[] = {};

struct plot {
	string root_get_member;
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
	int color;
	vector<int> systematic_groups;
};

vector<processed_tree> defineTrees() {
	vector<processed_tree> hft_trees;

	processed_tree tree_;

	vector<int> grouplist;
	grouplist.clear();

    
	// // FAKE trees
	// tree_.treename = "Fakes";
	// tree_.displayname = "Fakes";
	// tree_.tree_type = k_tree_fakes;
	// tree_.blinding_cuts = "";
	// tree_.special_cuts = "";
	// tree_.color = kGray;
	// grouplist.clear();
	// grouplist.push_back(1);// Fakes: group 1
	// tree_.systematic_groups = grouplist;
	// hft_trees.push_back(tree_);

	// monte-carlo trees
	// these apply to the group
	tree_.tree_type = k_tree_mc;
	tree_.blinding_cuts = "";
	tree_.special_cuts = "";

	// ZV
	tree_.treename = "ZV";
	tree_.displayname = "ZV";
	tree_.color = kSpring + 1;
	grouplist.clear();
	grouplist.push_back(0);
	grouplist.push_back(3);// XS
	tree_.systematic_groups = grouplist;
	hft_trees.push_back(tree_);

	// WW
	tree_.treename = "WW";
	tree_.displayname = "WW";
	tree_.color = kAzure + 4;
	grouplist.clear();
	grouplist.push_back(0);
	grouplist.push_back(4);// XS
	tree_.systematic_groups = grouplist;
	hft_trees.push_back(tree_);

	// Top
	tree_.treename = "Top";
	tree_.displayname = "tt+Wt";
	tree_.color = kRed + 1;
	grouplist.clear();
	grouplist.push_back(0);
	grouplist.push_back(5);// XS
	tree_.systematic_groups = grouplist;
	hft_trees.push_back(tree_);

	// Zjets
	tree_.treename = "Zjets";
	tree_.displayname = "Z+jets";
	tree_.color = kOrange;
	grouplist.clear();
	grouplist.push_back(0);
	grouplist.push_back(6);// XS
	tree_.systematic_groups = grouplist;
	hft_trees.push_back(tree_);

	// Higgs
	tree_.treename = "Higgs";
	tree_.displayname = "Higgs";
	tree_.color = kYellow - 9;
	grouplist.clear();
	grouplist.push_back(0);
	grouplist.push_back(7);// XS
	tree_.systematic_groups = grouplist;
	hft_trees.push_back(tree_);

	// DATA trees
	tree_.treename = "Data";
	tree_.displayname = "Data";
	tree_.tree_type = k_tree_data;
	tree_.blinding_cuts = "";//"(mlj>90000. || mljj>120000.)";
	tree_.special_cuts = "";
	tree_.color = kBlack;
	grouplist.clear();// no systematics on data
	tree_.systematic_groups = grouplist;
	hft_trees.push_back(tree_);

	return hft_trees;
}

vector<hft_systematic> defineSystematics() {
	vector<hft_systematic> hft_syst;

	hft_systematic syst;

	syst.adhoc_tcut2 = "";

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

    // WEIGHT SYSTEMATICS
    // these apply to all
    syst.sys_type = k_sys_weight;
    syst.up_name = "UP";
    syst.down_name = "DOWN";
    syst.preface = "syst_";
    syst.adhoc_tcut = "";
    syst.binding_group = 0;

    // BJET
    syst.basename = "BJET";
    hft_syst.push_back(syst);

    // CJET
    syst.basename = "CJET";
    hft_syst.push_back(syst);

    // BMISTAG
    syst.basename = "BMISTAG";
    hft_syst.push_back(syst);

    // ETRIGREW
    syst.basename = "ETRIGREW";
    hft_syst.push_back(syst);

    // MTRIGREW
    syst.basename = "MTRIGREW";
    hft_syst.push_back(syst);

    // ESF
    syst.basename = "ESF";
    hft_syst.push_back(syst);

    // MEFF
    syst.basename = "MEFF";
    hft_syst.push_back(syst);

    // PILEUP
    syst.basename = "PILEUP";
    hft_syst.push_back(syst);

    // MEFF
    syst.basename = "MEFF";
    hft_syst.push_back(syst);

    // SPLIT XS
    // these apply to all XS
    syst.sys_type = k_sys_weight;
    syst.up_name = "UP";
    syst.down_name = "DOWN";
    syst.preface = "syst_";
    syst.adhoc_tcut = "";

    // XS
    syst.basename = "XS";
    syst.binding_group = 3;
    hft_syst.push_back(syst);

    syst.basename = "XS";
    syst.binding_group = 4;
    hft_syst.push_back(syst);

    syst.basename = "XS";
    syst.binding_group = 5;
    hft_syst.push_back(syst);

    syst.basename = "XS";
    syst.binding_group = 7;
    hft_syst.push_back(syst);

    // XS on Z+jets
    syst.sys_type = k_sys_adhoc;
    syst.basename = "XS";
    syst.binding_group = 6;
    syst.preface = "syst_";
    syst.adhoc_tcut = "1.0014";
    syst.adhoc_tcut2 = "0.9986"; // note
    hft_syst.push_back(syst);

    /*
	// XS
	syst.basename = "XS";
	syst.binding_group = 3;
	hft_syst.push_back(syst);

	syst.basename = "XS";
	syst.binding_group = 4;
	hft_syst.push_back(syst);

	syst.basename = "XS";
	syst.binding_group = 5;
	hft_syst.push_back(syst);

	syst.basename = "XS";
	syst.binding_group = 6;
	hft_syst.push_back(syst);

	syst.basename = "XS";
	syst.binding_group = 7;
	hft_syst.push_back(syst);
    */

    /*
	// SPLIT GEN
	// these apply to all GEN
	syst.sys_type = k_sys_weight;
	syst.up_name = "UP";
	syst.down_name = "DOWN";
	syst.preface = "syst_";
	syst.adhoc_tcut = "";

	//GEN
	syst.basename = "GEN";
	syst.binding_group = 3;
	hft_syst.push_back(syst);

	syst.basename = "GEN";
	syst.binding_group = 4;
	hft_syst.push_back(syst);

	syst.basename = "GEN";
	syst.binding_group = 5;
	hft_syst.push_back(syst);

	syst.basename = "GEN";
	syst.binding_group = 6;
	hft_syst.push_back(syst);

	syst.basename = "GEN";
	syst.binding_group = 7;
	hft_syst.push_back(syst);

	//FAKE SYSTEMATICS (object)
	//these apply to all
	syst.sys_type = k_sys_object;
	syst.up_name = "UP";
	syst.down_name = "DOWN";
	syst.preface = "";
	syst.adhoc_tcut = "";
	syst.binding_group = 1;

	// ELFR
	syst.basename = "ELFR";
	hft_syst.push_back(syst);

	// MUFR
	syst.basename = "MUFR";
	hft_syst.push_back(syst);

	// ELRE
	syst.basename = "ELRE";
	hft_syst.push_back(syst);

	// MURE
	syst.basename = "MURE";
	hft_syst.push_back(syst);

	// ELFRAC
	syst.basename = "ELFRAC";
	hft_syst.push_back(syst);

	// MUFRAC
	syst.basename = "MUFRAC";
	hft_syst.push_back(syst);

    */

	// DONE //
	return hft_syst;
}

vector<plot> definePlots() {
	vector<plot> vect_plot;
	plot plot_;

	plot_.is_Log = false;
	plot_.reweight = "";

	// Section: VR_Etmiss_Zmm
	// Plot: table
	plot_.root_get_member = "1";
	plot_.root_get_factor = "";
	plot_.signal_region = "CR_el_and_mu";
	plot_.request_cuts = "";
	plot_.blind_data = "";
	// plot setup
	plot_.x_label = "";
	plot_.y_label = "";
	plot_.x_bin_width = 1.0;
	plot_.x_range_min = 0.0;
	plot_.x_range_max = 1.0;
	plot_.y_range_min = 0.0;
	plot_.y_range_max = 1.0;
	vect_plot.push_back(plot_);

    // Section: VR_Etmiss_Zee
    // Plot: table
    plot_.root_get_member = "1";
    plot_.root_get_factor = "";
    plot_.signal_region = "CR_EM";
    plot_.request_cuts = "";
    plot_.blind_data = "";
    // plot setup
    plot_.x_label = "";
    plot_.y_label = "";
    plot_.x_bin_width = 1.0;
    plot_.x_range_min = 0.0;
    plot_.x_range_max = 1.0;
    plot_.y_range_min = 0.0;
    plot_.y_range_max = 1.0;
    vect_plot.push_back(plot_);

    // Section: VR_Etmiss_Zmm
    // Plot: table
    plot_.root_get_member = "1";
    plot_.root_get_factor = "";
    plot_.signal_region = "CR_ME";
    plot_.request_cuts = "";
    plot_.blind_data = "";
    // plot setup
    plot_.x_label = "";
    plot_.y_label = "";
    plot_.x_bin_width = 1.0;
    plot_.x_range_min = 0.0;
    plot_.x_range_max = 1.0;
    plot_.y_range_min = 0.0;
    plot_.y_range_max = 1.0;
    vect_plot.push_back(plot_);

    // Section: VR_Etmiss_Zee
    // Plot: table
    plot_.root_get_member = "1";
    plot_.root_get_factor = "";
    plot_.signal_region = "CR_high_mcoll_EM";
    plot_.request_cuts = "";
    plot_.blind_data = "";
    // plot setup
    plot_.x_label = "";
    plot_.y_label = "";
    plot_.x_bin_width = 1.0;
    plot_.x_range_min = 0.0;
    plot_.x_range_max = 1.0;
    plot_.y_range_min = 0.0;
    plot_.y_range_max = 1.0;
    vect_plot.push_back(plot_);

    // Section: VR_Etmiss_Zmm
    // Plot: table
    plot_.root_get_member = "1";
    plot_.root_get_factor = "";
    plot_.signal_region = "CR_low_mcoll_EM";
    plot_.request_cuts = "";
    plot_.blind_data = "";
    // plot setup
    plot_.x_label = "";
    plot_.y_label = "";
    plot_.x_bin_width = 1.0;
    plot_.x_range_min = 0.0;
    plot_.x_range_max = 1.0;
    plot_.y_range_min = 0.0;
    plot_.y_range_max = 1.0;
    vect_plot.push_back(plot_);

    // Section: VR_Etmiss_Zee
    // Plot: table
    plot_.root_get_member = "1";
    plot_.root_get_factor = "";
    plot_.signal_region = "CR_high_mcoll_ME";
    plot_.request_cuts = "";
    plot_.blind_data = "";
    // plot setup
    plot_.x_label = "";
    plot_.y_label = "";
    plot_.x_bin_width = 1.0;
    plot_.x_range_min = 0.0;
    plot_.x_range_max = 1.0;
    plot_.y_range_min = 0.0;
    plot_.y_range_max = 1.0;
    vect_plot.push_back(plot_);

    // Section: VR_Etmiss_Zee
    // Plot: table
    plot_.root_get_member = "1";
    plot_.root_get_factor = "";
    plot_.signal_region = "CR_low_mcoll_ME";
    plot_.request_cuts = "";
    plot_.blind_data = "";
    // plot setup
    plot_.x_label = "";
    plot_.y_label = "";
    plot_.x_bin_width = 1.0;
    plot_.x_range_min = 0.0;
    plot_.x_range_max = 1.0;
    plot_.y_range_min = 0.0;
    plot_.y_range_max = 1.0;
    vect_plot.push_back(plot_);
	
	return vect_plot;
}

int main(int argc, char* argv[]) {

	vector<plot> request_plots = definePlots();
	vector<processed_tree> hft_trees = defineTrees();
	vector<hft_systematic> hft_syst = defineSystematics();

	// define the all-important TCut strings
	map <string, string> tcuts;
	stringstream tcuts_join;

    //VR_Etmiss_Zmm
    tcuts_join << "isMuMu && lept1Pt>35000 && lept2Pt>18000 && !(ll_M>(91200.-10000.) && ll_M<(91200.+10000.)) && nCentralBJets==0 && nForwardJets==0 && nCentralLightJets<2"; //  && MT2>40000
    tcuts["VR_Etmiss_Zmm"] = tcuts_join.str(); tcuts_join.str("");

    //VR_Etmiss_Zee
    tcuts_join << "isElEl && lept1Pt>35000 && lept2Pt>18000 && !(ll_M>(91200.-10000.) && ll_M<(91200.+10000.)) && nCentralBJets==0 && nForwardJets==0 && nCentralLightJets<2";
    tcuts["VR_Etmiss_Zee"] = tcuts_join.str(); tcuts_join.str("");

    //VR_LFV_Top
    tcuts_join << "(isEM || isME) && lept1Pt>30000 && lept2Pt>18000 && nCentralBJets>0 && nForwardJets==0";
    tcuts["VR_LFV_Top"] = tcuts_join.str(); tcuts_join.str("");

    //CR_el_and_mu
    tcuts_join << "isElMu && lept1Pt>35000. && lept2Pt>18000.";
    tcuts["CR_el_and_mu"] = tcuts_join.str(); tcuts_join.str("");

    //CR_EM
    tcuts_join << "isEM && lept1Pt>35000. && lept2Pt>18000.";
    tcuts["CR_EM"] = tcuts_join.str(); tcuts_join.str("");

    //CR_ME
    tcuts_join << "isME && lept1Pt>35000. && lept2Pt>18000.";
    tcuts["CR_ME"] = tcuts_join.str(); tcuts_join.str("");

    //CR_high_mcoll_EM                                              
    tcuts_join << "isEM && lept1Pt>35000. && lept2Pt>18000. && mcoll>150000. && nCentralBJets==0 && nForwardJets==0 ";
    tcuts["CR_high_mcoll_EM"] = tcuts_join.str(); tcuts_join.str("");

    //CR_high_mcoll_EM
    tcuts_join << "isEM && lept1Pt>35000. && lept2Pt>18000. && mcoll<100000. && nCentralBJets==0 && nForwardJets==0 ";
    tcuts["CR_low_mcoll_EM"] = tcuts_join.str(); tcuts_join.str("");

    //CR_high_mcoll_ME
    tcuts_join << "isME && lept1Pt>35000. && lept2Pt>18000. && mcoll>150000. && nCentralBJets==0 && nForwardJets==0 ";
    tcuts["CR_high_mcoll_ME"] = tcuts_join.str(); tcuts_join.str("");

    //CR_low_mcoll_ME
    tcuts_join << "isME && lept1Pt>35000. && lept2Pt>18000. && mcoll<100000. && nCentralBJets==0 && nForwardJets==0 ";
    tcuts["CR_low_mcoll_ME"] = tcuts_join.str(); tcuts_join.str("");

	tcuts["none"] = "1";

	TCanvas* tc = new TCanvas("canvas_1", "", 768, 768);

	vector < vector< vector<double> > > table_columns;
	vector<string> sr;

	stringstream filename_stream;
	filename_stream << DIR_0 << BG_FILE;

	string filename_ = filename_stream.str();

	for (int r_ = 0; r_ < request_plots.size(); r_++) {
		plot plot_ = request_plots[r_];
		sr.push_back(plot_.signal_region);
		vector< vector<double> > table_col_;// column is per process

		cout << setprecision(4);
		cout << std::fixed;
		cout << "+++++++++++++ " << plot_.signal_region << " +++++++++++++" << endl;

		// Complete/finalize TCut strings.
		stringstream signal_region_stream;
		signal_region_stream << "( ";
		signal_region_stream << tcuts[plot_.signal_region];
		if (plot_.request_cuts.compare("") != 0) signal_region_stream << " && " << plot_.request_cuts;
		signal_region_stream << " )";// Don't forget the close parentheses

		string signal_tcut = signal_region_stream.str();

		TH1F** hist_central = new TH1F*[hft_trees.size()];
		for (int t_ = 0; t_ < hft_trees.size(); t_++) {
			stringstream hist_;
			hist_ << "tree_" << plot_.signal_region << "_" << hft_trees[t_].treename;
			hist_central[t_] = new TH1F(hist_.str().data(), hist_.str().data(), 1, 0.0, 1.0);
			hist_central[t_]->Sumw2();
		}

		double** systematic_up = new double*[hft_trees.size()];
		for (int t_ = 0; t_ < hft_trees.size(); t_++) {
			systematic_up[t_] = new double[hft_syst.size()];
			for (int s_ = 0; s_ < hft_syst.size(); s_++) {
				systematic_up[t_][s_] = 0;
			}
		}

		double** systematic_down = new double*[hft_trees.size()];
		for (int t_ = 0; t_ < hft_trees.size(); t_++) {
			systematic_down[t_] = new double[hft_syst.size()];
			for (int s_ = 0; s_ < hft_syst.size(); s_++) {
				stringstream hist_;
				hist_ << "tree_" << plot_.signal_region << "_" << hft_trees[t_].treename << "_" << hft_syst[s_].basename << "_DOWN";
				systematic_down[t_][s_] = 0;
			}
		}

		vector<TH1F> total_ststematics;
		vector<double> tree_total;
		vector<double> tree_stat_error;

		// define central yields table
		for (int t_ = 0; t_ < hft_trees.size(); t_++) {
			stringstream chain_name;// select CENTRAL
			chain_name << hft_trees[t_].treename << "_" << m_central;

			TChain* chain_;
			chain_ = new TChain(chain_name.str().data());
			chain_->Add(filename_.data());

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

			stringstream draw_;
			draw_ << "1 >> +" << hist_central[t_]->GetName();

			chain_->Draw(draw_.str().data(), (sel + blind) * weight * tc_reweight);
			double error_ = 0;
			double integral_ = hist_central[t_]->IntegralAndError(0, hist_central[t_]->GetNbinsX() + 1, error_);

			tree_total.push_back(integral_);
			tree_stat_error.push_back(error_);

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
			// /}

			delete chain_;
		}

		// calculate all systematics
		for (int t_ = 0; t_ < hft_trees.size(); t_++) {
			if (!DEBUG_NTC) {
				cout << "<< " << hft_trees[t_].treename;
			}
			int lcount = 0;

			for (int s_ = 0; s_ < hft_syst.size(); s_++) {

				bool reweight_flag = false;

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
					reweight_flag = true;
				}
				else if (hft_syst[s_].sys_type == k_sys_adhoc) {
					chain_name_up << hft_trees[t_].treename << "_" << m_central;
					chain_name_down << hft_trees[t_].treename << "_" << m_central;
					reweight_flag = true;
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
				TCut blind = TCut(hft_trees[t_].blinding_cuts.data());

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
					syst_var_down = TCut(hft_syst[s_].adhoc_tcut.data());
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

				// temp histogram to get the variation
				TH1F* temp_hist_up = new TH1F("temp_hist_up", "temp_hist_up", 1, 0.0, 1.0);
				TH1F* temp_hist_down = new TH1F("temp_hist_down", "temp_hist_down", 1, 0.0, 1.0);

				stringstream draw_up;
				draw_up << "1.0 >> +" << "temp_hist_up";
				stringstream draw_down;
				draw_down << "1.0 >> +" << "temp_hist_down";

				if (hft_syst[s_].sys_type == k_sys_weight) {
					chain_up->Draw(draw_up.str().data(), sel * weight * syst_var_up * tc_reweight);
					chain_down->Draw(draw_down.str().data(), sel * weight * syst_var_down * tc_reweight);
				}
				else if (hft_syst[s_].sys_type == k_sys_adhoc) {
					chain_up->Draw(draw_up.str().data(), sel * weight * syst_var_up * tc_reweight);
					chain_down->Draw(draw_down.str().data(), sel * weight * syst_var_down * tc_reweight);
				}
				else {
					chain_up->Draw(draw_up.str().data(), sel * weight * tc_reweight);
					chain_down->Draw(draw_down.str().data(), sel * weight * tc_reweight);
				}

				double int_up = temp_hist_up->Integral(0, -1) - tree_total[t_];
				double int_down = temp_hist_down->Integral(0, -1) - tree_total[t_];

				// sign and sidedness adjustments
				if (hft_syst[s_].sys_type == k_sys_one_sided_object) {
					if (int_up > 0)	int_down = 0.0;
					else int_up = 0.0;
				}
				else if (hft_syst[s_].sys_type == k_sys_adhoc) {
					int_up = abs(int_up);
					int_down = -abs(int_down);
				}

				systematic_up[t_][s_] = int_up;
				systematic_down[t_][s_] = int_down;

				if (DEBUG_NTC) {
					// cout << hft_trees[t_].treename << " > " << hft_syst[s_].basename << " ";
					// cout << tree_total[t_] << ",+" << systematic_up[t_][s_] / tree_total[t_] << ",-" << systematic_down[t_][s_] / tree_total[t_] << endl;
					cout << hft_trees[t_].treename << " > " << hft_syst[s_].basename << " ";
					cout << tree_total[t_] << ",+" << systematic_up[t_][s_] << ",-" << systematic_down[t_][s_] << endl;
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

			if (!DEBUG_NTC) {
				cout << endl;
			}
		}

		// per grouping table
		for (int t_ = 0; t_ < hft_trees.size(); t_++) {
			vector<double> single_entry(4, 0.0);

			double syst_up = 0;
			double syst_down = 0;

			for (int s_ = 0; s_ < hft_syst.size(); s_++) {
				double up_var = systematic_up[t_][s_];
				double down_var = systematic_down[t_][s_];

				//syst_up += up_var * up_var;
				//syst_down += down_var * down_var;

				if (up_var > 0) syst_up += up_var * up_var;
				else syst_down += up_var * up_var;

				if (down_var < 0) syst_down += down_var * down_var;
				else syst_up += down_var * down_var;
			}

			syst_up = sqrt(syst_up);
			syst_down = sqrt(syst_down);

			cout << hft_trees[t_].treename << " > " << tree_total[t_] << " \xB1" << tree_stat_error[t_];
			cout << " + " << syst_up << " - " << syst_down << endl;

			single_entry[0] = tree_total[t_];
			single_entry[1] = tree_stat_error[t_];
			single_entry[2] = syst_up;
			single_entry[3] = syst_down;

			table_col_.push_back(single_entry);
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

		double total_stat_up = 0;
		double total_stat_down = 0;

		// totals
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

					double up_var = systematic_up[t_][s_];
					double down_var = systematic_down[t_][s_];

					syst_up += up_var;
					syst_down += down_var;
				}

				//total_stat_up += syst_up * syst_up;
				//total_stat_down += syst_down * syst_down;

				if (syst_up > 0) {
					total_stat_up += syst_up * syst_up;
				}
				else {
					total_stat_down += syst_up * syst_up;
				}

				if (syst_down < 0) {
					total_stat_down += syst_down * syst_down;
				}
				else {
					total_stat_up += syst_down * syst_down;
				}
			}
		}

		total_stat_up = sqrt(total_stat_up);
		total_stat_down = sqrt(total_stat_down);

		double realtotal = 0;
		double staterror = 0;
		for (int t_ = 0; t_ < hft_trees.size(); t_++) {
			if (hft_trees[t_].tree_type == k_tree_data) continue;
			realtotal += tree_total[t_];
			staterror += tree_stat_error[t_] * tree_stat_error[t_];
		}
		staterror = sqrt(staterror);

		cout << "Total" << " > " << realtotal << " \xB1" << staterror;
		cout << " + " << total_stat_up << " - " << total_stat_down << endl;

		vector<double> total_entry(4, 0.0);

		total_entry[0] = realtotal;
		total_entry[1] = staterror;
		total_entry[2] = total_stat_up;
		total_entry[3] = total_stat_down;

		table_col_.push_back(total_entry);

		table_columns.push_back(table_col_);

		for (int t_ = 0; t_ < hft_trees.size(); t_++) {
			delete hist_central[t_];
		}
		delete[] hist_central;

		for (int t_ = 0; t_ < hft_trees.size(); t_++) delete[] systematic_up[t_];
		delete[] systematic_up;

		for (int t_ = 0; t_ < hft_trees.size(); t_++) delete[] systematic_down[t_];
		delete[] systematic_down;

		tc->Clear();
	}

	delete tc;

	//vector < vector< vector<double> > > table_columns;
	//vector<string> sr;
	//vector<string> process_column;

	// Print nice tables
	vector<string> process_column;
	for (int t_ = 0; t_ < hft_trees.size(); t_++) {
		process_column.push_back(hft_trees[t_].displayname);
	}
	process_column.push_back("Total");

	cout << setprecision(2);
	cout << std::fixed;
	cout << endl;

	cout << "\t";
	for (int sr_ = 0; sr_ < sr.size(); sr_++) {
		cout << sr[sr_] << "  " << "\t\t\t";
	}
	cout << endl;

	for (int pr_ = 0; pr_ < process_column.size(); pr_++) {
		if (process_column[pr_].compare("Data") == 0) continue;
		if (process_column[pr_].compare("Total") == 0) 	cout << "---" << endl;
		cout << process_column[pr_] << "\t";
		for (int tc_ = 0; tc_ < table_columns.size(); tc_++) {
			cout << table_columns[tc_][pr_][0] << "\xB1" << table_columns[tc_][pr_][1];
			cout << "^+" << table_columns[tc_][pr_][2] << "_-" << table_columns[tc_][pr_][3] << "    \t";
		}
		cout << endl;
	}

	cout << "---" << endl;

	for (int pr_ = 0; pr_ < process_column.size(); pr_++) {
		if (process_column[pr_].compare("Data") != 0) continue;
		cout << process_column[pr_] << "\t";
		for (int tc_ = 0; tc_ < table_columns.size(); tc_++) {
			cout << table_columns[tc_][pr_][0] << "    \t\t\t";
		}
		cout << endl;
	}

	/////////////////////////
	cout << endl << endl;

	/////////////////////////
	cout << "Symmetrized" << endl;

	cout << "\t\t";
	for (int sr_ = 0; sr_ < sr.size(); sr_++) {
		cout << sr[sr_] << "  " << "\t";
	}
	cout << endl;

	for (int pr_ = 0; pr_ < process_column.size(); pr_++) {
		if (process_column[pr_].compare("Data") == 0) continue;
		if (process_column[pr_].compare("Total") == 0) 	cout << "---" << endl;
		cout << process_column[pr_] << "\t\t";
		for (int tc_ = 0; tc_ < table_columns.size(); tc_++) {
			cout << table_columns[tc_][pr_][0] << "\xB1";
			cout << sqrt(table_columns[tc_][pr_][1] * table_columns[tc_][pr_][1]
				+ pow(0.5 * (abs(table_columns[tc_][pr_][2]) + abs(table_columns[tc_][pr_][3])), 2.0)) << "\t";
		}
		cout << endl;
	}

	cout << "---" << endl;

	for (int pr_ = 0; pr_ < process_column.size(); pr_++) {
		if (process_column[pr_].compare("Data") != 0) continue;
		cout << process_column[pr_] << "\t\t";
		for (int tc_ = 0; tc_ < table_columns.size(); tc_++) {
			cout << table_columns[tc_][pr_][0] << "   \t";
		}
		cout << endl;
	}

	cout << "Done." << endl;
	return 0;
}