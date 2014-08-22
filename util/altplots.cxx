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
#include "TDirectory.h"
#include "TMath.h"

#include "TTree.h"
#include "TChain.h"
#include "TFile.h"
#include "TCut.h"

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
#include "Cintex/Cintex.h"

using namespace std;

#define DIR_0 "/gdata/atlas/suneetu/Documents/LFV_Higgs2014/output/R1_August_16/Processed/"
#define DIR_1 "/gdata/atlas/suneetu/Documents/LFV_Higgs2014/S2_Superflow/Superflow/plots/"

#define BG_FILE "LFVHBKG8TeV.root"
#define DATA_FILE "LFVHBKG8TeV.root"

const int n_process = 6;// added one for "Data"
const int n_mc_process = n_process - 2;
string process[] = { "WW", "Top", "Zjets", "ZV", "Higgs", "Data" };
//string process[] = { "Zjets", "Top", "Zjets", "ZV", "Higgs", "fake", "Data" };
int process_INSERT_ORDER[] = { 3, 1, 2, 0, 4, 5, 6 };
string process_name[] = { "WW", "tt+Wt", "Z+jets", "ZV", "Higgs", "Data" };
//string process_name[] = { "Z+jets", "tt+Wt", "Z+jets", "ZV", "Higgs", "Fake leptons", "Data" };
int stack_colors[] = { kAzure + 4, kRed + 1, kOrange, kSpring + 1, kYellow - 9, kMagenta - 7 };

const int n_systematics = 4;//11
string systematics[] = {
    "CENTRAL", "EESZ", "EER", "EESLOW",
    "EESMAT", "EESPS", "ID", "JES",
    "MS", "SCALEST", "TES" };

const int n_one_sided_syst = 3;
string one_sided_syst[] = { "CENTRAL", "JER", "RESOST" };

const int n_fakesyst = 3;// 7
string fake_syst[] = { "CENTRAL", "EESZ", "EER", "EESLOW", "EESMAT", "EESPS", "ID", "JES" };

const int n_weight_syst = 3;// 9
string weight_syst[] = { "BJET", "CJET", "BMISTAG", "ETRIGREW", "MTRIGREW", "ESF", "MEFF", "PILEUP" };

const int n_hists = 1;
string hists_entries[] = {};

const double pi = 3.141592653589793238462;

void myText(Double_t x, Double_t y, Color_t color, char *text) {

	Double_t tsize = 0.055;
	TLatex l; //l.SetTextAlign(12); 
	l.SetTextSize(tsize);
	l.SetNDC();
	l.SetTextColor(color);
	l.DrawLatex(x, y, text);
}

void convertErrorsToPoisson(TH1* histo, TGraphAsymmErrors* graph)
{
	// Needed variables
	double value = 0;
	double error_poisson_up = 0;
	double error_poisson_down = 0;
	double alpha = 0.158655, beta = 0.158655; // 68%

	// loop over bins and overwrite values
	for (int i = 1; i < histo->GetNbinsX() + 1; i++){
		value = histo->GetBinContent(i);
		if (value != 0) {
			error_poisson_up = 0.5*TMath::ChisquareQuantile(1 - beta, 2 * (value + 1)) - value;
			error_poisson_down = value - 0.5*TMath::ChisquareQuantile(alpha, 2 * value);
			graph->SetPoint(i - 1, histo->GetBinCenter(i), value);
			graph->SetPointError(i - 1, 0., 0., error_poisson_down, error_poisson_up);
		}
		else{
			graph->SetPoint(i - 1, histo->GetBinCenter(i), -10);
			graph->SetPointError(i - 1, 0., 0., 0., 0.);
		}
	}
}

struct plot {
	string root_get_member;
	string root_get_member_name;
	string root_get_factor;
	string signal_region;
	int n_jets;
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
};

vector<plot> definePlots() {
	vector<plot> vect_plot;
	plot plot_;

	plot_.is_Log = false;

	// chargeflip study
	plot_.root_get_member_name = "";
	plot_.reweight = "";

	plot_.is_Log = true;

    // VR_Etmiss_Zmm
    // VR_Etmiss_Zmm
    // VR_Etmiss_Zmm

    // VR_Etmiss_Zmm
    // Plot: metCorr
    plot_.root_get_member = "metCorr";
    plot_.root_get_factor = " / 1000.0";
    plot_.signal_region = "VR_Etmiss_Zmm";
    plot_.request_cuts = "";
    plot_.blind_data = "";
    plot_.x_label = "E_{T}^{miss-corr} [GeV]";
    plot_.y_label = "Events / 10 GeV";
    plot_.x_bin_width = 10.0;
    plot_.x_range_min = 0.0;
    plot_.x_range_max = 300.0;
    plot_.y_range_min = 1.0e-1;
    plot_.y_range_max = 1.0e9;
    plot_.long_name = "VR_Etmiss_Zmm_metCorr";
    vect_plot.push_back(plot_);

	return vect_plot;
}

int main(int argc, char* argv[]) {

	vector<plot> request_plots = definePlots();

	// define the all-important TCut strings
	map <string, string> tcuts;
	stringstream tcuts_join;

    //VR_Etmiss_Zmm
    tcuts_join << "isMuMu && lept1Pt>35000 && lept2Pt>18000 && !(ll_M>(91200.-10000.) && ll_M<(91200.+10000.)) && nCentralBJets==0 && nForwardJets==0 && nCentralLightJets<2";
    tcuts["VR_Etmiss_Zmm"] = tcuts_join.str(); tcuts_join.str("");// clear stringstream

    //VR_Etmiss_Zee
    tcuts_join << "isElEl && lept1Pt>35000 && lept2Pt>18000 && !(ll_M>(91200.-10000.) && ll_M<(91200.+10000.)) && nCentralBJets==0 && nForwardJets==0 && nCentralLightJets<2";
    tcuts["VR_Etmiss_Zee"] = tcuts_join.str(); tcuts_join.str("");//  && mlj<90000. 

    //VR_LFV_Top
    tcuts_join << "(isEM || isME) && lept1Pt>30000 && lept2Pt>18000 && nCentralBJets>0 && nForwardJets==0";
    tcuts["VR_LFV_Top"] = tcuts_join.str(); tcuts_join.str("");// clear stringstream

	tcuts["none"] = "";

	for (int r_ = 0; r_ < request_plots.size(); r_++) {
		plot plot_ = request_plots[r_];

		// Complete/finalize TCut strings.
		stringstream signal_region_stream;
		signal_region_stream << "( " << tcuts[plot_.signal_region];
		if (plot_.request_cuts.compare("") != 0) {
			if (plot_.signal_region.compare("none") == 0) {
				signal_region_stream << plot_.request_cuts;
			}
			else {
				signal_region_stream << " && " << plot_.request_cuts;
			}
		}
		signal_region_stream << " )";// Don't forget the close parenthesis

		string signal_tcut = signal_region_stream.str();

		cout << signal_tcut << endl;

		TH1F** new_hist = new TH1F*[n_process];
		for (int p_ = 0; p_ < n_process; p_++) {
			stringstream hist_;
			if (plot_.root_get_member_name.compare("") == 0) {
				hist_ << "h_" << plot_.signal_region << "_" << plot_.root_get_member << "_" << process[p_];
			}
			else {
				hist_ << "h_" << plot_.signal_region << "_" << plot_.root_get_member_name << "_" << process[p_];
			}
			// calculate the number of bins
			int n_bins = floor((plot_.x_range_max - plot_.x_range_min) / (plot_.x_bin_width) + 0.5);
			new_hist[p_] = new TH1F(hist_.str().data(), hist_.str().data(), n_bins, plot_.x_range_min, plot_.x_range_max);
			new_hist[p_]->Sumw2();
		}

		// WEIGHT DELTA UP
		double** weight_delta_up = new double*[n_weight_syst];
		for (int i = 0; i < n_weight_syst; i++){
			weight_delta_up[i] = new double[new_hist[0]->GetNbinsX()];
			for (int j = 0; j < new_hist[0]->GetNbinsX(); j++) weight_delta_up[i][j] = 0.0;
		}

		double* gen_weight_delta_up = new double[new_hist[0]->GetNbinsX()];
		for (int j = 0; j < new_hist[0]->GetNbinsX(); j++) gen_weight_delta_up[j] = 0.0;

		double* xs_weight_delta_up = new double[new_hist[0]->GetNbinsX()];
		for (int j = 0; j < new_hist[0]->GetNbinsX(); j++) xs_weight_delta_up[j] = 0.0;

		// WEIGHT DELTA DOWN
		double** weight_delta_down = new double*[n_weight_syst];
		for (int i = 0; i < n_weight_syst; i++){
			weight_delta_down[i] = new double[new_hist[0]->GetNbinsX()];
			for (int j = 0; j < new_hist[0]->GetNbinsX(); j++) weight_delta_down[i][j] = 0.0;
		}

		double* gen_weight_delta_down = new double[new_hist[0]->GetNbinsX()];
		for (int j = 0; j < new_hist[0]->GetNbinsX(); j++) gen_weight_delta_down[j] = 0.0;

		double* xs_weight_delta_down = new double[new_hist[0]->GetNbinsX()];
		for (int j = 0; j < new_hist[0]->GetNbinsX(); j++) xs_weight_delta_down[j] = 0.0;

		{
			int s_ = 0;// This section is for CENTRAL
			for (int p_ = 0; p_ < n_process; p_++) {//


				// open the file, and grab the right TChain
				stringstream filename_;
				filename_ << DIR_0 << BG_FILE;

				stringstream chain_name;
				chain_name << process[p_] << "_" << systematics[s_];

				TChain* chain_;
				chain_ = new TChain(chain_name.str().data());

				chain_->Add(filename_.str().data());

				if (chain_->IsZombie()) {
					cout << chain_name.str() << "not found." << endl;
					continue;
				}

				TCut sel = TCut("");// key TCut
				TCut weight = TCut("eventweight");// selects eventweight leaf
				stringstream signal_line;
				if (process[p_].compare("Data") == 0) {
					stringstream signal_region_stream;
					signal_region_stream << "( " << tcuts[plot_.signal_region];
					if (plot_.request_cuts.compare("") != 0) {
						if (plot_.signal_region.compare("none") == 0) {
							signal_region_stream << plot_.request_cuts;
						}
						else {
							signal_region_stream << " && " << plot_.request_cuts;
						}
					}
					if (plot_.signal_region.compare("none") != 0 && plot_.blind_data.compare("") != 0) {
						signal_region_stream << " && " << plot_.blind_data;
					}
					signal_region_stream << " )";// Don't forget the close parenthesis
					cout << "BLIND DATA" << endl;
					signal_line << signal_region_stream.str();
				}
				else {
					signal_line << signal_tcut.data();
				}
				sel += TCut(signal_line.str().data());// string.data()

				//weight.Print();
				//sel.Print();

				// START reweighting segment
				TCut tc_reweight = TCut("1.0");
				if (plot_.reweight.compare("") != 0) {
					tc_reweight = TCut(plot_.reweight.data());
					cout << "REWEIGHT!" << endl;
				}
				// END reweighting segment

				stringstream draw_;
				draw_ << plot_.root_get_member << plot_.root_get_factor << " >> +" << new_hist[p_]->GetName();

				//cout << "Draw: " << draw_.str() << endl;

				chain_->Draw(draw_.str().data(), sel * weight * tc_reweight);// reweight

				double integral_ = new_hist[p_]->Integral(0, -1);

				cout << new_hist[p_]->GetName() << ":" << integral_ << endl;

				delete chain_;
			}
		}


		// NOMINAL YIELD
		double* nominal_yield = new double[new_hist[0]->GetNbinsX()];
		for (int j = 0; j < new_hist[0]->GetNbinsX(); j++) nominal_yield[j] = 0.0;

		// calculate and store, per bin nominal_yield;
		TH1F* hist_prediction = (TH1F*)new_hist[0]->Clone();
		for (int p_ = 1; p_ < n_process - 2; p_++) {// add MC 
			hist_prediction->Add(new_hist[p_], 1.0);
		}
		for (int j = 1; j <= new_hist[0]->GetNbinsX(); j++) {// note bin numbering convention
			nominal_yield[j - 1] = hist_prediction->GetBinContent(j);
		}

		// calculate and store, per bin nominal_yield;
		TH1F* hist_prediction_wfakes = (TH1F*)new_hist[0]->Clone();
		for (int p_ = 1; p_ < n_process - 1; p_++) {// add MC 
			hist_prediction_wfakes->Add(new_hist[p_], 1.0);
		}

		// X BIN CENTER
		double* x_bin_center = new double[new_hist[0]->GetNbinsX()];
		for (int j = 1; j <= new_hist[0]->GetNbinsX(); j++) x_bin_center[j - 1] = new_hist[0]->GetBinCenter(j);
		// X BIN DOWN
		double* x_bin_down = new double[new_hist[0]->GetNbinsX()];
		for (int j = 0; j < new_hist[0]->GetNbinsX(); j++) x_bin_down[j] = new_hist[0]->GetBinWidth(j) / 2;
		// X BIN UP
		double* x_bin_up = new double[new_hist[0]->GetNbinsX()];
		for (int j = 0; j < new_hist[0]->GetNbinsX(); j++) x_bin_up[j] = new_hist[0]->GetBinWidth(j) / 2;
		// Y BIN DOWN
		double* y_bin_down = new double[new_hist[0]->GetNbinsX()];
		for (int j = 0; j < new_hist[0]->GetNbinsX(); j++) y_bin_down[j] = 0.0;
		// Y BIN UP
		double* y_bin_up = new double[new_hist[0]->GetNbinsX()];
		for (int j = 0; j < new_hist[0]->GetNbinsX(); j++) y_bin_up[j] = 0.0;

		// SYSTEMATIC DELTA UP
		double** systematic_delta_up = new double*[n_systematics - 1];
		for (int i = 0; i < n_systematics - 1; i++){
			systematic_delta_up[i] = new double[new_hist[0]->GetNbinsX()];
			for (int j = 0; j < new_hist[0]->GetNbinsX(); j++) systematic_delta_up[i][j] = 0.0;
		}


		// SYSTEMATIC DELTA DOWN
		double** systematic_delta_down = new double*[n_systematics - 1];
		for (int i = 0; i < n_systematics - 1; i++){
			systematic_delta_down[i] = new double[new_hist[0]->GetNbinsX()];
			for (int j = 0; j < new_hist[0]->GetNbinsX(); j++) systematic_delta_down[i][j] = 0.0;
		}

		// One Sided SYSTEMATIC DELTA UP
		double** one_sided_systematic_delta_up = new double*[n_systematics - 1];
		for (int i = 0; i < n_systematics - 1; i++){
			one_sided_systematic_delta_up[i] = new double[new_hist[0]->GetNbinsX()];
			for (int j = 0; j < new_hist[0]->GetNbinsX(); j++) one_sided_systematic_delta_up[i][j] = 0.0;
		}

		// One Sided SYSTEMATIC DELTA DOWN
		double** one_sided_systematic_delta_down = new double*[n_systematics - 1];
		for (int i = 0; i < n_systematics - 1; i++){
			one_sided_systematic_delta_down[i] = new double[new_hist[0]->GetNbinsX()];
			for (int j = 0; j < new_hist[0]->GetNbinsX(); j++) one_sided_systematic_delta_down[i][j] = 0.0;
		}

		// FAKE_SYSTEMATIC DELTA UP
		double** fake_systematic_delta_up = new double*[n_fakesyst - 1];
		for (int i = 0; i < n_fakesyst - 1; i++){
			fake_systematic_delta_up[i] = new double[new_hist[0]->GetNbinsX()];
			for (int j = 0; j < new_hist[0]->GetNbinsX(); j++) fake_systematic_delta_up[i][j] = 0.0;
		}

		// FAKE_SYSTEMATIC DELTA DOWN
		double** fake_systematic_delta_down = new double*[n_fakesyst - 1];
		for (int i = 0; i < n_fakesyst - 1; i++){
			fake_systematic_delta_down[i] = new double[new_hist[0]->GetNbinsX()];
			for (int j = 0; j < new_hist[0]->GetNbinsX(); j++) fake_systematic_delta_down[i][j] = 0.0;
		}

		// adhox_ZV
		double* adhoc_ZV = new double[new_hist[0]->GetNbinsX()];
		for (int j = 0; j < new_hist[0]->GetNbinsX(); j++) adhoc_ZV[j] = 0.0;

		// OBJECT SYSTEMATICS
		// OBJECT SYSTEMATICS
		// OBJECT SYSTEMATICS

		for (int s_ = 1; s_ < n_systematics; s_++) {
			cout << "OBJECT: " << systematics[s_] << endl;// Print 
			for (int p_ = 0; p_ < n_mc_process; p_++) {

				// open the file, and grab the right TChain
				stringstream filename_;
				filename_ << DIR_0 << BG_FILE;

				stringstream chain_name_up;
				chain_name_up << process[p_] << "_" << systematics[s_] << "UP";

				stringstream chain_name_down;
				chain_name_down << process[p_] << "_" << systematics[s_] << "DOWN";

				// Read in up chain
				TChain* chain_up;
				chain_up = new TChain(chain_name_up.str().data());

				chain_up->Add(filename_.str().data());

				if (chain_up->IsZombie()) {
					cout << chain_name_up.str() << "not found." << endl;
					continue;
				}

				// Rean in down chain
				TChain* chain_down;
				chain_down = new TChain(chain_name_down.str().data());

				chain_down->Add(filename_.str().data());

				if (chain_down->IsZombie()) {
					cout << chain_name_down.str() << "not found." << endl;
					continue;
				}

				TCut sel = TCut("");// key TCut
				TCut weight = TCut("eventweight");// selects eventweight leaf
				sel += TCut(signal_tcut.data());// string.data()

				// START reweighting segment
				TCut tc_reweight = TCut("1.0");
				if (plot_.reweight.compare("") != 0) {
					tc_reweight = TCut(plot_.reweight.data());
				}
				// END reweighting segment

				// temp histogram to get the variation
				TH1F* temp_hist_up = new TH1F("temp_hist_up", "temp_hist_up", new_hist[p_]->GetNbinsX(), plot_.x_range_min, plot_.x_range_max);
				TH1F* temp_hist_down = new TH1F("temp_hist_down", "temp_hist_down", new_hist[p_]->GetNbinsX(), plot_.x_range_min, plot_.x_range_max);

				stringstream draw_up;
				draw_up << plot_.root_get_member << plot_.root_get_factor << " >> +" << "temp_hist_up";
				stringstream draw_down;
				draw_down << plot_.root_get_member << plot_.root_get_factor << " >> +" << "temp_hist_down";

				// UP SECTION ----------------------------------------------------
				chain_up->Draw(draw_up.str().data(), sel * weight * tc_reweight);
				for (int j = 1; j <= new_hist[p_]->GetNbinsX(); j++) {// note bin numbering convention
					systematic_delta_up[s_ - 1][j - 1] += temp_hist_up->GetBinContent(j);
				}

				// DOWN SECTION --------------------------------------------------
				chain_down->Draw(draw_down.str().data(), sel * weight * tc_reweight);
				for (int j = 1; j <= new_hist[p_]->GetNbinsX(); j++) {// note bin numbering convention
					systematic_delta_down[s_ - 1][j - 1] += temp_hist_down->GetBinContent(j);
				}

				delete temp_hist_up;
				temp_hist_up = 0;
				delete temp_hist_down;
				temp_hist_down = 0;

				delete chain_up;
				delete chain_down;
			}
		}

		// ONE SIDED
		// ONE SIDED
		// ONE SIDED

		for (int s_ = 1; s_ < n_one_sided_syst; s_++) {

			cout << "ONE-SIDED: " << one_sided_syst[s_] << endl;// Print 

			for (int p_ = 0; p_ < n_mc_process; p_++) {

				// open the file, and grab the right TChain
				stringstream filename_;
				filename_ << DIR_0 << BG_FILE;

				stringstream chain_name_up;
				chain_name_up << process[p_] << "_" << one_sided_syst[s_];

				stringstream chain_name_down;
				chain_name_down << process[p_] << "_" << one_sided_syst[s_];

				// Read in up chain
				TChain* chain_up;
				chain_up = new TChain(chain_name_up.str().data());

				chain_up->Add(filename_.str().data());

				if (chain_up->IsZombie()) {
					cout << chain_name_up.str() << "not found." << endl;
					continue;
				}

				// Read in down chain
				TChain* chain_down;
				chain_down = new TChain(chain_name_down.str().data());

				chain_down->Add(filename_.str().data());

				if (chain_down->IsZombie()) {
					cout << chain_name_down.str() << "not found." << endl;
					continue;
				}

				TCut sel = TCut("");// key TCut
				TCut weight = TCut("eventweight");// selects eventweight leaf
				sel += TCut(signal_tcut.data());// string.data()

				// START reweighting segment
				TCut tc_reweight = TCut("1.0");
				if (plot_.reweight.compare("") != 0) {
					tc_reweight = TCut(plot_.reweight.data());
				}
				// END reweighting segment

				// temp histogram to get the variation
				TH1F* temp_hist_up = new TH1F("temp_hist_up", "temp_hist_up", new_hist[p_]->GetNbinsX(), plot_.x_range_min, plot_.x_range_max);
				TH1F* temp_hist_down = new TH1F("temp_hist_down", "temp_hist_down", new_hist[p_]->GetNbinsX(), plot_.x_range_min, plot_.x_range_max);

				stringstream draw_up;
				draw_up << plot_.root_get_member << plot_.root_get_factor << " >> +" << "temp_hist_up";
				stringstream draw_down;
				draw_down << plot_.root_get_member << plot_.root_get_factor << " >> +" << "temp_hist_down";

				// UP SECTION ----------------------------------------------------
				chain_up->Draw(draw_up.str().data(), sel * weight * tc_reweight);
				for (int j = 1; j <= new_hist[p_]->GetNbinsX(); j++) {// note bin numbering convention
					one_sided_systematic_delta_up[s_ - 1][j - 1] += temp_hist_up->GetBinContent(j);
				}

				// DOWN SECTION --------------------------------------------------
				chain_down->Draw(draw_down.str().data(), sel * weight * tc_reweight);
				for (int j = 1; j <= new_hist[p_]->GetNbinsX(); j++) {// note bin numbering convention
					one_sided_systematic_delta_down[s_ - 1][j - 1] += temp_hist_down->GetBinContent(j);
				}

				delete temp_hist_up;
				temp_hist_up = 0;
				delete temp_hist_down;
				temp_hist_down = 0;

				delete chain_up;
				delete chain_down;
			}
		}

		// WEIGHTS
		// WEIGHTS
		// WEIGHTS
		for (int s_ = 0; s_ < n_weight_syst; s_++) {

			cout << "WEIGHT: " << weight_syst[s_] << endl;// Print 
			for (int p_ = 0; p_ < n_mc_process; p_++) {

				// open the file, and grab the right TChain
				stringstream filename_;
				filename_ << DIR_0 << BG_FILE;

				stringstream chain_name_up;
				chain_name_up << process[p_] << "_" << systematics[0];

				stringstream chain_name_down;
				chain_name_down << process[p_] << "_" << systematics[0];

				// Read in up chain
				TChain* chain_up;
				chain_up = new TChain(chain_name_up.str().data());

				chain_up->Add(filename_.str().data());

				if (chain_up->IsZombie()) {
					cout << chain_name_up.str() << "not found." << endl;
					continue;
				}

				// Rean in down chain
				TChain* chain_down;
				chain_down = new TChain(chain_name_down.str().data());

				chain_down->Add(filename_.str().data());

				if (chain_down->IsZombie()) {
					cout << chain_name_down.str() << "not found." << endl;
					continue;
				}

				stringstream weight_name_up;
				weight_name_up << "syst_" << weight_syst[s_] << "UP";
				stringstream weight_name_down;
				weight_name_down << "syst_" << weight_syst[s_] << "DOWN";

				TCut sel = TCut("");// key TCut
				TCut weight = TCut("eventweight");// selects eventweight leaf
				TCut syst_var_up = TCut(weight_name_up.str().data());
				TCut syst_var_down = TCut(weight_name_down.str().data());
				sel += TCut(signal_tcut.data());// string.data()

				// START reweighting segment
				TCut tc_reweight = TCut("1.0");
				if (plot_.reweight.compare("") != 0) {
					tc_reweight = TCut(plot_.reweight.data());
				}
				// END reweighting segment


				// temp histogram to get the variation
				TH1F* temp_hist_up = new TH1F("temp_hist_up", "temp_hist_up", new_hist[p_]->GetNbinsX(), plot_.x_range_min, plot_.x_range_max);
				TH1F* temp_hist_down = new TH1F("temp_hist_down", "temp_hist_down", new_hist[p_]->GetNbinsX(), plot_.x_range_min, plot_.x_range_max);

				stringstream draw_up;
				draw_up << plot_.root_get_member << plot_.root_get_factor << " >> +" << "temp_hist_up";
				stringstream draw_down;
				draw_down << plot_.root_get_member << plot_.root_get_factor << " >> +" << "temp_hist_down";

				// UP SECTION ----------------------------------------------------
				chain_up->Draw(draw_up.str().data(), sel * weight * syst_var_up * tc_reweight);
				for (int j = 1; j <= new_hist[p_]->GetNbinsX(); j++) {// note bin numbering convention
					if (weight_syst[s_].compare("GEN") == 0) {
						gen_weight_delta_up[j - 1] += temp_hist_up->GetBinContent(j);
					}
					else if (weight_syst[s_].compare("XS") == 0) {
						xs_weight_delta_up[j - 1] += temp_hist_up->GetBinContent(j);
					}
					else {
						weight_delta_up[s_][j - 1] += temp_hist_up->GetBinContent(j);
					}
				}

				// DOWN SECTION --------------------------------------------------
				chain_down->Draw(draw_down.str().data(), sel * weight * syst_var_down * tc_reweight);
				for (int j = 1; j <= new_hist[p_]->GetNbinsX(); j++) {// note bin numbering convention
					if (weight_syst[s_].compare("GEN") == 0) {
						gen_weight_delta_down[j - 1] += temp_hist_down->GetBinContent(j);
					}
					else if (weight_syst[s_].compare("XS") == 0) {
						xs_weight_delta_down[j - 1] += temp_hist_down->GetBinContent(j);
					}
					else {
						weight_delta_down[s_][j - 1] += temp_hist_down->GetBinContent(j);
					}
				}

				delete temp_hist_up;
				temp_hist_up = 0;
				delete temp_hist_down;
				temp_hist_down = 0;

				delete chain_up;
				delete chain_down;
			}
		}


		// ADHOC ZV
		// ADHOC_ZV
		// ADHOC_ZV

		bool zv_sanity_once = true;

		for (int p_ = 0; p_ < n_mc_process; p_++) {

			cout << process_name[p_] << endl;// Print 

			// open the file, and grab the right TChain
			stringstream filename_;
			filename_ << DIR_0 << BG_FILE;

			stringstream chain_name_up;
			chain_name_up << process[p_] << "_" << systematics[0];

			// Read in up chain
			TChain* chain_up;
			chain_up = new TChain(chain_name_up.str().data());

			chain_up->Add(filename_.str().data());

			if (chain_up->IsZombie()) {
				cout << chain_name_up.str() << "not found." << endl;
				continue;
			}

			TCut syst_var_up;//

			syst_var_up = TCut("1.0");


			TCut sel = TCut(signal_tcut.data());// key TCut
			TCut weight = TCut("eventweight");// selects eventweight leaf

			// START reweighting segment
			TCut tc_reweight = TCut("1.0");
			if (plot_.reweight.compare("") != 0) {
				tc_reweight = TCut(plot_.reweight.data());
			}
			// END reweighting segment

			// temp histogram to get the variation
			TH1F* temp_hist_up = new TH1F("temp_hist_up", "temp_hist_up", new_hist[p_]->GetNbinsX(), plot_.x_range_min, plot_.x_range_max);

			stringstream draw_up;
			draw_up << plot_.root_get_member << plot_.root_get_factor << " >> +" << "temp_hist_up";

			// UP SECTION ----------------------------------------------------
			chain_up->Draw(draw_up.str().data(), sel * weight * syst_var_up * tc_reweight);
			for (int j = 1; j <= new_hist[p_]->GetNbinsX(); j++) {// note bin numbering convention
				adhoc_ZV[j - 1] += temp_hist_up->GetBinContent(j);
			}
			delete temp_hist_up;
			temp_hist_up = 0;

			delete chain_up;
		}

		{ // Calculate FAKE Systematic
			// Calculate FAKE Systematic
			// Calculate FAKE Systematic


			int p_ = n_mc_process;
			cout << process_name[p_] << endl;// Print 
			for (int s_ = 1; s_ < n_fakesyst; s_++) {
				// open the file, and grab the right TChain
				stringstream filename_;
				filename_ << DIR_0 << BG_FILE;

				stringstream chain_name_up;
				chain_name_up << process[p_] << "_" << fake_syst[s_] << "UP";

				stringstream chain_name_down;
				chain_name_down << process[p_] << "_" << fake_syst[s_] << "DOWN";

				// Read in up chain
				TChain* chain_up;
				chain_up = new TChain(chain_name_up.str().data());

				chain_up->Add(filename_.str().data());

				if (chain_up->IsZombie()) {
					cout << chain_name_up.str() << "not found." << endl;
					continue;
				}

				// Rean in down chain
				TChain* chain_down;
				chain_down = new TChain(chain_name_down.str().data());

				chain_down->Add(filename_.str().data());

				if (chain_down->IsZombie()) {
					cout << chain_name_down.str() << "not found." << endl;
					continue;
				}

				TCut sel = TCut("");// key TCut
				TCut weight = TCut("eventweight");// selects eventweight leaf
				sel += TCut(signal_tcut.data());// string.data()

				// temp histogram to get the variation
				TH1F* temp_hist_up = new TH1F("temp_hist_up", "temp_hist_up", new_hist[p_]->GetNbinsX(), plot_.x_range_min, plot_.x_range_max);
				TH1F* temp_hist_down = new TH1F("temp_hist_down", "temp_hist_down", new_hist[p_]->GetNbinsX(), plot_.x_range_min, plot_.x_range_max);

				stringstream draw_up;
				draw_up << plot_.root_get_member << plot_.root_get_factor << " >> +" << "temp_hist_up";
				stringstream draw_down;
				draw_down << plot_.root_get_member << plot_.root_get_factor << " >> +" << "temp_hist_down";

				// UP SECTION ----------------------------------------------------
				chain_up->Draw(draw_up.str().data(), sel * weight);
				double integral_up = temp_hist_up->Integral(0, -1);

				// calculate and store, per bin syst_deltas;
				for (int j = 1; j <= new_hist[p_]->GetNbinsX(); j++) {// note bin numbering convention
					fake_systematic_delta_up[s_ - 1][j - 1] += temp_hist_up->GetBinContent(j);
				}

				// DOWN SECTION --------------------------------------------------
				chain_down->Draw(draw_down.str().data(), sel * weight);
				double integral_down = temp_hist_down->Integral(0, -1);

				// calculate and store, per bin syst_deltas;
				for (int j = 1; j <= new_hist[p_]->GetNbinsX(); j++) {// note bin numbering convention
					fake_systematic_delta_down[s_ - 1][j - 1] += temp_hist_down->GetBinContent(j);
				}

				cout << setprecision(4);

				cout << chain_name_up.str() << ": " << fake_systematic_delta_up[s_ - 1][0] << ", " << fake_systematic_delta_up[s_ - 1][1] << endl;
				cout << chain_name_down.str() << ": " << fake_systematic_delta_down[s_ - 1][0] << ", " << fake_systematic_delta_down[s_ - 1][1] << endl;

				delete temp_hist_up;
				temp_hist_up = 0;
				delete temp_hist_down;
				temp_hist_down = 0;

				delete chain_up;
				delete chain_down;
			}
		}

		// add up systematics
		// complication: we add these in quadrature. But the fakes add linearly!
		// we account this below.

		for (int j = 0; j < new_hist[0]->GetNbinsX(); j++) {
			double temp_total_up;
			double temp_total_down;

			// object systematics
			for (int i = 0; i < n_systematics - 1; i++){
				double delta_up = systematic_delta_up[i][j] - nominal_yield[j];
				double delta_down = systematic_delta_down[i][j] - nominal_yield[j];

				y_bin_up[j] += delta_up * delta_up;
				y_bin_down[j] += delta_down * delta_down;
			}

			// one-sided object systematics
			for (int i = 0; i < n_one_sided_syst - 1; i++){
				double delta_up = one_sided_systematic_delta_up[i][j] - nominal_yield[j];

				if (delta_up > 0) {
					y_bin_up[j] += delta_up * delta_up;
				}
				else {
					y_bin_down[j] += delta_up * delta_up;
				}
			}

			// weight systematics
			for (int i = 0; i < n_weight_syst; i++){
				if (weight_syst[i].compare("GEN") == 0 || weight_syst[i].compare("XS") == 0) continue;

				double delta_up = weight_delta_up[i][j] - nominal_yield[j];
				double delta_down = weight_delta_down[i][j] - nominal_yield[j];

				y_bin_up[j] += delta_up * delta_up;
				y_bin_down[j] += delta_down * delta_down;
			}

			// gen
			{
				double delta_up = gen_weight_delta_up[j] - nominal_yield[j];
				double delta_down = gen_weight_delta_down[j] - nominal_yield[j];

				y_bin_up[j] += delta_up * delta_up;
				y_bin_down[j] += delta_down * delta_down;
			}

			// xs
			{
				double delta_up = xs_weight_delta_up[j] - nominal_yield[j];
				double delta_down = xs_weight_delta_down[j] - nominal_yield[j];

				y_bin_up[j] += delta_up * delta_up;
				y_bin_down[j] += delta_down * delta_down;
			}
		}

		vector<double> y_bin_ratio_up;
		vector<double> y_bin_ratio_down;

		for (int j = 0; j < new_hist[0]->GetNbinsX(); j++) {
			double fake_cache_syst_up = 0;
			double fake_cache_syst_down = 0;
			for (int i = 0; i < n_fakesyst - 1; i++){
				double delta_up = fake_systematic_delta_up[i][j] - new_hist[n_process - 2]->GetBinContent(j + 1);
				double delta_down = fake_systematic_delta_down[i][j] - new_hist[n_process - 2]->GetBinContent(j + 1);

				y_bin_up[j] += delta_up * delta_up;
				y_bin_down[j] += delta_down * delta_down;
			}

			//// add adhoc ZV
			//double delta_adhocZV = adhoc_ZV[j] - nominal_yield[j];
			//y_bin_up[j] += delta_adhocZV * delta_adhocZV;
			//y_bin_down[j] += delta_adhocZV * delta_adhocZV;

			// add statistical error
			y_bin_up[j] += (hist_prediction->GetBinError(j + 1)) * (hist_prediction->GetBinError(j + 1));
			y_bin_down[j] += (hist_prediction->GetBinError(j + 1)) * (hist_prediction->GetBinError(j + 1));

			y_bin_up[j] += (hist_prediction_wfakes->GetBinError(j + 1)) * (hist_prediction_wfakes->GetBinError(j + 1));
			y_bin_down[j] += (hist_prediction_wfakes->GetBinError(j + 1)) * (hist_prediction_wfakes->GetBinError(j + 1));

			y_bin_up[j] = sqrt(y_bin_up[j]);
			y_bin_down[j] = sqrt(y_bin_down[j]);

			// for bottom pad
			y_bin_ratio_up.push_back(y_bin_up[j]);
			y_bin_ratio_down.push_back(y_bin_down[j]);

			// add the fakes to the nominal monte carlo
			nominal_yield[j] += new_hist[n_process - 2]->GetBinContent(j + 1);
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
		out_file << DIR_1 << plot_.long_name << ".eps";

		stringstream out_file_image;
		out_file_image << DIR_1 << plot_.long_name << ".png";

		TCanvas* tc = new TCanvas("canvas", "", 768, 768);
		TLegend* tl = new TLegend(0.68, 0.50, 0.88, 0.88);

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

		{	// FAKES on top
			int p_ = n_process - 2; // index of fakes
			new_hist[p_]->SetMarkerStyle(1);
			new_hist[p_]->SetLineColor(kBlack);
			new_hist[p_]->SetLineWidth(2);
			new_hist[p_]->SetFillColor(stack_colors[p_]);
			hist_stack->Add(new_hist[p_]);
		}

		for (int p_ = 0; p_ < n_process - 2; p_++) {
			new_hist[process_INSERT_ORDER[p_]]->SetMarkerStyle(1);
			new_hist[process_INSERT_ORDER[p_]]->SetLineColor(kBlack);
			new_hist[process_INSERT_ORDER[p_]]->SetLineWidth(2);
			new_hist[process_INSERT_ORDER[p_]]->SetFillColor(stack_colors[process_INSERT_ORDER[p_]]);
			hist_stack->Add(new_hist[process_INSERT_ORDER[p_]]);
		}

		for (int p_ = 0; p_ < n_process - 2; p_++) {
			tl->AddEntry(new_hist[p_], process_name[p_].data(), "f");
		}

		// Fakes last
		tl->AddEntry(new_hist[n_process - 2], process_name[n_process - 2].data(), "f");

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

		convertErrorsToPoisson(new_hist[n_process - 1], Data);

		TGraphAsymmErrors* asym_errors = new TGraphAsymmErrors(new_hist[0]->GetNbinsX(), x_bin_center, nominal_yield, x_bin_down, x_bin_up, y_bin_down, y_bin_up);
		asym_errors->SetFillStyle(3004);
		asym_errors->SetFillColor(kGray + 3);
		tl->AddEntry(asym_errors, "Bkg. Uncert.", "f");
		asym_errors->Draw("option same 02");

		// Decoration
		char annoyingLabel1[100] = "#bf{#it{ATLAS}} Internal", annoyingLabel2[100] = "#scale[0.6]{#int} L dt = 20.3 fb^{-1} #sqrt{s} = 8 TeV";
		myText(0.184, 0.83, kBlack, annoyingLabel1);
		myText(0.184, 0.75, kBlack, annoyingLabel2);
		tl->Draw();

		Data->SetLineWidth(2);
		Data->SetMarkerStyle(8);
		Data->SetMarkerSize(1.5);
		Data->SetLineColor(1);
		Data->Draw("option same pz");

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
		double* nominal_one = new double[new_hist[0]->GetNbinsX()];
		for (int j = 0; j < new_hist[0]->GetNbinsX(); j++) nominal_one[j] = 1.0;

		double* value_zero = new double[new_hist[0]->GetNbinsX()];
		for (int j = 0; j < new_hist[0]->GetNbinsX(); j++) value_zero[j] = 0.0;

		double* y_bin_down_ratio = new double[new_hist[0]->GetNbinsX()];
		for (int j = 0; j < new_hist[0]->GetNbinsX(); j++) {
			if (nominal_yield[j] != 0) {
				y_bin_down_ratio[j] = y_bin_ratio_down[j] / nominal_yield[j];
			}
			else {
				y_bin_down_ratio[j] = 0;
			}
		}

		double* y_bin_up_ratio = new double[new_hist[0]->GetNbinsX()];
		for (int j = 0; j < new_hist[0]->GetNbinsX(); j++) {
			if (nominal_yield[j] != 0) {
				y_bin_up_ratio[j] = y_bin_ratio_up[j] / nominal_yield[j];
			}
			else {
				y_bin_up_ratio[j] = 0;
			}
		}

		double* data_ratio = new double[new_hist[0]->GetNbinsX()];
		for (int j = 0; j < new_hist[0]->GetNbinsX(); j++) {
			if (nominal_yield[j] > 0 && new_hist[n_process - 1]->GetBinContent(j + 1) > 0) {
				data_ratio[j] = new_hist[n_process - 1]->GetBinContent(j + 1) / nominal_yield[j];
			}
			else {
				data_ratio[j] = 100.0;
			}
		}

		double* data_err_down_ratio = new double[new_hist[0]->GetNbinsX()];
		double* data_err_up_ratio = new double[new_hist[0]->GetNbinsX()];
		for (int j = 0; j < new_hist[0]->GetNbinsX(); j++) {
			if (nominal_yield[j] > 0) {

				data_err_down_ratio[j] = 1.0 - ((nominal_yield[j] - Data->GetErrorYlow(j)) / nominal_yield[j]);

				data_err_up_ratio[j] = ((nominal_yield[j] + Data->GetErrorYhigh(j)) / nominal_yield[j]) - 1.0;
			}
			else {
				data_err_down_ratio[j] = 0.0;
				data_err_up_ratio[j] = 0.0;
			}
		}

		TGraphAsymmErrors* asym_errors2 = new TGraphAsymmErrors(new_hist[0]->GetNbinsX(), x_bin_center, nominal_one, x_bin_down, x_bin_up, y_bin_down_ratio, y_bin_up_ratio);

		asym_errors2->SetFillStyle(3004);
		asym_errors2->SetFillColor(kGray + 3);

		TH1F* hist_ratio = (TH1F*)new_hist[0]->Clone();
		TGraphAsymmErrors* Data2 = new TGraphAsymmErrors(new_hist[0]->GetNbinsX(), x_bin_center, data_ratio, value_zero, value_zero, data_err_down_ratio, data_err_up_ratio);

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

		Data2->SetTitle("");
		Data2->SetLineWidth(2);
		Data2->SetMarkerStyle(8);
		Data2->SetMarkerSize(1.5);
		Data2->SetLineColor(1);
		Data2->Draw("option same 0pz");

		line_.Draw("same");
		line_up.Draw("same");
		line_down.Draw("same");

		asym_errors2->Draw("option same 02");

		stringstream out_file_root;
		out_file_root << DIR_1 << plot_.long_name << ".root";

		tc->Print(out_file.str().data());
		tc->Print(out_file_image.str().data());
		tc->SaveAs(out_file_root.str().data());
		tc->Close();

		delete asym_errors;
		delete asym_errors2;
		delete Data;
		delete Data2;

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

		for (int i = 0; i < n_weight_syst; i++) delete[] weight_delta_up[i];
		delete[] weight_delta_up;

		for (int i = 0; i < n_weight_syst; i++) delete[] weight_delta_down[i];
		delete[] weight_delta_down;

		delete[] gen_weight_delta_up;
		delete[] xs_weight_delta_up;
		delete[] gen_weight_delta_down;
		delete[] xs_weight_delta_down;

		for (int i = 0; i < n_systematics - 1; i++) delete[] systematic_delta_up[i];
		delete[] systematic_delta_up;

		for (int i = 0; i < n_systematics - 1; i++) delete[] systematic_delta_down[i];
		delete[] systematic_delta_down;

		for (int i = 0; i < n_one_sided_syst - 1; i++) delete[] one_sided_systematic_delta_up[i];
		delete[] one_sided_systematic_delta_up;

		for (int i = 0; i < n_one_sided_syst - 1; i++) delete[] one_sided_systematic_delta_down[i];
		delete[] one_sided_systematic_delta_down;

		for (int i = 0; i < n_fakesyst - 1; i++) delete[] fake_systematic_delta_up[i];
		delete[] fake_systematic_delta_up;

		for (int i = 0; i < n_fakesyst - 1; i++) delete[] fake_systematic_delta_down[i];
		delete[] fake_systematic_delta_down;

		delete hist_prediction_wfakes;
		delete hist_prediction;
		delete hist_ratio;
		for (int p_ = 0; p_ < n_process; p_++) {
			delete new_hist[p_];
		}
		delete[] new_hist;
		delete hist_stack;
		delete tl;
		delete tc;
	}
	cout << "Done." << endl;
	return 0;
}
