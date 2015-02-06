plot_.root_get_member_name = ""; // critical

plot_.is_Log = true;

// VR_Top
// VR_Top
// VR_Top

// VR_Top
// Plot: ll_M
plot_.root_get_member = "ll_M";
plot_.root_get_factor = "";
plot_.signal_region = "VR_Top";
plot_.request_cuts = "";
plot_.blind_data = "";
plot_.x_label = "m_{ll} [GeV]";
plot_.y_label = "Events / 25 GeV";
plot_.x_bin_width = 25.0;
plot_.x_range_min = 0.0;
plot_.x_range_max = 750.0;
plot_.y_range_min = 1.0e-1;
plot_.y_range_max = 1.0e9;
plot_.long_name = "VR_Top_ll_M";
vect_plot.push_back(plot_);

// VR_Top
// Plot: ll_Pt
plot_.root_get_member = "ll_Pt";
plot_.root_get_factor = "";
plot_.signal_region = "VR_Top";
plot_.request_cuts = "";
plot_.blind_data = "";
plot_.x_label = "p_{T} (ll) [GeV]";
plot_.y_label = "Events / 15 GeV";
plot_.x_bin_width = 15.0;
plot_.x_range_min = 0.0;
plot_.x_range_max = 500.0;
plot_.y_range_min = 1.0e-1;
plot_.y_range_max = 1.0e9;
plot_.long_name = "VR_Top_ll_Pt";
vect_plot.push_back(plot_);

// VR_Top
// Plot: met
plot_.root_get_member = "lept1Pt";
plot_.root_get_factor = "";
plot_.signal_region = "VR_Top";
plot_.request_cuts = "";
plot_.blind_data = "";
plot_.x_label = "p_{T, l1} [GeV]";
plot_.y_label = "Events / 20 GeV";
plot_.x_bin_width = 20.0;
plot_.x_range_min = 0.0;
plot_.x_range_max = 500.0;
plot_.y_range_min = 1.0e-1;
plot_.y_range_max = 1.0e9;
plot_.long_name = "VR_Top_lept1Pt";
vect_plot.push_back(plot_);

// VR_Top
// Plot: met
plot_.root_get_member = "lept2Pt";
plot_.root_get_factor = "";
plot_.signal_region = "VR_Top";
plot_.request_cuts = "";
plot_.blind_data = "";
plot_.x_label = "p_{T, l2} [GeV]";
plot_.y_label = "Events / 20 GeV";
plot_.x_bin_width = 20.0;
plot_.x_range_min = 0.0;
plot_.x_range_max = 500.0;
plot_.y_range_min = 1.0e-1;
plot_.y_range_max = 1.0e9;
plot_.long_name = "VR_Top_lept2Pt";
vect_plot.push_back(plot_);

// VR_Top
// Plot: met
plot_.root_get_member = "met";
plot_.root_get_factor = "";
plot_.signal_region = "VR_Top";
plot_.request_cuts = "";
plot_.blind_data = "";
plot_.x_label = "E_{T}^{miss} [GeV]";
plot_.y_label = "Events / 10 GeV";
plot_.x_bin_width = 10.0;
plot_.x_range_min = 0.0;
plot_.x_range_max = 300.0;
plot_.y_range_min = 1.0e-1;
plot_.y_range_max = 1.0e9;
plot_.long_name = "VR_Top_met";
vect_plot.push_back(plot_);

// VR_Top
// Plot: ll_deltaPhi
plot_.root_get_member = "acos(cos(lept1Phi - lept2Phi))";
plot_.root_get_member_name = "ll_deltaPhi";
plot_.root_get_factor = "";
plot_.signal_region = "VR_Top";
plot_.request_cuts = "";
plot_.blind_data = "";
plot_.x_label = "#Delta#phi (ll)";
plot_.y_label = "Events";
plot_.x_bin_width = Pi / 24;
plot_.x_range_min = 0.0;
plot_.x_range_max = 1.5 * Pi;
plot_.y_range_min = 1.0e-1;
plot_.y_range_max = 1.0e9;
plot_.long_name = "VR_Top_ll_deltaPhi";
vect_plot.push_back(plot_);
plot_.root_get_member_name = "";
plot_.root_get_member_name = "";
plot_.root_get_member_name = "";

// Plot: nCentralLightJets
// selection region
plot_.root_get_member = "nCentralLightJets";
plot_.root_get_factor = "";
plot_.signal_region = "VR_Top";
plot_.request_cuts = "";
plot_.blind_data = "";
// plot setup
plot_.x_label = "N_{jets}-central";
plot_.y_label = "Events";
plot_.x_bin_width = 1.0;
plot_.x_range_min = 0.0;
plot_.x_range_max = 10.0;
plot_.y_range_min = 1.0e-1;
plot_.y_range_max = 1.0e9;
plot_.long_name = "VR_Top_nCentralLightJets";
vect_plot.push_back(plot_);

// Plot: nForwardJets
// selection region
plot_.root_get_member = "nForwardJets";
plot_.root_get_factor = "";
plot_.signal_region = "VR_Top";
plot_.request_cuts = "";
plot_.blind_data = "";
// plot setup
plot_.x_label = "N_{for-jets}";
plot_.y_label = "Events";
plot_.x_bin_width = 1.0;
plot_.x_range_min = 0.0;
plot_.x_range_max = 10.0;
plot_.y_range_min = 1.0e-1;
plot_.y_range_max = 1.0e9;
plot_.long_name = "VR_Top_nForwardJets";
vect_plot.push_back(plot_);

// Plot: nCentralBJets
// selection region
plot_.root_get_member = "nCentralBJets";
plot_.root_get_factor = "";
plot_.signal_region = "VR_Top";
plot_.request_cuts = "";
plot_.blind_data = "";
// plot setup
plot_.x_label = "N_{b-jets}";
plot_.y_label = "Events";
plot_.x_bin_width = 1.0;
plot_.x_range_min = 0.0;
plot_.x_range_max = 10.0;
plot_.y_range_min = 1.0e-1;
plot_.y_range_max = 1.0e9;
plot_.long_name = "VR_Top_nCentralBJets";
vect_plot.push_back(plot_);

// VR_Top
// Plot: MT2
plot_.root_get_member = "MT2";
plot_.root_get_factor = "";
plot_.signal_region = "VR_Top";
plot_.request_cuts = "";
plot_.blind_data = "";
plot_.x_label = "m_{T2} [GeV]";
plot_.y_label = "Events / 10 GeV";
plot_.x_bin_width = 5.0;
plot_.x_range_min = 0.0;
plot_.x_range_max = 200.0;
plot_.y_range_min = 1.0e-1;
plot_.y_range_max = 1.0e9;
plot_.long_name = "VR_Top_MT2";
vect_plot.push_back(plot_);

// VR_Top
// Plot: m_eff
plot_.root_get_member = "lept1Pt + lept2Pt + met";
plot_.root_get_member_name = "m_eff";
plot_.root_get_factor = "";
plot_.signal_region = "VR_Top";
plot_.request_cuts = "";
plot_.blind_data = "";
plot_.x_label = "m_{lepeff} [GeV]";
plot_.y_label = "Events / 25 GeV";
plot_.x_bin_width = 25.0;
plot_.x_range_min = 0.0;
plot_.x_range_max = 1200.0;
plot_.y_range_min = 1.0e-1;
plot_.y_range_max = 1.0e9;
plot_.long_name = "VR_Top_m_eff";
vect_plot.push_back(plot_);

plot_.root_get_member_name = "";
plot_.root_get_member_name = "";
plot_.root_get_member_name = "";



























plot_.root_get_member_name = ""; // critical

plot_.is_Log = true;

// VR_WW
// VR_WW
// VR_WW

// VR_WW
// Plot: ll_M
plot_.root_get_member = "ll_M";
plot_.root_get_factor = "";
plot_.signal_region = "VR_WW";
plot_.request_cuts = "";
plot_.blind_data = "";
plot_.x_label = "m_{ll} [GeV]";
plot_.y_label = "Events / 25 GeV";
plot_.x_bin_width = 25.0;
plot_.x_range_min = 0.0;
plot_.x_range_max = 750.0;
plot_.y_range_min = 1.0e-1;
plot_.y_range_max = 1.0e9;
plot_.long_name = "VR_WW_ll_M";
vect_plot.push_back(plot_);

// VR_WW
// Plot: ll_Pt
plot_.root_get_member = "ll_Pt";
plot_.root_get_factor = "";
plot_.signal_region = "VR_WW";
plot_.request_cuts = "";
plot_.blind_data = "";
plot_.x_label = "p_{T} (ll) [GeV]";
plot_.y_label = "Events / 15 GeV";
plot_.x_bin_width = 15.0;
plot_.x_range_min = 0.0;
plot_.x_range_max = 500.0;
plot_.y_range_min = 1.0e-1;
plot_.y_range_max = 1.0e9;
plot_.long_name = "VR_WW_ll_Pt";
vect_plot.push_back(plot_);

// VR_WW
// Plot: met
plot_.root_get_member = "lept1Pt";
plot_.root_get_factor = "";
plot_.signal_region = "VR_WW";
plot_.request_cuts = "";
plot_.blind_data = "";
plot_.x_label = "p_{T, l1} [GeV]";
plot_.y_label = "Events / 20 GeV";
plot_.x_bin_width = 20.0;
plot_.x_range_min = 0.0;
plot_.x_range_max = 500.0;
plot_.y_range_min = 1.0e-1;
plot_.y_range_max = 1.0e9;
plot_.long_name = "VR_WW_lept1Pt";
vect_plot.push_back(plot_);

// VR_WW
// Plot: met
plot_.root_get_member = "lept2Pt";
plot_.root_get_factor = "";
plot_.signal_region = "VR_WW";
plot_.request_cuts = "";
plot_.blind_data = "";
plot_.x_label = "p_{T, l2} [GeV]";
plot_.y_label = "Events / 20 GeV";
plot_.x_bin_width = 20.0;
plot_.x_range_min = 0.0;
plot_.x_range_max = 500.0;
plot_.y_range_min = 1.0e-1;
plot_.y_range_max = 1.0e9;
plot_.long_name = "VR_WW_lept2Pt";
vect_plot.push_back(plot_);

// VR_WW
// Plot: met
plot_.root_get_member = "met";
plot_.root_get_factor = "";
plot_.signal_region = "VR_WW";
plot_.request_cuts = "";
plot_.blind_data = "";
plot_.x_label = "E_{T}^{miss} [GeV]";
plot_.y_label = "Events / 10 GeV";
plot_.x_bin_width = 10.0;
plot_.x_range_min = 0.0;
plot_.x_range_max = 300.0;
plot_.y_range_min = 1.0e-1;
plot_.y_range_max = 1.0e9;
plot_.long_name = "VR_WW_met";
vect_plot.push_back(plot_);

// VR_WW
// Plot: ll_deltaPhi
plot_.root_get_member = "acos(cos(lept1Phi - lept2Phi))";
plot_.root_get_member_name = "ll_deltaPhi";
plot_.root_get_factor = "";
plot_.signal_region = "VR_WW";
plot_.request_cuts = "";
plot_.blind_data = "";
plot_.x_label = "#Delta#phi (ll)";
plot_.y_label = "Events";
plot_.x_bin_width = Pi / 24;
plot_.x_range_min = 0.0;
plot_.x_range_max = 1.5 * Pi;
plot_.y_range_min = 1.0e-1;
plot_.y_range_max = 1.0e9;
plot_.long_name = "VR_WW_ll_deltaPhi";
vect_plot.push_back(plot_);
plot_.root_get_member_name = "";
plot_.root_get_member_name = "";
plot_.root_get_member_name = "";

// Plot: nCentralLightJets
// selection region
plot_.root_get_member = "nCentralLightJets";
plot_.root_get_factor = "";
plot_.signal_region = "VR_WW";
plot_.request_cuts = "";
plot_.blind_data = "";
// plot setup
plot_.x_label = "N_{jets}-central";
plot_.y_label = "Events";
plot_.x_bin_width = 1.0;
plot_.x_range_min = 0.0;
plot_.x_range_max = 10.0;
plot_.y_range_min = 1.0e-1;
plot_.y_range_max = 1.0e9;
plot_.long_name = "VR_WW_nCentralLightJets";
vect_plot.push_back(plot_);

// Plot: nForwardJets
// selection region
plot_.root_get_member = "nForwardJets";
plot_.root_get_factor = "";
plot_.signal_region = "VR_WW";
plot_.request_cuts = "";
plot_.blind_data = "";
// plot setup
plot_.x_label = "N_{for-jets}";
plot_.y_label = "Events";
plot_.x_bin_width = 1.0;
plot_.x_range_min = 0.0;
plot_.x_range_max = 10.0;
plot_.y_range_min = 1.0e-1;
plot_.y_range_max = 1.0e9;
plot_.long_name = "VR_WW_nForwardJets";
vect_plot.push_back(plot_);

// Plot: nCentralBJets
// selection region
plot_.root_get_member = "nCentralBJets";
plot_.root_get_factor = "";
plot_.signal_region = "VR_WW";
plot_.request_cuts = "";
plot_.blind_data = "";
// plot setup
plot_.x_label = "N_{b-jets}";
plot_.y_label = "Events";
plot_.x_bin_width = 1.0;
plot_.x_range_min = 0.0;
plot_.x_range_max = 10.0;
plot_.y_range_min = 1.0e-1;
plot_.y_range_max = 1.0e9;
plot_.long_name = "VR_WW_nCentralBJets";
vect_plot.push_back(plot_);

// VR_WW
// Plot: MT2
plot_.root_get_member = "MT2";
plot_.root_get_factor = "";
plot_.signal_region = "VR_WW";
plot_.request_cuts = "";
plot_.blind_data = "";
plot_.x_label = "m_{T2} [GeV]";
plot_.y_label = "Events / 10 GeV";
plot_.x_bin_width = 5.0;
plot_.x_range_min = 0.0;
plot_.x_range_max = 200.0;
plot_.y_range_min = 1.0e-1;
plot_.y_range_max = 1.0e9;
plot_.long_name = "VR_WW_MT2";
vect_plot.push_back(plot_);

// VR_WW
// Plot: m_eff
plot_.root_get_member = "lept1Pt + lept2Pt + met";
plot_.root_get_member_name = "m_eff";
plot_.root_get_factor = "";
plot_.signal_region = "VR_WW";
plot_.request_cuts = "";
plot_.blind_data = "";
plot_.x_label = "m_{lepeff} [GeV]";
plot_.y_label = "Events / 25 GeV";
plot_.x_bin_width = 25.0;
plot_.x_range_min = 0.0;
plot_.x_range_max = 1200.0;
plot_.y_range_min = 1.0e-1;
plot_.y_range_max = 1.0e9;
plot_.long_name = "VR_WW_m_eff";
vect_plot.push_back(plot_);

plot_.root_get_member_name = "";
plot_.root_get_member_name = "";
plot_.root_get_member_name = "";
























plot_.root_get_member_name = ""; // critical

plot_.is_Log = true; 

// VR_Zmm
// VR_Zmm
// VR_Zmm

// VR_Zmm
// Plot: ll_M
plot_.root_get_member = "ll_M";
plot_.root_get_factor = "";
plot_.signal_region = "VR_Zmm";
plot_.request_cuts = "";
plot_.blind_data = "";
plot_.x_label = "m_{ll} [GeV]";
plot_.y_label = "Events / 25 GeV";
plot_.x_bin_width = 25.0;
plot_.x_range_min = 0.0;
plot_.x_range_max = 750.0;
plot_.y_range_min = 1.0e-1;
plot_.y_range_max = 1.0e9;
plot_.long_name = "VR_Zmm_ll_M";
vect_plot.push_back(plot_);

// VR_Zmm
// Plot: ll_Pt
plot_.root_get_member = "ll_Pt";
plot_.root_get_factor = "";
plot_.signal_region = "VR_Zmm";
plot_.request_cuts = "";
plot_.blind_data = "";
plot_.x_label = "p_{T} (ll) [GeV]";
plot_.y_label = "Events / 15 GeV";
plot_.x_bin_width = 15.0;
plot_.x_range_min = 0.0;
plot_.x_range_max = 500.0;
plot_.y_range_min = 1.0e-1;
plot_.y_range_max = 1.0e9;
plot_.long_name = "VR_Zmm_ll_Pt";
vect_plot.push_back(plot_);

// VR_Zmm
// Plot: met
plot_.root_get_member = "lept1Pt";
plot_.root_get_factor = "";
plot_.signal_region = "VR_Zmm";
plot_.request_cuts = "";
plot_.blind_data = "";
plot_.x_label = "p_{T, l1} [GeV]";
plot_.y_label = "Events / 20 GeV";
plot_.x_bin_width = 20.0;
plot_.x_range_min = 0.0;
plot_.x_range_max = 500.0;
plot_.y_range_min = 1.0e-1;
plot_.y_range_max = 1.0e9;
plot_.long_name = "VR_Zmm_lept1Pt";
vect_plot.push_back(plot_);

// VR_Zmm
// Plot: met
plot_.root_get_member = "lept2Pt";
plot_.root_get_factor = "";
plot_.signal_region = "VR_Zmm";
plot_.request_cuts = "";
plot_.blind_data = "";
plot_.x_label = "p_{T, l2} [GeV]";
plot_.y_label = "Events / 20 GeV";
plot_.x_bin_width = 20.0;
plot_.x_range_min = 0.0;
plot_.x_range_max = 500.0;
plot_.y_range_min = 1.0e-1;
plot_.y_range_max = 1.0e9;
plot_.long_name = "VR_Zmm_lept2Pt";
vect_plot.push_back(plot_);

// VR_Zmm
// Plot: met
plot_.root_get_member = "met";
plot_.root_get_factor = "";
plot_.signal_region = "VR_Zmm";
plot_.request_cuts = "";
plot_.blind_data = "";
plot_.x_label = "E_{T}^{miss} [GeV]";
plot_.y_label = "Events / 10 GeV";
plot_.x_bin_width = 10.0;
plot_.x_range_min = 0.0;
plot_.x_range_max = 300.0;
plot_.y_range_min = 1.0e-1;
plot_.y_range_max = 1.0e9;
plot_.long_name = "VR_Zmm_met";
vect_plot.push_back(plot_);

// VR_Zmm
// Plot: ll_deltaPhi
plot_.root_get_member = "acos(cos(lept1Phi - lept2Phi))";
plot_.root_get_member_name = "ll_deltaPhi";
plot_.root_get_factor = "";
plot_.signal_region = "VR_Zmm";
plot_.request_cuts = "";
plot_.blind_data = "";
plot_.x_label = "#Delta#phi (ll)";
plot_.y_label = "Events";
plot_.x_bin_width = Pi / 24;
plot_.x_range_min = 0.0;
plot_.x_range_max = 1.5 * Pi;
plot_.y_range_min = 1.0e-1;
plot_.y_range_max = 1.0e9;
plot_.long_name = "VR_Zmm_ll_deltaPhi";
vect_plot.push_back(plot_);
plot_.root_get_member_name = "";
plot_.root_get_member_name = "";
plot_.root_get_member_name = "";

// Plot: nCentralLightJets
// selection region
plot_.root_get_member = "nCentralLightJets";
plot_.root_get_factor = "";
plot_.signal_region = "VR_Zmm";
plot_.request_cuts = "";
plot_.blind_data = "";
// plot setup
plot_.x_label = "N_{jets}-central";
plot_.y_label = "Events";
plot_.x_bin_width = 1.0;
plot_.x_range_min = 0.0;
plot_.x_range_max = 10.0;
plot_.y_range_min = 1.0e-1;
plot_.y_range_max = 1.0e9;
plot_.long_name = "VR_Zmm_nCentralLightJets";
vect_plot.push_back(plot_);

// Plot: nCentralBJets
// selection region
plot_.root_get_member = "nCentralBJets";
plot_.root_get_factor = "";
plot_.signal_region = "VR_Zmm";
plot_.request_cuts = "";
plot_.blind_data = "";
// plot setup
plot_.x_label = "N_{b-jets}";
plot_.y_label = "Events";
plot_.x_bin_width = 1.0;
plot_.x_range_min = 0.0;
plot_.x_range_max = 10.0;
plot_.y_range_min = 1.0e-1;
plot_.y_range_max = 1.0e9;
plot_.long_name = "VR_Zmm_nCentralBJets";
vect_plot.push_back(plot_);

// VR_Zmm
// Plot: MT2
plot_.root_get_member = "MT2";
plot_.root_get_factor = "";
plot_.signal_region = "VR_Zmm";
plot_.request_cuts = "";
plot_.blind_data = "";
plot_.x_label = "m_{T2} [GeV]";
plot_.y_label = "Events / 10 GeV";
plot_.x_bin_width = 5.0;
plot_.x_range_min = 0.0;
plot_.x_range_max = 200.0;
plot_.y_range_min = 1.0e-1;
plot_.y_range_max = 1.0e9;
plot_.long_name = "VR_Zmm_MT2";
vect_plot.push_back(plot_);

// VR_Zmm
// Plot: m_eff
plot_.root_get_member = "lept1Pt + lept2Pt + met";
plot_.root_get_member_name = "m_eff";
plot_.root_get_factor = "";
plot_.signal_region = "VR_Zmm";
plot_.request_cuts = "";
plot_.blind_data = "";
plot_.x_label = "m_{lepeff} [GeV]";
plot_.y_label = "Events / 25 GeV";
plot_.x_bin_width = 25.0;
plot_.x_range_min = 0.0;
plot_.x_range_max = 1200.0;
plot_.y_range_min = 1.0e-1;
plot_.y_range_max = 1.0e9;
plot_.long_name = "VR_Zmm_m_eff";
vect_plot.push_back(plot_);

plot_.root_get_member_name = "";
plot_.root_get_member_name = "";
plot_.root_get_member_name = "";
































plot_.root_get_member_name = ""; // critical

plot_.is_Log = true;

// VR_Zee
// VR_Zee
// VR_Zee

// VR_Zee
// Plot: ll_M
plot_.root_get_member = "ll_M";
plot_.root_get_factor = "";
plot_.signal_region = "VR_Zee";
plot_.request_cuts = "";
plot_.blind_data = "";
plot_.x_label = "m_{ll} [GeV]";
plot_.y_label = "Events / 25 GeV";
plot_.x_bin_width = 25.0;
plot_.x_range_min = 0.0;
plot_.x_range_max = 750.0;
plot_.y_range_min = 1.0e-1;
plot_.y_range_max = 1.0e9;
plot_.long_name = "VR_Zee_ll_M";
vect_plot.push_back(plot_);

// VR_Zee
// Plot: ll_Pt
plot_.root_get_member = "ll_Pt";
plot_.root_get_factor = "";
plot_.signal_region = "VR_Zee";
plot_.request_cuts = "";
plot_.blind_data = "";
plot_.x_label = "p_{T} (ll) [GeV]";
plot_.y_label = "Events / 15 GeV";
plot_.x_bin_width = 15.0;
plot_.x_range_min = 0.0;
plot_.x_range_max = 500.0;
plot_.y_range_min = 1.0e-1;
plot_.y_range_max = 1.0e9;
plot_.long_name = "VR_Zee_ll_Pt";
vect_plot.push_back(plot_);

// VR_Zee
// Plot: met
plot_.root_get_member = "lept1Pt";
plot_.root_get_factor = "";
plot_.signal_region = "VR_Zee";
plot_.request_cuts = "";
plot_.blind_data = "";
plot_.x_label = "p_{T, l1} [GeV]";
plot_.y_label = "Events / 20 GeV";
plot_.x_bin_width = 20.0;
plot_.x_range_min = 0.0;
plot_.x_range_max = 500.0;
plot_.y_range_min = 1.0e-1;
plot_.y_range_max = 1.0e9;
plot_.long_name = "VR_Zee_lept1Pt";
vect_plot.push_back(plot_);

// VR_Zee
// Plot: met
plot_.root_get_member = "lept2Pt";
plot_.root_get_factor = "";
plot_.signal_region = "VR_Zee";
plot_.request_cuts = "";
plot_.blind_data = "";
plot_.x_label = "p_{T, l2} [GeV]";
plot_.y_label = "Events / 20 GeV";
plot_.x_bin_width = 20.0;
plot_.x_range_min = 0.0;
plot_.x_range_max = 500.0;
plot_.y_range_min = 1.0e-1;
plot_.y_range_max = 1.0e9;
plot_.long_name = "VR_Zee_lept2Pt";
vect_plot.push_back(plot_);

// VR_Zee
// Plot: met
plot_.root_get_member = "met";
plot_.root_get_factor = "";
plot_.signal_region = "VR_Zee";
plot_.request_cuts = "";
plot_.blind_data = "";
plot_.x_label = "E_{T}^{miss} [GeV]";
plot_.y_label = "Events / 10 GeV";
plot_.x_bin_width = 10.0;
plot_.x_range_min = 0.0;
plot_.x_range_max = 300.0;
plot_.y_range_min = 1.0e-1;
plot_.y_range_max = 1.0e9;
plot_.long_name = "VR_Zee_met";
vect_plot.push_back(plot_);

// VR_Zee
// Plot: ll_deltaPhi
plot_.root_get_member = "acos(cos(lept1Phi - lept2Phi))";
plot_.root_get_member_name = "ll_deltaPhi";
plot_.root_get_factor = "";
plot_.signal_region = "VR_Zee";
plot_.request_cuts = "";
plot_.blind_data = "";
plot_.x_label = "#Delta#phi (ll)";
plot_.y_label = "Events";
plot_.x_bin_width = Pi / 24;
plot_.x_range_min = 0.0;
plot_.x_range_max = 1.5 * Pi;
plot_.y_range_min = 1.0e-1;
plot_.y_range_max = 1.0e9;
plot_.long_name = "VR_Zee_ll_deltaPhi";
vect_plot.push_back(plot_);
plot_.root_get_member_name = "";
plot_.root_get_member_name = "";
plot_.root_get_member_name = "";

// Plot: nCentralLightJets
// selection region
plot_.root_get_member = "nCentralLightJets";
plot_.root_get_factor = "";
plot_.signal_region = "VR_Zee";
plot_.request_cuts = "";
plot_.blind_data = "";
// plot setup
plot_.x_label = "N_{jets}-central";
plot_.y_label = "Events";
plot_.x_bin_width = 1.0;
plot_.x_range_min = 0.0;
plot_.x_range_max = 10.0;
plot_.y_range_min = 1.0e-1;
plot_.y_range_max = 1.0e9;
plot_.long_name = "VR_Zee_nCentralLightJets";
vect_plot.push_back(plot_);

// Plot: nForwardJets
// selection region
plot_.root_get_member = "nForwardJets";
plot_.root_get_factor = "";
plot_.signal_region = "VR_Zee";
plot_.request_cuts = "";
plot_.blind_data = "";
// plot setup
plot_.x_label = "N_{for-jets}";
plot_.y_label = "Events";
plot_.x_bin_width = 1.0;
plot_.x_range_min = 0.0;
plot_.x_range_max = 10.0;
plot_.y_range_min = 1.0e-1;
plot_.y_range_max = 1.0e9;
plot_.long_name = "VR_Zee_nForwardJets";
vect_plot.push_back(plot_);

// Plot: nCentralBJets
// selection region
plot_.root_get_member = "nCentralBJets";
plot_.root_get_factor = "";
plot_.signal_region = "VR_Zee";
plot_.request_cuts = "";
plot_.blind_data = "";
// plot setup
plot_.x_label = "N_{b-jets}";
plot_.y_label = "Events";
plot_.x_bin_width = 1.0;
plot_.x_range_min = 0.0;
plot_.x_range_max = 10.0;
plot_.y_range_min = 1.0e-1;
plot_.y_range_max = 1.0e9;
plot_.long_name = "VR_Zee_nCentralBJets";
vect_plot.push_back(plot_);

// VR_Zee
// Plot: MT2
plot_.root_get_member = "MT2";
plot_.root_get_factor = "";
plot_.signal_region = "VR_Zee";
plot_.request_cuts = "";
plot_.blind_data = "";
plot_.x_label = "m_{T2} [GeV]";
plot_.y_label = "Events / 10 GeV";
plot_.x_bin_width = 5.0;
plot_.x_range_min = 0.0;
plot_.x_range_max = 200.0;
plot_.y_range_min = 1.0e-1;
plot_.y_range_max = 1.0e9;
plot_.long_name = "VR_Zee_MT2";
vect_plot.push_back(plot_);

// VR_Zee
// Plot: m_eff
plot_.root_get_member = "lept1Pt + lept2Pt + met";
plot_.root_get_member_name = "m_eff";
plot_.root_get_factor = "";
plot_.signal_region = "VR_Zee";
plot_.request_cuts = "";
plot_.blind_data = "";
plot_.x_label = "m_{lepeff} [GeV]";
plot_.y_label = "Events / 25 GeV";
plot_.x_bin_width = 25.0;
plot_.x_range_min = 0.0;
plot_.x_range_max = 1200.0;
plot_.y_range_min = 1.0e-1;
plot_.y_range_max = 1.0e9;
plot_.long_name = "VR_Zee_m_eff";
vect_plot.push_back(plot_);

plot_.root_get_member_name = "";
plot_.root_get_member_name = "";
plot_.root_get_member_name = "";