plot_.root_get_member_name = ""; // critical

//non-log// plot_.is_Log = false;
//non-log// 
//non-log// // VR_MuMu
//non-log// // Plot: met
//non-log// plot_.root_get_member = "lept1ptcone30";
//non-log// plot_.root_get_factor = "";
//non-log// plot_.signal_region = "VR_MuMu";
//non-log// plot_.request_cuts = "";
//non-log// plot_.blind_data = "";
//non-log// plot_.x_label = "p_{T, cone}^{30, l1} [GeV]";
//non-log// plot_.y_label = "Events / 0.5 GeV";
//non-log// plot_.x_bin_width = 0.5;
//non-log// plot_.x_range_min = 0.0;
//non-log// plot_.x_range_max = 20.0;
//non-log// plot_.y_range_min = 0.0;
//non-log// plot_.y_range_max = 4.0e5;
//non-log// plot_.long_name = "VR_MuMu_lept1ptcone30";
//non-log// vect_plot.push_back(plot_);
//non-log// 
//non-log// 
//non-log// // VR_MuMu
//non-log// // Plot: met
//non-log// plot_.root_get_member = "lept2ptcone30";
//non-log// plot_.root_get_factor = "";
//non-log// plot_.signal_region = "VR_MuMu";
//non-log// plot_.request_cuts = "";
//non-log// plot_.blind_data = "";
//non-log// plot_.x_label = "p_{T, cone}^{30, l2} [GeV]";
//non-log// plot_.y_label = "Events / 0.5 GeV";
//non-log// plot_.x_bin_width = 0.5;
//non-log// plot_.x_range_min = 0.0;
//non-log// plot_.x_range_max = 20.0;
//non-log// plot_.y_range_min = 0.0;
//non-log// plot_.y_range_max = 4.0e5;
//non-log// plot_.long_name = "VR_MuMu_lept2ptcone30";
//non-log// vect_plot.push_back(plot_);


plot_.is_Log = true;

// VR_MuMu
// Plot: met
plot_.root_get_member = "lept1ptcone30";
plot_.root_get_factor = "";
plot_.signal_region = "VR_MuMu";
plot_.request_cuts = "";
plot_.blind_data = "";
plot_.x_label = "p_{T, cone}^{30, l1} [GeV]";
plot_.y_label = "Events / 50 GeV";
plot_.x_bin_width = 50.0;
plot_.x_range_min = 0.0;
plot_.x_range_max = 2000.0;
plot_.y_range_min = 0.01;
plot_.y_range_max = 1.0e8;
plot_.long_name = "VR_MuMu_lept1ptcone30_log";
vect_plot.push_back(plot_);


// VR_MuMu
// Plot: met
plot_.root_get_member = "lept2ptcone30";
plot_.root_get_factor = "";
plot_.signal_region = "VR_MuMu";
plot_.request_cuts = "";
plot_.blind_data = "";
plot_.x_label = "p_{T, cone}^{30, l2} [GeV]";
plot_.y_label = "Events / 50 GeV";
plot_.x_bin_width = 50.0;
plot_.x_range_min = 0.0;
plot_.x_range_max = 2000.0;
plot_.y_range_min = 0.01;
plot_.y_range_max = 1.0e8;
plot_.long_name = "VR_MuMu_lept2ptcone30_log";
vect_plot.push_back(plot_);