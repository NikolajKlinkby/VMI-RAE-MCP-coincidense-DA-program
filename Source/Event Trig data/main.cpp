#include <iostream>
#include <TFile.h>
#include <TTree.h>
#include <TObject.h>
#include <vector>

int main(int argc, char *argv[]){

    std::string root_file;
    std::vector<double_t> bin_size(4,0);
    std::vector<double> settings_time_vec(7,0);

    bin_size.at(0) = 0.5e-3; //event_time_bin_size
    bin_size.at(1) = 0.5e-8; //trig_time_bin_size
    bin_size.at(2) = 20; //peak_sum_bin_size
    bin_size.at(3) = 0.1e-4; //freq_bin_size

    settings_time_vec.at(0) = 0; // Trigger time lower
    settings_time_vec.at(1) = 1.e-4; // Trigger time upper
    settings_time_vec.at(2) = 0; // Event time lower
    settings_time_vec.at(3) = 1; // Event time upper
    settings_time_vec.at(4) = 0; // Max radius
    settings_time_vec.at(5) = 0; // Peak sum lower
    settings_time_vec.at(6) = 10000; // Peak sum upper

    bool event_time_set = false;
    bool peak_sum_set = false;
    bool trigger_time_changed = false;
    bool event_time_changed = false;
    bool peak_sum_changed = false;
    bool frequency_changed = false;

    char *ptr;

    if (argc > 1) {
        root_file = argv[1];
    }
    else {
        std::cout << "ERROR" << std::endl;
        std::cout << "Check you specified an argument" << std::endl;
        std::cout << "./Main <data_file_location>" << std::endl;
        assert(true);
    }

    auto * data_file = new TFile(root_file.c_str(), "update");

    // If settings already exist, load them.
    bool setting_exist = data_file->GetListOfKeys()->Contains("SETTINGS_BINS") && data_file->GetListOfKeys()->Contains("SETTINGS_TIMES");
    if (setting_exist){
        auto *settings_bin_tree = (TTree *) data_file->Get("SETTINGS_BINS");

        double settings_bins;

        settings_bin_tree->SetBranchAddress("settings_bin", &settings_bins);

        // Setting trigger time, and event time
        for (int i = 0; i < settings_bin_tree->GetEntries(); ++i) {
            settings_bin_tree->GetEntry(i);
            bin_size.at(i) = settings_bins;
        }

        auto *settings_time_tree = (TTree *) data_file->Get("SETTINGS_TIMES");

        double settings_time;

        settings_time_tree->SetBranchAddress("settings_time", &settings_time);

        // Setting trigger time, and event time
        for (int i = 0; i < settings_time_tree->GetEntries(); ++i) {
            settings_time_tree->GetEntry(i);
            settings_time_vec.at(i) = settings_time;
        }
        event_time_set = true;
        peak_sum_set = true;
    }

    // Checking the arguments given and checking if they changed from previous settings.
    argc--;
    argv++;
    int c;
    while ((c = getopt(argc, argv, "t:T:e:E:p:P:f:g:h:j:")) != -1) {
        switch (c) {
            case 't':
                settings_time_vec.at(0) = strtod(optarg, & ptr);
                break;
            case 'T':
                settings_time_vec.at(1) = strtod(optarg, & ptr);
                break;
            case 'e':
                settings_time_vec.at(2) = strtod(optarg, & ptr);
                break;
            case 'E':
                settings_time_vec.at(3) = strtod(optarg, & ptr);
                break;
            case 'p':
                settings_time_vec.at(5) = strtod(optarg, & ptr);
                break;
            case 'P':
                settings_time_vec.at(6) = strtod(optarg, & ptr);
                break;
            case 'f': //Freq bin
                bin_size.at(3) = strtod(optarg, & ptr);
                frequency_changed = true;
                break;
            case 'g': //event bin
                bin_size.at(0) = strtod(optarg, & ptr);
                event_time_changed = true;
                break;
            case 'h': //trigger bin
                bin_size.at(1) = strtod(optarg, & ptr);
                trigger_time_changed = true;
                break;
            case 'j': //Peak sum bin
                bin_size.at(2) = strtod(optarg, & ptr);
                peak_sum_changed = true;
                break;
            default:
                abort();
        }
    }

    // If settings don't exist, create it! or if settings where changed, change them!
    if (!setting_exist || trigger_time_changed || event_time_changed || peak_sum_changed || frequency_changed) {
        auto *settings_bin_tree = new TTree("SETTINGS_BINS","SETTINGS_BINS");

        double settings_bins;

        settings_bin_tree->Branch("settings_bin", &settings_bins, "settings_bin/D");

        // Setting trigger time, and event time
        for (auto &ent: bin_size) {
            settings_bins = ent;
            settings_bin_tree->Fill();
        }
        //Write
        settings_bin_tree->Write(nullptr, TObject::kOverwrite);
    }


    Int_t                           injection;
    Bool_t                          shutter;
    Double_t                        max_val;
    std::vector<Double_t>           * trig_time = nullptr;
    std::vector<Double_t>           * event_time = nullptr;
    std::vector<Double_t>           * peak_sum = nullptr;
    Double_t                        freq;

    TBranch                         * trig_time_branch = nullptr;
    TBranch                         * event_time_branch = nullptr;
    TBranch                         * peak_sum_branch = nullptr;
    TBranch                         * freq_branch = nullptr;
    TBranch                         *max_branch = nullptr;

    // Setting up root file
    auto * tree = (TTree*) data_file->Get("DATA");
    auto *max_tree = (TTree *) data_file->Get("MAXVALUES");
    auto * event_tree = new TTree("EVENT_TOT","EVENT_TOT");
    auto * trig_tree = new TTree("TRIG_TOT","TRIG_TOT");
    auto * peak_tree = new TTree("PEAK_TOT","PEAK_TOT");
    auto * freq_tree = new TTree("FREQ_TOT","FREQ_TOT");
    auto * freq_bin_tree = new TTree("FREQ_BIN","FREQ_BIN");
    auto * event_now_tree = new TTree("EVENT_NOW","EVENT_NOW");
    auto * trig_now_tree = new TTree("TRIG_NOW","TRIG_NOW");
    auto * peak_now_tree = new TTree("PEAK_NOW","PEAK_NOW");

    // Set up the branches of the root tree
    tree->SetBranchAddress("injection", & injection);
    tree->SetBranchAddress("shutter", & shutter);

    tree->SetBranchAddress("trig_time", & trig_time, & trig_time_branch);
    tree->SetBranchAddress("event_time", & event_time, & event_time_branch);
    tree->SetBranchAddress("peak_sum", & peak_sum, & peak_sum_branch);
    tree->SetBranchAddress("frequency", & freq, & freq_branch);
    max_tree->SetBranchAddress("max_values", &max_val, &max_branch); //x, y, r, event time

    std::vector<Double_t> event_bin = {0};
    std::vector<int> event_hist = {0};
    std::vector<Double_t> trig_bin = {0};
    std::vector<int> trig_hist = {0};
    std::vector<Double_t> peak_bin = {0};
    std::vector<int> peak_hist = {0};
    std::vector<int> event_now_hist = {0};
    std::vector<int> trig_now_hist = {0};
    std::vector<int> peak_now_hist = {0};

    int entries = tree->GetEntries();

    // Find a starting point
    std::vector<Double_t> freq_bin = {0};
    std::vector<int> freq_hist = {0};
    freq = -3;
    int count = 0;
    while (freq < 0){
        freq_branch->GetEntry(tree->LoadTree(count));
        count++;
        if (count == entries){
            break;
        }
        freq_bin = {freq};
    }

    double trig_time_lower = settings_time_vec.at(0);
    double trig_time_upper = settings_time_vec.at(1);
    double event_time_lower = settings_time_vec.at(2);
    double event_time_upper = settings_time_vec.at(3);
    double peak_sum_lower = settings_time_vec.at(5);
    double peak_sum_upper = settings_time_vec.at(6);
    // If event_time is not set, look for the maximum time
    if (!event_time_set) {
        max_branch->GetEntry(3);
        event_time_upper = max_val;
    }
    if (!peak_sum_set) {
        max_branch->GetEntry(6);
        peak_sum_upper = max_val;
    }

    Double_t trig_val;
    Double_t event_val;
    Double_t peak_val;
    Double_t freq_val;
    Double_t trig_now_val;
    Double_t event_now_val;
    Double_t peak_now_val;
    Double_t freq_bin_val;

    event_tree->Branch("event_hist", & event_val, "event_hist/D");
    trig_tree->Branch("trig_hist", & trig_val, "trig_hist/D");
    peak_tree->Branch("peak_hist", & peak_val, "peak_hist/D");
    freq_tree->Branch("freq_hist", & freq_val, "freq_hist/D");
    freq_bin_tree->Branch("freq_bin", & freq_bin_val, "freq_bin/D");
    event_now_tree->Branch("event_hist", & event_now_val, "event_hist/D");
    trig_now_tree->Branch("trig_hist", & trig_now_val, "trig_hist/D");
    peak_now_tree->Branch("peak_hist", & peak_now_val, "peak_hist/D");
    
    for(int i = 0; i < entries; i++){
        trig_time_branch->GetEntry(tree->LoadTree(i));
        event_time_branch->GetEntry(tree->LoadTree(i));
        peak_sum_branch->GetEntry(tree->LoadTree(i));
        for (int j = 0; j < trig_time->size(); j++){
            while (trig_time->at(j) >= trig_bin.back()){
                trig_bin.emplace_back(trig_bin.back()+bin_size.at(1));
                trig_hist.emplace_back(0);
                trig_now_hist.emplace_back(0);
            }
            while (event_time->at(j) >= event_bin.back()) {
                event_bin.emplace_back(event_bin.back() + bin_size.at(0));
                event_hist.emplace_back(0);
                event_now_hist.emplace_back(0);
            }
            while (peak_sum->at(j) >= peak_bin.back()){
                peak_bin.emplace_back(peak_bin.back()+bin_size.at(2));
                peak_hist.emplace_back(0);
                peak_now_hist.emplace_back(0);
            }
            auto lower_trig = std::lower_bound(trig_bin.begin(), trig_bin.end(),trig_time->at(j));
            auto lower_event = std::lower_bound(event_bin.begin(), event_bin.end(),event_time->at(j));
            auto lower_peak = std::lower_bound(peak_bin.begin(), peak_bin.end(),peak_sum->at(j));
            trig_hist.at(std::distance(trig_bin.begin(), lower_trig))++;
            event_hist.at(std::distance(event_bin.begin(), lower_event))++;
            peak_hist.at(std::distance(peak_bin.begin(), lower_peak))++;
            if (trig_time->at(j) < trig_time_upper && trig_time->at(j) > trig_time_lower &&
                event_time->at(j) < event_time_upper && event_time->at(j) > event_time_lower &&
                peak_sum->at(j) < peak_sum_upper && peak_sum->at(j) > peak_sum_lower) {
                trig_now_hist.at(std::distance(trig_bin.begin(), lower_trig))++;
                event_now_hist.at(std::distance(event_bin.begin(), lower_event))++;
                peak_now_hist.at(std::distance(peak_bin.begin(), lower_peak))++;
            }
        }
        freq_branch->GetEntry(tree->LoadTree(i));
        if (freq > 0) {
            while ((freq >= freq_bin.back()) && ((freq - freq_bin.back()) < 500*bin_size.at(3))) {
                freq_bin.emplace_back(freq_bin.back() + bin_size.at(3));
                freq_hist.emplace_back(0);
            }
            while ((freq <= freq_bin.at(0)) && (freq_bin.at(0) - freq < 500*bin_size.at(3))) {
                freq_bin.insert(freq_bin.begin(),freq_bin.at(0) - bin_size.at(3));
                freq_hist.insert(freq_hist.begin(),0);
            }
            if ((freq <= freq_bin.back()) && (freq >= freq_bin.at(0))) {
                auto lower_freq = std::lower_bound(freq_bin.begin(), freq_bin.end(), freq);
                freq_hist.at(std::distance(freq_bin.begin(), lower_freq))++;
            }
        }
    }

    for (auto & ent : event_hist) {
        event_val = ent;
        event_tree->Fill();
    }
    for (auto & ent : event_now_hist) {
        event_now_val = ent;
        event_now_tree->Fill();
    }
    for (auto & ent : trig_hist) {
        trig_val = ent;
        trig_tree->Fill();
    }
    for (auto & ent : trig_now_hist) {
        trig_now_val = ent;
        trig_now_tree->Fill();
    }
    for (auto & ent : peak_hist) {
        peak_val = ent;
        peak_tree->Fill();
    }
    for (auto & ent : peak_now_hist) {
        peak_now_val = ent;
        peak_now_tree->Fill();
    }
    for (auto & ent : freq_hist) {
        freq_val = ent;
        freq_tree->Fill();
    }
    for (auto & ent : freq_bin) {
        freq_bin_val = ent;
        freq_bin_tree->Fill();
    }
    event_tree->Write(nullptr, TObject::kOverwrite);
    event_now_tree->Write(nullptr, TObject::kOverwrite);
    trig_tree->Write(nullptr, TObject::kOverwrite);
    peak_tree->Write(nullptr, TObject::kOverwrite);
    trig_now_tree->Write(nullptr, TObject::kOverwrite);
    peak_now_tree->Write(nullptr, TObject::kOverwrite);
    freq_tree->Write(nullptr, TObject::kOverwrite);
    freq_bin_tree->Write(nullptr, TObject::kOverwrite);
    data_file->Close();
    return 0;
}
