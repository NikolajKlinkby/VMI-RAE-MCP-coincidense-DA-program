#include <iostream>
#include <TFile.h>
#include <TTree.h>
#include <TObject.h>
#include <cstring>
#include <vector>
#include <fstream>

int main(int argc, char *argv[]){

    std::string doc_name;

    if (argc > 1) {
         doc_name = argv[1];
    }
    else {
        std::cout << "ERROR" << std::endl;
        std::cout << "Check you specified an argument" << std::endl;
        std::cout << "./Main <data_file_location>" << std::endl;
        assert(true);
    }

    std::ifstream file; file.open(doc_name);

    Double_t                        mm_conversion = 48.979591837;
    Int_t                           no_coin_count;
    Int_t                           injection;
    Bool_t                          chirp;
    Double_t                        frequency;
    Double_t                        power;
    Bool_t                          shutter;
    std::vector<Double_t>           X;
    std::vector<Double_t>           Y;
    std::vector<Double_t>           trig_time;
    std::vector<Double_t>           event_time;
    std::vector<Double_t>           peak_sum;

    Double_t                        max_val;

    // Setting up root file

    auto root_file = doc_name; root_file += ".root";
    auto * data_file = new TFile(root_file.c_str(), "recreate");
    std::cout << data_file->IsOpen() << std::endl;
    auto * tree = new TTree("DATA", "DATA FRAME");
    auto * max_tree = new TTree("MAXVALUES", "MAX VALUES");
    tree->SetAutoSave(0);

    // Set up the branches of the root tree
    tree->Branch("no_coin_count", & no_coin_count);
    tree->Branch("injection", & injection);
    tree->Branch("chirp", & chirp);
    tree->Branch("frequency",& frequency);
    tree->Branch("power", & power);
    tree->Branch("shutter", & shutter);
    tree->Branch("X", & X);
    tree->Branch("Y", & Y);
    tree->Branch("trig_time", & trig_time);
    tree->Branch("event_time", & event_time);
    tree->Branch("peak_sum", & peak_sum);
    max_tree->Branch("max_values", & max_val);

    // Checks for maximum values x,y,r, event time, min_x, min_y, peak_sum
    std::vector<double> max(7,0); max.at(4) = 100; max.at(5) = 100;
    double curr_r;

    // Load the data into file
    bool injection_data = false;
    std::string line;
    std::string line_sub;
    std::string data_delim = "\t";
    size_t pos;
    int count;
    while (getline(file,line)) {

        // Search for no coin count
        if (line.find("non coin count: ", 0) != std::string::npos){
            // If no coin count is reached, end of previous data range
            injection_data = false;

            // Check if it is the first occurrence of no coin count
            if (!X.empty()) {
                // Fill the data to the tree
                tree->Fill();
                // Reset data addresses
                X.clear();
                Y.clear();
                trig_time.clear();
                event_time.clear();
                peak_sum.clear();
            }
            //Set no coin count
            line.erase(0,16);
            no_coin_count = std::stoi(line);
        }

        // If data section is reached, fill the data.
        else if (injection_data) {
            if (line.empty()) {
                // EOF
                break;
            } else {
                count = 0;
                //Dividing the data and filling it in vectors
                while ((pos = line.find(data_delim)) != std::string::npos) {
                    if (count == 0) {
                        line_sub = line.substr(0, pos);
                        X.emplace_back(std::stod(line_sub)*mm_conversion);
                        count++;
                    }
                    else if (count == 1) {
                        line_sub = line.substr(0, pos);
                        Y.emplace_back(std::stod(line_sub)*mm_conversion);
                        count++;
                        // Assuming both Y and X now has a values,
                        if ( X.at(X.size()-1) > max.at(0)){
                            max.at(0) = X.at(X.size()-1);
                        }
                        if ( X.at(X.size()-1) != 0 && X.at(X.size()-1) < max.at(4)){
                            max.at(4) = X.at(X.size()-1);
                        }
                        if ( Y.at(Y.size()-1) > max.at(1)){
                            max.at(1) = Y.at(Y.size()-1);
                        }
                        if ( Y.at(Y.size()-1) != 0 && Y.at(Y.size()-1) < max.at(5)){
                            max.at(5) = Y.at(Y.size()-1);
                        }
                        curr_r = pow(X.at(X.size()-1)*X.at(X.size()-1)+Y.at(Y.size()-1)*Y.at(Y.size()-1),1./2.);
                        if ( curr_r > max.at(2)){
                            max.at(2) = curr_r;
                        }
                    }
                    else if (count == 2) {
                        line_sub = line.substr(0, pos);
                        trig_time.emplace_back(std::stod(line_sub));
                        count++;
                    }
                    else if ( count == 3) {
                        line_sub = line.substr(0, pos);
                        event_time.emplace_back(std::stod(line_sub));
                        if (std::stod(line_sub) > max.at(3)) {
                            max.at(3) = std::stod(line_sub);
                        }
                    }
                    line.erase(0,pos+data_delim.length());
                }
                peak_sum.emplace_back(std::stod(line));
                if (std::stod(line) > max.at(6)) {
                    max.at(6) = std::stod(line);
                }
            }
        }

        // Search for injection number
        else if (line.find("injection #:", 0) != std::string::npos) {
            line.erase(0, 12);
            injection = std::stoi(line);
        }

        // Search for chirp on/off
        else if (line.find("chirp on/off:", 0) != std::string::npos) {
            line.erase(0, 13);
            chirp = (bool) std::stoi(line);
        }

        // Search for laser frequency and power
        else if (line.find("Frequency (THz): ", 0) != std::string::npos) {
            while ((pos = line.find("    ")) != std::string::npos) {
                line_sub = line.substr(0, pos);
                line_sub.erase(0, 17);
                frequency = std::stod(line_sub);
                line.erase(0, pos + data_delim.length());
            }
            line.erase(0, 16);
            power = std::stod(line);
        }

        // Search for shutter status
        else if (line.find("Shutter status:", 0) != std::string::npos) {
                line.erase(0, 16);
                shutter = (bool) std::stoi(line);
        }

        // If we reached the DATA section set bool to true
        else if (line.find("DATA:", 0) != std::string::npos) {
            injection_data = true;
        }
    }

    for (auto & ent : max) {
        max_val = ent;
        max_tree->Fill();
    }

    // Close down
    //tree->Write(nullptr,TObject::kOverwrite);
    //max_tree->Write(nullptr,TObject::kOverwrite);
    data_file->Write();
    data_file->Close();
    file.close();

    return 1;
}