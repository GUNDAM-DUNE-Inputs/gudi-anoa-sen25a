#include <filesystem>
#include <regex>
#include <iostream>
#include <string>
#include <vector>
#include <map>

#include <TFile.h>
#include <TTree.h>
#include <TGraph.h>
#include <TCanvas.h>

namespace fs = std::filesystem;


// extract the parameter name, e.g. "sin12"
std::string extractParamName(const std::string& fn) {
    const std::string marker = "_With_";
    auto pos = fn.find(marker);
    if (pos == std::string::npos) {
        throw std::runtime_error("No '_With_' in filename: " + fn);
    }
    size_t start = pos + marker.size();
    size_t end = start;
    // consume letters and digits only
    while (end < fn.size() && std::isalnum(static_cast<unsigned char>(fn[end]))) {
        ++end;
    }
    if (end == start) {
        throw std::runtime_error("Empty parameter name in filename: " + fn);
    }
    return fn.substr(start, end - start);
}

// Extract the underscore-encoded prior, e.g. "0_00261" -> 0.00261,
// "7_41e_05" -> 7.41e-05, or "__2_433" -> -2.433
// Recognizes a leading double underscore as a negative sign.
double extractPrior(const std::string& fn) {
    // Regex captures the delimiter underscores and the numeric part before "_fixed"
    static const std::regex re(R"(^.*_With_[a-z0-9]+(_+)([0-9eE_]+)_toysForProfiling_Asimov_ToyFit_0_InjToyPar_parameterInjector.*$)");
    std::smatch m;
    if (!std::regex_match(fn, m, re) || m.size() < 3) {
        throw std::runtime_error("Bad filename for prior: " + fn);
    }
    std::string delim = m[1];       // one or more underscores
    std::string core  = m[2];       // e.g. "0_00261", "7_41e_05", or "2_433"

    bool negative = (delim.size() > 1);

    // handle scientific notation if present
    size_t posE = core.find_first_of("eE");
    double val;
    if (posE != std::string::npos) {
        std::string mant = core.substr(0, posE);
        std::string expPart = core.substr(posE + 1);
        std::replace(mant.begin(), mant.end(), '_', '.');
        if (!expPart.empty() && expPart[0] == '_') {
            expPart = '-' + expPart.substr(1);
        }
        std::string combined = mant + 'e' + expPart;
        val = std::stod(combined);
    } else {
        std::replace(core.begin(), core.end(), '_', '.');
        val = std::stod(core);
    }
    return negative ? -val : val;
}

int plotOApar() {
    const std::string rootDir = "../../outputs/profiling/";
    std::map<std::string, std::vector<double>> chi2map, priormap;

    // 1) scan directory
    for (auto& entry : fs::directory_iterator(rootDir)) {
        if (!entry.is_regular_file()) continue;
        std::string fn = entry.path().filename().string();
        //std::cout<<fn<<endl;

            auto param = extractParamName(fn);   // e.g. "sin12"
            auto prior = extractPrior(fn);       // e.g. 0.00241

            //std::cout<<"param: "<<param<<" prior: "<<prior<<endl;

            // open root file
            TFile* f = TFile::Open((rootDir+fn).c_str());
            //std::cout<<"Opening: "<<(rootDir+fn).c_str()<<endl;

            if (!f || f->IsZombie()) { f->Close(); continue; }
            auto* tFit = (TTree*)f->Get("FitterEngine/postFit/bestFitStats");
            if (!tFit) { f->Close(); continue; }

            tFit->GetEntry(0);
            int fitOK = tFit->GetLeaf("fitStatusCode")->GetValue(0);
            int covOK = tFit->GetLeaf("covStatusCode")->GetValue(0);
            if (fitOK!=0 || covOK!=3) { f->Close(); continue; }

            double chi2 = tFit->GetLeaf("totalLikelihoodAtBestFit")->GetValue(0);
            f->Close();

            // bucket it
            priormap[param].push_back(prior);
            chi2map[param].push_back(chi2);
            std::cout<<param<<"- prior: "<<prior<<" chi2: "<<chi2<<std::endl;
    }
    //sin12, sin13, sin23, m21,m32, dcp
    //0.303, 0.02225, 0.452, 7.41E-5, 2.51E-3, -2.233 //priors
    //0.013, 0.0007, 0.021, 0.0007, 1.8e-6, 3.5e-5, 1.17      // sigma

    double minX[6] = {0.24, 0.01875, 0.347, 6.51e-5, 0.002335, -3}; 
    double maxX[6] = {0.38, 0.02575, 0.557,8.31e-5, 0.002685, 3};

    int count = -1;
    // Plot each parameter into its own PDF
    for (const auto& kv : priormap) {
        count++;
        const std::string& param = kv.first;
        auto priors = kv.second;               // make a mutable copy
        auto chi2s  = chi2map[param];          // mutable copy
        int n = priors.size();
        if (n == 0) continue;

        // Zip priors and chi2s, sort by priors, then unzip
        std::vector<std::pair<double,double>> data;
        data.reserve(n);
        for (int i = 0; i < n; ++i) data.emplace_back(priors[i], chi2s[i]);
        std::sort(data.begin(), data.end(), [](auto &a, auto &b){ return a.first < b.first; });
        for (int i = 0; i < n; ++i) {
            priors[i] = data[i].first;
            chi2s[i]  = data[i].second;
            std::cout<<"prior2: "<<priors[i]<<" chi22: "<<chi2s[i]<<std::endl;
        }

        // Now plot sorted data
        TGraph graph(n, priors.data(), chi2s.data());
        graph.SetTitle(param.c_str());
        //graph.GetXaxis()->SetRangeUser(minX[count],maxX[count]);
        //if (count == 2) graph.GetYaxis()->SetRangeUser(0, 10);
        graph.GetYaxis()->SetTitle("Chi2");
        graph.SetMarkerStyle(20);

        TCanvas canvas(("c_" + param).c_str(), param.c_str(), 800, 600);
        graph.Draw("APL");
        canvas.SaveAs((param + ".png").c_str());
    }
    return 0;
}