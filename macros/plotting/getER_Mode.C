#include <TFile.h>
#include <TDirectory.h>
#include <TKey.h>
#include <TH1D.h>
#include <TClass.h>
#include <iostream>

void getER_Mode()
{   
    const char* inputFileName = "../outputs/gundamFitter_config_DUNE_Asimov_DryRun.root";
    // A lookup table from interaction mode number to descriptive text
    // (You can rename or adjust titles to match exactly what you want printed)
    static std::map<int, std::string> modeTitle = {
        {-1, "UNDEFINED"},
        {0,  "QE"},
        {1,  "Resonant"},
        {2,  "DIS"},
        {3,  "Coherent"},
        {4,  "Coherent Elastic"},
        {5,  "Electron Scattering"},
        {6,  "Inverse Muon Decay Annihliation"},
        {7,  "Inverse Beta Decay"},
        {8,  "Glasgow Resonance"},
        {9,  "Atmospheric Muon Nu Gamma"},
        {10, "MEC aka 2p2h"},
        {11, "Diffractive"},
        {12, "kEM"},
        {13, "kWeakMix"}
    };

    string modeTitle_alt[14] = {
        "QE",
        "RES",
        "DIS",
        "COH",
        "COHEL",
        "NuEl",
        "IMDAnn",
        "IBD",
        "GlRES",
        "AnuGam",
        "MEC",
        "DIFF",
        "kEM",
        "kWeakMix"
    };

    // Read Files
    TFile* f = TFile::Open(inputFileName, "READ");
    if (!f || f->IsZombie()) {
        std::cerr << "Error opening file!" << std::endl;
        return;
    }



    // Navigate to the directory: "FitterEngine/preFit/plots/histograms"
    TDirectory* histDir = dynamic_cast<TDirectory*>(
        f->Get("FitterEngine/preFit/model")
    );
    if (!histDir) {
        std::cerr << "Could not find directory: FitterEngine/preFit/model/" << std::endl;
        f->Close();
        return;
    }

    static const int NMODES = 14;
    double weightSumCC_e [NMODES] = {0.0};  // CC sums for all "#nu_{e} Selec" subDirs
    double weightSumNC_e [NMODES] = {0.0};  // NC sums for all "#nu_{e} Selec" subDirs
    double weightSumCC_mu[NMODES] = {0.0};  // CC sums for all "#nu_{#mu} Selec" subDirs
    double weightSumNC_mu[NMODES] = {0.0};  // NC sums for all "#nu_{#mu} Selec" subDirs


    // Loop over the eight subfolders (e.g., "FHC #nu_{#mu}_like CC", "FHC #nu_{#mu}_like NC", etc.)
    TIter nextSubfolder(histDir->GetListOfKeys());
    TKey* keySubfolder;
    while ((keySubfolder = (TKey*)nextSubfolder())) {

        // Only process directories
        if (strcmp(keySubfolder->GetClassName(), "TDirectoryFile") != 0) {
            continue;
        }

        // subDir is something like "FHC #nu_{#mu}_like CC"
        TDirectory* subDir = dynamic_cast<TDirectory*>(histDir->Get(keySubfolder->GetName()));
        if (!subDir){
            std::cout<<"No subDir"<<std::endl;
            continue;
        } 

        std::string subDirName = subDir->GetName();

        bool isNuE  = (subDirName.find("#nu_{e} Selec")  != std::string::npos);
        bool isNuMu = (subDirName.find("#nu_{#mu} Selec") != std::string::npos);

        bool isNonOsc_e = (subDirName.find("#nu_{e} x #nu_{e}") != std::string::npos);
        bool isNonOsc_ae = (subDirName.find("#bar{#nu_{e}} x #bar{#nu_{e}}") != std::string::npos);
        bool isNonOsc_mu = (subDirName.find("#nu_{#mu} x #nu_{#mu}") != std::string::npos);
        bool isNonOsc_amu = (subDirName.find("#bar{#nu_{#mu}} x #bar{#nu_{#mu}}") != std::string::npos);

        bool isNonOsc = 0;
        if (isNonOsc_e || isNonOsc_ae || isNonOsc_mu || isNonOsc_amu)  isNonOsc = 1;

        if (!isNuE && !isNuMu) {
            // Skip any other category
            std::cout<<"no isNuE or isNuMu found"<<std::endl;
            continue;
        }

        //std::cout<<"isNuE: "<<isNuE<<" isNuMu:"<<isNuMu<<std::endl;

        TTree* tree = dynamic_cast<TTree*>(subDir->Get("events_TTree"));
        if (!tree) {
            std::cout << "No TTree named \"events_TTree\" in directory: " 
                      << subDir->GetName() << std::endl;
            continue;
        }

        Int_t    mcIsCC        = -999;
        Int_t    mcMode      = 0;
        Float_t mcGenWeight = 0.0;

        // (Use the bare leaf names for split branches.)
        TLeaf* leafMode       = tree->GetLeaf("mcMode");
        TLeaf* leafIsCC       = tree->GetLeaf("mcIsCC");
        //TLeaf* leafGenWeight  = tree->GetLeaf("mcGenWeight");
        TLeaf* leafGenWeight  = tree->GetLeaf("eventWeight");


        double weightSumCC_dir [NMODES] = {0.0};
        double weightSumNC_dir [NMODES] = {0.0};
        double allmodes = 0;
        Long64_t nEntries = tree->GetEntries();
        for (Long64_t i = 0; i < nEntries; ++i) {
            tree->GetEntry(i);

            // You can read leaf values by casting leaf->GetValue(0) to the correct type:
            Int_t   mode       = static_cast<Int_t>(leafMode->GetValue(0));
            Int_t   isCC       = static_cast<Int_t>(leafIsCC->GetValue(0));
            Float_t weight     = static_cast<Float_t>(leafGenWeight->GetValue(0));
            //Int_t weight     = 1;


            // only accumulate when mcIsCC == 1 and idcc == 0
            if (mode < 0 || mode >= NMODES) continue;
            //std::cout<<"mode: "<<mode<<" isCC: "<<isCC<<" weight: "<<weight<<std::endl;

            // use when raw events or unoscillated
            //if (!isNonOsc) continue;

            if (isCC == 1 && isNuE == 1) {
                weightSumCC_e[mode] += weight;
            } else if (isCC == 1 && isNuMu == 1) {
                weightSumCC_mu[mode] += weight;
            } else if (isCC == 0 && isNuE == 1){
                weightSumNC_e[mode] += weight;
            } else if (isCC == 0 && isNuMu == 1){
                weightSumNC_mu[mode] += weight; 
            }

            if (isCC == 1)
                weightSumCC_dir[mode] += weight;
            else
                weightSumNC_dir[mode] += weight;

            allmodes += weight; 

        }

            for (int modeVal = 0; modeVal < NMODES; ++modeVal) {

                    std::cout << "Directory: " << subDirName
                      << " , Mode: " << modeVal << std::setprecision(10)
                      //<< ", " << outModeTitle
                      << ", Integral CC: " << weightSumCC_dir[modeVal]
                      << ", Integral NC: " << weightSumNC_dir[modeVal]
                      << std::endl;
            }
            
            std::cout << "Integral all modes: " << allmodes << std::endl;

    }

    f->Close();

    std::cout << "\n=== Combined sums for all \"#nu_{e} Selec\" subdirectories ===\n";
    for (int modeVal = 0; modeVal < NMODES; ++modeVal) {
        std::cout << "Mode " << modeVal << " "<< modeTitle_alt [modeVal] << std::setprecision(10)
                  << " : sum(mcGenWeight)[CC] = " << weightSumCC_e[modeVal]
                  << " , sum(mcGenWeight)[NC] = " << weightSumNC_e[modeVal]
                  << "\n";
    }

    std::cout << "\n=== Combined sums for all \"#nu_{#mu} Selec\" subdirectories ===\n";
    for (int modeVal = 0; modeVal < NMODES; ++modeVal) {
        std::cout << "Mode " << modeVal  << " "<< modeTitle_alt [modeVal] << std::setprecision(10)
                  << " : sum(mcGenWeight)[CC] = " << weightSumCC_mu[modeVal]
                  << " , sum(mcGenWeight)[NC] = " << weightSumNC_mu[modeVal]
                  << "\n";
    }

}
