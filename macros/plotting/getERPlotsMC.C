#include <TFile.h>
#include <TDirectory.h>
#include <TKey.h>
#include <TH1D.h>
#include <TClass.h>
#include <iostream>

void getERPlotsMC()
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

    //--- 1D histograms (same as before) ---
    TH1D* h_mcEnu_e_CC  = new TH1D("h_mcEnu_e_CC",  "E_{#nu} (#nu_{e} CC);E [GeV];Events", 50, 0, 10);
    TH1D* h_mcEnu_e_NC  = new TH1D("h_mcEnu_e_NC",  "E_{#nu} (#nu_{e} NC);E [GeV];Events", 50, 0, 10);
    TH1D* h_mcEnu_e  = new TH1D("h_mcEnu_e",  "E_{#nu} (#nu_{e});E [GeV];Events", 50, 0, 10);
    TH1D* h_mcEnu_mu_CC = new TH1D("h_mcEnu_mu_CC", "E_{#nu} (#nu_{#mu} CC);E [GeV];Events", 50, 0, 10);
    TH1D* h_mcEnu_mu_NC = new TH1D("h_mcEnu_mu_NC", "E_{#nu} (#nu_{#mu} NC);E [GeV];Events", 50, 0, 10);
    TH1D* h_mcEnu_mu = new TH1D("h_mcEnu_mu", "E_{#nu} (#nu_{#mu});E [GeV];Events", 50, 0, 10);

    TH1D* h_cosZ_e_CC  = new TH1D("h_cosZ_e_CC",  "cos#theta_Z (#nu_{e} CC);cos#theta_Z;Events", 20, -1, 1);
    TH1D* h_cosZ_e_NC  = new TH1D("h_cosZ_e_NC",  "cos#theta_Z (#nu_{e} NC);cos#theta_Z;Events", 20, -1, 1);
    TH1D* h_cosZ_e  = new TH1D("h_cosZ_e",  "cos#theta_Z (#nu_{e});cos#theta_Z;Events", 20, -1, 1);
    TH1D* h_cosZ_mu_CC = new TH1D("h_cosZ_mu_CC", "cos#theta_Z (#nu_{#mu} CC);cos#theta_Z;Events", 20, -1, 1);
    TH1D* h_cosZ_mu_NC = new TH1D("h_cosZ_mu_NC", "cos#theta_Z (#nu_{#mu} NC);cos#theta_Z;Events", 20, -1, 1);
    TH1D* h_cosZ_mu = new TH1D("h_cosZ_mu", "cos#theta_Z (#nu_{#mu});cos#theta_Z;Events", 20, -1, 1);

    //--- 2D histograms (cos#theta_Z vs E) ---
    TH2D* h2_e_CC  = new TH2D("h2_e_CC",  "E vs cos#theta_Z (#nu_{e} CC);E [GeV];cos#theta_Z", 50, 0, 10, 20, -1, 1);
    TH2D* h2_e_NC  = new TH2D("h2_e_NC",  "E vs cos#theta_Z (#nu_{e} NC);E [GeV];cos#theta_Z", 50, 0, 10, 20, -1, 1);
    TH2D* h2_mu_CC = new TH2D("h2_mu_CC", "E vs cos#theta_Z (#nu_{#mu} CC);E [GeV];cos#theta_Z", 50, 0, 10, 20, -1, 1);
    TH2D* h2_mu_NC = new TH2D("h2_mu_NC", "E vs cos#theta_Z (#nu_{#mu} NC);E [GeV];cos#theta_Z", 50, 0, 10, 20, -1, 1);
    TH2D* h2_mu = new TH2D("h2_mu", "E vs cos#theta_Z (#nu_{#mu});E [GeV];cos#theta_Z", 50, 0, 10, 20, -1, 1);
    TH2D* h2_e = new TH2D("h2_e", "E vs cos#theta_Z (#nu_{e});E [GeV];cos#theta_Z", 50, 0, 10, 20, -1, 1);


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
    //double weightSumCC_e [NMODES] = {0.0};  // CC sums for all "#nu_{e} Selec" subDirs
    //double weightSumNC_e [NMODES] = {0.0};  // NC sums for all "#nu_{e} Selec" subDirs
    //double weightSumCC_mu[NMODES] = {0.0};  // CC sums for all "#nu_{#mu} Selec" subDirs
    //double weightSumNC_mu[NMODES] = {0.0};  // NC sums for all "#nu_{#mu} Selec" subDirs


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

        //bool isNonOsc_e = (subDirName.find("#nu_{e} x #nu_{e}") != std::string::npos);
        //bool isNonOsc_ae = (subDirName.find("#bar{#nu_{e}} x #bar{#nu_{e}}") != std::string::npos);
        //bool isNonOsc_mu = (subDirName.find("#nu_{#mu} x #nu_{#mu}") != std::string::npos);
        //bool isNonOsc_amu = (subDirName.find("#bar{#nu_{#mu}} x #bar{#nu_{#mu}}") != std::string::npos);

        //bool isNonOsc = 0;
        //if (isNonOsc_e || isNonOsc_ae || isNonOsc_mu || isNonOsc_amu)  isNonOsc = 1;

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
        Float_t mcEnu = 0.0;
        Float_t mcthetaZ = 0.0;


        // (Use the bare leaf names for split branches.)
        TLeaf* leafMode       = tree->GetLeaf("mcMode");
        TLeaf* leafIsCC       = tree->GetLeaf("mcIsCC");
        //TLeaf* leafGenWeight  = tree->GetLeaf("mcGenWeight");
        TLeaf* leafGenWeight  = tree->GetLeaf("eventWeight");
        TLeaf* leafEnu       = tree->GetLeaf("mcEnu");
        TLeaf* leafCosZ      = tree->GetLeaf("mcthetaZ");


        //double weightSumCC_dir [NMODES] = {0.0};
        //double weightSumNC_dir [NMODES] = {0.0};
        //double allmodes = 0;
        Long64_t nEntries = tree->GetEntries();
        for (Long64_t i = 0; i < nEntries; ++i) {
            tree->GetEntry(i);

            int    isCC  = leafIsCC->GetValue(0);
            double E     = leafEnu ->GetValue(0);
            double cz    = leafCosZ->GetValue(0);
            double w     = leafGenWeight->GetValue(0);

            if (isCC==1 && isNuE) {
                h_mcEnu_e_CC ->Fill(E,w);
                h_cosZ_e_CC  ->Fill(cz,w);
                h2_e_CC      ->Fill(E,cz,w);
            }
            else if(isCC==0 && isNuE) {
                h_mcEnu_e_NC ->Fill(E,w);
                h_cosZ_e_NC  ->Fill(cz,w);
                h2_e_NC      ->Fill(E,cz,w);
            }
            else if(isCC==1 && isNuMu) {
                h_mcEnu_mu_CC->Fill(E,w);
                h_cosZ_mu_CC ->Fill(cz,w);
                h2_mu_CC     ->Fill(E,cz,w);
            }
            else if(isCC==0 && isNuMu) {
                h_mcEnu_mu_NC->Fill(E,w);
                h_cosZ_mu_NC ->Fill(cz,w);
                h2_mu_NC     ->Fill(E,cz,w);
            }

            if (isNuE) {
        
                h_mcEnu_e ->Fill(E,w);
                h_cosZ_e  ->Fill(cz,w);
                h2_e      ->Fill(E,cz,w);

            }

            if (isNuMu) {
        
                h_mcEnu_mu->Fill(E,w);
                h_cosZ_mu ->Fill(cz,w);
                h2_mu     ->Fill(E,cz,w);

            }

        }

    } 
    //--- Produce PDF plots ---
    TCanvas c("c","c",800,600);
    const std::vector<TH1*> h1s = {
        h_mcEnu_e_CC, h_mcEnu_e_NC, h_mcEnu_e, h_mcEnu_mu_CC, h_mcEnu_mu_NC, h_mcEnu_mu,
        h_cosZ_e_CC,  h_cosZ_e_NC, h_cosZ_e,  h_cosZ_mu_CC,  h_cosZ_mu_NC, h_cosZ_mu
    };
    const std::vector<TH2*> h2s = {
        h2_e_CC, h2_e_NC, h2_e, h2_mu_CC, h2_mu_NC, h2_mu
    };

    // open multi-page pdf
    c.Print("all_plots_Osc_MaCh3_fine_binning_v1.pdf(");
    int count = 0;
    for (auto h : h1s) {
        c.cd();
        h->SetLineColor(kRed);
        if (count > 5) h->GetYaxis()->SetRangeUser(700, 2300);

        h->Draw();
        h->Draw("hist,sames");
        c.Print("all_plots_Osc_MaCh3_fine_binning_v1.pdf");
        count++;
    }
    for (auto h2 : h2s) {
        c.cd();
        h2->Draw("COLZ");
        c.Print("all_plots_Osc_MaCh3_fine_binning_v1.pdf");
    }
    // close multi-page pdf
    c.Print("all_plots_Osc_MaCh3_fine_binning_v1.pdf)");

    f->Close();
}

