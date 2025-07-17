#include <TFile.h>
#include <TDirectory.h>
#include <TKey.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TTree.h>
#include <TLeaf.h>
#include <iostream>
#include <string>

void getERPlots_compMC3() {
  // --- 1) Open the 2D-template file and grab the two histograms ---
  TFile* fHist = TFile::Open("ev_nufit_binned_noprodhav.root","READ");
  if(!fHist || fHist->IsZombie()) {
    std::cerr<<"Cannot open in.root\n"; return;
  }
  TH2D* h2_e = (TH2D*)fHist->Get("hTrueCosineZ_TrueNeutrinoEnergynueselec");
  TH2D* h2_mu= (TH2D*)fHist->Get("hTrueCosineZ_TrueNeutrinoEnergynumuselec");
  if(!h2_e||!h2_mu){
    std::cerr<<"Missing one of the 2D hists\n"; return;
  }

  // --- 2) Make projections (they inherit the exact binning) ---
  TH1D* projE_e  = h2_e->ProjectionY("projE_e");  projE_e->SetTitle("E_{#nu} (#nu_{e} select)");
  TH1D* projZ_e  = h2_e->ProjectionX("projZ_e");  projZ_e->SetTitle("cos#theta_Z (#nu_{e} select)");
  TH1D* projE_mu = h2_mu->ProjectionY("projE_mu"); projE_mu->SetTitle("E_{#nu} (#nu_{#mu} select)");
  TH1D* projZ_mu = h2_mu->ProjectionX("projZ_mu"); projZ_mu->SetTitle("cos#theta_Z (#nu_{#mu} select)");

  // --- 3) Open the TTree file ---
  TFile* fTree = TFile::Open(
    "../outputs/gundamFitter_config_DUNE_Asimov_DryRun.root","READ"
  );
  if(!fTree||fTree->IsZombie()){
    std::cerr<<"Cannot open Asimov file\n"; return;
  }

  // --- 4) Clone the projections to make empty tree‐filled histos (same binning) ---
  TH1D* hE_e_tree  = (TH1D*)projE_e->Clone("hE_e_tree");   hE_e_tree->SetTitle("E_{#nu} (Tree #nu_{e})");  hE_e_tree->Reset();
  TH1D* hZ_e_tree  = (TH1D*)projZ_e->Clone("hZ_e_tree");   hZ_e_tree->SetTitle("cos#theta_Z (Tree #nu_{e})"); hZ_e_tree->Reset();
  TH1D* hE_mu_tree = (TH1D*)projE_mu->Clone("hE_mu_tree"); hE_mu_tree->SetTitle("E_{#nu} (Tree #nu_{#mu})"); hE_mu_tree->Reset();
  TH1D* hZ_mu_tree = (TH1D*)projZ_mu->Clone("hZ_mu_tree"); hZ_mu_tree->SetTitle("cos#theta_Z (Tree #nu_{#mu})");hZ_mu_tree->Reset();

  TDirectory* histDir = dynamic_cast<TDirectory*>(
        fTree->Get("FitterEngine/preFit/model")
    );
    if (!histDir) {
        std::cerr << "Could not find directory: FitterEngine/preFit/model/" << std::endl;
        fTree->Close();
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

      // fill only CC+NC for νe or νμ as appropriate
      if(isNuE){
        hE_e_tree ->Fill(E,w);
        hZ_e_tree ->Fill(cz,w);
      }
      if(isNuMu){
        hE_mu_tree->Fill(E,w);
        hZ_mu_tree->Fill(cz,w);
      }
    }
  }

  // --- 7) Draw overlays and save to PDF ---
  /*TCanvas c("c","c",800,600);
  c.Print("compare.pdf[");

  auto drawOverlay = [&](TH1D* h1, TH1D* h2, bool range){
    h1->SetLineColor(kRed);   h1->SetLineWidth(2);
    h2->SetLineColor(kBlue);  h2->SetLineWidth(2);
    gPad->Clear();
    if (range){
      h1->GetXaxis()->SetRangeUser(0,5);
      //h2->GetXaxis()->SetRangeUser(0,10);
      //h1->GetXaxis()->LogX();
      //gPad->SetLogx(1);
      //h2->GetXaxis()->SetRangeUser(0.2,10);
    }
    h1->Draw("hist");
    h2->Draw("hist same");

    auto leg = new TLegend(0.75, 0.55, 0.98, 0.38);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->AddEntry(h1, "MaCh3", "l");
    leg->AddEntry(h2, "GUNDAM", "l");
    leg->Draw();

    c.Print("compare.pdf");
  };

  c.cd();
  drawOverlay(projZ_e,  hZ_e_tree, 0);
  drawOverlay(projZ_mu, hZ_mu_tree, 0);

  //projE_e->GetXaxis()->SetRangeUser(0,10);
  //projE_mu->GetXaxis()->SetRangeUser(0,10);

  drawOverlay(projE_e,  hE_e_tree, 1);
  drawOverlay(projE_mu, hE_mu_tree, 1);

  c.Print("compare.pdf]");*/

    // --- Prepare multi‐page PDF and big canvas ---
  TCanvas c("c","c",800,800);
  c.Print("compare_with_ratio.pdf[");

  // --- Helper: draw overlay + ratio ---
  auto drawOverlayWithRatio = [&](TH1D* hRef, TH1D* hTest, bool applyRange){
    // clear previous pads
    c.Clear();
    gStyle->SetOptStat(0);

    // --- Top pad: overlays ---
    TPad* pad1 = new TPad("pad1","",0,0.3,1,1);
    pad1->SetBottomMargin(0.02);
    pad1->Draw();
    pad1->cd();

    if (applyRange) {
      //hRef ->GetXaxis()->SetRangeUser(0,5);
      //hTest->GetXaxis()->SetRangeUser(0,5);
      gPad->SetLogx(1);
    }

    hRef ->SetLineColor(kRed);    hRef ->SetLineWidth(2);
    hTest->SetLineColor(kBlue);   hTest->SetLineWidth(2);

    hRef ->Draw("hist");
    hTest->Draw("hist same");

    TLegend leg(0.65,0.7,0.88,0.88);
    leg.SetBorderSize(0);
    leg.SetFillStyle(0);
    leg.AddEntry(hRef,  "MaCh3", "l");
    leg.AddEntry(hTest, "GUNDAM","l");
    leg.Draw();

    // --- Bottom pad: ratio ---
    c.cd();
    TPad* pad2 = new TPad("pad2","",0,0.0,1,0.3);
    pad2->SetTopMargin(0.03);
    pad2->SetBottomMargin(0.3);
    pad2->Draw();
    pad2->cd();

    // build ratio
    TH1D* hRatio = (TH1D*)hTest->Clone("hRatio");
    hRatio->SetTitle("");            // ← remove the title banner
    hRatio->Divide(hRef);

    // style ratio plot
    hRatio->SetMarkerStyle(21);
    hRatio->SetMarkerSize(0.8);
    hRatio->GetYaxis()->SetTitle("GUNDAM/​MaCh3");
    hRatio->GetYaxis()->SetNdivisions(505);
    hRatio->GetYaxis()->SetTitleSize(0.12);
    hRatio->GetYaxis()->SetTitleOffset(0.4);
    hRatio->GetYaxis()->SetLabelSize(0.08);
    hRatio->GetXaxis()->SetTitleSize(0.12);
    hRatio->GetXaxis()->SetLabelSize(0.10);

    if (applyRange) {
      //hRef ->GetXaxis()->SetRangeUser(0,5);
      //hTest->GetXaxis()->SetRangeUser(0,5);
      gPad->SetLogx(1);
    }

    // unity line
    TF1 one("one","1",hRatio->GetXaxis()->GetXmin(),hRatio->GetXaxis()->GetXmax());
    one.SetLineStyle(2);
    one.SetLineColor(kBlack);

    hRatio->Draw("ep");
    one.Draw("same");

    // return to main canvas and print
    c.cd();
    c.Print("compare_with_ratio.pdf");

    // cleanup
    delete pad1;
    delete pad2;
    delete hRatio;
  };

  // --- Draw each pair ---
  drawOverlayWithRatio(projZ_e,  hZ_e_tree,  false);
  drawOverlayWithRatio(projZ_mu, hZ_mu_tree, false);
  drawOverlayWithRatio(projE_e,  hE_e_tree,  true);
  drawOverlayWithRatio(projE_mu, hE_mu_tree, true);

  // --- Close PDF ---
  c.Print("compare_with_ratio.pdf]");

  

// --- Dump template binning to the console ---
  std::cout<<"\nνₑ template cosΘₙ (X-axis) binning:\n";
  {
    TAxis* ax = h2_e->GetXaxis();
    int    n  = ax->GetNbins();
    for(int i=1;i<=n;++i){
      std::cout<<"  ["<<ax->GetBinLowEdge(i)<<", "<<ax->GetBinUpEdge(i)<<"]\n";
    }
  }
  std::cout<<"\nνₑ template Eₙᵤ (Y-axis) binning:\n";
  {
    TAxis* ay = h2_e->GetYaxis();
    int    m  = ay->GetNbins();
    for(int j=1;j<=m;++j){
      std::cout<<"  ["<<ay->GetBinLowEdge(j)<<", "<<ay->GetBinUpEdge(j)<<"]\n";
    }
  }

  std::cout<<"\nν_μ template cosΘₙ (X-axis) binning:\n";
  {
    TAxis* ax = h2_mu->GetXaxis();
    int    n  = ax->GetNbins();
    for(int i=1;i<=n;++i){
      std::cout<<"  ["<<ax->GetBinLowEdge(i)<<", "<<ax->GetBinUpEdge(i)<<"]\n";
    }
  }
  std::cout<<"\nν_μ template Eₙᵤ (Y-axis) binning:\n";
  {
    TAxis* ay = h2_mu->GetYaxis();
    int    m  = ay->GetNbins();
    for(int j=1;j<=m;++j){
      std::cout<<"  ["<<ay->GetBinLowEdge(j)<<", "<<ay->GetBinUpEdge(j)<<"]\n";
    }
  }

  // cleanup
  fHist->Close();
  fTree->Close();
}
