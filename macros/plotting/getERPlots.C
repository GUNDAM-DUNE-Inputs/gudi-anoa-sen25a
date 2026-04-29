#include <TFile.h>
#include <TDirectory.h>
#include <TKey.h>
#include <TTree.h>
#include <TLeaf.h>
#include <TH1D.h>
#include <TH2D.h>
#include <THStack.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TString.h>
#include <TColor.h>

#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <cmath>
#include <cstring>

namespace {

    // =========================================================
    // Configuration
    // =========================================================
    static const int NMODES = 14;

    // Set to true to use custom reco-energy binning separate for numu and nue selections.
    static const bool kUseCustomRecoBinning = true;

    enum class SelectionType {
        NuMu,
        NuE
    };

    // ---------------------------------------------------------
    // Custom reco-energy bin edges for numu selection
    // Edit these as needed.
    // ---------------------------------------------------------
    const std::vector<double> kRecoBinEdgesNuMu = {
        0.0, 0.12, 0.22, 0.37, 0.63,
        1.5, 4.0, 11.0, 19.0,
        32.0, 100.0
    };

    // ---------------------------------------------------------
    // Custom reco-energy bin edges for nue selection
    // Edit these as needed.
    // ---------------------------------------------------------
    const std::vector<double> kRecoBinEdgesNuE = {
        0.0, 0.07, 0.1, 0.15, 0.22,
        0.33, 0.47, 0.7, 1.1,
        1.7, 2.9, 5.0, 9.0, 16.0, 32.0, 48.0, 100.0
    };

    // Fallback uniform binning if kUseCustomRecoBinning = false
    static const int    kRecoNBinsUniform = 50;
    static const double kRecoMinUniform   = 0.0;
    static const double kRecoMaxUniform   = 10.0;

    // Truth-energy binning
    static const int    kTrueNBins = 50;
    static const double kTrueMin   = 0.0;
    static const double kTrueMax   = 10.0;

    // cosZ binning
    static const int    kCosNBins = 100;
    static const double kCosMin   = -1.0;
    static const double kCosMax   =  1.0;

    const char* modeTitle[NMODES] = {
        "QE", "RES", "DIS", "COH", "COHEL", "NuEl", "IMDAnn",
        "IBD", "GlRES", "AnuGam", "MEC", "DIFF", "kEM", "kWeakMix"
    };

    int modeColor[NMODES] = {
        kAzure+1, kOrange+7, kGreen+2, kMagenta+1, kCyan+1, kYellow+2,
        kRed+1, kBlue-7, kSpring+5, kViolet+1, kPink+7, kTeal+3,
        kGray+2, kOrange-3
    };

    // =========================================================
    // Category container
    // =========================================================
    struct Category {
        std::string name;
        std::string label;

        SelectionType selectionType;

        // MC stacked by mode
        std::vector<TH1D*> recoE_mc;
        std::vector<TH1D*> recoCos_mc;
        std::vector<TH1D*> trueE_mc;
        std::vector<TH1D*> trueCos_mc;

        // Data inclusive
        TH1D* recoE_data   = nullptr;
        TH1D* recoCos_data = nullptr;
        TH1D* trueE_data   = nullptr;
        TH1D* trueCos_data = nullptr;

        // 2D histograms
        TH2D* reco2D_mc    = nullptr;   // model reco: cosZ vs E_reco
        TH2D* true2D_mc    = nullptr;   // model true: cosZ vs E_true
        TH2D* reco2D_data  = nullptr;   // data reco:  cosZ vs E_reco

        // Inclusive weighted event rates
        double rate_mc   = 0.0;
        double rate_data = 0.0;

        // Inclusive raw entry counts
        Long64_t n_mc   = 0;
        Long64_t n_data = 0;

        // Per-mode weighted event rates
        std::vector<double> rate_mc_mode;
        std::vector<double> rate_data_mode;

        // Per-mode raw entry counts
        std::vector<Long64_t> n_mc_mode;
        std::vector<Long64_t> n_data_mode;
    };

    // =========================================================
    // Helpers
    // =========================================================
    void styleData(TH1D* h) {
        if (!h) return;

        h->SetMarkerStyle(20);
        h->SetMarkerSize(0.9);
        h->SetMarkerColor(kBlack);
        h->SetLineColor(kBlack);
        h->SetLineWidth(2);
    }

    const std::vector<double>& getRecoBinEdges(SelectionType selectionType) {
        if (selectionType == SelectionType::NuMu) {
            return kRecoBinEdgesNuMu;
        }

        return kRecoBinEdgesNuE;
    }

    bool hasCustomRecoBinning(SelectionType selectionType) {
        const auto& edges = getRecoBinEdges(selectionType);
        return kUseCustomRecoBinning && edges.size() >= 2;
    }

    TH1D* makeRecoHist1D(const TString& name,
                         const TString& title,
                         SelectionType selectionType)
    {
        if (hasCustomRecoBinning(selectionType)) {
            const auto& edges = getRecoBinEdges(selectionType);

            return new TH1D(
                name,
                title,
                static_cast<int>(edges.size()) - 1,
                edges.data()
            );
        }

        return new TH1D(
            name,
            title,
            kRecoNBinsUniform,
            kRecoMinUniform,
            kRecoMaxUniform
        );
    }

    TH2D* makeRecoHist2D(const TString& name,
                         const TString& title,
                         SelectionType selectionType)
    {
        if (hasCustomRecoBinning(selectionType)) {
            const auto& edges = getRecoBinEdges(selectionType);

            return new TH2D(
                name,
                title,
                static_cast<int>(edges.size()) - 1,
                edges.data(),
                kCosNBins,
                kCosMin,
                kCosMax
            );
        }

        return new TH2D(
            name,
            title,
            kRecoNBinsUniform,
            kRecoMinUniform,
            kRecoMaxUniform,
            kCosNBins,
            kCosMin,
            kCosMax
        );
    }

    TH1D* makeTrueHist1D(const TString& name, const TString& title) {
        return new TH1D(
            name,
            title,
            kTrueNBins,
            kTrueMin,
            kTrueMax
        );
    }

    TH2D* makeTrueHist2D(const TString& name, const TString& title) {
        return new TH2D(
            name,
            title,
            kTrueNBins,
            kTrueMin,
            kTrueMax,
            kCosNBins,
            kCosMin,
            kCosMax
        );
    }

    void initCategory(Category& c,
                      const std::string& name,
                      const std::string& label,
                      SelectionType selectionType)
    {
        c.name  = name;
        c.label = label;
        c.selectionType = selectionType;

        c.recoE_mc.resize(NMODES, nullptr);
        c.recoCos_mc.resize(NMODES, nullptr);
        c.trueE_mc.resize(NMODES, nullptr);
        c.trueCos_mc.resize(NMODES, nullptr);

        c.rate_mc_mode.assign(NMODES, 0.0);
        c.rate_data_mode.assign(NMODES, 0.0);

        c.n_mc_mode.assign(NMODES, 0);
        c.n_data_mode.assign(NMODES, 0);

        for (int m = 0; m < NMODES; ++m) {
            c.recoE_mc[m] = makeRecoHist1D(
                Form("hRecoE_%s_mc_mode%d", name.c_str(), m),
                Form("Reco E_{#nu} by mode (%s);Reco E_{#nu} [GeV];Events", label.c_str()),
                c.selectionType
            );

            c.recoCos_mc[m] = new TH1D(
                Form("hRecoCos_%s_mc_mode%d", name.c_str(), m),
                Form("Reco cos#theta_{Z} by mode (%s);Reco cos#theta_{Z};Events", label.c_str()),
                kCosNBins,
                kCosMin,
                kCosMax
            );

            c.trueE_mc[m] = makeTrueHist1D(
                Form("hTrueE_%s_mc_mode%d", name.c_str(), m),
                Form("True E_{#nu} by mode (%s);True E_{#nu} [GeV];Events", label.c_str())
            );

            c.trueCos_mc[m] = new TH1D(
                Form("hTrueCos_%s_mc_mode%d", name.c_str(), m),
                Form("True cos#theta_{Z} by mode (%s);True cos#theta_{Z};Events", label.c_str()),
                kCosNBins,
                kCosMin,
                kCosMax
            );

            c.recoE_mc[m]->Sumw2();
            c.recoCos_mc[m]->Sumw2();
            c.trueE_mc[m]->Sumw2();
            c.trueCos_mc[m]->Sumw2();
        }

        c.recoE_data = makeRecoHist1D(
            Form("hRecoE_%s_data", name.c_str()),
            Form("Reco E_{#nu} (%s);Reco E_{#nu} [GeV];Events", label.c_str()),
            c.selectionType
        );

        c.recoCos_data = new TH1D(
            Form("hRecoCos_%s_data", name.c_str()),
            Form("Reco cos#theta_{Z} (%s);Reco cos#theta_{Z};Events", label.c_str()),
            kCosNBins,
            kCosMin,
            kCosMax
        );

        c.trueE_data = makeTrueHist1D(
            Form("hTrueE_%s_data", name.c_str()),
            Form("True E_{#nu} (%s);True E_{#nu} [GeV];Events", label.c_str())
        );

        c.trueCos_data = new TH1D(
            Form("hTrueCos_%s_data", name.c_str()),
            Form("True cos#theta_{Z} (%s);True cos#theta_{Z};Events", label.c_str()),
            kCosNBins,
            kCosMin,
            kCosMax
        );

        c.recoE_data->Sumw2();
        c.recoCos_data->Sumw2();
        c.trueE_data->Sumw2();
        c.trueCos_data->Sumw2();

        styleData(c.recoE_data);
        styleData(c.recoCos_data);
        styleData(c.trueE_data);
        styleData(c.trueCos_data);

        c.reco2D_mc = makeRecoHist2D(
            Form("hReco2D_%s_mc", name.c_str()),
            Form("Model reco: cos#theta_{Z} vs E_{#nu} (%s);Reco E_{#nu} [GeV];Reco cos#theta_{Z}", label.c_str()),
            c.selectionType
        );

        c.true2D_mc = makeTrueHist2D(
            Form("hTrue2D_%s_mc", name.c_str()),
            Form("Model true: cos#theta_{Z} vs E_{#nu} (%s);True E_{#nu} [GeV];True cos#theta_{Z}", label.c_str())
        );

        c.reco2D_data = makeRecoHist2D(
            Form("hReco2D_%s_data", name.c_str()),
            Form("Data reco: cos#theta_{Z} vs E_{#nu} (%s);Reco E_{#nu} [GeV];Reco cos#theta_{Z}", label.c_str()),
            c.selectionType
        );

        c.reco2D_mc->Sumw2();
        c.true2D_mc->Sumw2();
        c.reco2D_data->Sumw2();
    }

    THStack* makeStack(const std::vector<TH1D*>& hs,
                       const char* name,
                       TLegend*& leg)
    {
        THStack* st = new THStack(name, "");

        leg = new TLegend(0.64, 0.38, 0.90, 0.90);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->SetTextSize(0.03);

        for (int m = 0; m < NMODES; ++m) {
            if (!hs[m]) continue;
            if (hs[m]->Integral() <= 0.0) continue;

            hs[m]->SetFillColor(modeColor[m]);
            hs[m]->SetLineColor(kBlack);
            hs[m]->SetLineWidth(1);

            st->Add(hs[m]);
            leg->AddEntry(hs[m], modeTitle[m], "f");
        }

        return st;
    }

    void drawStackWithData(TCanvas& can,
                           const std::string& pdf,
                           const std::vector<TH1D*>& mcHists,
                           TH1D* dataHist,
                           const TString& title,
                           bool drawData = true)
    {
        TLegend* leg = nullptr;
        THStack* st = makeStack(mcHists, Form("st_%s", dataHist->GetName()), leg);

        can.Clear();
        can.SetLogz(false);

        st->SetTitle(title);
        st->Draw("hist");

        double maxY = st->GetMaximum();

        if (drawData && dataHist && dataHist->Integral() > 0.0) {
            if (dataHist->GetMaximum() > maxY) {
                maxY = dataHist->GetMaximum();
            }
        }

        if (maxY > 0.0) {
            st->SetMaximum(1.25 * maxY);
        }

        if (drawData && dataHist && dataHist->Integral() > 0.0) {
            dataHist->Draw("E1 same");
            leg->AddEntry(dataHist, "Data", "lep");
        }

        leg->Draw();
        can.Print(pdf.c_str());

        delete st;
        delete leg;
    }

    void draw2DPage(TCanvas& can,
                    const std::string& pdf,
                    TH2D* h2,
                    const TString& title,
                    bool logz = true)
    {
        if (!h2) return;

        can.Clear();
        can.SetRightMargin(0.14);
        can.SetLogz(logz);

        h2->SetTitle(title);
        h2->Draw("COLZ");

        can.Print(pdf.c_str());

        can.SetLogz(false);
        can.SetRightMargin(0.10);
    }

    void fillOneCategory(Category& cat,
                         bool isData,
                         int mode,
                         double w,
                         double Ereco,
                         double czreco,
                         double Etrue,
                         double cztrue)
    {
        // -----------------------------------------------------
        // Inclusive event-rate counters
        // -----------------------------------------------------
        if (!isData) {
            cat.rate_mc += w;
            cat.n_mc    += 1;
        } else {
            cat.rate_data += w;
            cat.n_data    += 1;
        }

        // -----------------------------------------------------
        // Per-mode event-rate counters
        //
        // For data/Asimov, this works only if mcMode exists.
        // If the data tree has no mcMode leaf, inclusive data rates
        // are still printed, but data per-mode rates remain zero.
        // -----------------------------------------------------
        if (mode >= 0 && mode < NMODES) {
            if (!isData) {
                cat.rate_mc_mode[mode] += w;
                cat.n_mc_mode[mode]    += 1;
            } else {
                cat.rate_data_mode[mode] += w;
                cat.n_data_mode[mode]    += 1;
            }
        }

        // -----------------------------------------------------
        // Fill histograms
        // -----------------------------------------------------
        if (!isData) {
            if (mode < 0 || mode >= NMODES) {
                return;
            }

            if (std::isfinite(Ereco)) {
                cat.recoE_mc[mode]->Fill(Ereco, w);
            }

            if (std::isfinite(czreco)) {
                cat.recoCos_mc[mode]->Fill(czreco, w);
            }

            if (std::isfinite(Etrue)) {
                cat.trueE_mc[mode]->Fill(Etrue, w);
            }

            if (std::isfinite(cztrue)) {
                cat.trueCos_mc[mode]->Fill(cztrue, w);
            }

            if (std::isfinite(Ereco) && std::isfinite(czreco)) {
                cat.reco2D_mc->Fill(Ereco, czreco, w);
            }

            if (std::isfinite(Etrue) && std::isfinite(cztrue)) {
                cat.true2D_mc->Fill(Etrue, cztrue, w);
            }
        }
        else {
            if (std::isfinite(Ereco)) {
                cat.recoE_data->Fill(Ereco, w);
            }

            if (std::isfinite(czreco)) {
                cat.recoCos_data->Fill(czreco, w);
            }

            // Keep only if preFit/data contains truth branches.
            if (std::isfinite(Etrue)) {
                cat.trueE_data->Fill(Etrue, w);
            }

            if (std::isfinite(cztrue)) {
                cat.trueCos_data->Fill(cztrue, w);
            }

            if (std::isfinite(Ereco) && std::isfinite(czreco)) {
                cat.reco2D_data->Fill(Ereco, czreco, w);
            }
        }
    }

    void fillFromBaseDir(TDirectory* baseDir,
                         Category& numuSel_CC,
                         Category& nueSel_CC,
                         Category& numuSel_NC,
                         Category& nueSel_NC,
                         Category& numuSel_All,
                         Category& nueSel_All,
                         bool isData)
    {
        if (!baseDir) return;

        TIter nextSubfolder(baseDir->GetListOfKeys());
        TKey* keySubfolder = nullptr;

        while ((keySubfolder = static_cast<TKey*>(nextSubfolder()))) {
            if (std::strcmp(keySubfolder->GetClassName(), "TDirectoryFile") != 0) {
                continue;
            }

            TDirectory* subDir = dynamic_cast<TDirectory*>(baseDir->Get(keySubfolder->GetName()));
            if (!subDir) continue;

            std::string subDirName = subDir->GetName();

            bool isNuE  = (subDirName.find("#nu_{e} Selec")   != std::string::npos);
            bool isNuMu = (subDirName.find("#nu_{#mu} Selec") != std::string::npos);

            if (!isNuE && !isNuMu) continue;

            TTree* tree = dynamic_cast<TTree*>(subDir->Get("events_TTree"));
            if (!tree) continue;

            TLeaf* leafMode        = tree->GetLeaf("mcMode");
            TLeaf* leafIsCC        = tree->GetLeaf("mcIsCC");
            TLeaf* leafGenWeight   = tree->GetLeaf("eventWeight");

            TLeaf* leafEnuRecoEle  = tree->GetLeaf("recoEnuEcalo");
            TLeaf* leafCosZRecoEle = tree->GetLeaf("recoEle_thetaZ");

            TLeaf* leafEnuRecoMu   = tree->GetLeaf("recoEnuLepCalo");
            TLeaf* leafCosZRecoMu  = tree->GetLeaf("recoMu_thetaZ");

            TLeaf* leafEnuTruth    = tree->GetLeaf("mcEnu");
            TLeaf* leafCosZTruth   = tree->GetLeaf("mcthetaZ");

            if (!leafIsCC || !leafGenWeight) {
                std::cerr << "Skipping tree in directory " << subDirName
                          << " because mcIsCC or eventWeight is missing.\n";
                continue;
            }

            const Long64_t nEntries = tree->GetEntries();
            //const Long64_t nEntries = 100;

            for (Long64_t i = 0; i < nEntries; ++i) {
                tree->GetEntry(i);

                int isCC = static_cast<int>(leafIsCC->GetValue(0));

                int mode = -1;
                if (leafMode) {
                    mode = static_cast<int>(leafMode->GetValue(0));
                }

                double w = static_cast<double>(leafGenWeight->GetValue(0));

                double Etrue = NAN;
                if (leafEnuTruth) {
                    Etrue = static_cast<double>(leafEnuTruth->GetValue(0));
                }

                double cztrue = NAN;
                if (leafCosZTruth) {
                    cztrue = static_cast<double>(leafCosZTruth->GetValue(0));
                }

                Category* catCCNC = nullptr;
                Category* catAll  = nullptr;

                double Ereco  = NAN;
                double czreco = NAN;

                // -----------------------------------------------------
                // Select the reconstructed variables and categories.
                //
                // catAll receives both CC and NC events.
                // catCCNC preserves the original CC/NC-separated output.
                // -----------------------------------------------------
                if (isNuMu) {
                    catAll = &numuSel_All;

                    if (leafEnuRecoMu) {
                        Ereco = static_cast<double>(leafEnuRecoMu->GetValue(0));
                    }

                    if (leafCosZRecoMu) {
                        czreco = static_cast<double>(leafCosZRecoMu->GetValue(0));
                    }

                    if (isCC == 1) {
                        catCCNC = &numuSel_CC;
                    }
                    else if (isCC == 0) {
                        catCCNC = &numuSel_NC;
                    }
                }
                else if (isNuE) {
                    catAll = &nueSel_All;

                    if (leafEnuRecoEle) {
                        Ereco = static_cast<double>(leafEnuRecoEle->GetValue(0));
                    }

                    if (leafCosZRecoEle) {
                        czreco = static_cast<double>(leafCosZRecoEle->GetValue(0));
                    }

                    if (isCC == 1) {
                        catCCNC = &nueSel_CC;
                    }
                    else if (isCC == 0) {
                        catCCNC = &nueSel_NC;
                    }
                }

                if (!catAll) continue;

                // Fill inclusive CC+NC category.
                fillOneCategory(
                    *catAll,
                    isData,
                    mode,
                    w,
                    Ereco,
                    czreco,
                    Etrue,
                    cztrue
                );

                // Fill original CC/NC-separated category.
                if (catCCNC) {
                    fillOneCategory(
                        *catCCNC,
                        isData,
                        mode,
                        w,
                        Ereco,
                        czreco,
                        Etrue,
                        cztrue
                    );
                }
            }
        }
    }

    void printCategoryPDF(const Category& c) {
        TCanvas can("can", "can", 1000, 700);

        std::string pdf = c.name + ".pdf";

        can.Print((pdf + "[").c_str());

        drawStackWithData(
            can,
            pdf,
            c.recoE_mc,
            c.recoE_data,
            Form("Reco E_{#nu} (%s);Reco E_{#nu} [GeV];Events", c.label.c_str()),
            true
        );

        drawStackWithData(
            can,
            pdf,
            c.recoCos_mc,
            c.recoCos_data,
            Form("Reco cos#theta_{Z} (%s);Reco cos#theta_{Z};Events", c.label.c_str()),
            true
        );

        bool hasTrueData =
            (c.trueE_data && c.trueE_data->Integral() > 0.0) ||
            (c.trueCos_data && c.trueCos_data->Integral() > 0.0);

        drawStackWithData(
            can,
            pdf,
            c.trueE_mc,
            c.trueE_data,
            Form("True E_{#nu} (%s);True E_{#nu} [GeV];Events", c.label.c_str()),
            hasTrueData
        );

        drawStackWithData(
            can,
            pdf,
            c.trueCos_mc,
            c.trueCos_data,
            Form("True cos#theta_{Z} (%s);True cos#theta_{Z};Events", c.label.c_str()),
            hasTrueData
        );

        draw2DPage(
            can,
            pdf,
            c.reco2D_mc,
            Form("Model reco: cos#theta_{Z} vs E_{#nu} (%s);Reco E_{#nu} [GeV];Reco cos#theta_{Z}", c.label.c_str()),
            true
        );

        draw2DPage(
            can,
            pdf,
            c.true2D_mc,
            Form("Model true: cos#theta_{Z} vs E_{#nu} (%s);True E_{#nu} [GeV];True cos#theta_{Z}", c.label.c_str()),
            true
        );

        draw2DPage(
            can,
            pdf,
            c.reco2D_data,
            Form("Data reco: cos#theta_{Z} vs E_{#nu} (%s);Reco E_{#nu} [GeV];Reco cos#theta_{Z}", c.label.c_str()),
            true
        );

        can.Print((pdf + "]").c_str());
    }

    void printEventRatesPerMode(const Category& c)
    {
        std::cout << "\n------------------------------------------------------------\n";
        std::cout << "Event rates per interaction mode for " << c.name
                  << " [" << c.label << "]\n";
        std::cout << "------------------------------------------------------------\n";

        std::cout << std::fixed << std::setprecision(6);

        std::cout << std::setw(14) << std::left << "Mode"
                  << " | " << std::setw(16) << "Model rate"
                  << " | " << std::setw(16) << "Data rate"
                  << " | " << std::setw(12) << "MC entries"
                  << " | " << std::setw(12) << "Data entries"
                  << "\n";

        std::cout << "------------------------------------------------------------\n";

        double sumModel = 0.0;
        double sumData  = 0.0;

        Long64_t sumModelEntries = 0;
        Long64_t sumDataEntries  = 0;

        for (int m = 0; m < NMODES; ++m) {
            const double rMC   = c.rate_mc_mode[m];
            const double rData = c.rate_data_mode[m];

            const Long64_t nMC   = c.n_mc_mode[m];
            const Long64_t nData = c.n_data_mode[m];

            if (rMC == 0.0 && rData == 0.0 && nMC == 0 && nData == 0) {
                continue;
            }

            sumModel += rMC;
            sumData  += rData;

            sumModelEntries += nMC;
            sumDataEntries  += nData;

            std::cout << std::setw(14) << std::left << modeTitle[m]
                      << " | " << std::setw(16) << rMC
                      << " | " << std::setw(16) << rData
                      << " | " << std::setw(12) << nMC
                      << " | " << std::setw(12) << nData
                      << "\n";
        }

        std::cout << "------------------------------------------------------------\n";
        std::cout << std::setw(14) << std::left << "TOTAL MODES"
                  << " | " << std::setw(16) << sumModel
                  << " | " << std::setw(16) << sumData
                  << " | " << std::setw(12) << sumModelEntries
                  << " | " << std::setw(12) << sumDataEntries
                  << "\n";

        std::cout << std::setw(14) << std::left << "INCLUSIVE"
                  << " | " << std::setw(16) << c.rate_mc
                  << " | " << std::setw(16) << c.rate_data
                  << " | " << std::setw(12) << c.n_mc
                  << " | " << std::setw(12) << c.n_data
                  << "\n";

        if (std::fabs(sumModel - c.rate_mc) > 1e-8 ||
            std::fabs(sumData  - c.rate_data) > 1e-8) {
            std::cout << "  Note: TOTAL MODES may differ from INCLUSIVE if some entries "
                      << "have missing or invalid mcMode.\n";
        }
    }

    void printEventRates(const Category& numuSel_CC,
                         const Category& nueSel_CC,
                         const Category& numuSel_NC,
                         const Category& nueSel_NC,
                         const Category& numuSel_All,
                         const Category& nueSel_All)
    {
        auto printOne = [](const Category& c) {
            std::cout << std::fixed << std::setprecision(6)
                      << "  " << std::setw(22) << std::left << c.name
                      << " | model rate = " << std::setw(16) << c.rate_mc
                      << " | data rate = "  << std::setw(16) << c.rate_data
                      << " | model entries = " << std::setw(10) << c.n_mc
                      << " | data entries = "  << c.n_data
                      << "\n";
        };

        const double numuSel_model = numuSel_CC.rate_mc + numuSel_NC.rate_mc;
        const double numuSel_data  = numuSel_CC.rate_data + numuSel_NC.rate_data;

        const double nueSel_model  = nueSel_CC.rate_mc + nueSel_NC.rate_mc;
        const double nueSel_data   = nueSel_CC.rate_data + nueSel_NC.rate_data;

        const double nc_model      = numuSel_NC.rate_mc + nueSel_NC.rate_mc;
        const double nc_data       = numuSel_NC.rate_data + nueSel_NC.rate_data;

        const double cc_model      = numuSel_CC.rate_mc + nueSel_CC.rate_mc;
        const double cc_data       = numuSel_CC.rate_data + nueSel_CC.rate_data;

        const double total_model   = cc_model + nc_model;
        const double total_data    = cc_data + nc_data;

        std::cout << "\n============================================================\n";
        std::cout << "Detailed category event rates\n";
        std::cout << "============================================================\n";

        printOne(numuSel_CC);
        printOne(nueSel_CC);
        printOne(numuSel_NC);
        printOne(nueSel_NC);
        printOne(numuSel_All);
        printOne(nueSel_All);

        std::cout << "\n============================================================\n";
        std::cout << "Summary event rates\n";
        std::cout << "============================================================\n";

        std::cout << std::fixed << std::setprecision(6);

        std::cout << "  numu selection total from CC+NC categories : model = " << numuSel_model
                  << " , data = " << numuSel_data << "\n";

        std::cout << "  numu selection total from inclusive category: model = " << numuSel_All.rate_mc
                  << " , data = " << numuSel_All.rate_data << "\n";

        std::cout << "  nue  selection total from CC+NC categories : model = " << nueSel_model
                  << " , data = " << nueSel_data << "\n";

        std::cout << "  nue  selection total from inclusive category: model = " << nueSel_All.rate_mc
                  << " , data = " << nueSel_All.rate_data << "\n";

        std::cout << "  NC total             : model = " << nc_model
                  << " , data = " << nc_data << "\n";

        std::cout << "  CC total             : model = " << cc_model
                  << " , data = " << cc_data << "\n";

        std::cout << "  ALL selected         : model = " << total_model
                  << " , data = " << total_data << "\n";

        std::cout << "\n============================================================\n";
        std::cout << "Per-mode event rates\n";
        std::cout << "============================================================\n";

        printEventRatesPerMode(numuSel_CC);
        printEventRatesPerMode(nueSel_CC);
        printEventRatesPerMode(numuSel_NC);
        printEventRatesPerMode(nueSel_NC);
        printEventRatesPerMode(numuSel_All);
        printEventRatesPerMode(nueSel_All);

        std::cout << "\n============================================================\n\n";
    }

    void printRecoBinning()
    {
        auto printEdges = [](const std::string& label,
                             SelectionType selectionType)
        {
            std::cout << "  " << label << ": ";

            if (hasCustomRecoBinning(selectionType)) {
                const auto& edges = getRecoBinEdges(selectionType);

                std::cout << "custom binning with edges:\n    ";

                for (size_t i = 0; i < edges.size(); ++i) {
                    std::cout << edges[i];

                    if (i + 1 < edges.size()) {
                        std::cout << ", ";
                    }
                }

                std::cout << "\n";
            } else {
                std::cout << "uniform binning: "
                          << kRecoNBinsUniform << " bins from "
                          << kRecoMinUniform << " to "
                          << kRecoMaxUniform << "\n";
            }
        };

        std::cout << "\nReco binning configuration:\n";

        printEdges("numu selection", SelectionType::NuMu);
        printEdges("nue  selection", SelectionType::NuE);

        std::cout << "\n";
    }

} // namespace

// =========================================================
// Main entry point
// =========================================================
void getERPlots()
{
    gStyle->SetOptStat(0);

    const char* inputFileName =
        "/gpfs/projects/McGrewGroup/uyevarou/Work/atm/gudi-anoa-sen25a/outputs/"
        "gundamFitter_config_DUNE_Asimov_DryRun.root";

    printRecoBinning();

    TFile* f = TFile::Open(inputFileName, "READ");

    if (!f || f->IsZombie()) {
        std::cerr << "Error opening file: " << inputFileName << std::endl;
        return;
    }

    // MC model
    TDirectory* modelDir =
        dynamic_cast<TDirectory*>(f->Get("FitterEngine/preFit/model"));

    if (!modelDir) {
        std::cerr << "Could not find FitterEngine/preFit/model" << std::endl;
        f->Close();
        return;
    }

    // Data / Asimov data
    TDirectory* dataDir =
        dynamic_cast<TDirectory*>(f->Get("FitterEngine/preFit/data"));

    if (!dataDir) {
        std::cerr << "Could not find FitterEngine/preFit/data" << std::endl;
        f->Close();
        return;
    }

    Category numuSel_CC;
    Category nueSel_CC;
    Category numuSel_NC;
    Category nueSel_NC;

    // New inclusive CC+NC categories.
    Category numuSel_All;
    Category nueSel_All;

    initCategory(
        numuSel_CC,
        "modeStacks_numuSel_CC",
        "#nu_{#mu} selection, CC",
        SelectionType::NuMu
    );

    initCategory(
        nueSel_CC,
        "modeStacks_nueSel_CC",
        "#nu_{e} selection, CC",
        SelectionType::NuE
    );

    initCategory(
        numuSel_NC,
        "modeStacks_numuSel_NC",
        "#nu_{#mu} selection, NC",
        SelectionType::NuMu
    );

    initCategory(
        nueSel_NC,
        "modeStacks_nueSel_NC",
        "#nu_{e} selection, NC",
        SelectionType::NuE
    );

    initCategory(
        numuSel_All,
        "modeStacks_numuSel_All",
        "#nu_{#mu} selection, CC+NC",
        SelectionType::NuMu
    );

    initCategory(
        nueSel_All,
        "modeStacks_nueSel_All",
        "#nu_{e} selection, CC+NC",
        SelectionType::NuE
    );

    fillFromBaseDir(
        modelDir,
        numuSel_CC,
        nueSel_CC,
        numuSel_NC,
        nueSel_NC,
        numuSel_All,
        nueSel_All,
        false
    );

    fillFromBaseDir(
        dataDir,
        numuSel_CC,
        nueSel_CC,
        numuSel_NC,
        nueSel_NC,
        numuSel_All,
        nueSel_All,
        true
    );

    printEventRates(
        numuSel_CC,
        nueSel_CC,
        numuSel_NC,
        nueSel_NC,
        numuSel_All,
        nueSel_All
    );

    // Original CC/NC-separated outputs.
    printCategoryPDF(numuSel_CC);
    printCategoryPDF(nueSel_CC);
    printCategoryPDF(numuSel_NC);
    printCategoryPDF(nueSel_NC);

    // New inclusive CC+NC outputs.
    printCategoryPDF(numuSel_All);
    printCategoryPDF(nueSel_All);

    f->Close();
}
