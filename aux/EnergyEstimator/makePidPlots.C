void MakePlot(TTree *tree, const std::string &variable, const std::string &cut, const unsigned int nBins, const float min, const float max, const std::string &name)
{
    TCanvas c;
    gStyle->SetOptStat(0);

    // Set up the histograms
    TH1F *hPionG = new TH1F("hPionG", "", nBins, min, max);
    TH1F *hPionNG = new TH1F("hPionNG", "", nBins, min, max);
    TH1F *hMuon = new TH1F("hMuon", "", nBins, min, max);
    TH1F *hProton = new TH1F("hProton", "", nBins, min, max);
    TH1F *hExt = new TH1F("hExt", "", nBins, min, max);

    hPionG->SetLineColor(kGreen + 1);
    hPionNG->SetLineColor(kGreen + 1);
    hMuon->SetLineColor(kAzure - 2);
    hProton->SetLineColor(kOrange - 3);
    hExt->SetLineColor(kGray + 3);

    hPionG->SetLineWidth(2);
    hPionNG->SetLineWidth(2);
    hMuon->SetLineWidth(2);
    hProton->SetLineWidth(2);
    hExt->SetLineWidth(2);

    hPionNG->SetLineStyle(2); // Dashed

    // Add them to a vector along with the relevent selection string
    std::vector< std::tuple<TH1F *, std::string, std::string> > histograms = {
        {hPionG,  "hPionG",  "truePdgCode == 211 && isGolden"},
        {hPionNG, "hPionNG", "truePdgCode == 211 && !isGolden"},
        {hMuon,   "hMuon",   "truePdgCode == 13"},
        {hProton, "hProton", "truePdgCode == 2212"},
        {hExt,    "hExt",    "!hasMatchedMCParticle"}
    };

    // Fill the histograms
    float minY = std::numeric_limits<float>::max();
    float maxY = -std::numeric_limits<float>::max();
    for (const auto &entry : histograms)
    {
        const auto &pHist = std::get<0>(entry);
        const auto &histName = std::get<1>(entry);
        const auto &selection = std::get<2>(entry);

        // ATTN here we apply the range cut instead of a matching completeness cut
        tree->Draw((variable + " >> " + histName).c_str(), ("isSignal && isPrimary && range > 6 && (" + (cut.empty() ? "0==0" : cut ) + ") && (" + selection + ")").c_str());

        const auto nEntries = pHist->GetEntries();

        if (nEntries > std::numeric_limits<float>::epsilon())
        {
            pHist->Scale(1.f / pHist->GetEntries());
            minY = std::min(minY, static_cast<float>(pHist->GetMinimum()));
        }

        maxY = std::max(maxY, static_cast<float>(pHist->GetMaximum()));
    }

    // Draw all in one plot
    bool isFirst = true;
    for (const auto &entry : histograms)
    {
        const auto &pHist = std::get<0>(entry);

        if (isFirst)
        {
            pHist->Draw("hist");
            pHist->GetYaxis()->SetRangeUser(std::max(minY / 1.1f, std::numeric_limits<float>::epsilon()), maxY * 1.1f);
            isFirst = false;
        }
        else
        {
            pHist->Draw("hist same");
        }
    }

    // Save
    c.SaveAs((name + ".png").c_str());
    c.SaveAs((name + ".pdf").c_str());
    c.SaveAs((name + ".C").c_str());

    c.SetLogy();
    c.SaveAs((name + "_logy.png").c_str());
    c.SaveAs((name + "_logy.pdf").c_str());
    c.SaveAs((name + "_logy.C").c_str());

    // Clean up
    for (const auto &entry : histograms)
    {
        auto &pHist = std::get<0>(entry);
        delete pHist;
    }

    histograms.clear();
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void makePidPlots()
{
    MakePlot(particles, "range", "", 80, 0.f, 600.f, "range");

    MakePlot(particles, "thetaYZ", "", 50, -3.2f, 3.2f, "thetaYZ");
    MakePlot(particles, "thetaXZ", "", 50, -3.2f, 3.2f, "thetaXZ");
    MakePlot(particles, "thetaXY", "", 50, -3.2f, 3.2f, "thetaXY");

    MakePlot(particles, "muonLikelihoodW", "", 80, 0.f, 1.f, "muonLikelihoodW");
    MakePlot(particles, "pionLikelihoodW", "", 80, 0.f, 1.f, "pionLikelihoodW");
    MakePlot(particles, "protonLikelihoodW", "", 80, 0.f, 1.f, "protonLikelihoodW");
    MakePlot(particles, "mipLikelihoodW", "", 80, 0.f, 1.f, "mipLikelihoodW");
    
    MakePlot(particles, "log(protonLikelihoodW/mipLikelihoodW)", "protonLikelihoodW > -1 && mipLikelihoodW > -1", 80, -8.f, 8.f, "logProtonToMipLikelihoodW");
    MakePlot(particles, "log(pionLikelihoodW/mipLikelihoodW)", "pionLikelihoodW > -1 && mipLikelihoodW > -1", 80, -6.f, 6.f, "logPionToMipLikelihoodW");
    MakePlot(particles, "log(pionLikelihoodW/muonLikelihoodW)", "pionLikelihoodW > -1 && muonLikelihoodW > -1", 80, -2.f, 2.f, "logPionToMuonLikelihoodW");
    
    MakePlot(particles, "nDaughters", "", 5, 0, 5, "nDaughters");

    MakePlot(particles, "daughterCosOpeningAngle", "nDaughters > 0", 80, -1.f, 1.f, "daughterCosOpeningAngle");
}
