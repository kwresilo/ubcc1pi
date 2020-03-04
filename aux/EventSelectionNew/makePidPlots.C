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
        {hPionG,  "hPionG",  "t_pdgCode == 211 && t_isGolden"},
        {hPionNG, "hPionNG", "t_pdgCode == 211 && !t_isGolden"},
        {hMuon,   "hMuon",   "t_pdgCode == 13"},
        {hProton, "hProton", "t_pdgCode == 2212"},
        {hExt,    "hExt",    "!t_hasMatchedMCParticle"}
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
        tree->Draw((variable + " >> " + histName).c_str(), ("(isCC1PiEvent && r_areFeaturesAvailable && r_isContained && t_truthMatchCompleteness > 0.5) && (" + (cut.empty() ? "0==0" : cut ) + ") && (" + selection + ")").c_str());
        //tree->Draw((variable + " >> " + histName).c_str(), ("(r_isEventSelected && r_isMaxGoldenPionBDTResponse && r_areFeaturesAvailable && r_isContained) && (" + (cut.empty() ? "0==0" : cut ) + ") && (" + selection + ")").c_str());

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
    /*
    MakePlot(particles, "r_logBragg_pToMIP", "", 80, -8.f, 8.f, "r_logBragg_pToMIP");
    MakePlot(particles, "r_logBragg_piToMIP", "", 80, -4.f, 8.f, "r_logBragg_piToMIP");
    MakePlot(particles, "r_trackShower", "", 80, 0.f, 1.f, "r_trackShower");
    MakePlot(particles, "r_rmsTrackDeviation", "", 80, 0.f, 0.5f, "r_rmsTrackDeviation");
    MakePlot(particles, "r_nSpacePointsInSphere5", "", 80, 0, 160, "r_nSpacePointsInSphere5");
    MakePlot(particles, "r_nDownstreamHits", "", 80, 0, 400, "r_nDownstreamHits");
    MakePlot(particles, "r_nDownstreamHits", "", 80, 1, 401, "r_nDownstreamHits_noZero");
    MakePlot(particles, "r_nDescendents", "", 6, 0, 6, "r_nDescendents");
    */
    
    MakePlot(particles, "r_forwardLikelihood", "", 80, 0.4, 0.65, "r_forwardLikelihood");
    MakePlot(particles, "r_dEdxTruncMeanStart", "", 80, 0, 10, "r_dEdxTruncMeanStart");
    //MakePlot(particles, "r_goldenPionBDTResponse", "", 80, -0.75f, 0.6f, "r_goldenPionBDTResponse_goldenPionCandidates");
}
