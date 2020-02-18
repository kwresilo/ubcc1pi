/*
 *Br    0 :isSignal  : isSignal/O                                             *
*Entries :    57173 : Total  Size=      57824 bytes  File Size  =       9963 *
*Baskets :        2 : Basket Size=      32000 bytes  Compression=   5.75     *
*............................................................................*
 *Br    1 :isAvailable : isAvailable/O                                        *
*Entries :    57173 : Total  Size=      57842 bytes  File Size  =      12417 *
*Baskets :        2 : Basket Size=      32000 bytes  Compression=   4.62     *
*............................................................................*
 *Br    2 :shouldTrain : shouldTrain/O                                        *
*Entries :    57173 : Total  Size=      57842 bytes  File Size  =      11731 *
*Baskets :        2 : Basket Size=      32000 bytes  Compression=   4.89     *
*............................................................................*
 *Br    3 :trueIsGolden : trueIsGolden/O                                      *
*Entries :    57173 : Total  Size=      57848 bytes  File Size  =      12998 *
*Baskets :        2 : Basket Size=      32000 bytes  Compression=   4.41     *
*............................................................................*
 *Br    4 :truePdgCode : truePdgCode/I                                        *
*Entries :    57173 : Total  Size=     229865 bytes  File Size  =      37088 *
*Baskets :        8 : Basket Size=      32000 bytes  Compression=   6.18     *
*............................................................................*
 *Br    5 :truthMatchCompleteness : truthMatchCompleteness/F                  *
*Entries :    57173 : Total  Size=     229997 bytes  File Size  =     197251 *
*Baskets :        8 : Basket Size=      32000 bytes  Compression=   1.16     *
*............................................................................*
 *Br    6 :f_logBragg_pToMIP : f_logBragg_pToMIP/F                            *
*Entries :    57173 : Total  Size=     229937 bytes  File Size  =     201026 *
*Baskets :        8 : Basket Size=      32000 bytes  Compression=   1.14     *
*............................................................................*
 *Br    7 :f_logBragg_piToMIP : f_logBragg_piToMIP/F                          *
*Entries :    57173 : Total  Size=     229949 bytes  File Size  =     205136 *
*Baskets :        8 : Basket Size=      32000 bytes  Compression=   1.12     *
*............................................................................*
 *Br    8 :f_trackShower : f_trackShower/F                                    *
*Entries :    57173 : Total  Size=     229889 bytes  File Size  =     185962 *
*Baskets :        8 : Basket Size=      32000 bytes  Compression=   1.23     *
*............................................................................*
 *Br    9 :f_rmsTrackDeviation : f_rmsTrackDeviation/F                        *
*Entries :    57173 : Total  Size=     229961 bytes  File Size  =     184642 *
*Baskets :        8 : Basket Size=      32000 bytes  Compression=   1.24     *
*............................................................................*
 *Br   10 :f_rmsSequentialTrackDeviation : f_rmsSequentialTrackDeviation/F    *
*Entries :    57173 : Total  Size=     230081 bytes  File Size  =     187160 *
*Baskets :        8 : Basket Size=      32000 bytes  Compression=   1.23     *
*............................................................................*
 *Br   11 :f_offAxisSpacePointsRMSInSphere5 :                                 *
*         | f_offAxisSpacePointsRMSInSphere5/F                               *
*Entries :    57173 : Total  Size=     230117 bytes  File Size  =      99285 *
*Baskets :        8 : Basket Size=      32000 bytes  Compression=   2.31     *
*............................................................................*
 *Br   12 :f_offAxisSpacePointsRMSInSphere20 :                                *
*         | f_offAxisSpacePointsRMSInSphere20/F                              *
*Entries :    57173 : Total  Size=     230129 bytes  File Size  =     139722 *
*Baskets :        8 : Basket Size=      32000 bytes  Compression=   1.64     *
*............................................................................*
 *Br   13 :f_offAxisSpacePointsRMSInSphere30 :                                *
*         | f_offAxisSpacePointsRMSInSphere30/F                              *
*Entries :    57173 : Total  Size=     230129 bytes  File Size  =     149240 *
*Baskets :        8 : Basket Size=      32000 bytes  Compression=   1.54     *
*............................................................................*
 *Br   14 :f_nDownstreamHits : f_nDownstreamHits/I                            *
*Entries :    57173 : Total  Size=     229937 bytes  File Size  =      54579 *
*Baskets :        8 : Basket Size=      32000 bytes  Compression=   4.20     *
*............................................................................*
*/


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
        {hPionG,  "hPionG",  "truePdgCode == 211 && trueIsGolden"},
        {hPionNG, "hPionNG", "truePdgCode == 211 && !trueIsGolden"},
        {hMuon,   "hMuon",   "truePdgCode == 13"},
        {hProton, "hProton", "truePdgCode == 2212"},
        {hExt,    "hExt",    "truePdgCode < -999999"}
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
        tree->Draw((variable + " >> " + histName).c_str(), ("shouldTrain && (" + (cut.empty() ? "0==0" : cut ) + ") && (" + selection + ")").c_str());

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
    MakePlot(particles, "f_logBragg_pToMIP", "", 80, -8.f, 8.f, "f_logBragg_pToMIP");
    MakePlot(particles, "f_logBragg_piToMIP", "", 80, -4.f, 8.f, "f_logBragg_piToMIP");
    MakePlot(particles, "f_trackShower", "", 80, 0.f, 1.f, "f_trackShower");
    MakePlot(particles, "f_rmsTrackDeviation", "", 80, 0.f, 0.5f, "f_rmsTrackDeviation");
    MakePlot(particles, "f_rmsSequentialTrackDeviation", "", 80, 0.f, 0.05f, "f_rmsSequentialTrackDeviation");
    MakePlot(particles, "f_nSpacePointsInSphere5", "", 80, 0, 160, "f_nSpacePointsInSphere5");
    MakePlot(particles, "f_nOtherSpacePointsInSphere5", "", 80, 0, 80, "f_nOtherSpacePointsInSphere5");
    MakePlot(particles, "f_offAxisSpacePointsRMSInSphere5", "", 80, 0.f, 5.f, "f_offAxisSpacePointsRMSInSphere5");
    MakePlot(particles, "f_offAxisSpacePointsRMSInSphere20", "", 80, 0.f, 20.f, "f_offAxisSpacePointsRMSInSphere20");
    MakePlot(particles, "f_offAxisSpacePointsRMSInSphere30", "", 80, 0.f, 30.f, "f_offAxisSpacePointsRMSInSphere30");
    MakePlot(particles, "f_offAxisSpacePointsRMSInSphere5", "", 80, std::numeric_limits<float>::epsilon(), 5.f, "f_offAxisSpacePointsRMSInSphere5_noZero");
    MakePlot(particles, "f_offAxisSpacePointsRMSInSphere20", "", 80, std::numeric_limits<float>::epsilon(), 20.f, "f_offAxisSpacePointsRMSInSphere20_noZero");
    MakePlot(particles, "f_offAxisSpacePointsRMSInSphere30", "", 80, std::numeric_limits<float>::epsilon(), 30.f, "f_offAxisSpacePointsRMSInSphere30_noZero");
    MakePlot(particles, "f_nDownstreamHits", "", 80, 0, 400, "f_nDownstreamHits");
    MakePlot(particles, "f_nDownstreamHits", "", 80, 1, 401, "f_nDownstreamHits_noZero");
    MakePlot(particles, "f_nDescendents", "", 6, 0, 6, "f_nDescendents");
}
