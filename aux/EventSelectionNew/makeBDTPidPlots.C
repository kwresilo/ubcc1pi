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
        tree->Draw((variable + " >> " + histName).c_str(), ("(" + (cut.empty() ? "0==0" : cut ) + ") && (" + selection + ")").c_str());

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

int makeBDTPidPlots()
{
    const float minBDT = -0.75f;
    const float maxBDT = 0.6f;

    MakePlot(TestTree, "BDT", "", 60, minBDT, maxBDT, "BDTResponse");
        
    const std::string pionGString = "t_pdgCode == 211 && t_isGolden";
    const std::string pionNGString = "t_pdgCode == 211 && !t_isGolden";
    const std::string muonString = "t_pdgCode == 13";
    const std::string protonString = "t_pdgCode == 2212";

    const auto nPionG = TestTree->GetEntries(pionGString.c_str());
    const auto nPionNG = TestTree->GetEntries(pionNGString.c_str());
    const auto nMuon = TestTree->GetEntries(muonString.c_str());
    const auto nProton = TestTree->GetEntries(protonString.c_str());

    const size_t nSamples = 200;

    std::vector<float> fracPionGVect, fracPionNGVect, fracMuonVect, fracProtonVect, bdtCutVect;
    for (size_t i = 0; i < nSamples; ++i)
    {
        const auto bdtCut = minBDT + (maxBDT - minBDT) * (static_cast<float>(i) / static_cast<float>(nSamples - 1));
        const std::string cutStr = "BDT > " + std::to_string(bdtCut);

        const auto nPionGPassing = TestTree->GetEntries((cutStr + " && " + pionGString).c_str());
        const auto nPionNGPassing = TestTree->GetEntries((cutStr + " && " + pionNGString).c_str());
        const auto nMuonPassing = TestTree->GetEntries((cutStr + " && " + muonString).c_str());
        const auto nProtonPassing = TestTree->GetEntries((cutStr + " && " + protonString).c_str());
        
        const auto fracPionG = static_cast<float>(nPionGPassing) / static_cast<float>(nPionG);
        const auto fracPionNG = static_cast<float>(nPionNGPassing) / static_cast<float>(nPionNG);
        const auto fracMuon = static_cast<float>(nMuonPassing) / static_cast<float>(nMuon);
        const auto fracProton = static_cast<float>(nProtonPassing) / static_cast<float>(nProton);

        bdtCutVect.push_back(bdtCut);
        fracPionGVect.push_back(fracPionG);
        fracPionNGVect.push_back(fracPionNG);
        fracMuonVect.push_back(fracMuon);
        fracProtonVect.push_back(fracProton);
    }

    const auto bdtCutArray = bdtCutVect.data();
    const auto fracPionGArray = fracPionGVect.data();
    const auto fracPionNGArray = fracPionNGVect.data();
    const auto fracMuonArray = fracMuonVect.data();
    const auto fracProtonArray = fracProtonVect.data();

    TCanvas c;
    gStyle->SetOptStat(0);

    // Set up the histograms
    TGraph *gPionG = new TGraph(nSamples, bdtCutArray, fracPionGArray);
    TGraph *gPionNG = new TGraph(nSamples, bdtCutArray, fracPionNGArray);
    TGraph *gMuon = new TGraph(nSamples, bdtCutArray, fracMuonArray);
    TGraph *gProton = new TGraph(nSamples, bdtCutArray, fracProtonArray);

    gPionG->SetLineColor(kGreen + 1);
    gPionNG->SetLineColor(kGreen + 1);
    gMuon->SetLineColor(kAzure - 2);
    gProton->SetLineColor(kOrange - 3);

    gPionG->SetLineWidth(2);
    gPionNG->SetLineWidth(2);
    gMuon->SetLineWidth(2);
    gProton->SetLineWidth(2);

    gPionNG->SetLineStyle(2); // Dashed

    gPionG->Draw("AL");
    gPionNG->Draw("L");
    gMuon->Draw("L");
    gProton->Draw("L");
    
    // Save
    std::string name = "selectionEfficiencies";
    c.SaveAs((name + ".png").c_str());
    c.SaveAs((name + ".pdf").c_str());
    c.SaveAs((name + ".C").c_str());

    c.SetLogy();
    c.SaveAs((name + "_logy.png").c_str());
    c.SaveAs((name + "_logy.pdf").c_str());
    c.SaveAs((name + "_logy.C").c_str());
    c.SetLogy(0);

    // Make the ROC curves
    // Work out what is the signal
    std::vector<float> fracSignalVect;
    if (TestTree->GetEntries(("className == \"Signal\" && (" + pionGString + ")").c_str()) != 0)
    {
        fracSignalVect = fracPionGVect; 
    }
    else if (TestTree->GetEntries(("className == \"Signal\" && (" + pionNGString + ")").c_str()) != 0)
    {
        fracSignalVect = fracPionNGVect;
    }
    else if (TestTree->GetEntries(("className == \"Signal\" && (" + muonString + ")").c_str()) != 0)
    {
        fracSignalVect = fracMuonVect;
    }
    else if (TestTree->GetEntries(("className == \"Signal\" && (" + protonString + ")").c_str()) != 0)
    {
        fracSignalVect = fracProtonVect;
    }
    else
    {
        std::cout << "Can't work out which category is signal!" << std::endl;
        return 1;
    }

    const auto fracSignalArray = fracSignalVect.data();

    TGraph *gPionGRoc = new TGraph(nSamples, fracSignalArray, fracPionGArray);
    TGraph *gPionNGRoc = new TGraph(nSamples, fracSignalArray, fracPionNGArray);
    TGraph *gMuonRoc = new TGraph(nSamples, fracSignalArray, fracMuonArray);
    TGraph *gProtonRoc = new TGraph(nSamples, fracSignalArray, fracProtonArray);
    
    gPionGRoc->SetLineColor(kGreen + 1);
    gPionNGRoc->SetLineColor(kGreen + 1);
    gMuonRoc->SetLineColor(kAzure - 2);
    gProtonRoc->SetLineColor(kOrange - 3);

    gPionGRoc->SetLineWidth(2);
    gPionNGRoc->SetLineWidth(2);
    gMuonRoc->SetLineWidth(2);
    gProtonRoc->SetLineWidth(2);

    gPionNGRoc->SetLineStyle(2); // Dashed
    
    gPionNGRoc->Draw("AL");
    gMuonRoc->Draw("L");
    gProtonRoc->Draw("L");
    gPionGRoc->Draw("L");
    
    // Save
    name = "rocCurve";
    c.SaveAs((name + ".png").c_str());
    c.SaveAs((name + ".pdf").c_str());
    c.SaveAs((name + ".C").c_str());

    c.SetLogy();
    c.SaveAs((name + "_logy.png").c_str());
    c.SaveAs((name + "_logy.pdf").c_str());
    c.SaveAs((name + "_logy.C").c_str());
    c.SetLogy(0);
    
    // Calculate the integral of the ROC curve using trapezium integration
    // ATTN I think this method assumes a monatonically decreasing efficiency vs. bdt cut value for signal and background
    float rocSumPionG = 0.f;
    float rocSumPionNG = 0.f;
    float rocSumMuon = 0.f;
    float rocSumProton = 0.f;

    for (size_t i = 1; i < nSamples; ++i)
    {
        // Change in efficiencies in this interval
        const auto dFracSignal = fracSignalVect.at(i - 1) - fracSignalVect.at(i);
        const auto dFracPionG = fracPionGVect.at(i - 1) - fracPionGVect.at(i);
        const auto dFracPionNG = fracPionNGVect.at(i - 1) - fracPionNGVect.at(i);
        const auto dFracMuon = fracMuonVect.at(i - 1) - fracMuonVect.at(i);
        const auto dFracProton = fracProtonVect.at(i - 1) - fracProtonVect.at(i);

        rocSumPionG += (fracPionGVect.at(i) + 0.5f * dFracPionG) * dFracSignal;
        rocSumPionNG += (fracPionNGVect.at(i) + 0.5f * dFracPionNG) * dFracSignal;
        rocSumMuon += (fracMuonVect.at(i) + 0.5f * dFracMuon) * dFracSignal;
        rocSumProton += (fracProtonVect.at(i) + 0.5f * dFracProton) * dFracSignal;
    }

    std::cout << "FOM: Pion G  = " << rocSumPionG << std::endl;
    std::cout << "FOM: Pion NG = " << rocSumPionNG << std::endl;
    std::cout << "FOM: Muon    = " << rocSumMuon << std::endl;
    std::cout << "FOM: Proton  = " << rocSumProton << std::endl;

    delete gPionG;
    delete gPionNG;
    delete gMuon;
    delete gProton;
    
    delete gPionGRoc;
    delete gPionNGRoc;
    delete gMuonRoc;
    delete gProtonRoc;

    return 0;
}
