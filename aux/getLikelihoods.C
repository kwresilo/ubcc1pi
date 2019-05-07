// A hacky script to get the likelihood histograms for each particle type

void chooseBinning(const std::vector<float> &values, const float defaultBinWidth, const unsigned int minEntriesPerBin, std::vector<float> &binEdges)
{
    std::vector<float> sortedValues = values;
    std::sort(sortedValues.begin(), sortedValues.end());

    if (!binEdges.empty())
        throw std::invalid_argument("Input vector to hold bin edges isn't empty");

    // Always start with the first at 0 with the default bin width
    binEdges.push_back(-std::numeric_limits<float>::epsilon());
    float width = defaultBinWidth;

    // The number of entries in the current bin - start at zero
    unsigned int nEntries = 0;

    for (const auto &value : sortedValues)
    {
        // Skip the incalcuable values
        if (value < 0.f)
            continue;

        // The lower edge of the current bin
        const auto edgeLower = binEdges.back();
        const auto edgeUpper = edgeLower + width;

        if (value < edgeLower)
            throw std::invalid_argument("Found a value lower than the current bin low edge");

        // Check if the value is in the current bin
        if (value < edgeUpper)
        {
            nEntries++;
        }
        else
        {
            // If we have enough entries in the current bin, then great! Start a new bin
            if (nEntries >= minEntriesPerBin)
            {
                binEdges.push_back(edgeUpper);
                width = defaultBinWidth;
                nEntries = 1;
            }
            // Otherwise, we need to grow the bin to accommodate this point
            else
            {
                width = value - edgeLower;
                nEntries++;
            }
        }
    }

    // Include the largest value as the last bin edge if required
    const auto largestValue = sortedValues.back();
    if (largestValue > binEdges.back())
        binEdges.push_back(largestValue + std::numeric_limits<float>::epsilon());
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void scaleByBinWidth(TH1F *h, const std::vector<float> &binEdges)
{
    if (h->GetNbinsX() != binEdges.size() - 1)
        throw std::invalid_argument("Bin edges don't match the input histogram");

    for (unsigned int i = 1; i <= h->GetNbinsX(); ++i)
    {
        const auto width = binEdges.at(i) - binEdges.at(i - 1);
        const auto content = h->GetBinContent(i);
        const auto prob = content / width;

        h->SetBinContent(i, prob);
    }

    h->Scale(1.f / h->Integral("width"));
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void makeProbabilityHist(const int chosenPdgCode)
{
    // ---------------------------------------------------
    // Settings
    // ---------------------------------------------------
    const float minPurity = 0.0f;
    const float minCompleteness = 0.5f;
    const unsigned int minEntriesPer1DBin = 1000;
    // ---------------------------------------------------

    // Set up the branch addresses
    float trackShowerScore, protonMIPScore, muonPionScore, trueMatchPurity, trueMatchCompleteness;
    int truePdgCode;
    bool hasMatchedMCParticle;

    particles->SetBranchAddress("trackShowerScore", &trackShowerScore);
    particles->SetBranchAddress("protonMIPScore", &protonMIPScore);
    particles->SetBranchAddress("muonPionScore", &muonPionScore);

    particles->SetBranchAddress("hasMatchedMCParticle", &hasMatchedMCParticle);
    particles->SetBranchAddress("trueMatchPurity", &trueMatchPurity);
    particles->SetBranchAddress("trueMatchCompleteness", &trueMatchCompleteness);
    particles->SetBranchAddress("truePdgCode", &truePdgCode);

    // Extract the scores for the particles that pass the cuts
    std::vector<float> trackShowerScores, protonMIPScores, muonPionScores;
    for (int i = 0; i < particles->GetEntries(); ++i)
    {
        particles->GetEntry(i);

        // Only consider particles that have the chosen PDG and are sufficiently well reconstructed
        if (!hasMatchedMCParticle || trueMatchPurity < minPurity || trueMatchCompleteness < minCompleteness || truePdgCode != chosenPdgCode)
            continue;

        trackShowerScores.push_back(trackShowerScore);
        protonMIPScores.push_back(protonMIPScore);
        muonPionScores.push_back(muonPionScore);
    }

    // Choose the binning
    std::vector<float> trackShowerBins, protonMIPBins, muonPionBins;
    chooseBinning(trackShowerScores, 0.01, minEntriesPer1DBin, trackShowerBins);
    chooseBinning(protonMIPScores, 5, minEntriesPer1DBin, protonMIPBins);
    chooseBinning(muonPionScores, 5, minEntriesPer1DBin, muonPionBins);

    // Fill the bins
    TH1F *hTrackShower = new TH1F(("trackShower_" + std::to_string(chosenPdgCode)).c_str(), "", trackShowerBins.size() - 1, &trackShowerBins[0]);
    TH1F *hProtonMIP = new TH1F(("protonMIP_" + std::to_string(chosenPdgCode)).c_str(), "", protonMIPBins.size() - 1, &protonMIPBins[0]);
    TH1F *hMuonPion = new TH1F(("muonPion_" + std::to_string(chosenPdgCode)).c_str(), "", muonPionBins.size() - 1, &muonPionBins[0]);

    for (unsigned int i = 0, nParticles = trackShowerScores.size(); i < nParticles; ++i)
    {
        hTrackShower->Fill(trackShowerScores.at(i));
        hProtonMIP->Fill(protonMIPScores.at(i));
        hMuonPion->Fill(muonPionScores.at(i));
    }

    // Divide by the bin widths to get a probability distribution
    scaleByBinWidth(hTrackShower, trackShowerBins);
    scaleByBinWidth(hProtonMIP, protonMIPBins);
    scaleByBinWidth(hMuonPion, muonPionBins);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void getLikelihoods()
{
    makeProbabilityHist(13);
    makeProbabilityHist(211);
    makeProbabilityHist(2212);

    // Retrieve the objects
    TH1F *hTrackShower_muon;
    TH1F *hTrackShower_pion;
    TH1F *hTrackShower_proton;
    TH1F *hProtonMIP_muon;
    TH1F *hProtonMIP_pion;
    TH1F *hProtonMIP_proton;
    TH1F *hMuonPion_muon;
    TH1F *hMuonPion_pion;
    TH1F *hMuonPion_proton;

    gDirectory->GetObject("trackShower_13", hTrackShower_muon);
    gDirectory->GetObject("trackShower_211", hTrackShower_pion);
    gDirectory->GetObject("trackShower_2212", hTrackShower_proton);
    gDirectory->GetObject("protonMIP_13", hProtonMIP_muon);
    gDirectory->GetObject("protonMIP_211", hProtonMIP_pion);
    gDirectory->GetObject("protonMIP_2212", hProtonMIP_proton);
    gDirectory->GetObject("muonPion_13", hMuonPion_muon);
    gDirectory->GetObject("muonPion_211", hMuonPion_pion);
    gDirectory->GetObject("muonPion_2212", hMuonPion_proton);

    TFile *f = new TFile("ubcc1pi_pid_hist.root", "recreate");
    hTrackShower_muon->Write("hTrackShower_muon");
    hTrackShower_pion->Write("hTrackShower_pion");
    hTrackShower_proton->Write("hTrackShower_proton");
    hProtonMIP_muon->Write("hProtonMIP_muon");
    hProtonMIP_pion->Write("hProtonMIP_pion");
    hProtonMIP_proton->Write("hProtonMIP_proton");
    hMuonPion_muon->Write("hMuonPion_muon");
    hMuonPion_pion->Write("hMuonPion_pion");
    hMuonPion_proton->Write("hMuonPion_proton");

    f->Write();
}
