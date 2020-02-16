// ATTN this is a quick and very dirty way so I don't have to re-write this code. In the future this should be made better
const std::string momentumString = "(0.000782884 * std::pow(range, 1.12535) + 0.0623516 * std::pow(range, 0.292286))";


// -----------------------------------------------------------------------------------------------------------------------------------------

int GetEntriesInMomentumRange(TTree *pTree, const float minMomentum, const float maxMomentum)
{
    return pTree->GetEntries(("isGolden && isPrimary && truePdgCode == 211 && range > 6 && " + momentumString + " > " + std::to_string(minMomentum) + " && " + momentumString + " <= " + std::to_string(maxMomentum)).c_str());
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void SaveAs(TCanvas *c, const std::string &title)
{
    c->SaveAs((title + ".png").c_str());
    c->SaveAs((title + ".pdf").c_str());
    c->SaveAs((title + ".root").c_str());
    c->SaveAs((title + ".C").c_str());
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void chooseBinning()
{
    TCanvas *c = new TCanvas();
    c->SetMargin(0.15, 0.15, 0.15, 0.15);
    gStyle->SetOptStat(false);

    // Set the range of bins we should produce
    //   - The first bin will start at minMomentum
    //   - The last bin will end be between maxMomentum and cutoffMomentum
    //   - If we start trying to fine tune bin widths in order to reach the target entries per bin, then cut off if we are chaning the bin
    //     width in steps smaller than minDeltaMomentum
    const auto minMomentum = 0.1111f;  // ATTN this corresponds to the range threshold applied of 6cm
    const auto maxMomentum = 0.6f;
    const auto cutoffMomentum = 1.f;
    const auto minDeltaMomentum = 0.0001f;

    // Set the range of target number of entries per bin to test
    //   - Check entries between minEntriesPerBinTarget and maxEntriesPerBinTarget
    //   - Check in steps of entriesPerBinTargetStep
    const auto minEntriesPerBinTarget = 1200;
    const auto maxEntriesPerBinTarget = 3600;
    const auto entriesPerBinTargetStep = 400;

    // Do the loop over target entries per bin
    for (auto nEntriesPerBinTarget = minEntriesPerBinTarget; nEntriesPerBinTarget <= maxEntriesPerBinTarget; nEntriesPerBinTarget += entriesPerBinTargetStep)
    {
        // The initial momentum step to use when choosing bins, this is changed dynamically to converge on the bin size quickly
        const auto maxMomentumStep = (maxMomentum - minMomentum) * (static_cast<float>(nEntriesPerBinTarget) / static_cast<float>(GetEntriesInMomentumRange(particles, minMomentum, maxMomentum)));

        // Vectors to store output bin details
        std::vector<float> binLowerVector, binUpperVector, binWidthVector, nEntriesVector;

        // Define the starting bin
        auto momentumStep = maxMomentumStep;
        auto binLower = minMomentum;
        auto binUpper = binLower + momentumStep;
        auto nEntriesLast = 0;

        std::cout << "=====================================================" << std::endl;
        std::cout << "  - nEntriesPerBinTarget: " << nEntriesPerBinTarget << std::endl;
        std::cout << "-----------------------------------------------------" << std::endl;

        // Keep chaning the bin size and adding new bins until we have spanned the desired range
        while (binUpperVector.empty() || binUpperVector.back() < maxMomentum)
        {
            const auto nEntries = GetEntriesInMomentumRange(particles, binLower, binUpper);

            // Save the details of the current bin if it has the correct number of entries, or we have passed a cutoff threshold
            //   - std::abs(momentumStep) < minDeltaMomentum --> avoid fine tuning bins
            //   - binUpper >= cutoffMomentum             --> avoid an infinite final bin
            if (nEntries == nEntriesPerBinTarget || std::abs(momentumStep) < minDeltaMomentum || binUpper >= cutoffMomentum)
            {
                // Save the details
                binLowerVector.push_back(binLower);
                binUpperVector.push_back(binUpper);
                binWidthVector.push_back(binUpper - binLower);
                nEntriesVector.push_back(nEntries);

                // Print to screen
                std::cout << "    - bin: " << binLowerVector.size();
                std::cout << ", width: " << binWidthVector.back();
                std::cout << ", range: " << binLowerVector.back() << " -> " << binUpperVector.back();
                std::cout << ", entries: " << nEntriesVector.back() << std::endl;

                // Reset for the next bin
                momentumStep = maxMomentumStep;
                binLower = binUpper;
                binUpper = binLower + momentumStep;
                nEntriesLast = 0;

                // Bail if the last bin is getting huge
                if (binUpper >= cutoffMomentum)
                    break;
            }
            else
            {
                // Modify the step size (can be negative) to converge on the desired bin width
                //          current step size   * always step toward the target                * if we have overshot the target then reduce the step size, otherwise just keep stepping
                momentumStep = std::abs(momentumStep) * ((nEntries > nEntriesPerBinTarget) ? -1 : 1) * (((nEntries - nEntriesPerBinTarget) * (nEntriesLast - nEntriesPerBinTarget) > 0) ? 1.f : 0.5f);
            }

            // Grow (or shrink) the current bin
            binUpper = std::min(binUpper + momentumStep, cutoffMomentum);
            nEntriesLast = nEntries;
        }

        // Print a summary of the bins
        const auto nBins = binLowerVector.size();
        std::cout << "-----------------------------------------------------" << std::endl;
        std::cout << "  - nBins:                " << nBins << std::endl;
        std::cout << "  - minBinWidth:          " << *std::min_element(binWidthVector.begin(), binWidthVector.end()) << std::endl;
        std::cout << "  - maxBinWidth:          " << *std::max_element(binWidthVector.begin(), binWidthVector.end()) << std::endl;
        std::cout << "  - minEntries:           " << *std::min_element(nEntriesVector.begin(), nEntriesVector.end()) << std::endl;
        std::cout << "  - maxEntries:           " << *std::max_element(nEntriesVector.begin(), nEntriesVector.end()) << std::endl;

        // Convert the bins into the format that root wants to see
        auto finalBinsVector = binLowerVector;
        finalBinsVector.push_back(binUpperVector.back());
        const auto finalBins = finalBinsVector.data();

        // ==============================================================================
        // Make the 1D histograms for reconstructed range
        TH1F *h1 = new TH1F("h1", "", nBins, finalBins);
        TH1F *h1Scaled = new TH1F("h1Scaled", "", nBins, finalBins);
        
        h1->GetXaxis()->SetTitle("Reco momentum / GeV");
        h1->GetYaxis()->SetTitle("Number of reco particles");
        particles->Draw((momentumString + " >> h1").c_str(), "isGolden && isPrimary && truePdgCode == 211 && range > 6");
        SaveAs(c, "goldenPionMomentum_target" + std::to_string(nEntriesPerBinTarget) + "EntriesPerBin");
        
        // Scale the histogram bins by their widths
        for(unsigned int i = 1; i <= nBins; ++i)
        {
            const auto nEntries = h1->GetBinContent(i);
            const auto width = h1->GetBinWidth(i);
            h1Scaled->SetBinContent(i, nEntries / width);
        }

        h1Scaled->GetXaxis()->SetTitle("Reco momentum / GeV");
        h1Scaled->GetYaxis()->SetTitle("Number of reco particles / (bin width / GeV)");
        h1Scaled->Draw();
        SaveAs(c, "goldenPionMomentum_target" + std::to_string(nEntriesPerBinTarget) + "EntriesPerBin_scaled");

        // Divide the unscaled histogram to make the target entries per bin at 1 on the y-axis
        h1->GetYaxis()->SetTitle(("Number of reco particles / " + std::to_string(nEntriesPerBinTarget)).c_str());
        h1->Scale(1.f / static_cast<float>(nEntriesPerBinTarget));
        h1->Draw("hist");
        SaveAs(c, "goldenPionMomentum_target" + std::to_string(nEntriesPerBinTarget) + "EntriesPerBin_norm");

        // ==============================================================================
        // Make the 2D histogram showing reco-to-true range
        TH2F *h2 = new TH2F("h2", "", nBins, finalBins, nBins, finalBins);
        TH2F *h2Scaled = new TH2F("h2Scaled", "", nBins, finalBins, nBins, finalBins);
        TH2F *h2Equal = new TH2F("h2Equal", "", nBins, 0, nBins, nBins, 0, nBins);
        TH2F *h2EqualScaled = new TH2F("h2EqualScaled", "", nBins, 0, nBins, nBins, 0, nBins);
        TH2F *h2Migration = new TH2F("h2Migration", "", nBins, 0, nBins, nBins, 0, nBins);
        
        h2->GetXaxis()->SetTitle("True momentum / GeV");
        h2->GetYaxis()->SetTitle("Reco momentum / GeV");
        h2->GetZaxis()->SetTitle("Number of reco particles");
        h2->GetZaxis()->SetTitleOffset(1.6);
        particles->Draw((momentumString + ":trueMomentum >> h2").c_str(), "isGolden && isPrimary && truePdgCode == 211 && range > 6", "colz");
        SaveAs(c, "goldenPionMomentum_recoToTrue_target" + std::to_string(nEntriesPerBinTarget) + "EntriesPerBin");
        
        // ==============================================================================
        // Scale the histogram bins by their widths
        for(unsigned int i = 1; i <= nBins; ++i)
        {
            for(unsigned int j = 1; j <= nBins; ++j)
            {
                const auto nEntries = h2->GetBinContent(i, j);
                const auto widthX = h2->GetXaxis()->GetBinWidth(i);
                const auto widthY = h2->GetYaxis()->GetBinWidth(j);
                const auto area = widthX * widthY;
                const auto scaled = nEntries / area;
                
                h2Equal->SetBinContent(i, j, nEntries);

                h2Scaled->SetBinContent(i, j, scaled);
                h2EqualScaled->SetBinContent(i, j, scaled);
            }
        }

        h2Scaled->GetXaxis()->SetTitle("True momentum / GeV");
        h2Scaled->GetYaxis()->SetTitle("Reco momentum / GeV");
        h2Scaled->GetZaxis()->SetTitle("Number of reco particles / (bin area / GeV^2)");
        h2Scaled->GetZaxis()->SetTitleOffset(1.6);
        h2Scaled->Draw("colz");
        SaveAs(c, "goldenPionMomentum_recoToTrue_target" + std::to_string(nEntriesPerBinTarget) + "EntriesPerBin_scaled");
        
        // ==============================================================================
        
        // Setup the latex objects for the dynamic bin labels
        TLatex latexX;
        latexX.SetTextAlign(kVAlignCenter + kHAlignRight);
        latexX.SetTextSize(0.02);
        latexX.SetTextAngle(90);
        
        TLatex latexY;
        latexY.SetTextAlign(kVAlignCenter + kHAlignRight);
        latexY.SetTextSize(0.02);
        
        TLatex latexZ;
        latexZ.SetTextAlign(kVAlignCenter + kHAlignCenter);
        latexZ.SetTextSize(std::min(0.04, 0.01 * (37.f / static_cast<float>(nBins)))); // ATTN magic numbers scale the bin text size with the bin

        // Now draw the version with equal bin widths
        // Remove existing bin labels
        h2Equal->GetXaxis()->SetLabelOffset(999);
        h2Equal->GetXaxis()->SetLabelSize(0);
        h2Equal->GetYaxis()->SetLabelOffset(999);
        h2Equal->GetYaxis()->SetLabelSize(0);
        
        h2Equal->GetXaxis()->SetTitleOffset(1.9);
        h2Equal->GetYaxis()->SetTitleOffset(1.3);
        h2Equal->GetZaxis()->SetTitleOffset(1.6);

        h2Equal->GetXaxis()->SetTitle("True momentum / GeV");
        h2Equal->GetYaxis()->SetTitle("Reco momentum / GeV");
        h2Equal->GetZaxis()->SetTitle("Number of reco particles");

        h2Equal->Draw("colz");
        
        // Draw the bins
        for (unsigned int i = 0; i < finalBinsVector.size(); ++i)
        {
            std::stringstream ss;
            ss << std::fixed << std::setprecision(2) << finalBinsVector.at(i);
            const auto label = (ss.str() + "  ").c_str();

            latexX.DrawLatex(i, 0, label);
            latexY.DrawLatex(0, i, label);
        }
        
        SaveAs(c, "goldenPionMomentum_recoToTrue_target" + std::to_string(nEntriesPerBinTarget) + "EntriesPerBin_equal");
        
        // ==============================================================================
        
        // Remove existing bin labels
        h2EqualScaled->GetXaxis()->SetLabelOffset(999);
        h2EqualScaled->GetXaxis()->SetLabelSize(0);
        h2EqualScaled->GetYaxis()->SetLabelOffset(999);
        h2EqualScaled->GetYaxis()->SetLabelSize(0);
        
        h2EqualScaled->GetXaxis()->SetTitleOffset(1.9);
        h2EqualScaled->GetYaxis()->SetTitleOffset(1.3);
        h2EqualScaled->GetZaxis()->SetTitleOffset(1.6);

        h2EqualScaled->GetXaxis()->SetTitle("True momentum / GeV");
        h2EqualScaled->GetYaxis()->SetTitle("Reco momentum / GeV");
        h2EqualScaled->GetZaxis()->SetTitle("Number of reco particles / (bin area / GeV^2)");

        h2EqualScaled->Draw("colz");
        
        // Draw the bins
        for (unsigned int i = 0; i < finalBinsVector.size(); ++i)
        {
            std::stringstream ss;
            ss << std::fixed << std::setprecision(2) << finalBinsVector.at(i);
            const auto label = (ss.str() + "  ").c_str();

            latexX.DrawLatex(i, 0, label);
            latexY.DrawLatex(0, i, label);
        }
        
        SaveAs(c, "goldenPionMomentum_recoToTrue_target" + std::to_string(nEntriesPerBinTarget) + "EntriesPerBin_equalScaled");

        
        // ==============================================================================
        // Fill the migration matrix
        //   Migration (or smearing) matrix elements, Sji = P(reco value in bin j | true value in bin i)
        for(unsigned int i = 1; i <= nBins; ++i)
        {
            // Get the total entries in this true bin i, by summing over the reco bins
            auto totalTrue = 0.f;
            for (unsigned int j = 1; j <= nBins; ++j)
                totalTrue += h2Equal->GetBinContent(i, j);

            for(unsigned int j = 1; j <= nBins; ++j)
            {
                // Scale down by the total true entries to get a probability
                const auto migrationValue = h2Equal->GetBinContent(i, j) / totalTrue;
                h2Migration->SetBinContent(i, j, migrationValue);
            }
        }
        
        // Remove existing bin labels
        h2Migration->GetXaxis()->SetLabelOffset(999);
        h2Migration->GetXaxis()->SetLabelSize(0);
        h2Migration->GetYaxis()->SetLabelOffset(999);
        h2Migration->GetYaxis()->SetLabelSize(0);
        
        h2Migration->GetXaxis()->SetTitleOffset(1.9);
        h2Migration->GetYaxis()->SetTitleOffset(1.3);
        h2Migration->GetZaxis()->SetTitleOffset(1.6);

        h2Migration->GetXaxis()->SetTitle("True momentum / GeV");
        h2Migration->GetYaxis()->SetTitle("Reco momentum / GeV");
        h2Migration->GetZaxis()->SetTitle("P(reco | true)");

        h2Migration->Draw("colz");
        
        // Add the bin labels
        for (unsigned int i = 0; i < finalBinsVector.size(); ++i)
        {
            std::stringstream ss;
            ss << std::fixed << std::setprecision(2) << finalBinsVector.at(i);
            const auto label = (ss.str() + "  ").c_str();

            latexX.DrawLatex(i, 0, label);
            latexY.DrawLatex(0, i, label);
        }

        // Add the numerical values
        auto diagonalSum = 0.f;
        auto offDiagonalSum = 0.f;

        for(unsigned int i = 1; i <= nBins; ++i)
        {
            for(unsigned int j = 1; j <= nBins; ++j)
            {
                const auto value = h2Migration->GetBinContent(i, j);

                if (i == j)
                {
                    diagonalSum += value;
                }
                else
                {
                    offDiagonalSum += value;
                }

                std::stringstream ss;  
                ss << std::fixed << std::setprecision(2) << value;
                const auto valueStr = ss.str().c_str();

                latexZ.DrawLatex(i - 0.5f, j - 0.5, valueStr);
            }
        }

        SaveAs(c, "goldenPionMomentum_recoToTrue_target" + std::to_string(nEntriesPerBinTarget) + "EntriesPerBin_migration");
        
        std::cout << "  - maxEntries:           " << *std::max_element(nEntriesVector.begin(), nEntriesVector.end()) << std::endl;
        std::cout << "  - diagonal:             " << diagonalSum << std::endl;
        std::cout << "  - off-diagonal:         " << offDiagonalSum << std::endl;

        if (diagonalSum + offDiagonalSum > std::numeric_limits<float>::epsilon())
            std::cout << "  - diagonality:          " << diagonalSum / (diagonalSum + offDiagonalSum) << std::endl;

        std::cout << "-----------------------------------------------------" << std::endl;

        delete h1;
        delete h1Scaled;
        delete h2;
        delete h2Scaled;
        delete h2Equal;
        delete h2EqualScaled;
        delete h2Migration;
    }

    delete c;
}
