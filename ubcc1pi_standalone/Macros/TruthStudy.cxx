/**
 *  @file  ubcc1pi_standalone/Macros/TruthStudy.cxx
 *
 *  @brief The implementation file of the TruthStudy macro
 */

#include "ubcc1pi_standalone/Macros/Macros.h"

#include "ubcc1pi_standalone/Objects/FileReader.h"

#include "ubcc1pi_standalone/Helpers/AnalysisHelper.h"
#include "ubcc1pi_standalone/Helpers/FormattingHelper.h"
#include "ubcc1pi_standalone/Helpers/PlottingHelper.h"

#include <string>
#include <vector>
#include <utility>
#include <map>

using namespace ubcc1pi;

namespace ubcc1pi_macros
{

void TruthStudy(const Config &config)
{
    // Open the file
    FileReader reader(config.files.overlaysFileName);
    auto pEvent = reader.GetBoundEventAddress();
    const auto nEvents = reader.GetNumberOfEvents();

    // Setup the counters
    std::vector<std::pair< std::string, float > > interactionToEventCountVect;
    float totalCC1PiEventCount = 0.f;
    float totalCCInclusiveEventCount = 0.f;

    // Setup the plots
    TH1F *hPionMomentum = new TH1F("hPionMomentum", "", 100, 0.f, 1.2f);
    TH1F *hGoldenPionMomentum = new TH1F("hGoldenPionMomentum", "", 100, 0.f, 1.2f);
    TH1F *hNonGoldenPionMomentum = new TH1F("hNonGoldenPionMomentum", "", 100, 0.f, 1.2f);

    // Loop over the events
    for (unsigned int i = 0; i < nEvents; ++i)
    {
        AnalysisHelper::PrintLoadingBar(i, nEvents);
        reader.LoadEvent(i);

        // Ge the event weight
        const auto weight = AnalysisHelper::GetNominalEventWeight(pEvent);

        if (AnalysisHelper::IsTrueCCInclusive(pEvent, config.global.useAbsPdg))
            totalCCInclusiveEventCount += weight;

        // Only consider true CC1pi events
        if (!AnalysisHelper::IsTrueCC1Pi(pEvent, config.global.useAbsPdg))
            continue;

        const auto interactionString = pEvent->truth.interactionString();
        const auto nuEnergy = pEvent->truth.nuEnergy();

        // Count up the total weight
        totalCC1PiEventCount += weight;

        // Find the vector entry for this interaction string (if it exists)
        auto iter = std::find_if(interactionToEventCountVect.begin(), interactionToEventCountVect.end(), [&](const auto &x) {
            return x.first == interactionString;
        });

        // If we haven't seen this interaction before, then add it to the vetor
        if (iter == interactionToEventCountVect.end())
        {
            interactionToEventCountVect.emplace_back(interactionString, weight);
        }
        else
        {
            // Otherwise increment the weight
            iter->second += weight;
        }

        const auto analysisData = AnalysisHelper::GetTruthAnalysisData(pEvent->truth, config.global.useAbsPdg, config.global.protonMomentumThreshold);
        hPionMomentum->Fill(analysisData.pionMomentum, weight);

        if (analysisData.hasGoldenPion)
        {
            hGoldenPionMomentum->Fill(analysisData.pionMomentum, weight);
        }
        else
        {
            hNonGoldenPionMomentum->Fill(analysisData.pionMomentum, weight);
        }
    }

    // Sort the vector by event count
    std::sort(interactionToEventCountVect.begin(), interactionToEventCountVect.end(), [](const auto &a, const auto &b) {
        return a.second > b.second;
    });

    // Print the result
    FormattingHelper::Table table({"Interaction", "", "Events", "Fraction"});
    for (unsigned int i = 0; i < interactionToEventCountVect.size(); ++i)
    {
        const auto entry = interactionToEventCountVect.at(i);
        const auto interactionString = entry.first;
        const auto nEvents = entry.second;

        table.AddEmptyRow();
        table.SetEntry("Interaction", interactionString);
        table.SetEntry("Events", nEvents);
        table.SetEntry("Fraction", nEvents / totalCC1PiEventCount);
    }

    table.WriteToFile("truthStudy_interactionBreakdown.md");

    auto pCanvas = PlottingHelper::GetCanvas();
    PlottingHelper::SetLineStyle(hPionMomentum, PlottingHelper::Default);
    PlottingHelper::SetLineStyle(hGoldenPionMomentum, PlottingHelper::GoldenPion);
    PlottingHelper::SetLineStyle(hNonGoldenPionMomentum, PlottingHelper::NonGoldenPion);

    // Normalise the histograms
    const auto scaleFactor = 1.f / hPionMomentum->GetEntries();
    hPionMomentum->Scale(scaleFactor);
    hGoldenPionMomentum->Scale(scaleFactor);
    hNonGoldenPionMomentum->Scale(scaleFactor);

    hPionMomentum->Draw("hist");
    hGoldenPionMomentum->Draw("hist same");
    hNonGoldenPionMomentum->Draw("hist same");
    PlottingHelper::SaveCanvas(pCanvas, "truthStudy_pionMomentum");

    std::cout << "Total CC inclusive events: " << totalCCInclusiveEventCount << std::endl;
    std::cout << "Total CC1Pi events: " << totalCC1PiEventCount << std::endl;
}

} // namespace ubcc1pi_macros
