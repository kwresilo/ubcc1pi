/**
 *  @file  ubcc1pi_standalone/Macros/CCInclusiveMuonPIDStudy.cxx
 *
 *  @brief The implementation file of the CCInclusiveMuonPIDStudy macro
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

void CCInclusiveMuonPIDStudy(const Config &config)
{
    // Open the file
    FileReader reader(config.files.overlaysFileName);
    auto pEvent = reader.GetBoundEventAddress();
    const auto nEvents = reader.GetNumberOfEvents();
     
    // Setup the plots
    const std::string yLabel = "Fraction of events";
    PlottingHelper::MultiPlot momentumPlot("True muon momentum / GeV", yLabel, 40u, 0.f, 1.5f);
    PlottingHelper::MultiPlot cosThetaPlot("True muon cos(theta)", yLabel, 40u, -1.f, 1.f);

    // Setup the counters
    std::map<PlottingHelper::PlotStyle, float> trueParticleTypeToCountMap;
    float totalWeight = 0.f;

    // Loop over the events
    for (unsigned int i = 0; i < nEvents; ++i)
    {
        AnalysisHelper::PrintLoadingBar(i, nEvents);
        reader.LoadEvent(i);

        // Only consider true CC1pi events
        if (!AnalysisHelper::IsTrueCC1Pi(pEvent, config.global.useAbsPdg))
            continue;

        // Only consider events that pass the CC inclusive filter
        if (!pEvent->reco.passesCCInclusive())
            continue;

        const auto weight = AnalysisHelper::GetNominalEventWeight(pEvent);
        const auto analysisData = AnalysisHelper::GetTruthAnalysisData(pEvent->truth, config.global.useAbsPdg, config.global.protonMomentumThreshold);

        // Find the CC inclusive muon candidate
        bool foundCCInclusiveMuonCandidate = false;
        for (const auto &particle : pEvent->reco.particles)
        {
            if (!particle.isCCInclusiveMuonCandidate()) 
                continue;

            foundCCInclusiveMuonCandidate = true;
        
            // Get the plotting style based on the best matched MC particle (if it exists)
            const auto style = PlottingHelper::GetPlotStyle(particle, AnalysisHelper::Overlay, pEvent->truth.particles, false, config.global.useAbsPdg);

            auto iter = trueParticleTypeToCountMap.find(style);
            if (iter == trueParticleTypeToCountMap.end())
            {
                trueParticleTypeToCountMap.emplace(style, weight);
            }
            else
            {
                iter->second += weight;
            }

            totalWeight += weight;

            // Fill the plots
            momentumPlot.Fill(analysisData.muonMomentum, style, weight);
            cosThetaPlot.Fill(analysisData.muonCosTheta, style, weight);

            break;
        }

        if (!foundCCInclusiveMuonCandidate)
            throw std::logic_error("CCInclusiveMuonPIDStudy - no muon candidate found");
    }

    // Draw and save
    momentumPlot.SaveAsStacked("ccInclusiveMuonPID_momentum");
    cosThetaPlot.SaveAsStacked("ccInclusiveMuonPID_cosTheta");

    // Print the table
    FormattingHelper::Table table({"Particle type", "Events", "Fraction"});
    for (const auto &entry : trueParticleTypeToCountMap)
    {
        const auto style = entry.first;
        const auto nEvents = entry.second;

        std::string particleType;
        switch (style)
        {
            case PlottingHelper::Muon:
                particleType = "Muon";
                break;
            case PlottingHelper::Proton:
                particleType = "Proton";
                break;
            case PlottingHelper::GoldenPion:
                particleType = "Golden pion";
                break;
            case PlottingHelper::NonGoldenPion:
                particleType = "Non golden pion";
                break;
            case PlottingHelper::External:
                particleType = "External";
                break;
            default:
                particleType = "Other (" + std::to_string(style) + ")";
                break;
        }

        table.AddEmptyRow();
        table.SetEntry("Particle type", particleType);
        table.SetEntry("Events", nEvents);
        table.SetEntry("Fraction", nEvents / totalWeight);
    }

    table.WriteToFile("ccInclusiveMuonPID_breakdown.md");
}

} // namespace ubcc1pi_macros
