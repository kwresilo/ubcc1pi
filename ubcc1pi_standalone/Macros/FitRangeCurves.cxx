/**
 *  @file  ubcc1pi_standalone/Macros/FitRangeCurves.cxx
 *
 *  @brief The implementation file of the FitRangeCurves macro
 */

#include "ubcc1pi_standalone/Macros/Macros.h"

#include "ubcc1pi_standalone/Objects/FileReader.h"

#include "ubcc1pi_standalone/Helpers/AnalysisHelper.h"
#include "ubcc1pi_standalone/Helpers/PlottingHelper.h"
#include "ubcc1pi_standalone/Helpers/FormattingHelper.h"

#include <vector>
#include <map>

#include <TH2F.h>

using namespace ubcc1pi;

namespace ubcc1pi_macros
{

void FitRangeCurves(const Config &config)
{
    // Open the file
    FileReader reader(config.files.overlaysFileName);
    auto pEvent = reader.GetBoundEventAddress();
    const auto nEvents = reader.GetNumberOfEvents();

    // Setup the output data structure, index of true pdgCode
    std::map<int, std::vector<float> > rangeMap, momentumMap;
    std::map<int, std::shared_ptr<TH2F> > histMap;

    histMap.emplace(13, new TH2F("rangeCurve_muon", ";Range (cm); Momentum (GeV)", 60, 0.f, 300.f, 60, 0.f, 0.8f));
    histMap.emplace(211, new TH2F("rangeCurve_pion", ";Range (cm); Momentum (GeV)", 60, 0.f, 100.f, 60, 0.f, 0.4f));
    histMap.emplace(2212, new TH2F("rangeCurve_proton", ";Range (cm); Momentum (GeV)", 60, 0.f, 100.f, 60, 0.f, 1.2f));

    // Loop over the events
    for (unsigned int i = 0; i < nEvents; ++i)
    {
        AnalysisHelper::PrintLoadingBar(i, nEvents);
        reader.LoadEvent(i);

        // Only consider true CC1Pi events
        if (!AnalysisHelper::IsTrueCC1Pi(pEvent, config.global.useAbsPdg))
            continue;

        const auto weight = AnalysisHelper::GetNominalEventWeight(pEvent);
        const auto &truthParticles = pEvent->truth.particles;

        // Get the reco particles
        for (const auto &recoParticle : pEvent->reco.particles)
        {
            if (!AnalysisHelper::HasTrackFit(recoParticle))
                continue;

            if (!AnalysisHelper::IsContained(recoParticle))
                continue;

            bool hasGoodMatch = false;
            unsigned int truthParticleIndex = std::numeric_limits<unsigned int>::max();
            try
            {
                truthParticleIndex = AnalysisHelper::GetBestMatchedTruthParticleIndex(recoParticle, truthParticles);

                //const auto purity = recoParticle.truthMatchPurities().at(truthParticleIndex);
                const auto completeness = recoParticle.truthMatchCompletenesses().at(truthParticleIndex);
                hasGoodMatch = (completeness > 0.5f);
            }
            catch (const std::exception &)
            {
                hasGoodMatch = false;
            }

            if (!hasGoodMatch)
                continue;

            const auto &truthParticle = truthParticles.at(truthParticleIndex);
            const auto pdgCode = config.global.useAbsPdg ? std::abs(truthParticle.pdgCode()) : truthParticle.pdgCode();

            rangeMap[pdgCode].push_back(recoParticle.range());
            momentumMap[pdgCode].push_back(truthParticle.momentum());

            auto iter = histMap.find(pdgCode);
            if (iter != histMap.end())
            {
                iter->second->Fill(recoParticle.range(), truthParticle.momentum(), weight);
            }
        }
    }

    // Extract the used pdg codes
    std::vector<int> usedPdgCodes;
    for (const auto &entry : rangeMap)
    {
        if (std::find(usedPdgCodes.begin(), usedPdgCodes.end(), entry.first) == usedPdgCodes.end())
            usedPdgCodes.push_back(entry.first);
    }

    // Perform the fits
    auto pCanvas = PlottingHelper::GetCanvas();
    for (const auto &pdgCode : usedPdgCodes)
    {
        FormattingHelper::PrintLine();
        std::cout << "PdgCode = " << pdgCode << std::endl;
        FormattingHelper::PrintLine();

        const auto &ranges = rangeMap.at(pdgCode);
        const auto &momenta = momentumMap.at(pdgCode);

        AnalysisHelper::RangeToMomentumFitParameters params;
        if (pdgCode ==2212){
            AnalysisHelper::GetRangeToMomentumFitParameters(ranges, momenta, params);
        }
        else{
            AnalysisHelper::GetRangeToMomentumFitParameters(ranges, momenta, params);
        }

        FormattingHelper::Table table({"Parameter", "Value"});

        table.AddEmptyRow();
        table.SetEntry("Parameter", "a");
        table.SetEntry("Value", params.a);

        table.AddEmptyRow();
        table.SetEntry("Parameter", "b");
        table.SetEntry("Value", params.b);

        table.AddEmptyRow();
        table.SetEntry("Parameter", "c");
        table.SetEntry("Value", params.c);

        table.AddEmptyRow();
        table.SetEntry("Parameter", "d");
        table.SetEntry("Value", params.d);

        const std::string baseName = "fitRangeCurves_pdg" + std::to_string(pdgCode);
        table.WriteToFile(baseName + ".md");

        auto iter = histMap.find(pdgCode);
        if (iter != histMap.end())
        {
            iter->second->Draw("colz");

            auto pFunc = AnalysisHelper::GetRangeToMomentumFunction();
            AnalysisHelper::SetRangeToMomentumFunctionParameters(params, pFunc);

            pFunc->Draw("same");
            PlottingHelper::SaveCanvas(pCanvas, baseName + "_rangeVsMomentumFitted");
        }
    }
}

} // namespace ubcc1pi_macros
