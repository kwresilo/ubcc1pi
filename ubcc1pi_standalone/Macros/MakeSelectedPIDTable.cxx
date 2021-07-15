/**
 *  @file  ubcc1pi_standalone/Macros/MakeSelectedPIDTable.cxx
 *
 *  @brief The implementation file of the MakeSelectedPIDTable macro
 */

#include "ubcc1pi_standalone/Macros/Macros.h"

#include "ubcc1pi_standalone/Objects/FileReader.h"

#include "ubcc1pi_standalone/Helpers/AnalysisHelper.h"
#include "ubcc1pi_standalone/Helpers/SelectionHelper.h"
#include "ubcc1pi_standalone/Helpers/FormattingHelper.h"
#include "ubcc1pi_standalone/Helpers/NormalisationHelper.h"

using namespace ubcc1pi;

namespace ubcc1pi_macros
{

void MakeSelectedPIDTable(const Config &config)
{
    //
    // Setup the input files
    //
    std::vector< std::tuple<AnalysisHelper::SampleType, std::string, float> > inputData;

    inputData.emplace_back(AnalysisHelper::Overlay, config.files.overlaysFileName, NormalisationHelper::GetOverlaysNormalisation(config));
    inputData.emplace_back(AnalysisHelper::Dirt,    config.files.dirtFileName,     NormalisationHelper::GetDirtNormalisation(config));
    inputData.emplace_back(AnalysisHelper::DataEXT, config.files.dataEXTFileName,  NormalisationHelper::GetDataEXTNormalisation(config));

    //
    // Get the selection
    //
    auto selection = SelectionHelper::GetSelection(config.global.selection);
    const auto allCuts = selection.GetCuts();
    const auto lastCut = config.global.selection=="Default"&&config.makeSelectedPIDTable.useGenericSelection ? config.global.lastCutGeneric : allCuts.back();

    std::cout << "Getting PID table after cut: " << lastCut << std::endl;

    if (std::find(allCuts.begin(), allCuts.end(), lastCut) == allCuts.end())
        throw std::invalid_argument("MakeSelectedPIDTable - chosen cut \"" + lastCut + "\" isn't known to the selection");

    //  Counter with index [recoPdgCode][truePdgCode][isSignalOnly]
    std::unordered_map< int, std::unordered_map< int, std::unordered_map<bool, float > > > recoToTruePdgMap;

    // Counter with index [recoPdgCode][truePdgCode] for dirt events
    std::unordered_map< int, std::unordered_map< int, float > > recoToTruePdgMapDirt;

    for (const auto [sampleType, fileName, normalisation] : inputData)
    {
        std::cout << "Reading input file: " << fileName << std::endl;

        FileReader reader(fileName);
        auto pEvent = reader.GetBoundEventAddress();

        const bool isOverlay = (sampleType == AnalysisHelper::Overlay);
        const bool isDirt = (sampleType == AnalysisHelper::Dirt);

        const auto nEvents = reader.GetNumberOfEvents();
        for (unsigned int i = 0; i < nEvents; ++i)
        {
            AnalysisHelper::PrintLoadingBar(i, nEvents);

            reader.LoadEvent(i);

            const auto weight = AnalysisHelper::GetNominalEventWeight(pEvent) * normalisation;

            // For speed skip events that didn't even pass the pre-selection
            if (!pEvent->reco.passesCCInclusive())
                continue;

            const auto truthParticles = pEvent->truth.particles; // ATTN this is empty for data events
            const auto recoParticles = pEvent->reco.particles;

            // Determine if this is a signal event
            const auto nGoldenPions = AnalysisHelper::CountGoldenParticlesWithPdgCode(AnalysisHelper::SelectVisibleParticles(truthParticles), 211, config.global.useAbsPdg);
            bool isTrueSignal = isOverlay && AnalysisHelper::IsTrueCC1Pi(pEvent, config.global.useAbsPdg) &&
                                (config.makeSelectedPIDTable.goldenPionIsSignal ? (nGoldenPions != 0) : true);

            //// BEGIN TEST
            if (isTrueSignal && config.makeSelectedPIDTable.onlyLowMomentumPions)
            {
                // Get the true pion momentum
                const auto truthData = AnalysisHelper::GetTruthAnalysisData(pEvent->truth, config.global.useAbsPdg, config.global.protonMomentumThreshold);

                // Check if this pion has sufficiently low momentum
                isTrueSignal = (truthData.pionMomentum < config.makeSelectedPIDTable.pionMomentumThreshold);
            }
            //// END TEST

            // Run the event selection and store which cuts are passed
            const auto &[isSelectedTotal, cutsPassed, assignedPdgCodes] = selection.Execute(pEvent);
            const auto isSelectedAtChosenCut = (std::find(cutsPassed.begin(), cutsPassed.end(), lastCut) != cutsPassed.end());

            // Don't bother with events that weren't selected
            if (!isSelectedAtChosenCut)
                continue;

            if (assignedPdgCodes.size() != recoParticles.size())
                throw std::logic_error("MakeSelectedPIDTable - The output particle PDG codes is the wrong size");

            for (unsigned int index = 0; index < recoParticles.size(); ++index)
            {
                const auto &particle = recoParticles.at(index);

                // Get the true PDG code (using 0 for EXTs)
                int truePdgCode = 0;
                if (isOverlay || isDirt)
                {
                    try
                    {
                        const auto truthParticle = AnalysisHelper::GetBestMatchedTruthParticle(particle, truthParticles, true);
                        truePdgCode = config.global.useAbsPdg ? std::abs(truthParticle.pdgCode()) : truthParticle.pdgCode();
                    }
                    catch (std::exception &)
                    {
                    }
                }

                const auto recoPdgCode = assignedPdgCodes.at(index);
                if (recoPdgCode != 13 && recoPdgCode != 211 && recoPdgCode != 2212)
                    throw std::logic_error("MakeSelectedPIDTable - Unknown assigned PDG code of " + std::to_string(recoPdgCode));

                // Get the mapping from true PDG to count for this reco pdg code, making it if it doesn't exist
                auto &truePdgToCountMap = recoToTruePdgMap[recoPdgCode];

                // Increment the count
                auto iter = truePdgToCountMap.find(truePdgCode);
                if (iter == truePdgToCountMap.end())
                {
                    std::unordered_map<bool, float> newEntry;
                    newEntry.emplace(true, isTrueSignal ? weight : 0.f);   // Only count "true" key for signal events
                    newEntry.emplace(false, weight);                       // Count "false" key for all events

                    truePdgToCountMap.emplace(truePdgCode, newEntry);
                }
                else
                {
                    if (isTrueSignal)
                        iter->second.at(true) += weight;

                    iter->second.at(false) += weight;
                }

                // Do the same for dirt events only
                if (sampleType == AnalysisHelper::Dirt)
                {
                    auto &truePdgToCountMapDirt = recoToTruePdgMapDirt[recoPdgCode];
                    auto iterDirt = truePdgToCountMapDirt.find(truePdgCode);
                    if (iterDirt == truePdgToCountMapDirt.end())
                    {
                        truePdgToCountMapDirt.emplace(truePdgCode, weight);
                    }
                    else
                    {
                        iterDirt->second += weight;
                    }
                }
            }
        }
    }

    // Extract the PDG codes we want to print
    std::vector<int> recoPdgs, truePdgs;
    std::unordered_map<int, float> recoPdgToTotalWeightMap, recoPdgToTotalWeightMapSignal, recoPdgToTotalWeightMapDirt;

    for (const auto &recoEntry : recoToTruePdgMap)
    {
        const auto &recoPdgCode = recoEntry.first;

        // Store this PDG code if not yet seen
        if (std::find(recoPdgs.begin(), recoPdgs.end(), recoPdgCode) == recoPdgs.end())
        {
            recoPdgs.push_back(recoPdgCode);
            recoPdgToTotalWeightMap.emplace(recoPdgCode, 0.f);
            recoPdgToTotalWeightMapSignal.emplace(recoPdgCode, 0.f);
            recoPdgToTotalWeightMapDirt.emplace(recoPdgCode, 0.f);
        }

        for (const auto &trueEntry : recoEntry.second)
        {
            const auto &truePdgCode = trueEntry.first;

            recoPdgToTotalWeightMap.at(recoPdgCode) += trueEntry.second.at(false);
            recoPdgToTotalWeightMapSignal.at(recoPdgCode) += trueEntry.second.at(true);

            // Store this PDG code if not yet seen
            if (std::find(truePdgs.begin(), truePdgs.end(), truePdgCode) == truePdgs.end())
                truePdgs.push_back(truePdgCode);
        }
    }

    for (const auto &recoEntry : recoToTruePdgMapDirt)
    {
        const auto &recoPdgCode = recoEntry.first;

        // Store this PDG code if not yet seen
        if (std::find(recoPdgs.begin(), recoPdgs.end(), recoPdgCode) == recoPdgs.end())
        {
            recoPdgs.push_back(recoPdgCode);
            recoPdgToTotalWeightMap.emplace(recoPdgCode, 0.f);
            recoPdgToTotalWeightMapSignal.emplace(recoPdgCode, 0.f);
            recoPdgToTotalWeightMapDirt.emplace(recoPdgCode, 0.f);
        }

        for (const auto &trueEntry : recoEntry.second)
        {
            const auto &truePdgCode = trueEntry.first;
            const auto weight = trueEntry.second;

            recoPdgToTotalWeightMapDirt.at(recoPdgCode) += weight;

            // Store this PDG code if not yet seen
            if (std::find(truePdgs.begin(), truePdgs.end(), truePdgCode) == truePdgs.end())
                truePdgs.push_back(truePdgCode);
        }
    }

    std::sort(recoPdgs.begin(), recoPdgs.end());
    std::sort(truePdgs.begin(), truePdgs.end());

    // Now fill the output tables
    auto tableHeaders = std::vector<std::string>({"True PDG", ""});
    for (const auto &pdg : recoPdgs)
        tableHeaders.push_back(std::to_string(pdg));

    FormattingHelper::Table table(tableHeaders);
    FormattingHelper::Table signalTable(tableHeaders);
    FormattingHelper::Table dirtTable(tableHeaders);

    // Add the sum row
    table.AddEmptyRow();
    signalTable.AddEmptyRow();
    dirtTable.AddEmptyRow();

    table.SetEntry("True PDG", "all");
    signalTable.SetEntry("True PDG", "all");
    dirtTable.SetEntry("True PDG", "all");

    for (const auto &recoPdg : recoPdgs)
    {
        const auto count = recoPdgToTotalWeightMap.at(recoPdg);
        const auto signalCount = recoPdgToTotalWeightMapSignal.at(recoPdg);
        const auto dirtCount = recoPdgToTotalWeightMapDirt.at(recoPdg);

        table.SetEntry(std::to_string(recoPdg), count);
        signalTable.SetEntry(std::to_string(recoPdg), signalCount);
        dirtTable.SetEntry(std::to_string(recoPdg), dirtCount);
    }

    for (const auto &truePdg : truePdgs)
    {
        table.AddEmptyRow();
        signalTable.AddEmptyRow();
        dirtTable.AddEmptyRow();

        table.SetEntry("True PDG", truePdg);
        signalTable.SetEntry("True PDG", truePdg);
        dirtTable.SetEntry("True PDG", truePdg);

        for (const auto &recoPdg : recoPdgs)
        {
            const auto countTotal = recoPdgToTotalWeightMap.at(recoPdg);
            const auto signalCountTotal = recoPdgToTotalWeightMapSignal.at(recoPdg);
            const auto dirtCountTotal = recoPdgToTotalWeightMapDirt.at(recoPdg);

            float count = 0.f;
            float signalCount = 0.f;
            float dirtCount = 0.f;

            const auto &entry = recoToTruePdgMap.at(recoPdg);
            const auto iter = entry.find(truePdg);

            if (iter != entry.end())
            {
                count = iter->second.at(false);
                signalCount = iter->second.at(true);
            }

            table.SetEntry(std::to_string(recoPdg), count / countTotal);
            signalTable.SetEntry(std::to_string(recoPdg), signalCount / signalCountTotal);

            const auto &entryDirt = recoToTruePdgMapDirt.at(recoPdg);
            const auto iterDirt = entryDirt.find(truePdg);

            if (iterDirt != entryDirt.end())
            {
                dirtCount = iterDirt->second;
            }

            dirtTable.SetEntry(std::to_string(recoPdg), dirtCount / dirtCountTotal);
        }
    }

    // Print the tables
    const auto prefix = std::string("pidTable") + (config.makeSelectedPIDTable.goldenPionIsSignal ? "_goldenPionIsSignal" : "") + std::string("_atCut-") + lastCut;

    FormattingHelper::PrintLine();
    std::cout << "All events" << std::endl;
    FormattingHelper::PrintLine();
    table.WriteToFile(prefix + "_allEvents.md");

    std::cout << std::endl;
    FormattingHelper::PrintLine();
    std::cout << "Signal events" << std::endl;
    FormattingHelper::PrintLine();
    signalTable.WriteToFile(prefix + "_signalEvents.md");

    std::cout << std::endl;
    FormattingHelper::PrintLine();
    std::cout << "Dirt events" << std::endl;
    FormattingHelper::PrintLine();
    dirtTable.WriteToFile(prefix + "_dirtEvents.md");
}

} // namespace ubcc1pi_macros
