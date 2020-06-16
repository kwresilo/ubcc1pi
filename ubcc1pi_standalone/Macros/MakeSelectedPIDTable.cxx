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
    auto selection = SelectionHelper::GetDefaultSelection();
    const auto allCuts = selection.GetCuts();
    const auto lastCut = config.makeSelectedPIDTable.useGenericSelection ? config.global.lastCutGeneric : allCuts.back();

    std::cout << "Getting PID table after cut: " << lastCut << std::endl;

    if (std::find(allCuts.begin(), allCuts.end(), lastCut) == allCuts.end())
        throw std::invalid_argument("MakeSelectedPIDTable - chosen cut \"" + lastCut + "\" isn't known to the selection");

    //  Counter with index [recoPdgCode][truePdgCode][isSignalOnly]
    std::unordered_map< int, std::unordered_map< int, std::unordered_map<bool, float > > > recoToTruePdgMap;
                        
    float cc1piMuonCount = 0.f;
    float ccInclusiveMuonCount = 0.f;
    float ccInclusiveAnd1piMuonCount = 0.f;
    float cc1piMuonCorrectCount = 0.f;
    float ccInclusiveMuonCorrectCount = 0.f;
    float ccInclusiveAnd1piMuonCorrectCount = 0.f;

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
            const auto isTrueSignal = isOverlay && AnalysisHelper::IsTrueCC1Pi(pEvent, config.global.useAbsPdg) &&
                                      (config.makeSelectedPIDTable.goldenPionIsSignal ? (nGoldenPions != 0) : true);

            // Run the event selection and store which cuts are passed
            std::vector<std::string> cutsPassed;
            std::vector<int> assignedPdgCodes;
            selection.Execute(pEvent, cutsPassed, assignedPdgCodes);
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

                //// BEGIN DEBUG
                const auto isTrueMuon = (truePdgCode == 13);
                const auto isCC1PiMuonCandidate = (recoPdgCode == 13);
                const auto isCCInclusiveMuonCandidate = particle.isCCInclusiveMuonCandidate();

                if (isCC1PiMuonCandidate)
                {
                    cc1piMuonCount += weight;
                    cc1piMuonCorrectCount += isTrueMuon ? weight : 0.f;
                }

                if (isCCInclusiveMuonCandidate)
                {
                    ccInclusiveMuonCount += weight;
                    ccInclusiveMuonCorrectCount += isTrueMuon ? weight : 0.f;
                }

                if (isCC1PiMuonCandidate && isCCInclusiveMuonCandidate)
                {
                    ccInclusiveAnd1piMuonCount += weight;
                    ccInclusiveAnd1piMuonCorrectCount += isTrueMuon ? weight : 0.f;
                }
                //// END DEBUG

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
            }
        }
    }
    
    // Extract the PDG codes we want to print
    std::vector<int> recoPdgs, truePdgs;
    std::unordered_map<int, float> recoPdgToTotalWeightMap, recoPdgToTotalWeightMapSignal;

    for (const auto &recoEntry : recoToTruePdgMap)
    {
        const auto &recoPdgCode = recoEntry.first;

        // Store this PDG code if not yet seen
        if (std::find(recoPdgs.begin(), recoPdgs.end(), recoPdgCode) == recoPdgs.end())
        {
            recoPdgs.push_back(recoPdgCode);
            recoPdgToTotalWeightMap.emplace(recoPdgCode, 0.f);
            recoPdgToTotalWeightMapSignal.emplace(recoPdgCode, 0.f);
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

    std::sort(recoPdgs.begin(), recoPdgs.end());
    std::sort(truePdgs.begin(), truePdgs.end());

    // Now fill the output tables
    auto tableHeaders = std::vector<std::string>({"True PDG", ""});
    for (const auto &pdg : recoPdgs)
        tableHeaders.push_back(std::to_string(pdg));

    FormattingHelper::Table table(tableHeaders);
    FormattingHelper::Table signalTable(tableHeaders);

    // Add the sum row
    table.AddEmptyRow();
    signalTable.AddEmptyRow();

    table.SetEntry("True PDG", "all");
    signalTable.SetEntry("True PDG", "all");

    for (const auto &recoPdg : recoPdgs)
    {
        const auto count = recoPdgToTotalWeightMap.at(recoPdg);
        const auto signalCount = recoPdgToTotalWeightMapSignal.at(recoPdg);
            
        table.SetEntry(std::to_string(recoPdg), count);
        signalTable.SetEntry(std::to_string(recoPdg), signalCount);
    }

    for (const auto &truePdg : truePdgs)
    {
        table.AddEmptyRow();
        signalTable.AddEmptyRow();

        table.SetEntry("True PDG", truePdg);
        signalTable.SetEntry("True PDG", truePdg);

        for (const auto &recoPdg : recoPdgs)
        {
            const auto countTotal = recoPdgToTotalWeightMap.at(recoPdg);
            const auto signalCountTotal = recoPdgToTotalWeightMapSignal.at(recoPdg);
            

            float count = 0.f;
            float signalCount = 0.f;

            const auto &entry = recoToTruePdgMap.at(recoPdg);
            const auto iter = entry.find(truePdg);

            if (iter != entry.end())
            {
                count = iter->second.at(false);
                signalCount = iter->second.at(true);
            }

            table.SetEntry(std::to_string(recoPdg), count / countTotal);
            signalTable.SetEntry(std::to_string(recoPdg), signalCount / signalCountTotal);
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
    std::cout << "Muons" << std::endl;
    std::cout << "  - CC1Pi muons = " << cc1piMuonCount << std::endl;
    std::cout << "  - CC1Pi correct muons = " << cc1piMuonCorrectCount << std::endl;
    std::cout << "  - CCInc muons = " << ccInclusiveMuonCount << std::endl;
    std::cout << "  - CCInc correct muons = " << ccInclusiveMuonCorrectCount << std::endl;
    std::cout << "  - CC 1Pi & Inc muons = " << ccInclusiveAnd1piMuonCount << std::endl;
    std::cout << "  - CC 1Pi & Inc correct muons = " << ccInclusiveAnd1piMuonCorrectCount << std::endl;
}

} // namespace ubcc1pi_macros
