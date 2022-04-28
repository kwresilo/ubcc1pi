/**
 *  @file  ubcc1pi_standalone/Macros/CountPOT.cxx
 *
 *  @brief The implementation file of the CountPOT macro
 */

#include "ubcc1pi_standalone/Macros/Macros.h"

#include "ubcc1pi_standalone/Objects/FileReader.h"

#include "ubcc1pi_standalone/Helpers/AnalysisHelper.h"

using namespace ubcc1pi;

namespace ubcc1pi_macros
{

void CountPOT(const Config &config)
{
    // Get the input file names
    std::vector<std::string> fileNames;

    for (const auto run: config.global.runs)
    {
        switch (run)
        {
            case 1: 
            {
                if (config.countPOT.useOverlays) fileNames.push_back(config.filesRun1.overlaysFileName);
                if (config.countPOT.useDirt) fileNames.push_back(config.filesRun1.dirtFileName);
                if (config.countPOT.useNuWro) fileNames.push_back(config.filesRun1.nuWroFileName);
                if (config.countPOT.useDetectorVariations)
                {
                    for (const auto &[name, fileName] : config.filesRun1.detVarFiles)
                    {
                        fileNames.push_back(fileName);
                    }
                }
                break;
            }
            case 2:
            { 
                if (config.countPOT.useOverlays) fileNames.push_back(config.filesRun2.overlaysFileName);
                if (config.countPOT.useDirt) fileNames.push_back(config.filesRun2.dirtFileName);
                if (config.countPOT.useNuWro) fileNames.push_back(config.filesRun2.nuWroFileName);
                if (config.countPOT.useDetectorVariations)
                {
                    for (const auto &[name, fileName] : config.filesRun2.detVarFiles)
                    {
                        fileNames.push_back(fileName);
                    }
                }
                break;
            }
            case 3:
            { 
                if (config.countPOT.useOverlays) fileNames.push_back(config.filesRun3.overlaysFileName);
                if (config.countPOT.useDirt) fileNames.push_back(config.filesRun3.dirtFileName);
                if (config.countPOT.useNuWro) fileNames.push_back(config.filesRun3.nuWroFileName);
                if (config.countPOT.useDetectorVariations)
                {
                    for (const auto &[name, fileName] : config.filesRun3.detVarFiles)
                    {
                        fileNames.push_back(fileName);
                    }
                }
                break;
            }        
            default:
            {
                throw std::logic_error("TrainBDTs - Invalid run number");
            }
        }
    }

    // For each of the input files
    for (const auto &inputFileName : fileNames)
    {
        // Open the file
        FileReader reader(inputFileName);
        auto pSubrun = reader.GetBoundSubrunAddress();
        const auto nSubruns = reader.GetNumberOfSubruns();

        auto pEvent = reader.GetBoundEventAddress();
        const auto nEvents = reader.GetNumberOfEvents();

        // Count the total POT
        std::cout << "Processing file: " << inputFileName << std::endl;
        std::cout << "  - Getting total POT over " << nSubruns << " sub-runs." << std::endl;

        float totalPOT = 0.f;
        for (unsigned int i = 0; i < nSubruns; ++i)
        {
            reader.LoadSubrun(i);

            if (pSubrun->totalPOT.IsSet())
                totalPOT += pSubrun->totalPOT();
        }

        std::cout << "  - POT = " << totalPOT << std::endl;

        // Count the total number of events
        std::cout << "  - Getting event counts over " << nEvents << " events." << std::endl;
        float nEventsPassingCCInc = 0.f;
        float nEventsTrueCCInc = 0.f;
        float nEventsWeighted = 0.f;
        float nEventsPassingCCIncWeighted = 0.f;
        float nEventsTrueCCIncWeighted = 0.f;

        for (unsigned int i = 0; i < nEvents; ++i)
        {
            reader.LoadEvent(i);

            const auto weight = AnalysisHelper::GetNominalEventWeight(pEvent);
            const auto passesCCInc = pEvent->reco.passesCCInclusive.IsSet() && pEvent->reco.passesCCInclusive();
            const auto trueCCInc = AnalysisHelper::IsTrueCCInclusive(pEvent, true);

            nEventsPassingCCInc += (passesCCInc ? 1 : 0);
            nEventsTrueCCInc += (trueCCInc ? 1 : 0);
            nEventsWeighted += weight;
            nEventsPassingCCIncWeighted += (passesCCInc ? weight : 0);
            nEventsTrueCCIncWeighted += (trueCCInc ? weight : 0);
        }

        std::cout << "  - nEvents (unweighted) = " << nEvents << std::endl;
        std::cout << "  - nEvents passing CCInc. (unweighted) = " << nEventsPassingCCInc << std::endl;
        std::cout << "  - nEvents true CCInc. (unweighted) = " << nEventsTrueCCInc << std::endl;
        std::cout << "  - nEvents (nominal weight) = " << nEventsWeighted << std::endl;
        std::cout << "  - nEvents passing CCInc. (nominal weight) = " << nEventsPassingCCIncWeighted << std::endl;
        std::cout << "  - nEvents true CCInc. (nominal weight) = " << nEventsTrueCCIncWeighted << std::endl;
    }
}

} // namespace ubcc1pi_macros
