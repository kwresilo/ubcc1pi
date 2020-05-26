#include "ubcc1pi_standalone/Objects/FileReader.h"

#include "ubcc1pi_standalone/Helpers/AnalysisHelper.h"

using namespace ubcc1pi;

int CountPOT(const std::string &inputFileName)
{
    FileReader reader(inputFileName);
    auto pSubrun = reader.GetBoundSubrunAddress();
    auto pEvent = reader.GetBoundEventAddress();
    
    const auto nSubruns = reader.GetNumberOfSubruns();
    std::cout << "Getting total POT over " << nSubruns << " sub-runs." << std::endl;
    
    float totalPOT = 0.f;
    for (unsigned int i = 0; i < nSubruns; ++i)
    {
        reader.LoadSubrun(i);

        if (pSubrun->totalPOT.IsSet())
            totalPOT += pSubrun->totalPOT();
    }

    std::cout << totalPOT << std::endl;
   

    const auto nEvents = reader.GetNumberOfEvents();
    std::cout << "Corresponds to " << nEvents << " events." << std::endl;
    float totalEventWeight = 0.f;
    float totalEventWeightPassingCCInclusive = 0.f;
    for (unsigned int i = 0; i < nEvents; ++i)
    {
        AnalysisHelper::PrintLoadingBar(i, nEvents);

        reader.LoadEvent(i);

        const auto weight = AnalysisHelper::GetNominalEventWeight(pEvent);
        totalEventWeight += weight;
        
        if (pEvent->reco.passesCCInclusive())
            totalEventWeightPassingCCInclusive += weight;
    }

    std::cout << "   - Nominal event weight = " << totalEventWeight << std::endl;
    std::cout << "   - ... passing CC inclusive = " << totalEventWeightPassingCCInclusive << std::endl;

    return 0;
}
