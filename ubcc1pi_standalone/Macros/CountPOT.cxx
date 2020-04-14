#include "ubcc1pi_standalone/Objects/FileReader.h"

using namespace ubcc1pi;

int CountPOT(const std::string &inputFileName)
{
    FileReader reader(inputFileName);
    auto pSubrun = reader.GetBoundSubrunAddress();
    
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

    return 0;
}
