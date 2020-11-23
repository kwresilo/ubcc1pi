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

    if (config.countPOT.useOverlays)
        fileNames.push_back(config.files.overlaysFileName);
    
    if (config.countPOT.useDirt)
        fileNames.push_back(config.files.dirtFileName);

    if (config.countPOT.useDetectorVariations)
    {
        for (const auto &[run, name, fileName] : config.files.detVarFiles)
        {
            fileNames.push_back(fileName);
        }
    }

    // For each of the input files
    for (const auto &inputFileName : fileNames)
    {
        // Open the file
        FileReader reader(inputFileName);
        auto pSubrun = reader.GetBoundSubrunAddress();
        const auto nSubruns = reader.GetNumberOfSubruns();
        
        // Count the total POT
        std::cout << "Processing file: " << inputFileName << std::endl;
        std::cout << "  - Getting total POT over " << nSubruns << " sub-runs." << std::endl;
        
        float totalPOT = 0.f;
        for (unsigned int i = 0; i < nSubruns; ++i)
        {
            reader.LoadSubrun(i);
            totalPOT += pSubrun->totalPOT();
        }

        std::cout << "  - POT = " << totalPOT << std::endl;
    }
}

} // namespace ubcc1pi_macros
