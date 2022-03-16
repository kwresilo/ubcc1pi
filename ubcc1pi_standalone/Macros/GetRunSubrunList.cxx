/**
 *  @file  ubcc1pi_standalone/Macros/GetRunSubrunList.cxx
 *
 *  @brief The implementation file of the GetRunSubrunList macro
 */

#include "ubcc1pi_standalone/Macros/Macros.h"

#include "ubcc1pi_standalone/Objects/FileReader.h"

#include <stdexcept>
#include <fstream>

using namespace ubcc1pi;

namespace ubcc1pi_macros
{

void GetRunSubrunList(const Config &config)
{
    // Get the input and output file names
    std::vector< std::pair<std::string, std::string> > fileNames;

    for (const auto run: config.global.runs)
    {
        switch (run)
        {
            case 1:
            {
                if (config.getRunSubrunList.useDataEXT) fileNames.emplace_back(config.filesRun1.dataEXTFileName, "run1_runSubrunList_dataEXT.txt");
                if (config.getRunSubrunList.useDataBNB) fileNames.emplace_back(config.filesRun1.dataBNBFileName, "run1_runSubrunList_dataBNB.txt");
                if (config.getRunSubrunList.useNuWro) fileNames.emplace_back(config.filesRun1.nuWroFileName, "run1_runSubrunList_NuWro.txt");
                break;
            }
            case 2: 
            {
                if (config.getRunSubrunList.useDataEXT) fileNames.emplace_back(config.filesRun2.dataEXTFileName, "run2_runSubrunList_dataEXT.txt");
                if (config.getRunSubrunList.useDataBNB) fileNames.emplace_back(config.filesRun2.dataBNBFileName, "run2_runSubrunList_dataBNB.txt");
                if (config.getRunSubrunList.useNuWro) fileNames.emplace_back(config.filesRun2.nuWroFileName, "run2_runSubrunList_NuWro.txt");
                break;
            }
            case 3: 
            {
                if (config.getRunSubrunList.useDataEXT) fileNames.emplace_back(config.filesRun3.dataEXTFileName, "run3_runSubrunList_dataEXT.txt");
                if (config.getRunSubrunList.useDataBNB) fileNames.emplace_back(config.filesRun3.dataBNBFileName, "run3_runSubrunList_dataBNB.txt");
                if (config.getRunSubrunList.useNuWro) fileNames.emplace_back(config.filesRun3.nuWroFileName, "run3_runSubrunList_NuWro.txt");
                break;
            }        
            default:
                throw std::logic_error("GetRunSubrunList - Invalid run number");
        }
    }

    // For each of the input files
    for (const auto &ioFileName : fileNames)
    {
        const auto inputFileName = ioFileName.first;
        const auto outputFileName = ioFileName.second;

        // Open the output file
        ofstream outFile;
        outFile.open(outputFileName);

        // Read the subruns
        FileReader reader(inputFileName);
        auto pSubrun = reader.GetBoundSubrunAddress();

        const auto nSubruns = reader.GetNumberOfSubruns();
        std::cout << "Writing " << nSubruns << " sub-runs to file: " << outputFileName << std::endl;

        for (unsigned int i = 0; i < nSubruns; ++i)
        {
            reader.LoadSubrun(i);

            if (!pSubrun->run.IsSet() || !pSubrun->subRun.IsSet())
            {
                outFile.close();
                std::logic_error("Error while reading subrun, metadata isn't set!");
            }

            outFile << pSubrun->run() << " " << pSubrun->subRun() << std::endl;
        }

        outFile.close();
    }
}

} // namespace ubcc1pi_macros
