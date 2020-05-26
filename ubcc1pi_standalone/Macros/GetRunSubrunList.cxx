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

    if (config.getRunSubrunList.useDataEXT)
        fileNames.emplace_back(config.files.dataEXTFileName, "runSubrunList_dataEXT.txt");
    
    if (config.getRunSubrunList.useDataBNB)
        fileNames.emplace_back(config.files.dataBNBFileName, "runSubrunList_dataBNB.txt");
    
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
