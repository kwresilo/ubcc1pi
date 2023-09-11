/**
 *  @file  ubcc1pi_standalone/Macros/GetRunSubrunList.cxx
 *
 *  @brief The implementation file of the GetRunSubrunList macro
 */

#include "ubcc1pi_standalone/Macros/Macros.h"
#include "ubcc1pi_standalone/Objects/FileReader.h"
#include "ubcc1pi_standalone/Helpers/AnalysisHelper.h"

#include <stdexcept>
#include <fstream>

using namespace ubcc1pi;

namespace ubcc1pi_macros
{

void GetRunSubrunList(const Config &config)
{
    // const bool zarkosToolExists = std::filesystem::exists("/uboone/app/users/zarko/getDataInfo.py"); // C++ 17 only
    ifstream f("/uboone/app/users/zarko/getDataInfo.py");
    const bool zarkosToolExists = f.good(); 

    // For each of the input files
    for (const auto &[sampleType, useThisFile, inputFilePath] : config.inputFiles)
    {
        if(sampleType!=AnalysisHelper::DataBNB && sampleType != AnalysisHelper::DataEXT) continue;
        // std::filesystem::path pathObj(inputFilePath);
        // std::string outputFilePath = "./" + pathObj.stem().string() + ".list"; // C++ 17 only
        std::string outputFilePath = "./";
        std::string fileName = inputFilePath.substr(inputFilePath.find_last_of("/\\") + 1);
        fileName = fileName.substr(0, fileName.find_last_of("."));
        outputFilePath += fileName + ".list";
        std::cout << outputFilePath << std::endl;

        // Open the output file
        ofstream outFile;
        outFile.open(outputFilePath);

        // Read the subruns
        FileReader<EventPeLEE, SubrunPeLEE> reader(inputFilePath, false);
        auto pSubrun = reader.GetBoundSubrunAddress();

        const auto nSubruns = reader.GetNumberOfSubruns();
        std::cout << "Writing " << nSubruns << " sub-runs to file: " << outputFilePath << std::endl;

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

        if(zarkosToolExists)
        {
            std::system(std::string("/uboone/app/users/zarko/getDataInfo.py -v2 --run-subrun-list ").append(outputFilePath).c_str());
        }
    }
}

} // namespace ubcc1pi_macros
