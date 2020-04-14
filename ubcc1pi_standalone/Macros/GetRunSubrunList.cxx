#include <fstream>

#include "ubcc1pi_standalone/Objects/FileReader.h"

using namespace ubcc1pi;

int GetRunSubrunList(const std::string &inputFileName, const std::string &outputFileName = "runSubrunList.txt")
{
    ofstream outFile;
    outFile.open(outputFileName);

    FileReader reader(inputFileName);
    auto pSubrun = reader.GetBoundSubrunAddress();
 
    const auto nSubruns = reader.GetNumberOfSubruns();
    std::cout << "Writing " << nSubruns << " sub-runs to file: " << outputFileName << std::endl;
    
    for (unsigned int i = 0; i < nSubruns; ++i)
    {
        reader.LoadSubrun(i);

        if (!pSubrun->run.IsSet() || !pSubrun->subRun.IsSet())
        {
            std::cerr << "Error while reading subrun, metadata isn't set!" << std::endl;
            outFile.close();
            return 1;
        }

        outFile << pSubrun->run() << " " << pSubrun->subRun() << std::endl;
    }

    outFile.close();

    return 0;
}
