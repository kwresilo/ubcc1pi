#include "ubcc1pi_standalone/Objects/FileReader.h"

using namespace ubcc1pi;

int Test()
{
    //FileReader reader("/uboone/data/users/asmith/ubcc1pi/samples/smallTest/ubcc1piAnalysis_overlay.root");
    FileReader reader("/uboone/data/users/asmith/ubcc1pi/samples/smallTest/ubcc1piAnalysis_dataBNB.root");
    auto pEvent = reader.GetBoundEventAddress();
    auto pSubrun = reader.GetBoundSubrunAddress();

    // Print the events
    const auto nEvents = reader.GetNumberOfEvents();
    for (unsigned int i = 0; i < nEvents; ++i)
    {
        reader.LoadEvent(i);

        std::cout << "EVENT " << i << std::endl;
        pEvent->Print();
    }
    
    // Print the subruns
    const auto nSubruns = reader.GetNumberOfSubruns();
    for (unsigned int i = 0; i < nSubruns; ++i)
    {
        reader.LoadSubrun(i);

        std::cout << "SUBRUN " << i << std::endl;
        pSubrun->Print();
    }

    return 0;
}
