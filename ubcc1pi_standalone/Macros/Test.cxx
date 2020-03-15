#include "ubcc1pi_standalone/Objects/FileReader.h"

using namespace ubcc1pi;

int Test()
{
    FileReader reader("/uboone/app/users/asmith/cc1pi/dev/test/sandbox/ubcc1piAnalysis.root");
    auto pEvent = reader.GetBoundEventAddress();

    const auto nEvents = reader.GetNumberOfEvents();
    for (unsigned int i = 0; i < nEvents; ++i)
    {
        reader.LoadEvent(i);

        std::cout << "LOADED EVENT " << i << std::endl;
        pEvent->Print();
    }

    return 0;
}
