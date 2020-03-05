/**
 *  @file  ubcc1pi/Objects/EventFactory.cxx
 *
 *  @brief The implementation of the event factory class
 */

#include "ubcc1pi/Objects/EventFactory.h"

namespace ubcc1pi
{

void EventFactory::PopulateEvent(const art::Event &event, Event *pOutputEvent)
{
    pOutputEvent->Reset();

    // Populate the metadata
    auto &metadata = pOutputEvent->metadata;
    EventFactory::PopulateEventMetadata(event, metadata);

    // Populate the truth info
    //auto &truth = pOutputEvent->truth;
    //EventFactory::PopulateEventTruthInfo(event, truth);

    // Populate the reco info
    //auto &reco = pOutputEvent->reco;
    //EventFactory::PopulateEventRecoInfo(event, reco);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void EventFactory::PopulateEventMetadata(const art::Event &event, Event::Metadata &metadata)
{
    metadata.run.Set(event.run());
    metadata.subRun.Set(event.subRun());
    metadata.event.Set(event.event());
}

} // namespace ubcc1pi
