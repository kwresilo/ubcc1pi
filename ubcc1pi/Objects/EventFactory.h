/**
 *  @file  ubcc1pi/Objects/EventFactory.h
 *
 *  @brief The header file for the event factory class
 */

#ifndef UBCC1PI_OBJECTS_EVENT_FACTORY
#define UBCC1PI_OBJECTS_EVENT_FACTORY

#include "ubcc1pi/Objects/Event.h"
#include "art/Framework/Principal/Event.h"

namespace ubcc1pi
{

/**
 *  @brief  The event factory class is used to populate a ubcc1pi::Event from an art::Event 
 *          This implementation is separated from the Event class itself so that one can use the Event without any coupling to LArSoft
 */
class EventFactory
{
    public:

        /**
         *  @brief  Populate the output event with information from the input event
         *
         *  @param  event the input event
         *  @param  pOutputEvent the output event to populate
         */
        static void PopulateEvent(const art::Event &event, Event *pOutputEvent);

    private:

        /**
         *  @brief  Populate the metadata of the event
         *
         *  @param  event the input event
         *  @param  metadata the output metadata
         */
        static void PopulateEventMetadata(const art::Event &event, Event::Metadata &metadata);
};

} // namespace ubcc1pi

#endif
