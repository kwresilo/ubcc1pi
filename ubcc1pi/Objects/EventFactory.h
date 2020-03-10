/**
 *  @file  ubcc1pi/Objects/EventFactory.h
 *
 *  @brief The header file for the event factory class
 */

#ifndef UBCC1PI_OBJECTS_EVENT_FACTORY
#define UBCC1PI_OBJECTS_EVENT_FACTORY

#include "ubcc1pi/Objects/Event.h"
#include "art/Framework/Principal/Event.h"

#include "fhiclcpp/types/Atom.h"
#include "canvas/Utilities/InputTag.h"

#include "ubcc1pi/Helpers/CollectionHelper.h"

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
         *  @brief  The configuration required to build a ubcc1pi::Event
         */
        struct Config
        {
            fhicl::Atom<bool> HasTruthInfo
            {
                fhicl::Name("HasTruthInfo"),
                fhicl::Comment("If the input events have truth info that we should look for")
            };

            fhicl::Atom<art::InputTag> MCTruthLabel
            {
                fhicl::Name("MCTruthLabel"),
                fhicl::Comment("The label for the neutrino MCTruth producer")
            };
            
            fhicl::Atom<art::InputTag> MCParticleLabel
            {
                fhicl::Name("MCParticleLabel"),
                fhicl::Comment("The label for the MCParticle producer")
            };
            
            fhicl::Atom<art::InputTag> BacktrackerLabel
            {
                fhicl::Name("BacktrackerLabel"),
                fhicl::Comment("The label for the MCParticle to hit backtracker producer")
            };
        };

        /**
         *  @brief  Populate the output event with information from the input event
         *
         *  @param  event the input event
         *  @param  config the configuration options
         *  @param  pOutputEvent the output event to populate
         */
        static void PopulateEvent(const art::Event &event, const Config &config, Event *pOutputEvent);

    private:

        /**
         *  @brief  Populate the metadata of the event
         *
         *  @param  event the input event
         *  @param  config the configuration options
         *  @param  metadata the output metadata
         */
        static void PopulateEventMetadata(const art::Event &event, const Config &config, Event::Metadata &metadata);

        /**
         *  @brief  Populate the truth information of the event
         *
         *  @param  event the input event
         *  @param  config the configuration options
         *  @param  truth the output truth information
         */
        static void PopulateEventTruthInfo(const art::Event &event, const Config &config, Event::Truth &truth);

        /**
         *  @brief  Populate the truth particle information
         *
         *  @param  mcParticle the input MCParticle
         *  @param  mcParticleToHitWeights
         *  @param  particle the output particle
         */
        static void PopulateEventTruthParticleInfo(const art::Ptr<simb::MCParticle> &mcParticle, Event::Truth::Particle &particle);
};

} // namespace ubcc1pi

#endif
