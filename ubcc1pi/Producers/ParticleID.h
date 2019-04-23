/**
 *  @file  ubcc1pi/Producers/ParticleID.h
 *
 *  @brief The header file for the particle ID producer.
 */

#ifndef UBCC1PI_PRODUCERS_PARTICLE_ID
#define UBCC1PI_PRODUCERS_PARTICLE_ID

#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDProducer.h"

#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Table.h"
#include "canvas/Utilities/InputTag.h"

namespace ubcc1pi
{

/**
 *  @brief  The particle ID producer class
 */
class ParticleID : public art::EDProducer
{
    public:
        /**
         *  @brief  The configuration structure
         */
        struct Config
        {
            fhicl::Atom<art::InputTag> PFParticleLabel
            {
                fhicl::Name("PFParticleLabel"),
                fhicl::Comment("The label for the PFParticle producer (Pandora)")
            };
            
            fhicl::Atom<art::InputTag> HitLabel
            {
                fhicl::Name("HitLabel"),
                fhicl::Comment("The label for the Hit producer")
            };
        };

        /**
         *  @brief  Constructor
         *
         *  @param  config the set of input fhicl parameters
         */
        ParticleID(const art::EDProducer::Table<Config> &config);

        /**
         *  @brief  Produce collections in the event record
         *
         *  @param  event the event to add to
         */
        void produce(art::Event &event);

    private:
            
        art::EDProducer::Table<Config>  m_config;          ///< The FHiCL configuration options
};

} // namespace ubcc1pi

DEFINE_ART_MODULE(ubcc1pi::ParticleID)

#endif
