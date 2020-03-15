/**
 *  @file  ubcc1pi/Producers/HitSlicer.h
 *
 *  @brief The header file for the hit slicer producer, this module makes a copy of only the hits in the neutrino slice
 */

#ifndef UBCC1PI_PRODUCERS_HIT_SLICER
#define UBCC1PI_PRODUCERS_HIT_SLICER

#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDProducer.h"

#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Table.h"
#include "canvas/Utilities/InputTag.h"

namespace ubcc1pi
{

/**
 *  @brief  The pfparticle hierarchy producer class
 */
class HitSlicer : public art::EDProducer
{
    public:
        /**
         *  @brief  The configuration structure
         */
        struct Config
        {
            /**
             *  @brief  The slice label
             */
            fhicl::Atom<art::InputTag> SliceLabel
            {
                fhicl::Name("SliceLabel"),
                fhicl::Comment("The label for the slice producer")
            };
            
            /**
             *  @brief  The slice to hits label
             */
            fhicl::Atom<art::InputTag> SliceToHitsLabel
            {
                fhicl::Name("SliceToHitsLabel"),
                fhicl::Comment("The label for the slice to hits association producer")
            };
            
            /**
             *  @brief  The slice to PFParticles label
             */
            fhicl::Atom<art::InputTag> SliceToPFParticlesLabel
            {
                fhicl::Name("SliceToPFParticlesLabel"),
                fhicl::Comment("The label for the slice to pfparticle association producer")
            };
            
            /**
             *  @brief  The MCParticle label
             */
            fhicl::Atom<art::InputTag> MCParticleLabel
            {
                fhicl::Name("MCParticleLabel"),
                fhicl::Comment("The label for the MCParticle producer label")
            };
            
            /**
             *  @brief  The backtracker label
             */
            fhicl::Atom<art::InputTag> BacktrackerLabel
            {
                fhicl::Name("BacktrackerLabel"),
                fhicl::Comment("The label for the MCParticle to hits producer label")
            };
        };
        
        /**
         *  @brief  Constructor
         *
         *  @param  config the set of input fhicl parameters
         */
        HitSlicer(const art::EDProducer::Table<Config> &config);

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

DEFINE_ART_MODULE(ubcc1pi::HitSlicer)

#endif
