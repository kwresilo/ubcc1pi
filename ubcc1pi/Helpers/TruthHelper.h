/**
 *  @file  ubcc1pi/Helpers/TruthHelper.h
 *
 *  @brief The header file for the truth helper class
 */

#ifndef UBCC1PI_HELPERS_TRUTH_HELPER
#define UBCC1PI_HELPERS_TRUTH_HELPER

#include "art/Framework/Principal/Event.h"
#include "canvas/Utilities/InputTag.h"
#include "nusimdata/SimulationBase/MCTruth.h"

#include "ubcc1pi/Helpers/CollectionHelper.h"

namespace ubcc1pi
{

/**
 *  @brief  The truth helper class
 */
class TruthHelper
{
    public:
        /**
         *  @brief  Get the MCTruth with the most energetic MCNeutrino
         *
         *  @param  event the art event
         *  @param  mcTruthLabel the label of the MCTruth producer
         *
         *  @return the MCTruth
         */
        static art::Ptr<simb::MCTruth> GetNeutrinoMCTruth(const art::Event &event, const art::InputTag &mcTruthLabel);

        /**
         *  @brief  Get the MCParticles associated to the input MCTruth
         *
         *  @param  event the art event
         *  @param  mcTruthLabel the label of the MCTruth producer
         *  @param  mcParticleLabel the label of the MCParticle producer
         *  @param  mcTruth the MCTruth object
         *
         *  @return the associated MCParticles
         */
        static MCParticleVector GetMCParticlesFromMCTruth(const art::Event &event, const art::InputTag &mcTruthLabel, const art::InputTag &mcParticleLabel, const art::Ptr<simb::MCTruth> &mcTruth);

        /**
         *  @brief  Get the primary MCParticles
         *
         *  @param  allMCParticles the input MCParticles
         *
         *  @return the primary MCParticles in the input list
         */
        static MCParticleVector GetPrimaryMCParticles(const MCParticleVector &allMCParticles);
};

} // namespace ubcc1pi

#endif
