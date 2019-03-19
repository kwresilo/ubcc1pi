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

#include <unordered_map>

namespace ubcc1pi
{

/**
 *  @brief  A map intended to hold the mapping from MCParticle.TrackId() -> MCParticle - for navigation of the MCParticle hierarchy 
 */
typedef std::unordered_map<int, art::Ptr<simb::MCParticle> > MCParticleMap;

/**
 *  @brief  The truth helper class
 */
class TruthHelper
{
    public:
        /**
         *  @brief  The interaction class - represents a neutrino interaction
         */
        class Interaction
        {
            public:
                // No default constructor
                Interaction() = delete;

                /**
                 *  @brief  Constructor
                 *
                 *  @param  event the art event
                 *  @param  mcTruthLabel the label of the MCTruth producer
                 *  @param  mcParticleLabel the label of the MCParticle producer
                 */
                Interaction(const art::Event &event, const art::InputTag &mcTruthLabel, const art::InputTag &mcParticleLabel);

                /**
                 *  @brief  Get the MCParticle of the neutrino
                 */
                simb::MCParticle GetNeutrino() const;

                /**
                 *  @brief  Get the current type (charged or neutral)
                 */
                simb::curr_type_ GetCCNC() const;

                /**
                 *  @brief  Get the interaction mode (QE, RES, DIS, ...)
                 */
                simb::int_type_ GetIteractionMode() const;

                /**
                 *  @brief  Get the full list of MCParticles associated with this interaction
                 */
                MCParticleVector GetAllMCParticles() const;

                /**
                 *  @brief  Print the details of the interaction
                 */
                //void PrintInfo() const;

            private:
                simb::MCParticle  m_neutrino;         ///< The neutrino MCParticle
                MCParticleVector  m_allMCParticles;   ///< The MCParticles associated with the interaction

                simb::curr_type_  m_ccnc;             ///< If the iteraction is CC or NC
                simb::int_type_   m_mode;             ///< The interaction mode: QE, RES, DIS...
        };

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
         *  @brief  Get the primary MCParticles - the immediate products of the neutrino interaction leaving the nucleus
         *
         *  @param  allMCParticles the input MCParticles
         *
         *  @return the primary MCParticles in the input list
         */
        static MCParticleVector GetPrimaryMCParticles(const MCParticleVector &allMCParticles);

        /**
         *  @brief  Check if an MCParticle is primary
         *
         *  @param  particle the particle to check
         */
        static bool IsPrimary(const art::Ptr<simb::MCParticle> &particle);

        /**
         *  @brief  Get the mapping from MCParticle.TrackId() -> MCParticle
         *
         *  @param  allMCParticles the input MCParticles
         *
         *  @return the output MCParticleMap map
         */
        static MCParticleMap GetMCParticleMap(const MCParticleVector &allMCParticles);

        /**
         *  @brief  Get the mother MCParticle of the input particle
         *
         *  @param  particle the particle
         *  @param  mcParticleMap the mapping from TrackId -> MCParticle
         *
         *  @return the mother MCParticle
         */
        static art::Ptr<simb::MCParticle> GetMother(const art::Ptr<simb::MCParticle> &particle, const MCParticleMap &mcParticleMap);
        
        /**
         *  @brief  Get the daughter MCParticles of the input particle
         *
         *  @param  particle the particle
         *  @param  mcParticleMap the mapping from TrackId -> MCParticle
         *
         *  @return the daughter MCParticles
         */
        static MCParticleVector GetDaughters(const art::Ptr<simb::MCParticle> &particle, const MCParticleMap &mcParticleMap);
};

} // namespace ubcc1pi

#endif
