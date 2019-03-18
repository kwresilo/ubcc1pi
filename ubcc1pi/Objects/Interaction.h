/**
 *  @file  ubcc1pi/Objects/Interaction.h
 *
 *  @brief The header file for the interaction class
 */

#ifndef UBCC1PI_OBJECTS_INTERACTION
#define UBCC1PI_OBJECTS_INTERACTION

#include "ubcc1pi/Helpers/CollectionHelper.h"
#include "ubcc1pi/Helpers/TruthHelper.h"

#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"
#include "nusimdata/SimulationBase/MCParticle.h"

namespace ubcc1pi
{

/**
 *  @brief  The interaction class - represents a neutrino interaction
 */
class Interaction
{
    public:
        // No default constructor
        Interaction() = delete;

        /**
         *  @breif  Constructor
         *
         *  @param  event the art event
         *  @param  mcTruthLabel the label of the MCTruth producer
         *  @param  mcParticle the label of the MCParticle producer
         */
        Interaction(const art::Event &event, const art::InputTag &mcTruthLabel, const art::InputTag &mcParticleLabel);

        /**
         *  @breif  Get the MCParticle of the neutrino
         */
        simb::MCParticle GetNeutrino() const;
        
        /**
         *  @breif  Get the current type (charged or neutral)
         */
        simb::curr_type_ GetCCNC() const;
        
        /**
         *  @breif  Get the interaction mode (QE, RES, DIS, ...)
         */
        simb::int_type_ GetIteractionMode() const;

        /**
         *  @breif  Get the full list of MCParticles associated with this interaction
         */
        MCParticleVector GetAllMCParticles() const;

        /**
         *  @breif  Print the details of the interaction
         */
        void PrintInfo() const;

    private:
        simb::MCParticle  m_neutrino;         ///< The neutrino MCParticle
        MCParticleVector  m_allMCParticles;   ///< The MCParticles associated with the interaction

        simb::curr_type_  m_ccnc;             ///< If the iteraction is CC or NC
        simb::int_type_   m_mode;             ///< The interaction mode: QE, RES, DIS...
};

} // namespace ubcc1pi

#endif
