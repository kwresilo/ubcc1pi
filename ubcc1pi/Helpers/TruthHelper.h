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

                /**
                 *  @brief  No default constructor
                 */
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
                 *  @brief  Get the neutrino MCParticle
                 *
                 *  @return neutrino MCParticle
                 */
                simb::MCParticle GetNeutrino() const;

                /**
                 *  @brief  Get the interaction type CC vs. NC (charged current vs. neutral current)
                 *
                 *  @return the interaction type
                 */
                simb::curr_type_ GetCCNC() const;

                /**
                 *  @brief  Get the interaction mode (QE, RES, DIS, ...)
                 *
                 *  @return the interaction mode
                 */
                simb::int_type_ GetInteractionMode() const;

                /**
                 *  @brief  Get the full list of MCParticles associated with this interaction
                 *
                 *  @return all MCParticles
                 */
                MCParticleVector GetAllMCParticles() const;

            private:
                simb::MCParticle  m_neutrino;         ///< The neutrino MCParticle
                MCParticleVector  m_allMCParticles;   ///< The MCParticles associated with the interaction

                simb::curr_type_  m_ccnc;             ///< If the iteraction is CC or NC
                simb::int_type_   m_mode;             ///< The interaction mode: QE, RES, DIS...
        };

        /**
         *  @brief  Class describing the scatter of a particle
         */
        class Scatter
        {
            public:
                /**
                 *  @brief  Default constructor
                 *
                 *  @param  isElastic if the scatter was elastic
                 *  @param  initialMomentum the momentum of the incident particle just before the scatter in the lab frame
                 *  @param  finalMomentum the momentum of the incident particle just after the scatter in the lab frame
                 *  @param  finalParticle the incident particle tracked through to the final state
                 *  @param  products the interaction products of the scatter
                 */
                Scatter(const bool isElastic, const TVector3 initialMomentum, const TVector3 finalMomentum, const art::Ptr<simb::MCParticle> finalParticle, const MCParticleVector products);

                /**
                 *  @brief  Get the cosine of the scattering angle
                 *
                 *  @return cos(theta)
                 */
                float GetScatteringCosTheta() const;

                /**
                 *  @brief  Get the fractional momentum lost (Pi - Pf) / Pi
                 *
                 *  @return fractional momentum lost
                 */
                float GetMomentumFractionLost() const;

                bool                          m_isElastic;           ///< If the scatter was elastic
                TVector3                      m_initialMomentum;     ///< The momentum of the incident particle just before the scatter in the lab frame
                TVector3                      m_finalMomentum;       ///< The momentum of the incident particle just after the scatter in the lab frame
                art::Ptr<simb::MCParticle>    m_finalParticle;       ///< The incident particle tracked through to the final state
                MCParticleVector              m_products;            ///< The interaction products of the scatter
        };

        typedef std::vector<Scatter> ScatterVector; ///< A vector of scatter details

        /**
         *  @brief  Class describing the end-state of a particle
         */
        class EndState
        {
            public:
                // Forward declaration
                enum Type : unsigned int;

                /**
                 *  @brief  Default constructor
                 *
                 *  @param  type the end state type
                 *  @param  finalParticleMomentum the momentum of the particle just before the end state interaction in the lab frame
                 *  @param  products the final state products of the interaction
                 */
                EndState(const Type type, const TVector3 &finalParticleMomentum, const MCParticleVector &products);

                /**
                 *  @brief  The type of end state enumeration
                 */
                enum Type : unsigned int
                {
                    None,                  ///< The particle just stops after running out of energy
                    DecayToMuon,           ///< The particle decays to a muon (which could then in turn decay to a Michel)
                    InelasticAbsorption,   ///< The particle interacts inelastically with a nucleus and is absorbed
                    Pi0ChargeExchange,     ///< The particle interacts inelastically with a nucleus and a Pi0 is in the final state
                    Other                  ///< None of the above
                };

                Type                m_type;                     ///< The type of end state
                TVector3            m_finalParticleMomentum;    ///< The momentum of the particle just before the end state interaction in the lab frame
                MCParticleVector    m_products;                 ///< The final state products of the interaction
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
         *
         *  @return if the MCParticle is primary
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

        /**
         *  @brief  Get the MCParticles downstream of the input particle
         *
         *  @param  particle the input MCParticle
         *  @param  mcParticleMap the mapping from TrackId -> MCParticle
         *
         *  @return the downstream MCParticles
         */
        static MCParticleVector GetDownstreamParticles(const art::Ptr<simb::MCParticle> &particle, const MCParticleMap &mcParticleMap);

        /**
         *  @brief  Get the MCParticles downstream of the input particle
         *
         *  @param  particle the input MCParticle
         *  @param  mcParticleMap the mapping from TrackId -> MCParticle
         *  @param  downstreamParticles the output vector of downstream particles
         */
        static void GetDownstreamParticles(const art::Ptr<simb::MCParticle> &particle, const MCParticleMap &mcParticleMap, MCParticleVector &downstreamParticles);

        /**
         *  @brief  Navigate the particle heirarchy to follow it's elastic and inelastic scatters
         *
         *  @param  particle the particle to follow
         *  @param  mcParticleMap the mapping from MCParticle ID to art Ptr
         *  @param  scatters the output vector of all scatters
         *  @param  scatteredParticle the MCParticle after following all scatters, this is the one that will stop, interact or decay
         */
        static void FollowScatters(const art::Ptr<simb::MCParticle> &particle, const MCParticleMap &mcParticleMap, ScatterVector &scatters, art::Ptr<simb::MCParticle> &scatteredParticle);

        /**
         *  @brief  Get the details of an elastic scatter from the incident and target particle
         *
         *  @param  incident the incident particle that has scattered elastically
         *  @param  target the target MCParticle on which the particle elastically scattered
         *
         *  @return The scatter details
         */
        static Scatter GetElasticScatter(const art::Ptr<simb::MCParticle> &incident, const art::Ptr<simb::MCParticle> &target);

        /**
         *  @brief  Get the end state of an MCParticle
         *
         *  @param  particle the MCParticle after all scatters
         *  @param  mcParticleMap the mapping from MCParticle ID to art Ptr
         *
         *  @return the end state
         */
        static EndState GetEndState(const art::Ptr<simb::MCParticle> &particle, const MCParticleMap &mcParticleMap);
};

} // namespace ubcc1pi

#endif
