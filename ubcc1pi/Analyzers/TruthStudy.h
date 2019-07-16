/**
 *  @file  ubcc1pi/Analyzers/TruthStudy.h
 *
 *  @brief The header file for the truth study analyzer.
 */

#ifndef UBCC1PI_ANALYZERS_TRUTH_STUDY
#define UBCC1PI_ANALYZERS_TRUTH_STUDY

#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDAnalyzer.h"

#include "art/Framework/Services/Optional/TFileService.h"
#include <TTree.h>

#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Table.h"
#include "canvas/Utilities/InputTag.h"

#include "ubcc1pi/Helpers/TruthHelper.h"
#include "ubcc1pi/Helpers/CollectionHelper.h"

namespace ubcc1pi
{

/**
 *  @brief  The truth study class
 */
class TruthStudy : public art::EDAnalyzer
{
    public:
        /**
         *  @brief  The configuration structure
         */
        struct Config
        {
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
         *  @brief  The interaction output structure
         */
        struct InteractionOutput
        {
            // Event metadata
            int          m_run;          ///< The run number
            int          m_subRun;       ///< The subrun number
            int          m_event;        ///< The event number
            // Interaction information
            std::string  m_interaction;  ///< The truth interaction type
            bool         m_isNuFiducial; ///< If the neutrino is fiducial
            bool         m_isSignal;     ///< If the event is a CC1Pi signal event
            int          m_nMuMinus;     ///< The number of mu-
            int          m_nMuPlus;      ///< The number of mu+
            int          m_nPiPlus;      ///< The number of pi+
            int          m_nPiMinus;     ///< The number of pi-
            int          m_nKPlus;       ///< The number of K+
            int          m_nKMinus;      ///< The number of K-
            int          m_nProton;      ///< The number of p
            int          m_nNeutron;     ///< The number of n
            int          m_nPhoton;      ///< The number of gamma
            int          m_nTotal;       ///< The total number of particles
            // Kinematics
            float        m_nuE;          ///< The neutrino energy
        };
        
        /**
         *  @brief  The signal output structure
         */
        struct SignalOutput
        {
            // Event metadata
            int          m_run;                                 ///< The run number
            int          m_subRun;                              ///< The subrun number
            int          m_event;                               ///< The event number
            // Interaction information                          
            std::string  m_interaction;                         ///< The truth interaction type
            int          m_nProton;                             ///< The number of p
            // Kinematics                                       
            float        m_nuE;                                 ///< The neutrino energy
            float        m_muMomX;                              ///< The muon momentum in x
            float        m_muMomY;                              ///< The muon momentum in y
            float        m_muMomZ;                              ///< The muon momentum in z
            float        m_piMomX;                              ///< The pion momentum in x
            float        m_piMomY;                              ///< The pion momentum in y
            float        m_piMomZ;                              ///< The pion momentum in z
            // Pion interactions
            int                 m_piNElasticScatters;           ///< The number of elastic scatters
            int                 m_piNInelasticScatters;         ///< The number of elastic scatters
            int                 m_piNScatters;                  ///< The total number of scatters
            std::vector<bool>   m_piScatterIsElasticVect;       ///< Vector of bools one entry per scatter: was the scatter elastic
            std::vector<float>  m_piScatterCosThetaVect;        ///< Vector of cosine-scattering angles, one entry per scatter
            std::vector<float>  m_piScatterMomFracLostVect;     ///< Vector of fraction momentum lost by the pion, one entry per scatter
            float               m_piInitalMom;                  ///< The total initial pion momentum
            float               m_piScatteredMom;               ///< The total pion momentum after the final scatter
            float               m_piFinalMom;                   ///< The total final pion momentum just before the end-state interaction
            int                 m_piEndState;                   ///< The end-state interaction enumeration
            float               m_piEndStateProductsHitWeightU; ///< The total hit weight of all end-state products - U
            float               m_piEndStateProductsHitWeightV; ///< The total hit weight of all end-state products - V
            float               m_piEndStateProductsHitWeightW; ///< The total hit weight of all end-state products - W
        };

        /**
         *  @brief  Constructor
         *
         *  @param  config the set of input fhicl parameters
         */
        TruthStudy(const art::EDAnalyzer::Table<Config> &config);

        /**
         *  @brief  Analyze an event
         *
         *  @param  event the event to analyze
         */
        void analyze(const art::Event &event);

    private:

        /**
         *  @brief  Class describing the scatter of a pion
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
         *  @brief  Class describing the end-state of a pion
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
                 *  @param  finalPionMomentum the momentum of the pion just before the end state interaction in the lab frame
                 *  @param  products the final state products of the interaction
                 */
                EndState(const Type type, const TVector3 &finalPionMomentum, const MCParticleVector &products);

                /**
                 *  @brief  The type of end state enumeration
                 */
                enum Type : unsigned int
                {
                    None,                  ///< The pion just stops after running out of energy
                    DecayToMuon,           ///< The pion decays to a muon (which could then in turn decay to a Michel)
                    InelasticAbsorption,   ///< The pion interacts inelastically with a nucleus and is absorbed
                    Pi0ChargeExchange,     ///< The pion interacts inelastically with a nucleus and a Pi0 is in the final state
                    Other                  ///< None of the above
                };

                Type                m_type;                 ///< The type of end state
                TVector3            m_finalPionMomentum;    ///< The momentum of the pion just before the end state interaction in the lab frame
                MCParticleVector    m_products;             ///< The final state products of the interaction 
        };

        /**
         *  @brief  Navigate the pion's heirarchy to follow it's elastic and inelastic scatters 
         *
         *  @param  pion the pion to follow
         *  @param  mcParticleMap the mapping from MCParticle ID to art Ptr
         *  @param  scatters the output vector of all scatters
         *  @param  scatteredPion the final pion MCParticle after following all scatters, this is the one that will stop, interact or decay
         */
        void FollowScatters(const art::Ptr<simb::MCParticle> &pion, const MCParticleMap &mcParticleMap, ScatterVector &scatters, art::Ptr<simb::MCParticle> &scatteredPion) const;

        /**
         *  @brief  Get the details of an elastic scatter from the incident and target particle
         *
         *  @param  incident the incident pion that has scattered elastically
         *  @param  target the target MCParticle on which the pion elastically scattered
         *
         *  @return The scatter details
         */
        Scatter GetElasticScatter(const art::Ptr<simb::MCParticle> &incident, const art::Ptr<simb::MCParticle> &target) const;

        /**
         *  @brief  Get the end state of the pion
         *
         *  @param  pion the pion MCParticle after all scatters
         *  @param  mcParticleMap the mapping from MCParticle ID to art Ptr
         *
         *  @return the end state
         */
        EndState GetEndState(const art::Ptr<simb::MCParticle> &pion, const MCParticleMap &mcParticleMap) const;

        /**
         *  @brief  Print information about the hierarchy of MCParticles that starts with the input
         *
         *  @param  particle the seed MCParticle
         *  @param  mcParticleMap the mapping from MCParticle ID to art Ptr
         *  @param  event the art event
         */
        void PrintHierarchy(const art::Ptr<simb::MCParticle> &particle, const MCParticleMap &mcParticleMap, const art::Event &event) const;
            
        art::EDAnalyzer::Table<Config>  m_config;            ///< The FHiCL configuration options
        TTree                          *m_pInteractionTree;  ///< The output tree for all interactions
        TTree                          *m_pSignalTree;       ///< The output tree for signal interactions

        InteractionOutput               m_interactionOutput; ///< The output objects for all interactions
        SignalOutput                    m_signalOutput;      ///< The output objects for signal interactions
};

} // namespace ubcc1pi

DEFINE_ART_MODULE(ubcc1pi::TruthStudy)

#endif
