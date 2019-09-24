/**
 *  @file  ubcc1pi/Analyzers/EnergyEstimator.h
 *
 *  @brief The header file for the energy estimator study analyzer.
 */

#ifndef UBCC1PI_ANALYZERS_ENERGY_ESTIMATOR_STUDY
#define UBCC1PI_ANALYZERS_ENERGY_ESTIMATOR_STUDY

#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDAnalyzer.h"

#include "art/Framework/Services/Optional/TFileService.h"
#include <TTree.h>

#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Table.h"
#include "canvas/Utilities/InputTag.h"

#include "ubcc1pi/Helpers/CollectionHelper.h"
#include "ubcc1pi/Helpers/BacktrackHelper.h"
#include "ubcc1pi/Helpers/RecoHelper.h"

namespace ubcc1pi
{

/**
 *  @brief  The energy estimator study class
 */
class EnergyEstimator : public art::EDAnalyzer
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
                fhicl::Comment("The label for the MCParticle to Hit backtracker producer")
            };
            
            fhicl::Atom<art::InputTag> PFParticleLabel
            {
                fhicl::Name("PFParticleLabel"),
                fhicl::Comment("The label for the PFParticle producer (Pandora)")
            };
            
            fhicl::Atom<art::InputTag> TrackLabel
            {
                fhicl::Name("TrackLabel"),
                fhicl::Comment("The label for the Track producer")
            };
            
            fhicl::Atom<art::InputTag> PIDLabel
            {
                fhicl::Name("PIDLabel"),
                fhicl::Comment("The label for the PID producer")
            };
            
            fhicl::Atom<art::InputTag> CalorimetryLabel
            {
                fhicl::Name("CalorimetryLabel"),
                fhicl::Comment("The label for the calorimetry producer")
            };
        };

        /**
         *  @brief  The output particle level structure
         */
        struct OutputParticle
        {
            // Event metadata
            int          m_run;                   ///< The run number
            int          m_subRun;                ///< The subrun number
            int          m_event;                 ///< The event number
            bool         m_isSignal;              ///< If the event is a true CC1Pi signal
        
            // Matched True Particle details
            bool         m_hasMatchedMCParticle;  ///< If the PFParticle matches to any neutrino MCParticles
            float        m_trueMatchPurity;       ///< Match purity to the true particle
            float        m_trueMatchCompleteness; ///< Match completeness to the true particle

            // MCParticle details
            int          m_truePdgCode;           ///< The particle PDG code
            float        m_trueMomentum;          ///< The particle momentum
            float        m_trueEndMomentum;       ///< The particle end momentum
            float        m_trueKE;                ///< The particle kinetic energy
            bool         m_isStopping;            ///< If the particle is stopping (zero end momentum)
            TVector3     m_trueStart;             ///< The true start position
            TVector3     m_trueEnd;               ///< The true end position
            bool         m_isContained;           ///< If the particle is contained
            int          m_nScatters;             ///< The number of elasic scatters
            bool         m_isGolden;              ///< If the particle is golden = contained, stopping & no scatters
            float        m_trueRange;             ///< The true range of the particle
            float        m_trueEndPointDist;      ///< The true start-end distance

            // Reco particle variables
            bool         m_hasTrack;              ///< Has an associated track
            TVector3     m_start;                 ///< The start position (SCE corrected)
            TVector3     m_end;                   ///< The end position (SCE corrected)
            float        m_length;                ///< The track-length (Not SCE corrected)
            float        m_range;                 ///< The range of the particle (SCE corrected)
            float        m_endPointDist;          ///< The start-end distance (SCE corrected)
        };

        /**
         *  @brief  Constructor
         *
         *  @param  config the set of input fhicl parameters
         */
        EnergyEstimator(const art::EDAnalyzer::Table<Config> &config);

        /**
         *  @brief  Analyze an event
         *
         *  @param  event the event to analyze
         */
        void analyze(const art::Event &event);

        /**
         *  @brief  Get the range of an MCParticle by following all of it's trajectory points
         *
         *  @param  mcParticle the MCParticle who's range we wish to retrieve
         *
         *  @return the range of the MCParticle
         */
        float GetRange(const art::Ptr<simb::MCParticle> &mcParticle) const;

        /**
         *  @brief  Get the straight line distance from the MCParticle's start and end points
         *
         *  @param  mcParticle the MCParticle who's endpoint distance we wish to retrieve
         *
         *  @return the endpoint distance
         */
        float GetEndpointDistance(const art::Ptr<simb::MCParticle> &mcParticle) const;

        /**
         *  @brief  Get the range of a track by following it's trajectory points
         *
         *  @param  track the track
         *  @param  pSpaceChargeService the spacecharge service
         *
         *  @return the range
         */
        float GetRange(const art::Ptr<recob::Track> &track, const spacecharge::SpaceChargeService::provider_type *const pSpaceChargeService) const;

        /**
         *  @brief  Get the straight line distance from a track's by start and end trajectory points
         *
         *  @param  track the track
         *  @param  pSpaceChargeService the spacecharge service
         *
         *  @return the endpoint distance
         */
        float GetEndpointDistance(const art::Ptr<recob::Track> &track, const spacecharge::SpaceChargeService::provider_type *const pSpaceChargeService) const;

        /**
         *  @brief  Get the valid points in a track's trajectory
         *
         *  @param  track the track
         *
         *  @return the vector of valid points
         */
        std::vector<size_t> GetValidPoints(const art::Ptr<recob::Track> &track) const;

        /**
         *  @brief  Reset the output object to dummy values, ready for the next PFParticle to fill
         */
        void ResetParticleTree();


    private:
        art::EDAnalyzer::Table<Config>  m_config;         ///< The FHiCL configuration options
        
        TTree                          *m_pParticleTree;  ///< The output tree for all particle level data
        OutputParticle                  m_outputParticle; ///< The output particle-level object
};

} // namespace ubcc1pi

DEFINE_ART_MODULE(ubcc1pi::EnergyEstimator)

#endif
