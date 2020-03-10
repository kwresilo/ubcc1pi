/**
 *  @file  ubcc1pi/Helpers/AnalysisHelper.h
 *
 *  @brief The header file for the analysis helper class
 */

#ifndef UBCC1PI_HELPERS_ANALYSIS_HELPER
#define UBCC1PI_HELPERS_ANALYSIS_HELPER

#include <TVector3.h>

#include "ubcc1pi/Helpers/CollectionHelper.h"
#include "ubcc1pi/Helpers/TruthHelper.h"
#include "ubcc1pi/Helpers/RecoHelper.h"
#include "ubcc1pi/Helpers/BacktrackHelper.h"

namespace ubcc1pi
{

/**
 *  @brief  The analysis helper class
 *
 *          Generic functions that are relevent for many analysis tasks, e.g. check if a point is in the fiducial volume
 */
class AnalysisHelper
{
    public:
        /**
         *  @brief  Check if a position is contained within the TPC with some border margin
         *
         *  @param  point the point to check
         *  @param  border the size of the border around the TPC volume
         */
        static bool IsContained(const TVector3 &point, const float border);

        /**
         *  @brief  Check if a position lies in the fiducial volume
         *
         *  @param  point the point to check
         */
        static bool IsFiducial(const TVector3 &point);

        /**
         *  @brief  Check if a MCParticle passes the momentum threshold to be considered by the analysis
         *
         *  @param  particle the particle to check
         */
        static bool PassesMomentumThreshold(const art::Ptr<simb::MCParticle> &particle);

        /**
         *  @brief  Check if a MCParticle is considered visible to the detector
         *
         *  @param  particle the particle to check
         */
        static bool IsVisible(const art::Ptr<simb::MCParticle> &particle);

        /**
         *  @brief  Get the primary MCParticles from the input interaction 
         *
         *  @param  interaction the interaction
         *
         *  @return primary MCParticles
         */
        static MCParticleVector GetPrimaryParticles(const TruthHelper::Interaction &interaction);

        // TODO fill in final doxygen!
        /**
         *  @brief  
         *
         *  @param  interaction
         *
         *  @return 
         */
        static MCParticleVector GetVisibleFinalStates(const TruthHelper::Interaction &interaction);

        /**
         *  @brief  
         *
         *  @param  interaction
         *
         *  @return 
         */
        static MCParticleVector GetReconstructableFinalStates(const TruthHelper::Interaction &interaction);

        /**
         *  @brief  
         *
         *  @param  interaction
         *
         *  @return 
         */
        static bool IsNeutrinoVertexFiducial(const TruthHelper::Interaction &interaction);

        /**
         *  @brief  
         *
         *  @param  interaction
         *
         *  @return 
         */
        static bool IsCC1PiSignal(const TruthHelper::Interaction &interaction);

        /**
         *  @brief  
         *
         *  @param  particles
         *  @param  pdg
         *
         *  @return 
         */
        static unsigned int CountParticlesWithPDG(const MCParticleVector &particles, const int pdg);
        
        /**
         *  @brief  
         *
         *  @param  event
         *  @param  mcTruthLabel
         *  @param  mcParticleLabel
         *  @param  backtrackerLabel
         *  @param  pfParticleLabel
         *
         *  @return 
         */
        static BacktrackHelper::BacktrackerData GetBacktrackerData(const art::Event &event, const art::InputTag &mcTruthLabel, const art::InputTag &mcParticleLabel, const art::InputTag &backtrackerLabel, const art::InputTag &pfParticleLabel);

        /**
         *  @brief  
         *
         *  @param  event
         *  @param  mcTruthLabel
         *  @param  mcParticleLabel
         *  @param  backtrackerLabel
         *  @param  pfParticleLabel
         *  @param  sliceLabel
         *  @param  hitLabel
         *
         *  @return 
         */
        static BacktrackHelper::SliceMetadata GetSliceMetadata(const art::Event &event, const art::InputTag &mcTruthLabel, const art::InputTag &mcParticleLabel, const art::InputTag &backtrackerLabel, const art::InputTag &pfParticleLabel, const art::InputTag &sliceLabel, const art::InputTag &hitLabel);

    private:
        /**
         *  @brief  Get the first visible MCParticles in the hierarchy - these are the final state particles we expect to produce tracks and showers
         *          E.g. Suppose a neutrino iteraction produced a pi0, this function would return the decay photons (since the pi0 itself doesn't produce any tracks / showers)
         *
         *  @param  allMCParticles the input MCParticles
         *
         *  @return the visible final state MCParticles
         */
        static MCParticleVector GetVisibleFinalStates(const MCParticleVector &allMCParticles);

        /**
         *  @brief  Get the first visible MCParticles downstream of the input particle
         *
         *  @param  particle the input particle
         *  @param  mcParticleMap the mapping from TrackId -> MCParticle
         *  @param  finalStates the output vector of final state particles to populate
         */
        static void GetVisibleFinalStates(const art::Ptr<simb::MCParticle> &particle, const MCParticleMap &mcParticleMap, MCParticleVector &finalStates);
};

} // namespace ubcc1pi

#endif
