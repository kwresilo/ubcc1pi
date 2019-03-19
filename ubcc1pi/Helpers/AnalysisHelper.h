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

        // TODO add doxygen comments
        MCParticleVector GetPrimaryParticles(const TruthHelper::Interaction &interaction);
        MCParticleVector GetVisibleFinalStates(const TruthHelper::Interaction &interaction);
        MCParticleVector GetReconstructableFinalStates(const TruthHelper::Interaction &interaction);
        bool IsNeutrinoVertexFiducial(const TruthHelper::Interaction &interaction);
        bool IsCC1PiSignal(const TruthHelper::Interaction &interaction);

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
