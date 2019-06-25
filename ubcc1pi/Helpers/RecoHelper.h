/**
 *  @file  ubcc1pi/Helpers/RecoHelper.h
 *
 *  @brief The header file for the reco helper class
 */

#ifndef UBCC1PI_HELPERS_RECO_HELPER
#define UBCC1PI_HELPERS_RECO_HELPER

#include "ubcc1pi/Helpers/CollectionHelper.h"

#include "larevt/SpaceChargeServices/SpaceChargeService.h"

#include <unordered_map>

namespace ubcc1pi
{

/**
 *  @brief  A map intended to hold the mapping from PFParticle.TrackId() -> PFParticle - for navigation of the PFParticle hierarchy 
 */
typedef std::unordered_map<int, art::Ptr<recob::PFParticle> > PFParticleMap;

/**
 *  @brief  The reco helper class
 */
class RecoHelper
{
    public:
        // TODO doxygen comments
        static PFParticleMap GetPFParticleMap(const PFParticleVector &allPFParticles);
        static PFParticleVector GetNeutrinos(const PFParticleVector &allPFParticles);
        static art::Ptr<recob::PFParticle> GetNeutrino(const PFParticleVector &allPFParticles);

        /**
         *  @brief  Get the PFParticles that have been deemed the neutrino final states from the PFParticle hierarchy
         *
         *  @param  allPFParticles the input vector of all PFParticles
         */
        static PFParticleVector GetNeutrinoFinalStates(const PFParticleVector &allPFParticles);

        static TVector3 GetRecoNeutrinoVertex(const art::Event &event, const PFParticleVector &allPFParticles, const art::InputTag &vertexLabel);

        static art::Ptr<recob::PFParticle> GetParent(const art::Ptr<recob::PFParticle> &particle, const PFParticleMap &pfParticleMap);
        static PFParticleVector GetDaughters(const art::Ptr<recob::PFParticle> &particle, const PFParticleMap &pfParticleMap);
        static PFParticleVector GetDownstreamParticles(const art::Ptr<recob::PFParticle> &particle, const PFParticleMap &pfParticleMap);
        static void GetDownstreamParticles(const art::Ptr<recob::PFParticle> &particle, const PFParticleMap &pfParticleMap, PFParticleVector &downstreamParticles);
        static bool IsSliceSelectedAsNu(const art::Ptr<recob::Slice> &slice, const SlicesToPFParticles &slicesToPFParticles);

        static bool HasMetadataValue(const art::Ptr<larpandoraobj::PFParticleMetadata> &metadata, const std::string &property);
        
        static float GetMetadataFloat(const art::Ptr<larpandoraobj::PFParticleMetadata> &metadata, const std::string &property);
        static int GetMetadataInt(const art::Ptr<larpandoraobj::PFParticleMetadata> &metadata, const std::string &property);
        static bool GetMetadataBool(const art::Ptr<larpandoraobj::PFParticleMetadata> &metadata, const std::string &property);

        static float GetTrackScore(const art::Ptr<larpandoraobj::PFParticleMetadata> &metadata);
        static float GetTopologicalScore(const art::Ptr<larpandoraobj::PFParticleMetadata> &metadata);
        static unsigned int CountHitsInView(const HitVector &hits, const geo::View_t &view);

        static const spacecharge::SpaceChargeService::provider_type *const GetSpaceChargeService();
        static TVector3 CorrectForSpaceCharge(const TVector3 &position, const spacecharge::SpaceChargeService::provider_type *const pSpaceChargeService);
};

} // namespace ubcc1pi

#endif
