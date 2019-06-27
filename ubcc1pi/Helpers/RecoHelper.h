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
        // TODO fill out comments
        /**
         *  @brief  
         *
         *  @param  allPFParticles
         *
         *  @return 
         */
        static PFParticleMap GetPFParticleMap(const PFParticleVector &allPFParticles);

        /**
         *  @brief  
         *
         *  @param  allPFParticles
         *
         *  @return 
         */
        static PFParticleVector GetNeutrinos(const PFParticleVector &allPFParticles);

        /**
         *  @brief  
         *
         *  @param  allPFParticles
         *
         *  @return 
         */
        static art::Ptr<recob::PFParticle> GetNeutrino(const PFParticleVector &allPFParticles);

        /**
         *  @brief  Get the PFParticles that have been deemed the neutrino final states from the PFParticle hierarchy
         *
         *  @param  allPFParticles the input vector of all PFParticles
         */
        static PFParticleVector GetNeutrinoFinalStates(const PFParticleVector &allPFParticles);

        /**
         *  @brief  
         *
         *  @param  event
         *  @param  allPFParticles
         *  @param  vertexLabel
         *
         *  @return 
         */
        static TVector3 GetRecoNeutrinoVertex(const art::Event &event, const PFParticleVector &allPFParticles, const art::InputTag &vertexLabel);

        /**
         *  @brief  
         *
         *  @param  particle
         *  @param  pfParticleMap
         *
         *  @return 
         */
        static art::Ptr<recob::PFParticle> GetParent(const art::Ptr<recob::PFParticle> &particle, const PFParticleMap &pfParticleMap);

        /**
         *  @brief  
         *
         *  @param  particle
         *  @param  pfParticleMap
         *
         *  @return 
         */
        static PFParticleVector GetDaughters(const art::Ptr<recob::PFParticle> &particle, const PFParticleMap &pfParticleMap);

        /**
         *  @brief  
         *
         *  @param  particle
         *  @param  pfParticleMap
         *
         *  @return 
         */
        static PFParticleVector GetDownstreamParticles(const art::Ptr<recob::PFParticle> &particle, const PFParticleMap &pfParticleMap);

        /**
         *  @brief  
         *
         *  @param  particle
         *  @param  pfParticleMap
         *  @param  downstreamParticles
         */
        static void GetDownstreamParticles(const art::Ptr<recob::PFParticle> &particle, const PFParticleMap &pfParticleMap, PFParticleVector &downstreamParticles);

        /**
         *  @brief  
         *
         *  @param  slice
         *  @param  slicesToPFParticles
         *
         *  @return 
         */
        static bool IsSliceSelectedAsNu(const art::Ptr<recob::Slice> &slice, const SlicesToPFParticles &slicesToPFParticles);

        /**
         *  @brief  
         *
         *  @param  metadata
         *  @param  property
         *
         *  @return 
         */
        static bool HasMetadataValue(const art::Ptr<larpandoraobj::PFParticleMetadata> &metadata, const std::string &property);

        /**
         *  @brief  
         *
         *  @param  metadata
         *  @param  property
         *
         *  @return 
         */
        static float GetMetadataFloat(const art::Ptr<larpandoraobj::PFParticleMetadata> &metadata, const std::string &property);

        /**
         *  @brief  
         *
         *  @param  metadata
         *  @param  property
         *
         *  @return 
         */
        static int GetMetadataInt(const art::Ptr<larpandoraobj::PFParticleMetadata> &metadata, const std::string &property);

        /**
         *  @brief  
         *
         *  @param  metadata
         *  @param  property
         *
         *  @return 
         */
        static bool GetMetadataBool(const art::Ptr<larpandoraobj::PFParticleMetadata> &metadata, const std::string &property);

        /**
         *  @brief  
         *
         *  @param  metadata
         *
         *  @return 
         */
        static float GetTrackScore(const art::Ptr<larpandoraobj::PFParticleMetadata> &metadata);

        /**
         *  @brief  
         *
         *  @param  metadata
         *
         *  @return 
         */
        static float GetTopologicalScore(const art::Ptr<larpandoraobj::PFParticleMetadata> &metadata);

        /**
         *  @brief  
         *
         *  @param  hits
         *  @param  view
         *
         *  @return 
         */
        static unsigned int CountHitsInView(const HitVector &hits, const geo::View_t &view);

        /**
         *  @brief  
         *
         *  @return 
         */
        static const spacecharge::SpaceChargeService::provider_type *const GetSpaceChargeService();

        /**
         *  @brief  
         *
         *  @param  position
         *  @param  pSpaceChargeService
         *
         *  @return 
         */
        static TVector3 CorrectForSpaceCharge(const TVector3 &position, const spacecharge::SpaceChargeService::provider_type *const pSpaceChargeService);
};

} // namespace ubcc1pi

#endif
