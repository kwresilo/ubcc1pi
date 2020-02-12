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
#include <functional>

namespace ubcc1pi
{

/**
 *  @brief  A map intended to hold the mapping from PFParticle.Self() -> PFParticle - for navigation of the PFParticle hierarchy 
 */
typedef std::unordered_map<int, art::Ptr<recob::PFParticle> > PFParticleMap;

/**
 *  @brief  The reco helper class
 */
class RecoHelper
{
    public:
        /**
         *  @brief  Get the mapping from PFParticle.Self() -> PFParticle
         *
         *  @param  allPFParticles input vector of all PFParticles in the event
         *
         *  @return the PFParticle map
         */
        static PFParticleMap GetPFParticleMap(const PFParticleVector &allPFParticles);

        /**
         *  @brief  Check if the input PFParticle is a reconstructed neutrino 
         *
         *  @param  pfParticle the input PFParticle
         *
         *  @return if the particle is a neutrino
         */
        static bool IsNeutrino(const art::Ptr<recob::PFParticle> &pfParticle);

        /**
         *  @brief  Get all neutrino PFParticles in the input list
         *
         *  @param  allPFParticles the input vector of PFParticles
         *
         *  @return the neutrino PFParticles
         */
        static PFParticleVector GetNeutrinos(const PFParticleVector &allPFParticles);

        /**
         *  @brief  Get a single neutrino from the input list
         *
         *  @param  allPFParticles the input vector of all PFParticles
         *
         *  @return the neutrino PFParticle
         */
        static art::Ptr<recob::PFParticle> GetNeutrino(const PFParticleVector &allPFParticles);

        /**
         *  @brief  Get the PFParticles that have been deemed the neutrino final states from the PFParticle hierarchy
         *
         *  @param  allPFParticles the input vector of all PFParticles
         *
         *  @return the neutrino final state PFParticles
         */
        static PFParticleVector GetNeutrinoFinalStates(const PFParticleVector &allPFParticles);

        /**
         *  @brief  Get the mapping from new PFParticles to old PFParticles via their common ID
         *
         *  @param  pfParticlesOld the input vector of old PFParticles (produced by Pandora)
         *  @param  pfParticlesNew the input vector of new PFParticles (with IDs that map to the old PFParticles)
         *
         *  @return the new to old mapping
         */
        static PFParticleToPFParticles GetNewToOldPFParticlesMap(const PFParticleVector &pfParticlesOld, const PFParticleVector &pfParticlesNew);

        /**
         *  @brief  Get the reconstructed neutrino interaction vertex position
         *
         *  @param  event the art event record
         *  @param  allPFParticles the input vector of all PFParticles
         *  @param  vertexLabel the label for the vertex producer module
         *
         *  @return the 3D position of the reconstructed interaction vertex
         */
        static TVector3 GetRecoNeutrinoVertex(const art::Event &event, const PFParticleVector &allPFParticles, const art::InputTag &vertexLabel);

        /**
         *  @brief  Get the parent of the input PFParticle
         *
         *  @param  particle the child PFParticle
         *  @param  pfParticleMap the PFParticle mapping 
         *
         *  @return the partent PFParticle
         */
        static art::Ptr<recob::PFParticle> GetParent(const art::Ptr<recob::PFParticle> &particle, const PFParticleMap &pfParticleMap);

        /**
         *  @brief  Get the daughters of the input PFParticle
         *
         *  @param  particle the parent PFParticle
         *  @param  pfParticleMap the PFParticle mapping
         *
         *  @return the daughter PFParticles
         */
        static PFParticleVector GetDaughters(const art::Ptr<recob::PFParticle> &particle, const PFParticleMap &pfParticleMap);

        /**
         *  @brief  Collect all particles downstream of the input PFParticle (including the particle itself)
         *
         *  @param  particle the input PFParticle
         *  @param  pfParticleMap the PFParticle mapping
         *
         *  @return the downstream PFParticles
         */
        static PFParticleVector GetDownstreamParticles(const art::Ptr<recob::PFParticle> &particle, const PFParticleMap &pfParticleMap);

        /**
         *  @brief  Collect all particles downstream of the input PFParticle (including the particle itself)
         *
         *  @param  particle the input PFParticle
         *  @param  pfParticleMap the PFParticle mapping
         *  @param  downstreamParticles the output downstream PFParticles
         */
        static void GetDownstreamParticles(const art::Ptr<recob::PFParticle> &particle, const PFParticleMap &pfParticleMap, PFParticleVector &downstreamParticles);

        /**
         *  @brief  Get the generation of a given PFParticle in the hierarchy
         *
         *  @param  particle the input PFParticle
         *  @param  pfParticleMap the PFParticle mapping
         *
         *  @return the generation
         */
        static unsigned int GetGeneration(const art::Ptr<recob::PFParticle> &particle, const PFParticleMap &pfParticleMap);

        /**
         *  @brief  Determine if a given slice has been selected as a neutrino by looking at the PdgCodes of its PFParticles
         *
         *  @param  slice the slice to check
         *  @param  slicesToPFParticles the mapping from slices to PFParticles
         *
         *  @return if the slice is selected as a neutrino
         */
        static bool IsSliceSelectedAsNu(const art::Ptr<recob::Slice> &slice, const SlicesToPFParticles &slicesToPFParticles);

        /**
         *  @brief  Determine if a given property exists in the input metadata object
         *
         *  @param  metadata the metadata object to check
         *  @param  property the property value to look up
         *
         *  @return if the proerty exists in the metadata
         */
        static bool HasMetadataValue(const art::Ptr<larpandoraobj::PFParticleMetadata> &metadata, const std::string &property);

        /**
         *  @brief  Get a metadata value with floating point type
         *
         *  @param  metadata the metadata
         *  @param  property the property name
         *
         *  @return the value
         */
        static float GetMetadataFloat(const art::Ptr<larpandoraobj::PFParticleMetadata> &metadata, const std::string &property);

        /**
         *  @brief  Get a metadata value converted to integer type
         *
         *  @param  metadata the metadata
         *  @param  property the property name
         *
         *  @return the value
         */
        static int GetMetadataInt(const art::Ptr<larpandoraobj::PFParticleMetadata> &metadata, const std::string &property);

        /**
         *  @brief  Get a metadata value converted to boolean type
         *
         *  @param  metadata the metadata
         *  @param  property the property name
         *
         *  @return the value
         */
        static bool GetMetadataBool(const art::Ptr<larpandoraobj::PFParticleMetadata> &metadata, const std::string &property);

        /**
         *  @brief  Get the track vs. shower score assigned by pandora from a metadata object
         *
         *  @param  metadata the metadata
         *
         *  @return the track score
         */
        static float GetTrackScore(const art::Ptr<larpandoraobj::PFParticleMetadata> &metadata);

        /**
         *  @brief  Get the neutrino vs. cosmic topological score from pandora from a metadata object
         *
         *  @param  metadata the metadata
         *
         *  @return the neutrino vs. cosmic topological score
         */
        static float GetTopologicalScore(const art::Ptr<larpandoraobj::PFParticleMetadata> &metadata);

        /**
         *  @brief  Count the number of hits in a given view
         *
         *  @param  hits the full vector of hits
         *  @param  view the view to use
         *
         *  @return the number of hits
         */
        static unsigned int CountHitsInView(const HitVector &hits, const geo::View_t &view);

        /**
         *  @brief  Get the space charge service for later use
         *
         *  @return the space charge service
         */
        static const spacecharge::SpaceChargeService::provider_type *const GetSpaceChargeService();

        /**
         *  @brief  Correct a given position for space charge (reco -> true)
         *
         *  @param  position the input position
         *  @param  pSpaceChargeService the space charge service
         *
         *  @return the space charge corrected 3D position
         */
        static TVector3 CorrectForSpaceCharge(const TVector3 &position, const spacecharge::SpaceChargeService::provider_type *const pSpaceChargeService);

        /**
         *  @brief  Get the score associated with a given PID criteria
         *
         *  @param  pid the input pid object
         *  @param  fCriteria a function that returns true if a given algorithm score passes some criteria
         *
         *  @return the score
         */
        static float GetPidScore(const art::Ptr<anab::ParticleID> &pid, const std::function<bool(const anab::sParticleIDAlgScores &)> &fCriteria);

        /**
         *  @brief  Convert a plane mask from a PID algorithm into a human readable view
         *
         *  @param  planeMask the input plane mask
         *
         *  @return the view
         */
        static geo::View_t GetView(const std::bitset<8> &planeMask);

        /**
         *  @brief  Get the bragg likelihood from a PID object
         *
         *  @param  pid the input pid
         *  @param  pdg the assumed pdg
         *  @param  view the view
         *
         *  @return the bragg likelihood
         */
        static float GetBraggLikelihood(const art::Ptr<anab::ParticleID> &pid, const int &pdg, const geo::View_t &view);
};

} // namespace ubcc1pi

#endif
