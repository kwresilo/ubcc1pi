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
         *  @brief  Get the degrees of freedom associated with a given PID criteria
         *
         *  @param  pid the input pid object
         *  @param  fCriteria a function that returns true if a given algorithm score passes some criteria
         *
         *  @return the degrees of freedom
         */
        static float GetPidDegreesOfFreedom(const art::Ptr<anab::ParticleID> &pid, const std::function<bool(const anab::sParticleIDAlgScores &)> &fCriteria);

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
         *  @param  dir if the fit is forward or backward
         *
         *  @return the bragg likelihood
         */
        static float GetBraggLikelihood(const art::Ptr<anab::ParticleID> &pid, const int &pdg, const geo::View_t &view, const anab::kTrackDir &dir = anab::kForward);

        /**
         *  @brief  Get the bragg likelihood from a PID object using collection plane and falling back on induction planes if not available
         *
         *  @param  pid the input pid
         *  @param  pdg the assumed pdg
         *  @param  dir if the fit is forward or backward
         *  @param  yzAngle the angle of the input track in the YZ plane
         *  @param  sin2AngleThreshold the threshold within which the track musn't points along the W wire direction to use the info on that plane
         *
         *  @return the bragg likelihood over all planes
         */
        static float GetBraggLikelihood(const art::Ptr<anab::ParticleID> &pid, const int &pdg, const anab::kTrackDir &dir, const float yzAngle, const float sin2AngleThreshold);

        /**
         *  @brief  Get the number of degrees of freedom associated with a given bragg likelihood
         *
         *  @param  pid the input pid
         *  @param  pdg the assumed pdg
         *  @param  view the view
         *  @param  dir if the fit is forward or backward
         *
         *  @return the number of degrees of freedom
         */
        static float GetBraggLikelihoodDegreesOfFreedom(const art::Ptr<anab::ParticleID> &pid, const int &pdg, const geo::View_t &view, const anab::kTrackDir &dir);

        /**
         *  @brief  Get the indices of the trajectory points that are valid in an input track
         *
         *  @param  track the input track
         *
         *  @return the valid trajectory point indices
         */
        static std::vector<size_t> GetValidPoints(const art::Ptr<recob::Track> &track);

        /**
         *  @brief  Get the range of a track by integrating the distances between trajectory points
         *
         *  @param  track the input track
         *  @param  pSpaceChargeService the input space charge service
         *
         *  @return the range of the track
         */
        static float GetRange(const art::Ptr<recob::Track> &track, const spacecharge::SpaceChargeService::provider_type *const pSpaceChargeService);

        /**
         *  @brief  Get the distance from the start of a track to a given point in longitudinal and transverse coordinates based on the track direction
         *
         *  @param  track the track
         *  @param  point the point
         *  @param  pSpaceChargeService the space charge service
         *  @param  transverseDist the output transverse distance
         *  @param  longitudinalDist the output longitudinal distance
         */
        static void GetDistanceToPoint(const art::Ptr<recob::Track> &track, const TVector3 &point, const spacecharge::SpaceChargeService::provider_type *const pSpaceChargeService, float &transverseDist, float &longitudinalDist);

        /**
         *  @brief  Get the standard deviation of the angular differences between sequential track segments
         *
         *  @param  track the input track
         *
         *  @return the wiggliness
         */
        static float GetTrackWiggliness(const art::Ptr<recob::Track> &track);

        /**
         *  @brief  Count the input spacepoints which are within a given distance of the input track end
         *
         *  @param  track the input track
         *  @param  spacePoints the input spacepoints
         *  @param  distance the threshold distance
         *  @param  pSpaceChargeService the input space charge service
         *
         *  @return number of spacepoints near track end
         */
        static unsigned int CountSpacePointsNearTrackEnd(const art::Ptr<recob::Track> &track, const SpacePointVector &spacePoints, const float distance, const spacecharge::SpaceChargeService::provider_type *const pSpaceChargeService);

        /**
         *  @brief  Get the mean dEdx at the start of the track excluding any hits greater than one standard deviation from the median dEdx
         *
         *  @param  dedxPerHit the input dEdx per hit
         *  @param  residualRangePerHit the input residual ranges per hit
         *  @param  nHitsToSkip the number of hits to skip at the start of the track
         *  @param  lengthFraction the fraction of the total range to use to isolate the start of the track
         *
         *  @return the truncated mean dEdx at the start of the track
         */
        static float GetTruncatedMeandEdxAtTrackStart(const std::vector<float> dedxPerHit, const std::vector<float> &residualRangePerHit, const unsigned int nHitsToSkip, const float lengthFraction);
};

} // namespace ubcc1pi

#endif
