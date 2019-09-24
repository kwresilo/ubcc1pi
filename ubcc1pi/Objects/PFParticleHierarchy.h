/**
 *  @file  ubcc1pi/Objects/PFParticleHierarchy.h
 *
 *  @brief The header file for the pfparticle heirarchy object
 */

#ifndef UBCC1PI_OBJECTS_PFPARTICLE_HIERARCHY
#define UBCC1PI_OBJECTS_PFPARTICLE_HIERARCHY

#include "ubcc1pi/Helpers/CollectionHelper.h"

namespace ubcc1pi
{

class PFParticleHierarchy
{
    public:
        /**
         *  @brief  Constructor
         *
         *  @param  event the art event
         *  @param  pfParticles the full list of neutrino induced PFParticles to arrange into a hierarchy
         *  @param  pfParticlesToTracks the mapping from PFParticles to tracks
         *  @param  vertex the reconstructed vertex position
         */
        PFParticleHierarchy(const PFParticleVector &pfParticles, const PFParticleToTracks &pfParticleToTracks, const TVector3 &vertex);

        /**
         *  @brief  Get the final state particles, daughters of the neutrino
         *
         *  @return vector of final state particles
         */
        PFParticleVector GetFinalStates() const;

        /**
         *  @brief  Check if a given particle is a final state of the neutrino
         *
         *  @param  pfParticle the PFParticle to check
         *
         *  @return if the particle is a final state
         */
        bool IsFinalState(const art::Ptr<recob::PFParticle> &pfParticle) const;

        /**
         *  @brief  Get the parent of the input PFParticle (NB. the neutrino has no parent)
         *
         *  @param  pfParticle the PFParticle to get the parent of
         *
         *  @return the parent PFParticle
         */
        art::Ptr<recob::PFParticle> GetParent(const art::Ptr<recob::PFParticle> &pfParticle) const;

        /**
         *  @brief  Get the daughters of the input PFParticle
         *
         *  @param  pfParticle the PFParticle to get the daughters of
         *
         *  @return the daughters
         */
        PFParticleVector GetDaughters(const art::Ptr<recob::PFParticle> &pfParticle) const;

    private:
        TVector3 GetUpstreamEndpoint(const art::Ptr<recob::PFParticle> &pfParticle) const;
        TVector3 GetDownstreamEndpoint(const art::Ptr<recob::PFParticle> &pfParticle) const;

        PFParticleToPFParticles  m_pfParticleToParentsMap; ///< The mapping from PFParticles to their parent 
};

} // namespace ubc1pi

