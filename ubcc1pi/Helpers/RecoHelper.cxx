/**
 *  @file  ubcc1pi/Helpers/RecoHelper.cxx
 *
 *  @brief The implementation file for the reco helper class
 */

#include "ubcc1pi/Helpers/RecoHelper.h"

namespace ubcc1pi
{

PFParticleMap RecoHelper::GetPFParticleMap(const PFParticleVector &allPFParticles)
{
    PFParticleMap pfParticleMap;

    for (const auto &pfParticle : allPFParticles)
    {
        if (!pfParticleMap.emplace(pfParticle->Self(), pfParticle).second)
            throw cet::exception("RecoHelper::GetPFParticleMap") << " - Found repeated PFParticle with Self = " << pfParticle->Self() << "." << std::endl;
    }

    return pfParticleMap;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

art::Ptr<recob::PFParticle> RecoHelper::GetNeutrino(const PFParticleVector &allPFParticles)
{
    bool foundNu = false;
    art::Ptr<recob::PFParticle> neutrino;

    for (const auto &particle : allPFParticles)
    {
        const auto pdg = particle->PdgCode();
        if (pdg == 12 || // Nue
            pdg == 14 || // Numu
            pdg == 16 )  // Nutau
        {
            if (foundNu)
                throw cet::exception("RecoHelper::GetNeutrino") << " - Found multiple neutrino PFParticles." << std::endl;

            foundNu = true;
            neutrino = particle;
        }
    }
            
    if (!foundNu)
        throw cet::exception("RecoHelper::GetNeutrino") << " - Didn't find a neutrino PFParticle." << std::endl;

    return neutrino;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

PFParticleVector RecoHelper::GetNeutrinoFinalStates(const PFParticleVector &allPFParticles)
{
    // Get the neutrino PFParticle - if there is one!
    art::Ptr<recob::PFParticle> neutrino;
    try
    {
        neutrino = RecoHelper::GetNeutrino(allPFParticles);
    }
    catch (const cet::exception &)
    {
        // No neutrino found
        return PFParticleVector();
    }
    
    const auto pfParticleMap = RecoHelper::GetPFParticleMap(allPFParticles);
    return RecoHelper::GetDaughters(neutrino, pfParticleMap);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

art::Ptr<recob::PFParticle> RecoHelper::GetParent(const art::Ptr<recob::PFParticle> &particle, const PFParticleMap &pfParticleMap)
{
    if (particle->IsPrimary())
        throw cet::exception("RecoHelper::GetParent") << " - PFParticle is primary, so doesn't have a parent." << std::endl;

    const auto parentIter = pfParticleMap.find(particle->Parent());
    if (parentIter == pfParticleMap.end())
        throw cet::exception("RecoHelper::GetParent") << " - Couldn't find parent PFParticle in hierarchy." << std::endl;

    return parentIter->second;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

PFParticleVector RecoHelper::GetDaughters(const art::Ptr<recob::PFParticle> &particle, const PFParticleMap &pfParticleMap)
{
    PFParticleVector daughters;

    for (int i = 0; i < particle->NumDaughters(); ++i)
    {
        const auto daughterIter = pfParticleMap.find(particle->Daughter(i));
        if (daughterIter == pfParticleMap.end())
            throw cet::exception("RecoHelper::GetDaughter") << " - Couldn't find daughter PFParticle in hierarchy." << std::endl;

        daughters.push_back(daughterIter->second);
    }

    return daughters;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

PFParticleVector RecoHelper::GetDownstreamParticles(const art::Ptr<recob::PFParticle> &particle, const PFParticleMap &pfParticleMap)
{
    PFParticleVector downstreamParticles;
    RecoHelper::GetDownstreamParticles(particle, pfParticleMap, downstreamParticles);

    return downstreamParticles;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void RecoHelper::GetDownstreamParticles(const art::Ptr<recob::PFParticle> &particle, const PFParticleMap &pfParticleMap, PFParticleVector &downstreamParticles)
{
    downstreamParticles.push_back(particle);

    for (const auto &daughter : RecoHelper::GetDaughters(particle, pfParticleMap))
        RecoHelper::GetDownstreamParticles(daughter, pfParticleMap, downstreamParticles);
}

} // namespace ubcc1pi
