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

PFParticleVector RecoHelper::GetNeutrinos(const PFParticleVector &allPFParticles)
{
    PFParticleVector neutrinos;

    for (const auto &particle : allPFParticles)
    {
        const auto pdg = particle->PdgCode();
        if (pdg == 12 || // Nue
            pdg == 14 || // Numu
            pdg == 16 )  // Nutau
        {
            neutrinos.push_back(particle);
        }
    }
            
    return neutrinos;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

art::Ptr<recob::PFParticle> RecoHelper::GetNeutrino(const PFParticleVector &allPFParticles)
{
    const auto &neutrinos = RecoHelper::GetNeutrinos(allPFParticles);

    if (neutrinos.empty())
        throw cet::exception("RecoHelper::GetNeutrino") << " - Didn't find a neutrino PFParticle." << std::endl;

    if (neutrinos.size() > 1)
        throw cet::exception("RecoHelper::GetNeutrino") << " - Found multiple neutrino PFParticles." << std::endl;

    return neutrinos.front();
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

TVector3 RecoHelper::GetRecoNeutrinoVertex(const art::Event &event, const PFParticleVector &allPFParticles, const art::InputTag &vertexLabel)
{
    const auto pfpToVertex = CollectionHelper::GetAssociation<recob::PFParticle, recob::Vertex>(event, vertexLabel);

    try
    {
        const auto neutrino = RecoHelper::GetNeutrino(allPFParticles);
        const auto vertex = CollectionHelper::GetSingleAssociated(neutrino, pfpToVertex);

        return TVector3(vertex->position().X(), vertex->position().Y(), vertex->position().Z());
    }
    catch (const cet::exception &)
    {
        return TVector3(-std::numeric_limits<float>::max(), -std::numeric_limits<float>::max(), -std::numeric_limits<float>::max());
    }
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

// -----------------------------------------------------------------------------------------------------------------------------------------
        
bool RecoHelper::IsSliceSelectedAsNu(const art::Ptr<recob::Slice> &slice, const SlicesToPFParticles &slicesToPFParticles)
{
    const auto neutrinos = RecoHelper::GetNeutrinos(CollectionHelper::GetManyAssociated(slice, slicesToPFParticles));
    
    if (neutrinos.size() > 1)
        throw cet::exception("RecoHelper::IsSliceSelectedAsNu") << " - Found multiple neutrino PFParticles associated to slice." << std::endl;
 
    return (!neutrinos.empty());
}

// -----------------------------------------------------------------------------------------------------------------------------------------
        
bool RecoHelper::HasMetadataValue(const art::Ptr<larpandoraobj::PFParticleMetadata> &metadata, const std::string &property)
{
    const auto &properties = metadata->GetPropertiesMap();
    return (properties.find(property) != properties.end());
}

// -----------------------------------------------------------------------------------------------------------------------------------------

float RecoHelper::GetMetadataFloat(const art::Ptr<larpandoraobj::PFParticleMetadata> &metadata, const std::string &property)
{
    if (!RecoHelper::HasMetadataValue(metadata, property))
        throw cet::exception("RecoHelper::GetMetadataValue") << " - Can't find metadata with key: " << property << std::endl;

    return metadata->GetPropertiesMap().at(property);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

int RecoHelper::GetMetadataInt(const art::Ptr<larpandoraobj::PFParticleMetadata> &metadata, const std::string &property)
{
    return static_cast<int>(std::round(RecoHelper::GetMetadataFloat(metadata, property)));
}

// -----------------------------------------------------------------------------------------------------------------------------------------

bool RecoHelper::GetMetadataBool(const art::Ptr<larpandoraobj::PFParticleMetadata> &metadata, const std::string &property)
{
    switch (RecoHelper::GetMetadataInt(metadata, property))
    {
        case 0:
            return false;
        case 1:
            return true;
        default:
            throw cet::exception("RecoHelper::GetMetadataValue") << " - Can't interpret metadata: " << property << " as boolean." << std::endl;
    }
}

// -----------------------------------------------------------------------------------------------------------------------------------------

float RecoHelper::GetTrackScore(const art::Ptr<larpandoraobj::PFParticleMetadata> &metadata)
{
    return RecoHelper::GetMetadataFloat(metadata, "TrackScore");
}

// -----------------------------------------------------------------------------------------------------------------------------------------

unsigned int RecoHelper::CountHitsInView(const HitVector &hits, const geo::View_t &view)
{
    unsigned int count = 0;
    for (const auto &hit : hits)
    {
        if (hit->View() == view)
            count++;
    }

    return count;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

const spacecharge::SpaceChargeService::provider_type *const RecoHelper::GetSpaceChargeService()
{
    return lar::providerFrom<spacecharge::SpaceChargeService>();
}

// -----------------------------------------------------------------------------------------------------------------------------------------

TVector3 RecoHelper::CorrectForSpaceCharge(const TVector3 &position, const spacecharge::SpaceChargeService::provider_type *const pSpaceChargeService)
{
    const auto offset = pSpaceChargeService->GetCalPosOffsets(geo::Point_t(position.X(), position.Y(), position.Z()));
    return TVector3(position.X() - offset.X(), position.Y() + offset.Y(), position.Z() + offset.Z());
}

} // namespace ubcc1pi
