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
        
bool RecoHelper::IsNeutrino(const art::Ptr<recob::PFParticle> &pfParticle)
{
    const auto pdg = std::abs(pfParticle->PdgCode());
    return (pdg == 12 ||  // Nue
            pdg == 14 ||  // Numu
            pdg == 16 );  // Nutau
}

// -----------------------------------------------------------------------------------------------------------------------------------------

PFParticleVector RecoHelper::GetNeutrinos(const PFParticleVector &allPFParticles)
{
    PFParticleVector neutrinos;

    for (const auto &particle : allPFParticles)
    {
        if (RecoHelper::IsNeutrino(particle))
            neutrinos.push_back(particle);
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

PFParticleToPFParticles RecoHelper::GetNewToOldPFParticlesMap(const PFParticleVector &pfParticlesOld, const PFParticleVector &pfParticlesNew)
{
    PFParticleToPFParticles outputMap;

    for (const auto &newPFParticle : pfParticlesNew)
    {
        outputMap.emplace(newPFParticle, PFParticleVector());
        for (const auto &oldPFParticle : pfParticlesOld)
        {
            if (newPFParticle->Self() == oldPFParticle->Self())
            {
                outputMap.at(newPFParticle).push_back(oldPFParticle);
            }
        }    
    }

    return outputMap;
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
        
unsigned int RecoHelper::GetGeneration(const art::Ptr<recob::PFParticle> &particle, const PFParticleMap &pfParticleMap)
{
    unsigned int generation = 0;

    art::Ptr<recob::PFParticle> nextParticle = particle;
    while (!nextParticle->IsPrimary())
    {
        nextParticle = RecoHelper::GetParent(nextParticle, pfParticleMap);
        generation++;
    }

    return generation;
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

float RecoHelper::GetTopologicalScore(const art::Ptr<larpandoraobj::PFParticleMetadata> &metadata)
{
    return RecoHelper::GetMetadataFloat(metadata, "NuScore");
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

// -----------------------------------------------------------------------------------------------------------------------------------------
            
float RecoHelper::GetPidScore(const art::Ptr<anab::ParticleID> &pid, const std::function<bool(const anab::sParticleIDAlgScores &)> &fCriteria)
{
    bool found = false;
    float score = -std::numeric_limits<float>::max();

    for (const auto &algo : pid->ParticleIDAlgScores())
    {
        if (!fCriteria(algo))
            continue;

        if (found)
            throw cet::exception("RecoHelper::GetPidScore") << " - Ambiguous criteria supplied." << std::endl;

        found = true;
        score = algo.fValue;
    }

    return score;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

geo::View_t RecoHelper::GetView(const std::bitset<8> &planeMask)
{
    // Here is a hack to get around a bug in the PID code. Some algorithms call W = 0, U = 1, V = 2. But others call W = 7, U = 6, V = 5
    const bool usesW = planeMask.test(0) || planeMask.test(7);
    const bool usesU = planeMask.test(1) || planeMask.test(6);
    const bool usesV = planeMask.test(2) || planeMask.test(5);
    
    if (usesW && !usesU && !usesV)
        return geo::kW;
    
    if (!usesW && usesU && !usesV)
        return geo::kU;
    
    if (!usesW && !usesU && usesV)
        return geo::kV;

    return geo::kUnknown;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

float RecoHelper::GetBraggLikelihood(const art::Ptr<anab::ParticleID> &pid, const int &pdg, const geo::View_t &view)
{
    return RecoHelper::GetPidScore(pid, [&](const anab::sParticleIDAlgScores &algo) -> bool {
        return (algo.fAlgName == "BraggPeakLLH"              &&
                algo.fTrackDir == anab::kForward             &&
                algo.fAssumedPdg == pdg                      &&
                RecoHelper::GetView(algo.fPlaneMask) == view );
    });
}

} // namespace ubcc1pi
