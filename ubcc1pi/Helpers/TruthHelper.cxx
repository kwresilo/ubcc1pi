/**
 *  @file  ubcc1pi/Helpers/TruthHelper.cxx
 *
 *  @brief The implementation of the truth helper class
 */

#include "ubcc1pi/Helpers/TruthHelper.h"
#include "ubcc1pi/Helpers/AnalysisHelper.h"

#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"

namespace ubcc1pi
{

art::Ptr<simb::MCTruth> TruthHelper::GetNeutrinoMCTruth(const art::Event &event, const art::InputTag &mcTruthLabel)
{
    art::Ptr<simb::MCTruth> selectedTruth;
    bool foundNeutrino = false;
    float maxEnergy = -std::numeric_limits<float>::max();

    for (const auto &mcTruth : CollectionHelper::GetCollection<simb::MCTruth>(event, mcTruthLabel))
    {
        // We only care about beam neutrinos
        if (mcTruth->Origin() != simb::kBeamNeutrino)
            continue;

        // We only care about the most energetic neutrino
        const auto energy = mcTruth->GetNeutrino().Nu().E();
        if (energy < maxEnergy)
            continue;

        // Save the details of this MCTruth
        foundNeutrino = true;
        maxEnergy = energy;
        selectedTruth = mcTruth;
    }

    if (!foundNeutrino)
        throw cet::exception("TruthHelper::GetNeutrinoMCTruth") << " - No beam neutrino MCTruth objects found." << std::endl;

    return selectedTruth;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

MCParticleVector TruthHelper::GetMCParticlesFromMCTruth(const art::Event &event, const art::InputTag &mcTruthLabel, const art::InputTag &mcParticleLabel, const art::Ptr<simb::MCTruth> &mcTruth)
{
    MCParticleVector outputMCParticles;

    const auto mcTruthToMCParticle = CollectionHelper::GetAssociation<simb::MCTruth, simb::MCParticle>(event, mcTruthLabel, mcParticleLabel);
    const auto mcParticleToMCTruth = CollectionHelper::GetReversedAssociation(mcTruthToMCParticle);

    for (const auto &mcParticle : CollectionHelper::GetCollection<simb::MCParticle>(event, mcParticleLabel))
    {
        const auto associatedMCTruth = CollectionHelper::GetSingleAssociated(mcParticle, mcParticleToMCTruth);

        // Insist that the MCParticle originated from the desired MCTruth
        if (associatedMCTruth != mcTruth)
            continue;

        outputMCParticles.push_back(mcParticle);
    }

    return outputMCParticles;
}

// -----------------------------------------------------------------------------------------------------------------------------------------
        
MCParticleVector TruthHelper::GetPrimaryMCParticles(const MCParticleVector &allMCParticles)
{
    MCParticleVector primaryMCParticles;

    for (const auto &mcParticle : allMCParticles)
    {
        if (TruthHelper::IsPrimary(mcParticle))
            primaryMCParticles.push_back(mcParticle);
    }

    return primaryMCParticles;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

bool TruthHelper::IsPrimary(const art::Ptr<simb::MCParticle> &particle)
{
    return (particle->Process() == "primary");    
}

// -----------------------------------------------------------------------------------------------------------------------------------------

MCParticleMap TruthHelper::GetMCParticleMap(const MCParticleVector &allMCParticles)
{
    MCParticleMap mcParticleMap;

    for (const auto &mcParticle : allMCParticles)
    {
        if (!mcParticleMap.emplace(mcParticle->TrackId(), mcParticle).second)
            throw cet::exception("TruthHelper::GetMCParticleMap") << " - Found repeated MCParticle with TrackId = " << mcParticle->TrackId() << "." << std::endl;
    }

    return mcParticleMap;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

art::Ptr<simb::MCParticle> TruthHelper::GetMother(const art::Ptr<simb::MCParticle> &particle, const MCParticleMap &mcParticleMap)
{
    const auto motherIter = mcParticleMap.find(particle->Mother());
    if (motherIter == mcParticleMap.end())
        throw cet::exception("TruthHelper::GetMother") << " - Couldn't find mother MCParticle in hierarchy. Are you trying to get the mother of a primary?" << std::endl;

    return motherIter->second;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

MCParticleVector TruthHelper::GetDaughters(const art::Ptr<simb::MCParticle> &particle, const MCParticleMap &mcParticleMap)
{
    MCParticleVector daughters;

    for (int i = 0; i < particle->NumberDaughters(); ++i)
    {
        const auto daughterIter = mcParticleMap.find(particle->Daughter(i));
        if (daughterIter == mcParticleMap.end())
        {
            // ATTN this can happen - MCParticles outside of the cryostat are dropped 
            std::cout << "WARNING - Can't find daughter MCParticle in hierarchy, likely out of cryostat - skipping" << std::endl;
            continue;
        }

        daughters.push_back(daughterIter->second);
    }

    return daughters;
}

// -----------------------------------------------------------------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------------------------------------------------------------

TruthHelper::Interaction::Interaction(const art::Event &event, const art::InputTag &mcTruthLabel, const art::InputTag &mcParticleLabel)
{
    const auto nuMCTruth = TruthHelper::GetNeutrinoMCTruth(event, mcTruthLabel);

    if (!nuMCTruth->NeutrinoSet())
        throw cet::exception("Interaction::Interaction") << " - Beam neutrino MCTruth block doesn't have it's neutrino information filled." << std::endl;

    const auto nu = nuMCTruth->GetNeutrino();
    m_neutrino = nu.Nu();
    m_ccnc = static_cast<simb::curr_type_>(nu.CCNC());
    m_mode = static_cast<simb::int_type_>(nu.Mode());
    
    m_allMCParticles = TruthHelper::GetMCParticlesFromMCTruth(event, mcTruthLabel, mcParticleLabel, nuMCTruth);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

simb::MCParticle TruthHelper::Interaction::GetNeutrino() const
{
    return m_neutrino;
}

// -----------------------------------------------------------------------------------------------------------------------------------------
        
simb::curr_type_ TruthHelper::Interaction::GetCCNC() const
{
    return m_ccnc;
}

// -----------------------------------------------------------------------------------------------------------------------------------------
        
simb::int_type_ TruthHelper::Interaction::GetInteractionMode() const
{
    return m_mode;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

MCParticleVector TruthHelper::Interaction::GetAllMCParticles() const
{
    return m_allMCParticles;
}

} // namespace ubcc1pi
