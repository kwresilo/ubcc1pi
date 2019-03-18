/**
 *  @file  ubcc1pi/Helpers/TruthHelper.cxx
 *
 *  @brief The implementation of the truth helper class
 */

#include "ubcc1pi/Helpers/TruthHelper.h"

#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"

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
        if (mcParticle->Process() == "primary")
            primaryMCParticles.push_back(mcParticle);
    }

    return primaryMCParticles;
}

} // namespace ubcc1pi
