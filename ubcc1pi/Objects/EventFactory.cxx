/**
 *  @file  ubcc1pi/Objects/EventFactory.cxx
 *
 *  @brief The implementation of the event factory class
 */

#include "ubcc1pi/Objects/EventFactory.h"

#include "ubcc1pi/Helpers/TruthHelper.h"
#include "ubcc1pi/Helpers/AnalysisHelper.h"
#include "ubcc1pi/Helpers/DebugHelper.h"

namespace ubcc1pi
{

void EventFactory::PopulateEvent(const art::Event &event, const Config &config, Event *pOutputEvent)
{
    pOutputEvent->Reset();

    // Populate the metadata
    auto &metadata = pOutputEvent->metadata;
    EventFactory::PopulateEventMetadata(event, config, metadata);

    // Populate the truth info
    if (config.HasTruthInfo())
    {
        auto &truth = pOutputEvent->truth;
        EventFactory::PopulateEventTruthInfo(event, config, truth);
    }

    // Populate the reco info
    //auto &reco = pOutputEvent->reco;
    //EventFactory::PopulateEventRecoInfo(event, reco);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void EventFactory::PopulateEventMetadata(const art::Event &event, const Config &config, Event::Metadata &metadata)
{
    metadata.run.Set(event.run());
    metadata.subRun.Set(event.subRun());
    metadata.event.Set(event.event());
    metadata.hasTruthInfo.Set(config.HasTruthInfo());
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void EventFactory::PopulateEventTruthInfo(const art::Event &event, const Config &config, Event::Truth &truth)
{
    const TruthHelper::Interaction interaction(event, config.MCTruthLabel(), config.MCParticleLabel());

    // Fill the interaction type info
    truth.isCC.Set(interaction.GetCCNC() == simb::kCC);
    truth.interactionMode.Set(interaction.GetInteractionMode());
    truth.interactionString.Set(DebugHelper::GetInteractionString(interaction, true));

    // Fill the neutrino info
    const auto neutrino = interaction.GetNeutrino();
    truth.nuPdgCode.Set(neutrino.PdgCode());
    truth.nuEnergy.Set(neutrino.E());
    truth.nuVertex.Set(TVector3(neutrino.Vx(), neutrino.Vy(), neutrino.Vz()));

    // Fill the particle info
    // Here we only persist the final state particles that are visible (could possibly, but need not actually, produce hits)
    const auto visibleFinalStates = AnalysisHelper::GetVisibleFinalStates(interaction);
    for (const auto &mcParticle : visibleFinalStates)
    {
        Event::Truth::Particle particle;
        EventFactory::PopulateEventTruthParticleInfo(mcParticle, particle);
        truth.particles.push_back(particle);
    }
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void EventFactory::PopulateEventTruthParticleInfo(const art::Ptr<simb::MCParticle> &mcParticle, Event::Truth::Particle &particle)
{
    particle.pdgCode.Set(mcParticle->PdgCode());
    particle.momentumX.Set(mcParticle->Px());
    particle.momentumY.Set(mcParticle->Py());
    particle.momentumZ.Set(mcParticle->Pz());
    particle.momentum.Set(mcParticle->P());
    particle.energy.Set(mcParticle->E());
    particle.mass.Set(mcParticle->Mass());
}

} // namespace ubcc1pi
