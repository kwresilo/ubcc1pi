/**
 *  @file  ubcc1pi/Producers/ParticleID_module.cc
 *
 *  @brief The implementation file for the particle id producer.
 */

#include "ubcc1pi/Producers/ParticleID.h"

#include "ubcc1pi/Helpers/CollectionHelper.h"
#include "ubcc1pi/Helpers/RecoHelper.h"

namespace ubcc1pi
{

ParticleID::ParticleID(const art::EDProducer::Table<Config> &config) :
    art::EDProducer(config),
    m_config(config)
{
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void ParticleID::produce(art::Event &event)
{
    const auto pfParticleLabel = m_config().PFParticleLabel();
    const auto hitLabel = m_config().HitLabel();

    const auto allPFParticles = CollectionHelper::GetCollection<recob::PFParticle>(event, pfParticleLabel);
    const auto finalStates = RecoHelper::GetNeutrinoFinalStates(allPFParticles);
    const auto pfpToHits = CollectionHelper::GetAssociationViaCollection<recob::PFParticle, recob::Cluster, recob::Hit>(event, pfParticleLabel, pfParticleLabel, pfParticleLabel);

    for (const auto &pfParticle : finalStates)
    {
        std::cout << "PFParticle " << pfParticle->Self() << std::endl;

        const auto hits = CollectionHelper::GetManyAssociated(pfParticle, pfpToHits);
        std::cout << "  - nHits : " << hits.size() << std::endl;
    }
}

} // namespace ubcc1pi
