/**
 *  @file  ubcc1pi/Helpers/TruthHelper.cxx
 *
 *  @brief The implementation of the truth helper class
 */

#include "ubcc1pi/Helpers/TruthHelper.h"

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
        {
            std::cout<<"DEBUG TruthHelper::GetMCParticlesFromMCTruth associatedMCTruth != mcTruth"<<std::endl;
            continue;
        }

        outputMCParticles.push_back(mcParticle);
    }

    std::cout<<"DEBUG TruthHelper::GetMCParticlesFromMCTruth outputMCParticles.size() "<<outputMCParticles.size()<<std::endl;
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
    return (particle->Process() == "primary" && particle->StatusCode() == 1);
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
            // ATTN this can happen - MCParticles outside of the cryostat are dropped and the hierarchy gets truncated
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

// -----------------------------------------------------------------------------------------------------------------------------------------

MCParticleVector TruthHelper::GetDownstreamParticles(const art::Ptr<simb::MCParticle> &particle, const MCParticleMap &mcParticleMap)
{
    MCParticleVector downstreamParticles;
    TruthHelper::GetDownstreamParticles(particle, mcParticleMap, downstreamParticles);

    return downstreamParticles;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void TruthHelper::GetDownstreamParticles(const art::Ptr<simb::MCParticle> &particle, const MCParticleMap &mcParticleMap, MCParticleVector &downstreamParticles)
{
    downstreamParticles.push_back(particle);

    for (const auto &daughter : TruthHelper::GetDaughters(particle, mcParticleMap))
        TruthHelper::GetDownstreamParticles(daughter, mcParticleMap, downstreamParticles);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void TruthHelper::FollowScatters(const art::Ptr<simb::MCParticle> &particle, const MCParticleMap &mcParticleMap, ScatterVector &scatters, art::Ptr<simb::MCParticle> &scatteredParticle)
{
    scatteredParticle = particle;

    // First look for elastic scatters
    for (const auto &daughter : TruthHelper::GetDaughters(particle, mcParticleMap))
    {
        if (daughter->Process() == "hadElastic")
        {
            try
            {
                scatters.push_back(TruthHelper::GetElasticScatter(particle, daughter));
            }
            catch (const cet::exception &ex)
            {
                std::cout << "TruthHelper::FollowScatters. WARNING - " << ex.explain_self() << std::endl;
            }
        }
    }

    // Now follow inelastic scatters
    if (particle->EndProcess().find("Inelastic") == std::string::npos)
        return;

    MCParticleVector inelasticProducts;
    art::Ptr<simb::MCParticle> finalStateParticle;
    bool foundInelasticScatter = false;

    for (const auto &daughter : TruthHelper::GetDaughters(particle, mcParticleMap))
    {
        if (daughter->Process().find("Inelastic") != std::string::npos)
        {
            // If the daughter has the same PDG as the incident particle, then assume it's the same particle before & after the scatter
            if (daughter->PdgCode() != particle->PdgCode())
            {
                inelasticProducts.push_back(daughter);
                continue;
            }

            // If there are multiple particles from the inelastic collision, then don't treat it as a scatter
            if (foundInelasticScatter)
            {
                foundInelasticScatter = false;
                break;
            }

            finalStateParticle = daughter;
            foundInelasticScatter = true;
        }
    }

    if (!foundInelasticScatter)
        return;

    // Get the last momentum point before the scatter
    TVector3 initialMomentum(-std::numeric_limits<float>::max(), -std::numeric_limits<float>::max(), -std::numeric_limits<float>::max());
    for (unsigned int i = 0; i < particle->NumberTrajectoryPoints(); ++i)
    {
        if (particle->T(i) >= finalStateParticle->T())
            break;

        initialMomentum = particle->Momentum(i).Vect();
    }

    scatters.emplace_back(false, initialMomentum, finalStateParticle->Momentum().Vect(), finalStateParticle, inelasticProducts);
    TruthHelper::FollowScatters(finalStateParticle, mcParticleMap, scatters, scatteredParticle);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

TruthHelper::Scatter TruthHelper::GetElasticScatter(const art::Ptr<simb::MCParticle> &incident, const art::Ptr<simb::MCParticle> &target)
{
    if (target->Process() != "hadElastic")
        throw cet::exception("TruthHelper::GetElasticScatter") << " - Input target particle doesn't look to be from an elastic scatter" << std::endl;

    if (target->Mother() != incident->TrackId())
        throw cet::exception("TruthHelper::GetElasticScatter") << " - Input target particle isn't a daughter of the input incident particle" << std::endl;

    const unsigned int nTrajectoryPoints = incident->NumberTrajectoryPoints();
    if (nTrajectoryPoints < 2)
        throw cet::exception("TruthHelper::GetElasticScatter") << " - Incident particle has less than 2 trajectory points!" << std::endl;

    // Find the trajectory points just before and just after the scatter
    const float scatterTime(target->T());
    unsigned int postScatterIndex = std::numeric_limits<unsigned int>::max();
    bool foundScatterPoint = false;

    for (unsigned int i = 1; i < incident->NumberTrajectoryPoints(); ++i)
    {
        if (incident->T(i - 1) <= scatterTime && incident->T(i) >= scatterTime && incident->T(i - 1) < incident->T(i))
        {
            foundScatterPoint = true;
            postScatterIndex = i;
        }
    }
    unsigned int preScatterIndex = postScatterIndex - 1;

    if (!foundScatterPoint)
        throw cet::exception("TruthHelper::GetElasticScatter") << " - Couldn't find trajectory points before and after the scatter" << std::endl;

    const auto initialMomentum(incident->Momentum(preScatterIndex).Vect());
    const auto finalMomentum(incident->Momentum(postScatterIndex).Vect());

    return TruthHelper::Scatter(true, initialMomentum, finalMomentum, incident, { target });
}

// -----------------------------------------------------------------------------------------------------------------------------------------

TruthHelper::EndState TruthHelper::GetEndState(const art::Ptr<simb::MCParticle> &particle, const MCParticleMap &mcParticleMap)
{
    EndState::Type type = EndState::Type::Other;
    const auto finalParticleMomentum = particle->Momentum(std::max(static_cast<unsigned int>(0), particle->NumberTrajectoryPoints() - 2)).Vect();
    MCParticleVector products;

    // ATTN for now, just return OTHER for all non-pions, we only use the end-state of the pions
    if (particle->PdgCode() != 211)
        return TruthHelper::EndState(type, finalParticleMomentum, products);

    // Follow the scatters to make sure we are really looking at the correct MCParticle
    ScatterVector scatters;
    art::Ptr<simb::MCParticle> scatteredParticle;
    TruthHelper::FollowScatters(particle, mcParticleMap, scatters, scatteredParticle);
    if (particle != scatteredParticle)
        throw cet::exception("TruthHelper::GetEndState") << " - Input MCParticle undergos a scatter before reaching it's end-state" << std::endl;

    bool hasPi0 = false;
    bool hasDecayMuon = false;
    bool hasDecayMuonNeutrino = false;

    for (const auto &daughter : TruthHelper::GetDaughters(particle, mcParticleMap))
    {
        // Ignore ionisation electrons
        if (daughter->PdgCode() == 11 && daughter->Process() == "hIoni")
            continue;

        // Ignore nuclei from elastic scatters
        if (daughter->Process() == "hadElastic")
            continue;

        // Treat everything else as an end-state interaction product
        products.push_back(daughter);

        if (daughter->PdgCode() == 111 && daughter->Process() == "pi+Inelastic")
            hasPi0 = true;

        if (daughter->PdgCode() == -13 && daughter->Process() == "Decay")
            hasDecayMuon = true;

        if (daughter->PdgCode() == 14 && daughter->Process() == "Decay")
            hasDecayMuonNeutrino = true;
    }

    // Work out the end-state type
    if (products.empty())
    {
        type = EndState::Type::None;
    }
    else if (hasDecayMuon && hasDecayMuonNeutrino && products.size() == 2)
    {
        type = EndState::Type::DecayToMuon;
    }
    else if (particle->EndProcess() == "pi+Inelastic")
    {
        type = hasPi0 ? EndState::Type::Pi0ChargeExchange : EndState::Type::InelasticAbsorption;
    }

    return TruthHelper::EndState(type, finalParticleMomentum, products);
}

// -----------------------------------------------------------------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------------------------------------------------------------

TruthHelper::Scatter::Scatter(const bool isElastic, const TVector3 initialMomentum, const TVector3 finalMomentum, const art::Ptr<simb::MCParticle> finalParticle, const MCParticleVector products) :
    m_isElastic(isElastic),
    m_initialMomentum(initialMomentum),
    m_finalMomentum(finalMomentum),
    m_finalParticle(finalParticle),
    m_products(products)
{
}

// -----------------------------------------------------------------------------------------------------------------------------------------

float TruthHelper::Scatter::GetScatteringCosTheta() const
{
    if (m_initialMomentum.Mag() * m_finalMomentum.Mag() < std::numeric_limits<float>::epsilon())
        return -std::numeric_limits<float>::max();

    return m_initialMomentum.Dot(m_finalMomentum) / (m_initialMomentum.Mag() * m_finalMomentum.Mag());
}

// -----------------------------------------------------------------------------------------------------------------------------------------

float TruthHelper::Scatter::GetMomentumFractionLost() const
{
    if (m_initialMomentum.Mag() * m_finalMomentum.Mag() < std::numeric_limits<float>::epsilon())
        return -std::numeric_limits<float>::max();

    return (m_initialMomentum.Mag() - m_finalMomentum.Mag()) / m_initialMomentum.Mag();
}

// -----------------------------------------------------------------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------------------------------------------------------------

TruthHelper::EndState::EndState(const Type type, const TVector3 &finalParticleMomentum, const MCParticleVector &products) :
    m_type(type),
    m_finalParticleMomentum(finalParticleMomentum),
    m_products(products)
{
}


} // namespace ubcc1pi
