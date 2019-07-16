/**
 *  @file  ubcc1pi/Analyzers/TruthStudy_module.cc
 *
 *  @brief The implementation file for the truth study analyzer.
 */

#include "ubcc1pi/Analyzers/TruthStudy.h"

#include "ubcc1pi/Helpers/TruthHelper.h"
#include "ubcc1pi/Helpers/AnalysisHelper.h"
#include "ubcc1pi/Helpers/DebugHelper.h"
#include "ubcc1pi/Helpers/BacktrackHelper.h"

namespace ubcc1pi
{

TruthStudy::TruthStudy(const art::EDAnalyzer::Table<Config> &config) :
    art::EDAnalyzer(config),
    m_config(config)
{
    // Setup the output trees
    art::ServiceHandle<art::TFileService> fileService;
    
    m_pInteractionTree = fileService->make<TTree>("interactions", "");
    m_pInteractionTree->Branch("run",          &m_interactionOutput.m_run);
    m_pInteractionTree->Branch("subRun",       &m_interactionOutput.m_subRun);
    m_pInteractionTree->Branch("event",        &m_interactionOutput.m_event);
    m_pInteractionTree->Branch("interaction",  &m_interactionOutput.m_interaction);
    m_pInteractionTree->Branch("isNuFiducial", &m_interactionOutput.m_isNuFiducial);
    m_pInteractionTree->Branch("isSignal",     &m_interactionOutput.m_isSignal);
    m_pInteractionTree->Branch("nMuMinus",     &m_interactionOutput.m_nMuMinus);
    m_pInteractionTree->Branch("nMuPlus",      &m_interactionOutput.m_nMuPlus);
    m_pInteractionTree->Branch("nPiPlus",      &m_interactionOutput.m_nPiPlus);
    m_pInteractionTree->Branch("nPiMinus",     &m_interactionOutput.m_nPiMinus);
    m_pInteractionTree->Branch("nKPlus",       &m_interactionOutput.m_nKPlus);
    m_pInteractionTree->Branch("nKMinus",      &m_interactionOutput.m_nKMinus);
    m_pInteractionTree->Branch("nProton",      &m_interactionOutput.m_nProton);
    m_pInteractionTree->Branch("nNeutron",     &m_interactionOutput.m_nNeutron);
    m_pInteractionTree->Branch("nPhoton",      &m_interactionOutput.m_nPhoton);
    m_pInteractionTree->Branch("nTotal",       &m_interactionOutput.m_nTotal);
    m_pInteractionTree->Branch("nuE",          &m_interactionOutput.m_nuE);
    
    m_pSignalTree = fileService->make<TTree>("signalInteractions", "");
    m_pSignalTree->Branch("run",                          &m_signalOutput.m_run);
    m_pSignalTree->Branch("subRun",                       &m_signalOutput.m_subRun);
    m_pSignalTree->Branch("event",                        &m_signalOutput.m_event);
    m_pSignalTree->Branch("interaction",                  &m_signalOutput.m_interaction);
    m_pSignalTree->Branch("nProton",                      &m_signalOutput.m_nProton);
    m_pSignalTree->Branch("nuE",                          &m_signalOutput.m_nuE);
    m_pSignalTree->Branch("muMomX",                       &m_signalOutput.m_muMomX);
    m_pSignalTree->Branch("muMomY",                       &m_signalOutput.m_muMomY);
    m_pSignalTree->Branch("muMomZ",                       &m_signalOutput.m_muMomZ);
    m_pSignalTree->Branch("piMomX",                       &m_signalOutput.m_piMomX);
    m_pSignalTree->Branch("piMomY",                       &m_signalOutput.m_piMomY);
    m_pSignalTree->Branch("piMomZ",                       &m_signalOutput.m_piMomZ);
    m_pSignalTree->Branch("piNElasticScatters",           &m_signalOutput.m_piNElasticScatters);
    m_pSignalTree->Branch("piNInelasticScatters",         &m_signalOutput.m_piNInelasticScatters);
    m_pSignalTree->Branch("piNScatters",                  &m_signalOutput.m_piNScatters);
    m_pSignalTree->Branch("piScatterIsElasticVect",       &m_signalOutput.m_piScatterIsElasticVect);
    m_pSignalTree->Branch("piScatterCosThetaVect",        &m_signalOutput.m_piScatterCosThetaVect);
    m_pSignalTree->Branch("piScatterMomFracLostVect",     &m_signalOutput.m_piScatterMomFracLostVect);
    m_pSignalTree->Branch("piInitalMom",                  &m_signalOutput.m_piInitalMom);
    m_pSignalTree->Branch("piScatteredMom",               &m_signalOutput.m_piScatteredMom);
    m_pSignalTree->Branch("piFinalMom",                   &m_signalOutput.m_piFinalMom);
    m_pSignalTree->Branch("piEndState",                   &m_signalOutput.m_piEndState);
    m_pSignalTree->Branch("piEndStateProductsHitWeightU", &m_signalOutput.m_piEndStateProductsHitWeightU);
    m_pSignalTree->Branch("piEndStateProductsHitWeightV", &m_signalOutput.m_piEndStateProductsHitWeightV);
    m_pSignalTree->Branch("piEndStateProductsHitWeightW", &m_signalOutput.m_piEndStateProductsHitWeightW);
}

// -----------------------------------------------------------------------------------------------------------------------------------------
 
void TruthStudy::analyze(const art::Event &event)
{
    // Extract the labels
    const auto mcTruthLabel = m_config().MCTruthLabel();
    const auto mcParticleLabel = m_config().MCParticleLabel();

    // Get the neutrino interaction from the event
    const TruthHelper::Interaction interaction(event, mcTruthLabel, mcParticleLabel);
    
    // Output the interaction information
    m_interactionOutput.m_run = event.run();
    m_interactionOutput.m_subRun = event.subRun();
    m_interactionOutput.m_event = event.event();
    m_interactionOutput.m_interaction = DebugHelper::GetInteractionString(interaction, true);
    m_interactionOutput.m_isNuFiducial = AnalysisHelper::IsNeutrinoVertexFiducial(interaction);
    m_interactionOutput.m_isSignal = AnalysisHelper::IsCC1PiSignal(interaction);
    m_interactionOutput.m_nuE = interaction.GetNeutrino().E();

    // Count the particles by PDG
    const auto particles = AnalysisHelper::GetReconstructableFinalStates(interaction);
    m_interactionOutput.m_nMuMinus = AnalysisHelper::CountParticlesWithPDG(particles, 13);
    m_interactionOutput.m_nMuPlus = AnalysisHelper::CountParticlesWithPDG(particles, -13);
    m_interactionOutput.m_nPiPlus = AnalysisHelper::CountParticlesWithPDG(particles, 211);
    m_interactionOutput.m_nPiMinus = AnalysisHelper::CountParticlesWithPDG(particles, -211);
    m_interactionOutput.m_nKPlus = AnalysisHelper::CountParticlesWithPDG(particles, 321);
    m_interactionOutput.m_nKMinus = AnalysisHelper::CountParticlesWithPDG(particles, -321);  
    m_interactionOutput.m_nProton = AnalysisHelper::CountParticlesWithPDG(particles, 2212); 
    m_interactionOutput.m_nNeutron = AnalysisHelper::CountParticlesWithPDG(particles, 2112);
    m_interactionOutput.m_nPhoton = AnalysisHelper::CountParticlesWithPDG(particles, 22);
    m_interactionOutput.m_nTotal = particles.size();
    
    m_pInteractionTree->Fill();

    // Only fill the signal tree if this is a signal event!
    if (!AnalysisHelper::IsCC1PiSignal(interaction))
        return;
    
    m_signalOutput.m_run = event.run();
    m_signalOutput.m_subRun = event.subRun();
    m_signalOutput.m_event = event.event();
    m_signalOutput.m_interaction = m_interactionOutput.m_interaction;
    m_signalOutput.m_nProton = m_interactionOutput.m_nProton;
    m_signalOutput.m_nuE = m_interactionOutput.m_nuE;

    // Get the momenta
    TVector3 muMom, piMom;
    art::Ptr<simb::MCParticle> pion, muon;

    for (const auto &particle : particles)
    {
        if (particle->PdgCode() == 13)
        {
            muMom = particle->Momentum().Vect();
            muon = particle;
        }
        
        if (particle->PdgCode() == 211)
        {
            piMom = particle->Momentum().Vect();
            pion = particle;
        }
    }

    m_signalOutput.m_muMomX = muMom.X();
    m_signalOutput.m_muMomY = muMom.Y();
    m_signalOutput.m_muMomZ = muMom.Z();
    
    m_signalOutput.m_piMomX = piMom.X();
    m_signalOutput.m_piMomY = piMom.Y();
    m_signalOutput.m_piMomZ = piMom.Z();
    
    // Work out what happened to the pion
    const auto mcParticleMap = TruthHelper::GetMCParticleMap(interaction.GetAllMCParticles());
    ScatterVector scatters;
    art::Ptr<simb::MCParticle> scatteredPion;
    this->FollowScatters(pion, mcParticleMap, scatters, scatteredPion);
    const auto endState = this->GetEndState(scatteredPion, mcParticleMap);
  
    // Output the information to the trees
    m_signalOutput.m_piNElasticScatters = 0;
    m_signalOutput.m_piNInelasticScatters = 0;
    m_signalOutput.m_piNScatters = scatters.size();
    m_signalOutput.m_piScatterIsElasticVect.clear();
    m_signalOutput.m_piScatterCosThetaVect.clear();
    m_signalOutput.m_piScatterMomFracLostVect.clear();

    for (const auto &scatter : scatters)
    {
        m_signalOutput.m_piNElasticScatters += scatter.m_isElastic ? 1 : 0;
        m_signalOutput.m_piNInelasticScatters += scatter.m_isElastic ? 0 : 1;

        m_signalOutput.m_piScatterIsElasticVect.push_back(scatter.m_isElastic);
        m_signalOutput.m_piScatterCosThetaVect.push_back(scatter.GetScatteringCosTheta());
        m_signalOutput.m_piScatterMomFracLostVect.push_back(scatter.GetMomentumFractionLost());
    }

    m_signalOutput.m_piInitalMom = pion->Momentum().Vect().Mag();
    m_signalOutput.m_piScatteredMom = scatteredPion->Momentum().Vect().Mag();
    m_signalOutput.m_piFinalMom = endState.m_finalPionMomentum.Mag();
    m_signalOutput.m_piEndState = endState.m_type;

    m_signalOutput.m_piEndStateProductsHitWeightU = 0.f;
    m_signalOutput.m_piEndStateProductsHitWeightV = 0.f;
    m_signalOutput.m_piEndStateProductsHitWeightW = 0.f;
    const auto mcParticleToHits = CollectionHelper::GetReversedAssociation(BacktrackHelper::GetHitToMCParticleWeightMap(event, m_config().MCParticleLabel(), m_config().BacktrackerLabel(), endState.m_products)); 
    for (const auto &particle : endState.m_products)
    {
        m_signalOutput.m_piEndStateProductsHitWeightU += BacktrackHelper::GetHitWeightInView(particle, mcParticleToHits, geo::kU);
        m_signalOutput.m_piEndStateProductsHitWeightV += BacktrackHelper::GetHitWeightInView(particle, mcParticleToHits, geo::kV);
        m_signalOutput.m_piEndStateProductsHitWeightW += BacktrackHelper::GetHitWeightInView(particle, mcParticleToHits, geo::kW);
    }

    // Fill the trees
    m_pSignalTree->Fill();
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void TruthStudy::FollowScatters(const art::Ptr<simb::MCParticle> &pion, const MCParticleMap &mcParticleMap, ScatterVector &scatters, art::Ptr<simb::MCParticle> &scatteredPion) const
{
    scatteredPion = pion;
    
    // First look for elastic scatters
    for (const auto &daughter : TruthHelper::GetDaughters(pion, mcParticleMap))
    {
        if (daughter->Process() == "hadElastic")
            scatters.push_back(this->GetElasticScatter(pion, daughter));
    }

    // Now follow inelastic scatters
    if (pion->EndProcess() != "pi+Inelastic")
        return;

    MCParticleVector inelasticProducts;
    art::Ptr<simb::MCParticle> finalStatePion;
    bool foundInelasticScatter = false;

    for (const auto &daughter : TruthHelper::GetDaughters(pion, mcParticleMap))
    {
        if (daughter->Process() == "pi+Inelastic")
        {
            if (daughter->PdgCode() != 211)
            {
                inelasticProducts.push_back(daughter);
                continue;
            }

            // If there are multiple pions from the inelastic collision, then don't treat it as a scatter
            if (foundInelasticScatter)
            {
                foundInelasticScatter = false;
                break;
            }
        
            finalStatePion = daughter;
            foundInelasticScatter = true;
        }
    }

    if (!foundInelasticScatter)
        return;

    // Get the last momentum point before the scatter
    TVector3 initialMomentum(-std::numeric_limits<float>::max(), -std::numeric_limits<float>::max(), -std::numeric_limits<float>::max());
    for (unsigned int i = 0; i < pion->NumberTrajectoryPoints(); ++i)
    {
        if (pion->T(i) >= finalStatePion->T())
            break;

        initialMomentum = pion->Momentum(i).Vect();
    }

    scatters.emplace_back(false, initialMomentum, finalStatePion->Momentum().Vect(), finalStatePion, inelasticProducts);
    this->FollowScatters(finalStatePion, mcParticleMap, scatters, scatteredPion);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

TruthStudy::Scatter TruthStudy::GetElasticScatter(const art::Ptr<simb::MCParticle> &incident, const art::Ptr<simb::MCParticle> &target) const
{
    if (target->Process() != "hadElastic")
        throw cet::exception("TruthStudy::GetElasticScatter") << " - Input target particle doesn't look to be from an elastic scatter" << std::endl;

    if (target->Mother() != incident->TrackId())
        throw cet::exception("TruthStudy::GetElasticScatter") << " - Input target particle isn't a daughter of the input incident particle" << std::endl;

    const unsigned int nTrajectoryPoints = incident->NumberTrajectoryPoints();
    if (nTrajectoryPoints < 2)
        throw cet::exception("TruthStudy::GetElasticScatter") << " - Incident particle has less than 2 trajectory points!" << std::endl;
    
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
        throw cet::exception("TruthStudy::GetElasticScatter") << " - Couldn't find trajectory points before and after the scatter" << std::endl;

    const auto initialMomentum(incident->Momentum(preScatterIndex).Vect());
    const auto finalMomentum(incident->Momentum(postScatterIndex).Vect());

    return TruthStudy::Scatter(true, initialMomentum, finalMomentum, incident, { target });
}

// -----------------------------------------------------------------------------------------------------------------------------------------

TruthStudy::EndState TruthStudy::GetEndState(const art::Ptr<simb::MCParticle> &pion, const MCParticleMap &mcParticleMap) const
{
    if (pion->PdgCode() != 211)
        throw cet::exception("TruthStudy::GetEndState") << " - Input MCParticle isn't a Pi+" << std::endl;

    // Follow the scatters to make sure we are really looking at the correct MCParticle
    ScatterVector scatters;
    art::Ptr<simb::MCParticle> scatteredPion;
    this->FollowScatters(pion, mcParticleMap, scatters, scatteredPion);
    if (pion != scatteredPion)
        throw cet::exception("TruthStudy::GetEndState") << " - Input MCParticle undergos a scatter before reaching it's end-state" << std::endl;

    bool hasPi0 = false;
    bool hasDecayMuon = false;
    bool hasDecayMuonNeutrino = false;

    MCParticleVector products;
    for (const auto &daughter : TruthHelper::GetDaughters(pion, mcParticleMap))
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
    EndState::Type type = EndState::Type::Other;

    if (products.empty())
    {
        type = EndState::Type::None;
    }
    else if (hasDecayMuon && hasDecayMuonNeutrino && products.size() == 2)
    {
        type = EndState::Type::DecayToMuon;
    }
    else if (pion->EndProcess() == "pi+Inelastic")
    {
        type = hasPi0 ? EndState::Type::Pi0ChargeExchange : EndState::Type::InelasticAbsorption;
    }
    
    const auto finalPionMomentum = pion->Momentum(std::max(static_cast<unsigned int>(0), pion->NumberTrajectoryPoints() - 2)).Vect();
    return TruthStudy::EndState(type, finalPionMomentum, products);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void TruthStudy::PrintHierarchy(const art::Ptr<simb::MCParticle> &particle, const MCParticleMap &mcParticleMap, const art::Event &event) const
{
    DebugHelper::Print(particle);
    
    // Get the hits associated to all of the daughters of the primary particle folding in the downstream particles in their hierarchies 
    const auto daughters = TruthHelper::GetDaughters(particle, mcParticleMap);
    const auto mcParticleToHits = CollectionHelper::GetReversedAssociation(BacktrackHelper::GetHitToMCParticleWeightMap(event, m_config().MCParticleLabel(), m_config().BacktrackerLabel(), daughters));

    for (const auto &daughter : TruthHelper::GetDaughters(particle, mcParticleMap))
    {
        DebugHelper::Print(daughter, 1);
        DebugHelper::PrintProperty("Hit weight U", BacktrackHelper::GetHitWeightInView(daughter, mcParticleToHits, geo::kU), 1);
        DebugHelper::PrintProperty("Hit weight V", BacktrackHelper::GetHitWeightInView(daughter, mcParticleToHits, geo::kV), 1);
        DebugHelper::PrintProperty("Hit weight W", BacktrackHelper::GetHitWeightInView(daughter, mcParticleToHits, geo::kW), 1);
    }
}

// -----------------------------------------------------------------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------------------------------------------------------------
                
TruthStudy::Scatter::Scatter(const bool isElastic, const TVector3 initialMomentum, const TVector3 finalMomentum, const art::Ptr<simb::MCParticle> finalParticle, const MCParticleVector products) : 
    m_isElastic(isElastic),
    m_initialMomentum(initialMomentum),
    m_finalMomentum(finalMomentum),
    m_finalParticle(finalParticle),
    m_products(products)
{
}

// -----------------------------------------------------------------------------------------------------------------------------------------
                
float TruthStudy::Scatter::GetScatteringCosTheta() const
{
    if (m_initialMomentum.Mag() * m_finalMomentum.Mag() < std::numeric_limits<float>::epsilon())
        return -std::numeric_limits<float>::max();

    return m_initialMomentum.Dot(m_finalMomentum) / (m_initialMomentum.Mag() * m_finalMomentum.Mag());
}

// -----------------------------------------------------------------------------------------------------------------------------------------
                
float TruthStudy::Scatter::GetMomentumFractionLost() const
{
    if (m_initialMomentum.Mag() * m_finalMomentum.Mag() < std::numeric_limits<float>::epsilon())
        return -std::numeric_limits<float>::max();

    return (m_initialMomentum.Mag() - m_finalMomentum.Mag()) / m_initialMomentum.Mag();
}

// -----------------------------------------------------------------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------------------------------------------------------------
                
TruthStudy::EndState::EndState(const Type type, const TVector3 &finalPionMomentum, const MCParticleVector &products) :
    m_type(type),
    m_finalPionMomentum(finalPionMomentum),
    m_products(products)
{
}

} // namespace ubcc1pi
