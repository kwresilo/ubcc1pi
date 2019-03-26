/**
 *  @file  ubcc1pi/Analyzers/TruthStudy_module.cc
 *
 *  @brief The implementation file for the truth study analyzer.
 */

#include "ubcc1pi/Analyzers/TruthStudy.h"

#include "ubcc1pi/Helpers/TruthHelper.h"
#include "ubcc1pi/Helpers/AnalysisHelper.h"
#include "ubcc1pi/Helpers/DebugHelper.h"

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
    m_pSignalTree->Branch("run",         &m_signalOutput.m_run);
    m_pSignalTree->Branch("subRun",      &m_signalOutput.m_subRun);
    m_pSignalTree->Branch("event",       &m_signalOutput.m_event);
    m_pSignalTree->Branch("interaction", &m_signalOutput.m_interaction);
    m_pSignalTree->Branch("nProton",     &m_signalOutput.m_nProton);
    m_pSignalTree->Branch("nuE",         &m_signalOutput.m_nuE);
    m_pSignalTree->Branch("muMomX",      &m_signalOutput.m_muMomX);
    m_pSignalTree->Branch("muMomY",      &m_signalOutput.m_muMomY);
    m_pSignalTree->Branch("muMomZ",      &m_signalOutput.m_muMomZ);
    m_pSignalTree->Branch("piMomX",      &m_signalOutput.m_piMomX);
    m_pSignalTree->Branch("piMomY",      &m_signalOutput.m_piMomY);
    m_pSignalTree->Branch("piMomZ",      &m_signalOutput.m_piMomZ);
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
    
    DebugHelper::Print(interaction);

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
    for (const auto &particle : particles)
    {
        if (particle->PdgCode() == 13)
            muMom = particle->Momentum().Vect();
        
        if (particle->PdgCode() == 211)
            piMom = particle->Momentum().Vect();
    }

    m_signalOutput.m_muMomX = muMom.X();
    m_signalOutput.m_muMomY = muMom.Y();
    m_signalOutput.m_muMomZ = muMom.Z();
    
    m_signalOutput.m_piMomX = piMom.X();
    m_signalOutput.m_piMomY = piMom.Y();
    m_signalOutput.m_piMomZ = piMom.Z();
    
    // Fill the trees
    m_pSignalTree->Fill();
}

} // namespace ubcc1pi
