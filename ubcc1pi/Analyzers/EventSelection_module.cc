/**
 *  @file  ubcc1pi/Analyzers/EventSelection_module.cc
 *
 *  @brief The implementation file for the event selection analyzer.
 */

#include "ubcc1pi/Analyzers/EventSelection.h"

#include "ubcc1pi/Helpers/TruthHelper.h"
#include "ubcc1pi/Helpers/AnalysisHelper.h"
#include "ubcc1pi/Helpers/DebugHelper.h"

namespace ubcc1pi
{

EventSelection::EventSelection(const art::EDAnalyzer::Table<Config> &config) :
    art::EDAnalyzer(config),
    m_config(config)
{
    // Setup the output trees
    art::ServiceHandle<art::TFileService> fileService;

    // Particle tree
    m_pParticleTree = fileService->make<TTree>("particles", "");
    m_pParticleTree->Branch("run", &m_outputParticle.m_run);
    m_pParticleTree->Branch("subRun", &m_outputParticle.m_subRun);
    m_pParticleTree->Branch("event", &m_outputParticle.m_event);
    m_pParticleTree->Branch("isSignal", &m_outputParticle.m_isSignal);
    
    m_pParticleTree->Branch("didSelectionFinish", &m_outputParticle.m_didSelectionFinish);
    m_pParticleTree->Branch("eventCC1PiScore", &m_outputParticle.m_eventCC1PiScore);
    m_pParticleTree->Branch("nProtonSelected", &m_outputParticle.m_nProtonSelected);

    m_pParticleTree->Branch("hasMatchedMCParticle", &m_outputParticle.m_hasMatchedMCParticle);
    m_pParticleTree->Branch("truePdgCode", &m_outputParticle.m_truePdgCode);
    m_pParticleTree->Branch("trueMomentum", &m_outputParticle.m_trueMomentum);
    m_pParticleTree->Branch("trueMatchPurity", &m_outputParticle.m_trueMatchPurity);
    m_pParticleTree->Branch("trueMatchCompleteness", &m_outputParticle.m_trueMatchCompleteness);

    m_pParticleTree->Branch("trackShowerScore", &m_outputParticle.m_trackShowerScore);
    m_pParticleTree->Branch("protonMIPScore", &m_outputParticle.m_protonMIPScore);
    m_pParticleTree->Branch("muonPionScore", &m_outputParticle.m_muonPionScore);

    m_pParticleTree->Branch("muonLikelihood", &m_outputParticle.m_muonLikelihood);
    m_pParticleTree->Branch("pionLikelihood", &m_outputParticle.m_pionLikelihood);
    m_pParticleTree->Branch("protonLikelihood", &m_outputParticle.m_protonLikelihood);
    
    m_pParticleTree->Branch("isMuonCandidate", &m_outputParticle.m_isMuonCandidate);
    m_pParticleTree->Branch("isPionCandidate", &m_outputParticle.m_isPionCandidate);
    m_pParticleTree->Branch("isProtonCandidate", &m_outputParticle.m_isProtonCandidate);

    // Event tree
    m_pEventTree = fileService->make<TTree>("events", "");

    m_pEventTree->Branch("run", &m_outputEvent.m_run);
    m_pEventTree->Branch("subRun", &m_outputEvent.m_subRun);
    m_pEventTree->Branch("event", &m_outputEvent.m_event);

    m_pEventTree->Branch("isSignal", &m_outputEvent.m_isSignal);
    m_pEventTree->Branch("interaction", &m_outputEvent.m_interaction);
    m_pEventTree->Branch("isNuFiducial", &m_outputEvent.m_isNuFiducial);
    m_pEventTree->Branch("nuE", &m_outputEvent.m_nuE);

    m_pEventTree->Branch("nMuMinus", &m_outputEvent.m_nMuMinus);
    m_pEventTree->Branch("nMuPlus", &m_outputEvent.m_nMuPlus);
    m_pEventTree->Branch("nPiPlus", &m_outputEvent.m_nPiPlus);
    m_pEventTree->Branch("nPiMinus", &m_outputEvent.m_nPiMinus);
    m_pEventTree->Branch("nKPlus", &m_outputEvent.m_nKPlus);
    m_pEventTree->Branch("nKMinus", &m_outputEvent.m_nKMinus);
    m_pEventTree->Branch("nProton", &m_outputEvent.m_nProton);
    m_pEventTree->Branch("nNeutron", &m_outputEvent.m_nNeutron);
    m_pEventTree->Branch("nPhoton", &m_outputEvent.m_nPhoton);
    m_pEventTree->Branch("nTotal", &m_outputEvent.m_nTotal);

    m_pEventTree->Branch("didSelectionFinish", &m_outputEvent.m_didSelectionFinish);
    m_pEventTree->Branch("cc1piScore", &m_outputEvent.m_cc1piScore);
    m_pEventTree->Branch("nProtonSelected", &m_outputEvent.m_nProtonSelected);
    m_pEventTree->Branch("nFinalStatePFParticles", &m_outputEvent.m_nFinalStatePFParticles);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void EventSelection::analyze(const art::Event &event)
{
    // Extract the labels
    const auto mcTruthLabel = m_config().MCTruthLabel();
    const auto mcParticleLabel = m_config().MCParticleLabel();
    const auto backtrackerLabel = m_config().BacktrackerLabel();
    const auto pfParticleLabel = m_config().PFParticleLabel();
    const auto trackLabel = m_config().TrackLabel();
    const auto pidLabel = m_config().PIDLabel();
    
    // Get the final state PFParticles
    const auto allPFParticles = CollectionHelper::GetCollection<recob::PFParticle>(event, pfParticleLabel);
    const auto finalStates = RecoHelper::GetNeutrinoFinalStates(allPFParticles);
    
    // Get the true neutrino interaction from the event
    const TruthHelper::Interaction interaction(event, mcTruthLabel, mcParticleLabel);
    const auto isTrueSignal = AnalysisHelper::IsCC1PiSignal(interaction);
   
    // Get the truth matching information for the PFParticles
    const auto backtrackerData = AnalysisHelper::GetBacktrackerData(event, mcTruthLabel, mcParticleLabel, backtrackerLabel, pfParticleLabel);

    // Make a CC1Pi event selector to do the legwork
    EventSelector selector(finalStates, event, pfParticleLabel, trackLabel, pidLabel);
    
    // Get the PIDs under the CC1Pi assumption
    float cc1piScore = -std::numeric_limits<float>::max();
    art::Ptr<recob::PFParticle> muonCandidate, pionCandidate;
    PFParticleVector protonCandidates;

    if (finalStates.size() >= 2)
        cc1piScore = selector.GetBestCC1PiScore(muonCandidate, pionCandidate, protonCandidates);

    // Output event level information
    m_outputEvent.m_run = event.run();
    m_outputEvent.m_subRun = event.subRun();
    m_outputEvent.m_event = event.event();
    m_outputEvent.m_isSignal = isTrueSignal;
    m_outputEvent.m_interaction = DebugHelper::GetInteractionString(interaction, true);
    m_outputEvent.m_isNuFiducial = AnalysisHelper::IsNeutrinoVertexFiducial(interaction);
    m_outputEvent.m_nuE = interaction.GetNeutrino().E();

    const auto mcParticles = AnalysisHelper::GetReconstructableFinalStates(interaction);
    m_outputEvent.m_nMuMinus = AnalysisHelper::CountParticlesWithPDG(mcParticles, 13);
    m_outputEvent.m_nMuPlus = AnalysisHelper::CountParticlesWithPDG(mcParticles, -13);
    m_outputEvent.m_nPiPlus = AnalysisHelper::CountParticlesWithPDG(mcParticles, 211);
    m_outputEvent.m_nPiMinus = AnalysisHelper::CountParticlesWithPDG(mcParticles, -211);
    m_outputEvent.m_nKPlus = AnalysisHelper::CountParticlesWithPDG(mcParticles, 321);
    m_outputEvent.m_nKMinus = AnalysisHelper::CountParticlesWithPDG(mcParticles, -321);  
    m_outputEvent.m_nProton = AnalysisHelper::CountParticlesWithPDG(mcParticles, 2212); 
    m_outputEvent.m_nNeutron = AnalysisHelper::CountParticlesWithPDG(mcParticles, 2112);
    m_outputEvent.m_nPhoton = AnalysisHelper::CountParticlesWithPDG(mcParticles, 22);
    m_outputEvent.m_nTotal = mcParticles.size();

    m_outputEvent.m_cc1piScore = cc1piScore;
    m_outputEvent.m_didSelectionFinish = (cc1piScore >= 0.f);
    m_outputEvent.m_nProtonSelected = protonCandidates.size();
    m_outputEvent.m_nFinalStatePFParticles = finalStates.size();

    // Output particles
    m_outputParticle.m_run = event.run();
    m_outputParticle.m_subRun = event.subRun();
    m_outputParticle.m_event = event.event();
    m_outputParticle.m_isSignal = isTrueSignal;

    for (const auto &pfParticle : finalStates)
    {
        // Event-level selection information
        m_outputParticle.m_eventCC1PiScore = cc1piScore;
        m_outputParticle.m_didSelectionFinish = (cc1piScore >= 0.f);
        m_outputParticle.m_nProtonSelected = protonCandidates.size();

        // Particle-level selection information
        m_outputParticle.m_trackShowerScore = selector.GetTrackShowerScore(pfParticle);
        m_outputParticle.m_protonMIPScore = selector.GetProtonMIPScore(pfParticle);
        m_outputParticle.m_muonPionScore = selector.GetMuonPionScore(pfParticle);

        m_outputParticle.m_muonLikelihood = selector.GetMuonLikelihood(pfParticle);
        m_outputParticle.m_pionLikelihood = selector.GetPionLikelihood(pfParticle);
        m_outputParticle.m_protonLikelihood = selector.GetProtonLikelihood(pfParticle);

        m_outputParticle.m_isMuonCandidate = false;
        m_outputParticle.m_isPionCandidate = false;
        m_outputParticle.m_isProtonCandidate = false;

        if (m_outputParticle.m_didSelectionFinish)
        {
            m_outputParticle.m_isMuonCandidate = (pfParticle == muonCandidate);
            m_outputParticle.m_isPionCandidate = (pfParticle == pionCandidate);

            for (const auto &protonCandidate : protonCandidates)
            {
                if (pfParticle != protonCandidate)
                    continue;

                m_outputParticle.m_isProtonCandidate = true;
                break;
            }
        }

        // MCParticle matching information
        m_outputParticle.m_hasMatchedMCParticle = false;
        m_outputParticle.m_truePdgCode = -std::numeric_limits<int>::max();
        m_outputParticle.m_trueMomentum = -std::numeric_limits<float>::max();
        m_outputParticle.m_trueMatchPurity = -std::numeric_limits<float>::max();
        m_outputParticle.m_trueMatchCompleteness = -std::numeric_limits<float>::max();
        
        try
        {
            const auto mcParticle = backtrackerData.GetBestMatchedMCParticle(pfParticle);
            const auto purity = backtrackerData.GetMatchPurity(pfParticle, mcParticle);
            const auto completeness = backtrackerData.GetMatchCompleteness(pfParticle, mcParticle);
        
            m_outputParticle.m_hasMatchedMCParticle = true;
            m_outputParticle.m_truePdgCode = mcParticle->PdgCode();
            m_outputParticle.m_trueMomentum = mcParticle->P();
            m_outputParticle.m_trueMatchPurity = purity;
            m_outputParticle.m_trueMatchCompleteness = completeness;
        }
        catch (const cet::exception &)
        {
        }

        m_pParticleTree->Fill();
    }
        
    m_pEventTree->Fill();
}

} // namespace ubcc1pi
