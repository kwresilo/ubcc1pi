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

    // Event tree
    m_pEventTree = fileService->make<TTree>("events", "");

    // Set up the output branches
    m_pEventTree->Branch("run", &m_outputEvent.m_run);
    m_pEventTree->Branch("subRun", &m_outputEvent.m_subRun);
    m_pEventTree->Branch("event", &m_outputEvent.m_event);

    m_pEventTree->Branch("isSignal", &m_outputEvent.m_isSignal);
    m_pEventTree->Branch("interaction", &m_outputEvent.m_interaction);
    m_pEventTree->Branch("trueNuE", &m_outputEvent.m_trueNuE);
    m_pEventTree->Branch("trueNuVtx", &m_outputEvent.m_trueNuVtx);
    m_pEventTree->Branch("isTrueNuFiducial", &m_outputEvent.m_isTrueNuFiducial);
    
    m_pEventTree->Branch("nMuMinus", &m_outputEvent.m_nMuMinus);
    m_pEventTree->Branch("nMuPlus", &m_outputEvent.m_nMuPlus);
    m_pEventTree->Branch("nPiPlus", &m_outputEvent.m_nPiPlus);
    m_pEventTree->Branch("nPiMinus", &m_outputEvent.m_nPiMinus);
    m_pEventTree->Branch("nKPlus", &m_outputEvent.m_nKPlus);
    m_pEventTree->Branch("nKMinus", &m_outputEvent.m_nKMinus);
    m_pEventTree->Branch("nProton", &m_outputEvent.m_nProton);
    m_pEventTree->Branch("nNeutron", &m_outputEvent.m_nNeutron);
    m_pEventTree->Branch("nPhoton", &m_outputEvent.m_nPhoton);
    m_pEventTree->Branch("nElectron", &m_outputEvent.m_nElectron);
    m_pEventTree->Branch("nPositron", &m_outputEvent.m_nPositron);
    m_pEventTree->Branch("nTotal", &m_outputEvent.m_nTotal);

    m_pEventTree->Branch("trueMuEnergy", &m_outputEvent.m_trueMuEnergy);
    m_pEventTree->Branch("trueMuTheta", &m_outputEvent.m_trueMuTheta);
    m_pEventTree->Branch("trueMuPhi", &m_outputEvent.m_trueMuPhi);
    m_pEventTree->Branch("truePiEnergy", &m_outputEvent.m_truePiEnergy);
    m_pEventTree->Branch("truePiTheta", &m_outputEvent.m_truePiTheta);
    m_pEventTree->Branch("truePiPhi", &m_outputEvent.m_truePiPhi);
    m_pEventTree->Branch("trueMuPiAngle", &m_outputEvent.m_trueMuPiAngle);
            
    m_pEventTree->Branch("nMCParticles", &m_outputEvent.m_nMCParticles);
    m_pEventTree->Branch("mcpIdVect", &m_outputEvent.m_mcpIdVect);
    m_pEventTree->Branch("mcpPDGVect", &m_outputEvent.m_mcpPDGVect);
    m_pEventTree->Branch("mcpIsTargetFinalStateVect", &m_outputEvent.m_mcpIsTargetFinalStateVect);
    m_pEventTree->Branch("mcpProcessVect", &m_outputEvent.m_mcpProcessVect);
    m_pEventTree->Branch("mcpMotherIdVect", &m_outputEvent.m_mcpMotherIdVect);
    m_pEventTree->Branch("mcpNDaughtersVect", &m_outputEvent.m_mcpNDaughtersVect);
    m_pEventTree->Branch("mcpDaughterIdsVect", &m_outputEvent.m_mcpDaughterIdsVect);
    m_pEventTree->Branch("mcpEnergyVect", &m_outputEvent.m_mcpEnergyVect);
    m_pEventTree->Branch("mcpMomentumVect", &m_outputEvent.m_mcpMomentumVect);
    m_pEventTree->Branch("mcpMomentumXVect", &m_outputEvent.m_mcpMomentumXVect);
    m_pEventTree->Branch("mcpMomentumYVect", &m_outputEvent.m_mcpMomentumYVect);
    m_pEventTree->Branch("mcpMomentumZVect", &m_outputEvent.m_mcpMomentumZVect);
    m_pEventTree->Branch("mcpNHitsUVect", &m_outputEvent.m_mcpNHitsUVect);
    m_pEventTree->Branch("mcpNHitsVVect", &m_outputEvent.m_mcpNHitsVVect);
    m_pEventTree->Branch("mcpNHitsWVect", &m_outputEvent.m_mcpNHitsWVect);
    m_pEventTree->Branch("mcpNGoodHitsUVect", &m_outputEvent.m_mcpNGoodHitsUVect);
    m_pEventTree->Branch("mcpNGoodHitsVVect", &m_outputEvent.m_mcpNGoodHitsVVect);
    m_pEventTree->Branch("mcpNGoodHitsWVect", &m_outputEvent.m_mcpNGoodHitsWVect);
    m_pEventTree->Branch("mcpHitWeightUVect", &m_outputEvent.m_mcpHitWeightUVect);
    m_pEventTree->Branch("mcpHitWeightVVect", &m_outputEvent.m_mcpHitWeightVVect);
    m_pEventTree->Branch("mcpHitWeightWVect", &m_outputEvent.m_mcpHitWeightWVect);

    m_pEventTree->Branch("hasRecoNeutrino", &m_outputEvent.m_hasRecoNeutrino);
    m_pEventTree->Branch("recoNuVtx", &m_outputEvent.m_recoNuVtx);
    m_pEventTree->Branch("isRecoNuFiducial", &m_outputEvent.m_isRecoNuFiducial);
    m_pEventTree->Branch("topologicalScore", &m_outputEvent.m_topologicalScore);
    m_pEventTree->Branch("nFinalStatePFPs", &m_outputEvent.m_nFinalStatePFPs);

    m_pEventTree->Branch("hasMatchedMCParticleVect", &m_outputEvent.m_hasMatchedMCParticleVect);
    m_pEventTree->Branch("matchedMCParticleIdVect", &m_outputEvent.m_matchedMCParticleIdVect);
    m_pEventTree->Branch("truePdgCodeVect", &m_outputEvent.m_truePdgCodeVect);
    m_pEventTree->Branch("truthMatchCompletenessVect", &m_outputEvent.m_truthMatchCompletenessVect);
    m_pEventTree->Branch("truthMatchPurityVect", &m_outputEvent.m_truthMatchPurityVect);
    m_pEventTree->Branch("trueEnergyVect", &m_outputEvent.m_trueEnergyVect);
    m_pEventTree->Branch("trueKEVect", &m_outputEvent.m_trueKEVect);
    m_pEventTree->Branch("trueMomentumXVect", &m_outputEvent.m_trueMomentumXVect);
    m_pEventTree->Branch("trueMomentumYVect", &m_outputEvent.m_trueMomentumYVect);
    m_pEventTree->Branch("trueMomentumZVect", &m_outputEvent.m_trueMomentumZVect);
    m_pEventTree->Branch("trueStartXVect", &m_outputEvent.m_trueStartXVect);
    m_pEventTree->Branch("trueStartYVect", &m_outputEvent.m_trueStartYVect);
    m_pEventTree->Branch("trueStartZVect", &m_outputEvent.m_trueStartZVect);
    
    m_pEventTree->Branch("hasTrackInfoVect", &m_outputEvent.m_hasTrackInfoVect);

    m_pEventTree->Branch("startXVect", &m_outputEvent.m_startXVect);
    m_pEventTree->Branch("startYVect", &m_outputEvent.m_startYVect);
    m_pEventTree->Branch("startZVect", &m_outputEvent.m_startZVect);
    m_pEventTree->Branch("endXVect", &m_outputEvent.m_endXVect);
    m_pEventTree->Branch("endYVect", &m_outputEvent.m_endYVect);
    m_pEventTree->Branch("endZVect", &m_outputEvent.m_endZVect);
    m_pEventTree->Branch("directionXVect", &m_outputEvent.m_directionXVect);
    m_pEventTree->Branch("directionYVect", &m_outputEvent.m_directionYVect);
    m_pEventTree->Branch("directionZVect", &m_outputEvent.m_directionZVect);
    m_pEventTree->Branch("thetaVect", &m_outputEvent.m_thetaVect);
    m_pEventTree->Branch("phiVect", &m_outputEvent.m_phiVect);
    m_pEventTree->Branch("yzAngleVect", &m_outputEvent.m_yzAngleVect);
    m_pEventTree->Branch("lengthVect", &m_outputEvent.m_lengthVect);
    m_pEventTree->Branch("isContainedVect", &m_outputEvent.m_isContainedVect);
    
    m_pEventTree->Branch("nHitsUVect", &m_outputEvent.m_nHitsUVect);
    m_pEventTree->Branch("nHitsVVect", &m_outputEvent.m_nHitsVVect);
    m_pEventTree->Branch("nHitsWVect", &m_outputEvent.m_nHitsWVect);
    m_pEventTree->Branch("trackShowerVect", &m_outputEvent.m_trackShowerVect);
            
    m_pEventTree->Branch("hasCalorimetryInfoVect", &m_outputEvent.m_hasCalorimetryInfoVect);
    m_pEventTree->Branch("dedxPerHitUVect", &m_outputEvent.m_dedxPerHitUVect);
    m_pEventTree->Branch("dedxPerHitVVect", &m_outputEvent.m_dedxPerHitVVect);
    m_pEventTree->Branch("dedxPerHitWVect", &m_outputEvent.m_dedxPerHitWVect);
    m_pEventTree->Branch("residualRangePerHitUVect", &m_outputEvent.m_residualRangePerHitUVect);
    m_pEventTree->Branch("residualRangePerHitVVect", &m_outputEvent.m_residualRangePerHitVVect);
    m_pEventTree->Branch("residualRangePerHitWVect", &m_outputEvent.m_residualRangePerHitWVect);

    m_pEventTree->Branch("hasPIDInfoVect", &m_outputEvent.m_hasPIDInfoVect);

    m_pEventTree->Branch("isWChi2pAvailableVect", &m_outputEvent.m_isWChi2pAvailableVect);
    m_pEventTree->Branch("chi2pWVect", &m_outputEvent.m_chi2pWVect);
    m_pEventTree->Branch("isUChi2pAvailableVect", &m_outputEvent.m_isUChi2pAvailableVect);
    m_pEventTree->Branch("chi2pUVect", &m_outputEvent.m_chi2pUVect);
    m_pEventTree->Branch("isVChi2pAvailableVect", &m_outputEvent.m_isVChi2pAvailableVect);
    m_pEventTree->Branch("chi2pVVect", &m_outputEvent.m_chi2pVVect);
    m_pEventTree->Branch("isUVChi2pAvailableVect", &m_outputEvent.m_isUVChi2pAvailableVect);
    m_pEventTree->Branch("chi2pUVVect", &m_outputEvent.m_chi2pUVVect);

    m_pEventTree->Branch("isWBraggpAvailableVect", &m_outputEvent.m_isWBraggpAvailableVect);
    m_pEventTree->Branch("braggpWVect", &m_outputEvent.m_braggpWVect);
    m_pEventTree->Branch("isUBraggpAvailableVect", &m_outputEvent.m_isUBraggpAvailableVect);
    m_pEventTree->Branch("braggpUVect", &m_outputEvent.m_braggpUVect);
    m_pEventTree->Branch("isVBraggpAvailableVect", &m_outputEvent.m_isVBraggpAvailableVect);
    m_pEventTree->Branch("braggpVVect", &m_outputEvent.m_braggpVVect);
    m_pEventTree->Branch("isUVBraggpAvailableVect", &m_outputEvent.m_isUVBraggpAvailableVect);
    m_pEventTree->Branch("braggpUVVect", &m_outputEvent.m_braggpUVVect);

    m_pEventTree->Branch("isWBraggpBackwardAvailableVect", &m_outputEvent.m_isWBraggpBackwardAvailableVect);
    m_pEventTree->Branch("braggpBackwardWVect", &m_outputEvent.m_braggpBackwardWVect);
    m_pEventTree->Branch("isUBraggpBackwardAvailableVect", &m_outputEvent.m_isUBraggpBackwardAvailableVect);
    m_pEventTree->Branch("braggpBackwardUVect", &m_outputEvent.m_braggpBackwardUVect);
    m_pEventTree->Branch("isVBraggpBackwardAvailableVect", &m_outputEvent.m_isVBraggpBackwardAvailableVect);
    m_pEventTree->Branch("braggpBackwardVVect", &m_outputEvent.m_braggpBackwardVVect);
    m_pEventTree->Branch("isUVBraggpBackwardAvailableVect", &m_outputEvent.m_isUVBraggpBackwardAvailableVect);
    m_pEventTree->Branch("braggpBackwardUVVect", &m_outputEvent.m_braggpBackwardUVVect);

    m_pEventTree->Branch("isWBraggMIPAvailableVect", &m_outputEvent.m_isWBraggMIPAvailableVect);
    m_pEventTree->Branch("braggMIPWVect", &m_outputEvent.m_braggMIPWVect);
    m_pEventTree->Branch("isUBraggMIPAvailableVect", &m_outputEvent.m_isUBraggMIPAvailableVect);
    m_pEventTree->Branch("braggMIPUVect", &m_outputEvent.m_braggMIPUVect);
    m_pEventTree->Branch("isVBraggMIPAvailableVect", &m_outputEvent.m_isVBraggMIPAvailableVect);
    m_pEventTree->Branch("braggMIPVVect", &m_outputEvent.m_braggMIPVVect);
    m_pEventTree->Branch("isUVBraggMIPAvailableVect", &m_outputEvent.m_isUVBraggMIPAvailableVect);
    m_pEventTree->Branch("braggMIPUVVect", &m_outputEvent.m_braggMIPUVVect);

    m_pEventTree->Branch("isWBraggMIPBackwardAvailableVect", &m_outputEvent.m_isWBraggMIPBackwardAvailableVect);
    m_pEventTree->Branch("braggMIPBackwardWVect", &m_outputEvent.m_braggMIPBackwardWVect);
    m_pEventTree->Branch("isUBraggMIPBackwardAvailableVect", &m_outputEvent.m_isUBraggMIPBackwardAvailableVect);
    m_pEventTree->Branch("braggMIPBackwardUVect", &m_outputEvent.m_braggMIPBackwardUVect);
    m_pEventTree->Branch("isVBraggMIPBackwardAvailableVect", &m_outputEvent.m_isVBraggMIPBackwardAvailableVect);
    m_pEventTree->Branch("braggMIPBackwardVVect", &m_outputEvent.m_braggMIPBackwardVVect);
    m_pEventTree->Branch("isUVBraggMIPBackwardAvailableVect", &m_outputEvent.m_isUVBraggMIPBackwardAvailableVect);
    m_pEventTree->Branch("braggMIPBackwardUVVect", &m_outputEvent.m_braggMIPBackwardUVVect);

    m_pEventTree->Branch("isWBraggRatioAvailableVect", &m_outputEvent.m_isWBraggRatioAvailableVect);
    m_pEventTree->Branch("braggRatioWVect", &m_outputEvent.m_braggRatioWVect);
    m_pEventTree->Branch("isUVBraggRatioAvailableVect", &m_outputEvent.m_isUVBraggRatioAvailableVect);
    m_pEventTree->Branch("braggRatioUVVect", &m_outputEvent.m_braggRatioUVVect);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void EventSelection::analyze(const art::Event &event)
{
    this->ResetEventTree();
    this->SetEventMetadata(event);
    this->SetEventTruthInfo(event);
    this->SetRecoInfo(event);

    m_pEventTree->Fill();
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void EventSelection::ResetEventTree()
{
    m_outputEvent.m_run = -std::numeric_limits<int>::max();
    m_outputEvent.m_subRun = -std::numeric_limits<int>::max();
    m_outputEvent.m_event = -std::numeric_limits<int>::max();
    m_outputEvent.m_isSignal = false;
    m_outputEvent.m_interaction = "";
    m_outputEvent.m_trueNuE = -std::numeric_limits<float>::max();
    m_outputEvent.m_trueNuVtx = TVector3(-std::numeric_limits<float>::max(), -std::numeric_limits<float>::max(), -std::numeric_limits<float>::max());
    m_outputEvent.m_isTrueNuFiducial = false;
    m_outputEvent.m_nMuMinus = -std::numeric_limits<int>::max();
    m_outputEvent.m_nMuPlus = -std::numeric_limits<int>::max();
    m_outputEvent.m_nPiPlus = -std::numeric_limits<int>::max();
    m_outputEvent.m_nPiMinus = -std::numeric_limits<int>::max();
    m_outputEvent.m_nKPlus = -std::numeric_limits<int>::max();
    m_outputEvent.m_nKMinus = -std::numeric_limits<int>::max();
    m_outputEvent.m_nProton = -std::numeric_limits<int>::max();
    m_outputEvent.m_nNeutron = -std::numeric_limits<int>::max();
    m_outputEvent.m_nPhoton = -std::numeric_limits<int>::max();
    m_outputEvent.m_nElectron = -std::numeric_limits<int>::max();
    m_outputEvent.m_nPositron = -std::numeric_limits<int>::max();
    m_outputEvent.m_nTotal = -std::numeric_limits<int>::max();
    m_outputEvent.m_trueMuEnergy = -std::numeric_limits<float>::max();
    m_outputEvent.m_trueMuTheta = -std::numeric_limits<float>::max();
    m_outputEvent.m_trueMuPhi = -std::numeric_limits<float>::max();
    m_outputEvent.m_truePiEnergy = -std::numeric_limits<float>::max();
    m_outputEvent.m_truePiTheta = -std::numeric_limits<float>::max();
    m_outputEvent.m_truePiPhi = -std::numeric_limits<float>::max();
    m_outputEvent.m_trueMuPiAngle = -std::numeric_limits<float>::max();
    m_outputEvent.m_nMCParticles = -std::numeric_limits<int>::max();
    m_outputEvent.m_mcpIdVect.clear();
    m_outputEvent.m_mcpPDGVect.clear();
    m_outputEvent.m_mcpIsTargetFinalStateVect.clear();
    m_outputEvent.m_mcpProcessVect.clear();
    m_outputEvent.m_mcpMotherIdVect.clear();
    m_outputEvent.m_mcpNDaughtersVect.clear();
    m_outputEvent.m_mcpDaughterIdsVect.clear();
    m_outputEvent.m_mcpEnergyVect.clear();
    m_outputEvent.m_mcpMomentumVect.clear();
    m_outputEvent.m_mcpMomentumXVect.clear();
    m_outputEvent.m_mcpMomentumYVect.clear();
    m_outputEvent.m_mcpMomentumZVect.clear();
    m_outputEvent.m_mcpNHitsUVect.clear();
    m_outputEvent.m_mcpNHitsVVect.clear();
    m_outputEvent.m_mcpNHitsWVect.clear();
    m_outputEvent.m_mcpNGoodHitsUVect.clear();
    m_outputEvent.m_mcpNGoodHitsVVect.clear();
    m_outputEvent.m_mcpNGoodHitsWVect.clear();
    m_outputEvent.m_mcpHitWeightUVect.clear();
    m_outputEvent.m_mcpHitWeightVVect.clear();
    m_outputEvent.m_mcpHitWeightWVect.clear();
    m_outputEvent.m_hasRecoNeutrino = false;
    m_outputEvent.m_recoNuVtx = TVector3(-std::numeric_limits<float>::max(), -std::numeric_limits<float>::max(), -std::numeric_limits<float>::max());
    m_outputEvent.m_isRecoNuFiducial = false;
    m_outputEvent.m_topologicalScore = -std::numeric_limits<float>::max();
    m_outputEvent.m_nFinalStatePFPs = -std::numeric_limits<int>::max();
    m_outputEvent.m_hasMatchedMCParticleVect.clear();
    m_outputEvent.m_matchedMCParticleIdVect.clear();
    m_outputEvent.m_truePdgCodeVect.clear();
    m_outputEvent.m_truthMatchCompletenessVect.clear();
    m_outputEvent.m_truthMatchPurityVect.clear();
    m_outputEvent.m_trueEnergyVect.clear();
    m_outputEvent.m_trueKEVect.clear();
    m_outputEvent.m_trueMomentumXVect.clear();
    m_outputEvent.m_trueMomentumYVect.clear();
    m_outputEvent.m_trueMomentumZVect.clear();
    m_outputEvent.m_trueStartXVect.clear();
    m_outputEvent.m_trueStartYVect.clear();
    m_outputEvent.m_trueStartZVect.clear();
    m_outputEvent.m_nHitsUVect.clear();
    m_outputEvent.m_nHitsVVect.clear();
    m_outputEvent.m_nHitsWVect.clear();
    m_outputEvent.m_trackShowerVect.clear();
    m_outputEvent.m_hasTrackInfoVect.clear();
    m_outputEvent.m_startXVect.clear();
    m_outputEvent.m_startYVect.clear();
    m_outputEvent.m_startZVect.clear();
    m_outputEvent.m_endXVect.clear();
    m_outputEvent.m_endYVect.clear();
    m_outputEvent.m_endZVect.clear();
    m_outputEvent.m_directionXVect.clear();
    m_outputEvent.m_directionYVect.clear();
    m_outputEvent.m_directionZVect.clear();
    m_outputEvent.m_thetaVect.clear();
    m_outputEvent.m_phiVect.clear();
    m_outputEvent.m_yzAngleVect.clear();
    m_outputEvent.m_lengthVect.clear();
    m_outputEvent.m_isContainedVect.clear();
    m_outputEvent.m_hasCalorimetryInfoVect.clear();
    m_outputEvent.m_dedxPerHitUVect.clear();
    m_outputEvent.m_dedxPerHitVVect.clear();
    m_outputEvent.m_dedxPerHitWVect.clear();
    m_outputEvent.m_residualRangePerHitUVect.clear();
    m_outputEvent.m_residualRangePerHitVVect.clear();
    m_outputEvent.m_residualRangePerHitWVect.clear();
    m_outputEvent.m_hasPIDInfoVect.clear();
    m_outputEvent.m_isWChi2pAvailableVect.clear();
    m_outputEvent.m_chi2pWVect.clear();
    m_outputEvent.m_isUChi2pAvailableVect.clear();
    m_outputEvent.m_chi2pUVect.clear();
    m_outputEvent.m_isVChi2pAvailableVect.clear();
    m_outputEvent.m_chi2pVVect.clear();
    m_outputEvent.m_isUVChi2pAvailableVect.clear();
    m_outputEvent.m_chi2pUVVect.clear();
    m_outputEvent.m_isWBraggpAvailableVect.clear();
    m_outputEvent.m_braggpWVect.clear();
    m_outputEvent.m_isUBraggpAvailableVect.clear();
    m_outputEvent.m_braggpUVect.clear();
    m_outputEvent.m_isVBraggpAvailableVect.clear();
    m_outputEvent.m_braggpVVect.clear();
    m_outputEvent.m_isUVBraggpAvailableVect.clear();
    m_outputEvent.m_braggpUVVect.clear();
    m_outputEvent.m_isWBraggpBackwardAvailableVect.clear();
    m_outputEvent.m_braggpBackwardWVect.clear();
    m_outputEvent.m_isUBraggpBackwardAvailableVect.clear();
    m_outputEvent.m_braggpBackwardUVect.clear();
    m_outputEvent.m_isVBraggpBackwardAvailableVect.clear();
    m_outputEvent.m_braggpBackwardVVect.clear();
    m_outputEvent.m_isUVBraggpBackwardAvailableVect.clear();
    m_outputEvent.m_braggpBackwardUVVect.clear();
    m_outputEvent.m_isWBraggMIPAvailableVect.clear();
    m_outputEvent.m_braggMIPWVect.clear();
    m_outputEvent.m_isUBraggMIPAvailableVect.clear();
    m_outputEvent.m_braggMIPUVect.clear();
    m_outputEvent.m_isVBraggMIPAvailableVect.clear();
    m_outputEvent.m_braggMIPVVect.clear();
    m_outputEvent.m_isUVBraggMIPAvailableVect.clear();
    m_outputEvent.m_braggMIPUVVect.clear();
    m_outputEvent.m_braggMIPBackwardWVect.clear();
    m_outputEvent.m_isUBraggMIPBackwardAvailableVect.clear();
    m_outputEvent.m_braggMIPBackwardUVect.clear();
    m_outputEvent.m_isVBraggMIPBackwardAvailableVect.clear();
    m_outputEvent.m_braggMIPBackwardVVect.clear();
    m_outputEvent.m_isUVBraggMIPBackwardAvailableVect.clear();
    m_outputEvent.m_braggMIPBackwardUVVect.clear();
    m_outputEvent.m_isWBraggRatioAvailableVect.clear();
    m_outputEvent.m_braggRatioWVect.clear();
    m_outputEvent.m_isUVBraggRatioAvailableVect.clear();
    m_outputEvent.m_braggRatioUVVect.clear();
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void EventSelection::SetEventMetadata(const art::Event &event)
{
    m_outputEvent.m_run = event.run();
    m_outputEvent.m_subRun = event.subRun();
    m_outputEvent.m_event = event.event();
}
    
// -----------------------------------------------------------------------------------------------------------------------------------------

void EventSelection::SetEventTruthInfo(const art::Event &event)
{
    const TruthHelper::Interaction interaction(event, m_config().MCTruthLabel(), m_config().MCParticleLabel());

    m_outputEvent.m_isSignal = AnalysisHelper::IsCC1PiSignal(interaction);
    m_outputEvent.m_interaction = DebugHelper::GetInteractionString(interaction, true);

    const auto mcNeutrino = interaction.GetNeutrino();
    m_outputEvent.m_trueNuE = mcNeutrino.E();
    m_outputEvent.m_trueNuVtx = TVector3(mcNeutrino.Vx(), mcNeutrino.Vy(), mcNeutrino.Vz());
    m_outputEvent.m_isTrueNuFiducial = AnalysisHelper::IsNeutrinoVertexFiducial(interaction);

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
    m_outputEvent.m_nElectron = AnalysisHelper::CountParticlesWithPDG(mcParticles, 11);
    m_outputEvent.m_nPositron = AnalysisHelper::CountParticlesWithPDG(mcParticles, -11);
    m_outputEvent.m_nTotal = mcParticles.size();

    const auto mcParticleToHits = CollectionHelper::GetAssociationWithData<simb::MCParticle, recob::Hit, anab::BackTrackerHitMatchingData>(event, m_config().MCParticleLabel(), m_config().BacktrackerLabel());
    this->SetMCParticleInfo(interaction.GetAllMCParticles(), mcParticles, mcParticleToHits);

    if (!m_outputEvent.m_isSignal)
        return;
    
    const auto muon = this->SelectMCParticleWithPdgCode(mcParticles, 13);
    m_outputEvent.m_trueMuEnergy = muon->E();
    m_outputEvent.m_trueMuTheta = this->GetTheta(muon);
    m_outputEvent.m_trueMuPhi = this->GetPhi(muon);
    
    const auto pion = this->SelectMCParticleWithPdgCode(mcParticles, 211);
    m_outputEvent.m_truePiEnergy = pion->E();
    m_outputEvent.m_truePiTheta = this->GetTheta(pion);
    m_outputEvent.m_truePiPhi = this->GetPhi(pion);

    m_outputEvent.m_trueMuPiAngle = this->GetOpeningAngle(muon, pion);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void EventSelection::SetMCParticleInfo(const MCParticleVector &allMCParticles, const MCParticleVector &reconstrutableFinalStates, const AssociationData<simb::MCParticle, recob::Hit, anab::BackTrackerHitMatchingData> &mcParticleToHits)
{
    m_outputEvent.m_nMCParticles = allMCParticles.size();

    for (const auto &mcParticle : allMCParticles)
    {
        m_outputEvent.m_mcpIdVect.push_back(mcParticle->TrackId());
        m_outputEvent.m_mcpPDGVect.push_back(mcParticle->PdgCode());

        const auto isTargetFinalState = (std::find(reconstrutableFinalStates.begin(), reconstrutableFinalStates.end(), mcParticle) != reconstrutableFinalStates.end());
        m_outputEvent.m_mcpIsTargetFinalStateVect.push_back(isTargetFinalState);

        m_outputEvent.m_mcpProcessVect.push_back(mcParticle->Process());
        m_outputEvent.m_mcpMotherIdVect.push_back(mcParticle->Mother());
        m_outputEvent.m_mcpNDaughtersVect.push_back(mcParticle->NumberDaughters());

        std::vector<int> daughterIds;
        for (int i = 0; i < mcParticle->NumberDaughters(); ++i)
        {
            daughterIds.push_back(mcParticle->Daughter(i));
        }
        m_outputEvent.m_mcpDaughterIdsVect.push_back(daughterIds);

        m_outputEvent.m_mcpEnergyVect.push_back(mcParticle->E());
        m_outputEvent.m_mcpMomentumVect.push_back(mcParticle->P());
        m_outputEvent.m_mcpMomentumXVect.push_back(mcParticle->Px());
        m_outputEvent.m_mcpMomentumYVect.push_back(mcParticle->Py());
        m_outputEvent.m_mcpMomentumZVect.push_back(mcParticle->Pz());
    
        m_outputEvent.m_mcpNHitsUVect.push_back(BacktrackHelper::CountHitsInView(mcParticle, mcParticleToHits, geo::kU));
        m_outputEvent.m_mcpNHitsVVect.push_back(BacktrackHelper::CountHitsInView(mcParticle, mcParticleToHits, geo::kV));
        m_outputEvent.m_mcpNHitsWVect.push_back(BacktrackHelper::CountHitsInView(mcParticle, mcParticleToHits, geo::kW));
        m_outputEvent.m_mcpNGoodHitsUVect.push_back(BacktrackHelper::CountGoodHitsInView(mcParticle, mcParticleToHits, geo::kU));
        m_outputEvent.m_mcpNGoodHitsVVect.push_back(BacktrackHelper::CountGoodHitsInView(mcParticle, mcParticleToHits, geo::kV));
        m_outputEvent.m_mcpNGoodHitsWVect.push_back(BacktrackHelper::CountGoodHitsInView(mcParticle, mcParticleToHits, geo::kW));
        m_outputEvent.m_mcpHitWeightUVect.push_back(BacktrackHelper::GetHitWeightInView(mcParticle, mcParticleToHits, geo::kU));
        m_outputEvent.m_mcpHitWeightVVect.push_back(BacktrackHelper::GetHitWeightInView(mcParticle, mcParticleToHits, geo::kV));
        m_outputEvent.m_mcpHitWeightWVect.push_back(BacktrackHelper::GetHitWeightInView(mcParticle, mcParticleToHits, geo::kW));
    }
}

// -----------------------------------------------------------------------------------------------------------------------------------------


// -----------------------------------------------------------------------------------------------------------------------------------------

art::Ptr<simb::MCParticle> EventSelection::SelectMCParticleWithPdgCode(const MCParticleVector &mcParticles, const int pdgCode) const
{
    bool foundParticle = false;
    art::Ptr<simb::MCParticle> outputParticle;

    for (const auto &mcParticle : mcParticles)
    {
        if (mcParticle->PdgCode() != pdgCode)
            continue;

        if (foundParticle)
            throw cet::exception("EventSelection::SelectParticleWithPdgCode") << " - Found multiple MCParticles with PDG code: " << pdgCode << std::endl;

        foundParticle = true;
        outputParticle = mcParticle;
    }

    if (!foundParticle)
        throw cet::exception("EventSelection::SelectParticleWithPdgCode") << " - Didn't find MCParticle with PDG code: " << pdgCode << std::endl;

    return outputParticle;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

float EventSelection::GetTheta(const art::Ptr<simb::MCParticle> &mcParticle) const
{
    return (mcParticle->P() < std::numeric_limits<float>::epsilon()) ? -std::numeric_limits<float>::max() : std::acos(mcParticle->Pz() / mcParticle->P());
}

// -----------------------------------------------------------------------------------------------------------------------------------------

float EventSelection::GetPhi(const art::Ptr<simb::MCParticle> &mcParticle) const
{
    return std::atan2(mcParticle->Px(), mcParticle->Py());
}

// -----------------------------------------------------------------------------------------------------------------------------------------

float EventSelection::GetOpeningAngle(const art::Ptr<simb::MCParticle> &muon, const art::Ptr<simb::MCParticle> &pion) const
{
    const auto muDir = TVector3(muon->Px(), muon->Py(), muon->Pz()).Unit();
    const auto piDir = TVector3(pion->Px(), pion->Py(), pion->Pz()).Unit();
    return std::acos(muDir.Dot(piDir));
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void EventSelection::SetRecoInfo(const art::Event &event)
{
    // Get the PFParticles and their associations
    const auto allPFParticles = CollectionHelper::GetCollection<recob::PFParticle>(event, m_config().PFParticleLabel());
    const auto backtrackerData = AnalysisHelper::GetBacktrackerData(event, m_config().MCTruthLabel(), m_config().MCParticleLabel(), m_config().BacktrackerLabel(), m_config().PFParticleLabel());
    const auto finalStates = RecoHelper::GetNeutrinoFinalStates(allPFParticles); 
    const auto pfpToTracks = CollectionHelper::GetAssociation<recob::PFParticle, recob::Track>(event, m_config().PFParticleLabel(), m_config().TrackLabel());
    const auto trackToPIDs = CollectionHelper::GetAssociation<recob::Track, anab::ParticleID>(event, m_config().TrackLabel(), m_config().PIDLabel());
    const auto trackToCalorimetries = CollectionHelper::GetAssociation<recob::Track, anab::Calorimetry>(event, m_config().TrackLabel(), m_config().CalorimetryLabel());
    const auto pfpToMetadata = CollectionHelper::GetAssociation<recob::PFParticle, larpandoraobj::PFParticleMetadata>(event, m_config().PFParticleLabel());
    const auto pSpaceChargeService = RecoHelper::GetSpaceChargeService();

    // Set the event level reco information
    this->SetEventRecoInfo(event, allPFParticles, finalStates, pfpToMetadata, pSpaceChargeService);

    // Set the particle level reco information
    for (unsigned int i = 0; i < finalStates.size(); ++i)

    {
        const auto &finalState = finalStates.at(i);
        this->SetPFParticleInfo(i, finalState, backtrackerData, pfpToTracks, trackToPIDs, trackToCalorimetries, pfpToMetadata, pSpaceChargeService);
    }
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void EventSelection::SetEventRecoInfo(const art::Event &event, const PFParticleVector &allPFParticles, const PFParticleVector &finalStates, const Association<recob::PFParticle, larpandoraobj::PFParticleMetadata> &pfpToMetadata, const spacecharge::SpaceChargeService::provider_type *const pSpaceChargeService)
{
    art::Ptr<recob::PFParticle> neutrino;
    m_outputEvent.m_hasRecoNeutrino = false;

    try
    {
        neutrino = RecoHelper::GetNeutrino(allPFParticles);
        m_outputEvent.m_hasRecoNeutrino = true;
    }
    catch (const cet::exception &)
    {
    }

    if (m_outputEvent.m_hasRecoNeutrino)
    {
        // Get the SCE corrected reconstructed neutrino vertex
        m_outputEvent.m_recoNuVtx = RecoHelper::CorrectForSpaceCharge(RecoHelper::GetRecoNeutrinoVertex(event, allPFParticles, m_config().PFParticleLabel()), pSpaceChargeService);
        m_outputEvent.m_isRecoNuFiducial = AnalysisHelper::IsFiducial(m_outputEvent.m_recoNuVtx);
        
        // Get the topological score for the selected neutrino slice
        const auto neutrinoMetadata = CollectionHelper::GetSingleAssociated(neutrino, pfpToMetadata);
        m_outputEvent.m_topologicalScore = RecoHelper::GetTopologicalScore(neutrinoMetadata); 

        // Count the final states
        m_outputEvent.m_nFinalStatePFPs = finalStates.size();
    }
    else
    {
        m_outputEvent.m_recoNuVtx = TVector3(-std::numeric_limits<float>::max(), -std::numeric_limits<float>::max(), -std::numeric_limits<float>::max());
        m_outputEvent.m_isRecoNuFiducial = false;
        m_outputEvent.m_topologicalScore = -std::numeric_limits<float>::max();
        m_outputEvent.m_nFinalStatePFPs = 0;
    }
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void EventSelection::SetPFParticleInfo(const unsigned int index, const art::Ptr<recob::PFParticle> &finalState, const BacktrackHelper::BacktrackerData &backtrackerData, const Association<recob::PFParticle, recob::Track> &pfpToTracks, const Association<recob::Track, anab::ParticleID> &trackToPIDs, const Association<recob::Track, anab::Calorimetry> &trackToCalorimetries, const Association<recob::PFParticle, larpandoraobj::PFParticleMetadata> &pfpToMetadata, const spacecharge::SpaceChargeService::provider_type *const pSpaceChargeService)
{
    this->SetPFParticleMCParticleMatchInfo(finalState, backtrackerData);
    this->SetPFParticlePandoraInfo(finalState, backtrackerData, pfpToMetadata);

    const auto nHitsU = m_outputEvent.m_nHitsUVect.back();
    const auto nHitsV = m_outputEvent.m_nHitsVVect.back();

    // Set the track info if it's available
    if (!CollectionHelper::HasAssociated(finalState, pfpToTracks))
    {
        this->SetDummyTrackInfo();
        this->SetDummyPIDInfo();
        this->SetDummyCalorimetryInfo();
    }
    else
    {
        const auto track = CollectionHelper::GetSingleAssociated(finalState, pfpToTracks);

        this->SetTrackInfo(track, pSpaceChargeService);
        const auto yzAngle = m_outputEvent.m_yzAngleVect.back();

        // Set the Calorimetry info if it's available
        if (!CollectionHelper::HasAssociated(track, trackToCalorimetries))
        {
            this->SetDummyCalorimetryInfo();
        }
        else
        {
            const auto calos = CollectionHelper::GetManyAssociated(track, trackToCalorimetries);
            this->SetCalorimetryInfo(calos);
        }
        
        // Set the PID info if it's available
        if (!CollectionHelper::HasAssociated(track, trackToPIDs))
        {
            this->SetDummyPIDInfo();
        }
        else
        {
            const auto pid = CollectionHelper::GetSingleAssociated(track, trackToPIDs);
            this->SetPIDInfo(pid, yzAngle, nHitsU, nHitsV);
        }
    }

    this->ValidateOutputVectorSizes(index);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void EventSelection::SetPFParticleMCParticleMatchInfo(const art::Ptr<recob::PFParticle> &finalState, const BacktrackHelper::BacktrackerData &backtrackerData)
{
    bool hasMatchedMCParticle = false;
    int matchedMCParticleId = -std::numeric_limits<int>::max();
    int truePdgCode = -std::numeric_limits<int>::max();
    float truthMatchCompleteness = -std::numeric_limits<float>::max();
    float truthMatchPurity = -std::numeric_limits<float>::max();
    float trueEnergy = -std::numeric_limits<float>::max();
    float trueKE = -std::numeric_limits<float>::max();
    TVector3 trueMomentum = TVector3(-std::numeric_limits<float>::max(), -std::numeric_limits<float>::max(), -std::numeric_limits<float>::max());
    TVector3 trueStart = TVector3(-std::numeric_limits<float>::max(), -std::numeric_limits<float>::max(), -std::numeric_limits<float>::max());

    try
    {
        // ATTN this is the line that would throw if the PFParticle doesn't match to any MCParticle
        const auto matchedMCP = backtrackerData.GetBestMatchedMCParticle(finalState);

        matchedMCParticleId = matchedMCP->TrackId();
        truePdgCode = matchedMCP->PdgCode();
        truthMatchCompleteness = backtrackerData.GetMatchCompleteness(finalState, matchedMCP);
        truthMatchPurity = backtrackerData.GetMatchPurity(finalState, matchedMCP);
        trueEnergy = matchedMCP->E();
        trueKE = matchedMCP->E() - matchedMCP->Mass();
        trueMomentum = TVector3(matchedMCP->Px(), matchedMCP->Py(), matchedMCP->Pz());
        trueStart = TVector3(matchedMCP->Vx(), matchedMCP->Vy(), matchedMCP->Vz());
        hasMatchedMCParticle = true;
    }
    catch (const cet::exception &)
    {
    }

    m_outputEvent.m_hasMatchedMCParticleVect.push_back(hasMatchedMCParticle);
    m_outputEvent.m_matchedMCParticleIdVect.push_back(matchedMCParticleId);
    m_outputEvent.m_truePdgCodeVect.push_back(truePdgCode);
    m_outputEvent.m_truthMatchCompletenessVect.push_back(truthMatchCompleteness);
    m_outputEvent.m_truthMatchPurityVect.push_back(truthMatchPurity);
    m_outputEvent.m_trueEnergyVect.push_back(trueEnergy);
    m_outputEvent.m_trueKEVect.push_back(trueKE);
    m_outputEvent.m_trueMomentumXVect.push_back(trueMomentum.X());
    m_outputEvent.m_trueMomentumYVect.push_back(trueMomentum.Y());
    m_outputEvent.m_trueMomentumZVect.push_back(trueMomentum.Z());
    m_outputEvent.m_trueStartXVect.push_back(trueStart.X());
    m_outputEvent.m_trueStartYVect.push_back(trueStart.Y());
    m_outputEvent.m_trueStartZVect.push_back(trueStart.Z());
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void EventSelection::SetPFParticlePandoraInfo(const art::Ptr<recob::PFParticle> &finalState, const BacktrackHelper::BacktrackerData &backtrackerData, const Association<recob::PFParticle, larpandoraobj::PFParticleMetadata> &pfpToMetadata)
{
    // ATTN here we are just using the backtracker data for convenience, for data we need to go manually go through the clusters
    const auto &hits = backtrackerData.GetHits(finalState);

    m_outputEvent.m_nHitsUVect.push_back(RecoHelper::CountHitsInView(hits, geo::kU));
    m_outputEvent.m_nHitsVVect.push_back(RecoHelper::CountHitsInView(hits, geo::kV));
    m_outputEvent.m_nHitsWVect.push_back(RecoHelper::CountHitsInView(hits, geo::kW));
    
    const auto metadata = CollectionHelper::GetSingleAssociated(finalState, pfpToMetadata);
    m_outputEvent.m_trackShowerVect.push_back(RecoHelper::GetTrackScore(metadata));
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void EventSelection::SetDummyTrackInfo()
{
    m_outputEvent.m_hasTrackInfoVect.push_back(false);
    m_outputEvent.m_startXVect.push_back(-std::numeric_limits<float>::max());
    m_outputEvent.m_startYVect.push_back(-std::numeric_limits<float>::max());
    m_outputEvent.m_startZVect.push_back(-std::numeric_limits<float>::max());
    m_outputEvent.m_endXVect.push_back(-std::numeric_limits<float>::max());
    m_outputEvent.m_endYVect.push_back(-std::numeric_limits<float>::max());
    m_outputEvent.m_endZVect.push_back(-std::numeric_limits<float>::max());
    m_outputEvent.m_directionXVect.push_back(-std::numeric_limits<float>::max());
    m_outputEvent.m_directionYVect.push_back(-std::numeric_limits<float>::max());
    m_outputEvent.m_directionZVect.push_back(-std::numeric_limits<float>::max());
    m_outputEvent.m_thetaVect.push_back(-std::numeric_limits<float>::max());
    m_outputEvent.m_phiVect.push_back(-std::numeric_limits<float>::max());
    m_outputEvent.m_yzAngleVect.push_back(-std::numeric_limits<float>::max());
    m_outputEvent.m_lengthVect.push_back(-std::numeric_limits<float>::max());
    m_outputEvent.m_isContainedVect.push_back(false);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void EventSelection::SetTrackInfo(const art::Ptr<recob::Track> &track, const spacecharge::SpaceChargeService::provider_type *const pSpaceChargeService)
{
    const auto start = RecoHelper::CorrectForSpaceCharge(TVector3(track->Start().X(), track->Start().Y(), track->Start().Z()), pSpaceChargeService);
    const auto end = RecoHelper::CorrectForSpaceCharge(TVector3(track->End().X(), track->End().Y(), track->End().Z()), pSpaceChargeService);
    const auto dir = TVector3(track->StartDirection().X(), track->StartDirection().Y(), track->StartDirection().Z());

    m_outputEvent.m_hasTrackInfoVect.push_back(true);
    m_outputEvent.m_startXVect.push_back(start.X());
    m_outputEvent.m_startYVect.push_back(start.Y());
    m_outputEvent.m_startZVect.push_back(start.Z());
    m_outputEvent.m_endXVect.push_back(end.X());
    m_outputEvent.m_endYVect.push_back(end.Y());
    m_outputEvent.m_endZVect.push_back(end.Z());
    m_outputEvent.m_directionXVect.push_back(dir.X());
    m_outputEvent.m_directionYVect.push_back(dir.Y());
    m_outputEvent.m_directionZVect.push_back(dir.Z());
    m_outputEvent.m_thetaVect.push_back(track->Theta());
    m_outputEvent.m_phiVect.push_back(track->Phi());
    m_outputEvent.m_yzAngleVect.push_back(this->GetYZAngle(dir));
    m_outputEvent.m_lengthVect.push_back(track->Length());
    m_outputEvent.m_isContainedVect.push_back(AnalysisHelper::IsContained(start, m_config().ContainmentBorder()) && AnalysisHelper::IsContained(end, m_config().ContainmentBorder()));
}

// -----------------------------------------------------------------------------------------------------------------------------------------

float EventSelection::GetYZAngle(const TVector3 &dir)
{
    return atan2(dir.Z(), dir.Y());
}

// -----------------------------------------------------------------------------------------------------------------------------------------
    
void EventSelection::SetDummyCalorimetryInfo()
{
    const std::vector<float> dummy;
    m_outputEvent.m_hasCalorimetryInfoVect.push_back(false);
    m_outputEvent.m_dedxPerHitUVect.push_back(dummy);
    m_outputEvent.m_dedxPerHitVVect.push_back(dummy);
    m_outputEvent.m_dedxPerHitWVect.push_back(dummy);
    m_outputEvent.m_residualRangePerHitUVect.push_back(dummy);
    m_outputEvent.m_residualRangePerHitVVect.push_back(dummy);
    m_outputEvent.m_residualRangePerHitWVect.push_back(dummy);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void EventSelection::SetDummyPIDInfo()
{
    m_outputEvent.m_hasPIDInfoVect.push_back(false);

    m_outputEvent.m_isWChi2pAvailableVect.push_back(false);
    m_outputEvent.m_chi2pWVect.push_back(-std::numeric_limits<float>::max());
    m_outputEvent.m_isUChi2pAvailableVect.push_back(false);
    m_outputEvent.m_chi2pUVect.push_back(-std::numeric_limits<float>::max());
    m_outputEvent.m_isVChi2pAvailableVect.push_back(false);
    m_outputEvent.m_chi2pVVect.push_back(-std::numeric_limits<float>::max());
    m_outputEvent.m_isUVChi2pAvailableVect.push_back(false);
    m_outputEvent.m_chi2pUVVect.push_back(-std::numeric_limits<float>::max());

    m_outputEvent.m_isWBraggpAvailableVect.push_back(false);
    m_outputEvent.m_braggpWVect.push_back(-std::numeric_limits<float>::max());
    m_outputEvent.m_isUBraggpAvailableVect.push_back(false);
    m_outputEvent.m_braggpUVect.push_back(-std::numeric_limits<float>::max());
    m_outputEvent.m_isVBraggpAvailableVect.push_back(false);
    m_outputEvent.m_braggpVVect.push_back(-std::numeric_limits<float>::max());
    m_outputEvent.m_isUVBraggpAvailableVect.push_back(false);
    m_outputEvent.m_braggpUVVect.push_back(-std::numeric_limits<float>::max());
    
    m_outputEvent.m_isWBraggpBackwardAvailableVect.push_back(false);
    m_outputEvent.m_braggpBackwardWVect.push_back(-std::numeric_limits<float>::max());
    m_outputEvent.m_isUBraggpBackwardAvailableVect.push_back(false);
    m_outputEvent.m_braggpBackwardUVect.push_back(-std::numeric_limits<float>::max());
    m_outputEvent.m_isVBraggpBackwardAvailableVect.push_back(false);
    m_outputEvent.m_braggpBackwardVVect.push_back(-std::numeric_limits<float>::max());
    m_outputEvent.m_isUVBraggpBackwardAvailableVect.push_back(false);
    m_outputEvent.m_braggpBackwardUVVect.push_back(-std::numeric_limits<float>::max());

    m_outputEvent.m_isWBraggMIPAvailableVect.push_back(false);
    m_outputEvent.m_braggMIPWVect.push_back(-std::numeric_limits<float>::max());
    m_outputEvent.m_isUBraggMIPAvailableVect.push_back(false);
    m_outputEvent.m_braggMIPUVect.push_back(-std::numeric_limits<float>::max());
    m_outputEvent.m_isVBraggMIPAvailableVect.push_back(false);
    m_outputEvent.m_braggMIPVVect.push_back(-std::numeric_limits<float>::max());
    m_outputEvent.m_isUVBraggMIPAvailableVect.push_back(false);
    m_outputEvent.m_braggMIPUVVect.push_back(-std::numeric_limits<float>::max());

    m_outputEvent.m_braggMIPBackwardWVect.push_back(-std::numeric_limits<float>::max());
    m_outputEvent.m_isUBraggMIPBackwardAvailableVect.push_back(false);
    m_outputEvent.m_braggMIPBackwardUVect.push_back(-std::numeric_limits<float>::max());
    m_outputEvent.m_isVBraggMIPBackwardAvailableVect.push_back(false);
    m_outputEvent.m_braggMIPBackwardVVect.push_back(-std::numeric_limits<float>::max());
    m_outputEvent.m_isUVBraggMIPBackwardAvailableVect.push_back(false);
    m_outputEvent.m_braggMIPBackwardUVVect.push_back(-std::numeric_limits<float>::max());
    
    m_outputEvent.m_isWBraggRatioAvailableVect.push_back(false);
    m_outputEvent.m_braggRatioWVect.push_back(-std::numeric_limits<float>::max());
    m_outputEvent.m_isUVBraggRatioAvailableVect.push_back(false);
    m_outputEvent.m_braggRatioUVVect.push_back(-std::numeric_limits<float>::max());
}

// -----------------------------------------------------------------------------------------------------------------------------------------
    
void EventSelection::SetCalorimetryInfo(const CalorimetryVector &calos)
{
    if (calos.empty())
        throw cet::exception("EventSelection::SetCalorimetryInfo") << " - Input vector of calorimetry objects is empty!" << std::endl;

    std::vector<float> dedxPerHitU, dedxPerHitV, dedxPerHitW;
    std::vector<float> residualRangePerHitU, residualRangePerHitV, residualRangePerHitW;

    for (const auto &calo : calos)
    {
        switch (calo->PlaneID().Plane)
        {
            case 0:
                if (!dedxPerHitU.empty() || !residualRangePerHitU.empty())
                    throw cet::exception("EventSelection::SetCalorimetryInfo") << " - Multiple input calorimetry objects in plane 0 (U)" << std::endl;

                dedxPerHitU = calo->dEdx();
                residualRangePerHitU = calo->ResidualRange();

                break;
            case 1:
                if (!dedxPerHitV.empty() || !residualRangePerHitV.empty())
                    throw cet::exception("EventSelection::SetCalorimetryInfo") << " - Multiple input calorimetry objects in plane 1 (V)" << std::endl;
                
                dedxPerHitV = calo->dEdx();
                residualRangePerHitV = calo->ResidualRange();

                break;
            case 2:
                if (!dedxPerHitW.empty() || !residualRangePerHitW.empty())
                    throw cet::exception("EventSelection::SetCalorimetryInfo") << " - Multiple input calorimetry objects in plane 2 (W)" << std::endl;
                
                dedxPerHitW = calo->dEdx();
                residualRangePerHitW = calo->ResidualRange();

                break;
            default: break;
        }
    }

    if (dedxPerHitU.size() != residualRangePerHitU.size() || dedxPerHitV.size() != residualRangePerHitV.size() || dedxPerHitW.size() != residualRangePerHitW.size())
        throw cet::exception("EventSelection::SetCalorimetryInfo") << " - Invalid calorimetry object, different number of dedx and residual range points" << std::endl;

    m_outputEvent.m_hasCalorimetryInfoVect.push_back(true);
    m_outputEvent.m_dedxPerHitUVect.push_back(dedxPerHitU);
    m_outputEvent.m_dedxPerHitVVect.push_back(dedxPerHitV);
    m_outputEvent.m_dedxPerHitWVect.push_back(dedxPerHitW);
    m_outputEvent.m_residualRangePerHitUVect.push_back(residualRangePerHitU);
    m_outputEvent.m_residualRangePerHitVVect.push_back(residualRangePerHitV);
    m_outputEvent.m_residualRangePerHitWVect.push_back(residualRangePerHitW);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void EventSelection::SetPIDInfo(const art::Ptr<anab::ParticleID> &pid, const float yzAngle, const int nHitsU, const int nHitsV)
{
    // Set the defaults for the raw PID variables
    float chi2pW = -std::numeric_limits<float>::max();
    float chi2pU = -std::numeric_limits<float>::max();
    float chi2pV = -std::numeric_limits<float>::max();
    float chi2pUV = -std::numeric_limits<float>::max();

    float braggpW = -std::numeric_limits<float>::max();
    float braggpU = -std::numeric_limits<float>::max();
    float braggpV = -std::numeric_limits<float>::max();
    float braggpUV = -std::numeric_limits<float>::max();
    
    float braggpBackwardW = -std::numeric_limits<float>::max();
    float braggpBackwardU = -std::numeric_limits<float>::max();
    float braggpBackwardV = -std::numeric_limits<float>::max();
    float braggpBackwardUV = -std::numeric_limits<float>::max();

    float braggMIPW = -std::numeric_limits<float>::max();
    float braggMIPU = -std::numeric_limits<float>::max();
    float braggMIPV = -std::numeric_limits<float>::max();
    float braggMIPUV = -std::numeric_limits<float>::max();

    float braggMIPBackwardW = -std::numeric_limits<float>::max();
    float braggMIPBackwardU = -std::numeric_limits<float>::max();
    float braggMIPBackwardV = -std::numeric_limits<float>::max();
    float braggMIPBackwardUV = -std::numeric_limits<float>::max();
    
    bool isWChi2pAvailable = false;
    bool isUChi2pAvailable = false;
    bool isVChi2pAvailable = false;
    bool isUVChi2pAvailable = false;

    bool isWBraggpAvailable = false;
    bool isUBraggpAvailable = false;
    bool isVBraggpAvailable = false;
    bool isUVBraggpAvailable = false;
    
    bool isWBraggpBackwardAvailable = false;
    bool isUBraggpBackwardAvailable = false;
    bool isVBraggpBackwardAvailable = false;
    bool isUVBraggpBackwardAvailable = false;

    bool isWBraggMIPAvailable = false;
    bool isUBraggMIPAvailable = false;
    bool isVBraggMIPAvailable = false;
    bool isUVBraggMIPAvailable = false;

    bool isWBraggMIPBackwardAvailable = false;
    bool isUBraggMIPBackwardAvailable = false;
    bool isVBraggMIPBackwardAvailable = false;
    bool isUVBraggMIPBackwardAvailable = false;

    // Extract the raw PID values
    for (const auto &algo : pid->fParticleIDAlgScores)
    {
        const auto view = this->GetView(algo.fPlaneMask);

        // Chi2 algorithm under the proton hypothesis
        if (algo.fAlgName == "Chi2" && algo.fAssumedPdg == 2212)
            this->SetPIDVariables(view, algo.fValue, chi2pW, chi2pU, chi2pV, isWChi2pAvailable, isUChi2pAvailable, isVChi2pAvailable);

        // Bragg Likelihood algorithm under the proton hypothesis - forward going track
        if (algo.fAlgName == "BraggPeakLLH" && algo.fTrackDir == anab::kForward && algo.fAssumedPdg == 2212)
            this->SetPIDVariables(view, algo.fValue, braggpW, braggpU, braggpV, isWBraggpAvailable, isUBraggpAvailable, isVBraggpAvailable);
        
        // Bragg Likelihood algorithm under the proton hypothesis - backward going track
        if (algo.fAlgName == "BraggPeakLLH" && algo.fTrackDir == anab::kBackward && algo.fAssumedPdg == 2212)
            this->SetPIDVariables(view, algo.fValue, braggpBackwardW, braggpBackwardU, braggpBackwardV, isWBraggpBackwardAvailable, isUBraggpBackwardAvailable, isVBraggpBackwardAvailable);
        
        // Bragg Likelihood algorithm under the MIP hypothesis - forward going track
        if (algo.fAlgName == "BraggPeakLLH" && algo.fTrackDir == anab::kForward && algo.fAssumedPdg == 0)
            this->SetPIDVariables(view, algo.fValue, braggMIPW, braggMIPU, braggMIPV, isWBraggMIPAvailable, isUBraggMIPAvailable, isVBraggMIPAvailable);
        
        // Bragg Likelihood algorithm under the MIP hypothesis - backward going track
        if (algo.fAlgName == "BraggPeakLLH" && algo.fTrackDir == anab::kBackward && algo.fAssumedPdg == 0)
            this->SetPIDVariables(view, algo.fValue, braggMIPBackwardW, braggMIPBackwardU, braggMIPBackwardV, isWBraggMIPBackwardAvailable, isUBraggMIPBackwardAvailable, isVBraggMIPBackwardAvailable);
    }

    // Set the combined induction plane variables
    this->CombineInductionPlanes(yzAngle, nHitsU, nHitsV, chi2pU, chi2pV, isUChi2pAvailable, isVChi2pAvailable, chi2pUV, isUVChi2pAvailable);
    this->CombineInductionPlanes(yzAngle, nHitsU, nHitsV, braggpU, braggpV, isUBraggpAvailable, isVBraggpAvailable, braggpUV, isUVBraggpAvailable);
    this->CombineInductionPlanes(yzAngle, nHitsU, nHitsV, braggpBackwardU, braggpBackwardV, isUBraggpBackwardAvailable, isVBraggpBackwardAvailable, braggpBackwardUV, isUVBraggpBackwardAvailable);
    this->CombineInductionPlanes(yzAngle, nHitsU, nHitsV, braggMIPU, braggMIPV, isUBraggMIPAvailable, isVBraggMIPAvailable, braggMIPUV, isUVBraggMIPAvailable);
    this->CombineInductionPlanes(yzAngle, nHitsU, nHitsV, braggMIPBackwardU, braggMIPBackwardV, isUBraggMIPBackwardAvailable, isVBraggMIPBackwardAvailable, braggMIPBackwardUV, isUVBraggMIPBackwardAvailable);

    // Calculate the bragg peak ratio
    const bool isWBraggRatioAvailable = (isWBraggpAvailable && isWBraggMIPAvailable && braggMIPW >= std::numeric_limits<float>::epsilon());
    const bool isUVBraggRatioAvailable = (isUVBraggpAvailable && isUVBraggMIPAvailable && braggMIPUV >= std::numeric_limits<float>::epsilon());
    const float braggRatioW = isWBraggRatioAvailable ? (braggpW / braggMIPW) : -std::numeric_limits<float>::max();
    const float braggRatioUV = isUVBraggRatioAvailable ? (braggpUV / braggMIPUV) : -std::numeric_limits<float>::max();

    // Fill the vectors
    m_outputEvent.m_hasPIDInfoVect.push_back(true);

    m_outputEvent.m_isWChi2pAvailableVect.push_back(isWChi2pAvailable);
    m_outputEvent.m_chi2pWVect.push_back(chi2pW);
    m_outputEvent.m_isUChi2pAvailableVect.push_back(isUChi2pAvailable);
    m_outputEvent.m_chi2pUVect.push_back(chi2pU);
    m_outputEvent.m_isVChi2pAvailableVect.push_back(isVChi2pAvailable);
    m_outputEvent.m_chi2pVVect.push_back(chi2pV);
    m_outputEvent.m_isUVChi2pAvailableVect.push_back(isUVChi2pAvailable);
    m_outputEvent.m_chi2pUVVect.push_back(chi2pUV);

    m_outputEvent.m_isWBraggpAvailableVect.push_back(isWBraggpAvailable);
    m_outputEvent.m_braggpWVect.push_back(braggpW);
    m_outputEvent.m_isUBraggpAvailableVect.push_back(isUBraggpAvailable);
    m_outputEvent.m_braggpUVect.push_back(braggpU);
    m_outputEvent.m_isVBraggpAvailableVect.push_back(isVBraggpAvailable);
    m_outputEvent.m_braggpVVect.push_back(braggpV);
    m_outputEvent.m_isUVBraggpAvailableVect.push_back(isUVBraggpAvailable);
    m_outputEvent.m_braggpUVVect.push_back(braggpUV);
    
    m_outputEvent.m_isWBraggpBackwardAvailableVect.push_back(isWBraggpBackwardAvailable);
    m_outputEvent.m_braggpBackwardWVect.push_back(braggpBackwardW);
    m_outputEvent.m_isUBraggpBackwardAvailableVect.push_back(isUBraggpBackwardAvailable);
    m_outputEvent.m_braggpBackwardUVect.push_back(braggpBackwardU);
    m_outputEvent.m_isVBraggpBackwardAvailableVect.push_back(isVBraggpBackwardAvailable);
    m_outputEvent.m_braggpBackwardVVect.push_back(braggpBackwardV);
    m_outputEvent.m_isUVBraggpBackwardAvailableVect.push_back(isUVBraggpBackwardAvailable);
    m_outputEvent.m_braggpBackwardUVVect.push_back(braggpBackwardUV);

    m_outputEvent.m_isWBraggMIPAvailableVect.push_back(isWBraggMIPAvailable);
    m_outputEvent.m_braggMIPWVect.push_back(braggMIPW);
    m_outputEvent.m_isUBraggMIPAvailableVect.push_back(isUBraggMIPAvailable);
    m_outputEvent.m_braggMIPUVect.push_back(braggMIPU);
    m_outputEvent.m_isVBraggMIPAvailableVect.push_back(isVBraggMIPAvailable);
    m_outputEvent.m_braggMIPVVect.push_back(braggMIPV);
    m_outputEvent.m_isUVBraggMIPAvailableVect.push_back(isUVBraggMIPAvailable);
    m_outputEvent.m_braggMIPUVVect.push_back(braggMIPUV);

    m_outputEvent.m_isWBraggMIPBackwardAvailableVect.push_back(isWBraggMIPBackwardAvailable);
    m_outputEvent.m_braggMIPBackwardWVect.push_back(braggMIPBackwardW);
    m_outputEvent.m_isUBraggMIPBackwardAvailableVect.push_back(isUBraggMIPBackwardAvailable);
    m_outputEvent.m_braggMIPBackwardUVect.push_back(braggMIPBackwardU);
    m_outputEvent.m_isVBraggMIPBackwardAvailableVect.push_back(isVBraggMIPBackwardAvailable);
    m_outputEvent.m_braggMIPBackwardVVect.push_back(braggMIPBackwardV);
    m_outputEvent.m_isUVBraggMIPBackwardAvailableVect.push_back(isUVBraggMIPBackwardAvailable);
    m_outputEvent.m_braggMIPBackwardUVVect.push_back(braggMIPBackwardUV);
    
    m_outputEvent.m_isWBraggRatioAvailableVect.push_back(isWBraggRatioAvailable);
    m_outputEvent.m_braggRatioWVect.push_back(braggRatioW);
    m_outputEvent.m_isUVBraggRatioAvailableVect.push_back(isUVBraggRatioAvailable);
    m_outputEvent.m_braggRatioUVVect.push_back(braggRatioUV);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

geo::View_t EventSelection::GetView(const std::bitset<8> &planeMask) const
{
    // Here is a hack to get around a bug in the PID code. Different algorithms have different numbering conventions
    const bool usesW = planeMask.test(2) || planeMask.test(7);
    const bool usesU = planeMask.test(1) || planeMask.test(6);
    const bool usesV = planeMask.test(0) || planeMask.test(5);
    
    if (usesW && !usesU && !usesV)
        return geo::kW;
    
    if (!usesW && usesU && !usesV)
        return geo::kU;
    
    if (!usesW && !usesU && usesV)
        return geo::kV;

    return geo::kUnknown;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void EventSelection::SetPIDVariables(const geo::View_t &view, const float algoValue, float &wValue, float &uValue, float &vValue, bool &wSet, bool &uSet, bool &vSet) const
{
    switch (view)
    {
        case geo::kW:
            wValue = algoValue;
            wSet = true;
            break;
        case geo::kU:
            uValue = algoValue;
            uSet = true;
            break;
        case geo::kV:
            vValue = algoValue;
            vSet = true;
            break;
        default:
            throw cet::exception("EventSelection::SetPIDInfo") << " - Don't know how to interpret plane mask!" << std::endl;
    }
}

// -----------------------------------------------------------------------------------------------------------------------------------------
    
void EventSelection::CombineInductionPlanes(const float yzAngle, const int nHitsU, const int nHitsV, const float uValue, const float vValue, const bool uSet, const bool vSet, float &uvValue, bool &uvSet) const
{
    uvValue = -std::numeric_limits<float>::max();

    const float uAngle = -acos(0.5f);
    const float vAngle = +acos(0.5f);

    const bool uAngleCut = pow(sin(yzAngle - uAngle), 2) > m_config().Sin2YZAngleCut();
    const bool vAngleCut = pow(sin(yzAngle - vAngle), 2) > m_config().Sin2YZAngleCut();
    
    // Check if we have any information available
    const bool uAvailable = uSet && uAngleCut;
    const bool vAvailable = vSet && vAngleCut;

    uvSet = (uAvailable || vAvailable) && (nHitsU > 0 || nHitsV > 0);
    if (!uvSet)
        return;

    const float uWeight = static_cast<float>(uAvailable ? nHitsU : 0);
    const float vWeight = static_cast<float>(vAvailable ? nHitsV : 0);
    const float totalWeight = uWeight + vWeight;

    // ATTN here total weight must be non-zero due to nHits check above
    uvValue = (uWeight*uValue + vWeight*vValue) / totalWeight;
}

// -----------------------------------------------------------------------------------------------------------------------------------------
    
void EventSelection::ValidateOutputVectorSizes(const unsigned int index) const
{
    // At this point, we expect that every vector should be of the same size, and that size should be index + 1 (start counting from 0)
    // This function isn't strictly necessary, but it sanity checks that the vectors aren't going out of alignment due to any errors in the 
    // filling logic which could lead to a nightmare downstream!

    const auto expectedSize = index + 1;

    if (m_outputEvent.m_nFinalStatePFPs < 0)
        throw cet::exception("EventSelection::ValidateOutputVectorSizes") << " - Number of final state PFParticles is negative: " << m_outputEvent.m_nFinalStatePFPs << std::endl;

    if (static_cast<unsigned int>(m_outputEvent.m_nFinalStatePFPs) < expectedSize)
        throw cet::exception("EventSelection::ValidateOutputVectorSizes") << " - Particle index out of range: " << expectedSize << " / " << m_outputEvent.m_nFinalStatePFPs << std::endl;

    if (m_outputEvent.m_hasMatchedMCParticleVect.size() != expectedSize ||
        m_outputEvent.m_matchedMCParticleIdVect.size() != expectedSize ||
        m_outputEvent.m_truePdgCodeVect.size() != expectedSize ||
        m_outputEvent.m_truthMatchCompletenessVect.size() != expectedSize ||
        m_outputEvent.m_truthMatchPurityVect.size() != expectedSize ||
        m_outputEvent.m_trueEnergyVect.size() != expectedSize ||
        m_outputEvent.m_trueKEVect.size() != expectedSize ||
        m_outputEvent.m_trueMomentumXVect.size() != expectedSize ||
        m_outputEvent.m_trueMomentumYVect.size() != expectedSize ||
        m_outputEvent.m_trueMomentumZVect.size() != expectedSize ||
        m_outputEvent.m_trueStartXVect.size() != expectedSize ||
        m_outputEvent.m_trueStartYVect.size() != expectedSize ||
        m_outputEvent.m_trueStartZVect.size() != expectedSize ||
        m_outputEvent.m_nHitsUVect.size() != expectedSize ||
        m_outputEvent.m_nHitsVVect.size() != expectedSize ||
        m_outputEvent.m_nHitsWVect.size() != expectedSize ||
        m_outputEvent.m_trackShowerVect.size() != expectedSize ||
        m_outputEvent.m_startXVect.size() != expectedSize ||
        m_outputEvent.m_startYVect.size() != expectedSize ||
        m_outputEvent.m_startZVect.size() != expectedSize ||
        m_outputEvent.m_endXVect.size() != expectedSize ||
        m_outputEvent.m_endYVect.size() != expectedSize ||
        m_outputEvent.m_endZVect.size() != expectedSize ||
        m_outputEvent.m_directionXVect.size() != expectedSize ||
        m_outputEvent.m_directionYVect.size() != expectedSize ||
        m_outputEvent.m_directionZVect.size() != expectedSize ||
        m_outputEvent.m_thetaVect.size() != expectedSize ||
        m_outputEvent.m_phiVect.size() != expectedSize ||
        m_outputEvent.m_yzAngleVect.size() != expectedSize ||
        m_outputEvent.m_lengthVect.size() != expectedSize ||
        m_outputEvent.m_isContainedVect.size() != expectedSize ||
        m_outputEvent.m_hasCalorimetryInfoVect.size() != expectedSize ||
        m_outputEvent.m_dedxPerHitUVect.size() != expectedSize ||
        m_outputEvent.m_dedxPerHitVVect.size() != expectedSize ||
        m_outputEvent.m_dedxPerHitWVect.size() != expectedSize ||
        m_outputEvent.m_residualRangePerHitUVect.size() != expectedSize ||
        m_outputEvent.m_residualRangePerHitVVect.size() != expectedSize ||
        m_outputEvent.m_residualRangePerHitWVect.size() != expectedSize ||
        m_outputEvent.m_hasPIDInfoVect.size() != expectedSize ||
        m_outputEvent.m_isWChi2pAvailableVect.size() != expectedSize ||
        m_outputEvent.m_chi2pWVect.size() != expectedSize ||
        m_outputEvent.m_isUChi2pAvailableVect.size() != expectedSize ||
        m_outputEvent.m_chi2pUVect.size() != expectedSize ||
        m_outputEvent.m_isVChi2pAvailableVect.size() != expectedSize ||
        m_outputEvent.m_chi2pVVect.size() != expectedSize ||
        m_outputEvent.m_isUVChi2pAvailableVect.size() != expectedSize ||
        m_outputEvent.m_chi2pUVVect.size() != expectedSize ||
        m_outputEvent.m_isWBraggpAvailableVect.size() != expectedSize ||
        m_outputEvent.m_braggpWVect.size() != expectedSize ||
        m_outputEvent.m_isUBraggpAvailableVect.size() != expectedSize ||
        m_outputEvent.m_braggpUVect.size() != expectedSize ||
        m_outputEvent.m_isVBraggpAvailableVect.size() != expectedSize ||
        m_outputEvent.m_braggpVVect.size() != expectedSize ||
        m_outputEvent.m_isUVBraggpAvailableVect.size() != expectedSize ||
        m_outputEvent.m_braggpUVVect.size() != expectedSize ||
        m_outputEvent.m_isWBraggpBackwardAvailableVect.size() != expectedSize ||
        m_outputEvent.m_braggpBackwardWVect.size() != expectedSize ||
        m_outputEvent.m_isUBraggpBackwardAvailableVect.size() != expectedSize ||
        m_outputEvent.m_braggpBackwardUVect.size() != expectedSize ||
        m_outputEvent.m_isVBraggpBackwardAvailableVect.size() != expectedSize ||
        m_outputEvent.m_braggpBackwardVVect.size() != expectedSize ||
        m_outputEvent.m_isUVBraggpBackwardAvailableVect.size() != expectedSize ||
        m_outputEvent.m_braggpBackwardUVVect.size() != expectedSize ||
        m_outputEvent.m_isWBraggMIPAvailableVect.size() != expectedSize ||
        m_outputEvent.m_braggMIPWVect.size() != expectedSize ||
        m_outputEvent.m_isUBraggMIPAvailableVect.size() != expectedSize ||
        m_outputEvent.m_braggMIPUVect.size() != expectedSize ||
        m_outputEvent.m_isVBraggMIPAvailableVect.size() != expectedSize ||
        m_outputEvent.m_braggMIPVVect.size() != expectedSize ||
        m_outputEvent.m_isUVBraggMIPAvailableVect.size() != expectedSize ||
        m_outputEvent.m_braggMIPUVVect.size() != expectedSize ||
        m_outputEvent.m_isWBraggMIPBackwardAvailableVect.size() != expectedSize ||
        m_outputEvent.m_braggMIPBackwardWVect.size() != expectedSize ||
        m_outputEvent.m_isUBraggMIPBackwardAvailableVect.size() != expectedSize ||
        m_outputEvent.m_braggMIPBackwardUVect.size() != expectedSize ||
        m_outputEvent.m_isVBraggMIPBackwardAvailableVect.size() != expectedSize ||
        m_outputEvent.m_braggMIPBackwardVVect.size() != expectedSize ||
        m_outputEvent.m_isUVBraggMIPBackwardAvailableVect.size() != expectedSize ||
        m_outputEvent.m_braggMIPBackwardUVVect.size() != expectedSize ||
        m_outputEvent.m_isWBraggRatioAvailableVect.size() != expectedSize ||
        m_outputEvent.m_braggRatioWVect.size() != expectedSize ||
        m_outputEvent.m_isUVBraggRatioAvailableVect.size() != expectedSize ||
        m_outputEvent.m_braggRatioUVVect.size() != expectedSize)
    {
        throw cet::exception("EventSelection::ValidateOutputVectorSizes") << " - Output particle vectors are out of sync!" << std::endl;
    }
}

} // namespace ubcc1pi
