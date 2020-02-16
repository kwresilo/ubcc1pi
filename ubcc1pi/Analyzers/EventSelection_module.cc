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
//    m_pEventTree->Branch("topologicalScore", &m_outputEvent.m_topologicalScore);
    m_pEventTree->Branch("nFinalStatePFPs", &m_outputEvent.m_nFinalStatePFPs);

    m_pEventTree->Branch("hasMatchedMCParticleVect", &m_outputEvent.m_hasMatchedMCParticleVect);
    m_pEventTree->Branch("matchedMCParticleIdVect", &m_outputEvent.m_matchedMCParticleIdVect);
    m_pEventTree->Branch("truePdgCodeVect", &m_outputEvent.m_truePdgCodeVect);
    m_pEventTree->Branch("truthMatchCompletenessVect", &m_outputEvent.m_truthMatchCompletenessVect);
    m_pEventTree->Branch("truthMatchPurityVect", &m_outputEvent.m_truthMatchPurityVect);
    m_pEventTree->Branch("trueEnergyVect", &m_outputEvent.m_trueEnergyVect);
    m_pEventTree->Branch("trueKEVect", &m_outputEvent.m_trueKEVect);
    m_pEventTree->Branch("trueRangeVect", &m_outputEvent.m_trueRangeVect);
    m_pEventTree->Branch("trueMomentumXVect", &m_outputEvent.m_trueMomentumXVect);
    m_pEventTree->Branch("trueMomentumYVect", &m_outputEvent.m_trueMomentumYVect);
    m_pEventTree->Branch("trueMomentumZVect", &m_outputEvent.m_trueMomentumZVect);
    m_pEventTree->Branch("trueStartXVect", &m_outputEvent.m_trueStartXVect);
    m_pEventTree->Branch("trueStartYVect", &m_outputEvent.m_trueStartYVect);
    m_pEventTree->Branch("trueStartZVect", &m_outputEvent.m_trueStartZVect);
    m_pEventTree->Branch("trueIsContainedVect", &m_outputEvent.m_trueIsContainedVect);
    m_pEventTree->Branch("trueIsStoppingVect", &m_outputEvent.m_trueIsStoppingVect);
    m_pEventTree->Branch("trueNScattersVect", &m_outputEvent.m_trueNScattersVect);
    m_pEventTree->Branch("trueIsGoldenVect", &m_outputEvent.m_trueIsGoldenVect);
    
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
    m_pEventTree->Branch("xyAngleVect", &m_outputEvent.m_xyAngleVect);
    m_pEventTree->Branch("xzAngleVect", &m_outputEvent.m_xzAngleVect);
    m_pEventTree->Branch("lengthVect", &m_outputEvent.m_lengthVect);
    m_pEventTree->Branch("rangeVect", &m_outputEvent.m_rangeVect);
    m_pEventTree->Branch("isContainedVect", &m_outputEvent.m_isContainedVect);
    
    m_pEventTree->Branch("nHitsUVect", &m_outputEvent.m_nHitsUVect);
    m_pEventTree->Branch("nHitsVVect", &m_outputEvent.m_nHitsVVect);
    m_pEventTree->Branch("nHitsWVect", &m_outputEvent.m_nHitsWVect);
    m_pEventTree->Branch("nDescendentHitsUVect", &m_outputEvent.m_nDescendentHitsUVect);
    m_pEventTree->Branch("nDescendentHitsVVect", &m_outputEvent.m_nDescendentHitsVVect);
    m_pEventTree->Branch("nDescendentHitsWVect", &m_outputEvent.m_nDescendentHitsWVect);
    m_pEventTree->Branch("trackShowerVect", &m_outputEvent.m_trackShowerVect);
    m_pEventTree->Branch("nDaughtersVect", &m_outputEvent.m_nDaughtersVect);
    m_pEventTree->Branch("nDescendentsVect", &m_outputEvent.m_nDescendentsVect);
            
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
    
    m_pEventTree->Branch("isWChi2muAvailableVect", &m_outputEvent.m_isWChi2muAvailableVect);
    m_pEventTree->Branch("chi2muWVect", &m_outputEvent.m_chi2muWVect);
    m_pEventTree->Branch("isUChi2muAvailableVect", &m_outputEvent.m_isUChi2muAvailableVect);
    m_pEventTree->Branch("chi2muUVect", &m_outputEvent.m_chi2muUVect);
    m_pEventTree->Branch("isVChi2muAvailableVect", &m_outputEvent.m_isVChi2muAvailableVect);
    m_pEventTree->Branch("chi2muVVect", &m_outputEvent.m_chi2muVVect);
    m_pEventTree->Branch("isUVChi2muAvailableVect", &m_outputEvent.m_isUVChi2muAvailableVect);
    m_pEventTree->Branch("chi2muUVVect", &m_outputEvent.m_chi2muUVVect);
    
    m_pEventTree->Branch("isWChi2piAvailableVect", &m_outputEvent.m_isWChi2piAvailableVect);
    m_pEventTree->Branch("chi2piWVect", &m_outputEvent.m_chi2piWVect);
    m_pEventTree->Branch("isUChi2piAvailableVect", &m_outputEvent.m_isUChi2piAvailableVect);
    m_pEventTree->Branch("chi2piUVect", &m_outputEvent.m_chi2piUVect);
    m_pEventTree->Branch("isVChi2piAvailableVect", &m_outputEvent.m_isVChi2piAvailableVect);
    m_pEventTree->Branch("chi2piVVect", &m_outputEvent.m_chi2piVVect);
    m_pEventTree->Branch("isUVChi2piAvailableVect", &m_outputEvent.m_isUVChi2piAvailableVect);
    m_pEventTree->Branch("chi2piUVVect", &m_outputEvent.m_chi2piUVVect);
    
    m_pEventTree->Branch("isWChi2MIPAvailableVect", &m_outputEvent.m_isWChi2MIPAvailableVect);
    m_pEventTree->Branch("chi2MIPWVect", &m_outputEvent.m_chi2MIPWVect);
    m_pEventTree->Branch("isUChi2MIPAvailableVect", &m_outputEvent.m_isUChi2MIPAvailableVect);
    m_pEventTree->Branch("chi2MIPUVect", &m_outputEvent.m_chi2MIPUVect);
    m_pEventTree->Branch("isVChi2MIPAvailableVect", &m_outputEvent.m_isVChi2MIPAvailableVect);
    m_pEventTree->Branch("chi2MIPVVect", &m_outputEvent.m_chi2MIPVVect);
    m_pEventTree->Branch("isUVChi2MIPAvailableVect", &m_outputEvent.m_isUVChi2MIPAvailableVect);
    m_pEventTree->Branch("chi2MIPUVVect", &m_outputEvent.m_chi2MIPUVVect);

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
    
    m_pEventTree->Branch("isWBraggmuAvailableVect", &m_outputEvent.m_isWBraggmuAvailableVect);
    m_pEventTree->Branch("braggmuWVect", &m_outputEvent.m_braggmuWVect);
    m_pEventTree->Branch("isUBraggmuAvailableVect", &m_outputEvent.m_isUBraggmuAvailableVect);
    m_pEventTree->Branch("braggmuUVect", &m_outputEvent.m_braggmuUVect);
    m_pEventTree->Branch("isVBraggmuAvailableVect", &m_outputEvent.m_isVBraggmuAvailableVect);
    m_pEventTree->Branch("braggmuVVect", &m_outputEvent.m_braggmuVVect);
    m_pEventTree->Branch("isUVBraggmuAvailableVect", &m_outputEvent.m_isUVBraggmuAvailableVect);
    m_pEventTree->Branch("braggmuUVVect", &m_outputEvent.m_braggmuUVVect);
    
    m_pEventTree->Branch("isWBraggpiAvailableVect", &m_outputEvent.m_isWBraggpiAvailableVect);
    m_pEventTree->Branch("braggpiWVect", &m_outputEvent.m_braggpiWVect);
    m_pEventTree->Branch("isUBraggpiAvailableVect", &m_outputEvent.m_isUBraggpiAvailableVect);
    m_pEventTree->Branch("braggpiUVect", &m_outputEvent.m_braggpiUVect);
    m_pEventTree->Branch("isVBraggpiAvailableVect", &m_outputEvent.m_isVBraggpiAvailableVect);
    m_pEventTree->Branch("braggpiVVect", &m_outputEvent.m_braggpiVVect);
    m_pEventTree->Branch("isUVBraggpiAvailableVect", &m_outputEvent.m_isUVBraggpiAvailableVect);
    m_pEventTree->Branch("braggpiUVVect", &m_outputEvent.m_braggpiUVVect);
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
//    m_outputEvent.m_topologicalScore = -std::numeric_limits<float>::max();
    m_outputEvent.m_nFinalStatePFPs = -std::numeric_limits<int>::max();
    m_outputEvent.m_hasMatchedMCParticleVect.clear();
    m_outputEvent.m_matchedMCParticleIdVect.clear();
    m_outputEvent.m_truePdgCodeVect.clear();
    m_outputEvent.m_truthMatchCompletenessVect.clear();
    m_outputEvent.m_truthMatchPurityVect.clear();
    m_outputEvent.m_trueEnergyVect.clear();
    m_outputEvent.m_trueKEVect.clear();
    m_outputEvent.m_trueRangeVect.clear();
    m_outputEvent.m_trueMomentumXVect.clear();
    m_outputEvent.m_trueMomentumYVect.clear();
    m_outputEvent.m_trueMomentumZVect.clear();
    m_outputEvent.m_trueStartXVect.clear();
    m_outputEvent.m_trueStartYVect.clear();
    m_outputEvent.m_trueStartZVect.clear();
    m_outputEvent.m_trueIsContainedVect.clear();
    m_outputEvent.m_trueIsStoppingVect.clear();
    m_outputEvent.m_trueNScattersVect.clear();
    m_outputEvent.m_trueIsGoldenVect.clear();
    m_outputEvent.m_nHitsUVect.clear();
    m_outputEvent.m_nHitsVVect.clear();
    m_outputEvent.m_nHitsWVect.clear();
    m_outputEvent.m_nDescendentHitsUVect.clear();
    m_outputEvent.m_nDescendentHitsVVect.clear();
    m_outputEvent.m_nDescendentHitsWVect.clear();
    m_outputEvent.m_trackShowerVect.clear();
    m_outputEvent.m_nDaughtersVect.clear();
    m_outputEvent.m_nDescendentsVect.clear();
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
    m_outputEvent.m_xyAngleVect.clear();
    m_outputEvent.m_xzAngleVect.clear();
    m_outputEvent.m_lengthVect.clear();
    m_outputEvent.m_rangeVect.clear();
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
    m_outputEvent.m_isWChi2muAvailableVect.clear();
    m_outputEvent.m_chi2muWVect.clear();
    m_outputEvent.m_isUChi2muAvailableVect.clear();
    m_outputEvent.m_chi2muUVect.clear();
    m_outputEvent.m_isVChi2muAvailableVect.clear();
    m_outputEvent.m_chi2muVVect.clear();
    m_outputEvent.m_isUVChi2muAvailableVect.clear();
    m_outputEvent.m_chi2muUVVect.clear();
    m_outputEvent.m_isWChi2piAvailableVect.clear();
    m_outputEvent.m_chi2piWVect.clear();
    m_outputEvent.m_isUChi2piAvailableVect.clear();
    m_outputEvent.m_chi2piUVect.clear();
    m_outputEvent.m_isVChi2piAvailableVect.clear();
    m_outputEvent.m_chi2piVVect.clear();
    m_outputEvent.m_isUVChi2piAvailableVect.clear();
    m_outputEvent.m_chi2piUVVect.clear();
    m_outputEvent.m_isWChi2MIPAvailableVect.clear();
    m_outputEvent.m_chi2MIPWVect.clear();
    m_outputEvent.m_isUChi2MIPAvailableVect.clear();
    m_outputEvent.m_chi2MIPUVect.clear();
    m_outputEvent.m_isVChi2MIPAvailableVect.clear();
    m_outputEvent.m_chi2MIPVVect.clear();
    m_outputEvent.m_isUVChi2MIPAvailableVect.clear();
    m_outputEvent.m_chi2MIPUVVect.clear();
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
    m_outputEvent.m_isWBraggMIPBackwardAvailableVect.clear();
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
    m_outputEvent.m_isWBraggmuAvailableVect.clear();
    m_outputEvent.m_braggmuWVect.clear();
    m_outputEvent.m_isUBraggmuAvailableVect.clear();
    m_outputEvent.m_braggmuUVect.clear();
    m_outputEvent.m_isVBraggmuAvailableVect.clear();
    m_outputEvent.m_braggmuVVect.clear();
    m_outputEvent.m_isUVBraggmuAvailableVect.clear();
    m_outputEvent.m_braggmuUVVect.clear();
    m_outputEvent.m_isWBraggpiAvailableVect.clear();
    m_outputEvent.m_braggpiWVect.clear();
    m_outputEvent.m_isUBraggpiAvailableVect.clear();
    m_outputEvent.m_braggpiUVect.clear();
    m_outputEvent.m_isVBraggpiAvailableVect.clear();
    m_outputEvent.m_braggpiVVect.clear();
    m_outputEvent.m_isUVBraggpiAvailableVect.clear();
    m_outputEvent.m_braggpiUVVect.clear();
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
    const auto pfParticleMap = RecoHelper::GetPFParticleMap(allPFParticles);
    const auto backtrackerData = AnalysisHelper::GetBacktrackerData(event, m_config().MCTruthLabel(), m_config().MCParticleLabel(), m_config().BacktrackerLabel(), m_config().PFParticleLabel());
    const auto finalStates = RecoHelper::GetNeutrinoFinalStates(allPFParticles); 
    const auto pfpToTracks = CollectionHelper::GetAssociation<recob::PFParticle, recob::Track>(event, m_config().PFParticleLabel(), m_config().TrackLabel());
    const auto pfpToHits = CollectionHelper::GetAssociationViaCollection<recob::PFParticle, recob::Cluster, recob::Hit>(event,  m_config().PFParticleLabel(), m_config().PFParticleLabel(), m_config().PFParticleLabel());
    const auto trackToPIDs = CollectionHelper::GetAssociation<recob::Track, anab::ParticleID>(event, m_config().TrackLabel(), m_config().PIDLabel());
    const auto trackToCalorimetries = CollectionHelper::GetAssociation<recob::Track, anab::Calorimetry>(event, m_config().TrackLabel(), m_config().CalorimetryLabel());
    const auto pfpToMetadata = CollectionHelper::GetAssociation<recob::PFParticle, larpandoraobj::PFParticleMetadata>(event, m_config().PFParticleLabel());
    const auto pSpaceChargeService = RecoHelper::GetSpaceChargeService();

    const TruthHelper::Interaction interaction(event, m_config().MCTruthLabel(), m_config().MCParticleLabel());
    const auto mcParticleMap = TruthHelper::GetMCParticleMap(interaction.GetAllMCParticles());

    // Set the event level reco information
    this->SetEventRecoInfo(event, allPFParticles, finalStates, pfpToMetadata, pSpaceChargeService);

    // Set the particle level reco information
    for (unsigned int i = 0; i < finalStates.size(); ++i)

    {
        const auto &finalState = finalStates.at(i);
        this->SetPFParticleInfo(i, finalState, backtrackerData, mcParticleMap, pfParticleMap, pfpToHits, pfpToTracks, trackToPIDs, trackToCalorimetries, pfpToMetadata, pSpaceChargeService);
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
//        const auto neutrinoMetadata = CollectionHelper::GetSingleAssociated(neutrino, pfpToMetadata);
//        m_outputEvent.m_topologicalScore = RecoHelper::GetTopologicalScore(neutrinoMetadata); 

        // Count the final states
        m_outputEvent.m_nFinalStatePFPs = finalStates.size();
    }
    else
    {
        m_outputEvent.m_recoNuVtx = TVector3(-std::numeric_limits<float>::max(), -std::numeric_limits<float>::max(), -std::numeric_limits<float>::max());
        m_outputEvent.m_isRecoNuFiducial = false;
//        m_outputEvent.m_topologicalScore = -std::numeric_limits<float>::max();
        m_outputEvent.m_nFinalStatePFPs = 0;
    }
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void EventSelection::SetPFParticleInfo(const unsigned int index, const art::Ptr<recob::PFParticle> &finalState, const BacktrackHelper::BacktrackerData &backtrackerData, const MCParticleMap &mcParticleMap, const PFParticleMap &pfParticleMap, const Association<recob::PFParticle, recob::Hit> &pfpToHits, const Association<recob::PFParticle, recob::Track> &pfpToTracks, const Association<recob::Track, anab::ParticleID> &trackToPIDs, const Association<recob::Track, anab::Calorimetry> &trackToCalorimetries, const Association<recob::PFParticle, larpandoraobj::PFParticleMetadata> &pfpToMetadata, const spacecharge::SpaceChargeService::provider_type *const pSpaceChargeService)
{
    this->SetPFParticleMCParticleMatchInfo(finalState, backtrackerData, mcParticleMap);
    this->SetPFParticlePandoraInfo(finalState, pfpToHits, pfParticleMap, pfpToMetadata);

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

void EventSelection::SetPFParticleMCParticleMatchInfo(const art::Ptr<recob::PFParticle> &finalState, const BacktrackHelper::BacktrackerData &backtrackerData, const MCParticleMap &mcParticleMap)
{
    bool hasMatchedMCParticle = false;
    int matchedMCParticleId = -std::numeric_limits<int>::max();
    int truePdgCode = -std::numeric_limits<int>::max();
    float truthMatchCompleteness = -std::numeric_limits<float>::max();
    float truthMatchPurity = -std::numeric_limits<float>::max();
    float trueEnergy = -std::numeric_limits<float>::max();
    float trueKE = -std::numeric_limits<float>::max();
    float trueRange = -std::numeric_limits<float>::max();
    TVector3 trueMomentum = TVector3(-std::numeric_limits<float>::max(), -std::numeric_limits<float>::max(), -std::numeric_limits<float>::max());
    TVector3 trueStart = TVector3(-std::numeric_limits<float>::max(), -std::numeric_limits<float>::max(), -std::numeric_limits<float>::max());
    bool trueIsContained = false;
    bool trueIsStopping = false;
    int trueNScatters = -std::numeric_limits<int>::max();
    bool isGolden = false;

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
        trueRange = this->GetRange(matchedMCP);
        trueMomentum = TVector3(matchedMCP->Px(), matchedMCP->Py(), matchedMCP->Pz());
        trueStart = TVector3(matchedMCP->Vx(), matchedMCP->Vy(), matchedMCP->Vz());
        hasMatchedMCParticle = true;
        
        const auto trueEndMomentum = matchedMCP->Momentum(std::max(static_cast<unsigned int>(0), matchedMCP->NumberTrajectoryPoints() - 2)).Vect().Mag();
        trueIsStopping = (trueEndMomentum <= std::numeric_limits<float>::epsilon());
   
        const auto trueEnd = TVector3(matchedMCP->EndPosition().Vect());
        trueIsContained = AnalysisHelper::IsContained(trueStart, m_config().ContainmentBorder()) && AnalysisHelper::IsContained(trueEnd, m_config().ContainmentBorder());

        trueNScatters = 0;
        for (const auto &daughter : TruthHelper::GetDaughters(matchedMCP, mcParticleMap))
        {
            if (daughter->Process() == "hadElastic")
                trueNScatters++;
        }
            
        isGolden = trueIsStopping && (trueNScatters == 0) && trueIsContained;
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
    m_outputEvent.m_trueRangeVect.push_back(trueRange);
    m_outputEvent.m_trueMomentumXVect.push_back(trueMomentum.X());
    m_outputEvent.m_trueMomentumYVect.push_back(trueMomentum.Y());
    m_outputEvent.m_trueMomentumZVect.push_back(trueMomentum.Z());
    m_outputEvent.m_trueStartXVect.push_back(trueStart.X());
    m_outputEvent.m_trueStartYVect.push_back(trueStart.Y());
    m_outputEvent.m_trueStartZVect.push_back(trueStart.Z());
    m_outputEvent.m_trueIsContainedVect.push_back(trueIsContained);
    m_outputEvent.m_trueIsStoppingVect.push_back(trueIsStopping);
    m_outputEvent.m_trueNScattersVect.push_back(trueNScatters);
    m_outputEvent.m_trueIsGoldenVect.push_back(isGolden);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void EventSelection::SetPFParticlePandoraInfo(const art::Ptr<recob::PFParticle> &finalState, const Association<recob::PFParticle, recob::Hit> &pfpToHits, const PFParticleMap &pfParticleMap, const Association<recob::PFParticle, larpandoraobj::PFParticleMetadata> &pfpToMetadata)
{
    // Count the hits in the particle
    const auto &hits = CollectionHelper::GetManyAssociated(finalState, pfpToHits);
    m_outputEvent.m_nHitsUVect.push_back(RecoHelper::CountHitsInView(hits, geo::kU));
    m_outputEvent.m_nHitsVVect.push_back(RecoHelper::CountHitsInView(hits, geo::kV));
    m_outputEvent.m_nHitsWVect.push_back(RecoHelper::CountHitsInView(hits, geo::kW));
    
    // Count the descendent hits
    HitVector descendentHits;
    const auto &descendentParticles = RecoHelper::GetDownstreamParticles(finalState, pfParticleMap);
    for (const auto &descendentParticle : descendentParticles)
    {
        const auto hitsInParticle =  CollectionHelper::GetManyAssociated(descendentParticle, pfpToHits);
        descendentHits.insert(descendentHits.end(), hitsInParticle.begin(), hitsInParticle.end());
    }
    
    m_outputEvent.m_nDescendentHitsUVect.push_back(RecoHelper::CountHitsInView(descendentHits, geo::kU));
    m_outputEvent.m_nDescendentHitsVVect.push_back(RecoHelper::CountHitsInView(descendentHits, geo::kV));
    m_outputEvent.m_nDescendentHitsWVect.push_back(RecoHelper::CountHitsInView(descendentHits, geo::kW));

    // Get the track-shower score
    const auto metadata = CollectionHelper::GetSingleAssociated(finalState, pfpToMetadata);
    m_outputEvent.m_trackShowerVect.push_back(RecoHelper::GetTrackScore(metadata));

    // Count the daughters and descendents
    m_outputEvent.m_nDaughtersVect.push_back(finalState->NumDaughters());
    m_outputEvent.m_nDescendentsVect.push_back(descendentParticles.size() - 1);
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
    m_outputEvent.m_xyAngleVect.push_back(-std::numeric_limits<float>::max());
    m_outputEvent.m_xzAngleVect.push_back(-std::numeric_limits<float>::max());
    m_outputEvent.m_lengthVect.push_back(-std::numeric_limits<float>::max());
    m_outputEvent.m_rangeVect.push_back(-std::numeric_limits<float>::max());
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
    m_outputEvent.m_xyAngleVect.push_back(this->GetXYAngle(dir));
    m_outputEvent.m_xzAngleVect.push_back(this->GetXZAngle(dir));
    m_outputEvent.m_lengthVect.push_back(track->Length());
    m_outputEvent.m_rangeVect.push_back(this->GetRange(track, pSpaceChargeService));
    m_outputEvent.m_isContainedVect.push_back(AnalysisHelper::IsContained(start, m_config().ContainmentBorder()) && AnalysisHelper::IsContained(end, m_config().ContainmentBorder()));
}

// -----------------------------------------------------------------------------------------------------------------------------------------

float EventSelection::GetYZAngle(const TVector3 &dir)
{
    return atan2(dir.Z(), dir.Y());
}

// -----------------------------------------------------------------------------------------------------------------------------------------

float EventSelection::GetXYAngle(const TVector3 &dir)
{
    return atan2(dir.Y(), dir.X());
}

// -----------------------------------------------------------------------------------------------------------------------------------------

float EventSelection::GetXZAngle(const TVector3 &dir)
{
    return atan2(dir.Z(), dir.X());
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
    
    m_outputEvent.m_isWChi2muAvailableVect.push_back(false);
    m_outputEvent.m_chi2muWVect.push_back(-std::numeric_limits<float>::max());
    m_outputEvent.m_isUChi2muAvailableVect.push_back(false);
    m_outputEvent.m_chi2muUVect.push_back(-std::numeric_limits<float>::max());
    m_outputEvent.m_isVChi2muAvailableVect.push_back(false);
    m_outputEvent.m_chi2muVVect.push_back(-std::numeric_limits<float>::max());
    m_outputEvent.m_isUVChi2muAvailableVect.push_back(false);
    m_outputEvent.m_chi2muUVVect.push_back(-std::numeric_limits<float>::max());

    m_outputEvent.m_isWChi2piAvailableVect.push_back(false);
    m_outputEvent.m_chi2piWVect.push_back(-std::numeric_limits<float>::max());
    m_outputEvent.m_isUChi2piAvailableVect.push_back(false);
    m_outputEvent.m_chi2piUVect.push_back(-std::numeric_limits<float>::max());
    m_outputEvent.m_isVChi2piAvailableVect.push_back(false);
    m_outputEvent.m_chi2piVVect.push_back(-std::numeric_limits<float>::max());
    m_outputEvent.m_isUVChi2piAvailableVect.push_back(false);
    m_outputEvent.m_chi2piUVVect.push_back(-std::numeric_limits<float>::max());
    
    m_outputEvent.m_isWChi2MIPAvailableVect.push_back(false);
    m_outputEvent.m_chi2MIPWVect.push_back(-std::numeric_limits<float>::max());
    m_outputEvent.m_isUChi2MIPAvailableVect.push_back(false);
    m_outputEvent.m_chi2MIPUVect.push_back(-std::numeric_limits<float>::max());
    m_outputEvent.m_isVChi2MIPAvailableVect.push_back(false);
    m_outputEvent.m_chi2MIPVVect.push_back(-std::numeric_limits<float>::max());
    m_outputEvent.m_isUVChi2MIPAvailableVect.push_back(false);
    m_outputEvent.m_chi2MIPUVVect.push_back(-std::numeric_limits<float>::max());

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

    m_outputEvent.m_isWBraggMIPBackwardAvailableVect.push_back(false);
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
    
    m_outputEvent.m_isWBraggmuAvailableVect.push_back(false);
    m_outputEvent.m_braggmuWVect.push_back(-std::numeric_limits<float>::max());
    m_outputEvent.m_isUBraggmuAvailableVect.push_back(false);
    m_outputEvent.m_braggmuUVect.push_back(-std::numeric_limits<float>::max());
    m_outputEvent.m_isVBraggmuAvailableVect.push_back(false);
    m_outputEvent.m_braggmuVVect.push_back(-std::numeric_limits<float>::max());
    m_outputEvent.m_isUVBraggmuAvailableVect.push_back(false);
    m_outputEvent.m_braggmuUVVect.push_back(-std::numeric_limits<float>::max());
    
    m_outputEvent.m_isWBraggpiAvailableVect.push_back(false);
    m_outputEvent.m_braggpiWVect.push_back(-std::numeric_limits<float>::max());
    m_outputEvent.m_isUBraggpiAvailableVect.push_back(false);
    m_outputEvent.m_braggpiUVect.push_back(-std::numeric_limits<float>::max());
    m_outputEvent.m_isVBraggpiAvailableVect.push_back(false);
    m_outputEvent.m_braggpiVVect.push_back(-std::numeric_limits<float>::max());
    m_outputEvent.m_isUVBraggpiAvailableVect.push_back(false);
    m_outputEvent.m_braggpiUVVect.push_back(-std::numeric_limits<float>::max());
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
    
    float chi2muW = -std::numeric_limits<float>::max();
    float chi2muU = -std::numeric_limits<float>::max();
    float chi2muV = -std::numeric_limits<float>::max();
    float chi2muUV = -std::numeric_limits<float>::max();
    
    float chi2piW = -std::numeric_limits<float>::max();
    float chi2piU = -std::numeric_limits<float>::max();
    float chi2piV = -std::numeric_limits<float>::max();
    float chi2piUV = -std::numeric_limits<float>::max();
    
    float chi2MIPW = -std::numeric_limits<float>::max();
    float chi2MIPU = -std::numeric_limits<float>::max();
    float chi2MIPV = -std::numeric_limits<float>::max();
    float chi2MIPUV = -std::numeric_limits<float>::max();

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
    
    float braggmuW = -std::numeric_limits<float>::max();
    float braggmuU = -std::numeric_limits<float>::max();
    float braggmuV = -std::numeric_limits<float>::max();
    float braggmuUV = -std::numeric_limits<float>::max();
    
    float braggpiW = -std::numeric_limits<float>::max();
    float braggpiU = -std::numeric_limits<float>::max();
    float braggpiV = -std::numeric_limits<float>::max();
    float braggpiUV = -std::numeric_limits<float>::max();
    
    bool isWChi2pAvailable = false;
    bool isUChi2pAvailable = false;
    bool isVChi2pAvailable = false;
    bool isUVChi2pAvailable = false;
    
    bool isWChi2muAvailable = false;
    bool isUChi2muAvailable = false;
    bool isVChi2muAvailable = false;
    bool isUVChi2muAvailable = false;
    
    bool isWChi2piAvailable = false;
    bool isUChi2piAvailable = false;
    bool isVChi2piAvailable = false;
    bool isUVChi2piAvailable = false;
    
    bool isWChi2MIPAvailable = false;
    bool isUChi2MIPAvailable = false;
    bool isVChi2MIPAvailable = false;
    bool isUVChi2MIPAvailable = false;

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
    
    bool isWBraggmuAvailable = false;
    bool isUBraggmuAvailable = false;
    bool isVBraggmuAvailable = false;
    bool isUVBraggmuAvailable = false;
    
    bool isWBraggpiAvailable = false;
    bool isUBraggpiAvailable = false;
    bool isVBraggpiAvailable = false;
    bool isUVBraggpiAvailable = false;

    // Extract the raw PID values
    for (const auto &algo : pid->fParticleIDAlgScores)
    {
        const auto view = this->GetView(algo.fPlaneMask);

        // Chi2 algorithm under the proton hypothesis
        if (algo.fAlgName == "Chi2" && algo.fAssumedPdg == 2212)
            this->SetPIDVariables(view, algo.fValue, chi2pW, chi2pU, chi2pV, isWChi2pAvailable, isUChi2pAvailable, isVChi2pAvailable);
        
        // Chi2 algorithm under the muon hypothesis
        if (algo.fAlgName == "Chi2" && algo.fAssumedPdg == 13)
            this->SetPIDVariables(view, algo.fValue, chi2muW, chi2muU, chi2muV, isWChi2muAvailable, isUChi2muAvailable, isVChi2muAvailable);
        
        // Chi2 algorithm under the pion hypothesis
        if (algo.fAlgName == "Chi2" && algo.fAssumedPdg == 211)
            this->SetPIDVariables(view, algo.fValue, chi2piW, chi2piU, chi2piV, isWChi2piAvailable, isUChi2piAvailable, isVChi2piAvailable);
        
        // Chi2 algorithm under the MIP hypothesis
        if (algo.fAlgName == "Chi2" && algo.fAssumedPdg == 0)
            this->SetPIDVariables(view, algo.fValue, chi2MIPW, chi2MIPU, chi2MIPV, isWChi2MIPAvailable, isUChi2MIPAvailable, isVChi2MIPAvailable);

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
        
        if (algo.fAlgName == "BraggPeakLLH" && algo.fTrackDir == anab::kForward && algo.fAssumedPdg == 13)
            this->SetPIDVariables(view, algo.fValue, braggmuW, braggmuU, braggmuV, isWBraggmuAvailable, isUBraggmuAvailable, isVBraggmuAvailable);
        
        if (algo.fAlgName == "BraggPeakLLH" && algo.fTrackDir == anab::kForward && algo.fAssumedPdg == 211)
            this->SetPIDVariables(view, algo.fValue, braggpiW, braggpiU, braggpiV, isWBraggpiAvailable, isUBraggpiAvailable, isVBraggpiAvailable);
    }

    // Set the combined induction plane variables
    this->CombineInductionPlanes(yzAngle, nHitsU, nHitsV, chi2pU, chi2pV, isUChi2pAvailable, isVChi2pAvailable, chi2pUV, isUVChi2pAvailable);
    this->CombineInductionPlanes(yzAngle, nHitsU, nHitsV, chi2muU, chi2muV, isUChi2muAvailable, isVChi2muAvailable, chi2muUV, isUVChi2muAvailable);
    this->CombineInductionPlanes(yzAngle, nHitsU, nHitsV, chi2piU, chi2piV, isUChi2piAvailable, isVChi2piAvailable, chi2piUV, isUVChi2piAvailable);
    this->CombineInductionPlanes(yzAngle, nHitsU, nHitsV, chi2MIPU, chi2MIPV, isUChi2MIPAvailable, isVChi2MIPAvailable, chi2MIPUV, isUVChi2MIPAvailable);

    this->CombineInductionPlanes(yzAngle, nHitsU, nHitsV, braggpU, braggpV, isUBraggpAvailable, isVBraggpAvailable, braggpUV, isUVBraggpAvailable);
    this->CombineInductionPlanes(yzAngle, nHitsU, nHitsV, braggpBackwardU, braggpBackwardV, isUBraggpBackwardAvailable, isVBraggpBackwardAvailable, braggpBackwardUV, isUVBraggpBackwardAvailable);
    this->CombineInductionPlanes(yzAngle, nHitsU, nHitsV, braggMIPU, braggMIPV, isUBraggMIPAvailable, isVBraggMIPAvailable, braggMIPUV, isUVBraggMIPAvailable);
    this->CombineInductionPlanes(yzAngle, nHitsU, nHitsV, braggMIPBackwardU, braggMIPBackwardV, isUBraggMIPBackwardAvailable, isVBraggMIPBackwardAvailable, braggMIPBackwardUV, isUVBraggMIPBackwardAvailable);
    this->CombineInductionPlanes(yzAngle, nHitsU, nHitsV, braggmuU, braggmuV, isUBraggmuAvailable, isVBraggmuAvailable, braggmuUV, isUVBraggmuAvailable);
    this->CombineInductionPlanes(yzAngle, nHitsU, nHitsV, braggpiU, braggpiV, isUBraggpiAvailable, isVBraggpiAvailable, braggpiUV, isUVBraggpiAvailable);

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
    
    m_outputEvent.m_isWChi2muAvailableVect.push_back(isWChi2muAvailable);
    m_outputEvent.m_chi2muWVect.push_back(chi2muW);
    m_outputEvent.m_isUChi2muAvailableVect.push_back(isUChi2muAvailable);
    m_outputEvent.m_chi2muUVect.push_back(chi2muU);
    m_outputEvent.m_isVChi2muAvailableVect.push_back(isVChi2muAvailable);
    m_outputEvent.m_chi2muVVect.push_back(chi2muV);
    m_outputEvent.m_isUVChi2muAvailableVect.push_back(isUVChi2muAvailable);
    m_outputEvent.m_chi2muUVVect.push_back(chi2muUV);
    
    m_outputEvent.m_isWChi2piAvailableVect.push_back(isWChi2piAvailable);
    m_outputEvent.m_chi2piWVect.push_back(chi2piW);
    m_outputEvent.m_isUChi2piAvailableVect.push_back(isUChi2piAvailable);
    m_outputEvent.m_chi2piUVect.push_back(chi2piU);
    m_outputEvent.m_isVChi2piAvailableVect.push_back(isVChi2piAvailable);
    m_outputEvent.m_chi2piVVect.push_back(chi2piV);
    m_outputEvent.m_isUVChi2piAvailableVect.push_back(isUVChi2piAvailable);
    m_outputEvent.m_chi2piUVVect.push_back(chi2piUV);
    
    m_outputEvent.m_isWChi2MIPAvailableVect.push_back(isWChi2MIPAvailable);
    m_outputEvent.m_chi2MIPWVect.push_back(chi2MIPW);
    m_outputEvent.m_isUChi2MIPAvailableVect.push_back(isUChi2MIPAvailable);
    m_outputEvent.m_chi2MIPUVect.push_back(chi2MIPU);
    m_outputEvent.m_isVChi2MIPAvailableVect.push_back(isVChi2MIPAvailable);
    m_outputEvent.m_chi2MIPVVect.push_back(chi2MIPV);
    m_outputEvent.m_isUVChi2MIPAvailableVect.push_back(isUVChi2MIPAvailable);
    m_outputEvent.m_chi2MIPUVVect.push_back(chi2MIPUV);

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
    
    m_outputEvent.m_isWBraggmuAvailableVect.push_back(isWBraggmuAvailable);
    m_outputEvent.m_braggmuWVect.push_back(braggmuW);
    m_outputEvent.m_isUBraggmuAvailableVect.push_back(isUBraggmuAvailable);
    m_outputEvent.m_braggmuUVect.push_back(braggmuU);
    m_outputEvent.m_isVBraggmuAvailableVect.push_back(isVBraggmuAvailable);
    m_outputEvent.m_braggmuVVect.push_back(braggmuV);
    m_outputEvent.m_isUVBraggmuAvailableVect.push_back(isUVBraggmuAvailable);
    m_outputEvent.m_braggmuUVVect.push_back(braggmuUV);
    
    m_outputEvent.m_isWBraggpiAvailableVect.push_back(isWBraggpiAvailable);
    m_outputEvent.m_braggpiWVect.push_back(braggpiW);
    m_outputEvent.m_isUBraggpiAvailableVect.push_back(isUBraggpiAvailable);
    m_outputEvent.m_braggpiUVect.push_back(braggpiU);
    m_outputEvent.m_isVBraggpiAvailableVect.push_back(isVBraggpiAvailable);
    m_outputEvent.m_braggpiVVect.push_back(braggpiV);
    m_outputEvent.m_isUVBraggpiAvailableVect.push_back(isUVBraggpiAvailable);
    m_outputEvent.m_braggpiUVVect.push_back(braggpiUV);
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

float EventSelection::GetRange(const art::Ptr<simb::MCParticle> &mcParticle) const
{
    float range = 0.f;

    for (unsigned int i = 1; i < mcParticle->NumberTrajectoryPoints(); ++i)
    {
        range += (mcParticle->Position(i).Vect() - mcParticle->Position(i - 1).Vect()).Mag();
    }

    return range;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

std::vector<size_t> EventSelection::GetValidPoints(const art::Ptr<recob::Track> &track) const
{
    std::vector<size_t> validPoints;

    const auto firstValidPoint = track->FirstValidPoint();
    validPoints.push_back(firstValidPoint);

    auto nextValidPoint = track->NextValidPoint(firstValidPoint + 1);
    while (nextValidPoint != recob::TrackTrajectory::InvalidIndex)
    {
        validPoints.push_back(nextValidPoint);
        nextValidPoint = track->NextValidPoint(nextValidPoint + 1);
    }

    return validPoints;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

float EventSelection::GetRange(const art::Ptr<recob::Track> &track, const spacecharge::SpaceChargeService::provider_type *const pSpaceChargeService) const
{
    const auto validPoints = this->GetValidPoints(track);
    if (validPoints.size() < 2)
        return 0.f;

    float range = 0.f;
    for (unsigned int i = 1; i < validPoints.size(); ++i)
    {
        const auto pos = track->LocationAtPoint(validPoints.at(i));
        const auto posPrev = track->LocationAtPoint(validPoints.at(i - 1));

        const auto posVect = TVector3(pos.X(), pos.Y(), pos.Z());
        const auto posPrevVect = TVector3(posPrev.X(), posPrev.Y(), posPrev.Z());

        const auto posVectSCE = RecoHelper::CorrectForSpaceCharge(posVect, pSpaceChargeService);
        const auto posPrevVectSCE = RecoHelper::CorrectForSpaceCharge(posPrevVect, pSpaceChargeService);

        range += (posVectSCE - posPrevVectSCE).Mag();
    }

    return range;
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

    if (m_outputEvent.m_hasMatchedMCParticleVect.size() != expectedSize)
        throw cet::exception("EventSelection::ValidateOutputVectorSizes") << "m_hasMatchedMCParticleVect - is out of sync" << std::endl;

    if (m_outputEvent.m_matchedMCParticleIdVect.size() != expectedSize)
        throw cet::exception("EventSelection::ValidateOutputVectorSizes") << "m_matchedMCParticleIdVect - is out of sync" << std::endl;

    if (m_outputEvent.m_truePdgCodeVect.size() != expectedSize)
        throw cet::exception("EventSelection::ValidateOutputVectorSizes") << "m_truePdgCodeVect - is out of sync" << std::endl;

    if (m_outputEvent.m_truthMatchCompletenessVect.size() != expectedSize)
        throw cet::exception("EventSelection::ValidateOutputVectorSizes") << "m_truthMatchCompletenessVect - is out of sync" << std::endl;

    if (m_outputEvent.m_truthMatchPurityVect.size() != expectedSize)
        throw cet::exception("EventSelection::ValidateOutputVectorSizes") << "m_truthMatchPurityVect - is out of sync" << std::endl;

    if (m_outputEvent.m_trueEnergyVect.size() != expectedSize)
        throw cet::exception("EventSelection::ValidateOutputVectorSizes") << "m_trueEnergyVect - is out of sync" << std::endl;

    if (m_outputEvent.m_trueKEVect.size() != expectedSize)
        throw cet::exception("EventSelection::ValidateOutputVectorSizes") << "m_trueKEVect - is out of sync" << std::endl;
    
    if (m_outputEvent.m_trueRangeVect.size() != expectedSize)
        throw cet::exception("EventSelection::ValidateOutputVectorSizes") << "m_trueRangeVect - is out of sync" << std::endl;

    if (m_outputEvent.m_trueMomentumXVect.size() != expectedSize)
        throw cet::exception("EventSelection::ValidateOutputVectorSizes") << "m_trueMomentumXVect - is out of sync" << std::endl;

    if (m_outputEvent.m_trueMomentumYVect.size() != expectedSize)
        throw cet::exception("EventSelection::ValidateOutputVectorSizes") << "m_trueMomentumYVect - is out of sync" << std::endl;

    if (m_outputEvent.m_trueMomentumZVect.size() != expectedSize)
        throw cet::exception("EventSelection::ValidateOutputVectorSizes") << "m_trueMomentumZVect - is out of sync" << std::endl;

    if (m_outputEvent.m_trueStartXVect.size() != expectedSize)
        throw cet::exception("EventSelection::ValidateOutputVectorSizes") << "m_trueStartXVect - is out of sync" << std::endl;

    if (m_outputEvent.m_trueStartYVect.size() != expectedSize)
        throw cet::exception("EventSelection::ValidateOutputVectorSizes") << "m_trueStartYVect - is out of sync" << std::endl;

    if (m_outputEvent.m_trueStartZVect.size() != expectedSize)
        throw cet::exception("EventSelection::ValidateOutputVectorSizes") << "m_trueStartZVect - is out of sync" << std::endl;
    
    if (m_outputEvent.m_trueIsContainedVect.size() != expectedSize)
        throw cet::exception("EventSelection::ValidateOutputVectorSizes") << "m_trueIsContainedVect - is out of sync" << std::endl;
    
    if (m_outputEvent.m_trueIsStoppingVect.size() != expectedSize)
        throw cet::exception("EventSelection::ValidateOutputVectorSizes") << "m_trueIsStoppingVect - is out of sync" << std::endl;
    
    if (m_outputEvent.m_trueNScattersVect.size() != expectedSize)
        throw cet::exception("EventSelection::ValidateOutputVectorSizes") << "m_trueNScattersVect - is out of sync" << std::endl;
    
    if (m_outputEvent.m_trueIsGoldenVect.size() != expectedSize)
        throw cet::exception("EventSelection::ValidateOutputVectorSizes") << "m_trueIsGoldenVect - is out of sync" << std::endl;

    if (m_outputEvent.m_nHitsUVect.size() != expectedSize)
        throw cet::exception("EventSelection::ValidateOutputVectorSizes") << "m_nHitsUVect - is out of sync" << std::endl;

    if (m_outputEvent.m_nHitsVVect.size() != expectedSize)
        throw cet::exception("EventSelection::ValidateOutputVectorSizes") << "m_nHitsVVect - is out of sync" << std::endl;

    if (m_outputEvent.m_nHitsWVect.size() != expectedSize)
        throw cet::exception("EventSelection::ValidateOutputVectorSizes") << "m_nHitsWVect - is out of sync" << std::endl;
    
    if (m_outputEvent.m_nDescendentHitsUVect.size() != expectedSize)
        throw cet::exception("EventSelection::ValidateOutputVectorSizes") << "m_nDescendentHitsUVect - is out of sync" << std::endl;

    if (m_outputEvent.m_nDescendentHitsVVect.size() != expectedSize)
        throw cet::exception("EventSelection::ValidateOutputVectorSizes") << "m_nDescendentHitsVVect - is out of sync" << std::endl;

    if (m_outputEvent.m_nDescendentHitsWVect.size() != expectedSize)
        throw cet::exception("EventSelection::ValidateOutputVectorSizes") << "m_nDescendentHitsWVect - is out of sync" << std::endl;

    if (m_outputEvent.m_trackShowerVect.size() != expectedSize)
        throw cet::exception("EventSelection::ValidateOutputVectorSizes") << "m_trackShowerVect - is out of sync" << std::endl;
    
    if (m_outputEvent.m_nDaughtersVect.size() != expectedSize)
        throw cet::exception("EventSelection::ValidateOutputVectorSizes") << "m_nDaughtersVect - is out of sync" << std::endl;
    
    if (m_outputEvent.m_nDescendentsVect.size() != expectedSize)
        throw cet::exception("EventSelection::ValidateOutputVectorSizes") << "m_nDescendentsVect - is out of sync" << std::endl;

    if (m_outputEvent.m_startXVect.size() != expectedSize)
        throw cet::exception("EventSelection::ValidateOutputVectorSizes") << "m_startXVect - is out of sync" << std::endl;

    if (m_outputEvent.m_startYVect.size() != expectedSize)
        throw cet::exception("EventSelection::ValidateOutputVectorSizes") << "m_startYVect - is out of sync" << std::endl;

    if (m_outputEvent.m_startZVect.size() != expectedSize)
        throw cet::exception("EventSelection::ValidateOutputVectorSizes") << "m_startZVect - is out of sync" << std::endl;

    if (m_outputEvent.m_endXVect.size() != expectedSize)
        throw cet::exception("EventSelection::ValidateOutputVectorSizes") << "m_endXVect - is out of sync" << std::endl;

    if (m_outputEvent.m_endYVect.size() != expectedSize)
        throw cet::exception("EventSelection::ValidateOutputVectorSizes") << "m_endYVect - is out of sync" << std::endl;

    if (m_outputEvent.m_endZVect.size() != expectedSize)
        throw cet::exception("EventSelection::ValidateOutputVectorSizes") << "m_endZVect - is out of sync" << std::endl;

    if (m_outputEvent.m_directionXVect.size() != expectedSize)
        throw cet::exception("EventSelection::ValidateOutputVectorSizes") << "m_directionXVect - is out of sync" << std::endl;

    if (m_outputEvent.m_directionYVect.size() != expectedSize)
        throw cet::exception("EventSelection::ValidateOutputVectorSizes") << "m_directionYVect - is out of sync" << std::endl;

    if (m_outputEvent.m_directionZVect.size() != expectedSize)
        throw cet::exception("EventSelection::ValidateOutputVectorSizes") << "m_directionZVect - is out of sync" << std::endl;

    if (m_outputEvent.m_thetaVect.size() != expectedSize)
        throw cet::exception("EventSelection::ValidateOutputVectorSizes") << "m_thetaVect - is out of sync" << std::endl;

    if (m_outputEvent.m_phiVect.size() != expectedSize)
        throw cet::exception("EventSelection::ValidateOutputVectorSizes") << "m_phiVect - is out of sync" << std::endl;

    if (m_outputEvent.m_yzAngleVect.size() != expectedSize)
        throw cet::exception("EventSelection::ValidateOutputVectorSizes") << "m_yzAngleVect - is out of sync" << std::endl;
    
    if (m_outputEvent.m_xyAngleVect.size() != expectedSize)
        throw cet::exception("EventSelection::ValidateOutputVectorSizes") << "m_xyAngleVect - is out of sync" << std::endl;
    
    if (m_outputEvent.m_xzAngleVect.size() != expectedSize)
        throw cet::exception("EventSelection::ValidateOutputVectorSizes") << "m_xzAngleVect - is out of sync" << std::endl;

    if (m_outputEvent.m_lengthVect.size() != expectedSize)
        throw cet::exception("EventSelection::ValidateOutputVectorSizes") << "m_lengthVect - is out of sync" << std::endl;
    
    if (m_outputEvent.m_rangeVect.size() != expectedSize)
        throw cet::exception("EventSelection::ValidateOutputVectorSizes") << "m_rangeVect - is out of sync" << std::endl;

    if (m_outputEvent.m_isContainedVect.size() != expectedSize)
        throw cet::exception("EventSelection::ValidateOutputVectorSizes") << "m_isContainedVect - is out of sync" << std::endl;

    if (m_outputEvent.m_hasCalorimetryInfoVect.size() != expectedSize)
        throw cet::exception("EventSelection::ValidateOutputVectorSizes") << "m_hasCalorimetryInfoVect - is out of sync" << std::endl;

    if (m_outputEvent.m_dedxPerHitUVect.size() != expectedSize)
        throw cet::exception("EventSelection::ValidateOutputVectorSizes") << "m_dedxPerHitUVect - is out of sync" << std::endl;

    if (m_outputEvent.m_dedxPerHitVVect.size() != expectedSize)
        throw cet::exception("EventSelection::ValidateOutputVectorSizes") << "m_dedxPerHitVVect - is out of sync" << std::endl;

    if (m_outputEvent.m_dedxPerHitWVect.size() != expectedSize)
        throw cet::exception("EventSelection::ValidateOutputVectorSizes") << "m_dedxPerHitWVect - is out of sync" << std::endl;

    if (m_outputEvent.m_residualRangePerHitUVect.size() != expectedSize)
        throw cet::exception("EventSelection::ValidateOutputVectorSizes") << "m_residualRangePerHitUVect - is out of sync" << std::endl;

    if (m_outputEvent.m_residualRangePerHitVVect.size() != expectedSize)
        throw cet::exception("EventSelection::ValidateOutputVectorSizes") << "m_residualRangePerHitVVect - is out of sync" << std::endl;

    if (m_outputEvent.m_residualRangePerHitWVect.size() != expectedSize)
        throw cet::exception("EventSelection::ValidateOutputVectorSizes") << "m_residualRangePerHitWVect - is out of sync" << std::endl;

    if (m_outputEvent.m_hasPIDInfoVect.size() != expectedSize)
        throw cet::exception("EventSelection::ValidateOutputVectorSizes") << "m_hasPIDInfoVect - is out of sync" << std::endl;

    if (m_outputEvent.m_isWChi2pAvailableVect.size() != expectedSize)
        throw cet::exception("EventSelection::ValidateOutputVectorSizes") << "m_isWChi2pAvailableVect - is out of sync" << std::endl;

    if (m_outputEvent.m_chi2pWVect.size() != expectedSize)
        throw cet::exception("EventSelection::ValidateOutputVectorSizes") << "m_chi2pWVect - is out of sync" << std::endl;

    if (m_outputEvent.m_isUChi2pAvailableVect.size() != expectedSize)
        throw cet::exception("EventSelection::ValidateOutputVectorSizes") << "m_isUChi2pAvailableVect - is out of sync" << std::endl;

    if (m_outputEvent.m_chi2pUVect.size() != expectedSize)
        throw cet::exception("EventSelection::ValidateOutputVectorSizes") << "m_chi2pUVect - is out of sync" << std::endl;

    if (m_outputEvent.m_isVChi2pAvailableVect.size() != expectedSize)
        throw cet::exception("EventSelection::ValidateOutputVectorSizes") << "m_isVChi2pAvailableVect - is out of sync" << std::endl;

    if (m_outputEvent.m_chi2pVVect.size() != expectedSize)
        throw cet::exception("EventSelection::ValidateOutputVectorSizes") << "m_chi2pVVect - is out of sync" << std::endl;

    if (m_outputEvent.m_isUVChi2pAvailableVect.size() != expectedSize)
        throw cet::exception("EventSelection::ValidateOutputVectorSizes") << "m_isUVChi2pAvailableVect - is out of sync" << std::endl;

    if (m_outputEvent.m_chi2pUVVect.size() != expectedSize)
        throw cet::exception("EventSelection::ValidateOutputVectorSizes") << "m_chi2pUVVect - is out of sync" << std::endl;
    
    if (m_outputEvent.m_isWChi2muAvailableVect.size() != expectedSize)
        throw cet::exception("EventSelection::ValidateOutputVectorSizes") << "m_isWChi2muAvailableVect - is out of sync" << std::endl;

    if (m_outputEvent.m_chi2muWVect.size() != expectedSize)
        throw cet::exception("EventSelection::ValidateOutputVectorSizes") << "m_chi2muWVect - is out of sync" << std::endl;

    if (m_outputEvent.m_isUChi2muAvailableVect.size() != expectedSize)
        throw cet::exception("EventSelection::ValidateOutputVectorSizes") << "m_isUChi2muAvailableVect - is out of sync" << std::endl;

    if (m_outputEvent.m_chi2muUVect.size() != expectedSize)
        throw cet::exception("EventSelection::ValidateOutputVectorSizes") << "m_chi2muUVect - is out of sync" << std::endl;

    if (m_outputEvent.m_isVChi2muAvailableVect.size() != expectedSize)
        throw cet::exception("EventSelection::ValidateOutputVectorSizes") << "m_isVChi2muAvailableVect - is out of sync" << std::endl;

    if (m_outputEvent.m_chi2muVVect.size() != expectedSize)
        throw cet::exception("EventSelection::ValidateOutputVectorSizes") << "m_chi2muVVect - is out of sync" << std::endl;

    if (m_outputEvent.m_isUVChi2muAvailableVect.size() != expectedSize)
        throw cet::exception("EventSelection::ValidateOutputVectorSizes") << "m_isUVChi2muAvailableVect - is out of sync" << std::endl;

    if (m_outputEvent.m_chi2muUVVect.size() != expectedSize)
        throw cet::exception("EventSelection::ValidateOutputVectorSizes") << "m_chi2muUVVect - is out of sync" << std::endl;
    
    if (m_outputEvent.m_isWChi2piAvailableVect.size() != expectedSize)
        throw cet::exception("EventSelection::ValidateOutputVectorSizes") << "m_isWChi2piAvailableVect - is out of sync" << std::endl;

    if (m_outputEvent.m_chi2piWVect.size() != expectedSize)
        throw cet::exception("EventSelection::ValidateOutputVectorSizes") << "m_chi2piWVect - is out of sync" << std::endl;

    if (m_outputEvent.m_isUChi2piAvailableVect.size() != expectedSize)
        throw cet::exception("EventSelection::ValidateOutputVectorSizes") << "m_isUChi2piAvailableVect - is out of sync" << std::endl;

    if (m_outputEvent.m_chi2piUVect.size() != expectedSize)
        throw cet::exception("EventSelection::ValidateOutputVectorSizes") << "m_chi2piUVect - is out of sync" << std::endl;

    if (m_outputEvent.m_isVChi2piAvailableVect.size() != expectedSize)
        throw cet::exception("EventSelection::ValidateOutputVectorSizes") << "m_isVChi2piAvailableVect - is out of sync" << std::endl;

    if (m_outputEvent.m_chi2piVVect.size() != expectedSize)
        throw cet::exception("EventSelection::ValidateOutputVectorSizes") << "m_chi2piVVect - is out of sync" << std::endl;

    if (m_outputEvent.m_isUVChi2piAvailableVect.size() != expectedSize)
        throw cet::exception("EventSelection::ValidateOutputVectorSizes") << "m_isUVChi2piAvailableVect - is out of sync" << std::endl;

    if (m_outputEvent.m_chi2piUVVect.size() != expectedSize)
        throw cet::exception("EventSelection::ValidateOutputVectorSizes") << "m_chi2piUVVect - is out of sync" << std::endl;
    
    if (m_outputEvent.m_isWChi2MIPAvailableVect.size() != expectedSize)
        throw cet::exception("EventSelection::ValidateOutputVectorSizes") << "m_isWChi2MIPAvailableVect - is out of sync" << std::endl;

    if (m_outputEvent.m_chi2MIPWVect.size() != expectedSize)
        throw cet::exception("EventSelection::ValidateOutputVectorSizes") << "m_chi2MIPWVect - is out of sync" << std::endl;

    if (m_outputEvent.m_isUChi2MIPAvailableVect.size() != expectedSize)
        throw cet::exception("EventSelection::ValidateOutputVectorSizes") << "m_isUChi2MIPAvailableVect - is out of sync" << std::endl;

    if (m_outputEvent.m_chi2MIPUVect.size() != expectedSize)
        throw cet::exception("EventSelection::ValidateOutputVectorSizes") << "m_chi2MIPUVect - is out of sync" << std::endl;

    if (m_outputEvent.m_isVChi2MIPAvailableVect.size() != expectedSize)
        throw cet::exception("EventSelection::ValidateOutputVectorSizes") << "m_isVChi2MIPAvailableVect - is out of sync" << std::endl;

    if (m_outputEvent.m_chi2MIPVVect.size() != expectedSize)
        throw cet::exception("EventSelection::ValidateOutputVectorSizes") << "m_chi2MIPVVect - is out of sync" << std::endl;

    if (m_outputEvent.m_isUVChi2MIPAvailableVect.size() != expectedSize)
        throw cet::exception("EventSelection::ValidateOutputVectorSizes") << "m_isUVChi2MIPAvailableVect - is out of sync" << std::endl;

    if (m_outputEvent.m_chi2MIPUVVect.size() != expectedSize)
        throw cet::exception("EventSelection::ValidateOutputVectorSizes") << "m_chi2MIPUVVect - is out of sync" << std::endl;

    if (m_outputEvent.m_isWBraggpAvailableVect.size() != expectedSize)
        throw cet::exception("EventSelection::ValidateOutputVectorSizes") << "m_isWBraggpAvailableVect - is out of sync" << std::endl;

    if (m_outputEvent.m_braggpWVect.size() != expectedSize)
        throw cet::exception("EventSelection::ValidateOutputVectorSizes") << "m_braggpWVect - is out of sync" << std::endl;

    if (m_outputEvent.m_isUBraggpAvailableVect.size() != expectedSize)
        throw cet::exception("EventSelection::ValidateOutputVectorSizes") << "m_isUBraggpAvailableVect - is out of sync" << std::endl;

    if (m_outputEvent.m_braggpUVect.size() != expectedSize)
        throw cet::exception("EventSelection::ValidateOutputVectorSizes") << "m_braggpUVect - is out of sync" << std::endl;

    if (m_outputEvent.m_isVBraggpAvailableVect.size() != expectedSize)
        throw cet::exception("EventSelection::ValidateOutputVectorSizes") << "m_isVBraggpAvailableVect - is out of sync" << std::endl;

    if (m_outputEvent.m_braggpVVect.size() != expectedSize)
        throw cet::exception("EventSelection::ValidateOutputVectorSizes") << "m_braggpVVect - is out of sync" << std::endl;

    if (m_outputEvent.m_isUVBraggpAvailableVect.size() != expectedSize)
        throw cet::exception("EventSelection::ValidateOutputVectorSizes") << "m_isUVBraggpAvailableVect - is out of sync" << std::endl;

    if (m_outputEvent.m_braggpUVVect.size() != expectedSize)
        throw cet::exception("EventSelection::ValidateOutputVectorSizes") << "m_braggpUVVect - is out of sync" << std::endl;

    if (m_outputEvent.m_isWBraggpBackwardAvailableVect.size() != expectedSize)
        throw cet::exception("EventSelection::ValidateOutputVectorSizes") << "m_isWBraggpBackwardAvailableVect - is out of sync" << std::endl;

    if (m_outputEvent.m_braggpBackwardWVect.size() != expectedSize)
        throw cet::exception("EventSelection::ValidateOutputVectorSizes") << "m_braggpBackwardWVect - is out of sync" << std::endl;

    if (m_outputEvent.m_isUBraggpBackwardAvailableVect.size() != expectedSize)
        throw cet::exception("EventSelection::ValidateOutputVectorSizes") << "m_isUBraggpBackwardAvailableVect - is out of sync" << std::endl;

    if (m_outputEvent.m_braggpBackwardUVect.size() != expectedSize)
        throw cet::exception("EventSelection::ValidateOutputVectorSizes") << "m_braggpBackwardUVect - is out of sync" << std::endl;

    if (m_outputEvent.m_isVBraggpBackwardAvailableVect.size() != expectedSize)
        throw cet::exception("EventSelection::ValidateOutputVectorSizes") << "m_isVBraggpBackwardAvailableVect - is out of sync" << std::endl;

    if (m_outputEvent.m_braggpBackwardVVect.size() != expectedSize)
        throw cet::exception("EventSelection::ValidateOutputVectorSizes") << "m_braggpBackwardVVect - is out of sync" << std::endl;

    if (m_outputEvent.m_isUVBraggpBackwardAvailableVect.size() != expectedSize)
        throw cet::exception("EventSelection::ValidateOutputVectorSizes") << "m_isUVBraggpBackwardAvailableVect - is out of sync" << std::endl;

    if (m_outputEvent.m_braggpBackwardUVVect.size() != expectedSize)
        throw cet::exception("EventSelection::ValidateOutputVectorSizes") << "m_braggpBackwardUVVect - is out of sync" << std::endl;

    if (m_outputEvent.m_isWBraggMIPAvailableVect.size() != expectedSize)
        throw cet::exception("EventSelection::ValidateOutputVectorSizes") << "m_isWBraggMIPAvailableVect - is out of sync" << std::endl;

    if (m_outputEvent.m_braggMIPWVect.size() != expectedSize)
        throw cet::exception("EventSelection::ValidateOutputVectorSizes") << "m_braggMIPWVect - is out of sync" << std::endl;

    if (m_outputEvent.m_isUBraggMIPAvailableVect.size() != expectedSize)
        throw cet::exception("EventSelection::ValidateOutputVectorSizes") << "m_isUBraggMIPAvailableVect - is out of sync" << std::endl;

    if (m_outputEvent.m_braggMIPUVect.size() != expectedSize)
        throw cet::exception("EventSelection::ValidateOutputVectorSizes") << "m_braggMIPUVect - is out of sync" << std::endl;

    if (m_outputEvent.m_isVBraggMIPAvailableVect.size() != expectedSize)
        throw cet::exception("EventSelection::ValidateOutputVectorSizes") << "m_isVBraggMIPAvailableVect - is out of sync" << std::endl;

    if (m_outputEvent.m_braggMIPVVect.size() != expectedSize)
        throw cet::exception("EventSelection::ValidateOutputVectorSizes") << "m_braggMIPVVect - is out of sync" << std::endl;

    if (m_outputEvent.m_isUVBraggMIPAvailableVect.size() != expectedSize)
        throw cet::exception("EventSelection::ValidateOutputVectorSizes") << "m_isUVBraggMIPAvailableVect - is out of sync" << std::endl;

    if (m_outputEvent.m_braggMIPUVVect.size() != expectedSize)
        throw cet::exception("EventSelection::ValidateOutputVectorSizes") << "m_braggMIPUVVect - is out of sync" << std::endl;

    if (m_outputEvent.m_isWBraggMIPBackwardAvailableVect.size() != expectedSize)
        throw cet::exception("EventSelection::ValidateOutputVectorSizes") << "m_isWBraggMIPBackwardAvailableVect - is out of sync" << std::endl;

    if (m_outputEvent.m_braggMIPBackwardWVect.size() != expectedSize)
        throw cet::exception("EventSelection::ValidateOutputVectorSizes") << "m_braggMIPBackwardWVect - is out of sync" << std::endl;

    if (m_outputEvent.m_isUBraggMIPBackwardAvailableVect.size() != expectedSize)
        throw cet::exception("EventSelection::ValidateOutputVectorSizes") << "m_isUBraggMIPBackwardAvailableVect - is out of sync" << std::endl;

    if (m_outputEvent.m_braggMIPBackwardUVect.size() != expectedSize)
        throw cet::exception("EventSelection::ValidateOutputVectorSizes") << "m_braggMIPBackwardUVect - is out of sync" << std::endl;

    if (m_outputEvent.m_isVBraggMIPBackwardAvailableVect.size() != expectedSize)
        throw cet::exception("EventSelection::ValidateOutputVectorSizes") << "m_isVBraggMIPBackwardAvailableVect - is out of sync" << std::endl;

    if (m_outputEvent.m_braggMIPBackwardVVect.size() != expectedSize)
        throw cet::exception("EventSelection::ValidateOutputVectorSizes") << "m_braggMIPBackwardVVect - is out of sync" << std::endl;

    if (m_outputEvent.m_isUVBraggMIPBackwardAvailableVect.size() != expectedSize)
        throw cet::exception("EventSelection::ValidateOutputVectorSizes") << "m_isUVBraggMIPBackwardAvailableVect - is out of sync" << std::endl;

    if (m_outputEvent.m_braggMIPBackwardUVVect.size() != expectedSize)
        throw cet::exception("EventSelection::ValidateOutputVectorSizes") << "m_braggMIPBackwardUVVect - is out of sync" << std::endl;

    if (m_outputEvent.m_isWBraggRatioAvailableVect.size() != expectedSize)
        throw cet::exception("EventSelection::ValidateOutputVectorSizes") << "m_isWBraggRatioAvailableVect - is out of sync" << std::endl;

    if (m_outputEvent.m_braggRatioWVect.size() != expectedSize)
        throw cet::exception("EventSelection::ValidateOutputVectorSizes") << "m_braggRatioWVect - is out of sync" << std::endl;

    if (m_outputEvent.m_isUVBraggRatioAvailableVect.size() != expectedSize)
        throw cet::exception("EventSelection::ValidateOutputVectorSizes") << "m_isUVBraggRatioAvailableVect - is out of sync" << std::endl;

    if (m_outputEvent.m_braggRatioUVVect.size() != expectedSize)
        throw cet::exception("EventSelection::ValidateOutputVectorSizes") << "m_braggRatioUVVect - is out of sync" << std::endl;
    
    if (m_outputEvent.m_isWBraggmuAvailableVect.size() != expectedSize)
        throw cet::exception("EventSelection::ValidateOutputVectorSizes") << "m_isWBraggmuAvailableVect - is out of sync" << std::endl;

    if (m_outputEvent.m_braggmuWVect.size() != expectedSize)
        throw cet::exception("EventSelection::ValidateOutputVectorSizes") << "m_braggmuWVect - is out of sync" << std::endl;

    if (m_outputEvent.m_isUBraggmuAvailableVect.size() != expectedSize)
        throw cet::exception("EventSelection::ValidateOutputVectorSizes") << "m_isUBraggmuAvailableVect - is out of sync" << std::endl;

    if (m_outputEvent.m_braggmuUVect.size() != expectedSize)
        throw cet::exception("EventSelection::ValidateOutputVectorSizes") << "m_braggmuUVect - is out of sync" << std::endl;

    if (m_outputEvent.m_isVBraggmuAvailableVect.size() != expectedSize)
        throw cet::exception("EventSelection::ValidateOutputVectorSizes") << "m_isVBraggmuAvailableVect - is out of sync" << std::endl;

    if (m_outputEvent.m_braggmuVVect.size() != expectedSize)
        throw cet::exception("EventSelection::ValidateOutputVectorSizes") << "m_braggmuVVect - is out of sync" << std::endl;

    if (m_outputEvent.m_isUVBraggmuAvailableVect.size() != expectedSize)
        throw cet::exception("EventSelection::ValidateOutputVectorSizes") << "m_isUVBraggmuAvailableVect - is out of sync" << std::endl;

    if (m_outputEvent.m_braggmuUVVect.size() != expectedSize)
        throw cet::exception("EventSelection::ValidateOutputVectorSizes") << "m_braggmuUVVect - is out of sync" << std::endl;
    
    if (m_outputEvent.m_isWBraggpiAvailableVect.size() != expectedSize)
        throw cet::exception("EventSelection::ValidateOutputVectorSizes") << "m_isWBraggpiAvailableVect - is out of sync" << std::endl;

    if (m_outputEvent.m_braggpiWVect.size() != expectedSize)
        throw cet::exception("EventSelection::ValidateOutputVectorSizes") << "m_braggpiWVect - is out of sync" << std::endl;

    if (m_outputEvent.m_isUBraggpiAvailableVect.size() != expectedSize)
        throw cet::exception("EventSelection::ValidateOutputVectorSizes") << "m_isUBraggpiAvailableVect - is out of sync" << std::endl;

    if (m_outputEvent.m_braggpiUVect.size() != expectedSize)
        throw cet::exception("EventSelection::ValidateOutputVectorSizes") << "m_braggpiUVect - is out of sync" << std::endl;

    if (m_outputEvent.m_isVBraggpiAvailableVect.size() != expectedSize)
        throw cet::exception("EventSelection::ValidateOutputVectorSizes") << "m_isVBraggpiAvailableVect - is out of sync" << std::endl;

    if (m_outputEvent.m_braggpiVVect.size() != expectedSize)
        throw cet::exception("EventSelection::ValidateOutputVectorSizes") << "m_braggpiVVect - is out of sync" << std::endl;

    if (m_outputEvent.m_isUVBraggpiAvailableVect.size() != expectedSize)
        throw cet::exception("EventSelection::ValidateOutputVectorSizes") << "m_isUVBraggpiAvailableVect - is out of sync" << std::endl;

    if (m_outputEvent.m_braggpiUVVect.size() != expectedSize)
        throw cet::exception("EventSelection::ValidateOutputVectorSizes") << "m_braggpiUVVect - is out of sync" << std::endl;
}

} // namespace ubcc1pi
