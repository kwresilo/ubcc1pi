/**
 *  @file  ubcc1pi/Analyzers/PIDStudy_module.cc
 *
 *  @brief The implementation file for the PID study analyzer.
 */

#include "ubcc1pi/Analyzers/PIDStudy.h"

#include "ubcc1pi/Helpers/TruthHelper.h"
#include "ubcc1pi/Helpers/AnalysisHelper.h"
#include "ubcc1pi/Helpers/DebugHelper.h"


namespace ubcc1pi
{

PIDStudy::PIDStudy(const art::EDAnalyzer::Table<Config> &config) :
    art::EDAnalyzer(config),
    m_config(config)
{
    // Setup the output trees
    art::ServiceHandle<art::TFileService> fileService;
    
    m_pAlgorithmTree = fileService->make<TTree>("algorithms", "");
    m_pAlgorithmTree->Branch("run", &m_outputAlgorithm.m_run);
    m_pAlgorithmTree->Branch("subRun", &m_outputAlgorithm.m_subRun);
    m_pAlgorithmTree->Branch("event", &m_outputAlgorithm.m_event);
    m_pAlgorithmTree->Branch("isSignal", &m_outputAlgorithm.m_isSignal);
    m_pAlgorithmTree->Branch("nPFPHits", &m_outputAlgorithm.m_nPFPHits);
    m_pAlgorithmTree->Branch("truePdgCode", &m_outputAlgorithm.m_truePdgCode);
    m_pAlgorithmTree->Branch("trueMomentum", &m_outputAlgorithm.m_trueMomentum);
    m_pAlgorithmTree->Branch("trueMatchPurity", &m_outputAlgorithm.m_trueMatchPurity);
    m_pAlgorithmTree->Branch("trueMatchCompleteness", &m_outputAlgorithm.m_trueMatchCompleteness);
    m_pAlgorithmTree->Branch("name", &m_outputAlgorithm.m_name);
    m_pAlgorithmTree->Branch("variableType", &m_outputAlgorithm.m_variableType);
    m_pAlgorithmTree->Branch("trackDir", &m_outputAlgorithm.m_trackDir);
    m_pAlgorithmTree->Branch("nDOF", &m_outputAlgorithm.m_nDOF);
    m_pAlgorithmTree->Branch("assumedPdg", &m_outputAlgorithm.m_assumedPdg);
    m_pAlgorithmTree->Branch("planeWUsed", &m_outputAlgorithm.m_planeWUsed);
    m_pAlgorithmTree->Branch("planeUUsed", &m_outputAlgorithm.m_planeUUsed);
    m_pAlgorithmTree->Branch("planeVUsed", &m_outputAlgorithm.m_planeVUsed);
    m_pAlgorithmTree->Branch("value", &m_outputAlgorithm.m_value);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void PIDStudy::analyze(const art::Event &event)
{
    // Extract the labels
    const auto mcTruthLabel = m_config().MCTruthLabel();
    const auto mcParticleLabel = m_config().MCParticleLabel();
    const auto backtrackerLabel = m_config().BacktrackerLabel();
    const auto pfParticleLabel = m_config().PFParticleLabel();
    const auto trackLabel = m_config().TrackLabel();
    const auto pidLabel = m_config().PIDLabel();
    const auto sliceLabel = m_config().SliceLabel();
    const auto hitLabel = m_config().HitLabel();

    // Get the truth level information
    const TruthHelper::Interaction interaction(event, mcTruthLabel, mcParticleLabel);

    const auto pfpToTrack = CollectionHelper::GetAssociation<recob::PFParticle, recob::Track>(event, pfParticleLabel, trackLabel);
    const auto pfpToMetadata = CollectionHelper::GetAssociation<recob::PFParticle, larpandoraobj::PFParticleMetadata>(event, pfParticleLabel);
    const auto trackToPID = CollectionHelper::GetAssociation<recob::Track, anab::ParticleID>(event, trackLabel, pidLabel);
    
    // Get the reco-true matching information
    const auto backtrackerData = AnalysisHelper::GetBacktrackerData(event, mcTruthLabel, mcParticleLabel, backtrackerLabel, pfParticleLabel);
    const auto pfParticles = backtrackerData.GetPFParticles();

    // Get the maximum track length
    float maxTrackLength = -std::numeric_limits<float>::max();
    for (const auto &pfParticle : pfParticles)
    {
        try
        {
            const auto tracks = CollectionHelper::GetManyAssociated(pfParticle, pfpToTrack); 
            
            if (tracks.size() > 1)
                throw cet::exception("PIDStudy::analyze") << " - Multiple tracks associated to PFParticle!" << std::endl;
        
            if (!tracks.empty())
            {
                const auto track = tracks.front();
                if (track->Length() > maxTrackLength)
                {
                    maxTrackLength = track->Length();
                }
            }
        }
        catch (const cet::exception &)
        {
        }
    }

    // Store the event level information
    m_outputAlgorithm.m_run = event.run();
    m_outputAlgorithm.m_subRun = event.subRun();
    m_outputAlgorithm.m_event = event.event();
    m_outputAlgorithm.m_isSignal = AnalysisHelper::IsCC1PiSignal(interaction);

    for (const auto &pfParticle : pfParticles)
    {
        // Store the particle level information
        m_outputAlgorithm.m_nPFPHits = backtrackerData.GetNHits(pfParticle);
        try
        {
            const auto mcParticle = backtrackerData.GetBestMatchedMCParticle(pfParticle);

            m_outputAlgorithm.m_truePdgCode = mcParticle->PdgCode();
            m_outputAlgorithm.m_trueMomentum = mcParticle->P();
            m_outputAlgorithm.m_trueMatchPurity = backtrackerData.GetMatchPurity(pfParticle, mcParticle);
            m_outputAlgorithm.m_trueMatchCompleteness = backtrackerData.GetMatchCompleteness(pfParticle, mcParticle);
        }
        catch (const cet::exception &)
        {
            m_outputAlgorithm.m_truePdgCode = -std::numeric_limits<int>::max();
            m_outputAlgorithm.m_trueMomentum = -std::numeric_limits<float>::max();
            m_outputAlgorithm.m_trueMatchPurity = -std::numeric_limits<float>::max();
            m_outputAlgorithm.m_trueMatchCompleteness = -std::numeric_limits<float>::max();
        }
        
        // --------------------------------------------------------------------------

        // Store the Pandora track/shower ID algorithm output
        m_outputAlgorithm.m_name = "PandoraTrackShower";
        m_outputAlgorithm.m_variableType = static_cast<int>(anab::kScore);
        m_outputAlgorithm.m_trackDir = static_cast<int>(anab::kNoDirection);
        m_outputAlgorithm.m_nDOF = -std::numeric_limits<int>::max();
        m_outputAlgorithm.m_assumedPdg = -std::numeric_limits<int>::max();
        m_outputAlgorithm.m_planeWUsed = true;
        m_outputAlgorithm.m_planeUUsed = true;
        m_outputAlgorithm.m_planeVUsed = true;

        const auto metadata = CollectionHelper::GetSingleAssociated(pfParticle, pfpToMetadata);
        m_outputAlgorithm.m_value = RecoHelper::GetTrackScore(metadata);
        m_pAlgorithmTree->Fill();
        
        // --------------------------------------------------------------------------

        // If the PFParticle is track-like, then score the output of the PID algorithms on the track
        bool hasTrack = false;
        art::Ptr<recob::Track> track;
        try
        {
            const auto tracks = CollectionHelper::GetManyAssociated(pfParticle, pfpToTrack); 
            
            if (tracks.size() > 1)
                throw cet::exception("PIDStudy::analyze") << " - Multiple tracks associated to PFParticle!" << std::endl;
        
            if (!tracks.empty())
            {
                track = tracks.front();
                hasTrack = true;
            }
        }
        catch (const cet::exception &)
        {
        }

        if (!hasTrack)
            continue;
        
        // --------------------------------------------------------------------------
        
        // Store the track length
        m_outputAlgorithm.m_name = "TrackLength";
        m_outputAlgorithm.m_variableType = static_cast<int>(anab::kTrackLength);
        m_outputAlgorithm.m_trackDir = static_cast<int>(anab::kNoDirection);
        m_outputAlgorithm.m_nDOF = -std::numeric_limits<int>::max();
        m_outputAlgorithm.m_assumedPdg = -std::numeric_limits<int>::max();
        m_outputAlgorithm.m_planeWUsed = true;
        m_outputAlgorithm.m_planeUUsed = true;
        m_outputAlgorithm.m_planeVUsed = true;
        m_outputAlgorithm.m_value = track->Length();
        m_pAlgorithmTree->Fill();

        // Store the relative track length
        m_outputAlgorithm.m_name = "RelativeTrackLength";
        m_outputAlgorithm.m_variableType = static_cast<int>(anab::kTrackLength);
        m_outputAlgorithm.m_trackDir = static_cast<int>(anab::kNoDirection);
        m_outputAlgorithm.m_nDOF = -std::numeric_limits<int>::max();
        m_outputAlgorithm.m_assumedPdg = -std::numeric_limits<int>::max();
        m_outputAlgorithm.m_planeWUsed = true;
        m_outputAlgorithm.m_planeUUsed = true;
        m_outputAlgorithm.m_planeVUsed = true;
        m_outputAlgorithm.m_value = track->Length() / maxTrackLength;
        m_pAlgorithmTree->Fill();
        
        // --------------------------------------------------------------------------
      
        const auto pid = CollectionHelper::GetSingleAssociated(track, trackToPID);

        for (const auto &algo : pid->ParticleIDAlgScores())
        {
            m_outputAlgorithm.m_name = algo.fAlgName;
            m_outputAlgorithm.m_variableType = static_cast<int>(algo.fVariableType);
            m_outputAlgorithm.m_trackDir = static_cast<int>(algo.fTrackDir);
            m_outputAlgorithm.m_nDOF = algo.fNdf;
            m_outputAlgorithm.m_assumedPdg = algo.fAssumedPdg;
            m_outputAlgorithm.m_planeWUsed = algo.fPlaneMask[0];
            m_outputAlgorithm.m_planeUUsed = algo.fPlaneMask[1];
            m_outputAlgorithm.m_planeVUsed = algo.fPlaneMask[2];
            m_outputAlgorithm.m_value = algo.fValue;
            m_pAlgorithmTree->Fill();
        }
    }
        
}

} // namespace ubcc1pi
