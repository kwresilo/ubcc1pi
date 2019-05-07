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
    
    m_pParticleTree = fileService->make<TTree>("particles", "");
    m_pParticleTree->Branch("run", &m_outputParticle.m_run);
    m_pParticleTree->Branch("subRun", &m_outputParticle.m_subRun);
    m_pParticleTree->Branch("event", &m_outputParticle.m_event);
    m_pParticleTree->Branch("isSignal", &m_outputParticle.m_isSignal);

    m_pParticleTree->Branch("trackShowerScore", &m_outputParticle.m_trackShowerScore);
    m_pParticleTree->Branch("protonMIPScore", &m_outputParticle.m_protonMIPScore);
    m_pParticleTree->Branch("muonPionScore", &m_outputParticle.m_muonPionScore);

    m_pParticleTree->Branch("hasMatchedMCParticle", &m_outputParticle.m_hasMatchedMCParticle);
    m_pParticleTree->Branch("truePdgCode", &m_outputParticle.m_truePdgCode);
    m_pParticleTree->Branch("trueMomentum", &m_outputParticle.m_trueMomentum);
    m_pParticleTree->Branch("trueMatchPurity", &m_outputParticle.m_trueMatchPurity);
    m_pParticleTree->Branch("trueMatchCompleteness", &m_outputParticle.m_trueMatchCompleteness);
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

    // Dump the output to the terminal
    DebugHelper::Print(interaction);
    std::cout << "Is true CC1Pi+? : " << isTrueSignal << std::endl;

    // Event level information
    m_outputParticle.m_run = event.run();
    m_outputParticle.m_subRun = event.subRun();
    m_outputParticle.m_event = event.event();
    m_outputParticle.m_isSignal = isTrueSignal;

    for (const auto &pfParticle : finalStates)
    {
        const auto trackShowerScore = selector.GetTrackShowerScore(pfParticle);
        const auto protonMIPScore = selector.GetProtonMIPScore(pfParticle);
        const auto muonPionScore = selector.GetMuonPionScore(pfParticle);

        std::cout << "PFParticle " << pfParticle->Self() << std::endl;
        std::cout << "  - Track vs. Shower score : " << trackShowerScore << std::endl;
        std::cout << "  - Proton vs. MIP score   : " << protonMIPScore << std::endl;
        std::cout << "  - Muon vs. Pion score    : " << muonPionScore << std::endl;
       
        m_outputParticle.m_trackShowerScore = trackShowerScore;
        m_outputParticle.m_protonMIPScore = protonMIPScore;
        m_outputParticle.m_muonPionScore = muonPionScore;

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

            std::cout << "  - True MCParticle " << mcParticle->TrackId() << std::endl;
            std::cout << "    - PDG                  : " << mcParticle->PdgCode() << std::endl;
            std::cout << "    - Purity               : " << purity << std::endl;
            std::cout << "    - Completeness         : " << completeness << std::endl;
        
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
}

} // namespace ubcc1pi
