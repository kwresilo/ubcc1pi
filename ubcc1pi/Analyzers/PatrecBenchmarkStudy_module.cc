/**
 *  @file  ubcc1pi/Analyzers/PatrecBenchmarkStudy_module.cc
 *
 *  @brief The implementation of the patrec benchmark study analyzer.
 */

#include "ubcc1pi/Analyzers/PatrecBenchmarkStudy.h"

#include "ubcc1pi/Helpers/TruthHelper.h"
#include "ubcc1pi/Helpers/AnalysisHelper.h"
#include "ubcc1pi/Helpers/CollectionHelper.h"
#include "ubcc1pi/Helpers/DebugHelper.h"

namespace ubcc1pi
{

PatrecBenchmarkStudy::PatrecBenchmarkStudy(const art::EDAnalyzer::Table<Config> &config) :
    art::EDAnalyzer(config),
    m_config(config)
{
    // Setup the output trees
    art::ServiceHandle<art::TFileService> fileService;
    
    m_pEventTree = fileService->make<TTree>("events", "");
    m_pEventTree->Branch("run", &m_outputEvent.m_run);
    m_pEventTree->Branch("subRun", &m_outputEvent.m_subRun);
    m_pEventTree->Branch("event", &m_outputEvent.m_event);
            
    m_pEventTree->Branch("nuEnergy", &m_outputEvent.m_nuEnergy);
    m_pEventTree->Branch("hasSlice", &m_outputEvent.m_hasSlice);
    m_pEventTree->Branch("nNuHitsTotal", &m_outputEvent.m_nNuHitsTotal);
    m_pEventTree->Branch("hasNuHits", &m_outputEvent.m_hasNuHits);
    m_pEventTree->Branch("mcsPurity", &m_outputEvent.m_mcsPurity);
    m_pEventTree->Branch("mcsCompleteness", &m_outputEvent.m_mcsCompleteness);
    m_pEventTree->Branch("isMCSSelected", &m_outputEvent.m_isMCSSelected);

    m_pEventTree->Branch("mcp_pdg", &m_outputEvent.m_mcp_pdg);
    m_pEventTree->Branch("mcp_momentum", &m_outputEvent.m_mcp_momentum);
    m_pEventTree->Branch("mcp_hitWeight", &m_outputEvent.m_mcp_hitWeight);
    m_pEventTree->Branch("mcp_fracHitsInMCS", &m_outputEvent.m_mcp_fracHitsInMCS);
    m_pEventTree->Branch("mcp_nPFPMatches", &m_outputEvent.m_mcp_nPFPMatches);
    m_pEventTree->Branch("mcp_hasMatchedPFP", &m_outputEvent.m_mcp_hasMatchedPFP);
    m_pEventTree->Branch("mcp_bestMatchPurity", &m_outputEvent.m_mcp_bestMatchPurity);
    m_pEventTree->Branch("mcp_bestMatchCompleteness", &m_outputEvent.m_mcp_bestMatchCompleteness);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void PatrecBenchmarkStudy::analyze(const art::Event &event)
{
    // Get the truth level information
    const TruthHelper::Interaction interaction(event, m_config().MCTruthLabel(), m_config().MCParticleLabel());

    // Only use signal events
    if (!AnalysisHelper::IsCC1PiSignal(interaction))
        return;

    // Clear any previous values, ready for this event
    this->ResetOutputEvent(event);

    // Get the metadata on the slices
    const auto sliceMetadata = AnalysisHelper::GetSliceMetadata(event, m_config().MCTruthLabel(), m_config().MCParticleLabel(),
                                                                m_config().BacktrackerLabel(), m_config().PFParticleLabel(),
                                                                m_config().SliceLabel(), m_config().HitLabel());
    /* BEGIN DEV */
    m_outputEvent.m_nuEnergy = interaction.GetNeutrino().E();
    m_outputEvent.m_hasSlice = !sliceMetadata.GetSlices().empty(); 
    m_outputEvent.m_nNuHitsTotal = sliceMetadata.GetTotalNumberOfNuInducedHits();
    m_outputEvent.m_hasNuHits = (m_outputEvent.m_nNuHitsTotal > 0);

    art::Ptr<recob::Slice> mostCompleteSlice;
    if (m_outputEvent.m_hasNuHits && m_outputEvent.m_hasSlice)
    {
        mostCompleteSlice = sliceMetadata.GetMostCompleteSlice();
        m_outputEvent.m_mcsPurity = sliceMetadata.GetPurity(mostCompleteSlice);
        m_outputEvent.m_mcsCompleteness = sliceMetadata.GetCompleteness(mostCompleteSlice);
        m_outputEvent.m_isMCSSelected = sliceMetadata.IsMostCompleteSliceSelected();
    }

    // Get the reco-true matching information folding in downstream particles for PFParticle & MCParticles
    const auto backtrackerData = AnalysisHelper::GetBacktrackerData(event, m_config().MCTruthLabel(), m_config().MCParticleLabel(),
                                                                    m_config().BacktrackerLabel(), m_config().PFParticleLabel());

    for (const auto &mcParticle : backtrackerData.GetMCParticles())
    {
        // Get the MCParticle details
        const int pdg = mcParticle->PdgCode();
        const float momentum = mcParticle->P();
        const float hitWeight = backtrackerData.GetWeight(mcParticle);
 
        // Get the fraction of hits that are in the most complete slice
        float fracHitsInMCS = -std::numeric_limits<float>::max();
        if (m_outputEvent.m_hasNuHits && m_outputEvent.m_hasSlice && hitWeight > std::numeric_limits<float>::epsilon())
        {
            const auto sharedHits = CollectionHelper::GetIntersection(backtrackerData.GetHits(mcParticle), sliceMetadata.GetHits(mostCompleteSlice));
            const auto sharedHitWeight = backtrackerData.GetWeight(sharedHits, mcParticle);
            fracHitsInMCS = sharedHitWeight / hitWeight;
        }

        // Get the PFParticles which have this MCParticle as their strongest match
        const auto matchedPFParticles = backtrackerData.GetBestMatchedPFParticles(mcParticle);
        const int nPFPMatches = matchedPFParticles.size();
        const bool hasMatchedPFP = (nPFPMatches > 0);

        // Get the matching to the PFParticle with the largest completeness
        float bestMatchPurity = -std::numeric_limits<float>::max();
        float bestMatchCompleteness = -std::numeric_limits<float>::max();
        for (const auto &pfParticle : matchedPFParticles)
        {
            const float purity = backtrackerData.GetMatchPurity(pfParticle, mcParticle);
            const float completeness = backtrackerData.GetMatchCompleteness(pfParticle, mcParticle);

            if (completeness > bestMatchCompleteness)
            {
                bestMatchPurity = purity;
                bestMatchCompleteness = completeness;
            }
        }
    
        // Fill the output structure
        m_outputEvent.m_mcp_pdg.push_back(pdg);
        m_outputEvent.m_mcp_momentum.push_back(momentum);
        m_outputEvent.m_mcp_hitWeight.push_back(hitWeight);
        m_outputEvent.m_mcp_fracHitsInMCS.push_back(fracHitsInMCS);
        m_outputEvent.m_mcp_nPFPMatches.push_back(nPFPMatches);
        m_outputEvent.m_mcp_hasMatchedPFP.push_back(hasMatchedPFP);
        m_outputEvent.m_mcp_bestMatchPurity.push_back(bestMatchPurity);
        m_outputEvent.m_mcp_bestMatchCompleteness.push_back(bestMatchCompleteness);
    }

    m_pEventTree->Fill();
}
    
// -----------------------------------------------------------------------------------------------------------------------------------------

void PatrecBenchmarkStudy::ResetOutputEvent(const art::Event &event)
{
    m_outputEvent.m_run = event.run();
    m_outputEvent.m_subRun = event.subRun();
    m_outputEvent.m_event = event.event();

    m_outputEvent.m_nuEnergy = -std::numeric_limits<float>::max();
    m_outputEvent.m_hasSlice = false;
    m_outputEvent.m_nNuHitsTotal = -std::numeric_limits<int>::max();
    m_outputEvent.m_hasNuHits = false;

    m_outputEvent.m_mcsPurity = -std::numeric_limits<float>::max();
    m_outputEvent.m_mcsCompleteness = -std::numeric_limits<float>::max();
    m_outputEvent.m_isMCSSelected = false;

    m_outputEvent.m_mcp_pdg.clear();
    m_outputEvent.m_mcp_momentum.clear();
    m_outputEvent.m_mcp_hitWeight.clear();
    m_outputEvent.m_mcp_fracHitsInMCS.clear();
    m_outputEvent.m_mcp_nPFPMatches.clear();
    m_outputEvent.m_mcp_hasMatchedPFP.clear();
    m_outputEvent.m_mcp_bestMatchPurity.clear();
    m_outputEvent.m_mcp_bestMatchCompleteness.clear();
}

} // namespace ubcc1pi
