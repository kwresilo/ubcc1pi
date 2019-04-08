/**
 *  @file  ubcc1pi/Analyzers/RecoStudy_module.cc
 *
 *  @brief The implementation file for the reco study analyzer.
 */

#include "ubcc1pi/Analyzers/RecoStudy.h"

#include "ubcc1pi/Helpers/TruthHelper.h"
#include "ubcc1pi/Helpers/AnalysisHelper.h"
#include "ubcc1pi/Helpers/DebugHelper.h"


namespace ubcc1pi
{

RecoStudy::RecoStudy(const art::EDAnalyzer::Table<Config> &config) :
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
    m_pEventTree->Branch("nProtons", &m_outputEvent.m_nProtons);
    m_pEventTree->Branch("nPFPs", &m_outputEvent.m_nPFPs);
    m_pEventTree->Branch("nMuMatches", &m_outputEvent.m_nMuMatches);
    m_pEventTree->Branch("nPiMatches", &m_outputEvent.m_nPiMatches);
    m_pEventTree->Branch("nProtonsMatched", &m_outputEvent.m_nProtonsMatched);
    m_pEventTree->Branch("nProtonsMatchedOnce", &m_outputEvent.m_nProtonsMatchedOnce);
    m_pEventTree->Branch("nUnmatchedPFPs", &m_outputEvent.m_nUnmatchedPFPs);
    m_pEventTree->Branch("nNuHitsTotal", &m_outputEvent.m_nNuHitsTotal);
    m_pEventTree->Branch("hasChosenSlice", &m_outputEvent.m_hasChosenSlice);
    m_pEventTree->Branch("chosenSlicePurity", &m_outputEvent.m_chosenSlicePurity);
    m_pEventTree->Branch("chosenSliceCompleteness", &m_outputEvent.m_chosenSliceCompleteness);
    m_pEventTree->Branch("isChosenSliceMostComplete", &m_outputEvent.m_isChosenSliceMostComplete);

    m_pParticleTree = fileService->make<TTree>("particles", "");
    m_pParticleTree->Branch("run", &m_outputParticle.m_run);
    m_pParticleTree->Branch("subRun", &m_outputParticle.m_subRun);
    m_pParticleTree->Branch("event", &m_outputParticle.m_event);
    m_pParticleTree->Branch("nPFPsInEvent", &m_outputParticle.m_nPFPsInEvent);
    m_pParticleTree->Branch("eventHasChosenSlice", &m_outputParticle.m_eventHasChosenSlice);
    m_pParticleTree->Branch("chosenSlicePurity", &m_outputParticle.m_chosenSlicePurity);
    m_pParticleTree->Branch("chosenSliceCompleteness", &m_outputParticle.m_chosenSliceCompleteness);
    m_pParticleTree->Branch("isChosenSliceMostComplete", &m_outputParticle.m_isChosenSliceMostComplete);
    m_pParticleTree->Branch("nUnmatchedPFPsInEvent", &m_outputParticle.m_nUnmatchedPFPsInEvent);
    m_pParticleTree->Branch("pdgCode", &m_outputParticle.m_pdgCode);
    m_pParticleTree->Branch("momentum", &m_outputParticle.m_momentum);
    m_pParticleTree->Branch("mcHitWeight", &m_outputParticle.m_mcHitWeight);
    m_pParticleTree->Branch("mcHitWeightFromSlice", &m_outputParticle.m_mcHitWeightFromSlice);
    m_pParticleTree->Branch("nMatches", &m_outputParticle.m_nMatches);
    m_pParticleTree->Branch("isBestMatchTrack", &m_outputParticle.m_isBestMatchTrack);
    m_pParticleTree->Branch("bestMatchPurity", &m_outputParticle.m_bestMatchPurity);
    m_pParticleTree->Branch("bestMatchCompleteness", &m_outputParticle.m_bestMatchCompleteness);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

// TODO factorise out this function
void RecoStudy::analyze(const art::Event &event)
{
    // Extract the labels
    const auto mcTruthLabel = m_config().MCTruthLabel();
    const auto mcParticleLabel = m_config().MCParticleLabel();
    const auto backtrackerLabel = m_config().BacktrackerLabel();
    const auto pfParticleLabel = m_config().PFParticleLabel();
    const auto sliceLabel = m_config().SliceLabel();
    const auto hitLabel = m_config().HitLabel();

    // Get the truth level information
    const TruthHelper::Interaction interaction(event, mcTruthLabel, mcParticleLabel);

    // Only use signal events
    if (!AnalysisHelper::IsCC1PiSignal(interaction))
        return;

    // Get the metada on the slices
    const auto sliceMetadata = AnalysisHelper::GetSliceMetadata(event, mcTruthLabel, mcParticleLabel, backtrackerLabel, pfParticleLabel, sliceLabel, hitLabel);

    // Get the reco-true matching information
    const auto backtrackerData = AnalysisHelper::GetBacktrackerData(event, mcTruthLabel, mcParticleLabel, backtrackerLabel, pfParticleLabel);
    const auto mcParticles = backtrackerData.GetMCParticles();
    const auto pfParticles = backtrackerData.GetPFParticles();

    // Fill the event-level output structure
    m_outputEvent.m_run = event.run();
    m_outputEvent.m_subRun = event.subRun();
    m_outputEvent.m_event = event.event();
    m_outputEvent.m_nuEnergy = interaction.GetNeutrino().E();
    m_outputEvent.m_nProtons = AnalysisHelper::CountParticlesWithPDG(mcParticles, 2212);
    m_outputEvent.m_nPFPs = pfParticles.size();
    m_outputEvent.m_nNuHitsTotal = sliceMetadata.GetTotalNumberOfNuInducedHits();

    const auto chosenSlices = sliceMetadata.GetSelectedNeutrinoSlices();
    if (chosenSlices.size() > 1)
        throw cet::exception("RecoStudy::analyze") << " - Multiple slices were selected as the neutrino." << std::endl;

    m_outputEvent.m_hasChosenSlice = false;
    m_outputEvent.m_chosenSlicePurity = -std::numeric_limits<float>::max();
    m_outputEvent.m_chosenSliceCompleteness = -std::numeric_limits<float>::max();
    m_outputEvent.m_isChosenSliceMostComplete = false;

    if (!chosenSlices.empty())
    {
        const auto chosenSlice = chosenSlices.front();
        m_outputEvent.m_hasChosenSlice = true;
        m_outputEvent.m_chosenSlicePurity = sliceMetadata.GetPurity(chosenSlice);
        m_outputEvent.m_chosenSliceCompleteness = sliceMetadata.GetCompleteness(chosenSlice);
        m_outputEvent.m_isChosenSliceMostComplete = sliceMetadata.IsMostCompleteSliceSelected(); 
    }

    m_outputEvent.m_nProtonsMatched = 0;
    m_outputEvent.m_nProtonsMatchedOnce = 0;
    
    bool foundMu = false;
    bool foundPi = false;
    for (const auto &mcParticle : mcParticles)
    {
        const auto pdg = mcParticle->PdgCode();
        const auto nMatches = backtrackerData.GetBestMatchedPFParticles(mcParticle).size();

        // Fill the event level output structure
        if (!foundMu && pdg == 13) // Muon
        {
            foundMu = true;
            m_outputEvent.m_nMuMatches = nMatches;
        }
        else if (!foundPi && pdg == 211) // Pi+
        {
            foundPi = true;
            m_outputEvent.m_nPiMatches = nMatches;
        }
        else if (pdg == 2212) // Proton
        {
            m_outputEvent.m_nProtonsMatched += ((nMatches > 0) ? 1 : 0);
            m_outputEvent.m_nProtonsMatchedOnce += ((nMatches == 1) ? 1 : 0);
        }
        else
        {
            throw cet::exception("RecoStudy::analyze") << " - The event is not CC1Pi!" << std::endl;
        }
    }
    
    // Count the PFParticles that aren't matched to any MCParticles
    m_outputEvent.m_nUnmatchedPFPs = 0;
    for (const auto &pfParticle : pfParticles)
    {
        try
        {
            backtrackerData.GetBestMatchedMCParticle(pfParticle);
        }
        catch(const cet::exception &)
        {
            m_outputEvent.m_nUnmatchedPFPs++;
        }
    }
    
    // Finished with this event
    m_pEventTree->Fill();

    // -----------------------

    // Fill the output particle structure
    for (const auto &mcParticle : mcParticles)
    {
        // Event information
        m_outputParticle.m_run = event.run();
        m_outputParticle.m_subRun = event.subRun();
        m_outputParticle.m_event = event.event();
        m_outputParticle.m_nPFPsInEvent = m_outputEvent.m_nPFPs;
        m_outputParticle.m_nUnmatchedPFPsInEvent = m_outputEvent.m_nUnmatchedPFPs;
        
        m_outputParticle.m_eventHasChosenSlice = m_outputEvent.m_hasChosenSlice;
        m_outputParticle.m_chosenSlicePurity = m_outputEvent.m_chosenSlicePurity;
        m_outputParticle.m_chosenSliceCompleteness = m_outputEvent.m_chosenSliceCompleteness;
        m_outputParticle.m_isChosenSliceMostComplete = m_outputEvent.m_isChosenSliceMostComplete;

        // MCParticle information
        const auto pdg = mcParticle->PdgCode();
        m_outputParticle.m_pdgCode = pdg;
        m_outputParticle.m_momentum = mcParticle->P();
        m_outputParticle.m_mcHitWeight = backtrackerData.GetWeight(mcParticle);
        m_outputParticle.m_mcHitWeightFromSlice = 0.f;

        if (!chosenSlices.empty())
        {
            // Get the hits from the MCParticle that are also in the chosen slice
            const auto chosenSlice = chosenSlices.front();
            const auto sharedHits = CollectionHelper::GetIntersection(backtrackerData.GetHits(mcParticle), sliceMetadata.GetHits(chosenSlice));
            m_outputParticle.m_mcHitWeightFromSlice = backtrackerData.GetWeight(sharedHits, mcParticle);

            std::cout << "DEBUG" << std::endl;
            std::cout << "ID                = " << mcParticle->TrackId() << std::endl;
            std::cout << "PDG               = " << pdg << std::endl;
            std::cout << "MCParticle hits   = " << backtrackerData.GetHits(mcParticle).size() << std::endl;
            std::cout << "MCParticle weight = " << backtrackerData.GetWeight(backtrackerData.GetHits(mcParticle), mcParticle) << std::endl;
            std::cout << "... checksum ^^^  = " << backtrackerData.GetWeight(mcParticle) << std::endl;
            std::cout << "Slice hits        = " << sliceMetadata.GetHits(chosenSlice).size() << std::endl;
            std::cout << "Shared hits       = " << sharedHits.size() << std::endl;
            std::cout << "Shared weight     = " << m_outputParticle.m_mcHitWeightFromSlice << std::endl;
        }

        // Match information
        const auto matchedPFParticles = backtrackerData.GetBestMatchedPFParticles(mcParticle);
        m_outputParticle.m_nMatches = matchedPFParticles.size();
        
        // Best match information
        m_outputParticle.m_isBestMatchTrack = false;
        m_outputParticle.m_bestMatchPurity = -std::numeric_limits<float>::max();
        m_outputParticle.m_bestMatchCompleteness = -std::numeric_limits<float>::max();

        if (!matchedPFParticles.empty())
        {
            // Find the PFParticle match with the highest purity
            auto bestMatchedPFParticle = matchedPFParticles.front();
            auto bestMatchPurity = -std::numeric_limits<float>::max();

            for (const auto matchedPFParticle : matchedPFParticles)
            {
                const auto purity = backtrackerData.GetMatchPurity(matchedPFParticle, mcParticle);
                if (purity < bestMatchPurity)
                    continue;

                bestMatchedPFParticle = matchedPFParticle;
                bestMatchPurity = purity;
            }

            const auto absBestMatchPdg = std::abs(bestMatchedPFParticle->PdgCode());
            m_outputParticle.m_isBestMatchTrack = (absBestMatchPdg != 22 && absBestMatchPdg != 11);  // Not an electron or a photon
            m_outputParticle.m_bestMatchPurity = bestMatchPurity;
            m_outputParticle.m_bestMatchCompleteness = backtrackerData.GetMatchCompleteness(bestMatchedPFParticle, mcParticle);
        }
        
        // Finished with this particle
        m_pParticleTree->Fill();
    }

    DebugHelper::Print(interaction);
    DebugHelper::Print(sliceMetadata);
    DebugHelper::Print(backtrackerData);
}

} // namespace ubcc1pi
