/**
 *  @file  ubcc1pi/Objects/EventFactory.cxx
 *
 *  @brief The implementation of the event factory class
 */

#include "ubcc1pi/Objects/EventFactory.h"

#include "ubcc1pi/Helpers/TruthHelper.h"
#include "ubcc1pi/Helpers/RecoHelper.h"
#include "ubcc1pi/Helpers/TruthHelper.h"
#include "ubcc1pi/Helpers/BacktrackHelper.h"
#include "ubcc1pi/Helpers/DebugHelper.h"

namespace ubcc1pi
{

void EventFactory::PopulateEvent(const art::Event &event, const Config &config, Event *pOutputEvent)
{
    pOutputEvent->Reset();
   
    // Get the output data structures
    auto &metadata = pOutputEvent->metadata;
    auto &truth = pOutputEvent->truth;
    auto &reco = pOutputEvent->reco;

    // Populate the structures with information from the input event
    EventFactory::PopulateEventMetadata(event, config, metadata);
    
    MCParticleVector finalStateMCParticles;
    EventFactory::PopulateEventTruthInfo(event, config, truth, finalStateMCParticles);
    EventFactory::PopulateEventRecoInfo(event, config, finalStateMCParticles, reco);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void EventFactory::PopulateEventMetadata(const art::Event &event, const Config &config, Event::Metadata &metadata)
{
    metadata.run.Set(event.run());
    metadata.subRun.Set(event.subRun());
    metadata.event.Set(event.event());
    metadata.hasTruthInfo.Set(config.HasTruthInfo());
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void EventFactory::PopulateEventTruthInfo(const art::Event &event, const Config &config, Event::Truth &truth, MCParticleVector &finalStateMCParticles)
{
    if (!config.HasTruthInfo())
        return;
    
    // Fill the slice information
    EventFactory::PopulateEventTruthSliceInfo(event, config, truth);

    // Get the generator level interaction
    const TruthHelper::Interaction interaction(event, config.MCTruthLabel(), config.MCParticleLabel());
    
    // Fill the interaction type info
    truth.isCC.Set(interaction.GetCCNC() == simb::kCC);
    truth.interactionMode.Set(interaction.GetInteractionMode());
    truth.interactionString.Set(DebugHelper::GetInteractionString(interaction, true));

    // Fill the neutrino info
    const auto neutrino = interaction.GetNeutrino();
    truth.nuPdgCode.Set(neutrino.PdgCode());
    truth.nuEnergy.Set(neutrino.E());
    truth.nuVertex.Set(TVector3(neutrino.Vx(), neutrino.Vy(), neutrino.Vz()));
    
    finalStateMCParticles = TruthHelper::GetPrimaryMCParticles(interaction.GetAllMCParticles());
    truth.nFinalStates.Set(finalStateMCParticles.size());

    // Fill the truth particle info
    const auto hitToMCParticleWeights = BacktrackHelper::GetHitToMCParticleWeightMap(event, config.MCParticleLabel(), config.BacktrackerLabel(), finalStateMCParticles);
    const auto mcParticleToHitWeights = CollectionHelper::GetReversedAssociation(hitToMCParticleWeights);

    for (const auto &mcParticle : finalStateMCParticles)
    {
        Event::Truth::Particle particle;
        EventFactory::PopulateEventTruthParticleInfo(mcParticle, mcParticleToHitWeights, particle);
        truth.particles.push_back(particle);
    }
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void EventFactory::PopulateEventTruthParticleInfo(const art::Ptr<simb::MCParticle> &mcParticle, const MCParticlesToHitWeights &mcParticleToHitWeights, Event::Truth::Particle &particle)
{
    particle.pdgCode.Set(mcParticle->PdgCode());
    particle.momentumX.Set(mcParticle->Px());
    particle.momentumY.Set(mcParticle->Py());
    particle.momentumZ.Set(mcParticle->Pz());
    particle.momentum.Set(mcParticle->P());
    particle.energy.Set(mcParticle->E());
    particle.mass.Set(mcParticle->Mass());

    particle.hitWeightU.Set(BacktrackHelper::GetHitWeightInView(mcParticle, mcParticleToHitWeights, geo::kU));
    particle.hitWeightV.Set(BacktrackHelper::GetHitWeightInView(mcParticle, mcParticleToHitWeights, geo::kV));
    particle.hitWeightW.Set(BacktrackHelper::GetHitWeightInView(mcParticle, mcParticleToHitWeights, geo::kW));
}

// -----------------------------------------------------------------------------------------------------------------------------------------
    
void EventFactory::PopulateEventRecoInfo(const art::Event &event, const Config &config, const MCParticleVector &finalStateMCParticles, Event::Reco &reco)
{
    // Fill the slice information
    EventFactory::PopulateEventRecoSliceInfo(event, config, reco);

    // Find the neutrino PFParticle
    const auto allPFParticles = CollectionHelper::GetCollection<recob::PFParticle>(event, config.PFParticleLabel());
    const auto neutrinos = RecoHelper::GetNeutrinos(allPFParticles);

    if (neutrinos.size() > 1)
        throw cet::exception("EventFactory::PopulateEventRecoInfo") << " - Found multiple neutrino PFParticles." << std::endl;

    const auto hasNeutrino = !neutrinos.empty();
    reco.hasNeutrino.Set(hasNeutrino);

    if (!hasNeutrino)
        return;
    
    // Set the neutrino information
    const auto pSpaceChargeService = RecoHelper::GetSpaceChargeService();
    const auto nuVertex = RecoHelper::GetRecoNeutrinoVertex(event, neutrinos, config.VertexLabel());
    const auto nuVertexSCC = RecoHelper::CorrectForSpaceCharge(nuVertex, pSpaceChargeService);

    reco.nuVertex.Set(nuVertexSCC);
    reco.nuPdgCode.Set(neutrinos.front()->PdgCode());
    
    const auto finalStatePFParticles = RecoHelper::GetNeutrinoFinalStates(allPFParticles);
    reco.nFinalStates.Set(finalStatePFParticles.size());

    // If we have truth info, then run the reco-true-matching
    std::shared_ptr<BacktrackHelper::BacktrackerData> pBacktrackerData;
    if (config.HasTruthInfo())
    {
        const auto hitsToPfps = BacktrackHelper::GetHitToPFParticleMap(event, config.PFParticleLabel(), finalStatePFParticles, true);
        const auto hitsToMcps = BacktrackHelper::GetHitToMCParticleWeightMap(event, config.MCParticleLabel(), config.BacktrackerLabel(), finalStateMCParticles);

        pBacktrackerData = std::make_shared<BacktrackHelper::BacktrackerData>(finalStatePFParticles, finalStateMCParticles, hitsToPfps, hitsToMcps);
    }

    // Create the output reco particles
    for (const auto &pfParticle : finalStatePFParticles)
    {
        Event::Reco::Particle particle;
        EventFactory::PopulateEventRecoParticleInfo(config, pfParticle, finalStateMCParticles, pBacktrackerData, particle);
        reco.particles.push_back(particle);
    }
}   

// -----------------------------------------------------------------------------------------------------------------------------------------

void EventFactory::PopulateEventRecoParticleInfo(const Config &config, const art::Ptr<recob::PFParticle> &pfParticle, const MCParticleVector &finalStateMCParticles, const std::shared_ptr<BacktrackHelper::BacktrackerData> &pBacktrackerData, Event::Reco::Particle &particle)
{
    particle.pdgCode.Set(pfParticle->PdgCode());

    // Reco-true matching information
    if (!config.HasTruthInfo())
        return;

    if (!pBacktrackerData)
        throw cet::exception("EventFactory::PopulateEventRecoParticleInfo") << " - Invalid backtracker data." << std::endl;
    
    // Get the MCParticle with the strongest match 
    art::Ptr<simb::MCParticle> bestMatchedMCParticle;
    try
    {
        bestMatchedMCParticle = pBacktrackerData->GetBestMatchedMCParticle(pfParticle);
        particle.hasMatchedMCParticle.Set(true);
    }
    catch (const cet::exception &)
    {
        particle.hasMatchedMCParticle.Set(false);
    }
    
    // Fill the reco-true matches for this PFParticle
    std::vector<float> truthMatchPurities, truthMatchCompletenesses;
    for (unsigned int i = 0; i < finalStateMCParticles.size(); ++i)
    {
        const auto mcParticle = finalStateMCParticles.at(i);

        truthMatchPurities.push_back(pBacktrackerData->GetMatchPurity(pfParticle, mcParticle));
        truthMatchCompletenesses.push_back(pBacktrackerData->GetMatchCompleteness(pfParticle, mcParticle));
        
        if (particle.hasMatchedMCParticle() && mcParticle == bestMatchedMCParticle)
            particle.bestMatchedMCParticleIndex.Set(i);
    }

    particle.truthMatchPurities.Set(truthMatchPurities);
    particle.truthMatchCompletenesses.Set(truthMatchCompletenesses);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void EventFactory::PopulateEventRecoSliceInfo(const art::Event &event, const Config &config, Event::Reco &reco)
{
    // Reco slice information
    const auto slices = CollectionHelper::GetCollection<recob::Slice>(event, config.SliceLabel());
    const auto slicesToHits = CollectionHelper::GetAssociation<recob::Slice, recob::Hit>(event, config.SliceLabel());
    const auto pfParticlesToSlices = CollectionHelper::GetAssociation<recob::PFParticle, recob::Slice>(event, config.PFParticleLabel(), config.SliceLabel());
    const auto slicesToPFParticles = CollectionHelper::GetReversedAssociation(pfParticlesToSlices);
    const auto pfParticlesToMetadata = CollectionHelper::GetAssociation<recob::PFParticle, larpandoraobj::PFParticleMetadata>(event, config.PFParticleLabel());
    
    // Work out which slice (if any) was selected as the neutrino
    SlicesToBool sliceToIsSelectedAsNu;
    art::Ptr<recob::Slice> selectedSlice;
    for (const auto &slice : slices)
    {
        const auto isSelectedAsNu = RecoHelper::IsSliceSelectedAsNu(slice, slicesToPFParticles);

        if (isSelectedAsNu)
        {
            if (selectedSlice.isNonnull())
                throw cet::exception("AnalysisHelper::PopulateEventRecoSliceInfo") << " - Multiple slices selected as neutrino." << std::endl;

            selectedSlice = slice;
        }
    }

    const auto hasSelectedSlice = selectedSlice.isNonnull();

    // Get the topological scores of the slices
    std::vector<float> sliceTopologicalScores;
    for (const auto &slice : slices)
    {
        float topologicalScore = -1.f;

        // The score of the slice is held in the metadata of the primary PFParticles in that slice
        for (const auto &pfParticle : CollectionHelper::GetManyAssociated(slice, slicesToPFParticles))
        {
            if (pfParticle->IsPrimary())
            {
                const auto metadata = CollectionHelper::GetSingleAssociated(pfParticle, pfParticlesToMetadata);

                try
                {
                    topologicalScore = RecoHelper::GetTopologicalScore(metadata);
                    break;
                }
                catch (const cet::exception &) {}
            }
        }

        // ATTN here we choose to possibly push a dummy value to keep the vectors in sync
        sliceTopologicalScores.push_back(topologicalScore);

        if (slice == selectedSlice)
            reco.selectedTopologicalScore.Set(topologicalScore);
    }
    
    reco.nSlices.Set(slices.size());
    reco.hasSelectedSlice.Set(hasSelectedSlice);
    reco.sliceTopologicalScores.Set(sliceTopologicalScores);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void EventFactory::PopulateEventTruthSliceInfo(const art::Event &event, const Config &config, Event::Truth &truth)
{
    // ATTN I could probably do this much more efficiently, for now I'm reusing an old implementation
    const auto slices = CollectionHelper::GetCollection<recob::Slice>(event, config.SliceLabel());
    const auto slicesToHits = CollectionHelper::GetAssociation<recob::Slice, recob::Hit>(event, config.SliceLabel());
    const auto pfParticlesToSlices = CollectionHelper::GetAssociation<recob::PFParticle, recob::Slice>(event, config.PFParticleLabel(), config.SliceLabel());
    const auto slicesToPFParticles = CollectionHelper::GetReversedAssociation(pfParticlesToSlices);
    const auto nuMCTruth = TruthHelper::GetNeutrinoMCTruth(event, config.MCTruthLabel());
    const auto nuMCParticles = TruthHelper::GetMCParticlesFromMCTruth(event, config.MCTruthLabel(), config.MCParticleLabel(), nuMCTruth);
    const auto hitsToIsNuInduced = BacktrackHelper::GetHitsToIsNuInducedMap(event, config.HitLabel(), config.MCParticleLabel(), config.BacktrackerLabel(), nuMCParticles);

    // Get the required mapping from slice to "is selected as neutrino"
    SlicesToBool sliceToIsSelectedAsNu;
    for (const auto &slice : slices)
    {
        const auto isSelectedAsNu = RecoHelper::IsSliceSelectedAsNu(slice, slicesToPFParticles);

        if (!sliceToIsSelectedAsNu.emplace(slice, isSelectedAsNu).second)
            throw cet::exception("AnalysisHelper::PopulateEventRecoSliceInfo") << " - Repeated slices." << std::endl;
    }

    const auto sliceMetadata = BacktrackHelper::SliceMetadata(slices, sliceToIsSelectedAsNu, slicesToHits, hitsToIsNuInduced);
    const auto nNeutrinoHits = sliceMetadata.GetTotalNumberOfNuInducedHits();

    // Get the purities and completenesses of the slices
    std::vector<float> slicePurities, sliceCompletenesses;
    std::vector<bool> sliceIsSelectedAsNu;
    for (const auto &slice : slices)
    {
        const auto nHitsInSlice = sliceMetadata.GetNumberOfHits(slice);

        slicePurities.push_back(nHitsInSlice == 0 ? -std::numeric_limits<float>::max() : sliceMetadata.GetPurity(slice));
        sliceCompletenesses.push_back(nNeutrinoHits == 0 ? -std::numeric_limits<float>::max() : sliceMetadata.GetCompleteness(slice));
        sliceIsSelectedAsNu.push_back(sliceToIsSelectedAsNu.at(slice));
    }

    truth.slicePurities.Set(slicePurities);
    truth.sliceCompletenesses.Set(sliceCompletenesses);
    truth.sliceIsSelectedAsNu.Set(sliceIsSelectedAsNu);
}

} // namespace ubcc1pi
