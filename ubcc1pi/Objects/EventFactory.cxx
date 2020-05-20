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

    // Get the event weight information
    if (config.GetEventWeights())
    {
        const auto mcEventWeights = CollectionHelper::GetCollection<evwgh::MCEventWeight>(event, config.MCEventWeightLabel());

        if (!mcEventWeights.empty())
        {
            const auto &mcEventWeight = mcEventWeights.front();

            for (const auto &entry : mcEventWeight->fWeight)
            {
                const auto &label = entry.first;
                const auto &weights = entry.second;

                if (label == "splines_general_Spline")
                {
                    if (truth.splineEventWeight.IsSet())
                        throw cet::exception("EventFactory::PopulateEventTruthInfo") << " - Found multiple spline event weight entries." << std::endl;

                    if (weights.size() > 1)
                        throw cet::exception("EventFactory::PopulateEventTruthInfo") << " - Found multiple spline event weights." << std::endl;

                    if (!weights.empty())
                        truth.splineEventWeight.Set(weights.front());
                }
                
                if (label == "TunedCentralValue_Genie")
                {
                    if (truth.genieTuneEventWeight.IsSet())
                        throw cet::exception("EventFactory::PopulateEventTruthInfo") << " - Found multiple genie tune event weight entries." << std::endl;

                    if (weights.size() > 1)
                        throw cet::exception("EventFactory::PopulateEventTruthInfo") << " - Found multiple genie tune event weights." << std::endl;

                    if (!weights.empty())
                        truth.genieTuneEventWeight.Set(weights.front());
                }
            }
        }
    }

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
    const auto mcParticleMap = TruthHelper::GetMCParticleMap(interaction.GetAllMCParticles());

    for (const auto &mcParticle : finalStateMCParticles)
    {
        truth.particles.emplace_back();
        auto &particle = truth.particles.back();

        EventFactory::PopulateEventTruthParticleInfo(event, config, mcParticle, mcParticleToHitWeights, mcParticleMap, particle);
    }
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void EventFactory::PopulateEventTruthParticleInfo(const art::Event &event, const Config &config, const art::Ptr<simb::MCParticle> &mcParticle, const MCParticlesToHitWeights &mcParticleToHitWeights, const MCParticleMap &mcParticleMap, Event::Truth::Particle &particle)
{
    particle.pdgCode.Set(mcParticle->PdgCode());

    const auto start = mcParticle->Position().Vect();
    particle.startX.Set(start.X());
    particle.startY.Set(start.Y());
    particle.startZ.Set(start.Z());
    
    const auto end = mcParticle->EndPosition().Vect();
    particle.endX.Set(end.X());
    particle.endY.Set(end.Y());
    particle.endZ.Set(end.Z());

    particle.momentumX.Set(mcParticle->Px());
    particle.momentumY.Set(mcParticle->Py());
    particle.momentumZ.Set(mcParticle->Pz());
    particle.momentum.Set(mcParticle->P());
    particle.energy.Set(mcParticle->E());
    particle.mass.Set(mcParticle->Mass());

    particle.hitWeightU.Set(BacktrackHelper::GetHitWeightInView(mcParticle, mcParticleToHitWeights, geo::kU));
    particle.hitWeightV.Set(BacktrackHelper::GetHitWeightInView(mcParticle, mcParticleToHitWeights, geo::kV));
    particle.hitWeightW.Set(BacktrackHelper::GetHitWeightInView(mcParticle, mcParticleToHitWeights, geo::kW));
    
    // Follow the hierarchy information
    TruthHelper::ScatterVector scatters;
    art::Ptr<simb::MCParticle> scatteredParticle;
    TruthHelper::FollowScatters(mcParticle, mcParticleMap, scatters, scatteredParticle);
    const auto endState = TruthHelper::GetEndState(scatteredParticle, mcParticleMap);

    // Output the scattering info
    unsigned int nElasticScatters = 0;
    unsigned int nInelasticScatters = 0;
    
    std::vector<float> scatterCosThetas, scatterMomentumFracsLost;
    std::vector<int> scatterIsElastic;
    
    for (const auto &scatter : scatters)
    {
        nElasticScatters += scatter.m_isElastic ? 1 : 0;
        nInelasticScatters += !scatter.m_isElastic ? 1 : 0;

        scatterCosThetas.push_back(scatter.GetScatteringCosTheta());
        scatterMomentumFracsLost.push_back(scatter.GetMomentumFractionLost());
        scatterIsElastic.push_back(scatter.m_isElastic ? 1 : 0);
    }

    particle.nElasticScatters.Set(nElasticScatters);
    particle.nInelasticScatters.Set(nInelasticScatters);
    particle.scatterCosThetas.Set(scatterCosThetas);
    particle.scatterMomentumFracsLost.Set(scatterMomentumFracsLost);
    particle.scatterIsElastic.Set(scatterIsElastic);

    // The momentum of the particle after all scatters
    particle.scatteredMomentum.Set(scatteredParticle->Momentum().Vect().Mag());

    // The momentum of the particle before the end-state
    const auto endMomentum = endState.m_finalParticleMomentum.Mag();
    particle.endMomentum.Set(endMomentum);
    particle.isStopping.Set(endMomentum <= std::numeric_limits<float>::epsilon());

    // The end-state
    particle.endState.Set(endState.m_type);
   
    // The hit weight of the end-state
    const auto endStateMCParticleToHitWeights = CollectionHelper::GetReversedAssociation(BacktrackHelper::GetHitToMCParticleWeightMap(event, config.MCParticleLabel(), config.BacktrackerLabel(), endState.m_products));

    float endStateProductsHitWeightU = 0;
    float endStateProductsHitWeightV = 0;
    float endStateProductsHitWeightW = 0;
    for (const auto &particle : endState.m_products)
    {
        endStateProductsHitWeightU += BacktrackHelper::GetHitWeightInView(particle, endStateMCParticleToHitWeights, geo::kU);
        endStateProductsHitWeightV += BacktrackHelper::GetHitWeightInView(particle, endStateMCParticleToHitWeights, geo::kV);
        endStateProductsHitWeightW += BacktrackHelper::GetHitWeightInView(particle, endStateMCParticleToHitWeights, geo::kW);
    }

    particle.endStateProductsHitWeightU.Set(endStateProductsHitWeightU);
    particle.endStateProductsHitWeightV.Set(endStateProductsHitWeightV);
    particle.endStateProductsHitWeightW.Set(endStateProductsHitWeightW);
}

// -----------------------------------------------------------------------------------------------------------------------------------------
    
void EventFactory::PopulateEventRecoInfo(const art::Event &event, const Config &config, const MCParticleVector &finalStateMCParticles, Event::Reco &reco)
{
    // Fill the slice information
    EventFactory::PopulateEventRecoSliceInfo(event, config, reco);

    // Start by assuming the event doesn't pass the selection
    reco.passesCCInclusive.Set(false);

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
    const auto pfParticleMap = RecoHelper::GetPFParticleMap(allPFParticles);
    const auto pfParticleToT0s = CollectionHelper::GetAssociation<recob::PFParticle, anab::T0>(event, config.PFParticleLabel(), config.CCInclusiveLabel());
    const auto pfParticleToMetadata = CollectionHelper::GetAssociation<recob::PFParticle, larpandoraobj::PFParticleMetadata>(event, config.PFParticleLabel(), config.PFParticleLabel());
    const auto pfParticleToTracks = CollectionHelper::GetAssociation<recob::PFParticle, recob::Track>(event, config.PFParticleLabel(), config.TrackLabel());
    const auto pfParticleToHits = CollectionHelper::GetAssociationViaCollection<recob::PFParticle, recob::Cluster, recob::Hit>(event, config.PFParticleLabel(), config.PFParticleLabel(), config.PFParticleLabel());
    const auto trackToPIDs = CollectionHelper::GetAssociation<recob::Track, anab::ParticleID>(event, config.TrackLabel(), config.PIDLabel());
    const auto trackToCalorimetries = CollectionHelper::GetAssociation<recob::Track, anab::Calorimetry>(event, config.TrackLabel(), config.CalorimetryLabel());
    const auto spacePoints = CollectionHelper::GetCollection<recob::SpacePoint>(event, config.SpacePointLabel());

    for (const auto &pfParticle : finalStatePFParticles)
    {
        reco.particles.emplace_back();
        auto &particle = reco.particles.back();

        EventFactory::PopulateEventRecoParticleInfo(config, pfParticle, pfParticleMap, pfParticleToT0s, pfParticleToMetadata, pfParticleToHits, pfParticleToTracks, trackToPIDs, trackToCalorimetries, spacePoints, finalStateMCParticles, pBacktrackerData, nuVertex, particle);

        // If we have a muon candidate from the CC inclusive selection, then this event passed
        if (particle.isCCInclusiveMuonCandidate())
            reco.passesCCInclusive.Set(true);
    }
}   

// -----------------------------------------------------------------------------------------------------------------------------------------

void EventFactory::PopulateEventRecoParticleInfo(const Config &config, const art::Ptr<recob::PFParticle> &pfParticle, const PFParticleMap &pfParticleMap, const PFParticleToT0s &pfParticleToT0s, const PFParticleToMetadata &pfParticleToMetadata, const PFParticleToHits &pfParticleToHits, const PFParticleToTracks &pfParticleToTracks, const TrackToPIDs &trackToPIDs, const TrackToCalorimetries &trackToCalorimetries, const SpacePointVector &spacePoints, const MCParticleVector &finalStateMCParticles, const std::shared_ptr<BacktrackHelper::BacktrackerData> &pBacktrackerData, const TVector3 &nuVertex, Event::Reco::Particle &particle)
{
    if (!config.HasTruthInfo() && pBacktrackerData)
        throw cet::exception("EventFactory::PopulateEventRecoParticleInfo") << " - Supplied with backtracker data even though there should be no truth info." << std::endl;
  
    particle.isCCInclusiveMuonCandidate.Set(EventFactory::IsCCInclusiveMuonCandidate(pfParticle, pfParticleToT0s));

    EventFactory::PopulateEventRecoParticlePatRecInfo(pfParticle, pfParticleMap, pfParticleToMetadata, pfParticleToHits, particle);

    try
    {
        const auto track = CollectionHelper::GetSingleAssociated(pfParticle, pfParticleToTracks);
        EventFactory::PopulateEventRecoParticleTrackInfo(config, track, spacePoints, nuVertex, particle);

        try
        {
            const auto pid = CollectionHelper::GetSingleAssociated(track, trackToPIDs);
            EventFactory::PopulateEventRecoParticlePIDInfo(config, pid, particle);
        }
        catch (const cet::exception &) {}
        
        const auto calos = CollectionHelper::GetManyAssociated(track, trackToCalorimetries);
        EventFactory::PopulateEventRecoParticleCalorimetryInfo(config, calos, particle);
    }
    catch (const cet::exception &) {}

    // TODO consider refactoring below into separate function
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

bool EventFactory::IsCCInclusiveMuonCandidate(const art::Ptr<recob::PFParticle> &pfParticle, const PFParticleToT0s &pfParticleToT0s)
{
    // ATTN the CC inclusve producer will associate a T0 to the muon candidate PFParticle if the event has passed the selection
    return CollectionHelper::HasAssociated(pfParticle, pfParticleToT0s);
}

// -----------------------------------------------------------------------------------------------------------------------------------------
    
void EventFactory::PopulateEventRecoParticlePatRecInfo(const art::Ptr<recob::PFParticle> &pfParticle, const PFParticleMap &pfParticleMap, const PFParticleToMetadata &pfParticleToMetadata, const PFParticleToHits &pfParticleToHits, Event::Reco::Particle &particle)
{
    particle.pdgCode.Set(pfParticle->PdgCode());

    const auto hits = CollectionHelper::GetManyAssociated(pfParticle, pfParticleToHits);
    particle.nHitsU.Set(RecoHelper::CountHitsInView(hits, geo::kU));
    particle.nHitsV.Set(RecoHelper::CountHitsInView(hits, geo::kV));
    particle.nHitsW.Set(RecoHelper::CountHitsInView(hits, geo::kW));

    particle.nDaughters.Set(pfParticle->NumDaughters());

    // Get the descendent PFParticles (not including the primary)
    PFParticleVector descendents;
    for (const auto &downstreamParticle : RecoHelper::GetDownstreamParticles(pfParticle, pfParticleMap))
        if (downstreamParticle != pfParticle) descendents.push_back(downstreamParticle);

    particle.nDescendents.Set(descendents.size());

    unsigned int nDescendentHitsU = 0;
    unsigned int nDescendentHitsV = 0;
    unsigned int nDescendentHitsW = 0;
    unsigned int nHitsInLargestDescendent = 0;
    for (const auto &descendent : descendents)
    {
        const auto descendentHits = CollectionHelper::GetManyAssociated(descendent, pfParticleToHits);

        const auto nHitsU = RecoHelper::CountHitsInView(descendentHits, geo::kU);
        const auto nHitsV = RecoHelper::CountHitsInView(descendentHits, geo::kV);
        const auto nHitsW = RecoHelper::CountHitsInView(descendentHits, geo::kW);

        const auto nHitsTot = nHitsU + nHitsV + nHitsW;
        nHitsInLargestDescendent = std::max(nHitsInLargestDescendent, nHitsTot);

        nDescendentHitsU += nHitsU;
        nDescendentHitsV += nHitsV;
        nDescendentHitsW += nHitsW;
    }

    particle.nDescendentHitsU.Set(nDescendentHitsU);
    particle.nDescendentHitsV.Set(nDescendentHitsV);
    particle.nDescendentHitsW.Set(nDescendentHitsW);
    particle.nHitsInLargestDescendent.Set(nHitsInLargestDescendent);
    
    const auto metadata = CollectionHelper::GetSingleAssociated(pfParticle, pfParticleToMetadata);
    const auto trackScore = RecoHelper::GetTrackScore(metadata);

    // ATTN track score -1.f means that it wasn't able to be calculated
    if (trackScore > -0.5f)
        particle.trackScore.Set(trackScore);
}

// -----------------------------------------------------------------------------------------------------------------------------------------
        
void EventFactory::PopulateEventRecoParticleTrackInfo(const Config &config, const art::Ptr<recob::Track> &track, const SpacePointVector &spacePoints, const TVector3 &nuVertex, Event::Reco::Particle &particle)
{
    const auto pSpaceChargeService = RecoHelper::GetSpaceChargeService();

    const auto start = RecoHelper::CorrectForSpaceCharge(TVector3(track->Start().X(), track->Start().Y(), track->Start().Z()), pSpaceChargeService);
    particle.startX.Set(start.X());
    particle.startY.Set(start.Y());
    particle.startZ.Set(start.Z());
    
    const auto end = RecoHelper::CorrectForSpaceCharge(TVector3(track->End().X(), track->End().Y(), track->End().Z()), pSpaceChargeService);
    particle.endX.Set(end.X());
    particle.endY.Set(end.Y());
    particle.endZ.Set(end.Z());

    const auto direction = TVector3(track->StartDirection().X(), track->StartDirection().Y(), track->StartDirection().Z());
    particle.directionX.Set(direction.X());
    particle.directionY.Set(direction.Y());
    particle.directionZ.Set(direction.Z());

    particle.yzAngle.Set(std::atan2(direction.Z(), direction.Y()));
    particle.xyAngle.Set(std::atan2(direction.Y(), direction.X()));
    particle.xzAngle.Set(std::atan2(direction.Z(), direction.X()));

    particle.length.Set(track->Length());
    particle.range.Set(RecoHelper::GetRange(track, pSpaceChargeService));

    float transverseVertexDist = -std::numeric_limits<float>::max();
    float longitudinalVertexDist = -std::numeric_limits<float>::max();
    RecoHelper::GetDistanceToPoint(track, nuVertex, pSpaceChargeService, transverseVertexDist, longitudinalVertexDist);

    particle.transverseVertexDist.Set(transverseVertexDist);
    particle.longitudinalVertexDist.Set(longitudinalVertexDist);
    particle.wiggliness.Set(RecoHelper::GetTrackWiggliness(track));
    
    particle.nSpacePointsNearEnd.Set(RecoHelper::CountSpacePointsNearTrackEnd(track, spacePoints, config.TrackEndSpacePointDistance(), pSpaceChargeService));
}

// -----------------------------------------------------------------------------------------------------------------------------------------
            
void EventFactory::PopulateEventRecoParticlePIDInfo(const Config &config, const art::Ptr<anab::ParticleID> &pid, Event::Reco::Particle &particle)
{
    EventFactory::SetBraggLikelihood(pid, 13, geo::kU, anab::kForward, particle.likelihoodForwardMuonU);
    EventFactory::SetBraggLikelihood(pid, 13, geo::kV, anab::kForward, particle.likelihoodForwardMuonV);
    EventFactory::SetBraggLikelihood(pid, 13, geo::kW, anab::kForward, particle.likelihoodForwardMuonW);
    EventFactory::SetBraggLikelihood(pid, 13, geo::kU, anab::kBackward, particle.likelihoodBackwardMuonU);
    EventFactory::SetBraggLikelihood(pid, 13, geo::kV, anab::kBackward, particle.likelihoodBackwardMuonV);
    EventFactory::SetBraggLikelihood(pid, 13, geo::kW, anab::kBackward, particle.likelihoodBackwardMuonW);
    
    EventFactory::SetBraggLikelihood(pid, 211, geo::kU, anab::kForward, particle.likelihoodForwardPionU);
    EventFactory::SetBraggLikelihood(pid, 211, geo::kV, anab::kForward, particle.likelihoodForwardPionV);
    EventFactory::SetBraggLikelihood(pid, 211, geo::kW, anab::kForward, particle.likelihoodForwardPionW);
    EventFactory::SetBraggLikelihood(pid, 211, geo::kU, anab::kBackward, particle.likelihoodBackwardPionU);
    EventFactory::SetBraggLikelihood(pid, 211, geo::kV, anab::kBackward, particle.likelihoodBackwardPionV);
    EventFactory::SetBraggLikelihood(pid, 211, geo::kW, anab::kBackward, particle.likelihoodBackwardPionW);
    
    EventFactory::SetBraggLikelihood(pid, 2212, geo::kU, anab::kForward, particle.likelihoodForwardProtonU);
    EventFactory::SetBraggLikelihood(pid, 2212, geo::kV, anab::kForward, particle.likelihoodForwardProtonV);
    EventFactory::SetBraggLikelihood(pid, 2212, geo::kW, anab::kForward, particle.likelihoodForwardProtonW);
    EventFactory::SetBraggLikelihood(pid, 2212, geo::kU, anab::kBackward, particle.likelihoodBackwardProtonU);
    EventFactory::SetBraggLikelihood(pid, 2212, geo::kV, anab::kBackward, particle.likelihoodBackwardProtonV);
    EventFactory::SetBraggLikelihood(pid, 2212, geo::kW, anab::kBackward, particle.likelihoodBackwardProtonW);
    
    EventFactory::SetBraggLikelihood(pid, 0, geo::kU, anab::kForward, particle.likelihoodMIPU);
    EventFactory::SetBraggLikelihood(pid, 0, geo::kV, anab::kForward, particle.likelihoodMIPV);
    EventFactory::SetBraggLikelihood(pid, 0, geo::kW, anab::kForward, particle.likelihoodMIPW);
    
    // We need the yz angle to do the rest
    if (!particle.yzAngle.IsSet())
        return;

    const auto yzAngle = particle.yzAngle();
    const auto sin2AngleThreshold = config.Sin2AngleThreshold();
    const auto nHitsU = particle.nHitsU();
    const auto nHitsV = particle.nHitsV();

    EventFactory::SetBraggLikelihood(pid, 13, anab::kForward, yzAngle, sin2AngleThreshold, nHitsU, nHitsV, particle.likelihoodForwardMuon);
    EventFactory::SetBraggLikelihood(pid, 13, anab::kBackward, yzAngle, sin2AngleThreshold, nHitsU, nHitsV, particle.likelihoodBackwardMuon);
    
    EventFactory::SetBraggLikelihood(pid, 211, anab::kForward, yzAngle, sin2AngleThreshold, nHitsU, nHitsV, particle.likelihoodForwardPion);
    EventFactory::SetBraggLikelihood(pid, 211, anab::kBackward, yzAngle, sin2AngleThreshold, nHitsU, nHitsV, particle.likelihoodBackwardPion);
    
    EventFactory::SetBraggLikelihood(pid, 2212, anab::kForward, yzAngle, sin2AngleThreshold, nHitsU, nHitsV, particle.likelihoodForwardProton);
    EventFactory::SetBraggLikelihood(pid, 2212, anab::kBackward, yzAngle, sin2AngleThreshold, nHitsU, nHitsV, particle.likelihoodBackwardProton);
    
    EventFactory::SetBraggLikelihood(pid, 0, anab::kForward, yzAngle, sin2AngleThreshold, nHitsU, nHitsV, particle.likelihoodMIP);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void EventFactory::SetBraggLikelihood(const art::Ptr<anab::ParticleID> &pid, const int &pdg, const geo::View_t &view, const anab::kTrackDir &dir, Member<float> &member)
{
    const auto likelihood = RecoHelper::GetBraggLikelihood(pid, pdg, view, dir);

    // When likelihood isn't available the default is -floatmax
    if (likelihood > -1.f)
        member.Set(likelihood);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void EventFactory::SetBraggLikelihood(const art::Ptr<anab::ParticleID> &pid, const int &pdg, const anab::kTrackDir &dir, const float yzAngle, const float sin2AngleThreshold, const unsigned int nHitsU, const unsigned int nHitsV, Member<float> &member)
{
    const auto likelihood = RecoHelper::GetBraggLikelihood(pid, pdg, dir, yzAngle, sin2AngleThreshold, nHitsU, nHitsV);

    // When likelihood isn't available the default is -floatmax
    if (likelihood > -1.f)
        member.Set(likelihood);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void EventFactory::PopulateEventRecoParticleCalorimetryInfo(const Config &config, const CalorimetryVector &calos, Event::Reco::Particle &particle)
{
    float uWeight = 0.f;
    float vWeight = 0.f;
    for (const auto &calo : calos)
    {        
        const auto dedxPerHit = calo->dEdx();
        const auto residualRangePerHit = calo->ResidualRange();
        const auto truncatedMeandEdx = RecoHelper::GetTruncatedMeandEdxAtTrackStart(dedxPerHit, residualRangePerHit, config.NHitsToSkip(), config.LengthFraction());

        Member<float> *pMember = nullptr;

        switch (calo->PlaneID().Plane)
        {
            case geo::kU:
                pMember = &particle.truncatedMeandEdxU;
                uWeight = std::max(0.f, (dedxPerHit.size() * config.LengthFraction()) - config.NHitsToSkip());
                break;
            case geo::kV:
                pMember = &particle.truncatedMeandEdxV;
                vWeight = std::max(0.f, (dedxPerHit.size() * config.LengthFraction()) - config.NHitsToSkip());
                break;
            case geo::kW:
                pMember = &particle.truncatedMeandEdxW;
                break;
            default: break;
        }

        if (pMember && truncatedMeandEdx > -1.f)
            pMember->Set(truncatedMeandEdx);
    }

    const auto yzAngle = particle.yzAngle();
    const bool isTrackAlongWWire = (std::pow(std::sin(yzAngle), 2) < config.Sin2AngleThreshold());
    if (!isTrackAlongWWire && particle.truncatedMeandEdxW.IsSet())
    {
        particle.truncatedMeandEdx.Set(particle.truncatedMeandEdxW());
    }
    else
    {
        const auto hasU = particle.truncatedMeandEdxU.IsSet() && uWeight > 0.f;
        const auto hasV = particle.truncatedMeandEdxV.IsSet() && vWeight > 0.f;

        if (hasU || hasV)
        {
            float truncatedMeandEdx = 0.f;
            if (hasU) truncatedMeandEdx += particle.truncatedMeandEdxU() * uWeight;
            if (hasV) truncatedMeandEdx += particle.truncatedMeandEdxV() * vWeight;

            truncatedMeandEdx /= (uWeight + vWeight);
            particle.truncatedMeandEdx.Set(truncatedMeandEdx);
        }
    }
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
    art::Ptr<recob::Slice> selectedSlice;
    std::vector<bool> sliceIsSelectedAsNu;
    for (const auto &slice : slices)
    {
        const auto isSelectedAsNu = RecoHelper::IsSliceSelectedAsNu(slice, slicesToPFParticles);
        sliceIsSelectedAsNu.push_back(isSelectedAsNu);
        
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
        // Define the dummy value for slices without a topological score
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

        // ATTN here we can possibly push a dummy value to keep the vectors in sync
        sliceTopologicalScores.push_back(topologicalScore);

        if (slice == selectedSlice)
            reco.selectedTopologicalScore.Set(topologicalScore);
    }
    
    reco.nSlices.Set(slices.size());
    reco.hasSelectedSlice.Set(hasSelectedSlice);
    reco.sliceTopologicalScores.Set(sliceTopologicalScores);
    reco.sliceIsSelectedAsNu.Set(sliceIsSelectedAsNu);
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
    for (const auto &slice : slices)
    {
        const auto nHitsInSlice = sliceMetadata.GetNumberOfHits(slice);

        slicePurities.push_back(nHitsInSlice == 0 ? -std::numeric_limits<float>::max() : sliceMetadata.GetPurity(slice));
        sliceCompletenesses.push_back(nNeutrinoHits == 0 ? -std::numeric_limits<float>::max() : sliceMetadata.GetCompleteness(slice));
    }

    truth.slicePurities.Set(slicePurities);
    truth.sliceCompletenesses.Set(sliceCompletenesses);
}

} // namespace ubcc1pi
