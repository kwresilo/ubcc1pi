/**
 *  @file  ubcc1pi/Analyzers/EnergyEstimator_module.cc
 *
 *  @brief The implementation file for the PID study analyzer.
 */

#include "ubcc1pi/Analyzers/EnergyEstimator.h"

#include "ubcc1pi/Helpers/TruthHelper.h"
#include "ubcc1pi/Helpers/AnalysisHelper.h"
#include "ubcc1pi/Helpers/RecoHelper.h"
#include "ubcc1pi/Helpers/DebugHelper.h"


namespace ubcc1pi
{

EnergyEstimator::EnergyEstimator(const art::EDAnalyzer::Table<Config> &config) :
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

    m_pParticleTree->Branch("hasMatchedMCParticle", &m_outputParticle.m_hasMatchedMCParticle);
    m_pParticleTree->Branch("trueMatchPurity", &m_outputParticle.m_trueMatchPurity);
    m_pParticleTree->Branch("trueMatchCompleteness", &m_outputParticle.m_trueMatchCompleteness);
    m_pParticleTree->Branch("truePdgCode", &m_outputParticle.m_truePdgCode);
    m_pParticleTree->Branch("trueMomentum", &m_outputParticle.m_trueMomentum);
    m_pParticleTree->Branch("trueEndMomentum", &m_outputParticle.m_trueEndMomentum);
    m_pParticleTree->Branch("trueKE", &m_outputParticle.m_trueKE);
    m_pParticleTree->Branch("isStopping", &m_outputParticle.m_isStopping);
    m_pParticleTree->Branch("trueStart", &m_outputParticle.m_trueStart);
    m_pParticleTree->Branch("trueEnd", &m_outputParticle.m_trueEnd);
    m_pParticleTree->Branch("isContained", &m_outputParticle.m_isContained);
    m_pParticleTree->Branch("nScatters", &m_outputParticle.m_nScatters);
    m_pParticleTree->Branch("isGolden", &m_outputParticle.m_isGolden);
    m_pParticleTree->Branch("trueRange", &m_outputParticle.m_trueRange);
    m_pParticleTree->Branch("trueEndPointDist", &m_outputParticle.m_trueEndPointDist);

    m_pParticleTree->Branch("hasTrack", &m_outputParticle.m_hasTrack);
    m_pParticleTree->Branch("start", &m_outputParticle.m_start);
    m_pParticleTree->Branch("end", &m_outputParticle.m_end);
    m_pParticleTree->Branch("length", &m_outputParticle.m_length);
    m_pParticleTree->Branch("range", &m_outputParticle.m_range);
    m_pParticleTree->Branch("endPointDist", &m_outputParticle.m_endPointDist);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void EnergyEstimator::analyze(const art::Event &event)
{
    // Extract the labels
    const auto mcTruthLabel = m_config().MCTruthLabel();
    const auto mcParticleLabel = m_config().MCParticleLabel();
    const auto backtrackerLabel = m_config().BacktrackerLabel();
    const auto pfParticleLabel = m_config().PFParticleLabel();
    const auto trackLabel = m_config().TrackLabel();
    const auto pidLabel = m_config().PIDLabel();
    const auto calorimetryLabel = m_config().CalorimetryLabel();
    
    // Get the neutrino final state PFParticles
    const auto allPFParticles = CollectionHelper::GetCollection<recob::PFParticle>(event, pfParticleLabel);
    const auto pfParticleMap = RecoHelper::GetPFParticleMap(allPFParticles);
    const auto pfParticles = RecoHelper::GetNeutrinoFinalStates(allPFParticles);

    // Get the associations from PFParticles to tracks
    const auto pfpToTrack = CollectionHelper::GetAssociation<recob::PFParticle, recob::Track>(event, pfParticleLabel, trackLabel);
    
    // Get the reco-true matching information
    const auto backtrackerData = AnalysisHelper::GetBacktrackerData(event, mcTruthLabel, mcParticleLabel, backtrackerLabel, pfParticleLabel);

    const TruthHelper::Interaction interaction(event, mcTruthLabel, mcParticleLabel);
    const auto mcParticleMap = TruthHelper::GetMCParticleMap(interaction.GetAllMCParticles());

    // Get the space-charge service
    const auto pSpaceChargeService = RecoHelper::GetSpaceChargeService();

    m_outputParticle.m_run = event.run();
    m_outputParticle.m_subRun = event.subRun();
    m_outputParticle.m_event = event.event();
    m_outputParticle.m_isSignal = AnalysisHelper::IsCC1PiSignal(interaction);

    for (const auto &pfParticle : pfParticles)
    {
        this->ResetParticleTree();

        // Get the matched MCParticle
        art::Ptr<simb::MCParticle> matchedMCParticle;
        try
        {
            matchedMCParticle = backtrackerData.GetBestMatchedMCParticle(pfParticle); 
        }
        catch (const cet::exception &)
        {
            continue;
        }
    
        m_outputParticle.m_hasMatchedMCParticle = true;
        m_outputParticle.m_trueMatchPurity = backtrackerData.GetMatchPurity(pfParticle, matchedMCParticle);
        m_outputParticle.m_trueMatchCompleteness = backtrackerData.GetMatchCompleteness(pfParticle, matchedMCParticle);
        m_outputParticle.m_truePdgCode = matchedMCParticle->PdgCode();

        m_outputParticle.m_trueMomentum = matchedMCParticle->Momentum().Vect().Mag();
        m_outputParticle.m_trueEndMomentum = matchedMCParticle->Momentum(std::max(static_cast<unsigned int>(0), matchedMCParticle->NumberTrajectoryPoints() - 2)).Vect().Mag();
        m_outputParticle.m_trueKE = matchedMCParticle->E() - matchedMCParticle->Mass();
        m_outputParticle.m_isStopping = (m_outputParticle.m_trueEndMomentum <= std::numeric_limits<float>::epsilon());

        m_outputParticle.m_trueStart = matchedMCParticle->Position().Vect();
        m_outputParticle.m_trueEnd = matchedMCParticle->EndPosition().Vect();
        m_outputParticle.m_isContained = AnalysisHelper::IsContained(m_outputParticle.m_trueStart, 5.f) && AnalysisHelper::IsContained(m_outputParticle.m_trueEnd, 5.f);

        m_outputParticle.m_nScatters = 0;
        for (const auto &daughter : TruthHelper::GetDaughters(matchedMCParticle, mcParticleMap))
        {
            if (daughter->Process() == "hadElastic")
                m_outputParticle.m_nScatters++;
        }

        m_outputParticle.m_isGolden = m_outputParticle.m_isStopping && (m_outputParticle.m_nScatters == 0) && m_outputParticle.m_isContained;
        m_outputParticle.m_trueRange = this->GetRange(matchedMCParticle);
        m_outputParticle.m_trueEndPointDist = this->GetEndpointDistance(matchedMCParticle);


        // Now get the reconstructed track information
        art::Ptr<recob::Track> track;
        try
        {
            track = CollectionHelper::GetSingleAssociated(pfParticle, pfpToTrack);
        }
        catch (const cet::exception &)
        {
            continue;
        }

        m_outputParticle.m_hasTrack = true;
        
        m_outputParticle.m_start = RecoHelper::CorrectForSpaceCharge(TVector3(track->Start().X(), track->Start().Y(), track->Start().Z()), pSpaceChargeService);
        m_outputParticle.m_end = RecoHelper::CorrectForSpaceCharge(TVector3(track->End().X(), track->End().Y(), track->End().Z()), pSpaceChargeService);
                        
        m_outputParticle.m_length = track->Length();
        m_outputParticle.m_range = this->GetRange(track, pSpaceChargeService);
        m_outputParticle.m_endPointDist = this->GetEndpointDistance(track, pSpaceChargeService);

        m_pParticleTree->Fill();
    }
}

// -----------------------------------------------------------------------------------------------------------------------------------------

float EnergyEstimator::GetRange(const art::Ptr<simb::MCParticle> &mcParticle) const
{
    float range = 0.f;

    for (unsigned int i = 1; i < mcParticle->NumberTrajectoryPoints(); ++i)
    {
        range += (mcParticle->Position(i).Vect() - mcParticle->Position(i - 1).Vect()).Mag();
    }

    return range;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

float EnergyEstimator::GetEndpointDistance(const art::Ptr<simb::MCParticle> &mcParticle) const
{
    return (mcParticle->Position(mcParticle->NumberTrajectoryPoints() - 1).Vect() - mcParticle->Position(0).Vect()).Mag();
}

// -----------------------------------------------------------------------------------------------------------------------------------------

std::vector<size_t> EnergyEstimator::GetValidPoints(const art::Ptr<recob::Track> &track) const
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

float EnergyEstimator::GetRange(const art::Ptr<recob::Track> &track, const spacecharge::SpaceChargeService::provider_type *const pSpaceChargeService) const
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

float EnergyEstimator::GetEndpointDistance(const art::Ptr<recob::Track> &track, const spacecharge::SpaceChargeService::provider_type *const pSpaceChargeService) const
{
    const auto validPoints = this->GetValidPoints(track);
    if (validPoints.size() < 2)
        return 0.f;

    const auto start = track->LocationAtPoint(validPoints.front());
    const auto end = track->LocationAtPoint(validPoints.back());

    const auto startVect = TVector3(start.X(), start.Y(), start.Z());
    const auto endVect = TVector3(end.X(), end.Y(), end.Z());

    const auto startVectSCE = RecoHelper::CorrectForSpaceCharge(startVect, pSpaceChargeService);
    const auto endVectSCE = RecoHelper::CorrectForSpaceCharge(endVect, pSpaceChargeService);

    return (endVectSCE - startVectSCE).Mag();
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void EnergyEstimator::ResetParticleTree()
{
    m_outputParticle.m_hasMatchedMCParticle = false;
    m_outputParticle.m_trueMatchPurity = -std::numeric_limits<float>::max();
    m_outputParticle.m_trueMatchCompleteness = -std::numeric_limits<float>::max();
    m_outputParticle.m_truePdgCode = -std::numeric_limits<int>::max();
    m_outputParticle.m_trueMomentum = -std::numeric_limits<float>::max();
    m_outputParticle.m_trueEndMomentum = -std::numeric_limits<float>::max();
    m_outputParticle.m_trueKE = -std::numeric_limits<float>::max();
    m_outputParticle.m_isStopping = false;
    m_outputParticle.m_trueStart = TVector3(-std::numeric_limits<float>::max(), -std::numeric_limits<float>::max(), -std::numeric_limits<float>::max());
    m_outputParticle.m_trueEnd = TVector3(-std::numeric_limits<float>::max(), -std::numeric_limits<float>::max(), -std::numeric_limits<float>::max());
    m_outputParticle.m_isContained = false;
    m_outputParticle.m_nScatters = -std::numeric_limits<int>::max();
    m_outputParticle.m_isGolden = false;
    m_outputParticle.m_trueRange = -std::numeric_limits<float>::max();
    m_outputParticle.m_trueEndPointDist = -std::numeric_limits<float>::max();
    m_outputParticle.m_hasTrack = false;
    m_outputParticle.m_start = TVector3(-std::numeric_limits<float>::max(), -std::numeric_limits<float>::max(), -std::numeric_limits<float>::max());
    m_outputParticle.m_end = TVector3(-std::numeric_limits<float>::max(), -std::numeric_limits<float>::max(), -std::numeric_limits<float>::max());
    m_outputParticle.m_length = -std::numeric_limits<float>::max();
    m_outputParticle.m_range = -std::numeric_limits<float>::max();
    m_outputParticle.m_endPointDist = -std::numeric_limits<float>::max();
}

} // namespace ubcc1pi
