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
           
    m_pParticleTree->Branch("generation", &m_outputParticle.m_generation);
    m_pParticleTree->Branch("isPrimary", &m_outputParticle.m_isPrimary);
    m_pParticleTree->Branch("nDaughters", &m_outputParticle.m_nDaughters);
    m_pParticleTree->Branch("integratedRange", &m_outputParticle.m_integratedRange);
    m_pParticleTree->Branch("hasTrack", &m_outputParticle.m_hasTrack);
    m_pParticleTree->Branch("start", &m_outputParticle.m_start);
    m_pParticleTree->Branch("end", &m_outputParticle.m_end);
    m_pParticleTree->Branch("thetaYZ", &m_outputParticle.m_thetaYZ);
    m_pParticleTree->Branch("thetaXZ", &m_outputParticle.m_thetaXZ);
    m_pParticleTree->Branch("thetaXY", &m_outputParticle.m_thetaXY);
    m_pParticleTree->Branch("length", &m_outputParticle.m_length);
    m_pParticleTree->Branch("range", &m_outputParticle.m_range);
    m_pParticleTree->Branch("endPointDist", &m_outputParticle.m_endPointDist);
    m_pParticleTree->Branch("parentCosOpeningAngle", &m_outputParticle.m_parentCosOpeningAngle);
    m_pParticleTree->Branch("daughterCosOpeningAngle", &m_outputParticle.m_daughterCosOpeningAngle);
    
    m_pParticleTree->Branch("hasPid", &m_outputParticle.m_hasPid);
    m_pParticleTree->Branch("muonLikelihoodU", &m_outputParticle.m_muonLikelihoodU);
    m_pParticleTree->Branch("muonLikelihoodV", &m_outputParticle.m_muonLikelihoodV);
    m_pParticleTree->Branch("muonLikelihoodW", &m_outputParticle.m_muonLikelihoodW);
    m_pParticleTree->Branch("pionLikelihoodU", &m_outputParticle.m_pionLikelihoodU);
    m_pParticleTree->Branch("pionLikelihoodV", &m_outputParticle.m_pionLikelihoodV);
    m_pParticleTree->Branch("pionLikelihoodW", &m_outputParticle.m_pionLikelihoodW);
    m_pParticleTree->Branch("protonLikelihoodU", &m_outputParticle.m_protonLikelihoodU);
    m_pParticleTree->Branch("protonLikelihoodV", &m_outputParticle.m_protonLikelihoodV);
    m_pParticleTree->Branch("protonLikelihoodW", &m_outputParticle.m_protonLikelihoodW);
    m_pParticleTree->Branch("mipLikelihoodU", &m_outputParticle.m_mipLikelihoodU);
    m_pParticleTree->Branch("mipLikelihoodV", &m_outputParticle.m_mipLikelihoodV);
    m_pParticleTree->Branch("mipLikelihoodW", &m_outputParticle.m_mipLikelihoodW);
    
    m_pParticleTree->Branch("shouldUseMuon", &m_outputParticle.m_shouldUseMuon);
    m_pParticleTree->Branch("shouldUsePion", &m_outputParticle.m_shouldUsePion);
    m_pParticleTree->Branch("shouldUseProton", &m_outputParticle.m_shouldUseProton);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void EnergyEstimator::analyze(const art::Event &event)
{
    std::cout << "RUN:    " << event.run() << std::endl;
    std::cout << "SUBRUN: " << event.subRun() << std::endl;
    std::cout << "EVENT:  " << event.event() << std::endl;

    // Extract the labels
    const auto mcTruthLabel = m_config().MCTruthLabel();
    const auto mcParticleLabel = m_config().MCParticleLabel();
    const auto backtrackerLabel = m_config().BacktrackerLabel();
    const auto pfParticleLabel = m_config().PFParticleLabel();
    const auto trackLabel = m_config().TrackLabel();
    const auto pidLabel = m_config().PIDLabel();
    //const auto calorimetryLabel = m_config().CalorimetryLabel();
   
    std::cout << "I got the labels" << std::endl;

    // Get the neutrino final state PFParticles
    const auto allPFParticles = CollectionHelper::GetCollection<recob::PFParticle>(event, pfParticleLabel);
    const auto pfParticleMap = RecoHelper::GetPFParticleMap(allPFParticles);
    const auto finalStatePFParticles = RecoHelper::GetNeutrinoFinalStates(allPFParticles);
    PFParticleVector pfParticles;

    for (const auto &pfParticle : finalStatePFParticles)
        RecoHelper::GetDownstreamParticles(pfParticle, pfParticleMap, pfParticles);
    
    std::cout << "I have the PFParticles" << std::endl;

    // Get the associations from PFParticles to tracks and pid
    const auto pfpToTrack = CollectionHelper::GetAssociation<recob::PFParticle, recob::Track>(event, pfParticleLabel, trackLabel);
    const auto trackToPID = CollectionHelper::GetAssociation<recob::Track, anab::ParticleID>(event, trackLabel, pidLabel);
    
    // Get the reco-true matching information
    const auto backtrackerData = AnalysisHelper::GetBacktrackerData(event, mcTruthLabel, mcParticleLabel, backtrackerLabel, pfParticleLabel);

    const TruthHelper::Interaction interaction(event, mcTruthLabel, mcParticleLabel);
    const auto mcParticleMap = TruthHelper::GetMCParticleMap(interaction.GetAllMCParticles());

    // Get the space-charge service
    const auto pSpaceChargeService = RecoHelper::GetSpaceChargeService();
    
    std::cout << "I have the associations" << std::endl;

    m_outputParticle.m_run = event.run();
    m_outputParticle.m_subRun = event.subRun();
    m_outputParticle.m_event = event.event();
    m_outputParticle.m_isSignal = AnalysisHelper::IsCC1PiSignal(interaction);

    std::cout << "Going to output the particles now" << std::endl;
    std::cout << "nPFParticles : " << pfParticles.size() << std::endl;

    for (const auto &pfParticle : pfParticles)
    {
        std::cout << "Resetting particle tree" << std::endl;
        this->ResetParticleTree();
        
        std::cout << "PFParticle -----------------------------------------------" << std::endl;
        std::cout << "  generation:   " << RecoHelper::GetGeneration(pfParticle, pfParticleMap) << std::endl;

        // Get the matched MCParticle
        art::Ptr<simb::MCParticle> matchedMCParticle;
        art::Ptr<recob::PFParticle> topLevelPFParticle = pfParticle;
        try
        {
            // ATTN here we match the top-level PFParticle as that's all we can do right now without extra work
            while (std::find(finalStatePFParticles.begin(), finalStatePFParticles.end(), topLevelPFParticle) == finalStatePFParticles.end())
            {
                topLevelPFParticle = RecoHelper::GetParent(topLevelPFParticle, pfParticleMap);
            }
                
            matchedMCParticle = backtrackerData.GetBestMatchedMCParticle(topLevelPFParticle); 
            m_outputParticle.m_hasMatchedMCParticle = true;
        }
        catch (const cet::exception &)
        {
            m_outputParticle.m_hasMatchedMCParticle = false;
        }
    
        if (m_outputParticle.m_hasMatchedMCParticle)
        {
            m_outputParticle.m_trueMatchPurity = backtrackerData.GetMatchPurity(topLevelPFParticle, matchedMCParticle);
            m_outputParticle.m_trueMatchCompleteness = backtrackerData.GetMatchCompleteness(topLevelPFParticle, matchedMCParticle);
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

            std::cout << "  PDG:               " << m_outputParticle.m_truePdgCode << std::endl;
            std::cout << "  fullHeirarchyMatch:" << std::endl;
            std::cout << "    Purity:          " << m_outputParticle.m_trueMatchPurity << std::endl;
            std::cout << "    Completeness:    " << m_outputParticle.m_trueMatchCompleteness << std::endl;
            std::cout << "  trueRange:         " << m_outputParticle.m_trueRange << std::endl;
            std::cout << "  trueContained:     " << m_outputParticle.m_isContained << std::endl;
        }

        // Get the PFParticle hierarchy info
        m_outputParticle.m_generation = RecoHelper::GetGeneration(pfParticle, pfParticleMap);
        m_outputParticle.m_isPrimary = (m_outputParticle.m_generation == 1);
        m_outputParticle.m_nDaughters = pfParticle->NumDaughters();
        m_outputParticle.m_integratedRange = this->GetIntegratedRange(pfParticle, pfParticleMap, pfpToTrack, pSpaceChargeService);

        // Now get the reconstructed track information
        art::Ptr<recob::Track> track;
        m_outputParticle.m_hasTrack = this->GetTrack(pfParticle, pfpToTrack, track);

        if (m_outputParticle.m_hasTrack)
        {
            m_outputParticle.m_start = RecoHelper::CorrectForSpaceCharge(TVector3(track->Start().X(), track->Start().Y(), track->Start().Z()), pSpaceChargeService);
            m_outputParticle.m_end = RecoHelper::CorrectForSpaceCharge(TVector3(track->End().X(), track->End().Y(), track->End().Z()), pSpaceChargeService);

            const auto dir = track->StartDirection();
            const auto dirVect = TVector3(dir.X(), dir.Y(), dir.Z()).Unit();

            m_outputParticle.m_thetaYZ = std::atan2(dirVect.Y(), dirVect.Z());
            m_outputParticle.m_thetaXZ = std::atan2(dirVect.X(), dirVect.Z());
            m_outputParticle.m_thetaXY = std::atan2(dirVect.X(), dirVect.Y());

            m_outputParticle.m_length = track->Length();
            m_outputParticle.m_range = this->GetRange(track, pSpaceChargeService);
            m_outputParticle.m_endPointDist = this->GetEndpointDistance(track, pSpaceChargeService);

            std::cout << "  recoRange:         " << m_outputParticle.m_range << std::endl;
            std::cout << "    (r-t)/t:         " << (m_outputParticle.m_range - m_outputParticle.m_trueRange) / m_outputParticle.m_trueRange << std::endl;
            
            std::cout << "  recoRangeInt:      " << m_outputParticle.m_integratedRange << std::endl;
            std::cout << "    (r-t)/t:         " << (m_outputParticle.m_integratedRange - m_outputParticle.m_trueRange) / m_outputParticle.m_trueRange << std::endl;

            // Get angle to parent
            if (!m_outputParticle.m_isPrimary)
            {
                const auto parentPFParticle = RecoHelper::GetParent(pfParticle, pfParticleMap);
                art::Ptr<recob::Track> parentTrack;
                if (this->GetTrack(parentPFParticle, pfpToTrack, parentTrack))
                {
                    m_outputParticle.m_parentCosOpeningAngle = this->GetCosOpeningAngle(track, parentTrack);
                }
            }

            std::cout << "  nDaughters:        " << pfParticle->NumDaughters() << std::endl;

            // Get angle to largest daughter
            if (pfParticle->NumDaughters() >= 1)
            {
                art::Ptr<recob::Track> longestDaughterTrack;

                const auto daughterPFParticles = RecoHelper::GetDaughters(pfParticle, pfParticleMap);
                for (const auto &daughter : daughterPFParticles)
                {
                    art::Ptr<recob::Track> daughterTrack;
                    if (this->GetTrack(daughter, pfpToTrack, daughterTrack))
                    {
                        if (longestDaughterTrack.isNull() || daughterTrack->Length() > longestDaughterTrack->Length())
                            longestDaughterTrack = daughterTrack;
                    }
                }

                if (longestDaughterTrack.isNonnull())
                {
                    m_outputParticle.m_daughterCosOpeningAngle = this->GetCosOpeningAngle(longestDaughterTrack, track);
                    std::cout << "  cosDauAngle:       " << m_outputParticle.m_daughterCosOpeningAngle << std::endl;
                }
            }

            if (m_outputParticle.m_isPrimary &&
                    m_outputParticle.m_isContained &&
                    m_outputParticle.m_truePdgCode == 13 && 
                    (m_outputParticle.m_range - m_outputParticle.m_trueRange) / m_outputParticle.m_trueRange < -0.8 &&
                    m_outputParticle.m_trueMatchCompleteness > 0.5f)
            {
                std::cout << "Broken muon" << std::endl;
            }
            
            // Get the PID info
            art::Ptr<anab::ParticleID> pid;
            try
            {
                pid = CollectionHelper::GetSingleAssociated(track, trackToPID);
                m_outputParticle.m_hasPid = true;
            }
            catch (const cet::exception &)
            {
                m_outputParticle.m_hasPid = false;
            }

            if (m_outputParticle.m_hasPid)
            {
                m_outputParticle.m_muonLikelihoodU = RecoHelper::GetBraggLikelihood(pid, 13, geo::kU);
                m_outputParticle.m_muonLikelihoodV = RecoHelper::GetBraggLikelihood(pid, 13, geo::kV);
                m_outputParticle.m_muonLikelihoodW = RecoHelper::GetBraggLikelihood(pid, 13, geo::kW);
                
                m_outputParticle.m_pionLikelihoodU = RecoHelper::GetBraggLikelihood(pid, 211, geo::kU);
                m_outputParticle.m_pionLikelihoodV = RecoHelper::GetBraggLikelihood(pid, 211, geo::kV);
                m_outputParticle.m_pionLikelihoodW = RecoHelper::GetBraggLikelihood(pid, 211, geo::kW);
                
                m_outputParticle.m_protonLikelihoodU = RecoHelper::GetBraggLikelihood(pid, 2212, geo::kU);
                m_outputParticle.m_protonLikelihoodV = RecoHelper::GetBraggLikelihood(pid, 2212, geo::kV);
                m_outputParticle.m_protonLikelihoodW = RecoHelper::GetBraggLikelihood(pid, 2212, geo::kW);
                
                m_outputParticle.m_mipLikelihoodU = RecoHelper::GetBraggLikelihood(pid, 0, geo::kU);
                m_outputParticle.m_mipLikelihoodV = RecoHelper::GetBraggLikelihood(pid, 0, geo::kV);
                m_outputParticle.m_mipLikelihoodW = RecoHelper::GetBraggLikelihood(pid, 0, geo::kW);
            }

            m_outputParticle.m_shouldUseMuon   = this->CheckIfShouldUseForCalculation(pfParticle, pfParticleMap, pfpToTrack, trackToPID, 13);
            m_outputParticle.m_shouldUsePion   = this->CheckIfShouldUseForCalculation(pfParticle, pfParticleMap, pfpToTrack, trackToPID, 211);
            m_outputParticle.m_shouldUseProton = this->CheckIfShouldUseForCalculation(pfParticle, pfParticleMap, pfpToTrack, trackToPID, 2212);
        }
        
        m_pParticleTree->Fill();
    }


    std::cout << "Beans" << std::endl;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

bool EnergyEstimator::GetTrack(const art::Ptr<recob::PFParticle> &pfParticle, const PFParticleToTracks &pfpToTrack, art::Ptr<recob::Track> &track) const
{
    try
    {
        track = CollectionHelper::GetSingleAssociated(pfParticle, pfpToTrack);
    }
    catch (const cet::exception &)
    {
        return false;
    }

    return true;
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
        
float EnergyEstimator::GetIntegratedRange(const art::Ptr<recob::PFParticle> &pfParticle, const PFParticleMap &pfParticleMap, const PFParticleToTracks &pfpToTrack, const spacecharge::SpaceChargeService::provider_type *const pSpaceChargeService) const
{
    float integratedRange = 0.f;

    art::Ptr<recob::PFParticle> nextParticle = pfParticle;
    while (true)
    {
        art::Ptr<recob::Track> track;
        if (this->GetTrack(nextParticle, pfpToTrack, track))
            integratedRange += this->GetRange(track, pSpaceChargeService);

        if (nextParticle->IsPrimary())
            break;

        nextParticle = RecoHelper::GetParent(nextParticle, pfParticleMap);
    }

    return integratedRange;
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

float EnergyEstimator::GetCosOpeningAngle(const art::Ptr<recob::Track> &daughterTrack, const art::Ptr<recob::Track> &parentTrack) const
{
    const auto daughterDir = daughterTrack->StartDirection();
    const auto daughterDirVect = TVector3(daughterDir.X(), daughterDir.Y(), daughterDir.Z()).Unit();

    const auto parentDir = parentTrack->EndDirection();
    const auto parentDirVect = TVector3(parentDir.X(), parentDir.Y(), parentDir.Z()).Unit();

    return daughterDirVect.Dot(parentDirVect);
}

// -----------------------------------------------------------------------------------------------------------------------------------------
            
bool EnergyEstimator::CheckIfShouldUseForCalculation(const art::Ptr<recob::PFParticle> &pfParticle, const PFParticleMap &pfParticleMap, const PFParticleToTracks &pfpToTrack, const TrackToPIDs &trackToPID, unsigned int pdgCode) const
{
    // Get the upstream PFParticles (not including the neutrino) 
    PFParticleVector pfParticlesInHierarchy;
    art::Ptr<recob::PFParticle> nextPFParticle = pfParticle;
    while (!nextPFParticle->IsPrimary())
    {
        if (pfParticle != nextPFParticle)
            pfParticlesInHierarchy.push_back(nextPFParticle);

        nextPFParticle = RecoHelper::GetParent(nextPFParticle, pfParticleMap);
    }

    // Get the downstream PFParticles following all entries in the hierarchy
    RecoHelper::GetDownstreamParticles(pfParticle, pfParticleMap, pfParticlesInHierarchy);

    // Filter out broken tracks
    PFParticleVector filteredCandidates;
    for (const auto &candidate : pfParticlesInHierarchy)
    {
        const auto daughters = RecoHelper::GetDaughters(candidate, pfParticleMap);

        // Daughterless particles pass
        if (daughters.empty())
        {
            filteredCandidates.push_back(candidate);
            continue;    
        }
        
        // Must have a track associated
        art::Ptr<recob::Track> track;
        if (!this->GetTrack(candidate, pfpToTrack, track))
            continue;

        // Get the opening angle with each daughter
        bool passesAngleCut = true;
        for (const auto &daughter : daughters)
        {
            art::Ptr<recob::Track> daughterTrack;
            if (!this->GetTrack(daughter, pfpToTrack, daughterTrack))
                continue;

            // If the particle is very co-linear with any of it's daughters, then don't consider it
            if (this->GetCosOpeningAngle(daughterTrack, track) >= 0.995)
            {
                passesAngleCut = false;
                break;
            }
        }

        if (!passesAngleCut)
            continue;

        art::Ptr<anab::ParticleID> pid;
        try
        {
            pid = CollectionHelper::GetSingleAssociated(track, trackToPID);
        }
        catch (const cet::exception &)
        {
            continue;
        }

        // Cut out any particles that look more like a MIP than the particle they are supposed to be
        const auto likelihoodW = RecoHelper::GetBraggLikelihood(pid, pdgCode, geo::kW);
        const auto mipLikelihoodW = RecoHelper::GetBraggLikelihood(pid, 0, geo::kW);
        if (likelihoodW >= 0.f && mipLikelihoodW >= 0.f)
        {
            if (likelihoodW > mipLikelihoodW)
                filteredCandidates.push_back(candidate);

            continue;
        }
        
        const auto likelihoodV = RecoHelper::GetBraggLikelihood(pid, pdgCode, geo::kV);
        const auto mipLikelihoodV = RecoHelper::GetBraggLikelihood(pid, 0, geo::kV);
        if (likelihoodV >= 0.f && mipLikelihoodV >= 0.f)
        {
            if (likelihoodV > mipLikelihoodV)
                filteredCandidates.push_back(candidate);

            continue;
        }
        
        const auto likelihoodU = RecoHelper::GetBraggLikelihood(pid, pdgCode, geo::kU);
        const auto mipLikelihoodU = RecoHelper::GetBraggLikelihood(pid, 0, geo::kU);
        if (likelihoodU >= 0.f && mipLikelihoodU >= 0.f)
        {
            if (likelihoodU > mipLikelihoodU)
                filteredCandidates.push_back(candidate);

            continue;
        }

        // No PID available
    }

    if (filteredCandidates.empty())
        return false;

    // Select the most upstream candidates where possible
    unsigned int minGeneration = std::numeric_limits<unsigned int>::max();
    PFParticleVector finalCandidates;
    for (const auto &candidate : filteredCandidates)
    {
        const auto generation = RecoHelper::GetGeneration(candidate, pfParticleMap);

        if (generation < minGeneration)
        {
            minGeneration = generation;
            finalCandidates.clear();
        }

        if (generation == minGeneration)
            finalCandidates.push_back(candidate);
    }

    if (finalCandidates.empty())
        throw cet::exception("EnergyEstimator::CheckIfShouldUseForCalculation") << " - Logic error, sanity check failed." << std::endl;

    // Now if there are multiple options of the same generation, choose the longest as a tie-breaker
    float bestScore = -std::numeric_limits<float>::max();
    art::Ptr<recob::PFParticle> bestCandidate;

    for (const auto &candidate : finalCandidates)
    {
        // This shouldn't fail as we insisted a track was available above
        const auto track = CollectionHelper::GetSingleAssociated(pfParticle, pfpToTrack);
        const auto length = track->Length();

        if (length > bestScore)
        {
            bestScore = length;
            bestCandidate = candidate;
        }
    }

    if (bestCandidate.isNull())
        throw cet::exception("EnergyEstimator::CheckIfShouldUseForCalculation") << " - All particles passing cut have no length." << std::endl;

    return (pfParticle == bestCandidate);
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
    m_outputParticle.m_generation = -std::numeric_limits<int>::max();
    m_outputParticle.m_isPrimary = false;
    m_outputParticle.m_nDaughters = -std::numeric_limits<int>::max();
    m_outputParticle.m_integratedRange = -std::numeric_limits<float>::max();
    m_outputParticle.m_hasTrack = false;
    m_outputParticle.m_start = TVector3(-std::numeric_limits<float>::max(), -std::numeric_limits<float>::max(), -std::numeric_limits<float>::max());
    m_outputParticle.m_end = TVector3(-std::numeric_limits<float>::max(), -std::numeric_limits<float>::max(), -std::numeric_limits<float>::max());
    m_outputParticle.m_thetaYZ = -std::numeric_limits<float>::max(); 
    m_outputParticle.m_thetaXZ = -std::numeric_limits<float>::max(); 
    m_outputParticle.m_thetaXY = -std::numeric_limits<float>::max(); 
    m_outputParticle.m_length = -std::numeric_limits<float>::max();
    m_outputParticle.m_range = -std::numeric_limits<float>::max();
    m_outputParticle.m_endPointDist = -std::numeric_limits<float>::max();        
    m_outputParticle.m_parentCosOpeningAngle = -std::numeric_limits<float>::max();
    m_outputParticle.m_daughterCosOpeningAngle = -std::numeric_limits<float>::max();
    m_outputParticle.m_hasPid = false;
    m_outputParticle.m_muonLikelihoodU = -std::numeric_limits<float>::max(); 
    m_outputParticle.m_muonLikelihoodV = -std::numeric_limits<float>::max(); 
    m_outputParticle.m_muonLikelihoodW = -std::numeric_limits<float>::max(); 
    m_outputParticle.m_pionLikelihoodU = -std::numeric_limits<float>::max(); 
    m_outputParticle.m_pionLikelihoodV = -std::numeric_limits<float>::max(); 
    m_outputParticle.m_pionLikelihoodW = -std::numeric_limits<float>::max(); 
    m_outputParticle.m_protonLikelihoodU = -std::numeric_limits<float>::max(); 
    m_outputParticle.m_protonLikelihoodV = -std::numeric_limits<float>::max(); 
    m_outputParticle.m_protonLikelihoodW = -std::numeric_limits<float>::max(); 
    m_outputParticle.m_mipLikelihoodU = -std::numeric_limits<float>::max(); 
    m_outputParticle.m_mipLikelihoodV = -std::numeric_limits<float>::max(); 
    m_outputParticle.m_mipLikelihoodW = -std::numeric_limits<float>::max(); 
    m_outputParticle.m_shouldUseMuon = false;
    m_outputParticle.m_shouldUsePion = false;
    m_outputParticle.m_shouldUseProton = false;
}

} // namespace ubcc1pi
