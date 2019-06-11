/**
 *  @file  ubcc1pi/Analyzers/PIDStudy_module.cc
 *
 *  @brief The implementation file for the PID study analyzer.
 */

#include "ubcc1pi/Analyzers/PIDStudy.h"

#include "ubcc1pi/Helpers/TruthHelper.h"
#include "ubcc1pi/Helpers/AnalysisHelper.h"
#include "ubcc1pi/Helpers/RecoHelper.h"
#include "ubcc1pi/Helpers/DebugHelper.h"


namespace ubcc1pi
{

PIDStudy::PIDStudy(const art::EDAnalyzer::Table<Config> &config) :
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
    
    m_pParticleTree->Branch("nuVertex", &m_outputParticle.m_nuVertex);
    m_pParticleTree->Branch("nuVertexCorrected", &m_outputParticle.m_nuVertexCorrected);

    m_pParticleTree->Branch("hasMatchedMCParticle", &m_outputParticle.m_hasMatchedMCParticle);
    m_pParticleTree->Branch("truePdgCode", &m_outputParticle.m_truePdgCode);
    m_pParticleTree->Branch("trueMomentum", &m_outputParticle.m_trueMomentum);
    m_pParticleTree->Branch("trueMatchPurity", &m_outputParticle.m_trueMatchPurity);
    m_pParticleTree->Branch("trueMatchCompleteness", &m_outputParticle.m_trueMatchCompleteness);

    m_pParticleTree->Branch("nHitsU", &m_outputParticle.m_nHitsU);
    m_pParticleTree->Branch("nHitsV", &m_outputParticle.m_nHitsV);
    m_pParticleTree->Branch("nHitsW", &m_outputParticle.m_nHitsW);

    m_pParticleTree->Branch("start", &m_outputParticle.m_start);
    m_pParticleTree->Branch("end", &m_outputParticle.m_end);
    m_pParticleTree->Branch("startCorrected", &m_outputParticle.m_startCorrected);
    m_pParticleTree->Branch("endCorrected", &m_outputParticle.m_endCorrected);
    
    m_pParticleTree->Branch("length", &m_outputParticle.m_length);

    m_pParticleTree->Branch("trackShower", &m_outputParticle.m_trackShower);
    m_pParticleTree->Branch("mipFraction", &m_outputParticle.m_mipFraction);
    m_pParticleTree->Branch("primaryFraction", &m_outputParticle.m_primaryFraction);

    m_pParticleTree->Branch("chi2_mu_U", &m_outputParticle.m_chi2_mu_U);
    m_pParticleTree->Branch("chi2_pi_U", &m_outputParticle.m_chi2_pi_U);
    m_pParticleTree->Branch("chi2_p_U", &m_outputParticle.m_chi2_p_U);
    m_pParticleTree->Branch("chi2_mu_V", &m_outputParticle.m_chi2_mu_V);
    m_pParticleTree->Branch("chi2_pi_V", &m_outputParticle.m_chi2_pi_V);
    m_pParticleTree->Branch("chi2_p_V", &m_outputParticle.m_chi2_p_V);
    m_pParticleTree->Branch("chi2_mu_W", &m_outputParticle.m_chi2_mu_W);
    m_pParticleTree->Branch("chi2_pi_W", &m_outputParticle.m_chi2_pi_W);
    m_pParticleTree->Branch("chi2_p_W", &m_outputParticle.m_chi2_p_W);

    m_pParticleTree->Branch("braggPeakLLH_mu_U", &m_outputParticle.m_braggPeakLLH_mu_U);
    m_pParticleTree->Branch("braggPeakLLH_pi_U", &m_outputParticle.m_braggPeakLLH_pi_U);
    m_pParticleTree->Branch("braggPeakLLH_p_U", &m_outputParticle.m_braggPeakLLH_p_U);
    m_pParticleTree->Branch("braggPeakLLH_MIP_U", &m_outputParticle.m_braggPeakLLH_MIP_U);
    m_pParticleTree->Branch("braggPeakLLH_mu_V", &m_outputParticle.m_braggPeakLLH_mu_V);
    m_pParticleTree->Branch("braggPeakLLH_pi_V", &m_outputParticle.m_braggPeakLLH_pi_V);
    m_pParticleTree->Branch("braggPeakLLH_p_V", &m_outputParticle.m_braggPeakLLH_p_V);
    m_pParticleTree->Branch("braggPeakLLH_MIP_V", &m_outputParticle.m_braggPeakLLH_MIP_V);
    m_pParticleTree->Branch("braggPeakLLH_mu_W", &m_outputParticle.m_braggPeakLLH_mu_W);
    m_pParticleTree->Branch("braggPeakLLH_pi_W", &m_outputParticle.m_braggPeakLLH_pi_W);
    m_pParticleTree->Branch("braggPeakLLH_p_W", &m_outputParticle.m_braggPeakLLH_p_W);
    m_pParticleTree->Branch("braggPeakLLH_MIP_W", &m_outputParticle.m_braggPeakLLH_MIP_W);
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
    const auto calorimetryLabel = m_config().CalorimetryLabel();
    
    // Get the neutrino final state PFParticles
    const auto allPFParticles = CollectionHelper::GetCollection<recob::PFParticle>(event, pfParticleLabel);
    const auto pfParticleMap = RecoHelper::GetPFParticleMap(allPFParticles);
    const auto pfParticles = RecoHelper::GetNeutrinoFinalStates(allPFParticles);

    // Get the associations from PFParticles to metadata, tracks and PID
    const auto pfpToTrack = CollectionHelper::GetAssociation<recob::PFParticle, recob::Track>(event, pfParticleLabel, trackLabel);
    const auto pfpToMetadata = CollectionHelper::GetAssociation<recob::PFParticle, larpandoraobj::PFParticleMetadata>(event, pfParticleLabel);
    const auto pfpToHits = CollectionHelper::GetAssociationViaCollection<recob::PFParticle, recob::Cluster, recob::Hit>(event, pfParticleLabel, pfParticleLabel, pfParticleLabel);
    const auto trackToPID = CollectionHelper::GetAssociation<recob::Track, anab::ParticleID>(event, trackLabel, pidLabel);
    const auto trackToCalorimetry = CollectionHelper::GetAssociation<recob::Track, anab::Calorimetry>(event, trackLabel, calorimetryLabel);
    
    // Get the reco-true matching information
    const auto backtrackerData = AnalysisHelper::GetBacktrackerData(event, mcTruthLabel, mcParticleLabel, backtrackerLabel, pfParticleLabel);

    // Get the space-charge service
    const auto pSpaceChargeService = RecoHelper::GetSpaceChargeService();

    // Store the event level information
    this->SetEventInfo(event, allPFParticles);

    for (const auto &pfParticle : pfParticles)
    {
        this->ResetParticleInfo();
        this->SetMatchedMCParticleInfo(pfParticle, backtrackerData);
       
        // Number of hits
        const auto &hits = backtrackerData.GetHits(pfParticle);
        m_outputParticle.m_nHitsU = RecoHelper::CountHitsInView(hits, geo::kU);
        m_outputParticle.m_nHitsV = RecoHelper::CountHitsInView(hits, geo::kV);
        m_outputParticle.m_nHitsW = RecoHelper::CountHitsInView(hits, geo::kW);

        // Track-shower score 
        const auto metadata = CollectionHelper::GetSingleAssociated(pfParticle, pfpToMetadata);
        m_outputParticle.m_trackShower = RecoHelper::GetTrackScore(metadata);
        
        art::Ptr<recob::Track> track;
        if (!this->GetTrack(pfParticle, pfpToTrack, track))
            continue;
   
        // Get the track start and end points, and correct them for space charge
        m_outputParticle.m_start = TVector3(track->Start().X(), track->Start().Y(), track->Start().Z());
        m_outputParticle.m_end = TVector3(track->End().X(), track->End().Y(), track->End().Z());
        m_outputParticle.m_startCorrected = RecoHelper::CorrectForSpaceCharge(m_outputParticle.m_start, pSpaceChargeService);
        m_outputParticle.m_endCorrected = RecoHelper::CorrectForSpaceCharge(m_outputParticle.m_end, pSpaceChargeService);

        // Track length
        m_outputParticle.m_length = track->Length();

        /* BEGIN TEST */
        m_outputParticle.m_mipFraction = -std::numeric_limits<float>::max();
        const auto calos = CollectionHelper::GetManyAssociated(track, trackToCalorimetry);
        float nPoints = 0.f;
        float nMIPPoints = 0.f;
        for (const auto &calo : calos)
        {
            for (const auto &dEdx : calo->dEdx())
            {
                nPoints += 1.f;
                if (dEdx < 2.5)
                {
                    nMIPPoints += 1.f;
                }
            }
        }

        if (nPoints >= std::numeric_limits<float>::epsilon())
        {
            m_outputParticle.m_mipFraction = nMIPPoints / nPoints;
        }
        /* END TEST */

        /* BEGIN TEST */
        m_outputParticle.m_primaryFraction = -std::numeric_limits<float>::max();

        unsigned int nHitsTotal = 0;
        for (const auto &downstreamPFParticle : RecoHelper::GetDownstreamParticles(pfParticle, pfParticleMap))
        {
            nHitsTotal += CollectionHelper::GetManyAssociated(downstreamPFParticle, pfpToHits).size();
        }
        const auto nPrimaryHits = CollectionHelper::GetManyAssociated(pfParticle, pfpToHits).size();

        if (nHitsTotal > 0)
        {
            m_outputParticle.m_primaryFraction = static_cast<float>(nPrimaryHits) / static_cast<float>(nHitsTotal);
        }
        /* END TEST */

        // PID algorithm outputs
        const auto pid = CollectionHelper::GetSingleAssociated(track, trackToPID);
        for (const auto &algo : pid->ParticleIDAlgScores())
        {
            const auto view = this->GetView(algo.fPlaneMask);

            // Chi2 algorithms
            if (algo.fAlgName == "Chi2")
            {
                switch (algo.fAssumedPdg)
                {
                    case 13:
                        switch (view)
                        {
                            case geo::kW:
                                m_outputParticle.m_chi2_mu_W = algo.fValue;
                                break;
                            case geo::kU:
                                m_outputParticle.m_chi2_mu_U = algo.fValue;
                                break;
                            case geo::kV:
                                m_outputParticle.m_chi2_mu_V = algo.fValue;
                                break;
                            default: break;
                        }
                        break;
                    case 211:
                        switch (view)
                        {
                            case geo::kW:
                                m_outputParticle.m_chi2_pi_W = algo.fValue;
                                break;
                            case geo::kU:
                                m_outputParticle.m_chi2_pi_U = algo.fValue;
                                break;
                            case geo::kV:
                                m_outputParticle.m_chi2_pi_V = algo.fValue;
                                break;
                            default: break;
                        }
                        break;
                    case 2212:
                        switch (view)
                        {
                            case geo::kW:
                                m_outputParticle.m_chi2_p_W = algo.fValue;
                                break;
                            case geo::kU:
                                m_outputParticle.m_chi2_p_U = algo.fValue;
                                break;
                            case geo::kV:
                                m_outputParticle.m_chi2_p_V = algo.fValue;
                                break;
                            default: break;
                        }
                        break;
                    default: break;
                }
            }
            
            // Bragg peak algorithms
            if (algo.fAlgName == "BraggPeakLLH" && algo.fTrackDir == anab::kForward)
            {
                switch (algo.fAssumedPdg)
                {
                    case 13:
                        switch (view)
                        {
                            case geo::kW:
                                m_outputParticle.m_braggPeakLLH_mu_W = algo.fValue;
                                break;
                            case geo::kU:
                                m_outputParticle.m_braggPeakLLH_mu_U = algo.fValue;
                                break;
                            case geo::kV:
                                m_outputParticle.m_braggPeakLLH_mu_V = algo.fValue;
                                break;
                            default: break;
                        }
                        break;
                    case 211:
                        switch (view)
                        {
                            case geo::kW:
                                m_outputParticle.m_braggPeakLLH_pi_W = algo.fValue;
                                break;
                            case geo::kU:
                                m_outputParticle.m_braggPeakLLH_pi_U = algo.fValue;
                                break;
                            case geo::kV:
                                m_outputParticle.m_braggPeakLLH_pi_V = algo.fValue;
                                break;
                            default: break;
                        }
                        break;
                    case 2212:
                        switch (view)
                        {
                            case geo::kW:
                                m_outputParticle.m_braggPeakLLH_p_W = algo.fValue;
                                break;
                            case geo::kU:
                                m_outputParticle.m_braggPeakLLH_p_U = algo.fValue;
                                break;
                            case geo::kV:
                                m_outputParticle.m_braggPeakLLH_p_V = algo.fValue;
                                break;
                            default: break;
                        }
                        break;
                    case 0:
                        switch (view)
                        {
                            case geo::kW:
                                m_outputParticle.m_braggPeakLLH_MIP_W = algo.fValue;
                                break;
                            case geo::kU:
                                m_outputParticle.m_braggPeakLLH_MIP_U = algo.fValue;
                                break;
                            case geo::kV:
                                m_outputParticle.m_braggPeakLLH_MIP_V = algo.fValue;
                                break;
                            default: break;
                        }
                        break;
                    default: break;
                }
            }
        }
        
        m_pParticleTree->Fill();
    }
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void PIDStudy::SetEventInfo(const art::Event &event, const PFParticleVector &allPFParticles)
{    
    m_outputParticle.m_run = event.run();
    m_outputParticle.m_subRun = event.subRun();
    m_outputParticle.m_event = event.event();
    
    const TruthHelper::Interaction interaction(event, m_config().MCTruthLabel(), m_config().MCParticleLabel());
    m_outputParticle.m_isSignal = AnalysisHelper::IsCC1PiSignal(interaction);

    m_outputParticle.m_nuVertex = RecoHelper::GetRecoNeutrinoVertex(event, allPFParticles, m_config().PFParticleLabel());
    m_outputParticle.m_nuVertexCorrected = RecoHelper::CorrectForSpaceCharge(m_outputParticle.m_nuVertex, RecoHelper::GetSpaceChargeService());
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void PIDStudy::ResetParticleInfo()
{
    m_outputParticle.m_hasMatchedMCParticle = false;
    m_outputParticle.m_truePdgCode = -std::numeric_limits<int>::max();
    m_outputParticle.m_trueMomentum = -std::numeric_limits<float>::max();
    m_outputParticle.m_trueMatchPurity = -std::numeric_limits<float>::max();
    m_outputParticle.m_trueMatchCompleteness = -std::numeric_limits<float>::max();
    m_outputParticle.m_nHitsU = -std::numeric_limits<int>::max();
    m_outputParticle.m_nHitsV = -std::numeric_limits<int>::max();
    m_outputParticle.m_nHitsW = -std::numeric_limits<int>::max();
    m_outputParticle.m_length = -std::numeric_limits<float>::max();
    m_outputParticle.m_trackShower = -std::numeric_limits<float>::max();
    m_outputParticle.m_chi2_mu_U = -std::numeric_limits<float>::max();
    m_outputParticle.m_chi2_pi_U = -std::numeric_limits<float>::max();
    m_outputParticle.m_chi2_p_U = -std::numeric_limits<float>::max();
    m_outputParticle.m_chi2_mu_V = -std::numeric_limits<float>::max();
    m_outputParticle.m_chi2_pi_V = -std::numeric_limits<float>::max();
    m_outputParticle.m_chi2_p_V = -std::numeric_limits<float>::max();
    m_outputParticle.m_chi2_mu_W = -std::numeric_limits<float>::max();
    m_outputParticle.m_chi2_pi_W = -std::numeric_limits<float>::max();
    m_outputParticle.m_chi2_p_W = -std::numeric_limits<float>::max();
    m_outputParticle.m_braggPeakLLH_mu_U = -std::numeric_limits<float>::max();
    m_outputParticle.m_braggPeakLLH_pi_U = -std::numeric_limits<float>::max();
    m_outputParticle.m_braggPeakLLH_p_U = -std::numeric_limits<float>::max();
    m_outputParticle.m_braggPeakLLH_MIP_U = -std::numeric_limits<float>::max();
    m_outputParticle.m_braggPeakLLH_mu_V = -std::numeric_limits<float>::max();
    m_outputParticle.m_braggPeakLLH_pi_V = -std::numeric_limits<float>::max();
    m_outputParticle.m_braggPeakLLH_p_V = -std::numeric_limits<float>::max();
    m_outputParticle.m_braggPeakLLH_MIP_V = -std::numeric_limits<float>::max();
    m_outputParticle.m_braggPeakLLH_mu_W = -std::numeric_limits<float>::max();
    m_outputParticle.m_braggPeakLLH_pi_W = -std::numeric_limits<float>::max();
    m_outputParticle.m_braggPeakLLH_p_W = -std::numeric_limits<float>::max();
    m_outputParticle.m_braggPeakLLH_MIP_W = -std::numeric_limits<float>::max();
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void PIDStudy::SetMatchedMCParticleInfo(const art::Ptr<recob::PFParticle> &pfParticle, const BacktrackHelper::BacktrackerData &backtrackerData)
{
    try
    {
        const auto mcParticle = backtrackerData.GetBestMatchedMCParticle(pfParticle);

        m_outputParticle.m_truePdgCode = mcParticle->PdgCode();
        m_outputParticle.m_trueMomentum = mcParticle->P();
        m_outputParticle.m_trueMatchPurity = backtrackerData.GetMatchPurity(pfParticle, mcParticle);
        m_outputParticle.m_trueMatchCompleteness = backtrackerData.GetMatchCompleteness(pfParticle, mcParticle);
        m_outputParticle.m_hasMatchedMCParticle = true;
    }
    catch (const cet::exception &)
    {
        m_outputParticle.m_truePdgCode = -std::numeric_limits<int>::max();
        m_outputParticle.m_trueMomentum = -std::numeric_limits<float>::max();
        m_outputParticle.m_trueMatchPurity = -std::numeric_limits<float>::max();
        m_outputParticle.m_trueMatchCompleteness = -std::numeric_limits<float>::max();
        m_outputParticle.m_hasMatchedMCParticle = false;
    }
}

// -----------------------------------------------------------------------------------------------------------------------------------------
    
bool PIDStudy::GetTrack(const art::Ptr<recob::PFParticle> &pfParticle, const Association<recob::PFParticle, recob::Track> &pfpToTrack, art::Ptr<recob::Track> &track)
{
    try
    {
        const auto tracks = CollectionHelper::GetManyAssociated(pfParticle, pfpToTrack); 
        
        if (tracks.empty())
            return false;

        if (tracks.size() > 1)
            throw cet::exception("PIDStudy::GetTrack") << " - Multiple tracks associated to PFParticle!" << std::endl;

        track = tracks.front();
        return true;
    }
    catch (const cet::exception &)
    {
        return false;
    }
}       

// -----------------------------------------------------------------------------------------------------------------------------------------

geo::View_t PIDStudy::GetView(const std::bitset<8> &planeMask) const
{
    // Here is a hack to get around a bug in the PID code. Some algorithms call W = 0, U = 1, V = 2. But others call W = 7, U = 6, V = 5
    const bool usesW = planeMask.test(0) || planeMask.test(7);
    const bool usesU = planeMask.test(1) || planeMask.test(6);
    const bool usesV = planeMask.test(2) || planeMask.test(5);
    
    if (usesW && !usesU && !usesV)
        return geo::kW;
    
    if (!usesW && usesU && !usesV)
        return geo::kU;
    
    if (!usesW && !usesU && usesV)
        return geo::kV;

    return geo::kUnknown;

}

} // namespace ubcc1pi
