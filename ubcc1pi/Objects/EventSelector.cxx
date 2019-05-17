/**
 *  @file  ubcc1pi/Objects/EventSelector.h
 *
 *  @brief The implementation file of the event selector class
 */

#include "ubcc1pi/Objects/EventSelector.h"

#include "ubcc1pi/Helpers/RecoHelper.h"

namespace ubcc1pi
{

EventSelector::EventSelector(const PFParticleVector &finalStates, const art::Event &event, const art::InputTag &pfParticleLabel, const art::InputTag &trackLabel, const art::InputTag &pidLabel) :
    m_finalStates(finalStates),
    m_pEvent(&event)
{
    m_pfParticleToMetadata = CollectionHelper::GetAssociation<recob::PFParticle, larpandoraobj::PFParticleMetadata>(*m_pEvent, pfParticleLabel);
    m_pfParticleToTrack = CollectionHelper::GetAssociation<recob::PFParticle, recob::Track>(*m_pEvent, pfParticleLabel, trackLabel);
    m_trackToPID = CollectionHelper::GetAssociation<recob::Track, anab::ParticleID>(*m_pEvent, trackLabel, pidLabel);

    // TODO make this work with FW_SEARCH_PATHS
    TFile file(this->GetPIDHistFileName().c_str());

    file.GetObject("hTrackShower_muon", m_hTrackShower_muon);
    file.GetObject("hTrackShower_pion", m_hTrackShower_pion);
    file.GetObject("hTrackShower_proton", m_hTrackShower_proton);
    file.GetObject("hProtonMIP_muon", m_hProtonMIP_muon);
    file.GetObject("hProtonMIP_pion", m_hProtonMIP_pion);
    file.GetObject("hProtonMIP_proton", m_hProtonMIP_proton);
    file.GetObject("hMuonPion_muon", m_hMuonPion_muon);
    file.GetObject("hMuonPion_pion", m_hMuonPion_pion);
    file.GetObject("hMuonPion_proton", m_hMuonPion_proton);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

std::string EventSelector::GetPIDHistFileName() const
{
    const std::string pidHistFileName = "ubcc1pi_pid_hist.root";
    std::string fullPIDHistFileName;
    cet::search_path sp("FW_SEARCH_PATH");

    if (!sp.find_file(pidHistFileName, fullPIDHistFileName))
        throw cet::exception("EventSelector::EventSelector") << " - Failed to find file: " << pidHistFileName << " in FW_SEARCH_PATH." << std::endl;

    return fullPIDHistFileName;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

float EventSelector::GetTrackShowerScore(const art::Ptr<recob::PFParticle> &pfParticle) const
{
    const auto metadata = CollectionHelper::GetSingleAssociated(pfParticle, m_pfParticleToMetadata);
    return RecoHelper::GetTrackScore(metadata);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

float EventSelector::GetProtonMIPScore(const art::Ptr<recob::PFParticle> &pfParticle) const
{
    art::Ptr<anab::ParticleID> pid;

    try
    {
        const auto track = CollectionHelper::GetSingleAssociated(pfParticle, m_pfParticleToTrack);
        pid = CollectionHelper::GetSingleAssociated(track, m_trackToPID);
    }
    catch (const cet::exception &)
    {
        return -1.f;
    }

    // Find the score from the required algorithm
    for (const auto &algorithm : pid->ParticleIDAlgScores())
    {
        if (algorithm.fAlgName == "Chi2" && algorithm.fAssumedPdg == 2212 && algorithm.fPlaneMask[0])
            return algorithm.fValue;
    }

    return -1.f;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

float EventSelector::GetMuonPionScore(const art::Ptr<recob::PFParticle> &pfParticle) const
{
    try
    {
        const auto track = CollectionHelper::GetSingleAssociated(pfParticle, m_pfParticleToTrack);
        return track->Length();
    }
    catch (const cet::exception &)
    {
    }

    return -1.f;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

float EventSelector::GetLikelihood(const float value, const TH1F *pHistogram) const
{
    const auto min = pHistogram->GetBinLowEdge(1);
    const auto max = pHistogram->GetBinLowEdge(pHistogram->GetNbinsX()) + pHistogram->GetBinWidth(pHistogram->GetNbinsX());

    if (value < min || value > max)
        return -1.f;

    return pHistogram->GetBinContent(pHistogram->FindFixBin(value));
}

// -----------------------------------------------------------------------------------------------------------------------------------------

float EventSelector::GetLikelihood(const art::Ptr<recob::PFParticle> &pfParticle, const TH1F *pTrackShower, const TH1F *pProtonMIP, const TH1F *pMuonPion) const
{
    const auto trackShowerScore = this->GetTrackShowerScore(pfParticle);
    const auto protonMIPScore = this->GetProtonMIPScore(pfParticle);
    const auto muonPionScore = this->GetMuonPionScore(pfParticle);

    if (trackShowerScore < 0.f || protonMIPScore < 0.f || muonPionScore < 0.f)
        return -1.f;

    const auto trackShowerLikelihood = this->GetLikelihood(trackShowerScore, pTrackShower);
    const auto protonMIPLikelihood = this->GetLikelihood(protonMIPScore, pProtonMIP);
    const auto muonPionLikelihood = this->GetLikelihood(muonPionScore, pMuonPion);
    
    if (trackShowerLikelihood < 0.f || protonMIPLikelihood < 0.f || muonPionLikelihood < 0.f)
        return -1.f;

    return trackShowerLikelihood * protonMIPLikelihood * muonPionLikelihood;
}

// -----------------------------------------------------------------------------------------------------------------------------------------
        
float EventSelector::GetMuonLikelihood(const art::Ptr<recob::PFParticle> &pfParticle) const
{
    return this->GetLikelihood(pfParticle, m_hTrackShower_muon, m_hProtonMIP_muon, m_hMuonPion_muon);
}

// -----------------------------------------------------------------------------------------------------------------------------------------
        
float EventSelector::GetPionLikelihood(const art::Ptr<recob::PFParticle> &pfParticle) const
{
    return this->GetLikelihood(pfParticle, m_hTrackShower_pion, m_hProtonMIP_pion, m_hMuonPion_pion);
}

// -----------------------------------------------------------------------------------------------------------------------------------------
        
float EventSelector::GetProtonLikelihood(const art::Ptr<recob::PFParticle> &pfParticle) const
{
    return this->GetLikelihood(pfParticle, m_hTrackShower_proton, m_hProtonMIP_proton, m_hMuonPion_proton);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

float EventSelector::GetBestCC1PiScore(art::Ptr<recob::PFParticle> &bestMuon, art::Ptr<recob::PFParticle> &bestPion, PFParticleVector &bestProtons) const
{
    if (m_finalStates.size() < 2)
        throw cet::exception("EventSelector::GetBestCC1PiScore") << " - A CC1Pi+ event must have at least 2 particles, " << m_finalStates.size() << " were provided" << std::endl;

    if (!bestProtons.empty())
        throw cet::exception("EventSelector::GetBestCC1PiScore") << " - The input vetor of best protons wasn't empty" << std::endl;

    // Filter out particles that have an incalculable likelihood under any hypothesisi
    PFParticleVector usableParticles;
    for (const auto &particle : m_finalStates)
    {
        if (this->GetMuonLikelihood(particle) < 0.f || this->GetPionLikelihood(particle) < 0.f || this->GetProtonLikelihood(particle) < 0.f)
            continue;

        usableParticles.push_back(particle);
    }

    if (usableParticles.size() < 2)
        return -1.f;

    // Find the best PID for each particle assuming the event is CC1Pi
    float bestTotalLikelihood = -std::numeric_limits<float>::max();
    for (const auto &muonCandidate : usableParticles)
    {
        for (const auto &pionCandidate : usableParticles)
        {
            if (pionCandidate == muonCandidate)
                continue;

            const auto muonLikelihood = this->GetMuonLikelihood(muonCandidate);
            const auto pionLikelihood = this->GetPionLikelihood(pionCandidate);

            float totalLikelihood = muonLikelihood * pionLikelihood;
            for (const auto &protonCandidate : usableParticles)
            {
                if (protonCandidate == muonCandidate || protonCandidate == pionCandidate)
                    continue;

                const auto protonLikelihood = this->GetProtonLikelihood(protonCandidate);
                totalLikelihood *= protonLikelihood;
            }

            if (totalLikelihood < bestTotalLikelihood)
                continue;

            bestTotalLikelihood = totalLikelihood;
            bestMuon = muonCandidate;
            bestPion = pionCandidate;
        }
    }

    // Collect the protons
    for (const auto &particle : usableParticles)
    {
        if (particle != bestMuon && particle != bestPion)
            bestProtons.push_back(particle);
    }

    return bestTotalLikelihood;
}

} // namespace ubcc1pi
