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

    std::cout << "Number of tracks = " << CollectionHelper::GetCollection<recob::Track>(*m_pEvent, trackLabel).size() << std::endl;
    std::cout << "Track to PID associations = " << m_trackToPID.size() << std::endl;
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
        std::cout << "Looking for track" << std::endl;
        const auto track = CollectionHelper::GetSingleAssociated(pfParticle, m_pfParticleToTrack);
        
        std::cout << "... looking for PID" << std::endl;
        pid = CollectionHelper::GetSingleAssociated(track, m_trackToPID);
    }
    catch (const cet::exception &)
    {
        std::cout << "... Failed!" << std::endl;
        return -1.f;
    }

    // Find the score from the required algorithm
    std::cout << "... looking for Chi2" << std::endl;
    for (const auto &algorithm : pid->ParticleIDAlgScores())
    {
        if (algorithm.fAlgName == "Chi2" && algorithm.fAssumedPdg == 2212 && algorithm.fPlaneMask[0])
            return algorithm.fValue;
    }
   
    std::cout << "... Failed!" << std::endl;

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
        
float EventSelector::GetMuonLikelihood(const art::Ptr<recob::PFParticle> &pfParticle) const
{
    /*
    const auto trackShowerScore = this->GetTrackShowerScore(pfParticle);
    const auto protonMIPScore = this->GetProtonMIPScore(pfParticle);
    const auto muonPionScore = this->GetMuonPionScore(pfParticle);
    */

    return 0.5f;
}

// -----------------------------------------------------------------------------------------------------------------------------------------
        
float EventSelector::GetPionLikelihood(const art::Ptr<recob::PFParticle> &pfParticle) const
{
    /*
    const auto trackShowerScore = this->GetTrackShowerScore(pfParticle);
    const auto protonMIPScore = this->GetProtonMIPScore(pfParticle);
    const auto muonPionScore = this->GetMuonPionScore(pfParticle);
    */

    return 0.5f;
}

// -----------------------------------------------------------------------------------------------------------------------------------------
        
float EventSelector::GetProtonLikelihood(const art::Ptr<recob::PFParticle> &pfParticle) const
{
    /*
    const auto trackShowerScore = this->GetTrackShowerScore(pfParticle);
    const auto protonMIPScore = this->GetProtonMIPScore(pfParticle);
    const auto muonPionScore = this->GetMuonPionScore(pfParticle);
    */

    return 0.5f;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

float EventSelector::GetBestCC1PiScore(art::Ptr<recob::PFParticle> &bestMuon, art::Ptr<recob::PFParticle> &bestPion) const
{
    if (m_finalStates.size() < 2)
        throw cet::exception("EventSelector::GetBestCC1PiScore") << " - A CC1Pi+ event must have at least 2 particles, " << m_finalStates.size() << " were provided" << std::endl;

    float bestTotalLikelihood = -std::numeric_limits<float>::max();

    for (const auto &muonCandidate : m_finalStates)
    {
        for (const auto &pionCandidate : m_finalStates)
        {
            if (pionCandidate == muonCandidate)
                continue;

            const auto muonLikelihood = this->GetMuonLikelihood(muonCandidate);
            const auto pionLikelihood = this->GetPionLikelihood(pionCandidate);

            std::cout << muonCandidate->Self() << " = Muon : " << muonLikelihood << std::endl; 
            std::cout << pionCandidate->Self() << " = Pion : " << pionLikelihood << std::endl; 

            float totalLikelihood = muonLikelihood * pionLikelihood;
            for (const auto &protonCandidate : m_finalStates)
            {
                if (protonCandidate == muonCandidate || protonCandidate == pionCandidate)
                    continue;

                const auto protonLikelihood = this->GetProtonLikelihood(protonCandidate);
                std::cout << protonCandidate->Self() << " = Proton : " << protonLikelihood << std::endl;

                totalLikelihood *= protonLikelihood;
            }

            std::cout << "Total score = " << totalLikelihood << std::endl;
            std::cout << "---------------------------------------------------------" << std::endl;

            if (totalLikelihood < bestTotalLikelihood)
                continue;

            bestTotalLikelihood = totalLikelihood;
            bestMuon = muonCandidate;
            bestPion = pionCandidate;
        }
    }

    return bestTotalLikelihood;
}

} // namespace ubcc1pi
