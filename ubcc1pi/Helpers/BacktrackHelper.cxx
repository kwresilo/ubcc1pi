/**
 *  @file  ubcc1pi/Helpers/BacktrackHelper.cxx
 *
 *  @brief The implementation file for the back tracking helper class
 */

#include "ubcc1pi/Helpers/BacktrackHelper.h"

#include "ubcc1pi/Helpers/CollectionHelper.h"
#include "ubcc1pi/Helpers/RecoHelper.h"
#include "ubcc1pi/Helpers/AnalysisHelper.h"

namespace ubcc1pi
{

HitsToPFParticles BacktrackHelper::GetHitToPFParticleMap(const art::Event &event, const art::InputTag &pfParticleLabel, const PFParticleVector &finalStates)
{
    HitsToPFParticles outputMap;

    const auto pfParticles = CollectionHelper::GetCollection<recob::PFParticle>(event, pfParticleLabel);
    const auto pfpToCluster = CollectionHelper::GetAssociation<recob::PFParticle, recob::Cluster>(event, pfParticleLabel);
    const auto clusterToHits = CollectionHelper::GetAssociation<recob::Cluster, recob::Hit>(event, pfParticleLabel);

    for (const auto &finalState : finalStates)
    {
        for (const auto &particle : RecoHelper::GetDownstreamParticles(finalState, RecoHelper::GetPFParticleMap(pfParticles)))
        {
            for (const auto &cluster : CollectionHelper::GetManyAssociated(particle, pfpToCluster))
            {
                for (const auto &hit : CollectionHelper::GetManyAssociated(cluster, clusterToHits))
                {
                    outputMap[hit].push_back(finalState);
                }
            }
        }
    }

    return outputMap;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

HitsToMCParticleWeights BacktrackHelper::GetHitToMCParticleWeightMap(const art::Event &event, const art::InputTag &mcParticleLabel, const art::InputTag &backtrackerLabel, const MCParticleVector &finalStates)
{
    HitsToMCParticleWeights outputMap;

    const auto mcParticles = CollectionHelper::GetCollection<simb::MCParticle>(event, mcParticleLabel);
    const auto mcpToHitMap = CollectionHelper::GetAssociationWithData<simb::MCParticle, recob::Hit, anab::BackTrackerHitMatchingData>(event, mcParticleLabel, backtrackerLabel);

    for (const auto &finalState : finalStates)
    {
        for (const auto &particle : TruthHelper::GetDownstreamParticles(finalState, TruthHelper::GetMCParticleMap(mcParticles)))
        {
            CollectionData<recob::Hit, anab::BackTrackerHitMatchingData> hits;
            try
            {
                hits = CollectionHelper::GetManyAssociatedWithData(particle, mcpToHitMap);
            }
            catch (const cet::exception &)
            {
                // An MCParticle can have no hits
                continue;
            }

            for (const auto &entry : hits)
            {
                const auto weight = entry.second.ideFraction;

                auto &mcParticleWeightVector = outputMap[entry.first];
                auto iter = std::find_if(mcParticleWeightVector.begin(), mcParticleWeightVector.end(), [&](const auto &x){return x.first == finalState;});

                if (iter == mcParticleWeightVector.end())
                {
                    mcParticleWeightVector.emplace_back(finalState, weight);
                    continue;
                }

                iter->second += weight;
            }
        }
    }

    return outputMap;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

HitsToBool BacktrackHelper::GetHitsToIsNuInducedMap(const art::Event &event, const art::InputTag &hitLabel, const art::InputTag &mcParticleLabel, const art::InputTag &backtrackerLabel, const MCParticleVector &nuMCParticles)
{
    const auto hits = CollectionHelper::GetCollection<recob::Hit>(event, hitLabel);
    const auto mcpToHitMap = CollectionHelper::GetAssociation<simb::MCParticle, recob::Hit>(event, mcParticleLabel, backtrackerLabel);
    const auto hitToMcpMap = CollectionHelper::GetReversedAssociation(mcpToHitMap);
    
    HitsToBool outputMap;

    for (const auto &hit : hits)
    {
        // Start by assuming that the hit isn't neutrino induced, and change this later if required
        if (!outputMap.emplace(hit, false).second)
            throw cet::exception("BacktrackHelper::GetHitsToIsNuInducedMap") << " - Repeated input hits." << std::endl;

        // Check the hit is associated to an MCParticle
        const auto iter = hitToMcpMap.find(hit);
        if (iter == hitToMcpMap.end())
            continue;

        // Check the hit is associated to a neutrino induced MCParticle
        bool isNuInduced = false;
        for (const auto &mcParticle : iter->second)
        {
            if (std::find(nuMCParticles.begin(), nuMCParticles.end(), mcParticle) != nuMCParticles.end())
            {
                isNuInduced = true;
                break;
            }
        }

        if (isNuInduced)
            outputMap.at(hit) = true;
    }

    return outputMap;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

BacktrackHelper::BacktrackerData::BacktrackerData(const PFParticleVector &pfParticles, const MCParticleVector &mcParticles, const HitsToPFParticles &hitsToPfps, const HitsToMCParticleWeights &hitsToMcps) :
    m_mcParticles(mcParticles),
    m_pfParticles(pfParticles),
    m_hitsToMcps(hitsToMcps)
{
    this->CheckConsistency(m_mcParticles, hitsToMcps);
    this->CheckConsistency(m_pfParticles, hitsToPfps);

    const auto allHits = this->CollectHits(hitsToPfps, hitsToMcps);

    for (const auto &hit : this->CollectHits(hitsToPfps, hitsToMcps))
    {
        art::Ptr<recob::PFParticle> pfParticle;
        const bool hasPFParticle = this->CollectPFParticle(hit, hitsToPfps, pfParticle);
        
        CollectionData<simb::MCParticle, float> mcParticleWeights;
        const bool hasMCParticles = this->CollectMCParticleWeights(hit, hitsToMcps, mcParticleWeights);

        // Fill the PFParticle weight map
        if (hasPFParticle)
        {
            // Check if we already have an entry for this particle
            if (m_pfParticleWeightMap.find(pfParticle) == m_pfParticleWeightMap.end())
                m_pfParticleWeightMap.emplace(pfParticle, 0.f);

            m_pfParticleWeightMap.at(pfParticle) += 1.f;
        }
        
        // Fill the MCParticle weight map
        if (hasMCParticles)
        {
            for (const auto &entry : mcParticleWeights)
            {
                const auto mcParticle = entry.first;
                const auto weight = entry.second;
            
                // Check if we already have an entry for this particle
                if (m_mcParticleWeightMap.find(mcParticle) == m_mcParticleWeightMap.end())
                    m_mcParticleWeightMap.emplace(mcParticle, 0.f);
    
                m_mcParticleWeightMap.at(mcParticle) += weight;
                
                // Fill the PFParticle -> MCParticle matched weight map
                if (hasPFParticle)
                {
                    // Get the map entry for this PFParticle, creating it if it doesn't exist
                    auto &entry = m_matchMap[pfParticle];

                    auto iter = std::find_if(entry.begin(), entry.end(), [&](const auto &x){return x.first == mcParticle;});
                    if (iter == entry.end())
                    {
                        entry.emplace_back(mcParticle, weight);
                    }
                    else
                    {
                        iter->second += weight;
                    }
                }
            }
        }
    }
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void BacktrackHelper::BacktrackerData::CheckConsistency(const MCParticleVector &mcParticles, const HitsToMCParticleWeights &hitsToMcps) const
{
    for (const auto &entry : hitsToMcps)
    {
        for (const auto &mcpWeightPair : entry.second)
        {
            if (std::find(mcParticles.begin(), mcParticles.end(), mcpWeightPair.first) == mcParticles.end())
                throw cet::exception("BacktrackerData::CheckConsistency") << " - MCParticle from input map isn't in the input vector." << std::endl;
        }
    }
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void BacktrackHelper::BacktrackerData::CheckConsistency(const PFParticleVector &pfParticles, const HitsToPFParticles &hitsToPfps) const
{
    for (const auto &entry : hitsToPfps)
    {
        for (const auto &pfParticle : entry.second)
        {
            if (std::find(pfParticles.begin(), pfParticles.end(), pfParticle) == pfParticles.end())
                throw cet::exception("BacktrackerData::CheckConsistency") << " - PFParticle from input map isn't in the input vector." << std::endl;
        }
    }
}

// -----------------------------------------------------------------------------------------------------------------------------------------

HitVector BacktrackHelper::BacktrackerData::CollectHits(const HitsToPFParticles &hitsToPfps, const HitsToMCParticleWeights &hitsToMcps) const
{
    HitVector hits;
    
    for (const auto &entry : hitsToPfps)
    {
        if (std::find(hits.begin(), hits.end(), entry.first) == hits.end())
            hits.push_back(entry.first);
    }
    
    for (const auto &entry : hitsToMcps)
    {
        if (std::find(hits.begin(), hits.end(), entry.first) == hits.end())
            hits.push_back(entry.first);
    }

    // Sort the hits to ensure reproducibility
    std::sort(hits.begin(), hits.end(), [](const auto &hit1, const auto &hit2){return hit1->PeakTime() < hit2->PeakTime();});

    return hits;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

bool BacktrackHelper::BacktrackerData::CollectPFParticle(const art::Ptr<recob::Hit> &hit, const HitsToPFParticles &hitsToPfps, art::Ptr<recob::PFParticle> &outputParticle) const
{
    try
    {
        outputParticle = CollectionHelper::GetSingleAssociated(hit, hitsToPfps);
        return true;
    }
    catch (const cet::exception &)
    {
        return false;
    }
}

// -----------------------------------------------------------------------------------------------------------------------------------------

bool BacktrackHelper::BacktrackerData::CollectMCParticleWeights(const art::Ptr<recob::Hit> &hit, const HitsToMCParticleWeights &hitsToMcps, CollectionData<simb::MCParticle, float> &outputParticleWeights) const
{
    try
    {
        outputParticleWeights = CollectionHelper::GetManyAssociatedWithData(hit, hitsToMcps);
        return true;
    }
    catch (const cet::exception &)
    {
        return false;
    }
}

// -----------------------------------------------------------------------------------------------------------------------------------------
                
MCParticleVector BacktrackHelper::BacktrackerData::GetMCParticles() const
{
    return m_mcParticles;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

PFParticleVector BacktrackHelper::BacktrackerData::GetPFParticles() const
{    
    return m_pfParticles;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

unsigned int BacktrackHelper::BacktrackerData::GetNHits(const art::Ptr<recob::PFParticle> &pfParticle) const
{
    return static_cast<unsigned int>(std::round(this->GetWeight(pfParticle)));
}

// -----------------------------------------------------------------------------------------------------------------------------------------

float BacktrackHelper::BacktrackerData::GetWeight(const art::Ptr<recob::PFParticle> &pfParticle) const
{
    const auto iter = m_pfParticleWeightMap.find(pfParticle);

    if (iter == m_pfParticleWeightMap.end())
        return 0.f;

    return iter->second;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

float BacktrackHelper::BacktrackerData::GetWeight(const art::Ptr<simb::MCParticle> &mcParticle) const
{
    const auto iter = m_mcParticleWeightMap.find(mcParticle);

    if (iter == m_mcParticleWeightMap.end())
        return 0.f;

    return iter->second;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

HitVector BacktrackHelper::BacktrackerData::GetHits(const art::Ptr<simb::MCParticle> &mcParticle) const
{
    HitVector hits;

    for (const auto &entry : m_hitsToMcps)
    {
        for (const auto &pair : entry.second)
        {
            if (pair.first == mcParticle)
            {
                if (std::find(hits.begin(), hits.end(), entry.first) != hits.end())
                    throw cet::exception("BacktrackerData::GetHits") << " - Found repeated hits associated to the input MCParticle" << std::endl;

                hits.push_back(entry.first);
            }
        }
    }

    return hits;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

float BacktrackHelper::BacktrackerData::GetMatchWeight(const art::Ptr<recob::PFParticle> &pfParticle, const art::Ptr<simb::MCParticle> &mcParticle) const
{
    const auto iter = m_matchMap.find(pfParticle);

    // There were no entries for this PFParticle - MCParticle pair, so the match weight is zero
    if (iter == m_matchMap.end())
        return 0.f;

    for (const auto &entry : iter->second)
    {
        if (entry.first == mcParticle)
            return entry.second;
    }

    // There were no entries for this PFParticle - MCParticle pair, so the match weight is zero
    return 0.f;
}

// -----------------------------------------------------------------------------------------------------------------------------------------
                
float BacktrackHelper::BacktrackerData::GetMatchPurity(const art::Ptr<recob::PFParticle> &pfParticle, const art::Ptr<simb::MCParticle> &mcParticle) const
{
    const auto pfWeight = this->GetWeight(pfParticle);

    if (pfWeight < std::numeric_limits<float>::epsilon())
    {
        // Input PFParticle has no associated hits
        return -1.f;
    }

    return this->GetMatchWeight(pfParticle, mcParticle) / pfWeight;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

float BacktrackHelper::BacktrackerData::GetMatchCompleteness(const art::Ptr<recob::PFParticle> &pfParticle, const art::Ptr<simb::MCParticle> &mcParticle) const
{
    const auto mcWeight = this->GetWeight(mcParticle);

    if (mcWeight < std::numeric_limits<float>::epsilon())
    {
        // Input MCParticle has no associated hits
        return -1.f;
    }

    return this->GetMatchWeight(pfParticle, mcParticle) / mcWeight;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

art::Ptr<simb::MCParticle> BacktrackHelper::BacktrackerData::GetBestMatchedMCParticle(const art::Ptr<recob::PFParticle> &pfParticle) const
{
    art::Ptr<simb::MCParticle> bestMCParticle;
    float bestWeight = -std::numeric_limits<float>::max();

    for (const auto &mcParticle : m_mcParticles)
    {
        const auto weight = this->GetMatchWeight(pfParticle, mcParticle);

        if (weight < bestWeight)
            continue;

        bestWeight = weight;
        bestMCParticle = mcParticle;
    }

    if (bestWeight < std::numeric_limits<float>::epsilon())
        throw cet::exception("BacktrackerData::GetBestMatchedMCParticle") << " - Input PFParticle matches to no MCParticles - comes from noise!" << std::endl;

    return bestMCParticle;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

PFParticleVector BacktrackHelper::BacktrackerData::GetBestMatchedPFParticles(const art::Ptr<simb::MCParticle> &mcParticle) const
{
    PFParticleVector pfParticles;

    for (const auto &pfParticle : m_pfParticles)
    {
        try
        {
            if (this->GetBestMatchedMCParticle(pfParticle) == mcParticle)
                pfParticles.push_back(pfParticle);
        }
        catch (const cet::exception &)
        {
        }
    }

    return pfParticles;
}

// -----------------------------------------------------------------------------------------------------------------------------------------
                
float BacktrackHelper::BacktrackerData::GetWeight(const art::Ptr<recob::Hit> &hit, const art::Ptr<simb::MCParticle> &mcParticle) const
{
    CollectionData<simb::MCParticle, float> mcParticleWeights;
    const bool hasMCParticles = this->CollectMCParticleWeights(hit, m_hitsToMcps, mcParticleWeights);

    if (!hasMCParticles)
        return 0.f;

    float weight = 0.f;
    for (const auto &entry : mcParticleWeights)
    {
        if (mcParticle == entry.first)
            weight += entry.second;
    }

    return weight;
}

// -----------------------------------------------------------------------------------------------------------------------------------------
                
float BacktrackHelper::BacktrackerData::GetWeight(const HitVector &hits, const art::Ptr<simb::MCParticle> &mcParticle) const
{
    float weight = 0.f;
    for (const auto &hit : hits)
    {
        weight += this->GetWeight(hit, mcParticle);
    }

    return weight;
}

// -----------------------------------------------------------------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------------------------------------------------------------
                
BacktrackHelper::SliceMetadata::SliceMetadata(const SliceVector &slices, const SlicesToBool &sliceToIsSelectedAsNu, const SlicesToHits &slicesToHits, const HitsToBool &hitsToIsNuInduced) :
    m_slices(slices),
    m_sliceToIsSelectedAsNu(sliceToIsSelectedAsNu),
    m_slicesToHits(slicesToHits),
    m_hitsToIsNuInduced(hitsToIsNuInduced)
{
}

// -----------------------------------------------------------------------------------------------------------------------------------------

SliceVector BacktrackHelper::SliceMetadata::GetSlices() const
{
    return m_slices;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

SliceVector BacktrackHelper::SliceMetadata::GetSelectedNeutrinoSlices() const
{
    SliceVector neutrinos;
    for (const auto &slice : m_slices)
    {
        const auto iter = m_sliceToIsSelectedAsNu.find(slice);
        if (iter == m_sliceToIsSelectedAsNu.end())
            throw cet::exception("SliceMetadata::GetSelectedNeutrinoSlices") << " - map entry for slice is missing." << std::endl;

        if (iter->second)
            neutrinos.push_back(slice);
    }

    return neutrinos;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

float BacktrackHelper::SliceMetadata::GetPurity(const art::Ptr<recob::Slice> &slice) const
{
    const auto nHits = this->GetNumberOfHits(slice);

    if (nHits == 0)
        throw cet::exception("SliceMetadata::GetPurity") << " - slice has no hits." << std::endl;

    return (static_cast<float>(this->GetNumberOfNuInducedHits(slice)) / static_cast<float>(nHits));
}

// -----------------------------------------------------------------------------------------------------------------------------------------
                
float BacktrackHelper::SliceMetadata::GetCompleteness(const art::Ptr<recob::Slice> &slice) const
{
    const auto nNuHits = this->GetTotalNumberOfNuInducedHits();

    if (nNuHits == 0)
        throw cet::exception("SliceMetadata::GetCompleteness") << " - there are no neutrino hits." << std::endl;

    return (static_cast<float>(this->GetNumberOfNuInducedHits(slice)) / static_cast<float>(nNuHits));
}

// -----------------------------------------------------------------------------------------------------------------------------------------

unsigned int BacktrackHelper::SliceMetadata::GetNumberOfHits(const art::Ptr<recob::Slice> &slice) const
{
    return CollectionHelper::GetManyAssociated(slice, m_slicesToHits).size();
}

// -----------------------------------------------------------------------------------------------------------------------------------------

HitVector BacktrackHelper::SliceMetadata::GetHits(const art::Ptr<recob::Slice> &slice) const
{
    return CollectionHelper::GetManyAssociated(slice, m_slicesToHits);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

unsigned int BacktrackHelper::SliceMetadata::GetNumberOfNuInducedHits(const art::Ptr<recob::Slice> &slice) const
{
    return this->CountNuInducedHits(CollectionHelper::GetManyAssociated(slice, m_slicesToHits));
}

// -----------------------------------------------------------------------------------------------------------------------------------------
                
unsigned int BacktrackHelper::SliceMetadata::GetTotalNumberOfNuInducedHits() const
{
    return this->CountNuInducedHits(this->GetAllHits());
}

// -----------------------------------------------------------------------------------------------------------------------------------------

art::Ptr<recob::Slice> BacktrackHelper::SliceMetadata::GetMostCompleteSlice() const
{
    if (this->GetTotalNumberOfNuInducedHits() == 0)
        throw cet::exception("SliceMetadata::GetMostCompleteSlice") << " - there are no neutrino hits! No slice is more complete than any other" << std::endl;

    if (m_slices.empty())
        throw cet::exception("SliceMetadata::GetMostCompleteSlice") << " - there are no slices to choose from." << std::endl;

    float bestCompleteness = -std::numeric_limits<float>::max();
    art::Ptr<recob::Slice> mostCompleteSlice;
    bool foundMostCompleteSlice = false;

    for (const auto &slice : m_slices)
    {
        const auto completeness = this->GetCompleteness(slice);
        if (completeness < bestCompleteness)
            continue;

        bestCompleteness = completeness;
        mostCompleteSlice = slice;
        foundMostCompleteSlice = true;
    }

    if (!foundMostCompleteSlice)
        throw cet::exception("SliceMetadata::GetMostCompleteSlice") << " - couldn't find the most complete slice" << std::endl;

    return mostCompleteSlice;
}

// -----------------------------------------------------------------------------------------------------------------------------------------
                
bool BacktrackHelper::SliceMetadata::IsMostCompleteSliceSelected() const
{
    const auto iter = m_sliceToIsSelectedAsNu.find(this->GetMostCompleteSlice());
    if (iter == m_sliceToIsSelectedAsNu.end())
        throw cet::exception("SliceMetadata::IsMostCompleteSliceSelected") << " - map entry for slice is missing." << std::endl;

    return iter->second;
}

// -----------------------------------------------------------------------------------------------------------------------------------------
                
HitVector BacktrackHelper::SliceMetadata::GetAllHits() const
{
    HitVector outputHits;

    for (const auto &slice : m_slices)
    {
        for (const auto &hit : CollectionHelper::GetManyAssociated(slice, m_slicesToHits))
        {
            outputHits.push_back(hit);
        }
    }

    return outputHits;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

unsigned int BacktrackHelper::SliceMetadata::CountNuInducedHits(const HitVector &hits) const
{
    unsigned int nNuInducedHits = 0;

    for (const auto &hit : hits)
    {
        const auto iter = m_hitsToIsNuInduced.find(hit);
        if (iter == m_hitsToIsNuInduced.end())
            throw cet::exception("SliceMetadata::CountNuInducedHits") << " - map entry for hit is missing." << std::endl;

        if (iter->second)
            nNuInducedHits++;
    }

    return nNuInducedHits;
}

} // namespace ubcc1pi
