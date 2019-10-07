/**
 *  @file  ubcc1pi/Producers/PFParticleHierarchy_module.cc
 *
 *  @brief The implementation of the PFParticle hierarchy producer
 */

#include "ubcc1pi/Producers/PFParticleHierarchy.h"

#include "ubcc1pi/Helpers/RecoHelper.h"
#include "ubcc1pi/Helpers/DebugHelper.h"

#include "art/Persistency/Common/PtrMaker.h"

namespace ubcc1pi
{

PFParticleHierarchy::PFParticleHierarchy(const art::EDProducer::Table<Config> &config) :
    m_config(config)
{
    // Produce a new set of PFParticles which are identical to those from Pandora but with updated parent/daughter IDs
    produces< std::vector<recob::PFParticle> >();

    // Produce a vertex for each PFParticle at the SCE corrected spacepoint closest to the match point at which the parent/daughter link was made
    produces< std::vector<recob::Vertex> >();
    produces< art::Assns<recob::PFParticle, recob::Vertex> >();

    // Produce a new metadata collection to store useful info about the parent-daughter matches
    //   - distToParent = closest SCE corrected spacepoint distance between this PFParticle and it's new parent 
    //   - distToOldParent = same as above but to the default parent assigned by Pandora
    produces< std::vector<larpandoraobj::PFParticleMetadata> >();
    produces< art::Assns<recob::PFParticle, larpandoraobj::PFParticleMetadata> >();
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void PFParticleHierarchy::produce(art::Event &event)
{
    // Create the output collections
    std::unique_ptr< std::vector<recob::PFParticle> > outputPFParticles( new std::vector<recob::PFParticle> );
    std::unique_ptr< std::vector<recob::Vertex> > outputVertices( new std::vector<recob::Vertex> );
    std::unique_ptr< art::Assns<recob::PFParticle, recob::Vertex> > outputParticlesToVertices( new art::Assns<recob::PFParticle, recob::Vertex> );
    std::unique_ptr< std::vector<larpandoraobj::PFParticleMetadata> > outputMetadata( new std::vector<larpandoraobj::PFParticleMetadata> );
    std::unique_ptr< art::Assns<recob::PFParticle, larpandoraobj::PFParticleMetadata> > outputParticlesToMetadata( new art::Assns<recob::PFParticle, larpandoraobj::PFParticleMetadata> );

    // Re-interpret the hierarchy
    const auto pfParticles = this->GetNeutrinoInducedPFParticles(event);
    if (!pfParticles.empty())
    {
        const auto pfParticleToPointsMap = this->GetPFParticlesTo3DPointsMap(event, pfParticles);

        PFParticlePairSeparationMap pfParticlePairSeparationMap;
        PFParticlePairMatchPositionMap pfParticlePairMatchPositionMap; 
        this->GetPFParticleMatches(pfParticles, pfParticleToPointsMap, pfParticlePairSeparationMap, pfParticlePairMatchPositionMap);

        // Initialize the match points and the orphans
        MatchPointVector matchPoints;
        auto orphans = pfParticles;
        this->SeedMatchPoints(pfParticleToPointsMap, orphans, matchPoints);

        // Now iteratively make matches
        while (!orphans.empty())
        {
            // Add orphans to existing match points if they are within the threshold
            this->AddOrphansToExistingMatchPoints(orphans, matchPoints);
            if (orphans.empty())
                break;

            // With the remaining orphans, choose the one with the next strongest match (either to an existing match point, or another
            // non-orphan PFParticle). Then make this match, adding a new match point if required
            this->MakeBestMatch(pfParticlePairSeparationMap, pfParticlePairMatchPositionMap, pfParticleToPointsMap, orphans, matchPoints);
        }

        // Get hierarchy from the match points and the new PFParticles and Vertices
        PFParticleToMatchPointMap pfParticleToMatchPointMap;
        IdToIdMap parentMap;
        this->GetParentMap(pfParticles, matchPoints, pfParticleToMatchPointMap, parentMap);

        // Fill the output vectors
        PFParticleToIndexMap pfParticleToNewIndicesMap;
        this->GetNewPFParticles(pfParticles, parentMap, outputPFParticles, pfParticleToNewIndicesMap);
        this->GetNewVertices(pfParticles, pfParticleToNewIndicesMap, pfParticleToMatchPointMap, pfParticleToPointsMap, outputVertices, outputParticlesToVertices, event);
        this->GetNewMetadata(pfParticles, pfParticleToNewIndicesMap, pfParticlePairSeparationMap, parentMap, outputMetadata, outputParticlesToMetadata, event);
    }

    // Add the output collections to the art event
    event.put(std::move(outputPFParticles));
    event.put(std::move(outputVertices));
    event.put(std::move(outputParticlesToVertices));
    event.put(std::move(outputMetadata));
    event.put(std::move(outputParticlesToMetadata));
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void PFParticleHierarchy::SeedMatchPoints(const PFParticlesToPoints &pfParticleToPointsMap, PFParticleVector &orphans, MatchPointVector &matchPoints) const
{
    // Add a seed match point at the neutrino vertex
    const auto neutrino = RecoHelper::GetNeutrino(orphans);
    const auto nuVertex = pfParticleToPointsMap.at(neutrino).front();

    const auto seedMatchPoint = std::make_shared<MatchPoint>(nuVertex, neutrino);
    seedMatchPoint->CachePFParticleDistances(orphans, pfParticleToPointsMap);
    matchPoints.push_back(seedMatchPoint);

    // Remove the neutrino from the orphan vector
    const auto iter = std::find(orphans.begin(), orphans.end(), neutrino);
    if (iter == orphans.end())
        throw cet::exception("PFParticleHierarchy::MakeBestMatch - scrambled list of orphans");

    orphans.erase(iter);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

PFParticleVector PFParticleHierarchy::GetNeutrinoInducedPFParticles(const art::Event &event) const
{
    const auto allPFParticles = CollectionHelper::GetCollection<recob::PFParticle>(event, m_config().PFParticleLabel());

    art::Ptr<recob::PFParticle> neutrino;
    try
    {
        neutrino = RecoHelper::GetNeutrino(allPFParticles);
    }
    catch (const cet::exception &)
    {
        return PFParticleVector();
    }

    const auto pfParticleMap = RecoHelper::GetPFParticleMap(allPFParticles);
    return RecoHelper::GetDownstreamParticles(neutrino, pfParticleMap);
}

// -----------------------------------------------------------------------------------------------------------------------------------------
        
PFParticleHierarchy::PFParticlesToPoints PFParticleHierarchy::GetPFParticlesTo3DPointsMap(const art::Event &event, const PFParticleVector &pfParticles) const
{
    PFParticlesToPoints pfParticleToPointsMap;

    // Nothing more to do if there aren't any PFParticles
    if (pfParticles.empty())
        return pfParticleToPointsMap;

    // Get the associated space-points and vertices
    const auto pfParticlesToSpacePoints = CollectionHelper::GetAssociation<recob::PFParticle, recob::SpacePoint>(event, m_config().PFParticleLabel(), m_config().SpacePointLabel());
    const auto pfParticlesToVertices = CollectionHelper::GetAssociation<recob::PFParticle, recob::Vertex>(event, m_config().PFParticleLabel(), m_config().VertexLabel());
    const auto pSpaceChargeService = RecoHelper::GetSpaceChargeService();
    
    const auto neutrino = RecoHelper::GetNeutrino(pfParticles);

    for (const auto &pfParticle : pfParticles)
    {
        auto points = PFParticlesToPoints::mapped_type();

        // For the neutrino PFParticle use the vertex as it's only point
        if (pfParticle == neutrino)
        {
            const auto vertex = CollectionHelper::GetSingleAssociated(pfParticle, pfParticlesToVertices);
            const auto correctedPoint = RecoHelper::CorrectForSpaceCharge(TVector3(vertex->position().X(), vertex->position().Y(), vertex->position().Z()), pSpaceChargeService);
            points.push_back(correctedPoint);
        }

        // Space-charge correct the spacepoints
        const auto spacePoints = CollectionHelper::GetManyAssociated(pfParticle, pfParticlesToSpacePoints);
        for (const auto &spacePoint : spacePoints)
        {
            const auto correctedPoint = RecoHelper::CorrectForSpaceCharge(TVector3(spacePoint->XYZ()[0], spacePoint->XYZ()[1], spacePoint->XYZ()[2]), pSpaceChargeService);
            points.push_back(correctedPoint);
        }
        
        if (!pfParticleToPointsMap.emplace(pfParticle, points).second)
        {
            throw cet::exception("PFParticleHierarchy::GetPFParticlesTo3DPointsMap - repeated input PFParticles");
        }
    }

    return pfParticleToPointsMap;
}

// -----------------------------------------------------------------------------------------------------------------------------------------
        
void PFParticleHierarchy::GetPFParticleMatches(const PFParticleVector &pfParticles, const PFParticlesToPoints &pfParticleToPointsMap, PFParticlePairSeparationMap &pfParticlePairSeparationMap, PFParticlePairMatchPositionMap &pfParticlePairMatchPositionMap) const
{
    for (unsigned int i = 0; i < pfParticles.size(); ++i)
    {
        const auto iPoints = pfParticleToPointsMap.at(pfParticles.at(i));
        for (unsigned int j = i + 1; j < pfParticles.size(); ++j)
        {
            const auto jPoints = pfParticleToPointsMap.at(pfParticles.at(j));

            // Find the closest two points
            float minSqrDist = std::numeric_limits<float>::max();
            unsigned int aBest = std::numeric_limits<unsigned int>::max();
            unsigned int bBest = std::numeric_limits<unsigned int>::max();

            for (unsigned int a = 0; a < iPoints.size(); ++a)
            {
                for (unsigned int b = 0; b < jPoints.size(); ++b)
                {
                    const float sqrDist = (iPoints.at(a) - jPoints.at(b)).Mag2();
                    if (sqrDist < minSqrDist)
                    {
                        minSqrDist = sqrDist;
                        aBest = a;
                        bBest = b;
                    }
                }
            }

            // Find the matching position
            const float minDist = std::pow(minSqrDist, 0.5f);
            TVector3 bestPosition(-std::numeric_limits<float>::max(), -std::numeric_limits<float>::max(), -std::numeric_limits<float>::max());
            if (!iPoints.empty() && !jPoints.empty())
            {
                bestPosition = (iPoints.at(aBest) + jPoints.at(bBest)) * 0.5f;
            }

            // Save the outcome
            pfParticlePairSeparationMap[pfParticles.at(i)].emplace_back(pfParticles.at(j), minDist);
            pfParticlePairSeparationMap[pfParticles.at(j)].emplace_back(pfParticles.at(i), minDist);
            
            pfParticlePairMatchPositionMap[pfParticles.at(i)].emplace_back(pfParticles.at(j), bestPosition);
            pfParticlePairMatchPositionMap[pfParticles.at(j)].emplace_back(pfParticles.at(i), bestPosition);
        }
    }
}

// -----------------------------------------------------------------------------------------------------------------------------------------

art::Ptr<recob::PFParticle> PFParticleHierarchy::GetClosestPFParticle(const art::Ptr<recob::PFParticle> &orphan, const PFParticleVector &orphans, const PFParticlePairSeparationMap &pfParticlePairSeparationMap, float &minDist) const
{
    minDist = std::numeric_limits<float>::max();
    art::Ptr<recob::PFParticle> closestPFParticle;
    bool foundParticle = false;

    for (const auto &entry : pfParticlePairSeparationMap.at(orphan))
    {
        // Don't find an orphan PFParticle
        if (std::find(orphans.begin(), orphans.end(), entry.first) != orphans.end())
            continue;

        if (entry.second < minDist)
        {
            minDist = entry.second;
            closestPFParticle = entry.first;
            foundParticle = true;
        }
    }

    if (!foundParticle)
        throw cet::exception("PFParticleHierarchy::GetClosestPFParticle - there are no connected particles in the input map");

    return closestPFParticle;
}

// -----------------------------------------------------------------------------------------------------------------------------------------
        
std::shared_ptr<PFParticleHierarchy::MatchPoint> PFParticleHierarchy::GetClosestMatchPoint(const art::Ptr<recob::PFParticle> &orphan, const MatchPointVector &matchPoints) const
{
    if (matchPoints.empty())
        throw cet::exception("PFParticleHierarchy::GetClosestMatchPoint - input vector of match points is empty");

    float minDist = std::numeric_limits<float>::max();
    std::shared_ptr<MatchPoint> bestPoint;

    for (const auto &matchPoint : matchPoints)
    {
        const float dist = matchPoint->GetDistance(orphan);
        if (dist < minDist)
        {
            minDist = dist;
            bestPoint = matchPoint;
        }
    }

    return bestPoint;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void PFParticleHierarchy::AddOrphansToExistingMatchPoints(PFParticleVector &orphans, MatchPointVector &matchPoints) const
{
    const float thresh = m_config().DistanceThreshold();

    PFParticleVector newConnections;
    for (const auto &orphan : orphans)
    {
        const auto closestPoint = this->GetClosestMatchPoint(orphan, matchPoints);
        const auto dist = closestPoint->GetDistance(orphan);

        if (dist > thresh)
            continue;

        closestPoint->AddDaughter(orphan);
        newConnections.push_back(orphan);
    }

    // Remove the new connections from the list of orphans
    for (const auto &pfParticle : newConnections)
    {
        const auto iter = std::find(orphans.begin(), orphans.end(), pfParticle);
        if (iter == orphans.end())
            throw cet::exception("PFParticleHierarchy::AddOrphansToExistingMatchPoints - scrambled list of orphans");

        orphans.erase(iter);
    }
}
        
// -----------------------------------------------------------------------------------------------------------------------------------------

void PFParticleHierarchy::MakeBestMatch(const PFParticlePairSeparationMap &pfParticlePairSeparationMap, const PFParticlePairMatchPositionMap &pfParticlePairMatchPositionMap, const PFParticlesToPoints &pfParticleToPointsMap, PFParticleVector &orphans, MatchPointVector &matchPoints) const
{
    // Check if there is anything to do
    if (orphans.empty())
        throw cet::exception("PFParticleHierarchy::MakeBestMatch - Input vector of orphans is empty");

    art::Ptr<recob::PFParticle> bestOrphan, bestParent;
    float minDist = std::numeric_limits<float>::max();
    bool shouldMakeNewMatchPoint = false;
    TVector3 newMatchPointPosition(-std::numeric_limits<float>::max(), -std::numeric_limits<float>::max(), -std::numeric_limits<float>::max());
    std::shared_ptr<MatchPoint> bestMatchPoint;

    // Find the orphan which is closest to either an existing match point, or another non-orphan PFParticle
    for (const auto &orphan : orphans)
    {
        // Get the closest match point
        const auto closestMatchPoint = this->GetClosestMatchPoint(orphan, matchPoints);
        const auto matchPointDist = closestMatchPoint->GetDistance(orphan);

        // Get the closest PFParticle
        float pfParticleDist = std::numeric_limits<float>::max();
        const auto closestPFParticle = this->GetClosestPFParticle(orphan, orphans, pfParticlePairSeparationMap, pfParticleDist);

        // If this isn't the smallest distance yet, then continue to the next orphan
        if (matchPointDist > minDist && pfParticleDist > minDist)
            continue;

        // This is the best orphan so far
        bestOrphan = orphan;

        // Determine if it's better to link to an existing match point or make a new one
        shouldMakeNewMatchPoint = (pfParticleDist < matchPointDist);

        if (shouldMakeNewMatchPoint)
        {
            // It's best to make a new match point, so get the required details
            minDist = pfParticleDist;
            bestParent = closestPFParticle;

            // Look up the new match point position in the map
            bool foundNewMatchPointPosition = false;
            for (const auto &entry : pfParticlePairMatchPositionMap.at(orphan))
            {
                if (entry.first != closestPFParticle)
                    continue;

                newMatchPointPosition = entry.second;
                foundNewMatchPointPosition = true;
                break;
            }

            if (!foundNewMatchPointPosition)
                throw cet::exception("PFParticleHierarchy::MakeBestMatch - Can't find new match point position in input map");
        }
        else
        {
            // It's best to link to an existing match point
            minDist = matchPointDist;
            bestMatchPoint = closestMatchPoint;
        }
    }

    // Make a new match point if required
    if (shouldMakeNewMatchPoint)
    {
        bestMatchPoint = std::make_shared<MatchPoint>(newMatchPointPosition, bestParent);
        bestMatchPoint->CachePFParticleDistances(orphans, pfParticleToPointsMap);

        matchPoints.push_back(bestMatchPoint);
    }
    
    // Add the orphan to the best match point
    bestMatchPoint->AddDaughter(bestOrphan);

    // Remove the orphan from the vector
    const auto iter = std::find(orphans.begin(), orphans.end(), bestOrphan);
    if (iter == orphans.end())
        throw cet::exception("PFParticleHierarchy::MakeBestMatch - scrambled list of orphans");

    orphans.erase(iter);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void PFParticleHierarchy::GetParentMap(const PFParticleVector &pfParticles, const MatchPointVector &matchPoints, PFParticleToMatchPointMap &pfParticleToMatchPointMap, IdToIdMap &parentMap) const
{
    // Extract the daughter -> parent links from the match points
    std::shared_ptr<MatchPoint> neutrinoMatchPoint;
    bool hasNeutrinoMatchPoint = false;

    for (const auto &matchPoint : matchPoints)
    {
        const auto parent = matchPoint->GetParent();
        const auto parentID = parent->Self();

        // Find the neutrino match point
        if (RecoHelper::IsNeutrino(parent))
        {
            neutrinoMatchPoint = matchPoint;
            hasNeutrinoMatchPoint = true;
        }

        // Add the entries in the maps for the daughters
        for (const auto &daughter : matchPoint->GetDaughters())
        {
            const auto daughterID = daughter->Self();

            if (!parentMap.emplace(daughterID, parentID).second)
                throw cet::exception("PFParticleHierarchy::GetParentMap - bad hierarchy, found particle with multiple parents");
            
            if (!pfParticleToMatchPointMap.emplace(daughter, matchPoint).second)
                throw cet::exception("PFParticleHierarchy::GetParentMap - bad hierarchy, found particle with multiple parents");
        }
    }

    if (!hasNeutrinoMatchPoint)
        throw cet::exception("PFParticleHierarchy::GetParentMap - no match point found with the neutrino a the parent");

    // Sanity check that every PFParticle has been included, and add the neutrino
    for (const auto &pfParticle : pfParticles)
    {
        const bool isInMap = (parentMap.find(pfParticle->Self()) != parentMap.end());
        const bool isNeutrino = RecoHelper::IsNeutrino(pfParticle);

        if (isNeutrino)
        {
            if (isInMap)
                throw cet::exception("PFParticleHierarchy::GetParentMap - the neutrino PFParticle has been added as a daughter!");

            parentMap.emplace(pfParticle->Self(), recob::PFParticle::kPFParticlePrimary);
            pfParticleToMatchPointMap.emplace(pfParticle, neutrinoMatchPoint);
        }
        else if (!isInMap)
        {
            throw cet::exception("PFParticleHierarchy::GetParentMap - found PFParticle that's not been added to the hierarchy");
        }
    }
}

// -----------------------------------------------------------------------------------------------------------------------------------------
        
void PFParticleHierarchy::GetNewPFParticles(const PFParticleVector &pfParticles, const IdToIdMap &parentMap, std::unique_ptr< std::vector<recob::PFParticle> > &outputPFParticles, PFParticleToIndexMap &pfParticleToNewIndicesMap) const
{
    // Get the reverse daughter mapping
    std::unordered_map<size_t, std::vector<size_t> > daughterMap;
    for (const auto &entry : parentMap)
        daughterMap[entry.second].push_back(entry.first);

    size_t newId = 0;
    for (const auto &pfParticle : pfParticles)
    {
        const auto pdg = pfParticle->PdgCode();
        const auto self = pfParticle->Self();
        const auto parent = parentMap.at(self);
        auto daughters = std::vector<size_t>();

        const auto daughtersIter = daughterMap.find(self);
        if (daughtersIter != daughterMap.end())
            daughters = daughtersIter->second;

        outputPFParticles->emplace_back(pdg, self, parent, daughters);
        pfParticleToNewIndicesMap.emplace(pfParticle, newId);
        ++newId;
    }
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void PFParticleHierarchy::GetNewVertices(const PFParticleVector &pfParticles, const PFParticleToIndexMap &pfParticleToNewIndicesMap, const PFParticleToMatchPointMap &pfParticleToMatchPointMap, const PFParticlesToPoints &pfParticleToPointsMap, std::unique_ptr< std::vector<recob::Vertex> > &outputVertices, std::unique_ptr< art::Assns<recob::PFParticle, recob::Vertex> > &outputParticlesToVertices, art::Event &event) const
{
    for (const auto &pfParticle : pfParticles)
    {
        const auto newPFParticleIndex = pfParticleToNewIndicesMap.at(pfParticle);
        const auto matchPoint = pfParticleToMatchPointMap.at(pfParticle);

        // Find the spacepoint that is closest to the match point
        float minSqrDist = std::numeric_limits<float>::max();
        TVector3 closestPoint(-std::numeric_limits<float>::max(), -std::numeric_limits<float>::max(), -std::numeric_limits<float>::max());

        for (const auto &point : pfParticleToPointsMap.at(pfParticle))
        {
            const float sqrDist = (matchPoint->GetPosition() - point).Mag2();
            if (sqrDist > minSqrDist)
                continue;

            minSqrDist = sqrDist;
            closestPoint = point;
        }

        // Make the vertex
        double xyz[3] = {closestPoint.X(), closestPoint.Y(), closestPoint.Z()};
        const int id = static_cast<int>(outputVertices->size());
        outputVertices->emplace_back(xyz, id);

        // Make the association
        const art::PtrMaker<recob::PFParticle> makePtrPFP(event);
        const art::PtrMaker<recob::Vertex> makePtrVertex(event);
        outputParticlesToVertices->addSingle(makePtrPFP(newPFParticleIndex), makePtrVertex(id));
    }
}

// -----------------------------------------------------------------------------------------------------------------------------------------
        
void PFParticleHierarchy::GetNewMetadata(const PFParticleVector &pfParticles, const PFParticleToIndexMap &pfParticleToNewIndicesMap, const PFParticlePairSeparationMap &pfParticlePairSeparationMap, const IdToIdMap &parentMap, std::unique_ptr< std::vector<larpandoraobj::PFParticleMetadata> > &outputMetadata, std::unique_ptr< art::Assns<recob::PFParticle, larpandoraobj::PFParticleMetadata> > &outputParticlesToMetadata, art::Event &event) const
{
    const auto pfParticleMap = RecoHelper::GetPFParticleMap(pfParticles);
    for (const auto &pfParticle : pfParticles)
    {
        const auto newPFParticleIndex = pfParticleToNewIndicesMap.at(pfParticle);
        larpandoraobj::PFParticleMetadata::PropertiesMap properties;

        if (!RecoHelper::IsNeutrino(pfParticle))
        {
            const auto parent = pfParticleMap.at(parentMap.at(pfParticle->Self()));
            const auto oldParent = pfParticleMap.at(pfParticle->Parent());
            float distToParent = -std::numeric_limits<float>::max();
            float distToOldParent = -std::numeric_limits<float>::max();

            for (const auto &entry : pfParticlePairSeparationMap.at(pfParticle))
            {
                if (entry.first == parent)
                    distToParent = entry.second;

                if (entry.first == oldParent)
                    distToOldParent = entry.second;
            }

            properties.emplace("distToParent", distToParent);
            properties.emplace("distToOldParent", distToOldParent);
        }


        const art::PtrMaker<recob::PFParticle> makePtrPFP(event);
        const art::PtrMaker<larpandoraobj::PFParticleMetadata> makePtrMetadata(event);
        outputParticlesToMetadata->addSingle(makePtrPFP(newPFParticleIndex), makePtrMetadata(outputMetadata->size()));
        outputMetadata->emplace_back(properties);
    }
}

// -----------------------------------------------------------------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------------------------------------------------------------

PFParticleHierarchy::MatchPoint::MatchPoint(const TVector3 &position, const art::Ptr<recob::PFParticle> &parent) : 
    m_position(position),
    m_parent(parent)
{
}

// -----------------------------------------------------------------------------------------------------------------------------------------
                
void PFParticleHierarchy::MatchPoint::AddDaughter(const art::Ptr<recob::PFParticle> &daughter)
{
    m_daughters.push_back(daughter);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

art::Ptr<recob::PFParticle> PFParticleHierarchy::MatchPoint::GetParent() const
{
    return m_parent;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

PFParticleVector PFParticleHierarchy::MatchPoint::GetDaughters() const
{
    return m_daughters;
}

// -----------------------------------------------------------------------------------------------------------------------------------------
                
void PFParticleHierarchy::MatchPoint::CachePFParticleDistances(const PFParticleVector &orphans, const PFParticlesToPoints &pfParticleToPointsMap)
{
    for (const auto &pfParticle : orphans)
    {
        // Get the minimum distance between the PFParticle and this point
        float minSqrDist = std::numeric_limits<float>::max();

        for (const auto &point : pfParticleToPointsMap.at(pfParticle))
        {
            minSqrDist = std::min(minSqrDist, static_cast<float>((point - m_position).Mag2()));
        }

        const float minDist = std::pow(minSqrDist, 0.5f);
        m_pfParticleDistances.emplace(pfParticle, minDist);
    }
}

// -----------------------------------------------------------------------------------------------------------------------------------------
            
float PFParticleHierarchy::MatchPoint::GetDistance(const art::Ptr<recob::PFParticle> &orphan) const
{
    const auto iter = m_pfParticleDistances.find(orphan);
    if (iter == m_pfParticleDistances.end())
        throw cet::exception("PFParticleHierarchy::MatchPoint::GetDistance - distanc to input PFParticle isn't cahced");

    return iter->second;
}

// -----------------------------------------------------------------------------------------------------------------------------------------
                
TVector3 PFParticleHierarchy::MatchPoint::GetPosition() const
{
    return m_position;
}

} // namespace ubcc1pi
