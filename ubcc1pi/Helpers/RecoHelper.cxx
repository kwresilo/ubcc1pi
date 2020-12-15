/**
 *  @file  ubcc1pi/Helpers/RecoHelper.cxx
 *
 *  @brief The implementation file for the reco helper class
 */

#include "ubcc1pi/Helpers/RecoHelper.h"

namespace ubcc1pi
{

PFParticleMap RecoHelper::GetPFParticleMap(const PFParticleVector &allPFParticles)
{
    PFParticleMap pfParticleMap;

    for (const auto &pfParticle : allPFParticles)
    {
        if (!pfParticleMap.emplace(pfParticle->Self(), pfParticle).second)
            throw cet::exception("RecoHelper::GetPFParticleMap") << " - Found repeated PFParticle with Self = " << pfParticle->Self() << "." << std::endl;
    }

    return pfParticleMap;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

bool RecoHelper::IsNeutrino(const art::Ptr<recob::PFParticle> &pfParticle)
{
    const auto pdg = std::abs(pfParticle->PdgCode());
    return (pdg == 12 ||  // Nue
            pdg == 14 ||  // Numu
            pdg == 16 );  // Nutau
}

// -----------------------------------------------------------------------------------------------------------------------------------------

PFParticleVector RecoHelper::GetNeutrinos(const PFParticleVector &allPFParticles)
{
    PFParticleVector neutrinos;

    for (const auto &particle : allPFParticles)
    {
        if (RecoHelper::IsNeutrino(particle))
            neutrinos.push_back(particle);
    }

    return neutrinos;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

art::Ptr<recob::PFParticle> RecoHelper::GetNeutrino(const PFParticleVector &allPFParticles)
{
    const auto &neutrinos = RecoHelper::GetNeutrinos(allPFParticles);

    if (neutrinos.empty())
        throw cet::exception("RecoHelper::GetNeutrino") << " - Didn't find a neutrino PFParticle." << std::endl;

    if (neutrinos.size() > 1)
        throw cet::exception("RecoHelper::GetNeutrino") << " - Found multiple neutrino PFParticles." << std::endl;

    return neutrinos.front();
}

// -----------------------------------------------------------------------------------------------------------------------------------------

PFParticleVector RecoHelper::GetNeutrinoFinalStates(const PFParticleVector &allPFParticles)
{
    // Get the neutrino PFParticle - if there is one!
    art::Ptr<recob::PFParticle> neutrino;
    try
    {
        neutrino = RecoHelper::GetNeutrino(allPFParticles);
    }
    catch (const cet::exception &)
    {
        // No neutrino found
        return PFParticleVector();
    }

    const auto pfParticleMap = RecoHelper::GetPFParticleMap(allPFParticles);
    return RecoHelper::GetDaughters(neutrino, pfParticleMap);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

TVector3 RecoHelper::GetRecoNeutrinoVertex(const art::Event &event, const PFParticleVector &allPFParticles, const art::InputTag &vertexLabel)
{
    const auto pfpToVertex = CollectionHelper::GetAssociation<recob::PFParticle, recob::Vertex>(event, vertexLabel);

    try
    {
        const auto neutrino = RecoHelper::GetNeutrino(allPFParticles);
        const auto vertex = CollectionHelper::GetSingleAssociated(neutrino, pfpToVertex);

        return TVector3(vertex->position().X(), vertex->position().Y(), vertex->position().Z());
    }
    catch (const cet::exception &)
    {
        return TVector3(-std::numeric_limits<float>::max(), -std::numeric_limits<float>::max(), -std::numeric_limits<float>::max());
    }
}

// -----------------------------------------------------------------------------------------------------------------------------------------

float RecoHelper::GetRecoFlashChi2(const art::Event &event, const PFParticleVector &allPFParticles, const art::InputTag &pfparticleLabel, const art::InputTag &flashmatchLabel)
{
    const auto nuFlashScoreAssoc = CollectionHelper::GetAssociation<recob::PFParticle, anab::T0>(event, pfparticleLabel, flashmatchLabel);

    try
    {
        const auto neutrino = RecoHelper::GetNeutrino(allPFParticles);

        const auto T0_flashchi = CollectionHelper::GetSingleAssociated(neutrino, nuFlashScoreAssoc);

        return (float)T0_flashchi->TriggerConfidence();
    }
    catch (const cet::exception &)
    {
        return -std::numeric_limits<float>::max();
    }
}

// -----------------------------------------------------------------------------------------------------------------------------------------

art::Ptr<recob::OpFlash> RecoHelper::GetLargestFlash(const art::Event &event, const art::InputTag &flashLabel)
{
    float flash_pe = -std::numeric_limits<float>::max();
    int f_largest = -1;

    try
    {
        const auto flashes = CollectionHelper::GetCollection<recob::OpFlash>(event, flashLabel);

        for (size_t f=0; f< flashes.size(); f++){
            auto flash = flashes->at(f);

            if (flash.TotalPE() > flash_pe){
                flash_pe = flash.TotalPE();
                f_largest = (int)f;
            }
        } // end loop over f in flashes

        return flashes.at(f_largest);
    }
    catch (const cet::exception &)
    {
        throw cet::exception("RecoHelper::GetLargestFlash") << " - Didn't find a largest flash." << std::endl;
    }
}

// -----------------------------------------------------------------------------------------------------------------------------------------

art::Ptr<recob::PFParticle> RecoHelper::GetParent(const art::Ptr<recob::PFParticle> &particle, const PFParticleMap &pfParticleMap)
{
    if (particle->IsPrimary())
        throw cet::exception("RecoHelper::GetParent") << " - PFParticle is primary, so doesn't have a parent." << std::endl;

    const auto parentIter = pfParticleMap.find(particle->Parent());
    if (parentIter == pfParticleMap.end())
        throw cet::exception("RecoHelper::GetParent") << " - Couldn't find parent PFParticle in hierarchy." << std::endl;

    return parentIter->second;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

PFParticleVector RecoHelper::GetDaughters(const art::Ptr<recob::PFParticle> &particle, const PFParticleMap &pfParticleMap)
{
    PFParticleVector daughters;

    for (int i = 0; i < particle->NumDaughters(); ++i)
    {
        const auto daughterIter = pfParticleMap.find(particle->Daughter(i));
        if (daughterIter == pfParticleMap.end())
            throw cet::exception("RecoHelper::GetDaughter") << " - Couldn't find daughter PFParticle in hierarchy." << std::endl;

        daughters.push_back(daughterIter->second);
    }

    return daughters;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

PFParticleVector RecoHelper::GetDownstreamParticles(const art::Ptr<recob::PFParticle> &particle, const PFParticleMap &pfParticleMap)
{
    PFParticleVector downstreamParticles;
    RecoHelper::GetDownstreamParticles(particle, pfParticleMap, downstreamParticles);

    return downstreamParticles;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void RecoHelper::GetDownstreamParticles(const art::Ptr<recob::PFParticle> &particle, const PFParticleMap &pfParticleMap, PFParticleVector &downstreamParticles)
{
    downstreamParticles.push_back(particle);

    for (const auto &daughter : RecoHelper::GetDaughters(particle, pfParticleMap))
        RecoHelper::GetDownstreamParticles(daughter, pfParticleMap, downstreamParticles);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

unsigned int RecoHelper::GetGeneration(const art::Ptr<recob::PFParticle> &particle, const PFParticleMap &pfParticleMap)
{
    unsigned int generation = 0;

    art::Ptr<recob::PFParticle> nextParticle = particle;
    while (!nextParticle->IsPrimary())
    {
        nextParticle = RecoHelper::GetParent(nextParticle, pfParticleMap);
        generation++;
    }

    return generation;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

bool RecoHelper::IsSliceSelectedAsNu(const art::Ptr<recob::Slice> &slice, const SlicesToPFParticles &slicesToPFParticles)
{
    const auto neutrinos = RecoHelper::GetNeutrinos(CollectionHelper::GetManyAssociated(slice, slicesToPFParticles));

    if (neutrinos.size() > 1)
        throw cet::exception("RecoHelper::IsSliceSelectedAsNu") << " - Found multiple neutrino PFParticles associated to slice." << std::endl;

    return (!neutrinos.empty());
}

// -----------------------------------------------------------------------------------------------------------------------------------------

bool RecoHelper::HasMetadataValue(const art::Ptr<larpandoraobj::PFParticleMetadata> &metadata, const std::string &property)
{
    const auto &properties = metadata->GetPropertiesMap();
    return (properties.find(property) != properties.end());
}

// -----------------------------------------------------------------------------------------------------------------------------------------

float RecoHelper::GetMetadataFloat(const art::Ptr<larpandoraobj::PFParticleMetadata> &metadata, const std::string &property)
{
    if (!RecoHelper::HasMetadataValue(metadata, property))
        throw cet::exception("RecoHelper::GetMetadataValue") << " - Can't find metadata with key: " << property << std::endl;

    return metadata->GetPropertiesMap().at(property);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

int RecoHelper::GetMetadataInt(const art::Ptr<larpandoraobj::PFParticleMetadata> &metadata, const std::string &property)
{
    return static_cast<int>(std::round(RecoHelper::GetMetadataFloat(metadata, property)));
}

// -----------------------------------------------------------------------------------------------------------------------------------------

bool RecoHelper::GetMetadataBool(const art::Ptr<larpandoraobj::PFParticleMetadata> &metadata, const std::string &property)
{
    switch (RecoHelper::GetMetadataInt(metadata, property))
    {
        case 0:
            return false;
        case 1:
            return true;
        default:
            throw cet::exception("RecoHelper::GetMetadataValue") << " - Can't interpret metadata: " << property << " as boolean." << std::endl;
    }
}

// -----------------------------------------------------------------------------------------------------------------------------------------

float RecoHelper::GetTrackScore(const art::Ptr<larpandoraobj::PFParticleMetadata> &metadata)
{
    return RecoHelper::GetMetadataFloat(metadata, "TrackScore");
}

// -----------------------------------------------------------------------------------------------------------------------------------------

float RecoHelper::GetTopologicalScore(const art::Ptr<larpandoraobj::PFParticleMetadata> &metadata)
{
    return RecoHelper::GetMetadataFloat(metadata, "NuScore");
}

// -----------------------------------------------------------------------------------------------------------------------------------------

unsigned int RecoHelper::CountHitsInView(const HitVector &hits, const geo::View_t &view)
{
    unsigned int count = 0;
    for (const auto &hit : hits)
    {
        if (hit->View() == view)
            count++;
    }

    return count;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

const spacecharge::SpaceChargeService::provider_type *const RecoHelper::GetSpaceChargeService()
{
    return lar::providerFrom<spacecharge::SpaceChargeService>();
}

// -----------------------------------------------------------------------------------------------------------------------------------------

TVector3 RecoHelper::CorrectForSpaceCharge(const TVector3 &position, const spacecharge::SpaceChargeService::provider_type *const pSpaceChargeService)
{
    const auto offset = pSpaceChargeService->GetCalPosOffsets(geo::Point_t(position.X(), position.Y(), position.Z()));
    return TVector3(position.X() - offset.X(), position.Y() + offset.Y(), position.Z() + offset.Z());
}

// -----------------------------------------------------------------------------------------------------------------------------------------

float RecoHelper::GetPidScore(const art::Ptr<anab::ParticleID> &pid, const std::function<bool(const anab::sParticleIDAlgScores &)> &fCriteria)
{
    bool found = false;
    float score = -std::numeric_limits<float>::max();

    for (const auto &algo : pid->ParticleIDAlgScores())
    {
        if (!fCriteria(algo))
            continue;

        if (found)
            throw cet::exception("RecoHelper::GetPidScore") << " - Ambiguous criteria supplied." << std::endl;

        found = true;
        score = algo.fValue;
    }

    return score;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

geo::View_t RecoHelper::GetView(const std::bitset<8> &planeMask)
{
    const bool usesW = planeMask.test(2);
    const bool usesU = planeMask.test(0);
    const bool usesV = planeMask.test(1);

    if (usesW && !usesU && !usesV)
        return geo::kW;

    if (!usesW && usesU && !usesV)
        return geo::kU;

    if (!usesW && !usesU && usesV)
        return geo::kV;

    return geo::kUnknown;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

float RecoHelper::GetBraggLikelihood(const art::Ptr<anab::ParticleID> &pid, const int &pdg, const geo::View_t &view, const anab::kTrackDir &dir)
{
    return RecoHelper::GetPidScore(pid, [&](const anab::sParticleIDAlgScores &algo) -> bool {
        return (algo.fAlgName == "BraggPeakLLH"              &&
                algo.fTrackDir == dir                        &&
                algo.fAssumedPdg == pdg                      &&
                RecoHelper::GetView(algo.fPlaneMask) == view );
    });
}

// -----------------------------------------------------------------------------------------------------------------------------------------

float RecoHelper::GetBraggLikelihood(const art::Ptr<anab::ParticleID> &pid, const int &pdg, const anab::kTrackDir &dir, const float yzAngle, const float sin2AngleThreshold, const unsigned int nHitsU, const unsigned int nHitsV)
{
    const auto piBy3 = std::acos(0.5f);

    const auto likelihoodW = RecoHelper::GetBraggLikelihood(pid, pdg, geo::kW, dir);
    const bool isTrackAlongWWire = (std::pow(std::sin(yzAngle), 2) < sin2AngleThreshold);
    const auto hasW = (likelihoodW > -1.f && !isTrackAlongWWire);

    if (hasW)
        return likelihoodW;

    // Otherwise the track is along the W wire direction, so just average the other two planes weighted by the number of degrees of freedom
    const auto likelihoodU = RecoHelper::GetBraggLikelihood(pid, pdg, geo::kU, dir);
    const bool isTrackAlongUWire = (std::pow(std::sin(yzAngle - piBy3), 2) < sin2AngleThreshold);
    const auto hasU = (likelihoodU > -1.f && !isTrackAlongUWire);

    const auto likelihoodV = RecoHelper::GetBraggLikelihood(pid, pdg, geo::kV, dir);
    const bool isTrackAlongVWire = (std::pow(std::sin(yzAngle + piBy3), 2) < sin2AngleThreshold);
    const auto hasV = (likelihoodV > -1.f && !isTrackAlongVWire);

    const auto dofU = hasU ? static_cast<float>(nHitsU) : 0.f;
    const auto dofV = hasV ? static_cast<float>(nHitsV) : 0.f;
    const auto dofUV = dofU + dofV;

    const auto weightU = hasU ? (likelihoodU * dofU) : 0.f;
    const auto weightV = hasV ? (likelihoodV * dofV) : 0.f;

    if ((!hasU && !hasV) || (nHitsU == 0 && nHitsV == 0) || dofUV <= std::numeric_limits<float>::epsilon())
        return -std::numeric_limits<float>::max();

    return (weightU + weightV) / dofUV;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

std::vector<size_t> RecoHelper::GetValidPoints(const art::Ptr<recob::Track> &track)
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

float RecoHelper::GetRange(const art::Ptr<recob::Track> &track, const spacecharge::SpaceChargeService::provider_type *const pSpaceChargeService)
{
    const auto validPoints = RecoHelper::GetValidPoints(track);
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

void RecoHelper::GetDistanceToPoint(const art::Ptr<recob::Track> &track, const TVector3 &point, const spacecharge::SpaceChargeService::provider_type *const pSpaceChargeService, float &transverseDist, float &longitudinalDist)
{
    const auto start = RecoHelper::CorrectForSpaceCharge(TVector3(track->Start().X(), track->Start().Y(), track->Start().Z()), pSpaceChargeService);
    const auto direction = TVector3(track->StartDirection().X(), track->StartDirection().Y(), track->StartDirection().Z());
    const auto position = RecoHelper::CorrectForSpaceCharge(point, pSpaceChargeService);
    const auto dist2 = (start - position).Mag2();

    longitudinalDist = direction.Dot(start - position);
    transverseDist = std::sqrt(dist2 - std::pow(longitudinalDist, 2));
}

// -----------------------------------------------------------------------------------------------------------------------------------------

float RecoHelper::GetTrackWiggliness(const art::Ptr<recob::Track> &track)
{
    const auto validPoints = RecoHelper::GetValidPoints(track);
    if (validPoints.size() < 3)
        return 0.f;

    std::vector<float> thetaVector;
    float thetaSum = 0.f;
    for (unsigned int i = 1; i < validPoints.size(); ++i)
    {
        const auto dir = track->DirectionAtPoint(validPoints.at(i));
        const auto dirPrev = track->DirectionAtPoint(validPoints.at(i - 1));

        // Bind between -1 and 1 at floating precision to avoid issues with cast from double
        const auto cosTheta = std::min(1.f, std::max(-1.f, static_cast<float>(dir.Dot(dirPrev))));
        const auto theta = std::acos(cosTheta);

        thetaSum += theta;
        thetaVector.push_back(theta);
    }

    const auto thetaMean = thetaSum / static_cast<float>(thetaVector.size());

    float thetaDiffSum = 0.f;
    for (const auto &theta : thetaVector)
    {
        thetaDiffSum += std::pow(theta - thetaMean, 2);
    }

    const auto variance = thetaDiffSum / static_cast<float>(thetaVector.size() - 1);

    return std::sqrt(variance);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

unsigned int RecoHelper::CountSpacePointsNearTrackEnd(const art::Ptr<recob::Track> &track, const SpacePointVector &spacePoints, const float distance, const spacecharge::SpaceChargeService::provider_type *const pSpaceChargeService)
{
    const auto distSquared = distance * distance;
    const auto end = RecoHelper::CorrectForSpaceCharge(TVector3(track->End().X(), track->End().Y(), track->End().Z()), pSpaceChargeService);

    unsigned int nPoints = 0;
    for (const auto &spacePoint : spacePoints)
    {
        const auto spacePointCorrected = RecoHelper::CorrectForSpaceCharge(TVector3(spacePoint->XYZ()[0], spacePoint->XYZ()[1], spacePoint->XYZ()[2]), pSpaceChargeService);

        if ((end - spacePointCorrected).Mag2() < distSquared)
            nPoints++;
    }

    return nPoints;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

float RecoHelper::GetTruncatedMeandEdxAtTrackStart(const std::vector<float> dedxPerHit, const std::vector<float> &residualRangePerHit, const unsigned int nHitsToSkip, const float lengthFraction)
{
    if (dedxPerHit.size() != residualRangePerHit.size())
        throw cet::exception("RecoHelper::GetTruncatedMeandEdxAtTrackStart") << " - dEdx per hit and residual range vectors have different sizes" << std::endl;

    // Check if the variable is calculable
    if (residualRangePerHit.size() <= nHitsToSkip)
        return -std::numeric_limits<float>::max();

    // Make make the vector of pairs to keep track of the incides, and find the maximum residual range
    std::vector<std::pair<float, unsigned int> > residualRangeIndices;
    float maxResidualRange = -std::numeric_limits<float>::max();
    for (unsigned int i = 0; i < residualRangePerHit.size(); ++i)
    {
        const auto residualRange = residualRangePerHit.at(i);
        maxResidualRange = std::max(maxResidualRange, residualRange);

        residualRangeIndices.emplace_back(residualRange, i);
    }

    const auto residualRangeCutoff = maxResidualRange * lengthFraction;

    // Sort the residual ranges such that the largest residual range (closest to the start of the track) is first
    std::sort(residualRangeIndices.begin(), residualRangeIndices.end(), [](auto &a, auto &b) {
        return a.first > b.first;
    });

    // Get the dEdx of the hits at the start of the track
    std::vector<float> dedxPerHitAtStart;
    for (unsigned int i = nHitsToSkip; i < residualRangeIndices.size(); ++i)
    {
        const auto entry = residualRangeIndices.at(i);
        const auto residualRange = entry.first;
        const auto hitIndex = entry.second;

        // ATTN small residual ranges are at the start of the track
        if (residualRange < residualRangeCutoff)
            continue;

        dedxPerHitAtStart.push_back(dedxPerHit.at(hitIndex));
    }

    const auto nHits = dedxPerHitAtStart.size();
    if (nHits == 0)
        return -std::numeric_limits<float>::max();

    // Sort the dEdx so we can find the median
    std::sort(dedxPerHitAtStart.begin(), dedxPerHitAtStart.end());
    const auto median = dedxPerHitAtStart.at(nHits / 2);

    // Now find the mean
    float total = 0.f;
    for (const auto &dEdx : dedxPerHitAtStart)
        total += dEdx;

    const auto mean = total / static_cast<float>(nHits);

    // Now find the variance
    float squareSum = 0.f;
    for (const auto &dEdx : dedxPerHitAtStart)
        squareSum += std::pow(dEdx - mean, 2);

    const auto variance = squareSum / static_cast<float>(nHits);

    // Get the mean dEdx of the hits within one standard deviation of the median
    float truncatedTotal = 0.f;
    unsigned int nTruncatedHits = 0;
    for (const auto &dEdx : dedxPerHitAtStart)
    {
        if (std::pow(dEdx - median, 2) > variance)
            continue;

        truncatedTotal += dEdx;
        nTruncatedHits++;
    }

    if (nTruncatedHits == 0)
        return -std::numeric_limits<float>::max();

    const auto truncatedMean = truncatedTotal / static_cast<float>(nTruncatedHits);

    return truncatedMean;
}

} // namespace ubcc1pi
