#include "ubcc1pi_standalone/Helpers/AnalysisHelper.h"

#include "ubcc1pi_standalone/Helpers/GeometryHelper.h"

#include <stdexcept>

namespace ubcc1pi
{

bool AnalysisHelper::IsTrueCC1Pi(const std::shared_ptr<Event> &pEvent)
{
    if (!pEvent->metadata.hasTruthInfo())
        throw std::invalid_argument("AnalysisHelper::IsTrueCC1Pi - Input event doesn't have truth information!");

    const auto truth = pEvent->truth;

    // Insist the true neutrino is fiducial
    if (!AnalysisHelper::IsFiducial(truth.nuVertex()))
        return false;
    
    // Check the topology
    unsigned int nMu = 0;
    unsigned int nProton = 0;
    unsigned int nPion = 0;
    unsigned int nOther = 0;

    for (const auto &particle : pEvent->truth.particles)
    {
        if (!AnalysisHelper::PassesVisibilityThreshold(particle))
            continue;

        switch (particle.pdgCode())
        {
            case 13:
                nMu++;
                break;
            case 2212:
                nProton++;
                break;
            case 211:
                nPion++;
                break;
            default:
                nOther++;
                break;
        }
    }

    // Insist the CC1Pi topology
    return (nMu == 1 && nPion == 1 && nOther == 0);
}

// -----------------------------------------------------------------------------------------------------------------------------------------
    
bool AnalysisHelper::PassesVisibilityThreshold(const Event::Truth::Particle &particle)
{
    // First only consider specific particle types that can lead to hits
    const auto absPDG = std::abs(particle.pdgCode());
    const auto isVisible = (absPDG == 11   || // Electron
                            absPDG == 13   || // Muon
                            absPDG == 2212 || // Proton
                            absPDG == 2112 || // Neutron
                            absPDG == 22   || // Photon
                            absPDG == 211  || // Charged pion
                            absPDG == 321  || // Charged kaon
                            absPDG == 3112 || // Sigma minus
                            absPDG == 3222 || // Sigma plus
                            absPDG == 3312 ); // Hyperon minus

    if (!isVisible)
        return false;

    // Apply the momentum threhsolds
    float thresholdMomentum = -std::numeric_limits<float>::max();
    switch (absPDG)
    {
        case 2112: // Neutron
            thresholdMomentum = std::numeric_limits<float>::max(); // ATTN infinite threshold = never let a neutron pass
            break;
        default: break;
    }

    return (particle.momentum() >= thresholdMomentum);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

bool AnalysisHelper::IsPointWithinMargins(const TVector3 &point, const float lowXMargin, const float highXMargin, const float lowYMargin, const float highYMargin, const float lowZMargin, const float highZMargin)
{
    return ((point.x() > GeometryHelper::lowX + lowXMargin)   && 
            (point.x() < GeometryHelper::highX - highXMargin) &&
            (point.y() > GeometryHelper::lowY + lowYMargin)   &&
            (point.y() < GeometryHelper::highY - highYMargin) &&
            (point.z() > GeometryHelper::lowZ + lowZMargin)   &&
            (point.z() < GeometryHelper::highZ - highZMargin) );
}
        
// -----------------------------------------------------------------------------------------------------------------------------------------

bool AnalysisHelper::IsFiducial(const TVector3 &point)
{
    return AnalysisHelper::IsPointWithinMargins(point, 10.f, 10.f, 10.f, 10.f, 10.f, 50.f);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

bool AnalysisHelper::IsContained(const TVector3 &point)
{
    return AnalysisHelper::IsPointWithinMargins(point, 5.f, 5.f, 5.f, 5.f, 5.f, 5.f);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

bool AnalysisHelper::IsContained(const Event::Truth::Particle &particle)
{
    const auto start = TVector3(particle.startX(), particle.startY(), particle.startZ());
    const auto end = TVector3(particle.endX(), particle.endY(), particle.endZ());

    return (AnalysisHelper::IsContained(start) && AnalysisHelper::IsContained(end));
}

// -----------------------------------------------------------------------------------------------------------------------------------------

bool AnalysisHelper::IsContained(const Event::Reco::Particle &particle)
{
    if (!particle.startX.IsSet() || !particle.startY.IsSet() || !particle.startZ.IsSet() || !particle.endX.IsSet() || !particle.endY.IsSet() || !particle.endZ.IsSet())
        throw std::invalid_argument("AnalysisHelper::IsContained - input reco particle doesn't have a fitted track start-end points");

    const auto start = TVector3(particle.startX(), particle.startY(), particle.startZ());
    const auto end = TVector3(particle.endX(), particle.endY(), particle.endZ());

    return (AnalysisHelper::IsContained(start) && AnalysisHelper::IsContained(end));
}

// -----------------------------------------------------------------------------------------------------------------------------------------

unsigned int AnalysisHelper::GetBestMatchedTruthParticleIndex(const Event::Reco::Particle &recoParticle, const std::vector<Event::Truth::Particle> &truthParticles, const bool applyVisibilityThreshold)
{
    if (!recoParticle.hasMatchedMCParticle.IsSet())
        throw std::invalid_argument("AnalysisHelper::GetBestMatchedTruthParticleIndex - input reco particle is from an event without truth info");

    if (!recoParticle.hasMatchedMCParticle())
        throw std::invalid_argument("AnalysisHelper::GetBestMatchedTruthParticleIndex - input reco particle has no matched truth particle");

    if (truthParticles.size() != recoParticle.truthMatchPurities().size())
        throw std::invalid_argument("AnalysisHelper::GetBestMatchedTruthParticleIndex - wrong number of input truth particles");

    if (truthParticles.empty())
        throw std::invalid_argument("AnalysisHelper::GetBestMatchedTruthParticleIndex - no truth particles supplied but reco particle has a match!");


    unsigned int bestMatchedParticleId = std::numeric_limits<unsigned int>::max();
    float bestMatchScore = -std::numeric_limits<float>::max();
    bool foundMatch = false;

    for (unsigned int i = 0; i < truthParticles.size(); ++i)
    {
        if (!AnalysisHelper::PassesVisibilityThreshold(truthParticles.at(i)))
            continue;

        const auto score = recoParticle.truthMatchPurities().at(i);
        if (score < bestMatchScore)
            continue;

        bestMatchScore = score;
        bestMatchedParticleId = i;
        foundMatch = true;
    }

    if (!foundMatch)
        throw std::logic_error("AnalysisHelper::GetBestMatchedTruthParticleIndex - input reco particle has no matched truth particle");
    
    return bestMatchedParticleId;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

Event::Truth::Particle AnalysisHelper::GetBestMatchedTruthParticle(const Event::Reco::Particle &recoParticle, const std::vector<Event::Truth::Particle> &truthParticles, const bool applyVisibilityThreshold)
{
    return truthParticles.at(AnalysisHelper::GetBestMatchedTruthParticleIndex(recoParticle, truthParticles, applyVisibilityThreshold));
}

// -----------------------------------------------------------------------------------------------------------------------------------------
            
bool AnalysisHelper::IsGolden(const Event::Truth::Particle &particle)
{
    return (particle.nElasticScatters() == 0 &&
            particle.nInelasticScatters() == 0 &&
            particle.isStopping() &&
            AnalysisHelper::IsContained(particle));
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void AnalysisHelper::PrintLoadingBar(const unsigned int numerator, const unsigned int denominator)
{
    if (denominator == 0)
        return;

    const auto width = 85u;
    const auto loadedWidth = static_cast<unsigned int>(std::floor(width * (static_cast<float>(numerator) / static_cast<float>(denominator))));
    const auto unloadedWidth = static_cast<unsigned int>(width - loadedWidth);
    
    // Work out if it's worth printing 
    const bool shouldPrint = (numerator == 0) || (static_cast<unsigned int>(std::floor(width * (static_cast<float>(numerator - 1) / static_cast<float>(denominator)))) != loadedWidth);

    if (!shouldPrint)
        return;

    if (numerator != 0)
        std::cout << "\033[F\033[F";

    std::cout << numerator << " / " << denominator << std::endl;
    std::cout << "|" << std::string(loadedWidth, '=') << std::string(unloadedWidth, ' ') << "|" << std::endl;
}

} // namespace ubcc1pi
