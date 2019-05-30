/**
 *  @file  ubcc1pi/Helpers/AnalysisHelper.cxx
 *
 *  @brief Implementation of the analysis helper class
 */

#include "ubcc1pi/Helpers/AnalysisHelper.h"

#include <vector>
#include <map>

namespace ubcc1pi
{

bool AnalysisHelper::IsFiducial(const TVector3 &point)
{
    // For readability
    const auto x = point.x();
    const auto y = point.y();
    const auto z = point.z();
   
    // Detector geometry
    const float xMin = 0.f;
    const float yMin = -116.5f;

    const float xMax = 256.35f;
    const float yMax = 116.5f;

    // Fiducial volume cuts
    const float xMargin = 12.f;
    const float yMargin = 35.f;
    const float zRange1Min = 25.f;
    const float zRange1Max = 675.f;
    const float zRange2Min = 775.f;
    const float zRange2Max = 951.8f;

    std::vector<bool> cuts = {
        (x > xMin + xMargin),
        (x < xMax - xMargin),

        (y > yMin + yMargin),
        (y < yMax - yMargin),

        (z > zRange1Min && z < zRange1Max) || (z > zRange2Min && z < zRange2Max)
    };

    return std::all_of(cuts.begin(), cuts.end(), [](const bool b){return b;});
}

// -----------------------------------------------------------------------------------------------------------------------------------------

bool AnalysisHelper::PassesMomentumThreshold(const art::Ptr<simb::MCParticle> &particle)
{
    // The mapping from PDG code to momentum threshold (GeV)
    std::map<int, float> pdgToThresholdMap;
    pdgToThresholdMap.emplace(2212, 0.3);
    pdgToThresholdMap.emplace(2112, std::numeric_limits<float>::max());

    // If we haven't specified this PDG code then assume the threshold is zero - the particle passes by default
    const auto iter = pdgToThresholdMap.find(particle->PdgCode());
    if (iter == pdgToThresholdMap.end())
        return true;

    // Otherwise insist that the particle has sufficient momentum to pass the threshold
    return (particle->P() >= iter->second);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

bool AnalysisHelper::IsVisible(const art::Ptr<simb::MCParticle> &particle)
{
    // Use the absolute PDG to include antiparticles
    const auto absPDG = std::abs(particle->PdgCode());

    return (absPDG == 11   || // Electron
            absPDG == 13   || // Muon
            absPDG == 2212 || // Proton
            absPDG == 2112 || // Neutron
            absPDG == 22   || // Photon
            absPDG == 211  || // Charged pion
            absPDG == 321  || // Charged kaon
            absPDG == 3112 || // Sigma minus
            absPDG == 3222 || // Sigma plus
            absPDG == 3312 ); // Hyperon minus
}

// -----------------------------------------------------------------------------------------------------------------------------------------

MCParticleVector AnalysisHelper::GetVisibleFinalStates(const MCParticleVector &allMCParticles)
{
    MCParticleVector finalStates;
    
    const auto mcParticleMap = TruthHelper::GetMCParticleMap(allMCParticles);
    const auto primaries = TruthHelper::GetPrimaryMCParticles(allMCParticles);

    // Iteratively collect up the visible final states
    for (const auto &primary : primaries)
        AnalysisHelper::GetVisibleFinalStates(primary, mcParticleMap, finalStates);

    return finalStates;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void AnalysisHelper::GetVisibleFinalStates(const art::Ptr<simb::MCParticle> &particle, const MCParticleMap &mcParticleMap, MCParticleVector &finalStates)
{
    // If this particle is visible, then use it!
    if (AnalysisHelper::IsVisible(particle))
    {
        finalStates.push_back(particle);
        return;
    }

    // Otherwise look to the daughters
    for (const auto &daughter : TruthHelper::GetDaughters(particle, mcParticleMap))
        AnalysisHelper::GetVisibleFinalStates(daughter, mcParticleMap, finalStates);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

MCParticleVector AnalysisHelper::GetPrimaryParticles(const TruthHelper::Interaction &interaction)
{
    return TruthHelper::GetPrimaryMCParticles(interaction.GetAllMCParticles());
}

// -----------------------------------------------------------------------------------------------------------------------------------------

MCParticleVector AnalysisHelper::GetVisibleFinalStates(const TruthHelper::Interaction &interaction)
{
    return AnalysisHelper::GetVisibleFinalStates(interaction.GetAllMCParticles());
}

// -----------------------------------------------------------------------------------------------------------------------------------------

MCParticleVector AnalysisHelper::GetReconstructableFinalStates(const TruthHelper::Interaction &interaction)
{
    MCParticleVector reconstructableParticles;

    for (const auto &finalState : AnalysisHelper::GetVisibleFinalStates(interaction))
    {
        if (AnalysisHelper::PassesMomentumThreshold(finalState))
            reconstructableParticles.push_back(finalState);
    }

    return reconstructableParticles;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

bool AnalysisHelper::IsNeutrinoVertexFiducial(const TruthHelper::Interaction &interaction)
{
    return AnalysisHelper::IsFiducial(interaction.GetNeutrino().Position().Vect());
}

// -----------------------------------------------------------------------------------------------------------------------------------------

bool AnalysisHelper::IsCC1PiSignal(const TruthHelper::Interaction &interaction)
{
    // The neutrino vertex must be fiducial
    if (!AnalysisHelper::IsNeutrinoVertexFiducial(interaction))
        return false;

    // Insist that we have exactly the correct number of particles passing the thresholds of each type
    unsigned int nMu = 0;
    unsigned int nPiPlus = 0;

    const auto particles = AnalysisHelper::GetReconstructableFinalStates(interaction);
    for (const auto &particle : particles)
    {
        switch (particle->PdgCode())
        {
            // Mu-
            case 13:
                nMu++;
                break;
            // Pi+
            case 211:
                nPiPlus++;
                break;
            // Protons
            case 2212:
                break;
            // All other particles
            default:
                return false;
        }
    }

    return (nMu == 1 && nPiPlus == 1);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

unsigned int AnalysisHelper::CountParticlesWithPDG(const MCParticleVector &particles, const int pdg)
{
    unsigned int count = 0;

    for (const auto &particle : particles)
    {
        if (particle->PdgCode() == pdg)
            count++;
    }

    return count;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

BacktrackHelper::BacktrackerData AnalysisHelper::GetBacktrackerData(const art::Event &event, const art::InputTag &mcTruthLabel, const art::InputTag &mcParticleLabel, const art::InputTag &backtrackerLabel, const art::InputTag &pfParticleLabel)
{
    const TruthHelper::Interaction interaction(event, mcTruthLabel, mcParticleLabel);
    const auto allPFParticles = CollectionHelper::GetCollection<recob::PFParticle>(event, pfParticleLabel);

    const auto finalStatePFParticles = RecoHelper::GetNeutrinoFinalStates(allPFParticles);
    const auto finalStateMCParticles = AnalysisHelper::GetReconstructableFinalStates(interaction);

    const auto hitsToPfps = BacktrackHelper::GetHitToPFParticleMap(event, pfParticleLabel, finalStatePFParticles);
    const auto hitsToMcps = BacktrackHelper::GetHitToMCParticleWeightMap(event, mcParticleLabel, backtrackerLabel, finalStateMCParticles);

    return BacktrackHelper::BacktrackerData(finalStatePFParticles, finalStateMCParticles, hitsToPfps, hitsToMcps);
}

// -----------------------------------------------------------------------------------------------------------------------------------------
        
BacktrackHelper::SliceMetadata AnalysisHelper::GetSliceMetadata(const art::Event &event, const art::InputTag &mcTruthLabel, const art::InputTag &mcParticleLabel, const art::InputTag &backtrackerLabel, const art::InputTag &pfParticleLabel, const art::InputTag &sliceLabel, const art::InputTag &hitLabel)
{
    const auto slices = CollectionHelper::GetCollection<recob::Slice>(event, sliceLabel);
    const auto slicesToHits = CollectionHelper::GetAssociation<recob::Slice, recob::Hit>(event, sliceLabel);
    const auto pfParticlesToSlices = CollectionHelper::GetAssociation<recob::PFParticle, recob::Slice>(event, pfParticleLabel, sliceLabel);
    const auto slicesToPFParticles = CollectionHelper::GetReversedAssociation(pfParticlesToSlices);    
    const auto nuMCTruth = TruthHelper::GetNeutrinoMCTruth(event, mcTruthLabel);
    const auto nuMCParticles = TruthHelper::GetMCParticlesFromMCTruth(event, mcTruthLabel, mcParticleLabel, nuMCTruth);
    const auto hitsToIsNuInduced = BacktrackHelper::GetHitsToIsNuInducedMap(event, hitLabel, mcParticleLabel, backtrackerLabel, nuMCParticles);
    
    SlicesToBool sliceToIsSelectedAsNu;
    for (const auto &slice : slices)
    {
        if (!sliceToIsSelectedAsNu.emplace(slice, RecoHelper::IsSliceSelectedAsNu(slice, slicesToPFParticles)).second)
            throw cet::exception("AnalysisHelper::GetSliceMetadata") << " - Repeated slices." << std::endl;
    }

    return BacktrackHelper::SliceMetadata(slices, sliceToIsSelectedAsNu, slicesToHits, hitsToIsNuInduced);
}

} // namespace ubcc1pi
