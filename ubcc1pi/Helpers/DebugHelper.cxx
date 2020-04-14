/**
 *  @file  ubcc1pi/Helpers/DebugHelper.cxx
 *
 *  @brief Implementation of the debug helper class
 */

#include "ubcc1pi/Helpers/DebugHelper.h"

#include <vector>
#include <unordered_map>

namespace ubcc1pi
{

void DebugHelper::Print(const art::Ptr<simb::MCParticle> &particle, const unsigned int depth)
{
    DebugHelper::PrintHeader("MCParticle", depth);
    DebugHelper::PrintProperty("ID", particle->TrackId(), depth);
    DebugHelper::PrintProperty("PDG", particle->PdgCode(), depth);
    DebugHelper::PrintProperty("Momentum", particle->P(), depth);
    DebugHelper::PrintProperty("Time", particle->T(), depth);
    DebugHelper::PrintProperty("Process", particle->Process(), depth);
    DebugHelper::PrintProperty("End Time", particle->EndT(), depth);
    DebugHelper::PrintProperty("End Process", particle->EndProcess(), depth);
    DebugHelper::PrintProperty("Lifetime", particle->EndT() - particle->T(), depth);
    DebugHelper::PrintProperty("Length", (particle->EndPosition() - particle->Position()).Vect().Mag(), depth);
    DebugHelper::PrintProperty("Depth", depth, depth);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void DebugHelper::Print(const MCParticleVector &particles, const unsigned int depth)
{
    for (const auto &particle : particles)
        DebugHelper::Print(particle, depth);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void DebugHelper::PrintSummary(const MCParticleVector &particles, const unsigned int depth)
{
    // Count the number of particles with each pdg code
    std::unordered_map<int, MCParticleVector> pdgToMCParticleMap;
    std::vector<int> usedPdgs;

    for (const auto &particle : particles)
    {
        const auto pdg = particle->PdgCode();
        pdgToMCParticleMap[pdg].push_back(particle);

        if (std::find(usedPdgs.begin(), usedPdgs.end(), pdg) == usedPdgs.end())
            usedPdgs.push_back(pdg);
    }

    // Print the summary
    std::sort(usedPdgs.begin(), usedPdgs.end());
    for (const auto &pdg : usedPdgs)
    {
        DebugHelper::PrintProperty("# (" + std::to_string(pdg) + ")", pdgToMCParticleMap.at(pdg).size(), depth);
    }
}

// -----------------------------------------------------------------------------------------------------------------------------------------

std::string DebugHelper::GetInteractionString(const TruthHelper::Interaction &interaction, const bool brief)
{
    // Convert the interaction type to a human readable string
    std::string nuType;
    const auto neutrino = interaction.GetNeutrino();
    switch (neutrino.PdgCode())
    {
        case 12:
            nuType = brief ? "nue" : "Electron Neutrino";
            break;
        case 14:
            nuType = brief ? "numu" : "Muon Neutrino";
            break;
        case 16:
            nuType = brief ? "nutau" : "Tau Neutrino";
            break;
        case -12:
            nuType = brief ? "nuebar" : "Anti-Electron Neutrino";
            break;
        case -14:
            nuType = brief ? "numubar" : "Anti-Muon Neutrino";
            break;
        case -16:
            nuType = brief ? "nutaubar" : "Anti-Tau Neutrino";
            break;
        default:
            throw cet::exception("Interaction::PrintInfo") << " - MC Neutrino has non-neutrino PDG code: " << neutrino.PdgCode() << std::endl;
    }

    std::string ccnc;
    switch (interaction.GetCCNC())
    {
        case simb::kCC:
            ccnc = brief ? "CC" : "Charged-Current";
            break;
        case simb::kNC:
            ccnc = brief ? "NC" : "Neutral-Current";
            break;
        default:
            throw cet::exception("Interaction::PrintInfo") << " - Interaction is not CC or NC!" << std::endl;
    }
    
    std::string mode;
    switch (interaction.GetInteractionMode())
    {
        case simb::kUnknownInteraction:
            mode = brief ? "Unknown" : "(Unknown Interaction Type)";
            break;
        case simb::kQE:
            mode = brief ? "QE" : "Quasi-Elastic";
            break;
        case simb::kRes:
            mode = brief ? "Res" : "Resonant";
            break;
        case simb::kDIS:
            mode = brief ? "DIS" : "Deep Inelastic Scattering";
            break;
        case simb::kCoh:
            mode = brief ? "Coh" : "Coherent";
            break;
        case simb::kCohElastic:
            mode = brief ? "CohElastic" : "Coherent Elastic";
            break;
        case simb::kElectronScattering:
            mode = brief ? "ElecScat" : "Electron Scattering";
            break;
        case simb::kIMDAnnihilation:
            mode = brief ? "IMDAnnihil" : "IMD Annihilation";
            break;
        case simb::kInverseBetaDecay:
            mode = brief ? "IBD" : "Inverse Beta Decay";
            break;
        case simb::kGlashowResonance:
            mode = brief ? "Glashow" : "Glashow Resonance";
            break;
        case simb::kAMNuGamma:
            mode = brief ? "AMNuGamma" : "AM Nu Gamma";
            break;
        case simb::kMEC:
            mode = brief ? "MEC" : "Meson Exchange Current";
            break;
        case simb::kDiffractive:
            mode = brief ? "Diff" : "Diffractive";
            break;
        case simb::kEM:
            mode = brief ? "EM" : "Electromagnetic";
            break;
        case simb::kWeakMix:
            mode = brief ? "WeakMix" : "Weak Mix";
            break;
        default:
            mode = "Other";
    }

    if (brief)
        return (nuType + "_" + ccnc + "_" + mode);
        
    return (nuType + ", " + ccnc + ", " + mode);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void DebugHelper::Print(const art::Ptr<recob::PFParticle> &particle, const unsigned int depth)
{
    DebugHelper::PrintHeader("PFParticle", depth);
    DebugHelper::PrintProperty("ID", particle->Self(), depth);
    DebugHelper::PrintProperty("Parent", particle->Parent(), depth);
    DebugHelper::PrintProperty("PDG", particle->PdgCode(), depth);
}

// -----------------------------------------------------------------------------------------------------------------------------------------
        
void DebugHelper::Print(const PFParticleVector &particles, const unsigned int depth)
{
    for (const auto &particle : particles)
        DebugHelper::Print(particle, depth);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void DebugHelper::Print(const BacktrackHelper::BacktrackerData &data, const unsigned int depth)
{
    DebugHelper::PrintHeader("Backtracker data", depth);

    DebugHelper::PrintHeader("PFParticles", depth);
    for (const auto &pfParticle : data.GetPFParticles())
    {
        DebugHelper::Print(pfParticle, depth + 1);
        DebugHelper::PrintProperty("# Hits", data.GetNHits(pfParticle), depth + 1);
        DebugHelper::PrintProperty("Best matched MCParticle", depth + 1);

        try
        {
            const auto mcParticle = data.GetBestMatchedMCParticle(pfParticle);
            DebugHelper::Print(mcParticle, depth + 2);
            DebugHelper::PrintProperty("Purity", data.GetMatchPurity(pfParticle, mcParticle), depth + 2);
            DebugHelper::PrintProperty("Completeness", data.GetMatchCompleteness(pfParticle, mcParticle), depth + 2);
        }
        catch (const cet::exception &)
        {
        } 
    }
    
    DebugHelper::PrintHeader("MCParticles", depth);
    for (const auto &mcParticle : data.GetMCParticles())
    {
        DebugHelper::Print(mcParticle, depth + 1);
        DebugHelper::PrintProperty("# Hits (weighted)", data.GetWeight(mcParticle), depth + 1);

        const auto pfParticles = data.GetBestMatchedPFParticles(mcParticle);
        DebugHelper::PrintProperty("Best matched PFParticles", pfParticles.size(), depth + 1);
            
        for (const auto &pfParticle : pfParticles)
        {
            DebugHelper::Print(pfParticle, depth + 2);
            DebugHelper::PrintProperty("# Hits", data.GetNHits(pfParticle), depth + 2);
            DebugHelper::PrintProperty("Purity", data.GetMatchPurity(pfParticle, mcParticle), depth + 2);
            DebugHelper::PrintProperty("Completeness", data.GetMatchCompleteness(pfParticle, mcParticle), depth + 2);
        }
    }
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void DebugHelper::Print(const BacktrackHelper::SliceMetadata &sliceMetadata, const unsigned int depth)
{
    const auto slices = sliceMetadata.GetSlices();
    const auto totalNuHits = sliceMetadata.GetTotalNumberOfNuInducedHits();
    const auto selectedSlices = sliceMetadata.GetSelectedNeutrinoSlices();

    DebugHelper::PrintHeader("Slice metadata", depth);
    DebugHelper::PrintProperty("# Slices", slices.size(), depth);
    DebugHelper::PrintProperty("# Slices selected", selectedSlices.size(), depth);
    DebugHelper::PrintProperty("# Total neutrino hits", totalNuHits, depth);

    art::Ptr<recob::Slice> mostCompleteSlice;
    if (totalNuHits > 0)
    {
        mostCompleteSlice = sliceMetadata.GetMostCompleteSlice();
        DebugHelper::PrintProperty("Most complete slice selected?", sliceMetadata.IsMostCompleteSliceSelected(), depth);
    }

    for (const auto &slice : slices)
    {
        DebugHelper::PrintHeader("Slice", depth + 1);

        const auto nHits = sliceMetadata.GetNumberOfHits(slice);
        DebugHelper::PrintProperty("# Hits", nHits, depth + 1);

        if (nHits > 0)
            DebugHelper::PrintProperty("Purity", sliceMetadata.GetPurity(slice), depth + 1);

        if (totalNuHits > 0)
        {
            DebugHelper::PrintProperty("Completeness", sliceMetadata.GetCompleteness(slice), depth + 1);
            DebugHelper::PrintProperty("Most complete?", slice == mostCompleteSlice, depth + 1);
        }

        const auto isSelected = (std::find(selectedSlices.begin(), selectedSlices.end(), slice) != selectedSlices.end());
        DebugHelper::PrintProperty("Is selected as nu?", isSelected, depth + 1);
    }
}

} // namespace ubcc1pi
