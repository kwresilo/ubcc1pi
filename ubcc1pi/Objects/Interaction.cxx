/**
 *  @file  ubcc1pi/Objects/Interaction.cxx
 *
 *  @brief Implementation of the interaction class
 */

#include "ubcc1pi/Objects/Interaction.h"

namespace ubcc1pi
{

Interaction::Interaction(const art::Event &event, const art::InputTag &mcTruthLabel, const art::InputTag &mcParticleLabel)
{
    const auto nuMCTruth = TruthHelper::GetNeutrinoMCTruth(event, mcTruthLabel);

    if (!nuMCTruth->NeutrinoSet())
        throw cet::exception("Interaction::Interaction") << " - Beam neutrino MCTruth block doesn't have it's neutrino information filled." << std::endl;

    const auto nu = nuMCTruth->GetNeutrino();
    m_neutrino = nu.Nu();
    m_ccnc = static_cast<simb::curr_type_>(nu.CCNC());
    m_mode = static_cast<simb::int_type_>(nu.Mode());
    
    m_allMCParticles = TruthHelper::GetMCParticlesFromMCTruth(event, mcTruthLabel, mcParticleLabel, nuMCTruth);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

simb::MCParticle Interaction::GetNeutrino() const
{
    return m_neutrino;
}

// -----------------------------------------------------------------------------------------------------------------------------------------
        
simb::curr_type_ Interaction::GetCCNC() const
{
    return m_ccnc;
}

// -----------------------------------------------------------------------------------------------------------------------------------------
        
simb::int_type_ Interaction::GetIteractionMode() const
{
    return m_mode;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

MCParticleVector Interaction::GetAllMCParticles() const
{
    return m_allMCParticles;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void Interaction::PrintInfo() const
{
    // Convert the interaction type to a human readable string
    std::string nuType;
    switch (m_neutrino.PdgCode())
    {
        case 12:
            nuType = "Electron Neutrino";
            break;
        case 14:
            nuType = "Muon Neutrino";
            break;
        case 16:
            nuType = "Tau Neutrino";
            break;
        case -12:
            nuType = "Anti-Electron Neutrino";
            break;
        case -14:
            nuType = "Anti-Muon Neutrino";
            break;
        case -16:
            nuType = "Anti-Tau Neutrino";
            break;
        default:
            throw cet::exception("Interaction::PrintInfo") << " - MC Neutrino has non-neutrino PDG code: " << m_neutrino.PdgCode() << std::endl;
    }

    std::string ccnc;
    switch (m_ccnc)
    {
        case simb::kCC:
            ccnc = "Charged-Current";
            break;
        case simb::kNC:
            ccnc = "Neutral-Current";
            break;
        default:
            throw cet::exception("Interaction::PrintInfo") << " - Interaction is not CC or NC!" << std::endl;
    }
    
    std::string mode;
    switch (m_mode)
    {
        case simb::kUnknownInteraction:
            mode = "(Unknown Interaction Type)";
            break;
        case simb::kQE:
            mode = "Quasi-Elastic";
            break;
        case simb::kRes:
            mode = "Resonant";
            break;
        case simb::kDIS:
            mode = "Deep Inelastic Scattering";
            break;
        case simb::kCoh:
            mode = "Coherent";
            break;
        case simb::kCohElastic:
            mode = "Coherent Elastic";
            break;
        case simb::kElectronScattering:
            mode = "Electron Scattering";
            break;
        case simb::kIMDAnnihilation:
            mode = "IMD Annihilation";
            break;
        case simb::kInverseBetaDecay:
            mode = "Inverse Beta Decay";
            break;
        case simb::kGlashowResonance:
            mode = "Glashow Resonance";
            break;
        case simb::kAMNuGamma:
            mode = "AM Nu Gamma";
            break;
        case simb::kMEC:
            mode = "Meson Exchange Current";
            break;
        case simb::kDiffractive:
            mode = "Diffractive";
            break;
        case simb::kEM:
            mode = "Electromagnetic";
            break;
        case simb::kWeakMix:
            mode = "Weak Mix";
            break;
        default:
            mode = "Other";
    }


    // Output the information about the neutrino interaction
    std::cout << "Neutrino interaction --------------------------------------------------------------------------------------" << std::endl;
    std::cout << "  - Neutrino energy       : " << m_neutrino.E() << std::endl;
    std::cout << "  - Interaction type      : " << nuType << " - " << ccnc << " - " << mode << std::endl;
    
    const auto finalStates = TruthHelper::GetPrimaryMCParticles(m_allMCParticles);

    // Count the final state particles
    std::map<int, int> pdgCountMap;
    std::vector<int> usedPdgs;
    for (const auto &finalState : finalStates)
    {
        const auto pdg = finalState->PdgCode();

        if (pdgCountMap.find(pdg) == pdgCountMap.end())
        {
            pdgCountMap.emplace(pdg, 1);
            usedPdgs.push_back(pdg);
        }
        else
        {
            pdgCountMap.at(pdg)++;
        }
    }
    
    std::cout << "  - Final state particles : ";
    for (const auto pdg : usedPdgs)
    {
        std::cout << " " << pdgCountMap.at(pdg) << " (" << pdg << ")  ";
    }

    std::cout << std::endl;
        
    // Output the final state particles
    for (const auto &finalState : finalStates)
    {
        std::cout << "    - Particle" << std::endl;
        std::cout << "      - PDG Code : " << finalState->PdgCode() << std::endl;
        std::cout << "      - Momentum : " << finalState->P() << std::endl;
    }

    std::cout << "-----------------------------------------------------------------------------------------------------------" << std::endl;
}

} // namespace ubcc1pi
