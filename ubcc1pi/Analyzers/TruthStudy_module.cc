/**
 *  @file  ubcc1pi/Analyzers/TruthStudy_module.cc
 *
 *  @brief The implementation file for the truth study analyzer.
 */

#include "ubcc1pi/Analyzers/TruthStudy.h"

#include "ubcc1pi/Helpers/TruthHelper.h"

namespace ubcc1pi
{

TruthStudy::TruthStudy(const art::EDAnalyzer::Table<Config> &config) :
    art::EDAnalyzer(config),
    m_config(config)
{
}

// -----------------------------------------------------------------------------------------------------------------------------------------
 
void TruthStudy::analyze(const art::Event &event)
{
    const auto mcTruthLabel = m_config().MCTruthLabel();
    const auto mcParticleLabel = m_config().MCParticleLabel();

    const TruthHelper::Interaction interaction(event, mcTruthLabel, mcParticleLabel);
    //interaction.PrintInfo();
}

} // namespace ubcc1pi
