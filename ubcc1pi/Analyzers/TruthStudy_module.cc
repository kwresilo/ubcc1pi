/**
 *  @file  ubcc1pi/Analyzers/TruthStudy_module.cc
 *
 *  @brief The implementation file for the truth study analyzer.
 */

#include "ubcc1pi/Analyzers/TruthStudy.h"

namespace ubcc1pi
{

TruthStudy::TruthStudy(const fhicl::ParameterSet &config) : art::EDAnalyzer(config)
{
}

// -----------------------------------------------------------------------------------------------------------------------------------------
 
void TruthStudy::analyze(const art::Event &event)
{
}

} // namespace ubcc1pi
