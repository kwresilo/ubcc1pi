/**
 *  @file  ubcc1pi/Analyzers/AnalysisFileWriter_module.cc
 *
 *  @brief The implementation of the test module
 */

#include "ubcc1pi/Analyzers/AnalysisFileWriter.h"

namespace ubcc1pi
{

AnalysisFileWriter::AnalysisFileWriter(const art::EDAnalyzer::Table<Config> &config) :
    art::EDAnalyzer(config),
    m_config(config().EventFactoryConfig()),
    m_writer(config().OutputFileName()),
    m_pEvent(m_writer.GetBoundEventAddress())
{
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void AnalysisFileWriter::analyze(const art::Event &event)
{
    EventFactory::PopulateEvent(event, m_config, m_pEvent);
    m_pEvent->Print();
    m_writer.FillEvent();
}

} // namespace ubcc1pi
