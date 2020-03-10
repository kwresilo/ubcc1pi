/**
 *  @file  ubcc1pi/Analyzers/Test_module.cc
 *
 *  @brief The implementation of the test module
 */

#include "ubcc1pi/Analyzers/Test.h"

namespace ubcc1pi
{

Test::Test(const art::EDAnalyzer::Table<EventFactory::Config> &config) :
    art::EDAnalyzer(config),
    m_config(config()),
    m_writer("testOutput.root"),
    m_pEvent(m_writer.GetBoundEventAddress())
{
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void Test::analyze(const art::Event &event)
{
    EventFactory::PopulateEvent(event, m_config, m_pEvent);
    m_pEvent->Print();
    m_writer.FillEvent();
}

} // namespace ubcc1pi
