/**
 *  @file  ubcc1pi/Analyzers/Test_module.cc
 *
 *  @brief The implementation of the test module
 */

#include "ubcc1pi/Analyzers/Test.h"

#include "ubcc1pi/Objects/EventFactory.h"

namespace ubcc1pi
{

Test::Test(const fhicl::ParameterSet &pset) : 
    art::EDAnalyzer(pset),
    m_writer("testOutput.root"),
    m_pEvent(m_writer.GetBoundEventAddress())
{
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void Test::analyze(const art::Event &event)
{
    EventFactory::PopulateEvent(event, m_pEvent);
    m_writer.FillEvent();
}

} // namespace ubcc1pi
