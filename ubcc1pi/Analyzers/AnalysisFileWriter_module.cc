/**
 *  @file  ubcc1pi/Analyzers/AnalysisFileWriter_module.cc
 *
 *  @brief The implementation of the test module
 */

#include "ubcc1pi/Analyzers/AnalysisFileWriter.h"
#include <ctime> // DEBUG

namespace ubcc1pi
{

AnalysisFileWriter::AnalysisFileWriter(const art::EDAnalyzer::Table<Config> &config) :
    art::EDAnalyzer(config),
    m_eventConfig(config().EventFactoryConfig()),
    m_subrunConfig(config().SubrunFactoryConfig()),
    m_writer(config().OutputFileName()),
    m_pEvent(m_writer.GetBoundEventAddress()),
    m_pSubrun(m_writer.GetBoundSubrunAddress())
{
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void AnalysisFileWriter::analyze(const art::Event &event)
{
    std::cout << "DEBUG - Populating event - Time: " << std::time(nullptr) << std::endl;
    EventFactory::PopulateEvent(event, m_eventConfig, m_pEvent);

    //std::cout << "DEBUG - Printing event - Time: " << std::time(nullptr) << std::endl;
    //m_pEvent->Print();

    std::cout << "DEBUG - Filling event before - Time: " << std::time(nullptr) << std::endl;
    m_writer.FillEvent();
    std::cout << "DEBUG - Filling event after - Time: " << std::time(nullptr) << std::endl;
}

// -----------------------------------------------------------------------------------------------------------------------------------------


void AnalysisFileWriter::endSubRun(const art::SubRun &subrun)
{
    std::cout << "DEBUG - Populating subrun - Time: " << std::time(nullptr) << std::endl;
    SubrunFactory::PopulateSubrun(subrun, m_subrunConfig, m_pSubrun);

    //std::cout << "DEBUG - Printing subrun - Time: " << std::time(nullptr) << std::endl;
    //m_pSubrun->Print();

    std::cout << "DEBUG - Filling subrun - Time: " << std::time(nullptr) << std::endl;
    m_writer.FillSubrun();
}

} // namespace ubcc1pi
