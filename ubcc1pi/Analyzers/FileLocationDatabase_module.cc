/**
 *  @file  ubcc1pi/Analyzers/FileLocationDatabase_module.cc
 *
 *  @brief The implementation of the file location database module
 */

#include "ubcc1pi/Analyzers/FileLocationDatabase.h"

#include "art/Framework/Core/FileBlock.h"

#include <stdexcept>

namespace ubcc1pi
{

FileLocationDatabase::FileLocationDatabase(const fhicl::ParameterSet &pset) :
    EDAnalyzer(pset),
    m_outputFileName("fileLocations.root"),
    m_pFile(new TFile(m_outputFileName.c_str(), "RECREATE")),
    m_pTree(new TTree("events", "events"))
{
    // Setup the branches
    m_pTree->Branch("run", &m_run);
    m_pTree->Branch("subRun", &m_subRun);
    m_pTree->Branch("event", &m_event);
    m_pTree->Branch("file", &m_currentFileName);
    m_pTree->Branch("skip", &m_skip);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

FileLocationDatabase::~FileLocationDatabase()
{
    m_pFile->Write();
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void FileLocationDatabase::analyze(const art::Event &event)
{
    if (m_currentFileName.empty())
        throw std::logic_error("FileLocationDatabase::analyze - current file name is unknown");

    // Set the current run-subrun-event for this event
    m_run = event.run();
    m_subRun = event.subRun();
    m_event = event.event();

    std::cout << m_run << ":" << m_subRun << ":" << m_event << ":" << m_skip << ":" << m_currentFileName << std::endl;

    // Fill the output tree and increment the skip
    m_pTree->Fill();
    m_skip++;
}

// -----------------------------------------------------------------------------------------------------------------------------------------
        
void FileLocationDatabase::respondToOpenInputFile(const art::FileBlock &fileBlock)
{
    m_currentFileName = fileBlock.fileName();
    m_skip = 0;
}

} // namespace ubcc1pi
