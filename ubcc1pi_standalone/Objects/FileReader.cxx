/**
 *  @file  ubcc1pi_standalone/Objects/FileReader.cxx
 *
 *  @brief The implementation of the file reader class
 */

#include "ubcc1pi_standalone/Objects/FileReader.h"

namespace ubcc1pi
{

FileReader::FileReader(const std::string &inputFile) :
    m_inputFile(inputFile),
    m_pFile(TFile::Open(m_inputFile.c_str())),
    m_pEventTree(static_cast<TTree*>(m_pFile->Get("events"))),
    m_pEvent(new Event())
{
    this->BindEventToTree();

    if (this->GetNumberOfEvents() == 0)
        std::cout << "ubcc1pi::FileReader: Warning, input file is empty" << std::endl;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

FileReader::~FileReader()
{
    // ATTN this function does the memory management for the pointers belonging to this class
    m_pFile->Close();

    delete m_pEvent;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void FileReader::BindEventToTree()
{
    m_pEvent->BindToInputTree(m_pEventTree);
}

// -----------------------------------------------------------------------------------------------------------------------------------------
        
Event * FileReader::GetBoundEventAddress()
{
    return m_pEvent;
}

// -----------------------------------------------------------------------------------------------------------------------------------------
        
unsigned int FileReader::GetNumberOfEvents() const
{
    return m_pEventTree->GetEntries();
}

// -----------------------------------------------------------------------------------------------------------------------------------------
        
void FileReader::LoadEvent(const unsigned int eventIndex)
{
    const auto nEvents = this->GetNumberOfEvents();
    if (eventIndex >= nEvents)
        throw std::range_error("Can't load event, input eventIndex is out of bounds");

    m_pEvent->Reset();
    m_pEventTree->GetEntry(eventIndex);
    m_pEvent->PrepareAfterTreeRead();
}

// -----------------------------------------------------------------------------------------------------------------------------------------

} // namespace ubcc1pi
