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
    m_pEvent(std::make_shared<Event>())
{
    this->BindEventToTree();

    std::cout << "Checking number of events isn't zero" << std::endl;
    if (this->GetNumberOfEvents() == 0)
        std::cout << "ubcc1pi::FileReader: Warning, input file is empty" << std::endl;
    
    std::cout << "    - we goooood" << std::endl;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

FileReader::~FileReader()
{
    // ATTN this function does the memory management for the pointers belonging to this class
    m_pFile->Close();
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void FileReader::BindEventToTree()
{
    m_pEvent->BindToInputTree(m_pEventTree);
}

// -----------------------------------------------------------------------------------------------------------------------------------------
        
std::shared_ptr<Event> FileReader::GetBoundEventAddress()
{
    std::cout << "Getting bound event address" << std::endl;
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

    std::cout << "Resetting event" << std::endl;
    m_pEvent->Reset();
    
    std::cout << "Getting the entry" << std::endl;
    m_pEventTree->GetEntry(eventIndex);
    
    std::cout << "Preparing after read" << std::endl;
    m_pEvent->PrepareAfterTreeRead();
}

// -----------------------------------------------------------------------------------------------------------------------------------------

} // namespace ubcc1pi
