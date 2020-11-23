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
    m_pSubrunTree(static_cast<TTree*>(m_pFile->Get("subruns"))),
    m_pEvent(std::make_shared<Event>()),
    m_pSubrun(std::make_shared<Subrun>())
{
    this->BindEventToTree();
    this->BindSubrunToTree();

    if (this->GetNumberOfEvents() == 0)
        std::cout << "ubcc1pi::FileReader: Warning, input file has no entries" << std::endl;

    if (this->GetNumberOfSubruns() == 0)
        std::cout << "ubcc1pi::FileReader: Warning, input file has no subruns" << std::endl;

    // By default, disable the branches that hold the event weights as these are large and only needed in a few cases
    this->DisableSystematicBranches();
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

void FileReader::BindSubrunToTree()
{
    m_pSubrun->BindToInputTree(m_pSubrunTree);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

std::shared_ptr<Event> FileReader::GetBoundEventAddress()
{
    return m_pEvent;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

std::shared_ptr<Subrun> FileReader::GetBoundSubrunAddress()
{
    return m_pSubrun;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void FileReader::DisableSystematicBranches()
{
    m_pEventTree->SetBranchStatus("truth_systParam*", false);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void FileReader::EnableSystematicBranches()
{
    m_pEventTree->SetBranchStatus("truth_systParam*", true);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

unsigned int FileReader::GetNumberOfEvents() const
{
    return m_pEventTree->GetEntries();
}

// -----------------------------------------------------------------------------------------------------------------------------------------

unsigned int FileReader::GetNumberOfSubruns() const
{
    return m_pSubrunTree->GetEntries();
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

void FileReader::LoadSubrun(const unsigned int subrunIndex)
{
    const auto nSubruns = this->GetNumberOfSubruns();
    if (subrunIndex >= nSubruns)
        throw std::range_error("Can't load subrun, input subrunIndex is out of bounds");

    m_pSubrun->Reset();
    m_pSubrunTree->GetEntry(subrunIndex);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

} // namespace ubcc1pi
