/**
 *  @file  ubcc1pi_standalone/Objects/FileReader.cxx
 *
 *  @brief The implementation of the file reader class
 */

#include "ubcc1pi_standalone/Objects/FileReader.h"

namespace ubcc1pi
{

template <class T, class U>
FileReader<T, U>::FileReader(const std::string &inputFile)
{
    m_inputFile = inputFile;
    m_pFile = TFile::Open(m_inputFile.c_str());
    std::string eventTreeName = "events";
    if constexpr (std::is_same_v<T, EventPeLEE>) eventTreeName = "nuselection/NeutrinoSelectionFilter";
    m_pEventTree = static_cast<TTree*>(m_pFile->Get(eventTreeName.c_str()));
    std::string subrunTreeName = "subruns";
    if constexpr (std::is_same_v<U, SubrunPeLEE>) subrunTreeName = "nuselection/SubRun";
    m_pSubrunTree = static_cast<TTree*>(m_pFile->Get(subrunTreeName.c_str()));
    m_pEvent = std::make_shared<T>();
    m_pSubrun = std::make_shared<U>();


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

template <class T, class U>
FileReader<T, U>::~FileReader()
{
    // ATTN this function does the memory management for the pointers belonging to this class
    m_pFile->Close();
}

// -----------------------------------------------------------------------------------------------------------------------------------------

template <class T, class U>
void FileReader<T, U>::BindEventToTree()
{
    m_pEvent->BindToInputTree(m_pEventTree);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

template <class T, class U>
void FileReader<T, U>::BindSubrunToTree()
{
    m_pSubrun->BindToInputTree(m_pSubrunTree);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

template <class T, class U>
std::shared_ptr<T> FileReader<T, U>::GetBoundEventAddress()
{
    return m_pEvent;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

template <class T, class U>
std::shared_ptr<U> FileReader<T, U>::GetBoundSubrunAddress()
{
    return m_pSubrun;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

template <class T, class U>
void FileReader<T, U>::DisableSystematicBranches()
{
    std::string weightBranches = "truth_systParam*";
    if constexpr (std::is_same_v<T, EventPeLEE>) weightBranches = "weights*";
    m_pEventTree->SetBranchStatus(weightBranches.c_str(), false);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

template <class T, class U>
void FileReader<T, U>::EnableSystematicBranches()
{
    std::string weightBranches = "truth_systParam*";
    if constexpr (std::is_same_v<T, EventPeLEE>) weightBranches = "weights*";
    m_pEventTree->SetBranchStatus(weightBranches.c_str(), true);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

template <class T, class U>
unsigned int FileReader<T, U>::GetNumberOfEvents() const
{
    return m_pEventTree->GetEntries();
}

// -----------------------------------------------------------------------------------------------------------------------------------------

template <class T, class U>
unsigned int FileReader<T, U>::GetNumberOfSubruns() const
{
    return m_pSubrunTree->GetEntries();
}

// -----------------------------------------------------------------------------------------------------------------------------------------

template <class T, class U>
void FileReader<T, U>::LoadEvent(const unsigned int eventIndex)
{
    const auto nEvents = this->GetNumberOfEvents();
    if (eventIndex >= nEvents)
        throw std::range_error("Can't load event, input eventIndex is out of bounds");


    m_pEvent->Reset();
    m_pEventTree->GetEntry(eventIndex);
    m_pEvent->PrepareAfterTreeRead();
}

// -----------------------------------------------------------------------------------------------------------------------------------------

template <class T, class U>
void FileReader<T, U>::LoadSubrun(const unsigned int subrunIndex)
{
    const auto nSubruns = this->GetNumberOfSubruns();
    if (subrunIndex >= nSubruns)
        throw std::range_error("Can't load subrun, input subrunIndex is out of bounds");

    m_pSubrun->Reset();
    m_pSubrunTree->GetEntry(subrunIndex);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

// Explicit instantiations
template class FileReader<Event, Subrun>;
template class FileReader<EventPeLEE, SubrunPeLEE>;

} // namespace ubcc1pi
