/**
 *  @file  ubcc1pi/Objects/FileWriter.cxx
 *
 *  @brief The implementation of the file writer class
 */

#include "ubcc1pi/Objects/FileWriter.h"

namespace ubcc1pi
{

FileWriter::FileWriter(const std::string &outputFile) :
    m_outputFile(outputFile),
    m_pFile(new TFile(m_outputFile.c_str(), "RECREATE")),
    m_pEventTree(new TTree("events", "events")),
    m_pSubrunTree(new TTree("subruns", "subruns")),
    m_pEvent(new Event()),
    m_pSubrun(new Subrun())
{
    this->BindEventToTree();
    this->BindSubrunToTree();
}

// -----------------------------------------------------------------------------------------------------------------------------------------

FileWriter::~FileWriter()
{
    std::cout << "BEGIN DEBUG ----------------------------" << std::endl;
    std::cout << "nEvents : " << m_pEventTree->GetEntries() << std::endl;
    std::cout << "nSubruns : " << m_pSubrunTree->GetEntries() << std::endl;
    std::cout << "END DEBUG ------------------------------" << std::endl;

    m_pFile->Write();

    // ATTN this function does the memory management for the pointers belonging to this class
    m_pFile->Close();

    delete m_pEvent;
    delete m_pSubrun;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void FileWriter::BindEventToTree()
{
    m_pEvent->BindToOutputTree(m_pEventTree);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void FileWriter::BindSubrunToTree()
{
    m_pSubrun->BindToOutputTree(m_pSubrunTree);
}

// -----------------------------------------------------------------------------------------------------------------------------------------
        
Event * FileWriter::GetBoundEventAddress()
{
    return m_pEvent;
}

// -----------------------------------------------------------------------------------------------------------------------------------------
        
Subrun * FileWriter::GetBoundSubrunAddress()
{
    return m_pSubrun;
}

// -----------------------------------------------------------------------------------------------------------------------------------------
        
void FileWriter::FillEvent()
{
    m_pEvent->PrepareForTreeFill();
    m_pEventTree->Fill();
}

// -----------------------------------------------------------------------------------------------------------------------------------------
        
void FileWriter::FillSubrun()
{
    m_pSubrunTree->Fill();
}

} // namespace ubcc1pi
