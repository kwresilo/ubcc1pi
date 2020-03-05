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
    m_pEvent(new Event())
{
    this->BindEventToTree();
}

// -----------------------------------------------------------------------------------------------------------------------------------------

FileWriter::~FileWriter()
{
    m_pFile->Write();

    // ATTN this function does the memory management for the pointers belonging to this class
    m_pFile->Close();

    delete m_pEvent;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void FileWriter::BindEventToTree()
{
    m_pEvent->BindToOutputTree(m_pEventTree);
}

// -----------------------------------------------------------------------------------------------------------------------------------------
        
Event * FileWriter::GetBoundEventAddress()
{
    return m_pEvent;
}

// -----------------------------------------------------------------------------------------------------------------------------------------
        
void FileWriter::FillEvent()
{
    m_pEvent->PrepareForTreeFill();
    m_pEventTree->Fill();
}

} // namespace ubcc1pi
