/**
 *  @file  ubcc1pi/Objects/FileWriter.h
 *
 *  @brief The header file for the event file writer class
 */

#ifndef UBCC1PI_OBJECTS_FILE_WRITER
#define UBCC1PI_OBJECTS_FILE_WRITER

#include "ubcc1pi_standalone/interface/Event.h"

#include <TFile.h>
#include <TTree.h>

namespace ubcc1pi
{

/**
 *  @brief  The File Writer class
 *          The file write own an Event which is bound to the output tree. This Event can be accessed and modified by the user.
 *          Once modified the file writer will fill the entries from the event into the output tree
 */
class FileWriter
{
    public:

        /**
         *  @brief  Constructor
         *
         *  @param  outputFile the path to the output file
         */
        FileWriter(const std::string &outputFile);

        /**
         *  @brief  Destructor
         */
        ~FileWriter();

        /**
         *  @brief  Get the event bound to the output tree
         *
         *  @return the address of the bound event
         */
        Event * GetBoundEventAddress();

        /**
         *  @brief  Fill the output tree with the event
         */
        void FillEvent();

    private:

        /**
         *  @brief  Bind the event member to the output tree
         */
        void BindEventToTree();

        std::string m_outputFile;  ///< The output file name

        TFile      *m_pFile;       ///< The output file
        TTree      *m_pEventTree;  ///< The output event tree

        Event      *m_pEvent;       ///< The output event
};

} // namespace ubcc1pi

#endif
