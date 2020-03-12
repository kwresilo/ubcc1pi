/**
 *  @file  ubcc1pi_standalone/Objects/FileReader.h
 *
 *  @brief The header file for the event file reader class
 */

#ifndef UBCC1PI_STANDALONE_OBJECTS_FILE_READER
#define UBCC1PI_STANDALONE_OBJECTS_FILE_READER

#include "ubcc1pi_standalone/Interface/Event.h"

#include <TFile.h>
#include <TTree.h>

namespace ubcc1pi
{

/**
 *  @brief  The File Reader class
 */
class FileReader
{
    public:

        /**
         *  @brief  Constructor
         *
         *  @param  inputFile the path to the input file
         */
        FileReader(const std::string &inputFile);

        /**
         *  @brief  Destructor
         */
        ~FileReader();
        
        /**
         *  @brief  Get the total number of events 
         *
         *  @return number of events
         */
        unsigned int GetNumberOfEvents() const;

        /**
         *  @brief  Get the event bound to the output tree
         *
         *  @return the address of the bound event
         */
        Event * GetBoundEventAddress();

        /**
         *  @brief  Load the event with the supplied index
         *
         *  @param  eventIndex the event to load
         */
        void LoadEvent(const unsigned int eventIndex);

    private:

        /**
         *  @brief  Bind the event member to the input tree
         */
        void BindEventToTree();

        std::string m_inputFile;  ///< The input file name

        TFile      *m_pFile;       ///< The input file
        TTree      *m_pEventTree;  ///< The input event tree

        Event      *m_pEvent;       ///< The input event
};

} // namespace ubcc1pi

#endif
