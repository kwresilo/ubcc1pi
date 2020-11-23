/**
 *  @file  ubcc1pi/Analyzers/FileLocationDatabase.h
 *
 *  @brief The header file for the file location database analyzer
 */

#ifndef UBCC1PI_ANALYZERS_FILE_LOCATION_DATABASE
#define UBCC1PI_ANALYZERS_FILE_LOCATION_DATABASE

#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDAnalyzer.h"

#include <TFile.h>
#include <TTree.h>

namespace ubcc1pi
{

/**
 *  @brief  The file location database class
 */
class FileLocationDatabase : public art::EDAnalyzer
{
    public:

        /**
         *  @brief  Constructor
         *
         *  @param  pset the input FHiCL parameters
         */
        FileLocationDatabase(const fhicl::ParameterSet &pset);

        /**
         *  @brief  Destructor
         */
        ~FileLocationDatabase();

        /**
         *  @brief  Analyze an event
         *
         *  @param  event the event to analyze
         */
        void analyze(const art::Event &event);

        /**
         *  @brief  Respond to the opening of an input file
         *
         *  @param  fileBlock the input file information
         */
        void respondToOpenInputFile(const art::FileBlock &fileBlock);

    private:

        int          m_run;             ///< The current run number
        int          m_subRun;          ///< The current subrun number
        int          m_event;           ///< The current event number
        int          m_skip;            ///< The number of events to skip in the input file to get to the current event
        std::string  m_currentFileName; ///< The name of the current input file

        std::string  m_outputFileName;  ///< The name of the output file
        TFile       *m_pFile;           ///< The output file
        TTree       *m_pTree;           ///< The output tree
};

} // namespace ubcc1pi

DEFINE_ART_MODULE(ubcc1pi::FileLocationDatabase)

#endif
