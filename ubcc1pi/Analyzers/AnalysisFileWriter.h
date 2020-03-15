/**
 *  @file  ubcc1pi/Analyzers/AnalysisFileWriter.h
 *
 *  @brief The header file for the analysis file writer
 */

#ifndef UBCC1PI_ANALYZERS_ANALYSIS_FILE_WRITER
#define UBCC1PI_ANALYZERS_ANALYSIS_FILE_WRITER

#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDAnalyzer.h"

#include "art/Framework/Services/Optional/TFileService.h"
#include <TTree.h>

#include "ubcc1pi/Objects/FileWriter.h"
#include "ubcc1pi/Objects/EventFactory.h"
#include "ubcc1pi/Objects/SubrunFactory.h"

namespace ubcc1pi
{

/**
 *  @brief  The analysis file writer class
 */
class AnalysisFileWriter : public art::EDAnalyzer
{
    public:

        /**
         *  @brief  The FHiCL configuration
         */
        struct Config
        {
            fhicl::Atom<std::string> OutputFileName
            {
                fhicl::Name("OutputFileName"),
                fhicl::Comment("The name of the output file")
            };

            fhicl::Table<EventFactory::Config> EventFactoryConfig
            {
                fhicl::Name("EventFactoryConfig"),
                fhicl::Comment("The configuration of the event factory")
            };
            
            fhicl::Table<SubrunFactory::Config> SubrunFactoryConfig
            {
                fhicl::Name("SubrunFactoryConfig"),
                fhicl::Comment("The configuration of the subrun factory")
            };
        };

        /**
         *  @brief  Constructor
         *
         *  @param  config the FHiCl configuration table
         */
        AnalysisFileWriter(const art::EDAnalyzer::Table<Config> &config);

        /**
         *  @brief  Analyze an event
         *
         *  @param  event the event to analyze
         */
        void analyze(const art::Event &event);

        /**
         *  @brief  Process the subruns
         *
         *  @param  subrun the subrun to process
         */
        void endSubRun(const art::SubRun &subrun);

    private:

        EventFactory::Config   m_eventConfig;    ///< The event factory configuration
        SubrunFactory::Config  m_subrunConfig;   ///< The subrun factory configuration
        FileWriter             m_writer;         ///< The file writer
        Event                 *m_pEvent;         ///< The event bound to the output tree
        Subrun                *m_pSubrun;        ///< The subrun bound to the output tree
};

} // namespace ubcc1pi

DEFINE_ART_MODULE(ubcc1pi::AnalysisFileWriter)

#endif
