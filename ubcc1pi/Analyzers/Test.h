/**
 *  @file  ubcc1pi/Analyzers/Test.h
 *
 *  @brief The header file for the test analyzer
 */

#ifndef UBCC1PI_ANALYZERS_TEST
#define UBCC1PI_ANALYZERS_TEST

#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDAnalyzer.h"

#include "art/Framework/Services/Optional/TFileService.h"
#include <TTree.h>

#include "ubcc1pi/Objects/FileWriter.h"

namespace ubcc1pi
{

/**
 *  @brief  The test analyzer class
 */
class Test : public art::EDAnalyzer
{
    public:
        /**
         *  @brief  Constructor
         */
        Test(const fhicl::ParameterSet &pset);

        /**
         *  @brief  Analyze an event
         *
         *  @param  event the event to analyze
         */
        void analyze(const art::Event &event);

    private:

        FileWriter  m_writer;   ///< The file writer
        Event      *m_pEvent;   ///< The event bound to the output tree
};

} // namespace ubcc1pi

DEFINE_ART_MODULE(ubcc1pi::Test)

#endif
