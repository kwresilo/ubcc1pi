/**
 *  @file  ubcc1pi/Analyzers/TruthStudy.h
 *
 *  @brief The header file for the truth study analyzer.
 */

#ifndef UBCC1PI_ANALYZERS_TRUTH_STUDY
#define UBCC1PI_ANALYZERS_TRUTH_STUDY

#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDAnalyzer.h"

namespace ubcc1pi
{

/**
 *  @brief  The truth study class
 */
class TruthStudy : public art::EDAnalyzer
{
    public:
        /**
         *  @brief  Constructor
         *
         *  @param  config the set of input fhicl parameters
         */
        TruthStudy(const fhicl::ParameterSet &config);

        /**
         *  @brief  Analyze an event
         *
         *  @param  event the event to analyze
         */
        void analyze(const art::Event &event);
};

} // namespace ubcc1pi

DEFINE_ART_MODULE(ubcc1pi::TruthStudy)

#endif
