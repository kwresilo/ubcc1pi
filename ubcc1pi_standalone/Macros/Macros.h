/**
 *  @file  ubcc1pi_standalone/Macros/Macros.h
 *
 *  @brief The header file for the macros
 */

#ifndef UBCC1PI_STANDALONE_MACROS_MACROS
#define UBCC1PI_STANDALONE_MACROS_MACROS

#include "ubcc1pi_standalone/Objects/Config.h"

using namespace ubcc1pi;

namespace ubcc1pi_macros
{

    /**
     *  @brief  Macro to print the configuration
     *
     *  @param  config the input configuration
     */
    void PrintConfig(const Config &config = Config());
    
    /**
     *  @brief  Count the POT for the overlays and dirt samples and print the results
     *
     *  @param  config the input configuration
     */
    void CountPOT(const Config &config = Config());

    /**
     *  @brief  Get the run-subrun list that's needed to run Zarco's POT / trigger / spill counting tool
     *
     *  @param  config the input configuration
     */
    void GetRunSubrunList(const Config &config = Config());
    
    /**
     *  @brief  Top-level macro that runs the full analysis chain
     *
     *  @param  config the input config
     */
    void RunFullAnalysis(const Config &config = Config());


} // namespace ubcc1pi_macros

#endif

