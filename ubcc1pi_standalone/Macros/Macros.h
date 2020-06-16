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
     *  @brief  Plot the variables that form the input to the particle ID BDTs after the CC inclusive selection
     *
     *  @param  config the input configuration
     */
    void PlotInputVariables(const Config &config = Config());

    /**
     *  @brief  Train the BDTs
     *
     *  @param  config the input configuration
     */
    void TrainBDTs(const Config &config = Config());

    /**
     *  @brief  Run the N-1 BDT study
     *
     *  @param  config the input configuration
     */
    void NMinusOneBDTStudy(const Config &config = Config());
    
    /**
     *  @brief  Make the event selection table
     *
     *  @param  config the input configuration
     */
    void MakeEventSelectionTable(const Config &config = Config());
    
    /**
     *  @brief  Make the PID table
     *
     *  @param  config the input configuration
     */
    void MakeSelectedPIDTable(const Config &config = Config());

    /**
     *  @brief  Make the selection efficiency plots
     *
     *  @param  config the input configuration
     */
    void MakeEventSelectionEfficiencyPlots(const Config &config = Config());

    /**
     *  @brief  Plot the reconstructed variables for the muon of the CC inclusive selection
     *
     *  @param  config the input configuration
     */
    void PlotCCInclusiveMuonRecoVariables(const Config &config = Config());

    /**
     *  @brief  Plot the recontructed variables
     *
     *  @param  config the input configuration
     */
    void PlotReconstructedVariables(const Config &config = Config());

    /**
     *  @brief  Extract the cross-sections
     *
     *  @param  config the input configuration
     */
    void ExtractXSecs(const Config &config = Config());

    /**
     *  @brief  Top-level macro that runs the full analysis chain
     *
     *  @param  config the input config
     */
    void RunFullAnalysis(const Config &config = Config());


} // namespace ubcc1pi_macros

#endif

