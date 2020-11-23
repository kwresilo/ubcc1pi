/**
 *  @file  ubcc1pi_standalone/Macros/RunFullAnalysis.cxx
 *
 *  @brief The implementation file of the RunFullAnalysis macro
 */

#include "ubcc1pi_standalone/Macros/Macros.h"

#include "ubcc1pi_standalone/Objects/FileReader.h"

#include "ubcc1pi_standalone/Helpers/NormalisationHelper.h"

using namespace ubcc1pi;

namespace ubcc1pi_macros
{

void RunFullAnalysis(const Config &config)
{
    PrintConfig(config);
    std::cout << "TODO - Still need to put all the steps of the analysis together into this RunFullAnalysis macro" << std::endl;

    /*
    TruthStudy(config);
    CCInclusiveMuonPIDStudy(config);
    MultiPlanePIDDemo(config);
    PlotInputVariables(config);
    GetCorrelationPlots(config);
    TrainBDTs(config);

    MakeEventSelectionTable(config);

    // Run this macro twice using the generic and the golden selection
    auto customConfig = config;
    customConfig.makeSelectedPIDTable.useGenericSelection = true;
    MakeSelectedPIDTable(customConfig);

    customConfig.makeSelectedPIDTable.useGenericSelection = false;
    MakeSelectedPIDTable(customConfig);

    MakeEventSelectionEfficiencyPlots(config);
    PlotReconstructedVariables(config);
    ExtractXSecs(config);
    */
}

}

