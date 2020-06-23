/**
 *  @file  ubcc1pi_standalone/Macros/MakeEventSelectionTable.cxx
 *
 *  @brief The implementation file of the MakeEventSelectionTable macro
 */

#include "ubcc1pi_standalone/Macros/Macros.h"

#include "ubcc1pi_standalone/Objects/FileReader.h"

#include "ubcc1pi_standalone/Helpers/SelectionHelper.h"
#include "ubcc1pi_standalone/Helpers/AnalysisHelper.h"
#include "ubcc1pi_standalone/Helpers/BDTHelper.h"
#include "ubcc1pi_standalone/Helpers/NormalisationHelper.h"

using namespace ubcc1pi;

namespace ubcc1pi_macros
{

//int MakeEventSelectionTable(const std::string &overlayFileName, const float overlayWeight, const std::string &dataEXTFileName, const float extWeight, const std::string &dataBNBFileName, const bool shouldOptimize = false, const unsigned int nScanPoints = 20u, const float processFraction = 0.2f)
void MakeEventSelectionTable(const Config &config) 
{
    // Set up the selection
    auto selection = SelectionHelper::GetDefaultSelection();
    
    // Set the parameters to be optimized
    if (config.makeEventSelectionTable.shouldOptimize)
    {
        selection.EnableOptimization("2NonProtons", false, -0.4f, 0.4f);
        selection.EnableOptimization("openingAngle", true, 2.3f, 3.14f);
        selection.EnableOptimization("topologicalScore", true, 0.f, 0.8f);
        selection.EnableOptimization("noShowers", true, 0.f, 0.2f);
        selection.EnableOptimization("likelyGoldenPion", true, -0.2f, 0.2f, "S G"); // This is optimized for golden signal events
    }
    
    // Get bodge factors
    const auto overlayWeight = NormalisationHelper::GetOverlaysNormalisation(config);
    const auto dataEXTWeight = NormalisationHelper::GetDataEXTNormalisation(config);
    const auto dirtWeight = NormalisationHelper::GetDirtNormalisation(config);
    
    // Optimize the cuts
    if (config.makeEventSelectionTable.shouldOptimize)
        selection.Optimize(config.files.dataBNBFileName, config.files.overlaysFileName, overlayWeight, config.files.dataEXTFileName, dataEXTWeight, config.files.dirtFileName, dirtWeight, config.makeEventSelectionTable.nScanPoints, config.makeEventSelectionTable.processFraction);

    // Run the selection
    selection.Execute(config.files.dataBNBFileName, config.files.overlaysFileName, overlayWeight, config.files.dataEXTFileName, dataEXTWeight, config.files.dirtFileName, dirtWeight, true, 1.f, 10u);
}

} // namespace ubcc1pi_macros
