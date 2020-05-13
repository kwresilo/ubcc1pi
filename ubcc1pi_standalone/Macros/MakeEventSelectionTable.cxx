#include "ubcc1pi_standalone/Objects/FileReader.h"

#include "ubcc1pi_standalone/Helpers/SelectionHelper.h"
#include "ubcc1pi_standalone/Helpers/AnalysisHelper.h"
#include "ubcc1pi_standalone/Helpers/BDTHelper.h"

using namespace ubcc1pi;

int MakeEventSelectionTable(const std::string &overlayFileName, const float overlayWeight, const std::string &dataEXTFileName, const float extWeight, const std::string &dataBNBFileName, const bool shouldOptimize = false, const unsigned int nScanPoints = 20u, const float processFraction = 0.2f)
{
    // Set up the selection
    auto selection = SelectionHelper::GetDefaultSelection();
    
    // Set the parameters to be optimized
    if (shouldOptimize)
    {
        selection.EnableOptimization("2NonProtons", false, -0.4f, 0.4f); // This is optimized for golden signal events
        selection.EnableOptimization("openingAngle", true, 2.3f, 3.14f);
        selection.EnableOptimization("topologicalScore", true, 0.f, 0.8f);
        selection.EnableOptimization("noShowers", true, 0.f, 0.2f);
        selection.EnableOptimization("likelyGoldenPion", true, -0.2f, 0.2f, "S G"); // This is optimized for golden signal events
    }
    
    // Get bodge factors
    const float bodgeFactor = 28678.f / 22534.f; // ATTN this factor is a normalisation added so we can compare the shape of the distributions, can't exist in the final result!
    const auto overlayBodgeWeight = overlayWeight * bodgeFactor;
    const auto dataEXTBodgeWeight = extWeight * bodgeFactor;
    
    // Optimize the cuts
    if (shouldOptimize)
        selection.Optimize(dataBNBFileName, overlayFileName, overlayBodgeWeight, dataEXTFileName, dataEXTBodgeWeight, "", 0.f, nScanPoints, processFraction);

    // Run the selection
    selection.Execute(dataBNBFileName, overlayFileName, overlayBodgeWeight, dataEXTFileName, dataEXTBodgeWeight, "", 0.f, true, 1.f, 10u);

    return 0;
}
