/**
 *  @file  ubcc1pi_standalone/Macros/ExtractXSecs.cxx
 *
 *  @brief The implementation file of the ExtractXSecs macro
 */

#include "ubcc1pi_standalone/Macros/Macros.h"
#include "ubcc1pi_standalone/Objects/FileReader.h"

#include "ubcc1pi_standalone/Helpers/AnalysisHelper.h"
#include "ubcc1pi_standalone/Helpers/NormalisationHelper.h"
#include "ubcc1pi_standalone/Helpers/CrossSectionHelper.h"

using namespace ubcc1pi;

namespace ubcc1pi_macros
{

void ExtractXSecs(const Config &config)
{
    //
    // Setup the input files
    //
    std::vector< std::tuple<AnalysisHelper::SampleType, std::string, float> > inputData;
    inputData.emplace_back(AnalysisHelper::Overlay, config.files.overlaysFileName, NormalisationHelper::GetOverlaysNormalisation(config));

    // Setup an object that holds the details of the systematic parameters to apply
    CrossSectionHelper::CrossSection::SystParams systParams;
    systParams.nBootstrapUniverses = config.extractXSecs.nBootstrapUniverses;
    systParams.fluxDimensions = config.extractXSecs.fluxDimensions;
    systParams.xsecDimensions = config.extractXSecs.xsecDimensions;
    systParams.detVarDimensions = config.extractXSecs.detVarDimensions;

    // Setup the flux reweightor
    CrossSectionHelper::FluxReweightor fluxReweightor(config.flux.binEdges, config.flux.energyBins, systParams.fluxDimensions);
    std::cout << "Integrated flux: " << fluxReweightor.GetIntegratedNominalFlux() << std::endl;

    // Loop over the events
    for (const auto &[sampleType, fileName, normalisation] : inputData)
    {
        std::cout << "Reading input file: " << fileName << std::endl;
        const auto isOverlay = (sampleType == AnalysisHelper::Overlay);

        // Open the input file for reading and enable the branches with systematic event weights (if required)
        FileReader reader(fileName);
        if (isOverlay) reader.EnableSystematicBranches();
        auto pEvent = reader.GetBoundEventAddress();

        const auto nEvents = reader.GetNumberOfEvents();
        for (unsigned int i = 0; i < nEvents; ++i)
        {
            AnalysisHelper::PrintLoadingBar(i, nEvents);
            reader.LoadEvent(i);

            // Get the nominal event weight, scaled by the sample normalisation
            const auto weight = AnalysisHelper::GetNominalEventWeight(pEvent) * normalisation;

            // Get the flux weights
            const auto fluxWeights = CrossSectionHelper::GetWeightsMap(pEvent->truth, systParams.fluxDimensions, config.extractXSecs.mutuallyExclusiveDimensions);

            // Add the true neutrino energy to the flux reweightor
            fluxReweightor.AddEvent(pEvent->truth.nuEnergy(), weight, fluxWeights);
        }
    }
}

} // namespace ubcc1pi_macros
