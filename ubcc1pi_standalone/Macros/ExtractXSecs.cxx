/**
 *  @file  ubcc1pi_standalone/Macros/ExtractXSecs.cxx
 *
 *  @brief The implementation file of the ExtractXSecs macro
 */

#include "ubcc1pi_standalone/Macros/Macros.h"
#include "ubcc1pi_standalone/Objects/FileReader.h"

#include "ubcc1pi_standalone/Helpers/AnalysisHelper.h"
#include "ubcc1pi_standalone/Helpers/NormalisationHelper.h"
#include "ubcc1pi_standalone/Helpers/SelectionHelper.h"
#include "ubcc1pi_standalone/Helpers/CrossSectionHelper.h"

using namespace ubcc1pi;

namespace ubcc1pi_macros
{

void ExtractXSecs(const Config &config)
{
    // -------------------------------------------------------------------------------------------------------------------------------------
    // Setup the input files
    // -------------------------------------------------------------------------------------------------------------------------------------
    std::vector< std::tuple<AnalysisHelper::SampleType, std::string, std::string, float> > inputData;
    inputData.emplace_back(AnalysisHelper::Overlay, "", config.files.overlaysFileName, NormalisationHelper::GetOverlaysNormalisation(config));
    inputData.emplace_back(AnalysisHelper::Dirt,    "", config.files.dirtFileName, NormalisationHelper::GetDirtNormalisation(config));
    inputData.emplace_back(AnalysisHelper::DataEXT, "", config.files.dataEXTFileName, NormalisationHelper::GetDataEXTNormalisation(config));
    inputData.emplace_back(AnalysisHelper::DataBNB, "", config.files.dataBNBFileName, 1.f);

    // Add the detector variation files
    for (const auto &[name, fileName] : config.files.detVarFiles)
        inputData.emplace_back(AnalysisHelper::DetectorVariation, name, fileName, NormalisationHelper::GetDetectorVariationNormalisation(config, name));

    // -------------------------------------------------------------------------------------------------------------------------------------
    // Setup an object that holds the details of the systematic parameters to apply
    // -------------------------------------------------------------------------------------------------------------------------------------
    CrossSectionHelper::CrossSection::SystParams systParams;
    systParams.nBootstrapUniverses = config.extractXSecs.nBootstrapUniverses;
    systParams.fluxDimensions = config.extractXSecs.fluxDimensions;
    systParams.xsecDimensions = config.extractXSecs.xsecDimensions;
    systParams.detVarDimensions = config.extractXSecs.detVarDimensions;

    // -------------------------------------------------------------------------------------------------------------------------------------
    // Setup an object that holds the information about how we should scale an event rate to obtain a cross-section
    // -------------------------------------------------------------------------------------------------------------------------------------
    // - Flux             [10^-10 cm^-2 POT^-1]
    // - Exposure POT     [10^20 POT]               (stored in config as [POT])
    // - Target density   [10^31 nucleons/cm^3]     (stored in config as [10^23 nucleons/cm^3])
    // - Fiducial volume  [cm^3]
    //
    // - Cross-section    [Flux * Exposure * Target density * Fiducial volume]^-1 = [10^-41 cm^2 / nucleon]
    //
    // ATTN here we use a FluxReweightor to specify the flux. This is used to get the reweighted flux in each systematic universe.
    CrossSectionHelper::CrossSection::ScalingData scalingData;
    scalingData.pFluxReweightor = std::make_shared<CrossSectionHelper::FluxReweightor>(config.flux.binEdges, config.flux.energyBins, systParams.fluxDimensions);
    scalingData.exposurePOT = config.norms.dataBNBTor875WCut / (1e20);
    scalingData.nTargets = config.global.targetDensity * (1e-8) * AnalysisHelper::GetFiducialVolume();

    // -------------------------------------------------------------------------------------------------------------------------------------
    // Setup the event selection
    // -------------------------------------------------------------------------------------------------------------------------------------
    auto selection = SelectionHelper::GetDefaultSelection();

    // -------------------------------------------------------------------------------------------------------------------------------------
    // Setup the cross-section objects
    // -------------------------------------------------------------------------------------------------------------------------------------

    // ATTN here we use the machinary for a differential cross-section, and treat the total cross-section as a single-bin measurement.
    // The "kinematic quantity" in this case is just a dummy parameter. Here we define a single bin with edges arbitrarily chosen to be
    // (-1 -> +1), and we request that the cross-section object does not apply bin-width scaling. When we fill this object, we will use the
    // dummy kinematic quantity with a value of 0 for all events. Again this is arbitraty, as long as it's within the bin edges we chose.
    // In this way the single bin contains all events. It's just a trick to avoid implementing extra logic for the total cross-section.
    const auto dummyValue = 0.f;

    // Here we make a map from a name of the cross-section to the cross-section object itself. This way we can iterate through the
    // cross-section objects as desired.
    std::map<std::string, CrossSectionHelper::CrossSection> xsecMap;
    xsecMap.emplace("total", CrossSectionHelper::CrossSection(systParams, {-1.f, 1.f}, false, false, false));

    // -------------------------------------------------------------------------------------------------------------------------------------
    // Count the events
    // -------------------------------------------------------------------------------------------------------------------------------------
    // Loop over the files
    for (const auto &[sampleType, sampleName, fileName, normalisation] : inputData)
    {
        std::cout << "Reading input file: " << fileName << std::endl;

        const auto isOverlay = (sampleType == AnalysisHelper::Overlay);
        const auto isDirt    = (sampleType == AnalysisHelper::Dirt);
        const auto isDataBNB = (sampleType == AnalysisHelper::DataBNB);
        const auto isDetVar  = (sampleType == AnalysisHelper::DetectorVariation);

        // Open the input file for reading and enable the branches with systematic event weights (if required)
        FileReader reader(fileName);

        if (isOverlay)
            reader.EnableSystematicBranches();

        auto pEvent = reader.GetBoundEventAddress();

        // Loop over the events in the file
        const auto nEvents = reader.GetNumberOfEvents();
        for (unsigned int i = 0; i < nEvents; ++i)
        {
            AnalysisHelper::PrintLoadingBar(i, nEvents);
            reader.LoadEvent(i);

            // Run the selection
            const auto &[passedGoldenSelection, cutsPassed, assignedPdgCodes] = selection.Execute(pEvent);
            const bool passedGenericSelection = SelectionHelper::IsCutPassed(cutsPassed, config.global.lastCutGeneric);

            // Get the reco analysis data (if available, otherwise set to dummy values)
            const auto recoData = (
                passedGenericSelection
                    ? AnalysisHelper::GetRecoAnalysisData(pEvent->reco, assignedPdgCodes, passedGoldenSelection)
                    : AnalysisHelper::GetDummyAnalysisData()
            );

            // TODO apply reco-level phase-space restrictions
            const bool passesPhaseSpaceReco = true;
            const bool isSelectedGolden = passedGoldenSelection && passesPhaseSpaceReco;
            const bool isSelectedGeneric = passedGenericSelection && passesPhaseSpaceReco;

            // Handle BNB data
            if (isDataBNB)
            {
                // Count the events that are selected
                if (isSelectedGeneric)
                {
                    xsecMap.at("total").AddSelectedBNBDataEvent(dummyValue);
                }

                // For BNB data that's all we need to do!
                continue;
            }

            // Get the nominal event weight, scaled by the sample normalisation
            const auto weight = AnalysisHelper::GetNominalEventWeight(pEvent) * normalisation;

            // Determine if this is truly a CC1Pi event
            const auto isTrueCC1Pi = (isOverlay || isDetVar) && AnalysisHelper::IsTrueCC1Pi(pEvent, config.global.useAbsPdg);

            // Get the truth analysis data (if available, otherwise set to dummy values)
            const auto truthData = (
                isTrueCC1Pi
                    ? AnalysisHelper::GetTruthAnalysisData(pEvent->truth, config.global.useAbsPdg, config.global.protonMomentumThreshold)
                    : AnalysisHelper::GetDummyAnalysisData()
            );

            // TODO apply truth-level phase-space restrictions
            const bool passesPhaseSpaceTruth = true;
            const bool isSignal = isTrueCC1Pi && passesPhaseSpaceTruth;

            // Handle the detector variation samples (unisims)
            if (isDetVar)
            {
                // Handle signal events
                if (isSignal)
                {
                    xsecMap.at("total").AddSignalEventDetVar(dummyValue, dummyValue, isSelectedGeneric, weight, sampleName);
                }
                // Handle selected background events
                else if (isSelectedGeneric)
                {
                    xsecMap.at("total").AddSelectedBackgroundEventDetVar(dummyValue, weight, sampleName);
                }

                // For detector variation samples, that's all we need to do!
                continue;
            }

            // Get the flux weights
            const auto fluxWeights = (
                isOverlay
                    ? CrossSectionHelper::GetWeightsMap(pEvent->truth, systParams.fluxDimensions, config.extractXSecs.mutuallyExclusiveDimensions)
                    : CrossSectionHelper::GetUnitWeightsMap(systParams.fluxDimensions)
            );

            // Fill the flux-reweightor with all overlay events in the fiducial volume
            if (isOverlay && AnalysisHelper::IsFiducial(pEvent->truth.nuVertex()))
            {
                scalingData.pFluxReweightor->AddEvent(pEvent->truth.nuEnergy(), weight, fluxWeights);
            }

            // Get the cross-section weights
            // ATTN here we optionally scale the cross-section weights down by the genieTuneEventWeight - this is done so we don't
            // double count this weight (once one in the nominal event weight, and once in the xsec systematic event weights)
            const auto xsecWeightsScaleFactor = (isOverlay && config.extractXSecs.scaleXSecWeights) ? pEvent->truth.genieTuneEventWeight() : 1.f;
            const auto xsecWeights = (
                isOverlay
                    ? CrossSectionHelper::ScaleWeightsMap(CrossSectionHelper::GetWeightsMap(pEvent->truth, systParams.xsecDimensions, config.extractXSecs.mutuallyExclusiveDimensions), xsecWeightsScaleFactor)
                    : CrossSectionHelper::GetUnitWeightsMap(systParams.xsecDimensions)
            );

            // Handle signal events
            if (isSignal)
            {
                xsecMap.at("total").AddSignalEvent(dummyValue, dummyValue, isSelectedGeneric, weight, fluxWeights, xsecWeights);
            }
            // Handle selected background events
            else if (isSelectedGeneric)
            {
                xsecMap.at("total").AddSelectedBackgroundEvent(dummyValue, isDirt, weight, fluxWeights, xsecWeights);
            }
        }
    }

    // -------------------------------------------------------------------------------------------------------------------------------------
    // Calculate the cross-section
    // -------------------------------------------------------------------------------------------------------------------------------------

    // Loop over all cross-section objects
    for (const auto &[name, xsec] : xsecMap)
    {
        std::cout << "Processing cross-section: " << name << std::endl;

        // Get the cross-section as measured with BNB data along with it's uncertainties
        const auto data = xsec.GetBNBDataCrossSection(scalingData);
        const auto dataStatUncertainties = xsec.GetBNBDataCrossSectionStatUncertainty(scalingData);
        const auto dataSystBiasCovariances = xsec.GetBNBDataCrossSectionSystUncertainties(scalingData);

        // Get the predicted cross-section along with it's MC stat uncertainty
        const auto prediction = xsec.GetPredictedCrossSection(scalingData);
        const auto predictionStatBiasCovariance = xsec.GetPredictedCrossSectionStatUncertainty(scalingData);

        // TODO write this information out

        // Work out if this is a differential cross-section or a total cross-section. For the total cross-sections, we are done.
        const auto isDifferential = (xsec.GetMetadata().GetNBins() > 1);
        if (!isDifferential)
            continue;

        // Get the smearing matrix and along with it's uncertainties
        const auto smearingMatrix = xsec.GetSmearingMatrix();
        const auto smearingMatrixStatyBiasCovariances = xsec.GetSmearingMatrixSystUncertainties();

        // TODO write this information out
    }
}

} // namespace ubcc1pi_macros
