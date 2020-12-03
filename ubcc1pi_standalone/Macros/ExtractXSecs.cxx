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
#include "ubcc1pi_standalone/Helpers/FormattingHelper.h"

using namespace ubcc1pi;

namespace ubcc1pi_macros
{

void ExtractXSecs(const Config &config)
{
    // -------------------------------------------------------------------------------------------------------------------------------------
    // Setup the input files
    // -------------------------------------------------------------------------------------------------------------------------------------
    std::cout << "Setting up input files" << std::endl;

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
    std::cout << "Setting up systematic parameters" << std::endl;

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
    std::cout << "Setting up scaling data" << std::endl;

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
    std::cout << "Setting up cross-section objects" << std::endl;

    // Here we make a map from a name of the cross-section to the cross-section object itself.
    // This way we can iterate through the cross-section objects and reduce code-bloat.
    // The first index is an identifier for the selection that's applied (generic or goldlen), the second index is an identifier for the
    // kinematic quantity that's relevant for the cross-section (e.g. muonMomentum), and the mapped type is the cross-section object
    std::map<std::string, std::map<std::string, CrossSectionHelper::CrossSection> > xsecMap;

    for (const auto &selectionName : {"generic", "golden"})
    {
        // Add the differential cross-sections
        for (const auto &[name, binning, scaleByBinWidth] : std::vector< std::tuple<std::string, Config::Global::Binning, bool> > {

            // The names of the cross-section kinematic parameters, and their binning information. The third (boolean) parameter indicates
            // if the cross-section bins should be scaled by their width
            { "muonMomentum", config.global.muonMomentum, true }

        }){
            // Add the cross-section object to the map using the binning from the input configuration
            const auto &[extendedBinEdges, hasUnderflow, hasOverflow] = CrossSectionHelper::GetExtendedBinEdges(binning.min, binning.max, binning.binEdges);
            xsecMap[selectionName].emplace(name, CrossSectionHelper::CrossSection(systParams, extendedBinEdges, hasUnderflow, hasOverflow, scaleByBinWidth));
        }

        // ATTN here we use the machinary for a differential cross-section, and treat the total cross-section as a single-bin measurement.
        // The "kinematic quantity" in this case is just a dummy parameter. Here we define a single bin with edges arbitrarily chosen to be
        // (-1 -> +1), and we request that the cross-section object does not apply bin-width scaling. When we fill this object, we will use
        // the dummy kinematic quantity with a value of 0 for all events. This is arbitrary, as long as it's within the bin edges we chose.
        // In this way the single bin contains all events. It's just a trick to avoid implementing extra logic for the total cross-section.
        xsecMap[selectionName].emplace("total", CrossSectionHelper::CrossSection(systParams, {-1.f, 1.f}, false, false, false));
    }

    // The dummy value that will be used as the "kinematic quantity" for the total cross-section
    const auto dummyValue = 0.f;

    // -------------------------------------------------------------------------------------------------------------------------------------
    // Setup the relevent "getters" for each cross-section
    // -------------------------------------------------------------------------------------------------------------------------------------
    // Here we define a map from the name of each cross-section to a function which pulls out the relevant kinematic quanitity from an input
    // analysis data object. Again this is done up-front to reduce code-bloat below.
    std::unordered_map< std::string, std::function<float(const AnalysisHelper::AnalysisData &)> > getValue;

    // Differential cross-section kinematic parameters
    getValue.emplace("muonMomentum", [](const auto &data) { return data.muonMomentum; });

    // ATTN as described above, for the total cross-section we don't have an associated kinematic quantity so we just return a dummy value
    getValue.emplace("total", [=](const auto &) { return dummyValue; });

    // -------------------------------------------------------------------------------------------------------------------------------------
    // Count the events
    // -------------------------------------------------------------------------------------------------------------------------------------
    std::cout << "Counting events" << std::endl;

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

            // -----------------------------------------------------------------------------------------------------------------------------
            // Work out if this event passed the selection and apply any additional phase-space cuts based on the input binning
            // -----------------------------------------------------------------------------------------------------------------------------

            // Run the selection
            const auto &[passedGoldenSelection, cutsPassed, assignedPdgCodes] = selection.Execute(pEvent);
            const auto passedGenericSelection = SelectionHelper::IsCutPassed(cutsPassed, config.global.lastCutGeneric);

            // Get the reco analysis data (if available, otherwise set to dummy values)
            const auto recoData = (
                passedGenericSelection
                    ? AnalysisHelper::GetRecoAnalysisData(pEvent->reco, assignedPdgCodes, passedGoldenSelection)
                    : AnalysisHelper::GetDummyAnalysisData()
            );

            // Here we apply reco-level phase-space restrictions
            // For any event that passes the generic selection, get the value of the kinematic quantity and check if it is outside of the
            // min/max values supplied in the binning. If so, then reject the event.
            bool passesPhaseSpaceReco = false;
            if (passedGenericSelection)
            {
                // Start by assuming the event passes the phase-space cuts
                passesPhaseSpaceReco = true;

                // Check the value of the kinematic quantity for each cross-section is within the supplied bin limits
                // Note that these bin limits include the underflow & overflow
                for (const auto &[selectionName, xsecs] : xsecMap)
                {
                    for (const auto &[name, xsec] : xsecs)
                    {
                        const auto &binEdges = xsec.GetBinEdges();
                        const auto value = getValue.at(name)(recoData);
                        if (value < binEdges.front() || value > binEdges.back())
                        {
                            passesPhaseSpaceReco = false;
                            break;
                        }
                    }

                    if (!passesPhaseSpaceReco)
                        break;
                }
            }

            const auto isSelectedGolden = passedGoldenSelection && passesPhaseSpaceReco;
            const auto isSelectedGeneric = passedGenericSelection && passesPhaseSpaceReco;

            // -----------------------------------------------------------------------------------------------------------------------------
            // Handle BNB data
            // -----------------------------------------------------------------------------------------------------------------------------
            if (isDataBNB)
            {
                for (auto &[selectionName, xsecs] : xsecMap)
                {
                    // Determine if we passed the relevant selection
                    const auto isSelected = (selectionName == "golden" ? isSelectedGolden : isSelectedGeneric);

                    // Only count events passing the selection
                    if (!isSelected)
                        continue;

                    for (auto &[name, xsec] : xsecs)
                    {
                        xsec.AddSelectedBNBDataEvent(getValue.at(name)(recoData));
                    }
                }

                // For BNB data that's all we need to do!
                continue;
            }

            // -----------------------------------------------------------------------------------------------------------------------------
            // Work out if this event is signal, and apply any phase-space restrictions based on the input binning
            // -----------------------------------------------------------------------------------------------------------------------------

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

            // Here we apply truth-level phase-space restrictions
            // For all true CC1Pi events, we check if the values of each kinematic variable are within the supplied limits. If not then the
            // event is not classed as "signal"
            bool passesPhaseSpaceTruth = false;
            if (isTrueCC1Pi)
            {
                // Start by assuming the event passes the phase-space cuts
                passesPhaseSpaceTruth = true;

                // Check the value of the kinematic quantity for each cross-section is within the supplied bin limits
                // Note that these bin limits include the underflow & overflow
                for (const auto &[selectionName, xsecs] : xsecMap)
                {
                    for (const auto &[name, xsec] : xsecs)
                    {
                        const auto &binEdges = xsec.GetBinEdges();
                        const auto value = getValue.at(name)(truthData);
                        if (value < binEdges.front() || value > binEdges.back())
                        {
                            passesPhaseSpaceTruth = false;
                            break;
                        }
                    }

                    if (!passesPhaseSpaceTruth)
                        break;
                }
            }

            const auto isSignal = isTrueCC1Pi && passesPhaseSpaceTruth;

            // -----------------------------------------------------------------------------------------------------------------------------
            // Handle the detector variation samples (unisims)
            // -----------------------------------------------------------------------------------------------------------------------------
            if (isDetVar)
            {
                // Handle signal events
                if (isSignal)
                {
                    for (auto &[selectionName, xsecs] : xsecMap)
                    {
                        // Determine if we passed the relevant selection
                        const auto isSelected = (selectionName == "golden" ? isSelectedGolden : isSelectedGeneric);

                        for (auto &[name, xsec] : xsecs)
                        {
                            const auto recoValue = getValue.at(name)(recoData);
                            const auto trueValue = getValue.at(name)(truthData);

                            xsec.AddSignalEventDetVar(recoValue, trueValue, isSelected, weight, sampleName);
                        }
                    }
                }
                // Handle selected background events
                else
                {
                    for (auto &[selectionName, xsecs] : xsecMap)
                    {
                        // Only use selected background events
                        const auto isSelected = (selectionName == "golden" ? isSelectedGolden : isSelectedGeneric);
                        if (!isSelected)
                            continue;

                        for (auto &[name, xsec] : xsecs)
                        {
                            const auto recoValue = getValue.at(name)(recoData);
                            xsec.AddSelectedBackgroundEventDetVar(recoValue, weight, sampleName);
                        }
                    }
                }

                // For detector variation samples, that's all we need to do!
                continue;
            }

            // -----------------------------------------------------------------------------------------------------------------------------
            // Handle all other events (i.e those from the nominal simulation): Overlays, dirt, EXT data
            // -----------------------------------------------------------------------------------------------------------------------------

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
            // double count this weight (once in the nominal event weight, and once in the xsec systematic event weights)
            const auto xsecWeightsScaleFactor = (isOverlay && config.extractXSecs.scaleXSecWeights) ? pEvent->truth.genieTuneEventWeight() : 1.f;
            const auto xsecWeights = (
                isOverlay
                    ? CrossSectionHelper::ScaleWeightsMap(CrossSectionHelper::GetWeightsMap(pEvent->truth, systParams.xsecDimensions, config.extractXSecs.mutuallyExclusiveDimensions), xsecWeightsScaleFactor)
                    : CrossSectionHelper::GetUnitWeightsMap(systParams.xsecDimensions)
            );

            // Handle signal events
            if (isSignal)
            {
                for (auto &[selectionName, xsecs] : xsecMap)
                {
                    // Determine if we passed the relevant selection
                    const auto isSelected = (selectionName == "golden" ? isSelectedGolden : isSelectedGeneric);

                    for (auto &[name, xsec] : xsecs)
                    {
                        const auto recoValue = getValue.at(name)(recoData);
                        const auto trueValue = getValue.at(name)(truthData);

                        xsec.AddSignalEvent(recoValue, trueValue, isSelected, weight, fluxWeights, xsecWeights);
                    }
                }
            }
            // Handle selected background events
            else
            {
                for (auto &[selectionName, xsecs] : xsecMap)
                {
                    // Only use selected background events
                    const auto isSelected = (selectionName == "golden" ? isSelectedGolden : isSelectedGeneric);
                    if (!isSelected)
                        continue;

                    for (auto &[name, xsec] : xsecs)
                    {
                        const auto recoValue = getValue.at(name)(recoData);
                        xsec.AddSelectedBackgroundEvent(recoValue, isDirt, weight, fluxWeights, xsecWeights);
                    }
                }
            }
        }
    }

    // -------------------------------------------------------------------------------------------------------------------------------------
    // Calculate the cross-sections
    // -------------------------------------------------------------------------------------------------------------------------------------

    // Loop over all cross-section objects
    for (const auto &[selectionName, xsecs] : xsecMap)
    {
        for (const auto &[name, xsec] : xsecs)
        {
            std::cout << "Processing cross-section: " << name << std::endl;

            // -----------------------------------------------------------------------------------------------------------------------------
            // Get the cross-section as measured with BNB data along with it's uncertainties
            // -----------------------------------------------------------------------------------------------------------------------------
            const auto data = xsec.GetBNBDataCrossSection(scalingData);
            std::cout << "BNB data cross-section (reco-space)" << std::endl;
            FormattingHelper::SaveMatrix(data, "xsec_" + selectionName + "_" + name + "_data.txt");

            const auto dataStatUncertainties = xsec.GetBNBDataCrossSectionStatUncertainty(scalingData);
            std::cout << "BNB data stat uncertainty" << std::endl;
            FormattingHelper::SaveMatrix(dataStatUncertainties, "xsec_" + selectionName + "_" + name + "_data_stat.txt");

            const auto dataSystBiasCovariances = xsec.GetBNBDataCrossSectionSystUncertainties(scalingData);
            for (const auto &[group, map] : dataSystBiasCovariances)
            {
                for (const auto &[paramName, biasCovariance] : map)
                {
                    const auto &[pBias, pCovariance] = biasCovariance;

                    std::cout << "BNB data syst uncertainty: " << group << " " << paramName << std::endl;
                    std::cout << "Bias vector" << std::endl;
                    FormattingHelper::SaveMatrix(*pBias, "xsec_" + selectionName + "_" + name + "_data_" + group + "_" + paramName + "_bias.txt");
                    std::cout << "Covariance matrix" << std::endl;
                    FormattingHelper::SaveMatrix(*pCovariance, "xsec_" + selectionName + "_" + name + "_data_" + group + "_" + paramName + "_covariance.txt");
                }
            }

            // -----------------------------------------------------------------------------------------------------------------------------
            // Get the predicted cross-section along with it's MC stat uncertainty
            // -----------------------------------------------------------------------------------------------------------------------------
            const auto prediction = xsec.GetPredictedCrossSection(scalingData);
            std::cout << "Predicted cross-section (truth-space)" << std::endl;
            FormattingHelper::SaveMatrix(prediction, "xsec_" + selectionName + "_" + name + "_prediction.txt");

            const auto &[pPredictionStatBias, pPredictionStatCovariance] = xsec.GetPredictedCrossSectionStatUncertainty(scalingData);
            std::cout << "Predicted cross-section stat uncertainty" << std::endl;
            std::cout << "Bias vector" << std::endl;
            FormattingHelper::SaveMatrix(*pPredictionStatBias, "xsec_" + selectionName + "_" + name + "_prediction_stat_bias.txt");
            std::cout << "Covariance matrix" << std::endl;
            FormattingHelper::SaveMatrix(*pPredictionStatCovariance, "xsec_" + selectionName + "_" + name + "_prediction_stat_covariance.txt");

            // -----------------------------------------------------------------------------------------------------------------------------
            // Get the smearing matrix and along with it's uncertainties
            // -----------------------------------------------------------------------------------------------------------------------------
            const auto smearingMatrix = xsec.GetSmearingMatrix();
            std::cout << "Smearing matrix (reco-space rows, truth-space columns)" << std::endl;
            FormattingHelper::SaveMatrix(smearingMatrix, "xsec_" + selectionName + "_" + name + "_smearingMatrix.txt");

            const auto smearingMatrixSystBiasCovariances = xsec.GetSmearingMatrixSystUncertainties();
            for (const auto &[group, map] : smearingMatrixSystBiasCovariances)
            {
                for (const auto &[paramName, biasCovariance] : map)
                {
                    const auto &[pBias, pCovariance] = biasCovariance;

                    std::cout << "Smearing matrix syst uncertainty: " << group << " " << paramName << std::endl;
                    std::cout << "Bias vector" << std::endl;
                    FormattingHelper::SaveMatrix(*pBias, "xsec_" + selectionName + "_" + name + "_smearingMatrix_" + group + "_" + paramName + "_bias.txt");
                    std::cout << "Covariance matrix" << std::endl;
                    FormattingHelper::SaveMatrix(*pCovariance, "xsec_" + selectionName + "_" + name + "_smearingMatrix_" + group + "_" + paramName + "_covariance.txt");
                }
            }
        }
    }
}

} // namespace ubcc1pi_macros
