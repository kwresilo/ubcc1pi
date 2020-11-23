/**
 *  @file  ubcc1pi_standalone/Macros/ExtractCCInclusiveXSecs.cxx
 *
 *  @brief The implementation file of the ExtractCCInclusiveXSecs macro
 */

#include "ubcc1pi_standalone/Macros/Macros.h"
#include "ubcc1pi_standalone/Objects/FileReader.h"

#include "ubcc1pi_standalone/Helpers/AnalysisHelper.h"
#include "ubcc1pi_standalone/Helpers/SelectionHelper.h"
#include "ubcc1pi_standalone/Helpers/CrossSectionHelper.h"
#include "ubcc1pi_standalone/Helpers/NormalisationHelper.h"
#include "ubcc1pi_standalone/Helpers/FormattingHelper.h"

using namespace ubcc1pi;

namespace ubcc1pi_macros
{

void ExtractCCInclusiveXSecs(const Config &config)
{
    //
    // Setup the input files
    //
    std::vector< std::tuple<AnalysisHelper::SampleType, std::string, float, std::string, std::string> > inputData;

    inputData.emplace_back(AnalysisHelper::Overlay, config.files.overlaysFileName, NormalisationHelper::GetOverlaysNormalisation(config), "", "");
    inputData.emplace_back(AnalysisHelper::Dirt,    config.files.dirtFileName,     NormalisationHelper::GetDirtNormalisation(config), "", "");
    inputData.emplace_back(AnalysisHelper::DataEXT, config.files.dataEXTFileName,  NormalisationHelper::GetDataEXTNormalisation(config), "", "");
    inputData.emplace_back(AnalysisHelper::DataBNB, config.files.dataBNBFileName,  1.f, "", "");

    // Add the detector variation files
    for (const auto &[runId, detVarParam, fileName] : config.files.detVarFiles)
    {
        inputData.emplace_back(AnalysisHelper::DetectorVariation, fileName, NormalisationHelper::GetDetectorVariationNormalisation(config, runId, detVarParam), runId, detVarParam);
    }

    // Print the input data
    FormattingHelper::Table inputDataTable({"Type", "File", "Norm", "Events", "Events*Norm", "Run", "Detector parameter"});
    for (const auto &[sampleType, fileName, normalisation, runId, detVarParam] : inputData)
    {
        inputDataTable.AddEmptyRow();
        inputDataTable.SetEntry("Type", AnalysisHelper::GetSampleTypeName(sampleType));
        inputDataTable.SetEntry("File", fileName);
        inputDataTable.SetEntry("Norm", normalisation);
        inputDataTable.SetEntry("Run", runId.empty() ? "N/A" : runId);
        inputDataTable.SetEntry("Detector parameter", detVarParam.empty() ? "N/A" : detVarParam);

        FileReader reader(fileName);
        const auto nEvents = reader.GetNumberOfEvents();
        const auto normEvents = normalisation * nEvents;
        inputDataTable.SetEntry("Events", nEvents);
        inputDataTable.SetEntry("Events*Norm", normEvents);
    }
    inputDataTable.Print();

    // Get the selection
    auto selection = SelectionHelper::GetDefaultSelection();

    // Get the common input data and do any required unit conversions
    //
    // - Flux:            [10^-10 cm^-2 POT^-1]
    // - Exposure:        [10^20 POT]                 (stored in config as [POT])
    // - Target density:  [10^31 nucleons/cm^3]       (stored in config as [10^23 nucleons/cm^3])
    // - Fiducial volume: [cm^3]
    //
    // - Cross-section:   [Flux * Exposure * Target density * Fiducial volume]^-1 = [10^-41 cm^2 / nucleon]
    //
    CrossSectionHelper::CrossSection::InputData inputXSecData;
    inputXSecData.m_flux = CrossSectionHelper::GetFluxHistogram(config.flux);
    inputXSecData.m_exposurePOT = config.norms.dataBNBTor875WCut / (1e20);
    inputXSecData.m_nTargets = config.global.targetDensity * (1e-8) * AnalysisHelper::GetFiducialVolume();
    inputXSecData.m_fluxParameters = config.extractXSecs.fluxParams;

    // Print out the common input data
    std::cout << "flux = " << CrossSectionHelper::GetTotalFlux(inputXSecData.m_flux) << " 10^-10 cm^-2 POT^-1" << std::endl;
    std::cout << "exposure = " << inputXSecData.m_exposurePOT << " 10^20 POT" << std::endl;
    std::cout << "FV = " << AnalysisHelper::GetFiducialVolume() << " cm^3" << std::endl;
    std::cout << "nTargets = " << inputXSecData.m_nTargets << " 10^31 nucleons" << std::endl;

    // Build the systematic universe sizes map - this defines up-front the systematic parameters we want to apply and the number of
    // universes the we expect for each
    auto &systUniverseSizesMap = inputXSecData.m_systUniverseSizesMap;

    // Add the bootstrap universes (for the MC stat uncertainty)
    systUniverseSizesMap.emplace("bootstrap", config.extractXSecs.nBootstrapUniverses);

    // Add the independent systematic parameters
    for (const auto &[param, nUniverses] : config.extractXSecs.systematicParams)
        systUniverseSizesMap.emplace(param, nUniverses);

    // Add the mutually exclusive systematic parameters
    for (const auto &[param, nUniverses, exclusiveParams] : config.extractXSecs.mutuallyExclusiveSystematicParams)
        systUniverseSizesMap.emplace(param, nUniverses);

    // Build the list of detector variation parameters
    for (const auto &[runId, detVarParam, fileName] : config.files.detVarFiles)
    {
        // Store this as a parameter that the cross-section should expect
        inputXSecData.m_detVarParameters.emplace_back(runId, detVarParam);
    }

    //
    // Setup the cross-section objects
    //

    // Total cross-section - here we re-use the machinary of the differential cross-section by making a single bin measurement. The bin
    // edges we supply are just dummy values as we aren't measuring a kinematic paramter, here chosen to be [-1 -> 1]. We also request that
    // the cross-section isn't scaled by bin width.
    CrossSectionHelper::CrossSection xsec_total({-1.f, 1.f}, false, false, false, inputXSecData);

    // Loop over the events
    for (const auto &[sampleType, fileName, normalisation, runId, detVarParam] : inputData)
    {
        std::cout << "Reading input file: " << fileName << std::endl;
        const auto isOverlay = (sampleType == AnalysisHelper::Overlay);
        const auto isDataBNB = (sampleType == AnalysisHelper::DataBNB);
        const auto isDetVar = (sampleType == AnalysisHelper::DetectorVariation);

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

            // Get the systematic event weights non BNB data events
            CrossSectionHelper::SystematicWeightsMap systWeightsMap;
            if (!isDataBNB && !isDetVar)
            {
                if (isOverlay)
                {
                    // For overlay events get the event weights from the input file
                    CrossSectionHelper::AddSystematicWeights(pEvent->truth, config.extractXSecs.systematicParams, systWeightsMap);
                    CrossSectionHelper::AddMutuallyExclusiveSystematicWeights(pEvent->truth, config.extractXSecs.mutuallyExclusiveSystematicParams, systWeightsMap);

                    // Generate random boostrap weights
                    CrossSectionHelper::AddBootstrapWeights(config.extractXSecs.nBootstrapUniverses, systWeightsMap);

                    // When we get the nominal event weight, this includes a parameter to make our MC look like the MCC9 tune of GENIE
                    // This parameter is stored under event.truth.genieTuneEventWeight. In general we want to apply this weight to every
                    // event. But, the systematic weights relating to GENIE are calculated using the MCC9 tune - so if we apply those
                    // weights AND genieTuneEventWeight we will "double count". To account for this, here we scale down the GENIE systematic
                    // weights by genieTuneEventWeight so they can behave like all other weights.
                    const auto genieTuneEventWeight = pEvent->truth.splineEventWeight.IsSet() ? pEvent->truth.genieTuneEventWeight() : 0.f;
                    if (pEvent->truth.splineEventWeight.IsSet() && std::abs(genieTuneEventWeight) > std::numeric_limits<float>::epsilon())
                    {
                        for (const auto &genieParam : config.extractXSecs.genieParams)
                        {
                            // Get the weights for this genie parameter 
                            auto genieParamIter = systWeightsMap.find(genieParam);
                            if (genieParamIter == systWeightsMap.end())
                                throw std::logic_error("ExtractCCInclusiveXSecs - Supplied GENIE parameter: " + genieParam + " not found");

                            // Divide each weight down by the GENIE tune event weight
                            for (auto &weight : genieParamIter->second)
                                weight /= genieTuneEventWeight; 
                        }
                    }

                    // Add all overlay events in the fiducial volume to the cross-section objects so that we can use to total true neutrino
                    // energy spectrum in each universe to estimate the impact of the flux systematic variations on the total flux.
                    if (AnalysisHelper::IsFiducial(pEvent->truth.nuVertex()))
                    {
                        xsec_total.AddNeutrinoEvent(pEvent->truth.nuEnergy(), weight, systWeightsMap);
                    }
                }
                else
                {
                    // For non-overlay events just use unit weights
                    CrossSectionHelper::AddUnitWeights(config.extractXSecs.systematicParams, systWeightsMap);
                    CrossSectionHelper::AddUnitWeights(config.extractXSecs.mutuallyExclusiveSystematicParams, systWeightsMap);

                    // Generate random boostrap weights
                    CrossSectionHelper::AddBootstrapWeights(config.extractXSecs.nBootstrapUniverses, systWeightsMap);
                }
            }


            // Check if this event passed the CCInclusive selection
            auto isSelected = pEvent->reco.passesCCInclusive();

            // Determine if the event is truly a signal event
            const auto isSignal = (isOverlay || isDetVar) && AnalysisHelper::IsTrueCCInclusive(pEvent, config.global.useAbsPdg);

            // We only need to count signal events, or any event that passes the generic selection - everything else we can safely skip
            const bool shouldCount = isSignal || isSelected;
            if (!shouldCount)
                continue;

            // Count selected BNB data events
            if (isDataBNB)
            {
                if (isSelected)
                {
                    xsec_total.AddSelectedBNBDataEvent(0.f); // ATTN here a dummy value is used of 0.f
                }
            }
            else if (isDetVar)
            {
                if (isSignal)
                {
                    // Count selected signal events
                    xsec_total.AddSignalEventDetVar(0.f, 0.f, isSelected, weight, runId, detVarParam); // ATTN here a dummy value is used of 0.f
                }
                else
                {
                    // Count selected background events
                    if (isSelected)
                    {
                        xsec_total.AddSelectedBackgroundOverlayEventDetVar(0.f, weight, runId, detVarParam); // ATTN here a dummy value is used of 0.f
                    }
                }
            }
            else
            {
                if (isSignal)
                {
                    // Count selected signal events
                    xsec_total.AddSignalEvent(0.f, 0.f, isSelected, weight, systWeightsMap); // ATTN here a dummy value is used of 0.f
                }
                else
                {
                    // Count selected background events
                    if (isSelected)
                    {
                        xsec_total.AddSelectedBackgroundEvent(0.f, weight, systWeightsMap, isOverlay); // ATTN here a dummy value is used of 0.f
                    }
                }
            }
        }
    }

    // Get the flux variations for each parameter and save them
    auto pCanvas = PlottingHelper::GetCanvas();
    for (const auto &[paramName, pHist] : xsec_total.GetFluxVariations())
    {
        pHist->Draw("colz");
        inputXSecData.m_flux->Draw("hist same");
        PlottingHelper::SaveCanvas(pCanvas, "xsec_ccinc_fluxVariations_" + paramName);
    }

    // Calculate the cross-section for each kinematic parameter and print it out!
    for (auto &[namePrefix, xsec] : std::vector<std::pair<std::string, CrossSectionHelper::CrossSection> >({
        {"total", xsec_total}
    }))
    {
        FormattingHelper::PrintLine();
        std::cout << namePrefix << std::endl;
        FormattingHelper::PrintLine();

        // Get the cross-section itself
        std::cout << "Cross-section" << std::endl;
        const auto &xsec_plot = xsec.GetCrossSection();
        FormattingHelper::SaveHistAsTable(xsec_plot, "xsec_ccinc_" + namePrefix + ".md");

        // Get the stat uncertainties on the cross-seciton
        const auto &xsec_statErr = xsec.GetCrossSectionStatUncertainty();
        FormattingHelper::SaveHistAsTable(xsec_statErr, "xsec_ccinc_" + namePrefix + "_statUncertainty.md");

        // Get the covariance matrices and bias vectors for each systematic parameter
        std::cout << "Getting covariance & biases" << std::endl;
        const auto &xsec_covariances = xsec.GetCrossSectionCovarianceMatrices();
        for (const auto &[paramName, covarianceBiasPair] : xsec_covariances)
        {
            const auto &covarianceMatrix = covarianceBiasPair.first;
            const auto &biasVector = covarianceBiasPair.second;

            std::cout << "Bias vector" << std::endl;
            FormattingHelper::SaveHistAsTable(biasVector, "xsec_ccinc_" + namePrefix + "_" + paramName + "_bias.md");
            PlottingHelper::SaveBiasVector(biasVector, xsec_plot, xsec.HasUnderflowBin(), xsec.HasOverflowBin(), "xsec_ccinc_" + namePrefix + "_" + paramName);

            std::cout << "Covariance matrix" << std::endl;
            FormattingHelper::SaveHistAsTable(covarianceMatrix, "xsec_ccinc" + namePrefix + "_"+ paramName + "_covariance.md");
            PlottingHelper::SaveCovarianceMatrix(covarianceMatrix, xsec_plot, xsec.HasUnderflowBin(), xsec.HasOverflowBin(), "xsec_ccinc_" + namePrefix + "_" + paramName);
        }

        // Get the total covariance matrix, accounting for the biases
        std::cout << "Getting total covariance matrix" << std::endl;
        const auto &xsec_totalCovariance = CrossSectionHelper::GetTotalCovarianceMatrix(xsec_covariances);
        FormattingHelper::SaveHistAsTable(xsec_totalCovariance, "xsec_ccinc_" + namePrefix + "_allParams_covariance.md");
        PlottingHelper::SaveCovarianceMatrix(xsec_totalCovariance, xsec_plot, xsec.HasUnderflowBin(), xsec.HasOverflowBin(), "xsec_ccinc_" + namePrefix + "_allParams");

        // Save the cross-section plot with error bars
        PlottingHelper::SaveCrossSection(xsec_plot, xsec_statErr, xsec_totalCovariance, xsec.HasUnderflowBin(), xsec.HasOverflowBin(), "xsec_ccinc_" + namePrefix);

        // For performance clear the cached cross-sections from memory
        xsec.ClearCache();
    }
}

} // namespace ubcc1pi_macros
