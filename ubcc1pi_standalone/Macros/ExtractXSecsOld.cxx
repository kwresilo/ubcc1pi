/**
 *  @file  ubcc1pi_standalone/Macros/ExtractXSecsOld.cxx
 *
 *  @brief The implementation file of the ExtractXSecsOld macro
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

void ExtractXSecsOld(const Config &config)
{
    // -------------------------------------------------------------------------------------------------------------------------------------
    // Setup the input files
    // -------------------------------------------------------------------------------------------------------------------------------------
    // Here we define a vector of tuples with 4 entries
    //   - First, the sample type (e.g. overlay)
    //   - Second, a string which is used to identify a given detector variation sample (for other sample type, this is unused)
    //   - Third, the path to the input file
    //   - Fourth, the normalisation factor to apply to all events in that file
    std::vector< std::tuple<AnalysisHelper::SampleType, std::string, std::string, float> > inputData;
    inputData.emplace_back(AnalysisHelper::Overlay, "", config.filesRun1.overlaysFileName, NormalisationHelper::GetOverlaysNormalisation(config, 1));
    inputData.emplace_back(AnalysisHelper::Dirt,    "", config.filesRun1.dirtFileName, NormalisationHelper::GetDirtNormalisation(config, 1));
    inputData.emplace_back(AnalysisHelper::DataEXT, "", config.filesRun1.dataEXTFileName, NormalisationHelper::GetDataEXTNormalisation(config, 1));
    inputData.emplace_back(AnalysisHelper::DataBNB, "", config.filesRun1.dataBNBFileName, 1.f);

    // Add the detector variation files //Todo uncomment !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // for (const auto &[name, fileName] : config.filesRun1.detVarFiles)
    //     inputData.emplace_back(AnalysisHelper::DetectorVariation, name, fileName, NormalisationHelper::GetDetectorVariationNormalisation(config, name, 1));

    // -------------------------------------------------------------------------------------------------------------------------------------
    // Setup an object that holds the details of the systematic parameters to apply
    // -------------------------------------------------------------------------------------------------------------------------------------
    // Here we read in the "dimensions" of the systematic parameters to apply. For multisim parameters (flux & xsec), this is a map from
    // each parameter name to the expected number of universes. For unisim parameters (detector variations), this is a map from the
    // identidier (i.e. name) of the detector variation sample, to the identifier of the corresponding central-value sample. Additionally
    // we set the number of bootstrap universes (for the MC stat uncertainty), the corresponding weights are generated for each event.
    CrossSectionHelper::CrossSection::SystParams systParams;
    systParams.nBootstrapUniverses = config.extractXSecs.nBootstrapUniverses;
    systParams.fluxDimensions = config.extractXSecs.fluxDimensions;
    systParams.xsecDimensions = config.extractXSecs.xsecDimensions;
    systParams.reintDimensions = config.extractXSecs.reintDimensions;
    systParams.detVarDimensions = config.extractXSecs.detVarDimensions;
    systParams.potFracUncertainty = config.extractXSecs.potFracUncertainty;

    // -------------------------------------------------------------------------------------------------------------------------------------
    // Setup an object that holds the information about how we should scale an event rate to obtain a cross-section
    // -------------------------------------------------------------------------------------------------------------------------------------
    // Here we specify the:
    // - Flux             [10^-10 cm^-2 POT^-1]
    // - Exposure POT     [10^20 POT]               (stored in config as [POT])
    // - Target density   [10^31 nucleons/cm^3]     (stored in config as [10^23 nucleons/cm^3])
    // - Fiducial volume  [cm^3]
    //
    // Hence the units of the eventual cross-seciton are:
    // - Cross-section    [Flux * Exposure * Target density * Fiducial volume]^-1 = [10^-41 cm^2 / nucleon]
    //
    // Here we use a FluxReweightor to specify the flux. For each universe of each flux systematic paramter, this is uses the ratio of the
    // total neutrino event rate (as a function of the true neutrino energy) to the same rate in the nominal simulation to reweight the
    // input neutrino flux distribution. The integrated flux is each universe is used to scale the selected event rate when calculating the
    // cross-section in that universe. For all non-flux parameters, the nominal integrated flux is used.
    CrossSectionHelper::CrossSection::ScalingData scalingData;

    // Read in the flux
    const auto fluxHistNames = CrossSectionHelper::GetNominalFluxHistNames(config.fluxRun1.nuPdgsSignal, config.fluxRun1.nuPdgToHistName, config.fluxRun1.nomHistPattern);
    const auto &[fluxBinEdges, fluxValues] = CrossSectionHelper::ReadNominalFlux(config.fluxRun1.fileName, fluxHistNames, config.fluxRun1.pot);
    std::cout << "Flux bin edges: " << std::endl;
    for(const auto &edges: fluxBinEdges)
        std::cout << edges << std::endl;
    scalingData.pFluxReweightor = std::make_shared<CrossSectionHelper::FluxReweightor>(fluxBinEdges, fluxValues, systParams.fluxDimensions);
    const auto debugFlux = scalingData.pFluxReweightor->GetIntegratedNominalFlux();
    std::cout << "\n- Flux:            " << debugFlux << " * 10^-10 cm^-2 POT^-1\n" << std::endl;
    scalingData.exposurePOT = config.normsRun1.dataBNBTor875WCut / (1e20);
    scalingData.nTargets = config.global.targetDensity * (1e-8) * AnalysisHelper::GetFiducialVolume();

    // std::cout << "- Flux:            " << scalingData.pFluxReweightor->GetIntegratedNominalFlux() << " * 10^-10 cm^-2 POT^-1" << std::endl;
    // std::cout << "- Exposure:        " << scalingData.exposurePOT << " * 10^20 POT" << std::endl;
    // std::cout << "- Target density:  " << config.global.targetDensity << " * 10^31 nucleons/cm^3" << std::endl;
    // std::cout << "- Fiducial volume: " << AnalysisHelper::GetFiducialVolume() << " * cm^3" << std::endl;
    // std::cout << "- nTargets:        " << scalingData.nTargets << " * 10^31 nucleons" << std::endl;

    // -------------------------------------------------------------------------------------------------------------------------------------
    // Setup the event selection
    // -------------------------------------------------------------------------------------------------------------------------------------
    // Here we are using the default CC1pi selection. The final cut of this selection (named "likelyGoldenPion") is used to identify events
    // in which the reconstructed & selected pion candidate is "golden". Here a golden pion is one that is: contained within the TPC,
    // doesn't undergo any scatters on the argon, and comes to rest (before any secondary interaction). If an event passes all the cuts of
    // this selection, it's said to pass the "golden" selection. If an event passes the penultimate cut, it is said to have passed the
    // "generic" selection. Events passing the generic selection are likely to be CC1pi, but may or may not have a golden pion. All events
    // that pass the golden selection also pass the generic selection - but not all events that pass the generic selection also pass the
    // golden selection.
    auto selection = SelectionHelper::GetDefaultSelection();

    // -------------------------------------------------------------------------------------------------------------------------------------
    // Setup the cross-section objects
    // -------------------------------------------------------------------------------------------------------------------------------------
    // Here we make a map from a name of the cross-section to the cross-section object itself. In this way, we can iterate through the
    // cross-section objects and reduce code-bloat. The first index is an identifier for the selection that's applied (generic or goldlen),
    // the second index is an identifier for the kinematic quantity that's relevant for the cross-section (e.g. muonMomentum), and the
    // mapped type is the cross-section object.
    std::map<std::string, std::map<std::string, CrossSectionHelper::CrossSection> > xsecMap;

    // We additionally make a map from each cross-section to the limits of the phase-space that we should consider. The key is the
    // identifier for the kinematic quantity and the mapped value is a pair containing the limits [min, max]
    std::map< std::string, std::pair<float, float> > phaseSpaceMap;

    // ATTN the configuration allows the user to enable or disable each cross-section. If a cross-section has been disabled, then it won't
    // be added the xSecMap. However, the phaseSpaceMap always includes all kinematic paramters!

    // Add the differential cross-sections
    for (const auto &[name, binning, scaleByBinWidth] : std::vector< std::tuple<std::string, Config::Global::Binning, bool> > {

        // The names of the cross-section kinematic parameters, and their binning information.
        // The third (boolean) parameter indicates if the cross-section bins should be scaled by their width
        { "muonCosTheta",  config.global.muonCosTheta,  true  },
        { "muonPhi",       config.global.muonPhi,       true  },
        { "muonMomentum",  config.global.muonMomentum,  true  },

        { "pionCosTheta",  config.global.pionCosTheta,  true  },
        { "pionPhi",       config.global.pionPhi,       true  },
        { "pionMomentum",  config.global.pionMomentum,  true  },

        { "muonPionAngle", config.global.muonPionAngle, true  },
        { "nProtons",      config.global.nProtons,      false }

    })
    {
        // Add to the phase-space map
        phaseSpaceMap.emplace(name, std::pair<float, float>({binning.min, binning.max}));

        // Here we calculate every cross-section using both the generic and golden selection. In the end we only use the golden selection for
        // the pion momentum, but we additionally apply it to other cross-sections as a cross-check.
        for (const auto &selectionName : {"generic", "golden"})
        {
            // Don't setup a cross-section object if it's been disabled in the configuration
            if (!config.extractXSecs.crossSectionIsEnabled.at(selectionName).at(name))
                continue;

            // Add the cross-section object to the map using the binning from the input configuration
            const auto &[extendedBinEdges, hasUnderflow, hasOverflow] = CrossSectionHelper::GetExtendedBinEdges(binning.min, binning.max, binning.binEdges);
            xsecMap[selectionName].emplace(name, CrossSectionHelper::CrossSection(systParams, extendedBinEdges, hasUnderflow, hasOverflow, scaleByBinWidth));
        }
    }

    // ATTN here we use the machinary for a differential cross-section, and treat the total cross-section as a single-bin measurement.
    // The "kinematic quantity" in this case is just a dummy parameter. Here we define a single bin with edges arbitrarily chosen to be
    // (-1 -> +1), and we request that the cross-section object does not apply bin-width scaling. When we fill this object, we will use
    // the dummy kinematic quantity with a value of 0 for all events. This is arbitrary, as long as it's within the bin edges we chose.
    // In this way the single bin contains all events. It's just a trick to avoid implementing extra logic for the total cross-section.
    for (const auto &selectionName : {"generic", "golden"})
    {
        // Don't setup a cross-section object if it's been disabled in the configuration
        if (!config.extractXSecs.crossSectionIsEnabled.at(selectionName).at("total"))
            continue;

        xsecMap[selectionName].emplace("total", CrossSectionHelper::CrossSection(systParams, {-1.f, 1.f}, false, false, false));
    }

    // The dummy value that will be used as the "kinematic quantity" for the total cross-section
    const auto dummyValue = 0.f;

    // Check to see if we have any cross-sections enabled
    if (xsecMap.empty())
    {
        std::cout << "All cross-sections have been disabled in the configuration! Nothing more to do" << std::endl;
        return;
    }

    // Print the names of the cross-sections we are going to extract
    std::cout << "The following cross-sections are enabled:" << std::endl;
    for (const auto &[selectionName, xsecs] : xsecMap)
    {
        std::cout << "  - Selection: " << selectionName << std::endl;
        for (const auto &entry : xsecs)
        {
            const auto &name = entry.first;
            std::cout << "    - " << name << std::endl;
        }
    }

    // -------------------------------------------------------------------------------------------------------------------------------------
    // Setup the relevent "getters" for each cross-section
    // -------------------------------------------------------------------------------------------------------------------------------------
    // Here we define a map from the name of each cross-section to a function which pulls out the relevant kinematic quanitity from an input
    // analysis data object. Again this is done up-front to reduce code-bloat below.
    std::unordered_map< std::string, std::function<float(const AnalysisHelper::AnalysisData &)> > getValue;

    // Differential cross-section kinematic parameters
    getValue.emplace("muonCosTheta",  [](const auto &data) { return data.muonCosTheta;  });
    getValue.emplace("muonPhi",       [](const auto &data) { return data.muonPhi;       });
    getValue.emplace("muonMomentum",  [](const auto &data) { return data.muonMomentum;  });

    getValue.emplace("pionCosTheta",  [](const auto &data) { return data.pionCosTheta;  });
    getValue.emplace("pionPhi",       [](const auto &data) { return data.pionPhi;       });
    getValue.emplace("pionMomentum",  [](const auto &data) { return data.pionMomentum;  });

    getValue.emplace("muonPionAngle", [](const auto &data) { return data.muonPionAngle; });
    getValue.emplace("nProtons",      [](const auto &data) { return data.nProtons;      });

    // ATTN as described above, for the total cross-section we don't have an associated kinematic quantity so we just return a dummy value
    getValue.emplace("total", [=](const auto &) { return dummyValue; });

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
        // std::cout<<"DEBUG Point 6"<<std::endl;

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

                // Check the value of the kinematic quantities are within the phase-space limits
                for (const auto &[name, minMax] : phaseSpaceMap)
                {
                    const auto &[min, max] = minMax;
                    const auto value = getValue.at(name)(recoData);

                    if (value < min || value > max)
                    {
                        passesPhaseSpaceReco = false;
                        break;
                    }
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
            // std::cout<<"DEBUG Point 7"<<std::endl;

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

                // Check the value of the kinematic quantities are within the phase-space limits
                for (const auto &[name, minMax] : phaseSpaceMap)
                {
                    const auto &[min, max] = minMax;
                    const auto value = getValue.at(name)(truthData);

                    if (value < min || value > max)
                    {
                        passesPhaseSpaceTruth = false;
                        break;
                    }
                }
            }

            // std::cout<<"DEBUG Point 8"<<std::endl;

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
            // std::cout<<"DEBUG Point 9"<<std::endl;

            // -----------------------------------------------------------------------------------------------------------------------------
            // Handle all other events (i.e those from the nominal simulation): Overlays, dirt, EXT data
            // -----------------------------------------------------------------------------------------------------------------------------

            // Get the flux weights
            const auto fluxWeights = (
                isOverlay
                    ? CrossSectionHelper::GetWeightsMap(pEvent->truth, systParams.fluxDimensions, config.extractXSecs.mutuallyExclusiveDimensions)
                    : CrossSectionHelper::GetUnitWeightsMap(systParams.fluxDimensions)
            );

            // Fill the flux-reweightor with all overlay events from a neutrinos of the desired flavour in the active volume
            if (isOverlay && AnalysisHelper::IsInActiveVolume(pEvent->truth.nuVertex()) &&
                std::find(config.fluxRun1.nuPdgsSignal.begin(), config.fluxRun1.nuPdgsSignal.end(), pEvent->truth.nuPdgCode()) != config.fluxRun1.nuPdgsSignal.end())
            {
                scalingData.pFluxReweightor->AddEvent(pEvent->truth.nuEnergy(), weight, fluxWeights);
            }

            // Get the cross-section weights
            // ATTN here we optionally scale the cross-section weights down by the genieTuneEventWeight - this is done so we don't
            // double count this weight (once in the nominal event weight, and once in the xsec systematic event weights)
            auto xsecWeightsScaleFactor = (isOverlay && config.extractXSecs.scaleXSecWeights) ? pEvent->truth.genieTuneEventWeight() : 1.f;
            xsecWeightsScaleFactor = std::max(xsecWeightsScaleFactor, 0.0001f);
            const auto xsecWeights = (
                isOverlay
                    ? CrossSectionHelper::ScaleWeightsMap(CrossSectionHelper::GetWeightsMap(pEvent->truth, systParams.xsecDimensions, config.extractXSecs.mutuallyExclusiveDimensions), xsecWeightsScaleFactor)
                    : CrossSectionHelper::GetUnitWeightsMap(systParams.xsecDimensions)
            );

            // Get the reinteraction weights
            const auto reintWeights = (
                isOverlay
                    ? CrossSectionHelper::GetWeightsMap(pEvent->truth, systParams.reintDimensions, config.extractXSecs.mutuallyExclusiveDimensions)
                    : CrossSectionHelper::GetUnitWeightsMap(systParams.reintDimensions)
            );

            // std::cout<<"DEBUG Point 10"<<std::endl;
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

                        xsec.AddSignalEvent(recoValue, trueValue, isSelected, weight, fluxWeights, xsecWeights, reintWeights, "Test");
                    }
                }
            }
            // Handle selected background events
            else
            {
                // std::cout<<"DEBUG Point 11"<<std::endl;
                for (auto &[selectionName, xsecs] : xsecMap)
                {
                    // Only use selected background events
                    const auto isSelected = (selectionName == "golden" ? isSelectedGolden : isSelectedGeneric);
                    if (!isSelected)
                        continue;

                    for (auto &[name, xsec] : xsecs)
                    {
                        const auto recoValue = getValue.at(name)(recoData);
                        std::vector<float> sidebandWeights;
                        xsec.AddSelectedBackgroundEvent(recoValue, isDirt, weight, fluxWeights, xsecWeights, reintWeights, sidebandWeights, "Test");
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
            // Get the event rates for BNB data, backgrounds, and signal
            // -----------------------------------------------------------------------------------------------------------------------------
            const auto selectedEventsData = xsec.GetSelectedBNBDataEvents();
            std::cout << "Selected BNB data events" << std::endl;
            FormattingHelper::SaveMatrix(selectedEventsData, "xsec_" + selectionName + "_" + name + "_data_selected_eventRate.txt");

            const auto selectedEventsBackground = xsec.GetSelectedBackgroundEvents();
            std::cout << "Selected background events" << std::endl;
            FormattingHelper::SaveMatrix(selectedEventsBackground, "xsec_" + selectionName + "_" + name + "_background_selected_eventRate.txt");

            const auto selectedEventsSignal = xsec.GetSelectedSignalEvents();
            std::cout << "Selected signal events" << std::endl;
            FormattingHelper::SaveMatrix(selectedEventsSignal, "xsec_" + selectionName + "_" + name + "_signal_selected_eventRate.txt");

            const auto allEventsSignal = xsec.GetSignalEvents();
            std::cout << "All signal events" << std::endl;
            FormattingHelper::SaveMatrix(allEventsSignal, "xsec_" + selectionName + "_" + name + "_signal_all_eventRate.txt");

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