/**
 *  @file  ubcc1pi_standalone/Macros/ExtractCC0PiNormalisation.cxx
 *
 *  @brief The implementation file of the ExtractCC0PiNormalisation macro
 */

#include "ubcc1pi_standalone/Macros/Macros.h"

#include "ubcc1pi_standalone/Objects/FileReader.h"

#include "ubcc1pi_standalone/Helpers/PlottingHelper.h"
#include "ubcc1pi_standalone/Helpers/AnalysisHelper.h"
#include "ubcc1pi_standalone/Helpers/SelectionHelper.h"
#include "ubcc1pi_standalone/Helpers/NormalisationHelper.h"
#include "ubcc1pi_standalone/Helpers/CrossSectionHelper.h"
#include "ubcc1pi_standalone/Helpers/FormattingHelper.h"

#include <TH2F.h>
#include <TMath.h>

using namespace ubcc1pi;

namespace ubcc1pi_macros
{

void ExtractCC0PiNormalisation(const Config &config)
{
    // std::cout<<"DEBUG - ExtractCC0PiNormalisation - Point 0"<<std::endl;
    // Setup the input files
    // std::vector< std::tuple<AnalysisHelper::SampleType, std::string, float> > inputData;
    // inputData.emplace_back(AnalysisHelper::Overlay, config.files.overlaysFileName, NormalisationHelper::GetOverlaysNormalisation(config));
    // inputData.emplace_back(AnalysisHelper::Dirt,    config.files.dirtFileName,     NormalisationHelper::GetDirtNormalisation(config));
    // inputData.emplace_back(AnalysisHelper::DataEXT, config.files.dataEXTFileName,  NormalisationHelper::GetDataEXTNormalisation(config));
    // inputData.emplace_back(AnalysisHelper::DataBNB, config.files.dataBNBFileName,  1.f);

    std::vector< std::tuple<AnalysisHelper::SampleType, std::string, std::string, float> > inputData;
    inputData.emplace_back(AnalysisHelper::Overlay, "", config.files.overlaysFileName, NormalisationHelper::GetOverlaysNormalisation(config));
    inputData.emplace_back(AnalysisHelper::Dirt,    "", config.files.dirtFileName, NormalisationHelper::GetDirtNormalisation(config));
    inputData.emplace_back(AnalysisHelper::DataEXT, "", config.files.dataEXTFileName, NormalisationHelper::GetDataEXTNormalisation(config));
    inputData.emplace_back(AnalysisHelper::DataBNB, "", config.files.dataBNBFileName, 1.f);

    // std::cout<<"DEBUG - ExtractCC0PiNormalisation - Point 1"<<std::endl;
    // Add the detector variation files
    for (const auto &[name, fileName] : config.files.detVarFiles)
        inputData.emplace_back(AnalysisHelper::DetectorVariation, name, fileName, NormalisationHelper::GetDetectorVariationNormalisation(config, name));

    // std::cout<<"DEBUG - ExtractCC0PiNormalisation - Point 2"<<std::endl;
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
    // std::cout<<"DEBUG - ExtractCC0PiNormalisation - Point 3"<<std::endl;
    // Read in the flux
    const auto fluxHistNames = CrossSectionHelper::GetNominalFluxHistNames(config.flux.nuPdgsSignal, config.flux.nuPdgToHistName, config.flux.nomHistPattern);
    const auto &[fluxBinEdges, fluxValues] = CrossSectionHelper::ReadNominalFlux(config.flux.fileName, fluxHistNames, config.flux.pot);
    scalingData.pFluxReweightor = std::make_shared<CrossSectionHelper::FluxReweightor>(fluxBinEdges, fluxValues, systParams.fluxDimensions);

    scalingData.exposurePOT = config.norms.dataBNBTor875WCut / (1e20);
    scalingData.nTargets = config.global.targetDensity * (1e-8) * AnalysisHelper::GetFiducialVolume();

    std::cout << "- Flux:            " << scalingData.pFluxReweightor->GetIntegratedNominalFlux() << " * 10^-10 cm^-2 POT^-1" << std::endl;
    std::cout << "- Exposure:        " << scalingData.exposurePOT << " * 10^20 POT" << std::endl;
    std::cout << "- Target density:  " << config.global.targetDensity << " * 10^31 nucleons/cm^3" << std::endl;
    std::cout << "- Fiducial volume: " << AnalysisHelper::GetFiducialVolume() << " * cm^3" << std::endl;
    std::cout << "- nTargets:        " << scalingData.nTargets << " * 10^31 nucleons" << std::endl;

    // Get the selections
    auto selection = SelectionHelper::GetSelection("CC0pi");

    // std::cout<<"DEBUG - ExtractCC0PiNormalisation - Point 4"<<std::endl;
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

    auto selectionName = "cc0pi";
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
        // std::cout<<"DEBUG - ExtractCC0PiNormalisation - Point 5"<<std::endl;
        // Add to the phase-space map
        phaseSpaceMap.emplace(name, std::pair<float, float>({binning.min, binning.max}));

        // Here we calculate every cross-section using both the generic and golden selection. In the end we only use the golden selection for
        // the pion momentum, but we additionally apply it to other cross-sections as a cross-check.
        // for (const auto &selectionName : {"generic", "golden"})
        // {
        //     // Don't setup a cross-section object if it's been disabled in the configuration
        //     if (!config.extractXSecs.crossSectionIsEnabled.at(selectionName).at(name))
        //         continue;

            // Add the cross-section object to the map using the binning from the input configuration
            const auto &[extendedBinEdges, hasUnderflow, hasOverflow] = CrossSectionHelper::GetExtendedBinEdges(binning.min, binning.max, binning.binEdges);
            xsecMap[selectionName].emplace(name, CrossSectionHelper::CrossSection(systParams, extendedBinEdges, hasUnderflow, hasOverflow, scaleByBinWidth));
        // }
    }

    // std::cout<<"DEBUG - ExtractCC0PiNormalisation - Point 6"<<std::endl;
    // ATTN here we use the machinary for a differential cross-section, and treat the total cross-section as a single-bin measurement.
    // The "kinematic quantity" in this case is just a dummy parameter. Here we define a single bin with edges arbitrarily chosen to be
    // (-1 -> +1), and we request that the cross-section object does not apply bin-width scaling. When we fill this object, we will use
    // the dummy kinematic quantity with a value of 0 for all events. This is arbitrary, as long as it's within the bin edges we chose.
    // In this way the single bin contains all events. It's just a trick to avoid implementing extra logic for the total cross-section.

    // // Don't setup a cross-section object if it's been disabled in the configuration
    // if (!config.extractXSecs.crossSectionIsEnabled.at(selectionName).at("total"))
    //     continue;
    xsecMap[selectionName].emplace("total", CrossSectionHelper::CrossSection(systParams, {-1.f, 1.f}, false, false, false));

    // The dummy value that will be used as the "kinematic quantity" for the total cross-section
    const auto dummyValue = 0.f;

    // Check to see if we have any cross-sections enabled
    if (xsecMap.empty())
    {
        std::cout << "All cross-sections have been disabled in the configuration! Nothing more to do" << std::endl;
        return;
    }

    // std::cout<<"DEBUG - ExtractCC0PiNormalisation - Point 7"<<std::endl;
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

    // std::cout<<"DEBUG - ExtractCC0PiNormalisation - Point 8"<<std::endl;
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

    getValue.emplace("pionCosTheta",  [](const auto &data) { return data.protonCosTheta;  }); // Getting proton instead of pion values
    getValue.emplace("pionPhi",       [](const auto &data) { return data.protonPhi;       }); // Getting proton instead of pion values
    getValue.emplace("pionMomentum",  [](const auto &data) { return data.protonMomentum;  }); // Getting proton instead of pion values

    getValue.emplace("muonPionAngle", [](const auto &data) { return data.muonProtonAngle; }); // Getting proton instead of pion values
    getValue.emplace("nProtons",      [](const auto &data) { return data.nProtons-1;      }); //Leading proton treated as pion in CC0pi analysis

    // ATTN as described above, for the total cross-section we don't have an associated kinematic quantity so we just return a dummy value
    getValue.emplace("total", [=](const auto &) { return dummyValue; });

    // std::cout<<"DEBUG - ExtractCC0PiNormalisation - Point 9"<<std::endl;
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

        // std::cout<<"DEBUG - ExtractCC0PiNormalisation - Point 10"<<std::endl;
        // Open the input file for reading and enable the branches with systematic event weights (if required)
        FileReader reader(fileName);

        // std::cout<<"DEBUG - ExtractCC0PiNormalisation - Point 11"<<std::endl;
        // if (isOverlay)
        //     reader.EnableSystematicBranches();

        if (isOverlay)
            reader.EnableSystematicBranches();

        auto pEvent = reader.GetBoundEventAddress();

        // Loop over the events in the file
        const auto nEvents = reader.GetNumberOfEvents();
        for (unsigned int i = 0; i < nEvents; ++i)
        {
            // std::cout<<"DEBUG - ExtractCC0PiNormalisation - Point 12"<<std::endl;
            AnalysisHelper::PrintLoadingBar(i, nEvents);
            // std::cout<<"DEBUG - ExtractCC0PiNormalisation - Point 13"<<std::endl;
            reader.LoadEvent(i);
            // std::cout<<"DEBUG - ExtractCC0PiNormalisation - Point 14"<<std::endl;
            // -----------------------------------------------------------------------------------------------------------------------------
            // Work out if this event passed the selection and apply any additional phase-space cuts based on the input binning
            // -----------------------------------------------------------------------------------------------------------------------------

            // Run the selection
            // const auto &[passedGoldenSelection, cutsPassed, assignedPdgCodes] = selection.Execute(pEvent);
            const auto &[passedCC0piSelection, cutsPassed, assignedPdgCodes] = selection.Execute(pEvent);
            // const auto passedGenericSelection = SelectionHelper::IsCutPassed(cutsPassed, config.global.lastCutGeneric);

            // std::cout<<"DEBUG - ExtractCC0PiNormalisation - Point 14.1 - passedCC0piSelection: "<<passedCC0piSelection<<std::endl;
            // Get the reco analysis data (if available, otherwise set to dummy values)
            const auto recoData = (
                passedCC0piSelection
                    ? AnalysisHelper::GetRecoAnalysisDataCC0Pi(pEvent->reco, assignedPdgCodes)
                    : AnalysisHelper::GetDummyAnalysisData()
            );

            // std::cout<<"DEBUG - ExtractCC0PiNormalisation - Point 15"<<std::endl;
            // Here we apply reco-level phase-space restrictions
            // For any event that passes the generic selection, get the value of the kinematic quantity and check if it is outside of the
            // min/max values supplied in the binning. If so, then reject the event.
            bool passesPhaseSpaceReco = false;
            if (passedCC0piSelection)
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
                        // std::cout<<"FAIL1!!!!!! DEBUG passesPhaseSpaceReco - name: "<<name<<" - value: "<<value<<" - min: "<<min<<" - max: "<<max<<std::endl;
                        passesPhaseSpaceReco = false;
                        break;
                    }
                }
            }
            // if(passesPhaseSpaceReco || passedCC0piSelection)
            //     std::cout<<"passedCC0piSelection: "<<passedCC0piSelection<<" - passesPhaseSpaceReco: "<<passesPhaseSpaceReco<<std::endl;

            const auto isSelected = passedCC0piSelection && passesPhaseSpaceReco;

            // std::cout<<"DEBUG - ExtractCC0PiNormalisation - Point 16"<<std::endl;
            // -----------------------------------------------------------------------------------------------------------------------------
            // Handle BNB data
            // -----------------------------------------------------------------------------------------------------------------------------
            if (isDataBNB && isSelected)
            {
                for (auto &[selectionName, xsecs] : xsecMap)
                {
                    for (auto &[name, xsec] : xsecs)
                    {
                        xsec.AddSelectedBNBDataEvent(getValue.at(name)(recoData));
                        // std::cout<<"Added event to xSec - Point 0"<<std::endl;
                    }
                }
                // For BNB data that's all we need to do!
                continue;
            }

            // -----------------------------------------------------------------------------------------------------------------------------
            // Work out if this event is signal, and apply any phase-space restrictions based on the input binning
            // -----------------------------------------------------------------------------------------------------------------------------

            // std::cout<<"DEBUG - ExtractCC0PiNormalisation - Point 17"<<std::endl;
            // Get the nominal event weight, scaled by the sample normalisation
            const auto weight = AnalysisHelper::GetNominalEventWeight(pEvent) * normalisation;

            // Determine if this is truly a CC1Pi event
            const auto IsTrueCC0Pi = (isOverlay || isDetVar) && AnalysisHelper::IsTrueCC0Pi(pEvent, config.global.useAbsPdg, config.global.protonMomentumThreshold);

            // std::cout<<"DEBUG - ExtractCC0PiNormalisation - Point 17.1"<<std::endl;
            // Get the truth analysis data (if available, otherwise set to dummy values)
            const auto truthData = (
                IsTrueCC0Pi
                    ? AnalysisHelper::GetTruthAnalysisDataCC0Pi(pEvent->truth, config.global.useAbsPdg, config.global.protonMomentumThreshold)
                    : AnalysisHelper::GetDummyAnalysisData()
            );

            // std::cout<<"DEBUG - ExtractCC0PiNormalisation - Point 17.2"<<std::endl;
            // Here we apply truth-level phase-space restrictions
            // For all true CC1Pi events, we check if the values of each kinematic variable are within the supplied limits. If not then the
            // event is not classed as "signal"
            bool passesPhaseSpaceTruth = false;
            if (IsTrueCC0Pi)
            {
                // Start by assuming the event passes the phase-space cuts
                passesPhaseSpaceTruth = true;
                // std::cout<<"DEBUG - ExtractCC0PiNormalisation - Point 17.3"<<std::endl;
                // Check the value of the kinematic quantities are within the phase-space limits
                for (const auto &[name, minMax] : phaseSpaceMap)
                {
                    // std::cout<<"DEBUG - ExtractCC0PiNormalisation - Point 17.4"<<std::endl;
                    const auto &[min, max] = minMax;
                    const auto value = getValue.at(name)(truthData);

                    if (value < min || value > max)
                    {
                        // std::cout<<"FAIL2!!!!!! DEBUG passesPhaseSpaceTruth - name: "<<name<<" - value: "<<value<<" - min: "<<min<<" - max: "<<max<<std::endl;
                        passesPhaseSpaceTruth = false;
                        break;
                    }
                }
            }
            // if(passesPhaseSpaceTruth)
            //     std::cout<<"passesPhaseSpaceTruth: "<<passesPhaseSpaceTruth<<std::endl;

            const auto isSignal = IsTrueCC0Pi && passesPhaseSpaceTruth;

            // std::cout<<"DEBUG - ExtractCC0PiNormalisation - Point 18"<<std::endl;
            // -----------------------------------------------------------------------------------------------------------------------------
            // Handle the detector variation samples (unisims)
            // -----------------------------------------------------------------------------------------------------------------------------
            if (isDetVar)
            {
                // Handle signal events
                if (isSelected)
                {
                    // std::cout<<"passedCC0piSelection: "<<passedCC0piSelection<<std::endl;
                    // std::cout<<"isSelected: "<<isSelected<<std::endl;
                    // std::cout<<"IsTrueCC0Pi: "<<IsTrueCC0Pi<<" isOverlay: "<<isOverlay<<" - isDetVar: "<<isDetVar<<std::endl;
                    // std::cout<<"isSignal: "<<isSignal<<std::endl;
                    if (isSignal)
                    {
                        for (auto &[selectionName, xsecs] : xsecMap)
                        {
                            // std::cout<<"DEBUG event adding - Point 0.1 - selectionName: "<<selectionName<<std::endl;
                            for (auto &[name, xsec] : xsecs)
                            {
                                // std::cout<<"DEBUG event adding - Point 0.2 - selectionName: "<<selectionName<<std::endl;
                                const auto recoValue = getValue.at(name)(recoData);
                                const auto trueValue = getValue.at(name)(truthData);
                                xsec.AddSignalEventDetVar(recoValue, trueValue, isSelected, weight, sampleName);
                                // std::cout<<"Added event to xSec - Point 1"<<std::endl;
                            }
                        }
                    }
                    // Handle selected background events
                    else
                    {
                        for (auto &[selectionName, xsecs] : xsecMap)
                        {
                            // std::cout<<"DEBUG event adding - Point 1.2 - selectionName: "<<selectionName<<std::endl;
                            for (auto &[name, xsec] : xsecs)
                            {
                                // std::cout<<"DEBUG event adding - Point 1.3 - name: "<<name<<std::endl;
                                const auto recoValue = getValue.at(name)(recoData);
                                xsec.AddSelectedBackgroundEventDetVar(recoValue, weight, sampleName);
                                // std::cout<<"Added event to xSec - Point 2"<<std::endl;
                            }
                        }
                    }
                }

                // For detector variation samples, that's all we need to do!
                continue;
            }

            // std::cout<<"DEBUG - ExtractCC0PiNormalisation - Point 19_000"<<std::endl;
            // -----------------------------------------------------------------------------------------------------------------------------
            // Handle all other events (i.e those from the nominal simulation): Overlays, dirt, EXT data
            // -----------------------------------------------------------------------------------------------------------------------------

            // Get the flux weights
            const auto fluxWeights = (
                isOverlay
                    ? CrossSectionHelper::GetWeightsMap(pEvent->truth, systParams.fluxDimensions, config.extractXSecs.mutuallyExclusiveDimensions)
                    : CrossSectionHelper::GetUnitWeightsMap(systParams.fluxDimensions)
            );

            // std::cout<<"DEBUG - ExtractCC0PiNormalisation - Point 19.1"<<std::endl;

            // std::cout<<"DEBUG - ExtractCC0PiNormalisation - Point 19.2"<<std::endl;
            // Fill the flux-reweightor with all overlay events from neutrinos of the desired flavour in the active volume
            if (isOverlay && AnalysisHelper::IsInActiveVolume(pEvent->truth.nuVertex()) &&
                std::find(config.flux.nuPdgsSignal.begin(), config.flux.nuPdgsSignal.end(), pEvent->truth.nuPdgCode()) != config.flux.nuPdgsSignal.end())
            {
                scalingData.pFluxReweightor->AddEvent(pEvent->truth.nuEnergy(), weight, fluxWeights);
            }

            // std::cout<<"DEBUG - ExtractCC0PiNormalisation - Point 19.3"<<std::endl;
            // Get the cross-section weights
            // ATTN here we optionally scale the cross-section weights down by the genieTuneEventWeight - this is done so we don't
            // double count this weight (once in the nominal event weight, and once in the xsec systematic event weights)

            const auto xsecWeightsScaleFactor = 1.f;
            // const auto xsecWeightsScaleFactor = (isOverlay && config.extractXSecs.scaleXSecWeights) ? pEvent->truth.genieTuneEventWeight() : 1.f;

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

            // std::cout<<"DEBUG - ExtractCC0PiNormalisation - Point 20"<<std::endl;
            // Handle signal events
            if (isSignal)
            {
                // std::cout<<"DEBUG event adding - Point 2.1 - isSignal: "<<isSignal<<std::endl;
                // std::cout<<"passedCC0piSelection: "<<passedCC0piSelection<<std::endl;
                // std::cout<<"isSelected: "<<isSelected<<std::endl;
                // std::cout<<"IsTrueCC0Pi: "<<IsTrueCC0Pi<<" isOverlay: "<<isOverlay<<" - isDetVar: "<<isDetVar<<std::endl;
                // std::cout<<"isSignal: "<<isSignal<<std::endl;
                for (auto &[selectionName, xsecs] : xsecMap)
                {
                    // std::cout<<"DEBUG event adding - Point 2.2 - selectionName: "<<selectionName<<std::endl;
                    // Determine if we passed the relevant selection

                    for (auto &[name, xsec] : xsecs)
                    {
                        // std::cout<<"DEBUG event adding - Point 2.3 - name: "<<name<<std::endl;
                        const auto recoValue = getValue.at(name)(recoData);
                        const auto trueValue = getValue.at(name)(truthData);

                        xsec.AddSignalEvent(recoValue, trueValue, isSelected, false, weight, fluxWeights, xsecWeights, reintWeights);
                        // std::cout<<"Added event to xSec - Point 3"<<std::endl;
                    }
                }
            }
            // Handle selected background events
            else
            {
                // Only use selected background events
                if(isSelected)
                {
                    // std::cout<<"DEBUG event adding - Point 3.1 - isSelected: "<<isSelected<<std::endl;
                    // std::cout<<"passedCC0piSelection: "<<passedCC0piSelection<<std::endl;
                    // std::cout<<"isSelected: "<<isSelected<<std::endl;
                    // std::cout<<"IsTrueCC0Pi: "<<IsTrueCC0Pi<<" isOverlay: "<<isOverlay<<" - isDetVar: "<<isDetVar<<std::endl;
                    // std::cout<<"isSignal: "<<isSignal<<std::endl;
                    for (auto &[selectionName, xsecs] : xsecMap)
                    {
                        // std::cout<<"DEBUG event adding - Point 3.2 - selectionName: "<<selectionName<<std::endl;
                        for (auto &[name, xsec] : xsecs)
                        {
                            // std::cout<<"DEBUG event adding - Point 3.3 - name: "<<name<<std::endl;
                            const auto recoValue = getValue.at(name)(recoData);
                            xsec.AddSelectedBackgroundEvent(recoValue, isDirt, false, weight, fluxWeights, xsecWeights, reintWeights);
                            // std::cout<<"Added event to xSec - Point 4"<<std::endl;
                        }
                    }
                }
            }
        }
    }

    // -------------------------------------------------------------------------------------------------------------------------------------
    // Save CC0pi counts for each cross-section
    // -------------------------------------------------------------------------------------------------------------------------------------

    // std::cout<<"DEBUG - ExtractCC0PiNormalisation - Point 21"<<std::endl;
    // Loop over all cross-section objects
    // TFile *fout = TFile::Open("TEST.root", "RECREATE");
    for (const auto &[selectionName, xsecs] : xsecMap)
    {
        for (const auto &[name, xsec] : xsecs)
        {
            std::cout << "Processing CC0Pi selection: " << name << std::endl;

            // -----------------------------------------------------------------------------------------------------------------------------
            // Get the event rates for BNB data, backgrounds, and signal
            // -----------------------------------------------------------------------------------------------------------------------------
            const auto selectedEventsData = xsec.GetSelectedBNBDataEvents();
            std::cout << "Selected BNB data events" << std::endl;
            FormattingHelper::SaveMatrix(selectedEventsData, "CC0pi_" + name + "_data_selected_eventRate.txt");


            const auto selectedEventsBackground = xsec.GetSelectedBackgroundEvents();
            std::cout << "Selected background events" << std::endl;
            FormattingHelper::SaveMatrix(selectedEventsBackground, "CC0pi_" + name + "_background_selected_eventRate.txt");

            const auto selectedEventsSignal = xsec.GetSelectedSignalEvents();
            std::cout << "Selected signal events" << std::endl;
            FormattingHelper::SaveMatrix(selectedEventsSignal, "CC0pi_" + name + "_signal_selected_eventRate.txt");

            const auto allEventsSignal = xsec.GetSignalEvents();
            std::cout << "All signal events" << std::endl;
            FormattingHelper::SaveMatrix(allEventsSignal, "CC0pi_" + name + "_signal_all_eventRate.txt");


            // std::vector<float> xx; // Fill the vector
            // for (int k = 0; k<3; k++) xx.push_back(k);
            // fout->WriteObject(&xx, (name+"_True").c_str());
            // fout->WriteObject(&xx, (name+"_Reco").c_str());
            // fout->WriteObject(&xx, (name+"_Ratio").c_str());
        }
    }
    // fout->Close();

}

} // namespace ubcc1pi_macros
