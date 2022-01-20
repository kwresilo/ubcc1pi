/**
 *  @file  ubcc1pi_standalone/Macros/ExtractSidebandFit.cxx
 *
 *  @brief The implementation file of the ExtractSidebandFit macro
 */

#include "ubcc1pi_standalone/Macros/Macros.h"
#include "ubcc1pi_standalone/Objects/FileReader.h"

#include "ubcc1pi_standalone/Helpers/AnalysisHelper.h"
#include "ubcc1pi_standalone/Helpers/NormalisationHelper.h"
#include "ubcc1pi_standalone/Helpers/SelectionHelper.h"
#include "ubcc1pi_standalone/Helpers/CrossSectionHelper.h"
#include "ubcc1pi_standalone/Helpers/FormattingHelper.h"
#include "ubcc1pi_standalone/Helpers/FittingHelper.h"
#include "ubsmear.h"

#include <fstream> // Todo: not use txt files

// Boost libraries
// #include "binary_iarchive.hpp"
#include "binary_oarchive.hpp"
#include "binary_object.hpp"
#include "map.hpp"
#include "vector.hpp"

using namespace ubcc1pi;

std::vector<float> x, y, errorY, S;

//______________________________________________________________________________

void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
    //calculate chi2
    if(x.size()!=y.size() || y.size()!=errorY.size() || x.size()*x.size()!=S.size())
        throw std::logic_error("Fitting function for ExtractSidebandFit - Incompatible input dimenstions.");

    Double_t chisq = 0;
    Double_t delta;
    Int_t nBins = x.size();
    
    auto xScaled = x; 
    for (Int_t i=0; i<nBins; i++)
    {
        xScaled[i]*=par[i];
    }
    
    std::vector<float> xSmeared(nBins, 0);

    for (Int_t i=0; i<nBins; i++)
    {
        for (Int_t j=0; j<nBins; j++)
        {
            xSmeared[i] +=  S[i+j*nBins]*xScaled[j];
        }
    }
    // const auto xSmeared = S*xScaled;

    for (Int_t i=0; i<nBins; i++) 
    {
        delta  = (y[i]-xSmeared[i])/errorY[i];
        chisq += delta*delta;
    }
    f = chisq;
}

namespace ubcc1pi_macros
{

void ExtractSidebandFit(const Config &config)
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
    // auto selection = SelectionHelper::GetDefaultSelection();
    auto sidebandSelection = SelectionHelper::GetSelection("CC0pi");

    // -------------------------------------------------------------------------------------------------------------------------------------
    // Setup the cross-section objects
    // -------------------------------------------------------------------------------------------------------------------------------------
    // Here we make a map from a name of the cross-section to the cross-section object itself. In this way, we can iterate through the
    // cross-section objects and reduce code-bloat. The first index is an identifier for the selection that's applied (generic or goldlen),
    // the second index is an identifier for the kinematic quantity that's relevant for the cross-section (e.g. muonMomentum), and the
    // mapped type is the cross-section object.
    // std::map<std::string, std::map<std::string, CrossSectionHelper::CrossSection> > xsecMap;
    std::map<std::string, std::map<std::string, CrossSectionHelper::CrossSection> > xsecMapSideband;
    // std::map<std::string, std::map<std::string, CrossSectionHelper::CrossSection> > xsecMapSideband2;

    // We additionally make a map from each cross-section to the limits of the phase-space that we should consider. The key is the
    // identifier for the kinematic quantity and the mapped value is a pair containing the limits [min, max]
    std::map< std::string, std::pair<float, float> > phaseSpaceMap;

    // ATTN the configuration allows the user to enable or disable each cross-section. If a cross-section has been disabled, then it won't
    // be added the xSecMap. However, the phaseSpaceMap always includes all kinematic parameters!

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

        // { "muonPionAngle", config.global.muonPionAngle, true  },
        { "nProtons",      config.global.nProtons,      false }

    })
    {
        // Add to the phase-space map
        phaseSpaceMap.emplace(name, std::pair<float, float>({binning.min, binning.max}));

        // Here we calculate every cross-section using both the generic and golden selection. In the end we only use the golden selection for
        // the pion momentum, but we additionally apply it to other cross-sections as a cross-check.
        // Add the cross-section object to the map using the binning from the input configuration
        const auto &[extendedBinEdges, hasUnderflow, hasOverflow] = CrossSectionHelper::GetExtendedBinEdges(binning.min, binning.max, binning.binEdges);
        for (const auto &selectionName : {"generic", "golden"})
        {
            // Don't setup a cross-section object if it's been disabled in the configuration
            if (config.extractXSecs.crossSectionIsEnabled.at(selectionName).at(name))
            {
                // xsecMap[selectionName].emplace(name, CrossSectionHelper::CrossSection(systParams, extendedBinEdges, hasUnderflow, hasOverflow, scaleByBinWidth));
                xsecMapSideband[selectionName].emplace(name, CrossSectionHelper::CrossSection(systParams, extendedBinEdges, hasUnderflow, hasOverflow, scaleByBinWidth));
                // xsecMapSideband2[selectionName].emplace(name, CrossSectionHelper::CrossSection(systParams, extendedBinEdges, hasUnderflow, hasOverflow, scaleByBinWidth));
            }
        }
        // if (config.extractXSecs.crossSectionIsEnabled.at("sideband").at(name))
        //     xsecMapSideband["sideband"].emplace(name, CrossSectionHelper::CrossSection(systParams, extendedBinEdges, hasUnderflow, hasOverflow, scaleByBinWidth));

    }

    // ATTN here we use the machinary for a differential cross-section, and treat the total cross-section as a single-bin measurement.
    // The "kinematic quantity" in this case is just a dummy parameter. Here we define a single bin with edges arbitrarily chosen to be
    // (-1 -> +1), and we request that the cross-section object does not apply bin-width scaling. When we fill this object, we will use
    // the dummy kinematic quantity with a value of 0 for all events. This is arbitrary, as long as it's within the bin edges we chose.
    // In this way the single bin contains all events. It's just a trick to avoid implementing extra logic for the total cross-section.
    for (const auto &selectionName : {"generic", "golden"})
    {
        // Don't setup a cross-section object if it's been disabled in the configuration
        if (config.extractXSecs.crossSectionIsEnabled.at(selectionName).at("total"))
        {
            // xsecMap[selectionName].emplace("total", CrossSectionHelper::CrossSection(systParams, {-1.f, 1.f}, false, false, false));
            xsecMapSideband[selectionName].emplace("total", CrossSectionHelper::CrossSection(systParams, {-1.f, 1.f}, false, false, false));
            // xsecMapSideband2[selectionName].emplace("total", CrossSectionHelper::CrossSection(systParams, {-1.f, 1.f}, false, false, false));
        }

    }
    // Don't setup a cross-section object if it's been disabled in the configuration
    // if (config.extractXSecs.crossSectionIsEnabled.at("sideband").at("total"))
    //     xsecMapSideband["sideband"].emplace("total", CrossSectionHelper::CrossSection(systParams, {-1.f, 1.f}, false, false, false));

    // The dummy value that will be used as the "kinematic quantity" for the total cross-section
    const auto dummyValue = 0.f;

    // // Check to see if we have any cross-sections enabled
    // if (xsecMap.empty())
    // {
    //     std::cout << "All cross-sections have been disabled in the configuration! Nothing more to do" << std::endl;
    //     return;
    // }

    // // Print the names of the cross-sections we are going to extract
    // std::cout << "The following cross-sections are enabled:" << std::endl;
    // for (const auto &[selectionName, xsecs] : xsecMap)
    // {
    //     std::cout << "  - Selection: " << selectionName << std::endl;
    //     for (const auto &entry : xsecs)
    //     {
    //         const auto &name = entry.first;
    //         std::cout << "    - " << name << std::endl;
    //     }
    // }

    // -------------------------------------------------------------------------------------------------------------------------------------
    // Setup the relevent "getters" for each the sideband
    // -------------------------------------------------------------------------------------------------------------------------------------
    // Here we define a map from the name of each cross-section to a function which pulls out the relevant kinematic quanitity from an input
    // analysis data object. Again this is done up-front to reduce code-bloat below.
    std::unordered_map< std::string, std::function<float(const AnalysisHelper::AnalysisData &)> > getSidebandValue;

    // Differential cross-section kinematic parameters
    getSidebandValue.emplace("muonCosTheta",  [](const auto &data) { return data.muonCosTheta;  });
    getSidebandValue.emplace("muonPhi",       [](const auto &data) { return data.muonPhi;       });
    getSidebandValue.emplace("muonMomentum",  [](const auto &data) { return data.muonMomentum;  });

    getSidebandValue.emplace("pionCosTheta",  [](const auto &data) { return data.protonCosTheta;  }); // Getting proton instead of pion values
    getSidebandValue.emplace("pionPhi",       [](const auto &data) { return data.protonPhi;       }); // Getting proton instead of pion values
    getSidebandValue.emplace("pionMomentum",  [](const auto &data) { return data.protonMomentum;  }); // Getting proton instead of pion values

    getSidebandValue.emplace("muonPionAngle", [](const auto &data) { return data.muonProtonAngle; }); // Getting proton instead of pion values
    getSidebandValue.emplace("nProtons",      [](const auto &data) { return data.nProtons-1;      }); //Leading proton treated as pion in CC0pi analysis

    // ATTN as described above, for the total cross-section we don't have an associated kinematic quantity so we just return a dummy value
    getSidebandValue.emplace("total", [=](const auto &) { return dummyValue; });



    // // -------------------------------------------------------------------------------------------------------------------------------------
    // // Generate CC0pi bin weights
    // // -------------------------------------------------------------------------------------------------------------------------------------
    // // if(config.global.useCC0piConstraint)
    // // {
    // std::map<std::string,std::vector<float>> cc0piConstraintMap;
    // for (auto &[selectionName, xsecs] : xsecMap)
    // {
    //     for (auto &[name, xsec] : xsecs)
    //     {
    //         std::vector<float> cc0piWeights;
    //         std::ifstream signalSelected(("CC0pi_" + name + "_signal_selected_eventRate.txt").c_str());
    //         std::ifstream dataSelected(("CC0pi_" + name + "_data_selected_eventRate.txt").c_str());
    //         std::ifstream backgroundSelected(("CC0pi_" + name + "_background_selected_eventRate.txt").c_str());
    //         std::string signalBinValue;
    //         std::string dataBinValue;
    //         std::string backgroundBinValue;
    //         // Read the next line from File untill it reaches the end.
    //         while (std::getline(signalSelected, signalBinValue) && std::getline(dataSelected, dataBinValue) && std::getline(backgroundSelected, backgroundBinValue))
    //         {
    //             cc0piWeights.push_back(std::stof(dataBinValue)/(std::stof(signalBinValue)+std::stof(backgroundBinValue)));
    //         }
    //         if(std::getline(signalSelected, signalBinValue) || std::getline(dataSelected, dataBinValue) || std::getline(backgroundSelected, backgroundBinValue))
    //             throw std::logic_error("ExtractSidebandFit - CC0pi weight files have different numbers of entries.");
    //         cc0piConstraintMap.emplace(name, cc0piWeights);
    //     }
    // }
    // // }
    // Loop over the files

    // for (const auto &[sampleType, sampleName, fileName, normalisation] : inputData)
    // {

    //     const auto isOverlay = (sampleType == AnalysisHelper::Overlay);
    //     const auto isDirt    = (sampleType == AnalysisHelper::Dirt);

    //     // Open the input file for reading and enable the branches with systematic event weights (if required)
    //     FileReader reader(fileName);

    //     auto pEvent = reader.GetBoundEventAddress();

    //     // Loop over the events in the file
    //     const auto nEvents = reader.GetNumberOfEvents();
    //     for (unsigned int i = 0; i < nEvents; ++i)
    //     {
    //                     if (isDataBNB)
    //         {
    //         for (auto &[selectionName, xsecs] : xsecMapSideband)
    //         {
    //             // Determine if we passed the relevant selection
    //             const auto isSelected = isSelectedMap.at(selectionName);
    //             // Only count events passing the selection
    //             if (!isSelected)
    //                 continue;

    //             for (auto &[name, xsec] : xsecs)
    //             {


    // -------------------------------------------------------------------------------------------------------------------------------------
    // -------------------------------------------------------------------------------------------------------------------------------------
    // -------------------------------------------------------------------------------------------------------------------------------------
    // -------------------------------------------------------------------------------------------------------------------------------------
    // -------------------------------------------------------------------------------------------------------------------------------------
    // -------------------------------------------------------------------------------------------------------------------------------------
    // Count the events for CC0pi
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

            // -----------------------------------------------------------------------------------------------------------------------------
            // Work out if this event passed the selection and apply any additional phase-space cuts based on the input binning
            // -----------------------------------------------------------------------------------------------------------------------------

            // Run the selection
            const auto &[passedSidebandSelection, sidebandCutsPassed, sidebandAssignedPdgCodes] = sidebandSelection.Execute(pEvent);

            // Get the reco analysis data (if available, otherwise set to dummy values)
            const auto recoData = (
                passedSidebandSelection
                    ? AnalysisHelper::GetRecoAnalysisDataCC0Pi(pEvent->reco, sidebandAssignedPdgCodes)
                    : AnalysisHelper::GetDummyAnalysisData()
            );

            // Here we apply reco-level phase-space restrictions
            // For any event that passes the generic selection, get the value of the kinematic quantity and check if it is outside of the
            // min/max values supplied in the binning. If so, then reject the event.
            bool passesPhaseSpaceReco = false;
            if (passedSidebandSelection)
            {
                // Start by assuming the event passes the phase-space cuts
                passesPhaseSpaceReco = true;

                // Check the value of the kinematic quantities are within the phase-space limits
                for (const auto &[name, minMax] : phaseSpaceMap)
                {
                    const auto &[min, max] = minMax;
                    const auto value = getSidebandValue.at(name)(recoData);

                    if (value < min || value > max)
                    {
                        passesPhaseSpaceReco = false;
                        std::cout << "Event failed reco phase-space cuts" << std::endl;
                        break;
                    }
                }
            }

            const auto isSelectedSideband = passedSidebandSelection  && passesPhaseSpaceReco;
            // std::map<std::string, bool> isSelectedMap = {{"sideband",isSelectedSideband}};
            std::map<std::string, bool> isSelectedMap = {{"generic",isSelectedSideband},{"golden",isSelectedSideband}};

            // -----------------------------------------------------------------------------------------------------------------------------
            // Handle BNB data
            // -----------------------------------------------------------------------------------------------------------------------------
            if (isDataBNB)
            {
                for (auto &[selectionName, xsecs] : xsecMapSideband)
                {
                    // Determine if we passed the relevant selection
                    const auto isSelected = isSelectedMap.at(selectionName);
                    // Only count events passing the selection
                    if (!isSelected)
                        continue;

                    for (auto &[name, xsec] : xsecs)
                    {
                        xsec.AddSelectedBNBDataEvent(getSidebandValue.at(name)(recoData));
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

            // Determine if this is truly a CC0Pi event
            const auto isTrueCC0Pi = (isOverlay || isDetVar) && AnalysisHelper::IsTrueCC0Pi(pEvent, config.global.useAbsPdg, config.global.protonMomentumThreshold);            

            // Get the truth analysis data (if available, otherwise set to dummy values)
            const auto truthData = (
                (isTrueCC0Pi)
                    ? AnalysisHelper::GetTruthAnalysisDataCC0Pi(pEvent->truth, config.global.useAbsPdg, config.global.protonMomentumThreshold)
                    : AnalysisHelper::GetDummyAnalysisData()
            );

            // Here we apply truth-level phase-space restrictions
            // For all true CC1Pi events, we check if the values of each kinematic variable are within the supplied limits. If not then the
            // event is not classed as "signal"
            bool passesPhaseSpaceTruth = false;
            if (isTrueCC0Pi)
            {
                // Start by assuming the event passes the phase-space cuts
                passesPhaseSpaceTruth = true;

                // Check the value of the kinematic quantities are within the phase-space limits
                for (const auto &[name, minMax] : phaseSpaceMap)
                {
                    const auto &[min, max] = minMax;
                    const auto value = getSidebandValue.at(name)(truthData);

                    if (value < min || value > max)
                    {
                        std::cout << "Event failed truth phase-space cuts - name: "<<name<<" - min: "<<min<<" - value: "<<value<<" - max: "<<max << std::endl;                        
                        passesPhaseSpaceTruth = false;
                        break;
                    }
                }
            }

            const auto isCC0PiSignal = isTrueCC0Pi && passesPhaseSpaceTruth;
            // std::map<std::string, bool> isSignalMap = {{"sideband",isCC0PiSignal}};
            std::map<std::string, bool> isSignalMap = {{"generic",isCC0PiSignal},{"golden",isCC0PiSignal}};

            // -----------------------------------------------------------------------------------------------------------------------------
            // Handle the detector variation samples (unisims)
            // -----------------------------------------------------------------------------------------------------------------------------
            if (isDetVar)
            {
                for (auto &[selectionName, xsecs] : xsecMapSideband)
                {
                    const auto isSignal = isSignalMap.at(selectionName);
                    // Handle signal events
                    if (isSignal)
                    {
                        // Determine if we passed the relevant selection
                        const auto isSelected = isSelectedMap.at(selectionName);

                        for (auto &[name, xsec] : xsecs)
                        {
                            const auto recoValue = getSidebandValue.at(name)(recoData);
                            const auto trueValue = getSidebandValue.at(name)(truthData);
                            xsec.AddSignalEventDetVar(recoValue, trueValue, isSelected, weight, sampleName);
                        }
                    }
                    // Handle selected background events
                    else
                    {
                        for (auto &[selectionName, xsecs] : xsecMapSideband)
                        {
                            // Only use selected background events
                            const auto isSelected = isSelectedMap.at(selectionName);
                            if (!isSelected)
                                continue;

                            for (auto &[name, xsec] : xsecs)
                            {
                                const auto recoValue = getSidebandValue.at(name)(recoData);
                                xsec.AddSelectedBackgroundEventDetVar(recoValue, weight, sampleName);
                            }
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

            // Fill the flux-reweightor with all overlay events from neutrinos of the desired flavour in the active volume
            if (isOverlay && AnalysisHelper::IsInActiveVolume(pEvent->truth.nuVertex()) &&
                std::find(config.flux.nuPdgsSignal.begin(), config.flux.nuPdgsSignal.end(), pEvent->truth.nuPdgCode()) != config.flux.nuPdgsSignal.end())
            {
                scalingData.pFluxReweightor->AddEvent(pEvent->truth.nuEnergy(), weight, fluxWeights);
            }
            // Get the cross-section weights
            // ATTN here we optionally scale the cross-section weights down by the genieTuneEventWeight - this is done so we don't
            // double count this weight (once in the nominal event weight, and once in the xsec systematic event weights)
            // const auto xsecWeightsScaleFactor = 1.f;
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
            for (auto &[selectionName, xsecs] : xsecMapSideband)
            {
                const auto isSignal = isSignalMap.at(selectionName);
                // Handle signal events
                if (isSignal)
                {
                    for (auto &[selectionName, xsecs] : xsecMapSideband)
                    {
                        // Determine if we passed the relevant selection
                        const auto isSelected = isSelectedMap.at(selectionName);

                        for (auto &[name, xsec] : xsecs)
                        {
                            const auto recoValue = getSidebandValue.at(name)(recoData);
                            const auto trueValue = getSidebandValue.at(name)(truthData);
                            const auto seedString =  selectionName + name + std::to_string(i);
                            xsec.AddSignalEvent(recoValue, trueValue, isSelected, weight, fluxWeights, xsecWeights, reintWeights, seedString);
                        }
                    }
                }
                // Handle selected background events
                else
                {
                    for (auto &[selectionName, xsecs] : xsecMapSideband)
                    {
                        // Only use selected background events
                        const auto isSelected = isSelectedMap.at(selectionName);
                        if (!isSelected)
                            continue;

                        for (auto &[name, xsec] : xsecs)
                        {
                            const auto recoValue = getSidebandValue.at(name)(recoData);
                            const auto seedString =  selectionName + name + std::to_string(i);
                            std::vector<float> bootstrapWeights; // Parameter only needed for CC1pi
                            xsec.AddSelectedBackgroundEvent(recoValue, isDirt, weight, fluxWeights, xsecWeights, reintWeights, bootstrapWeights, seedString);
                        }
                    }
                }
            }
        }
    }

    // -------------------------------------------------------------------------------------------------------------------------------------
    // Calculate the sideband weights
    // -------------------------------------------------------------------------------------------------------------------------------------
    // Loop over all cross-section objects
    typedef std::pair<std::vector<Double_t>,std::vector<Double_t>> paramAndErrorPair; // Todo: Improve code!
    //Parameters: selectionName name
    std::map<std::string,std::map<std::string, std::vector<float>>> cc0piCovarianceMap; 
    std::map<std::string, std::map<std::string, paramAndErrorPair>> cc0piNominalConstraintMap;
    //Parameters: selectionName name paramName (i.e. golden muonMomentum hadronProduction)
    std::map<std::string, std::map<std::string, std::map<std::string, std::vector<paramAndErrorPair>>>> cc0piUniverseConstraintMap;
    
    try
    {
        for (const auto &[selectionName, xsecs] : xsecMapSideband)
        {
            for (const auto &[name, xsec] : xsecs)
            {
                // -------------------------------------------------------------------------------------------------------------------------------------
                // Fit nominal
                // -------------------------------------------------------------------------------------------------------------------------------------

                // std::cout<<"_______________________Fitting: "<<selectionName<<" - "<<name<<"_______________________"<<std::endl;
                // const auto sidebandCVFit = xsec.GetPredictedCrossSection(scalingData);

                // std::cout<<"_______________________Fitting Point 0.1"<<std::endl;
                const auto selectedEventsData = xsec.GetSelectedBNBDataEvents();
                // std::cout<<"_______________________Fitting Point 0.2"<<std::endl;
                // Get the smearing matrix of selected events
                const auto smearingMatrixAllSelected = xsec.GetSmearingMatrixAllSelected();
                // std::cout<<"_______________________Fitting Point 0.3"<<std::endl;
                const auto selectedEventsBackgroundReco = xsec.GetSelectedBackgroundEvents();
                // std::cout<<"_______________________Fitting Point 0.4"<<std::endl;
                const auto selectedEventsSignalTruth = xsec.GetSelectedSignalEvents();
                // std::cout<<"_______________________Fitting Point 0.5"<<std::endl;
                auto signalData = selectedEventsData - selectedEventsBackgroundReco;

                FormattingHelper::SaveMatrix(selectedEventsSignalTruth, "SidebandFit_" + selectionName + "_" + name + "_selectedEventsSignalTruth.txt");
                FormattingHelper::SaveMatrix(selectedEventsData, "SidebandFit_" + selectionName + "_" + name + "_selectedEventsData.txt");
                FormattingHelper::SaveMatrix(selectedEventsBackgroundReco, "SidebandFit_" + selectionName + "_" + name + "_selectedEventsBackgroundReco.txt");
                FormattingHelper::SaveMatrix(smearingMatrixAllSelected, "SidebandFit_" + selectionName + "_" + name + "_smearingMatrixAllSelected.txt");


                for(unsigned int r = 0; r<signalData.GetRows(); r++) // DEBUG - TODO: REMOVE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                {
                    for(unsigned int c = 0; c<signalData.GetColumns(); c++) // DEBUG - TODO: REMOVE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    {
                        if(signalData.At(r, c) < 0)
                        {
                            std::cout<<"signalData.At("<<r<<","<<c<<") = "<<signalData.At(r, c)<<std::endl;
                            signalData.SetElement(r, c, std::max(signalData.At(r, c), 0.f));
                        }
                    }
                }

                // std::cout<<"_______________________Fitting Point 1"<<std::endl;
                // const auto &[pPredictionStatBias, pPredictionStatCovariance] = xsec.GetPredictedCrossSectionStatUncertainty(scalingData);
                // const auto predictionErrorMatrix = CrossSectionHelper::GetErrorMatrix(*pPredictionBiasVector, *pPredictionCovarianceMatrix);\

                std::vector<float> elements;
                const auto nBins = signalData.GetRows();
                // std::cout<<"_______________________Fitting Point 1.1"<<std::endl;
                for (unsigned int iBin = 0; iBin < nBins; ++iBin)
                {
                    // std::cout<<"_______________________Fitting Point 1.2"<<std::endl;
                    const auto value = signalData.At(iBin, 0);
                    if (value<0)
                    {
                        std::cout<<"ERROR: ExtractXSec - Background-removed signal data is negative."<<std::endl;
                        // throw std::logic_error("ERROR: ExtractXSec - Background-removed signal data is negative."); // TODO: Uncomment
                    }
                    elements.push_back(AnalysisHelper::GetCountUncertainty(std::max(value,0.f)));
                    // elements.push_back(AnalysisHelper::GetCountUncertainty(value));
                    // std::cout<<"_______________________Fitting Point 1.3"<<std::endl;
                }
                // std::cout<<"_______________________Fitting Point 2"<<std::endl;
                // const ubsmear::UBMatrix signalDataUncertainty(elements, nBins, 1);

                // std::cout<<"_______________________Fitting Point 3"<<std::endl;
                x = selectedEventsSignalTruth.GetValues();
                y = signalData.GetValues();
                errorY = elements;//signalDataUncertainty.GetValues();
                S = smearingMatrixAllSelected.GetValues();

                // std::cout<<"\nx (selectedEventsSignalTruth): \n";
                for (const auto &xValue : x)
                    std::cout<<xValue<<" ";

                // std::cout<<"\ny (signalData): \n";
                for (const auto &yValue : y)
                    std::cout<<yValue<<" ";

                // std::cout<<"\nerrorY (signalDataUncertainty): \n";
                for (const auto &errorYValue : errorY)
                    std::cout<<errorYValue<<" ";

                // std::cout<<"\nSmearing matrix: \n";
                for(unsigned int i = 0; i<S.size(); i++)
                {
                    if(i%nBins==0)
                        std::cout<<"\n";
                    std::cout<<S[i]<<" ";
                }
                

                // std::cout<<"_______________________Fitting Point 4"<<std::endl;
                // auto minimizer = FittingHelper(selectedEventsSignal, signalData, signalDataUncertainty, smearingMatrix);
                auto minimizer = FittingHelper(nBins);
                std::pair<std::vector<Double_t>, std::vector<Double_t>> result;
                
                std::vector<float> fitCovMatrixVector;
                minimizer.Fit(fcn, result, fitCovMatrixVector, 0);
                cc0piCovarianceMap[selectionName].emplace(name, fitCovMatrixVector);

                // std::cout<<"\nFitting covariance matrix: \n";
                // for(unsigned int i = 0; i<fitCovMatrixVector.size(); i++)
                // {
                //     if(i%nBins==0)
                //         std::cout<<"\n";
                //     std::cout<<fitCovMatrixVector[i]<<" ";
                // }
                // std::cout<<"_______________________Fitting Point 4.1"<<std::endl;
                // vector<float> fitCovMatrixFloat(fitCovMatrix.begin(), fitCovMatrix.end()); //Todo avoid this
                const ubsmear::UBMatrix sidebandCovMatrix(fitCovMatrixVector, nBins, nBins);
                // std::cout<<"_______________________Fitting Point 4.2"<<std::endl;
                FormattingHelper::SaveMatrix(sidebandCovMatrix, "SidebandFit_" + selectionName + "_" + name + "_sideband_stat_covariance.txt");
                // std::cout<<"_______________________Fitting Point 4.3"<<std::endl;
                // std::cout<<"_______________________Fitting Point 4.4"<<std::endl;
                vector<float> paramVector(result.first.begin(), result.first.end()); //Todo avoid this
                vector<float> paramErrorVector(result.second.begin(), result.second.end()); //Todo avoid this
                // std::cout<<"_______________________Fitting Point 4.5"<<std::endl;
                // std::cout<<"nBins :"<<nBins<<std::endl;
                // std::cout<<"\nerrorY (paramErrorVector Double_t): \n";
                // for (const auto &e : result.second)
                //     std::cout<<e<<" ";
                // std::cout<<"\nerrorY (paramErrorVector float): \n";
                // for (const auto &e : paramErrorVector)
                //     std::cout<<e<<" ";
                // std::cout<<std::endl;
                const ubsmear::UBMatrix sidebandParamVectorTruth(paramVector, nBins, 1);
                const ubsmear::UBMatrix sidebandErrorVectorTruth(paramErrorVector, nBins, 1);
                // std::cout<<"_______________________Fitting Point 4.6"<<std::endl;
                // const auto sidebandErrorVectorReco = smearingMatrixAllSelected*sidebandErrorVectorTruth; // Todo: check this multiplication is correctly computed 
                FormattingHelper::SaveMatrix(sidebandParamVectorTruth, "SidebandFit_" + selectionName + "_" + name + "_sideband_parameterVector.txt");
                FormattingHelper::SaveMatrix(sidebandErrorVectorTruth, "SidebandFit_" + selectionName + "_" + name + "_sideband_parameterErrorVector.txt");
                // const ubsmear::UBMatrix sidebandCovMatrix(fitCovMatrixVector, nBins, nBins);
                // FormattingHelper::SaveMatrix(sidebandCovMatrix, "SidebandFit_" + selectionName + "_" + name + "_sideband_covariance_matrix.txt");

                // std::cout<<"_______________________Fitting Point 5"<<std::endl;

                // std::cout<<"param: \n";
                // for (const auto &p : result.first)
                //     std::cout<<p<<" ";

                // std::cout<<"paramError: \n";
                // for (const auto &p : result.second)
                //     std::cout<<p<<" ";


                cc0piNominalConstraintMap[selectionName].emplace(name, result);
                // FittingHelper::Fit(selectedEventsSignal, signalData, signalDataUncertainty, smearingMatrix);

                // const auto sidebandWeights = xsec.GetSidebandWeights(scalingData);
                
                for (const auto &[paramName, nUniverses] : systParams.xsecDimensions)
                {
                    // -------------------------------------------------------------------------------------------------------------------------------------
                    // Fit each universe
                    // -------------------------------------------------------------------------------------------------------------------------------------

                    // const auto nUniverses = config.extractXSecs.nBootstrapUniverses;

                    const auto selectedSignalTruthUniverses = xsec.GetSelectedSignalRecoTruthMap().at("xsec").at(paramName);
                    const auto selectedBackgroundRecoUniverses = xsec.GetSelectedBackgroundRecoMap().at("xsec").at(paramName);

                    std::vector<paramAndErrorPair> resultVector;
                    for (unsigned int iUni = 0; iUni < nUniverses; ++iUni)
                    {
                        // AnalysisHelper::PrintLoadingBar(iUni, nUniverses);
                        const auto selectedSignalTruth = xsec.GetSignalSelectedTrue(selectedSignalTruthUniverses.at(iUni));
                        const auto selectedBackgoundReco = CrossSectionHelper::GetMatrixFromHist(selectedBackgroundRecoUniverses.at(iUni));
                        auto signalData = selectedEventsData - selectedBackgoundReco;

                        for(unsigned int r = 0; r<signalData.GetRows(); r++) // DEBUG - TODO: REMOVE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        {
                            for(unsigned int c = 0; c<signalData.GetColumns(); c++) // DEBUG - TODO: REMOVE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                            {
                                if(signalData.At(r, c) < 0)
                                {
                                    std::cout<<"universe signalData.At("<<r<<","<<c<<") = "<<signalData.At(r, c)<<std::endl;
                                    signalData.SetElement(r, c, std::max(signalData.At(r, c), 0.f));
                                }
                            }
                        }

                        
                        std::vector<float> elements;
                        const auto nBins = signalData.GetRows();
                        for (unsigned int iBin = 0; iBin < nBins; ++iBin)
                        {
                            const auto value = signalData.At(iBin, 0);
                            if (value<0)
                            {
                                std::cout<<"ERROR: ExtractXSec - Background-removed universe (!!!) signal data is negative."<<std::endl;
                                // throw std::logic_error("ERROR: ExtractXSec - Background-removed signal data is negative."); // TODO: Uncomment
                            }
                            elements.push_back(AnalysisHelper::GetCountUncertainty(std::max(value,0.f)));
                            // elements.push_back(AnalysisHelper::GetCountUncertainty(value));
                        }

                        x = selectedSignalTruth.GetValues();
                        y = signalData.GetValues();
                        errorY = elements;
                        std::pair<std::vector<Double_t>, std::vector<Double_t>> result;
                        std::vector<float> covMatrixInUniverse;

                        // std::cout<<"\nx (universe selectedEventsSignalTruth): \n";
                        // for (const auto &xValue : x)
                        //     std::cout<<xValue<<" ";

                        // std::cout<<"\ny (universe signalData): \n";
                        // for (const auto &yValue : y)
                        //     std::cout<<yValue<<" ";

                        // std::cout<<"\nerrorY (universe signalDataUncertainty): \n";
                        // for (const auto &errorYValue : errorY)
                        //     std::cout<<errorYValue<<" ";

                        minimizer.Fit(fcn, result, covMatrixInUniverse, 0);

                        // std::cout<<"\nParameter("<<iUni<<"):";
                        // for (const auto &r : result.first)
                        //     std::cout<<" "<<r;
                            
                        // std::cout<<"\nUncertainty:";
                        // for (const auto &r : result.second)
                        //     std::cout<<" "<<r;
                        
                        resultVector.push_back(result);
                    }
                    cc0piUniverseConstraintMap[selectionName][name].emplace(paramName, resultVector);
                }
            }
        }
    }
    catch(exception &e)
    {
        cout << "CC0pi - Caught exception Point 1: "<<e.what();
        return;
    }

	std::ofstream ofs1("cc0piCovarianceMap.bin", std::ios::binary);
    std::ofstream ofs2("cc0piNominalConstraintMap.bin", std::ios::binary);
    std::ofstream ofs3("cc0piUniverseConstraintMap.bin", std::ios::binary);
	

    boost::archive::binary_oarchive oarch1(ofs1);
    boost::archive::binary_oarchive oarch2(ofs2);
    boost::archive::binary_oarchive oarch3(ofs3);

    oarch1 << cc0piCovarianceMap;
    oarch2 << cc0piNominalConstraintMap;
    oarch3 << cc0piUniverseConstraintMap;
	
	ofs1.close();
    ofs2.close();
    ofs3.close();

    std::cout<<"------------- All Done -------------"<<std::endl;

    return;
}

} // namespace ubcc1pi_macros
