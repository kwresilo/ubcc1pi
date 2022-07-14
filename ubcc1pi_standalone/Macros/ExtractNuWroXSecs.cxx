/**
 *  @file  ubcc1pi_standalone/Macros/ExtractNuWroXSecs.cxx
 *
 *  @brief The implementation file of the ExtractNuWroXSecs macro
 */

#include "ubcc1pi_standalone/Macros/Macros.h"
#include "ubcc1pi_standalone/Objects/FileReader.h"

#include "ubcc1pi_standalone/Helpers/AnalysisHelper.h"
#include "ubcc1pi_standalone/Helpers/NormalisationHelper.h"
#include "ubcc1pi_standalone/Helpers/SelectionHelper.h"
#include "ubcc1pi_standalone/Helpers/CrossSectionHelper.h"
#include "ubcc1pi_standalone/Helpers/FormattingHelper.h"
// #include "ubcc1pi_standalone/Helpers/FittingHelper.h"
#include "ubsmear.h"

#include <fstream> // Todo: not use txt files
// Boost libraries
#include "binary_iarchive.hpp"
// #include "binary_oarchive.hpp"
#include "binary_object.hpp"
#include "map.hpp"
#include "vector.hpp"

using namespace ubcc1pi;

namespace ubcc1pi_macros
{

void ExtractNuWroXSecs(const Config &config)
{

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
    CrossSectionHelper::CrossSection::ScalingData scalingDataTruth;

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
    // auto sidebandSelection = SelectionHelper::GetSelection("CC0pi");

    // -------------------------------------------------------------------------------------------------------------------------------------
    // Setup the cross-section objects
    // -------------------------------------------------------------------------------------------------------------------------------------
    // Here we make a map from a name of the cross-section to the cross-section object itself. In this way, we can iterate through the
    // cross-section objects and reduce code-bloat. The first index is an identifier for the selection that's applied (generic or goldlen),
    // the second index is an identifier for the kinematic quantity that's relevant for the cross-section (e.g. muonMomentum), and the
    // mapped type is the cross-section object.
    std::map<std::string, std::map<std::string, CrossSectionHelper::CrossSection> > xsecMap;
    std::map<std::string, std::map<std::string, CrossSectionHelper::CrossSection> > xsecMapUnscaled;
    std::map<std::string, std::map<std::string, CrossSectionHelper::CrossSection> > xsecMapTrue; // True NuWro cross-sections
    std::map<std::string, std::map<std::string, CrossSectionHelper::CrossSection> > xsecMapOnlyCC0Pi; // Map is used to save both nuwro as data and as preiction - not usable to calculate a cross-section 
    std::map<std::string, std::map<std::string, CrossSectionHelper::CrossSection> > xsecMapOnlyCC0PiUnscaled; // Map is used to save both nuwro as data and as prediction - not usable to calculate a cross-section

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
        // Add the cross-section object to the map using the binning from the input configuration
        const auto &[extendedBinEdges, hasUnderflow, hasOverflow] = CrossSectionHelper::GetExtendedBinEdges(binning.min, binning.max, binning.binEdges);
        for (const auto &selectionName : {"generic", "golden"})
        {
            // Don't setup a cross-section object if it's been disabled in the configuration
            if (config.extractXSecs.crossSectionIsEnabled.at(selectionName).at(name))
            {
                xsecMap[selectionName].emplace(name, CrossSectionHelper::CrossSection(systParams, extendedBinEdges, hasUnderflow, hasOverflow, scaleByBinWidth));
                xsecMapUnscaled[selectionName].emplace(name, CrossSectionHelper::CrossSection(systParams, extendedBinEdges, hasUnderflow, hasOverflow, scaleByBinWidth));
                xsecMapTrue[selectionName].emplace(name, CrossSectionHelper::CrossSection(systParams, extendedBinEdges, hasUnderflow, hasOverflow, scaleByBinWidth));
                xsecMapOnlyCC0Pi[selectionName].emplace(name, CrossSectionHelper::CrossSection(systParams, extendedBinEdges, hasUnderflow, hasOverflow, scaleByBinWidth));
                xsecMapOnlyCC0PiUnscaled[selectionName].emplace(name, CrossSectionHelper::CrossSection(systParams, extendedBinEdges, hasUnderflow, hasOverflow, scaleByBinWidth));
            }
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
        if (config.extractXSecs.crossSectionIsEnabled.at(selectionName).at("total"))
        {
            xsecMap[selectionName].emplace("total", CrossSectionHelper::CrossSection(systParams, {-1.f, 1.f}, false, false, false));
            xsecMapUnscaled[selectionName].emplace("total", CrossSectionHelper::CrossSection(systParams, {-1.f, 1.f}, false, false, false));
            xsecMapTrue[selectionName].emplace("total", CrossSectionHelper::CrossSection(systParams, {-1.f, 1.f}, false, false, false));
            xsecMapOnlyCC0Pi[selectionName].emplace("total", CrossSectionHelper::CrossSection(systParams, {-1.f, 1.f}, false, false, false));
            xsecMapOnlyCC0PiUnscaled[selectionName].emplace("total", CrossSectionHelper::CrossSection(systParams, {-1.f, 1.f}, false, false, false));
        }
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

    // // -------------------------------------------------------------------------------------------------------------------------------------
    // // Load the sideband weights from file
    // // -------------------------------------------------------------------------------------------------------------------------------------
    // // Loop over all cross-section objects
    typedef std::pair<std::vector<Double_t>,std::vector<Double_t>> paramAndErrorPair; // Todo: Improve code!
    //Parameters: selectionName, name
    std::map<std::string,std::map<std::string, std::vector<float>>> cc0piCovarianceMap; 
    std::map<std::string, std::map<std::string, paramAndErrorPair>> cc0piNominalConstraintMap;
    //Parameters: selectionName, name, paramName (i.e. golden, muonMomentum, hadronProduction)
    std::map<std::string, std::map<std::string, std::map<std::string, std::vector<paramAndErrorPair>>>> cc0piUniverseConstraintMap;

    std::ifstream ifs1("cc0piCovarianceMapNuWro.bin", std::ios::binary);
    std::ifstream ifs2("cc0piNominalConstraintMapNuWro.bin", std::ios::binary);
    std::ifstream ifs3("cc0piUniverseConstraintMapNuWro.bin", std::ios::binary);
    
    boost::archive::binary_iarchive iarch1(ifs1);
    boost::archive::binary_iarchive iarch2(ifs2);
    boost::archive::binary_iarchive iarch3(ifs3);

    iarch1 >> cc0piCovarianceMap;
    iarch2 >> cc0piNominalConstraintMap;
    iarch3 >> cc0piUniverseConstraintMap;

    ifs1.close();
    ifs2.close();
    ifs3.close();

    // -------------------------------------------------------------------------------------------------------------------------------------
    // Setup the relevent "getters" for each cross-section and for the sideband
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

    std::unordered_map< std::string, std::function<float(const AnalysisHelper::AnalysisData &)> > getSidebandValue;
    // Differential cross-section kinematic parameters
    getSidebandValue.emplace("muonCosTheta",  [](const auto &data) { return data.muonCosTheta;  });
    getSidebandValue.emplace("muonPhi",       [](const auto &data) { return data.muonPhi;       });
    getSidebandValue.emplace("muonMomentum",  [](const auto &data) { return data.muonMomentum;  });
    getSidebandValue.emplace("pionCosTheta",  [](const auto &data) { return data.protonCosTheta;  }); // Getting proton instead of pion values
    getSidebandValue.emplace("pionPhi",       [](const auto &data) { return data.protonPhi;       }); // Getting proton instead of pion values
    getSidebandValue.emplace("pionMomentum",  [](const auto &data) { return data.protonMomentum;  }); // Getting proton instead of pion values
    getSidebandValue.emplace("muonPionAngle", [](const auto &data) { return data.muonProtonAngle; }); // Getting proton instead of pion values
    getSidebandValue.emplace("nProtons",      [](const auto &data) { return data.nProtons-1;      }); // Leading proton treated as pion in CC0pi analysis
    // ATTN as described above, for the total cross-section we don't have an associated kinematic quantity so we just return a dummy value
    getSidebandValue.emplace("total", [=](const auto &) { return dummyValue; });

    // -------------------------------------------------------------------------------------------------------------------------------------
    // Setup the input files
    // -------------------------------------------------------------------------------------------------------------------------------------
    // Here we define a vector of tuples with 4 entries
    //   - First, the sample type (e.g. overlay)
    //   - Second, a string which is used to identify a given detector variation sample (for other sample type, this is unused)
    //   - Third, the path to the input file
    //   - Fourth, the normalisation factor to apply to all events in that file
    const auto fluxHistNames = CrossSectionHelper::GetNominalFluxHistNames(config.flux.nuPdgsSignal, config.flux.nuPdgToHistName, config.flux.nomHistPattern);
    const auto &[fluxBinEdges, fluxValues] = CrossSectionHelper::ReadNominalFlux(config.flux.fileName, fluxHistNames, config.flux.pot);
    scalingData.pFluxReweightor = std::make_shared<CrossSectionHelper::FluxReweightor>(fluxBinEdges, fluxValues, systParams.fluxDimensions);
    scalingDataTruth.pFluxReweightor = std::make_shared<CrossSectionHelper::FluxReweightor>(fluxBinEdges, fluxValues, systParams.fluxDimensions);
    std::vector< std::tuple<AnalysisHelper::SampleType, std::string, std::string, float> > inputData;
    auto totalExposurePOT = 0.f;

    if(std::binary_search(config.global.runs.begin(), config.global.runs.end(), 1))
    {
        inputData.emplace_back(AnalysisHelper::Overlay, "", config.filesRun1.overlaysFileName, NormalisationHelper::GetOverlaysNormalisation(config, 1));
        // inputData.emplace_back(AnalysisHelper::Dirt,    "", config.filesRun1.dirtFileName, NormalisationHelper::GetDirtNormalisationToNuWro(config, 1));
        // inputData.emplace_back(AnalysisHelper::DataEXT, "", config.filesRun1.dataEXTFileName, NormalisationHelper::GetDataEXTNormalisationToNuWro(config, 1));
        inputData.emplace_back(AnalysisHelper::DataBNB, "", config.filesRun1.nuWroFileName, NormalisationHelper::GetNuWroNormalisation(config, 1));

        // Add the detector variation files
        for (const auto &[name, fileName] : config.filesRun1.detVarFiles)
        {
            inputData.emplace_back(AnalysisHelper::DetectorVariation, name, fileName, NormalisationHelper::GetDetectorVariationNormalisation(config, name, 1));
        }

        totalExposurePOT += config.normsRun1.nuWroPOT / (1e20);
    }
    if(std::binary_search(config.global.runs.begin(), config.global.runs.end(), 2))
    {
        inputData.emplace_back(AnalysisHelper::Overlay, "", config.filesRun2.overlaysFileName, NormalisationHelper::GetOverlaysNormalisation(config, 2));
        // inputData.emplace_back(AnalysisHelper::Dirt,    "", config.filesRun2.dirtFileName, NormalisationHelper::GetDirtNormalisationToNuWro(config, 2));
        // inputData.emplace_back(AnalysisHelper::DataEXT, "", config.filesRun2.dataEXTFileName, NormalisationHelper::GetDataEXTNormalisationToNuWro(config, 2));
        inputData.emplace_back(AnalysisHelper::DataBNB, "", config.filesRun2.nuWroFileName, NormalisationHelper::GetNuWroNormalisation(config, 2));

        // Add the detector variation files
        for (const auto &[name, fileName] : config.filesRun2.detVarFiles)
            inputData.emplace_back(AnalysisHelper::DetectorVariation, name, fileName, NormalisationHelper::GetDetectorVariationNormalisation(config, name, 2));

        totalExposurePOT += config.normsRun2.nuWroPOT / (1e20);
    }
    if(std::binary_search(config.global.runs.begin(), config.global.runs.end(), 3))
    {
        inputData.emplace_back(AnalysisHelper::Overlay, "", config.filesRun3.overlaysFileName, NormalisationHelper::GetOverlaysNormalisation(config, 3));
        // inputData.emplace_back(AnalysisHelper::Dirt,    "", config.filesRun3.dirtFileName, NormalisationHelper::GetDirtNormalisationToNuWro(config, 3));
        // inputData.emplace_back(AnalysisHelper::DataEXT, "", config.filesRun3.dataEXTFileName, NormalisationHelper::GetDataEXTNormalisationToNuWro(config, 3));
        inputData.emplace_back(AnalysisHelper::DataBNB, "", config.filesRun3.nuWroFileName, NormalisationHelper::GetNuWroNormalisation(config, 3));

        // Add the detector variation files
        for (const auto &[name, fileName] : config.filesRun3.detVarFiles)
            inputData.emplace_back(AnalysisHelper::DetectorVariation, name, fileName, NormalisationHelper::GetDetectorVariationNormalisation(config, name, 3));
        
        totalExposurePOT += config.normsRun3.nuWroPOT / (1e20);
    }

    scalingData.exposurePOT = totalExposurePOT;
    scalingData.nTargets = config.global.targetDensity * (1e-8) * AnalysisHelper::GetFiducialVolume();
    scalingDataTruth.exposurePOT = totalExposurePOT;
    scalingDataTruth.nTargets = config.global.targetDensity * (1e-8) * AnalysisHelper::GetFiducialVolume();

    // -------------------------------------------------------------------------------------------------------------------------------------
    // Count the events for CC1pi
    // -------------------------------------------------------------------------------------------------------------------------------------
    // Loop over the files
    for (const auto &[sampleType, sampleName, fileName, normalisation] : inputData)
    {
        std::cout << "Reading input file: " << fileName << std::endl;

        const auto isOverlay = (sampleType == AnalysisHelper::Overlay);
        const auto isDirt    = (sampleType == AnalysisHelper::Dirt);
        const auto isNuWro = (sampleType == AnalysisHelper::DataBNB);
        const auto isDetVar  = (sampleType == AnalysisHelper::DetectorVariation);

        // Open the input file for reading and enable the branches with systematic event weights (if required)
        FileReader reader(fileName);

        if (isOverlay || isNuWro)
            reader.EnableSystematicBranches();
        
        if (isDetVar) continue; // TODO: Remove this after debugging !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        auto pEvent = reader.GetBoundEventAddress();

        // Loop over the events in the file
        const auto nEvents = reader.GetNumberOfEvents();

        //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        for (unsigned int i = 0; i < nEvents; ++i)// TODO REMOVE Every 50th!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
            std::map<std::string, bool> isSelectedMap = {{"generic",isSelectedGeneric},{"golden",isSelectedGolden}};

            // Determine if this is truly a CC1Pi event
            const auto isTrueCC1Pi = (isOverlay || isDetVar || isNuWro) && AnalysisHelper::IsTrueCC1Pi(pEvent, config.global.useAbsPdg);
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

            const auto isCC1PiSignal = isTrueCC1Pi && passesPhaseSpaceTruth;
            const auto isSignal = isCC1PiSignal;

            // -----------------------------------------------------------------------------------------------------------------------------
            // Work out if this event is signal, and apply any phase-space restrictions based on the input binning
            // -----------------------------------------------------------------------------------------------------------------------------

            // Get the nominal event weight, scaled by the sample normalisation
            const auto weight = AnalysisHelper::GetNominalEventWeight(pEvent) * normalisation;
            
            // Determine if this is truly a CC0Pi event
            const auto isTrueCC0Pi = (isOverlay || isDetVar || isNuWro) && AnalysisHelper::IsTrueCC0Pi(pEvent, config.global.useAbsPdg, config.global.protonMomentumThreshold);

            // Get the truth analysis data (if available, otherwise set to dummy values)
            const auto sidebandTruthData = (
                isTrueCC0Pi
                    ? AnalysisHelper::GetTruthAnalysisDataCC0Pi(pEvent->truth, config.global.useAbsPdg, config.global.protonMomentumThreshold)
                    : AnalysisHelper::GetDummyAnalysisData()
            );

            bool passesSidebandPhaseSpaceTruth = false;
            if (isTrueCC0Pi)
            {
                // Start by assuming the event passes the phase-space cuts
                passesSidebandPhaseSpaceTruth = true;

                // Check the value of the kinematic quantities are within the phase-space limits
                for (const auto &[name, minMax] : phaseSpaceMap)
                {
                    const auto &[min, max] = minMax;
                    const auto value = getSidebandValue.at(name)(sidebandTruthData);

                    if (value < min || value > max)
                    {
                        if(name == "pionMomentum") continue; // Not really compatible with the pion momentum in the CC1pi selection
                        passesSidebandPhaseSpaceTruth = false;
                        break;
                    }
                }
            }

            const auto isCC0PiSignal = isTrueCC0Pi && passesSidebandPhaseSpaceTruth;
            
            // -----------------------------------------------------------------------------------------------------------------------------
            // Handle BNB data
            // -----------------------------------------------------------------------------------------------------------------------------
            if (isNuWro)
            {
                std::cout<<"isNuWro true - weight: "<<weight<<" - isCC0PiSignal: "<<isCC0PiSignal<<std::endl;
                for (auto &[selectionName, xsecs] : xsecMap)
                {
                    // Determine if we passed the relevant selection
                    const auto isSelected = isSelectedMap.at(selectionName);
                    // Only count events passing the selection
                    if (!isSelected)
                        continue;
                    std::cout<<"\tselectionName: "<<selectionName<<" selected"<<std::endl;
                    for (auto &[name, xsec] : xsecs)
                    {
                        xsec.AddWeightedSelectedBNBDataEvent(getValue.at(name)(recoData), weight);
                        xsecMapUnscaled.at(selectionName).at(name).AddWeightedSelectedBNBDataEvent(getValue.at(name)(recoData), weight);
                        if(isCC0PiSignal) xsecMapOnlyCC0Pi.at(selectionName).at(name).AddWeightedSelectedBNBDataEvent(getValue.at(name)(recoData), weight);
                    }
                }
                // For BNB data that's all we need to do!
                // continue;// Important: Nuwro files are also used as MC for xSecMapTrue
            }

            // -----------------------------------------------------------------------------------------------------------------------------
            // Handle the detector variation samples (unisims)
            // -----------------------------------------------------------------------------------------------------------------------------
            if (isDetVar)
            {
                // const auto isSignal = isSignalMap.at(selectionName);
                // Handle signal events
                if (isSignal)
                {
                    for (auto &[selectionName, xsecs] : xsecMap)
                    {
                        // Determine if we passed the relevant selection
                        const auto isSelected = isSelectedMap.at(selectionName);

                        for (auto &[name, xsec] : xsecs)
                        {
                            const auto recoValue = getValue.at(name)(recoData);
                            const auto trueValue = getValue.at(name)(truthData);
                            xsec.AddSignalEventDetVar(recoValue, trueValue, isSelected, weight, sampleName);
                            xsecMapUnscaled.at(selectionName).at(name).AddSignalEventDetVar(recoValue, trueValue, isSelected, weight, sampleName);
                            // xsecMapTrue.at(selectionName).at(name).AddSignalEventDetVar(recoValue, trueValue, isSelected, weight, sampleName);
                        }
                    }
                    
                }
                // Handle selected background events
                else
                {
                    for (auto &[selectionName, xsecs] : xsecMap)
                    {
                        // Only use selected background events
                        const auto isSelected = isSelectedMap.at(selectionName);
                        if (!isSelected)
                            continue;

                        for (auto &[name, xsec] : xsecs)
                        {
                            const auto recoValue = getValue.at(name)(recoData);
                            xsec.AddSelectedBackgroundEventDetVar(recoValue, weight, sampleName);
                            xsecMapUnscaled.at(selectionName).at(name).AddSelectedBackgroundEventDetVar(recoValue, weight, sampleName);
                            // xsecMapTrue.at(selectionName).at(name).AddSelectedBackgroundEventDetVar(recoValue, weight, sampleName);
                        }
                    }
                }

                // For detector variation samples, that's all we need to do!
                continue;
            }

            // -----------------------------------------------------------------------------------------------------------------------------
            // Handle all other events (i.e those from the nominal simulation): Overlays, dirt, EXT data
            // -----------------------------------------------------------------------------------------------------------------------------
            // std::cout<<"### Debug Point 0"<<std::endl;
            // Get the flux weights
            auto fluxWeights = (
                isOverlay
                    ? CrossSectionHelper::GetWeightsMap(pEvent->truth, systParams.fluxDimensions, config.extractXSecs.mutuallyExclusiveDimensions)
                    : CrossSectionHelper::GetUnitWeightsMap(systParams.fluxDimensions)
            );

            // Fill the flux-reweightor with all overlay events from a neutrinos of the desired flavour in the active volume
            if (isOverlay && AnalysisHelper::IsInActiveVolume(pEvent->truth.nuVertex()) &&
                std::find(config.flux.nuPdgsSignal.begin(), config.flux.nuPdgsSignal.end(), pEvent->truth.nuPdgCode()) != config.flux.nuPdgsSignal.end())
            {
                scalingData.pFluxReweightor->AddEvent(pEvent->truth.nuEnergy(), weight, fluxWeights);
            }

            // Fill the flux-reweightor with all overlay events from a neutrinos of the desired flavour in the active volume
            if (isNuWro && AnalysisHelper::IsInActiveVolume(pEvent->truth.nuVertex()) &&
                std::find(config.flux.nuPdgsSignal.begin(), config.flux.nuPdgsSignal.end(), pEvent->truth.nuPdgCode()) != config.flux.nuPdgsSignal.end())
            {
                scalingDataTruth.pFluxReweightor->AddEvent(pEvent->truth.nuEnergy(), weight, fluxWeights);
            }

            // std::cout<<"### Debug Point 2"<<std::endl;

            // Get the cross-section weights
            // ATTN here we optionally scale the cross-section weights down by the genieTuneEventWeight - this is done so we don't
            // double count this weight (once in the nominal event weight, and once in the xsec systematic event weights)
            // const auto xsecWeightsScaleFactor = 1.f;
            // std::cout<<"Before xsecWeightsScaleFactor"<<std::endl;
            auto xsecWeightsScaleFactor = (isOverlay && config.extractXSecs.scaleXSecWeights) ? pEvent->truth.genieTuneEventWeight() : 1.f;
            xsecWeightsScaleFactor = std::max(xsecWeightsScaleFactor, 0.0001f);
            // std::cout<<"After xsecWeightsScaleFactor"<<std::endl;

            // std::cout<<"xsecWeightsScaleFactor Debugging Point 1 - xsecWeightsScaleFactor: "<<xsecWeightsScaleFactor<<std::endl;

            auto xsecWeights = (
                isOverlay
                    ? CrossSectionHelper::ScaleWeightsMap(CrossSectionHelper::GetWeightsMap(pEvent->truth, systParams.xsecDimensions, config.extractXSecs.mutuallyExclusiveDimensions), xsecWeightsScaleFactor)
                    : CrossSectionHelper::GetUnitWeightsMap(systParams.xsecDimensions)
            );

            // std::cout<<"xsecWeightsScaleFactor Debugging Point 2"<<std::endl;

            // Get the reinteraction weights
            auto reintWeights = (
                isOverlay
                    ? CrossSectionHelper::GetWeightsMap(pEvent->truth, systParams.reintDimensions, config.extractXSecs.mutuallyExclusiveDimensions)
                    : CrossSectionHelper::GetUnitWeightsMap(systParams.reintDimensions)
            );

            // Handle signal events
            if (isSignal)
            {
                for (auto &[selectionName, xsecs] : xsecMap)
                {
                    // Determine if we passed the relevant selection
                    const auto isSelected = isSelectedMap.at(selectionName);

                    for (auto &[name, xsec] : xsecs)
                    {
                        const auto recoValue = getValue.at(name)(recoData);
                        const auto trueValue = getValue.at(name)(truthData);

                        const auto seedString =  selectionName + name + std::to_string(i);
                        if(!isNuWro)
                        {
                            xsec.AddSignalEvent(recoValue, trueValue, isSelected, weight, fluxWeights, xsecWeights, reintWeights, seedString);
                            xsecMapUnscaled.at(selectionName).at(name).AddSignalEvent(recoValue, trueValue, isSelected, weight, fluxWeights, xsecWeights, reintWeights, seedString);
                        } 
                        else
                        {
                            xsecMapTrue.at(selectionName).at(name).AddSignalEvent(recoValue, trueValue, isSelected, weight, fluxWeights, xsecWeights, reintWeights, seedString);
                        }
                    }
                }
            }
            // Handle selected background events
            else
            {
                // std::cout<<"xsecWeightsScaleFactor Debugging Point 4"<<std::endl;
                for (auto &[selectionName, xsecs] : xsecMap)
                {
                    // Only use selected background events
                    const auto isSelected = isSelectedMap.at(selectionName);
                    if (!isSelected)
                        continue;

                    for (auto &[name, xsec] : xsecs)
                    {
                        const auto recoValue = getValue.at(name)(recoData);
                        const auto trueValue = getValue.at(name)(truthData);
                        
                        // Add CC0pi constraint
                        const auto trueSidebandValue = getSidebandValue.at(name)(sidebandTruthData);
                        const auto cc0piNominalConstraintParam = cc0piNominalConstraintMap.at("generic").at(name).first;
                        const auto cc0piNominalConstraintParamError = cc0piNominalConstraintMap.at("generic").at(name).second;
                        const auto cc0piUniverseConstraintVector = cc0piUniverseConstraintMap.at("generic").at(name);

                        const auto scaledWeight = (
                            (isCC0PiSignal)
                                ? weight*xsec.GetSidebandScaling(trueSidebandValue, cc0piNominalConstraintParam)
                                : weight);

                        const auto scaledXSecWeights = (
                            (isCC0PiSignal)
                                ? xsec.GetSidebandUniverseScaling(xsecWeights, trueSidebandValue, cc0piNominalConstraintParam, cc0piUniverseConstraintVector)
                                : xsecWeights);

                        const auto scaledFluxWeights = (
                            (isCC0PiSignal)
                                ? xsec.GetSidebandUniverseScaling(fluxWeights, trueSidebandValue, cc0piNominalConstraintParam, cc0piUniverseConstraintVector)
                                : fluxWeights);

                        const auto scaledReintWeights = (
                            (isCC0PiSignal)
                                ? xsec.GetSidebandUniverseScaling(reintWeights, trueSidebandValue, cc0piNominalConstraintParam, cc0piUniverseConstraintVector)
                                : reintWeights);

                        const auto sidebandWeights = ( //Todo: arbitrary threshold in function. Remove??? // Also change from vector<float> to SystFloatMap for consistency
                            (isCC0PiSignal)
                                ? xsec.GetSidebandParameterWeights(trueSidebandValue, cc0piNominalConstraintParam, cc0piNominalConstraintParamError)
                                : std::vector<float>(systParams.nBootstrapUniverses, 1.0));


                        const auto seedString =  selectionName + name + std::to_string(i);

                        if(!isNuWro)
                        {
                            xsec.AddSelectedBackgroundEvent(recoValue, isDirt, scaledWeight, scaledFluxWeights, scaledXSecWeights, scaledReintWeights, sidebandWeights, seedString);
                            xsecMapUnscaled.at(selectionName).at(name).AddSelectedBackgroundEvent(recoValue, isDirt, weight, fluxWeights, xsecWeights, reintWeights, std::vector<float>(systParams.nBootstrapUniverses, 1.0), seedString);
                            if(isCC0PiSignal)
                            {
                                xsecMapOnlyCC0Pi.at(selectionName).at(name).AddSelectedBackgroundEvent(recoValue, isDirt, scaledWeight, scaledFluxWeights, scaledXSecWeights, scaledReintWeights, sidebandWeights, seedString);
                                xsecMapOnlyCC0PiUnscaled.at(selectionName).at(name).AddSelectedBackgroundEvent(recoValue, isDirt, weight, fluxWeights, xsecWeights, reintWeights, std::vector<float>(systParams.nBootstrapUniverses, 1.0), seedString);
                            }
                        }
                        else
                        {
                            xsecMapTrue.at(selectionName).at(name).AddSelectedBackgroundEvent(recoValue, isDirt, scaledWeight, scaledFluxWeights, scaledXSecWeights, scaledReintWeights, sidebandWeights, seedString);
                        }
                    }
                }
            }
        }
    }


    // Loop over all cross-section objects
    for(const auto&[postfix, map]: {std::make_pair("scaled",xsecMap),std::make_pair("unscaled", xsecMapUnscaled)})
    {
        for (const auto &[selectionName, xsecs] : map)
        {
            for (const auto &[name, xsec] : xsecs)
            {
                std::cout << "Processing cross-section: "<<selectionName<< " - " << name << std::endl;

                // -----------------------------------------------------------------------------------------------------------------------------
                // Get the event rates for BNB data, backgrounds, and signal
                // -----------------------------------------------------------------------------------------------------------------------------
                const auto selectedEventsData = xsec.GetSelectedBNBDataEvents();
                std::cout << "Selected BNB data events" << std::endl;
                FormattingHelper::SaveMatrix(selectedEventsData, "xsecNuWro_" + selectionName + "_" + name + "_data_selected_eventRate_"+postfix+".txt");

                const auto selectedEventsBackground = xsec.GetSelectedBackgroundEvents();
                std::cout << "Selected background events" << std::endl;
                FormattingHelper::SaveMatrix(selectedEventsBackground, "xsecNuWro_" + selectionName + "_" + name + "_background_selected_eventRate_"+postfix+".txt");

                const auto selectedEventsSignal = xsec.GetSelectedSignalEvents();
                std::cout << "Selected signal events" << std::endl;
                FormattingHelper::SaveMatrix(selectedEventsSignal, "xsecNuWro_" + selectionName + "_" + name + "_signal_selected_eventRate_"+postfix+".txt");

                const auto allEventsSignal = xsec.GetSignalEvents();
                std::cout << "All signal events" << std::endl;
                FormattingHelper::SaveMatrix(allEventsSignal, "xsecNuWro_" + selectionName + "_" + name + "_signal_all_eventRate_"+postfix+".txt");
                std::cout << "After all signal events" << std::endl;
                // -----------------------------------------------------------------------------------------------------------------------------
                // Get the cross-section as measured with BNB data along with it's uncertainties
                // -----------------------------------------------------------------------------------------------------------------------------
                const auto data = xsec.GetBNBDataCrossSection(scalingData);
                std::cout << "BNB data cross-section (reco-space)" << std::endl;
                FormattingHelper::SaveMatrix(data, "xsecNuWro_" + selectionName + "_" + name + "_data_"+postfix+".txt");

                const auto dataStatUncertainties = xsec.GetBNBDataCrossSectionStatUncertainty(scalingData);
                std::cout << "BNB data stat uncertainty" << std::endl;
                FormattingHelper::SaveMatrix(dataStatUncertainties, "xsecNuWro_" + selectionName + "_" + name + "_data_stat_"+postfix+".txt");

                const auto dataSystBiasCovariances = xsec.GetBNBDataCrossSectionSystUncertainties(scalingData);
                for (const auto &[group, map] : dataSystBiasCovariances)
                {
                    for (const auto &[paramName, biasCovariance] : map)
                    {
                        const auto &[pBias, pCovariance] = biasCovariance;

                        std::cout << "BNB data syst uncertainty: " << group << " " << paramName << std::endl;
                        std::cout << "Bias vector" << std::endl;
                        FormattingHelper::SaveMatrix(*pBias, "xsecNuWro_" + selectionName + "_" + name + "_data_" + group + "_" + paramName + "_bias_"+postfix+".txt");
                        std::cout << "Covariance matrix" << std::endl;
                        FormattingHelper::SaveMatrix(*pCovariance, "xsecNuWro_" + selectionName + "_" + name + "_data_" + group + "_" + paramName + "_covariance_"+postfix+".txt");
                    }
                }

                // -----------------------------------------------------------------------------------------------------------------------------
                // Get the predicted cross-section along with its MC stat uncertainty
                // -----------------------------------------------------------------------------------------------------------------------------
                const auto prediction = xsec.GetPredictedCrossSection(scalingData);
                std::cout << "Predicted cross-section (truth-space)" << std::endl;
                FormattingHelper::SaveMatrix(prediction, "xsecNuWro_" + selectionName + "_" + name + "_prediction_"+postfix+".txt");

                const auto &[pPredictionStatBias, pPredictionStatCovariance] = xsec.GetPredictedCrossSectionStatUncertainty(scalingData);
                std::cout << "Predicted cross-section stat uncertainty" << std::endl;
                std::cout << "Bias vector" << std::endl;
                FormattingHelper::SaveMatrix(*pPredictionStatBias, "xsecNuWro_" + selectionName + "_" + name + "_prediction_stat_bias_"+postfix+".txt");
                std::cout << "Covariance matrix" << std::endl;
                FormattingHelper::SaveMatrix(*pPredictionStatCovariance, "xsecNuWro_" + selectionName + "_" + name + "_prediction_stat_covariance_"+postfix+".txt");


                const auto &[pPredictionSidebandStatBias, pPredictionSidebandStatCovariance] = xsec.GetPredictedSidebandCrossSectionStatUncertainty(scalingData);
                std::cout << "Predicted cross-section stat uncertainty" << std::endl;
                std::cout << "Bias vector" << std::endl;
                FormattingHelper::SaveMatrix(*pPredictionSidebandStatBias, "xsecNuWro_" + selectionName + "_" + name + "_prediction_sideband_stat_bias_"+postfix+".txt");
                std::cout << "Covariance matrix" << std::endl;
                FormattingHelper::SaveMatrix(*pPredictionSidebandStatCovariance, "xsecNuWro_" + selectionName + "_" + name + "_prediction_sideband_stat_covariance_"+postfix+".txt");

                // -----------------------------------------------------------------------------------------------------------------------------
                // Get the smearing matrix along with its uncertainties
                // -----------------------------------------------------------------------------------------------------------------------------
                std::cout << "Smearing Matrix (reco-space rows, truth-space columns)" << std::endl;
                const auto smearingMatrix = xsec.GetSmearingMatrix();
                FormattingHelper::SaveMatrix(smearingMatrix, "xsecNuWro_" + selectionName + "_" + name + "_smearingMatrix_"+postfix+".txt");

                std::cout << "Smearing Matrix AllSelected" << std::endl;
                const auto smearingMatrixAllSelected = xsec.GetSmearingMatrixAllSelected();
                FormattingHelper::SaveMatrix(smearingMatrix, "xsecNuWro_" + selectionName + "_" + name + "_smearingMatrixAllSelected_"+postfix+".txt");

                std::cout << "Smearing Matrix SystBiasCovariances" << std::endl;
                const auto smearingMatrixSystBiasCovariances = xsec.GetSmearingMatrixSystUncertainties();
                std::cout << "Smearing Matrix SystBiasCovariances - After" << std::endl;
                for (const auto &[group, map] : smearingMatrixSystBiasCovariances)
                {
                    for (const auto &[paramName, biasCovariance] : map)
                    {
                        const auto &[pBias, pCovariance] = biasCovariance;

                        std::cout << "Smearing matrix syst uncertainty: " << group << " " << paramName << std::endl;
                        std::cout << "Bias vector" << std::endl;
                        FormattingHelper::SaveMatrix(*pBias, "xsecNuWro_" + selectionName + "_" + name + "_smearingMatrix_" + group + "_" + paramName + "_bias_"+postfix+".txt");
                        std::cout << "Covariance matrix" << std::endl;
                        FormattingHelper::SaveMatrix(*pCovariance, "xsecNuWro_" + selectionName + "_" + name + "_smearingMatrix_" + group + "_" + paramName + "_covariance_"+postfix+".txt");
                    }
                }
            }
        }
    }

    // Loop over all cross-section objects
    for (const auto &[selectionName, xsecs] : xsecMapTrue)
    {
        for (const auto &[name, xsec] : xsecs)
        {
            std::cout << "True: Processing cross-section: "<<selectionName<< " - " << name << std::endl;

            // -----------------------------------------------------------------------------------------------------------------------------
            // Get the event rates for BNB data, backgrounds, and signal
            // -----------------------------------------------------------------------------------------------------------------------------
            // const auto selectedEventsData = xsec.GetSelectedBNBDataEvents();
            // std::cout << "True: Selected BNB data events" << std::endl;
            // FormattingHelper::SaveMatrix(selectedEventsData, "xsecNuWro_true_" + selectionName + "_" + name + "_data_selected_eventRate.txt");

            // -----------------------------------------------------------------------------------------------------------------------------
            // Get the cross-section as measured with BNB data along with it's uncertainties
            // -----------------------------------------------------------------------------------------------------------------------------
            // const auto data = xsec.GetNuWroTrueCrossSection(scalingData);
            // std::cout << "True: BNB data cross-section (reco-space)" << std::endl;
            // FormattingHelper::SaveMatrix(data, "xsecNuWro_true_" + selectionName + "_" + name + "_data.txt");

            const auto prediction = xsec.GetPredictedCrossSection(scalingDataTruth);
            std::cout << "Predicted NuWro cross-section (truth-space)" << std::endl;
            FormattingHelper::SaveMatrix(prediction, "xsecNuWro_true_" + selectionName + "_" + name + "_prediction.txt");

            std::cout << "NuWro smearing Matrix (reco-space rows, truth-space columns)" << std::endl;
            const auto smearingMatrix = xsec.GetSmearingMatrix();
            FormattingHelper::SaveMatrix(smearingMatrix, "xsecNuWro_true_" + selectionName + "_" + name + "_smearingMatrix.txt");

            std::cout << "NuWro smearing Matrix all selected (reco-space rows, truth-space columns)" << std::endl;
            const auto smearingMatrixAllSelected = xsec.GetSmearingMatrixAllSelected();
            FormattingHelper::SaveMatrix(smearingMatrixAllSelected, "xsecNuWro_true_" + selectionName + "_" + name + "_smearingMatrixAllSelected.txt");
        }
    }


        // Loop over all cross-section objects
    for (const auto &[selectionName, xsecs] : xsecMapOnlyCC0Pi)
    {
        for (const auto &[name, xsec] : xsecs)
        {
            std::cout << "Only selected CC0pi: Processing cross-section: "<<selectionName<< " - " << name << std::endl;

            // -----------------------------------------------------------------------------------------------------------------------------
            // Get the event rates for BNB data, backgrounds, and signal
            // -----------------------------------------------------------------------------------------------------------------------------
            const auto selectedEventsData = xsec.GetSelectedBNBDataEvents();
            std::cout << "True: Selected BNB data events" << std::endl;
            FormattingHelper::SaveMatrix(selectedEventsData, "xsecNuWro_onlyCC0pi_" + selectionName + "_" + name + "_data_selected_eventRate.txt");

            const auto selectedEventsBackground = xsec.GetSelectedBackgroundEvents();
            std::cout << "Selected background events" << std::endl;
            FormattingHelper::SaveMatrix(selectedEventsBackground, "xsecNuWro_onlyCC0pi_" + selectionName + "_" + name + "_background_selected_eventRate_scaled.txt");

            const auto selectedEventsBackgroundUnscaled = xsecMapOnlyCC0PiUnscaled.at(selectionName).at(name).GetSelectedBackgroundEvents();
            std::cout << "Selected background events" << std::endl;
            FormattingHelper::SaveMatrix(selectedEventsBackgroundUnscaled, "xsecNuWro_onlyCC0pi_" + selectionName + "_" + name + "_background_selected_eventRate_unscaled.txt");
        }
    }

    std::cout<<"------------- All Done -------------"<<std::endl;
    return;
}

} // namespace ubcc1pi_macros
