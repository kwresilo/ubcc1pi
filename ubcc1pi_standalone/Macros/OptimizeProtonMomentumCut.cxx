// /**
//  *  @file  ubcc1pi_standalone/Macros/OptimizeProtonMomentumCut.cxx
//  *
//  *  @brief The implementation file of the OptimizeProtonMomentumCut macro
//  */

// #include "ubcc1pi_standalone/Macros/Macros.h"
// #include "ubcc1pi_standalone/Objects/FileReader.h"

// #include "ubcc1pi_standalone/Helpers/AnalysisHelper.h"
// #include "ubcc1pi_standalone/Helpers/NormalisationHelper.h"
// #include "ubcc1pi_standalone/Helpers/SelectionHelper.h"
// #include "ubcc1pi_standalone/Helpers/CrossSectionHelper.h"
// #include "ubcc1pi_standalone/Helpers/FormattingHelper.h"
// // #include "ubcc1pi_standalone/Helpers/FittingHelper.h"
// #include "ubsmear.h"

// #include <fstream> // Todo: not use txt files
// // Boost libraries
// #include "binary_iarchive.hpp"
// // #include "binary_oarchive.hpp"
// #include "binary_object.hpp"
// #include "map.hpp"
// #include "vector.hpp"

// using namespace ubcc1pi;

// namespace ubcc1pi_macros
// {

// void OptimizeProtonMomentumCut(const Config &config)
// {

//     // -------------------------------------------------------------------------------------------------------------------------------------
//     // Setup an object that holds the details of the systematic parameters to apply
//     // -------------------------------------------------------------------------------------------------------------------------------------
//     // Here we read in the "dimensions" of the systematic parameters to apply. For multisim parameters (flux & xsec), this is a map from
//     // each parameter name to the expected number of universes. For unisim parameters (detector variations), this is a map from the
//     // identidier (i.e. name) of the detector variation sample, to the identifier of the corresponding central-value sample. Additionally
//     // we set the number of bootstrap universes (for the MC stat uncertainty), the corresponding weights are generated for each event.
//     CrossSectionHelper::CrossSection::SystParams systParams;
//     systParams.nBootstrapUniverses = config.extractXSecs.nBootstrapUniverses;
//     systParams.fluxDimensions = config.extractXSecs.fluxDimensions;
//     systParams.xsecDimensions = config.extractXSecs.xsecDimensions;
//     systParams.reintDimensions = config.extractXSecs.reintDimensions;
//     systParams.detVarDimensions = config.extractXSecs.detVarDimensions;
//     systParams.potFracUncertainty = config.extractXSecs.potFracUncertainty;

//     // -------------------------------------------------------------------------------------------------------------------------------------
//     // Setup an object that holds the information about how we should scale an event rate to obtain a cross-section
//     // -------------------------------------------------------------------------------------------------------------------------------------
//     // Here we specify the:
//     // - Flux             [10^-10 cm^-2 POT^-1]
//     // - Exposure POT     [10^20 POT]               (stored in config as [POT])
//     // - Target density   [10^31 nucleons/cm^3]     (stored in config as [10^23 nucleons/cm^3])
//     // - Fiducial volume  [cm^3]
//     //
//     // Hence the units of the eventual cross-seciton are:
//     // - Cross-section    [Flux * Exposure * Target density * Fiducial volume]^-1 = [10^-41 cm^2 / nucleon]
//     //
//     // Here we use a FluxReweightor to specify the flux. For each universe of each flux systematic paramter, this is uses the ratio of the
//     // total neutrino event rate (as a function of the true neutrino energy) to the same rate in the nominal simulation to reweight the
//     // input neutrino flux distribution. The integrated flux is each universe is used to scale the selected event rate when calculating the
//     // cross-section in that universe. For all non-flux parameters, the nominal integrated flux is used.
//     CrossSectionHelper::CrossSection::ScalingData scalingData;
//     CrossSectionHelper::CrossSection::ScalingData scalingDataNuWroTruth;

//     // std::cout << "- Flux:            " << scalingData.pFluxReweightor->GetIntegratedNominalFlux() << " * 10^-10 cm^-2 POT^-1" << std::endl;
//     // std::cout << "- Exposure:        " << scalingData.exposurePOT << " * 10^20 POT" << std::endl;
//     // std::cout << "- Target density:  " << config.global.targetDensity << " * 10^31 nucleons/cm^3" << std::endl;
//     // std::cout << "- Fiducial volume: " << AnalysisHelper::GetFiducialVolume() << " * cm^3" << std::endl;
//     // std::cout << "- nTargets:        " << scalingData.nTargets << " * 10^31 nucleons" << std::endl;

//     // -------------------------------------------------------------------------------------------------------------------------------------
//     // Setup the event selection
//     // -------------------------------------------------------------------------------------------------------------------------------------
//     // Here we are using the default CC1pi selection. The final cut of this selection (named "likelyGoldenPion") is used to identify events
//     // in which the reconstructed & selected pion candidate is "golden". Here a golden pion is one that is: contained within the TPC,
//     // doesn't undergo any scatters on the argon, and comes to rest (before any secondary interaction). If an event passes all the cuts of
//     // this selection, it's said to pass the "golden" selection. If an event passes the penultimate cut, it is said to have passed the
//     // "generic" selection. Events passing the generic selection are likely to be CC1pi, but may or may not have a golden pion. All events
//     // that pass the golden selection also pass the generic selection - but not all events that pass the generic selection also pass the
//     // golden selection.
//     auto selection = SelectionHelper::GetDefaultSelection();


//     // The list of proton momentum thresholds to try 
//     const std::vector<float> protonMomentumThresholdVector = {0.25f, 0.26f, 0.27f, 0.28f, 0.29f, 0.3f, 0.31f, 0.32f, 0.33f, 0.34f, 0.35f};

//     // -------------------------------------------------------------------------------------------------------------------------------------
//     // Setup the cross-section objects
//     // -------------------------------------------------------------------------------------------------------------------------------------
//     // Here we make a map from a name of the cross-section to the cross-section object itself. In this way, we can iterate through the
//     // cross-section objects and reduce code-bloat. The first index is an identifier for the selection that's applied (generic or goldlen),
//     // the second index is an identifier for the kinematic quantity that's relevant for the cross-section (e.g. muonMomentum), and the
//     // mapped type is the cross-section object.
//     typedef std::map<std::string, std::map<std::string, CrossSectionHelper::CrossSection>> crossSectionMap;
//     // std::map<std::string, std::map<std::string, CrossSectionHelper::CrossSection>> xsecMapBNB;
//     // std::map<std::string, std::map<std::string, CrossSectionHelper::CrossSection>> xsecMapNuWro;
//     // std::map<std::string, std::map<std::string, CrossSectionHelper::CrossSection>> xsecMapNuWroTrue; // True NuWro cross-sections

//     std::vector<crossSectionMap> xsecMapBNBVector(protonMomentumThresholdVector.size());
//     std::vector<crossSectionMap> xsecMapNuWroVector(protonMomentumThresholdVector.size());
//     std::vector<crossSectionMap> xsecMapNuWroTrueVector(protonMomentumThresholdVector.size()); // True NuWro cross-sections

//     // We additionally make a map from each cross-section to the limits of the phase-space that we should consider. The key is the
//     // identifier for the kinematic quantity and the mapped value is a pair containing the limits [min, max]
//     std::map< std::string, std::pair<float, float> > phaseSpaceMap;

//     // ATTN the configuration allows the user to enable or disable each cross-section. If a cross-section has been disabled, then it won't
//     // be added the xsecMapNuWro. However, the phaseSpaceMap always includes all kinematic parameters!


//     // Add the differential cross-sections
//     for (const auto &[name, binning, scaleByBinWidth] : std::vector< std::tuple<std::string, Config::Global::Binning, bool> > {

//         // The names of the cross-section kinematic parameters, and their binning information.
//         // The third (boolean) parameter indicates if the cross-section bins should be scaled by their width
//         { "muonCosTheta",  config.global.muonCosTheta,  true  },
//         { "muonPhi",       config.global.muonPhi,       true  },
//         { "muonMomentum",  config.global.muonMomentum,  true  },

//         { "pionCosTheta",  config.global.pionCosTheta,  true  },
//         { "pionPhi",       config.global.pionPhi,       true  },
//         { "pionMomentum",  config.global.pionMomentum,  true  },

//         { "muonPionAngle", config.global.muonPionAngle, true  },
//         { "nProtons",      config.global.nProtons,      false }

//     })
//     {

//         // Add to the phase-space map
//         phaseSpaceMap.emplace(name, std::pair<float, float>({binning.min, binning.max}));


        
//         if(name != "nProtons") continue; // Only add the nProtons cross-section

//         // Here we calculate every cross-section using both the generic and golden selection. In the end we only use the golden selection for
//         // the pion momentum, but we additionally apply it to other cross-sections as a cross-check.
//         // Add the cross-section object to the map using the binning from the input configuration
//         const auto &[extendedBinEdges, hasUnderflow, hasOverflow] = CrossSectionHelper::GetExtendedBinEdges(binning.min, binning.max, binning.binEdges);
        

//         for (const auto &selectionName : {"generic", "golden"})
//         {
//             // Don't setup a cross-section object if it's been disabled in the configuration
//             if (config.extractXSecs.crossSectionIsEnabled.at(selectionName).at(name))
//             {

//                 for(unsigned int n=0; n<protonMomentumThresholdVector.size(); n++)
//                 {

//                     xsecMapBNBVector.at(n)[selectionName].emplace(name, CrossSectionHelper::CrossSection(systParams, extendedBinEdges, hasUnderflow, hasOverflow, scaleByBinWidth));
//                     xsecMapNuWroVector.at(n)[selectionName].emplace(name, CrossSectionHelper::CrossSection(systParams, extendedBinEdges, hasUnderflow, hasOverflow, scaleByBinWidth));
//                     xsecMapNuWroTrueVector.at(n)[selectionName].emplace(name, CrossSectionHelper::CrossSection(systParams, extendedBinEdges, hasUnderflow, hasOverflow, scaleByBinWidth));

//                 }
//             }
//         }
//     }


//     // ATTN here we use the machinary for a differential cross-section, and treat the total cross-section as a single-bin measurement.
//     // The "kinematic quantity" in this case is just a dummy parameter. Here we define a single bin with edges arbitrarily chosen to be
//     // (-1 -> +1), and we request that the cross-section object does not apply bin-width scaling. When we fill this object, we will use
//     // the dummy kinematic quantity with a value of 0 for all events. This is arbitrary, as long as it's within the bin edges we chose.
//     // In this way the single bin contains all events. It's just a trick to avoid implementing extra logic for the total cross-section.
//     // for (const auto &selectionName : {"generic", "golden"})
//     // {
//     //     // Don't setup a cross-section object if it's been disabled in the configuration
//     //     if (config.extractXSecs.crossSectionIsEnabled.at(selectionName).at("total"))
//     //     {
//     //         for(unsigned int n=0; n<protonMomentumThresholdVector.size(); n++)
//     //         {

//     //             xsecMapBNBVector.at(n).at(selectionName).emplace("total", CrossSectionHelper::CrossSection(systParams, {-1.f, 1.f}, false, false, false));
//     //             xsecMapNuWroVector.at(n).at(selectionName).emplace("total", CrossSectionHelper::CrossSection(systParams, {-1.f, 1.f}, false, false, false));
//     //             xsecMapNuWroTrueVector.at(n).at(selectionName).emplace("total", CrossSectionHelper::CrossSection(systParams, {-1.f, 1.f}, false, false, false));
//     //         }
//     //     }
//     // }


//     // The dummy value that will be used as the "kinematic quantity" for the total cross-section
//     const auto dummyValue = 0.f;

//     // Check to see if we have any cross-sections enabled
//     if (xsecMapNuWroVector.empty())
//     {
//         std::cout << "All cross-sections have been disabled in the configuration! Nothing more to do" << std::endl;
//         return;
//     }

//     // Print the names of the cross-sections we are going to extract
//     std::cout << "The following cross-sections are enabled:" << std::endl;
//     for (const auto &[selectionName, xsecs] : xsecMapNuWroVector.at(0))
//     {
//         std::cout << "  - Selection: " << selectionName << std::endl;
//         for (const auto &entry : xsecs)
//         {
//             const auto &name = entry.first;
//             std::cout << "    - " << name << std::endl;
//         }
//     }

//     // -------------------------------------------------------------------------------------------------------------------------------------
//     // Setup the relevent "getters" for each cross-section
//     // -------------------------------------------------------------------------------------------------------------------------------------
//     // Here we define a map from the name of each cross-section to a function which pulls out the relevant kinematic quanitity from an input
//     // analysis data object. Again this is done up-front to reduce code-bloat below.
//     std::unordered_map< std::string, std::function<float(const AnalysisHelper::AnalysisData &)> > getValue;
//     // Differential cross-section kinematic parameters
//     getValue.emplace("muonCosTheta",  [](const auto &data) { return data.muonCosTheta;  });
//     getValue.emplace("muonPhi",       [](const auto &data) { return data.muonPhi;       });
//     getValue.emplace("muonMomentum",  [](const auto &data) { return data.muonMomentum;  });
//     getValue.emplace("pionCosTheta",  [](const auto &data) { return data.pionCosTheta;  });
//     getValue.emplace("pionPhi",       [](const auto &data) { return data.pionPhi;       });
//     getValue.emplace("pionMomentum",  [](const auto &data) { return data.pionMomentum;  });
//     getValue.emplace("muonPionAngle", [](const auto &data) { return data.muonPionAngle; });
//     getValue.emplace("nProtons",      [](const auto &data) { return data.nProtons;      });

//     // -------------------------------------------------------------------------------------------------------------------------------------
//     // Setup the input files
//     // -------------------------------------------------------------------------------------------------------------------------------------
//     // Here we define a vector of tuples with 4 entries
//     //   - First, the sample type (e.g. overlay)
//     //   - Second, a string which is used to identify a given detector variation sample (for other sample type, this is unused)
//     //   - Third, the path to the input file
//     //   - Fourth, the normalisation factor to apply to all events in that file
//     const auto fluxHistNames = CrossSectionHelper::GetNominalFluxHistNames(config.flux.nuPdgsSignal, config.flux.nuPdgToHistName, config.flux.nomHistPattern);
//     const auto &[fluxBinEdges, fluxValues] = CrossSectionHelper::ReadNominalFlux(config.flux.fileName, fluxHistNames, config.flux.pot);
//     scalingData.pFluxReweightor = std::make_shared<CrossSectionHelper::FluxReweightor>(fluxBinEdges, fluxValues, systParams.fluxDimensions);
//     scalingDataNuWroTruth.pFluxReweightor = std::make_shared<CrossSectionHelper::FluxReweightor>(fluxBinEdges, fluxValues, systParams.fluxDimensions);
//     std::vector< std::tuple<AnalysisHelper::SampleType, std::string, std::string, float> > inputData;
//     auto totalExposurePOT = 0.f;

//     if(std::binary_search(config.global.runs.begin(), config.global.runs.end(), 1))
//     {
//         inputData.emplace_back(AnalysisHelper::Overlay, "", config.filesRun1.overlaysFileName, NormalisationHelper::GetOverlaysNormalisation(config, 1));
//         inputData.emplace_back(AnalysisHelper::Dirt,    "", config.filesRun1.dirtFileName, NormalisationHelper::GetDirtNormalisation(config, 1));
//         inputData.emplace_back(AnalysisHelper::DataEXT, "", config.filesRun1.dataEXTFileName, NormalisationHelper::GetDataEXTNormalisation(config, 1));
//         inputData.emplace_back(AnalysisHelper::NuWro,   "", config.filesRun1.nuWroFileName, NormalisationHelper::GetNuWroNormalisation(config, 1));
//         inputData.emplace_back(AnalysisHelper::DataBNB, "", config.filesRun1.dataBNBFileName, 1.f);
//         // Add the detector variation files
//         for (const auto &[name, fileName] : config.filesRun1.detVarFiles)
//         {
//             inputData.emplace_back(AnalysisHelper::DetectorVariation, name, fileName, NormalisationHelper::GetDetectorVariationNormalisation(config, name, 1));
//         }
//         totalExposurePOT += config.normsRun1.dataBNBTor875WCut / (1e20);
//     }
//     if(std::binary_search(config.global.runs.begin(), config.global.runs.end(), 2))
//     {
//         inputData.emplace_back(AnalysisHelper::Overlay, "", config.filesRun2.overlaysFileName, NormalisationHelper::GetOverlaysNormalisation(config, 2));
//         inputData.emplace_back(AnalysisHelper::Dirt,    "", config.filesRun2.dirtFileName, NormalisationHelper::GetDirtNormalisation(config, 2));
//         inputData.emplace_back(AnalysisHelper::DataEXT, "", config.filesRun2.dataEXTFileName, NormalisationHelper::GetDataEXTNormalisation(config, 2));
//         inputData.emplace_back(AnalysisHelper::NuWro,   "", config.filesRun2.nuWroFileName, NormalisationHelper::GetNuWroNormalisation(config, 2));
//         inputData.emplace_back(AnalysisHelper::DataBNB, "", config.filesRun2.dataBNBFileName, 1.f);
//         // Add the detector variation files
//         for (const auto &[name, fileName] : config.filesRun2.detVarFiles)
//         {
//             inputData.emplace_back(AnalysisHelper::DetectorVariation, name, fileName, NormalisationHelper::GetDetectorVariationNormalisation(config, name, 2));
//         }
//         totalExposurePOT += config.normsRun2.dataBNBTor875WCut / (1e20);
//     }
//     if(std::binary_search(config.global.runs.begin(), config.global.runs.end(), 3))
//     {
//         inputData.emplace_back(AnalysisHelper::Overlay, "", config.filesRun3.overlaysFileName, NormalisationHelper::GetOverlaysNormalisation(config, 3));
//         inputData.emplace_back(AnalysisHelper::Dirt,    "", config.filesRun3.dirtFileName, NormalisationHelper::GetDirtNormalisation(config, 3));
//         inputData.emplace_back(AnalysisHelper::DataEXT, "", config.filesRun3.dataEXTFileName, NormalisationHelper::GetDataEXTNormalisation(config, 3));
//         inputData.emplace_back(AnalysisHelper::NuWro,   "", config.filesRun3.nuWroFileName, NormalisationHelper::GetNuWroNormalisation(config, 3));
//         inputData.emplace_back(AnalysisHelper::DataBNB, "", config.filesRun3.dataBNBFileName, 1.f);

//         // Add the detector variation files
//         for (const auto &[name, fileName] : config.filesRun3.detVarFiles)
//         {
//             inputData.emplace_back(AnalysisHelper::DetectorVariation, name, fileName, NormalisationHelper::GetDetectorVariationNormalisation(config, name, 3));
//         }
//         totalExposurePOT += config.normsRun3.dataBNBTor875WCut / (1e20);
//     }

//     scalingData.exposurePOT = totalExposurePOT;
//     scalingData.nTargets = config.global.targetDensity * (1e-8) * AnalysisHelper::GetFiducialVolume();
//     scalingDataNuWroTruth.exposurePOT = totalExposurePOT; //todo check this is correct and shouldn't use nuwro values
//     scalingDataNuWroTruth.nTargets = config.global.targetDensity * (1e-8) * AnalysisHelper::GetFiducialVolume();

//     std::cout<<"Testing proton momentum thresholds: ";
//     for (const auto &protonMomentumThreshold : protonMomentumThresholdVector) std::cout<<protonMomentumThreshold<<" ";
//     std::cout<<std::endl;

//     // -------------------------------------------------------------------------------------------------------------------------------------
//     // Count the events for CC1pi
//     // -------------------------------------------------------------------------------------------------------------------------------------
//     // Loop over the files
//     for (const auto &[sampleType, sampleName, fileName, normalisation] : inputData)
//     {
//         std::cout << "Reading input file: " << fileName << std::endl;

//         const auto isOverlay = (sampleType == AnalysisHelper::Overlay);
//         const auto isDirt    = (sampleType == AnalysisHelper::Dirt);
//         const auto isNuWro   = (sampleType == AnalysisHelper::NuWro);
//         const auto isDataBNB = (sampleType == AnalysisHelper::DataBNB);
//         const auto isDetVar  = (sampleType == AnalysisHelper::DetectorVariation);
//         const auto isDataEXT = (sampleType == AnalysisHelper::DataEXT);

//         if(isDataBNB || isDataEXT || isDirt) continue; // todo completely remove bnb // bnb and nuwro should give same values as this cut depends on MC (overlay) only

//         // Open the input file for reading and enable the branches with systematic event weights (if required)
//         FileReader reader(fileName);

//         if (isOverlay || isNuWro)
//             reader.EnableSystematicBranches();
        
//         // if (isDetVar) continue; // TODO: Remove this after debugging !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//         // if(isDataBNB || isDirt || isDataEXT) continue; // Todo: Remove !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

//         auto pEvent = reader.GetBoundEventAddress();

//         // Loop over the events in the file
//         const auto nEvents = reader.GetNumberOfEvents();

//         //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//         //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//         //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//         //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//         //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//         //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//         // std::cout<<"\n##############\nOnly counting every 30th event!\n##############"<<std::endl;
//         for (unsigned int i = 0; i < nEvents; ++i) // todo change back to every event!!!!!!!!!!!!!!!!!!!!!!
//         {
//             AnalysisHelper::PrintLoadingBar(i, nEvents);
//             reader.LoadEvent(i);
//             // -----------------------------------------------------------------------------------------------------------------------------
//             // Work out if this event passed the selection and apply any additional phase-space cuts based on the input binning
//             // -----------------------------------------------------------------------------------------------------------------------------
//             // Run the selection
//             const auto &[passedGoldenSelection, cutsPassed, assignedPdgCodes] = selection.Execute(pEvent);
//             const auto passedGenericSelection = SelectionHelper::IsCutPassed(cutsPassed, config.global.lastCutGeneric);

//             // Get the reco analysis data (if available, otherwise set to dummy values)
//             const auto recoData = (
//                 passedGenericSelection
//                     ? AnalysisHelper::GetRecoAnalysisData(pEvent->reco, assignedPdgCodes, passedGoldenSelection)
//                     : AnalysisHelper::GetDummyAnalysisData()
//             );


//             // Here we apply reco-level phase-space restrictions
//             // For any event that passes the generic selection, get the value of the kinematic quantity and check if it is outside of the
//             // min/max values supplied in the binning. If so, then reject the event.
//             bool passesPhaseSpaceReco = false;
//             if (passedGenericSelection)
//             {
//                 // Start by assuming the event passes the phase-space cuts
//                 passesPhaseSpaceReco = true;

//                 // Check the value of the kinematic quantities are within the phase-space limits
//                 for (const auto &[name, minMax] : phaseSpaceMap)
//                 {
//                     const auto &[min, max] = minMax;
//                     const auto value = getValue.at(name)(recoData);

//                     if (value < min || value > max)
//                     {
//                         passesPhaseSpaceReco = false;
//                         break;
//                     }
//                 }
//             }


//             const auto isSelectedGolden = passedGoldenSelection && passesPhaseSpaceReco;
//             const auto isSelectedGeneric = passedGenericSelection && passesPhaseSpaceReco;
//             std::map<std::string, bool> isSelectedMap = {{"generic",isSelectedGeneric},{"golden",isSelectedGolden}};

//             // Determine if this is truly a CC1Pi event
//             const auto isTrueCC1Pi = (isOverlay || isDetVar || isNuWro) && AnalysisHelper::IsTrueCC1Pi(pEvent, config.global.useAbsPdg);
//             // Get the truth analysis data (if available, otherwise set to dummy values)
//             std::vector<AnalysisHelper::AnalysisData> truthDataVector;
//             for(const auto &protonMomentumThreshold : protonMomentumThresholdVector)
//             {
//                 const auto truthData = (
//                     isTrueCC1Pi
//                         ? AnalysisHelper::GetTruthAnalysisData(pEvent->truth, config.global.useAbsPdg, protonMomentumThreshold)
//                         : AnalysisHelper::GetDummyAnalysisData()
//                 );
//                 truthDataVector.push_back(truthData);
//             }

//             // Here we apply truth-level phase-space restrictions
//             // For all true CC1Pi events, we check if the values of each kinematic variable are within the supplied limits. If not then the
//             // event is not classed as "signal"

//             vector<bool> passesPhaseSpaceTruthVector(protonMomentumThresholdVector.size(), false);
//             if (isTrueCC1Pi)
//             {

//                 for(unsigned int n = 0; n < protonMomentumThresholdVector.size(); ++n)
//                 {

//                     // Check the value of the kinematic quantities are within the phase-space limits
//                     for (const auto &[name, minMax] : phaseSpaceMap)
//                     {

//                         const auto &[min, max] = minMax;
//                         // Start by assuming the event passes the phase-space cuts

//                         passesPhaseSpaceTruthVector.at(n) = true;

//                         const auto value = getValue.at(name)(truthDataVector.at(n));

//                         if (value < min || value > max)
//                         {

//                             passesPhaseSpaceTruthVector.at(n) = false;

//                             break;
//                         }

//                     }

//                 }

//             }


//             auto isCC1PiSignalVector = passesPhaseSpaceTruthVector;
//             for (unsigned int n = 0; n < protonMomentumThresholdVector.size(); ++n)
//             {

//                 isCC1PiSignalVector.at(n) = isCC1PiSignalVector.at(n) && isTrueCC1Pi;
//             }
//             const auto isSignalVector = isCC1PiSignalVector;
//             // const auto isCC1PiSignal = isTrueCC1Pi && passesPhaseSpaceTruth;
//             // const auto isSignal = isCC1PiSignal;



//             // Get the nominal event weight, scaled by the sample normalisation
//             const auto weight = AnalysisHelper::GetNominalEventWeight(pEvent) * normalisation;

//             // -----------------------------------------------------------------------------------------------------------------------------
//             // Handle 'data' (NuWro and BNB)
//             // -----------------------------------------------------------------------------------------------------------------------------

//             if (isDataBNB)
//             {
//                 for(unsigned int n=0; n<xsecMapBNBVector.size(); n++)
//                 {
//                     for (auto &[selectionName, xsecs] : xsecMapBNBVector.at(n))
//                     {
//                         // Determine if we passed the relevant selection
//                         const auto isSelected = isSelectedMap.at(selectionName);
//                         // Only count events passing the selection
//                         if (!isSelected) continue;

//                         for (auto &[name, xsec] : xsecs)
//                             xsec.AddSelectedBNBDataEvent(getValue.at(name)(recoData));
//                     }
//                 }

//                 // For BNB data that's all we need to do!
//                 continue;
//             }

//             if (isNuWro)
//             {
//                 for(unsigned int n=0; n<xsecMapBNBVector.size(); n++)
//                 {
//                     for (auto &[selectionName, xsecs] : xsecMapNuWroVector.at(n))
//                     {
//                         // Determine if we passed the relevant selection
//                         const auto isSelected = isSelectedMap.at(selectionName);
//                         // Only count events passing the selection
//                         if (!isSelected)
//                             continue;

//                         for (auto &[name, xsec] : xsecs)
//                         {
//                             xsec.AddWeightedSelectedBNBDataEvent(getValue.at(name)(recoData), weight);
//                         }
//                     }
//                 }
//                 // For BNB data that's all we need to do!
//                 // continue;// Important: Nuwro files are also used as MC for xsecMapNuWroTrue
//             }

//             // -----------------------------------------------------------------------------------------------------------------------------
//             // Handle the detector variation samples (unisims)
//             // -----------------------------------------------------------------------------------------------------------------------------

//             if (isDetVar)
//             {
//                 // const auto isSignal = isSignalMap.at(selectionName);
//                 // Handle signal events
//                 for(unsigned int n=0; n<xsecMapBNBVector.size(); n++)
//                 {
//                     if (isSignalVector.at(n))
//                     {
//                         for (auto &[selectionName, xsecs] : xsecMapBNBVector.at(n))
//                         {
//                             // Determine if we passed the relevant selection
//                             const auto isSelected = isSelectedMap.at(selectionName);

//                             for (auto &[name, xsec] : xsecs)
//                             {
//                                 const auto recoValue = getValue.at(name)(recoData);
//                                 const auto trueValue = getValue.at(name)(truthDataVector.at(n));
//                                 xsec.AddSignalEventDetVar(recoValue, trueValue, isSelected, weight, sampleName);
//                             }
//                         }
//                     }
//                     // Handle selected background events
//                     else
//                     {
//                         for (auto &[selectionName, xsecs] : xsecMapBNBVector.at(n))
//                         {
//                             // Only use selected background events
//                             const auto isSelected = isSelectedMap.at(selectionName);
//                             if (!isSelected)
//                                 continue;

//                             for (auto &[name, xsec] : xsecs)
//                             {
//                                 const auto recoValue = getValue.at(name)(recoData);
//                                 xsec.AddSelectedBackgroundEventDetVar(recoValue, weight, sampleName);
//                             }
//                         }
//                     }
//                 }

//                 // For detector variation samples, that's all we need to do!
//                 continue;
//             }


//             // -----------------------------------------------------------------------------------------------------------------------------
//             // Handle all other events (i.e those from the nominal simulation): Overlays, dirt, EXT data
//             // -----------------------------------------------------------------------------------------------------------------------------
//             // std::cout<<"###  DP 0";
//             // Get the flux weights
//             auto fluxWeights = (
//                 isOverlay
//                     ? CrossSectionHelper::GetWeightsMap(pEvent->truth, systParams.fluxDimensions, config.extractXSecs.mutuallyExclusiveDimensions)
//                     : CrossSectionHelper::GetUnitWeightsMap(systParams.fluxDimensions)
//             );

//             // Fill the flux-reweightor with all overlay events from a neutrinos of the desired flavour in the active volume
//             if (isOverlay && AnalysisHelper::IsInActiveVolume(pEvent->truth.nuVertex()) &&
//                 std::find(config.flux.nuPdgsSignal.begin(), config.flux.nuPdgsSignal.end(), pEvent->truth.nuPdgCode()) != config.flux.nuPdgsSignal.end())
//             {
//                 scalingData.pFluxReweightor->AddEvent(pEvent->truth.nuEnergy(), weight, fluxWeights);
//             }

//             // Fill the flux-reweightor with all overlay events from a neutrinos of the desired flavour in the active volume
//             if (isNuWro && AnalysisHelper::IsInActiveVolume(pEvent->truth.nuVertex()) &&
//                 std::find(config.flux.nuPdgsSignal.begin(), config.flux.nuPdgsSignal.end(), pEvent->truth.nuPdgCode()) != config.flux.nuPdgsSignal.end())
//             {
//                 scalingDataNuWroTruth.pFluxReweightor->AddEvent(pEvent->truth.nuEnergy(), weight, fluxWeights);
//             }

//             // Get the cross-section weights
//             // ATTN here we optionally scale the cross-section weights down by the genieTuneEventWeight - this is done so we don't
//             // double count this weight (once in the nominal event weight, and once in the xsec systematic event weights)
//             // const auto xsecWeightsScaleFactor = 1.f;
//             // std::cout<<"Before xsecWeightsScaleFactor"<<std::endl;
//             auto xsecWeightsScaleFactor = (isOverlay && config.extractXSecs.scaleXSecWeights) ? pEvent->truth.genieTuneEventWeight() : 1.f;
//             xsecWeightsScaleFactor = std::max(xsecWeightsScaleFactor, 0.0001f);
//             // std::cout<<"After xsecWeightsScaleFactor"<<std::endl;

//             // std::cout<<"xsecWeightsScaleFactor Debugging Point 1 - xsecWeightsScaleFactor: "<<xsecWeightsScaleFactor<<std::endl;

//             auto xsecWeights = (
//                 isOverlay
//                     ? CrossSectionHelper::ScaleWeightsMap(CrossSectionHelper::GetWeightsMap(pEvent->truth, systParams.xsecDimensions, config.extractXSecs.mutuallyExclusiveDimensions), xsecWeightsScaleFactor)
//                     : CrossSectionHelper::GetUnitWeightsMap(systParams.xsecDimensions)
//             );

//             // std::cout<<"xsecWeightsScaleFactor Debugging Point 2"<<std::endl;

//             // Get the reinteraction weights
//             auto reintWeights = (
//                 isOverlay
//                     ? CrossSectionHelper::GetWeightsMap(pEvent->truth, systParams.reintDimensions, config.extractXSecs.mutuallyExclusiveDimensions)
//                     : CrossSectionHelper::GetUnitWeightsMap(systParams.reintDimensions)
//             );

//             // Handle signal events

//             for(unsigned int n=0; n<xsecMapBNBVector.size(); n++)
//             {

//                 if (isSignalVector.at(n))
//                 {

//                     for (auto &[selectionName, xsecs] : xsecMapNuWroVector.at(n))
//                     {

//                         // Determine if we passed the relevant selection
//                         const auto isSelected = isSelectedMap.at(selectionName);


//                         for (auto &[name, xsec] : xsecs)
//                         {

//                             const auto recoValue = getValue.at(name)(recoData);
//                             const auto trueValue = getValue.at(name)(truthDataVector.at(n));

//                             const auto seedString =  selectionName + name + std::to_string(i);


//                             if(isNuWro)
//                             {

//                                 xsecMapNuWroTrueVector.at(n).at(selectionName).at(name).AddSignalEvent(recoValue, trueValue, isSelected, weight, fluxWeights, xsecWeights, reintWeights, seedString);
//                             } else
//                             {

//                                 if(isOverlay) // Genie MC for NuWro (no EXT or dirt events in NuWro)
//                                 {

//                                     xsec.AddSignalEvent(recoValue, trueValue, isSelected, weight, fluxWeights, xsecWeights, reintWeights, seedString);
//                                 }
//                                 if(isOverlay || isDataEXT || isDirt) // Genie MC for Data
//                                 {

//                                     xsecMapBNBVector.at(n).at(selectionName).at(name).AddSignalEvent(recoValue, trueValue, isSelected, weight, fluxWeights, xsecWeights, reintWeights, seedString);
//                                 }
//                             }

//                         }
//                     }

//                 }
//                 // Handle selected background events
//                 else
//                 {

//                     // std::cout<<"xsecWeightsScaleFactor Debugging Point 4"<<std::endl;
//                     for (auto &[selectionName, xsecs] : xsecMapNuWroVector.at(n))
//                     {
//                         // Only use selected background events
//                         const auto isSelected = isSelectedMap.at(selectionName);
//                         if (!isSelected)
//                             continue;

//                         for (auto &[name, xsec] : xsecs)
//                         {
//                             const auto recoValue = getValue.at(name)(recoData);
//                             const auto seedString =  selectionName + name + std::to_string(i);

//                             if(isOverlay || isDataEXT || isDirt)
//                             {
//                                 const auto trueValue = getValue.at(name)(truthDataVector.at(n));
//                                 if(config.global.useNuWro && isOverlay)
//                                 {
//                                     const auto scaledWeightNuWro = weight;
//                                     const auto scaledXSecWeightsNuWro = xsecWeights;
//                                     const auto scaledFluxWeightsNuWro = fluxWeights;
//                                     const auto scaledReintWeightsNuWro = reintWeights;
//                                     const auto sidebandWeightsNuWro = std::vector<float>(systParams.nBootstrapUniverses, 1.0);

//                                     xsec.AddSelectedBackgroundEvent(recoValue, isDirt, weight, fluxWeights, xsecWeights, reintWeights, std::vector<float>(systParams.nBootstrapUniverses, 1.0), seedString);
//                                 }
//                                 if(config.global.useBNB && (isOverlay || isDataEXT || isDirt))
//                                 {
//                                     const auto scaledWeightBNB = weight;
//                                     const auto scaledXSecWeightsBNB = xsecWeights;
//                                     const auto scaledFluxWeightsBNB = fluxWeights;
//                                     const auto scaledReintWeightsBNB = reintWeights;
//                                     const auto sidebandWeightsBNB = std::vector<float>(systParams.nBootstrapUniverses, 1.0);

//                                     xsecMapBNBVector.at(n).at(selectionName).at(name).AddSelectedBackgroundEvent(recoValue, isDirt, scaledWeightBNB, scaledFluxWeightsBNB, scaledXSecWeightsBNB, scaledReintWeightsBNB, sidebandWeightsBNB, seedString);
//                                 }

//                             } 
//                             else if(isNuWro)
//                             {
//                                 // The multisim weights here are meaningless
//                                 xsecMapNuWroTrueVector.at(n).at(selectionName).at(name).AddSelectedBackgroundEvent(recoValue, false, weight, fluxWeights, xsecWeights, reintWeights, std::vector<float>(systParams.nBootstrapUniverses, 1.0), seedString);
//                             }
//                             else
//                             {
//                                 std::cout<<"Datatype of signal event is unaccounted for."<<std::endl;
//                                 throw std::logic_error("Datatype of signal event is unaccounted for.");
//                             }
//                         }
//                     }

//                 }

//             }


//         }
//     }

//     std::vector<float> traceVector; 
//     // Loop over all cross-section objects
//     for(unsigned int n=0; n<xsecMapBNBVector.size(); n++)
//     {
//         for(const auto&[postfix, map]: {std::make_pair("NuWro", xsecMapNuWroVector.at(n))})//, std::make_pair("BNB", xsecMapBNBVector.at(n))})
//         {

//             for (const auto &[selectionName, xsecs] : map)
//             {

//                 for (const auto &[name, xsec] : xsecs)
//                 {
//                     std::cout << "Processing cross-section: "<<selectionName<< " - " << name << std::endl;

//                     // -----------------------------------------------------------------------------------------------------------------------------
//                     // Get the event rates for BNB data, backgrounds, and signal
//                     // -----------------------------------------------------------------------------------------------------------------------------
//                     const auto selectedEventsData = xsec.GetSelectedBNBDataEvents();
//                     std::cout << "Selected BNB data events" << std::endl;
//                     FormattingHelper::SaveMatrix(selectedEventsData, "xsec_" + selectionName + "_" + name + "_data_selected_eventRate_"+postfix+"_"+std::to_string(n)+".txt");

//                     const auto selectedEventsBackground = xsec.GetSelectedBackgroundEvents();
//                     std::cout << "Selected background events" << std::endl;
//                     FormattingHelper::SaveMatrix(selectedEventsBackground, "xsec_" + selectionName + "_" + name + "_background_selected_eventRate_"+postfix+"_"+std::to_string(n)+".txt");

//                     const auto selectedEventsSignal = xsec.GetSelectedSignalEvents();
//                     std::cout << "Selected signal events" << std::endl;
//                     FormattingHelper::SaveMatrix(selectedEventsSignal, "xsec_" + selectionName + "_" + name + "_signal_selected_eventRate_"+postfix+"_"+std::to_string(n)+".txt");

//                     const auto allEventsSignal = xsec.GetSignalEvents();
//                     std::cout << "All signal events" << std::endl;
//                     FormattingHelper::SaveMatrix(allEventsSignal, "xsec_" + selectionName + "_" + name + "_signal_all_eventRate_"+postfix+"_"+std::to_string(n)+".txt");
//                     std::cout << "After all signal events" << std::endl;
//                     // -----------------------------------------------------------------------------------------------------------------------------
//                     // Get the cross-section as measured with BNB data along with it's uncertainties
//                     // -----------------------------------------------------------------------------------------------------------------------------
//                     const auto data = xsec.GetBNBDataCrossSection(scalingData);
//                     std::cout << "BNB data cross-section (reco-space)" << std::endl;
//                     FormattingHelper::SaveMatrix(data, "xsec_" + selectionName + "_" + name + "_data_"+postfix+"_"+std::to_string(n)+".txt");

//                     const auto dataStatUncertainties = xsec.GetBNBDataCrossSectionStatUncertainty(scalingData);
//                     std::cout << "BNB data stat uncertainty" << std::endl;
//                     FormattingHelper::SaveMatrix(dataStatUncertainties, "xsec_" + selectionName + "_" + name + "_data_stat_"+postfix+"_"+std::to_string(n)+".txt");

//                     const auto dataSystBiasCovariances = xsec.GetBNBDataCrossSectionSystUncertainties(scalingData);
//                     for (const auto &[group, map] : dataSystBiasCovariances)
//                     {
//                         for (const auto &[paramName, biasCovariance] : map)
//                         {
//                             const auto &[pBias, pCovariance] = biasCovariance;

//                             std::cout << "BNB data syst uncertainty: " << group << " " << paramName << std::endl;
//                             std::cout << "Bias vector" << std::endl;
//                             FormattingHelper::SaveMatrix(*pBias, "xsec_" + selectionName + "_" + name + "_data_" + group + "_" + paramName + "_bias_"+postfix+"_"+std::to_string(n)+".txt");
//                             std::cout << "Covariance matrix" << std::endl;
//                             FormattingHelper::SaveMatrix(*pCovariance, "xsec_" + selectionName + "_" + name + "_data_" + group + "_" + paramName + "_covariance_"+postfix+"_"+std::to_string(n)+".txt");
//                         }
//                     }

//                     // -----------------------------------------------------------------------------------------------------------------------------
//                     // Get the predicted cross-section along with its MC stat uncertainty
//                     // -----------------------------------------------------------------------------------------------------------------------------
//                     const auto prediction = xsec.GetPredictedCrossSection(scalingData);
//                     std::cout << "Predicted cross-section (truth-space)" << std::endl;
//                     FormattingHelper::SaveMatrix(prediction, "xsec_" + selectionName + "_" + name + "_prediction_"+postfix+"_"+std::to_string(n)+".txt");

//                     const auto &[pPredictionStatBias, pPredictionStatCovariance] = xsec.GetPredictedCrossSectionStatUncertainty(scalingData);
//                     std::cout << "Predicted cross-section stat uncertainty" << std::endl;
//                     std::cout << "Bias vector" << std::endl;
//                     FormattingHelper::SaveMatrix(*pPredictionStatBias, "xsec_" + selectionName + "_" + name + "_prediction_stat_bias_"+postfix+"_"+std::to_string(n)+".txt");
//                     std::cout << "Covariance matrix" << std::endl;
//                     FormattingHelper::SaveMatrix(*pPredictionStatCovariance, "xsec_" + selectionName + "_" + name + "_prediction_stat_covariance_"+postfix+"_"+std::to_string(n)+".txt");


//                     // const auto &[pPredictionSidebandStatBias, pPredictionSidebandStatCovariance] = xsec.GetPredictedSidebandCrossSectionStatUncertainty(scalingData);
//                     // std::cout << "Predicted cross-section stat uncertainty" << std::endl;
//                     // std::cout << "Bias vector" << std::endl;
//                     // FormattingHelper::SaveMatrix(*pPredictionSidebandStatBias, "xsec_" + selectionName + "_" + name + "_prediction_sideband_stat_bias_"+postfix+".txt");
//                     // std::cout << "Covariance matrix" << std::endl;
//                     // FormattingHelper::SaveMatrix(*pPredictionSidebandStatCovariance, "xsec_" + selectionName + "_" + name + "_prediction_sideband_stat_covariance_"+postfix+".txt");

//                     // -----------------------------------------------------------------------------------------------------------------------------
//                     // Get the smearing matrix along with its uncertainties
//                     // -----------------------------------------------------------------------------------------------------------------------------
//                     std::cout << "Smearing Matrix (reco-space rows, truth-space columns)" << std::endl;
//                     const auto smearingMatrix = xsec.GetSmearingMatrix();
//                     FormattingHelper::SaveMatrix(smearingMatrix, "xsec_" + selectionName + "_" + name + "_smearingMatrix_"+postfix+"_"+std::to_string(n)+".txt");
//                     auto trace = 0.f;
//                     for(unsigned int k=0; k<smearingMatrix.GetRows(); k++) trace += smearingMatrix.At(k,k);
//                     traceVector.push_back(trace);

//                     std::cout << "Smearing Matrix AllSelected" << std::endl;
//                     const auto smearingMatrixAllSelected = xsec.GetSmearingMatrixAllSelected();
//                     FormattingHelper::SaveMatrix(smearingMatrix, "xsec_" + selectionName + "_" + name + "_smearingMatrixAllSelected_"+postfix+"_"+std::to_string(n)+".txt");

//                     std::cout << "Smearing Matrix SystBiasCovariances" << std::endl;
//                     const auto smearingMatrixSystBiasCovariances = xsec.GetSmearingMatrixSystUncertainties();
//                     std::cout << "Smearing Matrix SystBiasCovariances - After" << std::endl;
//                     for (const auto &[group, map] : smearingMatrixSystBiasCovariances)
//                     {
//                         for (const auto &[paramName, biasCovariance] : map)
//                         {
//                             const auto &[pBias, pCovariance] = biasCovariance;

//                             std::cout << "Smearing matrix syst uncertainty: " << group << " " << paramName << std::endl;
//                             std::cout << "Bias vector" << std::endl;
//                             FormattingHelper::SaveMatrix(*pBias, "xsec_" + selectionName + "_" + name + "_smearingMatrix_" + group + "_" + paramName + "_bias_"+postfix+"_"+std::to_string(n)+".txt");
//                             std::cout << "Covariance matrix" << std::endl;
//                             FormattingHelper::SaveMatrix(*pCovariance, "xsec_" + selectionName + "_" + name + "_smearingMatrix_" + group + "_" + paramName + "_covariance_"+postfix+"_"+std::to_string(n)+".txt");
//                         }
//                     }
//                 }
//             }
//         }
//     }

//     // Loop over all cross-section objects
//     for(unsigned int n=0; n<xsecMapNuWroTrueVector.size(); n++)
//     {
//         for (const auto &[selectionName, xsecs] : xsecMapNuWroTrueVector.at(n))
//         {
//             for (const auto &[name, xsec] : xsecs)
//             {
//                 std::cout << "True: Processing cross-section: "<<selectionName<< " - " << name << std::endl;

//                 // -----------------------------------------------------------------------------------------------------------------------------
//                 // Get the event rates for BNB data, backgrounds, and signal
//                 // -----------------------------------------------------------------------------------------------------------------------------
//                 // const auto selectedEventsData = xsec.GetSelectedBNBDataEvents();
//                 // std::cout << "True: Selected BNB data events" << std::endl;
//                 // FormattingHelper::SaveMatrix(selectedEventsData, "xsec_" + selectionName + "_" + name + "_data_selected_eventRate_NuWroTruth.txt");

//                 // -----------------------------------------------------------------------------------------------------------------------------
//                 // Get the cross-section as measured with BNB data along with it's uncertainties
//                 // -----------------------------------------------------------------------------------------------------------------------------
//                 // const auto data = xsec.GetNuWroTrueCrossSection(scalingData);
//                 // std::cout << "True: BNB data cross-section (reco-space)" << std::endl;
//                 // FormattingHelper::SaveMatrix(data, "xsec_" + selectionName + "_" + name + "_data_NuWroTruth.txt");

//                 const auto prediction = xsec.GetPredictedCrossSection(scalingDataNuWroTruth);
//                 std::cout << "Predicted NuWro cross-section (truth-space)" << std::endl;
//                 FormattingHelper::SaveMatrix(prediction, "xsec_" + selectionName + "_" + name + "_prediction_NuWroTruth._"+std::to_string(n)+"txt");

//                 std::cout << "NuWro smearing Matrix (reco-space rows, truth-space columns)" << std::endl;
//                 const auto smearingMatrix = xsec.GetSmearingMatrix();
//                 FormattingHelper::SaveMatrix(smearingMatrix, "xsec_" + selectionName + "_" + name + "_smearingMatrix_NuWroTruth_"+std::to_string(n)+".txt");

//                 std::cout << "NuWro smearing Matrix all selected (reco-space rows, truth-space columns)" << std::endl;
//                 const auto smearingMatrixAllSelected = xsec.GetSmearingMatrixAllSelected();
//                 FormattingHelper::SaveMatrix(smearingMatrixAllSelected, "xsec_" + selectionName + "_" + name + "_smearingMatrixAllSelected_NuWroTruth_"+std::to_string(n)+".txt");
//             }
//         }
//     }

//     std::cout<<"------------- Done saving results to file -------------"<<std::endl;



//     std::cout<<"------------- Getting smearing matrix error matrix -------------"<<std::endl;
//     // std::map<unsigned int, ubsmear::UBMatrix> totalErrorMatrixMap;
//     for(unsigned int n=0; n<xsecMapNuWroTrueVector.size(); n++)
//     {
//         for (const auto &[selectionName, xsecs] : xsecMapNuWroTrueVector.at(n))
//         {
//             for (const auto &[name, xsec] : xsecs)
//             {
//                 // Define a lambda function to get a matrix from a file for brevity
//                 const auto getMatrix = [&, selectionName = selectionName, name = name](const std::string &identifier) -> ubsmear::UBMatrix {
//                     return ubsmear::UBFileHelper::ReadMatrix("xsec_" + selectionName + "_" + name + "_" + identifier + "_" + std::to_string(n) + ".txt");
//                 };

//                 const auto getMatrixTrue = [&, selectionName = selectionName, name = name](const std::string &identifier) -> ubsmear::UBMatrix {
//                     return ubsmear::UBFileHelper::ReadMatrix("xsec_" + selectionName + "_" + name + "_" + identifier + "_NuWroTruth_" + std::to_string(n) + ".txt");
//                 };

//                 // Define a lambda function to get a matrix from a file an trim any overflow / underflow bins
//                 const auto &metadata = xsec.GetMetadata();
//                 const auto getTrimmedMatrix = [&] (const std::string &identifier) -> ubsmear::UBMatrix {
//                     return ubsmear::UBSmearingHelper::TrimUnderOverflowBins(getMatrix(identifier), metadata);
//                 };

//                 // Define a lambda function to get a matrix from a file an trim any overflow / underflow bins
//                 // const auto &metadataTrue = metadataMap.at(name);
//                 const auto getTrimmedMatrixTrue = [&] (const std::string &identifier) -> ubsmear::UBMatrix {
//                     return ubsmear::UBSmearingHelper::TrimUnderOverflowBins(getMatrixTrue(identifier), metadata);
//                 };

//                 // Get the data cross-section
//                 const auto data = getTrimmedMatrix("data_NuWro");
//                 // const auto dataTrue = getTrimmedMatrixTrue("data"); //getTrimmedMatrixTrue("data");

//                 // Get the smearing matrix
//                 const auto smearingMatrix = getMatrix("smearingMatrix_NuWro");
//                 const auto smearingMatrixNuwro = getMatrixTrue("smearingMatrix");

//                 // Build the total error matrix for each group of systematic parameters
//                 // Here errorMatrixMap is indexed by [quantity][group], where quantity = "data" or "smearingMatrix"
//                 // Here totalErrorMatrixMap is index by [quantity] and contains the sum (over all groups) of the entries in errorMatrixMap
//                 std::map<std::string, std::map<std::string, ubsmear::UBMatrix> > errorMatrixMap;
//                 // std::map<std::string, ubsmear::UBMatrix> totalErrorMatrixMap;

//                 // -----------------------------------------------------------------------------------------------------------------------------
//                 // Read the matrices from disk
//                 // -----------------------------------------------------------------------------------------------------------------------------
//                 // Loop over the two quantities that have systematic uncertainties
//                 for (const std::string &quantity : {"data", "smearingMatrix"})
//                 {
//                     if(quantity=="data") continue; // Only interested in smearing matrix uncertainties // todo: implement correctly 

//                     // Get the number of bins (for the smearing matrix there are N^2 bins when flattened)
//                     const auto nBins = (quantity == "data" ? data.GetRows() : std::pow(smearingMatrix.GetRows(), 2));

//                     // Define the function to use when getting a matrix (i.e. should we trim the under/overflow bins or not?)
//                     const auto &getMatrixFunction = (quantity == "data"
//                         ? std::function<ubsmear::UBMatrix(const std::string&)>(getTrimmedMatrix)
//                         : std::function<ubsmear::UBMatrix(const std::string&)>(getMatrix));

//                     // Setup an empty error matrix for this quantity
//                     auto errorMatrixTotalSum =  ubsmear::UBMatrixHelper::GetZeroMatrix(nBins, nBins);

//                     // For the data cross-section, we also have a stat uncertainty
//                     // For the smearing matrix, just use a zero vector
//                     const auto statUncertainties = (quantity == "data") ? getMatrixFunction("data_stat_NuWro") : ubsmear::UBMatrixHelper::GetZeroMatrix(nBins, 1);

//                     // To convert these stat uncertainties into an error matrix, we produce a diagonal matrix whose diagonal entries contain the
//                     // variances (i.e. square of the stat uncerainties)
//                     const auto statVariances = ubsmear::ElementWiseOperation(statUncertainties, statUncertainties, [](const auto &l, const auto &r) { return l * r; });
//                     const auto statErrorMatrix = ubsmear::UBMatrixHelper::GetDiagonalMatrix(statVariances);

//                     // Store this in the map
//                     errorMatrixMap[quantity].emplace("stat", statErrorMatrix);

//                     // Add this to the grand total error matrix
//                     errorMatrixTotalSum = errorMatrixTotalSum + statErrorMatrix;

//                     // Handle the multisim parameters
//                     for (const auto &[group, dimensions] : std::map<std::string, CrossSectionHelper::SystDimensionsMap>(
//                         {
//                             {"flux", config.extractXSecs.fluxDimensions},
//                             {"xsec", config.extractXSecs.xsecDimensions},
//                             {"reint", config.extractXSecs.reintDimensions},
//                             {"misc", {
//                                 {"bootstrap", config.extractXSecs.nBootstrapUniverses},
//                                 // {"sidebandWeights", config.extractXSecs.nBootstrapUniverses},
//                                 {"dirt", 2},
//                                 {"POT", 0}
//                             }}
//                         }))
//                     {
//                         // Setup an empty error matrix for this group
//                         auto errorMatrixTotal =  ubsmear::UBMatrixHelper::GetZeroMatrix(nBins, nBins);

//                         // Loop over the parameters in this group
//                         for (const auto &[paramName, nUniverses] : dimensions)
//                         {
//                             // Get the bias vector and covariance matrix
//                             const auto biasVector = getMatrixFunction(quantity + "_" + group + "_" + paramName + "_bias_NuWro");
//                             const auto covarianceMatrix = getMatrixFunction(quantity + "_" + group + "_" + paramName + "_covariance_NuWro");

//                             // Get the total error matrix from the bias and covariance
//                             const auto errorMatrix = CrossSectionHelper::GetErrorMatrix(biasVector, covarianceMatrix);

//                             // Add this error matrix to the total
//                             errorMatrixTotal = errorMatrixTotal + errorMatrix;
//                         }

//                         // Store this total in the map
//                         errorMatrixMap[quantity].emplace(group, errorMatrixTotal);

//                         // Add this total to the grand total error matrix
//                         errorMatrixTotalSum = errorMatrixTotalSum + errorMatrixTotal;
//                     }

//                     // Handle the unisim parameters
//                     for (const auto &[group, dimensions] : std::map<std::string, CrossSectionHelper::SystUnisimDimensionsMap>(
//                         {
//                             {"detector", config.extractXSecs.detVarDimensions}
//                         }))
//                     {
//                         // Setup an empty error matrix for this group
//                         auto errorMatrixTotal = ubsmear::UBMatrixHelper::GetZeroMatrix(nBins, nBins);

//                         // Loop over the parameters in this group
//                         for (const auto &[paramName, cvName] : dimensions)
//                         {
//                             // Get the bias vector and (dummy) covariance matrix

//                             const auto biasVector = getMatrixFunction(quantity + "_" + group + "_" + paramName + "_bias_NuWro");
//                             const auto covarianceMatrix = getMatrixFunction(quantity + "_" + group + "_" + paramName + "_covariance_NuWro");

//                             // Get the total error matrix from the bias and covariance
//                             const auto errorMatrix = CrossSectionHelper::GetErrorMatrix(biasVector, covarianceMatrix);

//                             // Add this error matrix to the total
//                             errorMatrixTotal = errorMatrixTotal + errorMatrix;
//                         }

//                         // Store this total in the map
//                         errorMatrixMap[quantity].emplace(group, errorMatrixTotal);

//                         // Add this total to the grand total error matrix
//                         errorMatrixTotalSum = errorMatrixTotalSum + errorMatrixTotal;
//                     }

//                     // Add the grand total to the map
//                     // totalErrorMatrixMap.emplace(n, errorMatrixTotalSum);

//                     std::cout<<"Trace of smearing matrix for proton momentum threshold "<<protonMomentumThresholdVector.at(n)<<" is "<<traceVector.at(n)<<std::endl;
//                     std::cout << "Smearing matrix error matrix: ";
//                     auto trace = 0.f;
//                     for(unsigned int k=0; k<errorMatrixTotalSum.GetRows(); k++)
//                     {
//                         std::cout<<errorMatrixTotalSum.At(k,k)<<" ";
//                         trace += errorMatrixTotalSum.At(k,k);
//                     }
//                     std::cout<<" | trace: " << trace << std::endl;
//                     FormattingHelper::SaveMatrix(errorMatrixTotalSum, "xsec_" + selectionName + "_" + name + "_smearingMatrix_errorMatrix_NuWro_" + std::to_string(n) + ".txt");
//                 }
//             }
//         }
//     }


//     // std::vector<std::string> protonMomentumString;
//     // for(const auto &v: protonMomentumThresholdVector)
//     // {
//     //     protonMomentumString.push_back(std::to_string(v));
//     // }

//     // // protonMomentumString.insert(protonMomentumString.begin(), "0");

//     // FormattingHelper::Table tableScaled({"Selection", "Cross-section", "nBins", "", "Chi2", "DoF", "Chi2/DoF", "p-value"});
//     // tableScaled.AddEmptyRow();
//     // tableScaled.SetEntry("Selection", selectionName);
//     // tableScaled.SetEntry("Cross-section", xsecName);
//     // tableScaled.SetEntry("nBins", dataScaled.GetRows());
//     // tableScaled.SetEntry("Chi2", chi2Scaled);
//     // tableScaled.SetEntry("DoF", degreesOfFreedomScaled);
//     // tableScaled.SetEntry("Chi2/DoF", chi2Scaled / static_cast<float>(degreesOfFreedomScaled));
//     // tableScaled.SetEntry("p-value", pValueScaled);
//     // tableScaled.WriteToFile("XSecPlots_goodnessOfFitStatistics_NuWroScaled.md");


//     std::cout<<"------------- All Done -------------"<<std::endl;
//     // return;
// }

// } // namespace ubcc1pi_macros
