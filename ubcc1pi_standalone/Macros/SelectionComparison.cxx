/**
 *  @file  ubcc1pi_standalone/Macros/SelectionComparison.cxx
 *
 *  @brief The implementation file of the SelectionComparison macro
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

// #include <fstream> // Todo: not use txt files
// // Boost libraries
// #include "binary_iarchive.hpp"
// // #include "binary_oarchive.hpp"
// #include "binary_object.hpp"
// #include "map.hpp"
// #include "vector.hpp"

using namespace ubcc1pi;

namespace ubcc1pi_macros
{

void SelectionComparison(const Config &config)
{

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
    auto sidebandSelection = SelectionHelper::GetSelection("CC0pi");

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
    }

    // -------------------------------------------------------------------------------------------------------------------------------------
    // Setup the relevent "getters" for each cross-section
    // -------------------------------------------------------------------------------------------------------------------------------------
    // Here we define a map from the name of each cross-section to a function which pulls out the relevant kinematic quanitity from an input
    // analysis data object. Again this is done up-front to reduce code-bloat below.
    std::unordered_map< std::string, std::function<float(const AnalysisHelper::AnalysisData &)> > getValue;

    // The dummy value that will be used as the "kinematic quantity" for the total cross-section
    const auto dummyValue = 0.f;

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
    // Setup the relevent "getters" for the sideband
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


    // -------------------------------------------------------------------------------------------------------------------------------------
    // Setup the input files
    // -------------------------------------------------------------------------------------------------------------------------------------
    // Here we define a vector of tuples with 4 entries
    //   - First, the sample type (e.g. overlay)
    //   - Second, a string which is used to identify a given detector variation sample (for other sample type, this is unused)
    //   - Third, the path to the input file
    //   - Fourth, the normalisation factor to apply to all events in that file
    std::vector< std::tuple<AnalysisHelper::SampleType, std::string, std::string, float> > inputData;
    for (const auto run: config.global.runs)
    {
        if(run == 1)
        {
            inputData.emplace_back(AnalysisHelper::Overlay, "", config.filesRun1.overlaysFileName, NormalisationHelper::GetOverlaysNormalisationToNuWro(config, 1));
            inputData.emplace_back(AnalysisHelper::Overlay, "", config.filesRun1.nuWroFileName, 1.f);
        }
        else if(run == 2)
        {
            inputData.emplace_back(AnalysisHelper::Overlay, "", config.filesRun2.overlaysFileName, NormalisationHelper::GetOverlaysNormalisation(config, 2));
        }
        else if(run == 3)
        {
            inputData.emplace_back(AnalysisHelper::Overlay, "", config.filesRun3.overlaysFileName, NormalisationHelper::GetOverlaysNormalisation(config, 3));
        }
        else throw std::logic_error("ExtractSidebandFit - Invalid run number");
    }

    // std::cout<<"### Debug Point -4"<<std::endl;
    // -------------------------------------------------------------------------------------------------------------------------------------
    // Count the events for CC1pi
    // -------------------------------------------------------------------------------------------------------------------------------------
    // Loop over the files
    for (const auto &[sampleType, sampleName, fileName, normalisation] : inputData)
    {
        std::vector<std::vector<float>> signalSelectionMatrix = {{0.f,0.f,0.f,0.f},{0.f,0.f,0.f,0.f},{0.f,0.f,0.f,0.f},{0.f,0.f,0.f,0.f},{0.f,0.f,0.f,0.f}};
        // std::vector<std::vector<float>> signalSelectionSidebandMatrix = {{0.f,0.f},{0.f,0.f}};
        // std::vector<std::vector<float>> signalSidebandSelectionMatrix = {{0.f,0.f},{0.f,0.f}};
        // std::vector<std::vector<float>> signalSidebandSelectionSidebandMatrix = {{0.f,0.f},{0.f,0.f}};
        // std::vector<std::vector<float>> signalSidebandSelectionSelectionSidebandMatrix = {{0.f,0.f},{0.f,0.f}};

        std::cout << "Reading input file: " << fileName << std::endl;
        const auto isOverlay = (sampleType == AnalysisHelper::Overlay);
        const auto isDetVar  = (sampleType == AnalysisHelper::DetectorVariation);

        // Open the input file for reading and enable the branches with systematic event weights (if required)
        FileReader reader(fileName);

        if (isOverlay)
            reader.EnableSystematicBranches();

        auto pEvent = reader.GetBoundEventAddress();

        // Loop over the events in the file
        const auto nEvents = reader.GetNumberOfEvents();
        // std::cout<<"!!!!!Only using every 3th event."<<std::endl;
        for (unsigned int i = 0; i < nEvents; i+=1)
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
            
            // std::cout<<"### Debug Point -3"<<std::endl;
            const auto isSelectedGolden = passedGoldenSelection && passesPhaseSpaceReco;
            const auto isSelectedGeneric = passedGenericSelection && passesPhaseSpaceReco;
            std::map<std::string, bool> isSelectedMap = {{"generic",isSelectedGeneric},{"golden",isSelectedGolden}};


            // -----------------------------------------------------------------------------------------------------------------------------
            // Work out if this event passed the selection and apply any additional phase-space cuts based on the input binning
            // -----------------------------------------------------------------------------------------------------------------------------

            // Run the selection
            const auto &[passedSidebandSelection, sidebandCutsPassed, sidebandAssignedPdgCodes] = sidebandSelection.Execute(pEvent);

            // Get the reco analysis data (if available, otherwise set to dummy values)
            const auto recoDataSideband = (
                passedSidebandSelection
                    ? AnalysisHelper::GetRecoAnalysisDataCC0Pi(pEvent->reco, sidebandAssignedPdgCodes)
                    : AnalysisHelper::GetDummyAnalysisData()
            );

            // Here we apply reco-level phase-space restrictions
            // For any event that passes the generic selection, get the value of the kinematic quantity and check if it is outside of the
            // min/max values supplied in the binning. If so, then reject the event.
            bool passesPhaseSpaceRecoSideband = false;
            if (passedSidebandSelection)
            {
                // Start by assuming the event passes the phase-space cuts
                passesPhaseSpaceRecoSideband = true;

                // Check the value of the kinematic quantities are within the phase-space limits
                for (const auto &[name, minMax] : phaseSpaceMap)
                {
                    const auto &[min, max] = minMax;
                    const auto value = getSidebandValue.at(name)(recoDataSideband);

                    if (value < min || value > max)
                    {
                        passesPhaseSpaceRecoSideband = false;
                        // std::cout << "Event failed reco phase-space cuts" << std::endl;
                        break;
                    }
                }
            }

            const auto isSelectedSideband = passedSidebandSelection && passesPhaseSpaceRecoSideband;
            // std::map<std::string, bool> isSelectedMap = {{"sideband",isSelectedSideband}};
            std::map<std::string, bool> isSelectedSidebandMap = {{"generic",isSelectedSideband},{"golden",isSelectedSideband}};

            // -----------------------------------------------------------------------------------------------------------------------------
            // Work out if this event is signal, and apply any phase-space restrictions based on the input binning
            // -----------------------------------------------------------------------------------------------------------------------------

            // Get the nominal event weight, scaled by the sample normalisation
            const auto weight = AnalysisHelper::GetNominalEventWeight(pEvent) * normalisation;

            // Determine if this is truly a CC1Pi event
            const auto isTrueCC1Pi = (isOverlay || isDetVar) && AnalysisHelper::IsTrueCC1Pi(pEvent, config.global.useAbsPdg);
            
            // Determine if this is truly a CC0Pi event
            const auto isTrueCC0Pi = (isOverlay || isDetVar) && AnalysisHelper::IsTrueCC0Pi(pEvent, config.global.useAbsPdg, config.global.protonMomentumThreshold);

            // Get the truth analysis data (if available, otherwise set to dummy values)
            const auto truthData = (
                isTrueCC1Pi
                    ? AnalysisHelper::GetTruthAnalysisData(pEvent->truth, config.global.useAbsPdg, config.global.protonMomentumThreshold)
                    : AnalysisHelper::GetDummyAnalysisData()
            );

            // Get the truth analysis data (if available, otherwise set to dummy values)
            const auto sidebandTruthData = (
                isTrueCC0Pi
                    ? AnalysisHelper::GetTruthAnalysisDataCC0Pi(pEvent->truth, config.global.useAbsPdg, config.global.protonMomentumThreshold)
                    : AnalysisHelper::GetDummyAnalysisData()
            );

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

            bool passesSidebandPhaseSpaceTruth = false;
            if (isTrueCC0Pi)
            {
                // Start by assuming the event passes the phase-space cuts
                passesSidebandPhaseSpaceTruth = true;

                // Check the value of the kinematic quantities are within the phase-space limits
                for (const auto &[name, minMax] : phaseSpaceMap)
                {
                    if(name == "pionMomentum") continue; // Not really compatible with the pion momentum in the CC1pi selection
                    const auto &[min, max] = minMax;
                    const auto value = getSidebandValue.at(name)(sidebandTruthData);

                    if (value < min || value > max)
                    {
                        // std::cout<<"passesSidebandPhaseSpaceTruth failed due to: "<<name<<" - min: "<<min<<" - max: "<<max<<" - value: "<<value<<std::endl;
                        passesSidebandPhaseSpaceTruth = false;
                        break;
                    }
                }
            }

            const auto isCC1PiSignal = isTrueCC1Pi && passesPhaseSpaceTruth;
            const auto isCC0PiSignal = isTrueCC0Pi && passesSidebandPhaseSpaceTruth;
            const auto isSignal = isCC1PiSignal;

            const auto selectionName = "generic";
            const auto isSelected = isSelectedMap.at(selectionName);
            // const auto isSelectedSideband = isSelectedSidebandMap.at(selectionName);

            const std::vector<bool> signalBool              = {true, isCC1PiSignal, isCC0PiSignal, isTrueCC0Pi, isTrueCC1Pi};
            const std::vector<bool> selectionBool           = {true, isSelected, isSelectedSideband, isSelected && isSelectedSideband}; 
            for (int n = 0; n<5; ++n)
            {
                for (int m = 0; m<4; ++m)
                {
                    signalSelectionMatrix[n][m] += weight*(signalBool[n] && selectionBool[m]);
                    // signalSelectionSidebandMatrix[n][m] += weight*(signalBool[n] && selectionBoolSideband[m]);
                    // signalSidebandSelectionMatrix[n][m] += weight*(signalBoolSideband[n] && selectionBool[m]);
                    // signalSidebandSelectionSidebandMatrix[n][m] += weight*(signalBoolSideband[n] && selectionBoolSideband[m]);
                    // signalSidebandSelectionSelectionSidebandMatrix[n][m] += weight*(signalBoolSideband[n] && selectionBool[m] && selectionBoolSideband[m]);
                }
            }
            // std::cout<<weight<<","<<isCC0PiSignal<<","<<isCC1PiSignal<<","<<isSelected<<","<<isSelectedSideband<<","<<isTrueCC0Pi<<","<<isTrueCC1Pi<<","
            // <<passesSidebandPhaseSpaceTruth<<","<<passesPhaseSpaceTruth<<","<<passedSidebandSelection<<","<<passesPhaseSpaceRecoSideband<<","<<passedGenericSelection<<","<<passesPhaseSpaceReco<<std::endl;
        }

        std::cout<<"\n\nsignalSelectionMatrix: "<<std::endl;
        for (int n = 0; n<5; ++n)
        {
            for (int m = 0; m<4; ++m) std::cout<<signalSelectionMatrix[n][m]<<" ";
            std::cout<<std::endl;
        }

        // std::cout<<"signalSelectionSidebandMatrix: "<<std::endl;
        // for (int n = 0; n<2; ++n)
        // {
        //     for (int m = 0; m<2; ++m) std::cout<<signalSelectionSidebandMatrix[n][m]<<" ";
        //     std::cout<<std::endl;
        // }

        // std::cout<<"signalSidebandSelectionMatrix: "<<std::endl;
        // for (int n = 0; n<2; ++n)
        // {
        //     for (int m = 0; m<2; ++m) std::cout<<signalSidebandSelectionMatrix[n][m]<<" ";
        //     std::cout<<std::endl;
        // }

        // std::cout<<"signalSidebandSelectionSidebandMatrix: "<<std::endl;
        // for (int n = 0; n<2; ++n)
        // {
        //     for (int m = 0; m<2; ++m) std::cout<<signalSidebandSelectionSidebandMatrix[n][m]<<" ";
        //     std::cout<<std::endl;
        // }

        // std::cout<<"signalSidebandSelectionSelectionSidebandMatrix: "<<std::endl;
        // for (int n = 0; n<2; ++n)
        // {
        //     for (int m = 0; m<2; ++m) std::cout<<signalSidebandSelectionSelectionSidebandMatrix[n][m]<<" ";
        //     std::cout<<std::endl;
        // }

        // std::cout<<"\n\n";
    }
}

} // namespace ubcc1pi_macros
