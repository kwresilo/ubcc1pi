/**
 *  @file  ubcc1pi_standalone/Macros/Analyzer.cxx
 *
 *  @brief The implementation file of the Analyzer macro
 */

#include "ubcc1pi_standalone/Macros/Macros.h"
#include "ubcc1pi_standalone/Objects/FileReader.h"

#include "ubcc1pi_standalone/Helpers/AnalysisHelper.h"
// #include "ubcc1pi_standalone/Helpers/NormalisationHelper.h"
#include "ubcc1pi_standalone/Helpers/SelectionHelper.h"
#include "ubcc1pi_standalone/Helpers/ExtractionHelper.h"
#include "ubcc1pi_standalone/Objects/TreeWriter.h"

// #include <ctime>
// #include <chrono>
// #include <iomanip>

// // Boost libraries
// // #include "binary_iarchive.hpp"
// #include "binary_oarchive.hpp"
// #include "binary_object.hpp"
// #include "map.hpp"
// #include "vector.hpp"

using namespace ubcc1pi;

namespace ubcc1pi_macros
{

void Analyzer(const Config &config)
{

    // A map from each cross-section to the limits of the phase-space that should be considered. The key is the
    //     // identifier for the kinematic quantity and the mapped value is a pair containing the limits [min, max]
    std::map< std::string, std::pair<float, float> > phaseSpaceMap;
    for (const auto &[name, binning, scaleByBinWidth] : std::vector< std::tuple<std::string, Config::Global::Binning, bool> > {
        { "muonMomentum", config.global.muonMomentum, true },

        { "pionMomentum", config.global.pionMomentum, true },

        { "protonMomentum", config.global.pionMomentum, true }
    })

    {
        // Add to the phase-space map
        phaseSpaceMap.emplace(name, std::pair<float, float>({binning.min, binning.max}));
    }



    // LOOK INTO THESE SELECTIONS, GET DEFAULT DELECTION APPLIES BDT PART
    auto selectionCC1Pi = SelectionHelper::GetCCInclusiveSelection();

    const auto currentTime = std::chrono::system_clock::now();
    const auto currentTimeT = std::chrono::system_clock::to_time_t(currentTime);
    std::stringstream currentDateTime;
    currentDateTime << std::put_time(std::localtime(&currentTimeT), "%Y-%m-%d-%H-%M"); // Convert std::tm to a string

    // -------------------------------------------------------------------------------------------------------------------------------------
    // Setup the output file
    // -------------------------------------------------------------------------------------------------------------------------------------
    // std::vector<std::string> inputFiles;
    // inputFiles.push_back("/uboone/data/users/jdetje/ubcc1piVSpelee/pelee/neutrinoselection_filt_0_4k.root");
    // TreeWriter treeWriter(inputFiles, config.global.outputFile.c_str());

    // Loop over the files
    //
    for (const auto &[sampleType, useThisFile, filePath] : config.inputFiles)
    {
        if(!useThisFile) continue;

        // Find the last occurrence of '/'
        size_t lastSlash = filePath.find_last_of('/');
        // Find the last occurrence of '.'
        size_t lastDot = filePath.find_last_of('.');
        // Extract the substring between the last '/' and the last '.'fileList
        std::string fileName = filePath.substr(lastSlash + 1, lastDot - lastSlash - 1);
        std::cout << "Reading input file: " << fileName << std::endl;

        const auto outputFilePath = config.global.outputPath + fileName + "_ubcc1pi_" + currentDateTime.str() + ".root";
        // Todo: Change the code to combine all input files of the same type across different runs into one output file
        // WRITES NEW TREE, PROBABLY NEED TO MODIIFY
        // CLASSIFIES TYPE OF DATA FOR SELECTION 
        TreeWriter treeWriter(std::vector<std::string>({filePath}), outputFilePath.c_str());
        const auto isOverlay = (sampleType == AnalysisHelper::Overlay);
        const auto isDirt    = (sampleType == AnalysisHelper::Dirt);
        const auto isNuWro   = (sampleType == AnalysisHelper::NuWro);
        const auto isDataBNB = (sampleType == AnalysisHelper::DataBNB);
        const auto isDetVar  = (sampleType == AnalysisHelper::DetectorVariation);
        const auto isDataEXT = (sampleType == AnalysisHelper::DataEXT);
        const auto isMC = isOverlay;

        // if(isDetVar) continue; // Todo: Remove!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        // if(isDataBNB || isDirt || isDataEXT) continue; // Todo: Remove!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        // Open the input file for reading and enable the branches with systematic event weights (if required)
        // TEMPLATED: 2 VERSIONS ONE WITH BRACKETS (NEEDED TO PELEE NTUPLES)
        FileReader<EventPeLEE, SubrunPeLEE> readerPeLEE(filePath, isMC);

        if (isMC) readerPeLEE.EnableSystematicBranches(); // Todo: Is this correct/optimal?

	
        // Loop over the events in the file
        const auto nEvents = readerPeLEE.GetNumberOfEvents();
        const auto pEventPeLEE = readerPeLEE.GetBoundEventAddress();
        std::cout << "### Only processing every 10th event ###" << std::endl;
        for (unsigned int i = 0; i < 1000; i++) //nEvents; i++) // Todo: Remove!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        {

	    //LOAD EVENTS
            AnalysisHelper::PrintLoadingBar(i, nEvents);
            readerPeLEE.LoadEvent(i);
            Event event(*pEventPeLEE, isMC);
            const auto pEvent = std::make_shared<Event>(event);

            // ################################################################
            // #################### Reco CC1pi ####################
            // ################################################################
            const auto &[passedGoldenSelection, cutsPassedCC1Pi, assignedPdgCodesCC1Pi] = selectionCC1Pi.Execute(pEvent);
            const auto passedGenericSelection = SelectionHelper::IsCutPassed(cutsPassedCC1Pi, config.global.lastCutGeneric);

            // Get the reco analysis data (if available, otherwise set to dummy values)
            const auto recoDataCC1Pi = (
                passedGenericSelection
                    ? AnalysisHelper::GetRecoAnalysisData(pEvent->reco, assignedPdgCodesCC1Pi, passedGoldenSelection)
                    : AnalysisHelper::GetDummyAnalysisData()
            );

            // Here we apply reco-level phase-space restrictions
            // For any event that passes the generic selection, get the value of the kinematic quantity and check if it is outside of the
            // min/max values supplied in the binning. If so, then reject the event.
            auto passesPhaseSpaceRecoCC1Pi = false;
            if (passedGenericSelection)
            {
                // Start by assuming the event passes the phase-space cuts
                passesPhaseSpaceRecoCC1Pi = true;

                // Check the value of the kinematic quantities are within the phase-space limits
                for (const auto &[name, minMax] : phaseSpaceMap)
                {
                    const auto &[min, max] = minMax;
                    const auto value = getValueCC1Pi.at(name)(recoDataCC1Pi);

                    if (value < min || value > max)
                    {
                        passesPhaseSpaceRecoCC1Pi = false;
                        break;
                    }
                }
            }

	    
            const auto isSelectedGenericCC1Pi = passedGenericSelection && passesPhaseSpaceRecoCC1Pi;
            const auto isSelectedGoldenCC1Pi = passedGoldenSelection && passesPhaseSpaceRecoCC1Pi;

            // if(isSelectedGenericCC0Pi || isSelectedGenericCC1Pi) {std::cout << "\nDEBUG isSelectedGenericCC0Pi: " << isSelectedGenericCC0Pi << " isSelectedGenericCC1Pi: " << isSelectedGenericCC1Pi << "\n" << std::endl;}
            // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            // Reco variables
            // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
	    // ADD BRANCHES TO EMPTY NTUPLES?
            treeWriter.SetOutputBranchAddress("cc1pi_selected_generic", (void*)&isSelectedGenericCC1Pi, "cc1pi_selected_generic/O" );
            treeWriter.SetOutputBranchAddress("cc1pi_selected_golden", (void*)&isSelectedGoldenCC1Pi, "cc1pi_selected_golden/O" );

            treeWriter.SetOutputBranchAddress("cc1pi_reco_muonMomentum", (void*)&recoDataCC1Pi.muonMomentum, "cc1pi_reco_muonMomentum/F" );
            treeWriter.SetOutputBranchAddress("cc1pi_reco_muonCosTheta", (void*)&recoDataCC1Pi.muonCosTheta, "cc1pi_reco_muonCosTheta/F" );
            treeWriter.SetOutputBranchAddress("cc1pi_reco_muonPhi", (void*)&recoDataCC1Pi.muonPhi, "cc1pi_reco_muonPhi/F" );
            treeWriter.SetOutputBranchAddress("cc1pi_reco_pionMomentum", (void*)&recoDataCC1Pi.pionMomentum, "cc1pi_reco_pionMomentum/F" );
            treeWriter.SetOutputBranchAddress("cc1pi_reco_pionCosTheta", (void*)&recoDataCC1Pi.pionCosTheta, "cc1pi_reco_pionCosTheta/F" );
            treeWriter.SetOutputBranchAddress("cc1pi_reco_pionPhi", (void*)&recoDataCC1Pi.pionPhi, "cc1pi_reco_pionPhi/F" );
            treeWriter.SetOutputBranchAddress("cc1pi_reco_muonPionAngle", (void*)&recoDataCC1Pi.muonPionAngle, "cc1pi_reco_muonPionAngle/F" );
            treeWriter.SetOutputBranchAddress("cc1pi_reco_nProtons", (void*)&recoDataCC1Pi.nProtons, "cc1pi_reco_nProtons/I" );

            if (!isDataBNB) // For BNB data that's all we need to do!
            {

                // #####################################################
                // #################### Truth CC1pi ####################
                // #####################################################
                const auto isTrueCC1Pi = (isOverlay || isDetVar || isNuWro) && AnalysisHelper::IsTrueCC1Pi(pEvent, config.global.useAbsPdg);

                // Get the truth analysis data (if available, otherwise set to dummy values)
                const auto truthDataCC1Pi = (
                    isTrueCC1Pi
                        ? AnalysisHelper::GetTruthAnalysisData(pEvent->truth, config.global.useAbsPdg, config.global.protonMomentumThreshold)
                        : AnalysisHelper::GetDummyAnalysisData()
                );

                // Here we apply truth-level phase-space restrictions
                // For all true CC1Pi events, we check if the values of each kinematic variable are within the supplied limits. If not then the
                // event is not classed as "signal"
                bool passesPhaseSpaceTruthCC1Pi = false;
                if (isTrueCC1Pi)
                {
                    // Start by assuming the event passes the phase-space cuts
                    passesPhaseSpaceTruthCC1Pi = true;

                    // Check the value of the kinematic quantities are within the phase-space limits
                    for (const auto &[name, minMax] : phaseSpaceMap)
                    {
                        const auto &[min, max] = minMax;
                        const auto value = getValueCC1Pi.at(name)(truthDataCC1Pi);

                        if (value < min || value > max)
                        {
                            passesPhaseSpaceTruthCC1Pi = false;
                            break;
                        }
                    }
                }
                const auto isCC1PiSignal = isTrueCC1Pi && passesPhaseSpaceTruthCC1Pi;

                // if(isCC0PiSignal || isCC1PiSignal) {std::cout << "\nDEBUG isCC0PiSignal: " << isCC0PiSignal << " isCC1PiSignal: " << isCC1PiSignal << "\n" << std::endl;}

                // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                // Truth variables
                // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                treeWriter.SetOutputBranchAddress("cc0pi_signal", (void*)&isCC0PiSignal, "cc0pi_signal/O" );
                treeWriter.SetOutputBranchAddress("cc1pi_signal", (void*)&isCC1PiSignal, "cc1pi_signal/O" );

                treeWriter.SetOutputBranchAddress("cc1pi_truth_muonMomentum", (void*)&truthDataCC1Pi.muonMomentum, "cc1pi_truth_muonMomentum/F" );
                treeWriter.SetOutputBranchAddress("cc1pi_truth_muonCosTheta", (void*)&truthDataCC1Pi.muonCosTheta, "cc1pi_truth_muonCosTheta/F" );
                treeWriter.SetOutputBranchAddress("cc1pi_truth_muonPhi", (void*)&truthDataCC1Pi.muonPhi, "cc1pi_truth_muonPhi/F" );
                treeWriter.SetOutputBranchAddress("cc1pi_truth_pionMomentum", (void*)&truthDataCC1Pi.pionMomentum, "cc1pi_truth_pionMomentum/F" );
                treeWriter.SetOutputBranchAddress("cc1pi_truth_pionCosTheta", (void*)&truthDataCC1Pi.pionCosTheta, "cc1pi_truth_pionCosTheta/F" );
                treeWriter.SetOutputBranchAddress("cc1pi_truth_pionPhi", (void*)&truthDataCC1Pi.pionPhi, "cc1pi_truth_pionPhi/F" );
                treeWriter.SetOutputBranchAddress("cc1pi_truth_muonPionAngle", (void*)&truthDataCC1Pi.muonPionAngle, "cc1pi_truth_muonPionAngle/F" );
                treeWriter.SetOutputBranchAddress("cc1pi_truth_nProtons", (void*)&truthDataCC1Pi.nProtons, "cc1pi_truth_nProtons/I" );

                treeWriter.SetOutputBranchAddress("cc0pi_truth_muonMomentum", (void*)&truthDataCC0Pi.muonMomentum, "cc0pi_truth_muonMomentum/F" );
                treeWriter.SetOutputBranchAddress("cc0pi_truth_muonCosTheta", (void*)&truthDataCC0Pi.muonCosTheta, "cc0pi_truth_muonCosTheta/F" );
                treeWriter.SetOutputBranchAddress("cc0pi_truth_muonPhi", (void*)&truthDataCC0Pi.muonPhi, "cc0pi_truth_muonPhi/F" );
                treeWriter.SetOutputBranchAddress("cc0pi_truth_pionMomentum", (void*)&truthDataCC0Pi.pionMomentum, "cc0pi_truth_pionMomentum/F" );
                treeWriter.SetOutputBranchAddress("cc0pi_truth_pionCosTheta", (void*)&truthDataCC0Pi.pionCosTheta, "cc0pi_truth_pionCosTheta/F" );
                treeWriter.SetOutputBranchAddress("cc0pi_truth_pionPhi", (void*)&truthDataCC0Pi.pionPhi, "cc0pi_truth_pionPhi/F" );
                treeWriter.SetOutputBranchAddress("cc0pi_truth_muonPionAngle", (void*)&truthDataCC0Pi.muonPionAngle, "cc0pi_truth_muonPionAngle/F" );
                treeWriter.SetOutputBranchAddress("cc0pi_truth_nProtons", (void*)&truthDataCC0Pi.nProtons, "cc0pi_truth_nProtons/I" );

                // Other reco variables
                // ################# Signal definition flags #################
                // MC event category
                // Event weights

                // Backtracked neutrino purity and completeness
                // Number of neutrino slices identified by the SliceID
                // CC1pi cuts
                // CC0pi cuts

                // ################# MC event properties #################
                treeWriter.SetOutputBranchAddress("mc_nu_pdg", (void*)pEventPeLEE->truth.nu_pdg.GetAddress(), "mc_nu_pdg/I");
                // ################# Reco event properties #################
            }
            else
            {
                // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                // Dummy truth variables
                // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                const auto dummyBool = false;
                const auto dummyFloat = -std::numeric_limits<float>::max();
                const auto dummyInt = -std::numeric_limits<int>::max();

                treeWriter.SetOutputBranchAddress("cc0pi_signal", (void*)&dummyBool, "cc0pi_signal/O" );
                treeWriter.SetOutputBranchAddress("cc1pi_signal", (void*)&dummyBool, "cc1pi_signal/O" );

                treeWriter.SetOutputBranchAddress("cc1pi_truth_muonMomentum", (void*)&dummyFloat, "cc1pi_truth_muonMomentum/F" );
                treeWriter.SetOutputBranchAddress("cc1pi_truth_muonCosTheta", (void*)&dummyFloat, "cc1pi_truth_muonCosTheta/F" );
                treeWriter.SetOutputBranchAddress("cc1pi_truth_muonPhi", (void*)&dummyFloat, "cc1pi_truth_muonPhi/F" );
                treeWriter.SetOutputBranchAddress("cc1pi_truth_pionMomentum", (void*)&dummyFloat, "cc1pi_truth_pionMomentum/F" );
                treeWriter.SetOutputBranchAddress("cc1pi_truth_pionCosTheta", (void*)&dummyFloat, "cc1pi_truth_pionCosTheta/F" );
                treeWriter.SetOutputBranchAddress("cc1pi_truth_pionPhi", (void*)&dummyFloat, "cc1pi_truth_pionPhi/F" );
                treeWriter.SetOutputBranchAddress("cc1pi_truth_muonPionAngle", (void*)&dummyFloat, "cc1pi_truth_muonPionAngle/F" );
                treeWriter.SetOutputBranchAddress("cc1pi_truth_nProtons", (void*)&dummyInt, "cc1pi_truth_nProtons/I" );

                treeWriter.SetOutputBranchAddress("cc0pi_truth_muonMomentum", (void*)&dummyFloat, "cc0pi_truth_muonMomentum/F" );
                treeWriter.SetOutputBranchAddress("cc0pi_truth_muonCosTheta", (void*)&dummyFloat, "cc0pi_truth_muonCosTheta/F" );
                treeWriter.SetOutputBranchAddress("cc0pi_truth_muonPhi", (void*)&dummyFloat, "cc0pi_truth_muonPhi/F" );
                treeWriter.SetOutputBranchAddress("cc0pi_truth_pionMomentum", (void*)&dummyFloat, "cc0pi_truth_pionMomentum/F" );
                treeWriter.SetOutputBranchAddress("cc0pi_truth_pionCosTheta", (void*)&dummyFloat, "cc0pi_truth_pionCosTheta/F" );
                treeWriter.SetOutputBranchAddress("cc0pi_truth_pionPhi", (void*)&dummyFloat, "cc0pi_truth_pionPhi/F" );
                treeWriter.SetOutputBranchAddress("cc0pi_truth_muonPionAngle", (void*)&dummyFloat, "cc0pi_truth_muonPionAngle/F" );
                treeWriter.SetOutputBranchAddress("cc0pi_truth_nProtons", (void*)&dummyInt, "cc0pi_truth_nProtons/I" );

                treeWriter.SetOutputBranchAddress("mc_nu_pdg", (void*)&dummyInt, "mc_nu_pdg/I");
            }

            // // std::cout<<"DEBUG Point 6"<<std::endl;

            // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            // File-level variables
            // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            treeWriter.SetOutputBranchAddress("is_mc", (void*)&isMC, "is_mc/O");
            treeWriter.SetOutputBranchAddress("is_Overlay", (void*)&isOverlay, "is_Overlay/O");
            treeWriter.SetOutputBranchAddress("is_Dirt", (void*)&isDirt, "is_Dirt/O");
            treeWriter.SetOutputBranchAddress("is_NuWro", (void*)&isNuWro, "is_NuWro/O");
            treeWriter.SetOutputBranchAddress("is_DataBNB", (void*)&isDataBNB, "is_DataBNB/O");
            treeWriter.SetOutputBranchAddress("is_DetVar", (void*)&isDetVar, "is_DetVar/O");
            treeWriter.SetOutputBranchAddress("is_DataEXT", (void*)&isDataEXT, "is_DataEXT/O");
            // // std::cout<<"DEBUG Point 7"<<std::endl;

            treeWriter.CreateNoNewBranches(); // Only creates branches for the first run through the loop
            // // std::cout<<"DEBUG Point 8"<<std::endl;
            treeWriter.Fill();
            // // std::cout<<"DEBUG Point 9"<<std::endl;
        } // End of event-level iteration

    } // End of file-level iterration
    std::cout<<"-----------------Done-----------------"<<std::endl;


}

} // namespace ubcc1pi_macros
