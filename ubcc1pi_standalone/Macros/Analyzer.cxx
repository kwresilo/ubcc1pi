/**
 *  @file  ubcc1pi_standalone/Macros/Analyzer.cxx
 *
 *  @brief The implementation file of the Analyzer macro
 */

#include "ubcc1pi_standalone/Macros/Macros.h"
#include "ubcc1pi_standalone/Objects/FileReader.h"

#include "ubcc1pi_standalone/Helpers/AnalysisHelper.h"
#include "ubcc1pi_standalone/Helpers/NormalisationHelper.h"
#include "ubcc1pi_standalone/Helpers/SelectionHelper.h"
// #include "ubcc1pi_standalone/Helpers/CrossSectionHelper.h"
// #include "ubcc1pi_standalone/Helpers/FormattingHelper.h"
// #include "ubcc1pi_standalone/Helpers/FittingHelper.h"
#include "ubcc1pi_standalone/Helpers/ExtractionHelper.h"
// #include "ubcc1pi_standalone/ubsmear/inc/ubsmear/Helpers/UBSmearingHelper.h"
#include "ubcc1pi_standalone/Objects/TreeWriter.h"

// Boost libraries
// #include "binary_iarchive.hpp"
#include "binary_oarchive.hpp"
#include "binary_object.hpp"
#include "map.hpp"
#include "vector.hpp"

using namespace ubcc1pi;

namespace ubcc1pi_macros
{

void Analyzer(const Config &config)
{

    // We additionally make a map from each cross-section to the limits of the phase-space that we should consider. The key is the
    // identifier for the kinematic quantity and the mapped value is a pair containing the limits [min, max]
    std::map< std::string, std::pair<float, float> > phaseSpaceMap;
    for (const auto &[name, binning, scaleByBinWidth] : std::vector< std::tuple<std::string, Config::Global::Binning, bool> > {

        // The names of the cross-section kinematic parameters, and their binning information.
        // The third (boolean) parameter indicates if the cross-section bins should be scaled by their width
        { "muonCosTheta", config.global.muonCosTheta, true }, 
        { "muonPhi", config.global.muonPhi, true }, 
        { "muonMomentum", config.global.muonMomentum, true }, 

        { "pionCosTheta", config.global.pionCosTheta, true }, 
        { "pionPhi", config.global.pionPhi, true  }, 
        { "pionMomentum", config.global.pionMomentum, true }, 

        { "muonPionAngle", config.global.muonPionAngle, true }, 
        { "nProtons", config.global.nProtons, false }

    })
    {
        // Add to the phase-space map
        phaseSpaceMap.emplace(name, std::pair<float, float>({binning.min, binning.max}));
    }

    // -------------------------------------------------------------------------------------------------------------------------------------
    // Setup the relevent "getters" for each cross-section and for the CC0pi selection
    // -------------------------------------------------------------------------------------------------------------------------------------
    ExtractionHelper::AnalysisValueMap getValueCC1Pi;
    ExtractionHelper::AnalysisValueMap getValueCC0Pi;
    ExtractionHelper::PopulateAnalysisValueMap(getValueCC1Pi, false);
    ExtractionHelper::PopulateAnalysisValueMap(getValueCC0Pi, true);

    std::cout<<"..........................................\nUSING Modified CC0pi Selection: muonLikeProtonValue=-0.48f, barelyResemblingProtonValue=0.12f\n.........................................."<<std::endl;
    auto selectionCC0Pi = SelectionHelper::GetCC0piSelectionModified(-0.48f, 0.12f);
    auto selectionCC1Pi = SelectionHelper::GetDefaultSelection();

    ExtractionHelper::InputFileList inputData;
    typedef std::vector< std::tuple<AnalysisHelper::SampleType, std::string, std::string, float> > InputFileList;
    // inputData.emplace_back(AnalysisHelper::Overlay, "", "/uboone/data/users/jdetje/ubcc1piVSpelee/pelee/neutrinoselection_filt_0_4k.root", 1);
    // inputData.emplace_back(AnalysisHelper::Overlay, "", "/uboone/app/users/jdetje/searchingfornues/files/steps4/neutrinoselection_filt_upodated5.root", 1);
    inputData.emplace_back(AnalysisHelper::DataBNB, "", "/uboone/data/users/jdetje/pelee_v08_00_00_70/bnb_beam_on_peleeTuple_uboone_v08_00_00_70_run2_E1_head5.root", 1);


    // -------------------------------------------------------------------------------------------------------------------------------------
    // Setup the output file
    // -------------------------------------------------------------------------------------------------------------------------------------
    // std::vector<std::string> inputFiles;
    // inputFiles.push_back("/uboone/data/users/jdetje/ubcc1piVSpelee/pelee/neutrinoselection_filt_0_4k.root");
    // TreeWriter treeWriter(inputFiles, config.global.outputFile.c_str());

    // Loop over the files
    for (const auto &[sampleType, sampleName, fileName, normalisation] : inputData)
    {
        std::cout << "Reading input file: " << fileName << std::endl;
        // Todo: Change the code to combine all input files of the same type across different runs into one output file
        TreeWriter treeWriter(std::vector<std::string>({fileName}), config.global.outputFile.c_str());
        const auto isOverlay = (sampleType == AnalysisHelper::Overlay);
        const auto isDirt    = (sampleType == AnalysisHelper::Dirt);
        const auto isNuWro   = (sampleType == AnalysisHelper::NuWro);
        const auto isDataBNB = (sampleType == AnalysisHelper::DataBNB);
        const auto isDetVar  = (sampleType == AnalysisHelper::DetectorVariation);
        const auto isDataEXT = (sampleType == AnalysisHelper::DataEXT);
        const auto isMC = (!isDataBNB || !isDataEXT);

        // if(isDetVar) continue; // Todo: Remove!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        // if(isDataBNB || isDirt || isDataEXT) continue; // Todo: Remove!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        // Open the input file for reading and enable the branches with systematic event weights (if required)
        FileReader<EventPeLEE, SubrunPeLEE> readerPeLEE(fileName);

        if (isMC) readerPeLEE.EnableSystematicBranches(); // Todo: Is this correct/optimal?

        // Loop over the events in the file
        const auto nEvents = readerPeLEE.GetNumberOfEvents();
        const auto pEventPeLEE = readerPeLEE.GetBoundEventAddress();
        for (unsigned int i = 0; i < nEvents; ++i)
        {
            AnalysisHelper::PrintLoadingBar(i, nEvents);
            readerPeLEE.LoadEvent(i);
            Event event(*pEventPeLEE, isMC);
            const auto pEvent = std::make_shared<Event>(event);

            // ####################################################################
            // #################### Reco CC0pi ####################
            // ####################################################################
            const auto &[passedGoldenselectionCC0Pi, cutsPassedCC0Pi, assignedPdgCodesCC0Pi] = selectionCC0Pi.Execute(pEvent);
            const auto passedGenericselectionCC0Pi = SelectionHelper::IsCutPassed(cutsPassedCC0Pi, config.global.lastCutGeneric);

            // Get the reco analysis data (if available, otherwise set to dummy values)
            const auto recoDataCC0Pi = (
                passedGenericselectionCC0Pi
                    ? AnalysisHelper::GetRecoAnalysisDataCC0Pi(pEvent->reco, assignedPdgCodesCC0Pi, passedGoldenselectionCC0Pi)
                    : AnalysisHelper::GetDummyAnalysisData()
            );

            // Here we apply reco-level phase-space restrictions
            // For any event that passes the generic selection, get the value of the kinematic quantity and check if it is outside of the
            // min/max values supplied in the binning. If so, then reject the event.
            auto passesPhaseSpaceRecoCC0Pi = false;
            if (passedGenericselectionCC0Pi)
            {
                // Start by assuming the event passes the phase-space cuts
                passesPhaseSpaceRecoCC0Pi = true;

                // Check the value of the kinematic quantities are within the phase-space limits
                for (const auto &[name, minMax] : phaseSpaceMap)
                {
                    const auto &[min, max] = minMax;
                    const auto value = getValueCC0Pi.at(name)(recoDataCC0Pi);

                    if (value < min || value > max)
                    {
                        passesPhaseSpaceRecoCC0Pi = false;
                        break;
                    }
                }
            }
            const auto isSelectedGenericCC0Pi = passedGenericselectionCC0Pi && passesPhaseSpaceRecoCC0Pi;
            const auto isSelectedGoldenCC0Pi = passedGoldenselectionCC0Pi && passesPhaseSpaceRecoCC0Pi;

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

            // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            // Reco variables
            // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            treeWriter.SetOutputBranchAddress("cc0pi_selected_generic", (void*)&isSelectedGenericCC0Pi, "cc0pi_selected_generic/O" );
            treeWriter.SetOutputBranchAddress("cc0pi_selected_golden", (void*)&isSelectedGoldenCC0Pi, "cc0pi_selected_golden/O" );
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

            treeWriter.SetOutputBranchAddress("cc0pi_reco_muonMomentum", (void*)&recoDataCC0Pi.muonMomentum, "cc0pi_reco_muonMomentum/F" );
            treeWriter.SetOutputBranchAddress("cc0pi_reco_muonCosTheta", (void*)&recoDataCC0Pi.muonCosTheta, "cc0pi_reco_muonCosTheta/F" );
            treeWriter.SetOutputBranchAddress("cc0pi_reco_muonPhi", (void*)&recoDataCC0Pi.muonPhi, "cc0pi_reco_muonPhi/F" );
            treeWriter.SetOutputBranchAddress("cc0pi_reco_pionMomentum", (void*)&recoDataCC0Pi.pionMomentum, "cc0pi_reco_pionMomentum/F" );
            treeWriter.SetOutputBranchAddress("cc0pi_reco_pionCosTheta", (void*)&recoDataCC0Pi.pionCosTheta, "cc0pi_reco_pionCosTheta/F" );
            treeWriter.SetOutputBranchAddress("cc0pi_reco_pionPhi", (void*)&recoDataCC0Pi.pionPhi, "cc0pi_reco_pionPhi/F" );
            treeWriter.SetOutputBranchAddress("cc0pi_reco_muonPionAngle", (void*)&recoDataCC0Pi.muonPionAngle, "cc0pi_reco_muonPionAngle/F" );
            treeWriter.SetOutputBranchAddress("cc0pi_reco_nProtons", (void*)&recoDataCC0Pi.nProtons, "cc0pi_reco_nProtons/I" );

            // treeWriter.SetOutputBranchAddress("muon_candidate_idx", &pEventPeLEE->muon_candidate_idx_, "muon_candidate_idx/I" );
            // treeWriter.SetOutputBranchAddress("pion_candidate_idx", &pEventPeLEE->pion_candidate_idx_, "pion_candidate_idx/I" );
            // treeWriter.SetOutputBranchAddress("lead_p_candidate_idx", &pEventPeLEE->lead_p_candidate_idx_, "lead_p_candidate_idx/I" );

            if (!isDataBNB) // For BNB data that's all we need to do!
            {
                // #####################################################
                // #################### Truth CC0pi ####################
                // #####################################################
                // Determine if this is truly a CC0Pi event
                const auto isTrueCC0Pi = (isOverlay || isDetVar || isNuWro) && AnalysisHelper::IsTrueCC0Pi(pEvent, config.global.useAbsPdg, config.global.protonMomentumThreshold);

                // Get the truth analysis data (if available, otherwise set to dummy values)
                const auto truthDataCC0Pi = (
                    (isTrueCC0Pi)
                        ? AnalysisHelper::GetTruthAnalysisDataCC0Pi(pEvent->truth, config.global.useAbsPdg, config.global.protonMomentumThreshold)
                        : AnalysisHelper::GetDummyAnalysisData()
                );

                // Here we apply truth-level phase-space restrictions
                // For all true CC0Pi events, we check if the values of each kinematic variable are within the supplied limits. If not then the
                // event is not classed as "signal"
                bool passesPhaseSpaceTruthCC0Pi = false;
                if (isTrueCC0Pi)
                {
                    // Start by assuming the event passes the phase-space cuts
                    passesPhaseSpaceTruthCC0Pi = true;

                    // Check the value of the kinematic quantities are within the phase-space limits
                    for (const auto &[name, minMax] : phaseSpaceMap)
                    {
                        // if(name == "pionMomentum") continue; // Not really compatible with the pion momentum in the CC1pi selection // todo check this
                        const auto &[min, max] = minMax;
                        const auto value = getValueCC0Pi.at(name)(truthDataCC0Pi);

                        if (value < min || value > max)
                        {
                            passesPhaseSpaceTruthCC0Pi = false;
                            break;
                        }
                    }
                }
                const auto isCC0PiSignal = isTrueCC0Pi && passesPhaseSpaceTruthCC0Pi;

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
                treeWriter.SetOutputBranchAddress("cc0pi_signal", (void*)&dummyBool, "cc0pi_signal/O" );
                treeWriter.SetOutputBranchAddress("cc1pi_signal", (void*)&dummyBool, "cc1pi_signal/O" );
            }

            // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            // Weights
            // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            std::cout<<"DEBUG Point 0"<<std::endl;
            std::map<std::string, const std::vector<double>*> weightPtrMap;
            if(isMC)
            {
                std::cout<<"DEBUG Point 1"<<std::endl;
                // General systematic weights
                for (auto &weights : *pEventPeLEE->truth.weights.GetAddress())
                {
                    weightPtrMap[weights.first] = &weights.second;
                    treeWriter.SetObjectOutputBranchAddress<std::vector<double>>("weight_" + weights.first, weightPtrMap.at(weights.first));
                }

                std::cout<<"DEBUG Point 2"<<std::endl;
                // GENIE weights
                treeWriter.SetOutputBranchAddress("spline_weight", (void*)pEventPeLEE->truth.weightSpline.GetAddress(), "spline_weight/F");
                treeWriter.SetOutputBranchAddress("tuned_cv_weight", (void*)pEventPeLEE->truth.weightTune.GetAddress(), "tuned_cv_weight/F");
                std::cout<<"DEBUG Point 3"<<std::endl;
            }
            else
            {
                std::cout<<"DEBUG Point 4"<<std::endl;
                // Non-mc weights are set to 1
                const auto defaultWeight = 1.f;
                treeWriter.SetOutputBranchAddress("spline_weight", (void*)&defaultWeight, "spline_weight/F");
                treeWriter.SetOutputBranchAddress("tuned_cv_weight", (void*)&defaultWeight, "tuned_cv_weight/F");
                std::cout<<"DEBUG Point 5"<<std::endl;
            }
            std::cout<<"DEBUG Point 6"<<std::endl;

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
            std::cout<<"DEBUG Point 7"<<std::endl;

            treeWriter.CreateNoNewBranches(); // Only creates branches for the first run through the loop
            std::cout<<"DEBUG Point 8"<<std::endl;
            treeWriter.Fill();
            std::cout<<"DEBUG Point 9"<<std::endl;
        } // End of event-level iteration

    } // End of file-level iterration
    std::cout<<"-----------------Done-----------------"<<std::endl;


}

} // namespace ubcc1pi_macros
