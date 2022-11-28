/**
 *  @file  ubcc1pi_standalone/Macros/ExtractNuWroXSecs2.cxx
 *
 *  @brief The implementation file of the ExtractNuWroXSecs2 macro
 */

#include "ubcc1pi_standalone/Macros/Macros.h"
#include "ubcc1pi_standalone/Objects/FileReader.h"

#include "ubcc1pi_standalone/Helpers/AnalysisHelper.h"
#include "ubcc1pi_standalone/Helpers/NormalisationHelper.h"
// #include "ubcc1pi_standalone/Helpers/SelectionHelper.h"
// #include "ubcc1pi_standalone/Helpers/CrossSectionHelper.h"
// #include "ubcc1pi_standalone/Helpers/FormattingHelper.h"
// #include "ubcc1pi_standalone/Helpers/FittingHelper.h"

using namespace ubcc1pi;

namespace ubcc1pi_macros
{

void QuerryWeights(const Config &config)
{
    CrossSectionHelper::CrossSection::SystParams systParams;
    systParams.nBootstrapUniverses = config.extractXSecs.nBootstrapUniverses;
    systParams.fluxDimensions = config.extractXSecs.fluxDimensions;
    systParams.xsecDimensions = config.extractXSecs.xsecDimensions;
    systParams.reintDimensions = config.extractXSecs.reintDimensions;
    systParams.detVarDimensions = config.extractXSecs.detVarDimensions;
    systParams.potFracUncertainty = config.extractXSecs.potFracUncertainty;


    std::vector< std::tuple<AnalysisHelper::SampleType, std::string, std::string, float> > inputData;
    if(std::binary_search(config.global.runs.begin(), config.global.runs.end(), 1))
    {
        inputData.emplace_back(AnalysisHelper::Overlay, "", config.filesRun1.overlaysFileName, NormalisationHelper::GetOverlaysNormalisation(config, 1));
    }
    if(std::binary_search(config.global.runs.begin(), config.global.runs.end(), 2))
    {
        inputData.emplace_back(AnalysisHelper::Overlay, "", config.filesRun2.overlaysFileName, NormalisationHelper::GetOverlaysNormalisation(config, 2));
    }
    if(std::binary_search(config.global.runs.begin(), config.global.runs.end(), 3))
    {
        inputData.emplace_back(AnalysisHelper::Overlay, "", config.filesRun3.overlaysFileName, NormalisationHelper::GetOverlaysNormalisation(config, 3));
    }


    // -------------------------------------------------------------------------------------------------------------------------------------
    // Count the events for CC1pi
    // -------------------------------------------------------------------------------------------------------------------------------------
    // Loop over the files
    for (const auto &[sampleType, sampleName, fileName, normalisation] : inputData)
    {
        std::cout << "Reading input file: " << fileName << std::endl;

        const auto isOverlay = (sampleType == AnalysisHelper::Overlay);
        // Open the input file for reading and enable the branches with systematic event weights (if required)
        FileReader reader(fileName);

        if (isOverlay)// || isNuWro)
            reader.EnableSystematicBranches();
        
        // if (isDetVar) continue; // TODO: Remove this after debugging !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        // if(isDataBNB || isDirt || isDataEXT) continue; // Todo: Remove !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        auto pEvent = reader.GetBoundEventAddress();

        // Loop over the events in the file
        const auto nEvents = reader.GetNumberOfEvents();
        std::vector<float> summedValues(2, 0.f);

        for (unsigned int i = 0; i < nEvents; ++i) // todo change back to every event!!!!!!!!!!!!!!!!!!!!!!
        {
            // AnalysisHelper::PrintLoadingBar(i, nEvents);
            reader.LoadEvent(i);
            // -----------------------------------------------------------------------------------------------------------------------------
            // Work out if this event passed the selection and apply any additional phase-space cuts based on the input binning
            // -----------------------------------------------------------------------------------------------------------------------------
            // Run the selection

            // const auto &[passedGoldenSelection, cutsPassed, assignedPdgCodes] = selection.Execute(pEvent);
            // const auto passedGenericSelection = SelectionHelper::IsCutPassed(cutsPassed, config.global.lastCutGeneric);

            // // Get the reco analysis data (if available, otherwise set to dummy values)
            // const auto recoData = (
            //     passedGenericSelection
            //         ? AnalysisHelper::GetRecoAnalysisData(pEvent->reco, assignedPdgCodes, passedGoldenSelection)
            //         : AnalysisHelper::GetDummyAnalysisData()
            // );

            // // Here we apply reco-level phase-space restrictions
            // // For any event that passes the generic selection, get the value of the kinematic quantity and check if it is outside of the
            // // min/max values supplied in the binning. If so, then reject the event.
            // bool passesPhaseSpaceReco = false;
            // if (passedGenericSelection)
            // {
            //     // Start by assuming the event passes the phase-space cuts
            //     passesPhaseSpaceReco = true;

            //     // Check the value of the kinematic quantities are within the phase-space limits
            //     for (const auto &[name, minMax] : phaseSpaceMap)
            //     {
            //         const auto &[min, max] = minMax;
            //         const auto value = getValue.at(name)(recoData);

            //         if (value < min || value > max)
            //         {
            //             passesPhaseSpaceReco = false;
            //             break;
            //         }
            //     }
            // }

            // const auto isSelectedGolden = passedGoldenSelection && passesPhaseSpaceReco;
            // const auto isSelectedGeneric = passedGenericSelection && passesPhaseSpaceReco;
            // std::map<std::string, bool> isSelectedMap = {{"generic",isSelectedGeneric},{"golden",isSelectedGolden}};

            // // Determine if this is truly a CC1Pi event
            // const auto isTrueCC1Pi = (isOverlay || isDetVar || isNuWro) && AnalysisHelper::IsTrueCC1Pi(pEvent, config.global.useAbsPdg);
            // // Get the truth analysis data (if available, otherwise set to dummy values)
            // const auto truthData = (
            //     isTrueCC1Pi
            //         ? AnalysisHelper::GetTruthAnalysisData(pEvent->truth, config.global.useAbsPdg, config.global.protonMomentumThreshold)
            //         : AnalysisHelper::GetDummyAnalysisData()
            // );

            // // Here we apply truth-level phase-space restrictions
            // // For all true CC1Pi events, we check if the values of each kinematic variable are within the supplied limits. If not then the
            // // event is not classed as "signal"
            // bool passesPhaseSpaceTruth = false;
            // if (isTrueCC1Pi)
            // {
            //     // Start by assuming the event passes the phase-space cuts
            //     passesPhaseSpaceTruth = true;

            //     // Check the value of the kinematic quantities are within the phase-space limits
            //     for (const auto &[name, minMax] : phaseSpaceMap)
            //     {
            //         const auto &[min, max] = minMax;
            //         const auto value = getValue.at(name)(truthData);

            //         if (value < min || value > max)
            //         {
            //             passesPhaseSpaceTruth = false;
            //             break;
            //         }
            //     }
            // }

            // const auto isCC1PiSignal = isTrueCC1Pi && passesPhaseSpaceTruth;
            // const auto isSignal = isCC1PiSignal;

            // // -----------------------------------------------------------------------------------------------------------------------------
            // // Work out if this event is signal, and apply any phase-space restrictions based on the input binning
            // // -----------------------------------------------------------------------------------------------------------------------------

            // // Get the nominal event weight, scaled by the sample normalisation
            // const auto weight = AnalysisHelper::GetNominalEventWeight(pEvent) * normalisation;
            
            // // Determine if this is truly a CC0Pi event
            // const auto isTrueCC0Pi = (isOverlay || isDetVar || isNuWro) && AnalysisHelper::IsTrueCC0Pi(pEvent, config.global.useAbsPdg, config.global.protonMomentumThreshold);

            // // Get the truth analysis data (if available, otherwise set to dummy values)
            // const auto sidebandTruthData = (
            //     isTrueCC0Pi
            //         ? AnalysisHelper::GetTruthAnalysisDataCC0Pi(pEvent->truth, config.global.useAbsPdg, config.global.protonMomentumThreshold)
            //         : AnalysisHelper::GetDummyAnalysisData()
            // );

            // bool passesSidebandPhaseSpaceTruth = false;
            // if (isTrueCC0Pi)
            // {
            //     // Start by assuming the event passes the phase-space cuts
            //     passesSidebandPhaseSpaceTruth = true;

            //     // Check the value of the kinematic quantities are within the phase-space limits
            //     for (const auto &[name, minMax] : phaseSpaceMap)
            //     {
            //         const auto &[min, max] = minMax;
            //         const auto value = getSidebandValue.at(name)(sidebandTruthData);

            //         if (value < min || value > max)
            //         {
            //             if(name == "pionMomentum") continue; // Not really compatible with the pion momentum in the CC1pi selection
            //             passesSidebandPhaseSpaceTruth = false;
            //             break;
            //         }
            //     }
            // }

            // const auto isCC0PiSignal = isTrueCC0Pi && passesSidebandPhaseSpaceTruth;

            // -----------------------------------------------------------------------------------------------------------------------------
            // Handle all other events (i.e those from the nominal simulation): Overlays, dirt, EXT data
            // -----------------------------------------------------------------------------------------------------------------------------
            // std::cout<<"### Debug Point 0"<<std::endl;
            // Get the flux weights
            // auto fluxWeights = (
            //     isOverlay
            //         ? CrossSectionHelper::GetWeightsMap(pEvent->truth, systParams.fluxDimensions, config.extractXSecs.mutuallyExclusiveDimensions)
            //         : CrossSectionHelper::GetUnitWeightsMap(systParams.fluxDimensions)
            // );

            // // Fill the flux-reweightor with all overlay events from a neutrinos of the desired flavour in the active volume
            // if (isOverlay && AnalysisHelper::IsInActiveVolume(pEvent->truth.nuVertex()) &&
            //     std::find(config.flux.nuPdgsSignal.begin(), config.flux.nuPdgsSignal.end(), pEvent->truth.nuPdgCode()) != config.flux.nuPdgsSignal.end())
            // {
            //     scalingData.pFluxReweightor->AddEvent(pEvent->truth.nuEnergy(), weight, fluxWeights);
            // }

            // // Fill the flux-reweightor with all overlay events from a neutrinos of the desired flavour in the active volume
            // if (isNuWro && AnalysisHelper::IsInActiveVolume(pEvent->truth.nuVertex()) &&
            //     std::find(config.flux.nuPdgsSignal.begin(), config.flux.nuPdgsSignal.end(), pEvent->truth.nuPdgCode()) != config.flux.nuPdgsSignal.end())
            // {
            //     scalingDataNuWroTruth.pFluxReweightor->AddEvent(pEvent->truth.nuEnergy(), weight, fluxWeights);
            // }

            // Get the cross-section weights
            // ATTN here we optionally scale the cross-section weights down by the genieTuneEventWeight - this is done so we don't
            // double count this weight (once in the nominal event weight, and once in the xsec systematic event weights)
            // const auto xsecWeightsScaleFactor = 1.f;
            // std::cout<<"Before xsecWeightsScaleFactor"<<std::endl;




            const auto isTrueCC0Pi = AnalysisHelper::IsTrueCC0Pi(pEvent, config.global.useAbsPdg, config.global.protonMomentumThreshold);
            // if(!isTrueCC0Pi) continue;

            // auto xsecWeightsScaleFactor = (isOverlay && config.extractXSecs.scaleXSecWeights) ? pEvent->truth.genieTuneEventWeight() : 1.f;
            // xsecWeightsScaleFactor = std::max(xsecWeightsScaleFactor, 0.0001f);
            auto xsecWeightsScaleFactor = 1.f;
            // std::cout<<"After xsecWeightsScaleFactor"<<std::endl;

            // std::cout<<"xsecWeightsScaleFactor Debugging Point 1 - xsecWeightsScaleFactor: "<<xsecWeightsScaleFactor<<std::endl;

            auto xsecWeights = (
                isOverlay
                    ? CrossSectionHelper::ScaleWeightsMap(CrossSectionHelper::GetWeightsMap(pEvent->truth, systParams.xsecDimensions, config.extractXSecs.mutuallyExclusiveDimensions), xsecWeightsScaleFactor)
                    : CrossSectionHelper::GetUnitWeightsMap(systParams.xsecDimensions)
            );

            const auto weightVector = xsecWeights.at("DecayAngMEC_UBGenie");


            auto largeWeight = false;
            for(const auto& weight : weightVector)
            {
                if(weight>10) largeWeight = true;
            }
            if(largeWeight)
            {
                for(const auto& weight : weightVector)
                {
                    if(weight) std::cout<<weight<<",";
                }
                std::cout<<std::endl;
            }

            // auto show = false;
            // // for(const auto &weight : weightVector) if(weight > 2.f || weight < 0.5f) show = true;
            // // if(xsecWeightsScaleFactor > 2.f || xsecWeightsScaleFactor < 0.5f) show = true;
            // if(std::abs(weightVector.at(0)-weightVector.at(1))>0.01 && isTrueCC0Pi) show = true;
            // if(show)
            // {
            //     std::cout<<i<<":  isTrueCC0Pi: "<<isTrueCC0Pi<<" xsecWeightsScaleFactor: "<<xsecWeightsScaleFactor<<" VecFFCCQEshape_UBGenie: ";
            //     for(const auto &weight : weightVector) std::cout<<weight<<" ";
            //     std::cout<<std::endl;
            //     summedValues[0] += weightVector.at(0);
            //     summedValues[1] += weightVector.at(1);
            // }
        }
        std::cout<<"FileName: "<<fileName<<" summedValues[0]: "<<summedValues.at(0)<<" summedValues[1]: "<<summedValues.at(1)<<std::endl;
    }

    std::cout<<"------------- All Done -------------"<<std::endl;
    // return;
}

} // namespace ubcc1pi_macros
