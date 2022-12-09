/**
 *  @file  ubcc1pi_standalone/Macros/OptimizeSidebandBinEdges.cxx
 *
 *  @brief The implementation file of the OptimizeSidebandBinEdges macro
 */

#include "ubcc1pi_standalone/Macros/Macros.h"
#include "ubcc1pi_standalone/Objects/FileReader.h"

#include "ubcc1pi_standalone/Helpers/AnalysisHelper.h"
#include "ubcc1pi_standalone/Helpers/NormalisationHelper.h"
#include "ubcc1pi_standalone/Helpers/SelectionHelper.h"
// #include "ubcc1pi_standalone/Helpers/CrossSectionHelper.h"
// #include "ubcc1pi_standalone/Helpers/FormattingHelper.h"
#include "ubcc1pi_standalone/Helpers/ExtractionHelper.h"
#include "ubsmear.h"

using namespace ubcc1pi;

namespace ubcc1pi_macros
{

void OptimizeSidebandBinEdges(const Config &config)
{

    // -------------------------------------------------------------------------------------------------------------------------------------
    // Setup the relevent "getters" for each cross-section
    // -------------------------------------------------------------------------------------------------------------------------------------
    // Here we define a map from the name of each cross-section to a function which pulls out the relevant kinematic quanitity from an input
    // analysis data object. Again this is done up-front to reduce code-bloat below.
    ExtractionHelper::AnalysisValueMap getSidebandValue;
    ExtractionHelper::PopulateAnalysisValueMap(getSidebandValue, true); // true: creates sideband getters


    // We additionally make a map from each cross-section to the limits of the phase-space that we should consider. The key is the
    // identifier for the kinematic quantity and the mapped value is a pair containing the limits [min, max]
    std::map< std::string, std::pair<float, float> > phaseSpaceMap;

    // ATTN the configuration allows the user to enable or disable each cross-section. If a cross-section has been disabled, then it won't
    // be added the xSecMap. However, the phaseSpaceMap always includes all kinematic parameters!

    // Add the differential cross-sections
    for (const auto &[name, binning, scaleByBinWidth] : std::vector< std::tuple<std::string, Config::Global::Binning, bool> > {

        // The names of the cross-section kinematic parameters, and their binning information.
        // The third (boolean) parameter indicates if the cross-section bins should be scaled by their width
        { "muonCosTheta", config.global.muonCosTheta, false  }, // scaleByBinWidth not relevant for fit and incompatible with overflow bin fit.
        { "muonPhi", config.global.muonPhi, false  }, 
        { "muonMomentum", config.global.sidebandMuonMomentum, false  }, // modified version to treat overflow as regular bins

        { "pionCosTheta", config.global.pionCosTheta, false  }, 
        { "pionPhi", config.global.pionPhi, false  }, 
        { "pionMomentum", config.global.sidebandPionMomentum, false  }, // modified version to with different bins

        { "muonPionAngle", config.global.muonPionAngle, false  }, 
        { "nProtons", config.global.nProtons, false }

    })
    {
        // Add to the phase-space map
        phaseSpaceMap.emplace(name, std::pair<float, float>({binning.min, binning.max}));
    }

    std::cout<<"..........................................\nUSING Modified CC0pi Selection: muonLikeProtonValue=-0.48f, barelyResemblingProtonValue=0.12f\n.........................................."<<std::endl;
    auto sidebandSelection = SelectionHelper::GetCC0piSelectionModified(-0.48f, 0.12f);

    // -------------------------------------------------------------------------------------------------------------------------------------
    // Setup the input files
    // -------------------------------------------------------------------------------------------------------------------------------------
    // Here we define a vector of tuples with 4 entries
    //   - First, the sample type (e.g. overlay)
    //   - Second, a string which is used to identify a given detector variation sample (for other sample type, this is unused)
    //   - Third, the path to the input file
    //   - Fourth, the normalisation factor to apply to all events in that file
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
    std::map<std::string, std::vector<std::pair<float, float>>> eventTruthValueMap;
    float eventWeightTotal = 0;
    // std::vector<std::string> names = {"muonCosTheta", "muonPhi", "muonMomentum", "pionCosTheta", "pionPhi", "pionMomentum", "muonPionAngle"};//, "nProtons"};
    std::vector<std::string> names = {"pionMomentum"};//, "nProtons"};
    for(const auto &name: names) eventTruthValueMap.emplace(name, std::vector<std::pair<float, float>>());
    for (const auto &[sampleType, sampleName, fileName, normalisation] : inputData)
    {
        // Open the input file for reading and enable the branches with systematic event weights (if required)
        FileReader reader(fileName);
        auto pEvent = reader.GetBoundEventAddress();
        // Loop over the events in the file
        const auto nEvents = reader.GetNumberOfEvents();

        std::cout<<"############################\nUsing 1/1 of events!\n############################"<<std::endl;
        for (unsigned int i = 0; i < nEvents; ++i)
        {
            AnalysisHelper::PrintLoadingBar(i, nEvents);
            reader.LoadEvent(i);
            //std::cout<<"DEBUG - K-2"<<std::endl;

            // -----------------------------------------------------------------------------------------------------------------------------
            // Work out if this event passed the selection and apply any additional phase-space cuts based on the input binning
            // -----------------------------------------------------------------------------------------------------------------------------

            // Run the selection
            const auto &[passedGoldenSidebandSelection, sidebandCutsPassed, sidebandAssignedPdgCodes] = sidebandSelection.Execute(pEvent);

            if(!passedGoldenSidebandSelection) continue;

            //std::cout<<"DEBUG - K-1.3"<<std::endl;
            const auto passedGenericSidebandSelection = SelectionHelper::IsCutPassed(sidebandCutsPassed, config.global.lastCutGeneric);
            
            // if(!passedGenericSidebandSelection) continue;

            //std::cout<<"DEBUG - K-1.2"<<std::endl;
            // Get the reco analysis data (if available, otherwise set to dummy values)
            const auto recoData = (
                passedGenericSidebandSelection
                    ? AnalysisHelper::GetRecoAnalysisDataCC0Pi(pEvent->reco, sidebandAssignedPdgCodes, passedGoldenSidebandSelection)
                    : AnalysisHelper::GetDummyAnalysisData()
            );

            //std::cout<<"DEBUG - K-1.1"<<std::endl;
            //std::cout<<"DEBUG - K-1"<<std::endl;

            // Here we apply reco-level phase-space restrictions
            // For any event that passes the generic selection, get the value of the kinematic quantity and check if it is outside of the
            // min/max values supplied in the binning. If so, then reject the event.
            bool passesPhaseSpaceReco = false;
            if (passedGenericSidebandSelection)
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
                        break;
                    }
                }
            }

            // const auto isSelectedSidebandGeneric = passedGenericSidebandSelection && passesPhaseSpaceReco;
            // const auto isSelectedSidebandGolden = passedGoldenSidebandSelection && passesPhaseSpaceReco;

            // // std::map<std::string, bool> isSelectedMap = {{"sideband", isSelectedSideband}};
            // std::map<std::string, bool> isSelectedMap = {{"generic", isSelectedSidebandGeneric}, {"golden", isSelectedSidebandGolden}};
            //std::cout<<"DEBUG - K1"<<std::endl;
            if(!passesPhaseSpaceReco) continue;

            // -----------------------------------------------------------------------------------------------------------------------------
            // Work out if this event is signal, and apply any phase-space restrictions based on the input binning
            // -----------------------------------------------------------------------------------------------------------------------------

            // Determine if this is truly a CC0Pi event
            const auto isTrueCC0Pi = AnalysisHelper::IsTrueCC0Pi(pEvent, config.global.useAbsPdg, config.global.protonMomentumThreshold);

            if(!isTrueCC0Pi) continue;

            // Get the truth analysis data (if available, otherwise set to dummy values)
            const auto truthData = (
                (isTrueCC0Pi)
                    ? AnalysisHelper::GetTruthAnalysisDataCC0Pi(pEvent->truth, config.global.useAbsPdg, config.global.protonMomentumThreshold)
                    : AnalysisHelper::GetDummyAnalysisData()
            );

            // Here we apply truth-level phase-space restrictions
            // For all true CC1Pi events, we check if the values of each kinematic variable are within the supplied limits. If not then the
            // event is not classed as "signal"
            //std::cout<<"DEBUG - K2"<<std::endl;
            bool passesPhaseSpaceTruth = false;
            if (isTrueCC0Pi)
            {
                // Start by assuming the event passes the phase-space cuts
                passesPhaseSpaceTruth = true;

                // Check the value of the kinematic quantities are within the phase-space limits
                for (const auto &[name, minMax] : phaseSpaceMap)
                {
                    // if(name == "pionMomentum") continue; // Not really compatible with the pion momentum in the CC1pi selection // todo check this
                    const auto &[min, max] = minMax;
                    const auto value = getSidebandValue.at(name)(truthData);

                    if (value < min || value > max)
                    {
                        passesPhaseSpaceTruth = false;
                        break;
                    }
                }
            }

            if(!passesPhaseSpaceTruth) continue;

            // const auto recoValue = getSidebandValue.at(name)(recoData);
            //std::cout<<"DEBUG - K3"<<std::endl;
            const auto weight = AnalysisHelper::GetNominalEventWeight(pEvent) * normalisation;
            for(const auto &name : names)
            {
                const auto trueValue = getSidebandValue.at(name)(truthData);
                eventTruthValueMap[name].push_back(std::make_pair(trueValue, weight));
            }
            //std::cout<<"DEBUG - K4"<<std::endl;
            eventWeightTotal += weight;
        }
    }

    std::cout<<"eventWeightTotal: "<<eventWeightTotal<<std::endl;

    //std::cout<<"DEBUG - K5"<<std::endl;
    std::map<std::string, Config::Global::Binning> binningMap{
        // The names of the cross-section kinematic parameters, and their binning information.
        // { "muonCosTheta", config.global.muonCosTheta },
        // { "muonPhi", config.global.muonPhi}, 
        // { "muonMomentum", config.global.sidebandMuonMomentum },
        // { "pionCosTheta", config.global.pionCosTheta }, 
        // { "pionPhi", config.global.pionPhi }, 
        { "pionMomentum", config.global.sidebandPionMomentum },
        // { "muonPionAngle", config.global.muonPionAngle }
    };
    // { "nProtons", config.global.nProtons }

    //std::cout<<"DEBUG - K6"<<std::endl;
    for (const auto &[name, binning]: binningMap)
    {
        //std::cout<<"DEBUG - K7"<<std::endl;
        const auto minEdge = binning.min;
        const auto maxEdge = binning.max;
        const auto diff = maxEdge - minEdge;
        const auto nBins = binning.binEdges.size() - 1;// ATTN: Does not consider under/overflow bins
        //std::cout<<"DEBUG - K8"<<std::endl;
        // const auto binEdges = PlottingHelper::GenerateUniformBinEdges(diff/0.01, minEdge, maxEdge);
        //std::cout<<"DEBUG - K9"<<std::endl;
        auto lowerEdge = minEdge;
        auto upperEdge = minEdge+0.01;
        std::cout<<name<<" "<<minEdge;
        // //std::cout<<"DEBUG - K10"<<std::endl;
        unsigned int edgeNum = 1;
        while(edgeNum!=nBins)
        {
            // //std::cout<<"DEBUG - K11 - "<<upperEdge<<"/"<<maxEdge<<std::endl;
            auto weightSum = 0.f;
            for(const auto &[value, weight]: eventTruthValueMap[name])
            {
                if(value >= lowerEdge && value < upperEdge) weightSum += weight;
            }
            // std::cout<<weightSum<<"/"<<eventWeightTotal<<" - "<<lowerEdge<<" -> "<<upperEdge<<std::endl;
            if(weightSum>=eventWeightTotal/nBins)
            {
                std::cout<<" "<<upperEdge<<"("<<weightSum<<")"<<std::endl;
                edgeNum += 1;
                lowerEdge = upperEdge;
            }
            // //std::cout<<"DEBUG - K12"<<std::endl;
            upperEdge += 0.01;
        }
        std::cout<<" "<<maxEdge<<std::endl;
    }

    std::cout<<"------------- All Done -------------"<<std::endl;
    // return;
}

} // namespace ubcc1pi_macros
