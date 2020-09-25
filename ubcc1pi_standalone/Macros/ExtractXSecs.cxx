/**
 *  @file  ubcc1pi_standalone/Macros/ExtractXSecs.cxx
 *
 *  @brief The implementation file of the ExtractXSecs macro
 */

#include "ubcc1pi_standalone/Macros/Macros.h"
#include "ubcc1pi_standalone/Objects/FileReader.h"

#include "ubcc1pi_standalone/Helpers/AnalysisHelper.h"
#include "ubcc1pi_standalone/Helpers/SelectionHelper.h"
#include "ubcc1pi_standalone/Helpers/CrossSectionHelper.h"
#include "ubcc1pi_standalone/Helpers/NormalisationHelper.h"
#include "ubcc1pi_standalone/Helpers/FormattingHelper.h"

using namespace ubcc1pi;

namespace ubcc1pi_macros
{

void ExtractXSecs(const Config &config)
{
    //
    // Setup the input files
    // 
    std::vector< std::tuple<AnalysisHelper::SampleType, std::string, float> > inputData;
    
    inputData.emplace_back(AnalysisHelper::Overlay, config.files.overlaysFileName, NormalisationHelper::GetOverlaysNormalisation(config)); 
    inputData.emplace_back(AnalysisHelper::Dirt,    config.files.dirtFileName,     NormalisationHelper::GetDirtNormalisation(config)); 
    inputData.emplace_back(AnalysisHelper::DataEXT, config.files.dataEXTFileName,  NormalisationHelper::GetDataEXTNormalisation(config)); 
    inputData.emplace_back(AnalysisHelper::DataBNB, config.files.dataBNBFileName,  1.f);


    //
    // Load the first overlay event so we can store the names of the systematic parameters to apply and the numbers of universes
    // 
    FileReader overlayReader(config.files.overlaysFileName);
    if (overlayReader.GetNumberOfEvents() == 0)
        throw std::invalid_argument("ExtractXSecs - No overlay events!");

    overlayReader.EnableSystematicBranches();
    auto pOverlayEvent = overlayReader.GetBoundEventAddress();
    overlayReader.LoadEvent(0);

    // Get the dimensions of the systematic weights (without the boostrap weights added)
    auto systWeightsMap_firstEvent = CrossSectionHelper::GetSystematicWeightsMap(pOverlayEvent->truth);
    const auto systUniverseSizesMapNoBootstrap = CrossSectionHelper::GetSystematicUniverseSizesMap(systWeightsMap_firstEvent);
   
    // Add the boostrap weights, and get the dimensions again
    // TODO make the number of boostrap universes configurable
    const unsigned int nBootstrapUniverses = 1000u;
    CrossSectionHelper::AddBootstrapWeights(nBootstrapUniverses, systWeightsMap_firstEvent);
    const auto systUniverseSizesMap = CrossSectionHelper::GetSystematicUniverseSizesMap(systWeightsMap_firstEvent);

    // Print the systematic parameters we are going to apply
    FormattingHelper::Table systParamsTable({"Parameter", "Universes"});
    for (const auto &[name, nUniverses] : systUniverseSizesMap)
    {
        systParamsTable.AddEmptyRow();
        systParamsTable.SetEntry("Parameter", name);
        systParamsTable.SetEntry("Universes", nUniverses);
    }
    systParamsTable.Print();


    //
    // Setup the cross-section objects
    //
    
    // Muon momentum
    std::vector<float> muonMomentum_extendedBinEdges;
    bool muonMomentum_hasUnderflow, muonMomentum_hasOverflow;
    CrossSectionHelper::GetExtendedBinEdges(
        config.global.muonMomentum.min, config.global.muonMomentum.max, config.global.muonMomentum.binEdges,
        muonMomentum_extendedBinEdges, muonMomentum_hasUnderflow, muonMomentum_hasOverflow);

    CrossSectionHelper::CrossSection xsec_muonMomentum(muonMomentum_extendedBinEdges, muonMomentum_hasUnderflow, muonMomentum_hasOverflow, true, systUniverseSizesMap);

    // Get the selection
    auto selection = SelectionHelper::GetDefaultSelection();
    
    // Loop over the events
    for (const auto [sampleType, fileName, normalisation] : inputData)
    {
        std::cout << "Reading input file: " << fileName << std::endl;

        FileReader reader(fileName);
        reader.EnableSystematicBranches();
        auto pEvent = reader.GetBoundEventAddress();
        const auto isOverlay = (sampleType == AnalysisHelper::Overlay);
        const auto isDataBNB = (sampleType == AnalysisHelper::DataBNB);

        const auto nEvents = reader.GetNumberOfEvents();
        for (unsigned int i = 0; i < nEvents; ++i)
        {
            AnalysisHelper::PrintLoadingBar(i, nEvents);
            reader.LoadEvent(i);
            
            // Get the nominal event weight, scaled by the sample normalisation
            const auto weight = AnalysisHelper::GetNominalEventWeight(pEvent) * normalisation;

            // Run the event selection and store which cuts are passed
            std::vector<std::string> cutsPassed;
            std::vector<int> assignedPdgCodes;
            auto passedGoldenSelection = selection.Execute(pEvent, cutsPassed, assignedPdgCodes);
            auto passedGenericSelection = (std::find(cutsPassed.begin(), cutsPassed.end(), config.global.lastCutGeneric) != cutsPassed.end());
            
            // Determine if the phase-space restrictions are met in reco
            auto recoData = AnalysisHelper::GetDummyAnalysisData();
            bool passesPhaseSpaceReco = true;
            if (passedGenericSelection)
            {
                // Get the reconstructed kinematic parameters
                recoData = AnalysisHelper::GetRecoAnalysisData(pEvent->reco, assignedPdgCodes, passedGoldenSelection);

                if (recoData.muonMomentum < config.global.muonMomentum.min || recoData.muonMomentum > config.global.muonMomentum.max)
                    passesPhaseSpaceReco = false;
            }

            // Apply the phase-space restrictions in reco as an additional event selection cut
            const bool isSelectedGeneric = passedGenericSelection && passesPhaseSpaceReco;
            const bool isSelectedGolden = passedGoldenSelection && passesPhaseSpaceReco;
 
            // Determine if the event is truly a signal event
            const auto isTrueCC1Pi = isOverlay && AnalysisHelper::IsTrueCC1Pi(pEvent, config.global.useAbsPdg);

            // Determine of the phase-space restrictions are met in truth
            auto truthData = AnalysisHelper::GetDummyAnalysisData();
            bool passesPhaseSpaceTruth = true;
            if (isTrueCC1Pi)
            {
                // Get the true kinematic parameters
                truthData = AnalysisHelper::GetTruthAnalysisData(pEvent->truth, config.global.useAbsPdg, config.global.protonMomentumThreshold);

                if (truthData.muonMomentum < config.global.muonMomentum.min || truthData.muonMomentum > config.global.muonMomentum.max)
                    passesPhaseSpaceTruth = false;
            }

            // A signal event is a true CC1Pi event that passes the phase space restrictions
            const bool isSignal = isTrueCC1Pi && passesPhaseSpaceTruth;

            // We only need to count signal events, or any event that passes the generic selection - everything else we can savely skip
            const bool shouldCount = isSignal || isSelectedGeneric;
            if (!shouldCount)
                continue;

            // Count the BNB data event
            if (isDataBNB)
            {
                if (isSelectedGeneric)
                {
                    xsec_muonMomentum.AddSelectedBNBDataEvent(recoData.muonMomentum);
                }
            }
            else
            {
                // Count the signal event
                if (isSignal)
                {
                    // Get the systematic event weights
                    auto systWeightsMap = CrossSectionHelper::GetSystematicWeightsMap(pOverlayEvent->truth);
                    CrossSectionHelper::AddBootstrapWeights(nBootstrapUniverses, systWeightsMap);

                    xsec_muonMomentum.AddSignalEvent(truthData.muonMomentum, recoData.muonMomentum, isSelectedGeneric, weight, systWeightsMap);

                }
                // Count the background event
                else
                {
                    // Get the systematic event weights
                    auto systWeightsMap = CrossSectionHelper::GetUnitSystematicWeightsMap(systUniverseSizesMapNoBootstrap);
                    CrossSectionHelper::AddBootstrapWeights(nBootstrapUniverses, systWeightsMap);

                    if (isSelectedGeneric)
                    {
                        xsec_muonMomentum.AddSelectedBackgroundEvent(recoData.muonMomentum, weight, systWeightsMap);
                    }
                }
            }
        }
    }

    // Calculate the cross-section and print it out!
    const auto xsec_muonMomentum_plot = xsec_muonMomentum.GetCrossSection();
    for (unsigned int iBin = 1u; iBin <= static_cast<unsigned int>(xsec_muonMomentum_plot->GetXaxis()->GetNbins()); ++iBin)
    {
        std::cout << iBin << " : " << xsec_muonMomentum_plot->GetBinContent(iBin) << std::endl;
    }
}

} // namespace ubcc1pi_macros
