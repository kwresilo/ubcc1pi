/**
 *  @file  ubcc1pi_standalone/Macros/PlotSystematicWeights.cxx
 *
 *  @brief The implementation file of the PlotSystematicWeights macro
 */

#include "ubcc1pi_standalone/Macros/Macros.h"
#include "ubcc1pi_standalone/Objects/FileReader.h"

#include "ubcc1pi_standalone/Helpers/AnalysisHelper.h"
#include "ubcc1pi_standalone/Helpers/SelectionHelper.h"
#include "ubcc1pi_standalone/Helpers/CrossSectionHelper.h"

using namespace ubcc1pi;

namespace ubcc1pi_macros
{

void PlotSystematicWeights(const Config &config)
{
    // Open the input file for reading
    FileReader reader(config.files.overlaysFileName);
    reader.EnableSystematicBranches();

    auto pEvent = reader.GetBoundEventAddress();
    const auto nEvents = reader.GetNumberOfEvents();

    // Get the selection
    auto selection = SelectionHelper::GetDefaultSelection();

    // Setup a mapping with index of [systParamName] to the min/max of the weight in that universe
    std::map<std::string, std::pair<float, float> > minMaxMap;

    // Setup a mapping from [systParamName] to the sum of square of the difference of the weights from unity
    std::map<std::string, float> sqrDiffSumMap;

    // Setup a mapping from [systParamName] to the total number of weights used
    std::map<std::string, float> totalMap;

    // Get the min/max values for the plot limits
    unsigned int nBadWeights = 0u;
    std::cout << "Getting extremal weights" << std::endl;
    std::vector<unsigned int> fiducialEventIndices;

    for (unsigned int i = 0; i < nEvents; ++i)
    {
        AnalysisHelper::PrintLoadingBar(i, nEvents);
        reader.LoadEvent(i);

        // Insist the true neutrino is fiducial
        if (!AnalysisHelper::IsFiducial(pEvent->truth.nuVertex()))
            continue;

        fiducialEventIndices.push_back(i);

        // Get the systematic event weights
        CrossSectionHelper::SystematicWeightsMap systWeightsMap;
        CrossSectionHelper::AddSystematicWeights(pEvent->truth, config.extractXSecs.systematicParams, systWeightsMap);
        CrossSectionHelper::AddMutuallyExclusiveSystematicWeights(pEvent->truth, config.extractXSecs.mutuallyExclusiveSystematicParams, systWeightsMap);
        CrossSectionHelper::AddBootstrapWeights(config.extractXSecs.nBootstrapUniverses, systWeightsMap);

        // Add the nominal event weight
        const auto weight = AnalysisHelper::GetNominalEventWeight(pEvent);
        systWeightsMap.emplace("nominal", std::vector<float>({weight}));

        // Print out the weights from the first event
        if (fiducialEventIndices.size() == 1)
        {
            std::cout << "Systematic weights: " << systWeightsMap.size() << std::endl;
            for (const auto &[param, weights] : systWeightsMap)
                std::cout << param << " - " << weights.size() << std::endl;
        }

        // Get the minimum and maximum weights for each parameter
        for (const auto &[param, weights] : systWeightsMap)
        {
            if (weights.empty())
                continue;

            // Filter the weights to remove unphysical values
            std::vector<float> filteredWeights(weights.size());
            auto itFilter = std::copy_if(weights.begin(), weights.end(), filteredWeights.begin(), [](const float &x) {
                return (x < (std::numeric_limits<float>::max() - std::numeric_limits<float>::epsilon()));
            });
            filteredWeights.resize(std::distance(filteredWeights.begin(), itFilter));

            // Count the number of bad weights
            nBadWeights += (weights.size() - filteredWeights.size());

            const auto minWeight = *std::min_element(filteredWeights.begin(), filteredWeights.end());
            const auto maxWeight = *std::max_element(filteredWeights.begin(), filteredWeights.end());

            float sqrDiffSum = 0.f;
            float totalSum = 0.f;
            for (const auto &weight : filteredWeights)
            {
                // Skip weights equal to exactly one - we want the spread of the distribution
                if (std::abs(weight - 1.f) <= std::numeric_limits<float>::epsilon())
                    continue;

                sqrDiffSum += std::pow(1 - weight, 2);
                totalSum += 1.f;
            }

            auto iter = minMaxMap.find(param);
            if (iter == minMaxMap.end())
            {
                minMaxMap.emplace(param, std::pair<float, float>(minWeight, maxWeight));
                sqrDiffSumMap.emplace(param, sqrDiffSum);
                totalMap.emplace(param, totalSum);
            }
            else
            {
                auto &minMaxPair = iter->second;
                auto &min = minMaxPair.first;
                auto &max = minMaxPair.second;

                min = std::min(min, minWeight);
                max = std::max(max, maxWeight);

                sqrDiffSumMap.at(param) += sqrDiffSum;
                totalMap.at(param) += totalSum;
            }
        }
    }

    std::cout << "Number of bad weights (== inf): " << nBadWeights << std::endl;

    // Setup a mapping with index of [systParamName] to plot for that parameter
    std::map<std::string, std::shared_ptr<PlottingHelper::MultiPlot> > plotMap, noOnesPlotMap, selectedPlotMap;
    for (const auto &[param, minMaxPair] : minMaxMap)
    {
        // Get the RMS difference of the weights from unity
        const auto rms = std::pow(sqrDiffSumMap.at(param) / totalMap.at(param), 0.5f);

        std::cout << param << " : " << minMaxPair.first << " -> " << minMaxPair.second << " (" << rms << ")" << std::endl;

        const auto spread = 4.f * rms; // 4 is an abitrary choice to get reasonable looking plots
        const auto min = std::max(1 - spread, minMaxPair.first);
        const auto max = std::min(1 + spread, minMaxPair.second);

        plotMap.emplace(param, std::make_shared<PlottingHelper::MultiPlot>(param, "Number of weights", 100, min, max));
        noOnesPlotMap.emplace(param, std::make_shared<PlottingHelper::MultiPlot>(param, "Number of weights", 100, min, max));
        selectedPlotMap.emplace(param, std::make_shared<PlottingHelper::MultiPlot>(param, "Number of weights", 100, min, max));
    }

    // Loop over the events
    std::cout << "Filling plots" << std::endl;
    for (const auto i : fiducialEventIndices)
    {
        AnalysisHelper::PrintLoadingBar(i, nEvents);
        reader.LoadEvent(i);

        // Get the event type
        const auto style = PlottingHelper::GetPlotStyle(AnalysisHelper::Overlay, pEvent, config.global.useAbsPdg);

        // Run the event selection and store which cuts are passed
        std::vector<std::string> cutsPassed;
        std::vector<int> assignedPdgCodes;
        auto passedGoldenSelection = selection.Execute(pEvent, cutsPassed, assignedPdgCodes);
        auto passedGenericSelection = (std::find(cutsPassed.begin(), cutsPassed.end(), config.global.lastCutGeneric) != cutsPassed.end());

        // Get the nominal event weight
        const auto weight = AnalysisHelper::GetNominalEventWeight(pEvent);

        // Add this to the plot map
        plotMap.at("nominal")->Fill(weight, style);

        if (std::abs(weight - 1.f) > std::numeric_limits<float>::epsilon())
        {
            noOnesPlotMap.at("nominal")->Fill(weight, style);
        }

        if (passedGenericSelection)
        {
            selectedPlotMap.at("nominal")->Fill(weight, style);
        }

        // Get the systematic event weights
        CrossSectionHelper::SystematicWeightsMap systWeightsMap;
        CrossSectionHelper::AddSystematicWeights(pEvent->truth, config.extractXSecs.systematicParams, systWeightsMap);
        CrossSectionHelper::AddMutuallyExclusiveSystematicWeights(pEvent->truth, config.extractXSecs.mutuallyExclusiveSystematicParams, systWeightsMap);
        CrossSectionHelper::AddBootstrapWeights(config.extractXSecs.nBootstrapUniverses, systWeightsMap);

        // Add these to the total
        for (const auto &[param, weights] : systWeightsMap)
        {
            for (const auto &weight : weights)
            {
                plotMap.at(param)->Fill(weight, style);

                if (std::abs(weight - 1.f) > std::numeric_limits<float>::epsilon())
                {
                    noOnesPlotMap.at(param)->Fill(weight, style);
                }

                if (passedGenericSelection)
                {
                    selectedPlotMap.at(param)->Fill(weight, style);
                }
            }
        }
    }

    // Save the plots
    for (const auto &[param, pPlot] : plotMap)
    {
        pPlot->SaveAsStacked("plotSystematicWeights_" + param);
    }

    for (const auto &[param, pPlot] : noOnesPlotMap)
    {
        pPlot->SaveAsStacked("plotSystematicWeights_noOnes_" + param);
    }

    for (const auto &[param, pPlot] : selectedPlotMap)
    {
        pPlot->SaveAsStacked("plotSystematicWeights_selected_" + param);
    }
}

} // namespace ubcc1pi_macros
