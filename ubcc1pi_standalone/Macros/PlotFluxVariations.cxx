/**
 *  @file  ubcc1pi_standalone/Macros/PlotFluxVariations.cxx
 *
 *  @brief The implementation file of the PlotFluxVariations macro
 */

#include "ubcc1pi_standalone/Macros/Macros.h"
#include "ubcc1pi_standalone/Objects/FileReader.h"

#include "ubcc1pi_standalone/Helpers/AnalysisHelper.h"
#include "ubcc1pi_standalone/Helpers/NormalisationHelper.h"
#include "ubcc1pi_standalone/Helpers/GeometryHelper.h"
#include "ubcc1pi_standalone/Helpers/CrossSectionHelper.h"
#include "ubcc1pi_standalone/Helpers/FormattingHelper.h"
#include "ubcc1pi_standalone/Helpers/PlottingHelper.h"

#include <regex>

using namespace ubcc1pi;

namespace ubcc1pi_macros
{

void PlotFluxVariations(const Config &config)
{
    // Setup the dimensions
    CrossSectionHelper::SystDimensionsMap fluxDimensions;
    for (const auto &[paramName, nUniverses] : config.extractXSecs.fluxDimensions)
    {
        // If this is a parameter that's stored directly in the event weights, then just use it
        if (config.extractXSecs.mutuallyExclusiveDimensions.count(paramName) == 0)
        {
            fluxDimensions.emplace(paramName, nUniverses);
            continue;
        }

        // Otherwise, we have a parameter built from a set of mutually exclusive parameters
        // For the sake of these validation plots, we will use each parameter individually
        const auto &[mutuallyExclusiveParams, nUniversesCheck] = config.extractXSecs.mutuallyExclusiveDimensions.at(paramName);
        if (nUniverses != nUniversesCheck)
            throw std::invalid_argument("PlotFluxVariations - Wrong number of universes for mutually exclusive param: " + paramName);

        for (const auto &mutuallyExclusiveParamName : mutuallyExclusiveParams)
            fluxDimensions.emplace(mutuallyExclusiveParamName, nUniverses);
    }

    // Open the flux file for reading
    // Get the scaling factor to go from the event rate in the flux file, to the flux itself
    // Here we get the flux in neutrinos/POT/bin/cm2 by scaling the event rate (in the samples) down by POT and the cross-sectional area of the active volume
    // Here we also scale up the fluxes by 1e10 for the sake of comparison just so we are working with reasonable numbers
    const float unitsScaling = 1e10;
    const auto fluxScaleFactor = unitsScaling / (config.flux.pot * (GeometryHelper::highX - GeometryHelper::lowX) * (GeometryHelper::highY - GeometryHelper::lowY));

    // Get the flux bin edges and content
    const auto fluxHistNames = CrossSectionHelper::GetNominalFluxHistNames(config.flux.nuPdgsSignal, config.flux.nuPdgToHistName, config.flux.nomHistPattern);
    const auto &[fluxBinEdges, fluxBinValuesNominal] = CrossSectionHelper::ReadNominalFlux(config.flux.fileName, fluxHistNames, config.flux.pot);

    // Setup the flux reweightor
    CrossSectionHelper::FluxReweightor fluxReweightor(fluxBinEdges, fluxBinValuesNominal, fluxDimensions);
    std::cout << "Integrated nominal flux: " << fluxReweightor.GetIntegratedNominalFlux() << std::endl;

    std::cout << "Flux includes:" << std::endl;
    for (const auto &fluxHistName : fluxHistNames)
    {
        std::cout << "  - " << fluxHistName << std::endl;
    }

    // Open the flux file
    TFile *pFluxFile = new TFile(config.flux.fileName.c_str(), "READ");
    if (!pFluxFile->IsOpen())
        throw std::invalid_argument("PlotFluxVariations - Can't open flux file: " + config.flux.fileName);

    // Open the input file for reading and enable the branches with systematic event weights
    FileReader reader(config.files.overlaysFileName);
    reader.EnableSystematicBranches();

    // Get the bound event address
    auto pEvent = reader.GetBoundEventAddress();

    // Loop over the events
    const auto nEvents = reader.GetNumberOfEvents();
    for (unsigned int i = 0; i < nEvents; ++i)
    {
        AnalysisHelper::PrintLoadingBar(i, nEvents);
        reader.LoadEvent(i);

        // Only use the desired neutrino flavours
        if (std::find(config.flux.nuPdgsSignal.begin(), config.flux.nuPdgsSignal.end(), pEvent->truth.nuPdgCode()) == config.flux.nuPdgsSignal.end())
            continue;

        // Skip events that aren't in the active volume
        if (!AnalysisHelper::IsInActiveVolume(pEvent->truth.nuVertex()))
            continue;

        // Get the nominal event weight
        const auto weight = AnalysisHelper::GetNominalEventWeight(pEvent);

        // Get the flux weights
        const auto fluxWeights = CrossSectionHelper::GetWeightsMap(pEvent->truth, fluxDimensions);

        // Add this event to the flux reweightor
        fluxReweightor.AddEvent(pEvent->truth.nuEnergy(), weight, fluxWeights);
    }

    // Setup a table for the RMS difference between the fluxes in the files and using the reweighting method
    FormattingHelper::Table table({"Parameter", "True flux STD", "", "Frac diff RMS"});

    // Get the integrated flux in each universe
    for (const auto &[paramName, nUniverses] : fluxDimensions)
    {
        std::cout << "Parameter name: " << paramName << std::endl;

        // Get the directories correponding to this parameter
        std::unordered_map<int, TDirectoryFile *> nuPdgToDirMap;
        for (const auto &nuPdg : config.flux.nuPdgsSignal)
        {
            const auto nuNameIter = config.flux.nuPdgToHistName.find(nuPdg);
            if (nuNameIter == config.flux.nuPdgToHistName.end())
                throw std::logic_error("PlotFluxVariations - Can't find name for PDG code: " + std::to_string(nuPdg));

            const auto &nuName = nuNameIter->second;
            const auto dirName = std::regex_replace(std::regex_replace(config.flux.variationDirPattern, std::regex("NEUTRINO"), nuName), std::regex("PARAMNAME"), paramName);
            const auto pDir = static_cast<TDirectoryFile *>(pFluxFile->Get(dirName.c_str()));
            if (!pDir)
                throw std::logic_error("PlotFluxVariations - Can't get TDirectoryFile: " + dirName);

            nuPdgToDirMap.emplace(nuPdg, pDir);
        }

        // Store the true fluxes, and the fractional difference when using reweighting
        std::vector<float> trueFluxes, fractionalDiffs;

        // Keep track of the min and max fluxes seen
        float minFlux = +std::numeric_limits<float>::max();
        float maxFlux = -std::numeric_limits<float>::max();
        float maxFracDiff = -std::numeric_limits<float>::max();

        // Loop over the universes
        for (unsigned int iUni = 0; iUni < nUniverses; ++iUni)
        {
            // Get the reweighted flux
            const auto fluxReweighted = fluxReweightor.GetIntegratedFluxVariation(paramName, iUni);

            // Get the flux from the file for this universe
            float fluxUni = 0.f;
            for (const auto &nuPdg : config.flux.nuPdgsSignal)
            {
                // Get the name of the flux histogram
                const auto nuName = config.flux.nuPdgToHistName.at(nuPdg);
                const auto fluxName = std::regex_replace(std::regex_replace(std::regex_replace(config.flux.variationHistPattern, std::regex("PARAMNAME"), paramName), std::regex("UNIVERSEINDEX"), std::to_string(iUni)), std::regex("NEUTRINO"), nuName);
                const auto pFluxUni = static_cast<TH1F *>(nuPdgToDirMap.at(nuPdg)->Get(fluxName.c_str()));
                if (!pFluxUni)
                    throw std::logic_error("PlotFluxVariations - Can't get TH1F: " + fluxName);

                fluxUni += pFluxUni->Integral() * fluxScaleFactor;
            }

            float fracDiff = 0.f;
            if (std::abs(fluxUni) > std::numeric_limits<float>::epsilon())
            {
                fracDiff = (fluxReweighted - fluxUni) / fluxUni;
            }

            minFlux = std::min(minFlux, fluxUni);
            maxFlux = std::max(maxFlux, fluxUni);
            maxFracDiff = std::max(maxFracDiff, std::abs(fracDiff));

            // Store the result
            trueFluxes.push_back(fluxUni);
            fractionalDiffs.push_back(fracDiff);
        }

        // Make a histogram for the frational differences
        auto pHistTrueFlux = std::make_shared<TH1F>(("trueFlux_" + paramName).c_str(), "", 50, minFlux, maxFlux);
        auto pHistFracDiff = std::make_shared<TH1F>(("fracDiff_" + paramName).c_str(), "", 50, -maxFracDiff, maxFracDiff);
        auto pHist2D = std::make_shared<TH2F>(("fracDiff-vs-trueFlux_" + paramName).c_str(), "", 50, minFlux, maxFlux, 50, -maxFracDiff, maxFracDiff);
        for (unsigned int iUni = 0; iUni < nUniverses; ++iUni)
        {
            pHistTrueFlux->Fill(trueFluxes.at(iUni));
            pHistFracDiff->Fill(fractionalDiffs.at(iUni));
            pHist2D->Fill(trueFluxes.at(iUni), fractionalDiffs.at(iUni));
        }

        auto pCanvas = PlottingHelper::GetCanvas();
        pHistTrueFlux->Draw("hist");
        PlottingHelper::SaveCanvas(pCanvas, "plotFluxVariations_trueFlux_" + paramName);

        pHistFracDiff->Draw("hist");
        PlottingHelper::SaveCanvas(pCanvas, "plotFluxVariations_fracDiff_" + paramName);

        pHist2D->Draw("colz");
        PlottingHelper::SaveCanvas(pCanvas, "plotFluxVariations_fracDiff-vs-trueFlux_" + paramName);

        table.AddEmptyRow();
        table.SetEntry("Parameter", paramName);
        table.SetEntry("True flux STD", pHistTrueFlux->GetStdDev());
        table.SetEntry("Frac diff RMS", pHistFracDiff->GetRMS());
    }

    table.WriteToFile("plotFluxVariations_differences.md");
}

} // namespace ubcc1pi_macros
