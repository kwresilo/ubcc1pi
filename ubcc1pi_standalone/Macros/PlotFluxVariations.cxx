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
    TFile *pFluxFile = new TFile(config.flux.fileName.c_str(), "READ");
    if (!pFluxFile->IsOpen())
        throw std::invalid_argument("PlotFluxVariations - Can't open flux file: " + config.flux.fileName);

    // Get the nominal flux
    const auto pFluxHist = static_cast<TH1F *>(pFluxFile->Get(config.flux.nomHistName.c_str()));

    // Get the scaling factor to go from the event rate in the flux file, to the flux itself
    // Here we get the flux in neutrinos/POT/bin/cm2 by scaling the event rate (in the samples) down by POT and the cross-sectional area of the active volume
    // Here we also scale up the fluxes by 1e10 for the sake of comparison just so we are working with reasonable numbers
    const float unitsScaling = 1e10;
    const auto fluxScaleFactor = unitsScaling / (config.flux.pot * (GeometryHelper::highX - GeometryHelper::lowX) * (GeometryHelper::highY - GeometryHelper::lowY));

    // Get the flux bin edges and content
    std::vector<float> fluxBinEdges, fluxBinValuesNominal;
    const unsigned int nFluxBins = pFluxHist->GetNbinsX();
    for (unsigned int iBin = 1; iBin <= nFluxBins; ++iBin)
    {
        const auto lowEdge = pFluxHist->GetBinLowEdge(iBin);
        fluxBinEdges.push_back(lowEdge);

        // Add the upper edge of the last bin
        if (iBin == nFluxBins)
        {
            const auto binWidth = pFluxHist->GetBinWidth(iBin);
            fluxBinEdges.push_back(lowEdge + binWidth);
        }

        const auto nNeutrinos = pFluxHist->GetBinContent(iBin);
        fluxBinValuesNominal.push_back(nNeutrinos * fluxScaleFactor);
    }

    // Setup the flux reweightor
    CrossSectionHelper::FluxReweightor fluxReweightor(fluxBinEdges, fluxBinValuesNominal, fluxDimensions);
    std::cout << "Integrated nominal flux: " << fluxReweightor.GetIntegratedNominalFlux() << std::endl;

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

        // Only use numu events
        if (pEvent->truth.nuPdgCode() != 14)
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

        // Get the directory correponding to this parameter
        const auto dirName = std::regex_replace(config.flux.variationDirPattern, std::regex("PARAMNAME"), paramName);
        const auto pDir = static_cast<TDirectoryFile *>(pFluxFile->Get(dirName.c_str()));

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

            // Get the flux from the input file
            const auto fluxName = std::regex_replace(std::regex_replace(config.flux.variationHistPattern, std::regex("PARAMNAME"), paramName), std::regex("UNIVERSEINDEX"), std::to_string(iUni));

            const auto pFluxUni = static_cast<TH1F *>(pDir->Get(fluxName.c_str()));
            const float fluxUni = pFluxUni->Integral() * fluxScaleFactor;

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

            //// BEGIN DEBUG
            // Print the result
            std::cout << iUni << ": " << fluxUni << ", " << fluxReweighted << ", " << fracDiff << std::endl;
            //// END DEBUG
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
