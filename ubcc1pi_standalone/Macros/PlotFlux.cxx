/**
 *  @file  ubcc1pi_standalone/Macros/PlotFlux.cxx
 *
 *  @brief The implementation file of the PlotFlux macro
 */

#include "ubcc1pi_standalone/Macros/Macros.h"
#include "ubcc1pi_standalone/Objects/FileReader.h"

#include "ubcc1pi_standalone/Helpers/CrossSectionHelper.h"
#include "ubcc1pi_standalone/Helpers/PlottingHelper.h"

using namespace ubcc1pi;

namespace ubcc1pi_macros
{

void PlotFlux(const Config &config)
{
    // TODO make configurable
    const float maxNuEnergy = 6.8f;

    // Setup a map to hold the flux histograms
    std::map< int, std::shared_ptr<TH1F> > histMap;

    // The maximum & minimum flux bins
    float maxFlux = -std::numeric_limits<float>::max();
    float minFlux = +std::numeric_limits<float>::max();

    // Loop over the neutrino PDG codes
    for (const auto &pdg : {+14, -14, +12, -12})
    {
        // Get the flux for this neutrino flavour
        const auto fluxHistNames = CrossSectionHelper::GetNominalFluxHistNames({ pdg }, config.flux.nuPdgToHistName, config.flux.nomHistPattern);
        const auto &[fluxBinEdges, fluxValues] = CrossSectionHelper::ReadNominalFlux(config.flux.fileName, fluxHistNames, config.flux.pot);

        // Setup a histogram for this flux
        const auto name = fluxHistNames.front() + "_ubcc1pi";
        const auto nBins = fluxBinEdges.size() - 1;
        auto pHist = std::make_shared<TH1F>(name.c_str(), "", nBins, fluxBinEdges.data());

        // Fill the histogram
        for (unsigned int iBin = 1; iBin <= nBins; ++iBin)
        {
            const auto binContent = fluxValues.at(iBin - 1);
            pHist->SetBinContent(iBin, binContent);

            // Get the minimum and maximum values (provided this bin isn't zero)
            if (binContent >= std::numeric_limits<float>::epsilon() && fluxBinEdges.at(iBin) <= maxNuEnergy)
            {
                maxFlux = std::max(maxFlux, binContent);
                minFlux = std::min(minFlux, binContent);
            }
        }

        // Add the hisogram to the map
        histMap.emplace(pdg, pHist);
    }

    // Set the line colours
    PlottingHelper::SetLineStyle(histMap.at(+14), PlottingHelper::Primary);
    PlottingHelper::SetLineStyle(histMap.at(-14), PlottingHelper::Primary);
    PlottingHelper::SetLineStyle(histMap.at(+12), PlottingHelper::Secondary);
    PlottingHelper::SetLineStyle(histMap.at(-12), PlottingHelper::Secondary);

    // Used dashed lines for the antineutrinos
    histMap.at(-14)->SetLineStyle(2);
    histMap.at(-12)->SetLineStyle(2);

    // Add a bit of padding to the top of the plot
    maxFlux += std::exp(std::log(maxFlux - minFlux) * 0.05);

    // Make the plot
    auto pCanvas = PlottingHelper::GetCanvas();
    pCanvas->SetLogy();
    bool isFirst = true;
    for (const auto &[pdg, pHist] : histMap)
    {
        pHist->GetXaxis()->SetRangeUser(0.f, maxNuEnergy);
        pHist->GetYaxis()->SetRangeUser(minFlux, maxFlux);

        if (isFirst)
        {
            pHist->Draw("hist");
        }
        else
        {
            pHist->Draw("hist same");
        }

        isFirst = false;
    }

    PlottingHelper::SaveCanvas(pCanvas, "plotFlux_allFluxes");
}

} // namespace ubcc1pi_macros
