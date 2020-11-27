/**
 *  @file  ubcc1pi_standalone/Helpers/CrossSectionHelperNew.cxx
 *
 *  @brief The implementation of the cross section helper class
 */

#include "ubcc1pi_standalone/Helpers/CrossSectionHelperNew.h"

#include <stdexcept>

namespace ubcc1pi
{

CrossSectionHelperNew::FluxReweightor::FluxReweightor(const std::vector<float> &binEdges, const std::vector<float> &binValuesNominal, const SystDimensionsMap &fluxWeightsDimensions) :
    m_binEdges(binEdges),
    m_dimensions(fluxWeightsDimensions),
    m_pFluxNominal(CrossSectionHelperNew::GetTH1F(binEdges)),
    m_pSpectrumNominal(CrossSectionHelperNew::GetTH1F(binEdges)),
    m_spectrumVariations(CrossSectionHelperNew::GetSystTH1FMap(binEdges, fluxWeightsDimensions))
{
    // There should be one more bin edge than the number of bins
    if (binEdges.size() != binValuesNominal.size() + 1)
        throw std::invalid_argument("FluxReweightor::FluxReweightor - Number of bin edges doesn't match number of bins!");

    // Fill the nominal flux histogram
    for (unsigned int iBin = 1; iBin <= binValuesNominal.size(); ++iBin)
    {
        // ATTN here we shift the index by 1 as root enumerates bins from 1 and the vector enumerates from 0
        m_pFluxNominal->SetBinContent(iBin, binValuesNominal.at(iBin - 1));
    }
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void CrossSectionHelperNew::FluxReweightor::AddEvent(const float trueNuEnergy, const float nominalWeight, const SystFloatMap &fluxWeights)
{
    // Make sure we have a valid input map
    CrossSectionHelperNew::ValidateSystMap(fluxWeights, m_dimensions);

    // Fill the nominal spectrum
    m_pSpectrumNominal->Fill(trueNuEnergy, nominalWeight);

    // Fill the universes
    for (const auto &[paramName, weights] : fluxWeights)
    {
        auto &spectrumUniverses = m_spectrumVariations.at(paramName);

        for (unsigned int iUni = 0; iUni < weights.size(); ++iUni)
        {
            const auto universeWeight = nominalWeight * weights.at(iUni);
            spectrumUniverses.at(iUni)->Fill(trueNuEnergy, universeWeight);
        }
    }
}

// -----------------------------------------------------------------------------------------------------------------------------------------

std::shared_ptr<TH1F> CrossSectionHelperNew::FluxReweightor::GetNominalFlux() const
{
    return m_pFluxNominal;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

CrossSectionHelperNew::SystTH1FMap CrossSectionHelperNew::FluxReweightor::GetFluxVariations() const
{
    // Setup an empty map to fill
    const auto nBins = m_binEdges.size() - 1;
    auto map = CrossSectionHelperNew::GetSystTH1FMap(m_binEdges, m_dimensions);

    // Loop ove all paramteres and universes
    for (const auto &[paramName, universes] : m_spectrumVariations)
    {
        auto &reweightedFluxUniverses = map.at(paramName);

        for (unsigned int iUni = 0; iUni < universes.size(); ++iUni)
        {
            const auto pSpectrumVariation = universes.at(iUni);

            // Loop over each bin
            for (unsigned iBin = 1; iBin <= nBins; ++iBin)
            {
                // Get the nominal specturem and the specturm in this universe
                const auto nEventsNominal = m_pSpectrumNominal->GetBinContent(iBin);
                const auto nEventsVariation = pSpectrumVariation->GetBinContent(iBin);

                // Get the ratio of the variation to the nominal (if not possible use a unit weight)
                const auto weight = (nEventsNominal <= std::numeric_limits<float>::epsilon() ? 1.f : (nEventsVariation / nEventsNominal) );

                // Apply this weight to the nominal flux
                const auto nominalFlux = m_pFluxNominal->GetBinContent(iBin);
                const auto reweightedFlux = nominalFlux * weight;

                // Store the result
                reweightedFluxUniverses.at(iUni)->SetBinContent(iBin, reweightedFlux);
            }
        }
    }

    return map;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

float CrossSectionHelperNew::FluxReweightor::GetIntegratedFlux(const std::shared_ptr<TH1F> &pFlux) const
{
    float total = 0.f;

    for (unsigned int iBin = 1; iBin <= static_cast<unsigned int>(pFlux->GetNbinsX()); ++iBin)
    {
        // ATTN the input flux has already been scaled by bin width
        total += pFlux->GetBinContent(iBin);
    }

    return total;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

float CrossSectionHelperNew::FluxReweightor::GetIntegratedNominalFlux() const
{
    return this->GetIntegratedFlux(m_pFluxNominal);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

CrossSectionHelperNew::SystFloatMap CrossSectionHelperNew::FluxReweightor::GetIntegratedFluxVariations() const
{
    SystFloatMap map;

    for (const auto &[paramName, universes] : this->GetFluxVariations())
    {
        std::vector<float> integratedFluxes;
        for (const auto &pFlux : universes)
        {
            integratedFluxes.push_back(this->GetIntegratedFlux(pFlux));
        }

        map.emplace(paramName, integratedFluxes);
    }

    return map;
}

// -----------------------------------------------------------------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------------------------------------------------------------

std::shared_ptr<TH1F> CrossSectionHelperNew::GetTH1F(const std::vector<float> &binEdges)
{
    return std::make_shared<TH1F>(("xSecHist_" + std::to_string(m_histCount++)).c_str(), "", binEdges.size() - 1, binEdges.data());
}

// -----------------------------------------------------------------------------------------------------------------------------------------

CrossSectionHelperNew::SystTH1FMap CrossSectionHelperNew::GetSystTH1FMap(const std::vector<float> &binEdges, const SystDimensionsMap &dimensions)
{
    SystTH1FMap map;

    for (const auto &[paramName, nUniverses] : dimensions)
    {
        std::vector< std::shared_ptr<TH1F> > universes;
        for (unsigned int iUni = 0; iUni < nUniverses; ++iUni)
        {
            universes.push_back(CrossSectionHelperNew::GetTH1F(binEdges));
        }

        map.emplace(paramName, universes);
    }

    return map;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

CrossSectionHelperNew::SystFloatMap CrossSectionHelperNew::GetWeightsMap(const Event::Truth &truth, const SystDimensionsMap &dimensions)
{
    if (!truth.systParamNames.IsSet() || !truth.systParamFirstValueIndex.IsSet() || !truth.systParamValues.IsSet())
        throw std::invalid_argument("CrossSectionHelperNew::GetWeightsMap - Systematic parameters are not set in the input truth information");

    const auto &systParamNames = truth.systParamNames();
    const auto &systParamFirstValueIndex = truth.systParamFirstValueIndex();
    const auto &systParamValues = truth.systParamValues();

    // Get the total number of parameters available
    const auto nParameters = systParamNames.size();

    // Make a new map for the output
    SystFloatMap map;

    // Loop over all parameters in the input dimensions map
    for (const auto &[paramName, nUniverses] : dimensions)
    {
        // Find the desired name
        const auto iter = std::find(systParamNames.begin(), systParamNames.end(), paramName);
        if (iter == systParamNames.end())
            throw std::invalid_argument("CrossSectionHelperNew::GetWeightsMap - Unknown parameter: " + paramName);

        // Get the index of the requested parameter
        const unsigned int index = std::distance(systParamNames.begin(), iter);

        // Get the first and last value in the weights vector
        const unsigned int firstValueIndex = systParamFirstValueIndex.at(index);
        const unsigned int lastValueIndex = ((index == nParameters - 1u) ? systParamValues.size() : systParamFirstValueIndex.at(index + 1u));

        // Pick out the weights corresponding the the desired parameter
        const std::vector<float> weights(std::next(systParamValues.begin(), firstValueIndex), std::next(systParamValues.begin(), lastValueIndex));

        // Check that we have the right number of weights
        if (weights.size() != nUniverses)
            throw std::invalid_argument("CrossSectionHelperNew::GetWeightsMap - Number of weights for parameter " + paramName + " doesn't match the input dimensions");

        // Setup the output map entry
        auto &outputWeights = map[paramName];

        // Fill the output weights
        for (const auto &weight : weights)
        {
            // ATTN sometimes we have non-physical weights in the input. Here we decide what to do with them!

            // If the weight is non-negative and non-infinite then it's okay!
            if (weight >= 0.f && std::abs(std::abs(weight) - std::numeric_limits<float>::max()) > std::numeric_limits<float>::epsilon() && std::isfinite(weight))
            {
                outputWeights.push_back(weight);
                continue;
            }

            // Force negative weights to be zero and infinite weights to be unity
            outputWeights.push_back((weight < 0.f) ? 0.f : 1.f);
        }
    }

    return map;
}

} // namespace ubcc1pi
