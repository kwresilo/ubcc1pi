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
    m_dimensions(fluxWeightsDimensions),
    m_pFluxNominal(CrossSectionHelperNew::GetTH1F(binEdges)),
    m_pSpectrumNominal(CrossSectionHelperNew::GetTH1F(binEdges)),
    m_spectrumVariations(CrossSectionHelperNew::GetSystTH1FMap(binEdges, fluxWeightsDimensions)),

    // Define the cache-function that gets the integrated flux variations
    // ATTN it might look a bit odd to some to define a function in the constructor! To be clear, here we are initalizing the member
    // variable named m_getIntegratedFluxVariationFunction. This is a SystCacheFunction that takes a lambda function as it's first paramter
    // and stores the logic for later. When the user calls GetIntegratedFluxVariation() we ask the SystCacheFunction to check if it has the
    // value we are looking for in it's cache. If the answer is yes, then we just return the cached value. If the answer is no, then we call
    // the logic (defined in the lambda function below), store the result in the cache for later and return it.
    m_getIntegratedFluxVariationFunction([&](const std::string &paramName, const unsigned int universeIndex)
    {
        // Setup a new histogram to hold the reweighted flux
        auto pReweightedFlux = CrossSectionHelperNew::GetTH1F(binEdges);

        // Get the event rate spectrum in this universe
        const auto pSpectrumVariation = m_spectrumVariations.at(paramName).at(universeIndex);

        // Loop over each bin
        const auto nBins = binEdges.size() - 1;
        for (unsigned iBin = 1; iBin <= nBins; ++iBin)
        {
            // Get bin value for the nominal spectrum and the specturm in this universe
            const auto nEventsNominal = m_pSpectrumNominal->GetBinContent(iBin);
            const auto nEventsVariation = pSpectrumVariation->GetBinContent(iBin);

            // Get the ratio of the variation to the nominal (if not possible use a unit weight)
            const auto weight = (nEventsNominal <= std::numeric_limits<float>::epsilon() ? 1.f : (nEventsVariation / nEventsNominal) );

            // Apply this weight to the nominal flux
            const auto nominalFlux = m_pFluxNominal->GetBinContent(iBin);
            const auto reweightedFlux = nominalFlux * weight;

            // Store the result
            pReweightedFlux->SetBinContent(iBin, reweightedFlux);
        }

        // Return the integrated flux
        return this->GetIntegratedFlux(pReweightedFlux);

    }, m_dimensions)
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

    // Here we clear any cached values for the integrated flux variations because we are modifying the event rate spectra
    m_getIntegratedFluxVariationFunction.ClearCache();

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

float CrossSectionHelperNew::FluxReweightor::GetIntegratedFluxVariation(const std::string &paramName, const unsigned int universeIndex)
{
    // Call the cache-function to either calculate the return value or retrieve it from the cache
    return m_getIntegratedFluxVariationFunction(paramName, universeIndex);
}

// -----------------------------------------------------------------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------------------------------------------------------------

CrossSectionHelperNew::CrossSection::CrossSection(const SystParams &systParams, const std::vector<float> &binEdges, const bool hasUnderflow, const bool hasOverflow) :
    m_systParams(systParams),
    m_binEdges(binEdges),
    m_metadata(binEdges.size() - 1, hasUnderflow, hasOverflow),
    m_pSignal_true_nom(CrossSectionHelperNew::GetTH1F(binEdges)),
    m_signal_true_fluxes(CrossSectionHelperNew::GetSystTH1FMap(binEdges, systParams.fluxDimensions)),
    m_signal_true_xsecs(CrossSectionHelperNew::GetSystTH1FMap(binEdges, systParams.xsecDimensions)),
    m_signal_true_detVars(CrossSectionHelperNew::GetSystUnisimTH1FMap(binEdges, systParams.detVarDimensions)),
    m_signal_true_detCVs(CrossSectionHelperNew::GetCVSystUnisimTH1FMap(binEdges, systParams.detVarDimensions)),
    m_pSignal_selected_recoTrue_nom(CrossSectionHelperNew::GetTH2F(binEdges)),
    m_signal_selected_recoTrue_fluxes(CrossSectionHelperNew::GetSystTH2FMap(binEdges, systParams.fluxDimensions)),
    m_signal_selected_recoTrue_xsecs(CrossSectionHelperNew::GetSystTH2FMap(binEdges, systParams.xsecDimensions)),
    m_signal_selected_recoTrue_detVars(CrossSectionHelperNew::GetSystUnisimTH2FMap(binEdges, systParams.detVarDimensions)),
    m_signal_selected_recoTrue_detCVs(CrossSectionHelperNew::GetCVSystUnisimTH2FMap(binEdges, systParams.detVarDimensions)),
    m_pBackground_selected_reco_nom(CrossSectionHelperNew::GetTH1F(binEdges)),
    m_background_selected_reco_fluxes(CrossSectionHelperNew::GetSystTH1FMap(binEdges, systParams.fluxDimensions)),
    m_background_selected_reco_xsecs(CrossSectionHelperNew::GetSystTH1FMap(binEdges, systParams.xsecDimensions)),
    m_background_selected_reco_detVars(CrossSectionHelperNew::GetSystUnisimTH1FMap(binEdges, systParams.detVarDimensions)),
    m_background_selected_reco_detCVs(CrossSectionHelperNew::GetCVSystUnisimTH1FMap(binEdges, systParams.detVarDimensions)),
    m_pBNBData_selected_reco(CrossSectionHelperNew::GetTH1F(binEdges))
{
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void CrossSectionHelperNew::CrossSection::AddSignalEvent(const float trueValue, const float recoValue, const bool isSelected, const float nominalWeight, const SystFloatMap &fluxWeights, const SystFloatMap &xsecWeights)
{
    // TODO
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void CrossSectionHelperNew::CrossSection::AddSignalEventDetVar(const float trueValue, const float recoValue, const bool isSelected, const float nominalWeight, const std::string &paramName)
{
    // TODO
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void CrossSectionHelperNew::CrossSection::AddSignalEventDetVarCV(const float trueValue, const float recoValue, const bool isSelected, const float nominalWeight, const std::string &name)
{
    // TODO
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void CrossSectionHelperNew::CrossSection::AddSelectedBackgroundEvent(const float recoValue, const float nominalWeight, const SystFloatMap &fluxWeights, const SystFloatMap &xsecWeights)
{
    // TODO
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void CrossSectionHelperNew::CrossSection::AddSelectedBackgroundEventDetVar(const float recoValue, const float nominalWeight, const std::string &paramName)
{
    // TODO
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void CrossSectionHelperNew::CrossSection::AddSelectedBackgroundEventDetVarCV(const float recoValue, const float nominalWeight, const std::string &name)
{
    // TODO
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void CrossSectionHelperNew::CrossSection::AddSelectedBNBDataEvent(const float recoValue)
{
    // TODO
}

// -----------------------------------------------------------------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------------------------------------------------------------

std::shared_ptr<TH1F> CrossSectionHelperNew::GetTH1F(const std::vector<float> &binEdges)
{
    return std::make_shared<TH1F>(("xSecHist_" + std::to_string(m_histCount++)).c_str(), "", binEdges.size() - 1, binEdges.data());
}

// -----------------------------------------------------------------------------------------------------------------------------------------

std::shared_ptr<TH2F> CrossSectionHelperNew::GetTH2F(const std::vector<float> &binEdges)
{
    return std::make_shared<TH2F>(("xSecHist_" + std::to_string(m_histCount++)).c_str(), "", binEdges.size() - 1, binEdges.data(), binEdges.size() - 1, binEdges.data());
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

CrossSectionHelperNew::SystUnisimTH1FMap CrossSectionHelperNew::GetSystUnisimTH1FMap(const std::vector<float> &binEdges, const SystUnisimDimensionsMap &dimensions)
{
    SystUnisimTH1FMap map;

    for (const auto &[paramName, cvName] : dimensions)
    {
        map.emplace(paramName, CrossSectionHelperNew::GetTH1F(binEdges));
        (void) cvName; // Not required
    }

    return map;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

CrossSectionHelperNew::SystUnisimTH1FMap CrossSectionHelperNew::GetCVSystUnisimTH1FMap(const std::vector<float> &binEdges, const SystUnisimDimensionsMap &dimensions)
{
    SystUnisimTH1FMap map;

    for (const auto &[paramName, cvName] : dimensions)
    {
        if (map.find(cvName) != map.end())
            continue;

        map.emplace(cvName, CrossSectionHelperNew::GetTH1F(binEdges));
        (void) paramName; // Not required
    }

    return map;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

CrossSectionHelperNew::SystTH2FMap CrossSectionHelperNew::GetSystTH2FMap(const std::vector<float> &binEdges, const SystDimensionsMap &dimensions)
{
    SystTH2FMap map;

    for (const auto &[paramName, nUniverses] : dimensions)
    {
        std::vector< std::shared_ptr<TH2F> > universes;
        for (unsigned int iUni = 0; iUni < nUniverses; ++iUni)
        {
            universes.push_back(CrossSectionHelperNew::GetTH2F(binEdges));
        }

        map.emplace(paramName, universes);
    }

    return map;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

CrossSectionHelperNew::SystUnisimTH2FMap CrossSectionHelperNew::GetSystUnisimTH2FMap(const std::vector<float> &binEdges, const SystUnisimDimensionsMap &dimensions)
{
    SystUnisimTH2FMap map;

    for (const auto &[paramName, cvName] : dimensions)
    {
        map.emplace(paramName, CrossSectionHelperNew::GetTH2F(binEdges));
        (void) cvName; // Not required
    }

    return map;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

CrossSectionHelperNew::SystUnisimTH2FMap CrossSectionHelperNew::GetCVSystUnisimTH2FMap(const std::vector<float> &binEdges, const SystUnisimDimensionsMap &dimensions)
{
    SystUnisimTH2FMap map;

    for (const auto &[paramName, cvName] : dimensions)
    {
        if (map.find(cvName) != map.end())
            continue;

        map.emplace(cvName, CrossSectionHelperNew::GetTH2F(binEdges));
        (void) paramName; // Not required
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

// -----------------------------------------------------------------------------------------------------------------------------------------

CrossSectionHelperNew::SystFloatMap CrossSectionHelperNew::GetWeightsMap(const Event::Truth &truth, const SystDimensionsMap &dimensions, const SystMutuallyExclusiveDimensionsMap &mutuallyExclusiveDimensions)
{
    // Make a new map for the output
    SystFloatMap map;

    // Make new dimensions map to populate with the parameters that are available in the input truth information directly
    SystDimensionsMap availableParamDimensions;

    const auto &systParamNames = truth.systParamNames();

    // Find which parameters (if any) are built from mutually exclusive parameters
    for (const auto &[paramName, nUniverses] : dimensions)
    {
        // Check if this parameter name is available in the input truth information, if so just skip it for later
        if (std::find(systParamNames.begin(), systParamNames.end(), paramName) != systParamNames.end())
        {
            availableParamDimensions.emplace(paramName, nUniverses);
            continue;
        }

        // Insist that parameters not available in the input truth information are listed in the mutually exlusive dimensions
        if (mutuallyExclusiveDimensions.find(paramName) == mutuallyExclusiveDimensions.end())
            throw std::invalid_argument("CrossSectionHelperNew::GetWeightsMap - Parameter \"" + paramName + "\" isn't listed in the input truth information or in the supplied mutually exclusive parameter names");

        // We have a mutually exclusive parameter, so get it's weights
        const auto &[parameters, nUniversesCheck] = mutuallyExclusiveDimensions.at(paramName);
        if (nUniverses != nUniversesCheck)
            throw std::invalid_argument("CrossSectionHelperNew::GetWeightsMap - Inconsistent number of universes specified for parameter: \"" + paramName + "\"");

        const auto weights = CrossSectionHelperNew::GetMutuallyExclusiveWeights(truth, parameters, nUniverses);

        // Add the result to the output map
        map.emplace(paramName, weights);
    }

    // Now get the weights from the rest of the parameters
    auto availableParamWeights = CrossSectionHelperNew::GetWeightsMap(truth, availableParamDimensions);

    // Combine this together and return it!
    map.insert(availableParamWeights.begin(), availableParamWeights.end());
    return map;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

std::vector<float> CrossSectionHelperNew::GetMutuallyExclusiveWeights(const Event::Truth &truth, const std::vector<std::string> &parameters, const unsigned int nUniverses)
{
    // Build a dimensions map using the same number of universes for each parameter
    SystDimensionsMap dimensions;
    for (const auto &paramName : parameters)
        dimensions.emplace(paramName, nUniverses);

    // Get the weights map for these parameters individually
    const auto weightsMap = CrossSectionHelperNew::GetWeightsMap(truth, dimensions);

    // Extract the mutually exclusive weights
    std::vector<float> weights;
    for (unsigned int iUni = 0; iUni < nUniverses; ++iUni)
    {
        // Find the parameter name for which the weight is not exactly one
        float weight = 1.f;
        bool foundNonUnitWeight = false;

        for (const auto &paramName : parameters)
        {
            const auto universeWeight = weightsMap.at(paramName).at(iUni);
            if (std::abs(universeWeight - 1.f) > std::numeric_limits<float>::epsilon())
            {
                // This weight is something other than one!

                // If we find more than one weight that's not exactly one - then the parameters aren't mutually exclusive!
                if (foundNonUnitWeight)
                    throw std::logic_error("CrossSectionHelperNew::GetMutuallyExclusiveWeights - Input weights aren't mutually exclusive");

                // Store this weight
                weight = universeWeight;
                foundNonUnitWeight = true;
            }
        }

        // ATTN it's still posible that all universe weights are exactly one. This is okay - just use 1.f as the weight
        weights.push_back(weight);
    }

    return weights;
}

} // namespace ubcc1pi
