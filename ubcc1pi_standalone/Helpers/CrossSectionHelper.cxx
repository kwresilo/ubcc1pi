/**
 *  @file  ubcc1pi_standalone/Helpers/CrossSectionHelper.cxx
 *
 *  @brief The implementation of the cross section helper class
 */

#include "ubcc1pi_standalone/Helpers/CrossSectionHelper.h"

#include "ubcc1pi_standalone/Helpers/AnalysisHelper.h"
#include "ubcc1pi_standalone/Helpers/GeometryHelper.h"

#include <TFile.h>

#include <stdexcept>
#include <regex>

namespace ubcc1pi
{

CrossSectionHelper::FluxReweightor::FluxReweightor(const std::vector<float> &binEdges, const std::vector<float> &binValuesNominal, const SystDimensionsMap &fluxWeightsDimensions) :
    /// @cond Doxygen can't handle this initilizer list
    m_dimensions(fluxWeightsDimensions),
    m_pFluxNominal(CrossSectionHelper::GetTH1F(binEdges)),
    m_pSpectrumNominal(CrossSectionHelper::GetTH1F(binEdges)),
    m_spectrumVariations(CrossSectionHelper::GetSystTH1FMap(binEdges, fluxWeightsDimensions)),

    // Define the cache-function that gets the integrated flux variations
    // ATTN it might look a bit odd to some to define a function in the constructor! To be clear, here we are initalizing the member
    // variable named m_getIntegratedFluxVariationFunction. This is a SystCacheFunction that takes a lambda function as it's first paramter
    // and stores the logic for later. When the user calls GetIntegratedFluxVariation() we ask the SystCacheFunction to check if it has the
    // value we are looking for in it's cache. If the answer is yes, then we just return the cached value. If the answer is no, then we call
    // the logic (defined in the lambda function below), store the result in the cache for later and return it.
    m_getIntegratedFluxVariationFunction([&](const std::string &paramName, const unsigned int universeIndex)
    {
        // Setup a new histogram to hold the reweighted flux
        auto pReweightedFlux = CrossSectionHelper::GetTH1F(binEdges);

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
    /// @endcond
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

void CrossSectionHelper::FluxReweightor::AddEvent(const float trueNuEnergy, const float nominalWeight, const SystFloatMap &fluxWeights)
{
    // Make sure we have a valid input map
    CrossSectionHelper::ValidateSystMap(fluxWeights, m_dimensions);

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

std::shared_ptr<TH1F> CrossSectionHelper::FluxReweightor::GetNominalFlux() const
{
    return m_pFluxNominal;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

float CrossSectionHelper::FluxReweightor::GetIntegratedFlux(const std::shared_ptr<TH1F> &pFlux) const
{
    float total = 0.f;

    for (unsigned int iBin = 1; iBin <= static_cast<unsigned int>(pFlux->GetNbinsX()); ++iBin)
    {
        total += pFlux->GetBinContent(iBin);
    }

    return total;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

float CrossSectionHelper::FluxReweightor::GetIntegratedNominalFlux() const
{
    return this->GetIntegratedFlux(m_pFluxNominal);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

float CrossSectionHelper::FluxReweightor::GetIntegratedFluxVariation(const std::string &paramName, const unsigned int universeIndex)
{
    // Call the cache-function to either calculate the return value or retrieve it from the cache
    return m_getIntegratedFluxVariationFunction(paramName, universeIndex);
}

// -----------------------------------------------------------------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------------------------------------------------------------

CrossSectionHelper::CrossSection::CrossSection(const SystParams &systParams, const std::vector<float> &binEdges, const bool hasUnderflow, const bool hasOverflow, const bool scaleByBinWidth) :
    m_systParams(systParams),
    m_binEdges(binEdges),
    m_metadata(binEdges, hasUnderflow, hasOverflow, scaleByBinWidth),
    m_scaleByBinWidth(scaleByBinWidth),
    m_pSignal_true_nom(CrossSectionHelper::GetTH1F(binEdges)),
    m_pSignal_selected_recoTrue_nom(CrossSectionHelper::GetTH2F(binEdges)),
    m_pBackground_selected_reco_nom(CrossSectionHelper::GetTH1F(binEdges)),
    m_pBNBData_selected_reco(CrossSectionHelper::GetTH1F(binEdges))
{
    // Make sure the bin edges are valid
    const int nAnalysisBins = (binEdges.size() - 1) - (hasUnderflow ? 1 : 0) - (hasOverflow ? 1 : 0);
    if (nAnalysisBins <= 0)
        throw std::invalid_argument("CrossSection::CrossSection - Insufficient bin edges provided!");

    // Setup a dimensions map for the miscellaneous ("misc") multisim parameters.
    // These are added in ad-hoc (i.e. they don't come from the weights stored in the truth information)
    const SystDimensionsMap miscDimensions = {

        // Here we apply a weight to each event using the Poisson bootstrap method to estimate the MC stats uncertainty
        { "bootstrap", systParams.nBootstrapUniverses },

        // Here we apply a 2-universe multisim in which we weight the dirt up and down by 100%
        { "dirt", 2 }
    };

    // Setup the maps to hold the histograms for the multisim parameters
    m_signal_true_multisims.emplace("flux", CrossSectionHelper::GetSystTH1FMap(binEdges, systParams.fluxDimensions));
    m_signal_true_multisims.emplace("xsec", CrossSectionHelper::GetSystTH1FMap(binEdges, systParams.xsecDimensions));
    m_signal_true_multisims.emplace("misc", CrossSectionHelper::GetSystTH1FMap(binEdges, miscDimensions));

    m_signal_selected_recoTrue_multisims.emplace("flux", CrossSectionHelper::GetSystTH2FMap(binEdges, systParams.fluxDimensions));
    m_signal_selected_recoTrue_multisims.emplace("xsec", CrossSectionHelper::GetSystTH2FMap(binEdges, systParams.xsecDimensions));
    m_signal_selected_recoTrue_multisims.emplace("misc", CrossSectionHelper::GetSystTH2FMap(binEdges, miscDimensions));

    m_background_selected_reco_multisims.emplace("flux", CrossSectionHelper::GetSystTH1FMap(binEdges, systParams.fluxDimensions));
    m_background_selected_reco_multisims.emplace("xsec", CrossSectionHelper::GetSystTH1FMap(binEdges, systParams.xsecDimensions));
    m_background_selected_reco_multisims.emplace("misc", CrossSectionHelper::GetSystTH1FMap(binEdges, miscDimensions));

    // Setup the maps to hold the histograms for the unisim parameters
    m_signal_true_unisims.emplace("detector", CrossSectionHelper::GetSystUnisimTH1FMap(binEdges, systParams.detVarDimensions));
    m_signal_selected_recoTrue_unisims.emplace("detector", CrossSectionHelper::GetSystUnisimTH2FMap(binEdges, systParams.detVarDimensions));
    m_background_selected_reco_unisims.emplace("detector", CrossSectionHelper::GetSystUnisimTH1FMap(binEdges, systParams.detVarDimensions));

    // ATTN we could in theory include other unisim parameters here using a different name. At the moment we only use "detector".
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void CrossSectionHelper::CrossSection::AddSignalEvent(const float recoValue, const float trueValue, const bool isSelected, const float nominalWeight, const SystFloatMap &fluxWeights, const SystFloatMap &xsecWeights)
{
    // Get the miscellaneous weights
    const SystFloatMap miscWeights = {

        // Generate the bootstrap weights
        { "bootstrap", CrossSectionHelper::GenerateBootstrapWeights(m_systParams.nBootstrapUniverses) },

        // A signal event is never dirt - so just use a unit weight
        { "dirt", {1.f, 1.f} }
    };

    m_pSignal_true_nom->Fill(trueValue, nominalWeight);
    CrossSectionHelper::FillSystTH1FMap(trueValue, nominalWeight, fluxWeights, m_signal_true_multisims.at("flux"));
    CrossSectionHelper::FillSystTH1FMap(trueValue, nominalWeight, xsecWeights, m_signal_true_multisims.at("xsec"));
    CrossSectionHelper::FillSystTH1FMap(trueValue, nominalWeight, miscWeights, m_signal_true_multisims.at("misc"));

    if (!isSelected)
        return;

    // ATTN the reco value is only used if the event is selected
    m_pSignal_selected_recoTrue_nom->Fill(recoValue, trueValue, nominalWeight);
    CrossSectionHelper::FillSystTH2FMap(recoValue, trueValue, nominalWeight, fluxWeights, m_signal_selected_recoTrue_multisims.at("flux"));
    CrossSectionHelper::FillSystTH2FMap(recoValue, trueValue, nominalWeight, xsecWeights, m_signal_selected_recoTrue_multisims.at("xsec"));
    CrossSectionHelper::FillSystTH2FMap(recoValue, trueValue, nominalWeight, miscWeights, m_signal_selected_recoTrue_multisims.at("misc"));
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void CrossSectionHelper::CrossSection::AddSignalEventDetVar(const float recoValue, const float trueValue, const bool isSelected, const float nominalWeight, const std::string &paramName)
{
    if (!static_cast<bool>(m_signal_true_unisims.at("detector").count(paramName)))
        throw std::invalid_argument("CrossSection::AddSignalEventDetVar - Unknown parameter or CV sample name \"" + paramName + "\"");

    m_signal_true_unisims.at("detector").at(paramName)->Fill(trueValue, nominalWeight);

    if (!isSelected)
        return;

    m_signal_selected_recoTrue_unisims.at("detector").at(paramName)->Fill(recoValue, trueValue, nominalWeight);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void CrossSectionHelper::CrossSection::AddSelectedBackgroundEvent(const float recoValue, const bool isDirt, const float nominalWeight, const SystFloatMap &fluxWeights, const SystFloatMap &xsecWeights)
{
    // The weight to apply as a result of dirt
    //   - For non-dirt events, just use zero (i.e. no change from the nominal simulation)
    //   - For dirt event we here use +-100%
    const auto dirtWeightDelta = isDirt ? 1.f : 0.f;

    // Get the miscellaneous weights
    const SystFloatMap miscWeights = {

        // Generate the bootstrap weights
        { "bootstrap", CrossSectionHelper::GenerateBootstrapWeights(m_systParams.nBootstrapUniverses) },

        // Apply the dirt weight
        { "dirt", {1.f - dirtWeightDelta, 1.f + dirtWeightDelta} }
    };

    m_pBackground_selected_reco_nom->Fill(recoValue, nominalWeight);
    CrossSectionHelper::FillSystTH1FMap(recoValue, nominalWeight, fluxWeights, m_background_selected_reco_multisims.at("flux"));
    CrossSectionHelper::FillSystTH1FMap(recoValue, nominalWeight, xsecWeights, m_background_selected_reco_multisims.at("xsec"));
    CrossSectionHelper::FillSystTH1FMap(recoValue, nominalWeight, miscWeights, m_background_selected_reco_multisims.at("misc"));
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void CrossSectionHelper::CrossSection::AddSelectedBackgroundEventDetVar(const float recoValue, const float nominalWeight, const std::string &paramName)
{
    if (!static_cast<bool>(m_background_selected_reco_unisims.at("detector").count(paramName)))
        throw std::invalid_argument("CrossSection::AddSelectedBackgroundEventDetVar - Unknown parameter or CV sample name \"" + paramName + "\"");

    m_background_selected_reco_unisims.at("detector").at(paramName)->Fill(recoValue, nominalWeight);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void CrossSectionHelper::CrossSection::AddSelectedBNBDataEvent(const float recoValue)
{
    m_pBNBData_selected_reco->Fill(recoValue);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

ubsmear::UBXSecMeta CrossSectionHelper::CrossSection::GetMetadata() const
{
    return m_metadata;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

std::vector<float> CrossSectionHelper::CrossSection::GetBinEdges() const
{
    return m_binEdges;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

ubsmear::UBMatrix CrossSectionHelper::CrossSection::GetBinWidths() const
{
    return m_metadata.GetBinWidths();
}

// -----------------------------------------------------------------------------------------------------------------------------------------

ubsmear::UBMatrix CrossSectionHelper::CrossSection::GetCrossSection(const ubsmear::UBMatrix &selected, const ubsmear::UBMatrix &backgrounds, const float integratedFlux, const float exposurePOT, const float nTargets) const
{
    // Check the input event rate has the right dimensions
    if (selected.GetColumns() != 1)
        throw std::invalid_argument("CrossSection::GetCrossSection - Input event rate isn't a column vector");

    if (selected.GetRows() != m_metadata.GetNBins())
        throw std::invalid_argument("CrossSection::GetCrossSection - Input event rate has the wrong number of bins");

    // Check the input backgrounds have the right dimensions
    if (backgrounds.GetColumns() != 1)
        throw std::invalid_argument("CrossSection::GetCrossSection - Input background rate isn't a column vector");

    if (backgrounds.GetRows() != m_metadata.GetNBins())
        throw std::invalid_argument("CrossSection::GetCrossSection - Input background rate has the wrong number of bins");

    // Get the normalisation factor
    const auto norm = integratedFlux * exposurePOT * nTargets;
    if (norm <= std::numeric_limits<float>::epsilon())
        throw std::invalid_argument("CrossSection::GetCrossSection - Product of flux, exposure and targets non-positive");

    // Get the column vector of bin widths
    const auto binWidths = this->GetBinWidths();

    // Now get the flux-integrated, forward-folded cross-section by
    //   - Subtracting the backgrounds
    //   - Scaling by the normalisation factor
    //   - Scaling (element wise) by the bin-widths
    return ubsmear::ElementWiseOperation(
        (selected - backgrounds),
        binWidths * norm,
        [](const float numerator, const float denominator) {

            if (denominator <= std::numeric_limits<float>::epsilon())
                throw std::logic_error("CrossSection::GetCrossSection - Found a bin in which the denominator (flux * POT * nTargets * binWidth) is <= 0");

            return numerator / denominator;
        }
    );
}

// -----------------------------------------------------------------------------------------------------------------------------------------

std::shared_ptr<ubsmear::UBMatrix> CrossSectionHelper::CrossSection::GetSmearingMatrix(const std::shared_ptr<TH1F> &pSignal_true, const std::shared_ptr<TH2F> &pSignal_selected_recoTrue) const
{
    // Get the input histograms as a matrix
    const auto allTrue = CrossSectionHelper::GetMatrixFromHist(pSignal_true);
    const auto selectedRecoTrue = CrossSectionHelper::GetMatrixFromHist(pSignal_selected_recoTrue);

    // Check theh matrices have the expected sizes
    if (!ubsmear::UBMatrixHelper::IsSquare(selectedRecoTrue))
        throw std::invalid_argument("CrossSection::GetSmearingMatrix - input reco-true matrix isn't square");

    const auto nBins = allTrue.GetRows();
    if (selectedRecoTrue.GetRows() != nBins)
        throw std::invalid_argument("CrossSection::GetSmearingMatrix - the input objects differ in number of bins");

    // Get the elements of the smearing matrix
    std::vector<float> elements;
    for (unsigned int iReco = 0; iReco < nBins; ++iReco)
    {
        for (unsigned int iTrue = 0; iTrue < nBins; ++iTrue)
        {
            // Get the total number of signal events in this true bin
            const auto denominator = allTrue.At(iTrue, 0);

            // ATTN. If we find a bin in which there are no signal events, then there's no sensible value for the smearing matrix element
            // Here we return null to indicate that the smearing matrix can't be found
            if (denominator <= std::numeric_limits<float>::epsilon())
                return std::shared_ptr<ubsmear::UBMatrix>(nullptr);

            elements.push_back(selectedRecoTrue.At(iReco, iTrue) / denominator);
        }
    }

    return std::make_shared<ubsmear::UBMatrix>(elements, nBins, nBins);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

ubsmear::UBMatrix CrossSectionHelper::CrossSection::GetSignalSelectedTrue(const std::shared_ptr<TH2F> &pSignal_selected_recoTrue) const
{
    // Get the input histogram as a matrix
    const auto recoTrueMatrix = CrossSectionHelper::GetMatrixFromHist(pSignal_selected_recoTrue);

    // Integrate over the reco bins
    std::vector<float> elements;
    for (unsigned int iTrue = 0; iTrue < recoTrueMatrix.GetColumns(); ++iTrue)
    {
        float sum = 0.f;
        for (unsigned int iReco = 0; iReco < recoTrueMatrix.GetRows(); ++iReco)
        {
            sum += recoTrueMatrix.At(iReco, iTrue);
        }

        elements.push_back(sum);
    }

    // Return the result as a column vector
    return ubsmear::UBMatrix(elements, recoTrueMatrix.GetColumns(), 1);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

ubsmear::UBMatrix CrossSectionHelper::CrossSection::GetSelectedBNBDataEvents() const
{
    return CrossSectionHelper::GetMatrixFromHist(m_pBNBData_selected_reco);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

ubsmear::UBMatrix CrossSectionHelper::CrossSection::GetSelectedBackgroundEvents() const
{
    return CrossSectionHelper::GetMatrixFromHist(m_pBackground_selected_reco_nom);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

ubsmear::UBMatrix CrossSectionHelper::CrossSection::GetSelectedSignalEvents() const
{
    return this->GetSignalSelectedTrue(m_pSignal_selected_recoTrue_nom);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

ubsmear::UBMatrix CrossSectionHelper::CrossSection::GetSignalEvents() const
{
    return CrossSectionHelper::GetMatrixFromHist(m_pSignal_true_nom);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

ubsmear::UBMatrix CrossSectionHelper::CrossSection::GetBNBDataCrossSection(const ScalingData &scalingData) const
{
    // Get the number of selected events in BNB data
    const auto selected = CrossSectionHelper::GetMatrixFromHist(m_pBNBData_selected_reco);

    // Get the number of predicted backgrounds in the nominal universe
    const auto backgrounds = CrossSectionHelper::GetMatrixFromHist(m_pBackground_selected_reco_nom);

    // Get the integrated flux in the nominal universe
    const auto integratedFlux = scalingData.pFluxReweightor->GetIntegratedNominalFlux();

    // Get the cross-section
    return this->GetCrossSection(selected, backgrounds, integratedFlux, scalingData.exposurePOT, scalingData.nTargets);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

ubsmear::UBMatrix CrossSectionHelper::CrossSection::GetPredictedCrossSection(const ScalingData &scalingData) const
{
    // Get the number of signal events in the nominal simulation in true bins
    const auto signal = CrossSectionHelper::GetMatrixFromHist(m_pSignal_true_nom);

    // Get the number of backgrounds is zero as we are effectively applying a "perfect" selection
    const auto zeroVector = ubsmear::UBMatrixHelper::GetZeroMatrix(signal.GetRows(), 1);

    // Get the integrated flux in the nominal universe
    const auto integratedFlux = scalingData.pFluxReweightor->GetIntegratedNominalFlux();

    // Get the cross-section
    return this->GetCrossSection(signal, zeroVector, integratedFlux, scalingData.exposurePOT, scalingData.nTargets);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

std::shared_ptr<ubsmear::UBMatrix> CrossSectionHelper::CrossSection::GetSmearingMatrixNominal() const
{
    // Get the smearing matrix in the nominal universe
    return this->GetSmearingMatrix(m_pSignal_true_nom, m_pSignal_selected_recoTrue_nom);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

ubsmear::UBMatrix CrossSectionHelper::CrossSection::GetSmearingMatrix() const
{
    // Get the smearing matrix in the nominal universe
    const auto pSmearingMatrix = this->GetSmearingMatrixNominal();

    // Check we were able to find the smearing matrix
    if (!pSmearingMatrix)
        throw std::logic_error("CrossSection::GetSmearingMatrix - There's a bin with no signal events! Can't get the smearing matrix");

    return *pSmearingMatrix;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

ubsmear::UBMatrix CrossSectionHelper::CrossSection::GetBNBDataCrossSectionStatUncertainty(const ScalingData &scalingData) const
{
    // Get the number of selected events in BNB data
    const auto selected = CrossSectionHelper::GetMatrixFromHist(m_pBNBData_selected_reco);

    // Get a vector containing the count uncertainties (Poisson) on each bin
    std::vector<float> elements;
    const auto nBins = selected.GetRows();
    for (unsigned int iBin = 0; iBin < nBins; ++iBin)
    {
        elements.push_back(AnalysisHelper::GetCountUncertainty(selected.At(iBin, 0)));
    }

    const ubsmear::UBMatrix selectedUncertainty(elements, nBins, 1);

    // Get the integrated flux in the nominal universe
    const auto integratedFlux = scalingData.pFluxReweightor->GetIntegratedNominalFlux();

    // Scale the uncertainty on the number of selected events by same factor as the cross-section.
    // ATTN to apply the scaling we reuse the GetCrossSection function but intead pass it the count uncertainties as "selected" a zero
    // vector as the "backgrounds"
    const auto zeroVector = ubsmear::UBMatrixHelper::GetZeroMatrix(nBins, 1);
    return this->GetCrossSection(selectedUncertainty, zeroVector, integratedFlux, scalingData.exposurePOT, scalingData.nTargets);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

CrossSectionHelper::SystBiasCovarianceMap CrossSectionHelper::CrossSection::GetBNBDataCrossSectionSystUncertainties(const ScalingData &scalingData) const
{
    // Setup the output map
    SystBiasCovarianceMap outputMap;

    // Add the multisim parameters
    // ATTN the use of m_signal_true_multisims is arbitrary, we could use any multisims map
    for (const auto &[group, map] : m_signal_true_multisims)
    {
        for (const auto &entry : map)
        {
            const auto &paramName = entry.first;
            outputMap[group].emplace(paramName, this->GetBNBDataCrossSectionDistributionParams(group, paramName, scalingData));
        }
    }

    // Add the unisim parameters
    // ATTN again the use of m_signal_true_unisims is arbitrary, we could use any unisims map
    for (const auto &[group, map] : m_signal_true_unisims)
    {
        // ATTN for now we only have detector parameters that are unisims, if others are added later then we would need to add a clause here
        // to get the relevent dimensions.
        if (group != "detector")
            throw std::logic_error("CrossSection::GetBNBDataCrossSectionSystUncertainties - Don't know dimensions of group: \"" + group + "\"");

        const auto dimensions = m_systParams.detVarDimensions;

        // Extract the central-value names from the dimensions object
        std::vector<std::string> cvNames;
        for (const auto &[paramName, cvName] : dimensions)
            cvNames.push_back(cvName);

        for (const auto &entry : map)
        {
            const auto &paramName = entry.first;

            // ATTN The paramName could either be a detector variation sample of the name of a CV sample.
            // Here we check if the paramName is for a CV sample and if so, we can skip it.
            const auto isCVName = (std::find(cvNames.begin(), cvNames.end(), paramName) != cvNames.end());
            if (isCVName)
                continue;

            // Get the central-value name that corresponds to this parameter
            const auto &cvName = dimensions.at(paramName);

            outputMap[group].emplace(paramName, this->GetBNBDataCrossSectionDistributionParamsUnisim(group, paramName, cvName, scalingData));
        }
    }

    // Get the special POT normalisation uncertainty
    const auto pXSecNom = std::make_shared<ubsmear::UBMatrix>( this->GetBNBDataCrossSection(scalingData) );
    outputMap["misc"].emplace("POT", this->GetDistributionParamsNormalisation(pXSecNom, m_systParams.potFracUncertainty));

    return outputMap;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

CrossSectionHelper::SystBiasCovarianceMap CrossSectionHelper::CrossSection::GetSmearingMatrixSystUncertainties() const
{
    // Setup the output map
    SystBiasCovarianceMap outputMap;

    // Add the multisim parameters
    // ATTN the use of m_signal_true_multisims is arbitrary, we could use any multisims map
    for (const auto &[group, map] : m_signal_true_multisims)
    {
        for (const auto &entry : map)
        {
            const auto &paramName = entry.first;
            outputMap[group].emplace(paramName, this->GetSmearingMatrixDistributionParams(group, paramName));
        }
    }

    // Add the unisim parameters
    // ATTN again the use of m_signal_true_unisims is arbitrary, we could use any unisims map
    for (const auto &[group, map] : m_signal_true_unisims)
    {
        // ATTN for now we only have detector parameters that are unisims, if others are added later then we would need to add a clause here
        // to get the relevent dimensions.
        if (group != "detector")
            throw std::logic_error("CrossSection::GetSmearingMatrixSystUncertainties - Don't know dimensions of group: \"" + group + "\"");

        const auto dimensions = m_systParams.detVarDimensions;

        // Extract the central-value names from the dimensions object
        std::vector<std::string> cvNames;
        for (const auto &[paramName, cvName] : dimensions)
            cvNames.push_back(cvName);

        for (const auto &entry : map)
        {
            const auto &paramName = entry.first;

            // ATTN The paramName could either be a detector variation sample of the name of a CV sample.
            // Here we check if the paramName is for a CV sample and if so, we can skip it.
            const auto isCVName = (std::find(cvNames.begin(), cvNames.end(), paramName) != cvNames.end());
            if (isCVName)
                continue;

            // Get the central-value name that corresponds to this parameter
            const auto &cvName = dimensions.at(paramName);

            outputMap[group].emplace(paramName, this->GetSmearingMatrixDistributionParamsUnisim(group, paramName, cvName));
        }
    }

    // Get the special POT normalisation uncertainty
    // ATTN no uncertainty is assigned to the smearing matrix due to POT normalisation (use a factor of 0.f)
    const auto pSmearingNom = CrossSectionHelper::FlattenMatrix(this->GetSmearingMatrixNominal());
    outputMap["misc"].emplace("POT", this->GetDistributionParamsNormalisation(pSmearingNom, 0.f));

    return outputMap;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

CrossSectionHelper::SystBiasCovariancePair CrossSectionHelper::CrossSection::GetPredictedCrossSectionStatUncertainty(const ScalingData &scalingData) const
{
    // Get the number of universes for this systematic parameter
    // ATTN this choice of m_signal_true_multisims here is arbitrary, any of the mutlisims maps would do
    const auto nUniverses = m_systParams.nBootstrapUniverses;

    // Get the nominal cross-section
    const auto pXSecNom = std::make_shared<ubsmear::UBMatrix>(this->GetPredictedCrossSection(scalingData));

    // Get the number of backgrounds is zero as we are effectively applying a "perfect" selection
    const auto zeroVector = ubsmear::UBMatrixHelper::GetZeroMatrix(pXSecNom->GetRows(), 1);

    // Get the integrated flux in the nominal universe
    const auto integratedFlux = scalingData.pFluxReweightor->GetIntegratedNominalFlux();

    // The bootstrap universes
    const auto universes = m_signal_true_multisims.at("misc").at("bootstrap");

    // Get the distribution parameters
    return this->GetDistributionParams([&](const unsigned int &iUni)
    {
        // This is the function we call to get the value in each universe

        // Get the number of signal events in the nominal simulation in true bins
        const auto signal = CrossSectionHelper::GetMatrixFromHist(universes.at(iUni));

        // Get the cross-section
        return std::make_shared<ubsmear::UBMatrix> ( this->GetCrossSection(signal, zeroVector, integratedFlux, scalingData.exposurePOT, scalingData.nTargets) );

    }, nUniverses, pXSecNom);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

ubsmear::UBMatrix CrossSectionHelper::CrossSection::GetBNBDataCrossSectionInUniverse(const std::string &group, const std::string &paramName, const unsigned int universeIndex, const ScalingData &scalingData) const
{
    // ATTN for speed, this function (and other similar functions) doesn't explicitly check if the input group, paramName, universeIndex is
    // valid. Instead we rely on the internal range checking of the .at() method for STL containers. This is for peformance reasons as this
    // function gets called once per universe!

    // Get the number of selected events in BNB data
    const auto selected = CrossSectionHelper::GetMatrixFromHist(m_pBNBData_selected_reco);

    // Get the number of predicted backgrounds in the supplied universe
    const auto backgrounds = CrossSectionHelper::GetMatrixFromHist(m_background_selected_reco_multisims.at(group).at(paramName).at(universeIndex));

    // Get the integrated flux in the supplied universe (if it's not a flux parameter, then use the nominal universe)
    const auto integratedFlux = (
        group == "flux"
            ? scalingData.pFluxReweightor->GetIntegratedFluxVariation(paramName, universeIndex)
            : scalingData.pFluxReweightor->GetIntegratedNominalFlux()
    );

    // Get the cross-section
    return this->GetCrossSection(selected, backgrounds, integratedFlux, scalingData.exposurePOT, scalingData.nTargets);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

ubsmear::UBMatrix CrossSectionHelper::CrossSection::GetBNBDataCrossSectionForUnisim(const std::string &group, const std::string &paramName, const ScalingData &scalingData) const
{
    // Get the number of selected events in BNB data
    const auto selected = CrossSectionHelper::GetMatrixFromHist(m_pBNBData_selected_reco);

    // Get the number of predicted backgrounds in the supplied universe
    const auto backgrounds = CrossSectionHelper::GetMatrixFromHist(m_background_selected_reco_unisims.at(group).at(paramName));

    // Get the integrated flux in the nominal universe
    const auto integratedFlux = scalingData.pFluxReweightor->GetIntegratedNominalFlux();

    // Get the cross-section
    return this->GetCrossSection(selected, backgrounds, integratedFlux, scalingData.exposurePOT, scalingData.nTargets);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

CrossSectionHelper::SystBiasCovariancePair CrossSectionHelper::CrossSection::GetDistributionParams(const std::function<std::shared_ptr<ubsmear::UBMatrix>(const unsigned int)> &func, const unsigned int nUniverses, const std::shared_ptr<ubsmear::UBMatrix> &pNominal) const
{
    if (nUniverses == 0)
        throw std::logic_error("CrossSection::GetDisributionParams - No universes supplied");

    // Insist the nominal input is a column vector
    if (pNominal->GetColumns() != 1)
        throw std::logic_error("CrossSection::GetDisributionParams - Input nominal is not a column vector");

    // Get the number of bins
    const auto nBins = pNominal->GetRows();

    // Setup some empty matrices to hold the mean vector and error matrices
    auto meanSum = ubsmear::UBMatrixHelper::GetZeroMatrix(nBins, 1);
    auto errorMatrixSum = ubsmear::UBMatrixHelper::GetZeroMatrix(nBins, nBins);

    // Loop over the universes
    unsigned int nValidUniverses = 0u;
    for (unsigned int iUni = 0; iUni < nUniverses; ++iUni)
    {
        //// BEGIN DEBUG
        AnalysisHelper::PrintLoadingBar(iUni, nUniverses);
        //// END DEBUG

        // Get the value in this universe
        const auto pUniverse = func(iUni);

        // Insist that we get a valid quanitiy in this universes
        if (!pUniverse)
            continue;

        // Count the number of valid universes
        nValidUniverses++;

        // Loop over the bins
        // ATTN for speed we don't expcitly check here if the universe has the same number of bins as the nominal. Instead we rely on the
        // internal range checking of the ubsmear::UBMatrix class
        for (unsigned int iBin = 0; iBin < nBins; ++iBin)
        {
            // Add up the universes
            meanSum.SetElement(iBin, 0, meanSum.At(iBin, 0) + pUniverse->At(iBin, 0));

            // Loop over the bins again
            // ATTN the error matrix is symmetic so here we only set one half of the off-diagonals within the universe loop and then copy
            // them over to the other half afterwards. This is done for performance reasons as it halves the number of inserts required
            const auto diffI = pUniverse->At(iBin, 0) - pNominal->At(iBin, 0);
            for (unsigned int jBin = 0; jBin <= iBin; ++jBin)
            {
                const auto diffJ = pUniverse->At(jBin, 0) - pNominal->At(jBin, 0);

                // Add up the error matrix elements
                errorMatrixSum.SetElement(iBin, jBin, errorMatrixSum.At(iBin, jBin) + diffI*diffJ);
            }
        }
    }

    std::cout << "DEBUG - Valid universes: " << nValidUniverses << " / " << nUniverses << std::endl;

    // Scale the sums by the number of universes
    if (nValidUniverses == 0)
        throw std::logic_error("CrossSection::GetDisributionParams - Desired quantity was invalid in all universes");

    const auto scaleFactor = 1.f / static_cast<float>(nValidUniverses);
    const auto mean = meanSum * scaleFactor;
    const auto errorMatrix = errorMatrixSum * scaleFactor;

    // Get the bias vector
    const auto pBias = std::make_shared<ubsmear::UBMatrix>(mean - *pNominal);

    // Get the covariance matrix
    auto pCovarianceMatrix = std::make_shared<ubsmear::UBMatrix>(ubsmear::UBMatrixHelper::GetZeroMatrix(nBins, nBins));
    for (unsigned int iBin = 0; iBin < nBins; ++iBin)
    {
        for (unsigned int jBin = 0; jBin <= iBin; ++jBin)
        {
            // ATTN here we use the fact that V_ij = E_ij - b_i*b_j
            //   - E_ij = error matrix element
            //   - V_ij = covariance matrix element
            //   - b_i = bias vector element
            //
            // This is done instead of calculating the covariance directly as this only requires one pass through the universe loop.
            const auto covariance = errorMatrix.At(iBin, jBin) - pBias->At(iBin, 0)*pBias->At(jBin, 0);

            pCovarianceMatrix->SetElement(iBin, jBin, covariance);
            pCovarianceMatrix->SetElement(jBin, iBin, covariance);
        }
    }

    return { pBias, pCovarianceMatrix };
}

// -----------------------------------------------------------------------------------------------------------------------------------------

CrossSectionHelper::SystBiasCovariancePair CrossSectionHelper::CrossSection::GetDistributionParamsUnisim(const std::shared_ptr<ubsmear::UBMatrix> &pVaried, const std::shared_ptr<ubsmear::UBMatrix> &pCentralValue, const std::shared_ptr<ubsmear::UBMatrix> &pNominal) const
{
    // Check the input matrices are valid
    if (!pVaried || !pCentralValue || !pNominal)
        throw std::invalid_argument("CrossSection::GetDistributionParamsUnisim - One or more of the inputs are not valid");

    // Check ths input matrices have the desired dimensions
    if (pVaried->GetColumns() != 1 || pCentralValue->GetColumns() != 1 || pNominal->GetColumns() != 1)
        throw std::invalid_argument("CrossSection::GetDistributionParamsUnisim - One or more of the inputs are not a column vector");

    const auto nBins = pVaried->GetRows();
    if (pCentralValue->GetRows() != nBins || pNominal->GetRows() != nBins)
        throw std::invalid_argument("CrossSection::GetDistributionParamsUnisim - Input column vectors have different sizes");

    // Get the fractional bias of the cross-section away from the central value
    const auto fracBias = ubsmear::ElementWiseOperation(
        (*pVaried - *pCentralValue), *pCentralValue,
        [](const float numerator, const float denominator)
        {
            // ATTN here we decide what to do if there's a bin of the central-value that's zero
            // In this case we don't have any information, so here we set the fractional bias to 0.
            if (std::abs(denominator) <= std::numeric_limits<float>::epsilon())
                return 0.f;

            return numerator / denominator;
        }
    );

    // Apply the fractional bias to the nominal cross-section
    const auto bias = ubsmear::ElementWiseOperation(
        fracBias, *pNominal,
        [](const float lhs, const float rhs) { return lhs * rhs; }
    );

    // Return the result (use a zero matrix for the covariance)
    return {
        std::make_shared<ubsmear::UBMatrix>(bias),
        std::make_shared<ubsmear::UBMatrix>(ubsmear::UBMatrixHelper::GetZeroMatrix(nBins, nBins))
    };
}

// -----------------------------------------------------------------------------------------------------------------------------------------

CrossSectionHelper::SystBiasCovariancePair CrossSectionHelper::CrossSection::GetDistributionParamsNormalisation(const std::shared_ptr<ubsmear::UBMatrix> &pNominal, const float fracUncertainty) const
{
    // Check the input matrix are valid
    if (!pNominal)
        throw std::invalid_argument("CrossSection::GetDistributionParamsNormalisation - The input nominal is not valid");

    // Check ths input matrix has the desired dimensions
    if (pNominal->GetColumns() != 1)
        throw std::invalid_argument("CrossSection::GetDistributionParamsNormalisation - The input is not a column vector");

    const auto nBins = pNominal->GetRows();

    // Scale the nominal value by the fractional uncertainty to get the bias, and use a zero matrix for the covariance
    return {
        std::make_shared<ubsmear::UBMatrix>((*pNominal) * fracUncertainty),
        std::make_shared<ubsmear::UBMatrix>(ubsmear::UBMatrixHelper::GetZeroMatrix(nBins, nBins))
    };
}

// -----------------------------------------------------------------------------------------------------------------------------------------

CrossSectionHelper::SystBiasCovariancePair CrossSectionHelper::CrossSection::GetBNBDataCrossSectionDistributionParams(const std::string &group, const std::string &paramName, const ScalingData &scalingData) const
{
    std::cout << "DEBUG - Processing parameter: " << group << ", " << paramName << std::endl;

    // Get the number of universes for this systematic parameter
    // ATTN this choice of m_signal_true_multisims here is arbitrary, any of the mutlisims maps would do
    const auto nUniverses = m_signal_true_multisims.at(group).at(paramName).size();

    // Get the nominal cross-section
    const auto pXSecNom = std::make_shared<ubsmear::UBMatrix>(this->GetBNBDataCrossSection(scalingData));

    // Get the distirbution parameters
    return this->GetDistributionParams([&](const unsigned int &iUni)
    {
        // This is the function we call to get the value in each universe
        return std::make_shared<ubsmear::UBMatrix>( this->GetBNBDataCrossSectionInUniverse(group, paramName, iUni, scalingData) );

    }, nUniverses, pXSecNom);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

CrossSectionHelper::SystBiasCovariancePair CrossSectionHelper::CrossSection::GetBNBDataCrossSectionDistributionParamsUnisim(const std::string &group, const std::string &paramName, const std::string &cvName, const ScalingData &scalingData) const
{
    std::cout << "DEBUG - Processing parameter: " << group << ", " << paramName << " (with central value: " << cvName << ")" << std::endl;

    // Get the cross-section using the unisim sample, in the central value sample and in the nominal universe
    const auto pXSec = std::make_shared<ubsmear::UBMatrix>( this->GetBNBDataCrossSectionForUnisim(group, paramName, scalingData) );
    const auto pXSecCV = std::make_shared<ubsmear::UBMatrix>( this->GetBNBDataCrossSectionForUnisim(group, cvName, scalingData) );
    const auto pXSecNom = std::make_shared<ubsmear::UBMatrix>( this->GetBNBDataCrossSection(scalingData) );

    // Get the parameters
    return this->GetDistributionParamsUnisim(pXSec, pXSecCV, pXSecNom);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

std::shared_ptr<ubsmear::UBMatrix> CrossSectionHelper::CrossSection::GetSmearingMatrixInUniverse(const std::string &group, const std::string &paramName, const unsigned int universeIndex) const
{
    return this->GetSmearingMatrix(m_signal_true_multisims.at(group).at(paramName).at(universeIndex), m_signal_selected_recoTrue_multisims.at(group).at(paramName).at(universeIndex));
}

// -----------------------------------------------------------------------------------------------------------------------------------------

std::shared_ptr<ubsmear::UBMatrix> CrossSectionHelper::CrossSection::GetSmearingMatrixForUnisim(const std::string &group, const std::string &paramName) const
{
    return this->GetSmearingMatrix(m_signal_true_unisims.at(group).at(paramName), m_signal_selected_recoTrue_unisims.at(group).at(paramName));
}

// -----------------------------------------------------------------------------------------------------------------------------------------

CrossSectionHelper::SystBiasCovariancePair CrossSectionHelper::CrossSection::GetSmearingMatrixDistributionParams(const std::string &group, const std::string &paramName) const
{
    std::cout << "DEBUG - Smearing matrix. Processing parameter: " << group << ", " << paramName << std::endl;

    // Get the number of universes for this systematic parameter
    // ATTN this choice of m_signal_true_multisims here is arbitrary, any of the mutlisims maps would do
    const auto nUniverses = m_signal_true_multisims.at(group).at(paramName).size();

    // Get the nominal smearing matrix
    const auto pSmearingNom = CrossSectionHelper::FlattenMatrix(this->GetSmearingMatrixNominal());
    if (!pSmearingNom)
        throw std::logic_error("CrossSection::GetSmearingMatrixDistributionParams - The nominal smearing matrix can't be calculated");

    // Get the distirbution parameters
    return this->GetDistributionParams([&](const unsigned int &iUni)
    {
        // This is the function we call to get the value in each universe
        return CrossSectionHelper::FlattenMatrix(this->GetSmearingMatrixInUniverse(group, paramName, iUni));

    }, nUniverses, pSmearingNom);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

CrossSectionHelper::SystBiasCovariancePair CrossSectionHelper::CrossSection::GetSmearingMatrixDistributionParamsUnisim(const std::string &group, const std::string &paramName, const std::string &cvName) const
{
    std::cout << "DEBUG - Smearing matrix. Processing parameter: " << group << ", " << paramName << " (with central value: " << cvName << ")" << std::endl;

    // Get the smearing matrix using the unisim sample, in the central value sample and in the nominal universe
    const auto pSmearing = CrossSectionHelper::FlattenMatrix(this->GetSmearingMatrixForUnisim(group, paramName));
    const auto pSmearingCV = CrossSectionHelper::FlattenMatrix(this->GetSmearingMatrixForUnisim(group, cvName));
    const auto pSmearingNom = CrossSectionHelper::FlattenMatrix(this->GetSmearingMatrixNominal());

    // Get the parameters
    return this->GetDistributionParamsUnisim(pSmearing, pSmearingCV, pSmearingNom);
}

// -----------------------------------------------------------------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------------------------------------------------------------

std::shared_ptr<TH1F> CrossSectionHelper::GetTH1F(const std::vector<float> &binEdges)
{
    return std::make_shared<TH1F>(("xSecHist_" + std::to_string(m_histCount++)).c_str(), "", binEdges.size() - 1, binEdges.data());
}

// -----------------------------------------------------------------------------------------------------------------------------------------

std::shared_ptr<TH2F> CrossSectionHelper::GetTH2F(const std::vector<float> &binEdges)
{
    return std::make_shared<TH2F>(("xSecHist_" + std::to_string(m_histCount++)).c_str(), "", binEdges.size() - 1, binEdges.data(), binEdges.size() - 1, binEdges.data());
}

// -----------------------------------------------------------------------------------------------------------------------------------------

CrossSectionHelper::SystTH1FMap CrossSectionHelper::GetSystTH1FMap(const std::vector<float> &binEdges, const SystDimensionsMap &dimensions)
{
    SystTH1FMap map;

    for (const auto &[paramName, nUniverses] : dimensions)
    {
        std::vector< std::shared_ptr<TH1F> > universes;
        for (unsigned int iUni = 0; iUni < nUniverses; ++iUni)
        {
            universes.push_back(CrossSectionHelper::GetTH1F(binEdges));
        }

        map.emplace(paramName, universes);
    }

    return map;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

CrossSectionHelper::SystUnisimTH1FMap CrossSectionHelper::GetSystUnisimTH1FMap(const std::vector<float> &binEdges, const SystUnisimDimensionsMap &dimensions)
{
    SystUnisimTH1FMap map;

    // Setup a vector to keep track of the used central-values sample names
    std::vector<std::string> allCVNames;

    for (const auto &[paramName, cvName] : dimensions)
    {
        // Include all parameter names
        map.emplace(paramName, CrossSectionHelper::GetTH1F(binEdges));

        // If we haven't seen this central value sample name before, then add it too!
        if (std::find(allCVNames.begin(), allCVNames.end(), cvName) == allCVNames.end())
        {
            if (map.find(cvName) != map.end())
                throw std::invalid_argument("CrossSectionHelper::GetSystUnisimTH1FMap - Name \"" + cvName + "\" is describes a parameter and a central-value sample");

            map.emplace(cvName, CrossSectionHelper::GetTH1F(binEdges));
            allCVNames.push_back(cvName);
        }
    }

    return map;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

CrossSectionHelper::SystTH2FMap CrossSectionHelper::GetSystTH2FMap(const std::vector<float> &binEdges, const SystDimensionsMap &dimensions)
{
    SystTH2FMap map;

    for (const auto &[paramName, nUniverses] : dimensions)
    {
        std::vector< std::shared_ptr<TH2F> > universes;
        for (unsigned int iUni = 0; iUni < nUniverses; ++iUni)
        {
            universes.push_back(CrossSectionHelper::GetTH2F(binEdges));
        }

        map.emplace(paramName, universes);
    }

    return map;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

CrossSectionHelper::SystUnisimTH2FMap CrossSectionHelper::GetSystUnisimTH2FMap(const std::vector<float> &binEdges, const SystUnisimDimensionsMap &dimensions)
{
    SystUnisimTH2FMap map;

    // Setup a vector to keep track of the used central-values sample names
    std::vector<std::string> allCVNames;

    for (const auto &[paramName, cvName] : dimensions)
    {
        // Include all parameter names
        map.emplace(paramName, CrossSectionHelper::GetTH2F(binEdges));

        // If we haven't seen this central value sample name before, then add it too!
        if (std::find(allCVNames.begin(), allCVNames.end(), cvName) == allCVNames.end())
        {
            if (map.find(cvName) != map.end())
                throw std::invalid_argument("CrossSectionHelper::GetSystUnisimTH2FMap - Name \"" + cvName + "\" is describes a parameter and a central-value sample");

            map.emplace(cvName, CrossSectionHelper::GetTH2F(binEdges));
            allCVNames.push_back(cvName);
        }
    }

    return map;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

CrossSectionHelper::SystFloatMap CrossSectionHelper::GetUnitWeightsMap(const SystDimensionsMap &dimensions)
{
    // Make a new map for the output
    SystFloatMap map;

    // Loop over all parameters in the input dimensions map
    for (const auto &[paramName, nUniverses] : dimensions)
    {
        // Add a weight of 1.f for each universe
        map.emplace(paramName, std::vector<float>(nUniverses, 1.f));
    }

    return map;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

CrossSectionHelper::SystFloatMap CrossSectionHelper::ScaleWeightsMap(const SystFloatMap &weightsMap, const float &divisor)
{
    // Check the input divisor isn't zero
    if (std::abs(divisor) <= std::numeric_limits<float>::epsilon())
        throw std::invalid_argument("CrossSectionHelper::ScaleWeightsMap - The input divisor is zero");

    auto map = weightsMap;
    for (auto &[paramName, weights] : map)
    {
        for (auto &weight : weights)
        {
            weight /= divisor;
        }
    }

    return map;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

CrossSectionHelper::SystFloatMap CrossSectionHelper::GetWeightsMap(const Event::Truth &truth, const SystDimensionsMap &dimensions)
{
    if (!truth.systParamNames.IsSet() || !truth.systParamFirstValueIndex.IsSet() || !truth.systParamValues.IsSet())
        throw std::invalid_argument("CrossSectionHelper::GetWeightsMap - Systematic parameters are not set in the input truth information");

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
            throw std::invalid_argument("CrossSectionHelper::GetWeightsMap - Unknown parameter: " + paramName);

        // Get the index of the requested parameter
        const unsigned int index = std::distance(systParamNames.begin(), iter);

        // Get the first and last value in the weights vector
        const unsigned int firstValueIndex = systParamFirstValueIndex.at(index);
        const unsigned int lastValueIndex = ((index == nParameters - 1u) ? systParamValues.size() : systParamFirstValueIndex.at(index + 1u));

        // Pick out the weights corresponding the the desired parameter
        const std::vector<float> weights(std::next(systParamValues.begin(), firstValueIndex), std::next(systParamValues.begin(), lastValueIndex));

        // Check that we have the right number of weights
        if (weights.size() != nUniverses)
            throw std::invalid_argument("CrossSectionHelper::GetWeightsMap - Number of weights for parameter " + paramName + " doesn't match the input dimensions");

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

CrossSectionHelper::SystFloatMap CrossSectionHelper::GetWeightsMap(const Event::Truth &truth, const SystDimensionsMap &dimensions, const SystMutuallyExclusiveDimensionsMap &mutuallyExclusiveDimensions)
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
            throw std::invalid_argument("CrossSectionHelper::GetWeightsMap - Parameter \"" + paramName + "\" isn't listed in the input truth information or in the supplied mutually exclusive parameter names");

        // We have a mutually exclusive parameter, so get it's weights
        const auto &[parameters, nUniversesCheck] = mutuallyExclusiveDimensions.at(paramName);
        if (nUniverses != nUniversesCheck)
            throw std::invalid_argument("CrossSectionHelper::GetWeightsMap - Inconsistent number of universes specified for parameter: \"" + paramName + "\"");

        const auto weights = CrossSectionHelper::GetMutuallyExclusiveWeights(truth, parameters, nUniverses);

        // Add the result to the output map
        map.emplace(paramName, weights);
    }

    // Now get the weights from the rest of the parameters
    auto availableParamWeights = CrossSectionHelper::GetWeightsMap(truth, availableParamDimensions);

    // Combine this together and return it!
    map.insert(availableParamWeights.begin(), availableParamWeights.end());
    return map;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

std::vector<float> CrossSectionHelper::GetMutuallyExclusiveWeights(const Event::Truth &truth, const std::vector<std::string> &parameters, const unsigned int nUniverses)
{
    // Build a dimensions map using the same number of universes for each parameter
    SystDimensionsMap dimensions;
    for (const auto &paramName : parameters)
        dimensions.emplace(paramName, nUniverses);

    // Get the weights map for these parameters individually
    const auto weightsMap = CrossSectionHelper::GetWeightsMap(truth, dimensions);

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
                    throw std::logic_error("CrossSectionHelper::GetMutuallyExclusiveWeights - Input weights aren't mutually exclusive");

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

// -----------------------------------------------------------------------------------------------------------------------------------------

std::vector<float> CrossSectionHelper::GenerateBootstrapWeights(const unsigned int nUniverses)
{
    std::poisson_distribution<int> poisson(1.f);

    std::vector<float> weights;
    for (unsigned int iUni = 0; iUni < nUniverses; ++iUni)
    {
        weights.push_back(static_cast<float>(poisson(m_generator)));
    }

    return weights;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void CrossSectionHelper::FillSystTH1FMap(const float value, const float nominalWeight, const SystFloatMap &weights, SystTH1FMap &map)
{
    // Validate the dimensions of the input SystMaps
    CrossSectionHelper::ValidateSystMap(weights, CrossSectionHelper::GetSystMapDimensions(map));

    // Fill the histograms in each systematic universe
    for (const auto &[paramName, universeWeights] : weights)
    {
        const auto nUniverses = universeWeights.size();
        auto &histVector = map.at(paramName);

        // ATTN this should never happen because of the above ValidateSystMap check, but including here for sanity
        if (histVector.size() != nUniverses)
            throw std::logic_error("CrossSectionHelper::FillSystTH1FMap - Input weights map and SystTH1FMap have different dimensions");

        for (unsigned int iUni = 0; iUni < nUniverses; ++iUni)
        {
            // ATTN we use the product of the nominal weight and the universe weight
            const auto weight = nominalWeight * universeWeights.at(iUni);

            histVector.at(iUni)->Fill(value, weight);
        }
    }
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void CrossSectionHelper::FillSystTH2FMap(const float xValue, const float yValue, const float nominalWeight, const SystFloatMap &weights, SystTH2FMap &map)
{
    // Validate the dimensions of the input SystMaps
    CrossSectionHelper::ValidateSystMap(weights, CrossSectionHelper::GetSystMapDimensions(map));

    // Fill the histograms in each systematic universe
    for (const auto &[paramName, universeWeights] : weights)
    {
        const auto nUniverses = universeWeights.size();
        auto &histVector = map.at(paramName);

        // ATTN this should never happen because of the above ValidateSystMap check, but including here for sanity
        if (histVector.size() != nUniverses)
            throw std::logic_error("CrossSectionHelper::FillSystTH1FMap - Input weights map and SystTH1FMap have different dimensions");

        for (unsigned int iUni = 0; iUni < nUniverses; ++iUni)
        {
            // ATTN we use the product of the nominal weight and the universe weight
            const auto weight = nominalWeight * universeWeights.at(iUni);

            histVector.at(iUni)->Fill(xValue, yValue, weight);
        }
    }
}

// -----------------------------------------------------------------------------------------------------------------------------------------

ubsmear::UBMatrix CrossSectionHelper::GetMatrixFromHist(const std::shared_ptr<TH1F> &pHist)
{
    std::vector<float> elements;

    const unsigned int nBins = pHist->GetNbinsX();
    for (unsigned int iBin = 1; iBin <= nBins; ++iBin)
    {
        elements.push_back(pHist->GetBinContent(iBin));
    }

    return ubsmear::UBMatrix(elements, nBins, 1);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

ubsmear::UBMatrix CrossSectionHelper::GetMatrixFromHist(const std::shared_ptr<TH2F> &pHist)
{
    std::vector< std::vector<float> > elements;

    const unsigned int nBinsX = pHist->GetNbinsX();
    const unsigned int nBinsY = pHist->GetNbinsY();

    for (unsigned int iBinX = 1; iBinX <= nBinsX; ++iBinX)
    {
        // Make a new row
        elements.emplace_back();
        auto &row = elements.back();

        for (unsigned int iBinY = 1; iBinY <= nBinsY; ++iBinY)
        {
            row.push_back(pHist->GetBinContent(iBinX, iBinY));
        }
    }

    return ubsmear::UBMatrix(elements);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

std::shared_ptr<ubsmear::UBMatrix> CrossSectionHelper::FlattenMatrix(const std::shared_ptr<ubsmear::UBMatrix> &pMatrix)
{
    if (!pMatrix)
        throw std::invalid_argument("CrossSectionHelper::FlattenMatrix - Input matrix pointer is null");

    return std::make_shared<ubsmear::UBMatrix>(ubsmear::UBSmearingHelper::Flatten(*pMatrix));
}

// -----------------------------------------------------------------------------------------------------------------------------------------

std::tuple< std::vector<float>, bool, bool > CrossSectionHelper::GetExtendedBinEdges(const float min, const float max, const std::vector<float> &binEdges)
{
    // Ensure that we have at least one input bin
    if (binEdges.size() < 2)
        throw std::invalid_argument("CrossSectionHelper::GetExtendedBinEdges - Fewer than two input bin edges were supplied");

    // Ensure that the input bin edges are sorted
    if (!std::is_sorted(binEdges.begin(), binEdges.end()))
        throw std::invalid_argument("CrossSectionHelper::GetExtendedBinEdges - Input bin edges are not sorted into ascending order");

    // Get the first and last bin edges supplied
    const auto firstEdge = binEdges.front();
    const auto lastEdge = binEdges.back();

    // Check if the first/last edge is already at the min/max
    const auto isFirstEdgeAtMin = (std::abs(firstEdge - min) <= std::numeric_limits<float>::epsilon());
    const auto isLastEdgeAtMax = (std::abs(lastEdge - max) <= std::numeric_limits<float>::epsilon());

    // If the first/last edge isn't at min/max, then make sure that it's not below/above
    if (!isFirstEdgeAtMin && firstEdge < min)
        throw std::invalid_argument("CrossSectionHelper::GetExtendedBinEdges - Lowest input bin edge is smaller than the supplied minimum value");

    if (!isLastEdgeAtMax && lastEdge > max)
        throw std::invalid_argument("CrossSectionHelper::GetExtendedBinEdges - Uppermost input bin edge is larget than the supplied maximum value");

    // Determine if we need to extend to an extra underflow/overflow bin edge
    const auto hasUnderflow = !isFirstEdgeAtMin;
    const auto hasOverflow = !isLastEdgeAtMax;

    // Setup the output bin edges
    std::vector<float> extendedBinEdges;

    // Add the underflow bin if required
    if (hasUnderflow)
        extendedBinEdges.push_back(min);

    // Add the supplied bins
    extendedBinEdges.insert(extendedBinEdges.end(), binEdges.begin(), binEdges.end());

    // Add the overflow bin if required
    if (hasOverflow)
        extendedBinEdges.push_back(max);

    return {extendedBinEdges, hasUnderflow, hasOverflow};
}

// -----------------------------------------------------------------------------------------------------------------------------------------

std::vector<std::string> CrossSectionHelper::GetNominalFluxHistNames(const std::vector<int> &nuPdgsSignal, const std::map<int, std::string> &nuPdgToHistName, const std::string &nomHistPattern)
{
    std::vector<std::string> histNames;

    for (const auto &nuPdg : nuPdgsSignal)
    {
        const auto iter = nuPdgToHistName.find(nuPdg);
        if (iter == nuPdgToHistName.end())
            throw std::logic_error("CrossSectionHelper::GetNominalFluxHistNames - Can't find name corresponding to PDG code " + std::to_string(nuPdg));

        const auto &nuName = iter->second;
        const auto histName = std::regex_replace(nomHistPattern, std::regex("NEUTRINO"), nuName);
        histNames.push_back(histName);
    }

    return histNames;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

std::pair< std::vector<float>, std::vector<float> > CrossSectionHelper::ReadNominalFlux(const std::string &fileName, const std::vector<std::string> &histNames, const float pot)
{
    if (histNames.empty())
        throw std::invalid_argument("CrossSectionHelper::ReadNominalFlux - No histogram names supplied");

    // Open the flux file for reading
    TFile *pFluxFile = new TFile(fileName.c_str(), "READ");
    if (!pFluxFile->IsOpen())
        throw std::invalid_argument("CrossSectionHelper::ReadNominalFlux - Can't open flux file: " + fileName);

    // Get the scaling factor to go from the event rate in the flux file, to the flux itself
    // Here we get the flux in neutrinos/POT/bin/cm2 by scaling the event rate (in the samples) down by POT and the cross-sectional area of the active volume
    // Here we also scale up the fluxes by 1e10 for the sake of comparison just so we are working with reasonable numbers
    const float unitsScaling = 1e10;
    const auto fluxScaleFactor = unitsScaling / (pot * (GeometryHelper::highX - GeometryHelper::lowX) * (GeometryHelper::highY - GeometryHelper::lowY));

    bool isFirstHist = true;
    std::vector<float> fluxBinEdges, fluxBinValuesNominal;

    for (const auto &histName : histNames)
    {
        // Get the nominal flux
        const auto pFluxHist = static_cast<TH1F *>(pFluxFile->Get(histName.c_str()));
        if (!pFluxHist)
            throw std::logic_error("CrossSectionHelper::ReadNominalFlux - Input file doesn't contain histogram with name: " + histName);

        // Get the number of bins
        const unsigned int nFluxBins = pFluxHist->GetNbinsX();
        if (!isFirstHist && nFluxBins != (fluxBinEdges.size() - 1))
            throw std::logic_error("CrossSectionHelper::ReadNominalFlux - Supplied flux histograms have an inconsistent number of bins");

        // Get the flux bin edges and content
        for (unsigned int iBin = 1; iBin <= nFluxBins; ++iBin)
        {
            const auto flux = pFluxHist->GetBinContent(iBin) * fluxScaleFactor;
            const float lowEdge = pFluxHist->GetBinLowEdge(iBin);
            const float binWidth = pFluxHist->GetBinWidth(iBin);

            if (isFirstHist)
            {
                fluxBinEdges.push_back(lowEdge);

                // Add the upper edge of the last bin
                if (iBin == nFluxBins)
                {
                    fluxBinEdges.push_back(lowEdge + binWidth);
                }

                fluxBinValuesNominal.push_back(flux);
            }
            else
            {
                // If this isn't the first flux histogram, then add to the existing flux
                fluxBinValuesNominal.at(iBin - 1) += flux;

                // Check that the bin edges are consistent
                if (std::abs(fluxBinEdges.at(iBin - 1) - lowEdge) > std::numeric_limits<float>::epsilon())
                    throw std::logic_error("CrossSectionHelper::ReadNominalFlux - Supplied flux histograms have an inconsistent binning (low edge)");
            }
        }

        isFirstHist = false;
    }

    return {fluxBinEdges, fluxBinValuesNominal};
}

// -----------------------------------------------------------------------------------------------------------------------------------------

ubsmear::UBMatrix CrossSectionHelper::GetErrorMatrix(const ubsmear::UBMatrix &biasVector, const ubsmear::UBMatrix &covarianceMatrix)
{
    // Make sure the dimensions of the input are valid
    if (biasVector.GetColumns() != 1)
        throw std::invalid_argument("CrossSectionHelper::GetErrorMatrix - Input bias vector isn't a column vector");

    if (!ubsmear::UBMatrixHelper::IsSquare(covarianceMatrix))
        throw std::invalid_argument("CrossSectionHelper::GetErrorMatrix - Input covariance matrix isn't square");

    const auto nBins = covarianceMatrix.GetRows();
    if (biasVector.GetRows() != nBins)
        throw std::invalid_argument("CrossSectionHelper::GetErrorMatrix - Input bias vector and covariance matrix have incompatible dimensions");

    // Make a new matrix for the output
    auto errorMatrix = ubsmear::UBMatrixHelper::GetZeroMatrix(nBins, nBins);
    for (unsigned int iRow = 0; iRow < nBins; ++iRow)
    {
        for (unsigned int iCol = 0; iCol < nBins; ++iCol)
        {
            errorMatrix.SetElement(iRow, iCol, (biasVector.At(iRow, 0) * biasVector.At(iCol, 0)) + covarianceMatrix.At(iRow, iCol));
        }
    }

    return errorMatrix;
}

} // namespace ubcc1pi
