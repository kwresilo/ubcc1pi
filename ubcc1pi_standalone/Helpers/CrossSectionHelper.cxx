/**
 *  @file  ubcc1pi_standalone/Helpers/CrossSectionHelper.cxx
 *
 *  @brief The implementation of the cross section helper class
 */

#include "ubcc1pi_standalone/Helpers/CrossSectionHelper.h"

#include <stdexcept>

namespace ubcc1pi
{

CrossSectionHelper::CrossSection::CrossSection(const std::vector<float> &binEdges, const bool hasUnderflow, const bool hasOverflow, const bool scaleByBinWidth, const InputData &inputData) :
    m_binEdges(binEdges),
    m_hasUnderflow(hasUnderflow),
    m_hasOverflow(hasOverflow),
    m_scaleByBinWidth(scaleByBinWidth),
    m_inputData(inputData),
    m_allEventsNuEnergyNomTrue(this->GetEmptyFluxHist()),
    m_dataSelectedReco(this->GetEmptyHist1D()),
    m_signalSelectedNomRecoTrue(this->GetEmptyHist2D()),
    m_signalAllNomTrue(this->GetEmptyHist1D()),
    m_backgroundSelectedNomReco(this->GetEmptyHist1D()),
    m_shouldResetCache(true)
{
    // Setup the histograms for the universe variations
    for (const auto &[name, nUniverses] : m_inputData.m_systUniverseSizesMap)
    {
        auto &allEventsNuEnergyTrueVect = m_allEventsNuEnergyTrue[name];
        auto &signalSelectedRecoTrueVect = m_signalSelectedRecoTrue[name];
        auto &signalAllTrueVect = m_signalAllTrue[name];
        auto &backgroundSelectedRecoVect = m_backgroundSelectedReco[name];

        for (unsigned int iUni = 0; iUni < nUniverses; ++iUni)
        {
            allEventsNuEnergyTrueVect.push_back(this->GetEmptyFluxHist());
            signalSelectedRecoTrueVect.push_back(this->GetEmptyHist2D());
            signalAllTrueVect.push_back(this->GetEmptyHist1D());
            backgroundSelectedRecoVect.push_back(this->GetEmptyHist1D());
        }
    }
}

// -----------------------------------------------------------------------------------------------------------------------------------------
                
bool CrossSectionHelper::CrossSection::HasUnderflowBin() const
{
    return m_hasUnderflow;
}

// -----------------------------------------------------------------------------------------------------------------------------------------
                
bool CrossSectionHelper::CrossSection::HasOverflowBin() const
{
    return m_hasOverflow;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void CrossSectionHelper::CrossSection::AddNeutrinoEvent(const float trueNuEnergy, const float nominalWeight, const SystematicWeightsMap &systWeightsMap)
{
    // Check the validity of the input map if systematic weights
    if (!CrossSectionHelper::CheckSystematicWeightsMapDimensions(systWeightsMap, m_inputData.m_systUniverseSizesMap))
        throw std::invalid_argument("CrossSection::AddNeutrinoEvent - Input systematic weights map has the wrong dimensions!");

    // Add this event in the nominal universe
    m_allEventsNuEnergyNomTrue->Fill(trueNuEnergy, nominalWeight);
    
    // Add the event in each systematic universe
    for (const auto &[name, nUniverses] : m_inputData.m_systUniverseSizesMap)
    {
        auto &allEventsNuEnergyTrueVect = m_allEventsNuEnergyTrue.at(name);

        const auto &weights = systWeightsMap.at(name);
        for (unsigned int iUni = 0; iUni < nUniverses; ++iUni)
        {
            const auto universeWeight = weights.at(iUni);
            allEventsNuEnergyTrueVect.at(iUni)->Fill(trueNuEnergy, nominalWeight * universeWeight);
        }
    }

    // Mark that the inputs have changed and any cached values need to be re-calculated
    m_shouldResetCache = true;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void CrossSectionHelper::CrossSection::AddSignalEvent(const float trueValue, const float recoValue, const bool isSelected, const float nominalWeight, const SystematicWeightsMap &systWeightsMap)
{
    // Check the validity of the input map if systematic weights
    if (!CrossSectionHelper::CheckSystematicWeightsMapDimensions(systWeightsMap, m_inputData.m_systUniverseSizesMap))
        throw std::invalid_argument("CrossSection::AddSignalEvent - Input systematic weights map has the wrong dimensions!");

    // Add this event in the nominal universe
    m_signalAllNomTrue->Fill(trueValue, nominalWeight);

    if (isSelected)
        m_signalSelectedNomRecoTrue->Fill(trueValue, recoValue, nominalWeight);

    // Add the event in each systematic universe
    for (const auto &[name, nUniverses] : m_inputData.m_systUniverseSizesMap)
    {
        auto &signalAllTrueVect = m_signalAllTrue.at(name);
        auto &signalSelectedRecoTrueVect = m_signalSelectedRecoTrue.at(name);

        const auto &weights = systWeightsMap.at(name);
        for (unsigned int iUni = 0; iUni < nUniverses; ++iUni)
        {
            const auto universeWeight = weights.at(iUni);
            signalAllTrueVect.at(iUni)->Fill(trueValue, nominalWeight * universeWeight);
            
            if (isSelected)
                signalSelectedRecoTrueVect.at(iUni)->Fill(trueValue, recoValue, nominalWeight * universeWeight);
        }
    }
    
    // Mark that the inputs have changed and any cached values need to be re-calculated
    m_shouldResetCache = true;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void CrossSectionHelper::CrossSection::AddSelectedBackgroundEvent(const float recoValue, const float nominalWeight, const SystematicWeightsMap &systWeightsMap)
{
    // Check the validity of the input map if systematic weights
    if (!CrossSectionHelper::CheckSystematicWeightsMapDimensions(systWeightsMap, m_inputData.m_systUniverseSizesMap))
        throw std::invalid_argument("CrossSection::AddSignalEvent - Input systematic weights map has the wrong dimensions!");

    // Add this event in the nominal universe
    m_backgroundSelectedNomReco->Fill(recoValue, nominalWeight);

    // Add the event in each systematic universe
    for (const auto &[name, nUniverses] : m_inputData.m_systUniverseSizesMap)
    {
        auto &backgroundSelectedRecoVect = m_backgroundSelectedReco.at(name);

        const auto &weights = systWeightsMap.at(name);
        for (unsigned int iUni = 0; iUni < nUniverses; ++iUni)
        {
            const auto universeWeight = weights.at(iUni);
            backgroundSelectedRecoVect.at(iUni)->Fill(recoValue, nominalWeight * universeWeight);
        }
    }
    
    // Mark that the inputs have changed and any cached values need to be re-calculated
    m_shouldResetCache = true;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void CrossSectionHelper::CrossSection::AddSelectedBNBDataEvent(const float recoValue)
{
    m_dataSelectedReco->Fill(recoValue);
    
    // Mark that the inputs have changed and any cached values need to be re-calculated
    m_shouldResetCache = true;
}

// -----------------------------------------------------------------------------------------------------------------------------------------
                
std::shared_ptr<TH1F> CrossSectionHelper::CrossSection::GetEmptyFluxHist() const
{
    const std::string name = "xSecHist_" + std::to_string(m_histCount++);

    // Get the binning of the flux histogram
    const unsigned int nBins = m_inputData.m_flux->GetNbinsX();
    std::vector<float> edges;

    for (unsigned int iBin = 1; iBin <= nBins; ++iBin)
    {
        edges.push_back(m_inputData.m_flux->GetBinLowEdge(iBin));
    }

    // Add the last bin edge
    edges.push_back(m_inputData.m_flux->GetBinLowEdge(nBins) + m_inputData.m_flux->GetBinWidth(nBins));

    return std::make_shared<TH1F>(name.c_str(), "", edges.size() - 1, edges.data());
}

// -----------------------------------------------------------------------------------------------------------------------------------------
                
std::shared_ptr<TH1F> CrossSectionHelper::CrossSection::GetEmptyHist1D() const
{
    const std::string name = "xSecHist_" + std::to_string(m_histCount++);
    return std::make_shared<TH1F>(name.c_str(), "", m_binEdges.size() - 1, m_binEdges.data());
}

// -----------------------------------------------------------------------------------------------------------------------------------------
                
std::shared_ptr<TH2F> CrossSectionHelper::CrossSection::GetEmptyHist2D() const
{
    const std::string name = "xSecHist_" + std::to_string(m_histCount++);
    return std::make_shared<TH2F>(name.c_str(), "", m_binEdges.size() - 1, m_binEdges.data(), m_binEdges.size() - 1, m_binEdges.data());
}

// -----------------------------------------------------------------------------------------------------------------------------------------
                
bool CrossSectionHelper::CrossSection::IsUnderOverflowBin(const unsigned int binIndex) const
{
    const auto nBins = m_binEdges.size() - 1;
    if (binIndex > nBins)
        throw std::out_of_range("CrossSection::IsUnderOverflowBin - input bin index " + std::to_string(binIndex) + " is out of range");

    return ((binIndex == 1 && m_hasUnderflow) || (binIndex == nBins && m_hasOverflow));
}

// -----------------------------------------------------------------------------------------------------------------------------------------
                
std::shared_ptr<TH2F> CrossSectionHelper::CrossSection::GetSmearingMatrix(const std::shared_ptr<TH2F> &signalSelectedRecoTrue) const
{
    // Check the input matrix has the correct dimensions
    const auto nBins = m_binEdges.size() - 1;

    if (static_cast<unsigned int>(signalSelectedRecoTrue->GetXaxis()->GetNbins()) != nBins ||
        static_cast<unsigned int>(signalSelectedRecoTrue->GetYaxis()->GetNbins()) != nBins )
    {
        throw std::logic_error("CrossSection::GetSmearingMatrix - input matrix has the wrong dimensions");
    }

    auto smearingMatrix = this->GetEmptyHist2D();
    for (unsigned int iTrue = 1; iTrue <= nBins; ++iTrue)
    {
        // Get the sum of the reco values for this true bin
        float summedValue = 0.f;
        for (unsigned int iReco = 1; iReco <= nBins; ++iReco)
        {
            summedValue += signalSelectedRecoTrue->GetBinContent(iTrue, iReco);
        }

        if (summedValue <= std::numeric_limits<float>::epsilon())
        {
            // ATTN if a true bin has no entries, then we have no information to work with
            // Here we make the choice to give each reco bin an equal probability
            const auto value = 1.f / static_cast<float>(nBins);
            for (unsigned int iReco = 1; iReco <= nBins; ++iReco)
            {
                smearingMatrix->SetBinContent(iTrue, iReco, value);
            }
        }
        else
        {
            // Normalise the by a constant factor bins such that the columns add to unity
            const auto norm = 1.f / summedValue;

            for (unsigned int iReco = 1; iReco <= nBins; ++iReco)
            {
                const auto value = signalSelectedRecoTrue->GetBinContent(iTrue, iReco) * norm;
                smearingMatrix->SetBinContent(iTrue, iReco, value);
            }
        }
    }

    return smearingMatrix;
}

// -----------------------------------------------------------------------------------------------------------------------------------------
                
std::shared_ptr<TH1F> CrossSectionHelper::CrossSection::Smear(const std::shared_ptr<TH1F> &values, const std::shared_ptr<TH2F> &smearingMatrix) const
{
    // Check the input objects have the correct dimensions
    const auto nBins = m_binEdges.size() - 1;

    if (static_cast<unsigned int>(values->GetXaxis()->GetNbins()) != nBins)
        throw std::logic_error("CrossSection::Smear - input values has the wrong dimensions");

    if (static_cast<unsigned int>(smearingMatrix->GetXaxis()->GetNbins()) != nBins ||
        static_cast<unsigned int>(smearingMatrix->GetYaxis()->GetNbins()) != nBins)
    {
        throw std::logic_error("CrossSection::Smear - input matrix has the wrong dimensions");
    }


    auto smearedValues = this->GetEmptyHist1D();

    for (unsigned int iReco = 1; iReco <= nBins; ++iReco)
    {
        // Smear the true values in this reco bin
        float smearedValue = 0.f;
        for (unsigned int iTrue = 1; iTrue <= nBins; ++iTrue)
        {
            smearedValue += values->GetBinContent(iTrue) * smearingMatrix->GetBinContent(iTrue, iReco);
        }

        smearedValues->SetBinContent(iReco, smearedValue);
    }

    return smearedValues;
}

// -----------------------------------------------------------------------------------------------------------------------------------------
                
std::shared_ptr<TH1F> CrossSectionHelper::CrossSection::GetSmearedEfficiency(const std::shared_ptr<TH2F> &signalSelectedRecoTrue, const std::shared_ptr<TH1F> &signalAllTrue) const
{
    const auto smearingMatrix = this->GetSmearingMatrix(signalSelectedRecoTrue);
    return this->GetSmearedEfficiency(smearingMatrix, signalSelectedRecoTrue, signalAllTrue);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

std::shared_ptr<TH1F> CrossSectionHelper::CrossSection::GetSmearedEfficiency(const std::shared_ptr<TH2F> &smearingMatrix, const std::shared_ptr<TH2F> &signalSelectedRecoTrue, const std::shared_ptr<TH1F> &signalAllTrue) const
{
    // Check the input objects have the correct dimensions
    const auto nBins = m_binEdges.size() - 1;

    if (static_cast<unsigned int>(signalAllTrue->GetXaxis()->GetNbins()) != nBins)
        throw std::logic_error("CrossSection::GetSmearedEfficiency - input object signalAllTrue has the wrong dimensions");

    if (static_cast<unsigned int>(smearingMatrix->GetXaxis()->GetNbins()) != nBins ||
        static_cast<unsigned int>(smearingMatrix->GetYaxis()->GetNbins()) != nBins)
    {
        throw std::logic_error("CrossSection::GetSmearedEfficiency - input smearing matrix has the wrong dimensions");
    }
    
    if (static_cast<unsigned int>(signalSelectedRecoTrue->GetXaxis()->GetNbins()) != nBins ||
        static_cast<unsigned int>(signalSelectedRecoTrue->GetYaxis()->GetNbins()) != nBins)
    {
        throw std::logic_error("CrossSection::GetSmearedEfficiency - input signalSelectedRecoTrue matrix has the wrong dimensions");
    }

    // Get the total selected events per true bin
    auto signalSelectedTrue = this->GetEmptyHist1D();
    for (unsigned int iTrue = 1; iTrue <= nBins; ++iTrue)
    {
        // Get the sum of the reco values for this true bin
        float summedValue = 0.f;
        for (unsigned int iReco = 1; iReco <= nBins; ++iReco)
        {
            summedValue += signalSelectedRecoTrue->GetBinContent(iTrue, iReco);
        }

        signalSelectedTrue->SetBinContent(iTrue, summedValue);
    }

    // Now smear the numerator and denominators of the efficiency
    const auto signalSelectedSmearedReco = this->Smear(signalSelectedTrue, smearingMatrix);
    const auto signalAllSmearedReco = this->Smear(signalAllTrue, smearingMatrix);

    // Get the smeared efficiency
    auto effSmearedReco = this->GetEmptyHist1D();
    for (unsigned int iReco = 1; iReco <= nBins; ++iReco)
    {
        // Skip under/overflow bins
        if (this->IsUnderOverflowBin(iReco))
            continue;

        const auto numerator = signalSelectedSmearedReco->GetBinContent(iReco);
        const auto denominator = signalAllSmearedReco->GetBinContent(iReco);

        const auto isDenominatorZero = denominator <= std::numeric_limits<float>::epsilon();

        // ATTN here we check a zero denominator - in this we set a dummy efficiency of -inf
        const auto value = isDenominatorZero ? -std::numeric_limits<float>::max() : (numerator / denominator);
        effSmearedReco->SetBinContent(iReco, value);
    }

    return effSmearedReco;
}

// -----------------------------------------------------------------------------------------------------------------------------------------
                
void CrossSectionHelper::CrossSection::ClearCache()
{
    m_crossSectionCache.clear();
    m_smearingMatrixCache.clear();
    m_shouldResetCache = false;
}

// -----------------------------------------------------------------------------------------------------------------------------------------
                
std::shared_ptr<TH1F> CrossSectionHelper::CrossSection::GetCachedCrossSection(const std::string &systParameter, const unsigned int universeIndex)
{
    // Check that we have an entry for the requested systematic paramter
    const auto iter = m_inputData.m_systUniverseSizesMap.find(systParameter);
    if (iter == m_inputData.m_systUniverseSizesMap.end())
        throw std::invalid_argument("CrossSection::GetCachedCrossSection - unknown parameter: \"" + systParameter + "\"");

    // Get the number of universes
    const auto nUniverses = iter->second;
    if (universeIndex >= nUniverses)
        throw std::range_error("CrossSection::GetCachedCrossSection - Universe index for parameter: \"" + systParameter + "\" is out of range");

    // Reset the cached cross-section if required
    if (m_shouldResetCache)
        this->ClearCache();
        
    // Get the cached cross-section (if it exists)
    const auto paramIter = m_crossSectionCache.find(systParameter);
    if (paramIter != m_crossSectionCache.end())
    {
        const auto universeIter = paramIter->second.find(universeIndex);
        if (universeIter != paramIter->second.end())
            return universeIter->second;
    }

    // We don't have an entry for this cross-section, so make it!
    
    // Determine if this is a flux parameter and if so, reweight the total flux
    const auto isFluxParam = (std::find(m_inputData.m_fluxParameters.begin(), m_inputData.m_fluxParameters.end(), systParameter) != m_inputData.m_fluxParameters.end());

    const auto totalFlux = CrossSectionHelper::GetTotalFlux(isFluxParam
            ? this->GetReweightedFlux(m_allEventsNuEnergyTrue.at(systParameter).at(universeIndex))
            : m_inputData.m_flux);

    // Get the smeared efficiency distribution and cache the smearing matrix if required
    const auto smearingMatrix = this->GetCachedSmearingMatrix(systParameter, universeIndex);
    const auto smearedEff = this->GetSmearedEfficiency(smearingMatrix, m_signalSelectedRecoTrue.at(systParameter).at(universeIndex), m_signalAllTrue.at(systParameter).at(universeIndex));

    // Get the cross-section and cache it
    auto &cachedXSec = m_crossSectionCache[systParameter][universeIndex];
    cachedXSec = this->GetCrossSection(                                         // Get the total cross-section using:
        m_dataSelectedReco,                                                     //   - Selected BNB data
        m_backgroundSelectedReco.at(systParameter).at(universeIndex),           //   - Backgrounds in this universe
        smearedEff, totalFlux);                                                 //   - Smeared efficiency and flux in this universe

    return cachedXSec;
}

// -----------------------------------------------------------------------------------------------------------------------------------------
                
std::shared_ptr<TH2F> CrossSectionHelper::CrossSection::GetCachedSmearingMatrix(const std::string &systParameter, const unsigned int universeIndex)
{
    // Check that we have an entry for the requested systematic paramter
    const auto iter = m_inputData.m_systUniverseSizesMap.find(systParameter);
    if (iter == m_inputData.m_systUniverseSizesMap.end())
        throw std::invalid_argument("CrossSection::GetCachedSmearingMatrix - unknown parameter: \"" + systParameter + "\"");

    // Get the number of universes
    const auto nUniverses = iter->second;
    if (universeIndex >= nUniverses)
        throw std::range_error("CrossSection::GetCachedSmearingMatrix - Universe index for parameter: \"" + systParameter + "\" is out of range");

    // Reset the cached cross-section if required
    if (m_shouldResetCache)
        this->ClearCache();
        
    // Get the cached smearing matrix (if it exists)
    const auto paramIter = m_smearingMatrixCache.find(systParameter);
    if (paramIter != m_smearingMatrixCache.end())
    {
        const auto universeIter = paramIter->second.find(universeIndex);
        if (universeIter != paramIter->second.end())
            return universeIter->second;
    }

    // We don't have an entry for this smearing matrix, so make it!
    auto &cachedSmearingMatrix = m_smearingMatrixCache[systParameter][universeIndex];
    cachedSmearingMatrix = this->GetSmearingMatrix(m_signalSelectedRecoTrue.at(systParameter).at(universeIndex));

    return cachedSmearingMatrix;
}

// -----------------------------------------------------------------------------------------------------------------------------------------
                
std::shared_ptr<TH1F> CrossSectionHelper::CrossSection::GetCrossSection(const std::shared_ptr<TH1F> &selected, const std::shared_ptr<TH1F> &background, const std::shared_ptr<TH1F> &smearedEff, const float totalFlux) const
{
    const auto normFactor = 1.f / (totalFlux * m_inputData.m_exposurePOT * m_inputData.m_nTargets);
    
    // Check the input objects have the correct dimensions
    const auto nBins = m_binEdges.size() - 1;

    if (static_cast<unsigned int>(selected->GetXaxis()->GetNbins()) != nBins)
        throw std::logic_error("CrossSection::GetCrossSection - input selected distribution has the wrong dimensions");

    if (static_cast<unsigned int>(background->GetXaxis()->GetNbins()) != nBins)
        throw std::logic_error("CrossSection::GetCrossSection - input background distribution has the wrong dimensions");

    if (static_cast<unsigned int>(smearedEff->GetXaxis()->GetNbins()) != nBins)
        throw std::logic_error("CrossSection::GetCrossSection - input smeared efficiency distribution has the wrong dimensions");

    auto crossSection = this->GetEmptyHist1D();
    for (unsigned int iReco = 1; iReco <= nBins; ++iReco)
    {
        // Skip under/overflow bins
        if (this->IsUnderOverflowBin(iReco))
            continue;

        // In some universes, the efficiency may be incalculable (and so a negative dummy value is used) - here we use a dummy value again
        const auto efficiency = smearedEff->GetBinContent(iReco);
        if (std::abs(std::abs(efficiency) - std::numeric_limits<float>::max()) <= std::numeric_limits<float>::epsilon() ||
            efficiency <= std::numeric_limits<float>::epsilon())
        {
            crossSection->SetBinContent(iReco, -std::numeric_limits<float>::max());
            continue;
        }

        const auto binWidth = m_scaleByBinWidth ? (m_binEdges.at(iReco) - m_binEdges.at(iReco - 1)) : 1.f;
        const auto value = normFactor * (selected->GetBinContent(iReco) - background->GetBinContent(iReco)) / (efficiency * binWidth);

        crossSection->SetBinContent(iReco, value);
    }

    return crossSection;
}

// -----------------------------------------------------------------------------------------------------------------------------------------
                
std::shared_ptr<TH2F> CrossSectionHelper::CrossSection::GetSmearingMatrix() const
{
    return this->GetSmearingMatrix(m_signalSelectedNomRecoTrue);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

std::shared_ptr<TH1F> CrossSectionHelper::CrossSection::GetCrossSection() const
{
    const auto smearedEff = this->GetSmearedEfficiency(m_signalSelectedNomRecoTrue, m_signalAllNomTrue);
    const auto totalFlux = CrossSectionHelper::GetTotalFlux(m_inputData.m_flux);
    return this->GetCrossSection(m_dataSelectedReco, m_backgroundSelectedNomReco, smearedEff, totalFlux);
}

// -----------------------------------------------------------------------------------------------------------------------------------------
                
std::shared_ptr<TH1F> CrossSectionHelper::CrossSection::GetCrossSectionStatUncertainty() const
{
    // Setup the output histogram
    const auto nBins = m_binEdges.size() - 1;
    const auto nAnalysisBins = nBins - (m_hasUnderflow ? 1 : 0) - (m_hasOverflow ? 1 : 0);
    const std::string name = "xSecStatErr_" + std::to_string(m_histCount++);
    auto statErrHist = std::make_shared<TH1F>(name.c_str(), "", nAnalysisBins, 0, nAnalysisBins);

    const auto smearedEff = this->GetSmearedEfficiency(m_signalSelectedNomRecoTrue, m_signalAllNomTrue);
    const auto totalFlux = CrossSectionHelper::GetTotalFlux(m_inputData.m_flux);
    const auto normFactor = 1.f / (totalFlux * m_inputData.m_exposurePOT * m_inputData.m_nTargets);

    for (unsigned int iBin = 1u; iBin <= nBins; ++iBin)
    {
        if (this->IsUnderOverflowBin(iBin))
            continue;
        
        const auto binWidth = m_scaleByBinWidth ? (m_binEdges.at(iBin) - m_binEdges.at(iBin - 1)) : 1.f;
        const auto nDataEvents = m_dataSelectedReco->GetBinContent(iBin);
        const auto efficiency = smearedEff->GetBinContent(iBin);
        
        if (std::abs(std::abs(efficiency) - std::numeric_limits<float>::max()) <= std::numeric_limits<float>::epsilon() ||
            efficiency <= std::numeric_limits<float>::epsilon())
        {
            throw std::logic_error("CrossSection::GetCrossSectionStatUncertainty - Found bin in which smeared efficiency couldn't be calculated!");
        }

        // Get the uncertainty on the number of data events
        const auto nDataEventsErr = AnalysisHelper::GetCountUncertainty(nDataEvents);

        // Get the corresponding uncertainty on the cross-section
        // Cross-section = normFactor * (nDataEvents - nBackgroundsEvents) / (efficiency * binWidth)
        const auto uncertainty = nDataEventsErr * normFactor / (efficiency * binWidth);

        statErrHist->SetBinContent(iBin - (m_hasUnderflow ? 1 : 0), uncertainty);
    }

    return statErrHist;
}

// -----------------------------------------------------------------------------------------------------------------------------------------
                
std::map< std::string, std::vector< std::tuple<unsigned int, unsigned int, std::shared_ptr<TGraph> > > > CrossSectionHelper::CrossSection::GetCrossSectionBinScatterPlots()
{
    std::map< std::string, std::vector< std::tuple<unsigned int, unsigned int, std::shared_ptr<TGraph> > > > outputMap;

    for (const auto &[param, nUniverses] : m_inputData.m_systUniverseSizesMap)
    {
        outputMap.emplace(param, this->GetCrossSectionBinScatterPlots(param));
    }

    return outputMap;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

std::vector< std::tuple<unsigned int, unsigned int, std::shared_ptr<TGraph> > > CrossSectionHelper::CrossSection::GetCrossSectionBinScatterPlots(const std::string &systParameter)
{
    // Check that we have an entry for the requested systematic paramter
    const auto iter = m_inputData.m_systUniverseSizesMap.find(systParameter);
    if (iter == m_inputData.m_systUniverseSizesMap.end())
        throw std::invalid_argument("CrossSection::GetCrossSectionBinScatterPlots - unknown parameter: \"" + systParameter + "\"");
   
    const auto nUniverses = iter->second;

    // Setup the output vector
    std::vector< std::tuple<unsigned int, unsigned int, std::shared_ptr<TGraph> > > outputVect;

    // Loop over all pairs of bins
    const auto nBins = m_binEdges.size() - 1;
    for (unsigned int iBin = 1u; iBin <= nBins; ++iBin)
    {
        if (this->IsUnderOverflowBin(iBin))
            continue;
    
        for (unsigned int jBin = 1u; jBin <= iBin; ++jBin)
        {
            if (this->IsUnderOverflowBin(jBin))
                continue;

            // Get the cross-section in each universe from the cache - or add it if required
            std::vector<float> iValues, jValues;
            float minValue = std::numeric_limits<float>::max();
            float maxValue = -std::numeric_limits<float>::max();

            for (unsigned int iUni = 0u; iUni < nUniverses; ++iUni)
            {    
                // Get the values of the cross-section in the current pair of bins and make sure they are valid
                const auto crossSection = this->GetCachedCrossSection(systParameter, iUni);

                const float iValue = crossSection->GetBinContent(iBin);
                if (std::abs(std::abs(iValue) - std::numeric_limits<float>::max()) < std::numeric_limits<float>::epsilon())
                    continue;
                
                const float jValue = crossSection->GetBinContent(jBin);
                if (std::abs(std::abs(jValue) - std::numeric_limits<float>::max()) < std::numeric_limits<float>::epsilon())
                    continue;

                // Store these values
                iValues.push_back(iValue);
                jValues.push_back(jValue);

                // Get the minimum and maximum values
                minValue = std::min(minValue, std::min(iValue, jValue));
                maxValue = std::max(maxValue, std::max(iValue, jValue));
            }

            // Store the values in a scatter plot
            const auto nValues = iValues.size();
            auto pGraph = std::make_shared<TGraph>(nValues, iValues.data(), jValues.data());

            // Set the range of the graph to ensure it's square
            pGraph->GetXaxis()->SetLimits(minValue, maxValue);
            pGraph->GetYaxis()->SetLimits(minValue, maxValue);

            outputVect.emplace_back(iBin, jBin, pGraph);
        }
    }

    return outputVect;
}

// -----------------------------------------------------------------------------------------------------------------------------------------
                
CrossSectionHelper::CovarianceBiasPair CrossSectionHelper::CrossSection::GetCrossSectionCovarianceMatrix(const std::string &systParameter)
{
    // Check that we have an entry for the requested systematic paramter
    const auto iter = m_inputData.m_systUniverseSizesMap.find(systParameter);
    if (iter == m_inputData.m_systUniverseSizesMap.end())
        throw std::invalid_argument("CrossSection::GetCrossSectionCovarianceMatrix - unknown parameter: \"" + systParameter + "\"");
   
    // Get the number of universes
    const auto nUniverses = iter->second;
    if (nUniverses == 0)
        throw std::logic_error("CrossSection::GetCrossSectionCovarianceMatrix - No universes for parameter: \"" + systParameter + "\"");

    // Get the nominal cross section
    const auto crossSectionNom = this->GetCrossSection();
    
    // Get the cross-section in each universe and the mean cross-section over all universes
    std::vector< std::shared_ptr<TH1F> > crossSections;

    auto nUniversesPerBin = this->GetEmptyHist1D();
    auto summedCrossSection = this->GetEmptyHist1D();
    const auto nBins = m_binEdges.size() - 1;

    //// BEGIN DEBUG
    std::cout << "Getting covariance matrix for parameter: " << systParameter << std::endl;
    //// END DEBUG
    
    for (unsigned int iUni = 0; iUni < nUniverses; ++iUni)
    {
        //// BEGIN DEBUG
        AnalysisHelper::PrintLoadingBar(iUni, nUniverses);
        //// END DEBUG
        
        // Get the cross-section from the cache (or add it if we haven't calculated this cross section yet)
        const auto crossSection = this->GetCachedCrossSection(systParameter, iUni);
        crossSections.push_back(crossSection);

        // Add this to the sum for each bin
        for (unsigned int iBin = 1; iBin <= nBins; ++iBin)
        {
            if (this->IsUnderOverflowBin(iBin))
                continue;

            // Skip dummy bins in which the cross-section couldn't be extracted
            const auto binValue = crossSection->GetBinContent(iBin);
            if (std::abs(std::abs(binValue) - std::numeric_limits<float>::max()) < std::numeric_limits<float>::epsilon())
                continue;

            summedCrossSection->SetBinContent(iBin, summedCrossSection->GetBinContent(iBin) + binValue);
            nUniversesPerBin->SetBinContent(iBin, nUniversesPerBin->GetBinContent(iBin) + 1.f);
        }
    }

    // Get the mean cross-section
    auto meanCrossSection = this->GetEmptyHist1D();
    for (unsigned int iBin = 1; iBin <= nBins; ++iBin)
    {
        if (this->IsUnderOverflowBin(iBin))
            continue;

        const auto numerator = summedCrossSection->GetBinContent(iBin);
        const auto denominator = nUniversesPerBin->GetBinContent(iBin);

        // If the denominator is zero then use a dummy value
        const auto isDemoninatorZero = (denominator <= std::numeric_limits<float>::epsilon());
        const auto mean = isDemoninatorZero ? -std::numeric_limits<float>::max() : (numerator / denominator);
        meanCrossSection->SetBinContent(iBin, mean);
    }

    // Get the covairance matrix and bias vector
    const auto nAnalysisBins = nBins - (m_hasUnderflow ? 1 : 0) - (m_hasOverflow ? 1 : 0);
    const std::string covName = "xSecCov_" + std::to_string(m_histCount++);
    const std::string biasName = "xSecBias_" + std::to_string(m_histCount++);
    CovarianceBiasPair pair (
        std::make_shared<TH2F>(covName.c_str(), "", nAnalysisBins, 0, nAnalysisBins, nAnalysisBins, 0, nAnalysisBins),
        std::make_shared<TH1F>(biasName.c_str(), "", nAnalysisBins, 0, nAnalysisBins)
    );
    
    auto &covarianceMatrix = pair.first;
    auto &biasVector = pair.second;

    for (unsigned int iBin = 1; iBin <= nBins; ++iBin)
    {
        if (this->IsUnderOverflowBin(iBin))
            continue;

        const auto meanI = meanCrossSection->GetBinContent(iBin);
        const auto nomI = crossSectionNom->GetBinContent(iBin);

        // Check for dummy bins and warn if they exists
        if (std::abs(std::abs(nomI) - std::numeric_limits<float>::max()) <= std::numeric_limits<float>::epsilon())
        {
            std::cout << "WARNING - Bin: " << iBin << " of nominal cross-section is invalid" << std::endl;
            continue;
        }
        
        if (std::abs(std::abs(meanI) - std::numeric_limits<float>::max()) <= std::numeric_limits<float>::epsilon())
        {
            std::cout << "WARNING - Bin: " << iBin << " of has an incalculable mean value over all universe variations" << std::endl;
            continue;
        }

        // Set the bias vector element
        const auto bias = meanI - nomI;
        const auto iBinOutput = iBin - (m_hasUnderflow ? 1 : 0);
        biasVector->SetBinContent(iBinOutput, bias);

        for (unsigned int jBin = 1; jBin <= nBins; ++jBin)
        {
            if (this->IsUnderOverflowBin(jBin))
                continue;
        
            const auto meanJ = meanCrossSection->GetBinContent(jBin);
            const auto nomJ = crossSectionNom->GetBinContent(jBin);

            // Skip dummy bins
            if (std::abs(std::abs(nomJ) - std::numeric_limits<float>::max()) <= std::numeric_limits<float>::epsilon())
                continue;

            if (std::abs(std::abs(meanJ) - std::numeric_limits<float>::max()) <= std::numeric_limits<float>::epsilon())
                continue;

            float numerator = 0.f;
            float denominator = 0.f;
            for (unsigned int iUni = 0; iUni < nUniverses; ++iUni)
            {
                const auto &crossSection = crossSections.at(iUni);
                const auto sigmaI = crossSection->GetBinContent(iBin);
                const auto sigmaJ = crossSection->GetBinContent(jBin);

                // Skip dummy bins
                if (std::abs(std::abs(sigmaI) - std::numeric_limits<float>::max()) <= std::numeric_limits<float>::epsilon())
                    continue;
                
                if (std::abs(std::abs(sigmaJ) - std::numeric_limits<float>::max()) <= std::numeric_limits<float>::epsilon())
                    continue;
                
                numerator   += (sigmaI - meanI) * (sigmaJ - meanJ);
                denominator += 1.f;
            }

            if (denominator <= std::numeric_limits<float>::epsilon())
            {
                throw std::logic_error("CrossSection::GetCrossSectionCovarianceMatrix - For parameter: \"" + systParameter + "\", " + 
                    "can't calculate covariance in bin: " + std::to_string(iBin) + ", " + std::to_string(jBin));
            }

            // Set the covariance matrix element
            const auto covariance = numerator / denominator;
            const auto jBinOutput = jBin - (m_hasUnderflow ? 1 : 0);
            covarianceMatrix->SetBinContent(iBinOutput, jBinOutput, covariance);
        }
    }

    return pair;
}

// -----------------------------------------------------------------------------------------------------------------------------------------
                
CrossSectionHelper::CovarianceBiasPair CrossSectionHelper::CrossSection::GetSmearingMatrixCovarianceMatrix(const std::string &systParameter)
{
    // Check that we have an entry for the requested systematic paramter
    const auto iter = m_inputData.m_systUniverseSizesMap.find(systParameter);
    if (iter == m_inputData.m_systUniverseSizesMap.end())
        throw std::invalid_argument("CrossSection::GetSmearingMatrixCovarianceMatrix - unknown parameter: \"" + systParameter + "\"");
   
    // Get the number of universes
    const auto nUniverses = iter->second;
    if (nUniverses == 0)
        throw std::logic_error("CrossSection::GetSmearingMatrixCovarianceMatrix - No universes for parameter: \"" + systParameter + "\"");
    
    // Setup the covairance matrix and bias vector (here we flatten the bins of the smearing matrix into 1D)
    const auto nBins = m_binEdges.size() - 1;
    const auto nSmearingBins = nBins*nBins;

    const std::string covName = "xSecSmearingCov_" + std::to_string(m_histCount++);
    const std::string biasName = "xSecSmearingBias_" + std::to_string(m_histCount++);
    CovarianceBiasPair pair (
        std::make_shared<TH2F>(covName.c_str(), "", nSmearingBins, 0, nSmearingBins, nSmearingBins, 0, nSmearingBins),
        std::make_shared<TH1F>(biasName.c_str(), "", nSmearingBins, 0, nSmearingBins)
    );
    
    auto &covarianceMatrix = pair.first;
    auto &biasVector = pair.second;

    //// BEGIN DEBUG
    std::cout << "Getting mean smearing matrix for parameter: " << systParameter << std::endl;
    //// END DEBUG

    // Get the mean smearing matrix
    auto meanSmearingMatrix = this->GetEmptyHist2D();
    for (unsigned int iUni = 0; iUni < nUniverses; ++iUni)
    {
        //// BEGIN DEBUG
        AnalysisHelper::PrintLoadingBar(iUni, nUniverses);
        //// END DEBUG

        // Get the smearing matrix from the cache (or add it if we haven't calculated this smearing matrix yet)
        const auto smearingMatrix = this->GetCachedSmearingMatrix(systParameter, iUni);
        meanSmearingMatrix->Add(smearingMatrix.get());
    }
    const auto universeNormalisation = 1.f / static_cast<float>(nUniverses);
    meanSmearingMatrix->Scale(universeNormalisation);

    // Get the nominal smearing matrix
    const auto nomSmearingMatrix = this->GetSmearingMatrix();

    // Flatten the mean and nominal smearing matrices
    std::vector<float> meanSmearingMatrixFlat, nomSmearingMatrixFlat;
    for (unsigned int iTrue = 1u; iTrue <= nBins; ++iTrue)
    {
        for (unsigned int iReco = 1u; iReco <= nBins; ++iReco)
        {
            meanSmearingMatrixFlat.push_back(meanSmearingMatrix->GetBinContent(iTrue, iReco));
            nomSmearingMatrixFlat.push_back(nomSmearingMatrix->GetBinContent(iTrue, iReco));
        }
    }

    // Get the bias vector
    for (unsigned int iFlat = 0u; iFlat < nSmearingBins; ++iFlat)
    {
        const auto bias = meanSmearingMatrixFlat.at(iFlat) - nomSmearingMatrixFlat.at(iFlat);
        biasVector->SetBinContent(iFlat + 1, bias); // +1 to account for root enumerating bins from 1
    }
    
    //// BEGIN DEBUG
    std::cout << "Getting covariance matrix of smearing matrix for parameter: " << systParameter << std::endl;
    //// END DEBUG

    // Get the covariance matrix
    for (unsigned int iUni = 0u; iUni < nUniverses; ++iUni)
    {
        //// BEGIN DEBUG
        AnalysisHelper::PrintLoadingBar(iUni, nUniverses);
        //// END DEBUG
        
        // Get the smearing matrix from the cache
        const auto smearingMatrix = this->GetCachedSmearingMatrix(systParameter, iUni);

        // Flatten the smearing matrix
        std::vector<float> smearingMatrixFlat;
        for (unsigned int iTrue = 1u; iTrue <= nBins; ++iTrue)
        {
            for (unsigned int iReco = 1u; iReco <= nBins; ++iReco)
            {
                smearingMatrixFlat.push_back(smearingMatrix->GetBinContent(iTrue, iReco));
            }
        }

        // Loop over the flattened indices and add to the covariance matrix
        for (unsigned int iFlat = 0u; iFlat < nSmearingBins; ++iFlat)
        {
            for (unsigned int jFlat = 0u; jFlat <= iFlat; ++jFlat)
            {
                const auto currentValue = covarianceMatrix->GetBinContent(iFlat + 1, jFlat + 1); // +1 to account for root enumerating bins from 1
                const auto newValue = currentValue + (
                        (smearingMatrixFlat.at(iFlat) - meanSmearingMatrixFlat.at(iFlat)) * 
                        (smearingMatrixFlat.at(jFlat) - meanSmearingMatrixFlat.at(jFlat)));

                // The matrix is symmetric so we can switch iFlat -> jFlat
                covarianceMatrix->SetBinContent(iFlat + 1, jFlat + 1, newValue);
                covarianceMatrix->SetBinContent(jFlat + 1, iFlat + 1, newValue);
            }
        }
    }
    covarianceMatrix->Scale(universeNormalisation);

    return pair;
}

// -----------------------------------------------------------------------------------------------------------------------------------------
                
std::shared_ptr<TH1F> CrossSectionHelper::CrossSection::GetReweightedFlux(const std::shared_ptr<TH1F> &nuEnergyUniverse) const
{
    // Check the input histogram has the right number of bins
    const unsigned int nBins = m_inputData.m_flux->GetNbinsX();
    if (static_cast<unsigned int>(nuEnergyUniverse->GetNbinsX()) != nBins)
        throw std::invalid_argument("CrossSection::GetReweightedFlux - Input distribution has the wrong number of bins");

    // Get an empty histogram to fill with the re-weighted distribution
    auto reweightedFlux = this->GetEmptyFluxHist();

    // Reweight the flux
    for (unsigned int iBin = 1; iBin < nBins; ++iBin)
    {
        // Get the number of event in the nominal universe and the universe supplied
        const auto nEventsNom = m_allEventsNuEnergyNomTrue->GetBinContent(iBin);
        const auto nEventUni = nuEnergyUniverse->GetBinContent(iBin);

        // Get the ratio of the number of events - if not possible use a unit weight
        const auto weight = (nEventsNom <= std::numeric_limits<float>::epsilon()) ? 1.f : (nEventUni / nEventsNom);

        // Reweight the flux by this weight
        const auto fluxNom = m_inputData.m_flux->GetBinContent(iBin);
        const auto fluxUni = fluxNom * weight;

        reweightedFlux->SetBinContent(iBin, fluxUni);
    }

    return reweightedFlux;
}

// -----------------------------------------------------------------------------------------------------------------------------------------
                
std::map< std::string, std::shared_ptr<TH2F> > CrossSectionHelper::CrossSection::GetFluxVariations() const
{
    std::map< std::string, std::shared_ptr<TH2F> > outputMap;

    for (const auto &param : m_inputData.m_fluxParameters)
        outputMap.emplace(param, this->GetFluxVariations(param));

    return outputMap;
}

// -----------------------------------------------------------------------------------------------------------------------------------------
                
std::shared_ptr<TH2F> CrossSectionHelper::CrossSection::GetFluxVariations(const std::string &systParameter) const
{
    if (std::find(m_inputData.m_fluxParameters.begin(), m_inputData.m_fluxParameters.end(), systParameter) == m_inputData.m_fluxParameters.end())
        throw std::invalid_argument("CrossSection::GetFluxVariations - Input name \"" + systParameter + "\" is not one of the flux parameters supplied");

    const auto &allEventsNuEnergyTrueVect = m_allEventsNuEnergyTrue.at(systParameter);
    const auto nUniverses = allEventsNuEnergyTrueVect.size();

    // Get the reweighted flux in each universe
    float maxFlux = -std::numeric_limits<float>::max();
    float minFlux = +std::numeric_limits<float>::max();
    std::vector< std::shared_ptr<TH1F> > fluxVariations;
    for (unsigned int iUni = 0; iUni < nUniverses; ++iUni)
    {
        const auto flux = this->GetReweightedFlux(allEventsNuEnergyTrueVect.at(iUni));

        maxFlux = std::max(maxFlux, static_cast<float>(flux->GetMaximum()));
        minFlux = std::min(minFlux, static_cast<float>(flux->GetMinimum()));
        fluxVariations.push_back(flux);
    }

    // Setup the output histogram
    const std::string name = "xSecHist_" + std::to_string(m_histCount++);

    // Get the binning of the flux histogram
    const unsigned int nBinsX = m_inputData.m_flux->GetNbinsX();
    std::vector<float> edgesX;
    for (unsigned int iBin = 1; iBin <= nBinsX; ++iBin)
        edgesX.push_back(m_inputData.m_flux->GetBinLowEdge(iBin));

    // Add the last bin edge
    edgesX.push_back(m_inputData.m_flux->GetBinLowEdge(nBinsX) + m_inputData.m_flux->GetBinWidth(nBinsX));

    // Get the binning in Y
    const auto nBinsY = 100u;
    std::vector<float> edgesY;
    for (unsigned int iBin = 0; iBin < nBinsY; ++iBin)
        edgesY.push_back(minFlux + static_cast<float>(iBin) * (maxFlux - minFlux) / static_cast<float>(nBinsY - 1));

    auto outputHist = std::make_shared<TH2F>(name.c_str(), "", edgesX.size() - 1, edgesX.data(), edgesY.size() - 1, edgesY.data());
    
    // Fill the histogram with the flux variations
    for (const auto &flux : fluxVariations)
    {
        for (unsigned int iBinX = 1; iBinX <= nBinsX; ++iBinX)
        {
            const auto &fluxValue = flux->GetBinContent(iBinX);
            const auto nuEnergy = 0.5f * (edgesX.at(iBinX - 1) + edgesX.at(iBinX));

            outputHist->Fill(nuEnergy, fluxValue);
        }
    }

    return outputHist;
}

// -----------------------------------------------------------------------------------------------------------------------------------------
                
std::map< std::string, CrossSectionHelper::CovarianceBiasPair > CrossSectionHelper::CrossSection::GetCrossSectionCovarianceMatricies()
{
    std::map< std::string, CovarianceBiasPair> outputMap;

    for (const auto &[param, nUniverses] : m_inputData.m_systUniverseSizesMap)
    {
        const auto pair = this->GetCrossSectionCovarianceMatrix(param);
        outputMap.emplace(param, pair);
    }

    return outputMap;
}

// -----------------------------------------------------------------------------------------------------------------------------------------
                
std::map< std::string, CrossSectionHelper::CovarianceBiasPair > CrossSectionHelper::CrossSection::GetSmearingMatrixCovarianceMatricies()
{
    std::map< std::string, CovarianceBiasPair> outputMap;

    for (const auto &[param, nUniverses] : m_inputData.m_systUniverseSizesMap)
    {
        const auto pair = this->GetSmearingMatrixCovarianceMatrix(param);
        outputMap.emplace(param, pair);
    }

    return outputMap;
}

// -----------------------------------------------------------------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------------------------------------------------------------

std::vector<float> CrossSectionHelper::GetSystematicWeights(const Event::Truth &truth, const std::string &parameter)
{
    if (!truth.systParamNames.IsSet() || !truth.systParamFirstValueIndex.IsSet() || !truth.systParamValues.IsSet())
        throw std::invalid_argument("CrossSectionHelper::GetSystematicWeights - Systematic parameters are not set in the input truth information");

    const auto &systParamNames = truth.systParamNames();
    const auto &systParamFirstValueIndex = truth.systParamFirstValueIndex();
    const auto &systParamValues = truth.systParamValues();

    // Get the total number of parameters available
    if (systParamNames.size() != systParamFirstValueIndex.size())
        throw std::logic_error("CrossSectionHelper::GetSystematicWeights - Number of parameters is inconsistent in input truth information");
    
    const auto nParameters = systParamNames.size();

    // Find the requested parameter
    const auto iter = std::find(systParamNames.begin(), systParamNames.end(), parameter);
    if (iter == systParamNames.end())
        throw std::invalid_argument("CrossSectionHelper::GetSystematicWeights - Unknown parameter \"" + parameter + "\"");

    const unsigned int index = std::distance(systParamNames.begin(), iter);

    // Get the first and last index in the weights vector 
    const unsigned int firstValueIndex = systParamFirstValueIndex.at(index);
    const unsigned int lastValueIndex = ((index == nParameters - 1u) ? systParamValues.size() : systParamFirstValueIndex.at(index + 1u));

    if (firstValueIndex >= systParamValues.size() || lastValueIndex >= systParamValues.size())
        throw std::logic_error("CrossSectionHelper::GetSystematicWeights - Input systematic parameter first index vector is malformed");

    // Pick out the weights corresponding to the requested paramter
    std::vector<float> weights(std::next(systParamValues.begin(), firstValueIndex), std::next(systParamValues.begin(), lastValueIndex));
   
    // Ensure the the input weights are valid
    for (auto &weight : weights)
    {
        // If a weight is non-negative and non-infinite then it's fine!
        if (weight >= 0.f && std::abs(std::abs(weight) - std::numeric_limits<float>::max()) > std::numeric_limits<float>::epsilon() && std::isfinite(weight))
            continue;

        // Force negative weights to be zero and infinite weights to be unity
        weight = (weight < 0.f) ? 0.f : 1.f;
    }

    return weights;
}

// -----------------------------------------------------------------------------------------------------------------------------------------
        
void CrossSectionHelper::AddSystematicWeights(const Event::Truth &truth, const SystematicParamUniversesPairVector &params, SystematicWeightsMap &systWeightsMap)
{
    for (const auto &[name, nUniverses] : params)
    {
        const auto &weights = CrossSectionHelper::GetSystematicWeights(truth, name);
        
        if (weights.size() != nUniverses)
            throw std::logic_error("CrossSectionHelper::GetSystematicWeightsMap - Unexpected number of universes");

        if (!systWeightsMap.emplace(name, weights).second)
            throw std::invalid_argument("CrossSectionHelper::GetSystematicWeightsMap - Repeated systematic parameter: \"" + name + "\"");
    }
}

// -----------------------------------------------------------------------------------------------------------------------------------------
        
void CrossSectionHelper::AddMutuallyExclusiveSystematicWeights(const Event::Truth &truth, const MutuallyExclusiveParamVector &params, SystematicWeightsMap &systWeightsMap)
{
    for (const auto &[combinedName, nUniverses, paramsToCombine] : params)
    {
        // Build a SystematicParamUniversesPairVector for the parameters to combine
        SystematicParamUniversesPairVector paramUniverseVector;
        for (const auto &name : paramsToCombine)
            paramUniverseVector.emplace_back(name, nUniverses);
        
        // Get the weights for these parameters
        SystematicWeightsMap weightsMap;
        CrossSectionHelper::AddSystematicWeights(truth, paramUniverseVector, weightsMap);

        // Build the combined weights
        std::vector<float> combinedWeights;
        for (unsigned int iUni = 0u; iUni < nUniverses; ++iUni)
        {
            // Get the weights in this universe from each parameter that aren't exactly 1.f
            std::vector<float> variedWeights;
            std::vector<std::string> variedParams;
            for (const auto &name : paramsToCombine)
            {
                const auto &weight = weightsMap.at(name).at(iUni);

                if (std::abs(weight - 1.f) > std::numeric_limits<float>::epsilon())
                {
                    variedWeights.push_back(weight);
                    variedParams.push_back(name);
                }
            }

            // Check the parameters are truly mutually exclusive
            if (variedWeights.size() > 1)
            {
                std::string variedParamString = "";
                for (const auto &name : variedParams)
                    variedParamString += " " + name;

                throw std::logic_error("CrossSectionHelper::AddMutuallyExclusiveSystematicWeights - Input params are not mutually exclusive. \
                        Found a universe in which" + variedParamString + " each had values not equal to 1.f");
            }

            // Store the weight that was varied (or 1.f if there was no variation)
            combinedWeights.push_back(variedWeights.empty() ? 1.f : variedWeights.front());
        }
        
        // Add the combined weights to the output map
        if (!systWeightsMap.emplace(combinedName, combinedWeights).second)
            throw std::invalid_argument("CrossSectionHelper::AddMutuallyExclusiveSystematicWeights - Repeated systematic parameter: \"" + combinedName + "\"");
    }
}

// -----------------------------------------------------------------------------------------------------------------------------------------
        
void CrossSectionHelper::AddBootstrapWeights(const unsigned int nBootstrapUniverses, SystematicWeightsMap &systWeightsMap)
{
    // Add the entry in the input map for the bootstrap systematic parameter - intially with unit weights
    auto [iter, success] = systWeightsMap.emplace("bootstrap", std::vector<float>(nBootstrapUniverses, 1.f));

    if (!success)
        throw std::invalid_argument("CrossSectionHelper::AddBootstrapWeights - Bootstrap weights already exist in input map");

    // Populate the MC stats systematic weights using the Poisson bootstrap method
    std::poisson_distribution<int> poisson(1.f);
    auto &weights = iter->second;
    for (unsigned int iUni = 0; iUni < nBootstrapUniverses; ++iUni)
        weights.at(iUni) = static_cast<float>(poisson(m_generator));
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void CrossSectionHelper::AddUnitWeights(const SystematicParamUniversesPairVector &params, SystematicWeightsMap &systWeightsMap)
{
    for (const auto &[name, nUniverses] : params)
    {
        if (!systWeightsMap.emplace(name, std::vector<float>(nUniverses, 1.f)).second)
            throw std::logic_error("CrossSectionHelper::AddUnitWeights - Parameter \"" + name + "\" is already in the input map"); 
    }
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void CrossSectionHelper::AddUnitWeights(const MutuallyExclusiveParamVector &params, SystematicWeightsMap &systWeightsMap)
{
    for (const auto &[name, nUniverses, exclusiveNames] : params)
    {
        if (!systWeightsMap.emplace(name, std::vector<float>(nUniverses, 1.f)).second)
            throw std::logic_error("CrossSectionHelper::AddUnitWeights - Parameter \"" + name + "\" is already in the input map"); 
    }
}

// -----------------------------------------------------------------------------------------------------------------------------------------

bool CrossSectionHelper::CheckSystematicWeightsMapDimensions(const SystematicWeightsMap &systWeightsMap, const SystematicUniverseSizesMap &systUniverseSizesMap)
{
    if (systWeightsMap.size() != systUniverseSizesMap.size())
        return false;

    for (const auto &[name, nUniverses] : systUniverseSizesMap)
    {
        const auto iter = systWeightsMap.find(name);
        if (iter == systWeightsMap.end())
            return false;

        if (iter->second.size() != nUniverses)
            return false;
    }

    return true;
}

// -----------------------------------------------------------------------------------------------------------------------------------------
        
CrossSectionHelper::SystematicUniverseSizesMap CrossSectionHelper::GetSystematicUniverseSizesMap(const SystematicWeightsMap &systWeightsMap)
{
    // Sort the parameter names for reproducibility
    std::vector<std::string> names;
    for (const auto &[name, weights] : systWeightsMap)
        names.push_back(name);

    std::sort(names.begin(), names.end());

    // Store the number of universes for each parameter
    SystematicUniverseSizesMap outputMap;
    for (const auto &name : names)
    {
        const auto nUniverses = systWeightsMap.at(name).size();
        outputMap.emplace(name, nUniverses);
    }

    return outputMap;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void CrossSectionHelper::GetExtendedBinEdges(const float min, const float max, const std::vector<float> &binEdges, std::vector<float> &extendedBinEdges, bool &hasUnderflow, bool &hasOverflow)
{
    if (!extendedBinEdges.empty())
        throw std::invalid_argument("CrossSectionHelper::GetExtendedBinEdges - vector supplied to accept the output isn't empty");

    // Ensure the bin edges are sorted
    std::vector<float> sortedBinEdges(binEdges);
    std::sort(sortedBinEdges.begin(), sortedBinEdges.end());

    // Determine if underflow/overflow bins are required
    const auto firstBinEdge = binEdges.front();
    const auto lastBinEdge = binEdges.back();

    if (firstBinEdge < min)
        throw std::invalid_argument("CrossSectionHelper::GetExtendedBinEdges - lowest bin edge is less than the supplied minimum value");
    
    if (lastBinEdge > max)
        throw std::invalid_argument("CrossSectionHelper::GetExtendedBinEdges - highest bin edge is greater than the supplied maximum value");

    hasUnderflow = (firstBinEdge - min > std::numeric_limits<float>::epsilon());
    hasOverflow = (max - lastBinEdge > std::numeric_limits<float>::epsilon());

    // Extend the bins as required
    if (hasUnderflow)
        extendedBinEdges.push_back(min);
    
    extendedBinEdges.insert(extendedBinEdges.end(), binEdges.begin(), binEdges.end());

    if (hasOverflow)
        extendedBinEdges.push_back(max);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

std::shared_ptr<TH1F> CrossSectionHelper::GetFluxHistogram(const Config::Flux &flux)
{
    auto pHist = std::make_shared<TH1F>(("fluxHist_" + std::to_string(m_histCount++)).c_str(), "", flux.binEdges.size() - 1, flux.binEdges.data());
    const unsigned int nBins = pHist->GetNbinsX();

    if (flux.energyBins.size() != nBins)
        throw std::invalid_argument("CrossSectionHelper::GetFluxHistogram - Input number of input flux data points doesn't match the number of bins");

    for (unsigned int iBin = 1; iBin <= nBins; ++iBin)
    {
        // ATTN ROOT enumerates bins from 1, but we are reading data from a vector that enumerates from 0
        pHist->SetBinContent(iBin, flux.energyBins.at(iBin - 1));
    }

    return pHist;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

float CrossSectionHelper::GetTotalFlux(const std::shared_ptr<TH1F> &fluxHist)
{
    float total = 0.f;
    
    const unsigned int nBins = fluxHist->GetNbinsX();
    for (unsigned int iBin = 1; iBin <= nBins; ++iBin)
    {
        total += fluxHist->GetBinContent(iBin);
    }

    return total;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

std::pair<TVector2, TVector2> CrossSectionHelper::GetUncertaintyEigenVectors(const CovarianceBiasPair &covarianceBias, const unsigned int iBin, const unsigned int jBin)
{
    // Pick out the covariance matrix
    const auto &covarianceMatrix = covarianceBias.first;

    // Make sure we have square matricies of the same dimensions
    const auto &nBins = covarianceMatrix->GetNbinsX();
    if (covarianceMatrix->GetNbinsY() != nBins)
        throw std::invalid_argument("CrossSectionHelper::GetUncertaintyEigenVectors - Input covariance matrix is not square");

    // Make sure the requested bins are valid
    if (iBin == 0 || iBin > static_cast<unsigned int>(nBins))
    {
        throw std::range_error("CrossSectionHelper::GetUncertaintyEigenVectors - Supplied bin index iBin = " + std::to_string(iBin) +
            " is out of range of input matricies: 1 - " + std::to_string(nBins));
    }
    
    if (jBin == 0 || jBin > static_cast<unsigned int>(nBins))
    {
        throw std::range_error("CrossSectionHelper::GetUncertaintyEigenVectors - Supplied bin index jBin = " + std::to_string(jBin) +
            " is out of range of input matricies: 1 - " + std::to_string(nBins));
    }

    // Get the total error matrix elements for the bins we care about
    // Here we take the sum of the covariance and bias matrices as the total error
    const auto ii = covarianceMatrix->GetBinContent(iBin, iBin);
    const auto ij = covarianceMatrix->GetBinContent(iBin, jBin);
    const auto ji = covarianceMatrix->GetBinContent(jBin, iBin);
    const auto jj = covarianceMatrix->GetBinContent(jBin, jBin);

    // Get the discriminant of the eigen equation for this 2D covariance matrix
    const auto trace = ii + jj;
    const auto det = ii*jj - ij*ji;
    const auto disc = trace*trace - 4*det;

    if (disc < 0.f)
        throw std::logic_error("CrossSectionHelper::GetUncertaintyEigenVectors - Eigen equation has no solutions");

    // Get the eigenvalues
    const auto rootDisc = std::pow(disc, 0.5f);
    const auto lambdaMax = std::max(trace + rootDisc, trace - rootDisc) / 2.;
    const auto lambdaMin = std::min(trace + rootDisc, trace - rootDisc) / 2.;

    // Get the unit eigenvectors
    const auto vectMaxUnit = TVector2(ij, lambdaMax - ii).Unit();
    const auto vectMinUnit = TVector2(ij, lambdaMin - ii).Unit();

    // Get the eigen vectors which have length of the root of the eigen values
    const auto sigmaMax = std::pow(lambdaMax, 0.5f);
    const auto sigmaMin = std::pow(lambdaMin, 0.5f);

    const auto vectMax = sigmaMax * vectMaxUnit;
    const auto vectMin = sigmaMin * vectMinUnit;

    return std::pair<TVector2, TVector2>(vectMax, vectMin);
}

// -----------------------------------------------------------------------------------------------------------------------------------------
        
std::shared_ptr<TH2F> CrossSectionHelper::GetTotalCovarianceMatrix(const std::map< std::string, CovarianceBiasPair > &covarianceBiasPairs)
{
    std::shared_ptr<TH2F> totalCovarianceMatrix;

    for (const auto &[paramName, covarianceBiasPair] : covarianceBiasPairs)
    {
        const auto &covarianceMatrix = covarianceBiasPair.first;
        const auto &biasVector = covarianceBiasPair.second;

        // Check the binning is sensisble
        const unsigned int nBins = covarianceMatrix->GetNbinsX();
        if (static_cast<unsigned int>(covarianceMatrix->GetNbinsY()) != nBins)
            throw std::invalid_argument("CrossSectionHelper::GetTotalCovarianceMatrix - Input covariance matrix is not square");
        
        if (static_cast<unsigned int>(biasVector->GetNbinsX()) != nBins)
            throw std::invalid_argument("CrossSectionHelper::GetTotalCovarianceMatrix - Input bias vector has different dimensions to covariance matrix");

        // If this is our first entry, then setup the binning of the output covariance matrix
        if (!totalCovarianceMatrix)
        {
            totalCovarianceMatrix = std::make_shared<TH2F>(("totalCovariance_" + std::to_string(m_histCount++)).c_str(), "", nBins, 0, nBins, nBins, 0, nBins);

            // Initialize the total covariance matrix with zeros
            for (unsigned int iBin = 1; iBin <= nBins; ++iBin)
            {
                for (unsigned int jBin = 1; jBin <= nBins; ++jBin)
                {
                    totalCovarianceMatrix->SetBinContent(iBin, jBin, 0.f);
                }
            }
        }
        // Otherwise make sure that we have the expected number of bins
        else if (static_cast<unsigned int>(totalCovarianceMatrix->GetNbinsX()) != nBins)
        {
            throw std::invalid_argument("CrossSectionHelper::GetTotalCovarianceMatrix - Sizes of input covaraince matrices are not equal");
        }

        // Add to the total covariance matrix
        for (unsigned int iBin = 1u; iBin <= nBins; ++iBin)
        {
            for (unsigned int jBin = 1u; jBin <= nBins; ++jBin)
            {
                const auto covariance = covarianceMatrix->GetBinContent(iBin, jBin);
                const auto bias = biasVector->GetBinContent(iBin) * biasVector->GetBinContent(jBin);

                // Add the covariance and bias to the current total
                const auto currentValue = totalCovarianceMatrix->GetBinContent(iBin, jBin);
                const auto updatedValue = currentValue + covariance + bias;
                totalCovarianceMatrix->SetBinContent(iBin, jBin, updatedValue);
            }
        }
    }

    return totalCovarianceMatrix;
}

} // namespace ubcc1pi
