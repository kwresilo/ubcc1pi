/**
 *  @file  ubcc1pi_standalone/Helpers/CrossSectionHelper.cxx
 *
 *  @brief The implementation of the cross section helper class
 */

#include "ubcc1pi_standalone/Helpers/CrossSectionHelper.h"

#include <stdexcept>

namespace ubcc1pi
{

CrossSectionHelper::CrossSection::CrossSection(const std::vector<float> &binEdges, const bool hasUnderflow, const bool hasOverflow, const bool scaleByBinWidth, const SystematicUniverseSizesMap &systUniverseSizesMap) :
    m_binEdges(binEdges),
    m_hasUnderflow(hasUnderflow),
    m_hasOverflow(hasOverflow),
    m_scaleByBinWidth(scaleByBinWidth),
    m_systUniverseSizesMap(systUniverseSizesMap),
    m_dataSelectedReco(this->GetEmptyHist1D()),
    m_signalSelectedNomRecoTrue(this->GetEmptyHist2D()),
    m_signalAllNomTrue(this->GetEmptyHist1D()),
    m_backgroundSelectedNomReco(this->GetEmptyHist1D())
{
    // Setup the histograms for the universe variations
    for (const auto &[name, nUniverses] : m_systUniverseSizesMap)
    {
        auto &signalSelectedRecoTrueVect = m_signalSelectedRecoTrue[name];
        auto &signalAllTrueVect = m_signalAllTrue[name];
        auto &backgroundSelectedRecoVect = m_backgroundSelectedReco[name];

        for (unsigned int iUni = 0; iUni < nUniverses; ++iUni)
        {
            signalSelectedRecoTrueVect.push_back(this->GetEmptyHist2D());
            signalAllTrueVect.push_back(this->GetEmptyHist1D());
            backgroundSelectedRecoVect.push_back(this->GetEmptyHist1D());
        }
    }
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void CrossSectionHelper::CrossSection::AddSignalEvent(const float trueValue, const float recoValue, const bool isSelected, const float nominalWeight, const SystematicWeightsMap &systWeightsMap)
{
    // Check the validity of the input map if systematic weights
    if (!CrossSectionHelper::CheckSystematicWeightsMapDimensions(systWeightsMap, m_systUniverseSizesMap))
        throw std::invalid_argument("CrossSection::AddSignalEvent - Input systematic weights map has the wrong dimensions!");

    // Add this event in the nominal universe
    m_signalAllNomTrue->Fill(trueValue, nominalWeight);

    if (isSelected)
        m_signalSelectedNomRecoTrue->Fill(trueValue, recoValue, nominalWeight);

    // Add the event in each systematic universe
    for (const auto &[name, nUniverses] : m_systUniverseSizesMap)
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
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void CrossSectionHelper::CrossSection::AddSelectedBackgroundEvent(const float recoValue, const float nominalWeight, const SystematicWeightsMap &systWeightsMap)
{
    // Check the validity of the input map if systematic weights
    if (!CrossSectionHelper::CheckSystematicWeightsMapDimensions(systWeightsMap, m_systUniverseSizesMap))
        throw std::invalid_argument("CrossSection::AddSignalEvent - Input systematic weights map has the wrong dimensions!");

    // Add this event in the nominal universe
    m_backgroundSelectedNomReco->Fill(recoValue, nominalWeight);

    // Add the event in each systematic universe
    for (const auto &[name, nUniverses] : m_systUniverseSizesMap)
    {
        auto &backgroundSelectedRecoVect = m_backgroundSelectedReco.at(name);

        const auto &weights = systWeightsMap.at(name);
        for (unsigned int iUni = 0; iUni < nUniverses; ++iUni)
        {
            const auto universeWeight = weights.at(iUni);
            backgroundSelectedRecoVect.at(iUni)->Fill(recoValue, nominalWeight * universeWeight);
        }
    }
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void CrossSectionHelper::CrossSection::AddSelectedBNBDataEvent(const float recoValue)
{
    m_dataSelectedReco->Fill(recoValue);
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

    return ((binIndex == 0 && m_hasUnderflow) || (binIndex == nBins && m_hasOverflow));
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

        // Normalise the bins such that the columns add to unity
        if (summedValue <= std::numeric_limits<float>::epsilon())
            throw std::logic_error("CrossSection::GetSmearingMatrix - Found true bin with no selected signal events");

        const auto norm = 1.f / summedValue;

        for (unsigned int iReco = 1; iReco <= nBins; ++iReco)
        {
            const auto value = signalSelectedRecoTrue->GetBinContent(iTrue, iReco) * norm;
            smearingMatrix->SetBinContent(iTrue, iReco, value);
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
    // Check the input objects have the correct dimensions
    const auto nBins = m_binEdges.size() - 1;

    if (static_cast<unsigned int>(signalAllTrue->GetXaxis()->GetNbins()) != nBins)
        throw std::logic_error("CrossSection::GetSmearedEfficiency - input object signalAllTrue has the wrong dimensions");

    if (static_cast<unsigned int>(signalSelectedRecoTrue->GetXaxis()->GetNbins()) != nBins ||
        static_cast<unsigned int>(signalSelectedRecoTrue->GetYaxis()->GetNbins()) != nBins)
    {
        throw std::logic_error("CrossSection::GetSmearedEfficiency - input matrix has the wrong dimensions");
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
    const auto smearingMatrix = this->GetSmearingMatrix(signalSelectedRecoTrue);
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

        if (denominator <= std::numeric_limits<float>::epsilon())
            throw std::logic_error("CrossSection::GetSmearedEfficiency - denominator of smeared efficiency is zero!");

        const auto value = numerator / denominator;
        effSmearedReco->SetBinContent(iReco, value);
    }

    return effSmearedReco;
}

// -----------------------------------------------------------------------------------------------------------------------------------------
                
std::shared_ptr<TH1F> CrossSectionHelper::CrossSection::GetCrossSection(const std::shared_ptr<TH1F> &selected, const std::shared_ptr<TH1F> &background, const std::shared_ptr<TH1F> &smearedEff) const
{
    // TODO add in the flux, exposure and number of targets as configurable parameters ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    const float flux = 7.3f;       // e-10 / cm^2 / POT
    const float exposure = 1.455f; // e+20 POT
    const float targets = 4.1354f; // e+31 nucleons

    const auto normFactor = 1.f / (flux * exposure * targets);
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    // Check the input objects have the correct dimensions
    const auto nBins = m_binEdges.size() - 1;

    if (static_cast<unsigned int>(selected->GetXaxis()->GetNbins()) != nBins)
        throw std::logic_error("CrossSection::GetSmearedEfficiency - input selected distribution has the wrong dimensions");

    if (static_cast<unsigned int>(background->GetXaxis()->GetNbins()) != nBins)
        throw std::logic_error("CrossSection::GetSmearedEfficiency - input background distribution has the wrong dimensions");

    if (static_cast<unsigned int>(smearedEff->GetXaxis()->GetNbins()) != nBins)
        throw std::logic_error("CrossSection::GetSmearedEfficiency - input smeared efficiency distribution has the wrong dimensions");

    auto crossSection = this->GetEmptyHist1D();
    for (unsigned int iReco = 1; iReco <= nBins; ++iReco)
    {
        // Skip under/overflow bins
        if (this->IsUnderOverflowBin(iReco))
            continue;

        const auto binWidth = m_scaleByBinWidth ? (m_binEdges.at(iReco) - m_binEdges.at(iReco - 1)) : 1.f;
        const auto value = normFactor * (selected->GetBinContent(iReco) - background->GetBinContent(iReco)) / (smearedEff->GetBinContent(iReco) * binWidth);

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
    return this->GetCrossSection(m_dataSelectedReco, m_backgroundSelectedNomReco, smearedEff);
}

// -----------------------------------------------------------------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------------------------------------------------------------

CrossSectionHelper::SystematicWeightsMap CrossSectionHelper::GetSystematicWeightsMap(const Event::Truth &truth)
{
    if (!truth.systParamNames.IsSet() || !truth.systParamFirstValueIndex.IsSet() || !truth.systParamValues.IsSet())
        throw std::invalid_argument("CrossSectionHelper::GetSystematicWeightsMap - Systematic parameters are not set in the input truth information");

    const auto &systParamNames = truth.systParamNames();
    const auto &systParamFirstValueIndex = truth.systParamFirstValueIndex();
    const auto &systParamValues = truth.systParamValues();

    // Make the mapping of systematic weights
    if (systParamNames.size() != systParamFirstValueIndex.size())
        throw std::logic_error("CrossSectionHelper::GetSystematicWeightsMap - Number of parameters is inconsistent in input truth information");

    SystematicWeightsMap outputMap;

    // Add in the systematic weights from the input truth info
    const auto nParameters = systParamNames.size();
    for (unsigned int iParam = 0; iParam < nParameters; ++iParam)
    {
        const auto &name = systParamNames.at(iParam);

        // Get the weights
        const auto firstValueIndex = systParamFirstValueIndex.at(iParam);
        const auto lastValueIndex = ((iParam == nParameters - 1) ? systParamValues.size() : systParamFirstValueIndex.at(iParam + 1));
        const std::vector<float> values(std::next(systParamValues.begin(), firstValueIndex), std::next(systParamValues.begin(), lastValueIndex));

        // ATTN The nominal event weights are stored in this object too and have one value - we don't need these
        const auto nValues = values.size();
        if (nValues < 2)
            continue;
        
        // Looks good, so store the weights
        if (!outputMap.emplace(name, values).second)
            throw std::invalid_argument("CrossSectionHelper::GetSystematicWeightsMap - Found repeated parameter: \"" + name + "\"");
    }

    return outputMap;
}

// -----------------------------------------------------------------------------------------------------------------------------------------
        
void CrossSectionHelper::AddBootstrapWeights(const unsigned int nBootstrapUniverses, SystematicWeightsMap &systWeightsMap)
{
    // TODO actually add the weights
    systWeightsMap.emplace("bootstrap", std::vector<float>(nBootstrapUniverses, 1.f));
}

// -----------------------------------------------------------------------------------------------------------------------------------------
        
CrossSectionHelper::SystematicWeightsMap CrossSectionHelper::GetUnitSystematicWeightsMap(const SystematicUniverseSizesMap &systUniverseSizesMap)
{
    SystematicWeightsMap outputMap;
    for (const auto &[name, nUniverses] : systUniverseSizesMap)
        outputMap.emplace(name, std::vector<float>(nUniverses, 1.f));

    return outputMap;
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

} // namespace ubcc1pi
