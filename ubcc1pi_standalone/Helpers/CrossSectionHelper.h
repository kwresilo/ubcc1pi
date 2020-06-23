/**
 *  @file  ubcc1pi_standalone/Helpers/CrossSectionHelper.h
 *
 *  @brief The header file for the cross section helper class
 */

#ifndef UBCC1PI_STANDALONE_HELPERS_CROSS_SECTION_HELPER
#define UBCC1PI_STANDALONE_HELPERS_CROSS_SECTION_HELPER

#include "ubcc1pi_standalone/Helpers/AnalysisHelper.h"
#include "ubcc1pi_standalone/Helpers/FormattingHelper.h"
#include "ubcc1pi_standalone/Helpers/PlottingHelper.h"

#include <stdexcept>
#include <string>
#include <vector>

#include <TFile.h>
#include <TTree.h>
#include <TH2F.h>
#include <TLine.h>

namespace ubcc1pi
{

/**
 *  @brief  The cross section helper class
 */
class CrossSectionHelper
{
    public:
        
        /**
         *  @brief  The cross-section object
         */
        class XSec
        {
            public:
                /**
                 *  @brief  Constructor
                 *
                 *  @param  fileName the name of the output file
                 *  @param  min the minimum possible value that the quantity can take
                 *  @param  max the maximum possible value that the quantity can take
                 *  @param  useAbsPdg use absolute PDG code when classifying event types
                 *  @param  countProtonsInclusively count protons inclusively when classifying event types
                 */
                XSec(const std::string &fileName, const float &min, const float &max, const bool useAbsPdg, const bool countProtonsInclusively);

                /**
                 *  @brief  Destructor
                 */
                ~XSec();

                /**
                 *  @brief  Manually set the bins
                 *
                 *  @param  binEdges the bin edges - has size nBins + 1
                 */
                void SetBins(const std::vector<float> &binEdges);

                /**
                 *  @brief  Automatically choose the bins at sensible values
                 *
                 *  @param  lower the lower bin edge of the first bin
                 *  @param  upper the upper bin edge of the last bin
                 *  @param  targetFracError the target fractional uncertainty in each bin
                 *  @param  targetSmearingDiagonal the target entries on the diagonal of the smearing matrix
                 */
                void SetBinsAuto(const float lower, const float upper, const float targetFracError, const float targetSmearingDiagonal);

                /**
                 *  @brief  Add a signal event 
                 *
                 *  @param  pEvent the event
                 *  @param  isSelected if the event has been selected
                 *  @param  trueValue the true value of the quantity to measure
                 *  @param  recoValue the reco value of the quantity to meaasure
                 *  @param  weight the weight
                 *  @param  normalisation the sample normalisation
                 */
                void AddSignalEvent(const std::shared_ptr<Event> &pEvent, const bool &isSelected, const float &trueValue, const float &recoValue, const float weight, const float normalisation);

                /**
                 *  @brief  Add a background event that has been selected
                 *
                 *  @param  pEvent the event
                 *  @param  sampleType the sample type (Overlays, EXT, Dirt)
                 *  @param  recoValue the reco value of the quantity to meaasure
                 *  @param  weight the weight
                 *  @param  normalisation the sample normalisation
                 */
                void AddSelectedBackgroundEvent(const std::shared_ptr<Event> &pEvent, const AnalysisHelper::SampleType sampleType, const float &recoValue, const float weight, const float normalisation);
                
                /**
                 *  @brief  Add a BNB data event that has been selected
                 *
                 *  @param  pEvent the event
                 *  @param  recoValue the reco value of the quantity to meaasure
                 */
                void AddSelectedBNBDataEvent(const std::shared_ptr<Event> &pEvent, const float &recoValue);

                /**
                 *  @brief  Get the confusion matrix, M[reco][true]
                 *
                 *  @param  matrix the output confusion matrix
                 *  @param  uncertainties the output uncertainties on the confusion matrix elements
                 */
                void GetConfusionMatrix(std::vector< std::vector<float> > &matrix, std::vector< std::vector<float> > &uncertainties) const;
                
                /**
                 *  @brief  Get the smearing matrix, S[reco][true] = P(reco | true)
                 *
                 *  @param  matrix the output smearing matrix
                 *  @param  uncertainties the output uncertainties on the smearing matrix elements
                 */
                void GetSmearingMatrix(std::vector< std::vector<float> > &matrix, std::vector< std::vector<float> > &uncertainties) const;
                
                /**
                 *  @brief  Get the selection efficiencies in the true bins
                 *
                 *  @param  efficiencies the output efficiency per bin
                 *  @param  uncertainties the output uncertainty on the efficiency per bin
                 */
                void GetTrueEfficiencies(std::vector<float> &efficiencies, std::vector<float> &uncertainties) const;

                /**
                 *  @brief  Get the smeared selection efficiencies
                 *
                 *  @param  efficiencies the output smeared efficiency per bin
                 *  @param  uncertainties the output uncertainty on the efficiency per bin
                 */
                void GetSmearedEfficiencies(std::vector<float> &efficiencies, std::vector<float> &uncertainties) const;
                
                /**
                 *  @brief  Get the bin widths
                 *
                 *  @return the bin widths
                 */
                std::vector<float> GetBinWidths() const;

                /**
                 *  @brief  Get the cross-section from MC using true values
                 *
                 *  @param  bins the output cross-section values
                 *  @param  uncertainties the uncertainty on the bins
                 */
                void GetMCTrueCrossSection(std::vector<float> &bins, std::vector<float> &uncertainties) const;

                /**
                 *  @brief  Get the cross-section from MC using smeared true values
                 *
                 *  @param  bins the output cross-section values
                 *  @param  uncertainties the uncertainty on the bins
                 */
                void GetSmearedMCTrueCrossSection(std::vector<float> &bins, std::vector<float> &uncertainties) const;

                /**
                 *  @brief  Get the reco cross-section using either real BNB data or Overlays + EXT + Dirt
                 *
                 *  @param  useRealData if true, use real BNB data, if false use Overlays + EXT + Dirt
                 *  @param  bins the output cross-section values
                 *  @param  uncertainties the uncertainty on the bins
                 */
                void GetRecoCrossSection(const bool useRealData, std::vector<float> &bins, std::vector<float> &uncertainties) const;

                /**
                 *  @brief  Print the contents of the bins
                 *
                 *  @param  outputFileNamePrefix the prefix of the name of the output file we should save the tables to
                 */
                void PrintBinContents(const std::string &outputFileNamePrefix) const;
                
                /**
                 *  @brief  Make the plots
                 *
                 *  @param  fileNamePrefix the file name prefix
                 */
                void MakePlots(const std::string &fileNamePrefix) const;

            private:

                /**
                 *  @brief  Fill the output tree with an event
                 *
                 *  @param  sampleType the sample type
                 *  @param  sampleNorm the sample normalisation
                 *  @param  weight the event weight
                 *  @param  recoValue the reco value
                 *  @param  trueValue the true value
                 *  @param  isSignal if the event is signal
                 *  @param  isSelected if the event is selected
                 *  @param  classification the event classification
                 */
                void FillEvent(const AnalysisHelper::SampleType sampleType, const float sampleNorm, const float weight, const float &recoValue, const float &trueValue, const bool isSignal, const bool isSelected, const std::string &classification);

                /**
                 *  @brief  Get the initial bins within a given range with some minimum number of MC selected signal events per bin
                 *
                 *  @param  lower the lower bound of the bins
                 *  @param  upper the upper bound of the bins
                 *  @param  minEventsPerBin the minimum number of selected signal MC events per bin
                 *
                 *  @return the initial bins
                 */
                std::vector<float> GetInitialBins(const float lower, const float upper, const unsigned int minEventsPerBin) const;

                /**
                 *  @brief  Remove the bin with the largest statistical uncertainty
                 *
                 *  @param  targetFracError the target fractional statistical error, we wont remove a bin if the uncertainty is below this target value
                 *  @param  binEdges the bin edges to modify
                 *
                 *  @return boolean, true if a bin was removed - false otherwise
                 */
                bool RemoveLowestStatsBin(const float targetFracError, std::vector<float> &binEdges);

                /**
                 *  @brief  Get the factional uncertainties on the reconstructed cross-section 
                 *
                 *  @param  binEdges the input bin edges to use
                 *
                 *  @return the fractional uncertainty on the reconstructed cross-section in each bin
                 */
                std::vector<float> GetFractionalRecoCrossSectionErrors(const std::vector<float> &binEdges);

                /**
                 *  @brief  Remove the bin with the smallest diagonal entry of the smearing matrix
                 *
                 *  @param  targetSmearingDiagonal the target value of the diagonals of the smearing matrix, we wont remove a bin if the the diagonal is above this value
                 *  @param  binEdges the bin edges to modify
                 *
                 *  @return boolean, true if a bin was removed - false otherwise
                 */
                bool RemoveMostSmearedBin(const float targetSmearingDiagonal, std::vector<float> &binEdges);
                
                /**
                 *  @brief  Get the differential cross-section in a bin given some parameters
                 *
                 *  @param  nTotal the total number of events selected
                 *  @param  nBackgrounds the expected number of background events
                 *  @param  efficiency the efficiency of the selection
                 *  @param  binWidth the width of the bin
                 *
                 *  @return the cross-section value
                 */
                float GetCrossSection(const float nTotal, const float nBackgrounds, const float efficiency, const float binWidth) const;

                /**
                 *  @brief  Get the uncertainty on the differential cross-section in a bin
                 *
                 *  @param  nTotal the total number of events selected
                 *  @param  totalErr the error on the total count
                 *  @param  nBackgrounds the expected number of background events
                 *  @param  backgroundsErr the error on the background count
                 *  @param  efficiency the efficiency of the selection
                 *  @param  efficiencyErr the error on the efficiency
                 *  @param  binWidth the bin width
                 *
                 *  @return the error on the cross-section value
                 */
                float GetCrossSectionUncertainty(const float nTotal, const float totalErr, const float nBackgrounds, const float backgroundsErr, const float efficiency, const float efficiencyErr, const float binWidth) const;
                

                /**
                 *  @brief  Get the total number of signal events in each true bin
                 *
                 *  @param  mustBeSelected if we should insist that the events are selected
                 *
                 *  @return the true signal per bin
                 */
                std::vector<float> GetTrueSignalPerBin(const bool mustBeSelected) const;
                
                /**
                 *  @brief  Get the overlay normalisation 
                 *
                 *  @return the normalisation factor for overlays
                 */
                float GetOverlayNormalisation() const;
                
                /**
                 *  @brief  Get the selected BNB data in each bin
                 *
                 *  @param  bins the values per bin
                 *  @param  uncertainties the uncertainties
                 *  @param  useRealData if true we will use real BNB data, if false we use Overlays + EXT + Dirt as fake data
                 */
                void GetBNBDataPerBin(std::vector<float> &bins, std::vector<float> &uncertainties, const bool useRealData = true) const;

                /**
                 *  @brief  Get the expected number of backgrounds in each bin
                 *
                 *  @param  bins the values per bin
                 *  @param  uncertainties the uncertainties
                 */
                void GetBackgroundsPerBin(std::vector<float> &bins, std::vector<float> &uncertainties) const;

                /**
                 *  @brief  Get the index of a bin from the input value
                 *
                 *  @param  value the input value
                 *
                 *  @return the bin index
                 */
                unsigned int GetBinIndex(const float &value) const;

                /**
                 *  @brief  Get the name of the bin with the given index (can be numerical or underflow/overflow)
                 *
                 *  @param  index the bin index
                 *
                 *  @return the bin name
                 */
                std::string GetBinName(const unsigned int index) const;

                /**
                 *  @brief  Determine if a given bin is under or overflow
                 *
                 *  @param  index the bin index
                 *
                 *  @return boolean, true if underflow or overflow bin
                 */
                bool IsUnderOverFlowBin(const unsigned int index) const;
                
                /**
                 *  @brief  Smear the input bins
                 *
                 *  @param  inputBins the input bin values
                 *  @param  inputUncertainties the uncertainties on the bin values
                 *  @param  bins the output smeared bins
                 *  @param  uncertainties the output smeared uncertainties
                 */
                void SmearBins(const std::vector<float> &inputBins, const std::vector<float> &inputUncertainties, std::vector<float> &bins, std::vector<float> &uncertainties) const;

                /**
                 *  @brief  Print a matrix as a table
                 *
                 *  @param  matrix the matrix to print
                 *  @param  uncertainties the uncertainties on the values in the matrix
                 *  @param  outputFileName the output file name for the table
                 */
                void PrintMatrix(const std::vector< std::vector<float> > &matrix, const std::vector< std::vector<float> > &uncertainties, const std::string &outputFileName) const;

                /**
                 *  @brief  The output event structure
                 */
                struct OutputEvent
                {
                    int          m_sampleType;      ///< The sample type enumeration
                    float        m_sampleNorm;      ///< The sample normalisation factor   
                    float        m_weight;          ///< The event weight
                    float        m_recoValue;       ///< The reconstructed value
                    float        m_trueValue;       ///< The true value
                    bool         m_isSignal;        ///< If the event is signal
                    bool         m_isSelected;      ///< If the event is selected
                    std::string  m_classification;  ///< The event classification string
                };

                std::string    m_fileName; ///< The output file name
                TFile         *m_pFile;    ///< The output file
                TTree         *m_pTree;    ///< The output tree (we store in a tree so we don't use too much memory)

                float          m_min;      ///< The minimum possible value of the quantity
                float          m_max;      ///< The maximum possible value of the quantity

                std::vector<float> m_binEdges;     ///< The bin edges
                bool               m_hasUnderflow; ///< If we have an underflow bin
                bool               m_hasOverflow;  ///< If we have an overflow bin

                bool           m_useAbsPdg;                ///< Use absolute PDGs when classifying events
                bool           m_countProtonsInclusively;  ///< Count protons inclusively when classifying events

                float          m_nTargets;        ///< The number of target particles / 10^31
                float          m_integratedFlux;  ///< The integrated flux / 10^11 cm^-2

                OutputEvent    m_outputEvent;  ///< The output event struture to bind to the trees
            
        };
};

// -----------------------------------------------------------------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------------------------------------------------------------

CrossSectionHelper::XSec::XSec(const std::string &fileName, const float &min, const float &max, const bool useAbsPdg, const bool countProtonsInclusively) : 
    m_fileName(fileName),
    m_pFile(new TFile(m_fileName.c_str(), "RECREATE")),
    m_pTree(new TTree("events", "events")),
    m_min(min),
    m_max(max),
    m_hasUnderflow(false),
    m_hasOverflow(false),
    m_useAbsPdg(useAbsPdg),
    m_countProtonsInclusively(countProtonsInclusively),
    m_nTargets(4.1741), // TODO make configurable
    m_integratedFlux(1.26816) // TODO make configurable
{
    // Setup the output branches
    m_pTree->Branch("sampleType", &m_outputEvent.m_sampleType);
    m_pTree->Branch("sampleNorm", &m_outputEvent.m_sampleNorm);
    m_pTree->Branch("weight", &m_outputEvent.m_weight);
    m_pTree->Branch("recoValue", &m_outputEvent.m_recoValue);
    m_pTree->Branch("trueValue", &m_outputEvent.m_trueValue);
    m_pTree->Branch("isSignal", &m_outputEvent.m_isSignal);
    m_pTree->Branch("isSelected", &m_outputEvent.m_isSelected);
    m_pTree->Branch("classification", &m_outputEvent.m_classification);
}

// -----------------------------------------------------------------------------------------------------------------------------------------
                
CrossSectionHelper::XSec::~XSec()
{
    m_pFile->Write();
    m_pFile->Close();
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void CrossSectionHelper::XSec::SetBins(const std::vector<float> &binEdges)
{
    if (binEdges.size() < 2)
        throw std::invalid_argument("XSec::SetBins - There must be at least two bin edges");

    auto sortedBinEdges = binEdges;
    std::sort(sortedBinEdges.begin(), sortedBinEdges.end());

    if (sortedBinEdges.front() < m_min)
        throw std::invalid_argument("XSec::SetBins - Lowest bin edge is out if bounds");
    
    if (sortedBinEdges.back() > m_max)
        throw std::invalid_argument("XSec::SetBins - Highest bin edge is out if bounds");
   
    // Remove any binning we currently have
    m_binEdges.clear();

    // If required add underflow bin
    m_hasUnderflow = sortedBinEdges.front() > m_min;

    if (m_hasUnderflow)
        m_binEdges.push_back(m_min);

    // Add the main bins
    m_binEdges.insert(m_binEdges.end(), sortedBinEdges.begin(), sortedBinEdges.end());

    // If required add an overflow bin
    m_hasOverflow = sortedBinEdges.back() < m_max;

    if (m_hasOverflow)
        m_binEdges.push_back(m_max);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

std::vector<float> CrossSectionHelper::XSec::GetInitialBins(const float lower, const float upper, const unsigned int minEventsPerBin) const
{
    if (upper <= lower)
        throw std::invalid_argument("XSec::GetInitialBins - Invalid range: " + std::to_string(lower) + " -> " + std::to_string(upper));
    
    // Extract the reco values from the tree
    std::vector<float> recoValues;
    for (unsigned int i = 0; i < m_pTree->GetEntries(); ++i)
    {
        m_pTree->GetEntry(i);

        // We only care about overlays
        if (m_outputEvent.m_sampleType != AnalysisHelper::Overlay)
            continue;
        
        // Only use selected signal events
        if (!m_outputEvent.m_isSelected || !m_outputEvent.m_isSignal)
            continue;

        // Only use events with reconstructed values in the range we care about
        if (m_outputEvent.m_recoValue < lower || m_outputEvent.m_recoValue > upper)
            continue;

        recoValues.push_back(m_outputEvent.m_recoValue);
    }

    // Sort the reco values
    std::sort(recoValues.begin(), recoValues.end()); 

    // Work out how many bins we need to have the required number of events per bin
    const auto nBins = static_cast<unsigned int>(std::floor(static_cast<float>(recoValues.size()) / static_cast<float>(minEventsPerBin)));

    // Work out the average number of events we need to put in each of those bins
    const auto nEventsPerBin = static_cast<unsigned int>(std::floor(static_cast<float>(recoValues.size()) / static_cast<float>(nBins)));

    // Set the bin edges
    std::vector<float> binEdges;
    binEdges.push_back(lower);

    unsigned int recoValueIndex = nEventsPerBin - 1u;
    while (recoValueIndex < recoValues.size() - 1u)
    {
        const auto edge = recoValues.at(recoValueIndex);

        // This check in required in case we have loads of entries at the same floating point value (e.g. all zero)
        if (edge > binEdges.back())
        {
            binEdges.push_back(edge);
            recoValueIndex += nEventsPerBin;
            continue;
        }
        
        recoValueIndex++;
    }

    // Remove the last bin
    binEdges.erase(std::next(binEdges.end(), -1));
    
    if (upper > binEdges.back())
        binEdges.push_back(upper);

    return binEdges;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

bool CrossSectionHelper::XSec::RemoveMostSmearedBin(const float targetSmearingDiagonal, std::vector<float> &binEdges)
{
    // Use the current binning
    this->SetBins(binEdges);
        
    // Get the smearing matrix
    std::vector< std::vector<float> > smearingMatrix, smearingMatrixErr;
    this->GetSmearingMatrix(smearingMatrix, smearingMatrixErr);

    const auto nBins = smearingMatrix.size();
    
    // Find the bin with the smallest diagonal smearing matrix element
    float smallestSmearingElement = std::numeric_limits<float>::max();
    unsigned int mostSmearedBinIndex = 0u;

    // Loop over the bins
    for (unsigned int i = 0; i < nBins; ++i)
    {
        // Don't care about underflow or overflow bins
        if (this->IsUnderOverFlowBin(i))
            continue;

        // Get the diagonal entry corresponding to this bin
        const auto smearingElement = smearingMatrix.at(i).at(i);

        if (smearingElement < smallestSmearingElement)
        {
            smallestSmearingElement = smearingElement;
            mostSmearedBinIndex = i;
        }
    }

    // Check to see if the smearing is acceptable
    if (smallestSmearingElement > targetSmearingDiagonal)
        return false;

    // We have a bin with lots of smearing, so let's merge it with an adjacent bin to help the problem
    // We need to work out which edge of the bin we should remove
    bool removeLeftEdge = true;

    if (mostSmearedBinIndex - (m_hasUnderflow ? 1u : 0u) == 0u) // First bin
    {
        removeLeftEdge = false;
    }
    else if (mostSmearedBinIndex == nBins - 1u - (m_hasOverflow ? 1u : 0u)) // Last bin
    {
        removeLeftEdge = true;
    }
    else // A middle bin
    {
        // Merge with the adjacent bin with the smaller smearing element
        const auto smearingElementLeft = smearingMatrix.at(mostSmearedBinIndex - 1).at(mostSmearedBinIndex - 1);
        const auto smearingElementRight = smearingMatrix.at(mostSmearedBinIndex + 1).at(mostSmearedBinIndex + 1);

        removeLeftEdge = (smearingElementLeft < smearingElementRight);
    }

    // Remove the desired edge of the bad bin
    binEdges.erase(std::next(binEdges.begin(), mostSmearedBinIndex - (m_hasUnderflow ? 1u : 0u) + (removeLeftEdge ? 0u : 1u)));

    // We made a change
    return true;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

std::vector<float> CrossSectionHelper::XSec::GetFractionalRecoCrossSectionErrors(const std::vector<float> &binEdges)
{
    // Use the current binning
    this->SetBins(binEdges);

    // Get the cross-sections from MC (as fake data)
    std::vector<float> xSec, xSecErr;
    this->GetRecoCrossSection(false, xSec, xSecErr);
    const auto nBins = xSec.size();

    // Get the fractional errors in each bin
    std::vector<float> fracErrs;
    for (unsigned int i = 0; i < nBins; ++i)
    {
        // If there are no entries for the cross-section in a given bin then give it an "infinite" error
        if (xSec.at(i) <= std::numeric_limits<float>::epsilon())
        {
            fracErrs.push_back(std::numeric_limits<float>::max());
            continue;
        }

        fracErrs.push_back(xSecErr.at(i) / xSec.at(i));
    }

    return fracErrs;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

bool CrossSectionHelper::XSec::RemoveLowestStatsBin(const float targetFracError, std::vector<float> &binEdges)
{
    // Get the fractional errors on each bin of the reco cross-section (using MC)
    const auto fracErrs = this->GetFractionalRecoCrossSectionErrors(binEdges);

    // Now find the bin with the largest fractional error
    float largestFracErr = -std::numeric_limits<float>::max();
    unsigned int worstBinIndex = 0u;
    const auto nBins = fracErrs.size();
    for (unsigned int i = 0; i < nBins; ++i)
    {
        // Don't care about underflow or overflow bins
        if (this->IsUnderOverFlowBin(i))
            continue;

        const auto fracErr = fracErrs.at(i);

        if (fracErr > largestFracErr)
        {
            largestFracErr = fracErr;
            worstBinIndex = i;
        }
    }

    // Check to see if the fractional error is acceptable
    if (largestFracErr < targetFracError)
        return false;
    
    // We have a bin with a high fractional error
    // We need to work out which edge of the bin we should remove
    bool removeLeftEdge = true;

    if (worstBinIndex - (m_hasUnderflow ? 1u : 0u) == 0u) // First bin
    {
        removeLeftEdge = false;
    }
    else if (worstBinIndex == nBins - 1u - (m_hasOverflow ? 1u : 0u)) // Last bin
    {
        removeLeftEdge = true;
    }
    else // A middle bin
    {
        // Merge with the adjacent bin with the larger fractional uncertainty
        const auto fracErrorLeft = fracErrs.at(worstBinIndex - 1);
        const auto fracErrorRight = fracErrs.at(worstBinIndex + 1);

        removeLeftEdge = (fracErrorLeft > fracErrorRight);
    }

    // Remove the desired edge of the bad bin
    binEdges.erase(std::next(binEdges.begin(), worstBinIndex - (m_hasUnderflow ? 1u : 0u) + (removeLeftEdge ? 0u : 1u)));

    // We made a change
    return true;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void CrossSectionHelper::XSec::SetBinsAuto(const float lower, const float upper, const float targetFracError, const float targetSmearingDiagonal)
{
    // Get the initial bin edges
    // The aim here is to have as many bins as we can while still having reasonable statistics in each bin to calculate the fractional
    // uncertainties and smearing diagonal elements - the actual value doesn't matter that much
    const unsigned int minEntriesPerBin = 100u;
    auto binEdges = this->GetInitialBins(lower, upper, minEntriesPerBin);

    // Now iteratively merge / modify the bins
    while (true)
    {
        bool changeMade = false;

        if (this->RemoveMostSmearedBin(targetSmearingDiagonal, binEdges))
        {
            changeMade = true;
        } 
        else if (this->RemoveLowestStatsBin(targetFracError, binEdges))
        {
            changeMade = true;
        }

        const auto nBins = binEdges.size() - 1;

        // Break out if there is nothing left to do
        if (!changeMade || nBins == 1)
            break;
    }

    // Set the bin eges and we are done!
    this->SetBins(binEdges);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void CrossSectionHelper::XSec::FillEvent(const AnalysisHelper::SampleType sampleType, const float sampleNorm, const float weight, const float &recoValue, const float &trueValue, const bool isSignal, const bool isSelected, const std::string &classification)
{
    m_outputEvent.m_sampleType = sampleType;
    m_outputEvent.m_sampleNorm = sampleNorm;
    m_outputEvent.m_weight = weight;
    m_outputEvent.m_recoValue = recoValue;
    m_outputEvent.m_trueValue = trueValue;
    m_outputEvent.m_isSignal = isSignal;
    m_outputEvent.m_isSelected = isSelected;
    m_outputEvent.m_classification = classification;

    m_pTree->Fill();
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void CrossSectionHelper::XSec::AddSignalEvent(const std::shared_ptr<Event> &pEvent, const bool &isSelected, const float &trueValue, const float &recoValue, const float weight, const float normalisation)
{
    const auto classification = AnalysisHelper::GetClassificationString(pEvent, m_useAbsPdg, m_countProtonsInclusively);
    this->FillEvent(AnalysisHelper::Overlay, normalisation, weight, recoValue, trueValue, true, isSelected, classification);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void CrossSectionHelper::XSec::AddSelectedBackgroundEvent(const std::shared_ptr<Event> &pEvent, const AnalysisHelper::SampleType sampleType, const float &recoValue, const float weight, const float normalisation)
{
    if (sampleType == AnalysisHelper::DataBNB)
        throw std::invalid_argument("XSec::AddSelectedBackgroundEvent - You can't add BNB data as a known background");

    const auto classification = AnalysisHelper::GetClassificationString(pEvent, m_useAbsPdg, m_countProtonsInclusively);
    this->FillEvent(sampleType, normalisation, weight, recoValue, std::numeric_limits<float>::max(), false, true, classification);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void CrossSectionHelper::XSec::AddSelectedBNBDataEvent(const std::shared_ptr<Event> &pEvent, const float &recoValue)
{
    const auto classification = AnalysisHelper::GetClassificationString(pEvent, m_useAbsPdg, m_countProtonsInclusively);
    this->FillEvent(AnalysisHelper::DataBNB, 1.f, 1.f, recoValue, std::numeric_limits<float>::max(), false, true, classification);
}
                
// -----------------------------------------------------------------------------------------------------------------------------------------

unsigned int CrossSectionHelper::XSec::GetBinIndex(const float &value) const
{
    if (m_binEdges.size() < 2)
        throw std::logic_error("XSec::GetBinIndex - binning hasn't been set");

    const auto nBins = m_binEdges.size() - 1;
    for (unsigned int i = 0; i < nBins; ++i)
    {
        const auto binLower = m_binEdges.at(i);
        const auto binUpper = m_binEdges.at(i+1);

        if (value >= binLower && value <= binUpper)
            return i;
    }
   
    std::cout << "DEBUG - Value = " << value << std::endl;
    throw std::logic_error("XSec::GetBinIndex - input value is out of bounds");
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void CrossSectionHelper::XSec::GetConfusionMatrix(std::vector< std::vector<float> > &matrix, std::vector< std::vector<float> > &uncertainties) const
{
    if (m_binEdges.size() < 2)
        throw std::logic_error("XSec::GetConfusionMatrix - binning hasn't been set");
    
    if (!matrix.empty())
        throw std::logic_error("XSec::GetConfusionMatrix - input matrix isn't empty");
    
    if (!uncertainties.empty())
        throw std::logic_error("XSec::GetConfusionMatrix - input uncertainty matrix isn't empty");

    const auto nBins = m_binEdges.size() - 1;

    matrix = std::vector< std::vector<float> >(nBins, std::vector<float>(nBins, 0.f));
    uncertainties = std::vector< std::vector<float> >(nBins, std::vector<float>(nBins, 0.f));
    
    // Loop through the tree
    for (unsigned int i = 0; i < m_pTree->GetEntries(); ++i)
    {
        m_pTree->GetEntry(i);

        // We only care about overlays
        if (m_outputEvent.m_sampleType != AnalysisHelper::Overlay)
            continue;
        
        // Only use selected signal events
        if (!m_outputEvent.m_isSelected || !m_outputEvent.m_isSignal)
            continue;

        const auto trueBinIndex = this->GetBinIndex(m_outputEvent.m_trueValue);
        const auto recoBinIndex = this->GetBinIndex(m_outputEvent.m_recoValue);

        matrix.at(recoBinIndex).at(trueBinIndex) += m_outputEvent.m_weight;
    }

    for (unsigned int iReco = 0; iReco < nBins; ++iReco)
    {
        for (unsigned int iTrue = 0; iTrue < nBins; ++iTrue)
        {
            uncertainties.at(iReco).at(iTrue) = AnalysisHelper::GetCountUncertainty(matrix.at(iReco).at(iTrue));
        }
    }
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void CrossSectionHelper::XSec::GetSmearingMatrix(std::vector< std::vector<float> > &matrix, std::vector< std::vector<float> > &uncertainties) const
{
    if (!matrix.empty())
        throw std::logic_error("XSec::GetSmearingMatrix - input matrix isn't empty");
    
    if (!uncertainties.empty())
        throw std::logic_error("XSec::GetSmearingMatrix - input uncertainty matrix isn't empty");

    std::vector< std::vector<float> > confusionMatrix, confusionMatrixUncertainties;
    this->GetConfusionMatrix(confusionMatrix, confusionMatrixUncertainties);

    const auto nBins = confusionMatrix.size();
    matrix = std::vector< std::vector<float> >(nBins, std::vector<float>(nBins, 0.f));
    uncertainties = std::vector< std::vector<float> >(nBins, std::vector<float>(nBins, 0.f));

    // Get the total number of true entries in each bin
    auto truthWeights = std::vector<float>(nBins, 0.f);
    for (unsigned int iTrue = 0; iTrue < nBins; ++iTrue)
    {
        for (unsigned int iReco = 0; iReco < nBins; ++iReco)
        {
            truthWeights.at(iTrue) += confusionMatrix.at(iReco).at(iTrue);   
        }    
    }
    
    // Populate the smearing matrix
    for (unsigned int iTrue = 0; iTrue < nBins; ++iTrue)
    {
        const auto truthWeight = truthWeights.at(iTrue);
        if (truthWeight <= std::numeric_limits<float>::epsilon())
            throw std::invalid_argument("XSec::GetSmearingMatrix - Found truth bin with total weight: " + std::to_string(truthWeight));

        for (unsigned int iReco = 0; iReco < nBins; ++iReco)
        {
            matrix.at(iReco).at(iTrue) = confusionMatrix.at(iReco).at(iTrue) / truthWeight; 
            uncertainties.at(iReco).at(iTrue) = AnalysisHelper::GetEfficiencyUncertainty(confusionMatrix.at(iReco).at(iTrue), truthWeight);
        }    
    }
}

// -----------------------------------------------------------------------------------------------------------------------------------------

std::vector<float> CrossSectionHelper::XSec::GetTrueSignalPerBin(const bool mustBeSelected) const
{
    if (m_binEdges.size() < 2)
        throw std::logic_error("XSec::GetTrueSignalPerBin - binning hasn't been set");

    const auto nBins = m_binEdges.size() - 1;
    auto bins = std::vector<float>(nBins, 0.f);

    // Loop through the tree
    for (unsigned int i = 0; i < m_pTree->GetEntries(); ++i)
    {
        m_pTree->GetEntry(i);

        // We only care about overlays
        if (m_outputEvent.m_sampleType != AnalysisHelper::Overlay)
            continue;
        
        // Only use signal events
        if (!m_outputEvent.m_isSignal)
            continue;

        // Insist the event is selected if required
        if (mustBeSelected && !m_outputEvent.m_isSelected)
            continue;

        // Fill the bins
        const auto trueBinIndex = this->GetBinIndex(m_outputEvent.m_trueValue);
        bins.at(trueBinIndex) += m_outputEvent.m_weight;
    }

    return bins;
}

// -----------------------------------------------------------------------------------------------------------------------------------------
                
float CrossSectionHelper::XSec::GetOverlayNormalisation() const
{
    // Loop through the tree
    for (unsigned int i = 0; i < m_pTree->GetEntries(); ++i)
    {
        m_pTree->GetEntry(i);

        // ATTN here we are assume all events have same sample norm (which should be true, if it's not something is going wrong)
        if (m_outputEvent.m_sampleType == AnalysisHelper::Overlay)
            return m_outputEvent.m_sampleNorm;
    }

    throw std::logic_error("XSec::GetOverlayNormalisation - no overlays added");
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void CrossSectionHelper::XSec::GetBNBDataPerBin(std::vector<float> &bins, std::vector<float> &uncertainties, const bool useRealData) const
{
    if (m_binEdges.size() < 2)
        throw std::logic_error("XSec::GetBNBDataPerBin - binning hasn't been set");
    
    if (!bins.empty())
        throw std::invalid_argument("XSec::GetBNBDataPerBin - input bins vector isn't empty");
    
    if (!uncertainties.empty())
        throw std::invalid_argument("XSec::GetBNBDataPerBin - input uncertainties vector isn't empty");

    const auto nBins = m_binEdges.size() - 1;
    bins.resize(nBins, 0.f);
    uncertainties.resize(nBins, 0.f);
    
    // Loop through the tree
    for (unsigned int i = 0; i < m_pTree->GetEntries(); ++i)
    {
        m_pTree->GetEntry(i);

        // Check we have the type of data we want
        const auto isBNBData = m_outputEvent.m_sampleType == AnalysisHelper::DataBNB;
        if (useRealData && !isBNBData)
            continue;

        if (!useRealData && isBNBData)
            continue;
        
        // Insist the event is selected
        if (!m_outputEvent.m_isSelected)
            continue;

        // Fill the bins
        const auto recoBinIndex = this->GetBinIndex(m_outputEvent.m_recoValue);
        bins.at(recoBinIndex) += m_outputEvent.m_weight * m_outputEvent.m_sampleNorm;
    }
     
    // Fill the uncertainties
    for (unsigned int i = 0; i < nBins; ++i)
        uncertainties.at(i) = AnalysisHelper::GetCountUncertainty(bins.at(i));
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void CrossSectionHelper::XSec::GetBackgroundsPerBin(std::vector<float> &bins, std::vector<float> &uncertainties) const
{
    if (m_binEdges.size() < 2)
        throw std::logic_error("XSec::GetBackgroundsPerBin - binning hasn't been set");
    
    if (!bins.empty())
        throw std::invalid_argument("XSec::GetBackgroundsPerBin - input bins vector isn't empty");
    
    if (!uncertainties.empty())
        throw std::invalid_argument("XSec::GetBackgroundsPerBin - input uncertainties vector isn't empty");

    const auto nBins = m_binEdges.size() - 1;
    bins.resize(nBins, 0.f);
    uncertainties.resize(nBins, 0.f);

    std::vector< std::unordered_map<int, float> > weightsPerSampleType(nBins);
    std::unordered_map<int, float> sampleTypeNorms;
    std::vector<int> sampleTypes;

    // Loop through the tree
    for (unsigned int i = 0; i < m_pTree->GetEntries(); ++i)
    {
        m_pTree->GetEntry(i);

        // We only don't care about BNB data
        if (m_outputEvent.m_sampleType == AnalysisHelper::DataBNB)
            continue;
        
        // Insist the event is selected
        if (!m_outputEvent.m_isSelected)
            continue;

        // Insist the event isn't true signal
        if (m_outputEvent.m_sampleType == AnalysisHelper::Overlay && m_outputEvent.m_isSignal)
            continue;

        // Store the weights and normalisations
        const auto recoBinIndex = this->GetBinIndex(m_outputEvent.m_recoValue);

        auto &sampleWeightMap = weightsPerSampleType.at(recoBinIndex);
        if (sampleWeightMap.find(m_outputEvent.m_sampleType) == sampleWeightMap.end())
        {
            sampleWeightMap.emplace(m_outputEvent.m_sampleType, m_outputEvent.m_weight);
        }
        else
        {
            sampleWeightMap.at(m_outputEvent.m_sampleType) += m_outputEvent.m_weight;
        }

        sampleTypeNorms[m_outputEvent.m_sampleType] = m_outputEvent.m_sampleNorm;

        if (std::find(sampleTypes.begin(), sampleTypes.end(), m_outputEvent.m_sampleType) == sampleTypes.end())
            sampleTypes.push_back(m_outputEvent.m_sampleType);
    }
     
    // Fill the bins
    for (unsigned int i = 0; i < nBins; ++i)
    {
        float uncertaintySquared = 0.f;
        float binValue = 0.f;

        for (const auto &sampleType : sampleTypes)
        {
            const auto iter = weightsPerSampleType.at(i).find(sampleType);

            if (iter == weightsPerSampleType.at(i).end())
                continue;

            const auto weight = iter->second;
            const auto norm = sampleTypeNorms.at(sampleType);

            binValue += weight * norm;
            uncertaintySquared += std::pow(norm * AnalysisHelper::GetCountUncertainty(weight), 2);
        }

        uncertainties.at(i) = std::pow(uncertaintySquared, 0.5f);
        bins.at(i) = binValue;
    }
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void CrossSectionHelper::XSec::GetTrueEfficiencies(std::vector<float> &efficiencies, std::vector<float> &uncertainties) const
{
    if (m_binEdges.size() < 2)
        throw std::logic_error("XSec::GetTrueEfficiencies - binning hasn't been set");

    if (!efficiencies.empty())
        throw std::invalid_argument("XSec::GetTrueEfficiencies - input efficiencies vector isn't empty");
    
    if (!uncertainties.empty())
        throw std::invalid_argument("XSec::GetTrueEfficiencies - input uncertainties vector isn't empty");

    const auto nBins = m_binEdges.size() - 1;
    const auto binsAll = this->GetTrueSignalPerBin(false);
    const auto binsSelected = this->GetTrueSignalPerBin(true);

    efficiencies.resize(nBins, 0.f);
    uncertainties.resize(nBins, 0.f);
    
    // Calculate the efficiencies
    for (unsigned int i = 0; i < nBins; ++i)
    {
        const auto numerator = binsSelected.at(i);
        const auto denominator = binsAll.at(i);

        if (denominator <= std::numeric_limits<float>::epsilon())
            throw std::logic_error("XSec::GetTrueEfficiencies - found true bin with no entries");

        efficiencies.at(i) = numerator / denominator;
        uncertainties.at(i) = AnalysisHelper::GetEfficiencyUncertainty(numerator, denominator);
    }
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void CrossSectionHelper::XSec::GetSmearedEfficiencies(std::vector<float> &efficiencies, std::vector<float> &uncertainties) const
{
    if (m_binEdges.size() < 2)
        throw std::logic_error("XSec::GetSmearedEfficiencies - binning hasn't been set");

    if (!efficiencies.empty())
        throw std::invalid_argument("XSec::GetSmearedEfficiencies - input efficiencies vector isn't empty");
    
    if (!uncertainties.empty())
        throw std::invalid_argument("XSec::GetSmearedEfficiencies - input uncertainties vector isn't empty");
    
    const auto nBins = m_binEdges.size() - 1;
    efficiencies.resize(nBins, 0.f);
    uncertainties.resize(nBins, 0.f);

    // The the number of events before and after the selection in true bins
    const auto binsAll = this->GetTrueSignalPerBin(false);
    const auto binsSelected = this->GetTrueSignalPerBin(true);

    std::vector<float> binsAllErr, binsSelectedErr;
    for (unsigned int i = 0; i < nBins; ++i)
    {
        binsAllErr.push_back(AnalysisHelper::GetCountUncertainty(binsAll.at(i)));
        binsSelectedErr.push_back(AnalysisHelper::GetCountUncertainty(binsSelected.at(i)));
    }

    // Smear these counts
    std::vector<float> binsAllSmear, binsAllSmearErr, binsSelectedSmear, binsSelectedSmearErr;
    this->SmearBins(binsAll, binsAllErr, binsAllSmear, binsAllSmearErr);
    this->SmearBins(binsSelected, binsSelectedErr, binsSelectedSmear, binsSelectedSmearErr);
    
    // Calculate the smeared efficiencies
    for (unsigned int iReco = 0; iReco < nBins; ++iReco)
    {
        // Smear the numerator and denominator of the efficiency
        float allSmeared = binsAllSmear.at(iReco);
        float allSmearedErrSquared = std::pow(binsAllSmearErr.at(iReco), 2);

        float selectedSmeared = binsSelectedSmear.at(iReco);
        float selectedSmearedErrSquared = std::pow(binsSelectedSmearErr.at(iReco), 2);

        if (allSmeared <= std::numeric_limits<float>::epsilon())
            throw std::logic_error("XSec::GetSmearedEfficiencies - found smeared bin with no entries");
       
        if (selectedSmeared <= std::numeric_limits<float>::epsilon())
            throw std::logic_error("XSec::GetSmearedEfficiencies - found smeared bin with no selected entries");

        const auto efficiencySmeared = selectedSmeared / allSmeared;
        const auto uncertaintySmeared = efficiencySmeared * std::pow((selectedSmearedErrSquared / std::pow(selectedSmeared, 2)) + (allSmearedErrSquared / std::pow(allSmeared, 2)), 0.5f);

        efficiencies.at(iReco) = efficiencySmeared;
        uncertainties.at(iReco) = uncertaintySmeared;
    }
}

// -----------------------------------------------------------------------------------------------------------------------------------------

std::vector<float> CrossSectionHelper::XSec::GetBinWidths() const
{
    if (m_binEdges.size() < 2)
        throw std::logic_error("XSec::GetBinWidths - binning hasn't been set");
    
    const auto nBins = m_binEdges.size() - 1;
    
    std::vector<float> binWidths;
    for (unsigned int i = 0; i < nBins; ++i)
    {
        binWidths.push_back(m_binEdges.at(i+1) - m_binEdges.at(i));
    }

    return binWidths;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

float CrossSectionHelper::XSec::GetCrossSection(const float nTotal, const float nBackgrounds, const float efficiency, const float binWidth) const
{
    if (binWidth <= std::numeric_limits<float>::epsilon())
        throw std::invalid_argument("XSec::GetCrossSection - bin width is zero");
    
    if (efficiency <= std::numeric_limits<float>::epsilon())
        throw std::invalid_argument("XSec::GetCrossSection - efficiency is zero");

    return (nTotal - nBackgrounds) / (efficiency * binWidth * m_nTargets * m_integratedFlux);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

float CrossSectionHelper::XSec::GetCrossSectionUncertainty(const float nTotal, const float totalErr, const float nBackgrounds, const float backgroundsErr, const float efficiency, const float efficiencyErr, const float binWidth) const
{
    if (binWidth <= std::numeric_limits<float>::epsilon())
        throw std::invalid_argument("XSec::GetCrossSectionUncertainty - bin width is zero");
    
    if (efficiency <= std::numeric_limits<float>::epsilon())
        throw std::invalid_argument("XSec::GetCrossSectionUncertainty - efficiency is zero");

    const auto nSignalSelected = nTotal - nBackgrounds;
    const auto nSignalSelectedErr = std::pow(totalErr*totalErr + backgroundsErr*backgroundsErr, 0.5f);

    if (nSignalSelected <= std::numeric_limits<float>::epsilon())
        throw std::invalid_argument("XSec::GetCrossSectionUncertainty - found bin with no expected signal events");

    const auto nSignal = nSignalSelected / efficiency;
    const auto nSignalErr = nSignal * std::pow(std::pow(nSignalSelectedErr / nSignalSelected, 2) + std::pow(efficiencyErr / efficiency, 2), 0.5f);

    return nSignalErr / (binWidth * m_nTargets * m_integratedFlux);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void CrossSectionHelper::XSec::GetMCTrueCrossSection(std::vector<float> &bins, std::vector<float> &uncertainties) const
{
    if (m_binEdges.size() < 2)
        throw std::logic_error("XSec::GetMCTrueCrossSection - binning hasn't been set");
    
    if (!bins.empty())
        throw std::invalid_argument("XSec::GetMCTrueCrossSection - input bins vector isn't empty");
    
    if (!uncertainties.empty())
        throw std::invalid_argument("XSec::GetMCTrueCrossSection - input uncertainties vector isn't empty");

    const auto nBins = m_binEdges.size() - 1;
    bins.resize(nBins, 0.f);
    uncertainties.resize(nBins, 0.f);

    const auto binWidths = this->GetBinWidths();
    const auto nSignalPerBin = this->GetTrueSignalPerBin(false);
    const auto norm = this->GetOverlayNormalisation();

    for (unsigned int i = 0; i < nBins; ++i)
    {
        const auto nTotal = norm * nSignalPerBin.at(i);
        const auto totalErr = norm * AnalysisHelper::GetCountUncertainty(nTotal);
        const auto binWidth = binWidths.at(i);

        bins.at(i) = this->GetCrossSection(nTotal, 0.f, 1.f, binWidth);
        uncertainties.at(i) = this->GetCrossSectionUncertainty(nTotal, totalErr, 0.f, 0.f, 1.f, 0.f, binWidth);
    }
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void CrossSectionHelper::XSec::SmearBins(const std::vector<float> &inputBins, const std::vector<float> &inputUncertainties, std::vector<float> &bins, std::vector<float> &uncertainties) const
{
    if (m_binEdges.size() < 2)
        throw std::logic_error("XSec::SmearBins - binning hasn't been set");
    
    if (!bins.empty())
        throw std::invalid_argument("XSec::SmearBins - bins vector isn't empty");
    
    if (!uncertainties.empty())
        throw std::invalid_argument("XSec::SmearBins - uncertainties vector isn't empty");

    const auto nBins = m_binEdges.size() - 1;
    bins.resize(nBins, 0.f);
    uncertainties.resize(nBins, 0.f);

    if (inputBins.size() != nBins)
        throw std::invalid_argument("XSec::SmearBins - input bins has wrong size");
    
    if (inputUncertainties.size() != nBins)
        throw std::invalid_argument("XSec::SmearBins - input uncertainties has wrong size");

    std::vector< std::vector<float> > smearingMatrix, smearingMatrixUncertainties;
    this->GetSmearingMatrix(smearingMatrix, smearingMatrixUncertainties);

    // Apply the smearing matrix
    for (unsigned int iReco = 0; iReco < nBins; ++iReco)
    {
        float smearedValue = 0.f;
        float smearedErrSquared = 0.f;

        for (unsigned int iTrue = 0; iTrue < nBins; ++iTrue)
        {
            const auto smearingElement = smearingMatrix.at(iReco).at(iTrue);
            const auto smearingElementErr = smearingMatrixUncertainties.at(iReco).at(iTrue);
            const auto binValue = inputBins.at(iTrue);
            const auto binValueErr = inputUncertainties.at(iTrue);

            const auto contribution = smearingElement * binValue;
            smearedValue += contribution;
            
            if (smearingElement > std::numeric_limits<float>::epsilon() && binValue > std::numeric_limits<float>::epsilon())
            {
                smearedErrSquared += std::pow(contribution, 2) * (std::pow(smearingElementErr / smearingElement, 2) + std::pow(binValueErr / binValue, 2));
            }
        }

        bins.at(iReco) = smearedValue;
        uncertainties.at(iReco) = std::pow(smearedErrSquared, 0.5f);
    }
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void CrossSectionHelper::XSec::GetSmearedMCTrueCrossSection(std::vector<float> &bins, std::vector<float> &uncertainties) const
{
    if (m_binEdges.size() < 2)
        throw std::logic_error("XSec::GetSmearedMCTrueCrossSection - binning hasn't been set");
    
    if (!bins.empty())
        throw std::invalid_argument("XSec::GetSmearedMCTrueCrossSection - input bins vector isn't empty");
    
    if (!uncertainties.empty())
        throw std::invalid_argument("XSec::GetSmearedMCTrueCrossSection - input uncertainties vector isn't empty");

    const auto nBins = m_binEdges.size() - 1;
    bins.resize(nBins, 0.f);
    uncertainties.resize(nBins, 0.f);

    std::vector<float> mcTrueXSec, mcTrueXSecErr;
    this->GetMCTrueCrossSection(mcTrueXSec, mcTrueXSecErr);

    // Scale up by bin widths before smearing
    const auto binWidths = this->GetBinWidths();
    std::vector<float> mcTrueXSecScaled, mcTrueXSecErrScaled;
    for (unsigned int i = 0; i < nBins; ++i)
    {
        mcTrueXSecScaled.push_back(mcTrueXSec.at(i) * binWidths.at(i));
        mcTrueXSecErrScaled.push_back(mcTrueXSecErr.at(i) * binWidths.at(i));
    }

    // Now apply the smearing matrix
    std::vector<float> mcTrueXSecScaledSmeared, mcTrueXSecErrScaledSmeared;
    this->SmearBins(mcTrueXSecScaled, mcTrueXSecErrScaled, mcTrueXSecScaledSmeared, mcTrueXSecErrScaledSmeared);

    // Scale back down by the bin widths after smearing
    for (unsigned int i = 0; i < nBins; ++i)
    {
        const auto binWidth = binWidths.at(i);

        if (binWidth <= std::numeric_limits<float>::epsilon())
            throw std::invalid_argument("XSec::GetSmearedMCTrueCrossSection - bin width is zero");

        bins.at(i) = mcTrueXSecScaledSmeared.at(i) / binWidth;
        uncertainties.at(i) = mcTrueXSecErrScaledSmeared.at(i) / binWidth;
    }
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void CrossSectionHelper::XSec::GetRecoCrossSection(const bool useRealData, std::vector<float> &bins, std::vector<float> &uncertainties) const
{
    if (m_binEdges.size() < 2)
        throw std::logic_error("XSec::GetRecoCrossSection - binning hasn't been set");
    
    if (!bins.empty())
        throw std::invalid_argument("XSec::GetRecoCrossSection - input bins vector isn't empty");
    
    if (!uncertainties.empty())
        throw std::invalid_argument("XSec::GetRecoCrossSection - input uncertainties vector isn't empty");

    const auto nBins = m_binEdges.size() - 1;
    bins.resize(nBins, 0.f);
    uncertainties.resize(nBins, 0.f);

    std::vector<float> selected, selectedErr;
    this->GetBNBDataPerBin(selected, selectedErr, useRealData);
    
    std::vector<float> backgrounds, backgroundsErr;
    this->GetBackgroundsPerBin(backgrounds, backgroundsErr);

    std::vector<float> efficiencies, efficienciesErr;
    this->GetSmearedEfficiencies(efficiencies, efficienciesErr);

    const auto binWidths = this->GetBinWidths();

    for (unsigned int i = 0; i < nBins; ++i)
    {
        if (this->IsUnderOverFlowBin(i))
            continue;

        const auto xSec = this->GetCrossSection(selected.at(i), backgrounds.at(i), efficiencies.at(i), binWidths.at(i));
        const auto xSecErr = this->GetCrossSectionUncertainty(selected.at(i), selectedErr.at(i), backgrounds.at(i), backgroundsErr.at(i), efficiencies.at(i), efficienciesErr.at(i), binWidths.at(i));
        
        bins.at(i) = xSec;
        uncertainties.at(i) = xSecErr;
    }
}

// -----------------------------------------------------------------------------------------------------------------------------------------

std::string CrossSectionHelper::XSec::GetBinName(const unsigned int index) const
{
    if (m_binEdges.size() < 2)
        throw std::logic_error("XSec::GetBinName - binning hasn't been set");

    const auto nBins = m_binEdges.size() - 1;

    if (index >= nBins)
        throw std::invalid_argument("XSec::GetBinName - input index is out of bounds");

    if (index == 0 && m_hasUnderflow)
        return "Underflow";

    if (index == nBins - 1 && m_hasOverflow)
        return "Overflow";

    return std::to_string(index - (m_hasUnderflow ? 1 : 0));
}

// -----------------------------------------------------------------------------------------------------------------------------------------

bool CrossSectionHelper::XSec::IsUnderOverFlowBin(const unsigned int index) const
{
    if (m_binEdges.size() < 2)
        throw std::logic_error("XSec::IsUnderOverFlowBin - binning hasn't been set");

    const auto nBins = m_binEdges.size() - 1;

    if (index == 0 && m_hasUnderflow)
        return true;

    if (index == nBins - 1 && m_hasOverflow)
        return true;

    return false;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void CrossSectionHelper::XSec::PrintMatrix(const std::vector< std::vector<float> > &matrix, const std::vector< std::vector<float> > &uncertainties, const std::string &outputFileName) const
{
    // Fill the headers with the truth bins
    std::vector<std::string> headers = {"Bin"};
    for (unsigned int iTrue = 0; iTrue < matrix.size(); ++iTrue)
        headers.push_back("T " + this->GetBinName(iTrue));
    
    FormattingHelper::Table table(headers);

    for (unsigned int iReco = 0; iReco < matrix.size(); ++iReco)
    {
        table.AddEmptyRow();
        table.SetEntry("Bin", "R " + this->GetBinName(iReco));

        for (unsigned int iTrue = 0; iTrue < matrix.size(); ++iTrue)
        {
            const auto value = FormattingHelper::GetValueWithError(matrix.at(iReco).at(iTrue), uncertainties.at(iReco).at(iTrue));
            table.SetEntry("T " + this->GetBinName(iTrue), value);
        }
    }

    table.WriteToFile(outputFileName);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void CrossSectionHelper::XSec::PrintBinContents(const std::string &outputFileNamePrefix) const
{
    if (m_binEdges.size() < 2)
        throw std::logic_error("XSec::PrintBinContents - binning hasn't been set");

    const auto nBins = m_binEdges.size() - 1;
                
    // Get the bin values
    std::vector<float> bnbData, bnbDataErr, fakeBNBData, fakeBNBDataErr, backgrounds, backgroundsErr;
    this->GetBNBDataPerBin(bnbData, bnbDataErr, true);
    this->GetBNBDataPerBin(fakeBNBData, fakeBNBDataErr, false);
    this->GetBackgroundsPerBin(backgrounds, backgroundsErr);
    const auto trueSignal = this->GetTrueSignalPerBin(true);

    // Print the bins
    FormattingHelper::PrintLine();
    std::cout << "Bins" << std::endl;
    FormattingHelper::PrintLine();
    FormattingHelper::Table binTable({"Bin", "", "Min", "Max", "", "True MC signal", "", "BNB data", "Fake data (MC)", "data - MC", "", "Backgrounds MC"});
    for (unsigned int i = 0; i < nBins; ++i)
    {
        binTable.AddEmptyRow();
        binTable.SetEntry("Bin", this->GetBinName(i));
        binTable.SetEntry("Min", m_binEdges.at(i));
        binTable.SetEntry("Max", m_binEdges.at(i+1));
        binTable.SetEntry("True MC signal", trueSignal.at(i));
        binTable.SetEntry("BNB data", FormattingHelper::GetValueWithError(bnbData.at(i), bnbDataErr.at(i)));
        binTable.SetEntry("Fake data (MC)", FormattingHelper::GetValueWithError(fakeBNBData.at(i), fakeBNBDataErr.at(i)));

        const auto dataMinusMC = bnbData.at(i) - fakeBNBData.at(i);
        const auto dataMinusMCErr = std::pow(std::pow(bnbDataErr.at(i), 2) + std::pow(fakeBNBDataErr.at(i), 2), 0.5f);
        binTable.SetEntry("data - MC", FormattingHelper::GetValueWithError(dataMinusMC, dataMinusMCErr));

        binTable.SetEntry("Backgrounds MC", FormattingHelper::GetValueWithError(backgrounds.at(i), backgroundsErr.at(i)));
    }
    binTable.WriteToFile(outputFileNamePrefix + "_bins.md");
    
    // Print the confusion matrix
    FormattingHelper::PrintLine();
    std::cout << "Confusion matrix" << std::endl;
    FormattingHelper::PrintLine();

    std::vector< std::vector<float> > confusionMatrix, confusionMatrixUncertainties;
    this->GetConfusionMatrix(confusionMatrix, confusionMatrixUncertainties);
    this->PrintMatrix(confusionMatrix, confusionMatrixUncertainties, outputFileNamePrefix + "_confusionMatrix.md");
    
    // Print the smearing matrix
    FormattingHelper::PrintLine();
    std::cout << "Smearing matrix" << std::endl;
    FormattingHelper::PrintLine();

    std::vector< std::vector<float> > smearingMatrix, smearingMatrixUncertainties;
    this->GetSmearingMatrix(smearingMatrix, smearingMatrixUncertainties);
    this->PrintMatrix(smearingMatrix, smearingMatrixUncertainties, outputFileNamePrefix + "_smearingMatrix.md");
   
    // Print the efficiencies
    FormattingHelper::PrintLine();
    std::cout << "Selection efficiencies" << std::endl;
    FormattingHelper::PrintLine();
    FormattingHelper::Table efficiencyTable({"Bin", "True efficiency", "Smeared efficiency"});
    
    std::vector<float> trueEfficiencies, trueEfficiencyUncertainties;
    this->GetTrueEfficiencies(trueEfficiencies, trueEfficiencyUncertainties);
    
    std::vector<float> recoEfficiencies, recoEfficiencyUncertainties;
    this->GetSmearedEfficiencies(recoEfficiencies, recoEfficiencyUncertainties);

    for (unsigned int i = 0; i < nBins; ++i)
    {
        const auto trueEfficiency = FormattingHelper::GetValueWithError(trueEfficiencies.at(i), trueEfficiencyUncertainties.at(i));
        const auto recoEfficiency = FormattingHelper::GetValueWithError(recoEfficiencies.at(i), recoEfficiencyUncertainties.at(i));

        efficiencyTable.AddEmptyRow();
        efficiencyTable.SetEntry("Bin", this->GetBinName(i));
        efficiencyTable.SetEntry("True efficiency", trueEfficiency);
        efficiencyTable.SetEntry("Smeared efficiency", recoEfficiency);
    }
    efficiencyTable.WriteToFile(outputFileNamePrefix + "_efficiencies.md");

    // Print the extracted cross-sections
    FormattingHelper::PrintLine();
    std::cout << "Cross-sections" << std::endl;
    FormattingHelper::PrintLine();
    FormattingHelper::Table xSecTable({"Bin", "", "xSec MC true", "", "xSec MC smeared", "xSec MC reco", "xSec Data"});

    std::vector<float> xSecMCTrue, xSecMCTrueErr, xSecMCSmeared, xSecMCSmearedErr, xSecMCReco, xSecMCRecoErr, xSecData, xSecDataErr;
    this->GetMCTrueCrossSection(xSecMCTrue, xSecMCTrueErr);
    this->GetSmearedMCTrueCrossSection(xSecMCSmeared, xSecMCSmearedErr);
    this->GetRecoCrossSection(false, xSecMCReco, xSecMCRecoErr);
    this->GetRecoCrossSection(true, xSecData, xSecDataErr);

    for (unsigned int i = 0; i < nBins; ++i)
    {
        xSecTable.AddEmptyRow();
        xSecTable.SetEntry("Bin", this->GetBinName(i));
        xSecTable.SetEntry("xSec MC true", FormattingHelper::GetValueWithError(xSecMCTrue.at(i), xSecMCTrueErr.at(i)));
        xSecTable.SetEntry("xSec MC smeared", FormattingHelper::GetValueWithError(xSecMCSmeared.at(i), xSecMCSmearedErr.at(i)));
        xSecTable.SetEntry("xSec MC reco", FormattingHelper::GetValueWithError(xSecMCReco.at(i), xSecMCRecoErr.at(i)));
        xSecTable.SetEntry("xSec Data", FormattingHelper::GetValueWithError(xSecData.at(i), xSecDataErr.at(i)));
    }

    xSecTable.WriteToFile(outputFileNamePrefix + "_xSecs.md");

}

// -----------------------------------------------------------------------------------------------------------------------------------------

// TODO this function is vomit on toast.
void CrossSectionHelper::XSec::MakePlots(const std::string &fileNamePrefix) const
{
    if (m_binEdges.size() < 2)
        throw std::logic_error("XSec::MakePlotsPrint - binning hasn't been set");

    const auto nBins = m_binEdges.size() - 1;

    auto pCanvas = PlottingHelper::GetCanvas();

    // Make the smearing matrix plot
    std::vector< std::vector<float> > smearingMatrix, smearingMatrixUncertainties;
    this->GetSmearingMatrix(smearingMatrix, smearingMatrixUncertainties);

    TH2F *hSmearing = new TH2F(("hSmearing_" + fileNamePrefix).c_str(), "", nBins, 0, nBins, nBins, 0, nBins);
    TH2F *hSmearingErr = new TH2F(("hSmearingErr_" + fileNamePrefix).c_str(), "", nBins, 0, nBins, nBins, 0, nBins);
    
    for (unsigned int i = 1; i <= nBins; ++i)
    {
        const auto recoName = "R " + this->GetBinName(i-1);
        const auto trueName = "T " + this->GetBinName(i-1);

        hSmearing->GetXaxis()->SetBinLabel(i, recoName.c_str());
        hSmearing->GetYaxis()->SetBinLabel(i, trueName.c_str());
        
        hSmearingErr->GetXaxis()->SetBinLabel(i, recoName.c_str());
        hSmearingErr->GetYaxis()->SetBinLabel(i, trueName.c_str());
    }

    for (unsigned int iReco = 1; iReco <= nBins; ++iReco)
    {
        for (unsigned int iTrue = 1; iTrue <= nBins; ++iTrue)
        {
            hSmearing->SetBinContent(iReco, iTrue, smearingMatrix.at(iReco-1).at(iTrue-1));
            hSmearingErr->SetBinContent(iReco, iTrue, smearingMatrixUncertainties.at(iReco-1).at(iTrue-1));
        }
    }

    hSmearing->Draw("colz");
    PlottingHelper::SaveCanvas(pCanvas, fileNamePrefix + "_smearingMatrix");
    
    hSmearingErr->Draw("colz");
    PlottingHelper::SaveCanvas(pCanvas, fileNamePrefix + "_smearingMatrixUncertainties");

    // Now plot the cross-sections
    std::vector<float> xSecMCTrue, xSecMCTrueErr, xSecMCSmeared, xSecMCSmearedErr, xSecMCReco, xSecMCRecoErr, xSecData, xSecDataErr;
    this->GetMCTrueCrossSection(xSecMCTrue, xSecMCTrueErr);
    this->GetSmearedMCTrueCrossSection(xSecMCSmeared, xSecMCSmearedErr);
    this->GetRecoCrossSection(false, xSecMCReco, xSecMCRecoErr);
    this->GetRecoCrossSection(true, xSecData, xSecDataErr);
    
    std::vector<float> trueEfficiencies, trueEfficienciesErr, smearedEfficiencies, smearedEfficienciesErr;
    this->GetTrueEfficiencies(trueEfficiencies, trueEfficienciesErr);
    this->GetSmearedEfficiencies(smearedEfficiencies, smearedEfficienciesErr);


    const auto nNormalBins = nBins - (m_hasUnderflow ? 1u : 0u) - (m_hasOverflow ? 1u : 0u);
    std::vector<float> normalBinEdges;
    for (unsigned int i = 0; i < nBins; ++i)
    {
        if (this->IsUnderOverFlowBin(i))
            continue;

        // Add the lower edge of the bin
        normalBinEdges.push_back(m_binEdges.at(i));

        // If this is the last bin, then include the upper edge too
        if (i - (m_hasUnderflow ? 1u : 0u) == nNormalBins - 1)
            normalBinEdges.push_back(m_binEdges.at(i+1));
    }

    // Fill the histograms
    TH1F *hXSecMCTrue = new TH1F(("hXSecMCTrue_" + fileNamePrefix).c_str(), "", nNormalBins, normalBinEdges.data());
    TH1F *hXSecMCSmeared = new TH1F(("hXSecMCSmeared_" + fileNamePrefix).c_str(), "", nNormalBins, normalBinEdges.data());
    TH1F *hXSecMCReco = new TH1F(("hXSecMCReco_" + fileNamePrefix).c_str(), "", nNormalBins, normalBinEdges.data());
    TH1F *hXSecData = new TH1F(("hXSecData_" + fileNamePrefix).c_str(), "", nNormalBins, normalBinEdges.data());
    TH1F *hEff = new TH1F(("hEff_" + fileNamePrefix).c_str(), "", nNormalBins, normalBinEdges.data());
    TH1F *hEffSmeared = new TH1F(("hEffSmeared_" + fileNamePrefix).c_str(), "", nNormalBins, normalBinEdges.data());
    hXSecMCTrue->SetLineWidth(2);
    hXSecMCSmeared->SetLineWidth(2);
    hXSecMCReco->SetLineWidth(2);
    hXSecData->SetLineWidth(2);
    hEff->SetLineWidth(2);
    hEffSmeared->SetLineWidth(2);
    hEff->SetLineColor(kAzure - 2);
    hEffSmeared->SetLineColor(kOrange - 3);

    float maxValue = -std::numeric_limits<float>::max();
    float maxEfficiency = 0.1f;

    const int binOffset = (m_hasUnderflow ? 1 : 0) - 1; // Here we -1 as root bins are enumerated from 1 and ours are from 0
    for (unsigned int i = 1; i <= nNormalBins; ++i)
    {
        if (this->IsUnderOverFlowBin(i + binOffset))
            continue;

        hXSecMCTrue->SetBinContent(i, xSecMCTrue.at(i + binOffset));
        hXSecMCTrue->SetBinError(i, xSecMCTrueErr.at(i + binOffset));
        maxValue = std::max(maxValue, xSecMCTrue.at(i + binOffset) + xSecMCTrueErr.at(i + binOffset));
        
        hXSecMCSmeared->SetBinContent(i, xSecMCSmeared.at(i + binOffset));
        hXSecMCSmeared->SetBinError(i, xSecMCSmearedErr.at(i + binOffset));
        maxValue = std::max(maxValue, xSecMCSmeared.at(i + binOffset) + xSecMCSmearedErr.at(i + binOffset));
        
        hXSecMCReco->SetBinContent(i, xSecMCReco.at(i + binOffset));
        hXSecMCReco->SetBinError(i, xSecMCRecoErr.at(i + binOffset));
        maxValue = std::max(maxValue, xSecMCReco.at(i + binOffset) + xSecMCRecoErr.at(i + binOffset));
        
        hXSecData->SetBinContent(i, xSecData.at(i + binOffset));
        hXSecData->SetBinError(i, xSecDataErr.at(i + binOffset));
        maxValue = std::max(maxValue, xSecData.at(i + binOffset) + xSecDataErr.at(i + binOffset));

        hEff->SetBinContent(i, trueEfficiencies.at(i + binOffset));
        hEff->SetBinError(i, trueEfficienciesErr.at(i + binOffset));
        maxEfficiency = std::max(maxEfficiency, trueEfficiencies.at(i + binOffset) + trueEfficienciesErr.at(i + binOffset));

        hEffSmeared->SetBinContent(i, smearedEfficiencies.at(i + binOffset));
        hEffSmeared->SetBinError(i, smearedEfficienciesErr.at(i + binOffset));
        maxEfficiency = std::max(maxEfficiency, smearedEfficiencies.at(i + binOffset) + smearedEfficienciesErr.at(i + binOffset));
    }

    const std::vector< std::pair<TH1F *, std::string> > histos = {
        {hXSecMCTrue, "xSecMCTrue"},
        {hXSecMCSmeared, "xSecMCSmeared"},
        {hXSecMCReco, "xSecMCReco"},
        {hXSecData, "xSecData"}
    };

    for (const auto &entry : histos)
    {
        const auto &pHist = entry.first;
        const auto &name = entry.second;

        pHist->GetYaxis()->SetRangeUser(0.f, maxValue * 1.05f);

        auto pHistClone = static_cast<TH1F *>(pHist->Clone());
        const auto col = pHistClone->GetLineColor();
        pHistClone->SetFillStyle(1001);
        pHistClone->SetLineColorAlpha(col, 0.f);
        pHistClone->SetFillColorAlpha(col, 0.3f);

        pHistClone->Draw("e2");
        pHist->Draw("hist same");
    
        PlottingHelper::SaveCanvas(pCanvas, fileNamePrefix + "_" + name);
    }

    // Sanity check plot
    auto hXSecMCRecoClone = static_cast<TH1F *>(hXSecMCReco->Clone());
    auto hXSecMCSmearedClone = static_cast<TH1F *>(hXSecMCSmeared->Clone());
    hXSecMCRecoClone->SetLineColor(kBlack);
    hXSecMCRecoClone->SetLineStyle(2);
    hXSecMCSmearedClone->SetLineColor(kOrange - 3);

    hXSecMCSmearedClone->Draw("hist");
    hXSecMCRecoClone->Draw("hist same");
    PlottingHelper::SaveCanvas(pCanvas, fileNamePrefix + "_sanityCheck");



    // Now make the plot with all x-sec histograms
    // True
    /*
    auto hXSecMCTrueClone = static_cast<TH1F *>(hXSecMCTrue->Clone());
    hXSecMCTrue->SetLineColor(kAzure - 2);
    hXSecMCTrueClone->SetFillStyle(1001);
    hXSecMCTrueClone->SetLineColorAlpha(kAzure - 2, 0.f);
    hXSecMCTrueClone->SetFillColorAlpha(kAzure - 2, 0.3f);
    hXSecMCTrueClone->GetYaxis()->SetRangeUser(0.f, maxValue * 1.05f);
    
    hXSecMCTrueClone->Draw("e2");
    hXSecMCTrue->Draw("hist same");
    */

    // Smeared
//    auto hXSecMCSmearedClone = static_cast<TH1F *>(hXSecMCSmeared->Clone());
    hXSecMCSmeared->SetLineColor(kOrange - 3);
    hXSecMCSmearedClone->SetFillStyle(1001);
    hXSecMCSmearedClone->SetLineColorAlpha(kOrange - 3, 0.f);
    hXSecMCSmearedClone->SetFillColorAlpha(kOrange - 3, 0.3f);
    
    hXSecMCSmearedClone->Draw("e2");
    hXSecMCSmeared->Draw("hist same");

    // Data
    hXSecData->SetLineColor(kBlack);
    hXSecData->Draw("e1 same");

    PlottingHelper::SaveCanvas(pCanvas, fileNamePrefix + "_combined");

    // Make the efficiency plots
    hEff->GetYaxis()->SetRangeUser(0.f, maxEfficiency * 1.05);
    hEffSmeared->GetYaxis()->SetRangeUser(0.f, maxEfficiency * 1.05);

    auto hEffClone = static_cast<TH1F *>(hEff->Clone());
    hEffClone->SetFillStyle(1001);
    hEffClone->SetLineColorAlpha(hEffClone->GetLineColor(), 0.f);
    hEffClone->SetFillColorAlpha(hEffClone->GetLineColor(), 0.3f);
    
    auto hEffSmearedClone = static_cast<TH1F *>(hEffSmeared->Clone());
    hEffSmearedClone->SetFillStyle(1001);
    hEffSmearedClone->SetLineColorAlpha(hEffSmearedClone->GetLineColor(), 0.f);
    hEffSmearedClone->SetFillColorAlpha(hEffSmearedClone->GetLineColor(), 0.3f);

    hEffClone->Draw("e2");
    hEff->Draw("hist same");
    hEffSmearedClone->Draw("e2 same");
    hEffSmeared->Draw("hist same");

    PlottingHelper::SaveCanvas(pCanvas, fileNamePrefix + "_efficiencies");


    // Ratio plot
    auto pCanvasRatio = PlottingHelper::GetCanvas(960, 270);

    auto hXSecDataRatio = static_cast<TH1F *>(hXSecData->Clone());
    hXSecDataRatio->Sumw2();
    hXSecMCSmeared->Sumw2();
    hXSecDataRatio->Divide(hXSecMCSmeared);
    
    // Set the y-axis range
    float minRatio = 0.75f;
    float maxRatio = 1.25f;
    
    std::cout << "TEST : " << fileNamePrefix << std::endl;
    float chi2Sum = 0.f;
    float nDof = 0.f;
    for (unsigned int i = 1; i <= static_cast<unsigned int>(hXSecDataRatio->GetNbinsX()); ++i)
    {
        minRatio = std::min(minRatio, static_cast<float>(hXSecDataRatio->GetBinContent(i) - hXSecDataRatio->GetBinError(i)));
        maxRatio = std::max(maxRatio, static_cast<float>(hXSecDataRatio->GetBinContent(i) + hXSecDataRatio->GetBinError(i)));

        const auto diff = (hXSecDataRatio->GetBinContent(i) - 1.f) / hXSecDataRatio->GetBinError(i);
        std::cout << "  - " << i << " : " << hXSecDataRatio->GetBinContent(i) << ", " << hXSecDataRatio->GetBinError(i) << ", " << diff << std::endl;
        chi2Sum += std::pow(diff, 2);
        nDof += 1.f;
    }
    std::cout << "  - Chi2 = " << chi2Sum << " (" << nDof << ")" << std::endl;
    std::cout << "  - Chi2/dof = " << (chi2Sum / nDof) << std::endl;

    const auto padding = (maxRatio - minRatio) * 0.05f;
    hXSecDataRatio->GetYaxis()->SetRangeUser(std::max(0.f, minRatio - padding), maxRatio + padding);

    hXSecDataRatio->Draw("e1");
    pCanvasRatio->Update();

    // Add the lines at 0.8, 1.0 and 1.2
    TLine *l=new TLine(pCanvasRatio->GetUxmin(),1.0,pCanvasRatio->GetUxmax(),1.0);
    TLine *lPlus=new TLine(pCanvasRatio->GetUxmin(), 1.2f, pCanvasRatio->GetUxmax(), 1.2f);
    TLine *lMinus=new TLine(pCanvasRatio->GetUxmin(), 0.8f, pCanvasRatio->GetUxmax(), 0.8f);

    l->SetLineColor(kOrange - 3);
    lPlus->SetLineColor(kOrange - 3);
    lMinus->SetLineColor(kOrange - 3);

    lPlus->SetLineStyle(2);
    lMinus->SetLineStyle(2);

    l->Draw();
    lPlus->Draw();
    lMinus->Draw();

    PlottingHelper::SaveCanvas(pCanvasRatio, fileNamePrefix + "_ratio");
}

} // namespace ubcc1pi

#endif
