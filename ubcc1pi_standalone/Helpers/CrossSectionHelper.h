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
#include <random>
#include <functional>

#include <TFile.h>
#include <TTree.h>
#include <TH2F.h>
#include <TLine.h>
#include <TStyle.h>

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
            private:
                /**
                 *  @brief  The event structure
                 *          Holds the basic information needed about an event to calculate the cross-section
                 */
                struct OutputEvent
                {
                    int          m_sampleType;      ///< The sample type enumeration
                    float        m_sampleNorm;      ///< The sample normalisation factor   
                    float        m_weight;          ///< The nominal event weight
                    float        m_recoValue;       ///< The reconstructed value
                    float        m_trueValue;       ///< The true value
                    bool         m_isSignal;        ///< If the event is signal
                    bool         m_isSelected;      ///< If the event is selected
                    std::string  m_classification;  ///< The event classification string

                    std::map<std::string, std::vector<float> > m_systWeights; ///< The systematic event weights indexed by [systematicName][universeIndex]
                };

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
                 *  @param  targetEntriesPerBin the target number of entries in each bin
                 *  @param  targetSmearingDiagonal the target entries on the diagonal of the smearing matrix
                 */
                void SetBinsAuto(const float lower, const float upper, const unsigned int targetEntriesPerBin, const float targetSmearingDiagonal);

                /**
                 *  @brief  Add a signal event 
                 *
                 *  @param  pEvent the event
                 *  @param  isSelected if the event has been selected
                 *  @param  trueValue the true value of the quantity to measure
                 *  @param  recoValue the reco value of the quantity to meaasure
                 *  @param  weight the nominal weight
                 *  @param  normalisation the sample normalisation
                 */
                void AddSignalEvent(const std::shared_ptr<Event> &pEvent, const bool &isSelected, const float &trueValue, const float &recoValue, const float weight, const float normalisation);

                /**
                 *  @brief  Add a background event that has been selected
                 *
                 *  @param  pEvent the event
                 *  @param  sampleType the sample type (Overlays, EXT, Dirt)
                 *  @param  recoValue the reco value of the quantity to meaasure
                 *  @param  weight the nominal weight
                 *  @param  normalisation the sample normalisation
                 */
                void AddSelectedBackgroundEvent(const std::shared_ptr<Event> &pEvent, const AnalysisHelper::SampleType sampleType, const float &recoValue, const float weight, const float normalisation);
                
                /**
                 *  @brief  Add a BNB data event that has been selected
                 *
                 *  @param  pEvent the event
                 *  @param  recoValue the reco value of the quantity to measure
                 */
                void AddSelectedBNBDataEvent(const std::shared_ptr<Event> &pEvent, const float &recoValue);

                /**
                 *  @brief  Get the total weight of signal events per bin
                 *
                 *  @param  eventsPerBin the events per bin
                 */
//                void CountSignalEventsPerBin(std::vector<float> &eventsPerBin) const;
                
                /**
                 *  @brief  Count the selected BNB data events per reco bin
                 *
                 *  @return the bin contents
                 */
                std::vector<float> CountSelectedBNBDataEventsPerRecoBin();
                
                /**
                 *  @brief  Count the selected MC (as fake data) events per reco bin
                 *
                 *  @return the bin contents
                 */
                std::vector<float> CountSelectedMCEventsPerRecoBin();
                
                /**
                 *  @brief  Count the selected background events per reco bin
                 *
                 *  @return the bin contents
                 */
                std::vector<float> CountSelectedBackgroundEventsPerRecoBin();
                
                /**
                 *  @brief  Count the signal events per true bin
                 *
                 *  @return the bin contents
                 */
                std::vector<float> CountSignalEventsPerTrueBin();
                
                /**
                 *  @brief  Count the selected signal events per true bin
                 *
                 *  @return the bin contents
                 */
                std::vector<float> CountSelectedSignalEventsPerTrueBin();
                
                /**
                 *  @brief  Count the selected signal events per 2D reco-true bin
                 *
                 *  @return the bin contents
                 */
                std::vector< std::vector<float> > CountSelectedSignalEventsPerRecoTrueBin();

                /**
                 *  @brief  Get the efficiency of the selection per true bin
                 *
                 *  @return the bin efficiencies
                 */
                std::vector<float> GetEfficiencyPerTrueBin();
                
                /**
                 *  @brief  Get the smaering matrix [reco][true]
                 *
                 *  @return the smearing matrix
                 */
                std::vector< std::vector<float> > GetSmearingMatrix();

                /**
                 *  @brief  Get the smeared efficiency in reco-space
                 *
                 *  @return the reco bin efficiencies
                 */
                std::vector<float> GetSmearedEfficiencyPerRecoBin();

                /**
                 *  @brief  Get the cross-section from MC using true values
                 *
                 *  @param  includeUnderOverFlowBins if we should include underflow/overflow bins in the output (if false, dummy values are returned)
                 *
                 *  @return the MC true cross section
                 */
                std::vector<float> GetMCTrueCrossSectionPerTrueBin(const bool includeUnderOverFlowBins = false);

                /**
                 *  @brief  Get the cross-section from MC using smeared true values
                 *
                 *  @return the smeared MC true cross section
                 */
                std::vector<float> GetSmearedMCTrueCrossSectionPerRecoBin();

                /**
                 *  @brief  Get the MC stats uncertainty for the cross-section calculated from smeared true values
                 *
                 *  @return the smeared MC true cross section MC stat uncertainty
                 */
                std::vector<float> GetSmearedMCTrueCrossSectionMCStatUncertaintyPerRecoBin();
                
                /**
                 *  @brief  Get the forward folded cross-section in reco bins
                 *
                 *  @param  useRealData if we should use real BNB data (if false, MC is used as fake data)
                 *
                 *  @return the cross section per reco bin
                 */
                std::vector<float> GetCrossSectionPerRecoBin(const bool useRealData);
                
                /**
                 *  @brief  Get the stat uncertainty on the cross-section accounting only for the data (or fake data)
                 *
                 *  @param  useRealData if we should use real BNB data (if false, MC is used as fake data)
                 *
                 *  @return the cross section uncertainty per reco bin
                 */
                std::vector<float> GetCrossSectionStatUncertaintyPerRecoBin(const bool useRealData);
                
                /**
                 *  @brief  Get the stat uncertainty on the cross-section accounting only for the MC stats (via backgrounds, efficiency & smearing)
                 *          This function uses the bootstrap method to estiamte the uncertainty
                 *
                 *  @param  useRealData if we should use real BNB data (if false, MC is used as fake data)
                 *
                 *  @return the cross section uncertainty per reco bin
                 */
                std::vector<float> GetCrossSectionMCStatUncertaintyPerRecoBin(const bool useRealData);

                /**
                 *  @brief  Get the MC statistical uncertainty on the smearing matrix elements
                 *
                 *  @return the uncertainties on the smearing matrix elements
                 */
                std::vector< std::vector<float> > GetSmearingMatrixMCStatUncertainty();
                
                /**
                 *  @brief  Print the systematic universes
                 */
                void PrintSystematicUniverses() const;

                /**
                 *  @brief  Print the contents of the bins
                 *
                 *  @param  outputFileNamePrefix the prefix of the name of the output file we should save the tables to
                 */
                void PrintBinContents(const std::string &outputFileNamePrefix);

                /**
                 *  @brief  Make the relevant plots and save them to disk
                 *
                 *  @param  outputFileNamePrefix the prefix of the name of the output files
                 *  @param  includeSystematicUncertainties if we should run the systematic uncertainties (can be slow) and include them in the plots
                 */
                void MakePlots(const std::string &outputFileNamePrefix, const bool includeSystematicUncertainties);

                /**
                 *  @brief  Get a histogram filled from a set of bin values and uncertainties
                 *
                 *  @param  binValues the values of each bin
                 *  @param  binUncertainties the uncertainties of each bin
                 *  @param  includeUnderflow include the underflow bin (if it exists)
                 *  @param  includeOverflow include the overflow bin (if it exists)
                 *
                 *  @return the histogram filled with the specified values
                 */
                std::shared_ptr<TH1F> GetHistogram(const std::vector<float> &binValues, const std::vector<float> &binUncertainties, const bool includeUnderflow = false, const bool includeOverflow = false) const;

                /**
                 *  @brief  Get a 2D histogram filled with the supplied values
                 *
                 *  @param  binValues the bin values
                 *
                 *  @return the histogram
                 */
                std::shared_ptr<TH2F> GetHistogram(const std::vector< std::vector<float> > &binValues) const;
                
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
                 *  @brief  Setup branches in the output tree for a given systematic
                 *
                 *  @param  systematicParam the name of the systematic parameter
                 *  @param  nUniverses the number of universes
                 */
                void SetupSystematicBranches(const std::string &systematicParam, const unsigned int nUniverses);

                /**
                 *  @brief  Get the branch name for a given systematic parameter in a given universe
                 *
                 *  @param  systematicParam the parameter name
                 *  @param  universeIndex the universe index
                 *
                 *  @return the branch name
                 */
                std::string GetSystematicBranchName(const std::string &systematicParam, const unsigned int universeIndex) const;
                
                /**
                 *  @brief  Disable all systematic parameter branches for reading
                 */
                void DisableAllSystematicBranches();

                /**
                 *  @brief  Enable a branch for reading for a given systematic parameter
                 *
                 *  @param  systematicParam the parameter name
                 *  @param  universeIndex the universe index
                 */
                void EnableSystematicBranch(const std::string &systematicParam, const unsigned int universeIndex);

                /**
                 *  @brief  Get reco-true value pairs for selected "MC" events
                 *
                 *  @param  recoTrueValuePairs the output vector. The first element is a pair of reco & true values, the second element is the weight
                 */
                void GetSelectedMCRecoTrueValuePairs(std::vector< std::pair<std::pair<float, float>, float> > &recoTrueValuePairs);

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
                 *  @brief  Determine if a given bin is underflow
                 *
                 *  @param  index the bin index
                 *
                 *  @return boolean, true if underflow bin
                 */
                bool IsUnderflowBin(const unsigned int index) const;
                
                /**
                 *  @brief  Determine if a given bin is overflow
                 *
                 *  @param  index the bin index
                 *
                 *  @return boolean, true if overflow bin
                 */
                bool IsOverflowBin(const unsigned int index) const;
               
                /**
                 *  @brief  Get the bin widths
                 *
                 *  @return the bin widths
                 */
                std::vector<float> GetBinWidths() const;

                /**
                 *  @brief  Get the weight (accounting for sample normalisation) of an event after applying the supplied systematic parameters in the given universe
                 *
                 *  @param  event the event
                 *  @param  systematicParams the names of the systematic parameters to apply (can be empty, in which case the nominal universe is used)
                 *  @param  universeIndex the index of the universe to use
                 *
                 *  @return the event weight
                 */
                float GetEventWeight(const OutputEvent &event, const std::vector<std::string> &systematicParams, const unsigned int universeIndex) const;

                /**
                 *  @brief  Count the events per bin (using reconstructed values) passing the supplied critiera
                 *          The result is reweighted using the requested systematic parameters in the requested universe
                 *
                 *  @param  criteria the criteria to count an event
                 *  @param  systematicParams the names of the systematic parameters to apply (can be empty, in which case the nominal universe is used)
                 *  @param  universeIndex the index of the universe to use
                 *
                 *  @return the total weight of events passing the criteira per bin
                 */
                std::vector<float> CountEventsPerRecoBinWithCriteria(const std::function<bool(const OutputEvent &)> &criteria, const std::vector<std::string> &systematicParams, const unsigned int universeIndex);
                
                /**
                 *  @brief  Count the events per bin (using true values) passing the supplied critiera
                 *          The result is reweighted using the requested systematic parameters in the requested universe
                 *
                 *  @param  criteria the criteria to count an event
                 *  @param  systematicParams the names of the systematic parameters to apply (can be empty, in which case the nominal universe is used)
                 *  @param  universeIndex the index of the universe to use
                 *
                 *  @return the total weight of events passing the criteira per bin
                 */
                std::vector<float> CountEventsPerTrueBinWithCriteria(const std::function<bool(const OutputEvent &)> &criteria, const std::vector<std::string> &systematicParams, const unsigned int universeIndex);

                /**
                 *  @brief  Count the events per [reco][true] bin passing the supplied critiera
                 *          The result is reweighted using the requested systematic parameters in the requested universe
                 *
                 *  @param  criteria the criteria to count an event
                 *  @param  systematicParams the names of the systematic parameters to apply (can be empty, in which case the nominal universe is used)
                 *  @param  universeIndex the index of the universe to use
                 *
                 *  @return the total weight of events passing the criteira per bin
                 */
                std::vector< std::vector<float> > CountEventsPerRecoTrueBinWithCriteria(const std::function<bool(const OutputEvent &)> &criteria, const std::vector<std::string> &systematicParams, const unsigned int universeIndex);

                /**
                 *  @brief  Count the selected background events per reco bin - applying systematic parameters
                 *
                 *  @param  systematicParams the systematic parameters to apply
                 *  @param  universeIndex the universe to use
                 *
                 *  @return the bin contents
                 */
                std::vector<float> CountSelectedBackgroundEventsPerRecoBin(const std::vector<std::string> &systematicParams, const unsigned int universeIndex);
                
                /**
                 *  @brief  Count the signal events per true bin - applying systematic parameters
                 *
                 *  @param  systematicParams the systematic parameters to apply
                 *  @param  universeIndex the universe to use
                 *
                 *  @return the bin contents
                 */
                std::vector<float> CountSignalEventsPerTrueBin(const std::vector<std::string> &systematicParams, const unsigned int universeIndex);
                
                /**
                 *  @brief  Count the selected signal events per true bin - applying systematic parameters
                 *
                 *  @param  systematicParams the systematic parameters to apply
                 *  @param  universeIndex the universe to use
                 *
                 *  @return the bin contents
                 */
                std::vector<float> CountSelectedSignalEventsPerTrueBin(const std::vector<std::string> &systematicParams, const unsigned int universeIndex);
                
                /**
                 *  @brief  Count the selected signal events per 2D reco-true bin - applying systematic parameters
                 *
                 *  @param  systematicParams the systematic parameters to apply
                 *  @param  universeIndex the universe to use
                 *
                 *  @return the bin contents
                 */
                std::vector< std::vector<float> > CountSelectedSignalEventsPerRecoTrueBin(const std::vector<std::string> &systematicParams, const unsigned int universeIndex);
                
                /**
                 *  @brief  Get the efficiency of the selection per true bin - applying systematic parameters
                 *
                 *  @param  systematicParams the systematic parameters to apply
                 *  @param  universeIndex the universe to use
                 *
                 *  @return the bin efficiencies
                 */
                std::vector<float> GetEfficiencyPerTrueBin(const std::vector<std::string> &systematicParams, const unsigned int universeIndex);

                /**
                 *  @brief  Get the smearing matrix [reco][true] - applying systematic parameters
                 *
                 *  @param  systematicParams the systematic parameters to apply
                 *  @param  universeIndex the universe to use
                 *
                 *  @return the smearing matrix
                 */
                std::vector< std::vector<float> > GetSmearingMatrix(const std::vector<std::string> &systematicParams, const unsigned int universeIndex);

                /**
                 *  @brief  Smear a set of input bins (true-space) into reco-space using the supplied smearing matrix
                 *
                 *  @param  bins the input true-space bins
                 *  @param  smearingMatrix the smearing matrix
                 *
                 *  @return the output reco-space smeared bins
                 */
                std::vector<float> Smear(const std::vector<float> &bins, const std::vector< std::vector<float> > &smearingMatrix) const;
                
                /**
                 *  @brief  Get the smeared efficiencies in reco-space - applying systematic parameters
                 *
                 *  @param  systematicParams the systematic parameters to apply
                 *  @param  universeIndex the universe to use
                 *
                 *  @return the reco efficiency per bin
                 */
                std::vector<float> GetSmearedEfficiencyPerRecoBin(const std::vector<std::string> &systematicParams, const unsigned int universeIndex);
                
                /**
                 *  @brief  Get the cross-section from MC using true values - applying systematic parameters
                 *
                 *  @param  systematicParams the systematic parameters to apply
                 *  @param  universeIndex the universe to use
                 *  @param  includeUnderOverFlowBins if we should include underflow/overflow bins in the output (if false, dummy values are returned)
                 *
                 *  @return the MC true cross section
                 */
                std::vector<float> GetMCTrueCrossSectionPerTrueBin(const std::vector<std::string> &systematicParams, const unsigned int universeIndex, const bool includeUnderOverFlowBins = false);

                /**
                 *  @brief  Get the cross-section from MC using smeared true values - applying systematic parameters
                 *
                 *  @param  systematicParams the systematic parameters to apply
                 *  @param  universeIndex the universe to use
                 *
                 *  @return the smeared MC true cross section
                 */
                std::vector<float> GetSmearedMCTrueCrossSectionPerRecoBin(const std::vector<std::string> &systematicParams, const unsigned int universeIndex);

                /**
                 *  @brief  Get the forward folded cross-section in reco bins - applying systematic parameters
                 *
                 *  @param  systematicParams the systematic parameters to apply
                 *  @param  universeIndex the universe to use
                 *  @param  useRealData if we should use real BNB data (if false, MC is used as fake data)
                 *
                 *  @return the cross section per reco bin
                 */
                std::vector<float> GetCrossSectionPerRecoBin(const std::vector<std::string> &systematicParams, const unsigned int universeIndex, const bool useRealData);

                /**
                 *  @brief  Add two vectors in quadrature element-wise: sqrt(a^2 + b^2)
                 *
                 *  @param  a the first vector
                 *  @param  b the second vector
                 *
                 *  @return the sum in quadrature
                 */
                std::vector<float> AddInQuadrature(const std::vector<float> &a, const std::vector<float> &b) const;

                /**
                 *  @brief  Set the y-range of a vector of histograms
                 *
                 *  @param  hists the histograms to modify
                 *  @param  capAtZero is we should cap the lower y-value at zero
                 */
                void SetHistogramYRanges(std::vector< std::shared_ptr<TH1F> > &hists, const bool capAtZero) const;

                std::string        m_fileName; ///< The output file name
                TFile             *m_pFile;    ///< The output file (never explicitly written, but required by ROOT to fall back on disk if processing is too memory intensive)
                TTree             *m_pTree;    ///< The output tree 
                                  
                float              m_min;      ///< The minimum possible value of the quantity
                float              m_max;      ///< The maximum possible value of the quantity

                std::vector<float> m_binEdges;     ///< The bin edges
                bool               m_hasUnderflow; ///< If we have an underflow bin
                bool               m_hasOverflow;  ///< If we have an overflow bin

                bool               m_useAbsPdg;                ///< Use absolute PDGs when classifying events
                bool               m_countProtonsInclusively;  ///< Count protons inclusively when classifying events
                                  
                float              m_nTargets;        ///< The number of target particles / 10^31
                float              m_integratedFlux;  ///< The integrated flux / 10^11 cm^-2
                                  
                OutputEvent        m_outputEvent;         ///< The output event struture to bind to the trees

                std::default_random_engine m_generator;           ///< Random number generator
                unsigned int               m_nBootstrapUniverses; ///< The number of universes to use when finding the statistical uncertainty
                
                static unsigned int m_histIndex;   ///< An index used to ensure all histogram names are unique
        };
};

// -----------------------------------------------------------------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------------------------------------------------------------

unsigned int CrossSectionHelper::XSec::m_histIndex = 0u;

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
    m_integratedFlux(1.26816), // TODO make configurable
    m_nBootstrapUniverses(100u) // TODO make configurable
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

void CrossSectionHelper::XSec::PrintSystematicUniverses() const
{
    std::cout << "Systematic universes:" << std::endl;
    for (const auto &entry : m_outputEvent.m_systWeights)
    {
        std::cout << "  - " << entry.first << ", " << entry.second.size() << std::endl;
    }
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void CrossSectionHelper::XSec::DisableAllSystematicBranches()
{
    m_pTree->SetBranchStatus("systWeight_*", false);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void CrossSectionHelper::XSec::EnableSystematicBranch(const std::string &systematicParam, const unsigned int universeIndex)
{
    m_pTree->SetBranchStatus(this->GetSystematicBranchName(systematicParam, universeIndex).c_str(), true);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void CrossSectionHelper::XSec::GetSelectedMCRecoTrueValuePairs(std::vector< std::pair<std::pair<float, float>, float> > &recoTrueValuePairs)
{
    // For speed turn off the systematic parameters for reading
    this->DisableAllSystematicBranches();

    // Extract the reco values from the tree
    for (unsigned int i = 0; i < m_pTree->GetEntries(); ++i)
    {
        m_pTree->GetEntry(i);

        // We only care about "MC"
        if (m_outputEvent.m_sampleType == AnalysisHelper::DataBNB)
            continue;
        
        // Only use selected events
        if (!m_outputEvent.m_isSelected)
            continue;

        const auto weight = m_outputEvent.m_weight * m_outputEvent.m_sampleNorm;
        recoTrueValuePairs.emplace_back(std::pair<float, float>(m_outputEvent.m_recoValue, m_outputEvent.m_trueValue), weight);
    }
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void CrossSectionHelper::XSec::SetBinsAuto(const float lower, const float upper, const unsigned int targetEntriesPerBin, const float targetSmearingDiagonal)
{
    // Get the reco & true values
    std::vector< std::pair<std::pair<float, float>, float> > recoTrueValuePairs;
    this->GetSelectedMCRecoTrueValuePairs(recoTrueValuePairs);

    // Check we have at least 2 entries
    const auto nEntries = recoTrueValuePairs.size();
    if (nEntries < 2)
        throw std::logic_error("XSec::SetBinsAuto - there are fewer than 2 entries!");

    // Sort the entries by reco-range
    std::sort(recoTrueValuePairs.begin(), recoTrueValuePairs.end(), [] (const auto &a, const auto &b) {
        return a.first.first < b.first.first;
    });

    // Coun the number of entries with reco-range in the specified limits
    std::vector< std::pair<std::pair<float, float>, float> > entriesInRange;
    const auto nEntriesInRange = std::count_if(recoTrueValuePairs.begin(), recoTrueValuePairs.end(), [&] (const auto &x) {
        return x.first.first >= lower && x.first.first <= upper;
    });
    
    // Check we have enough entries
    if (nEntriesInRange < 2)
        throw std::logic_error("XSec::SetBinsAuto - there are fewer than 2 entries in the specified range");
    
    if (nEntriesInRange < targetEntriesPerBin)
        throw std::logic_error("XSec::SetBinsAuto - there aren't enough entries in the specified range to meet the target entries per bin");

    // Get the entries with reco-range in the specified limits
    entriesInRange.reserve(nEntriesInRange);
    std::copy_if(recoTrueValuePairs.begin(), recoTrueValuePairs.end(), std::back_inserter(entriesInRange), [&] (const auto &x) {
        return x.first.first >= lower && x.first.first <= upper;
    });

    // Now generate the bin edges
    std::vector<float> binEdges(1, lower);

    // Variables defining the current bin
    float lowerEdge = binEdges.back();
    unsigned int upperEdgeIndex = 0u;

    while (true)
    {
        // Increment the upper edge index and check if we have reached the end
        if (++upperEdgeIndex >= nEntriesInRange)
            break;

        const float upperEdge = entriesInRange.at(upperEdgeIndex).first.first;

        // Count the number of entries in the current bin
        auto nEntriesInBinReco = 0.f;
        auto nEntriesInBinTrue = 0.f;
        auto nEntriesInBinRecoAndTrue = 0.f;

        for (const auto &entry : entriesInRange)
        {
            const auto recoVal = entry.first.first;
            const auto trueVal = entry.first.second;
            const auto weight = entry.second;

            const auto inBinReco = (recoVal >= lowerEdge && recoVal < upperEdge);
            const auto inBinTrue = (trueVal >= lowerEdge && trueVal < upperEdge);

            nEntriesInBinReco += inBinReco ? weight : 0.f;
            nEntriesInBinTrue += inBinTrue ? weight : 0.f;
            nEntriesInBinRecoAndTrue += (inBinReco && inBinTrue) ? weight : 0.f;
        }

        // If we don't have enough entries then grow the bin
        if (nEntriesInBinReco < static_cast<float>(targetEntriesPerBin))
            continue;
        
        // If we don't have any true entries in the is bin, then we can't calculate the smearing matrix element
        if (nEntriesInBinTrue <= std::numeric_limits<float>::epsilon())
            continue;

        // Check that the amount of smearing is sufficiently low
        const auto smearingMatrixElement = nEntriesInBinRecoAndTrue / nEntriesInBinTrue;
        if (smearingMatrixElement < targetSmearingDiagonal)
            continue;

        // We have a valid bin!
        lowerEdge = upperEdge;
        binEdges.push_back(upperEdge);
    }

    if (binEdges.size() < 2)
        throw std::logic_error("XSec::SetBinsAuto - no bins automatically generated!");

    // Stretch the last bin to fit the desired range
    binEdges.at(binEdges.size() - 1) = upper;
    
    // Set the bin eges and we are done!
    this->SetBins(binEdges);    
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void CrossSectionHelper::XSec::FillEvent(const AnalysisHelper::SampleType sampleType, const float sampleNorm, const float weight, const float &recoValue, const float &trueValue, const bool isSignal, const bool isSelected, const std::string &classification)
{
    // Set the event parameters
    m_outputEvent.m_sampleType = sampleType;
    m_outputEvent.m_sampleNorm = sampleNorm;
    m_outputEvent.m_weight = weight;
    m_outputEvent.m_recoValue = recoValue;
    m_outputEvent.m_trueValue = trueValue;
    m_outputEvent.m_isSignal = isSignal;
    m_outputEvent.m_isSelected = isSelected;
    m_outputEvent.m_classification = classification;

    // Set the systematic weights
    // ATTN here we only have one "systematic" which is for the MC stats. In the future this can be extended by using other systematic event
    // weights. The m_systWeights object is a map from the parameter name to the weights of the event in each index, i.e. [paramName][index]
    // Here we extend the map for each new parameter & universe index, and also add a branch to the tree using a consistent naming
    // convention. In this way we can disable all branches and only load the weights for the parameters & universe that we care about - this
    // is done for performance reasons. To extend to other systematics, follow the procedure used here for the MC stats.
    
    // If this is the first event, then set up the branches
    // ATTN this is done here instead of the constructor so that future systmatics can be implemented without the need to know the names /
    // number of universes at compile time
    if (m_pTree->GetEntries() == 0)
    {
        this->SetupSystematicBranches("stat", m_nBootstrapUniverses);
    }

    // Populate the MC stats systematic weights using the Poisson bootstrap method
    std::poisson_distribution<int> poisson(1.f);
    auto &mcStatWeights = m_outputEvent.m_systWeights.at("stat");
    for (unsigned int iUni = 0; iUni < m_nBootstrapUniverses; ++iUni)
    {
        mcStatWeights.at(iUni) = static_cast<float>(poisson(m_generator));
    }

    // Store this event!
    m_pTree->Fill();
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void CrossSectionHelper::XSec::SetupSystematicBranches(const std::string &systematicParam, const unsigned int nUniverses)
{
    // Make map entries for this systematic parameter
    if (!m_outputEvent.m_systWeights.emplace(systematicParam, std::vector<float>(nUniverses, -std::numeric_limits<float>::max())).second)
    {
        throw std::logic_error("XSec::SetupSystematicBranches - parameter \"" + systematicParam + "\" already has been set up");
    }

    // Make a new branch for each universe and bind it to the relevant map entry
    for (unsigned int iUni = 0; iUni < nUniverses; ++iUni)
    {
        const auto branchName = this->GetSystematicBranchName(systematicParam, iUni);
        m_pTree->Branch(branchName.c_str(), &m_outputEvent.m_systWeights.at(systematicParam).at(iUni));
    }
}

// -----------------------------------------------------------------------------------------------------------------------------------------

std::string CrossSectionHelper::XSec::GetSystematicBranchName(const std::string &systematicParam, const unsigned int universeIndex) const
{
    return ("systWeight_" + systematicParam + "_" + std::to_string(universeIndex));
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
   
    throw std::logic_error("XSec::GetBinIndex - input value is out of bounds");
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
        return "UF";

    if (index == nBins - 1 && m_hasOverflow)
        return "OF";

    return std::to_string(index - (m_hasUnderflow ? 1 : 0));
}

// -----------------------------------------------------------------------------------------------------------------------------------------

bool CrossSectionHelper::XSec::IsUnderOverFlowBin(const unsigned int index) const
{
    return (this->IsUnderflowBin(index) || this->IsOverflowBin(index));
}

// -----------------------------------------------------------------------------------------------------------------------------------------

bool CrossSectionHelper::XSec::IsUnderflowBin(const unsigned int index) const
{
    if (m_binEdges.size() < 2)
        throw std::logic_error("XSec::IsUnderflowBin - binning hasn't been set");

    const auto nBins = m_binEdges.size() - 1;

    return (index == 0 && m_hasUnderflow);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

bool CrossSectionHelper::XSec::IsOverflowBin(const unsigned int index) const
{
    if (m_binEdges.size() < 2)
        throw std::logic_error("XSec::IsOverflowBin - binning hasn't been set");

    const auto nBins = m_binEdges.size() - 1;

    return (index == nBins - 1 && m_hasOverflow);
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

std::vector<float> CrossSectionHelper::XSec::CountEventsPerRecoBinWithCriteria(const std::function<bool(const OutputEvent &)> &criteria, const std::vector<std::string> &systematicParams, const unsigned int universeIndex)
{
    if (m_binEdges.size() < 2)
        throw std::logic_error("XSec::CountEventsPerRecoBinWithCriteria - binning hasn't been set");

    const auto nBins = m_binEdges.size() - 1;
    std::vector<float> binTotals(nBins, 0.f);
    
    // For speed turn off the systematic parameters we don't need for reading
    this->DisableAllSystematicBranches();
    for (const auto &paramName : systematicParams)
        this->EnableSystematicBranch(paramName, universeIndex);

    // Extract the reco values from the tree
    for (unsigned int i = 0; i < m_pTree->GetEntries(); ++i)
    {
        m_pTree->GetEntry(i);

        // Only count events passing the criteria
        if (!criteria(m_outputEvent))
            continue;

        const auto weight = this->GetEventWeight(m_outputEvent, systematicParams, universeIndex);
        const auto binIndex = this->GetBinIndex(m_outputEvent.m_recoValue);
        binTotals.at(binIndex) += weight;
    }

    return binTotals;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

std::vector<float> CrossSectionHelper::XSec::CountEventsPerTrueBinWithCriteria(const std::function<bool(const OutputEvent &)> &criteria, const std::vector<std::string> &systematicParams, const unsigned int universeIndex)
{
    if (m_binEdges.size() < 2)
        throw std::logic_error("XSec::CountEventsPerTrueBinWithCriteria - binning hasn't been set");

    const auto nBins = m_binEdges.size() - 1;
    std::vector<float> binTotals(nBins, 0.f);
    
    // For speed turn off the systematic parameters we don't need for reading
    this->DisableAllSystematicBranches();
    for (const auto &paramName : systematicParams)
        this->EnableSystematicBranch(paramName, universeIndex);

    // Extract the reco values from the tree
    for (unsigned int i = 0; i < m_pTree->GetEntries(); ++i)
    {
        m_pTree->GetEntry(i);

        // Only count events passing the criteria
        if (!criteria(m_outputEvent))
            continue;
        
        if (!m_outputEvent.m_isSignal)
            throw std::logic_error("XSec::CountEventsPerTrueBinWithCriteria - event passed the criteria but is not signal, no truth information available");

        const auto weight = this->GetEventWeight(m_outputEvent, systematicParams, universeIndex);
        const auto binIndex = this->GetBinIndex(m_outputEvent.m_trueValue);
        binTotals.at(binIndex) += weight;
    }

    return binTotals;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

std::vector< std::vector<float> > CrossSectionHelper::XSec::CountEventsPerRecoTrueBinWithCriteria(const std::function<bool(const OutputEvent &)> &criteria, const std::vector<std::string> &systematicParams, const unsigned int universeIndex)
{
    if (m_binEdges.size() < 2)
        throw std::logic_error("XSec::CountEventsPerRecoTrueBinWithCriteria - binning hasn't been set");

    const auto nBins = m_binEdges.size() - 1;

    // Setup the output matrix
    std::vector< std::vector<float> > binTotals;
    for (unsigned int i = 0; i < nBins; ++i)
        binTotals.emplace_back(nBins, 0.f);
    
    // For speed turn off the systematic parameters we don't need for reading
    this->DisableAllSystematicBranches();
    for (const auto &paramName : systematicParams)
        this->EnableSystematicBranch(paramName, universeIndex);

    // Extract the reco values from the tree
    for (unsigned int i = 0; i < m_pTree->GetEntries(); ++i)
    {
        m_pTree->GetEntry(i);

        // Only count events passing the criteria
        if (!criteria(m_outputEvent))
            continue;
        
        if (!m_outputEvent.m_isSignal)
            throw std::logic_error("XSec::CountEventsPerRecoTrueBinWithCriteria - event passed the criteria but is not signal, no truth information available");

        const auto weight = this->GetEventWeight(m_outputEvent, systematicParams, universeIndex);
        const auto trueBinIndex = this->GetBinIndex(m_outputEvent.m_trueValue);
        const auto recoBinIndex = this->GetBinIndex(m_outputEvent.m_recoValue);

        binTotals.at(recoBinIndex).at(trueBinIndex) += weight;
    }

    return binTotals;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

float CrossSectionHelper::XSec::GetEventWeight(const OutputEvent &event, const std::vector<std::string> &systematicParams, const unsigned int universeIndex) const
{
    // Get the nominal event weight
    const auto nominalWeight = event.m_weight * event.m_sampleNorm;

    // Re-weight the event using the supplied systematic parameters
    float reweightedValue = nominalWeight;

    const auto &systWeights = event.m_systWeights;
    for (const auto &param : systematicParams)
    {
        // Get the universe weights for this parameter
        const auto iter = systWeights.find(param);
        if (iter == systWeights.end())
            throw std::invalid_argument("XSec::GetEventWeight - Can't find requested systematic parameter: \"" + param + "\"");

        const auto &weights = iter->second;

        // Check we are trying to access a valid universe index
        const auto nUniverses = weights.size();
        if (universeIndex >= nUniverses)
            throw std::out_of_range("XSec::GetEventWeight - Systematic parameter: \"" + param + "\" has " + std::to_string(nUniverses) + " universes. Requested universe index, " + std::to_string(universeIndex) + " is out of bounds.");

        // Get the weight for this universe
        const auto universeWeight = weights.at(universeIndex);

        // Apply the weights
        reweightedValue *= universeWeight;
    }

    return reweightedValue;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

std::vector<float> CrossSectionHelper::XSec::CountSelectedBNBDataEventsPerRecoBin()
{
    return this->CountEventsPerRecoBinWithCriteria([](const OutputEvent &event) {

        return (event.m_sampleType == AnalysisHelper::DataBNB && event.m_isSelected);

    }, {}, 0u);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

std::vector<float> CrossSectionHelper::XSec::CountSelectedMCEventsPerRecoBin()
{
    return this->CountEventsPerRecoBinWithCriteria([](const OutputEvent &event) {

        return (event.m_sampleType != AnalysisHelper::DataBNB && event.m_isSelected);

    }, {}, 0u);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

std::vector<float> CrossSectionHelper::XSec::CountSelectedBackgroundEventsPerRecoBin()
{
    return this->CountSelectedBackgroundEventsPerRecoBin({}, 0u);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

std::vector<float> CrossSectionHelper::XSec::CountSelectedBackgroundEventsPerRecoBin(const std::vector<std::string> &systematicParams, const unsigned int universeIndex)
{
    return this->CountEventsPerRecoBinWithCriteria([](const OutputEvent &event) {

        return (event.m_sampleType != AnalysisHelper::DataBNB && !event.m_isSignal && event.m_isSelected);

    }, systematicParams, universeIndex);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

std::vector<float> CrossSectionHelper::XSec::CountSignalEventsPerTrueBin()
{
    return this->CountSignalEventsPerTrueBin({}, 0u);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

std::vector<float> CrossSectionHelper::XSec::CountSignalEventsPerTrueBin(const std::vector<std::string> &systematicParams, const unsigned int universeIndex)
{
    return this->CountEventsPerTrueBinWithCriteria([](const OutputEvent &event) {

        return event.m_isSignal;

    }, systematicParams, universeIndex);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

std::vector<float> CrossSectionHelper::XSec::CountSelectedSignalEventsPerTrueBin()
{
    return this->CountSelectedSignalEventsPerTrueBin({}, 0u);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

std::vector<float> CrossSectionHelper::XSec::CountSelectedSignalEventsPerTrueBin(const std::vector<std::string> &systematicParams, const unsigned int universeIndex)
{
    return this->CountEventsPerTrueBinWithCriteria([](const OutputEvent &event) {

        return (event.m_isSignal && event.m_isSelected);

    }, systematicParams, universeIndex);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

std::vector< std::vector<float> > CrossSectionHelper::XSec::CountSelectedSignalEventsPerRecoTrueBin()
{
    return this->CountSelectedSignalEventsPerRecoTrueBin({}, 0u);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

std::vector< std::vector<float> > CrossSectionHelper::XSec::CountSelectedSignalEventsPerRecoTrueBin(const std::vector<std::string> &systematicParams, const unsigned int universeIndex)
{
    return this->CountEventsPerRecoTrueBinWithCriteria([](const OutputEvent &event) {

        return (event.m_isSignal && event.m_isSelected);

    }, systematicParams, universeIndex);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

std::vector<float> CrossSectionHelper::XSec::GetEfficiencyPerTrueBin()
{
    return this->GetEfficiencyPerTrueBin({}, 0u);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

std::vector<float> CrossSectionHelper::XSec::GetEfficiencyPerTrueBin(const std::vector<std::string> &systematicParams, const unsigned int universeIndex)
{
    // Count the signal events per bin before and after the selection
    const auto numerators = this->CountSelectedSignalEventsPerTrueBin(systematicParams, universeIndex);
    const auto denominators = this->CountSignalEventsPerTrueBin(systematicParams, universeIndex);

    if (numerators.size() != denominators.size())
        throw std::logic_error("XSec::GetEfficiencyPerTrueBin - vector of numerators and denominators are of different sizes");

    const auto nBins = numerators.size();

    std::vector<float> efficiencies;
    for (unsigned int iBin = 0; iBin < nBins; ++iBin)
    {
        // We don't need the efficiency for the under/overflow bins so just push a dummy value
        if (this->IsUnderOverFlowBin(iBin))
        {
            efficiencies.push_back(-std::numeric_limits<float>::max());
            continue;
        }

        const auto denominator = denominators.at(iBin);
        if (denominator <= std::numeric_limits<float>::epsilon())
            throw std::logic_error("XSec::GetEfficiencyPerTrueBin - Found bin with no signal events, can't get efficiency");

        const auto numerator = numerators.at(iBin);

        efficiencies.push_back(numerator/denominator);
    }

    return efficiencies;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

std::vector< std::vector<float> > CrossSectionHelper::XSec::GetSmearingMatrix()
{
    return this->GetSmearingMatrix({}, 0u);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

std::vector< std::vector<float> > CrossSectionHelper::XSec::GetSmearingMatrix(const std::vector<std::string> &systematicParams, const unsigned int universeIndex)
{
    if (m_binEdges.size() < 2)
        throw std::logic_error("XSec::GetSmearingMatrix - binning hasn't been set");

    const auto nBins = m_binEdges.size() - 1;

    // Setup the empty matrix
    std::vector< std::vector<float> > smearingMatrix;
    for (unsigned int i = 0; i < nBins; ++i)
        smearingMatrix.emplace_back(nBins, 0.f);

 
    const auto confusionMatrix = this->CountSelectedSignalEventsPerRecoTrueBin(systematicParams, universeIndex);
    const auto denominators = this->CountSelectedSignalEventsPerTrueBin(systematicParams, universeIndex);

    for (unsigned int iReco = 0; iReco < denominators.size(); ++iReco)
    {
        for (unsigned int iTrue = 0; iTrue < denominators.size(); ++iTrue)
        {
            const auto denominator = denominators.at(iTrue);
            if (denominator <= std::numeric_limits<float>::epsilon())
                throw std::logic_error("XSec::GetSmearingMatrix - Found true bin with no selected signal events, can't get the smearing matrix element");

            smearingMatrix.at(iReco).at(iTrue) = confusionMatrix.at(iReco).at(iTrue) / denominator;
        }
    }

    return smearingMatrix;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

std::vector<float> CrossSectionHelper::XSec::Smear(const std::vector<float> &bins, const std::vector< std::vector<float> > &smearingMatrix) const
{
    if (bins.size() != smearingMatrix.size())
        throw std::invalid_argument("XSec::Smear - dimensions of input bins and smearing matrix do not match");

    std::vector<float> smearedBins;
    for (unsigned int iReco = 0; iReco < smearingMatrix.size(); ++iReco)
    {
        float val = 0.f;
            
        const auto &smearingValues = smearingMatrix.at(iReco);
        if (smearingValues.size() != smearingMatrix.size())
            throw std::invalid_argument("XSec::Smear - input smearing matrix is malformed");
    
        for (unsigned int iTrue = 0; iTrue < smearingMatrix.size(); ++iTrue)
        {
            val += bins.at(iTrue) * smearingValues.at(iTrue);
        }

        smearedBins.push_back(val);
    }

    return smearedBins;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

std::vector<float> CrossSectionHelper::XSec::GetSmearedEfficiencyPerRecoBin()
{
    return this->GetSmearedEfficiencyPerRecoBin({}, 0u);
}

// -----------------------------------------------------------------------------------------------------------------------------------------
                
std::vector<float> CrossSectionHelper::XSec::GetSmearedEfficiencyPerRecoBin(const std::vector<std::string> &systematicParams, const unsigned int universeIndex)
{
    // Count the signal events per bin before and after the selection
    const auto numerators = this->CountSelectedSignalEventsPerTrueBin(systematicParams, universeIndex);
    const auto denominators = this->CountSignalEventsPerTrueBin(systematicParams, universeIndex);

    if (numerators.size() != denominators.size())
        throw std::logic_error("XSec::GetSmearedEfficiencyPerRecoBin - vector of numerators and denominators are of different sizes");

    const auto nBins = numerators.size();

    // Smear the numerator and denominators into reco-space
    const auto smearingMatrix = this->GetSmearingMatrix(systematicParams, universeIndex);
    const auto smearedNumerators = this->Smear(numerators, smearingMatrix);
    const auto smearedDenominators = this->Smear(denominators, smearingMatrix);

    std::vector<float> efficiencies;
    for (unsigned int iBin = 0; iBin < nBins; ++iBin)
    {
        // We don't need the efficiency for the under/overflow bins so just push a dummy value
        if (this->IsUnderOverFlowBin(iBin))
        {
            efficiencies.push_back(-std::numeric_limits<float>::max());
            continue;
        }

        const auto denominator = smearedDenominators.at(iBin);
        if (denominator <= std::numeric_limits<float>::epsilon())
            throw std::logic_error("XSec::GetSmearedEfficiencyPerRecoBin - Found smeared bin with no events");

        const auto numerator = smearedNumerators.at(iBin);

        efficiencies.push_back(numerator/denominator);
    }

    return efficiencies;
}

// -----------------------------------------------------------------------------------------------------------------------------------------
                
std::vector<float> CrossSectionHelper::XSec::GetMCTrueCrossSectionPerTrueBin(const bool includeUnderOverFlowBins)
{
    return this->GetMCTrueCrossSectionPerTrueBin({}, 0u, includeUnderOverFlowBins);
}

// -----------------------------------------------------------------------------------------------------------------------------------------
                
std::vector<float> CrossSectionHelper::XSec::GetMCTrueCrossSectionPerTrueBin(const std::vector<std::string> &systematicParams, const unsigned int universeIndex, const bool includeUnderOverFlowBins)
{
    // Get the bin widths
    if (m_binEdges.size() < 2)
        throw std::logic_error("XSec::GetMCTrueCrossSectionPerTrueBin - binning hasn't been set");

    const auto nBins = m_binEdges.size() - 1;
    const auto binWidths = this->GetBinWidths();

    // Get the true number of signal events in each bin
    const auto signalEvents = this->CountSignalEventsPerTrueBin(systematicParams, universeIndex);
    
    std::vector<float> crossSections;
    for (unsigned int iTrue = 0; iTrue < nBins; ++iTrue)
    {
        if (this->IsUnderOverFlowBin(iTrue) && !includeUnderOverFlowBins)
        {
            crossSections.push_back(-std::numeric_limits<float>::max());
            continue;
        }

        const auto nSignal = signalEvents.at(iTrue);
        const auto binWidth = binWidths.at(iTrue);
    
        if (binWidth <= std::numeric_limits<float>::epsilon())
            throw std::invalid_argument("XSec::GetMCTrueCrossSectionPerTrueBin - bin width is zero");
    
        const auto crossSection = (nSignal) / (binWidth * m_nTargets * m_integratedFlux);
        crossSections.push_back(crossSection);
    }

    return crossSections;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

std::vector<float> CrossSectionHelper::XSec::GetSmearedMCTrueCrossSectionPerRecoBin()
{
    return this->GetSmearedMCTrueCrossSectionPerRecoBin({}, 0u);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

std::vector<float> CrossSectionHelper::XSec::GetSmearedMCTrueCrossSectionPerRecoBin(const std::vector<std::string> &systematicParams, const unsigned int universeIndex)
{
    // Get the bin widths
    if (m_binEdges.size() < 2)
        throw std::logic_error("XSec::GetSmearedMCTrueCrossSectionPerRecoBin - binning hasn't been set");

    const auto nBins = m_binEdges.size() - 1;
    const auto binWidths = this->GetBinWidths();
    
    // Get the MC true cross-section in true bins
    const auto trueCrossSections = this->GetMCTrueCrossSectionPerTrueBin(systematicParams, universeIndex, true);

    // Scale up by the bin-width to get quantities proportional to event rates
    std::vector<float> scaledTrueCrossSections;
    for (unsigned int iTrue = 0; iTrue < nBins; ++iTrue)
    {
        const auto trueCrossSection = trueCrossSections.at(iTrue);
        const auto binWidth = binWidths.at(iTrue);
        scaledTrueCrossSections.push_back(trueCrossSection * binWidth);
    }

    // Smear into reco-space
    const auto smearingMatrix = this->GetSmearingMatrix(systematicParams, universeIndex);
    const auto scaledSmearedCrossSections = this->Smear(scaledTrueCrossSections, smearingMatrix);

    // Scale back down by the bin widths
    std::vector<float> smearedCrossSections;
    for (unsigned int iReco = 0; iReco < nBins; ++iReco)
    {
        if (this->IsUnderOverFlowBin(iReco))
        {
            smearedCrossSections.push_back(-std::numeric_limits<float>::max());
            continue;
        }

        const auto scaledSmearedCrossSection = scaledSmearedCrossSections.at(iReco);
        const auto binWidth = binWidths.at(iReco);

        if (binWidth <= std::numeric_limits<float>::epsilon())
            throw std::invalid_argument("XSec::GetSmearedMCTrueCrossSectionPerRecoBin - bin width is zero");

        smearedCrossSections.push_back(scaledSmearedCrossSection / binWidth);
    }

    return smearedCrossSections;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

std::vector<float> CrossSectionHelper::XSec::GetSmearedMCTrueCrossSectionMCStatUncertaintyPerRecoBin()
{
    if (m_nBootstrapUniverses == 0)
        throw std::logic_error("XSec::GetSmearedMCTrueCrossSectionMCStatUncertaintyPerRecoBin - zero bootstrap universes have been generated");

    // Get the nominal cross-section
    const auto crossSectionsNom = this->GetSmearedMCTrueCrossSectionPerRecoBin();

    // Now get the cross-section in each universe of statistical variations
    const std::vector<std::string> systematicParams = {"stat"};
    std::vector<float> sumSquaredDiffs(crossSectionsNom.size(), 0.f);
    for (unsigned int iUni = 0; iUni < m_nBootstrapUniverses; ++iUni)
    {
        const auto crossSectionsUni = this->GetSmearedMCTrueCrossSectionPerRecoBin(systematicParams, iUni);

        for (unsigned int iReco = 0; iReco < crossSectionsNom.size(); ++iReco)
        {
            if (this->IsUnderOverFlowBin(iReco))
                continue;
            
            sumSquaredDiffs.at(iReco) += std::pow(crossSectionsUni.at(iReco) - crossSectionsNom.at(iReco), 2);
        }
    }
    
    std::vector<float> uncertainties;
    for (unsigned int iReco = 0; iReco < crossSectionsNom.size(); ++iReco)
    {
        if (this->IsUnderOverFlowBin(iReco))
        {
            uncertainties.push_back(-std::numeric_limits<float>::max());
            continue;
        }

        const auto standardDeviation = std::pow(sumSquaredDiffs.at(iReco) / static_cast<float>(m_nBootstrapUniverses), 0.5f); 
        uncertainties.push_back(standardDeviation);
    }

    return uncertainties;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

std::vector<float> CrossSectionHelper::XSec::GetCrossSectionPerRecoBin(const bool useRealData)
{
    return this->GetCrossSectionPerRecoBin({}, 0u, useRealData);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

std::vector<float> CrossSectionHelper::XSec::GetCrossSectionPerRecoBin(const std::vector<std::string> &systematicParams, const unsigned int universeIndex, const bool useRealData)
{
    // Get the bin widths
    if (m_binEdges.size() < 2)
        throw std::logic_error("XSec::GetCrossSectionPerRecoBin - binning hasn't been set");

    const auto nBins = m_binEdges.size() - 1;
    const auto binWidths = this->GetBinWidths();

    // Get the selected events
    const auto selectedEvents = useRealData ? this->CountSelectedBNBDataEventsPerRecoBin() : 
                                              this->CountSelectedMCEventsPerRecoBin();

    // Get the backgrounds
    const auto backgrounds = this->CountSelectedBackgroundEventsPerRecoBin(systematicParams, universeIndex);

    // Get the efficiencies
    const auto efficiencies = this->GetSmearedEfficiencyPerRecoBin(systematicParams, universeIndex);

    std::vector<float> crossSections;
    for (unsigned int iReco = 0; iReco < nBins; ++iReco)
    {
        // We don't need the cross-section in underflow/overflow bins so just use a dummy value
        if (this->IsUnderOverFlowBin(iReco))
        {
            crossSections.push_back(-std::numeric_limits<float>::max());
            continue;
        }

        const auto nTotal = selectedEvents.at(iReco);
        const auto nBackground = backgrounds.at(iReco);
        const auto binWidth = binWidths.at(iReco);
        const auto efficiency = efficiencies.at(iReco);
    
        if (binWidth <= std::numeric_limits<float>::epsilon())
            throw std::invalid_argument("XSec::GetCrossSectionPerRecoBin - bin width is zero");
    
        if (efficiency <= std::numeric_limits<float>::epsilon())
            throw std::invalid_argument("XSec::GetCrossSectionPerRecoBin - efficiency is zero");

        const auto crossSection = (nTotal - nBackground) / (efficiency * binWidth * m_nTargets * m_integratedFlux);
        crossSections.push_back(crossSection);
    }

    return crossSections;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

std::vector<float> CrossSectionHelper::XSec::GetCrossSectionStatUncertaintyPerRecoBin(const bool useRealData)
{
    // Get the bin widths
    if (m_binEdges.size() < 2)
        throw std::logic_error("XSec::GetCrossSectionStatUncertaintyPerRecoBin - binning hasn't been set");

    const auto nBins = m_binEdges.size() - 1;
    const auto binWidths = this->GetBinWidths();

    // Get the selected events
    const auto selectedEvents = useRealData ? this->CountSelectedBNBDataEventsPerRecoBin() : 
                                              this->CountSelectedMCEventsPerRecoBin();

    // Get the efficiencies in the nominal universe
    const auto efficiencies = this->GetSmearedEfficiencyPerRecoBin({}, 0u);

    std::vector<float> uncertainties;
    for (unsigned int iReco = 0; iReco < nBins; ++iReco)
    {
        if (this->IsUnderOverFlowBin(iReco))
        {
            uncertainties.push_back(-std::numeric_limits<float>::max());
            continue;
        }
        
        const auto nTotal = selectedEvents.at(iReco);
        const auto binWidth = binWidths.at(iReco);
        const auto efficiency = efficiencies.at(iReco);
    
        if (binWidth <= std::numeric_limits<float>::epsilon())
            throw std::invalid_argument("XSec::GetCrossSectionStatUncertaintyPerRecoBin - bin width is zero");
    
        if (efficiency <= std::numeric_limits<float>::epsilon())
            throw std::invalid_argument("XSec::GetCrossSectionStatUncertaintyPerRecoBin - efficiency is zero");

        const auto nTotalErr = AnalysisHelper::GetCountUncertainty(nTotal);
        const auto scaleFactor = 1.f / (efficiency * binWidth * m_nTargets * m_integratedFlux);

        const auto uncertainty = scaleFactor * nTotalErr;
        uncertainties.push_back(uncertainty);
    }

    return uncertainties;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

std::vector<float> CrossSectionHelper::XSec::GetCrossSectionMCStatUncertaintyPerRecoBin(const bool useRealData)
{
    if (m_nBootstrapUniverses == 0)
        throw std::logic_error("XSec::GetCrossSectionMCStatUncertaintyPerRecoBin - zero bootstrap universes have been generated");

    // Get the nominal cross-section
    const auto crossSectionsNom = this->GetCrossSectionPerRecoBin(useRealData);

    // Now get the cross-section in each universe of statistical variations
    const std::vector<std::string> systematicParams = {"stat"};
    std::vector<float> sumSquaredDiffs(crossSectionsNom.size(), 0.f);
    for (unsigned int iUni = 0; iUni < m_nBootstrapUniverses; ++iUni)
    {
        const auto crossSectionsUni = this->GetCrossSectionPerRecoBin(systematicParams, iUni, useRealData);

        for (unsigned int iReco = 0; iReco < crossSectionsNom.size(); ++iReco)
        {
            if (this->IsUnderOverFlowBin(iReco))
                continue;
            
            sumSquaredDiffs.at(iReco) += std::pow(crossSectionsUni.at(iReco) - crossSectionsNom.at(iReco), 2);
        }
    }
    
    std::vector<float> uncertainties;
    for (unsigned int iReco = 0; iReco < crossSectionsNom.size(); ++iReco)
    {
        if (this->IsUnderOverFlowBin(iReco))
        {
            uncertainties.push_back(-std::numeric_limits<float>::max());
            continue;
        }

        const auto standardDeviation = std::pow(sumSquaredDiffs.at(iReco) / static_cast<float>(m_nBootstrapUniverses), 0.5f); 
        uncertainties.push_back(standardDeviation);
    }

    return uncertainties;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

std::vector< std::vector<float> > CrossSectionHelper::XSec::GetSmearingMatrixMCStatUncertainty()
{
    if (m_nBootstrapUniverses == 0)
        throw std::logic_error("XSec::GetCrossSectionMCStatUncertaintyPerRecoBin - zero bootstrap universes have been generated");

    // Get the nominal smearing matrix
    const auto smearingMatrixNom = this->GetSmearingMatrix();
    const auto nBins = smearingMatrixNom.size();

    // Now get the cross-section in each universe of statistical variations
    const std::vector<std::string> systematicParams = {"stat"};

    // Setup an empty matrix
    std::vector< std::vector<float> > sumSquaredDiffs;
    for (unsigned int iBin = 0; iBin < nBins; ++iBin)
        sumSquaredDiffs.emplace_back(nBins, 0.f);

    for (unsigned int iUni = 0; iUni < m_nBootstrapUniverses; ++iUni)
    {
        const auto smearingMatrixUni = this->GetSmearingMatrix(systematicParams, iUni);

        for (unsigned int iReco = 0; iReco < nBins; ++iReco)
        {
            for (unsigned int iTrue = 0; iTrue < nBins; ++iTrue)
            {
                sumSquaredDiffs.at(iReco).at(iTrue) += std::pow(smearingMatrixUni.at(iReco).at(iTrue) - smearingMatrixNom.at(iReco).at(iTrue), 2);
            }
        }
    }
    
    std::vector< std::vector<float> > uncertainties;
    for (unsigned int iReco = 0; iReco < nBins; ++iReco)
    {
        std::vector<float> trueBins;
        for (unsigned int iTrue = 0; iTrue < nBins; ++iTrue)
        {
            const auto standardDeviation = std::pow(sumSquaredDiffs.at(iReco).at(iTrue) / static_cast<float>(m_nBootstrapUniverses), 0.5f); 
            trueBins.push_back(standardDeviation);
        }

        uncertainties.push_back(trueBins);
    }

    return uncertainties;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void CrossSectionHelper::XSec::PrintBinContents(const std::string &outputFileNamePrefix)
{
    if (m_binEdges.size() < 2)
        throw std::logic_error("XSec::PrintBinContents - binning hasn't been set");

    const auto nBins = m_binEdges.size() - 1;
    
    // Setup the table of the inputs to the cross-section
    FormattingHelper::Table inputsTable({
        "Bin", "",
        "Min", "Max", "Width", "",
        "Sig_t", "Eff_t", "",
        "SelMC_r", "SelData_r", "Bkg_r", "Eff_r"
    });
    const auto binWidths = this->GetBinWidths();
    const auto sig_t = this->CountSignalEventsPerTrueBin();
    const auto eff_t = this->GetEfficiencyPerTrueBin();
    const auto selMC_r = this->CountSelectedMCEventsPerRecoBin();
    const auto selData_r = this->CountSelectedBNBDataEventsPerRecoBin();
    const auto bkg_r = this->CountSelectedBackgroundEventsPerRecoBin();
    const auto eff_r = this->GetSmearedEfficiencyPerRecoBin();


    // Setup the table of cross-sections
    FormattingHelper::Table xSecTable({
        "Bin", "",
        "MCTrue_t", "MCSmear_r", "MCFake_r", "Data_r"
    });
    const auto xsecMCTrue_t = this->GetMCTrueCrossSectionPerTrueBin();
    const auto xsecMCSmear_r = this->GetSmearedMCTrueCrossSectionPerRecoBin();
    const auto xsecMCFake_r = this->GetCrossSectionPerRecoBin(false);
    const auto xsecData_r = this->GetCrossSectionPerRecoBin(true);

    // Setup the table for the smearing matrix
    std::vector<std::string> headers(1u, "Bin");
    for (unsigned int iReco = 0u; iReco < nBins; ++iReco)
        headers.push_back("r" + this->GetBinName(iReco));

    FormattingHelper::Table smearingTable(headers);
    const auto smearingMatrix = this->GetSmearingMatrix();


    // Populate the tables
    for (unsigned int iBin = 0u; iBin < nBins; ++iBin)
    {
        const auto binName = this->GetBinName(iBin);

        inputsTable.AddEmptyRow();
        inputsTable.SetEntry("Bin", binName);
        inputsTable.SetEntry("Min", m_binEdges.at(iBin));
        inputsTable.SetEntry("Max", m_binEdges.at(iBin + 1));
        inputsTable.SetEntry("Width", binWidths.at(iBin));
        inputsTable.SetEntry("Sig_t", sig_t.at(iBin));
        inputsTable.SetEntry("Eff_t", eff_t.at(iBin));
        inputsTable.SetEntry("SelMC_r", selMC_r.at(iBin));
        inputsTable.SetEntry("SelData_r", selData_r.at(iBin));
        inputsTable.SetEntry("Bkg_r", bkg_r.at(iBin));
        inputsTable.SetEntry("Eff_r", eff_r.at(iBin));
    
        xSecTable.AddEmptyRow();
        xSecTable.SetEntry("Bin", binName);
        xSecTable.SetEntry("MCTrue_t", xsecMCTrue_t.at(iBin));
        xSecTable.SetEntry("MCSmear_r", xsecMCSmear_r.at(iBin));
        xSecTable.SetEntry("MCFake_r", xsecMCFake_r.at(iBin));
        xSecTable.SetEntry("Data_r", xsecData_r.at(iBin));

        smearingTable.AddEmptyRow();
        smearingTable.SetEntry("Bin", "t" + binName);
        for (unsigned int iReco = 0u; iReco < nBins; ++iReco)
        {
            smearingTable.SetEntry(headers.at(iReco + 1), smearingMatrix.at(iReco).at(iBin));
        }
    }
        
    inputsTable.WriteToFile(outputFileNamePrefix + "_inputs.md");
    xSecTable.WriteToFile(outputFileNamePrefix + "_xSecs.md");
    smearingTable.WriteToFile(outputFileNamePrefix + "_smearing.md");
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void CrossSectionHelper::XSec::MakePlots(const std::string &outputFileNamePrefix, const bool includeSystematicUncertainties)
{
    // Get the data cross-section and it's uncertainties
    const auto xsecData = this->GetCrossSectionPerRecoBin(true);
    const auto xsecDataErr_stat = this->GetCrossSectionStatUncertaintyPerRecoBin(true);
    std::vector<float> xsecDataErr_syst(xsecData.size(), 0u);
    if (includeSystematicUncertainties)
    {
        xsecDataErr_syst = this->GetCrossSectionMCStatUncertaintyPerRecoBin(true);
    }
    const auto xsecDataErr = this->AddInQuadrature(xsecDataErr_stat, xsecDataErr_syst);

    // Get the smeared truth cross-section and it's uncertainties
    const auto xsecMCSmear = this->GetSmearedMCTrueCrossSectionPerRecoBin();
    std::vector<float> xsecMCSmearErr(xsecMCSmear.size(), 0u);
    if (includeSystematicUncertainties)
    {
        xsecMCSmearErr = this->GetSmearedMCTrueCrossSectionMCStatUncertaintyPerRecoBin();
    }

    // Make the plots
    auto hData = this->GetHistogram(xsecData, xsecDataErr);
    auto hMCSmear = this->GetHistogram(xsecMCSmear, xsecMCSmearErr);

    std::vector< shared_ptr<TH1F> > hists({hData, hMCSmear});
    this->SetHistogramYRanges(hists, true);

    PlottingHelper::SetLineStyle(hData, PlottingHelper::Primary);
    PlottingHelper::SetLineStyle(hMCSmear, PlottingHelper::Secondary);

    // Draw the histograms
    auto pCanvas = PlottingHelper::GetCanvas();

    // Draw the error band for MC
    auto hMCSmearClone = static_cast<TH1F *>(hMCSmear->Clone());
    const auto col = hMCSmearClone->GetLineColor();
    hMCSmearClone->SetFillStyle(1001);
    hMCSmearClone->SetLineColorAlpha(col, 0.f);
    hMCSmearClone->SetFillColorAlpha(col, 0.3f);

    hMCSmearClone->Draw("e2");
    hMCSmear->Draw("hist same");
    hData->Draw("e1 same");

    PlottingHelper::SaveCanvas(pCanvas, outputFileNamePrefix + "_dataMC");
    
    // Plot the smearing matrix
    const auto smearingMatrix = this->GetSmearingMatrix();
    auto hSmearingMatrix = this->GetHistogram(smearingMatrix);
    gStyle->SetPaintTextFormat("2.2f");
    hSmearingMatrix->SetMarkerSize(2.2); // Text size
    hSmearingMatrix->Draw("colz text");
    PlottingHelper::SaveCanvas(pCanvas, outputFileNamePrefix + "_smearingMatrix");
}

// -----------------------------------------------------------------------------------------------------------------------------------------

std::vector<float> CrossSectionHelper::XSec::AddInQuadrature(const std::vector<float> &a, const std::vector<float> &b) const
{
    if (a.size() != b.size())
        throw std::invalid_argument("XSec::AddInQuadrature - inputs are different size");

    std::vector<float> c;
    for (unsigned int i = 0; i < a.size(); ++i)
        c.push_back(std::pow(a.at(i)*a.at(i) + b.at(i)*b.at(i), 0.5f));

    return c;
}

// -----------------------------------------------------------------------------------------------------------------------------------------
                
std::shared_ptr<TH1F> CrossSectionHelper::XSec::GetHistogram(const std::vector<float> &binValues, const std::vector<float> &binUncertainties, const bool includeUnderflow, const bool includeOverflow) const
{
    if (m_binEdges.size() < 2)
        throw std::logic_error("XSec::GetHistogram - binning hasn't been set");

    const auto nBins = m_binEdges.size() - 1;

    if (binValues.size() != nBins)
        throw std::invalid_argument("XSec::GetHistogram - wrong number of input values for the bins already set");
    
    if (binUncertainties.size() != nBins)
        throw std::invalid_argument("XSec::GetHistogram - wrong number of input uncertainties for the bins already set");

    // Get the bin edges that we actually want to use
    std::vector<float> binEdges;
    const auto firstBinEdge = std::next(m_binEdges.begin(), (m_hasUnderflow && !includeUnderflow) ? 1u : 0u);
    const auto lastBinEdge = std::prev(m_binEdges.end(), (m_hasOverflow && !includeOverflow) ? 1u : 0u);
    binEdges.insert(binEdges.end(), firstBinEdge, lastBinEdge);
    
    // Setup the histogram
    const auto nBinsInUse = binEdges.size() - 1;
    std::shared_ptr<TH1F> pHist(new TH1F(("hXSec_" + std::to_string(m_histIndex++)).c_str(), "", nBinsInUse, binEdges.data()));
    
    // Fill the histogram with the desired values
    unsigned int histBinIndex = 1u;
    for (unsigned int iBin = 0; iBin < nBins; ++iBin)
    {
        if (this->IsUnderflowBin(iBin) && !includeUnderflow)
            continue;
        
        if (this->IsOverflowBin(iBin) && !includeOverflow)
            continue;

        pHist->SetBinContent(histBinIndex, binValues.at(iBin));
        pHist->SetBinError(histBinIndex, binUncertainties.at(iBin));

        histBinIndex++;
    }

    return pHist;
}

// -----------------------------------------------------------------------------------------------------------------------------------------
                
std::shared_ptr<TH2F> CrossSectionHelper::XSec::GetHistogram(const std::vector< std::vector<float> > &binValues) const
{
    if (m_binEdges.size() < 2)
        throw std::logic_error("XSec::GetHistogram - binning hasn't been set");

    const auto nBins = m_binEdges.size() - 1;

    if (binValues.size() != nBins)
        throw std::invalid_argument("XSec::GetHistogram - wrong number of input values for the bins already set");
    
    // Setup the histogram
    std::shared_ptr<TH2F> pHist(new TH2F(("hXSec_" + std::to_string(m_histIndex++)).c_str(), "", nBins, 0, nBins, nBins, 0, nBins));

    // Fill the histogram with the desired values
    for (unsigned int iBinX = 0; iBinX < nBins; ++iBinX)
    {
        const auto name = CrossSectionHelper::XSec::GetBinName(iBinX).c_str();
        pHist->GetXaxis()->SetBinLabel(iBinX+1, name);
        pHist->GetYaxis()->SetBinLabel(iBinX+1, name);

        for (unsigned int iBinY = 0; iBinY < nBins; ++iBinY)
        {
            // Put truth on x-axis & reco on y-axis
            pHist->SetBinContent(iBinY+1, iBinX+1, binValues.at(iBinX).at(iBinY));
        }
    }

    return pHist;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void CrossSectionHelper::XSec::SetHistogramYRanges(std::vector< std::shared_ptr<TH1F> > &hists, const bool capAtZero) const
{
    float min = std::numeric_limits<float>::max();
    float max = -std::numeric_limits<float>::max();

    for (const auto &pHist : hists)
    {
        std::cout << "Hist: " << pHist << std::endl;

        for (unsigned int iBin = 1; iBin <= static_cast<unsigned int>(pHist->GetNbinsX()); ++iBin)
        {    
            const auto content = static_cast<float>(pHist->GetBinContent(iBin));
            const auto error = static_cast<float>(pHist->GetBinError(iBin));

            std::cout << iBin << " : " << content << " +- " << error << std::endl;
            min = std::min(min, content - error);
            max = std::max(max, content + error);
        }
    }

    // Add some padding
    const auto padding = (max - min) * 0.05;
    min -= padding;
    max += padding;

    // Cap the lower end if required
    if (capAtZero)
        min = std::max(0.f, min);
    
    // Set the ranges
    for (auto &pHist : hists)
        pHist->GetYaxis()->SetRangeUser(min, max);
}

} // namespace ubcc1pi

#endif
