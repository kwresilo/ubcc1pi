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
                
                //std::vector<float> GetSelectedBNBDataEventsPerRecoBinUncertainty(const std::vector<std::string> &systematicParams, const unsigned int nUniverses) const;

                /**
                 *  @brief  Get the confusion matrix, M[reco][true]
                 *
                 *  @param  matrix the output confusion matrix
                 *  @param  uncertainties the output uncertainties on the confusion matrix elements
                 */
//                void GetConfusionMatrix(std::vector< std::vector<float> > &matrix, std::vector< std::vector<float> > &uncertainties) const;
                
                /**
                 *  @brief  Get the smearing matrix, S[reco][true] = P(reco | true)
                 *
                 *  @param  matrix the output smearing matrix
                 *  @param  uncertainties the output uncertainties on the smearing matrix elements
                 */
//                void GetSmearingMatrix(std::vector< std::vector<float> > &matrix, std::vector< std::vector<float> > &uncertainties) const;
                
                /**
                 *  @brief  Get the selection efficiencies in the true bins
                 *
                 *  @param  efficiencies the output efficiency per bin
                 *  @param  uncertainties the output uncertainty on the efficiency per bin
                 */
//                void GetTrueEfficiencies(std::vector<float> &efficiencies, std::vector<float> &uncertainties) const;

                /**
                 *  @brief  Get the smeared selection efficiencies
                 *
                 *  @param  efficiencies the output smeared efficiency per bin
                 *  @param  uncertainties the output uncertainty on the efficiency per bin
                 */
//                void GetSmearedEfficiencies(std::vector<float> &efficiencies, std::vector<float> &uncertainties) const;

                /**
                 *  @brief  Get the reco cross-section using either real BNB data or Overlays + EXT + Dirt
                 *
                 *  @param  useRealData if true, use real BNB data, if false use Overlays + EXT + Dirt
                 *  @param  bins the output cross-section values
                 *  @param  uncertainties the uncertainty on the bins
                 */
//                void GetRecoCrossSection(const bool useRealData, std::vector<float> &bins, std::vector<float> &uncertainties) const;
                
                /**
                 *  @brief  Make the plots
                 *
                 *  @param  fileNamePrefix the file name prefix
                 */
//                void MakePlots(const std::string &fileNamePrefix) const;

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
                 *  @brief  Get the differential cross-section in a bin given some parameters
                 *
                 *  @param  nTotal the total number of events selected
                 *  @param  nBackgrounds the expected number of background events
                 *  @param  efficiency the efficiency of the selection
                 *  @param  binWidth the width of the bin
                 *
                 *  @return the cross-section value
                 */
//                float GetCrossSection(const float nTotal, const float nBackgrounds, const float efficiency, const float binWidth) const;

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
//                float GetCrossSectionUncertainty(const float nTotal, const float totalErr, const float nBackgrounds, const float backgroundsErr, const float efficiency, const float efficiencyErr, const float binWidth) const;
                

                /**
                 *  @brief  Get the total number of signal events in each true bin
                 *
                 *  @param  mustBeSelected if we should insist that the events are selected
                 *
                 *  @return the true signal per bin
                 */
//                std::vector<float> GetTrueSignalPerBin(const bool mustBeSelected) const;
                
                /**
                 *  @brief  Get the overlay normalisation 
                 *
                 *  @return the normalisation factor for overlays
                 */
//                float GetOverlayNormalisation() const;
                
                /**
                 *  @brief  Get the selected BNB data in each bin
                 *
                 *  @param  bins the values per bin
                 *  @param  uncertainties the uncertainties
                 *  @param  useRealData if true we will use real BNB data, if false we use Overlays + EXT + Dirt as fake data
                 */
//                void GetBNBDataPerBin(std::vector<float> &bins, std::vector<float> &uncertainties, const bool useRealData = true) const;

                /**
                 *  @brief  Get the expected number of backgrounds in each bin
                 *
                 *  @param  bins the values per bin
                 *  @param  uncertainties the uncertainties
                 */
//                void GetBackgroundsPerBin(std::vector<float> &bins, std::vector<float> &uncertainties) const;
                
                /**
                 *  @brief  Get the purity in each bin
                 *
                 *  @param  bins the output purity values per bin
                 */
//                void GetPurityPerBin(std::vector<float> &bins) const;

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

                /**
                 *  @brief  Smear the input bins
                 *
                 *  @param  inputBins the input bin values
                 *  @param  inputUncertainties the uncertainties on the bin values
                 *  @param  bins the output smeared bins
                 *  @param  uncertainties the output smeared uncertainties
                 */
//                void SmearBins(const std::vector<float> &inputBins, const std::vector<float> &inputUncertainties, std::vector<float> &bins, std::vector<float> &uncertainties) const;

                /**
                 *  @brief  Print a matrix as a table
                 *
                 *  @param  matrix the matrix to print
                 *  @param  uncertainties the uncertainties on the values in the matrix
                 *  @param  outputFileName the output file name for the table
                 */
//                void PrintMatrix(const std::vector< std::vector<float> > &matrix, const std::vector< std::vector<float> > &uncertainties, const std::string &outputFileName) const;

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
        for (unsigned int iBinY = 0; iBinY < nBins; ++iBinY)
        {
            pHist->SetBinContent(iBinX, iBinY, binValues.at(iBinX).at(iBinY));
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

// -----------------------------------------------------------------------------------------------------------------------------------------

/*
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

void CrossSectionHelper::XSec::GetPurityPerBin(std::vector<float> &bins) const
{
    if (!bins.empty())
        throw std::invalid_argument("XSec::GetPurityPerBin - input bins vector isn't empty");
    
    const bool useRealData = false;
    std::vector<float> selected, selectedErr;
    this->GetBNBDataPerBin(selected, selectedErr, useRealData);
    
    std::vector<float> backgrounds, backgroundsErr;
    this->GetBackgroundsPerBin(backgrounds, backgroundsErr);

    const auto nBins = selected.size();
    if (backgrounds.size() != nBins)
        throw std::logic_error("XSec::GetPurityPerBin - Number of selected and backgrounds bins don't match!");

    for (unsigned int i = 0; i < nBins; ++i)
    {
        const auto selectedInBin = selected.at(i);
        if (selectedInBin <= std::numeric_limits<float>::epsilon())
            throw std::logic_error("XSec::GetPurityPerBin - Found empty bin - no selected event in MC");

        const auto signalInBin = selectedInBin - backgrounds.at(i);
        const auto purityInBin = signalInBin / selectedInBin;

        bins.push_back(purityInBin);
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

    std::vector<float> purities;
    this->GetPurityPerBin(purities);

    // Print the bins
    FormattingHelper::PrintLine();
    std::cout << "Bins" << std::endl;
    FormattingHelper::PrintLine();
    FormattingHelper::Table binTable({"Bin", "", "Min", "Max", "", "True MC signal", "", "BNB data", "Fake data (MC)", "data - MC", "", "Backgrounds MC", "Purity MC"});
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
        binTable.SetEntry("Purity MC", purities.at(i));
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

// TODO this function is vomit on toast - make it less so.
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

    gStyle->SetPaintTextFormat("2.2f");
    hSmearing->SetMarkerSize(2.2); // Text size
    hSmearing->Draw("colz text");
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
//    auto hXSecMCTrueClone = static_cast<TH1F *>(hXSecMCTrue->Clone());
//    hXSecMCTrue->SetLineColor(kAzure - 2);
//    hXSecMCTrueClone->SetFillStyle(1001);
//    hXSecMCTrueClone->SetLineColorAlpha(kAzure - 2, 0.f);
//    hXSecMCTrueClone->SetFillColorAlpha(kAzure - 2, 0.3f);
//    hXSecMCTrueClone->GetYaxis()->SetRangeUser(0.f, maxValue * 1.05f);
    
//    hXSecMCTrueClone->Draw("e2");
//    hXSecMCTrue->Draw("hist same");

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
*/
} // namespace ubcc1pi

#endif
