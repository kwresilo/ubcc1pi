/**
 *  @file  ubcc1pi_standalone/Helpers/CrossSectionHelper.h
 *
 *  @brief The header file for the cross section helper class
 */

#ifndef UBCC1PI_STANDALONE_HELPERS_CROSS_SECTION_HELPER
#define UBCC1PI_STANDALONE_HELPERS_CROSS_SECTION_HELPER

#include "ubcc1pi_standalone/Interface/Event.h"
#include "ubcc1pi_standalone/Objects/Config.h"

#include <string>
#include <unordered_map>
#include <map>
#include <vector>
#include <memory>
#include <random>
#include <TH1F.h>
#include <TH2F.h>
#include <TGraph.h>

namespace ubcc1pi
{

/**
 *  @brief  The cross section helper class
 */
class CrossSectionHelper
{
    public:
        
        /**
         *  @brief  A mapping from a systematic parameter name to the associated weights (one per universe)
         */
        typedef std::unordered_map<std::string, std::vector<float> > SystematicWeightsMap;

        /**
         *  @brief  A mapping from a systematic parameter name to the number of universes
         */
        typedef std::map<std::string, unsigned int> SystematicUniverseSizesMap;

        /**
         *  @brief  A pair of values, first = name of a systematic parameter, second = the number of universes
         */
        typedef std::pair<std::string, unsigned int> SystematicParamUniversesPair;

        /**
         *  @brief  A vector of SystematicParamUniversesPair
         */
        typedef std::vector<SystematicParamUniversesPair> SystematicParamUniversesPairVector;

        /**
         *  @brief  A set of mutually exclusive parameters, first = name for the combined parameters, second = number of universes, third = list of mutually exclusive params
         */
        typedef std::tuple<std::string, unsigned int, std::vector<std::string> > MutuallyExclusiveParam;

        /**
         *  @brief  A vector of MutuallyExclusiveParam
         */
        typedef std::vector<MutuallyExclusiveParam> MutuallyExclusiveParamVector;

        /**
         *  @brief  A pair of matricies - the first element is a covariance matrix, and the second if the bias vector
         */
        typedef std::pair<std::shared_ptr<TH2F>, std::shared_ptr<TH1F> > CovarianceBiasPair;

        /**
         *  @brief  The cross-section class
         */
        class CrossSection
        {
            private:
                /**
                 *  @brief  A map with indices [systematicParamterName][universeIndex] to an arbitrary type, T
                 *
                 *  @tparam T the mapped value type
                 */
                template <typename T>
                using SystMap = std::unordered_map< std::string, std::vector< T > >;
                
                /**
                 *  @brief  A map with indices [systematicParamterName][universeIndex] to an arbitrary type, T - used to cache values
                 *
                 *  @tparam T the mapped value type
                 */
                template <typename T>
                using SystCache = std::unordered_map<std::string, std::unordered_map<unsigned int, T> >;
            
            public:

                /**
                 *  @brief  The input information required by the cross-section class that is common to all cross-sections
                 */
                struct InputData
                {
                    std::shared_ptr<TH1F>        m_flux;                 ///< The flux distribution in bins of true neutrino energy
                    float                        m_exposurePOT;          ///< The total POT for the BNB data sample (units of e+20 POT)
                    float                        m_nTargets;             ///< The number of target nucleons (units e+31 nucleons)

                    SystematicUniverseSizesMap   m_systUniverseSizesMap; ///< The mapping from a systematic parameter name to the number of universes for that parameter
                    std::vector<std::string>     m_fluxParameters;       ///< The systematic parameters for which we should consider variations on the flux
                };

                /**
                 *  @brief  Constructor
                 *
                 *  @param  binEdges the bin edges of the cross-section
                 *  @param  hasUnderflow if the first bin is an underflow bin
                 *  @param  hasOverflow if the last bin is an overflow bin
                 *  @param  scaleByBinWidth if we should should divide by the bin-width when calculating the cross-section
                 *  @param  inputData the flux, exposure, number of targets, and systematic parameters to use for the cross-section calculation
                 */
                CrossSection(const std::vector<float> &binEdges, const bool hasUnderflow, const bool hasOverflow, const bool scaleByBinWidth, const InputData &inputData);

                /**
                 *  @brief  Check if this cross-section has an underflow bin
                 *
                 *  @return has underflow
                 */
                bool HasUnderflowBin() const;
                
                /**
                 *  @brief  Check if this cross-section has an overflow bin
                 *
                 *  @return has overflow
                 */
                bool HasOverflowBin() const;

                /**
                 *  @brief  Add a neutrino event - the event added isn't used explicitly in the cross-section calculation, but is needed to
                 *          be able to get the total flux in each universe.
                 *
                 *  @param  trueNuEnergy the true neutrino energy
                 *  @param  nominalWeight the nominal event weight
                 *  @param  systWeightsMap the systematic weights map
                 */
                void AddNeutrinoEvent(const float trueNuEnergy, const float nominalWeight, const SystematicWeightsMap &systWeightsMap);

                /**
                 *  @brief  Add a signal event
                 *
                 *  @param  trueValue the true value of the parameter
                 *  @param  recoValue the reconstructed value of the parameter
                 *  @param  isSelected if the event passes the selection
                 *  @param  nominalWeight the nominal event weight
                 *  @param  systWeightsMap the systematic weights map
                 */
                void AddSignalEvent(const float trueValue, const float recoValue, const bool isSelected, const float nominalWeight, const SystematicWeightsMap &systWeightsMap);
                
                /**
                 *  @brief  Add a selected background event with systematic weights
                 *
                 *  @param  recoValue the reconstructed value of the parameter
                 *  @param  nominalWeight the nominal event weight
                 *  @param  systWeightsMap the systematics weight map
                 */
                void AddSelectedBackgroundEvent(const float recoValue, const float nominalWeight, const SystematicWeightsMap &systWeightsMap);
                
                /**
                 *  @brief  Add a selected BNB data event
                 *
                 *  @param  recoValue the reconstructed value of the parameter
                 */
                void AddSelectedBNBDataEvent(const float recoValue);
                
                /**
                 *  @brief  Get the smearing matrix in the nominal universe
                 *
                 *  @return the smearing matrix (x-axis is true, y-axis is reco)
                 */
                std::shared_ptr<TH2F> GetSmearingMatrix() const;
                
                /**
                 *  @brief  Get the forward folded cross section in the nominal universe
                 *
                 *  @return the cross section in reco bins
                 */
                std::shared_ptr<TH1F> GetCrossSection() const;

                /**
                 *  @brief  Get the covariance matricies for the cross-section for each systematic paramter
                 *
                 *  @return a map from systematic paramter name to the covarianve matrix matrix / bias vector for that parameter
                 */
                std::map< std::string, CovarianceBiasPair > GetCrossSectionCovarianceMatricies();
                
                /**
                 *  @brief  Get the scatter plots of the value of the cross-sections in each pair of bins for each universe of each paramter
                 *
                 *  @return the scatter plots
                 */
                std::map< std::string, std::vector< std::tuple<unsigned int, unsigned int, std::shared_ptr<TGraph> > > > GetCrossSectionBinScatterPlots();
                
                /**
                 *  @brief  Get the flux variations for each flux systematic parameter
                 *
                 *  @return the mapping from a flux systematic parameter to the distribution of the variations of the flux
                 */
                std::map< std::string, std::shared_ptr<TH2F> > GetFluxVariations() const;
                
                /**
                 *  @brief  Check if a bin index is an underflow or an overflow bin
                 *
                 *  @param  binIndex the ROOT bin index (enumeratred from 1)
                 *
                 *  @return if the bin is under/overflow
                 */
                bool IsUnderOverflowBin(const unsigned int binIndex) const;

                /**
                 *  @brief  Clear the cached cross-sections to free up some memory
                 */
                void ClearCache();

            private:
                
                /**
                 *  @brief  Get a new empty 1D histogram with a unique name with the same binning as the flux
                 *
                 *  @return the histogram
                 */
                std::shared_ptr<TH1F> GetEmptyFluxHist() const;

                /**
                 *  @brief  Get a new empty 1D histogram with a unique name
                 *
                 *  @return the histogram
                 */
                std::shared_ptr<TH1F> GetEmptyHist1D() const;
                
                /**
                 *  @brief  Get a new empty 2D histogram with a unique name
                 *
                 *  @return the histogram
                 */
                std::shared_ptr<TH2F> GetEmptyHist2D() const;

                /**
                 *  @brief  Get the smearing matrix from the reco vs. true distribution
                 *
                 *  @param  signalSelectedRecoTrue the reco vs. true distribution of selected signal events
                 *
                 *  @return the smearing matrix
                 */
                std::shared_ptr<TH2F> GetSmearingMatrix(const std::shared_ptr<TH2F> &signalSelectedRecoTrue) const;

                /**
                 *  @brief  Smear the input values using the supplied smearing matrix
                 *
                 *  @param  values the input values to smear (in true bins)
                 *  @param  smearingMatrix the smearing matrix
                 *
                 *  @return the smeared values (in reco bins)
                 */
                std::shared_ptr<TH1F> Smear(const std::shared_ptr<TH1F> &values, const std::shared_ptr<TH2F> &smearingMatrix) const;

                /**
                 *  @brief  Get the smeared efficiency in reco bins
                 *
                 *  @param  signalSelectedRecoTrue the reco vs. true distribution of selected signal events
                 *  @param  signalAllTrue the true distribution of all signal events
                 *
                 *  @return the smeared efficiency per reco bin
                 */
                std::shared_ptr<TH1F> GetSmearedEfficiency(const std::shared_ptr<TH2F> &signalSelectedRecoTrue, const std::shared_ptr<TH1F> &signalAllTrue) const;

                /**
                 *  @brief  Get the forward-folded cross-section in reco bins
                 *
                 *  @param  selected the distribution of selected events (reco bins)
                 *  @param  background the distribution of background events (reco bins)
                 *  @param  smearedEff the smeared efficiency (reco bins)
                 *  @param  totalFlux the total flux (units e-10 / cm^2 / POT)
                 *
                 *  @return the cross-section
                 */
                std::shared_ptr<TH1F> GetCrossSection(const std::shared_ptr<TH1F> &selected, const std::shared_ptr<TH1F> &background, const std::shared_ptr<TH1F> &smearedEff, const float totalFlux) const;

                /**
                 *  @brief  Get the cached cross-section in a given systematic universe. If the cached value doesn't exist, then this
                 *          function will calculate the cross-section and cache it.
                 *
                 *  @param  systParameter the systematic paramter
                 *  @param  universeIndex the universe index
                 *
                 *  @return the cross-section
                 */
                std::shared_ptr<TH1F> GetCachedCrossSection(const std::string &systParameter, const unsigned int universeIndex);

                /**
                 *  @brief  Get the scatter plots showing the universe variations for each pair of bins
                 *
                 *  @param  systParameter the systematic paramter
                 *
                 *  @return a vector of tuples, first two elements are the indices of the bins and the last is the scatter plot
                 */
                std::vector< std::tuple<unsigned int, unsigned int, std::shared_ptr<TGraph> > > GetCrossSectionBinScatterPlots(const std::string &systParameter);
                
                /**
                 *  @brief  Get the cross-section covariance matrix for a given systematic paramter
                 *
                 *  @param  systParameter the systematic parameter to apply
                 *
                 *  @return a pair of matricies, first is the covariance matrix, second is the bias vector
                 */
                CovarianceBiasPair GetCrossSectionCovarianceMatrix(const std::string &systParameter);

                /**
                 *  @brief  Get a reweighted flux distribution according to the input neutrino energy spectrum (for the desired universe)
                 *          
                 *          In each true neutrino energy bin, this function finds the ratio of the total number of neutrino events in the
                 *          supplied universe to the nominal universe. The number of events in a given energy bin depends on the
                 *          cross-section and the flux at that energy. If the input universe differs only in the flux to the nominal
                 *          universe, then this ratio gives the fractional change in the flux. This is then used to reweight the nominal
                 *          flux distribution.
                 *
                 *  @param  nuEnergyUniverse the total neutrino energy spectrum in the desired universe
                 *
                 *  @return the reweighted flux distribution
                 */
                std::shared_ptr<TH1F> GetReweightedFlux(const std::shared_ptr<TH1F> &nuEnergyUniverse) const;
                
                /**
                 *  @brief  Get the distribution of reweighted fluxes for a given flux systematic parameter
                 *
                 *  @param  systParameter the systematic flux parameter
                 *
                 *  @return the flux universe variations 
                 */
                std::shared_ptr<TH2F> GetFluxVariations(const std::string &systParameter) const;

                // Binning information
                std::vector<float>                  m_binEdges;                   ///< The bin edges
                bool                                m_hasUnderflow;               ///< If the first bin is an underflow bin
                bool                                m_hasOverflow;                ///< If the last bin is an overflow bin
                bool                                m_scaleByBinWidth;            ///< If we should scale by bin width
                InputData                           m_inputData;                  ///< The extra information required to calculate the cross-section

                // The raw event counts in each systematic universe
                SystMap<std::shared_ptr<TH1F> >     m_allEventsNuEnergyTrue;      ///< The true neutrino energy spetrum of all events for each universe
                SystMap<std::shared_ptr<TH2F> >     m_signalSelectedRecoTrue;     ///< The 2D histograms of selected signal events with reco & true bin indices for each universe
                SystMap<std::shared_ptr<TH1F> >     m_signalAllTrue;              ///< The 1D histograms of all signal events with true bin indices for each universe
                SystMap<std::shared_ptr<TH1F> >     m_backgroundSelectedReco;     ///< The 1D histograms of selected background events with reco bin indices for each universe

                // The raw event counts in the nominal universe
                std::shared_ptr<TH1F>               m_allEventsNuEnergyNomTrue;   ///< The true neutrino energy spetrum of all events in the nominal universe
                std::shared_ptr<TH1F>               m_dataSelectedReco;           ///< The 1D histogram of selected BNB data events with reco bin indices for each universe
                std::shared_ptr<TH2F>               m_signalSelectedNomRecoTrue;  ///< The 1D histogram of selected BNB data events with reco bin indices in the nominal universe
                std::shared_ptr<TH1F>               m_signalAllNomTrue;           ///< The 1D histogram of all signal events with true bin indices in the nominal universe
                std::shared_ptr<TH1F>               m_backgroundSelectedNomReco;  ///< The 1D histogram of selected background events with reco bin indices in the nominal universe

                // The cached objects for speed
                SystCache<std::shared_ptr<TH1F> >   m_crossSectionCache;          ///< The cached cross-sections in each systematic universe to hold in memory
                bool                                m_shouldResetCache;           ///< If we need to re-calculate the cached values because something has changed
                
                static unsigned int                 m_histCount;                  ///< A counter for the total number of histograms - used to avoid name collisions
        };

        // ---------------------------------------------------------------------------------------------------------------------------------

        /**
         *  @brief  Get the event weights corresponding to a given systematic paramter (one per universe)
         *
         *  @param  truth the input event truth information
         *  @param  parameter the systematic parameter name
         *
         *  @return the event weights
         */
        static std::vector<float> GetSystematicWeights(const Event::Truth &truth, const std::string &parameter);
        
        /**
         *  @brief  Add the event weights for a given set of systematic parameters to the input systematic weights map
         *
         *  @param  truth the input truth information (this contains the weights)
         *  @param  params the parameters to add to the map - pairs of universe names and the number of universes expected
         *  @param  systWeightsMap the systematic weights map to update
         */
        static void AddSystematicWeights(const Event::Truth &truth, const SystematicParamUniversesPairVector &params, SystematicWeightsMap &systWeightsMap);

        /**
         *  @brief  add a set of mutually exclusive systematic parameters to the input systematics weights map and combine them as a single parameter
         *
         *  @param  truth the input truth information (this contains the weights)
         *  @param  params the input mutually exclusive parameters to combine 
         *  @param  systweightsmap the systematic weights map to update
         */
        static void AddMutuallyExclusiveSystematicWeights(const Event::Truth &truth, const MutuallyExclusiveParamVector &params, SystematicWeightsMap &systWeightsMap);

        /**
         *  @brief  Add the specified number of boostrap universes to the input systematic weights map
         *
         *  @param  nBootstrapUniverses the number of boostrap universes to use
         *  @param  systWeightsMap the systematic weights map to update
         */
        static void AddBootstrapWeights(const unsigned int nBootstrapUniverses, SystematicWeightsMap &systWeightsMap);

        /**
         *  @brief  Add weights of 1.f for each parameter in each universe supplied to the input systematic weights map
         *
         *  @param  params the input parameters 
         *  @param  systWeightsMap the systematic weights map to update
         */
        static void AddUnitWeights(const SystematicParamUniversesPairVector &params, SystematicWeightsMap &systWeightsMap);

        /**
         *  @brief  Add weights of 1.f for each combined parameter in each universe supplied to the input systematic weights map
         *
         *  @param  params the input paramters
         *  @param  systWeightsMap the systematic weights map to update
         */
        static void AddUnitWeights(const MutuallyExclusiveParamVector &params, SystematicWeightsMap &systWeightsMap);

        /**
         *  @brief  Get the mapping from systematic parameters name to the number of universe variations
         *
         *  @param systWeightsMap a mapping from systematic parameters to weights in each universe
         *
         *  @return the output systematic universe sizes map
         */
        static SystematicUniverseSizesMap GetSystematicUniverseSizesMap(const SystematicWeightsMap &systWeightsMap);

        /**
         *  @brief  Check that the input map of systematic weights has the supplied dimensions
         *
         *  @param  systWeightsMap the input mapping from systematic parameters to the weights in each universe
         *  @param  systUniverseSizesMap the input mappin from systematic parameters to the expected number of universes
         *
         *  @return boolean, true if the dimensions of the input objects match
         */
        static bool CheckSystematicWeightsMapDimensions(const SystematicWeightsMap &systWeightsMap, const SystematicUniverseSizesMap &systUniverseSizesMap);

        /**
         *  @brief  Extend the analysis bin edges to the supplied extremities
         *
         *  @param  min the minimum possible value
         *  @param  max the maximum possible value
         *  @param  binEdges the analysis bin edges
         *  @param  extendedBinEdges the output extended bin edges (from min -> max)
         *  @param  hasUnderflow if an underflow bin was added
         *  @param  hasOverflow if an overflow bin was added
         */
        static void GetExtendedBinEdges(const float min, const float max, const std::vector<float> &binEdges, std::vector<float> &extendedBinEdges, bool &hasUnderflow, bool &hasOverflow);

        /**
         *  @brief  Extend the analysis bin edges to the supplied extremities
         *
         *  @tparam T the type of the binning config
         *  @param  binningConfig the binning configuration, must have a min, max and binEdges property
         *
         *  @return a tuple containing the: extended bin edges, hasUnderflow, hasOverflow
         */
        template <typename T>
        static std::tuple<std::vector<float>, bool, bool> GetExtendedBinEdges(const T &binningConfig);
        
        /**
         *  @brief  Get the histogram corresponding to the input flux distribution
         *
         *  @param  flux the input flux distribution
         *
         *  @return the flux histogram
         */
        static std::shared_ptr<TH1F> GetFluxHistogram(const Config::Flux &flux);

        /**
         *  @brief  Get the total flux in the input flux histogram
         *
         *  @param  fluxHist the input flux histogram
         *
         *  @return the total flux
         */
        static float GetTotalFlux(const std::shared_ptr<TH1F> &fluxHist);

        /**
         *  @brief  Get the uncertainty for a given pair of bins in the input covariance as eigen vectors.
         *          This function returns the major & minor eigenvectors of the error sub-matrix for the supplied bins. The length of the
         *          vectors are the corresponding eigenvalue. In a scatter plot of values of the cross-section (in iBin vs. jBin) from the
         *          universes used to construct the covariance matrix... these vectors point along the primary & secondary axes
         *          of the 2D Gaussian distribution that we are using to model the uncertainty distribution. The length of the vectors is
         *          the standard deviation of the Gaussian along those axes.
         *
         *  @param  covarianceBias the input pair of corresponding covariance matric & bias vector
         *  @param  iBin a bin index in the inupt matrix
         *  @param  jBin another bin index in the input matrix
         *
         *  @return the uncertainty eigenvectors
         */
        static std::pair<TVector2, TVector2> GetUncertaintyEigenVectors(const CovarianceBiasPair &covarianceBias, const unsigned int iBin, const unsigned int jBin);
                
        /**
         *  @brief  Get the total covariance matrix by combining the input covarainces and biases for the systematic parameters used
         *
         *  @param  covarianceBiasPairs the input covariance matrix / bias vector pairs
         *
         *  @return the total covaraince matrix
         */
        static std::shared_ptr<TH2F> GetTotalCovarianceMatrix(const std::map< std::string, CovarianceBiasPair > &covarianceBiasPairs);
    
    private:
        static std::default_random_engine m_generator; ///< The random number generator
        static unsigned int               m_histCount; ///< A counter to ensure histogram names are unique
};

// -----------------------------------------------------------------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------------------------------------------------------------

unsigned int CrossSectionHelper::CrossSection::m_histCount = 0u;

// -----------------------------------------------------------------------------------------------------------------------------------------

unsigned int CrossSectionHelper::m_histCount = 0u;

// -----------------------------------------------------------------------------------------------------------------------------------------

std::default_random_engine CrossSectionHelper::m_generator = std::default_random_engine();
        
// -----------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
std::tuple<std::vector<float>, bool, bool> CrossSectionHelper::GetExtendedBinEdges(const T &binningConfig)
{
    std::vector<float> extendedBinEdges;
    bool hasUnderflow, hasOverflow;

    CrossSectionHelper::GetExtendedBinEdges(binningConfig.min, binningConfig.max, binningConfig.binEdges, extendedBinEdges, hasUnderflow, hasOverflow);

    return std::tuple<std::vector<float>, bool, bool>(extendedBinEdges, hasUnderflow, hasOverflow);
}

} // namespace ubcc1pi

#endif
