/**
 *  @file  ubcc1pi_standalone/Helpers/CrossSectionHelperNew.h
 *
 *  @brief The header file for the cross section helper class
 */

// TODO Don't forget to change the define statement below to remove NEW when required
#ifndef UBCC1PI_STANDALONE_HELPERS_CROSS_SECTION_NEW_HELPER
#define UBCC1PI_STANDALONE_HELPERS_CROSS_SECTION_NEW_HELPER

#include "ubcc1pi_standalone/Interface/Event.h"
#include "ubcc1pi_standalone/Objects/Config.h"

#include "ubsmear.h"

#include <functional>

namespace ubcc1pi
{

/**
 *  @brief  The cross section helper class
 */
class CrossSectionHelperNew
{
    public:
        /**
        *  @brief  A mapping from a systematic parameter name, to a vector of a quantity of arbitrary type T.
        *          Each element of the vector represents the value of a quantity in a different universe of the systematic paramter
        *
        *  @tparam T the mapped value type
        */
        template <typename T>
        using SystMap = std::unordered_map< std::string, std::vector<T> >;

        /**
        *  @brief  Like a SystMap but using an unordered_map instead of a vector as the mapped type.
        *          In this way it's possible for the map to have empty entries for a given universe index
        *
        *  @tparam T the mapped value type
        */
        template <typename T>
        using SystPartialMap = std::unordered_map< std::string, std::unordered_map<unsigned int, T> >;

        /**
        *  @brief  A mapping from a systematic parameter name, to a quantity of arbitrary type T.
        *
        *  @tparam T the mapped value type
        */
        template <typename T>
        using SystUnisimMap = std::unordered_map< std::string, T >;

        /**
        *  @brief  A mapping from a systematic parameter name to the number of universes
        */
        typedef std::unordered_map< std::string, unsigned int > SystDimensionsMap;

        /**
        *  @brief  A mapping from a systematic parameter name, to a vector of mutually exclusive systematic parameter names and the number
        *          universes. For example, if parameters {"a", "b", "c"} were all mutually exclusive (i.e. can not simultaneously vary),
        *          then this map might contain a new parameter "abc" as the key which maps to those individual parameter names. The new
        *          parameter is then treated like any other systematic parameter, with the number of universes specified. See the flux
        *          hadron production parameter as an example.
        */
        typedef std::unordered_map<std::string, std::pair< std::vector<std::string>, unsigned int> > SystMutuallyExclusiveDimensionsMap;

        /**
        *  @brief  A mapping from a unisim systematic parameter name, to the name of the central-value sample to which it should be compared
        */
        typedef std::unordered_map<std::string, std::string> SystUnisimDimensionsMap;

        /**
        *  @brief  A mapping from a systematic paramter name to a vector of floats, one per universe
        */
        typedef SystMap<float> SystFloatMap;

        /**
        *  @brief  A mapping from a systematic parameter name to a vector of 1D histograms, one per universe
        */
        typedef SystMap< std::shared_ptr<TH1F> > SystTH1FMap;

        /**
        *  @brief  A mapping from a systematic parameter name to a vector of 2D histograms, one per universe
        */
        typedef SystMap< std::shared_ptr<TH2F> > SystTH2FMap;

        /**
        *  @brief  A mapping from a systematic parameter name to a 1D histogram
        */
        typedef SystUnisimMap< std::shared_ptr<TH1F> > SystUnisimTH1FMap;

        /**
        *  @brief  A mapping from a systematic parameter name to a 2D histogram
        */
        typedef SystUnisimMap< std::shared_ptr<TH2F> > SystUnisimTH2FMap;

        /**
        *  @brief  A wrapper around an arbitrary function that caches the result of the function in a SystPartialMap.
        *          If the function has already been executed for a given systematic parameter / universe then the a call to the function
        *          will return the cached value instead. This is provides a helful way of storing the result of a function that we might
        *          want to use in multiple places. This is done for performance reasons, so please excuse the added complexity!
        *
        *  @tparam R the return type of the function
        *  @tparam Args the argument types of the function
        */
        template <typename R, typename... Args>
        class SystCacheFunction
        {
            public:
                /**
                *  @brief  The function type around which we are wrapping.
                *          The first parameter must be the systematic parameter name, the second must be the universe index, and the
                *          remaining paramters are arbitrary
                */
                using Function = std::function<R(const std::string &, const unsigned int, Args ... args)>;

                /**
                *  @brief  Constructor
                *
                *  @param  function the function definition
                *  @param  dimensions the dimensions of the SystPartialMap in which we should cache
                */
                SystCacheFunction(const Function &function, const SystDimensionsMap &dimensions);

                /**
                *  @brief  The () operator is used to call the function (or it's already been called retreive the cached value)
                *
                *  @param  paramName the systematic parameter name
                *  @param  universeIndex the universe index
                *  @param  args the arbitrary arguments to the function
                *
                *  @return the return value of the stored function
                */
                R operator()(const std::string &paramName, const unsigned int universeIndex, Args ... args);

                /**
                *  @brief  Clear the cached values
                */
                void ClearCache();

            private:

                Function           m_function;   ///< The function to execute
                SystDimensionsMap  m_dimensions; ///< The dimensions of the cache
                SystPartialMap<R>  m_cache;      ///< The cached return values
        };

        /**
        *  @brief  The flux reweightor class. Used to calculate the flux in a given systematic universe.
        *
        *          The flux is defined as a function of the true neutrino energy. Here we get the total event rate as a function of neutrino
        *          energy for the nominal universe and universes with a varied flux. The ratio of the event rate in each energy bin
        *          (varied/nominal) is used to scale the nominal flux to give the flux distribution in the varied universe.
        */
        class FluxReweightor
        {
            public:
                /**
                *  @brief  Constructor
                *
                *  @param  binEdges the bin edges in which the neutrino flux is defined
                *  @param  binValuesNominal the nominal neutrino flux in each energy bin
                *  @param  fluxWeightsDimensions the dimensions of the flux universes (mapping from flux paramter name to the number of universes)
                */
                FluxReweightor(const std::vector<float> &binEdges, const std::vector<float> &binValuesNominal, const SystDimensionsMap &fluxWeightsDimensions);

                /**
                *  @brief  Add an event to the total event rate spectrum in each systematic univerese
                *          Here the spectrum in the nominal universe is added to using the nominal weight, and the spectrum in each varied
                *          universe is added to using the product of the nominal & weight for that universe
                *
                *  @param  trueNuEnergy the true neutrino energy
                *  @param  nominalWeight the nominal event weight
                *  @param  fluxWeights the flux universe weights (mapping from flux parameter name to the weights in each universe)
                */
                void AddEvent(const float trueNuEnergy, const float nominalWeight, const SystFloatMap &fluxWeights);

                /**
                *  @brief  Get the flux distribution in the nominal universe
                *
                *  @return the nominal flux
                */
                std::shared_ptr<TH1F> GetNominalFlux() const;

                /**
                *  @brief  Get the integrated flux in the nominal universe
                *
                *  @return the integrated nominal flux
                */
                float GetIntegratedNominalFlux() const;

                /**
                *  @brief  Get the integrated flux in a given universe
                *          This function will try to retrieve the return value from a cache. If that's not possible, it will calculated the
                *          return value and store it in the cache for later
                *
                *  @param  paramName the flux parameter name
                *  @param  universeIndex the universe index
                *
                *  @return the integrated flux in the requested univese
                */
                float GetIntegratedFluxVariation(const std::string &paramName, const unsigned int universeIndex);

            private:

                /**
                *  @brief  Get the integrated flux
                *
                *  @param  pFlux the input flux (in energy bins)
                *
                *  @return the integrated flux
                */
                float GetIntegratedFlux(const std::shared_ptr<TH1F> &pFlux) const;

                SystDimensionsMap     m_dimensions;             ///< The systematic paramter dimensions
                std::shared_ptr<TH1F> m_pFluxNominal;           ///< The nominal flux distribution
                std::shared_ptr<TH1F> m_pSpectrumNominal;       ///< The event rate spectrum in the nominal universe
                SystTH1FMap           m_spectrumVariations;     ///< The event rate spectrum in each systematic universe

                SystCacheFunction<float> m_getIntegratedFluxVariationFunction;  ///< The cache function for getting the integrated flux in a given universe
        };

        /**
        *  @brief  The cross-section class. Describes a forward-folded differential cross-section.
        */
        class CrossSection
        {
            public:
                /**
                *  @brief  A structure storing the dimensions of the systematic parameters that are to be considered
                */
                struct SystParams
                {
                    unsigned int             nBootstrapUniverses; ///< The number of bootstrap universes to use for the MC stat uncertainty
                    SystDimensionsMap        fluxDimensions;      ///< The dimensions of the flux systematic parameters
                    SystDimensionsMap        xsecDimensions;      ///< The dimensions of the cross-section systematic parameters
                    SystUnisimDimensionsMap  detVarDimensions;    ///< The detector variation systematic parameters (and their central-value sample identifiers)
                };

                /**
                *  @brief  Constructor
                *
                *  @param  systParams the systematic parameters that will be applied
                *  @param  binEdges the bin edges for the kinematic quantity including any underflow and overflow bins
                *  @param  hasUnderflow if there is an underflow bin
                *  @param  hasOverflow if there is an overflow bin
                */
                CrossSection(const SystParams &systParams, const std::vector<float> &binEdges, const bool hasUnderflow, const bool hasOverflow);

                /**
                *  @brief  Add a simulated signal event with flux and cross-section universe weights
                *          The events added using this function are used to produce:
                *            - The prediction for the forward-folded cross-section in reco-space
                *            - The nominal smearing matrix - which maps truth-space to reco-space (including the selection efficiency)
                *            - The smearing matrix in each systematic universe
                *            - The covariance matrix for the smearing matrix variations due to flux & cross-section uncertainties
                *
                *  @param  trueValue the true value of the kinematic quantity
                *  @param  recoValue the reconstructed value of the kinematic quantity
                *  @param  isSelected if the event is selected
                *  @param  nominalWeight the nominal event weight
                *  @param  fluxWeights the flux systematic event weights
                *  @param  xsecWeights the cross-section systematic event weights
                */
                void AddSignalEvent(const float trueValue, const float recoValue, const bool isSelected, const float nominalWeight, const SystFloatMap &fluxWeights, const SystFloatMap &xsecWeights);

                /**
                *  @brief  Add a simulated signal event from a detector variation sample
                *          The events added using this function are used to produce:
                *            - The fractional difference in the smearing matrix with each detector variation applied (the corresponding central-value sample is used to get the fraction)
                *
                *  @param  trueValue the true value of the kinematic quantity
                *  @param  recoValue the reconstructed value of the kinematic quantity
                *  @param  isSelected if the event is selected
                *  @param  nominalWeight the nominal event weight
                *  @param  paramName the name of the detector variation parameter
                */
                void AddSignalEventDetVar(const float trueValue, const float recoValue, const bool isSelected, const float nominalWeight, const std::string &paramName);

                /**
                *  @brief  Add a simulated signal event from a detector variation central-value sample (i.e with no parameters varied)
                *
                *  @param  trueValue the true value of the kinematic quantity
                *  @param  recoValue the reconstructed value of the kinematic quantity
                *  @param  isSelected if the event is selected
                *  @param  nominalWeight the nominal event weight
                *  @param  name the name of the central value sample
                */
                void AddSignalEventDetVarCV(const float trueValue, const float recoValue, const bool isSelected, const float nominalWeight, const std::string &name);

                /**
                *  @brief  Add a simulated selected backgound event with flux and cross-section universe weights
                *          The events added using this function are used to produce:
                *            - The nominal background prediction (that get's subtracted away from the data)
                *            - The background prediction in each systematic universe
                *            - The covariance matrix for the background-subtracted data due to flux & cross-secion uncertainties
                *
                *  @param  recoValue the reconstructed value of the kinematic quantity
                *  @param  nominalWeight the nominal event weight
                *  @param  fluxWeights the flux systematic event weights
                *  @param  xsecWeights the cross-section systematic event weights
                */
                void AddSelectedBackgroundEvent(const float recoValue, const float nominalWeight, const SystFloatMap &fluxWeights, const SystFloatMap &xsecWeights);

                /**
                *  @brief  Add a simulated selected background event from a detector variation sample
                *          The events added using this function are used to produce:
                *            - The fractional difference in the background prediction with each detector variation applied (using the corresponding central-value sample)
                *
                *  @param  recoValue the reconstructed value of the kinematic quantity
                *  @param  nominalWeight the nominal event weight
                *  @param  paramName the name of the detector variation parameter
                */
                void AddSelectedBackgroundEventDetVar(const float recoValue, const float nominalWeight, const std::string &paramName);

                /**
                *  @brief  Add a simulated selected background event from a detector variation central-value sample
                *
                *  @param  recoValue the reconstructed value of the kinematic quantity
                *  @param  nominalWeight the nominal event weight
                *  @param  name the name of the central value sample
                */
                void AddSelectedBackgroundEventDetVarCV(const float recoValue, const float nominalWeight, const std::string &name);

                /**
                *  @brief  Add a selected event from real BNB data
                *          The events added using this function are used to produce:
                *            - The measured "forward-folded cross-section" (= background-substracted & flux normalised event rate)
                *
                *  @param  recoValue the reconstructed value of the kinematic quantity
                */
                void AddSelectedBNBDataEvent(const float recoValue);

                // TODO add getter functions for the desired quantities (e.g. smearing matrices etc.)
                // Consider using SystCacheFunction objects if something will be needed multiple times

            private:

                SystParams          m_systParams;                    ///< The systematic parameters
                std::vector<float>  m_binEdges;                      ///< The bin edges
                ubsmear::UBXSecMeta m_metadata;                      ///< The cross-section metadata

                // Signal event distributions in truth space (all signal events before any selection)
                std::shared_ptr<TH1F> m_pSignal_true_nom;    ///< The signal event distribution in truth space in the nominal universe
                SystTH1FMap           m_signal_true_fluxes;  ///< The signal event distribution in truth space for each flux universe
                SystTH1FMap           m_signal_true_xsecs;   ///< The signal event distribution in truth space for each cross-section universe
                SystUnisimTH1FMap     m_signal_true_detVars; ///< The signal event distribution in truth space for each detector variation sample
                SystUnisimTH1FMap     m_signal_true_detCVs;  ///< The signal event distribution in truth space for each detector variation central-value sample

                // Selected signal event distributions in reco-vs-truth space
                std::shared_ptr<TH2F> m_pSignal_selected_recoTrue_nom;    ///< The selected signal event distribution in reco-vs-truth space in the nominal universe
                SystTH2FMap           m_signal_selected_recoTrue_fluxes;  ///< The selected signal event distribution in reco-vs-truth space for each flux universe
                SystTH2FMap           m_signal_selected_recoTrue_xsecs;   ///< The selected signal event distribution in reco-vs-truth space for each cross-section universe
                SystUnisimTH2FMap     m_signal_selected_recoTrue_detVars; ///< The selected signal event distribution in reco-vs-truth space for each detector variation sample
                SystUnisimTH2FMap     m_signal_selected_recoTrue_detCVs;  ///< The selected signal event distribution in reco-vs-truth space for each detector variation central-value sample

                // Selected background event distributions in reco space
                std::shared_ptr<TH1F> m_pBackground_selected_reco_nom;    ///< The selected background event distribution in reco space in the nominal universe
                SystTH1FMap           m_background_selected_reco_fluxes;  ///< The selected background event distribution in reco space for each flux universe
                SystTH1FMap           m_background_selected_reco_xsecs;   ///< The selected background event distribution in reco space for each cross-section universe
                SystUnisimTH1FMap     m_background_selected_reco_detVars; ///< The selected background event distribution in reco space for each detector variation sample
                SystUnisimTH1FMap     m_background_selected_reco_detCVs;  ///< The selected background event distribution in reco space for each detector variation central-value sample

                std::shared_ptr<TH1F> m_pBNBData_selected_reco; ///< The selected BNB data event distribution in reco space
        };

        /**
        *  @brief  Validate an input SystMap by checking it has the same dimensions as the input dimensions map
        *          This function will raise an exception if the the input map is invalid
        *
        *  @tparam T the mapped type
        *  @param  systMap the input SystMap
        *  @param  dimensions the expect dimensions
        */
        template <typename T>
        static void ValidateSystMap(const SystMap<T> &systMap, const SystDimensionsMap &dimensions);

        /**
        *  @brief  Get a shared pointer to a new TH1F and give it a unique name
        *
        *  @param  binEdges the bin edges
        *
        *  @return the TH1F
        */
        static std::shared_ptr<TH1F> GetTH1F(const std::vector<float> &binEdges);

        /**
        *  @brief  Get a shared pointer to a new TH2F and give it a unique name
        *
        *  @param  binEdges the bin edges
        *
        *  @return the TH2F
        */
        static std::shared_ptr<TH2F> GetTH2F(const std::vector<float> &binEdges);

        /**
        *  @brief  Get a SystTH1FMap with the supplied dimensions containing new TH1F objects with the supplied binning
        *
        *  @param  binEdges the bin edges
        *  @param  dimensions the systematic dimensions map
        *
        *  @return the SystTH1FMap
        */
        static SystTH1FMap GetSystTH1FMap(const std::vector<float> &binEdges, const SystDimensionsMap &dimensions);

        /**
        *  @brief  Get a SystUnisimTH1FMap with the supplied dimensions containing new TH1F objects with the supplied binning
        *
        *  @param  binEdges the bin edges
        *  @param  dimensions the systematic dimensions map
        *
        *  @return the SystUnisimTH1FMap
        */
        static SystUnisimTH1FMap GetSystUnisimTH1FMap(const std::vector<float> &binEdges, const SystUnisimDimensionsMap &dimensions);

        /**
        *  @brief  Get a SystUnisimTH1FMap with the CV samples of the supplied dimensions containing new TH1F objects with the supplied binning
        *
        *  @param  binEdges the bin edges
        *  @param  dimensions the systematic dimensions map
        *
        *  @return the SystUnisimTH1FMap
        */
        static SystUnisimTH1FMap GetCVSystUnisimTH1FMap(const std::vector<float> &binEdges, const SystUnisimDimensionsMap &dimensions);

        /**
        *  @brief  Get a SystTH2FMap with the supplied dimensions containing new TH2F objects with the supplied binning
        *
        *  @param  binEdges the bin edges
        *  @param  dimensions the systematic dimensions map
        *
        *  @return the SystTH2FMap
        */
        static SystTH2FMap GetSystTH2FMap(const std::vector<float> &binEdges, const SystDimensionsMap &dimensions);

        /**
        *  @brief  Get a SystUnisimTH2FMap with the supplied dimensions containing new TH2F objects with the supplied binning
        *
        *  @param  binEdges the bin edges
        *  @param  dimensions the systematic dimensions map
        *
        *  @return the SystUnisimTH2FMap
        */
        static SystUnisimTH2FMap GetSystUnisimTH2FMap(const std::vector<float> &binEdges, const SystUnisimDimensionsMap &dimensions);

        /**
        *  @brief  Get a SystUnisimTH2FMap with the CV samples of the supplied dimensions containing new TH2F objects with the supplied binning
        *
        *  @param  binEdges the bin edges
        *  @param  dimensions the systematic dimensions map
        *
        *  @return the SystUnisimTH2FMap
        */
        static SystUnisimTH2FMap GetCVSystUnisimTH2FMap(const std::vector<float> &binEdges, const SystUnisimDimensionsMap &dimensions);

        /**
        *  @brief  Get the weights corresponding to the systematic parameters in the input dimensions object. The parameter names must
        *          correspond to those in the input truth information with matching number of universes.
        *
        *  @param  truth the input truth information
        *  @param  dimensions the dimensions of the weights map to obtain
        *
        *  @return the weights map
        */
        static SystFloatMap GetWeightsMap(const Event::Truth &truth, const SystDimensionsMap &dimensions);

        /**
        *  @brief  Get the weights corresponding to the input systematic paramters in the input dimensions object. The parameters must
        *          either correspon to those in the input truth information, or be one of the mutually exclusive parameters supplied. In
        *          every case the number of universes must match the input truth information.
        *
        *  @param  truth the input truth information
        *  @param  dimensions the dimensions of the weights map to obtain (can include musually exclusive parameters)
        *  @param  mutuallyExclusiveDimensions the mutualy exclusive parameter dimensions
        *
        *  @return the weights map
        */
        static SystFloatMap GetWeightsMap(const Event::Truth &truth, const SystDimensionsMap &dimensions, const SystMutuallyExclusiveDimensionsMap &mutuallyExclusiveDimensions);

        /**
        *  @brief  Get the combined weights from a set of mutually exclusive parameters (i.e. only one parameter at a time has a weight other than 1.f)
        *
        *  @param  truth the input truth information
        *  @param  parameters the mutually exclusive parameter names
        *  @param  nUniverses the number of universes (must be the same for all parameters)
        *
        *  @return the weights in each universe from the parameter that was not 1.f
        */
        static std::vector<float> GetMutuallyExclusiveWeights(const Event::Truth &truth, const std::vector<std::string> &parameters, const unsigned int nUniverses);

    private:

        static unsigned int m_histCount; ///< A counter to keep track of the number of histograms produced
};

// -----------------------------------------------------------------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------------------------------------------------------------

template <typename R, typename... Args>
inline CrossSectionHelperNew::SystCacheFunction<R, Args...>::SystCacheFunction(const Function &function, const SystDimensionsMap &dimensions) :
    m_function(function),
    m_dimensions(dimensions)
{
}

// -----------------------------------------------------------------------------------------------------------------------------------------

template <typename R, typename... Args>
inline R CrossSectionHelperNew::SystCacheFunction<R, Args...>::operator()(const std::string &paramName, const unsigned int universeIndex, Args ... args)
{
    // Check that we have a valid parameter name and universe index
    const auto iterDimensions = m_dimensions.find(paramName);
    if (iterDimensions == m_dimensions.end())
        throw std::invalid_argument("SystCacheFunction - Unknown paramName: " + paramName);

    const auto nUniverses = iterDimensions->second;
    if (universeIndex >= nUniverses)
        throw std::out_of_range("SystCacheFunction - Universe index is out of range");

    // Check if we have a cached entry for this parameter name
    const auto iterCacheName = m_cache.find(paramName);
    if (iterCacheName != m_cache.end())
    {
        // Check if we have a cached entry for this universe index
        const auto iterCacheUniverse = iterCacheName->second.find(universeIndex);
        if (iterCacheUniverse != iterCacheName->second.end())
        {
            // We have a cached entry!
            return iterCacheUniverse->second;
        }
    }

    // We don't have an entry in the cache, so make one
    const auto result = m_function(paramName, universeIndex, args...);
    m_cache[paramName][universeIndex] = result;

    // Return the result
    return result;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

template <typename R, typename... Args>
inline void CrossSectionHelperNew::SystCacheFunction<R, Args...>::ClearCache()
{
    m_cache.clear();
}

// -----------------------------------------------------------------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------------------------------------------------------------

unsigned int CrossSectionHelperNew::m_histCount = 0u;

// -----------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
inline void CrossSectionHelperNew::ValidateSystMap(const SystMap<T> &systMap, const SystDimensionsMap &dimensions)
{
    // Insist that the maps have the same number of entries
    if (systMap.size() != dimensions.size())
        throw std::invalid_argument("CrossSectionHelperNew::ValidateSystMap - Number of parameters in input map doesn't match the supplied dimensions");

    // Insist that each entry of the input map has the right number of universes
    for (const auto &[paramName, universes] : systMap)
    {
        const auto iter = dimensions.find(paramName);
        if (iter == dimensions.end())
            throw std::invalid_argument("CrossSectionHelperNew::ValidateSystMap - Unexpected parameter name: " + paramName);

        if (universes.size() != iter->second)
            throw std::invalid_argument("CrossSectionHelperNew::ValidateSystMap - Unexpected number of universes for parameter : " + paramName);
    }
}


} // namespace ubcc1pi

#endif
