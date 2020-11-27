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
        *  @brief  Like a SystMap but using a map instead of a vector as the mapped type.
        *          In this way it's possible for the map to have empty entries for a given universe index
        *
        *  @tparam T the mapped value type
        */
        template <typename T>
        using SystPartialMap = std::unordered_map< std::string, std::unordered_map<unsigned int, T> >;

        /**
        *  @brief  A mapping from a systematic parameter name to the number of universes
        */
        typedef std::unordered_map< std::string, unsigned int > SystDimensionsMap;

        /**
        *  @brief  A mapping from a systematic paramter name to a vector of floats, one per universe
        */
        typedef SystMap<float> SystFloatMap;

        /**
        *  @brief  A mapping from a systematic parameter name to a vector of 1D histograms, one per universe
        */
        typedef SystMap< std::shared_ptr<TH1F> > SystTH1FMap;

        /**
        *  @brief  A wrapper around an arbitrary function that caches the result of the function in a SystPartialMap.
        *          If the function has already been executed for a given systematic parameter / universe then the a call to the function
        *          will return the cached value instead. This is provides a helful way of storing the result of a function that we might
        *          want to use in multiple places.
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
                *  @brief  Get the flux distribution in each systematic universe
                *
                *  @return the varied fluxes
                */
                SystTH1FMap GetFluxVariations() const;

                /**
                *  @brief  Get the integrated flux in the nominal universe
                *
                *  @return the integrated nominal flux
                */
                float GetIntegratedNominalFlux() const;

                /**
                *  @brief  Get the integrated flux in each systematic universe
                *
                *  @return the varied integrated fluxes
                */
                SystFloatMap GetIntegratedFluxVariations() const;

            private:

                /**
                *  @brief  Get the integrated flux
                *
                *  @param  pFlux the input flux (in energy bins)
                *
                *  @return the integrated flux
                */
                float GetIntegratedFlux(const std::shared_ptr<TH1F> &pFlux) const;

                std::vector<float>    m_binEdges;               ///< The bin edges
                SystDimensionsMap     m_dimensions;             ///< The systematic paramter dimensions
                std::shared_ptr<TH1F> m_pFluxNominal;           ///< The nominal flux distribution
                std::shared_ptr<TH1F> m_pSpectrumNominal;       ///< The event rate spectrum in the nominal universe
                SystTH1FMap           m_spectrumVariations;     ///< The event rate spectrum in each systematic universe
        };

        /**
        *  @brief  The cross-section class. Describes a forward-folded differential cross-section.
        */
        class CrossSection
        {
            public:
                /**
                *  @brief  Constructor
                *
                *  @param  binEdges the bin edges for the kinematic quantity
                *  @param  hasUnderflow if there is an underflow bin
                *  @param  hasOverflow if there is an overflow bin
                */
                CrossSection(const std::vector<float> &binEdges, const bool hasUnderflow, const bool hasOverflow);

            private:

                std::vector<float>  m_binEdges; ///< The bin edges
                ubsmear::UBXSecMeta m_metadata; ///< The cross-section metadata
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
        *  @brief  Get a SystTH1FMap with the supplied dimensions containing new TH1F objects with the supplied binning
        *
        *  @param  binEdges the bin edges
        *  @param  dimensions the systematic dimensions map
        *
        *  @return the SystTH1FMap
        */
        static SystTH1FMap GetSystTH1FMap(const std::vector<float> &binEdges, const SystDimensionsMap &dimensions);

        /**
        *  @brief  Get the weights corresponding to the systematic parameters in the input dimensions object
        *
        *  @param  truth the input truth information
        *  @param  dimensions the dimensions of the weights map to obtain
        *
        *  @return the weights map
        */
        static SystFloatMap GetWeightsMap(const Event::Truth &truth, const SystDimensionsMap &dimensions);

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
