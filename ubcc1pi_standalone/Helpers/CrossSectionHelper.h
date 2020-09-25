/**
 *  @file  ubcc1pi_standalone/Helpers/CrossSectionHelper.h
 *
 *  @brief The header file for the cross section helper class
 */

#ifndef UBCC1PI_STANDALONE_HELPERS_CROSS_SECTION_HELPER
#define UBCC1PI_STANDALONE_HELPERS_CROSS_SECTION_HELPER

#include "ubcc1pi_standalone/Interface/Event.h"

#include <string>
#include <unordered_map>
#include <map>
#include <vector>
#include <memory>
#include <TH1F.h>
#include <TH2F.h>

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
            
            public:

                /**
                 *  @brief  Constructor
                 *
                 *  @param  binEdges the bin edges of the cross-section
                 *  @param  hasUnderflow if the first bin is an underflow bin
                 *  @param  hasOverflow if the last bin is an overflow bin
                 *  @param  scaleByBinWidth if we should should divide by the bin-width when calculating the cross-section
                 *  @param  systUniverseSizesMap the input mapping from the systematic parameters to the number of universes to expect
                 */
                CrossSection(const std::vector<float> &binEdges, const bool hasUnderflow, const bool hasOverflow, const bool scaleByBinWidth, const SystematicUniverseSizesMap &systUniverseSizesMap);

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

            private:

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
                 *  @brief  Check if a bin index is an underflow or an overflow bin
                 *
                 *  @param  binIndex the ROOT bin index (enumeratred from 1)
                 *
                 *  @return if the bin is under/overflow
                 */
                bool IsUnderOverflowBin(const unsigned int binIndex) const;

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
                 *
                 *  @return the cross-section
                 */
                std::shared_ptr<TH1F> GetCrossSection(const std::shared_ptr<TH1F> &selected, const std::shared_ptr<TH1F> &background, const std::shared_ptr<TH1F> &smearedEff) const;

                std::vector<float>                  m_binEdges;                   ///< The bin edges
                bool                                m_hasUnderflow;               ///< If the first bin is an underflow bin
                bool                                m_hasOverflow;                ///< If the last bin is an overflow bin
                bool                                m_scaleByBinWidth;            ///< If we should scale by bin width
                SystematicUniverseSizesMap          m_systUniverseSizesMap;       ///< The number of universes for each systematic paramter

                SystMap<std::shared_ptr<TH2F> >     m_signalSelectedRecoTrue;     ///< The 2D histograms of selected signal events with reco & true bin indices for each universe
                SystMap<std::shared_ptr<TH1F> >     m_signalAllTrue;              ///< The 1D histograms of all signal events with true bin indices for each universe
                SystMap<std::shared_ptr<TH1F> >     m_backgroundSelectedReco;     ///< The 1D histograms of selected background events with reco bin indices for each universe

                std::shared_ptr<TH1F>               m_dataSelectedReco;           ///< The 1D histogram of selected BNB data events with reco bin indices for each universe
                std::shared_ptr<TH2F>               m_signalSelectedNomRecoTrue;  ///< The 1D histogram of selected BNB data events with reco bin indices in the nominal universe
                std::shared_ptr<TH1F>               m_signalAllNomTrue;           ///< The 1D histogram of all signal events with true bin indices in the nominal universe
                std::shared_ptr<TH1F>               m_backgroundSelectedNomReco;  ///< The 1D histogram of selected background events with reco bin indices in the nominal universe

                static unsigned int                 m_histCount;                  ///< A counter for the total number of histograms - used to avoid name collisions
        };

        // ---------------------------------------------------------------------------------------------------------------------------------
        
        /**
         *  @brief  Get the systematic weights from the event truth information in a map from paramter name to the universe weights
         *
         *  @param  truth the event truth object
         *
         *  @return the systematic weights map
         */
        static SystematicWeightsMap GetSystematicWeightsMap(const Event::Truth &truth);

        /**
         *  @brief  Add the specified number of boostrap universes to the input systematic weights map
         *
         *  @param  nBootstrapUniverses the number of boostrap universes to use
         *  @param  systWeightsMap the systematic weights map to update
         */
        static void AddBootstrapWeights(const unsigned int nBootstrapUniverses, SystematicWeightsMap &systWeightsMap);

        /**
         *  @brief  Get a systematic weights map with the supplied dimensions with a value of one for each weight
         *
         *  @param  systUniverseSizesMap the mapping from systematic parameter to the number of universes
         *
         *  @return the systematic weights map
         */
        static SystematicWeightsMap GetUnitSystematicWeightsMap(const SystematicUniverseSizesMap &systUniverseSizesMap);

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
};

// -----------------------------------------------------------------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------------------------------------------------------------

unsigned int CrossSectionHelper::CrossSection::m_histCount = 0u;

} // namespace ubcc1pi

#endif
