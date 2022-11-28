/**
 *  @file  ubcc1pi_standalone/Helpers/ExtractionHelper.h
 *
 *  @brief The header file for the extraction helper class
 */

#ifndef UBCC1PI_STANDALONE_HELPERS_EXTRACTION_HELPER
#define UBCC1PI_STANDALONE_HELPERS_EXTRACTION_HELPER

#include "ubcc1pi_standalone/Objects/Config.h"

namespace ubcc1pi
{

    /**
     *  @brief  The extraction helper class
     */
    class ExtractionHelper
    {
        public:
            /**
             *  @brief  A vector of tuples with 4 entries
            *           - First, the sample type (e.g. overlay)
            *           - Second, a string which is used to identify a given detector variation sample (for other sample type, this is unused)
            *           - Third, the path to the input file
            *           - Fourth, the normalisation factor to apply to all events in that file
            */
            typedef std::vector< std::tuple<AnalysisHelper::SampleType, std::string, std::string, float> > InputFileList;

            /**
             *  @brief   A map from name of cross-section to function which gets relevant kinematic quanitity from an input analysis data object
            */
            typedef std::unordered_map< std::string, std::function<float(const AnalysisHelper::AnalysisData &)> > AnalysisValueMap;

            /**
             *  @brief  Get the normalisation factor for the dirt
             *
             *  @param  config the input configuration
             *  @param  inputData the vector of input files to be filled
             *  @param  totalExposurePOT the total exposure in POT is returned with this variable
             * 
             */
            static void PopulateInputFileList(const Config &config, InputFileList &inputData, float &totalExposurePOT);

            /**
             *  @brief  Setup the relevent "getters" for each cross-section by populating an AnalysisValueMap
             *
             *  @param  analysisValueMap the AnalysisValueMap to be populated
             *  @param  createSidebandVersion populate with the sideband (CC0pi) functions instead of normal CC1pi functions 
             */
            static void PopulateAnalysisValueMap(AnalysisValueMap &analysisValueMap, const bool createSidebandVersion=false);

            /**
             *  @brief  Get the selected events, cross-sections and uncertainties and save the results to txt files
             *
             *  @param  xsec the cross-section to be extracted
             *  @param  scalingData the scaling data to be used 
             *  @param  selectionName the name of the selection (generic/golden) to be used
             *  @param  name the name of the cross-section
             *  @param  postfix additional string to be added to end of txt file names 
             */
            static void SaveCrossSectionMatricies(const CrossSectionHelper::CrossSection &xsec, const CrossSectionHelper::CrossSection::ScalingData &scalingData, 
                const std::string &selectionName, const std::string &name, const std::string &postfix = std::string(""), const bool disableUncertainties=false);
    };

} // namespace ubcc1pi

#endif
