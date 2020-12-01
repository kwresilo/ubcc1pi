/**
 *  @file  ubcc1pi_standalone/Helpers/NormalisationHelper.h
 *
 *  @brief The header file for the normalisation helper class
 */

#ifndef UBCC1PI_STANDALONE_HELPERS_NORMALISATION_HELPER
#define UBCC1PI_STANDALONE_HELPERS_NORMALISATION_HELPER

#include "ubcc1pi_standalone/Objects/Config.h"

namespace ubcc1pi
{
    /**
     *  @brief  The normalisation helper class
     */
    class NormalisationHelper
    {
        public:
            /**
             *  @brief  Get the normalisation factor for the overlays
             *
             *  @param  config the input configuration
             *
             *  @return the normalisation factor
             */
            static float GetOverlaysNormalisation(const Config &config);

            /**
             *  @brief  Get the normalisation factor for the dirt
             *
             *  @param  config the input configuration
             *
             *  @return the normalisation factor
             */
            static float GetDirtNormalisation(const Config &config);

            /**
             *  @brief  Get the normalisation factor for the EXT data
             *
             *  @param  config the input configuration
             *
             *  @return the normalisation factor
             */
            static float GetDataEXTNormalisation(const Config &config);

            /**
             *  @brief  Get the normalisation factor for a detector variation file
             *
             *  @param  config the input configuration
             *  @param  paramName the varied detector parmater or central value sample name
             *
             *  @return the normalisation factor
             */
            static float GetDetectorVariationNormalisation(const Config &config, const std::string &paramName);
    };

} // namespace ubcc1pi

#endif
