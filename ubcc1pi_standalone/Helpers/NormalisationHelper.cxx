#include "ubcc1pi_standalone/Helpers/NormalisationHelper.h"

#include <stdexcept>

namespace ubcc1pi
{
        
float NormalisationHelper::GetOverlaysNormalisation(const Config &config)
{
    if (config.norms.overlaysPOT <= std::numeric_limits<float>::epsilon())
        throw std::invalid_argument("NormalisationHelper::GetOverlaysNormalisation - Overlay POT is invalid");

    return config.norms.dataBNBTor875WCut / config.norms.overlaysPOT;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

float NormalisationHelper::GetDirtNormalisation(const Config &config)
{
    if (config.norms.dirtPOT <= std::numeric_limits<float>::epsilon())
        throw std::invalid_argument("NormalisationHelper::GetDirtNormalisation - Dirt POT is invalid");
    
    return config.norms.dataBNBTor875WCut / config.norms.dirtPOT;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

float NormalisationHelper::GetDataEXTNormalisation(const Config &config)
{
    if (config.norms.dataEXTTriggers <= std::numeric_limits<float>::epsilon())
        throw std::invalid_argument("NormalisationHelper::GetDataEXTNormalisation - Data EXT trigger count is invalid");

    return config.norms.dataBNBE1DCNTWCut / config.norms.dataEXTTriggers;
}


}
