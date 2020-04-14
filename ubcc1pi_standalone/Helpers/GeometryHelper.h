/**
 *  @file  ubcc1pi_standalone/Helpers/GeometryHelper.h
 *
 *  @brief The header file for the geometry helper class
 */

#ifndef UBCC1PI_STANDALONE_HELPERS_GEOMETRY_HELPER
#define UBCC1PI_STANDALONE_HELPERS_GEOMETRY_HELPER

namespace ubcc1pi
{

/**
 *  @brief  The geometry helper class
 */
class GeometryHelper
{
    public:
        static float lowX;  ///< The low X border of the detector
        static float highX; ///< The low X border of the detector
        static float lowY;  ///< The low Y border of the detector
        static float highY; ///< The low Y border of the detector
        static float lowZ;  ///< The low Z border of the detector
        static float highZ; ///< The low Z border of the detector
};

// -----------------------------------------------------------------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------------------------------------------------------------

float GeometryHelper::lowX  = 0.f;
float GeometryHelper::highX = 256.35f;

float GeometryHelper::lowY  = -116.5f;
float GeometryHelper::highY = 116.5f;

float GeometryHelper::lowZ  = 0.f;
float GeometryHelper::highZ = 1036.8f;
    

} // namespace ubcc1pi

#endif
