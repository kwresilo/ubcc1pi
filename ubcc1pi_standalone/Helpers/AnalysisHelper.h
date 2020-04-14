/**
 *  @file  ubcc1pi_standalone/Helpers/AnalysisHelper.h
 *
 *  @brief The header file for the analysis helper class
 */

#ifndef UBCC1PI_STANDALONE_HELPERS_ANALYSIS_HELPER
#define UBCC1PI_STANDALONE_HELPERS_ANALYSIS_HELPER

#include "ubcc1pi_standalone/Interface/Event.h"

#include <memory>

#include <TVector3.h>

namespace ubcc1pi
{

/**
 *  @brief  The analysis helper class
 */
class AnalysisHelper
{
    public:
        /**
         *  @brief  Determine if the input event is truly a fiducial CC1Pi event
         *
         *  @param  pEvent the input event
         *
         *  @return boolean, true if CC1Pi
         */
        static bool IsTrueCC1Pi(const std::shared_ptr<Event> &pEvent);

        /**
         *  @brief  Determine if the input truth particle is deemed visible and possibly is above a momentum threshold
         *
         *  @param  particle the input truth particle
         *
         *  @return boolean, true if visible
         */
        static bool PassesVisibilityThreshold(const Event::Truth::Particle &particle);

        /**
         *  @brief  Determine if a given point is in the volume defined by the input margins
         *
         *  @param  point the point
         *  @param  lowXMargin the fiducial margin from the low X face
         *  @param  highXMargin the fiducial margin from the high X face
         *  @param  lowYMargin the fiducial margin from the low Y face
         *  @param  highYMargin the fiducial margin from the high Y face
         *  @param  lowZMargin the fiducial margin from the low Z face
         *  @param  highZMargin the fiducial margin from the high Z face
         *
         *  @return boolean, true if the input point is within the specified margins of the detector faces
         */
        static bool IsPointWithinMargins(const TVector3 &point, const float lowXMargin, const float highXMargin, const float lowYMargin, const float highYMargin, const float lowZMargin, const float highZMargin);

        /**
         *  @brief  Determine if a given point is in the fiducial volume
         *
         *  @param  point the point
         *
         *  @return boolean, true if in fiducial volume
         */
        static bool IsFiducial(const TVector3 &point);
        
        /**
         *  @brief  Determine if a given point is contained within the detector
         *
         *  @param  point the point
         *
         *  @return boolean, true if contained
         */
        static bool IsContained(const TVector3 &point);
        
        /**
         *  @brief  Determine if a given truth particle is contained 
         *
         *  @param  particle
         *
         *  @return boolean, true if contained
         */
        static bool IsContained(const Event::Truth::Particle &particle);
        
        /**
         *  @brief  Determine if a given reco particle is contained
         *
         *  @param  particle
         *
         *  @return boolean, true if contained
         */
        static bool IsContained(const Event::Reco::Particle &particle);
        
        /**
         *  @brief  Get the index of the truth particle that best matches to the input reco particle
         *
         *  @param  recoParticle the input reco particle
         *  @param  truthParticles the input list of all truth particles from the event
         *  @param  applyVisibilityThreshold whether to apply the visibility threshold when considering possible matches
         *
         *  @return the output truth particle index
         */
        static unsigned int GetBestMatchedTruthParticleIndex(const Event::Reco::Particle &recoParticle, const std::vector<Event::Truth::Particle> &truthParticles, const bool applyVisibilityThreshold = true);

        /**
         *  @brief  Get the truth particle that best matches to the input reco particle
         *
         *  @param  recoParticle the input reco particle
         *  @param  truthParticles the input list of all truth particles from the event
         *  @param  applyVisibilityThreshold whether to apply the visibility threshold when considering possible matches
         *
         *  @return the output truth particle
         */
        static Event::Truth::Particle GetBestMatchedTruthParticle(const Event::Reco::Particle &recoParticle, const std::vector<Event::Truth::Particle> &truthParticles, const bool applyVisibilityThreshold = true);

        /**
         *  @brief  Determin if the input particle is golden: no scatters, stopping and contained
         *
         *  @param  particle the input truth particle
         *
         *  @return boolean, true if golden
         */
        static bool IsGolden(const Event::Truth::Particle &particle);

        /**
         *  @brief  Print a loading bar showing the numerator/denominator ratio
         *
         *  @param  numerator the numerator
         *  @param  denominator the denominator
         */
        static void PrintLoadingBar(const unsigned int numerator, const unsigned int denominator);
};

// -----------------------------------------------------------------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------------------------------------------------------------

} // namespace ubcc1pi

#endif
