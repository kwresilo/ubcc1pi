/**
 *  @file  ubcc1pi_standalone/Helpers/AnalysisHelper.h
 *
 *  @brief The header file for the analysis helper class
 */

#ifndef UBCC1PI_STANDALONE_HELPERS_ANALYSIS_HELPER
#define UBCC1PI_STANDALONE_HELPERS_ANALYSIS_HELPER

#include "ubcc1pi_standalone/Interface/Event.h"

#include <memory>
#include <unordered_map>

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
         *  @brief  The sample type enumeration
         */
        enum SampleType
        {
            DataBNB,
            DataEXT,
            Overlay,
            Dirt
        };

        /**
         *  @brief  The vector of all sample types
         */
        static const std::vector<SampleType> AllSampleTypes;

        /**
         *  @brief  Get the sample type name string
         *
         *  @param  sampleType the input sample type
         *
         *  @return the sample type name
         */
        static std::string GetSampleTypeName(const SampleType &sampleType);

        /**
         *  @brief  The event counter class
         */
        class EventCounter
        {
            public:
                /**
                 *  @brief  Constructor
                 */
                EventCounter();

                /**
                 *  @brief  Count the input event and assign it a tag for reference
                 *
                 *  @param  tag the input tag, this can be anything - e.g. "passingPreselection" or "finalSelection". Use the special tag "all" for all input events
                 *  @param  sampleType the type of sample from which the event came
                 *  @param  pEvent the input event to count
                 *  @param  weight a weight to apply to the event
                 */
                void CountEvent(const std::string &tag, const SampleType &sampleType, const std::shared_ptr<Event> &pEvent, const float weight = 1.f);
                
                /**
                 *  @brief  Get the total weight for a given tag, sample and classification - if the entry doesn't exist a weight of zero is returned
                 *
                 *  @param  tag the tag
                 *  @param  sampleType the sample type
                 *  @param  classification the classification
                 *
                 *  @return the total weight
                 */
                float GetWeight(const std::string &tag, const SampleType &sampleType, const std::string &classification) const;

                /**
                 *  @brief  Get the total weight for a given tag, sample and classification
                 *
                 *  @param  tag the tag
                 *  @param  sampleType the sample type
                 *  @param  classification the classification
                 *  @param  weight the output weight
                 *
                 *  @return boolean, true if an entry exists
                 */
                bool GetWeight(const std::string &tag, const SampleType &sampleType, const std::string &classification, float &weight) const;
                
                /**
                 *  @brief  Get the total weight for a given tag over all non beam data sample types at a given tag
                 *
                 *  @param  tag the tag
                 *
                 *  @return the total weight
                 */
                float GetTotalMCWeight(const std::string &tag) const;

                /**
                 *  @brief  Get the total weight for BNB data with the given tag
                 *
                 *  @param  tag the tag
                 *
                 *  @return the weight
                 */
                float GetBNBDataWeight(const std::string &tag) const;
                
                /**
                 *  @brief  Get the total weight for signal events with the given tag
                 *
                 *  @param  tag the tag
                 *  @param  searchQuery additional string that must be found in the classification name to count as signal
                 *
                 *  @return the weight
                 */
                float GetSignalWeight(const std::string &tag, const std::string &searchQuery = "") const;
                
                /**
                 *  @brief  Get the total weight for background events with the given tag
                 *
                 *  @param  tag the tag
                 *  @param  searchQuery additional string that must be found in the classification name to count as signal
                 *
                 *  @return the weight
                 */
                float GetBackgroundWeight(const std::string &tag, const std::string &searchQuery = "") const;
                
                /**
                 *  @brief  Get the purity for signal events at a given tag
                 *
                 *  @param  tag the input tag
                 *  @param  searchQuery additional string that must be found in the classification name to count as signal
                 *
                 *  @return the purity
                 */
                float GetSignalPurity(const std::string &tag, const std::string &searchQuery = "") const;
                
                /**
                 *  @brief  Get the efficiency for signal events at the given tag
                 * 
                 *  @param  tag the input tag
                 *  @param  searchQuery additional string that must be found in the classification name to count as signal
                 *
                 *  @return the efficiency
                 */
                float GetSignalEfficiency(const std::string &tag, const std::string &searchQuery = "") const;

                /**
                 *  @brief  Get the list of signal classifications
                 *
                 *  @param  searchQuery additional string that must be found in the classification name to count as signal
                 *
                 *  @return classifications meeting criteria 
                 */
                std::vector<std::string> GetSignalClassifications(const std::string &searchQuery = "") const;

                /**
                 *  @brief  Determine if a given classification representes a signal event
                 *
                 *  @param  classification the input classification
                 *  @param  searchQuery additional string that must be found in the classification name to count as signal
                 *
                 *  @return boolea, true if signal
                 */
                bool IsSignalClassification(const std::string &classification, const std::string &searchQuery = "") const;

                /**
                 *  @brief  Get the purity for the specified subsample at the given tag
                 *
                 *  @param  tag the input tag
                 *  @param  sampleType the sample type
                 *  @param  classification the classification
                 *
                 *  @return the purity
                 */
                float GetPurity(const std::string &tag, const SampleType &sampleType, const std::string &classification) const;
                
                /**
                 *  @brief  Get the efficiency for the specified subsample at the given tag
                 *
                 *  @param  tag the input tag
                 *  @param  sampleType the sample type
                 *  @param  classification the classification
                 *
                 *  @return the efficiency
                 */
                float GetEfficiency(const std::string &tag, const SampleType &sampleType, const std::string &classification) const;

                /**
                 *  @brief  Get the tags
                 *
                 *  @return the tags
                 */
                std::vector<std::string> GetTags() const;

                /**
                 *  @brief  Print a breakdown of the counted events tag-by-tag broken down into signal and background
                 *
                 *  @param  outputFileName the name of the output file in which we should save the table
                 */
                void PrintBreakdownSummary(const std::string &outputFileName) const;
                
                /**
                 *  @brief  Print a breakdown of the counted events tag-by-tag broken down into individual classifications
                 *
                 *  @param  outputFileName the name of the output file in which we should save the table
                 *  @param  nEntries the number of backgrounds to print for each tag
                 */
                void PrintBreakdownDetails(const std::string &outputFileName, const unsigned int nEntries = 10u) const;

            private:

                /**
                 *  @brief  Define a 3-key map template for brevity
                 *          key 1 : tag string
                 *          key 2 : sample type enum
                 *          key 3 : event classification string (e.g. event topology)
                 *
                 *  @tparam T the mapped type
                 */
                template <typename T>
                using CounterMap = std::unordered_map<std::string, std::unordered_map<SampleType, std::unordered_map<std::string, T> > >; 

                std::vector<std::string> m_tags;            ///< The tags
                std::vector<std::string> m_classifications; ///< The classifications
                CounterMap<float>        m_eventWeightMap;  ///< The cumulative event weight
        };

        /**
         *  @brief  Get the nominal event weight
         *
         *  @param  pEvent the input event
         *
         *  @return the nominal event weight
         */
        static float GetNominalEventWeight(const std::shared_ptr<Event> &pEvent);

        /**
         *  @brief  Determine if the input event is truly a fiducial CC1Pi event
         *
         *  @param  pEvent the input event
         *  @param  useAbsPdg if we should use the absolute values of PDG codes
         *
         *  @return boolean, true if CC1Pi
         */
        static bool IsTrueCC1Pi(const std::shared_ptr<Event> &pEvent, const bool useAbsPdg);

        /**
         *  @brief  Determine if the input truth particle is deemed visible and possibly is above a momentum threshold
         *
         *  @param  particle the input truth particle
         *
         *  @return boolean, true if visible
         */
        static bool PassesVisibilityThreshold(const Event::Truth::Particle &particle);

        /**
         *  @brief  Select from the input truth particles those which pass the visibility threshold
         *
         *  @param  particles the input particles
         *
         *  @return the visible particles
         */
        static std::vector<Event::Truth::Particle> SelectVisibleParticles(const std::vector<Event::Truth::Particle> &particles);

        /**
         *  @brief  Count the number of particles in the input vector with the supplied PDG code
         *
         *  @param  particles the input particles
         *  @param  pdgCode the pdg code
         *  @param  useAbsPdg if we should group particles with the same absolute PDG code
         *
         *  @return the number of particles with the pdg code
         */
        static unsigned int CountParticlesWithPdgCode(const std::vector<Event::Truth::Particle> &particles, const int pdgCode, const bool useAbsPdg);
        
        /**
         *  @brief  Count the number of particles in the input vector with the supplied PDG code that are golden
         *
         *  @param  particles the input particles
         *  @param  pdgCode the pdg code
         *  @param  useAbsPdg if we should group particles with the same absolute PDG code
         *
         *  @return the number of particles with the pdg code that are golden
         */
        static unsigned int CountGoldenParticlesWithPdgCode(const std::vector<Event::Truth::Particle> &particles, const int pdgCode, const bool useAbsPdg);

        /**
         *  @brief  Get the mapping from PDG code to the number of particles with that PDG code in the input vector
         *
         *  @param  particles the input particles
         *  @param  useAbsPdg if we should group particles with the same absolute PDG code
         *  @param  foundPdgs the output vector of PDGs found (ordered)
         *  @param  pdgCodeCountMap the output map from PDG codes to counts
         */
        static void GetPdgCodeCountMap(const std::vector<Event::Truth::Particle> &particles, const bool useAbsPdg, std::vector<int> &foundPdgs, std::unordered_map<int, unsigned int> &pdgCodeCountMap);

        /**
         *  @brief  Get the topology string for the input particles, e.g. "1 Mu-  1 Pi+  X p"
         *
         *  @param  particles the input particles
         *  @param  useAbsPdg if we should group particles with the same absolute PDG code
         *  @param  countProtonsInclusively whether to count protons inclusively (with an X) or exclusively
         *
         *  @return the topology string
         */
        static std::string GetTopologyString(const std::vector<Event::Truth::Particle> &particles, const bool useAbsPdg, const bool countProtonsInclusively);

        /**
         *  @brief  Get a string which classifies the input event using truth information
         *
         *  @param  pEvent the input event
         *  @param  useAbsPdg if we should group particles with the same absolute PDG code
         *  @param  countProtonsInclusively whether to count protons inclusively (with an X) or exclusively
         *
         *  @return the classification string
         */
        static std::string GetClassificationString(const std::shared_ptr<Event> &pEvent, const bool useAbsPdg, const bool countProtonsInclusively);

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
         *  @brief  Determine if the input reco particle has a fitted track
         *
         *  @param  particle the particle
         *
         *  @return boolean, true if track fit information is available
         */
        static bool HasTrackFit(const Event::Reco::Particle &particle);
        
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
         *  @brief  Get the momentum of a pion from range
         *
         *  @param  range the input range
         *
         *  @return the momentum
         */
        static float GetPionMomentumFromRange(const float &range);
        
        /**
         *  @brief  Get the momentum of a muon from range
         *
         *  @param  range the input range
         *
         *  @return the momentum
         */
        static float GetMuonMomentumFromRange(const float &range);

        /**
         *  @brief  Get the ratio of two likelihoods
         *
         *  @param  numerator the numerator likelihood member 
         *  @param  denominator the denominator likelihood member
         *  @param  ratio the output ratio
         *
         *  @return boolean, true if the denominator was non zero and both numerator and denominator are available
         */
        static bool GetLogLikelihoodRatio(const Member<float> &numerator, const Member<float> &denominator, float &ratio);

        /**
         *  @brief  Get the softmax of two variables: exp(s) / (exp(s) + exp(b))
         *
         *  @param  signal the signal variable member
         *  @param  background the background variable member
         *  @param  softmax the output softmax
         *
         *  @return boolean, true if the value is calculable
         */
        static bool GetSoftmax(const Member<float> &signal, const Member<float> &background, float &softmax);

        /**
         *  @brief  Get the uncertainty on a count
         *          The uncertainty is calculated using a poisson distribution with uniform prior
         *
         *  @param  count the count
         *
         *  @return the uncertainty
         */
        static float GetCountUncertainty(const float &count);

        /**
         *  @brief  Get the uncertainty on an efficiency: (numerator, n) / (denominator, d)
         *          The uncertainty is calculated assuming a binomial distribution with a uniform prior
         *
         *  @param  numerator the numerator
         *  @param  denominator the denominator
         *
         *  @return the uncertainty
         */
        static float GetEfficiencyUncertainty(const float &numerator, const float &denominator);

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
        
const std::vector<AnalysisHelper::SampleType> AnalysisHelper::AllSampleTypes = {
    DataBNB,
    Overlay,
    DataEXT,
    Dirt
};

} // namespace ubcc1pi

#endif
