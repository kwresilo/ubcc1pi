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
                 *  @brief  Count the input event and assign it a tag for reference
                 *
                 *  @param  tag the input tag, this can be anything - e.g. "passingPreselection" or "finalSelection". Use the special tag "all" for all input events
                 *  @param  sampleType the type of sample from which the event came
                 *  @param  pEvent the input event to count
                 *  @param  weight a weight to apply to the event
                 */
                void CountEvent(const std::string &tag, const SampleType &sampleType, const std::shared_ptr<Event> &pEvent, const float weight = 1.f);

                /**
                 *  @brief  Print a breakdown of the counted events
                 *  
                 *  @param  nBackgrounds the number of backgrounds to print
                 */
                void PrintBreakdown(const unsigned int nBackgrounds = 10u) const;
                
                /**
                 *  @brief  Print a breakdown of the counted events with weights applied
                 *
                 *  @param  nBackgrounds the number of backgrounds to print
                 */
                void PrintBreakdownWithWeights(const unsigned int nBackgrounds = 10u) const;

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

                /**
                 *  @brief  Print the breakdown of the input counter map
                 *
                 *  @tparam T the mapped type
                 *  @param  counterMap the counter map to print
                 *  @param  nBackgrounds the number of backgrounds to print
                 */
                template <typename T>
                void PrintBreakdown(const CounterMap<T> &counterMap, const unsigned int nBackgrounds) const;

                std::vector<std::string> m_tags;           ///< The tags
                CounterMap<unsigned int> m_eventCountMap;  ///< The cumulative event count
                CounterMap<float>        m_eventWeightMap; ///< The cumulative event weight
        };

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
         *
         *  @return the number of particles with the pdg code
         */
        static unsigned int CountParticlesWithPdgCode(const std::vector<Event::Truth::Particle> &particles, const int pdgCode);
        
        /**
         *  @brief  Count the number of particles in the input vector with the supplied PDG code that are golden
         *
         *  @param  particles the input particles
         *  @param  pdgCode the pdg code
         *
         *  @return the number of particles with the pdg code that are golden
         */
        static unsigned int CountGoldenParticlesWithPdgCode(const std::vector<Event::Truth::Particle> &particles, const int pdgCode);

        /**
         *  @brief  Get the mapping from PDG code to the number of particles with that PDG code in the input vector
         *
         *  @param  particles the input particles
         *  @param  foundPdgs the output vector of PDGs found (ordered)
         *  @param  pdgCodeCountMap the output map from PDG codes to counts
         */
        static void GetPdgCodeCountMap(const std::vector<Event::Truth::Particle> &particles, std::vector<int> &foundPdgs, std::unordered_map<int, unsigned int> &pdgCodeCountMap);

        /**
         *  @brief  Get the topology string for the input particles, e.g. "1 Mu-  1 Pi+  X p"
         *
         *  @param  particles the input particles
         *  @param  countProtonsInclusively whether to count protons inclusively (with an X) or exclusively
         *
         *  @return the topology string
         */
        static std::string GetTopologyString(const std::vector<Event::Truth::Particle> &particles, const bool countProtonsInclusively = true);

        /**
         *  @brief  Get a string which classifies the input event using truth information
         *
         *  @param  pEvent the input event
         *  @param  countProtonsInclusively whether to count protons inclusively (with an X) or exclusively
         *
         *  @return the classification string
         */
        static std::string GetClassificationString(const std::shared_ptr<Event> &pEvent, const bool countProtonsInclusively = true);

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
         *  @brief  Get the ratio of two likelihoods
         *
         *  @param  numerator the numerator likelihood member 
         *  @param  denominator the denominator likelihood member
         *  @param  ratio the output ratio
         *
         *  @return boolean, true if the denominator was non zero and both numerator and denominator are available
         */
        static bool GetLikelihoodRatio(const Member<float> &numerator, const Member<float> &denominator, float &ratio);

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
    DataEXT,
    Overlay,
    Dirt
};

// -----------------------------------------------------------------------------------------------------------------------------------------

// TODO make this less repulsive
template <typename T>
void AnalysisHelper::EventCounter::PrintBreakdown(const CounterMap<T> &counterMap, const unsigned int nBackgrounds) const
{
    const auto lineSize = 87u;
    const auto alignWidthName = 40u;
    const auto alignWidthValue = 10u;

    // Check for the special "all" tag to use as a baseline
    const auto allIter = counterMap.find("all");
    const auto hasAllTag = (allIter != counterMap.end());

    // Reverse the tags to get chronological order
    std::vector<std::string> reverseTags = m_tags;
    std::reverse(reverseTags.begin(), reverseTags.end());
    for (unsigned int iTag = 0; iTag < reverseTags.size(); iTag++)
    {
        const auto &tag = reverseTags.at(iTag);

        // Print the header for this tag
        std::cout << std::string(lineSize, '=') << std::endl;
        std::cout << "Tag : " << tag << std::endl;
        std::cout << std::string(lineSize, '=') << std::endl;

        // Check if we have entries to print
        const auto iter = counterMap.find(tag);
        if (iter == counterMap.end())
        {
            std::cout << "No entries" << std::endl;
            continue;
        }

        const auto sampleTypeMap = iter->second;
        
        // Collect the total count over all sample types
        T totalCountMC = static_cast<T>(0);
        T totalCountDataBNB = static_cast<T>(0);
        std::string mcSamplesUsed = "";
        for (const auto &sampleType : AnalysisHelper::AllSampleTypes)
        {
            const auto iter2 = sampleTypeMap.find(sampleType);
            if (iter2 == sampleTypeMap.end())
                continue;
            
            const auto classificationToCountMap = iter2->second;

            if (classificationToCountMap.empty())
                continue;

            for (const auto &entry : classificationToCountMap)
            {
                if (sampleType == DataBNB)
                {
                    totalCountDataBNB += entry.second;
                }
                else
                {
                    totalCountMC += entry.second;
                }
            }
          
            if (sampleType != DataBNB)
                mcSamplesUsed += AnalysisHelper::GetSampleTypeName(sampleType) + "  ";
        }
        std::cout << "Total Data BNB: " << totalCountDataBNB << std::endl;
        std::cout << "Total \"MC\": " << totalCountMC << std::endl;
        std::cout << "  - This is the denominator of the purity for: " << mcSamplesUsed << std::endl;

        // Breakdown by sample type
        for (const auto &sampleType : AnalysisHelper::AllSampleTypes)
        {
            const auto iter2 = sampleTypeMap.find(sampleType);
            if (iter2 == sampleTypeMap.end())
                continue;
            
            const auto classificationToCountMap = iter2->second; 

            const auto sampleString = AnalysisHelper::GetSampleTypeName(sampleType);

            std::cout << std::string(lineSize, '-') << std::endl;
            std::cout << "Sample : " << sampleString << std::endl;
            std::cout << std::string(lineSize, '-') << std::endl;
            
            if (classificationToCountMap.empty())
            {
                std::cout << "No entries" << std::endl;
                continue;
            }

            // Separate the map entries into signal an backgrounds
            std::vector<std::pair<std::string, T> > signals, backgrounds;
            for (const auto &entry : classificationToCountMap)
            {
                const auto &classification = entry.first;
                const auto count = entry.second;

                // Check for signal and print it, we want to see all signal
                if (!classification.empty() && classification.at(0) == 'S')
                {
                    signals.emplace_back(classification, count);
                }
                else
                {
                    backgrounds.emplace_back(classification, count);
                }
            }

            // Sort by the count
            std::sort(signals.begin(), signals.end(), [](const auto &a, const auto &b){return a.second > b.second;});
            std::sort(backgrounds.begin(), backgrounds.end(), [](const auto &a, const auto &b){return a.second > b.second;});

            // Combine into a vector we want to print
            const auto nBackgroundsToPrint = std::min(static_cast<unsigned int>(backgrounds.size()), nBackgrounds);
            auto combined = signals;
            combined.insert(combined.end(), backgrounds.begin(), std::next(backgrounds.begin(), nBackgroundsToPrint));

            // Determine which columns we can print
            const auto purityDenominator = static_cast<float>((sampleType == DataBNB) ? totalCountDataBNB : totalCountMC);
            const auto shouldPrintPurity = purityDenominator > std::numeric_limits<float>::epsilon();
            const auto shouldPrintRate = hasAllTag;
            const auto shouldPrintRemoved = (iTag > 0);

            // Print the header
            std::cout << std::setw(alignWidthName) << std::left << "Topology" << " : " << std::setw(alignWidthValue) << "Count";
            if (shouldPrintPurity)
                std::cout << " : " << std::setw(alignWidthValue) << "Purity";

            if (shouldPrintRate)
                std::cout << " : " << std::setw(alignWidthValue) << "Rate";
            
            if (shouldPrintRemoved)
                std::cout << " : " << std::setw(alignWidthValue) << "Removed";

            std::cout << std::endl;

            // Print the values
            for (const auto &entry : combined)
            {
                std::cout << std::setw(alignWidthName) << std::left << entry.first << " : " << std::setw(alignWidthValue) << entry.second;
                if (shouldPrintPurity)
                {
                    const auto purity = static_cast<float>(entry.second) / purityDenominator;
                    std::cout << " : " << std::setw(alignWidthValue) << purity;
                }
            
                if (shouldPrintRate)
                {
                    bool printed = false;
                    const auto allIter2 = allIter->second.find(sampleType);

                    if (allIter2 != allIter->second.end())
                    {
                        const auto allIter3 = allIter2->second.find(entry.first);
                        
                        if (allIter3 != allIter2->second.end())
                        {
                            const auto allCount = allIter3->second;
                            const auto denom = static_cast<float>(allCount);
                            if (denom > std::numeric_limits<float>::epsilon())
                            {
                                const auto efficiency = static_cast<float>(entry.second) / denom;
                                std::cout << " : " << std::setw(alignWidthValue) << efficiency;
                                printed = true;
                            }
                        }
                    }

                    if (!printed)
                    {
                        std::cout << " : " << std::setw(alignWidthValue) << "?";
                    }
                }

                if (shouldPrintRemoved)
                {
                    bool printed = false;
                    const auto iterLast = counterMap.find(reverseTags.at(iTag - 1));

                    if (iterLast != counterMap.end())
                    {
                        const auto iterLast2 = iterLast->second.find(sampleType);
                        if (iterLast2 != iterLast->second.end())
                        {
                            const auto iterLast3 = iterLast2->second.find(entry.first);

                            if (iterLast3 != iterLast2->second.end())
                            {
                                const auto lastCount = iterLast3->second;
                                const auto denom = static_cast<float>(lastCount);

                                if (denom > std::numeric_limits<float>::epsilon())
                                {
                                    const auto fracRemoved = static_cast<float>(entry.second) / denom;
                                    std::cout << " : " << std::setw(alignWidthValue) << fracRemoved;
                                    printed = true;
                                }
                            }
                        }
                    }
                    
                    if (!printed)
                    {
                        std::cout << " : " << std::setw(alignWidthValue) << "?";
                    }
                }

                std::cout << std::endl;
            }
        }
    }
}

} // namespace ubcc1pi

#endif
