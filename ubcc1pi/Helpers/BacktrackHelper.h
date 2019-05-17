/**
 *  @file  ubcc1pi/Helpers/BacktrackHelper.h
 *
 *  @brief The header file for the back tracking helper class
 */

#ifndef UBCC1PI_HELPERS_BACKTRACK_HELPER
#define UBCC1PI_HELPERS_BACKTRACK_HELPER

#include "ubcc1pi/Helpers/CollectionHelper.h"

#include <unordered_map>
#include <tuple>

namespace ubcc1pi
{

/**
 *  @brief  The backtrack helper class
 *
 *          Helper functions for matching PFParticles to MCParticles
 */
class BacktrackHelper
{
    public:
        /**
         *  @brief  Class holding backtracking data
         */
        class BacktrackerData
        {
            public:
                typedef std::unordered_map<art::Ptr<simb::MCParticle>, float> MCParticleToFloatMap;
                typedef std::unordered_map<art::Ptr<recob::PFParticle>, float> PFParticleToFloatMap;
                typedef AssociationData<recob::PFParticle, simb::MCParticle, float> MatchMap;

                /**
                 *  @brief  Constructor
                 *
                 *  @param  pfParticles the input pfParticles
                 *  @param  mcParticles the input mcParticles
                 *  @param  hitsToPfps the input mapping from hits to PFParticles
                 *  @param  hitsToMcps the input mapping from hits to MCParticles
                 */
                BacktrackerData(const PFParticleVector &pfParticles, const MCParticleVector &mcParticles, const HitsToPFParticles &hitsToPfps, const HitsToMCParticleWeights &hitsToMcps);

                /**
                 *  @brief  Get the MCParticles
                 */
                MCParticleVector GetMCParticles() const;
                
                /**
                 *  @brief  Get the PFParticles
                 */
                PFParticleVector GetPFParticles() const;

                /**
                 *  @brief  Get the number of hits associated with a given PFParticle
                 *
                 *  @param  pfParticle the pfParticle
                 */
                unsigned int GetNHits(const art::Ptr<recob::PFParticle> &pfParticle) const;
                
                /**
                 *  @brief  Get the weight associated with a given PFParticle (= the number of hits)
                 *
                 *  @param  pfParticle the pfParticle
                 */
                float GetWeight(const art::Ptr<recob::PFParticle> &pfParticle) const;
                
                /**
                 *  @brief  Get the weight associated with a given MCParticle (= the weighted number of hits)
                 *
                 *  @param  mcParticle the mcParticle
                 */
                float GetWeight(const art::Ptr<simb::MCParticle> &mcParticle) const;

                /**
                 *  @brief  Get all hits that are associated the MCParticle with any weight
                 *
                 *  @param  mcParticle the mcParticle
                 */
                HitVector GetHits(const art::Ptr<simb::MCParticle> &mcParticle) const;
                
                /**
                 *  @brief  Get all hits that are associated the PFParticle
                 *
                 *  @param  pfParticle the pfParticle
                 */
                HitVector GetHits(const art::Ptr<recob::PFParticle> &pfParticle) const;
                
                /**
                 *  @brief  Get the weight of a given PFParticle-MCParticle pair
                 *
                 *  @param  pfParticle the pfParticle
                 *  @param  mcParticle the mcParticle
                 */
                float GetMatchWeight(const art::Ptr<recob::PFParticle> &pfParticle, const art::Ptr<simb::MCParticle> &mcParticle) const;

                /**
                 *  @brief  Get the purity of a given PFParticle-MCParticle pair = shared weight / PFParticle weight
                 *          Purity is the fraction of the PFParticle that represents the MCParticle
                 *
                 *  @param  pfParticle the pfParticle
                 *  @param  mcParticle the mcParticle
                 */
                float GetMatchPurity(const art::Ptr<recob::PFParticle> &pfParticle, const art::Ptr<simb::MCParticle> &mcParticle) const;
                
                /**
                 *  @brief  Get the completeness a given PFParticle-MCParticle pair = shared weight / MCParticle weight
                 *          Completeness is the fraction of the MCParticle that is represented by the PFParticle
                 *
                 *  @param  pfParticle the pfParticle
                 *  @param  mcParticle the mcParticle
                 */
                float GetMatchCompleteness(const art::Ptr<recob::PFParticle> &pfParticle, const art::Ptr<simb::MCParticle> &mcParticle) const;

                /**
                 *  @brief  Get the MCParticle with the highest weighted match to the given PFParticle
                 *
                 *  @param  pfParticle the pfParticle
                 */
                art::Ptr<simb::MCParticle> GetBestMatchedMCParticle(const art::Ptr<recob::PFParticle> &pfParticle) const;
                
                /**
                 *  @brief  Get the PFParticles which have the given MCParticle as their strongest match
                 *
                 *  @param  mcParticle the mcParticle
                 */
                PFParticleVector GetBestMatchedPFParticles(const art::Ptr<simb::MCParticle> &mcParticle) const;

                /**
                 *  @brief  Get the weight associated with the input Hit and MCParticle
                 *          weight = fraction of the energy of the hit that was provided by the supplied MCParticle
                 *
                 *  @param  hit the input hit
                 *  @param  mcParticle the input MCParticle
                 */
                float GetWeight(const art::Ptr<recob::Hit> &hit, const art::Ptr<simb::MCParticle> &mcParticle) const;

                /**
                 *  @brief  Get the summed weight associated with the input collection of Hits and MCParticle
                 *
                 *  @param  hits the input hits
                 *  @param  mcParticle the input MCParticle
                 */
                float GetWeight(const HitVector &hits, const art::Ptr<simb::MCParticle> &mcParticle) const;

            private:

                // TODO doxygen comments
                void CheckConsistency(const PFParticleVector &pfParticles, const HitsToPFParticles &hitsToPfps) const;
                void CheckConsistency(const MCParticleVector &mcParticles, const HitsToMCParticleWeights &hitsToMcps) const;
                HitVector CollectHits(const HitsToPFParticles &hitsToPfps, const HitsToMCParticleWeights &hitsToMcps) const;
                bool CollectPFParticle(const art::Ptr<recob::Hit> &hit, const HitsToPFParticles &hitsToPfps, art::Ptr<recob::PFParticle> &outputParticle) const;
                bool CollectMCParticleWeights(const art::Ptr<recob::Hit> &hit, const HitsToMCParticleWeights &hitsToMcps, CollectionData<simb::MCParticle, float> &outputParticleWeights) const;

                MCParticleVector        m_mcParticles;         ///< The MCParticles
                PFParticleVector        m_pfParticles;         ///< The PFParticles
                MCParticleToFloatMap    m_mcParticleWeightMap; ///< The mapping from MCParticle to the total weight (= weighted number of hits)
                PFParticleToFloatMap    m_pfParticleWeightMap; ///< The mapping from PFParticle to the total weight (= number of hits)
                MatchMap                m_matchMap;            ///< The matching between PFParticles and MCParticles along with the shared weight
                HitsToMCParticleWeights m_hitsToMcps;          ///< The mapping from hits to MCParticle along with the associated weights 
                HitsToPFParticles       m_hitsToPfps;          ///< The mapping from hits to PFParticles
        };

        /**
         *  @brief  Class holding the backtracking details of the slices
         */
        class SliceMetadata
        {
            public:
                /** 
                 *  @brief  Constructor
                 *
                 *  @param  slices the input slices
                 *  @param  sliceToIsSelectedAsNu the input map from a slice to whether it's been selected as the neutrino - There should only ever be 1 or zero slices selected
                 *  @param  slicesToHits the input mapping from hits to slices
                 *  @param  hitsToIsNuInduced the input mapping from a hit to whether it's neutrino induced or not
                 */
                SliceMetadata(const SliceVector &slices, const SlicesToBool &sliceToIsSelectedAsNu, const SlicesToHits &slicesToHits, const HitsToBool &hitsToIsNuInduced);

                /** 
                 *  @brief  Get the slices
                 */
                SliceVector GetSlices() const;

                /** 
                 *  @brief  Get the purity of a given slice
                 *          purity = fraction of the hits in the slice that are neutrino induced
                 *
                 *  @param  slice the input slice
                 */
                float GetPurity(const art::Ptr<recob::Slice> &slice) const;
                
                /** 
                 *  @brief  Get the completeness of a given slice
                 *          completeness = fraction of the neutrino induced hits in the event that are in the slice
                 *
                 *  @param  slice the input slice
                 */
                float GetCompleteness(const art::Ptr<recob::Slice> &slice) const;

                /** 
                 *  @brief  Get the number of hits in the slice
                 *
                 *  @param  slice the input slice
                 */
                unsigned int GetNumberOfHits(const art::Ptr<recob::Slice> &slice) const;
                
                /** 
                 *  @brief  Get the hits in the given slice
                 *
                 *  @param  slice the input slice
                 */
                HitVector GetHits(const art::Ptr<recob::Slice> &slice) const;
                
                /** 
                 *  @brief  Get the number of neutrino induced hits in the slice
                 *
                 *  @param  slice the input slice
                 */
                unsigned int GetNumberOfNuInducedHits(const art::Ptr<recob::Slice> &slice) const;
                
                /** 
                 *  @brief  Get the total number of neutrino induced hits
                 */
                unsigned int GetTotalNumberOfNuInducedHits() const;

                /** 
                 *  @brief  Get the slice with the largest completeness
                 */
                art::Ptr<recob::Slice> GetMostCompleteSlice() const;
                
                /** 
                 *  @brief  Get the slices that have been selected a neutrinos
                 */
                SliceVector GetSelectedNeutrinoSlices() const;
                
                /** 
                 *  @brief  Determine if the most complete slice has been selected as a neutrino
                 */
                bool IsMostCompleteSliceSelected() const;

            private:

                /** 
                 *  @brief  Get all hits in the event through the slices
                 */
                HitVector GetAllHits() const;

                /** 
                 *  @brief  Count the number of hits in the input vector that are neutrino induced
                 *
                 *  @param  hits the input hits
                 */
                unsigned int CountNuInducedHits(const HitVector &hits) const;

                SliceVector  m_slices;
                SlicesToBool m_sliceToIsSelectedAsNu;
                SlicesToHits m_slicesToHits;
                HitsToBool   m_hitsToIsNuInduced;
        };

    /**
     *  @brief  Get the mapping from neutrino final state PFParticles to hits, folding in all downstream PFParticles
     *
     *  @param  event the art event
     *  @param  pfParticleLabel the PFParticle producer label
     *  @param  finalStates the input vector of final state PFParticles
     */
    static HitsToPFParticles GetHitToPFParticleMap(const art::Event &event, const art::InputTag &pfParticleLabel, const PFParticleVector &finalStates);
    
    /**
     *  @brief  Get the mapping from neutrino final state MCParticles to hits, folding in all downstream particles
     *
     *  @param  event the art event
     *  @param  mcParticleLabel the MCParticle producer label
     *  @param  backtrackerLabel the MCParticle to hit producer label
     *  @param  finalStates the input vector of final state MCParticles
     */
    static HitsToMCParticleWeights GetHitToMCParticleWeightMap(const art::Event &event, const art::InputTag &mcParticleLabel, const art::InputTag &backtrackerLabel, const MCParticleVector &finalStates);

    /**
     *  @brief  Get the mapping from all hits, to whether the hit is induced (at least in part) by a neutrino induced MCParticle
     *
     *  @param  event the art event
     *  @param  hitLabel the Hit producer label
     *  @param  mcParticleLabel the MCParticle producer label
     *  @param  backtrackerLabel the MCParticle to hit producer label
     *  @param  nuMCParticles the input vector of neutrino induced MCParticles
     */
    static HitsToBool GetHitsToIsNuInducedMap(const art::Event &event, const art::InputTag &hitLabel, const art::InputTag &mcParticleLabel, const art::InputTag &backtrackerLabel, const MCParticleVector &nuMCParticles);
};

} // namespace ubcc1pi

#endif
