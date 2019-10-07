/**
 *  @file  ubcc1pi/Producers/PFParticleHierarchy.h
 *
 *  @brief The header file for the PFParticle hierarchy producer
 */

#ifndef UBCC1PI_PRODUCERS_PFPARTICLE_HIERARCHY
#define UBCC1PI_PRODUCERS_PFPARTICLE_HIERARCHY

#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDProducer.h"

#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Table.h"
#include "canvas/Utilities/InputTag.h"

#include "ubcc1pi/Helpers/CollectionHelper.h"

namespace ubcc1pi
{

/**
 *  @brief  The pfparticle hierarchy producer class
 */
class PFParticleHierarchy : public art::EDProducer
{
    public:
        /**
         *  @brief  The configuration structure
         */
        struct Config
        {
            fhicl::Atom<art::InputTag> PFParticleLabel
            {
                fhicl::Name("PFParticleLabel"),
                fhicl::Comment("The label for the PFParticle producer (Pandora)")
            };
            
            fhicl::Atom<art::InputTag> SpacePointLabel
            {
                fhicl::Name("SpacePointLabel"),
                fhicl::Comment("The label for the PFParticle -> SpacePoint producer (Pandora)")
            };
            
            fhicl::Atom<art::InputTag> VertexLabel
            {
                fhicl::Name("VertexLabel"),
                fhicl::Comment("The label for the PFParticle -> Vertex producer (Pandora)")
            };

            fhicl::Atom<float> DistanceThreshold
            {
                fhicl::Name("DistanceThreshold"),
                fhicl::Comment("The distance within which a PFParticle will automatically be associated to an existing match")
            };
        };
        
        /**
         *  @brief  Constructor
         *
         *  @param  config the set of input fhicl parameters
         */
        PFParticleHierarchy(const art::EDProducer::Table<Config> &config);

        /**
         *  @brief  Produce collections in the event record
         *
         *  @param  event the event to add to
         */
        void produce(art::Event &event);

    private:
        typedef std::unordered_map<art::Ptr<recob::PFParticle>, std::vector<TVector3> > PFParticlesToPoints;
        typedef AssociationData<recob::PFParticle, recob::PFParticle, float> PFParticlePairSeparationMap;
        typedef AssociationData<recob::PFParticle, recob::PFParticle, TVector3> PFParticlePairMatchPositionMap;
        typedef std::unordered_map<art::Ptr<recob::PFParticle>, float> PFParticleToFloatMap;
        typedef std::unordered_map<art::Ptr<recob::PFParticle>, size_t> PFParticleToIndexMap;
        typedef std::unordered_map<size_t, size_t> IdToIdMap;

        /**
         *  @brief  The match point class
         */
        class MatchPoint
        {
            public:
                /**
                 *  @brief  Constructor
                 *
                 *  @param  position the position of the matching point
                 *  @param  parent the parent PFParticle to which this match point relates
                 */
                MatchPoint(const TVector3 &position, const art::Ptr<recob::PFParticle> &parent);

                /**
                 *  @brief  Calculate the minimum distance between the input PFParticles and the match point for future use
                 *
                 *  @param  orphans the input vector of current orphan PFParticles
                 *  @param  pfParticleToPointsMap the input mapping from PFParticles to their point positions
                 */
                void CachePFParticleDistances(const PFParticleVector &orphans, const PFParticlesToPoints &pfParticleToPointsMap);

                /**
                 *  @brief  Get the distance between an orphan PFParticle and this match point, using the cached distances
                 *
                 *  @param  orphan the orphan PFParticle
                 *
                 *  @return the distance
                 */
                float GetDistance(const art::Ptr<recob::PFParticle> &orphan) const;

                /**
                 *  @brief  Add a parent daughter match at this point
                 *
                 *  @param  daughter the daughter PFParticle
                 */
                void AddDaughter(const art::Ptr<recob::PFParticle> &daughter);

                /**
                 *  @brief  Get the parent PFParticle
                 *
                 *  @return parent PFParticle
                 */
                art::Ptr<recob::PFParticle> GetParent() const;

                /**
                 *  @brief  Get the daughter PFParticles
                 *
                 *  @return the daughter PFParticles
                 */
                PFParticleVector GetDaughters() const;

                /**
                 *  @brief  Get the position
                 *
                 *  @return the position
                 */
                TVector3 GetPosition() const;
                
            private:
                TVector3                    m_position;            ///< The position of the match point
                art::Ptr<recob::PFParticle> m_parent;              ///< The parent PFParticle
                PFParticleVector            m_daughters;           ///< The daugher PFParticles
                PFParticleToFloatMap        m_pfParticleDistances; ///< The cached distances of PFParticles to this point
        };
        
        typedef std::vector<std::shared_ptr<MatchPoint> > MatchPointVector;
        typedef std::map<art::Ptr<recob::PFParticle>, std::shared_ptr<MatchPoint> > PFParticleToMatchPointMap;

        /**
         *  @brief  Collect the neutrino induced PFParticles from the event
         *  
         *  @param  event the art event
         *
         *  @return vector of neutrino induced PFParticles
         */
        PFParticleVector GetNeutrinoInducedPFParticles(const art::Event &event) const;

        /**
         *  @brief  Get the mapping from the input PFParticles to their space-charge corrected 3D spacepoints (as TVector3s)
         *
         *  @param  event the art event
         *  @param  pfParticles the input PFParticles
         *
         *  @return the output mapping from PFParticles to spacepoints
         */
        PFParticlesToPoints GetPFParticlesTo3DPointsMap(const art::Event &event, const PFParticleVector &pfParticles) const;

        /**
         *  @brief  Get the minimum separation and the position of minimum separation between each PFParticle pair
         *
         *  @param  pfParticles the input PFParticles
         *  @param  pfParticleToPointsMap the input mapping from PFParticles to spacepoints
         *  @param  pfParticlePairSeparationMap the output minimum separations between each PFParticle pair
         *  @param  pfParticlePairMatchPositionMap the output positions of minimum separation for each PFParticle pair
         */
        void GetPFParticleMatches(const PFParticleVector &pfParticles, const PFParticlesToPoints &pfParticleToPointsMap,
                                  PFParticlePairSeparationMap &pfParticlePairSeparationMap,
                                  PFParticlePairMatchPositionMap &pfParticlePairMatchPositionMap) const;

        /**
         *  @brief  Add a seed match point at the neutrino vertex
         *
         *  @param  pfParticleToPointsMap the input mapping from PFParticles to spacepoints
         *  @param  orphans the orphan PFParticles
         *  @param  matchPoints the match points
         */
        void SeedMatchPoints(const PFParticlesToPoints &pfParticleToPointsMap, PFParticleVector &orphans,
                             MatchPointVector &matchPoints) const;

        /**
         *  @brief  Get the closest PFParticle that isn't an orphan to the input orphan PFParticle
         *
         *  @param  orphan the orphan PFParticle with which to find the closest connected PFParticle
         *  @param  orphans the current list of orphans
         *  @param  pfParticlePairSeparationMap the minimum sepration between PFParticle pairs
         *  @param  minDist the output minimum distance to the closest PFParticle
         *
         *  @return the closest PFParticle
         */
        art::Ptr<recob::PFParticle> GetClosestPFParticle(const art::Ptr<recob::PFParticle> &orphan, 
                                                         const PFParticleVector &orphans,
                                                         const PFParticlePairSeparationMap &pfParticlePairSeparationMap,
                                                         float &minDist) const;

        /**
         *  @brief  Get the closest match point in the input vector to the input orphan PFParticle
         *
         *  @param  orphan the orphan PFParticle with which to find the closest match point
         *  @param  matchPoints the input vector of match points
         *
         *  @return the closest match point
         */
        std::shared_ptr<MatchPoint> GetClosestMatchPoint(const art::Ptr<recob::PFParticle> &orphan,
                                                         const MatchPointVector &matchPoints) const;

        /**
         *  @brief  Add the orphans to the input match points if they are within threshold
         *
         *  @param  orphans the input vector of orphan PFParticles
         *  @param  matchPoints the input match points
         */
        void AddOrphansToExistingMatchPoints(PFParticleVector &orphans, MatchPointVector &matchPoints) const;

        /**
         *  @brief  Make the best available match between a single orphan PFParticle and either an existing match point or a new one
         *
         *  @param  pfParticlePairSeparationMap the minimum separation distance between PFParticle pairs
         *  @param  pfParticlePairMatchPositionMap the position of minimum separation between PFParticle pairs
         *  @param  pfParticleToPointsMap the points for each PFParticle
         *  @param  orphans the input vector of orphan PFParticles
         *  @param  matchPoints the input match points
         */
        void MakeBestMatch(const PFParticlePairSeparationMap &pfParticlePairSeparationMap,
                           const PFParticlePairMatchPositionMap &pfParticlePairMatchPositionMap,
                           const PFParticlesToPoints &pfParticleToPointsMap,
                           PFParticleVector &orphans,
                           MatchPointVector &matchPoints) const;

        /**
         *  @brief  Get the mapping from PFParticle to it's parent
         *
         *  @param  pfParticles the input vector of pfparticles
         *  @param  matchPoints the input vector of match points
         *  @param  pfParticleToMatchPointMap the output mapping from PFParticles to their parent match points
         *  @param  parentMap the output mapping from PFParticle ID to parent ID
         *
         */
        void GetParentMap(const PFParticleVector &pfParticles, const MatchPointVector &matchPoints,
                          PFParticleToMatchPointMap &pfParticleToMatchPointMap, IdToIdMap &parentMap) const;

        /**
         *  @brief  Get the new PFParticles with the updated hierarchy information
         *
         *  @param  pfParticles the input original pfparticles
         *  @param  parentMap the input mapping from pfparticle ID to their parent IDs
         *  @param  outputPFParticles the output PFParticle vector to populate
         *  @param  pfParticleToNewIndicesMap the output mapping from old PFParticles to the new indices
         */
        void GetNewPFParticles(const PFParticleVector &pfParticles, const IdToIdMap &parentMap,
                               std::unique_ptr< std::vector<recob::PFParticle> > &outputPFParticles,
                               PFParticleToIndexMap &pfParticleToNewIndicesMap) const;

        /**
         *  @brief  Make new vertices for each PFParticle at their closest SCE corrected spacepoint to the parent match point
         *
         *  @param  pfParticles the input vector of PFParticles
         *  @param  pfParticleToNewIndicesMap the input mapping from old PFParticles to the index of the corresponding new PFParticle
         *  @param  pfParticleToMatchPointMap the input mapping from old PFParticles to their corresponding match point
         *  @param  pfParticleToPointsMap the input mapping from old PFParticles to their spacepoints
         *  @param  outputVertices output vector of vertices
         *  @param  outputParticlesToVertices the output association from PFParticles to vertices
         *  @param  event the art event
         */
        void GetNewVertices(const PFParticleVector &pfParticles, const PFParticleToIndexMap &pfParticleToNewIndicesMap,
                            const PFParticleToMatchPointMap &pfParticleToMatchPointMap, const PFParticlesToPoints &pfParticleToPointsMap,
                            std::unique_ptr< std::vector<recob::Vertex> > &outputVertices,
                            std::unique_ptr< art::Assns<recob::PFParticle, recob::Vertex> > &outputParticlesToVertices,
                            art::Event &event) const;

        /**
         *  @brief  Make new PFParticle metadata objets
         *
         *  @param  pfParticles the input vector of PFParticles
         *  @param  pfParticleToNewIndicesMap the input mapping from old PFParticles to the index of the corresponding new PFParticle
         *  @param  pfParticlePairSeparationMap the input mapping from PFParticle pairs to their closest separation
         *  @param  parentMap the input mapping from pfparticle ID to their new parent IDs
         *  @param  outputMetadata the output metadata vector
         *  @param  outputParticlesToMetadata the output PFParticle to metadata association
         *  @param  event the art event
         */
        void GetNewMetadata(const PFParticleVector &pfParticles, const PFParticleToIndexMap &pfParticleToNewIndicesMap,
                            const PFParticlePairSeparationMap &pfParticlePairSeparationMap, const IdToIdMap &parentMap,
                            std::unique_ptr< std::vector<larpandoraobj::PFParticleMetadata> > &outputMetadata,
                            std::unique_ptr< art::Assns<recob::PFParticle, larpandoraobj::PFParticleMetadata> > &outputParticlesToMetadata,
                            art::Event &event) const;
        
        art::EDProducer::Table<Config>  m_config;          ///< The FHiCL configuration options
};

} // namespace ubcc1pi

DEFINE_ART_MODULE(ubcc1pi::PFParticleHierarchy)

#endif
