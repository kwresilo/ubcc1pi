/**
 *  @file  ubcc1pi/Helpers/CollectionHelper.h
 *
 *  @brief The header file for the collection helper class
 */

#ifndef UBCC1PI_HELPERS_COLLECTION_HELPER
#define UBCC1PI_HELPERS_COLLECTION_HELPER

#include "art/Framework/Principal/Event.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/FindManyP.h"

#include <vector>
#include <unordered_map>

// ATTN could just forward declare these
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/SpacePoint.h"

#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"

namespace ubcc1pi
{

// Typedefs scoped to the namespace
template <typename T, typename D>
using ObjectData = std::pair<art::Ptr<T>, D>;

template <typename T>
using Collection = std::vector< art::Ptr<T> >;

template <typename T, typename D>
using CollectionData = std::vector<ObjectData<T, D> >;

template <typename L, typename R>
using Association = std::unordered_map< art::Ptr<L>, Collection<R> >;

template <typename L, typename R, typename D>
using AssociationData = std::unordered_map< art::Ptr<L>, CollectionData<R, D> >;

typedef Collection<simb::MCParticle> MCParticleVector;   ///< A collection of MCParticles
typedef Collection<recob::PFParticle> PFParticleVector;  ///< A collection of PFParticles
typedef Collection<recob::Hit> HitVector;                ///< A collection of Hits
typedef Collection<recob::Slice> SliceVector;            ///< A collection of Slices
typedef Collection<recob::Track> TrackVector;            ///< A colleciton of Tracks
typedef Collection<anab::Calorimetry> CalorimetryVector; ///< A collection of Calorimetry objects

typedef Association<recob::Hit, recob::PFParticle> HitsToPFParticles;                  ///< Association from Hit to PFParticle
typedef AssociationData<recob::Hit, simb::MCParticle, float> HitsToMCParticleWeights;  ///< Association from Hit to MCParticle along with the charge contribution weight
typedef AssociationData<simb::MCParticle, recob::Hit, float> MCParticlesToHitWeights;  ///< Association from MCParticle to Hits along with the charge contribution weight
typedef Association<recob::Slice, recob::Hit> SlicesToHits;                            ///< Association from Slice to Hit
typedef Association<recob::Slice, recob::PFParticle> SlicesToPFParticles;              ///< Association from Slice to PFParticle
typedef Association<recob::PFParticle, recob::Hit> PFParticleToHits;                   ///< Association from PFParticle to Hit
typedef Association<recob::PFParticle, recob::Track> PFParticleToTracks;               ///< Association from PFParticle to Tracks
typedef Association<recob::PFParticle, recob::PFParticle> PFParticleToPFParticles;     ///< Association from PFParticle to PFParticles
                
typedef std::unordered_map<art::Ptr<recob::Hit>, bool> HitsToBool;                     ///< Mapping from Hit to boolean
typedef std::unordered_map<art::Ptr<recob::Slice>, bool> SlicesToBool;                 ///< Mapping from Slice to boolean


/**
 *  @brief  The collection helper class
 */
class CollectionHelper
{
    public:
        /**
         *  @brief  Get a collection of objects in the desired format from the event
         *
         *  @param  event the art event
         *  @param  label the label of the collection producer
         */
        template <typename T>
        static Collection<T> GetCollection(const art::Event &event, const art::InputTag &label);
        
        /**
         *  @brief  Get an association between objects in the desired format from the event
         *
         *  @param  event the art event
         *  @param  collectionLabel the label of the producer which made the collection of type L
         *  @param  associationLabel the label of the association producer
         */
        template <typename L, typename R>
        static Association<L, R> GetAssociation(const art::Event &event, const art::InputTag &collectionLabel, const art::InputTag &associationLabel);
        
        /**
         *  @brief  Get an association between objects in the desired format from the event
         *
         *  @param  event the art event
         *  @param  label the label of the association and collection producer
         */
        template <typename L, typename R>
        static Association<L, R> GetAssociation(const art::Event &event, const art::InputTag &label);
        
        /**
         *  @brief  Get the reversed association (R -> L) from in input association (L -> R)
         *
         *  @param  forwardAssociation the forward L -> R association
         *
         *  @return the reversed association R -> L
         */
        template <typename L, typename R>
        static Association<R, L> GetReversedAssociation(const Association<L, R> &forwardAssociation);
        
        /**
         *  @brief  Get the objects associated to a given input object
         *
         *  @param  objectL the input object
         *  @param  association the association between objects of type L -> R
         *
         *  @return the objects of type R associated to the input object
         */
        template <typename L, typename R>
        static Collection<R> GetManyAssociated(const art::Ptr<L> &objectL, const Association<L, R> &association);
        
        /**
         *  @brief  Check if there is any association availble for the input object
         *
         *  @param  objectL the input object
         *  @param  association the association between objects of type L -> R
         *
         *  @return if the association is available
         */
        template <typename L, typename R>
        static bool HasAssociated(const art::Ptr<L> &objectL, const Association<L, R> &association);
        
        /**
         *  @brief  Get the single object associated to a given input object
         *
         *  @param  objectL the input object
         *  @param  association the association between objects of type L -> R
         *
         *  @return the object of type R associated to the input object
         */
        template <typename L, typename R>
        static art::Ptr<R> GetSingleAssociated(const art::Ptr<L> &objectL, const Association<L, R> &association);
        
        /**
         *  @brief  Get an association between objects in the desired format from the event with data
         *
         *  @param  event the art event
         *  @param  collectionLabel the label of the producer which made the collection of type L
         *  @param  associationLabel the label of the association producer
         */
        template <typename L, typename R, typename D>
        static AssociationData<L, R, D> GetAssociationWithData(const art::Event &event, const art::InputTag &collectionLabel, const art::InputTag &associationLabel);
        
        /**
         *  @brief  Get an association between objects in the desired format from the event with data
         *
         *  @param  event the art event
         *  @param  label the label of the association and collection producer
         */
        template <typename L, typename R, typename D>
        static AssociationData<L, R, D> GetAssociationWithData(const art::Event &event, const art::InputTag &label);
        
        /**
         *  @brief  Get the reversed association (R -> L) from in input association (L -> R) with data
         *
         *  @param  forwardAssociation the forward L -> R association
         *
         *  @return the reversed association R -> L
         */
        template <typename L, typename R, typename D>
        static AssociationData<R, L, D> GetReversedAssociation(const AssociationData<L, R, D> &forwardAssociation);
        
        /**
         *  @brief  Get the objects associated to a given input object
         *
         *  @param  objectL the input object
         *  @param  association the association between objects of type L -> R
         *
         *  @return the objects of type R associated to the input object
         */
        template <typename L, typename R, typename D>
        static Collection<R> GetManyAssociated(const art::Ptr<L> &objectL, const AssociationData<L, R, D> &association);
        
        /**
         *  @brief  Get the single object associated to a given input object
         *
         *  @param  objectL the input object
         *  @param  association the association between objects of type L -> R
         *
         *  @return the object of type R associated to the input object
         */
        template <typename L, typename R, typename D>
        static art::Ptr<R> GetSingleAssociated(const art::Ptr<L> &objectL, const AssociationData<L, R, D> &association);

        /**
         *  @brief  Get the objects associated to a given input object with data
         *
         *  @param  objectL the input object
         *  @param  association the association between objects of type L -> R
         *
         *  @return the objects of type R associated to the input object
         */
        template <typename L, typename R, typename D>
        static CollectionData<R, D> GetManyAssociatedWithData(const art::Ptr<L> &objectL, const AssociationData<L, R, D> &association);
        
        /**
         *  @brief  Get the single object associated to a given input object with data
         *
         *  @param  objectL the input object
         *  @param  association the association between objects of type L -> R
         *
         *  @return the object of type R associated to the input object
         */
        template <typename L, typename R, typename D>
        static ObjectData<R, D> GetSingleAssociatedWithData(const art::Ptr<L> &objectL, const AssociationData<L, R, D> &association);

        /**
         *  @brief  Get the objects that are in both collection A and B
         *
         *  @param  a the first collection A
         *  @param  b the second collection B
         */
        template <typename T>
        static Collection<T> GetIntersection(const Collection<T> &a, const Collection<T> &b);

        /**
         *  @brief  Get an association between object of type A and C, via an intermediate collection B
         *
         *  @param  event the art event
         *  @param  labelA the label of the A producer
         *  @param  labelB the label of the B producer
         *  @param  labelC the label of the B to C assoication producer
         */
        template <typename A, typename B, typename C>
        static Association<A, C> GetAssociationViaCollection(const art::Event &event, const art::InputTag &labelA, const art::InputTag &labelB, const art::InputTag &labelC);
};

// -----------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
inline Collection<T> CollectionHelper::GetCollection(const art::Event &event, const art::InputTag &label)
{
    Collection<T> outputCollection;

    const auto handle = event.getValidHandle< std::vector<T> >(label);
    for (size_t i = 0; i < handle->size(); ++i)
    {
        outputCollection.emplace_back(handle, i);    
    }

    return outputCollection;
}

// -----------------------------------------------------------------------------------------------------------------------------------------
        
template <typename L, typename R>
inline Association<L, R> CollectionHelper::GetAssociation(const art::Event &event, const art::InputTag &label)
{
    return CollectionHelper::GetAssociation<L, R>(event, label, label);  
}

// -----------------------------------------------------------------------------------------------------------------------------------------
        
template <typename L, typename R>
inline Association<L, R> CollectionHelper::GetAssociation(const art::Event &event, const art::InputTag &collectionLabel, const art::InputTag &associationLabel)
{
    Association<L, R> outputAssociation;

    const auto handle = event.getValidHandle< std::vector<L> >(collectionLabel);
    const art::FindManyP<R> assoc(handle, event, associationLabel);

    for (unsigned int i = 0; i < handle->size(); ++i)
    {
        const art::Ptr<L> objectL(handle, i);

        for (const auto &objectR : assoc.at(objectL.key()))
            outputAssociation[objectL].push_back(objectR);
    }

    return outputAssociation;
}

// -----------------------------------------------------------------------------------------------------------------------------------------
        
template <typename L, typename R>
inline Association<R, L> CollectionHelper::GetReversedAssociation(const Association<L, R> &forwardAssociation)
{
    Association<R, L> reverseAssociation;

    // ATTN does this mess up reproducibility?
    for (const auto &entry : forwardAssociation)
    {
        for (const auto &object : entry.second)
            reverseAssociation[object].push_back(entry.first);
    }

    return reverseAssociation;
}

// -----------------------------------------------------------------------------------------------------------------------------------------
        
template <typename L, typename R>
inline Collection<R> CollectionHelper::GetManyAssociated(const art::Ptr<L> &objectL, const Association<L, R> &association)
{
    const auto iter = association.find(objectL);

    // ATTN previously threw here. Changed, but need to test if this impacts any existing modules
    if (iter == association.end())
    {
        return Collection<R>();
        //throw cet::exception("CollectionHelper::GetManyAssociated") << " - No association entry found for the input object." << std::endl;
    }

    return iter->second;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

template <typename L, typename R>
inline bool CollectionHelper::HasAssociated(const art::Ptr<L> &objectL, const Association<L, R> &association)
{
    const auto iter = association.find(objectL);
    return (iter != association.end());
}

// -----------------------------------------------------------------------------------------------------------------------------------------
        
template <typename L, typename R>
inline art::Ptr<R> CollectionHelper::GetSingleAssociated(const art::Ptr<L> &objectL, const Association<L, R> &association)
{
    const auto objects = CollectionHelper::GetManyAssociated(objectL, association);

    if (objects.size() != 1)
        throw cet::exception("CollectionHelper::GetSingleAssociated") << " - Found " << objects.size() << " objects associated to the input object. Expected 1." << std::endl;

    return objects.front();
}

// -----------------------------------------------------------------------------------------------------------------------------------------
        
template <typename L, typename R, typename D>
inline AssociationData<L, R, D> CollectionHelper::GetAssociationWithData(const art::Event &event, const art::InputTag &label)
{
    return CollectionHelper::GetAssociationWithData<L, R, D>(event, label, label);
}

// -----------------------------------------------------------------------------------------------------------------------------------------
        
template <typename L, typename R, typename D>
inline AssociationData<L, R, D> CollectionHelper::GetAssociationWithData(const art::Event &event, const art::InputTag &collectionLabel, const art::InputTag &associationLabel)
{
    AssociationData<L, R, D> outputAssociation;

    const auto handle = event.getValidHandle< std::vector<L> >(collectionLabel);
    const art::FindManyP<R, D> assoc(handle, event, associationLabel);

    for (unsigned int i = 0; i < handle->size(); ++i)
    {
        const art::Ptr<L> objectL(handle, i);
        const auto key = objectL.key();
        const auto objects = assoc.at(key);
        const auto data = assoc.data(key);

        if (objects.size() != data.size())
            throw cet::exception("CollectionHelper::GetAssociationWithData") << " - Number of metadata doesn't match number of associated objects." << std::endl;

        for (unsigned int j = 0; j < objects.size(); ++j)
            outputAssociation[objectL].emplace_back(objects.at(j), *data.at(j));
    }

    return outputAssociation;
}

// -----------------------------------------------------------------------------------------------------------------------------------------
        
template <typename L, typename R, typename D>
inline AssociationData<R, L, D> CollectionHelper::GetReversedAssociation(const AssociationData<L, R, D> &forwardAssociation)
{
    AssociationData<R, L, D> reverseAssociation;

    for (const auto &entry : forwardAssociation)
    {
        for (const auto &object : entry.second)
            reverseAssociation[object.first].emplace_back(entry.first, object.second);
    }

    return reverseAssociation;
}

// -----------------------------------------------------------------------------------------------------------------------------------------
        
template <typename L, typename R, typename D>
inline Collection<R> CollectionHelper::GetManyAssociated(const art::Ptr<L> &objectL, const AssociationData<L, R, D> &association)
{
    const auto iter = association.find(objectL);

    if (iter == association.end())
        throw cet::exception("CollectionHelper::GetManyAssociated") << " - No association entry found for the input object." << std::endl;

    Collection<R> outputCollection;
    for (const auto &entry : iter->second)
        outputCollection.push_back(entry.first);

    return outputCollection;
}

// -----------------------------------------------------------------------------------------------------------------------------------------
        
template <typename L, typename R, typename D>
inline art::Ptr<R> CollectionHelper::GetSingleAssociated(const art::Ptr<L> &objectL, const AssociationData<L, R, D> &association)
{
    const auto objects = CollectionHelper::GetManyAssociated(objectL, association);

    if (objects.size() != 1)
        throw cet::exception("CollectionHelper::GetSingleAssociated") << " - Found " << objects.size() << " objects associated to the input object. Expected 1." << std::endl;

    return objects.front().first;
}

// -----------------------------------------------------------------------------------------------------------------------------------------
        
template <typename L, typename R, typename D>
inline CollectionData<R, D> CollectionHelper::GetManyAssociatedWithData(const art::Ptr<L> &objectL, const AssociationData<L, R, D> &association)
{
    const auto iter = association.find(objectL);

    // ATTN previously threw here. Changed, but need to test if this impacts any existing modules
    if (iter == association.end())
    {
        return CollectionData<R, D>();
        //throw cet::exception("CollectionHelper::GetManyAssociated") << " - No association entry found for the input object." << std::endl;
    }

    return iter->second;
}

// -----------------------------------------------------------------------------------------------------------------------------------------
        
template <typename L, typename R, typename D>
inline ObjectData<R, D> CollectionHelper::GetSingleAssociatedWithData(const art::Ptr<L> &objectL, const AssociationData<L, R, D> &association)
{
    const auto objects = CollectionHelper::GetManyAssociatedWithData(objectL, association);

    if (objects.size() != 1)
        throw cet::exception("CollectionHelper::GetSingleAssociated") << " - Found " << objects.size() << " objects associated to the input object. Expected 1." << std::endl;

    return objects.front();
}

// -----------------------------------------------------------------------------------------------------------------------------------------
        
template <typename T>
inline Collection<T> CollectionHelper::GetIntersection(const Collection<T> &a, const Collection<T> &b)
{
    Collection<T> intersection;

    for (const auto &objectA : a)
    {
        if (std::find(b.begin(), b.end(), objectA) != b.end())
            intersection.push_back(objectA);
    }
    
    return intersection;
}

// -----------------------------------------------------------------------------------------------------------------------------------------
        
template <typename A, typename B, typename C>
inline Association<A, C> CollectionHelper::GetAssociationViaCollection(const art::Event &event, const art::InputTag &labelA, const art::InputTag &labelB, const art::InputTag &labelC)
{
    const auto collectionA = CollectionHelper::GetCollection<A>(event, labelA);
    const auto assocAB = CollectionHelper::GetAssociation<A, B>(event, labelA, labelB);
    const auto assocBC = CollectionHelper::GetAssociation<B, C>(event, labelB, labelC);

    Association<A, C> outputAssociation;

    for (const auto &objectA : collectionA)
    {
        try
        {
            const auto collectionB = CollectionHelper::GetManyAssociated(objectA, assocAB);
            for (const auto &objectB : collectionB)
            {
                const auto collectionC = CollectionHelper::GetManyAssociated(objectB, assocBC);

                auto &outputCollection = outputAssociation[objectA];
                outputCollection.insert(outputCollection.end(), collectionC.begin(), collectionC.end());
            }
        }
        catch (const cet::exception &)
        {
        }
    }

    return outputAssociation;
}


} // namespace ubcc1pi

#endif
