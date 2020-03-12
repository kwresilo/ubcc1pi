/**
 *  @file  ubcc1pi_standalone/Interface/EventMembers.h
 *
 *  @brief The file defining the members of the event class
 *
 *  In this file we define all of the member variable that will go into the analysis files. The types and names of the variables are defined
 *  using pre-processor directives.
 *
 *  By defining the variables at the pre-processer level, all of the bookwork can be handled automatically. Each variable need only be
 *  defined in one place and will automatically:
 *      - Be added as a member variable of the ubcc1pi::Event class
 *      - Be bound to a branch if an analysis file is to be written or read using a consistent branch naming convention
 *      - Be assigned a default value and be checked such that you can't accidentally use a default variable at the analysis level 
 *      - Show up as a compiler error if used in C++ but not defined here
 *
 *  The Event class has three categories
 *      - Metadata
 *          - Defined by the UBCC1PI_MACRO_EVENT_METADATA_MEMBERS macro
 *          - Accessible in C++ via: event.metadata.variable()
 *          - Accessible in the root tree via: metadata_variable
 *      - Truth
 *          - Defined by the UBCC1PI_MACRO_EVENT_TRUTH_MEMBERS macro
 *          - Accessible in C++ via: event.truth.variable()
 *          - Accessible in the root tree via: truth_variable
 *          - Particle level information is defined by the UBCC1PI_MACRO_EVENT_TRUTH_PARTICLE_MEMBERS macro
 *              - Accessible in C++ via: event.truth.particles.at(i).variable()
 *              - Accessible in the root tree via: truth_particle_variable_vect
 *      - Reco
 *          - Defined by the UBCC1PI_MACRO_EVENT_RECO_MEMBERS macro
 *          - Accessible in C++ via: event.reco.variable()
 *          - Accessible in the root tree via: reco_variable
 *          - Particle level information is defined by the UBCC1PI_MACRO_EVENT_RECO_PARTICLE_MEMBERS macro
 *              - Accessible in C++ via: event.reco.particles.at(i).variable()
 *              - Accessible in the root tree via: reco_particle_variable_vect
 */


#ifndef UBCC1PI_STANDALONE_INTERFACE_EVENTMEMBERS
#define UBCC1PI_STANDALONE_INTERFACE_EVENTMEMBERS

#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <TVector3.h>

// The event metadata members
#define UBCC1PI_MACRO_EVENT_METADATA_MEMBERS(p, q, f)                                                                                      \
    f(p, q, int,         run)                                                                                                              \
    f(p, q, int,         subRun)                                                                                                           \
    f(p, q, int,         event)                                                                                                            \
    f(p, q, bool,        hasTruthInfo)                                                                                                     \

// The event truth information members
#define UBCC1PI_MACRO_EVENT_TRUTH_MEMBERS(p, q, f)                                                                                         \
    f(p, q, bool,                isCC)                                                                                                     \
    f(p, q, int,                 interactionMode)                                                                                          \
    f(p, q, std::string,         interactionString)                                                                                        \
    f(p, q, int,                 nuPdgCode)                                                                                                \
    f(p, q, float,               nuEnergy)                                                                                                 \
    f(p, q, TVector3,            nuVertex)                                                                                                 \
    f(p, q, int,                 nFinalStates)                                                                                             \
    f(p, q, std::vector<float>,  slicePurities)                                                                                            \
    f(p, q, std::vector<float>,  sliceCompletenesses)                                                                                      

// The event truth particle information members
#define UBCC1PI_MACRO_EVENT_TRUTH_PARTICLE_MEMBERS(p, q, f)                                                                                \
    f(p, q, int,                 pdgCode)                                                                                                  \
    f(p, q, float,               startX)                                                                                                   \
    f(p, q, float,               startY)                                                                                                   \
    f(p, q, float,               startZ)                                                                                                   \
    f(p, q, float,               endX)                                                                                                     \
    f(p, q, float,               endY)                                                                                                     \
    f(p, q, float,               endZ)                                                                                                     \
    f(p, q, float,               momentumX)                                                                                                \
    f(p, q, float,               momentumY)                                                                                                \
    f(p, q, float,               momentumZ)                                                                                                \
    f(p, q, float,               momentum)                                                                                                 \
    f(p, q, float,               energy)                                                                                                   \
    f(p, q, float,               mass)                                                                                                     \
    f(p, q, float,               hitWeightU)                                                                                               \
    f(p, q, float,               hitWeightV)                                                                                               \
    f(p, q, float,               hitWeightW)                                                                                               \
    f(p, q, int,                 nElasticScatters)                                                                                         \
    f(p, q, int,                 nInelasticScatters)                                                                                       \
    f(p, q, std::vector<float>,  scatterCosThetas)                                                                                         \
    f(p, q, std::vector<float>,  scatterMomentumFracsLost)                                                                                 \
    f(p, q, std::vector<bool>,   scatterIsElastic)                                                                                         \
    f(p, q, float,               scatteredMomentum)                                                                                        \
    f(p, q, float,               endMomentum)                                                                                              \
    f(p, q, bool,                isStopping)                                                                                               \
    f(p, q, int,                 endState)                                                                                                 \
    f(p, q, float,               endStateProductsHitWeightU)                                                                               \
    f(p, q, float,               endStateProductsHitWeightV)                                                                               \
    f(p, q, float,               endStateProductsHitWeightW)                                                                               

// The event reco information members
#define UBCC1PI_MACRO_EVENT_RECO_MEMBERS(p, q, f)                                                                                          \
    f(p, q, int,                 nSlices)                                                                                                  \
    f(p, q, bool,                hasSelectedSlice)                                                                                         \
    f(p, q, float,               selectedTopologicalScore)                                                                                 \
    f(p, q, std::vector<float>,  sliceTopologicalScores)                                                                                   \
    f(p, q, std::vector<bool>,   sliceIsSelectedAsNu)                                                                                      \
    f(p, q, bool,                hasNeutrino)                                                                                              \
    f(p, q, int,                 nuPdgCode)                                                                                                \
    f(p, q, TVector3,            nuVertex)                                                                                                 \
    f(p, q, int,                 nFinalStates)                                                                                                       

// The event reco particle information members
#define UBCC1PI_MACRO_EVENT_RECO_PARTICLE_MEMBERS(p, q, f)                                                                                 \
    f(p, q, int,                 pdgCode)                                                                                                  \
    f(p, q, int,                 nHitsU)                                                                                                   \
    f(p, q, int,                 nHitsV)                                                                                                   \
    f(p, q, int,                 nHitsW)                                                                                                   \
    f(p, q, int,                 nDaughters)                                                                                               \
    f(p, q, int,                 nDescendents)                                                                                             \
    f(p, q, int,                 nDescendentHitsU)                                                                                         \
    f(p, q, int,                 nDescendentHitsV)                                                                                         \
    f(p, q, int,                 nDescendentHitsW)                                                                                         \
    f(p, q, int,                 nHitsInLargestDescendent)                                                                                 \
    f(p, q, float,               trackScore)                                                                                               \
    f(p, q, float,               startX)                                                                                                   \
    f(p, q, float,               startY)                                                                                                   \
    f(p, q, float,               startZ)                                                                                                   \
    f(p, q, float,               endX)                                                                                                     \
    f(p, q, float,               endY)                                                                                                     \
    f(p, q, float,               endZ)                                                                                                     \
    f(p, q, float,               directionX)                                                                                               \
    f(p, q, float,               directionY)                                                                                               \
    f(p, q, float,               directionZ)                                                                                               \
    f(p, q, float,               yzAngle)                                                                                                  \
    f(p, q, float,               xyAngle)                                                                                                  \
    f(p, q, float,               xzAngle)                                                                                                  \
    f(p, q, float,               length)                                                                                                   \
    f(p, q, float,               range)                                                                                                    \
    f(p, q, float,               transverseVertexDist)                                                                                     \
    f(p, q, float,               longitudinalVertexDist)                                                                                   \
    f(p, q, float,               wiggliness)                                                                                               \
    f(p, q, int,                 nSpacePointsNearEnd)                                                                                      \
    f(p, q, float,               likelihoodForwardMuonU)                                                                                   \
    f(p, q, float,               likelihoodForwardMuonV)                                                                                   \
    f(p, q, float,               likelihoodForwardMuonW)                                                                                   \
    f(p, q, float,               likelihoodForwardMuon)                                                                                    \
    f(p, q, float,               likelihoodBackwardMuonU)                                                                                  \
    f(p, q, float,               likelihoodBackwardMuonV)                                                                                  \
    f(p, q, float,               likelihoodBackwardMuonW)                                                                                  \
    f(p, q, float,               likelihoodBackwardMuon)                                                                                   \
    f(p, q, float,               likelihoodForwardPionU)                                                                                   \
    f(p, q, float,               likelihoodForwardPionV)                                                                                   \
    f(p, q, float,               likelihoodForwardPionW)                                                                                   \
    f(p, q, float,               likelihoodForwardPion)                                                                                    \
    f(p, q, float,               likelihoodBackwardPionU)                                                                                  \
    f(p, q, float,               likelihoodBackwardPionV)                                                                                  \
    f(p, q, float,               likelihoodBackwardPionW)                                                                                  \
    f(p, q, float,               likelihoodBackwardPion)                                                                                   \
    f(p, q, float,               likelihoodForwardProtonU)                                                                                 \
    f(p, q, float,               likelihoodForwardProtonV)                                                                                 \
    f(p, q, float,               likelihoodForwardProtonW)                                                                                 \
    f(p, q, float,               likelihoodForwardProton)                                                                                  \
    f(p, q, float,               likelihoodBackwardProtonU)                                                                                \
    f(p, q, float,               likelihoodBackwardProtonV)                                                                                \
    f(p, q, float,               likelihoodBackwardProtonW)                                                                                \
    f(p, q, float,               likelihoodBackwardProton)                                                                                 \
    f(p, q, float,               likelihoodMIPU)                                                                                           \
    f(p, q, float,               likelihoodMIPV)                                                                                           \
    f(p, q, float,               likelihoodMIPW)                                                                                           \
    f(p, q, float,               likelihoodMIP)                                                                                            \
    f(p, q, float,               truncatedMeandEdxU)                                                                                       \
    f(p, q, float,               truncatedMeandEdxV)                                                                                       \
    f(p, q, float,               truncatedMeandEdxW)                                                                                       \
    f(p, q, float,               truncatedMeandEdx)                                                                                        \
    f(p, q, std::vector<float>,  truthMatchPurities)                                                                                       \
    f(p, q, std::vector<float>,  truthMatchCompletenesses)                                                                                 \
    f(p, q, bool,                hasMatchedMCParticle)                                                                                     \
    f(p, q, int,                 bestMatchedMCParticleIndex)                                                                                


// =========================================================================================================================================

// The macros that come below are all envokes with one of the above member lists
//     p = a prefix to a the variable name in the list (e.g. "metadata_run" is a branch added to the output trees)
//     q = an object that owns the member variable (e.g. truth.nuPdgCode)
//     t = the type of the member variable (e.g. int)
//     n = the name of the member variable (e.g. run)

// Define a macro that declares a member variable and an associated boolean to check if that variable has been set
#define UBCC1PI_MACRO_DECLARE_MEMBER(p, q, t, n)                                                                                           \
    Member<t> n;

// Define a macro that declares a member vector, used to bind particles to vectors in trees
#define UBCC1PI_MACRO_DECLARE_MEMBER_VECTOR(p, q, t, n)                                                                                    \
    std::vector<t> p##_##n##_vect;                                                                                                         \
    std::vector<bool> p##_##n##_isSet_vect;

// Define a macro that resets a member variable
#define UBCC1PI_MACRO_RESET_MEMBER(p, q, t, n)                                                                                             \
    q.n.Reset();

// Define a macro that resets a member vector
#define UBCC1PI_MACRO_RESET_MEMBER_VECTOR(p, q, t, n)                                                                                      \
    p##_##n##_vect.clear();                                                                                                                \
    p##_##n##_isSet_vect.clear();

// Define a macro that fills the output vectors
#define UBCC1PI_MACRO_FILL_MEMBER_VECTOR(p, q, t, n)                                                                                       \
    p##_##n##_vect.push_back(q.n.m_value);                                                                                                 \
    p##_##n##_isSet_vect.push_back(q.n.m_isSet);

// Define a macro that prints each of the member variables
#define UBCC1PI_MACRO_PRINT_MEMBER(p, q, t, n)                                                                                             \
    std::cout << std::setw(24) << "(" #t ")" << "  ";                                                                                      \
    std::cout << std::setw(40) << (#q "." #n "()") << "  ";                                                                                \
    std::cout << q.n.ToString() << std::endl;                                                                                              

// Define a macro to bind a member variable to an output branch
#define UBCC1PI_MACRO_BIND_OUTPUT_BRANCH(p, q, t, n)                                                                                       \
    pTree->Branch(#p "_" #n, &q.n.m_value);                                                                                                \
    pTree->Branch(#p "_" #n "_isSet", &q.n.m_isSet);

// Define a macro to bind a member variable to an output vector branch
#define UBCC1PI_MACRO_BIND_OUTPUT_VECTOR_BRANCH(p, q, t, n)                                                                                \
    pTree->Branch(#p "_" #n "_vect", &p##_##n##_vect);                                                                                     \
    pTree->Branch(#p "_" #n "_isSet_vect", &p##_##n##_isSet_vect);

#endif
