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
#include <map>
#include <string>
#include <memory>
#include <TVector3.h>

// In the definition of the members below the first two parameters (p, q) are used by later macros, the third parameter is a boolean which
// should be true if ROOT needs the address of a pointer when setting branch addresses (for non-standard types). The fourth is the type of
// the varaiable, and the fifth is the variable name itself

// The event metadata members
#define UBCC1PI_MACRO_EVENT_METADATA_MEMBERS(p, q, f)                                                                                      \
    f(p, q, false, int,         run)                                                                                                       \
    f(p, q, false, int,         subRun)                                                                                                    \
    f(p, q, false, int,         event)                                                                                                     \
    f(p, q, false, bool,        hasTruthInfo)                                                                                              \

// The event truth information members
//    f(p, q, false, float,                      splineEventWeight)                                                                          \
//    f(p, q, false, float,                      genieTuneEventWeight)                                                                       \

#define UBCC1PI_MACRO_EVENT_TRUTH_MEMBERS(p, q, f)                                                                                         \
    f(p, q, false, bool,                       isCC)                                                                                       \
    f(p, q, false, int,                        interactionMode)                                                                            \
    f(p, q, true,  std::string,                interactionString)                                                                          \
    f(p, q, false, int,                        nuPdgCode)                                                                                  \
    f(p, q, false, float,                      nuEnergy)                                                                                   \
    f(p, q, true,  TVector3,                   nuVertex)                                                                                   \
    f(p, q, false, int,                        nFinalStates)                                                                               \
    f(p, q, true,  std::vector<float>,         slicePurities)                                                                              \
    f(p, q, true,  std::vector<float>,         sliceCompletenesses)                                                                               

// The event truth particle information members
#define UBCC1PI_MACRO_EVENT_TRUTH_PARTICLE_MEMBERS(p, q, f)                                                                                \
    f(p, q, false, int,                 pdgCode)                                                                                           \
    f(p, q, false, float,               startX)                                                                                            \
    f(p, q, false, float,               startY)                                                                                            \
    f(p, q, false, float,               startZ)                                                                                            \
    f(p, q, false, float,               endX)                                                                                              \
    f(p, q, false, float,               endY)                                                                                              \
    f(p, q, false, float,               endZ)                                                                                              \
    f(p, q, false, float,               momentumX)                                                                                         \
    f(p, q, false, float,               momentumY)                                                                                         \
    f(p, q, false, float,               momentumZ)                                                                                         \
    f(p, q, false, float,               momentum)                                                                                          \
    f(p, q, false, float,               energy)                                                                                            \
    f(p, q, false, float,               mass)                                                                                              \
    f(p, q, false, float,               hitWeightU)                                                                                        \
    f(p, q, false, float,               hitWeightV)                                                                                        \
    f(p, q, false, float,               hitWeightW)                                                                                        \
    f(p, q, false, int,                 nElasticScatters)                                                                                  \
    f(p, q, false, int,                 nInelasticScatters)                                                                                \
    f(p, q, false, std::vector<float>,  scatterCosThetas)                                                                                  \
    f(p, q, false, std::vector<float>,  scatterMomentumFracsLost)                                                                          \
    f(p, q, false, std::vector<int>,    scatterIsElastic)                                                                                  \
    f(p, q, false, float,               scatteredMomentum)                                                                                 \
    f(p, q, false, float,               endMomentum)                                                                                       \
    f(p, q, false, bool,                isStopping)                                                                                        \
    f(p, q, false, int,                 endState)                                                                                          \
    f(p, q, false, float,               endStateProductsHitWeightU)                                                                        \
    f(p, q, false, float,               endStateProductsHitWeightV)                                                                        \
    f(p, q, false, float,               endStateProductsHitWeightW)                                                                        
    
// The event reco information members
#define UBCC1PI_MACRO_EVENT_RECO_MEMBERS(p, q, f)                                                                                          \
    f(p, q, false, bool,                passesCCInclusive)                                                                                 \
    f(p, q, false, int,                 nSlices)                                                                                           \
    f(p, q, false, bool,                hasSelectedSlice)                                                                                  \
    f(p, q, false, float,               selectedTopologicalScore)                                                                          \
    f(p, q, true,  std::vector<float>,  sliceTopologicalScores)                                                                            \
    f(p, q, true,  std::vector<bool>,   sliceIsSelectedAsNu)                                                                               \
    f(p, q, false, bool,                hasNeutrino)                                                                                       \
    f(p, q, false, int,                 nuPdgCode)                                                                                         \
    f(p, q, true,  TVector3,            nuVertex)                                                                                          \
    f(p, q, false, int,                 nFinalStates)                                                                                                

// The event reco particle information members
#define UBCC1PI_MACRO_EVENT_RECO_PARTICLE_MEMBERS(p, q, f)                                                                                 \
    f(p, q, false, bool,                isCCInclusiveMuonCandidate)                                                                        \
    f(p, q, false, int,                 pdgCode)                                                                                           \
    f(p, q, false, int,                 nHitsU)                                                                                            \
    f(p, q, false, int,                 nHitsV)                                                                                            \
    f(p, q, false, int,                 nHitsW)                                                                                            \
    f(p, q, false, int,                 nDaughters)                                                                                        \
    f(p, q, false, int,                 nDescendents)                                                                                      \
    f(p, q, false, int,                 nDescendentHitsU)                                                                                  \
    f(p, q, false, int,                 nDescendentHitsV)                                                                                  \
    f(p, q, false, int,                 nDescendentHitsW)                                                                                  \
    f(p, q, false, int,                 nHitsInLargestDescendent)                                                                          \
    f(p, q, false, float,               trackScore)                                                                                        \
    f(p, q, false, float,               startX)                                                                                            \
    f(p, q, false, float,               startY)                                                                                            \
    f(p, q, false, float,               startZ)                                                                                            \
    f(p, q, false, float,               endX)                                                                                              \
    f(p, q, false, float,               endY)                                                                                              \
    f(p, q, false, float,               endZ)                                                                                              \
    f(p, q, false, float,               directionX)                                                                                        \
    f(p, q, false, float,               directionY)                                                                                        \
    f(p, q, false, float,               directionZ)                                                                                        \
    f(p, q, false, float,               yzAngle)                                                                                           \
    f(p, q, false, float,               xyAngle)                                                                                           \
    f(p, q, false, float,               xzAngle)                                                                                           \
    f(p, q, false, float,               length)                                                                                            \
    f(p, q, false, float,               range)                                                                                             \
    f(p, q, false, float,               transverseVertexDist)                                                                              \
    f(p, q, false, float,               longitudinalVertexDist)                                                                            \
    f(p, q, false, float,               wiggliness)                                                                                        \
    f(p, q, false, int,                 nSpacePointsNearEnd)                                                                               \
    f(p, q, false, float,               likelihoodForwardMuonU)                                                                            \
    f(p, q, false, float,               likelihoodForwardMuonV)                                                                            \
    f(p, q, false, float,               likelihoodForwardMuonW)                                                                            \
    f(p, q, false, float,               likelihoodForwardMuon)                                                                             \
    f(p, q, false, float,               likelihoodBackwardMuonU)                                                                           \
    f(p, q, false, float,               likelihoodBackwardMuonV)                                                                           \
    f(p, q, false, float,               likelihoodBackwardMuonW)                                                                           \
    f(p, q, false, float,               likelihoodBackwardMuon)                                                                            \
    f(p, q, false, float,               likelihoodForwardPionU)                                                                            \
    f(p, q, false, float,               likelihoodForwardPionV)                                                                            \
    f(p, q, false, float,               likelihoodForwardPionW)                                                                            \
    f(p, q, false, float,               likelihoodForwardPion)                                                                             \
    f(p, q, false, float,               likelihoodBackwardPionU)                                                                           \
    f(p, q, false, float,               likelihoodBackwardPionV)                                                                           \
    f(p, q, false, float,               likelihoodBackwardPionW)                                                                           \
    f(p, q, false, float,               likelihoodBackwardPion)                                                                            \
    f(p, q, false, float,               likelihoodForwardProtonU)                                                                          \
    f(p, q, false, float,               likelihoodForwardProtonV)                                                                          \
    f(p, q, false, float,               likelihoodForwardProtonW)                                                                          \
    f(p, q, false, float,               likelihoodForwardProton)                                                                           \
    f(p, q, false, float,               likelihoodBackwardProtonU)                                                                         \
    f(p, q, false, float,               likelihoodBackwardProtonV)                                                                         \
    f(p, q, false, float,               likelihoodBackwardProtonW)                                                                         \
    f(p, q, false, float,               likelihoodBackwardProton)                                                                          \
    f(p, q, false, float,               likelihoodMIPU)                                                                                    \
    f(p, q, false, float,               likelihoodMIPV)                                                                                    \
    f(p, q, false, float,               likelihoodMIPW)                                                                                    \
    f(p, q, false, float,               likelihoodMIP)                                                                                     \
    f(p, q, false, float,               truncatedMeandEdxU)                                                                                \
    f(p, q, false, float,               truncatedMeandEdxV)                                                                                \
    f(p, q, false, float,               truncatedMeandEdxW)                                                                                \
    f(p, q, false, float,               truncatedMeandEdx)                                                                                 \
    f(p, q, false, std::vector<float>,  truthMatchPurities)                                                                                \
    f(p, q, false, std::vector<float>,  truthMatchCompletenesses)                                                                          \
    f(p, q, false, bool,                hasMatchedMCParticle)                                                                              \
    f(p, q, false, int,                 bestMatchedMCParticleIndex)                                                                         

#endif
