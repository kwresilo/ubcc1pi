/**
 *  @file  ubcc1pi/Objects/EventMembers.h
 *
 *  @brief The file defining the members of the event class
 *
 *  In this file we define all of the member variable that will go into the analysis files. The types and names of the variables are defined
 *  using pre-processor directives. Normally this would be bad practise by modern standards, but it is the best option for this use case.
 *
 *  By defining the variables at the pre-processer level, all of the bookwork can be handled automatically. Each variable need only be
 *  defined in one place and will automatically:
 *      - Be added as a member variable of the ubcc1pi::Event class
 *      - Be bound to a branch if an analysis file is to be written or read using a consistent branch naming convention
 *      - Be assigned a default value and be checked such that you can't accidentally use a default variable at the analysis level 
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


#ifndef UBCC1PI_OBJECTS_EVENTMEMBERS
#define UBCC1PI_OBJECTS_EVENTMEMBERS

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
    f(p, q, bool,        isCC)                                                                                                             \
    f(p, q, int,         interactionMode)                                                                                                  \
    f(p, q, std::string, interactionString)                                                                                                \
    f(p, q, int,         nuPdgCode)                                                                                                        \
    f(p, q, float,       nuEnergy)                                                                                                         \
    f(p, q, TVector3,    nuVertex)                                                                                                         

// The event truth particle information members
#define UBCC1PI_MACRO_EVENT_TRUTH_PARTICLE_MEMBERS(p, q, f)                                                                                \
    f(p, q, int,         pdgCode)                                                                                                          \
    f(p, q, float,       momentumX)                                                                                                        \
    f(p, q, float,       momentumY)                                                                                                        \
    f(p, q, float,       momentumZ)                                                                                                        \
    f(p, q, float,       momentum)                                                                                                         \
    f(p, q, float,       energy)                                                                                                           \
    f(p, q, float,       mass)

// The event reco information members
#define UBCC1PI_MACRO_EVENT_RECO_MEMBERS(p, q, f)                                                                                          \
    f(p, q, bool,        dummy)                                                                                                            

// The event reco particle information members
#define UBCC1PI_MACRO_EVENT_RECO_PARTICLE_MEMBERS(p, q, f)                                                                                 \
    f(p, q, int,         pdgCode)                                                                                                          

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
    std::cout << std::setw(20) << "(" #t ")" << "  ";                                                                                      \
    std::cout << (#q "." #n) << "() : ";                                                                                                   \
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
