/**
 *  @file  ubcc1pi_standalone/Interface/EventConversionMembers.h
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


#ifndef UBCC1PI_STANDALONE_INTERFACE_EVENTCONVERSIONMEMBERS
#define UBCC1PI_STANDALONE_INTERFACE_EVENTCONVERSIONMEMBERS

#include <iostream>
#include <iomanip>
#include <vector>
#include <map>
#include <string>
#include <memory>
#include <TVector3.h>

#include "ubcc1pi/Helpers/DebugHelper.h"

// In the definition of the members below the first two parameters (p, q) are used by later macros, the third parameter is a boolean which
// should be true if ROOT needs the address of a pointer when setting branch addresses (for non-standard types). The fourth is the type of
// the varaiable, and the fifth is the variable name itself

/** The event metadata members */
#define PELEE_TO_UBCC1PI_MACRO_EVENT_METADATA(p, q, f)                                                                                      \
    f(p.run(),    q.run,          false)                                                                                                    \
    f(p.sub(),    q.subRun,       false)                                                                                                    \
    f(p.evt(),    q.event,        false)
    //f(true,     q.hasTruthInfo, false) /*Placeholder*/

/** The event truth information members */
#define COMMA , // This is needed for the std::map macro below since it uses commas to separate arguments
#define PELEE_TO_UBCC1PI_MACRO_EVENT_TRUTH(p, q, f)                                                                                         \
    f(p.weightSpline.Get(),  q.splineEventWeight,   false)                                                                                                  \
    f(p.weightSplineTimesTune.Get(),  q.genieTuneEventWeight,   false)                                                                                                  \
    f( ,  q.systParamNames,   false)                                                                                                  \
    f( ,  q.systParamFirstValueIndex,   false)                                                                                                  \
    f( ,  q.systParamValues,   false)                                                                                                  \
    f(p.ccnc.Get() == 0,  q.isCC,   false)                                                                                                  \
    f(p.interaction.Get(),  q.interactionMode,   false)                                                                                                  \
    f(DebugHelper::GetInteractionString(p.interaction.Get() COMMA true),  q.interactionString,   false)                                                                                                  \
    f( ,  q.nuPdgCode,   false)                                                                                                  \
    f( ,  q.nuEnergy,   false)                                                                                                  \
    f( ,  q.nuVertex,   false)                                                                                                  \
    f( ,  q.nFinalStates,   false)                                                                                                  \
    f( ,  q.slicePurities,   false)                                                                                                  \
    f( ,  q.sliceCompletenesses,   false)                                                                                                  \



weights
weightsFlux
weightsGenie
weightsReint

interaction
nu_pdg
nu_e
true_nu_vtx_x
true_nu_vtx_y
true_nu_vtx_z
mc_purity
mc_completeness


/** The event truth particle information members */
#define PELEE_TO_UBCC1PI_MACRO_EVENT_TRUTH_PARTICLE(p, q, f)                                                                                \

/** The event reco information members */
#define PELEE_TO_UBCC1PI_MACRO_EVENT_RECO(p, q, f)                                                                                          \


/** The event reco particle information members */
#define PELEE_TO_UBCC1PI_MACRO_EVENT_RECO_PARTICLE(p, q, f)                                                                                 \


#endif
