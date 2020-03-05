/**
 *  @file  ubcc1pi/Objects/Event.cxx
 *
 *  @brief The implementation of the event class
 */

#include "ubcc1pi/Objects/Event.h"

#include <TFile.h>

namespace ubcc1pi
{
        
void Event::Print() const
{
    std::cout << std::string(80, '=');
    std::cout << "METADATA" << std::endl;
    UBCC1PI_MACRO_EVENT_METADATA_MEMBERS("", metadata, UBCC1PI_MACRO_PRINT_MEMBER)
    
    std::cout << "TRUTH" << std::endl;
    UBCC1PI_MACRO_EVENT_TRUTH_MEMBERS("", truth, UBCC1PI_MACRO_PRINT_MEMBER)

    for (const auto &particle : truth.particles)
    {
        std::cout << "TRUTH PARTICLE" << std::endl;
        UBCC1PI_MACRO_EVENT_TRUTH_PARTICLE_MEMBERS("", particle, UBCC1PI_MACRO_PRINT_MEMBER)
    }
}

// -----------------------------------------------------------------------------------------------------------------------------------------
        
void Event::BindToOutputTree(TTree * pTree)
{
    UBCC1PI_MACRO_EVENT_METADATA_MEMBERS(metadata, metadata, UBCC1PI_MACRO_BIND_OUTPUT_BRANCH)
    UBCC1PI_MACRO_EVENT_TRUTH_MEMBERS(truth, truth, UBCC1PI_MACRO_BIND_OUTPUT_BRANCH)
    UBCC1PI_MACRO_EVENT_TRUTH_PARTICLE_MEMBERS(truth_particle, "", UBCC1PI_MACRO_BIND_OUTPUT_VECTOR_BRANCH)
}

// -----------------------------------------------------------------------------------------------------------------------------------------
        
void Event::Reset()
{
    UBCC1PI_MACRO_EVENT_METADATA_MEMBERS("", metadata, UBCC1PI_MACRO_RESET_MEMBER)
    UBCC1PI_MACRO_EVENT_TRUTH_MEMBERS("", truth, UBCC1PI_MACRO_RESET_MEMBER)

    UBCC1PI_MACRO_EVENT_TRUTH_PARTICLE_MEMBERS(truth_particle, "", UBCC1PI_MACRO_RESET_MEMBER_VECTOR)
    truth.particles.clear();
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void Event::PrepareForTreeFill()
{
    for (const auto &particle : truth.particles)
    {
        UBCC1PI_MACRO_EVENT_TRUTH_PARTICLE_MEMBERS(truth_particle, particle, UBCC1PI_MACRO_FILL_MEMBER_VECTOR)
    }
}


} // namespace ubcc1pi
