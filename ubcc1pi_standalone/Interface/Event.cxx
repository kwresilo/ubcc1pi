/**
 *  @file  ubcc1pi_standalone/Interface/Event.cxx
 *
 *  @brief The implementation of the event class
 */

#include "ubcc1pi_standalone/Interface/Event.h"

namespace ubcc1pi
{

Event::Event()
{
    UBCC1PI_MACRO_EVENT_TRUTH_PARTICLE_MEMBERS(truth_particle, "", UBCC1PI_MACRO_INIT_MEMBER_VECTOR)
    UBCC1PI_MACRO_EVENT_RECO_PARTICLE_MEMBERS(reco_particle, "", UBCC1PI_MACRO_INIT_MEMBER_VECTOR)
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void Event::Print() const
{
    std::cout << std::string(80, '=') << std::endl;

    std::cout << std::string(80, '-') << std::endl;
    std::cout << "METADATA" << std::endl;
    std::cout << std::string(80, '-') << std::endl;

    UBCC1PI_MACRO_EVENT_METADATA_MEMBERS("", metadata, UBCC1PI_MACRO_PRINT_MEMBER)

    std::cout << std::string(80, '-') << std::endl;
    std::cout << "TRUTH" << std::endl;
    std::cout << std::string(80, '-') << std::endl;

    UBCC1PI_MACRO_EVENT_TRUTH_MEMBERS("", truth, UBCC1PI_MACRO_PRINT_MEMBER)

    for (const auto &particle : truth.particles)
    {
        std::cout << "TRUTH PARTICLE" << std::endl;
        UBCC1PI_MACRO_EVENT_TRUTH_PARTICLE_MEMBERS("", particle, UBCC1PI_MACRO_PRINT_MEMBER)
    }

    std::cout << std::string(80, '-') << std::endl;
    std::cout << "RECO" << std::endl;
    std::cout << std::string(80, '-') << std::endl;
    UBCC1PI_MACRO_EVENT_RECO_MEMBERS("", reco, UBCC1PI_MACRO_PRINT_MEMBER)

    for (const auto &particle : reco.particles)
    {
        std::cout << "RECO PARTICLE" << std::endl;
        UBCC1PI_MACRO_EVENT_RECO_PARTICLE_MEMBERS("", particle, UBCC1PI_MACRO_PRINT_MEMBER)
    }
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void Event::BindToOutputTree(TTree * pTree)
{
    UBCC1PI_MACRO_EVENT_METADATA_MEMBERS(metadata, metadata, UBCC1PI_MACRO_BIND_OUTPUT_BRANCH)
    UBCC1PI_MACRO_EVENT_TRUTH_MEMBERS(truth, truth, UBCC1PI_MACRO_BIND_OUTPUT_BRANCH)
    UBCC1PI_MACRO_EVENT_TRUTH_PARTICLE_MEMBERS(truth_particle, "", UBCC1PI_MACRO_BIND_OUTPUT_VECTOR_BRANCH)
    UBCC1PI_MACRO_EVENT_RECO_MEMBERS(reco, reco, UBCC1PI_MACRO_BIND_OUTPUT_BRANCH)
    UBCC1PI_MACRO_EVENT_RECO_PARTICLE_MEMBERS(reco_particle, "", UBCC1PI_MACRO_BIND_OUTPUT_VECTOR_BRANCH)
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void Event::BindToInputTree(TTree * pTree)
{
    UBCC1PI_MACRO_EVENT_METADATA_MEMBERS(metadata, metadata, UBCC1PI_MACRO_BIND_INPUT_BRANCH)
    UBCC1PI_MACRO_EVENT_TRUTH_MEMBERS(truth, truth, UBCC1PI_MACRO_BIND_INPUT_BRANCH)
    UBCC1PI_MACRO_EVENT_TRUTH_PARTICLE_MEMBERS(truth_particle, "", UBCC1PI_MACRO_BIND_INPUT_VECTOR_BRANCH)
    UBCC1PI_MACRO_EVENT_RECO_MEMBERS(reco, reco, UBCC1PI_MACRO_BIND_INPUT_BRANCH)
    UBCC1PI_MACRO_EVENT_RECO_PARTICLE_MEMBERS(reco_particle, "", UBCC1PI_MACRO_BIND_INPUT_VECTOR_BRANCH)
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void Event::Reset()
{
    UBCC1PI_MACRO_EVENT_METADATA_MEMBERS("", metadata, UBCC1PI_MACRO_RESET_MEMBER)

    UBCC1PI_MACRO_EVENT_TRUTH_MEMBERS("", truth, UBCC1PI_MACRO_RESET_MEMBER)
    UBCC1PI_MACRO_EVENT_TRUTH_PARTICLE_MEMBERS(truth_particle, "", UBCC1PI_MACRO_RESET_MEMBER_VECTOR)
    truth.particles.clear();

    UBCC1PI_MACRO_EVENT_RECO_MEMBERS("", reco, UBCC1PI_MACRO_RESET_MEMBER)
    UBCC1PI_MACRO_EVENT_RECO_PARTICLE_MEMBERS(reco_particle, "", UBCC1PI_MACRO_RESET_MEMBER_VECTOR)
    reco.particles.clear();
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void Event::PrepareForTreeFill()
{
    for (const auto &particle : truth.particles)
    {
        UBCC1PI_MACRO_EVENT_TRUTH_PARTICLE_MEMBERS(truth_particle, particle, UBCC1PI_MACRO_FILL_MEMBER_VECTOR)
    }

    for (const auto &particle : reco.particles)
    {
        UBCC1PI_MACRO_EVENT_RECO_PARTICLE_MEMBERS(reco_particle, particle, UBCC1PI_MACRO_FILL_MEMBER_VECTOR)
    }
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void Event::PrepareAfterTreeRead()
{
    unsigned int nTruthParticles;
    UBCC1PI_MACRO_EVENT_TRUTH_PARTICLE_MEMBERS(truth_particle, &nTruthParticles, UBCC1PI_MACRO_GET_MEMBER_VECTOR_SIZE)

    truth.particles.resize(nTruthParticles);
    UBCC1PI_MACRO_EVENT_TRUTH_PARTICLE_MEMBERS(truth_particle, truth.particles, UBCC1PI_MACRO_READ_MEMBER_VECTOR)

    unsigned int nRecoParticles;
    UBCC1PI_MACRO_EVENT_RECO_PARTICLE_MEMBERS(reco_particle, &nRecoParticles, UBCC1PI_MACRO_GET_MEMBER_VECTOR_SIZE)

    reco.particles.resize(nRecoParticles);
    UBCC1PI_MACRO_EVENT_RECO_PARTICLE_MEMBERS(reco_particle, reco.particles, UBCC1PI_MACRO_READ_MEMBER_VECTOR)
}


} // namespace ubcc1pi
