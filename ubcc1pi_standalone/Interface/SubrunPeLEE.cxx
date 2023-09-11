/**
 *  @file  ubcc1pi_standalone/Interface/Subrun.cxx
 *
 *  @brief The implementation of the subrun class
 */

#include "ubcc1pi_standalone/Interface/SubrunPeLEE.h"

namespace ubcc1pi
{

SubrunPeLEE::SubrunPeLEE(const bool hasTruthInfo) : hasTruthWeights(hasTruthInfo) {}

// -----------------------------------------------------------------------------------------------------------------------------------------

void SubrunPeLEE::Print() const
{
    std::cout << std::string(80, '=') << std::endl;

    std::cout << std::string(80, '-') << std::endl;
    std::cout << "SUBRUN" << std::endl;
    std::cout << std::string(80, '-') << std::endl;

    auto &subrun = *this;
    PELEE_MACRO_SUBRUN_MEMBERS("", subrun, PELEE_MACRO_PRINT_MEMBER)
    if(hasTruthWeights){PELEE_MACRO_SUBRUN_OPTIONAL_MEMBERS("", subrun, PELEE_MACRO_PRINT_MEMBER)}
}

// -----------------------------------------------------------------------------------------------------------------------------------------

bool SubrunPeLEE::HasTruthWeights() const
{
    return hasTruthWeights;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void SubrunPeLEE::BindToOutputTree(TTree * pTree)
{
    auto &subrun = *this;
    PELEE_MACRO_SUBRUN_MEMBERS(subrun, subrun, PELEE_MACRO_BIND_OUTPUT_BRANCH)
    if(hasTruthWeights){PELEE_MACRO_SUBRUN_OPTIONAL_MEMBERS(subrun, subrun, PELEE_MACRO_BIND_OUTPUT_BRANCH)}
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void SubrunPeLEE::BindToInputTree(TTree * pTree)
{
    auto &subrun = *this;
    PELEE_MACRO_SUBRUN_MEMBERS(subrun, subrun, PELEE_MACRO_BIND_INPUT_BRANCH)
    if(hasTruthWeights){PELEE_MACRO_SUBRUN_OPTIONAL_MEMBERS(subrun, subrun, PELEE_MACRO_BIND_INPUT_BRANCH)}
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void SubrunPeLEE::Reset()
{
    auto &subrun = *this;
    PELEE_MACRO_SUBRUN_MEMBERS("", subrun, PELEE_MACRO_RESET_MEMBER)
    if(hasTruthWeights){PELEE_MACRO_SUBRUN_OPTIONAL_MEMBERS("", subrun, PELEE_MACRO_RESET_MEMBER)}
}

} // namespace ubcc1pi
