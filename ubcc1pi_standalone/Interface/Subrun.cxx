/**
 *  @file  ubcc1pi_standalone/Interface/Subrun.cxx
 *
 *  @brief The implementation of the subrun class
 */

#include "ubcc1pi_standalone/Interface/Subrun.h"

namespace ubcc1pi
{

void Subrun::Print() const
{
    std::cout << std::string(80, '=') << std::endl;

    std::cout << std::string(80, '-') << std::endl;
    std::cout << "SUBRUN" << std::endl;
    std::cout << std::string(80, '-') << std::endl;

    auto &subrun = *this;
    UBCC1PI_MACRO_SUBRUN_MEMBERS("", subrun, UBCC1PI_MACRO_PRINT_MEMBER)
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void Subrun::BindToOutputTree(TTree * pTree)
{
    auto &subrun = *this;
    UBCC1PI_MACRO_SUBRUN_MEMBERS(subrun, subrun, UBCC1PI_MACRO_BIND_OUTPUT_BRANCH)
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void Subrun::BindToInputTree(TTree * pTree)
{
    auto &subrun = *this;
    UBCC1PI_MACRO_SUBRUN_MEMBERS(subrun, subrun, UBCC1PI_MACRO_BIND_INPUT_BRANCH)
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void Subrun::Reset()
{
    auto &subrun = *this;
    UBCC1PI_MACRO_SUBRUN_MEMBERS("", subrun, UBCC1PI_MACRO_RESET_MEMBER)
}

} // namespace ubcc1pi
