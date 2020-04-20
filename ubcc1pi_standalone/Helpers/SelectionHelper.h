/**
 *  @file  ubcc1pi_standalone/Helpers/SelectionHelper.h
 *
 *  @brief The header file for the selection helper class
 */

#ifndef UBCC1PI_STANDALONE_HELPERS_SELECTION_HELPER
#define UBCC1PI_STANDALONE_HELPERS_SELECTION_HELPER

#include "ubcc1pi_standalone/Interface/Event.h"
#include <memory>
#include <string>
#include <functional>

namespace ubcc1pi
{

/**
 *  @brief  The selection helper class
 */
class SelectionHelper
{
    public:
        /**
         *  @brief  The event selection class
         */
        class EventSelection
        {
            public:

                // Notes to self
                //
                // This class should own an event counter (maybe, move that class to the SelectionHelper).
                // The user declares some cuts at the top of their code and sets which can be optimised and how (cut ranges, can be disabled, etc.)
                // The user then writes an infinite loop and inside it applies their event selection using the cuts managed by this object
                // At the end of the selection, the user calls the method to update the cut parameters
                // They keep interating until the the cuts can no longer be updated and break out of the infinite loop
                // The user then extracts the optimal cut parameters
                //
                //
                // How to optimize?
                //
                // Optimize for efficiency x purity
                // 
                // Start by optimizing each parameter in turn (in the declared order), disable it, enable it and do a gradient descent from
                // the nominal value to get the optimal result, save the result and move on to the next cut.
                // Once all cuts are dealt with, iterate from the first cut applying gradient descent from the optimal position.
                // Repeat until no cut can be further optimized.


                /**
                 *  @brief  Declare a cut which can be optionally be disabled during optimisation
                 *
                 *  @param  name the name of the cut
                 *  @param  canDisable if the cut can be disabled
                 *  @param  shouldOptimize if this cut should be optimized
                 */
                void DeclareCut(const std::string &name, const bool canDisable, const bool shouldOptimize);
                
                /**
                 *  @brief  Declare a cut which can be optionally be disabled, and have the cut value varied during optimisation
                 *
                 *  @param  name the name of the cut
                 *  @param  canDisable if the cut can be disabled
                 *  @param  nominal the nominal cut value used when we aren't optimizing the parameter
                 *  @param  min the minimum cut value to test
                 *  @param  max the maximum cut value to test
                 *  @param  shouldOptimize if this cut should be optimized
                 */
                void DeclareCut(const std::string &name, const bool canDisable, const float &nominal, const float &min, const float &max, const bool shouldOptimize);

                /**
                 *  @brief  Perform an an optimisation step for the cut parameters
                 *
                 *  @param  shouldPrint whether to print what is currently being optimised
                 *
                 *  @return boolean, true if the parameters changed, false if they are optimal
                 */
                bool UpdateCutParameters(const bool shouldPrint = true);

                /**
                 *  @brief  Apply the cut with the given name using the current optimisation parameters
                 *
                 *  @param  name the cut to apply
                 *  @param  method the actual cut method, supplied as a lambda function which returns true if the cut passes
                 *
                 *  @return the result of the cut
                 */
                bool ApplyCut(const std::string &name, const std::function<bool()> &method) const;
                
                /**
                 *  @brief  Apply the cut with the given name using the current optimisation parameters
                 *
                 *  @param  name the cut to apply
                 *  @param  method the actual cut method, supplied as a lambda function which returns true if the cut passes. The parameter of the method is the cut value
                 *
                 *  @return the result of the cut
                 */
                bool ApplyCut(const std::string &name, const std::function<bool(const float &)> &method) const;

                /**
                 *  @brief  Determine if the cut is should be enabled for optimal performance 
                 *
                 *  @param  name the name of the cut
                 *
                 *  @return boolean, true if enabled
                 */
                bool ShouldBeEnabled(const std::string &name) const;
                
                /**
                 *  @brief  Get the optimal value of the cut
                 *
                 *  @param  name the name of the cut
                 *
                 *  @return the value
                 */
                float GetOptimalCutValue(const std::string &name) const;
        };
};

// -----------------------------------------------------------------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------------------------------------------------------------


} // namespace ubcc1pi

#endif
