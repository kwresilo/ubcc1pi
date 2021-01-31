/**
 *  @file  ubcc1pi_standalone/Macros/NMinusOneBDTStudyFull.cxx
 *
 *  @brief The implementation file of the NMinusOneBDTStudyFull macro
 */

#include "ubcc1pi_standalone/Macros/Macros.h"

using namespace ubcc1pi;

namespace ubcc1pi_macros
{

void NMinusOneBDTStudyFull(const Config &config)
{
    // ATTN here we run each iteration of the N-1 BDT study, and "hard-code" the feature to be removed in each step. If the BDT is changed
    // the NMinusOneBDTStudy should be executed manually step-by-step, and the most appropriate feature should be removed

    // Setup a custom config object based on the supplied configuration
    Config c = config;
    auto &featureNames = c.nMinusOneBDTStudy.featureNames;

    // Define a function to remove an element
    const auto Remove = [&](const std::string &feature)
    {
        std::cout << "Removing: " << feature << std::endl;
        featureNames.erase(std::remove(featureNames.begin(), featureNames.end(), feature), featureNames.end());
    };

    // Define a function to print the remaining elements
    const auto Print = [&]()
    {
        std::cout << "Remaining features: " << featureNames.size() << std::endl;

        for (const auto &feature : featureNames)
            std::cout << "  - " << feature << std::endl;
    };

    // Define the features to remove for each BDT in the order they should be removed
    std::map< PlottingHelper::PlotStyle, std::vector< std::pair<std::string, bool> > > featuresToRemove = {
        {
            PlottingHelper::Muon,
            {
                {"protonForward", false},
                {"muonForward", false},
                {"nSpacePointsNearEnd", false}
            }
        },
        {
            PlottingHelper::Proton,
            {
                {"muonForward", false},
                {"nDescendents", false},
                {"protonForward", true},
                {"nSpacePointsNearEnd", true}
            }
        },
        {
            PlottingHelper::GoldenPion,
            {
                {"protonForward", true},
                {"nSpacePointsNearEnd", true},
                {"muonForward", true}
            }
        }
    };

    // Iterate through the map and run the N-1 BDT study for each entry
    for (const auto &[signalType, features] : featuresToRemove)
    {
        c.nMinusOneBDTStudy.signalType = signalType;
        featureNames = config.nMinusOneBDTStudy.featureNames;

        for (const auto &[feature, shouldRunStudy] : features)
        {
            Print();

            if (shouldRunStudy)
            {
                NMinusOneBDTStudy(c);
            }
            else
            {
                std::cout << "Skipping." << std::endl;
            }

            Remove(feature);
        }

        // Do the final iteration
        Print();
        NMinusOneBDTStudy(c);
    }
}

} // namespace ubc1pi_macros
