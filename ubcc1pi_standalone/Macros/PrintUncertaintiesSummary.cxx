/**
 *  @file  ubcc1pi_standalone/Macros/PrintUncertaintiesSummary.cxx
 *
 *  @brief The implementation file of the PrintUncertaintiesSummary macro
 */

#include "ubcc1pi_standalone/Macros/Macros.h"
#include "ubcc1pi_standalone/Objects/FileReader.h"

#include "ubcc1pi_standalone/Helpers/FormattingHelper.h"
#include "ubcc1pi_standalone/Helpers/CrossSectionHelper.h"

#include "ubsmear.h"

using namespace ubcc1pi;

namespace ubcc1pi_macros
{

void PrintUncertaintiesSummary(const Config &config)
{
    // Define a lambda function to get a single value from file
    const auto getValue = [&](const std::string &identifier) -> float {
        const auto vect = ubsmear::UBFileHelper::ReadColumnVector("xsec_generic_total_" + identifier + ".txt");
        if (vect.GetRows() != 1)
            throw std::logic_error("PrintUncertaintiesSummary - Input file contains more than one value");

        return vect.At(0,0);
    };

    // Read in the values from the files produced by ExtractXSecs
    const auto xsecData = getValue("data");
    std::cout << "Flux-integrated, forward-folded cross section: " << xsecData << std::endl;

    FormattingHelper::Table table({"Group", "Parameter", "", "Universes", "", "Bias", "Std", "Total", "", "Frac"});

    // Keep track of the grand-total uncertainty
    float grandTotalSqr = 0.f;

    // Get the uncertainties for the multisims
    for (const auto &[group, dimensions] : std::map<std::string, CrossSectionHelper::SystDimensionsMap>({
        {"flux", config.extractXSecs.fluxDimensions},
        {"xsec", config.extractXSecs.xsecDimensions},
        {"misc", {
            {"bootstrap", config.extractXSecs.nBootstrapUniverses},
            {"dirt", 2}
        }}
    }))
    {
        float groupTotalSqr = 0.f;
        for (const auto &[paramName, nUniverses] : dimensions)
        {
            const auto bias = getValue("data_" + group + "_" + paramName + "_bias");
            const auto variance = getValue("data_" + group + "_" + paramName + "_covariance");
            const auto standardDeviation = std::pow(variance, 0.5f);

            const auto totalSqr = bias*bias + variance;
            const auto total = std::pow(totalSqr, 0.5f);

            const auto frac = total / xsecData;

            table.AddEmptyRow();
            table.SetEntry("Group", group);
            table.SetEntry("Parameter", paramName);
            table.SetEntry("Universes", nUniverses);
            table.SetEntry("Bias", bias);
            table.SetEntry("Std", standardDeviation);
            table.SetEntry("Total", total);
            table.SetEntry("Frac", frac);

            groupTotalSqr += totalSqr;
        }

        const auto groupTotal = std::pow(groupTotalSqr, 0.5f);
        const auto groupTotalFrac = groupTotal / xsecData;

        table.AddEmptyRow();
        table.SetEntry("Group", group);
        table.SetEntry("Total", groupTotal);
        table.SetEntry("Frac", groupTotalFrac);
        table.AddEmptyRow();

        grandTotalSqr += groupTotalSqr;
    }

    // Get the uncertainties for the unisims
    for (const auto &[group, dimensions] : std::map<std::string, CrossSectionHelper::SystUnisimDimensionsMap>({
        {"detector", config.extractXSecs.detVarDimensions}
    }))
    {
        float groupTotalSqr = 0.f;
        for (const auto &[paramName, cvName] : dimensions)
        {
            const auto bias = getValue("data_" + group + "_" + paramName + "_bias");

            const auto totalSqr = bias*bias;
            const auto total = std::pow(totalSqr, 0.5f);

            const auto frac = total / xsecData;

            table.AddEmptyRow();
            table.SetEntry("Group", group);
            table.SetEntry("Parameter", paramName);
            table.SetEntry("Bias", bias);
            table.SetEntry("Total", total);
            table.SetEntry("Frac", frac);

            groupTotalSqr += totalSqr;
        }

        const auto groupTotal = std::pow(groupTotalSqr, 0.5f);
        const auto groupTotalFrac = groupTotal / xsecData;

        table.AddEmptyRow();
        table.SetEntry("Group", group);
        table.SetEntry("Total", groupTotal);
        table.SetEntry("Frac", groupTotalFrac);
        table.AddEmptyRow();

        grandTotalSqr += groupTotalSqr;
    }

    // Add the grand total
    const auto grandTotal = std::pow(grandTotalSqr, 0.5f);
    const auto grandTotalFrac = grandTotal / xsecData;

    table.AddEmptyRow();
    table.SetEntry("Group", "all");
    table.SetEntry("Total", grandTotal);
    table.SetEntry("Frac", grandTotalFrac);

    table.WriteToFile("xsec_uncertaintiesSummary.md");
}

} // namespace ubcc1pi_macros
