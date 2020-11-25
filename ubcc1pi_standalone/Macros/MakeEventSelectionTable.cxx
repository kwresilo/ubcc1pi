/**
 *  @file  ubcc1pi_standalone/Macros/MakeEventSelectionTable.cxx
 *
 *  @brief The implementation file of the MakeEventSelectionTable macro
 */

#include "ubcc1pi_standalone/Macros/Macros.h"

#include "ubcc1pi_standalone/Objects/FileReader.h"

#include "ubcc1pi_standalone/Helpers/SelectionHelper.h"
#include "ubcc1pi_standalone/Helpers/AnalysisHelper.h"
#include "ubcc1pi_standalone/Helpers/NormalisationHelper.h"
#include "ubcc1pi_standalone/Helpers/FormattingHelper.h"

using namespace ubcc1pi;

namespace ubcc1pi_macros
{

void MakeEventSelectionTable(const Config &config)
{
    // Setup the input files
    std::vector< std::tuple<AnalysisHelper::SampleType, std::string, float> > inputData;

    inputData.emplace_back(AnalysisHelper::Overlay, config.files.overlaysFileName, NormalisationHelper::GetOverlaysNormalisation(config));
    inputData.emplace_back(AnalysisHelper::Dirt,    config.files.dirtFileName,     NormalisationHelper::GetDirtNormalisation(config));
    inputData.emplace_back(AnalysisHelper::DataEXT, config.files.dataEXTFileName,  NormalisationHelper::GetDataEXTNormalisation(config));
    inputData.emplace_back(AnalysisHelper::DataBNB, config.files.dataBNBFileName,  1.f);

    // Set up the selection
    auto selection = SelectionHelper::GetDefaultSelection();

    // Setup the event counter
    AnalysisHelper::EventCounter counter;

    // Loop over the input files
    for (const auto &[sampleType, fileName, normalisation] : inputData)
    {
        // Open the input file for reading
        FileReader reader(fileName);
        auto pEvent = reader.GetBoundEventAddress();
        std::cout << "Processing file - " << fileName << std::endl;

        // Loop over the events
        const auto nEvents = reader.GetNumberOfEvents();
        for (unsigned int i = 0; i < nEvents; ++i)
        {
            AnalysisHelper::PrintLoadingBar(i, nEvents);

            reader.LoadEvent(i);

            // Get the nominal event weight, scaled by the sample normalisation
            const auto weight = AnalysisHelper::GetNominalEventWeight(pEvent) * normalisation;

            // For overlay & dirt samples only, count all events before any selection cuts
            // ATTN The EXT and BNB data has been filtered through the CC inclusive selection, so we can't get the count before any cuts were applied
            if (sampleType == AnalysisHelper::Overlay || sampleType == AnalysisHelper::Dirt)
            {
                counter.CountEvent("all", sampleType, pEvent, weight);
            }

            // Run the selection
            const auto &[passedAllCuts, cutsPassed, assignedPdgCodes] = selection.Execute(pEvent);

            // Check which of the cuts the event passed
            for (const auto &cutName : selection.GetCuts())
            {
                // If the event didn't pass this cut, then skip it
                // TODO probably better to break here instead of continue - check this
                if (!SelectionHelper::IsCutPassed(cutsPassed, cutName))
                    continue;

                // Count the event, and use the cut name as the "tag"
                counter.CountEvent(cutName, sampleType, pEvent, weight);
            }
        }
    }

    // Print the cuts that were used
    FormattingHelper::PrintLine();
    std::cout << "Cuts" << std::endl;
    FormattingHelper::PrintLine();

    FormattingHelper::Table table({"Cut", "", "Value"});
    for (const auto &cutName : selection.GetCuts())
    {
        table.AddEmptyRow();
        table.SetEntry("Cut", cutName);

        if (selection.CutHasValue(cutName))
            table.SetEntry("Value", selection.GetCutValue(cutName));
    }

    table.WriteToFile("eventSelection_cuts.md");

    // Print the breakdown of the event counts
    FormattingHelper::PrintLine();
    std::cout << "Summary" << std::endl;
    FormattingHelper::PrintLine();
    counter.PrintBreakdownSummary("eventSelection_summary.md");

    FormattingHelper::PrintLine();
    std::cout << "Details" << std::endl;
    FormattingHelper::PrintLine();

    const auto nEntriesToPrint = 10u;
    counter.PrintBreakdownDetails("eventSelection_details.md", nEntriesToPrint);
}

} // namespace ubcc1pi_macros
