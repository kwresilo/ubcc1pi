/**
 *  @file  ubcc1pi_standalone/Macros/MakeSidebandEventSelectionTable.cxx
 *
 *  @brief The implementation file of the MakeSidebandEventSelectionTable macro
 */

#include "ubcc1pi_standalone/Macros/Macros.h"

#include "ubcc1pi_standalone/Objects/FileReader.h"

#include "ubcc1pi_standalone/Helpers/SelectionHelper.h"
#include "ubcc1pi_standalone/Helpers/AnalysisHelper.h"
#include "ubcc1pi_standalone/Helpers/NormalisationHelper.h"
#include "ubcc1pi_standalone/Helpers/FormattingHelper.h"

using namespace ubcc1pi;

const bool onlyTrueCC0piMode = false; 

namespace ubcc1pi_macros
{

void MakeSidebandEventSelectionTable(const Config &config)
{
    // Setup the input files
    std::vector< std::tuple<AnalysisHelper::SampleType, std::string, float> > inputData;

    // -------------------------------------------------------------------------------------------------------------------------------------
    // Setup the input files
    // -------------------------------------------------------------------------------------------------------------------------------------
    // Here we define a vector of tuples with 4 entries
    //   - First, the sample type (e.g. overlay)
    //   - Second, the path to the input file
    //   - Third, the normalisation factor to apply to all events in that file
    std::cout<<"##########################################\nUSING NUWRO AS DATA & Only CC0pi!\n##########################################"<<std::endl;
    for (const auto run: config.global.runs)
    {
        if(run == 1)
        {
            inputData.emplace_back(AnalysisHelper::Overlay, config.filesRun1.overlaysFileName, NormalisationHelper::GetOverlaysNormalisation(config, 1));
            inputData.emplace_back(AnalysisHelper::NuWro, config.filesRun1.nuWroFileName, NormalisationHelper::GetNuWroNormalisation(config, 1));
            if(!onlyTrueCC0piMode)
            {
                inputData.emplace_back(AnalysisHelper::Dirt,    config.filesRun1.dirtFileName, NormalisationHelper::GetDirtNormalisation(config, 1));
                inputData.emplace_back(AnalysisHelper::DataEXT, config.filesRun1.dataEXTFileName, NormalisationHelper::GetDataEXTNormalisation(config, 1));
                inputData.emplace_back(AnalysisHelper::DataBNB, config.filesRun1.dataBNBFileName, 1.f);
            }
        }
        else if(run == 2)
        {
            inputData.emplace_back(AnalysisHelper::Overlay, config.filesRun2.overlaysFileName, NormalisationHelper::GetOverlaysNormalisation(config, 2));
            inputData.emplace_back(AnalysisHelper::NuWro, config.filesRun2.nuWroFileName, NormalisationHelper::GetNuWroNormalisation(config, 1));
            if(!onlyTrueCC0piMode)
            {
                inputData.emplace_back(AnalysisHelper::Dirt,    config.filesRun2.dirtFileName, NormalisationHelper::GetDirtNormalisation(config, 2));
                inputData.emplace_back(AnalysisHelper::DataEXT, config.filesRun2.dataEXTFileName, NormalisationHelper::GetDataEXTNormalisation(config, 2));
                inputData.emplace_back(AnalysisHelper::DataBNB, config.filesRun2.dataBNBFileName, 1.f);
            }
        }
        else if(run == 3)
        {
            inputData.emplace_back(AnalysisHelper::Overlay, config.filesRun3.overlaysFileName, NormalisationHelper::GetOverlaysNormalisation(config, 3));
            inputData.emplace_back(AnalysisHelper::NuWro, config.filesRun3.nuWroFileName, NormalisationHelper::GetNuWroNormalisation(config, 1));
            if(!onlyTrueCC0piMode)
            {
                inputData.emplace_back(AnalysisHelper::Dirt,    config.filesRun3.dirtFileName, NormalisationHelper::GetDirtNormalisation(config, 3));
                inputData.emplace_back(AnalysisHelper::DataEXT, config.filesRun3.dataEXTFileName, NormalisationHelper::GetDataEXTNormalisation(config, 3));
                inputData.emplace_back(AnalysisHelper::DataBNB, config.filesRun3.dataBNBFileName, 1.f);
            }
        }
        else throw std::logic_error("PlotEventSelectionCuts - Invalid run number");
    }

    // for (const auto run: config.global.runs)
    // {
    //     if(run == 1)
    //     {
    //         inputData.emplace_back(AnalysisHelper::Overlay, config.filesRun1.overlaysFileName, NormalisationHelper::GetOverlaysNormalisation(config, 1));
    //         inputData.emplace_back(AnalysisHelper::Dirt,    config.filesRun1.dirtFileName, NormalisationHelper::GetDirtNormalisation(config, 1));
    //         inputData.emplace_back(AnalysisHelper::DataEXT, config.filesRun1.dataEXTFileName, NormalisationHelper::GetDataEXTNormalisation(config, 1));
    //         inputData.emplace_back(AnalysisHelper::DataBNB, config.filesRun1.dataBNBFileName, 1.f);
    //     }
    //     else if(run == 2)
    //     {
    //         inputData.emplace_back(AnalysisHelper::Overlay, config.filesRun2.overlaysFileName, NormalisationHelper::GetOverlaysNormalisation(config, 2));
    //         inputData.emplace_back(AnalysisHelper::Dirt,    config.filesRun2.dirtFileName, NormalisationHelper::GetDirtNormalisation(config, 2));
    //         inputData.emplace_back(AnalysisHelper::DataEXT, config.filesRun2.dataEXTFileName, NormalisationHelper::GetDataEXTNormalisation(config, 2));
    //         inputData.emplace_back(AnalysisHelper::DataBNB, config.filesRun2.dataBNBFileName, 1.f);
    //     }
    //     else if(run == 3)
    //     {
    //         inputData.emplace_back(AnalysisHelper::Overlay, config.filesRun3.overlaysFileName, NormalisationHelper::GetOverlaysNormalisation(config, 3));
    //         inputData.emplace_back(AnalysisHelper::Dirt,    config.filesRun3.dirtFileName, NormalisationHelper::GetDirtNormalisation(config, 3));
    //         inputData.emplace_back(AnalysisHelper::DataEXT, config.filesRun3.dataEXTFileName, NormalisationHelper::GetDataEXTNormalisation(config, 3));
    //         inputData.emplace_back(AnalysisHelper::DataBNB, config.filesRun3.dataBNBFileName, 1.f);
    //     }
    //     else throw std::logic_error("ExtractSidebandFit - Invalid run number");
    // }

    // Set up the selection
    // std::cout<<"..........................................\nUSING Modified CC0pi Selection: muonLikeProtonValue=-?.??f, barelyResemblingProtonValue=?.??f\n.........................................."<<std::endl;
    // auto selection = SelectionHelper::GetCC0piSelectionModified(-0.35f, 0.45f);
    // auto selection = SelectionHelper::GetCC0piSelectionModified(-0.25f, 0.09f);
    auto selection = SelectionHelper::GetSelection("CC0pi");

    // Setup the event counter
    AnalysisHelper::EventCounter counterBNB;
    AnalysisHelper::EventCounter counterNuWro;
    
    // Loop over the input files
    for (const auto &[sampleType, fileName, normalisation] : inputData)
    {
        // Open the input file for reading
        FileReader reader(fileName);
        auto pEvent = reader.GetBoundEventAddress();
        std::cout << "Processing file - " << fileName << std::endl;

        // Loop over the events
        const auto nEvents = reader.GetNumberOfEvents();
        
        std::cout<<"\n##############\nOnly counting every 25th event!\n##############"<<std::endl;
        for (unsigned int i = 0; i < nEvents/25; ++i)
        {
            AnalysisHelper::PrintLoadingBar(i, nEvents);
            reader.LoadEvent(i);

            if(onlyTrueCC0piMode)
            {
                const auto isTrueCC0Pi = AnalysisHelper::IsTrueCC0Pi(pEvent, config.global.useAbsPdg, config.global.protonMomentumThreshold); // todo remove this
                if(!isTrueCC0Pi) continue;
            }

            // Get the nominal event weight, scaled by the sample normalisation
            const auto weight = AnalysisHelper::GetNominalEventWeight(pEvent) * normalisation;

            // Run the selection
            const auto &[passedAllCuts, cutsPassed, assignedPdgCodes] = selection.Execute(pEvent);

            // For overlay & dirt samples only, count all events before any selection cuts
            // ATTN The EXT and BNB data has been filtered through the CC inclusive selection, so we can't get the count before any cuts were applied
            for(const auto &[counterType, pCounter]: {std::make_pair(AnalysisHelper::DataBNB, &counterBNB), std::make_pair(AnalysisHelper::NuWro, &counterNuWro)})
            {
                if(onlyTrueCC0piMode && counterType == AnalysisHelper::DataBNB) continue;
                if (sampleType != AnalysisHelper::Overlay
                    && ((counterType == AnalysisHelper::NuWro && sampleType != AnalysisHelper::NuWro) 
                    || (counterType == AnalysisHelper::DataBNB && !(sampleType == AnalysisHelper::DataBNB || sampleType == AnalysisHelper::Dirt || sampleType == AnalysisHelper::DataEXT))))
                    continue;

                if (sampleType == AnalysisHelper::Overlay || sampleType == AnalysisHelper::Dirt)
                {
                    pCounter->CountEvent("all", sampleType, pEvent, weight, true, config.global.protonMomentumThreshold);
                }

                const auto modifiedSampleType = (sampleType == AnalysisHelper::NuWro) ? AnalysisHelper::DataBNB : sampleType;
                // Check which of the cuts the event passed
                for (const auto &cutName : selection.GetCuts())
                {
                    // If the event didn't pass this cut, then move on
                    if (!SelectionHelper::IsCutPassed(cutsPassed, cutName))
                        break;

                    // Count the event, and use the cut name as the "tag"
                    pCounter->CountEvent(cutName, modifiedSampleType, pEvent, weight, true, config.global.protonMomentumThreshold);
                }
            }
        }
    }

    for(const auto &[postfix, pCounter]: {std::make_pair(std::string("BNB"), &counterBNB), std::make_pair(std::string("NuWro"), &counterNuWro)})
    {
        if(onlyTrueCC0piMode && postfix == "BNB") continue;
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

        table.WriteToFile("SidebandEventSelection_cuts_"+postfix+".md");

        // Print the breakdown of the event counts
        FormattingHelper::PrintLine();
        std::cout << "Summary" << std::endl;
        FormattingHelper::PrintLine();
        pCounter->PrintBreakdownSummary("SidebandEventSelection_summary_"+postfix+".md");

        FormattingHelper::PrintLine();
        std::cout << "Details" << std::endl;
        FormattingHelper::PrintLine();

        const auto nEntriesToPrint = 10u;
        pCounter->PrintBreakdownDetails("SidebandEventSelection_details_"+postfix+".md", nEntriesToPrint);
    }
}

} // namespace ubcc1pi_macros