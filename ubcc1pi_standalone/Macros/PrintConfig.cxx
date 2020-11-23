/**
 *  @file  ubcc1pi_standalone/Macros/PrintConfig.cxx
 *
 *  @brief The implementation file of the PrintConfig macro
 */

#include "ubcc1pi_standalone/Macros/Macros.h"

#include "ubcc1pi_standalone/Objects/FileReader.h"

#include "ubcc1pi_standalone/Helpers/FormattingHelper.h"

using namespace ubcc1pi;

namespace ubcc1pi_macros
{

void PrintConfig(const Config &config)
{
    FormattingHelper::Table table({"Variable", "Value"});

    // TODO make sure this is up to date!

    //
    // Files
    //
    table.AddEmptyRow();
    table.SetEntry("Variable", "files.overlaysFileName");
    table.SetEntry("Value", config.files.overlaysFileName);

    table.AddEmptyRow();
    table.SetEntry("Variable", "files.dirtFileName");
    table.SetEntry("Value", config.files.dirtFileName);

    table.AddEmptyRow();
    table.SetEntry("Variable", "files.dataEXTFileName");
    table.SetEntry("Value", config.files.dataEXTFileName);

    table.AddEmptyRow();
    table.SetEntry("Variable", "files.dataBNBFileName");
    table.SetEntry("Value", config.files.dataBNBFileName);


    //
    // Norms
    //
    table.AddEmptyRow();

    table.AddEmptyRow();
    table.SetEntry("Variable", "norms.overlaysPOT");
    table.SetEntry("Value", config.norms.overlaysPOT);

    table.AddEmptyRow();
    table.SetEntry("Variable", "norms.dirtPOT");
    table.SetEntry("Value", config.norms.dirtPOT);

    table.AddEmptyRow();
    table.SetEntry("Variable", "norms.dataEXTTriggers");
    table.SetEntry("Value", config.norms.dataEXTTriggers);

    table.AddEmptyRow();
    table.SetEntry("Variable", "norms.dataBNBTor875WCut");
    table.SetEntry("Value", config.norms.dataBNBTor875WCut);

    table.AddEmptyRow();
    table.SetEntry("Variable", "norms.dataBNBE1DCNTWCut");
    table.SetEntry("Value", config.norms.dataBNBE1DCNTWCut);

    //
    // Global
    //
    table.AddEmptyRow();

    table.AddEmptyRow();
    table.SetEntry("Variable", "global.useAbsPdg");
    table.SetEntry("Value", config.global.useAbsPdg);

    table.AddEmptyRow();
    table.SetEntry("Variable", "global.countProtonsInclusively");
    table.SetEntry("Value", config.global.countProtonsInclusively);

    table.AddEmptyRow();
    table.SetEntry("Variable", "global.lastCutGeneric");
    table.SetEntry("Value", config.global.lastCutGeneric);

    table.AddEmptyRow();
    table.SetEntry("Variable", "global.protonMomentumThreshold");
    table.SetEntry("Value", config.global.protonMomentumThreshold);

    //
    // CountPOT
    //
    table.AddEmptyRow();

    table.AddEmptyRow();
    table.SetEntry("Variable", "countPOT.useOverlays");
    table.SetEntry("Value", config.countPOT.useOverlays);

    table.AddEmptyRow();
    table.SetEntry("Variable", "countPOT.useDirt");
    table.SetEntry("Value", config.countPOT.useDirt);

    //
    // GetRunSubrunList
    //
    table.AddEmptyRow();

    table.AddEmptyRow();
    table.SetEntry("Variable", "getRunSubrunList.useDataEXT");
    table.SetEntry("Value", config.getRunSubrunList.useDataEXT);

    table.AddEmptyRow();
    table.SetEntry("Variable", "getRunSubrunList.useDataBNB");
    table.SetEntry("Value", config.getRunSubrunList.useDataBNB);

    //
    // MultiPlanePIDDemo
    //
    table.AddEmptyRow();

    table.AddEmptyRow();
    table.SetEntry("Variable", "multiPlanePIDDemo.sin2AngleThreshold");
    table.SetEntry("Value", config.multiPlanePIDDemo.sin2AngleThreshold);

    //
    // PlotInputVariables
    //
    table.AddEmptyRow();

    table.AddEmptyRow();
    table.SetEntry("Variable", "plotInputVariables.plotBDTResponses");
    table.SetEntry("Value", config.plotInputVariables.plotBDTResponses);

    //
    // NMinusOneBDTStudy
    //
    table.AddEmptyRow();

    table.AddEmptyRow();
    table.SetEntry("Variable", "nMinusOneBDTStudy.shouldTrainBDTs");
    table.SetEntry("Value", config.nMinusOneBDTStudy.shouldTrainBDTs);

    table.AddEmptyRow();
    table.SetEntry("Variable", "nMinusOneBDTStudy.signalType");
    table.SetEntry("Value", config.nMinusOneBDTStudy.signalType);

    for (unsigned int i = 0; i < config.nMinusOneBDTStudy.featureNames.size(); ++i)
    {
        table.AddEmptyRow();
        table.SetEntry("Variable", "nMinusOneBDTStudy.featureNames.at(" + std::to_string(i) + ")");
        table.SetEntry("Value", config.nMinusOneBDTStudy.featureNames.at(i));
    }

    table.AddEmptyRow();
    table.SetEntry("Variable", "nMinusOneBDTStudy.nSamplePoints");
    table.SetEntry("Value", config.nMinusOneBDTStudy.nSamplePoints);

    //
    // TrainBDTs
    //
    table.AddEmptyRow();

    table.AddEmptyRow();
    table.SetEntry("Variable", "trainBDTs.trainingFraction");
    table.SetEntry("Value", config.trainBDTs.trainingFraction);

    table.AddEmptyRow();
    table.SetEntry("Variable", "trainBDTs.onlyGoodTruthMatches");
    table.SetEntry("Value", config.trainBDTs.onlyGoodTruthMatches);

    table.AddEmptyRow();
    table.SetEntry("Variable", "trainBDTs.weightByCompleteness");
    table.SetEntry("Value", config.trainBDTs.weightByCompleteness);

    table.AddEmptyRow();
    table.SetEntry("Variable", "trainBDTs.shouldOptimize");
    table.SetEntry("Value", config.trainBDTs.shouldOptimize);

    table.AddEmptyRow();
    table.SetEntry("Variable", "trainBDTs.shouldMakePlots");
    table.SetEntry("Value", config.trainBDTs.shouldMakePlots);

    //
    // MakeEventSelectionTable
    //
    table.AddEmptyRow();

    table.AddEmptyRow();
    table.SetEntry("Variable", "makeEventSelectionTable.shouldOptimize");
    table.SetEntry("Value", config.makeEventSelectionTable.shouldOptimize);

    table.AddEmptyRow();
    table.SetEntry("Variable", "makeEventSelectionTable.nScanPoints");
    table.SetEntry("Value", config.makeEventSelectionTable.nScanPoints);

    table.AddEmptyRow();
    table.SetEntry("Variable", "makeEventSelectionTable.processFraction");
    table.SetEntry("Value", config.makeEventSelectionTable.processFraction);

    //
    // MakeSelectedPIDTable
    //
    table.AddEmptyRow();

    table.AddEmptyRow();
    table.SetEntry("Variable", "makeSelectedPIDTable.useGenericSelection");
    table.SetEntry("Value", config.makeSelectedPIDTable.useGenericSelection);

    table.AddEmptyRow();
    table.SetEntry("Variable", "makeSelectedPIDTable.goldenPionIsSignal");
    table.SetEntry("Value", config.makeSelectedPIDTable.goldenPionIsSignal);

    table.AddEmptyRow();
    table.SetEntry("Variable", "makeSelectedPIDTable.onlyLowMomentumPions");
    table.SetEntry("Value", config.makeSelectedPIDTable.onlyLowMomentumPions);

    table.AddEmptyRow();
    table.SetEntry("Variable", "makeSelectedPIDTable.pionMomentumThreshold");
    table.SetEntry("Value", config.makeSelectedPIDTable.pionMomentumThreshold);

    //
    // EfficiencyPlots
    //
    table.AddEmptyRow();

    table.AddEmptyRow();
    table.SetEntry("Variable", "efficiencyPlots.drawErrors");
    table.SetEntry("Value", config.efficiencyPlots.drawErrors);

    //
    // MakeBinningPlots
    //
    table.AddEmptyRow();

    table.AddEmptyRow();
    table.SetEntry("Variable", "makeBinningPlots.useFineBinEdges");
    table.SetEntry("Value", config.makeBinningPlots.useFineBinEdges);

    //
    // ExtractXSecs
    //
    table.AddEmptyRow();

    table.AddEmptyRow();
    table.SetEntry("Variable", "extractXSecs.nBootstrapUniverses");
    table.SetEntry("Value", config.extractXSecs.nBootstrapUniverses);

    for (unsigned int i = 0; i < config.extractXSecs.systematicParams.size(); ++i)
    {
        const auto &[paramName, nUniverses] = config.extractXSecs.systematicParams.at(i);
        table.AddEmptyRow();
        table.SetEntry("Variable", "config.extractXSecs.systematicParams.at(" + std::to_string(i) + ")");
        table.SetEntry("Value", paramName + " (" + std::to_string(nUniverses) + " universes)");
    }

    for (unsigned int i = 0; i < config.extractXSecs.mutuallyExclusiveSystematicParams.size(); ++i)
    {
        const auto &[combinedParamName, nUniverses, componentParams] = config.extractXSecs.mutuallyExclusiveSystematicParams.at(i);
        table.AddEmptyRow();
        table.SetEntry("Variable", "config.extractXSecs.mutuallyExclusiveSystematicParams.at(" + std::to_string(i) + ")");
        table.SetEntry("Value", combinedParamName + " (" + std::to_string(nUniverses) + " universes)");

        for (const auto &paramName : componentParams)
        {
            table.AddEmptyRow();
            table.SetEntry("Value", paramName);
        }
    }

    for (unsigned int i = 0; i < config.extractXSecs.fluxParams.size(); ++i)
    {
        const auto &paramName = config.extractXSecs.fluxParams.at(i);

        table.AddEmptyRow();
        table.SetEntry("Variable", "config.extractXSecs.fluxParams.at(" + std::to_string(i) + ")");
        table.SetEntry("Value", paramName);
    }

    for (unsigned int i = 0; i < config.extractXSecs.genieParams.size(); ++i)
    {
        const auto &paramName = config.extractXSecs.genieParams.at(i);

        table.AddEmptyRow();
        table.SetEntry("Variable", "config.extractXSecs.genieParams.at(" + std::to_string(i) + ")");
        table.SetEntry("Value", paramName);
    }

    // Print out the table
    table.Print();
}

} // namespace ubcc1pi_macros
