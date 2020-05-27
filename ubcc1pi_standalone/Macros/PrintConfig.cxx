#include "ubcc1pi_standalone/Macros/Macros.h"

#include "ubcc1pi_standalone/Objects/FileReader.h"

#include "ubcc1pi_standalone/Helpers/FormattingHelper.h"

using namespace ubcc1pi;

namespace ubcc1pi_macros
{

void PrintConfig(const Config &config)
{
    FormattingHelper::Table table({"Variable", "Value"});
    
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

    // Print out the table
    table.Print();
}

} // namespace ubcc1pi_macros
