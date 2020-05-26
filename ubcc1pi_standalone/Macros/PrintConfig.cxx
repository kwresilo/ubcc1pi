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

    table.Print();
}

}
