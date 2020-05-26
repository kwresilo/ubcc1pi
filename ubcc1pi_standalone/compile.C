#include <cstdlib>
#include <stdexcept>
#include <string>

{
    // Add ubcc1pi_standalone to the include path
    const auto ubcc1piStandaloneDir = std::getenv("UBCC1PI_STANDALONE_DIR");
    if (!ubcc1piStandaloneDir)
        throw std::logic_error("The environment variable UBCC1PI_STANDALONE_DIR hasn't been set!");

    std::cout << "Adding " << ubcc1piStandaloneDir << " to include path" << std::endl;
    gSystem->AddIncludePath(("  -I" + std::string(ubcc1piStandaloneDir) + "/../").c_str());


    // Compile the code using root
    for (const auto &file : std::vector<std::string>({

        // Interface
        "Interface/Event.cxx",
        "Interface/Subrun.cxx",

        // Objects
        "Objects/FileReader.cxx",

        // Helpers
        "Helpers/AnalysisHelper.cxx",
        "Helpers/BDTHelper.cxx",
        "Helpers/PlottingHelper.cxx",
        "Helpers/SelectionHelper.cxx",
        "Helpers/NormalisationHelper.cxx",

        // Macros
        "Macros/Test.cxx",
        "Macros/CountPOT.cxx",
        "Macros/GetRunSubrunList.cxx",
        "Macros/PlotInputVariables.cxx",
        "Macros/CCInclusiveTruthStudy.cxx",
        "Macros/TrainBDTs.cxx",
        "Macros/NMinusOneBDTStudy.cxx",
        "Macros/MakeEventSelectionTable.cxx",
        "Macros/MakeEventSelectionEfficiencyPlots.cxx",
        "Macros/MakeSelectedPIDTable.cxx",
        "Macros/PlotReconstructedVariables.cxx",
        "Macros/ExtractXSecs.cxx"

    }))
    {
        std::cout << "Compiling " << file << std::endl;
        gROOT->ProcessLine((".L " + std::string(ubcc1piStandaloneDir) + "/" + file + "+").c_str());
    }

    // Setting the global plotting environment
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);

    gStyle->SetTitleFont(133, "XYZ");
    gStyle->SetTitleSize(18, "XYZ");
    gStyle->SetTitleFontSize(18);

    gStyle->SetLabelFont(133, "XYZ");
    gStyle->SetLabelSize(16, "XYZ");

    gROOT->ForceStyle();
}
