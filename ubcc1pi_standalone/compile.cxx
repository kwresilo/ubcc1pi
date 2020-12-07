#include <cstdlib>
#include <stdexcept>
#include <string>

#include <vector>
#include <map>

{
    // Generate any required dictionaries
    gInterpreter->GenerateDictionary("std::map<std::string, std::vector<float> >", "map;string;vector");

    // Add ubcc1pi_standalone to the include path
    const auto ubcc1piStandaloneDir = std::getenv("UBCC1PI_STANDALONE_DIR");
    if (!ubcc1piStandaloneDir)
        throw std::logic_error("The environment variable UBCC1PI_STANDALONE_DIR hasn't been set!");

    std::cout << "Adding " << ubcc1piStandaloneDir << " to include path" << std::endl;
    gSystem->AddIncludePath(("  -I" + std::string(ubcc1piStandaloneDir) + "/../").c_str());

    std::cout << "Adding " << ubcc1piStandaloneDir << "/ubsmear/inc/ to include path" << std::endl;
    gSystem->AddIncludePath(("  -I" + std::string(ubcc1piStandaloneDir) + "/ubsmear/inc/").c_str());


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
        "Helpers/CrossSectionHelper.cxx",

        // Macros
        "Macros/PrintConfig.cxx",
        "Macros/CountPOT.cxx",
        "Macros/GetRunSubrunList.cxx",
        "Macros/TruthStudy.cxx",
        "Macros/FitRangeCurves.cxx",
        "Macros/SecondaryInteractionsStudy.cxx",
        "Macros/MomentumThresholdsStudy.cxx",
        "Macros/CCInclusiveMuonPIDStudy.cxx",
        "Macros/MultiPlanePIDDemo.cxx",
        "Macros/PlotInputVariables.cxx",
        "Macros/GetCorrelationPlots.cxx",
        "Macros/TrainBDTs.cxx",
        "Macros/NMinusOneBDTStudy.cxx",
        "Macros/MuonPIDStudy.cxx",
        "Macros/MakeEventSelectionTable.cxx",
        "Macros/MakeEventSelectionEfficiencyPlots.cxx",
        "Macros/PlotEventSelectionCuts.cxx",
        "Macros/MakeSelectedPIDTable.cxx",
        "Macros/PlotMuonRecoVariables.cxx",
        "Macros/PlotReconstructedVariables.cxx",
        "Macros/MakeBinningPlots.cxx",
        "Macros/PlotFluxVariations.cxx",
        "Macros/ExtractXSecs.cxx",

        "Macros/RunFullAnalysis.cxx"
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
