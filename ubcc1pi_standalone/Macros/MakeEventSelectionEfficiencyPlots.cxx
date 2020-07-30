/**
 *  @file  ubcc1pi_standalone/Macros/MakeEventSelectionEfficiencyPlots.cxx
 *
 *  @brief The implementation file of the MakeEventSelectionEfficiencyPlots macro
 */

#include "ubcc1pi_standalone/Macros/Macros.h"

#include "ubcc1pi_standalone/Objects/FileReader.h"

#include "ubcc1pi_standalone/Helpers/PlottingHelper.h"
#include "ubcc1pi_standalone/Helpers/AnalysisHelper.h"
#include "ubcc1pi_standalone/Helpers/SelectionHelper.h"
#include "ubcc1pi_standalone/Helpers/FormattingHelper.h"

using namespace ubcc1pi;

namespace ubcc1pi_macros
{

void MakeEventSelectionEfficiencyPlots(const Config &config)
{
    // Get the selection
    auto selection = SelectionHelper::GetDefaultSelection();
    const auto allCuts = selection.GetCuts();
    
    //// BEGIN TEST
    selection.SetCutNominalValue("muonRange", -std::numeric_limits<float>::max());
    selection.SetCutNominalValue("pionRange", -std::numeric_limits<float>::max());
    //// END TEST
    
    std::cout << "Making plots for cuts:" << std::endl;
    for (const auto &cut : allCuts)
        std::cout << " - " << cut << std::endl;

    // Read the input file
    FileReader reader(config.files.overlaysFileName);
    auto pEvent = reader.GetBoundEventAddress();

    // Set up the plots
    const auto drawErrors = config.efficiencyPlots.drawErrors;
    auto plot_nuEnergy = PlottingHelper::EfficiencyPlot("Neutrino energy / GeV", 40u, 0, 2.5f, allCuts, drawErrors);
    auto plot_nProtons = PlottingHelper::EfficiencyPlot("Proton multiplicity", 5u, 0, 5, allCuts, drawErrors);
    auto plot_muMomentum = PlottingHelper::EfficiencyPlot("True muon momentum / GeV", 40u, 0.f, 1.5f, allCuts, drawErrors);
    auto plot_muCosTheta = PlottingHelper::EfficiencyPlot("True muon cos(theta)", 40u, -1.f, 1.0f, allCuts, drawErrors);
    auto plot_muPhi = PlottingHelper::EfficiencyPlot("True muon phi / rad", 40u, -3.142f, 3.142f, allCuts, drawErrors);
    auto plot_muPiAngle = PlottingHelper::EfficiencyPlot("True muon-pion opening angle", 40u, 0.f, 3.142f, allCuts, drawErrors);
    auto plot_piMomentum = PlottingHelper::EfficiencyPlot("True pion momentum / GeV", 40u, 0.f, 1.0f, allCuts, drawErrors);
    auto plot_piCosTheta = PlottingHelper::EfficiencyPlot("True pion cos(theta)", 40u, -1.f, 1.0f, allCuts, drawErrors);
    auto plot_piPhi = PlottingHelper::EfficiencyPlot("True pion phi / rad", 40u, -3.142, 3.142f, allCuts, drawErrors);
    auto plot_piMomentumGolden = PlottingHelper::EfficiencyPlot("True golden pion momentum / GeV", 40u, 0.f, 0.4f, allCuts, drawErrors);
    auto plot_piCosThetaGolden = PlottingHelper::EfficiencyPlot("True golden pion cos(theta)", 40u, -1.f, 1.0f, allCuts, drawErrors);
    auto plot_piPhiGolden = PlottingHelper::EfficiencyPlot("True golden pion phi / rad", 40u, -3.142, 3.142f, allCuts, drawErrors);

    // Run the selection
    const auto nEvents = reader.GetNumberOfEvents();
    for (unsigned int i = 0; i < nEvents; ++i)
    {
        AnalysisHelper::PrintLoadingBar(i, nEvents);
        reader.LoadEvent(i);

        // Only care about signal events
        if (!AnalysisHelper::IsTrueCC1Pi(pEvent, config.global.useAbsPdg))
            continue;
        
        // Check which event selection cuts are passed by this event
        std::vector<std::string> cutsPassed;
        std::vector<int> assignedPdgCodes;
        const auto isSelected = selection.Execute(pEvent, cutsPassed, assignedPdgCodes);

        // Get the features we want to plot
        const auto analysisData = AnalysisHelper::GetTruthAnalysisData(pEvent->truth, config.global.useAbsPdg, config.global.protonMomentumThreshold);
        
        // Fill the plots at each step
        for (unsigned int i = 0; i < allCuts.size(); ++i)
        {
            const auto &cut = allCuts.at(i);
            const auto passedCut = (std::find(cutsPassed.begin(), cutsPassed.end(), cut) != cutsPassed.end());

            plot_nuEnergy.AddEvent(pEvent->truth.nuEnergy(), cut, passedCut);
            plot_nProtons.AddEvent(analysisData.nProtons, cut, passedCut);
            plot_muMomentum.AddEvent(analysisData.muonMomentum, cut, passedCut);
            plot_muCosTheta.AddEvent(analysisData.muonCosTheta, cut, passedCut);
            plot_muPhi.AddEvent(analysisData.muonPhi, cut, passedCut);
            plot_muPiAngle.AddEvent(analysisData.muonPionAngle, cut, passedCut);
            plot_piMomentum.AddEvent(analysisData.pionMomentum, cut, passedCut);
            plot_piCosTheta.AddEvent(analysisData.pionCosTheta, cut, passedCut);
            plot_piPhi.AddEvent(analysisData.pionPhi, cut, passedCut);

            if (!analysisData.hasGoldenPion)
                continue;

            plot_piMomentumGolden.AddEvent(analysisData.pionMomentum, cut, passedCut);
            plot_piCosThetaGolden.AddEvent(analysisData.pionCosTheta, cut, passedCut);
            plot_piPhiGolden.AddEvent(analysisData.pionPhi, cut, passedCut);
        }
    }

    plot_nuEnergy.SaveAs("efficiency_nuEnergy");
    plot_nProtons.SaveAs("efficiency_nProtons");
    plot_muMomentum.SaveAs("efficiency_muMomentum");
    plot_muCosTheta.SaveAs("efficiency_muCosTheta");
    plot_muPhi.SaveAs("efficiency_muPhi");
    plot_muPiAngle.SaveAs("efficiency_muPiAngle");
    plot_piMomentum.SaveAs("efficiency_piMomentum");
    plot_piCosTheta.SaveAs("efficiency_piCosTheta");
    plot_piPhi.SaveAs("efficiency_piPhi");
    plot_piMomentumGolden.SaveAs("efficiency_piMomentumGolden");
    plot_piCosThetaGolden.SaveAs("efficiency_piCosThetaGolden");
    plot_piPhiGolden.SaveAs("efficiency_piPhiGolden");
}

} // namespace ubcc1pi_macros
