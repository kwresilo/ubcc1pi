#include "ubcc1pi_standalone/Objects/FileReader.h"

#include "ubcc1pi_standalone/Helpers/PlottingHelper.h"
#include "ubcc1pi_standalone/Helpers/AnalysisHelper.h"
#include "ubcc1pi_standalone/Helpers/SelectionHelper.h"
#include "ubcc1pi_standalone/Helpers/FormattingHelper.h"

using namespace ubcc1pi;

int MakeEventSelectionEfficiencyPlots(const std::string &overlayFileName, const bool drawErrors = false, const bool useAbsPdg = true)
{
    // Get the selection
    auto selection = SelectionHelper::GetDefaultSelection();
    const auto allCuts = selection.GetCuts();

    std::cout << "Making plots for cuts:" << std::endl;
    for (const auto &cut : allCuts)
        std::cout << " - " << cut << std::endl;

    // Read the input file
    FileReader reader(overlayFileName);
    auto pEvent = reader.GetBoundEventAddress();

    // Set up the plots
    auto plot_nuEnergy = PlottingHelper::EfficiencyPlot("True neutrino energy / GeV", 40u, 0.3f, 2.8f, allCuts, drawErrors);
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
        if (!AnalysisHelper::IsTrueCC1Pi(pEvent, useAbsPdg))
            continue;
        
        // Check which event selection cuts are passed by this event
        std::vector<std::string> cutsPassed;
        std::vector<int> assignedPdgCodes;
        const auto isSelected = selection.Execute(pEvent, cutsPassed, assignedPdgCodes);

        // Get the features we want to plot
        const auto nuEnergy = pEvent->truth.nuEnergy();

        const auto visibleParticles = AnalysisHelper::SelectVisibleParticles(pEvent->truth.particles);
        const auto nProtons = AnalysisHelper::CountParticlesWithPdgCode(visibleParticles, 2212, false);

        float piMomentum = -std::numeric_limits<float>::max();
        float piCosTheta = -std::numeric_limits<float>::max();
        float piPhi = -std::numeric_limits<float>::max();

        float muMomentum = -std::numeric_limits<float>::max();
        float muCosTheta = -std::numeric_limits<float>::max();
        float muPhi = -std::numeric_limits<float>::max();

        TVector3 muDir(0.f, 0.f, 0.f);
        TVector3 piDir(0.f, 0.f, 0.f);

        bool hasGoldenPion = false;

        // Get the particle features
        for (const auto &particle : pEvent->truth.particles)
        {
            const auto absPdg = std::abs(particle.pdgCode());
            const auto dir = TVector3(particle.momentumX(), particle.momentumY(), particle.momentumZ()).Unit();
            const auto cosTheta = dir.Z();
            const auto phi = std::atan2(dir.Y(), dir.X());

            if (absPdg == 211)
            {
                piMomentum = particle.momentum();
                piCosTheta = cosTheta;
                piPhi = phi;
                piDir = dir;

                hasGoldenPion = AnalysisHelper::IsGolden(particle);
            }

            if (absPdg == 13)
            {
                muMomentum = particle.momentum();
                muCosTheta = cosTheta;
                muPhi = phi;
                muDir = dir;
            }
        }

        const auto muPiAngle = std::acos(muDir.Dot(piDir));
        
        // Fill the plots at each step
        for (unsigned int i = 0; i < allCuts.size(); ++i)
        {
            const auto &cut = allCuts.at(i);
            const auto passedCut = (std::find(cutsPassed.begin(), cutsPassed.end(), cut) != cutsPassed.end());

            plot_nuEnergy.AddEvent(nuEnergy, cut, passedCut);
            plot_nProtons.AddEvent(nProtons, cut, passedCut);
            plot_muMomentum.AddEvent(muMomentum, cut, passedCut);
            plot_muCosTheta.AddEvent(muCosTheta, cut, passedCut);
            plot_muPhi.AddEvent(muPhi, cut, passedCut);
            plot_muPiAngle.AddEvent(muPiAngle, cut, passedCut);
            plot_piMomentum.AddEvent(piMomentum, cut, passedCut);
            plot_piCosTheta.AddEvent(piCosTheta, cut, passedCut);
            plot_piPhi.AddEvent(piPhi, cut, passedCut);

            if (!hasGoldenPion)
                continue;

            plot_piMomentumGolden.AddEvent(piMomentum, cut, passedCut);
            plot_piCosThetaGolden.AddEvent(piCosTheta, cut, passedCut);
            plot_piPhiGolden.AddEvent(piPhi, cut, passedCut);
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

    return 0;
}
