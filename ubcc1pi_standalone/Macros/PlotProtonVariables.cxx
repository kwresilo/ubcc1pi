/**
 *  @file  ubcc1pi_standalone/Macros/PlotProtonVariables.cxx
 *
 *  @brief The implementation file of the PlotProtonVariables macro
 */

#include "ubcc1pi_standalone/Macros/Macros.h"

#include "ubcc1pi_standalone/Objects/FileReader.h"

#include "ubcc1pi_standalone/Helpers/PlottingHelper.h"
#include "ubcc1pi_standalone/Helpers/AnalysisHelper.h"
#include "ubcc1pi_standalone/Helpers/SelectionHelper.h"
#include "ubcc1pi_standalone/Helpers/NormalisationHelper.h"


using namespace ubcc1pi;

namespace ubcc1pi_macros
{

void PlotProtonVariables(const Config &config)
{
    //
    // Setup the input files
    //
    std::cout<<"PlotProtonVariables Point -1"<<std::endl;
    std::vector< std::tuple<AnalysisHelper::SampleType, std::string, float> > inputData;

    inputData.emplace_back(AnalysisHelper::Overlay, config.files.overlaysFileName, NormalisationHelper::GetOverlaysNormalisation(config));
    inputData.emplace_back(AnalysisHelper::Dirt,    config.files.dirtFileName,     NormalisationHelper::GetDirtNormalisation(config));
    inputData.emplace_back(AnalysisHelper::DataEXT, config.files.dataEXTFileName,  NormalisationHelper::GetDataEXTNormalisation(config));
    inputData.emplace_back(AnalysisHelper::DataBNB, config.files.dataBNBFileName,  1.f);

    //
    // Setup the plots
    //
    std::cout<<"PlotProtonVariables Point 0"<<std::endl;
    const std::string yLabel = "Number of particles";
    /*
    PlottingHelper::MultiPlot muonMomentumPlot("Muon momentum / GeV", yLabel, config.global.muonMomentum.binEdges);
    PlottingHelper::MultiPlot muonCosThetaPlot("Muon cos(theta)", yLabel, config.global.muonCosTheta.binEdges);
    PlottingHelper::MultiPlot muonPhiPlot("Muon phi / rad", yLabel, config.global.muonPhi.binEdges);

    PlottingHelper::MultiPlot muonMomentumParticlePlot("Muon momentum / GeV", yLabel, config.global.muonMomentum.binEdges);
    PlottingHelper::MultiPlot muonCosThetaParticlePlot("Muon cos(theta)", yLabel, config.global.muonCosTheta.binEdges);
    PlottingHelper::MultiPlot muonPhiParticlePlot("Muon phi / rad", yLabel, config.global.muonPhi.binEdges);

    PlottingHelper::MultiPlot pionMomentumPlot("Pion momentum / GeV", yLabel, config.global.pionMomentum.binEdges);
    PlottingHelper::MultiPlot pionCosThetaPlot("Pion cos(theta)", yLabel, config.global.pionCosTheta.binEdges);
    PlottingHelper::MultiPlot pionPhiPlot("Pion phi / rad", yLabel, config.global.pionPhi.binEdges);

    PlottingHelper::MultiPlot pionMomentumParticlePlot("Pion momentum / GeV", yLabel, config.global.pionMomentum.binEdges);
    PlottingHelper::MultiPlot pionCosThetaParticlePlot("Pion cos(theta)", yLabel, config.global.pionCosTheta.binEdges);
    PlottingHelper::MultiPlot pionPhiParticlePlot("Pion phi / rad", yLabel, config.global.pionPhi.binEdges);

    PlottingHelper::MultiPlot muonPionAnglePlot("Muon-pion opening angle / rad", yLabel, config.global.muonPionAngle.binEdges);
    PlottingHelper::MultiPlot nProtonsPlot("Proton multiplicity", yLabel, config.global.nProtons.binEdges);
    */

    // PlottingHelper::MultiPlot muonMomentumPlot("Muon momentum / GeV", yLabel, 50u, 0.f, 2.f, true, config.global.axisTitles);
    // PlottingHelper::MultiPlot muonCosThetaPlot("Muon cos(theta)", yLabel, 50u, config.global.muonCosTheta.min, config.global.muonCosTheta.max, true, config.global.axisTitles);
    // PlottingHelper::MultiPlot muonPhiPlot("Muon phi / rad", yLabel, 50u, config.global.muonPhi.min, config.global.muonPhi.max, true, config.global.axisTitles);

    // PlottingHelper::MultiPlot muonMomentumParticlePlot("Muon momentum / GeV", yLabel, 50u, 0.f, 2.f, true, config.global.axisTitles);
    // PlottingHelper::MultiPlot muonCosThetaParticlePlot("Muon cos(theta)", yLabel, 50u, config.global.muonCosTheta.min, config.global.muonCosTheta.max, true, config.global.axisTitles);
    // PlottingHelper::MultiPlot muonPhiParticlePlot("Muon phi / rad", yLabel, 50u, config.global.muonPhi.min, config.global.muonPhi.max, true, config.global.axisTitles);

    // PlottingHelper::MultiPlot pionMomentumPlot("Pion momentum / GeV", yLabel, 50u, 0.f, 0.8f, true, config.global.axisTitles);
    // PlottingHelper::MultiPlot pionCosThetaPlot("Pion cos(theta)", yLabel, 50u, config.global.pionCosTheta.min, config.global.pionCosTheta.max, true, config.global.axisTitles);
    // PlottingHelper::MultiPlot pionPhiPlot("Pion phi / rad", yLabel, 50u, config.global.pionPhi.min, config.global.pionPhi.max, true, config.global.axisTitles);

    // PlottingHelper::MultiPlot pionMomentumParticlePlot("Pion momentum / GeV", yLabel, 50u, 0.f, 0.8f, true, config.global.axisTitles);
    // PlottingHelper::MultiPlot pionCosThetaParticlePlot("Pion cos(theta)", yLabel, 50u, config.global.pionCosTheta.min, config.global.pionCosTheta.max, true, config.global.axisTitles);
    // PlottingHelper::MultiPlot pionPhiParticlePlot("Pion phi / rad", yLabel, 50u, config.global.pionPhi.min, config.global.pionPhi.max, true, config.global.axisTitles);

    // PlottingHelper::MultiPlot muonPionAnglePlot("Muon-pion opening angle / rad", yLabel, 50u, config.global.muonPionAngle.min, config.global.muonPionAngle.max, true, config.global.axisTitles);
    // PlottingHelper::MultiPlot nProtonsPlot("Proton multiplicity", yLabel, 5u, 0, 5, true, config.global.axisTitles);

    // The highest energy proton variables are plotted with three different selections  
    const std::vector<std::string> protonPlotNames{"nProtons>0", "nProtons==1", "nProtons>=2"};
    std::vector<PlottingHelper::MultiPlot> protonMomentumPlots, protonCosThetaPlots, protonPhiPlots;
    std::vector<PlottingHelper::MultiPlot> protonMomentumParticlePlots, protonCosThetaParticlePlots, protonPhiParticlePlots;
    std::vector<PlottingHelper::MultiPlot> protonPionAnglePlots, protonMuonAnglePlots;
    
    // Proton plots are generated for three different muliplicity ranges
    // const std::vector<std::pair<std::string, std::function<bool(unsigned int)>> {
    //     {"protons>0", [](unsigned int n){return n>0;}}, 
    //     {"protons=1", [](unsigned int n){return n=1;}},
    //     {"protons>=2", [](unsigned int n){return n>=2;}}
    // };
    
    std::cout<<"PlotProtonVariables Point 1"<<std::endl;

    for(auto const& plotName : protonPlotNames)
    {
        //TODO: Some of the values here are placeholders - make sure they are sensible
        protonMomentumPlots.emplace_back("Proton momentum / GeV (" + plotName + ")", yLabel, 50u, config.global.protonMomentum.min, config.global.protonMomentum.max, true, config.global.axisTitles);
        protonCosThetaPlots.emplace_back("Proton cos(theta) (" + plotName + ")", yLabel, 50u, config.global.protonCosTheta.min, config.global.protonCosTheta.max, true, config.global.axisTitles);
        protonPhiPlots.emplace_back("Proton phi / rad (" + plotName + ")", yLabel, 50u, config.global.protonPhi.min, config.global.protonPhi.max, true, config.global.axisTitles);

        protonMomentumParticlePlots.emplace_back("Proton momentum / GeV (" + plotName + ")", yLabel, 50u, config.global.protonMomentum.min, config.global.protonMomentum.max, true, config.global.axisTitles);
        protonCosThetaParticlePlots.emplace_back("Proton cos(theta) (" + plotName + ")", yLabel, 50u, config.global.protonCosTheta.min, config.global.protonCosTheta.max, true, config.global.axisTitles);
        protonPhiParticlePlots.emplace_back("Proton phi / rad (" + plotName + ")", yLabel, 50u, config.global.protonPhi.min, config.global.protonPhi.max, true, config.global.axisTitles);

        protonPionAnglePlots.emplace_back("Proton-pion opening angle / rad (" + plotName + ")", yLabel, 50u, config.global.protonPionAngle.min, config.global.protonPionAngle.max, true, config.global.axisTitles);
        protonMuonAnglePlots.emplace_back("Proton-muon opening angle / rad (" + plotName + ")", yLabel, 50u, config.global.protonMuonAngle.min, config.global.protonMuonAngle.max, true, config.global.axisTitles);
    }

    //
    // Get the selection
    //
    auto selection = SelectionHelper::GetDefaultSelection();

    // Loop over the events
    for (const auto [sampleType, fileName, normalisation] : inputData)
    {
        std::cout << "Reading input file: " << fileName << std::endl;

        FileReader reader(fileName);
        auto pEvent = reader.GetBoundEventAddress();

        const auto nEvents = reader.GetNumberOfEvents();
        for (unsigned int i = 0; i < nEvents; ++i)
        {
            AnalysisHelper::PrintLoadingBar(i, nEvents);
            reader.LoadEvent(i);

            // Run the event selection and store which cuts are passed
            const auto &[isSelectedGolden, cutsPassed, assignedPdgCodes] = selection.Execute(pEvent);
            const auto isSelectedGeneric = (std::find(cutsPassed.begin(), cutsPassed.end(), config.global.lastCutGeneric) != cutsPassed.end());

            // Only use events that at least pass the generic selection
            if (!isSelectedGeneric)
                continue;

            // Get the truth and reco analysis data
            const auto weight = AnalysisHelper::GetNominalEventWeight(pEvent) * normalisation;
            const auto recoData = AnalysisHelper::GetRecoAnalysisData(pEvent->reco, assignedPdgCodes, isSelectedGolden);
            // const auto truthData = AnalysisHelper::GetTruthAnalysisData(pEvent->truth, config.global.useAbsPdg, config.global.protonMomentumThreshold);

            const auto plotStyle = PlottingHelper::GetPlotStyle(sampleType, pEvent, config.global.useAbsPdg);

            // Get the true origin of the selected muon and pion candidates
            const auto &recoParticles = pEvent->reco.particles;
            const auto &truthParticles = pEvent->truth.particles;


            // const auto muonPlotStyle = PlottingHelper::GetPlotStyle(muon, sampleType, truthParticles, false, config.global.useAbsPdg);
            // const auto pionPlotStyle = PlottingHelper::GetPlotStyle(pion, sampleType, truthParticles, false, config.global.useAbsPdg);

            // Fill the plots
            // muonMomentumPlot.Fill(recoData.muonMomentum, plotStyle, weight);
            // muonCosThetaPlot.Fill(recoData.muonCosTheta, plotStyle, weight);
            // muonPhiPlot.Fill(recoData.muonPhi, plotStyle, weight);

            // muonMomentumParticlePlot.Fill(recoData.muonMomentum, muonPlotStyle, weight);
            // muonCosThetaParticlePlot.Fill(recoData.muonCosTheta, muonPlotStyle, weight);
            // muonPhiParticlePlot.Fill(recoData.muonPhi, muonPlotStyle, weight);

            // if (recoData.hasGoldenPion)
            // {
            //     pionMomentumPlot.Fill(recoData.pionMomentum, plotStyle, weight);
            //     pionMomentumParticlePlot.Fill(recoData.pionMomentum, pionPlotStyle, weight);
            // }

            // pionCosThetaPlot.Fill(recoData.pionCosTheta, plotStyle, weight);
            // pionPhiPlot.Fill(recoData.pionPhi, plotStyle, weight);

            // pionCosThetaParticlePlot.Fill(recoData.pionCosTheta, pionPlotStyle, weight);
            // pionPhiParticlePlot.Fill(recoData.pionPhi, pionPlotStyle, weight);

            // muonPionAnglePlot.Fill(recoData.muonPionAngle, plotStyle, weight);
            // nProtonsPlot.Fill(recoData.nProtons, plotStyle, weight);

            // Not all selected events have protons. Use only those that do to plot reconstructed proton variables.

            // Only get the reconstructed proton variables when a suitable candidate is present
            const auto protonIndex = SelectionHelper::GetLeadingProtonCandidateIndex(recoParticles, assignedPdgCodes);
            if (protonIndex!=std::numeric_limits<unsigned int>::max())
            {
                const auto &muon = recoParticles.at(AnalysisHelper::GetParticleIndexWithPdg(assignedPdgCodes, 13));
                const auto &pion = recoParticles.at(AnalysisHelper::GetParticleIndexWithPdg(assignedPdgCodes, 211));
                const auto &proton = recoParticles.at(protonIndex);
                const auto muonDir = TVector3(muon.directionX(), muon.directionY(), muon.directionZ()).Unit();
                const auto pionDir = TVector3(pion.directionX(), pion.directionY(), pion.directionZ()).Unit();
                const auto protonDir = TVector3(proton.directionX(), proton.directionY(), proton.directionZ()).Unit();
                
                const auto protonMomentum = AnalysisHelper::GetProtonMomentumFromRange(proton.range());
                const auto protonCosTheta = protonDir.Z();
                const auto protonPhi  = std::atan2(protonDir.Y(), protonDir.X());
                const auto protonPionAngle = std::acos(protonDir.Dot(pionDir));
                const auto protonMuonAngle = std::acos(protonDir.Dot(muonDir));

                const auto protonPlotStyle = PlottingHelper::GetPlotStyle(proton, sampleType, truthParticles, false, config.global.useAbsPdg);
                for (unsigned int plot = 0; plot < protonPlotNames.size(); ++plot)
                {
                    // Check proton multiplicity to match conditions described in protonPlotNames
                    // No need to to test for plot==0 due to: if (protonIndex!=std::numeric_limits<unsigned int>::max()) 
                    if(plot==1 && recoData.nProtons!=1) continue;
                    if(plot==2 && recoData.nProtons<2) continue;

                    protonMomentumPlots.at(plot).Fill(protonMomentum, plotStyle, weight);
                    protonCosThetaPlots.at(plot).Fill(protonCosTheta, plotStyle, weight);
                    protonPhiPlots.at(plot).Fill(protonPhi, plotStyle, weight);

                    protonMomentumParticlePlots.at(plot).Fill(protonMomentum, protonPlotStyle, weight);
                    protonCosThetaParticlePlots.at(plot).Fill(protonCosTheta, protonPlotStyle, weight);
                    protonPhiParticlePlots.at(plot).Fill(protonPhi, protonPlotStyle, weight);

                    protonPionAnglePlots.at(plot).Fill(protonPionAngle, plotStyle, weight);
                    protonMuonAnglePlots.at(plot).Fill(protonMuonAngle, plotStyle, weight);
                }
            }
        }
    }

    // muonMomentumPlot.SaveAsStacked("reco_muonMomentum",false,false,config.global.axisTitles);
    // muonCosThetaPlot.SaveAsStacked("reco_muonCosTheta",false,false,config.global.axisTitles);
    // muonPhiPlot.SaveAsStacked("reco_muonPhi",false,false,config.global.axisTitles);

    // muonMomentumParticlePlot.SaveAsStacked("reco_muonMomentum_particle",false,false,config.global.axisTitles);
    // muonCosThetaParticlePlot.SaveAsStacked("reco_muonCosTheta_particle",false,false,config.global.axisTitles);
    // muonPhiParticlePlot.SaveAsStacked("reco_muonPhi_particle",false,false,config.global.axisTitles);

    // pionMomentumPlot.SaveAsStacked("reco_pionMomentum",false,false,config.global.axisTitles);
    // pionCosThetaPlot.SaveAsStacked("reco_pionCosTheta",false,false,config.global.axisTitles);
    // pionPhiPlot.SaveAsStacked("reco_pionPhi",false,false,config.global.axisTitles);

    // pionMomentumParticlePlot.SaveAsStacked("reco_pionMomentum_particle",false,false,config.global.axisTitles);
    // pionCosThetaParticlePlot.SaveAsStacked("reco_pionCosTheta_particle",false,false,config.global.axisTitles);
    // pionPhiParticlePlot.SaveAsStacked("reco_pionPhi_particle",false,false,config.global.axisTitles);

    // muonPionAnglePlot.SaveAsStacked("reco_muonPionAngle",false,false,config.global.axisTitles);
    // nProtonsPlot.SaveAsStacked("reco_nProtons",false,false,config.global.axisTitles);

    for (unsigned int plot = 0; plot < protonPlotNames.size(); ++plot)
    {
        const auto plotName = protonPlotNames.at(plot);

        protonMomentumPlots.at(plot).SaveAsStacked("reco_protonMomentum_" + plotName,false,false,false,config.global.axisTitles);
        protonCosThetaPlots.at(plot).SaveAsStacked("reco_protonCosTheta_" + plotName,false,false,false,config.global.axisTitles);
        protonPhiPlots.at(plot).SaveAsStacked("reco_protonPhi_" + plotName,false,false,false,config.global.axisTitles);

        protonMomentumParticlePlots.at(plot).SaveAsStacked("reco_protonMomentum_particle_" + plotName,false,false,false,config.global.axisTitles);
        protonCosThetaParticlePlots.at(plot).SaveAsStacked("reco_protonCosTheta_particle_" + plotName,false,false,false,config.global.axisTitles);
        protonPhiParticlePlots.at(plot).SaveAsStacked("reco_protonPhi_particle_" + plotName,false,false,false,config.global.axisTitles);

        protonPionAnglePlots.at(plot).SaveAsStacked("reco_protonPionAngle_" + plotName,false,false,false,config.global.axisTitles);
        protonMuonAnglePlots.at(plot).SaveAsStacked("reco_protonMuonAngle_" + plotName,false,false,false,config.global.axisTitles);
    }
}

} // namespace ubcc1pi_macros
