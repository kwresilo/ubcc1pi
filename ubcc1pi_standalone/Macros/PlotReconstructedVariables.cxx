/**
 *  @file  ubcc1pi_standalone/Macros/PlotReconstructedVariables.cxx
 *
 *  @brief The implementation file of the PlotReconstructedVariables macro
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

void PlotReconstructedVariables(const Config &config)
{
    //
    // Setup the input files
    // 
    std::vector< std::tuple<AnalysisHelper::SampleType, std::string, float> > inputData;
    
    inputData.emplace_back(AnalysisHelper::Overlay, config.files.overlaysFileName, NormalisationHelper::GetOverlaysNormalisation(config)); 
    inputData.emplace_back(AnalysisHelper::Dirt,    config.files.dirtFileName,     NormalisationHelper::GetDirtNormalisation(config)); 
    inputData.emplace_back(AnalysisHelper::DataEXT, config.files.dataEXTFileName,  NormalisationHelper::GetDataEXTNormalisation(config)); 
    inputData.emplace_back(AnalysisHelper::DataBNB, config.files.dataBNBFileName,  1.f);

    //
    // Setup the plots
    //

    const std::string yLabel = "Number of particles";
    PlottingHelper::MultiPlot muonMomentumPlot("Muon momentum / GeV", yLabel, config.global.muonMomentum.binEdges);
    PlottingHelper::MultiPlot muonCosThetaPlot("Muon cos(theta)", yLabel, config.global.muonCosTheta.binEdges);
    PlottingHelper::MultiPlot muonPhiPlot("Muon phi / rad", yLabel, config.global.muonPhi.binEdges);

    PlottingHelper::MultiPlot pionMomentumPlot("Pion momentum / GeV", yLabel, config.global.pionMomentum.binEdges);
    PlottingHelper::MultiPlot pionCosThetaPlot("Pion cos(theta)", yLabel, config.global.pionCosTheta.binEdges);
    PlottingHelper::MultiPlot pionPhiPlot("Pion phi / rad", yLabel, config.global.pionPhi.binEdges);
    
    PlottingHelper::MultiPlot muonPionAnglePlot("Muon-pion opening angle / rad", yLabel, config.global.muonPionAngle.binEdges);
    PlottingHelper::MultiPlot nProtonsPlot("Proton multiplicity", yLabel, config.global.nProtons.binEdges);

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
            std::vector<std::string> cutsPassed;
            std::vector<int> assignedPdgCodes;
            const auto isSelectedGolden = selection.Execute(pEvent, cutsPassed, assignedPdgCodes);
            const auto isSelectedGeneric = (std::find(cutsPassed.begin(), cutsPassed.end(), config.global.lastCutGeneric) != cutsPassed.end());

            // Only use events that at least pass the generic selection
            if (!isSelectedGeneric)
                continue;
            
            // Get the truth and reco analysis data
            const auto weight = AnalysisHelper::GetNominalEventWeight(pEvent) * normalisation;
            const auto recoData = AnalysisHelper::GetRecoAnalysisData(pEvent->reco, assignedPdgCodes, isSelectedGolden);
    
            const auto plotStyle = PlottingHelper::GetPlotStyle(sampleType, pEvent, config.global.useAbsPdg);

            muonMomentumPlot.Fill(recoData.muonMomentum, plotStyle, weight);
            muonCosThetaPlot.Fill(recoData.muonCosTheta, plotStyle, weight);
            muonPhiPlot.Fill(recoData.muonPhi, plotStyle, weight);

            if (recoData.hasGoldenPion)
                pionMomentumPlot.Fill(recoData.pionMomentum, plotStyle, weight);

            pionCosThetaPlot.Fill(recoData.pionCosTheta, plotStyle, weight);
            pionPhiPlot.Fill(recoData.pionPhi, plotStyle, weight);
    
            muonPionAnglePlot.Fill(recoData.muonPionAngle, plotStyle, weight);
            nProtonsPlot.Fill(recoData.nProtons, plotStyle, weight);
        }
    }
            
    muonMomentumPlot.SaveAsStacked("reco_muonMomentum");
    muonCosThetaPlot.SaveAsStacked("reco_muonCosTheta");
    muonPhiPlot.SaveAsStacked("reco_muonPhi");

    pionMomentumPlot.SaveAsStacked("reco_pionMomentum");
    pionCosThetaPlot.SaveAsStacked("reco_pionCosTheta");
    pionPhiPlot.SaveAsStacked("reco_pionPhi");

    muonPionAnglePlot.SaveAsStacked("reco_muonPionAngle");
    nProtonsPlot.SaveAsStacked("reco_nProtons");
}

} // namespace ubcc1pi_macros
