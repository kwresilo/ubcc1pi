/**
 *  @file  ubcc1pi_standalone/Macros/PlotCCInclusiveMuonRecoVariables.cxx
 *
 *  @brief The implementation file of the PlotCCInclusiveMuonRecoVariables macro
 */

#include "ubcc1pi_standalone/Macros/Macros.h"

#include "ubcc1pi_standalone/Objects/FileReader.h"

#include "ubcc1pi_standalone/Helpers/PlottingHelper.h"
#include "ubcc1pi_standalone/Helpers/AnalysisHelper.h"
#include "ubcc1pi_standalone/Helpers/SelectionHelper.h"
#include "ubcc1pi_standalone/Helpers/NormalisationHelper.h"

#include <TH2F.h>

using namespace ubcc1pi;

namespace ubcc1pi_macros
{

void PlotCCInclusiveMuonRecoVariables(const Config &config)
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
    // Get the selection
    //
    auto selection = SelectionHelper::GetDefaultSelection();
    const auto cuts = selection.GetCuts();


    //
    // Setup the plots
    //

    const std::string yLabel = "Number of particles";
    std::vector<PlottingHelper::MultiPlot> muonMomentumPlots, muonCosThetaPlots, muonPhiPlots;
    for (const auto &cut : cuts)
    {    
        muonMomentumPlots.emplace_back("Muon momentum / GeV", yLabel, config.global.muonMomentum.binEdges);
        muonCosThetaPlots.emplace_back("Muon cos(theta)", yLabel, config.global.muonCosTheta.binEdges);
        muonPhiPlots.emplace_back("Muon phi / rad", yLabel, config.global.muonPhi.binEdges);
    }

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
            
            const auto recoParticles = pEvent->reco.particles;

            // Run the event selection and store which cuts are passed
            std::vector<std::string> cutsPassed;
            std::vector<int> assignedPdgCodes;
            const auto isSelected = selection.Execute(pEvent, cutsPassed, assignedPdgCodes);

            // Get the truth and reco analysis data
            const auto weight = AnalysisHelper::GetNominalEventWeight(pEvent) * normalisation;
            const auto plotStyle = PlottingHelper::GetPlotStyle(sampleType, pEvent, config.global.useAbsPdg);
    
            // Get the muon reco data
            bool foundRecoMuon = false;
            unsigned int muonIndex = std::numeric_limits<unsigned int>::max();
            for (unsigned int index = 0; index < recoParticles.size(); ++index)
            {
                if (recoParticles.at(index).isCCInclusiveMuonCandidate())
                {
                    if (foundRecoMuon)
                        throw std::logic_error("Found multiple muon candidates");

                    muonIndex = index;
                    foundRecoMuon = true;
                }
            }
            
            if (!foundRecoMuon)
                continue;

            const auto muon = recoParticles.at(muonIndex);
            const auto muonDir = TVector3(muon.directionX(), muon.directionY(), muon.directionZ()).Unit();
            const float muonCosTheta = muonDir.Z();
            const float muonPhi = std::atan2(muonDir.Y(), muonDir.X());
            const float muonMomentum = AnalysisHelper::GetMuonMomentum(muon);
               
            for (unsigned int iCut = 0; iCut < cuts.size(); ++iCut)
            {
                const auto passesCut = (std::find(cutsPassed.begin(), cutsPassed.end(), cuts.at(iCut)) != cutsPassed.end());

                if (!passesCut)
                    continue;

                muonMomentumPlots.at(iCut).Fill(muonMomentum, plotStyle, weight);
                muonCosThetaPlots.at(iCut).Fill(muonCosTheta, plotStyle, weight);
                muonPhiPlots.at(iCut).Fill(muonPhi, plotStyle, weight);
            }
        }
    }
    
    for (unsigned int iCut = 0; iCut < cuts.size(); ++iCut)
    {
        const auto &cut = cuts.at(iCut);

        const std::string suffix = std::to_string(iCut) + "_" + cut;
        muonMomentumPlots.at(iCut).SaveAsStacked("reco_ccinclusive_muonMomentum_" + suffix);
        muonCosThetaPlots.at(iCut).SaveAsStacked("reco_ccinclusive_muonCosTheta_" + suffix);
        muonPhiPlots.at(iCut).SaveAsStacked("reco_ccinclusive_muonPhi_" + suffix);
    }
}

} // namespace ubcc1pi_macros
