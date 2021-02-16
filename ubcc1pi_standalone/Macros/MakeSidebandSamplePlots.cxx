/**
 *  @file  ubcc1pi_standalone/Macros/PlotMuonRecoVariables.cxx
 *
 *  @brief The implementation file of the PlotMuonRecoVariables macro
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

void MakeSidebandSamplePlots(const Config &config)
{
    // Setup the input files
    std::vector< std::tuple<AnalysisHelper::SampleType, std::string, float> > inputData;

    inputData.emplace_back(AnalysisHelper::Overlay, config.files.overlaysFileName, NormalisationHelper::GetOverlaysNormalisation(config));
    inputData.emplace_back(AnalysisHelper::Dirt,    config.files.dirtFileName,     NormalisationHelper::GetDirtNormalisation(config));
    inputData.emplace_back(AnalysisHelper::DataEXT, config.files.dataEXTFileName,  NormalisationHelper::GetDataEXTNormalisation(config));
    inputData.emplace_back(AnalysisHelper::DataBNB, config.files.dataBNBFileName,  1.f);

    // Get the selection
    auto selection = SelectionHelper::GetSelection("CC0pi");
    const auto cuts = selection.GetCuts();

    // Get the muon BDT
    const auto muonFeatureNames = BDTHelper::MuonBDTFeatureNames;
    BDTHelper::BDT muonBDT("muon", muonFeatureNames);

    //
    // Setup the plots
    //
    const std::string yLabel = "Number of particles (norm. to bin width)";
    std::vector<PlottingHelper::MultiPlot> muonMomentumPlots, muonCosThetaPlots, muonPhiPlots;
    std::vector<PlottingHelper::MultiPlot> protonMomentumPlots, protonCosThetaPlots, protonPhiPlots;
    for (const auto &cut : cuts)
    {
        muonMomentumPlots.emplace_back("Muon momentum / GeV", yLabel, config.global.muonMomentum.binEdges, true, config.global.axisTitles);
        muonCosThetaPlots.emplace_back("Muon cos(theta)", yLabel, config.global.muonCosTheta.binEdges, true, config.global.axisTitles);
        muonPhiPlots.emplace_back("Muon phi / rad", yLabel, config.global.muonPhi.binEdges, true, config.global.axisTitles);

        protonMomentumPlots.emplace_back("Proton momentum / GeV", yLabel, config.global.pionMomentum.binEdges, true, config.global.axisTitles);
        protonCosThetaPlots.emplace_back("Proton cos(theta)", yLabel, config.global.pionCosTheta.binEdges, true, config.global.axisTitles);
        protonPhiPlots.emplace_back("Proton phi / rad", yLabel, config.global.pionPhi.binEdges, true, config.global.axisTitles);
    }

    // Loop over the events
    for (const auto [sampleType, fileName, normalisation] : inputData)
    {
        std::cout << "Reading input file: " << fileName << std::endl;

        FileReader reader(fileName);
        auto pEvent = reader.GetBoundEventAddress();

        const auto nEvents = reader.GetNumberOfEvents();
        for (unsigned int i = 0; i < /*nEvents*/100; ++i)
        {
            AnalysisHelper::PrintLoadingBar(i, nEvents);
            reader.LoadEvent(i);

            if (!pEvent->reco.passesCCInclusive())
                continue;

            const auto recoParticles = pEvent->reco.particles;

            // Run the event selection and store which cuts are passed
            const auto &[isSelected, cutsPassed, assignedPdgCodes] = selection.Execute(pEvent);

            // Get the truth and reco analysis data
            const auto weight = AnalysisHelper::GetNominalEventWeight(pEvent) * normalisation;
            const auto plotStyle = PlottingHelper::GetPlotStyle(sampleType, pEvent, config.global.useAbsPdg);

            // Get the muon reco data
            const auto muonIndex = SelectionHelper::GetMuonCandidateIndex(recoParticles, muonFeatureNames, muonBDT);
            const auto muon = recoParticles.at(muonIndex);
            const auto muonDir = TVector3(muon.directionX(), muon.directionY(), muon.directionZ()).Unit();
            const float muonCosTheta = muonDir.Z();
            const float muonPhi = std::atan2(muonDir.Y(), muonDir.X());
            const float muonMomentum = AnalysisHelper::GetMuonMomentum(muon);

            // Get the proton reco data
            const auto protonIndex = SelectionHelper::GetLeadingProtonCandidateIndex(recoParticles, muonFeatureNames, muonBDT);
            const auto proton = recoParticles.at(protonIndex);
            const auto protonDir = TVector3(proton.directionX(), proton.directionY(), proton.directionZ()).Unit();
            const float protonCosTheta = protonDir.Z();
            const float protonPhi = std::atan2(protonDir.Y(), protonDir.X());
            const float protonMomentum = AnalysisHelper::GetProtonMomentum(proton);

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
        muonMomentumPlots.at(iCut).SaveAsStacked("reco_cc1pi_muonMomentum_" + suffix, false, config.global.scaleByBinWidth, config.global.axisTitles);
        muonCosThetaPlots.at(iCut).SaveAsStacked("reco_cc1pi_muonCosTheta_" + suffix, false, config.global.scaleByBinWidth, config.global.axisTitles);
        muonPhiPlots.at(iCut).SaveAsStacked("reco_cc1pi_muonPhi_" + suffix, false, config.global.scaleByBinWidth, config.global.axisTitles);
    }
}

} // namespace ubcc1pi_macros
