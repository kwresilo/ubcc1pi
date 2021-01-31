/**
 *  @file  ubcc1pi_standalone/Macros/PlotInputVariables.cxx
 *
 *  @brief The implementation file of the PlotInputVariables macro
 */

#include "ubcc1pi_standalone/Macros/Macros.h"

#include "ubcc1pi_standalone/Objects/FileReader.h"

#include "ubcc1pi_standalone/Helpers/PlottingHelper.h"
#include "ubcc1pi_standalone/Helpers/BDTHelper.h"
#include "ubcc1pi_standalone/Helpers/AnalysisHelper.h"
#include "ubcc1pi_standalone/Helpers/NormalisationHelper.h"

#include <TH2F.h>

using namespace ubcc1pi;

namespace ubcc1pi_macros
{

void PlotInputVariables(const Config &config)
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
    // Set up the plots for each BDT feature
    //
    const auto featureNames = BDTHelper::ParticleBDTFeatureNames;
    const std::string yLabel = "Fraction of reco particles";
    std::vector< PlottingHelper::MultiPlot > plotVector, plotVectorSignal;

    for (const auto &featureName : featureNames)
    {
        if (featureName == "logBragg_pToMIP")
        {
            plotVector.emplace_back("log(L_p / L_MIP)", yLabel, 60, -9, 7);
            plotVectorSignal.emplace_back("log(L_p / L_MIP)", yLabel, 60, -9, 7);
            continue;
        }

        if (featureName == "logBragg_piToMIP")
        {
            plotVector.emplace_back("log(L_pi / L_MIP)", yLabel, 60, -4, 6);
            plotVectorSignal.emplace_back("log(L_pi / L_MIP)", yLabel, 60, -4, 6);
            continue;
        }

        if (featureName == "truncMeandEdx")
        {
            plotVector.emplace_back("Truncated Mean dEdx", yLabel, 60, 0, 10);
            plotVectorSignal.emplace_back("Truncated Mean dEdx", yLabel, 60, 0, 10);
            continue;
        }

        if (featureName == "protonForward")
        {
            plotVector.emplace_back("Proton forward likelihood", yLabel, 60, 0.42, 0.62);
            plotVectorSignal.emplace_back("Proton forward likelihood", yLabel, 60, 0.42, 0.62);
            continue;
        }

        if (featureName == "muonForward")
        {
            plotVector.emplace_back("Muon forward likelihood", yLabel, 60, 0.35, 0.65);
            plotVectorSignal.emplace_back("Muon forward likelihood", yLabel, 60, 0.35, 0.65);
            continue;
        }

        if (featureName == "nDescendents")
        {
            plotVector.emplace_back("Number of descendent particles", yLabel, 4, 0, 4);
            plotVectorSignal.emplace_back("Number of descendent particles", yLabel, 4, 0, 4);
            continue;
        }

        if (featureName == "nSpacePointsNearEnd")
        {
            plotVector.emplace_back("Number of spacepoints near track end", yLabel, 45, 0, 90);
            plotVectorSignal.emplace_back("Number of spacepoints near track end", yLabel, 45, 0, 90);
            continue;
        }

        if (featureName == "wiggliness")
        {
            // ATTN here we exclude zero to better show the distribution, but zero is used by the BDT
            plotVector.emplace_back("Wiggliness", yLabel, 30, 0, 0.004);
            plotVectorSignal.emplace_back("Wiggliness", yLabel, 30, 0, 0.004);
            continue;
        }

        if (featureName == "trackScore")
        {
            plotVector.emplace_back("Track score", yLabel, 30, 0, 1);
            plotVectorSignal.emplace_back("Track score", yLabel, 30, 0, 1);
            continue;
        }

        throw std::logic_error("PlotInputVariables - unknown feature: \"" + featureName + "\"");
    }

    // Setup the BDT outputs
    const auto goldenPionFeatureNames = BDTHelper::GoldenPionBDTFeatureNames;
    const auto protonFeatureNames = BDTHelper::ProtonBDTFeatureNames;
    const auto muonFeatureNames = BDTHelper::MuonBDTFeatureNames;

    PlottingHelper::MultiPlot muonBDTPlot("Muon BDT response", yLabel, 40, -0.85f, 0.50f);
    PlottingHelper::MultiPlot protonBDTPlot("Proton BDT response", yLabel, 40, -0.60f, 0.60f);
    PlottingHelper::MultiPlot goldenPionBDTPlot("Golden pion BDT response", yLabel, 40, -0.8f, 0.4f);

    std::shared_ptr<BDTHelper::BDT> pGoldenPionBDT, pProtonBDT, pMuonBDT;
    if (config.plotInputVariables.plotBDTResponses)
    {
        // Setup the BDTs
        pGoldenPionBDT = std::make_shared<BDTHelper::BDT>("goldenPion", goldenPionFeatureNames);
        pProtonBDT = std::make_shared<BDTHelper::BDT>("proton", protonFeatureNames);
        pMuonBDT = std::make_shared<BDTHelper::BDT>("muon", muonFeatureNames);
    }

    PlottingHelper::MultiPlot phiPlot("Phi / rad", yLabel, 50u, -3.142f, 3.142f);
    PlottingHelper::MultiPlot cosThetaPlot("cos(theta)", yLabel, 50u, -1.f, 1.f);

    TH2F *hPhiCosThetaData = new TH2F("hPhiCosThetaData", "", 100u, -3.142f, 3.142f, 100u, -1.f, 1.f);
    TH2F *hPhiCosThetaSim = new TH2F("hPhiCosThetaSim", "", 100u, -3.142f, 3.142f, 100u, -1.f, 1.f);

    //
    // Fill the plots
    //
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

            // Only use events passing the CC inclusive selection
            if (!pEvent->reco.passesCCInclusive())
                continue;

            const auto weight = normalisation * AnalysisHelper::GetNominalEventWeight(pEvent);
            const auto recoParticles = pEvent->reco.particles;

            const auto truthParticles = pEvent->truth.particles; // This will be empty for non MC events
            const auto isSignal = (sampleType == AnalysisHelper::Overlay && AnalysisHelper::IsTrueCC1Pi(pEvent, config.global.useAbsPdg));

            for (unsigned int index = 0; index < recoParticles.size(); ++index)
            {
                const auto &particle = recoParticles.at(index);

                // Get the plot style
                const auto particleStyle = PlottingHelper::GetPlotStyle(particle, sampleType, truthParticles, false, config.global.useAbsPdg);

                // Insist the particle has a fitted track
                if (!AnalysisHelper::HasTrackFit(particle))
                    continue;

                // Fill the angle plots
                const auto dir = TVector3(particle.directionX(), particle.directionY(), particle.directionZ()).Unit();
                const auto phi = std::atan2(dir.Y(), dir.X());
                const auto cosTheta = dir.Z();
                phiPlot.Fill(phi, particleStyle, weight);
                cosThetaPlot.Fill(cosTheta, particleStyle, weight);

                if (sampleType == AnalysisHelper::DataBNB)
                {
                    hPhiCosThetaData->Fill(phi, cosTheta, weight);
                }
                else
                {
                    hPhiCosThetaSim->Fill(phi, cosTheta, weight);
                }

                // Insist the particle is contained
                if (!AnalysisHelper::IsContained(particle))
                    continue;

                // Get the BDT features
                std::vector<float> features;
                if (!BDTHelper::GetBDTFeatures(particle, featureNames, features))
                    continue;

                // Fill the feature plots
                for (unsigned int iFeature = 0; iFeature < featureNames.size(); ++iFeature)
                {
                    plotVector.at(iFeature).Fill(features.at(iFeature), particleStyle, weight);

                    // ATTN there can be particles that are in a signal event (i.e. neutrons) but we don't want to plot, so here we check
                    // for the "Other" category to avoid including them
                    if (isSignal && particleStyle != PlottingHelper::Other)
                    {
                        if (particleStyle != PlottingHelper::Muon &&
                            particleStyle != PlottingHelper::Proton &&
                            particleStyle != PlottingHelper::NonGoldenPion &&
                            particleStyle != PlottingHelper::GoldenPion &&
                            particleStyle != PlottingHelper::External)
                        {
                            throw std::logic_error("Found signal event with reco particle matching to unexpected truth particle");
                        }


                        plotVectorSignal.at(iFeature).Fill(features.at(iFeature), particleStyle, weight);
                    }
                }

                // Fill the BDT plots
                if (config.plotInputVariables.plotBDTResponses)
                {
                    std::vector<float> goldenPionFeatures;
                    const auto areAllFeaturesAvailableGoldenPion = BDTHelper::GetBDTFeatures(particle, goldenPionFeatureNames, goldenPionFeatures);

                    std::vector<float> protonFeatures;
                    const auto areAllFeaturesAvailableProton = BDTHelper::GetBDTFeatures(particle, protonFeatureNames, protonFeatures);

                    std::vector<float> muonFeatures;
                    const auto areAllFeaturesAvailableMuon = BDTHelper::GetBDTFeatures(particle, muonFeatureNames, muonFeatures);

                    if (areAllFeaturesAvailableGoldenPion)
                    {
                        const auto goldenPionBDTResponse = pGoldenPionBDT->GetResponse(goldenPionFeatures);
                        goldenPionBDTPlot.Fill(goldenPionBDTResponse, particleStyle, weight);
                    }

                    if (areAllFeaturesAvailableProton)
                    {
                        const auto protonBDTResponse = pProtonBDT->GetResponse(protonFeatures);
                        protonBDTPlot.Fill(protonBDTResponse, particleStyle, weight);
                    }

                    if (areAllFeaturesAvailableMuon)
                    {
                        const auto muonBDTResponse = pMuonBDT->GetResponse(muonFeatures);
                        muonBDTPlot.Fill(muonBDTResponse, particleStyle, weight);
                    }
                }
            }
        }
    }

    // Save the plots
    for (unsigned int iFeature = 0; iFeature < featureNames.size(); ++iFeature)
    {
        const auto &featureName = featureNames.at(iFeature);

        const bool useLogY = (featureName == "trackScore" || featureName == "wiggliness");

        plotVector.at(iFeature).SaveAsStacked("inputVariables_" + featureName, false, false, useLogY);
        plotVectorSignal.at(iFeature).SaveAs("inputVariables_signal_" + featureName, false, false, 0, useLogY);
    }

    goldenPionBDTPlot.SaveAsStacked("inputVariables_goldenPionBDTResponse");
    protonBDTPlot.SaveAsStacked("inputVariables_protonBDTResponse");
    muonBDTPlot.SaveAsStacked("inputVariables_muonBDTResponse");

    phiPlot.SaveAsStacked("inputVariables_phi");
    cosThetaPlot.SaveAsStacked("inputVariables_cosTheta");

    phiPlot.SaveAs("inputVariables_phi_unstacked", false, false, 500u);
    cosThetaPlot.SaveAs("inputVariables_cosTheta_unstacked", false, false, 500u);

    auto pCanvas = PlottingHelper::GetCanvas();
    hPhiCosThetaSim->Draw("colz");
    PlottingHelper::SaveCanvas(pCanvas, "inputVariables_phi-cosTheta");

    hPhiCosThetaData->Draw("colz");
    PlottingHelper::SaveCanvas(pCanvas, "inputVariables_phi-cosTheta_data");
}

} // ubcc1pi macros
