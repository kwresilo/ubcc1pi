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
            plotVector.emplace_back("Wiggliness", yLabel, 60, 0, 0.004);
            plotVectorSignal.emplace_back("Wiggliness", yLabel, 60, 0, 0.004);
            continue;
        }

        if (featureName == "trackScore")
        {
            plotVector.emplace_back("Track score", yLabel, 60, 0, 1);
            plotVectorSignal.emplace_back("Track score", yLabel, 60, 0, 1);
            continue;
        }

        throw std::logic_error("PlotInputVariables - unknown feature: \"" + featureName + "\"");
    }

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

                bool isContained = false;
                try
                {
                    isContained = AnalysisHelper::IsContained(particle);
                }
                catch (const std::invalid_argument &) {}

                if (!isContained)
                    continue;

                // Get the plot style 
                const auto particleStyle = PlottingHelper::GetPlotStyle(particle, sampleType, truthParticles, false, config.global.useAbsPdg);

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
            }
        }
    }

    // Save the plots
    for (unsigned int iFeature = 0; iFeature < featureNames.size(); ++iFeature)
    {
        const auto &featureName = featureNames.at(iFeature);

        plotVector.at(iFeature).SaveAsStacked("inputVariables_" + featureName);
        plotVectorSignal.at(iFeature).SaveAs("inputVariables_signal_" + featureName);
    }
}

} // ubcc1pi macros
