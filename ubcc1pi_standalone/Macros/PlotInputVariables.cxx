#include "ubcc1pi_standalone/Objects/FileReader.h"

#include "ubcc1pi_standalone/Helpers/PlottingHelper.h"
#include "ubcc1pi_standalone/Helpers/AnalysisHelper.h"
#include "ubcc1pi_standalone/Helpers/SelectionHelper.h"
#include "ubcc1pi_standalone/Helpers/BDTHelper.h"

using namespace ubcc1pi;

int PlotInputVariables(const std::string &overlayFileName, const float overlayWeight, const std::string &dataEXTFileName, const float extWeight, const std::string &dataBNBFileName, const bool cleanSignalOnly = false, const bool useAbsPdg = true)
{
        
    const float bodgeFactor = 1.273f; // ATTN this factor is a normalisation added so we can compare the shape of the distributions, can't exist in the final result!
    
    // Get the selection
    auto selection = SelectionHelper::GetDefaultSelection();
    const auto allCuts = selection.GetCuts();

    // Set up the plots for each feature for every cut
    const auto featureNames = BDTHelper::ParticleBDTFeatureNames;
    const std::string yLabel = "Fraction of reco particles";
    std::vector< std::vector< PlottingHelper::MultiPlot > > plots; // Holds the input features plots, first index is cut second index is feature
    std::vector< PlottingHelper::MultiPlot > goldenPionBDTPlots;
    std::vector< PlottingHelper::MultiPlot > pionBDTPlots;
    std::vector< PlottingHelper::MultiPlot > protonBDTPlots;
    std::vector< PlottingHelper::MultiPlot > muonBDTPlots;
    std::vector< PlottingHelper::MultiPlot > goldenPionBDTOfPionPlots;

    for (unsigned int iCut = 0; iCut < allCuts.size(); ++iCut)
    {
        std::vector< PlottingHelper::MultiPlot > plotVector;

        for (const auto &featureName : featureNames)
        {
            if (featureName == "logBragg_pToMIP")
            {
                plotVector.emplace_back("log(L_p / L_MIP)", yLabel, 100, -8, 8);
                continue;
            }

            if (featureName == "logBragg_piToMIP")
            {
                plotVector.emplace_back("log(L_pi / L_MIP)", yLabel, 100, -4, 8);
                continue;
            }

            if (featureName == "truncMeandEdx")
            {
                plotVector.emplace_back("Truncated Mean dEdx", yLabel, 100, 0, 10);
                continue;
            }

            if (featureName == "protonForward")
            {
                plotVector.emplace_back("Proton forward likelihood", yLabel, 100, 0.3, 0.7);
                continue;
            }

            if (featureName == "muonForward")
            {
                plotVector.emplace_back("Muon forward likelihood", yLabel, 100, 0.3, 0.7);
                continue;
            }

            if (featureName == "nDescendents")
            {
                plotVector.emplace_back("Number of descendent particles", yLabel, 5, 0, 5);
                continue;
            }

            if (featureName == "nSpacePointsNearEnd")
            {
                plotVector.emplace_back("Number of spacepoints near track end", yLabel, 60, 0, 120);
                continue;
            }

            if (featureName == "wiggliness")
            {
                plotVector.emplace_back("Wiggliness", yLabel, 60, 0, 0.005);
                continue;
            }

            if (featureName == "trackScore")
            {
                plotVector.emplace_back("Track score", yLabel, 100, 0, 1);
                continue;
            }

            throw std::logic_error("PlotInputVariables - unknown feature: \"" + featureName + "\"");
        }

        plots.push_back(plotVector);
    
        goldenPionBDTPlots.emplace_back("Golden pion BDT response", yLabel, 50, -0.9f, 0.45f);
        pionBDTPlots.emplace_back("Pion BDT response", yLabel, 50, -0.8f, 0.7f);
        protonBDTPlots.emplace_back("Proton BDT response", yLabel, 50, -0.8f, 0.7f);
        muonBDTPlots.emplace_back("Muon BDT response", yLabel, 50, -0.9f, 0.7f);
    
        goldenPionBDTOfPionPlots.emplace_back("Golden pion BDT response of pion candidate", yLabel, 50, -0.9f, 0.45f);
    }
    
    // Set up the BDTs
    BDTHelper::BDT goldenPionBDT("goldenPion", featureNames); 
    BDTHelper::BDT pionBDT("pion", featureNames); 
    BDTHelper::BDT protonBDT("proton", featureNames); 
    BDTHelper::BDT muonBDT("muon", featureNames); 

    for (const auto fileName : {dataEXTFileName, dataBNBFileName, overlayFileName})
    {
        std::cout << "Reading input file: " << fileName << std::endl;

        const bool isBNBData = (fileName == dataBNBFileName);
        const bool isOverlay = (fileName == overlayFileName);
        const bool isEXTData = (fileName == dataEXTFileName);

        float weight = 1.f;
        if (isOverlay) weight = overlayWeight * bodgeFactor;
        if (isEXTData) weight = extWeight * bodgeFactor;

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
    
            if (cleanSignalOnly)
            {
                // Only use MC events 
                if (!pEvent->metadata.hasTruthInfo())
                    continue;

                // Only use signal events
                if (!AnalysisHelper::IsTrueCC1Pi(pEvent, useAbsPdg))
                    continue;
            }


            const auto truthParticles = pEvent->truth.particles;
            const auto recoParticles = pEvent->reco.particles;

            // Run the event selection and store which cuts are passed
            std::vector<std::string> cutsPassed;
            std::vector<int> assignedPdgCodes;
            const auto isSelected = selection.Execute(pEvent, cutsPassed, assignedPdgCodes);

            // Get the particles identified as the pion if it exists
            if (assignedPdgCodes.size() != recoParticles.size())
                throw std::logic_error("PlotInputVariables - The output particle PDG codes is the wrong size");

            unsigned int pionIndex = std::numeric_limits<unsigned int>::max();
            bool foundPion = false;
            for (unsigned int index = 0; index < assignedPdgCodes.size(); ++index)
            {
                const auto pdg = assignedPdgCodes.at(index);

                if (pdg == 211)
                {
                    if (foundPion)
                        throw std::logic_error("PlotInputVariables - Multiple pions found");

                    pionIndex = index;
                    foundPion = true;
                }
            }

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

                if (cleanSignalOnly)
                {
                    try
                    {
                        const auto truthParticleIndex = AnalysisHelper::GetBestMatchedTruthParticleIndex(particle, truthParticles);
                        const auto completeness = particle.truthMatchCompletenesses().at(truthParticleIndex);

                        if (completeness < 0.5f)
                            continue;
                    }
                    catch (const std::logic_error &)
                    {
                        continue;
                    }
                }

                const auto particleStyle = isBNBData ? PlottingHelper::BNBData : PlottingHelper::GetPlotStyle(particle, truthParticles);

                // Get the BDT features
                std::vector<float> features;
                if (!BDTHelper::GetBDTFeatures(particle, featureNames, features))
                    continue;
                   
                
                // Fill the plots
                for (unsigned int iCut = 0; iCut < allCuts.size(); ++iCut)
                {
                    // Only fill if we pass the cut
                    const auto &cut = allCuts.at(iCut);
                    if (std::find(cutsPassed.begin(), cutsPassed.end(), cut) == cutsPassed.end())
                        continue;

                    // Fill the feature plots
                    for (unsigned int iFeature = 0; iFeature < featureNames.size(); ++iFeature)
                        plots.at(iCut).at(iFeature).Fill(features.at(iFeature), particleStyle, weight);

                    // Fill the BDT response plots
                    const auto goldenPionBDTResponse = goldenPionBDT.GetResponse(features);
                    const auto pionBDTResponse = pionBDT.GetResponse(features);
                    const auto protonBDTResponse = protonBDT.GetResponse(features);
                    const auto muonBDTResponse = muonBDT.GetResponse(features);

                    goldenPionBDTPlots.at(iCut).Fill(goldenPionBDTResponse, particleStyle, weight);
                    pionBDTPlots.at(iCut).Fill(pionBDTResponse, particleStyle, weight);
                    protonBDTPlots.at(iCut).Fill(protonBDTResponse, particleStyle, weight);
                    muonBDTPlots.at(iCut).Fill(muonBDTResponse, particleStyle, weight);

                    if (foundPion && index == pionIndex)
                        goldenPionBDTOfPionPlots.at(iCut).Fill(goldenPionBDTResponse, particleStyle, weight);
                }
            }
        }
    }

    // Save the plots
    for (unsigned int iCut = 0; iCut < allCuts.size(); ++iCut)
    {
        const auto &cutName = allCuts.at(iCut);
        const auto suffix = "_" + std::to_string(iCut) + "_" + cutName + (cleanSignalOnly ? "_cleanSignal" : "");

        for (unsigned int iFeature = 0; iFeature < featureNames.size(); ++iFeature)
        {
            const auto &featureName = featureNames.at(iFeature);
            const auto name = featureName + suffix;
            plots.at(iCut).at(iFeature).SaveAs(name);

            if (!cleanSignalOnly)
                plots.at(iCut).at(iFeature).SaveAsStacked(name + "_stacked");
        }

        goldenPionBDTPlots.at(iCut).SaveAs("goldenPionBDTResponse" + suffix);
        pionBDTPlots.at(iCut).SaveAs("pionBDTResponse" + suffix);
        protonBDTPlots.at(iCut).SaveAs("protonBDTResponse" + suffix);
        muonBDTPlots.at(iCut).SaveAs("muonBDTResponse" + suffix);
    
        goldenPionBDTOfPionPlots.at(iCut).SaveAs("goldenPionBDTOfPionResponse" + suffix);
        
        if (!cleanSignalOnly)
        {
            goldenPionBDTPlots.at(iCut).SaveAsStacked("goldenPionBDTResponse" + suffix + "_stacked");
            pionBDTPlots.at(iCut).SaveAsStacked("pionBDTResponse" + suffix + "_stacked");
            protonBDTPlots.at(iCut).SaveAsStacked("protonBDTResponse" + suffix + "_stacked");
            muonBDTPlots.at(iCut).SaveAsStacked("muonBDTResponse" + suffix + "_stacked");
    
            goldenPionBDTOfPionPlots.at(iCut).SaveAsStacked("goldenPionBDTOfPionResponse" + suffix + "_stacked");
        }
    }

    return 0;
}
