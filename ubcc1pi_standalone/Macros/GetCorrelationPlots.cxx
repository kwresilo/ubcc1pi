/**
 *  @file  ubcc1pi_standalone/Macros/GetCorrelationPlots.cxx
 *
 *  @brief The implementation file of the GetCorrelationPlots macro
 */

#include "ubcc1pi_standalone/Macros/Macros.h"

#include "ubcc1pi_standalone/Objects/FileReader.h"

#include "ubcc1pi_standalone/Helpers/BDTHelper.h"
#include "ubcc1pi_standalone/Helpers/AnalysisHelper.h"
#include "ubcc1pi_standalone/Helpers/FormattingHelper.h"
#include "ubcc1pi_standalone/Helpers/PlottingHelper.h"

#include <TGraph.h>
#include <TH2F.h>
#include <TStyle.h>
#include <TVector3.h>

using namespace ubcc1pi;

namespace ubcc1pi_macros
{

void GetCorrelationPlots(const Config &config)
{
    // Setup the input file
    FileReader reader(config.files.overlaysFileName);
    auto pEvent = reader.GetBoundEventAddress();
    const auto nEvents = reader.GetNumberOfEvents();

    // Setup the data structures
    std::map<PlottingHelper::PlotStyle, std::map< std::string, std::vector<float> > > typeToFeatureNameToValuesMap;

    // Get the feature names
    const auto featureNames = BDTHelper::ParticleBDTFeatureNames;
    const auto nFeatures = featureNames.size();

    // Extract the values of the features
    for (unsigned int eventIndex = 0; eventIndex < nEvents; ++eventIndex)
    {
        AnalysisHelper::PrintLoadingBar(eventIndex, nEvents);
        reader.LoadEvent(eventIndex);

        // Event must be true CC1Pi
        if (!AnalysisHelper::IsTrueCC1Pi(pEvent, config.global.useAbsPdg))
            continue;

        // Event must pass the CCInclusive selection
        if (!pEvent->reco.passesCCInclusive())
            continue;

        for (const auto &particle : pEvent->reco.particles)
        {
            // Find the MC origin of this reco particle
            const auto style = PlottingHelper::GetPlotStyle(particle, AnalysisHelper::Overlay, pEvent->truth.particles, false, config.global.useAbsPdg);

            // Get the BDT features if available
            std::vector<float> features;
            if (!BDTHelper::GetBDTFeatures(particle, featureNames, features))
                continue;

            for (unsigned int iFeature = 0; iFeature < nFeatures; ++iFeature)
            {
                const auto featureName = featureNames.at(iFeature);
                const auto value = features.at(iFeature);

                // Store this value in the map
                typeToFeatureNameToValuesMap[style][featureName].push_back(value);
            }
        }
    }

    // Get the list of plot styles
    std::vector<PlottingHelper::PlotStyle> allStyles;
    for (const auto &entry : typeToFeatureNameToValuesMap)
        allStyles.push_back(entry.first);

    std::sort(allStyles.begin(), allStyles.end());

    // Setup the output correlation factor matricies
    std::map<PlottingHelper::PlotStyle, std::map<std::string, std::map< std::string, float > > > correlations;

    // Setup the canvas
    auto pCanvas = PlottingHelper::GetCanvas();

    // Now make the scatter plots for each pair of features
    for (unsigned int iFeature = 0; iFeature < nFeatures; ++iFeature)
    {
        const auto featureNameI = featureNames.at(iFeature);

        for (unsigned int jFeature = iFeature; jFeature < nFeatures; ++jFeature)
        {
            const auto featureNameJ = featureNames.at(jFeature);

            // Make a scatter graph for every style
            std::vector<TGraph *> graphVector;

            for (const auto &style : allStyles)
            {
                const auto &valuesI = typeToFeatureNameToValuesMap[style][featureNameI];
                const auto &valuesJ = typeToFeatureNameToValuesMap[style][featureNameJ];

                if (valuesI.size() != valuesJ.size())
                    throw std::logic_error("GetCorrelationPlots - found features with different number of entries!");

                const auto nValues = valuesI.size();

                // Make the graph
                TGraph *pGraph = new TGraph(nValues, valuesI.data(), valuesJ.data());
                pGraph->SetMarkerColor(PlottingHelper::GetColor(style));
                pGraph->Draw(style == allStyles.front() ? "AP" : "P");

                // Store the correlation
                const auto correlation = pGraph->GetCorrelationFactor();
                correlations[style][featureNameI][featureNameJ] = correlation;

                graphVector.push_back(pGraph);
            }

            // Save the graph
            PlottingHelper::SaveCanvas(pCanvas, "covariance_" + featureNameI + "-vs-" + featureNameJ);

            // Clean up the heap
            for (auto &pGraph : graphVector)
                delete pGraph;
        }
    }


    // Print the correlations
    FormattingHelper::Table table({"Particle type", "Feature I", "Feature J", "Correlation"});

    for (const auto &style : allStyles)
    {
        std::string particleType = "Other-" + std::to_string(style);
        switch (style)
        {
            case PlottingHelper::Muon:
                particleType = "muon";
                break;
            case PlottingHelper::Proton:
                particleType = "proton";
                break;
            case PlottingHelper::GoldenPion:
                particleType = "goldenPion";
                break;
            case PlottingHelper::NonGoldenPion:
                particleType = "nonGoldenPion";
                break;
            case PlottingHelper::External:
                particleType = "external";
                break;
            default: break;
        }

        // Set the number of decimal points to use when writing bin values
        gStyle->SetPaintTextFormat("2.2f");

        TH2F *pHist = new TH2F(("correlation_" + particleType).c_str(), "", nFeatures, 0, nFeatures, nFeatures, 0, nFeatures);
        pHist->SetMarkerSize(2.2); // Text size

        for (unsigned int iFeature = 0; iFeature < nFeatures; ++iFeature)
        {
            const auto featureNameI = featureNames.at(iFeature);

            // Set the bin labels (here we just use integers as a key, but could also use the actual feature names)
            const auto label = std::to_string(iFeature).c_str();
            pHist->GetXaxis()->SetBinLabel(iFeature + 1, label);
            pHist->GetYaxis()->SetBinLabel(iFeature + 1, label);

            for (unsigned int jFeature = iFeature; jFeature < nFeatures; ++jFeature)
            {
                const auto featureNameJ = featureNames.at(jFeature);

                const auto correlation = correlations.at(style).at(featureNameI).at(featureNameJ) * 100.f;

                table.AddEmptyRow();
                table.SetEntry("Particle type", particleType);
                table.SetEntry("Feature I", featureNameI);
                table.SetEntry("Feature J", featureNameJ);
                table.SetEntry("Correlation", correlation);

                // Set correlation plot bin values
                pHist->SetBinContent(iFeature + 1, jFeature + 1, correlation);
                pHist->SetBinContent(jFeature + 1, iFeature + 1, correlation); // Correlation is symmetric
            }
        }

        pHist->GetZaxis()->SetRangeUser(-100.f, 100.f);
        pHist->Draw("colz text");
        PlottingHelper::SaveCanvas(pCanvas, "correlationFactors_" + particleType);
    }

    table.WriteToFile("correlations.md");

}

} // namespace ubc1pi_macros
