/**
 *  @file  ubcc1pi_standalone/Macros/NMinusOneBDTDataMCStudy.cxx
 *
 *  @brief The implementation file of the NMinusOneBDTDataMCStudy macro
 */

#include "ubcc1pi_standalone/Macros/Macros.h"

#include "ubcc1pi_standalone/Objects/FileReader.h"

#include "ubcc1pi_standalone/Helpers/BDTHelper.h"
#include "ubcc1pi_standalone/Helpers/AnalysisHelper.h"
#include "ubcc1pi_standalone/Helpers/NormalisationHelper.h"

using namespace ubcc1pi;

namespace ubcc1pi_macros
{

void NMinusOneBDTDataMCStudy(const Config &config)
{
    //
    // Setup the input files
    //
    std::vector< std::tuple<AnalysisHelper::SampleType, std::string, float> > inputData;

    inputData.emplace_back(AnalysisHelper::Overlay, config.files.overlaysFileName, NormalisationHelper::GetOverlaysNormalisation(config));
    inputData.emplace_back(AnalysisHelper::Dirt,    config.files.dirtFileName,     NormalisationHelper::GetDirtNormalisation(config));
    inputData.emplace_back(AnalysisHelper::DataEXT, config.files.dataEXTFileName,  NormalisationHelper::GetDataEXTNormalisation(config));
    inputData.emplace_back(AnalysisHelper::DataBNB, config.files.dataBNBFileName,  1.f);


    // Define some strings that match the naming structure of BDT weights files produced using NMinusOneBDTStudy
    std::string nFeaturesString = std::to_string(config.nMinusOneBDTStudy.featureNames.size());
    for (const auto &feature : config.nMinusOneBDTStudy.featureNames)
    {
        nFeaturesString += "_" + feature;
    }

    std::string signalString = "";
    switch (config.nMinusOneBDTStudy.signalType)
    {
        case PlottingHelper::Muon:
            signalString = "muon";
            break;
        case PlottingHelper::Proton:
            signalString = "proton";
            break;
        case PlottingHelper::GoldenPion:
            signalString = "goldenPion";
            break;
        default:
            throw std::invalid_argument("NMinusOneBDTDataMCStudy - Input signal type must be: proton, muon or golden pion");
    }

    // Load the BDTs with each feature removed (must be pre-trained)
    std::vector< BDTHelper::BDT > bdtVector;
    std::vector< string> bdtNameVector;
    std::vector< std::vector<string> > bdtUsedFeatureNames;
    std::vector< PlottingHelper::MultiPlot > bdtPlotVector;

    for (unsigned int i = 0; i < config.nMinusOneBDTStudy.featureNames.size(); ++i)
    {
        const auto missingFeature = config.nMinusOneBDTStudy.featureNames.at(i);

        std::vector<std::string> usedFeatureNames;
        for (unsigned int j = 0; j < config.nMinusOneBDTStudy.featureNames.size(); ++j)
        {
            if (i == j) continue;

            usedFeatureNames.push_back(config.nMinusOneBDTStudy.featureNames.at(j));
        }

        const auto bdtName = signalString + "_N-" + nFeaturesString + "_minus_" + missingFeature;
        bdtNameVector.push_back(bdtName);
        bdtVector.emplace_back(bdtName, usedFeatureNames);
        bdtUsedFeatureNames.push_back(usedFeatureNames);

        // Setup the BDT plots
        switch (config.nMinusOneBDTStudy.signalType)
        {
            case PlottingHelper::Muon:
                bdtPlotVector.emplace_back("Muon BDT response", "Number of reco particles", 40, -0.85f, 0.50f);
                break;
            case PlottingHelper::Proton:
                bdtPlotVector.emplace_back("Proton BDT response", "Number of reco particles", 40, -0.60f, 0.60f);
                break;
            case PlottingHelper::GoldenPion:
                bdtPlotVector.emplace_back("Golden pion BDT response", "Number of reco particles", 40, -0.8f, 0.4f);
                break;
            default:
                throw std::invalid_argument("NMinusOneBDTDataMCStudy - Input signal type must be: proton, muon or golden pion");
        }
    }

    // Loop over all events and produce data/MC plots for each BDT response
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

                // Insist the particle is contained
                if (!AnalysisHelper::IsContained(particle))
                    continue;

                // Loop over the BDTs
                for (unsigned int iBDT = 0; iBDT < bdtVector.size(); ++iBDT)
                {
                    // Get the BDT features
                    std::vector<float> features;
                    if (!BDTHelper::GetBDTFeatures(particle, bdtUsedFeatureNames.at(iBDT), features))
                        continue;

                    // Get the BDT response
                    const auto bdtResponse = bdtVector.at(iBDT).GetResponse(features);

                    // Fill the plots
                    bdtPlotVector.at(iBDT).Fill(bdtResponse, particleStyle, weight);
                }
            }
        }
    }

    // Save the plots
    for (unsigned int iBDT = 0; iBDT < bdtVector.size(); ++iBDT)
    {
        bdtPlotVector.at(iBDT).SaveAsStacked("nMinusOneBDTDataMCStudy_" + bdtNameVector.at(iBDT));
    }
}

} // namespace ubc1pi_macros
