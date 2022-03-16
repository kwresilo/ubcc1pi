/**
 *  @file  ubcc1pi_standalone/Macros/TrainBDTs.cxx
 *
 *  @brief The implementation file of the TrainBDTs macro
 */

#include "ubcc1pi_standalone/Macros/Macros.h"

#include "ubcc1pi_standalone/Objects/FileReader.h"

#include "ubcc1pi_standalone/Helpers/BDTHelper.h"
#include "ubcc1pi_standalone/Helpers/AnalysisHelper.h"
#include "ubcc1pi_standalone/Helpers/PlottingHelper.h"
#include "ubcc1pi_standalone/Helpers/NormalisationHelper.h"

using namespace ubcc1pi;

namespace ubcc1pi_macros
{

void TrainBDTs(const Config &config)
{
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Extract the CC1Pi events tha pass the pre-selection
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    std::cout << "Finding CC1Pi events" << std::endl;

    std::vector< std::tuple<AnalysisHelper::SampleType, unsigned int, std::string, float> > inputData;
    for (const auto run: config.global.runs)
    {
        if(run == 1)
        {
            inputData.emplace_back(AnalysisHelper::Overlay, 1, config.filesRun1.overlaysFileName, NormalisationHelper::GetOverlaysNormalisation(config, 1));
        }
        else if(run == 2)
        {
            inputData.emplace_back(AnalysisHelper::Overlay, 2, config.filesRun2.overlaysFileName, NormalisationHelper::GetOverlaysNormalisation(config, 2));
        }
        else if(run == 3)
        {
            inputData.emplace_back(AnalysisHelper::Overlay, 3, config.filesRun3.overlaysFileName, NormalisationHelper::GetOverlaysNormalisation(config, 3));
        }
        else throw std::logic_error("ExtractSidebandFit - Invalid run number");
    }

    std::vector<std::pair<unsigned int, unsigned int>> cc1PiEventIndices;
    for (const auto &[sampleType, sampleRun, fileName, normalisation] : inputData)
    {
        // Read the input file
        FileReader reader(fileName);
        auto pEvent = reader.GetBoundEventAddress();

        for (unsigned int eventIndex = 0, nEvents = reader.GetNumberOfEvents(); eventIndex < nEvents; ++eventIndex)
        {
            AnalysisHelper::PrintLoadingBar(eventIndex, nEvents);
            reader.LoadEvent(eventIndex);

            // Event must be true CC1Pi
            if (!AnalysisHelper::IsTrueCC1Pi(pEvent, config.global.useAbsPdg))
                continue;

            // Event must pass the CCInclusive selection
            if (!pEvent->reco.passesCCInclusive())
                continue;

            cc1PiEventIndices.emplace_back(sampleRun, eventIndex);
        }
    }
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Randomly choose the training events
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    const auto nCC1PiEvents = cc1PiEventIndices.size();
    const auto nTrainingEvents = static_cast<unsigned int>(std::floor(static_cast<float>(nCC1PiEvents) * config.trainBDTs.trainingFraction));
    BDTHelper::EventShuffler shuffler(nCC1PiEvents, nTrainingEvents);
    std::cout << "Found " << nCC1PiEvents << " CC1Pi events passing CC inclusive seleciton. Using " << nTrainingEvents << " for training." << std::endl;

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Setup the BDTs to train
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    const auto goldenPionFeatureNames = BDTHelper::GoldenPionBDTFeatureNames;
    const auto protonFeatureNames = BDTHelper::ProtonBDTFeatureNames;
    const auto muonFeatureNames = BDTHelper::MuonBDTFeatureNames;
    BDTHelper::BDTFactory goldenPionBDTFactory("goldenPion", goldenPionFeatureNames);
    BDTHelper::BDTFactory protonBDTFactory("proton", protonFeatureNames);
    BDTHelper::BDTFactory muonBDTFactory("muon", muonFeatureNames);

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Fill the BDT training and testing entries
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    for (const auto &[sampleType, sampleRun, fileName, normalisation] : inputData)
    {
        // Read the input file
        FileReader reader(fileName);
        auto pEvent = reader.GetBoundEventAddress();

        std::cout << "Filling the BDT entries" << std::endl;
        for (unsigned int i = 0; i < nCC1PiEvents; ++i)
        {
            const auto run = cc1PiEventIndices.at(i).first;
            if(sampleRun != run) continue;

            AnalysisHelper::PrintLoadingBar(i, nCC1PiEvents);

            const auto eventIndex = cc1PiEventIndices.at(i).second;
            const auto isTrainingEvent = shuffler.IsTrainingEvent(i);
            reader.LoadEvent(eventIndex);

            const auto truthParticles = pEvent->truth.particles;
            const auto recoParticles = pEvent->reco.particles;
            const float eventWeight = AnalysisHelper::GetNominalEventWeight(pEvent)*normalisation;

            for (const auto &recoParticle : recoParticles)
            {
                // Only use contained particles for training
                if (!AnalysisHelper::HasTrackFit(recoParticle) || !AnalysisHelper::IsContained(recoParticle))
                    continue;

                // Determine the true origin of the reco particle
                bool isExternal = true;
                int truePdgCode = -std::numeric_limits<int>::max();
                bool trueIsGolden = false;
                float trueMomentum = -std::numeric_limits<float>::max();
                float completeness = -std::numeric_limits<float>::max();

                try
                {
                    const auto truthParticleIndex = AnalysisHelper::GetBestMatchedTruthParticleIndex(recoParticle, truthParticles);
                    const auto truthParticle = truthParticles.at(truthParticleIndex);

                    isExternal = false;
                    truePdgCode = truthParticle.pdgCode();
                    trueMomentum = truthParticle.momentum();
                    trueIsGolden = AnalysisHelper::IsGolden(truthParticle);
                    completeness = recoParticle.truthMatchCompletenesses().at(truthParticleIndex);
                }
                catch (const std::exception &) {}

                // Only use good matches for training
                if (config.trainBDTs.onlyGoodTruthMatches && (isExternal || completeness < 0.5f))
                    continue;

                // Extract the features
                std::vector<float> goldenPionFeatures;
                const auto areAllFeaturesAvailableGoldenPion = BDTHelper::GetBDTFeatures(recoParticle, goldenPionFeatureNames, goldenPionFeatures);

                std::vector<float> protonFeatures;
                const auto areAllFeaturesAvailableProton = BDTHelper::GetBDTFeatures(recoParticle, protonFeatureNames, protonFeatures);

                std::vector<float> muonFeatures;
                const auto areAllFeaturesAvailableMuon = BDTHelper::GetBDTFeatures(recoParticle, muonFeatureNames, muonFeatures);

                // Define the weight
                const auto completenessWeight = (config.trainBDTs.weightByCompleteness ? (isExternal ? 1.f : completeness) : 1.f);
                const auto weight = eventWeight * completenessWeight;

                if (areAllFeaturesAvailableGoldenPion)
                {
                    const bool isGoldenPion = !isExternal && truePdgCode == 211 && trueIsGolden;
                    goldenPionBDTFactory.AddEntry(goldenPionFeatures, isGoldenPion, isTrainingEvent, weight);
                }

                if (areAllFeaturesAvailableProton)
                {
                    const bool isProton = !isExternal && truePdgCode == 2212;
                    protonBDTFactory.AddEntry(protonFeatures, isProton, isTrainingEvent, weight);
                }

                if (areAllFeaturesAvailableMuon)
                {
                    const bool isMuon = !isExternal && truePdgCode == 13;
                    muonBDTFactory.AddEntry(muonFeatures, isMuon, isTrainingEvent, weight);
                }
            }
        }
    }

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Optimize the BDTs
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if (config.trainBDTs.shouldOptimize)
    {
        std::cout << "Optimizing golden pion BDT" << std::endl;
        goldenPionBDTFactory.OptimizeParameters();

        std::cout << "Optimizing proton BDT" << std::endl;
        protonBDTFactory.OptimizeParameters();

        std::cout << "Optimizing muon BDT" << std::endl;
        muonBDTFactory.OptimizeParameters();
    }

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Train and test the BDTs
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    std::cout << "Training and testing golden pion BDT" << std::endl;
    goldenPionBDTFactory.TrainAndTest();

    std::cout << "Training and testing proton BDT" << std::endl;
    protonBDTFactory.TrainAndTest();

    std::cout << "Training and testing muon BDT" << std::endl;
    muonBDTFactory.TrainAndTest();

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Make plots
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if (!config.trainBDTs.shouldMakePlots)
        return;

    const std::string yLabel = "Fraction of reco particles";
    PlottingHelper::MultiPlot goldenPionBDTPlot("Golden pion BDT response", yLabel, 50, -0.9f, 0.45f);
    PlottingHelper::MultiPlot protonBDTPlot("Proton BDT response", yLabel, 50, -0.8f, 0.7f);
    PlottingHelper::MultiPlot muonBDTPlot("Muon BDT response", yLabel, 50, -0.9f, 0.6f);

    // Using the newly trained BDT weight files, setup up a BDT for evaluation
    BDTHelper::BDT goldenPionBDT("goldenPion", goldenPionFeatureNames);
    BDTHelper::BDT protonBDT("proton", protonFeatureNames);
    BDTHelper::BDT muonBDT("muon", muonFeatureNames);

    std::cout << "Making BDT training plots" << std::endl;
    for (const auto &[sampleType, sampleRun, fileName, normalisation] : inputData)
    {
        // Read the input file
        FileReader reader(fileName);
        auto pEvent = reader.GetBoundEventAddress();

        std::cout << "Filling the BDT entries" << std::endl;

        for (unsigned int i = 0; i < nCC1PiEvents; ++i)
        {
            const auto run = cc1PiEventIndices.at(i).first;
            if(sampleRun != run) continue;

            AnalysisHelper::PrintLoadingBar(i, nCC1PiEvents);
            const auto eventIndex = cc1PiEventIndices.at(i).second;
            const auto isTrainingEvent = shuffler.IsTrainingEvent(i);
            reader.LoadEvent(eventIndex);

            const auto truthParticles = pEvent->truth.particles;
            const auto recoParticles = pEvent->reco.particles;
            const float eventWeight = AnalysisHelper::GetNominalEventWeight(pEvent)*normalisation;

            for (const auto &recoParticle : recoParticles)
            {
                // Only use contained particles for testing
                if (!AnalysisHelper::HasTrackFit(recoParticle) || !AnalysisHelper::IsContained(recoParticle))
                    continue;

                // Determine the true origin of the reco particle
                bool isExternal = true;
                int truePdgCode = -std::numeric_limits<int>::max();
                bool trueIsGolden = false;
                float completeness = -std::numeric_limits<float>::max();

                try
                {
                    const auto truthParticleIndex = AnalysisHelper::GetBestMatchedTruthParticleIndex(recoParticle, truthParticles);
                    const auto truthParticle = truthParticles.at(truthParticleIndex);

                    isExternal = false;
                    truePdgCode = config.global.useAbsPdg ? std::abs(truthParticle.pdgCode()) : truthParticle.pdgCode();
                    trueIsGolden = AnalysisHelper::IsGolden(truthParticle);
                    completeness = recoParticle.truthMatchCompletenesses().at(truthParticleIndex);
                }
                catch (const std::exception &) {}

                // Only use good matches for testing
                if (config.trainBDTs.onlyGoodTruthMatches && (isExternal || completeness < 0.5f))
                    continue;

                // Fill to the plots
                const auto style = PlottingHelper::GetPlotStyle(recoParticle, AnalysisHelper::Overlay, truthParticles, isTrainingEvent, config.global.useAbsPdg);

                // For these plots skip neutrons
                if (style == PlottingHelper::Other || style == PlottingHelper::OtherPoints)
                    continue;

                // Extract the features
                std::vector<float> goldenPionFeatures;
                const auto areAllFeaturesAvailableGoldenPion = BDTHelper::GetBDTFeatures(recoParticle, goldenPionFeatureNames, goldenPionFeatures);

                std::vector<float> protonFeatures;
                const auto areAllFeaturesAvailableProton = BDTHelper::GetBDTFeatures(recoParticle, protonFeatureNames, protonFeatures);

                std::vector<float> muonFeatures;
                const auto areAllFeaturesAvailableMuon = BDTHelper::GetBDTFeatures(recoParticle, muonFeatureNames, muonFeatures);

                if (areAllFeaturesAvailableGoldenPion)
                {
                    const auto goldenPionBDTResponse = goldenPionBDT.GetResponse(goldenPionFeatures);
                    goldenPionBDTPlot.Fill(goldenPionBDTResponse, style, eventWeight);
                }

                if (areAllFeaturesAvailableProton)
                {
                    const auto protonBDTResponse = protonBDT.GetResponse(protonFeatures);
                    protonBDTPlot.Fill(protonBDTResponse, style, eventWeight);
                }

                if (areAllFeaturesAvailableMuon)
                {
                    const auto muonBDTResponse = muonBDT.GetResponse(muonFeatures);
                    muonBDTPlot.Fill(muonBDTResponse, style, eventWeight);
                }
            }
        }
    }

    // Save the plots
    goldenPionBDTPlot.SaveAs("goldenPionBDTResponse");
    protonBDTPlot.SaveAs("protonBDTResponse");
    muonBDTPlot.SaveAs("muonBDTResponse");
}

} // namespace ubcc1pi_plots
