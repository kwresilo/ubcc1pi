#include "ubcc1pi_standalone/Macros/Macros.h"

#include "ubcc1pi_standalone/Objects/FileReader.h"

#include "ubcc1pi_standalone/Helpers/BDTHelper.h"
#include "ubcc1pi_standalone/Helpers/AnalysisHelper.h"
#include "ubcc1pi_standalone/Helpers/PlottingHelper.h"

using namespace ubcc1pi;

namespace ubcc1pi_macros
{

void TrainBDTs(const Config &config)
{
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Extract the CC1Pi events tha pass the pre-selection
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FileReader reader(config.files.overlaysFileName);
    auto pEvent = reader.GetBoundEventAddress();

    std::cout << "Finding CC1Pi events" << std::endl;
    std::vector<unsigned int> cc1PiEventIndices;

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

        cc1PiEventIndices.push_back(eventIndex);
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
    const auto featureNames = BDTHelper::ParticleBDTFeatureNames;
    BDTHelper::BDTFactory goldenPionBDTFactory("goldenPion", featureNames);
    BDTHelper::BDTFactory protonBDTFactory("proton", featureNames);

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Fill the BDT training and testing entries
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    std::cout << "Filling the BDT entries" << std::endl;
    for (unsigned int i = 0; i < nCC1PiEvents; ++i)
    {
        AnalysisHelper::PrintLoadingBar(i, nCC1PiEvents);

        const auto eventIndex = cc1PiEventIndices.at(i);
        const auto isTrainingEvent = shuffler.IsTrainingEvent(i);
        reader.LoadEvent(eventIndex);

        const auto truthParticles = pEvent->truth.particles;
        const auto recoParticles = pEvent->reco.particles;
        const float eventWeight = AnalysisHelper::GetNominalEventWeight(pEvent);

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
            std::vector<float> features;
            const auto areAllFeaturesAvailable = BDTHelper::GetBDTFeatures(recoParticle, featureNames, features);

            // Only use particles with all features available
            if (!areAllFeaturesAvailable)
                continue;

            // Define the weight
            const auto completenessWeight = (config.trainBDTs.weightByCompleteness ? (isExternal ? 1.f : completeness) : 1.f);
            const auto weight = eventWeight * completenessWeight;

            // Add the particle to the BDTs
            const bool isGoldenPion = !isExternal && truePdgCode == 211 && trueIsGolden;
            goldenPionBDTFactory.AddEntry(features, isGoldenPion, isTrainingEvent, weight);
            
            const bool isProton = !isExternal && truePdgCode == 2212;
            protonBDTFactory.AddEntry(features, isProton, isTrainingEvent, weight);
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
    }
    
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Train and test the BDTs
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    std::cout << "Training and testing golden pion BDT" << std::endl;
    goldenPionBDTFactory.TrainAndTest();

    std::cout << "Training and testing proton BDT" << std::endl;
    protonBDTFactory.TrainAndTest();

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Make plots
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if (!config.trainBDTs.shouldMakePlots)
        return;

    const std::string yLabel = "Fraction of reco particles";
    PlottingHelper::MultiPlot goldenPionBDTPlot("Golden pion BDT response", yLabel, 50, -0.9f, 0.45f);
    PlottingHelper::MultiPlot protonBDTPlot("Proton BDT response", yLabel, 50, -0.8f, 0.7f);

    // Using the newly trained BDT weight files, setup up a BDT for evaluation
    BDTHelper::BDT goldenPionBDT("goldenPion", featureNames); 
    BDTHelper::BDT protonBDT("proton", featureNames); 

    std::cout << "Making BDT training plots" << std::endl;
    for (unsigned int i = 0; i < nCC1PiEvents; ++i)
    {
        AnalysisHelper::PrintLoadingBar(i, nCC1PiEvents);

        const auto eventIndex = cc1PiEventIndices.at(i);
        const auto isTrainingEvent = shuffler.IsTrainingEvent(i);
        reader.LoadEvent(eventIndex);

        const auto truthParticles = pEvent->truth.particles;
        const auto recoParticles = pEvent->reco.particles;
        const float eventWeight = 1.f;

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
                truePdgCode = truthParticle.pdgCode();
                trueIsGolden = AnalysisHelper::IsGolden(truthParticle);
                completeness = recoParticle.truthMatchCompletenesses().at(truthParticleIndex);
            }
            catch (const std::exception &) {}

            // Only use good matches for testing
            if (config.trainBDTs.onlyGoodTruthMatches && (isExternal || completeness < 0.5f))
                continue;

            // Extract the features
            std::vector<float> features;
            const auto areAllFeaturesAvailable = BDTHelper::GetBDTFeatures(recoParticle, featureNames, features);

            // Only use particles with all features available
            if (!areAllFeaturesAvailable)
                continue;

            // Add the particle to the BDTs
            const auto goldenPionBDTResponse = goldenPionBDT.GetResponse(features);
            const auto protonBDTResponse = protonBDT.GetResponse(features);
           
            // Fill to the plots
            const auto style = PlottingHelper::GetPlotStyle(recoParticle, AnalysisHelper::Overlay, truthParticles, isTrainingEvent);

            // For these plots skip neutrons
            if (style == PlottingHelper::Other || style == PlottingHelper::OtherPoints)
                continue;

            goldenPionBDTPlot.Fill(goldenPionBDTResponse, style, eventWeight);
            protonBDTPlot.Fill(protonBDTResponse, style, eventWeight);

        }
    }
    
    // Save the plots
    goldenPionBDTPlot.SaveAs("goldenPionBDTResponse");
    protonBDTPlot.SaveAs("protonBDTResponse");
}

} // namespace ubcc1pi_plots
