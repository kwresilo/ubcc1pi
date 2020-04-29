#include "ubcc1pi_standalone/Objects/FileReader.h"

#include "ubcc1pi_standalone/Helpers/BDTHelper.h"
#include "ubcc1pi_standalone/Helpers/AnalysisHelper.h"
#include "ubcc1pi_standalone/Helpers/FormattingHelper.h"

using namespace ubcc1pi;

int NMinusOneBDTStudy(const std::string &overlayFileName, const bool useAbsPdg = true, const float trainingFraction = 0.5f, const bool onlyContained = true, const bool onlyGoodTruthMatches = false, const bool weightByCompleteness = true)
{
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Extract the CC1Pi events tha pass the pre-selection
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FileReader reader(overlayFileName);
    auto pEvent = reader.GetBoundEventAddress();

    std::cout << "Finding CC1Pi events" << std::endl;
    std::vector<unsigned int> cc1PiEventIndices;

    for (unsigned int eventIndex = 0, nEvents = reader.GetNumberOfEvents(); eventIndex < nEvents; ++eventIndex)
    {
        AnalysisHelper::PrintLoadingBar(eventIndex, nEvents);
        reader.LoadEvent(eventIndex);
        
        // Event must be true CC1Pi
        if (!AnalysisHelper::IsTrueCC1Pi(pEvent, useAbsPdg))
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
    const auto nTrainingEvents = static_cast<unsigned int>(std::floor(static_cast<float>(nCC1PiEvents) * trainingFraction));
    BDTHelper::EventShuffler shuffler(nCC1PiEvents, nTrainingEvents); 
    std::cout << "Found " << nCC1PiEvents << " CC1Pi events passing CC inclusive seleciton. Using " << nTrainingEvents << " for training." << std::endl;

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Setup the BDTs to train
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    const auto featureNames = BDTHelper::ParticleBDTFeatureNames;

    // Make copies of the list of features but remove one each time
    std::vector< std::vector<std::string> > reducedFeatureNames;
    std::vector< BDTHelper::BDTFactory > protonBDTFactoryVector;

    for (unsigned int i = 0; i < featureNames.size(); ++i)
    {
        std::vector<std::string> nameVector;
        for (unsigned int j = 0; j < featureNames.size(); ++j)
        {
            if (i == j) continue;

            nameVector.push_back(featureNames.at(j));
        }

        // Make the BDT and store the features
        protonBDTFactoryVector.emplace_back("proton_" + featureNames.at(i), nameVector);
        reducedFeatureNames.push_back(nameVector);
    }

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
        const float eventWeight = 1.f;

        for (const auto &recoParticle : recoParticles)
        {
            // Only use contained particles for training
            if (onlyContained && (!AnalysisHelper::HasTrackFit(recoParticle) || !AnalysisHelper::IsContained(recoParticle)))
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

            // Only use good matches for training
            if (onlyGoodTruthMatches && (isExternal || completeness < 0.5f))
                continue;

            // Define the weight
            const auto weight = eventWeight * (weightByCompleteness ? (isExternal ? 1.f : completeness) : 1.f);
            //const bool isGoldenPion = !isExternal && truePdgCode == 211 && trueIsGolden;
            //const bool isPion = !isExternal && truePdgCode == 211;
            const bool isProton = !isExternal && truePdgCode == 2212;
            //const bool isMuon = !isExternal && truePdgCode == 13;

            // Fill the BDTs 
            for (unsigned int iBDT = 0; iBDT < reducedFeatureNames.size(); ++iBDT)
            {
                const auto &usedFeatureNames = reducedFeatureNames.at(iBDT);
                auto &protonBDTFactory = protonBDTFactoryVector.at(iBDT);
            
                // Extract the features
                std::vector<float> features;
                const auto areAllFeaturesAvailable = BDTHelper::GetBDTFeatures(recoParticle, usedFeatureNames, features);
    
                // Only use particles with all features available
                if (!areAllFeaturesAvailable)
                    continue;
    
                // Add the particle to the BDTs
                protonBDTFactory.AddEntry(features, isProton, isTrainingEvent, weight);
            }
        }
    }
    
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Train and test the BDTs
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    std::cout << "Training and testing the BDTs" << std::endl;
    for (auto &bdtFactory : protonBDTFactoryVector)
    {
        bdtFactory.TrainAndTest();
    }
    
    
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Get the performance table
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    /*
    FormattingHelper::Table table({"Missing feature", "", });
    for (unsigned int iBDT = 0; iBDT < reducedFeatureNames.size(); ++iBDT)
    {
        // Get the trained BDT
        const auto missingFeature = featureNames.at(iBDT);
        const auto &usedFeatureNames = reducedFeatureNames.at(iBDT);
        BDTHelper::BDT protonBDT("proton_" + missingFeature, usedFeatureNames);


    }
    */
    return 0;
}
