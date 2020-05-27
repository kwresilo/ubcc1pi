#include "ubcc1pi_standalone/Macros/Macros.h"

#include "ubcc1pi_standalone/Objects/FileReader.h"

#include "ubcc1pi_standalone/Helpers/BDTHelper.h"
#include "ubcc1pi_standalone/Helpers/AnalysisHelper.h"
#include "ubcc1pi_standalone/Helpers/FormattingHelper.h"

using namespace ubcc1pi;

namespace ubcc1pi_macros
{

// TODO validate this whole macro before showing the results
void NMinusOneBDTStudy(const Config &config)
{
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Extract the CC1Pi events that pass the pre-selection
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

    // Make copies of the list of features but remove one each time
    std::vector< std::vector<std::string> > reducedFeatureNames;
    std::vector< BDTHelper::BDTFactory > protonBDTFactoryVector;
    std::vector< BDTHelper::BDTFactory > goldenPionBDTFactoryVector;

    for (unsigned int i = 0; i < featureNames.size(); ++i)
    {
        std::vector<std::string> nameVector;
        for (unsigned int j = 0; j < featureNames.size(); ++j)
        {
            if (i == j) continue;

            nameVector.push_back(featureNames.at(j));
        }

        // Make the BDT and store the features
        protonBDTFactoryVector.emplace_back("proton_minus_" + featureNames.at(i), nameVector);
        goldenPionBDTFactoryVector.emplace_back("goldenPion_minus_" + featureNames.at(i), nameVector);
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
            if (config.trainBDTs.onlyGoodTruthMatches && (isExternal || completeness < 0.5f))
                continue;

            // Define the weight
            const auto weight = eventWeight * (config.trainBDTs.weightByCompleteness ? (isExternal ? 1.f : completeness) : 1.f);
            const bool isGoldenPion = !isExternal && truePdgCode == 211 && trueIsGolden;
            const bool isProton = !isExternal && truePdgCode == 2212;

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
                protonBDTFactoryVector.at(iBDT).AddEntry(features, isProton, isTrainingEvent, weight);
                goldenPionBDTFactoryVector.at(iBDT).AddEntry(features, isGoldenPion, isTrainingEvent, weight);
            }
        }
    }
    
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Train and test the BDTs
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    std::cout << "Training and testing the proton BDTs" << std::endl;
    for (auto &bdtFactory : protonBDTFactoryVector)
        bdtFactory.TrainAndTest();
    
    std::cout << "Training and testing the golden pion BDTs" << std::endl;
    for (auto &bdtFactory : goldenPionBDTFactoryVector)
        bdtFactory.TrainAndTest();
    
    
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Get the performance table
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    std::cout << "Comparing performances of BDTs" << std::endl;

    // Get the trained BDTs
    std::vector< BDTHelper::BDT > protonBDTVector, goldenPionBDTVector;
    for (unsigned int iBDT = 0; iBDT < reducedFeatureNames.size(); ++iBDT)
    {
        // Get the trained BDT
        const auto missingFeature = featureNames.at(iBDT);

        const auto &usedFeatureNames = reducedFeatureNames.at(iBDT);
        protonBDTVector.emplace_back("proton_minus_" + missingFeature, usedFeatureNames);
        goldenPionBDTVector.emplace_back("goldenPion_minus_" + missingFeature, usedFeatureNames);
    }

    // Setup the counters (each element of the vector is a different BDT)
    const auto nBDTs = reducedFeatureNames.size();
    std::vector<float> trainingSignalTotalProton(nBDTs, 0.f);
    std::vector<float> trainingSignalTotalGoldenPion(nBDTs, 0.f);
    std::vector<float> trainingSignalSelectedProton(nBDTs, 0.f);
    std::vector<float> trainingSignalSelectedGoldenPion(nBDTs, 0.f);
    std::vector<float> trainingBackgroundSelectedProton(nBDTs, 0.f);
    std::vector<float> trainingBackgroundSelectedGoldenPion(nBDTs, 0.f);
    std::vector<float> testingSignalTotalProton(nBDTs, 0.f);
    std::vector<float> testingSignalTotalGoldenPion(nBDTs, 0.f);
    std::vector<float> testingSignalSelectedProton(nBDTs, 0.f);
    std::vector<float> testingSignalSelectedGoldenPion(nBDTs, 0.f);
    std::vector<float> testingBackgroundSelectedProton(nBDTs, 0.f);
    std::vector<float> testingBackgroundSelectedGoldenPion(nBDTs, 0.f);

    // Now test the BDTs
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
            
            const bool isGoldenPion = !isExternal && truePdgCode == 211 && trueIsGolden;
            const bool isProton = !isExternal && truePdgCode == 2212;

            // Only use good matches for testing
            if (config.trainBDTs.onlyGoodTruthMatches && (isExternal || completeness < 0.5f))
                continue;
            
            // Use the BDTs
            for (unsigned int iBDT = 0; iBDT < reducedFeatureNames.size(); ++iBDT)
            {
                const auto &usedFeatureNames = reducedFeatureNames.at(iBDT);

                // Extract the features
                std::vector<float> features;
                const auto areAllFeaturesAvailable = BDTHelper::GetBDTFeatures(recoParticle, usedFeatureNames, features);

                // Only use particles with all features available
                if (!areAllFeaturesAvailable)
                    continue;

                // Add the particle to the BDTs
                const auto protonBDTResponse = protonBDTVector.at(iBDT).GetResponse(features);
                const auto isSelectedProton = protonBDTResponse > 0.f;

                auto &signalTotalProton = (isTrainingEvent ? trainingSignalTotalProton.at(iBDT) : testingSignalTotalProton.at(iBDT));
                auto &signalSelectedProton = (isTrainingEvent ? trainingSignalSelectedProton.at(iBDT) : testingSignalSelectedProton.at(iBDT));
                auto &backgroundSelectedProton = (isTrainingEvent ? trainingBackgroundSelectedProton.at(iBDT) : testingBackgroundSelectedProton.at(iBDT));
            
                signalTotalProton += (isProton ? eventWeight : 0.f);
                signalSelectedProton += (( isProton && isSelectedProton ) ? eventWeight : 0.f);
                backgroundSelectedProton += (( !isProton && isSelectedProton ) ? eventWeight : 0.f);
                

                const auto goldenPionBDTResponse = goldenPionBDTVector.at(iBDT).GetResponse(features);
                const auto isSelectedGoldenPion = goldenPionBDTResponse > 0.f;
                
                auto &signalTotalGoldenPion = (isTrainingEvent ? trainingSignalTotalGoldenPion.at(iBDT) : testingSignalTotalGoldenPion.at(iBDT));
                auto &signalSelectedGoldenPion = (isTrainingEvent ? trainingSignalSelectedGoldenPion.at(iBDT) : testingSignalSelectedGoldenPion.at(iBDT));
                auto &backgroundSelectedGoldenPion = (isTrainingEvent ? trainingBackgroundSelectedGoldenPion.at(iBDT) : testingBackgroundSelectedGoldenPion.at(iBDT));
                
                signalTotalGoldenPion += (isGoldenPion ? eventWeight : 0.f);
                signalSelectedGoldenPion += (( isGoldenPion && isSelectedGoldenPion ) ? eventWeight : 0.f);
                backgroundSelectedGoldenPion += (( !isGoldenPion && isSelectedGoldenPion ) ? eventWeight : 0.f);
            }
        }
    }
    

    // Put the table data is a sensible format
    std::vector< std::tuple< std::string, std::vector<float>, std::vector<float>, std::vector<float> > > tableData;
    tableData.emplace_back("nMinusOne_proton_training.md", trainingSignalTotalProton, trainingSignalSelectedProton, trainingBackgroundSelectedProton);
    tableData.emplace_back("nMinusOne_proton_testing.md", testingSignalTotalProton, testingSignalSelectedProton, testingBackgroundSelectedProton);
    tableData.emplace_back("nMinusOne_goldenPion_training.md", trainingSignalTotalGoldenPion, trainingSignalSelectedGoldenPion, trainingBackgroundSelectedGoldenPion);
    tableData.emplace_back("nMinusOne_goldenPion_testing.md", testingSignalTotalGoldenPion, testingSignalSelectedGoldenPion, testingBackgroundSelectedGoldenPion);

    for (const auto [tableFileName, signalTotalVect, signalSelectedVect, backgroundSelectedVect] : tableData)
    {
        std::cout << "Writing table: " << tableFileName << std::endl;

        FormattingHelper::Table table({"Missing feature", "", "Signal total", "Signal selected", "Background selected", "", "Efficiency", "Purity", "ExP"});

        for (unsigned int iBDT = 0; iBDT < reducedFeatureNames.size(); ++iBDT)
        {
            // Get the trained BDT
            const auto missingFeature = featureNames.at(iBDT);

            const auto signalTotal =  signalTotalVect.at(iBDT);
            const auto signalSelected =  signalSelectedVect.at(iBDT);
            const auto backgroundSelected =  backgroundSelectedVect.at(iBDT);
            const auto allSelected = signalSelected + backgroundSelected;

            if (signalTotal <= std::numeric_limits<float>::epsilon())
                throw std::logic_error("There are no signal events!");

            if (allSelected <= std::numeric_limits<float>::epsilon())
                throw std::logic_error("There are no selected events when removing: " + missingFeature);

            const auto efficiency = signalSelected / signalTotal;
            const auto purity = signalSelected / allSelected;
            const auto efficiencyTimesPurity = efficiency * purity;

            table.AddEmptyRow();
            table.SetEntry("Missing feature", missingFeature);
            table.SetEntry("Signal total", signalTotal);
            table.SetEntry("Signal selected", signalSelected);
            table.SetEntry("Background selected", backgroundSelected);
            table.SetEntry("Efficiency", efficiency);
            table.SetEntry("Purity", purity);
            table.SetEntry("ExP", efficiencyTimesPurity);
        }

        table.WriteToFile(tableFileName);
    }   
}

} // namespace ubc1pi_macros
