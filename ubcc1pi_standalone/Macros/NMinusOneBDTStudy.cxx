/**
 *  @file  ubcc1pi_standalone/Macros/NMinusOneBDTStudy.cxx
 *
 *  @brief The implementation file of the NMinusOneBDTStudy macro
 */

#include "ubcc1pi_standalone/Macros/Macros.h"

#include "ubcc1pi_standalone/Objects/FileReader.h"

#include "ubcc1pi_standalone/Helpers/BDTHelper.h"
#include "ubcc1pi_standalone/Helpers/AnalysisHelper.h"
#include "ubcc1pi_standalone/Helpers/FormattingHelper.h"

#include <TGraph.h>
#include <TLine.h>
#include <TH1F.h>

using namespace ubcc1pi;

namespace ubcc1pi_macros
{

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
            throw std::invalid_argument("NMinusOneBDTStudy - Input signal type must be: proton, muon or golden pion");
    }

    // Make copies of the list of features but remove one each time
    std::vector< std::vector<std::string> > reducedFeatureNames;
    std::vector< BDTHelper::BDTFactory > bdtFactoryVector;

    for (unsigned int i = 0; i < config.nMinusOneBDTStudy.featureNames.size(); ++i)
    {
        std::vector<std::string> nameVector;
        for (unsigned int j = 0; j < config.nMinusOneBDTStudy.featureNames.size(); ++j)
        {
            if (i == j) continue;

            nameVector.push_back(config.nMinusOneBDTStudy.featureNames.at(j));
        }

        // Make the BDT and store the features
        if (config.nMinusOneBDTStudy.shouldTrainBDTs)
        {
            bdtFactoryVector.emplace_back(signalString + "_N-" + nFeaturesString + "_minus_" + config.nMinusOneBDTStudy.featureNames.at(i), nameVector);
        }
            
        reducedFeatureNames.push_back(nameVector);
    }

    // Now add the BDT using all features
    if (config.nMinusOneBDTStudy.shouldTrainBDTs)
    {
        bdtFactoryVector.emplace_back(signalString + "_N-" + nFeaturesString + "_minus_none", config.nMinusOneBDTStudy.featureNames);
    }

    reducedFeatureNames.push_back(config.nMinusOneBDTStudy.featureNames);


    // Train the BDTs if required
    if (config.nMinusOneBDTStudy.shouldTrainBDTs)
    {

        // Fill the BDT factories with training examples
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
                    truePdgCode = config.global.useAbsPdg ? std::abs(truthParticle.pdgCode()) : truthParticle.pdgCode();
                    trueIsGolden = AnalysisHelper::IsGolden(truthParticle);
                    completeness = recoParticle.truthMatchCompletenesses().at(truthParticleIndex);
                }
                catch (const std::exception &) {}

                // Only use good matches for training
                if (config.trainBDTs.onlyGoodTruthMatches && (isExternal || completeness < 0.5f))
                    continue;

                // Define the weight
                const auto weight = eventWeight * (config.trainBDTs.weightByCompleteness ? (isExternal ? 1.f : completeness) : 1.f);

                bool isSignal = false;
                switch (config.nMinusOneBDTStudy.signalType)
                {
                    case PlottingHelper::Muon:
                        isSignal = !isExternal && truePdgCode == 13;
                        break;
                    case PlottingHelper::Proton:
                        isSignal = !isExternal && truePdgCode == 2212;
                        break;
                    case PlottingHelper::GoldenPion:
                        isSignal = !isExternal && truePdgCode == 211 && trueIsGolden;
                        break;
                    default:
                        throw std::invalid_argument("NMinusOneBDTStudy - Input signal type must be: proton, muon or golden pion");
                }

                // Fill the BDTs 
                for (unsigned int iBDT = 0; iBDT < reducedFeatureNames.size(); ++iBDT)
                {
                    const auto &usedFeatureNames = reducedFeatureNames.at(iBDT);

                    // Extract the features
                    std::vector<float> features;
                    const auto areAllFeaturesAvailable = BDTHelper::GetBDTFeatures(recoParticle, usedFeatureNames, features);

                    // Only use particles with all features available
                    if (!areAllFeaturesAvailable)
                        continue;

                    // Add the particle to the BDT
                    bdtFactoryVector.at(iBDT).AddEntry(features, isSignal, isTrainingEvent, weight);
                }
            }
        }

        // Train the BDTs (also runs the standard TMVA testing for overtraining etc)
        std::cout << "Training and testing the BDTs" << std::endl;
        for (auto &bdtFactory : bdtFactoryVector)
            bdtFactory.TrainAndTest();

    }
    
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Get the performance table / plots
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    std::cout << "Comparing performances of BDTs" << std::endl;

    // Get the trained BDTs
    std::vector< BDTHelper::BDT > bdtVector;
    for (unsigned int iBDT = 0; iBDT < reducedFeatureNames.size(); ++iBDT)
    {
        // Get the trained BDT
        const auto missingFeature = (iBDT < config.nMinusOneBDTStudy.featureNames.size()) ? config.nMinusOneBDTStudy.featureNames.at(iBDT) : "none";

        const auto &usedFeatureNames = reducedFeatureNames.at(iBDT);
        bdtVector.emplace_back(signalString + "_N-" + nFeaturesString + "_minus_" + missingFeature, usedFeatureNames);
    }

    // Setup the structure to hold the data
    BDTHelper::BDTResponseMap bdtResponseMap;
    
    // Now test the BDTs
    for (unsigned int i = 0; i < nCC1PiEvents; ++i)
    {
        AnalysisHelper::PrintLoadingBar(i, nCC1PiEvents);

        const auto eventIndex = cc1PiEventIndices.at(i);
        const auto isTrainingEvent = shuffler.IsTrainingEvent(i);
        if (isTrainingEvent)
            continue;

        reader.LoadEvent(eventIndex);

        const auto truthParticles = pEvent->truth.particles;
        const auto recoParticles = pEvent->reco.particles;
        const float eventWeight = AnalysisHelper::GetNominalEventWeight(pEvent);
        
        for (const auto &recoParticle : recoParticles)
        {
            // Only use contained particles for testing
            if (!AnalysisHelper::HasTrackFit(recoParticle) || !AnalysisHelper::IsContained(recoParticle))
                continue;

            // Find the MC origin of this reco particle
            const auto particleType = PlottingHelper::GetPlotStyle(recoParticle, AnalysisHelper::Overlay, truthParticles, false, config.global.useAbsPdg);

            // Loop over all BDTs
            for (unsigned int iBDT = 0; iBDT < reducedFeatureNames.size(); ++iBDT)
            {
                const auto missingFeature = (iBDT < config.nMinusOneBDTStudy.featureNames.size()) ? config.nMinusOneBDTStudy.featureNames.at(iBDT) : "none";
                const auto &usedFeatureNames = reducedFeatureNames.at(iBDT);

                // Get the features if they are available
                std::vector<float> features;
                if (!BDTHelper::GetBDTFeatures(recoParticle, usedFeatureNames, features))
                    continue;

                // Get the BDT responses
                const auto bdtResponse = bdtVector.at(iBDT).GetResponse(features);
            
                // Store the result
                bdtResponseMap[missingFeature][particleType].emplace_back(bdtResponse, eventWeight);

                // Add to the special category that sums all backgrounds
                if (particleType != config.nMinusOneBDTStudy.signalType)
                    bdtResponseMap[missingFeature][PlottingHelper::Other].emplace_back(bdtResponse, eventWeight);
            }
        }
    }

    // Define the particle types we care about
    const std::vector<PlottingHelper::PlotStyle> particleTypes = {
        PlottingHelper::Other, // All backgrounds together
        PlottingHelper::External,
        PlottingHelper::Proton,
        PlottingHelper::Muon,
        PlottingHelper::NonGoldenPion,
        PlottingHelper::GoldenPion
    };


    // Make the N-1 plots
    auto pCanvas = PlottingHelper::GetCanvas();
    std::map<PlottingHelper::PlotStyle, TH1F *> histMap;
    const auto allFeautreNames = BDTHelper::ParticleBDTFeatureNames;
    for (const auto &type : particleTypes)
    {
        if (type == config.nMinusOneBDTStudy.signalType)
            continue;

        auto pHist = new TH1F(("pHist_" + std::to_string(type)).c_str(), "", reducedFeatureNames.size(), 0, reducedFeatureNames.size());
        auto pAxis = pHist->GetXaxis();

        PlottingHelper::SetLineStyle(pHist, type);

        // Use a dashed line for the other category
        if (type == PlottingHelper::Other)
            pHist->SetLineStyle(2);

        for (unsigned int iBDT = 0; iBDT < reducedFeatureNames.size(); ++iBDT)
        {
            if (iBDT < config.nMinusOneBDTStudy.featureNames.size())
            {
                const auto missingFeature = config.nMinusOneBDTStudy.featureNames.at(iBDT);
                const auto iter = std::find(allFeautreNames.begin(), allFeautreNames.end(), missingFeature);

                if (iter == allFeautreNames.end())
                    throw std::logic_error("NMinusOneBDTStudy - Unknown feature: " + missingFeature);

                const auto featureIndex = std::distance(allFeautreNames.begin(), iter);
                pAxis->SetBinLabel(iBDT + 1, std::to_string(featureIndex).c_str());
            }
            else
            {
                pAxis->SetBinLabel(iBDT + 1, "none");
            }
        }

        histMap.emplace(type, pHist);
    }

    // Fill the plots
    float minROCIntegral = std::numeric_limits<float>::max();
    float maxROCIntegral = -std::numeric_limits<float>::max();

    for (unsigned int iBDT = 0; iBDT < reducedFeatureNames.size(); ++iBDT)
    {
        const auto missingFeature = (iBDT < config.nMinusOneBDTStudy.featureNames.size()) ? config.nMinusOneBDTStudy.featureNames.at(iBDT) : "none";

        const auto &bdtResponses = bdtResponseMap.at(missingFeature);
        const auto &signalResponses = bdtResponses.at(config.nMinusOneBDTStudy.signalType);

        bool isFirst = true;
        std::vector<TGraph *> curves;
        for (const auto &background : particleTypes)
        {
            if (background == config.nMinusOneBDTStudy.signalType)
                continue;

            const auto &backgroundResponses = bdtResponses.at(background);

            // Get the ROC curve
            std::vector<float> signalPassingRates, backgroundRejectionRates;
            std::vector<float> signalPassingRateErrs, backgroundRejectionRateErrs;
            BDTHelper::GetROCCurve(signalResponses, backgroundResponses, config.nMinusOneBDTStudy.nSamplePoints, signalPassingRates, backgroundRejectionRates, signalPassingRateErrs, backgroundRejectionRateErrs);

            // Get the ROC integral
            const auto rocIntegral = BDTHelper::GetROCIntegral(signalPassingRates, backgroundRejectionRates);
            const auto rocIntegralErr = BDTHelper::GetROCIntegralError(signalPassingRates, backgroundRejectionRates, signalPassingRateErrs, backgroundRejectionRateErrs);
            histMap.at(background)->SetBinContent(iBDT + 1, rocIntegral);
            histMap.at(background)->SetBinError(iBDT + 1, rocIntegralErr);

            minROCIntegral = std::min(minROCIntegral, rocIntegral);
            maxROCIntegral = std::max(maxROCIntegral, rocIntegral);

            TGraph *pROCCurve = new TGraph(config.nMinusOneBDTStudy.nSamplePoints, signalPassingRates.data(), backgroundRejectionRates.data());
            PlottingHelper::SetLineStyle(pROCCurve, PlottingHelper::GetColor(background));

            pROCCurve->Draw(isFirst ? "AL" : "L");
            isFirst = false;
        }

        PlottingHelper::SaveCanvas(pCanvas, "rocCurve_N-" + nFeaturesString + "_" + signalString + "_minus_" + missingFeature);

        // Clean up the heap
        for (auto &pGraph : curves)
            delete pGraph;
    }
        
    // Make the ROC integral plot
    float padding = (maxROCIntegral - minROCIntegral) * 0.1;
    minROCIntegral -= padding;
    maxROCIntegral += padding;

    bool isFirstHist = true;
    for (const auto &background : particleTypes)
    {
        if (background == config.nMinusOneBDTStudy.signalType)
            continue;

        auto pHist = histMap.at(background);
        pHist->GetYaxis()->SetRangeUser(minROCIntegral, maxROCIntegral);
       
        // Draw the histogram with it's uncertainties
        auto pHistClone = static_cast<TH1F *>(pHist->Clone());
        pHistClone->SetFillStyle(1001);
        pHistClone->SetLineColorAlpha(pHist->GetLineColor(), 0.f);
        pHistClone->SetFillColorAlpha(pHist->GetLineColor(), 0.3f);
        pHistClone->Draw(isFirstHist ? "e2" : "e2 same");
        pHist->Draw("hist same");
    
        // Draw the line at the none value
        const auto nBins = pHist->GetNbinsX();
        const auto noneROCIntegral = pHist->GetBinContent(nBins);
        TLine *pLine = new TLine(0, noneROCIntegral, nBins, noneROCIntegral);
        pLine->SetLineColor(pHist->GetLineColor());
        pLine->SetLineStyle(3);
        pLine->Draw();

        isFirstHist = false;
    }
    PlottingHelper::SaveCanvas(pCanvas, "rocCurveIntegrals_N-" + nFeaturesString + "_" + signalString);

    // Make the output table
    FormattingHelper::Table table({"Feature removed", "ID", "", "ROC integral", "Uncertainty", "", "Difference", "Sigma", "", "Should remove"});

    const auto pHist = histMap.at(PlottingHelper::Other);
    const auto noneROCIntegral = pHist->GetBinContent(pHist->GetNbinsX());
    const auto noneROCIntegralErr = pHist->GetBinError(pHist->GetNbinsX());
    
    // Store the features that aren't doing much
    std::vector< std::pair<unsigned int, float> > badFeatures;

    for (unsigned int iBDT = 0; iBDT < reducedFeatureNames.size(); ++iBDT)
    {
        const auto missingFeature = (iBDT < config.nMinusOneBDTStudy.featureNames.size()) ? config.nMinusOneBDTStudy.featureNames.at(iBDT) : "none";

        table.AddEmptyRow();
        table.SetEntry("Feature removed", missingFeature);

        if (missingFeature != "none")
        {
            const auto iter = std::find(allFeautreNames.begin(), allFeautreNames.end(), missingFeature);
            if (iter == allFeautreNames.end())
                throw std::logic_error("NMinusOneBDTStudy - Unknown feature: " + missingFeature);

            const auto featureIndex = std::distance(allFeautreNames.begin(), iter);
            table.SetEntry("ID", featureIndex);
        }


        // Get the ROC integral for this feature
        const auto rocIntegral = pHist->GetBinContent(iBDT + 1);
        const auto rocIntegralErr = pHist->GetBinError(iBDT + 1);
        table.SetEntry("ROC integral", rocIntegral);
        table.SetEntry("Uncertainty", rocIntegralErr);

        // Get the drop wrt. no features removed
        const auto diff = noneROCIntegral - rocIntegral;
        table.SetEntry("Difference", diff);

        // Get the sigma for this difference
        const auto combinedErr = std::pow(rocIntegralErr*rocIntegralErr + noneROCIntegralErr*noneROCIntegralErr, 0.5f);
        const auto sigma = diff / combinedErr;
        table.SetEntry("Sigma", sigma);

        if (missingFeature == "none")
            continue;
        
        // Assume the feature shouldn't be removed for now
        table.SetEntry("Should remove", "No");

        // These features clearly do some good
        if (sigma > 1.f)
            continue;

        badFeatures.emplace_back(iBDT, diff);
    }

    if (!badFeatures.empty())
    {
        const auto worstFeatureIndex = std::min_element(badFeatures.begin(), badFeatures.end(), [](const auto &a, const auto &b){return a.second < b.second;})->first;
        table.SetEntry("Should remove", worstFeatureIndex, "Yes");
    }

    table.WriteToFile("rocCurveIntegralsTable_N-" + nFeaturesString + "_" + signalString + ".md");
        
    // Clean up the heap
    for (auto &entry : histMap)
        delete entry.second;

}

} // namespace ubc1pi_macros
