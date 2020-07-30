/**
 *  @file  ubcc1pi_standalone/Macros/MomentumThresholdsStudy.cxx
 *
 *  @brief The implementation file of the MomentumThresholdsStudy macro
 */

#include "ubcc1pi_standalone/Macros/Macros.h"

#include "ubcc1pi_standalone/Objects/FileReader.h"

#include "ubcc1pi_standalone/Helpers/AnalysisHelper.h"
#include "ubcc1pi_standalone/Helpers/PlottingHelper.h"
#include "ubcc1pi_standalone/Helpers/SelectionHelper.h"
#include "ubcc1pi_standalone/Helpers/FormattingHelper.h"

#include <map>

#include <TH2F.h>

using namespace ubcc1pi;

namespace ubcc1pi_macros
{

void MomentumThresholdsStudy(const Config &config)
{
    // Open the file
    FileReader reader(config.files.overlaysFileName);
    auto pEvent = reader.GetBoundEventAddress();
    const auto nEvents = reader.GetNumberOfEvents();

    // Get the selection and switch off the momentum threshold cuts
    auto selection = SelectionHelper::GetDefaultSelection();
    selection.SetCutNominalValue("muonRange", -std::numeric_limits<float>::max());
    selection.SetCutNominalValue("pionRange", -std::numeric_limits<float>::max());

    // Setup the plots
    std::shared_ptr<TH2F> pHistErrorMuonTruth(new TH2F("histErrorMuonTruth", "", 30u, 0.f, 0.8f, 30u, -1.f, 0.6f));
    std::shared_ptr<TH2F> pHistErrorMuonReco(new TH2F("histErrorMuonReco", "", 30u, 0.f, 0.8f, 30u, -1.f, 0.6f));
    std::shared_ptr<TH2F> pHistErrorPionTruth(new TH2F("histErrorPionTruth", "", 30u, 0.f, 0.4f, 30u, -1.f, 1.f));
    std::shared_ptr<TH2F> pHistErrorPionReco(new TH2F("histErrorPionReco", "", 30u, 0.f, 0.4f, 30u, -1.f, 1.f));
    
    std::shared_ptr<TH2F> pHistMuonTruthReco(new TH2F("histMuonTruthReco", "", 30u, 0.f, 0.8f, 30u, 0.f, 0.8f));
    std::shared_ptr<TH2F> pHistPionTruthReco(new TH2F("histPionTruthReco", "", 30u, 0.f, 0.4f, 30u, 0.f, 0.4f));
    
    // Loop over the events
    for (unsigned int i = 0; i < nEvents; ++i)
    {
        AnalysisHelper::PrintLoadingBar(i, nEvents);
        reader.LoadEvent(i);

        // Only consider true CC1Pi events
        if (!AnalysisHelper::IsTrueCC1Pi(pEvent, config.global.useAbsPdg))
            continue;
            
        // Insist we pass the CC inclusive
        if (!pEvent->reco.passesCCInclusive())
            continue;
    
        // Insist we pass the generic selection
        std::vector<std::string> cutsPassed;
        std::vector<int> assignedPdgCodes;
        const auto passedGoldenSelection = selection.Execute(pEvent, cutsPassed, assignedPdgCodes);
        const auto passedGenericSelection = (std::find(cutsPassed.begin(), cutsPassed.end(), config.global.lastCutGeneric) != cutsPassed.end());
        
        if (!passedGenericSelection)
            continue;

        // Get the analysis data
        const auto truthData = AnalysisHelper::GetTruthAnalysisData(pEvent->truth, config.global.useAbsPdg, config.global.protonMomentumThreshold);
        const auto recoData = AnalysisHelper::GetRecoAnalysisData(pEvent->reco, assignedPdgCodes, passedGoldenSelection);
        const auto weight = AnalysisHelper::GetNominalEventWeight(pEvent);

        // Fill the muon plots
        const auto muonMomentumError = (recoData.muonMomentum - truthData.muonMomentum) / truthData.muonMomentum;
        pHistErrorMuonTruth->Fill(truthData.muonMomentum, muonMomentumError, weight); 
        pHistErrorMuonReco->Fill(recoData.muonMomentum, muonMomentumError, weight); 
        pHistMuonTruthReco->Fill(truthData.muonMomentum, recoData.muonMomentum);
        
        if (!recoData.hasGoldenPion)
            continue;

        // Fill the pion plots
        const auto pionMomentumError = (recoData.pionMomentum - truthData.pionMomentum) / truthData.pionMomentum;
        pHistErrorPionTruth->Fill(truthData.pionMomentum, pionMomentumError, weight);
        pHistErrorPionReco->Fill(recoData.pionMomentum, pionMomentumError, weight);
        pHistPionTruthReco->Fill(truthData.pionMomentum, recoData.pionMomentum);
    }

    // Draw the plots
    auto pCanvas = PlottingHelper::GetCanvas();

    pHistErrorMuonTruth->Draw("colz");
    PlottingHelper::SaveCanvas(pCanvas, "momentumThresholds_genericSelection_muon_momError-trueMomentum");
    
    pHistErrorMuonReco->Draw("colz");
    PlottingHelper::SaveCanvas(pCanvas, "momentumThresholds_genericSelection_muon_momError-recoMomentum");
    
    pHistMuonTruthReco->Draw("colz");
    PlottingHelper::SaveCanvas(pCanvas, "momentumThresholds_genericSelection_muon_recoMomentum-truthMomentum");
    
    
    pHistErrorPionTruth->Draw("colz");
    PlottingHelper::SaveCanvas(pCanvas, "momentumThresholds_goldenSelection_pion_momError-trueMomentum");
    
    pHistErrorPionReco->Draw("colz");
    PlottingHelper::SaveCanvas(pCanvas, "momentumThresholds_goldenSelection_pion_momError-recoMomentum");
    
    pHistPionTruthReco->Draw("colz");
    PlottingHelper::SaveCanvas(pCanvas, "momentumThresholds_goldenSelection_pion_recoMomentum-truthMomentum");

    /*
    // Setup the plots (indexed by histMaps[descriptor][pdg code])
    using PDGToHistMap = std::map<int, std::shared_ptr<TH1F> >;

    std::map<std::string, PDGToHistMap > histMaps;
    histMaps.emplace("all", PDGToHistMap());
    histMaps.emplace("withHit", PDGToHistMap());
    histMaps.emplace("withMatch", PDGToHistMap());
    histMaps.emplace("withMatchAccurate", PDGToHistMap());

    for (unsigned int i = 0; i < nEvents; ++i)
    {
        AnalysisHelper::PrintLoadingBar(i, nEvents);
        reader.LoadEvent(i);

        // Only consider true CC1Pi events
        if (!AnalysisHelper::IsTrueCC1Pi(pEvent, config.global.useAbsPdg))
            continue;
            
        // Insist we pass the CC inclusive
        if (!pEvent->reco.passesCCInclusive())
            continue;
    
        // Insist we pass the generic selection
        std::vector<std::string> cutsPassed;
        std::vector<int> assignedPdgCodes;
        const auto passedGoldenSelection = selection.Execute(pEvent, cutsPassed, assignedPdgCodes);
        const auto passedGenericSelection = (std::find(cutsPassed.begin(), cutsPassed.end(), config.global.lastCutGeneric) != cutsPassed.end());
        
        if (!passedGenericSelection)
            continue;    
      
        const auto weight = AnalysisHelper::GetNominalEventWeight(pEvent);
        const auto &truthParticles = pEvent->truth.particles;
        const auto &recoParticles = pEvent->reco.particles;

        // Get the mapping from truth particles to reco particles
        std::map<unsigned int, std::vector<unsigned int> > trueToRecoIndicesMap;
        std::map<unsigned int, std::vector<float> > trueToRecoCompletenessesMap;

        for (unsigned int recoIndex = 0; recoIndex < recoParticles.size(); ++recoIndex)
        {
            const auto &recoParticle = recoParticles.at(recoIndex);

            try
            {
                const auto truthIndex = AnalysisHelper::GetBestMatchedTruthParticleIndex(recoParticle, truthParticles);
                const auto completeness = recoParticle.truthMatchCompletenesses().at(truthIndex);

                trueToRecoIndicesMap[truthIndex].push_back(recoIndex);
                trueToRecoCompletenessesMap[truthIndex].push_back(completeness);
            }
            catch (const std::exception &)
            {
            }
        }

        // Now consider all truth particles
        for (unsigned int truthIndex = 0; truthIndex < truthParticles.size(); ++truthIndex)
        {
            const auto &truthParticle = truthParticles.at(truthIndex);
            const auto pdgCode = config.global.useAbsPdg ? std::abs(truthParticle.pdgCode()) : truthParticle.pdgCode();

            // Ensure this is a particle we care about
            if (!AnalysisHelper::PassesVisibilityThreshold(truthParticle))
                continue;

            // Apply the thresholds we care about for measuring momentum by range
            if (!AnalysisHelper::IsContained(truthParticle))
                continue;

            // Only use golden pions
            if (pdgCode == 211 && !AnalysisHelper::IsGolden(truthParticle))
                continue;
            
            // Get the true momentum
            const auto momentum = truthParticle.momentum();

            // Get the reco matches
            const auto recoMatchIndices = trueToRecoIndicesMap[truthIndex];
            const auto hasRecoMatch = !recoMatchIndices.empty();

            // Setup the plots if this is the first time we have seen this PDG code
            for (auto &histMapEntry : histMaps)
            {
                const auto &descriptor = histMapEntry.first;
                auto &histMap = histMapEntry.second;

                if (histMap.find(pdgCode) != histMap.end())
                    continue;
                
                switch (std::abs(pdgCode))
                {
                    // Protons
                    case 2212:
                        histMap.emplace(pdgCode, new TH1F(("thresholds_" + std::to_string(pdgCode) + "_" + descriptor).c_str(), "", 60u, 0.f, 1.2f));
                        break;
                    // Pions
                    case 211:
                        histMap.emplace(pdgCode, new TH1F(("thresholds_" + std::to_string(pdgCode) + "_" + descriptor).c_str(), "", 40u, 0.f, 0.3f));
                        break;
                    // Muons / Default
                    default:
                        histMap.emplace(pdgCode, new TH1F(("thresholds_" + std::to_string(pdgCode) + "_" + descriptor).c_str(), "", 30u, 0.f, 0.3f));
                }
            }

            // Fill the plots
            histMaps.at("all").at(pdgCode)->Fill(momentum, weight);
        
            // Count the hits
            // ATTN we actually store hit weights (each hit that matches to the particle is given a weight = fraction of charge of the hit
            // contributed by the particle). Here we take the summed floor of these weights as a measure of the number of hits.
            const auto nHits = static_cast<unsigned int>(std::floor(
                        std::floor(truthParticle.hitWeightU()) + 
                        std::floor(truthParticle.hitWeightV()) + 
                        std::floor(truthParticle.hitWeightW())));

            if (nHits > 0u)
            {
                histMaps.at("withHit").at(pdgCode)->Fill(momentum, weight);
            }

            if (!hasRecoMatch)
                continue;
            
            histMaps.at("withMatch").at(pdgCode)->Fill(momentum, weight);

    
            // Only attempt to measure the momentum in reco of pions and muons
            if (pdgCode != 13 && pdgCode != 211)
                continue;

            // Get the best matched reco particles (most complete)
            const auto recoMatchCompletenesses = trueToRecoCompletenessesMap.at(truthIndex);
            const auto bestRecoParticle = recoParticles.at(recoMatchIndices.at(std::distance(recoMatchCompletenesses.begin(), std::max_element(recoMatchCompletenesses.begin(), recoMatchCompletenesses.end()))));

            if (!bestRecoParticle.range.IsSet())
                continue;
            
            const auto recoMomentum = (pdgCode == 13 ? AnalysisHelper::GetMuonMomentum(bestRecoParticle) : AnalysisHelper::GetPionMomentumFromRange(bestRecoParticle.range()));
            const auto fracMomentumError = std::abs((recoMomentum - momentum) / momentum);

            if (fracMomentumError > 0.2f)
                continue;
            
            histMaps.at("withMatchAccurate").at(pdgCode)->Fill(momentum, weight);
        }
    }

    // Save the plots
    auto pCanvas = PlottingHelper::GetCanvas();

    for (const auto &pdgCode : {13, 211, 2212})
    {
        PlottingHelper::SetLineStyle(histMaps.at("all").at(pdgCode), PlottingHelper::Primary);
        PlottingHelper::SetLineStyle(histMaps.at("withHit").at(pdgCode), PlottingHelper::Secondary);
        PlottingHelper::SetLineStyle(histMaps.at("withMatch").at(pdgCode), PlottingHelper::Tertiary);
        PlottingHelper::SetLineStyle(histMaps.at("withMatchAccurate").at(pdgCode), PlottingHelper::Quaternary);
        
        histMaps.at("all").at(pdgCode)->Draw("hist");
        histMaps.at("withHit").at(pdgCode)->Draw("hist same");
        histMaps.at("withMatch").at(pdgCode)->Draw("hist same");
        histMaps.at("withMatchAccurate").at(pdgCode)->Draw("hist same");
        PlottingHelper::SaveCanvas(pCanvas, "momentumThresholds_momentum_" + std::to_string(pdgCode));
        
        histMaps.at("withHit").at(pdgCode)->Divide(histMaps.at("all").at(pdgCode).get());
        histMaps.at("withMatch").at(pdgCode)->Divide(histMaps.at("all").at(pdgCode).get());
        histMaps.at("withMatchAccurate").at(pdgCode)->Divide(histMaps.at("all").at(pdgCode).get());
        histMaps.at("withHit").at(pdgCode)->Draw("hist");
        histMaps.at("withHit").at(pdgCode)->GetYaxis()->SetRangeUser(0.f, 1.f);
        histMaps.at("withMatch").at(pdgCode)->Draw("hist same");
        histMaps.at("withMatchAccurate").at(pdgCode)->Draw("hist same");
        PlottingHelper::SaveCanvas(pCanvas, "momentumThresholds_momentum_ratio_" + std::to_string(pdgCode));
    }
    */
}

} // namespace ubcc1pi_macros
