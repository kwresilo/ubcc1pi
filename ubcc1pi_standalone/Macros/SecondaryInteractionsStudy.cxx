/**
 *  @file  ubcc1pi_standalone/Macros/SecondaryInteractionsStudy.cxx
 *
 *  @brief The implementation file of the SecondaryInteractionsStudy macro
 */

#include "ubcc1pi_standalone/Macros/Macros.h"

#include "ubcc1pi_standalone/Objects/FileReader.h"

#include "ubcc1pi_standalone/Helpers/AnalysisHelper.h"
#include "ubcc1pi_standalone/Helpers/PlottingHelper.h"
#include "ubcc1pi_standalone/Helpers/FormattingHelper.h"

#include <string>
#include <vector>
#include <utility>
#include <map>

#include <TH1F.h>
#include <TH2F.h>
#include <TStyle.h>
#include <TGraph.h>

using namespace ubcc1pi;

namespace ubcc1pi_macros
{

void SecondaryInteractionsStudy(const Config &config)
{
    // Open the file
    // const std::string prefix = "NuWro";
    // FileReader reader(config.filesRun1.nuWroFileName);
    const std::string prefix = "Overlay";
    FileReader reader(config.filesRun1.overlaysFileName);
    auto pEvent = reader.GetBoundEventAddress();
    const auto nEvents = reader.GetNumberOfEvents();

    //
    // Setup the data structures for the table
    //
    std::map<int, float> pionEndStateMap;
    std::map<int, float> stoppingPionEndStateMap;
    float totalWeight = 0.f;

    //
    // Setup the plots
    //

    // Pion momentum plot
    PlottingHelper::MultiPlot pionMomentumPlot_endState("Pion momentum / GeV", "Number of particles", 50, 0.f, 1.2f, false);
    PlottingHelper::MultiPlot pionEndMomentumPlot_endState("Pion end momentum / GeV", "Number of particles", 50, 0.f, 1.2f, false);

    // nScatters plot
    std::shared_ptr<TH2F> pNScattersHist(new TH2F("nScatters", "", 3, 0, 3, 3, 0, 3));

    // Momentum lost in scatter plots
    std::shared_ptr<TH1F> pMomLostElastic(new TH1F("momLostElastic", "", 50, 0.f, 1.f));
    std::shared_ptr<TH1F> pMomLostInelastic(new TH1F("momLostInelastic", "", 50, 0.f, 1.f));

    // Reco-true momentum resolution plots (index of map is a descriptor, eg. "contained")
    std::map<std::string, std::shared_ptr<TH1F> > pionMomentumResolutionPlots, muonMomentumResolutionPlots, protonMomentumResolutionPlots;

    pionMomentumResolutionPlots.emplace("all", new TH1F("pionRes_all", "", 50, -1.f, 1.5f));
    pionMomentumResolutionPlots.emplace("contained", new TH1F("pionRes_contained", "", 50, -1.f, 1.5f));
    pionMomentumResolutionPlots.emplace("contained-stopping", new TH1F("pionRes_contained-stopping", "", 50, -1.f, 1.5f));
    pionMomentumResolutionPlots.emplace("contained-stopping-noScatters", new TH1F("pionRes_contained-stopping-noScatters", "", 50, -1.f, 1.5f));

    muonMomentumResolutionPlots.emplace("all_range", new TH1F("muonRes_all_range", "", 50, -1.f, 0.7f));
    muonMomentumResolutionPlots.emplace("contained_range", new TH1F("muonRes_contained_range", "", 50, -1.f, 0.7f));
    muonMomentumResolutionPlots.emplace("all_MCS", new TH1F("muonRes_all_MCS", "", 50, -1.f, 0.7f));
    muonMomentumResolutionPlots.emplace("contained_MCS", new TH1F("muonRes_contained_MCS", "", 50, -1.f, 0.7f));

    protonMomentumResolutionPlots.emplace("all_range", new TH1F("protonRes_all_range", ";(Reco-True)/True Momentum;Number of Particles", 50, -1.f, 0.7f));
    protonMomentumResolutionPlots.emplace("contained_range", new TH1F("protonRes_contained_range", ";(Reco-True)/True Momentum;Number of Particles", 50, -1.f, 0.7f));
    protonMomentumResolutionPlots.emplace("all_range_larsoft", new TH1F("protonRes_all_range_larsoft", ";(Reco-True)/True Momentum;Number of Particles", 50, -1.f, 0.7f));
    protonMomentumResolutionPlots.emplace("contained_range_larsoft", new TH1F("protonRes_contained_range_larsoft", ";(Reco-True)/True Momentum;Number of Particles", 50, -1.f, 0.7f));

    // Reco-true momentum resolution vs. reco range
    std::shared_ptr<TH2F> pPionMomVsRange(new TH2F("pionMomVsRange", "", 60, 0.f, 100.f, 60, 0.f, 0.4f));
    std::shared_ptr<TH2F> pPionResVsRange(new TH2F("pionResVsRange", "", 60, 0.f, 100.f, 60, -1.f, 0.7f));

    std::shared_ptr<TH2F> pMuonMomVsRange(new TH2F("muonMomVsRange", "", 60, 0.f, 300.f, 60, 0.f, 0.8f));
    std::shared_ptr<TH2F> pMuonResVsRange(new TH2F("muonResVsRange", "", 60, 0.f, 300.f, 60, -1.f, 0.7f));

    std::shared_ptr<TH2F> pProtonMomVsRange(new TH2F("protonMomVsRange", ";Range (cm);Momentum (GeV)", 60, 0.f, 100.f, 60, 0.f, 1.2f));
    std::shared_ptr<TH2F> pProtonResVsRange(new TH2F("protonResVsRange", ";Range (cm);(Reco-True)/True Momentum", 60, 0.f, 100.f, 60, -1.f, 0.7f));
    std::shared_ptr<TH2F> pProtonResVsRangeLarsoft(new TH2F("protonResVsRangeLarsoft", ";Range (cm);(Reco-True)/True Momentum", 60, 0.f, 100.f, 60, -1.f, 0.7f));

    // Loop over the events
    for (unsigned int i = 0; i < nEvents; ++i)
    {
        AnalysisHelper::PrintLoadingBar(i, nEvents);
        reader.LoadEvent(i);

        // Only consider true CC1Pi events
        if (!AnalysisHelper::IsTrueCC1Pi(pEvent, config.global.useAbsPdg))
            continue;

        const auto weight = AnalysisHelper::GetNominalEventWeight(pEvent);
        const auto &truthParticles = pEvent->truth.particles;

        // Get the pion & muon
        bool foundPion = false;
        bool foundMuon = false;
        unsigned int pionIndex = std::numeric_limits<unsigned int>::max();
        unsigned int muonIndex = std::numeric_limits<unsigned int>::max();

        for (unsigned int index = 0; index < truthParticles.size(); ++index)
        {
            const auto &particle = truthParticles.at(index);
            const auto pdgCode = config.global.useAbsPdg ? std::abs(particle.pdgCode()) : particle.pdgCode();

            if (pdgCode == 211)
            {
                if (foundPion)
                    throw std::logic_error("SecondaryInteractionsStudy - Found multiple pions in a signal event!");

                pionIndex = index;
                foundPion = true;
            }

            if (pdgCode == 13)
            {
                if (foundMuon)
                    throw std::logic_error("SecondaryInteractionsStudy - Found multiple muons in a signal event!");

                muonIndex = index;
                foundMuon = true;
            }
        }

        if (!foundPion)
            throw std::logic_error("SecondaryInteractionsStudy - Couldn't find the pion in a signal event!");

        if (!foundMuon)
            throw std::logic_error("SecondaryInteractionsStudy - Couldn't find the muon in a signal event!");

        // Get the relevant information
        const auto &pion = truthParticles.at(pionIndex);
        const auto &muon = truthParticles.at(muonIndex);

        const auto pionMomentum = pion.momentum();
        const auto pionEndMomentum = pion.endMomentum();
        const auto pionEndState = pion.endState();

        const auto muonMomentum = muon.momentum();

        // Get the reco particles
        for (const auto &recoParticle : pEvent->reco.particles)
        {
            if (!AnalysisHelper::HasTrackFit(recoParticle))
                continue;

            const auto style = PlottingHelper::GetPlotStyle(recoParticle, AnalysisHelper::Overlay, truthParticles, false, config.global.useAbsPdg);

            if (style == PlottingHelper::GoldenPion || style == PlottingHelper::NonGoldenPion)
            {
                const auto pionRecoMomentum = AnalysisHelper::GetPionMomentumFromRange(recoParticle.range());
                const auto pionResolution = (pionRecoMomentum - pionMomentum) / pionMomentum;

                // Fill the pion momentum resolution plots
                pionMomentumResolutionPlots.at("all")->Fill(pionResolution, weight);

                if (AnalysisHelper::IsContained(pion))
                {
                    pionMomentumResolutionPlots.at("contained")->Fill(pionResolution, weight);

                    if (pion.isStopping())
                    {
                        pionMomentumResolutionPlots.at("contained-stopping")->Fill(pionResolution, weight);

                        if (pion.nElasticScatters() == 0 && pion.nInelasticScatters() == 0)
                        {
                            pionMomentumResolutionPlots.at("contained-stopping-noScatters")->Fill(pionResolution, weight);
                        }
                    }
                }

                if (AnalysisHelper::IsGolden(pion))
                {
                    pPionMomVsRange->Fill(recoParticle.range(), pionMomentum, weight);
                    pPionResVsRange->Fill(recoParticle.range(), pionResolution, weight);
                }
            }
            else if (style == PlottingHelper::Muon)
            {
                const auto muonRecoMomentumRange = AnalysisHelper::GetMuonMomentumFromRange(recoParticle.range());
                const auto muonRecoMomentumMCS = AnalysisHelper::GetMuonMomentumFromMCS(recoParticle);
                const auto muonResolutionRange = (muonRecoMomentumRange - muonMomentum) / muonMomentum;
                const auto muonResolutionMCS = (muonRecoMomentumMCS - muonMomentum) / muonMomentum;

                // Fill the muon momentum resolution plots
                muonMomentumResolutionPlots.at("all_range")->Fill(muonResolutionRange, weight);
                muonMomentumResolutionPlots.at("all_MCS")->Fill(muonResolutionMCS, weight);

                if (AnalysisHelper::IsContained(muon))
                {
                    muonMomentumResolutionPlots.at("contained_range")->Fill(muonResolutionRange, weight);
                    muonMomentumResolutionPlots.at("contained_MCS")->Fill(muonResolutionMCS, weight);

                    pMuonMomVsRange->Fill(recoParticle.range(), muonMomentum, weight);
                    pMuonResVsRange->Fill(recoParticle.range(), muonResolutionRange, weight);
                }
            }
            else if (style == PlottingHelper::Proton)
            {
                // Check for a truth-matched true proton.
                // If not matched to a true proton, skip this particle
                const auto &proton = AnalysisHelper::GetBestMatchedTruthParticle(recoParticle, truthParticles);
                if (proton.pdgCode()!=2212) continue;
                const auto protonMomentum = proton.momentum();

                const auto protonRecoMomentumRange = AnalysisHelper::GetProtonMomentumFromRange(recoParticle.range());
                const auto protonResolutionRange = (protonRecoMomentumRange - protonMomentum) / protonMomentum;

                const auto protonRecoMomentumRangeLarsoft = AnalysisHelper::GetProtonMomentumFromRangeLarsoft(recoParticle.range());
                const auto protonResolutionRangeLarsoft = (protonRecoMomentumRangeLarsoft - protonMomentum) / protonMomentum;

                // Fill the muon momentum resolution plots
                protonMomentumResolutionPlots.at("all_range")->Fill(protonResolutionRange, weight);

                protonMomentumResolutionPlots.at("all_range_larsoft")->Fill(protonResolutionRangeLarsoft, weight);

                if (AnalysisHelper::IsContained(proton))
                {
                    protonMomentumResolutionPlots.at("contained_range")->Fill(protonResolutionRange, weight);
                    protonMomentumResolutionPlots.at("contained_range_larsoft")->Fill(protonResolutionRangeLarsoft, weight);

                    pProtonMomVsRange->Fill(recoParticle.range(), protonMomentum, weight);
                    pProtonResVsRange->Fill(recoParticle.range(), protonResolutionRange, weight);

                    pProtonResVsRangeLarsoft->Fill(recoParticle.range(), protonResolutionRangeLarsoft, weight);
                }
            }
        }

        // ATTN in the backend code, we return "Other" as the end state for anything that's not a pi+... But here we can have pi- in the
        // signal samples - for now, this check restricts these plots to only pi+, but if we were to change the backend and re-run we could
        // make it for pi- as well.
        if (pion.pdgCode() != 211)
            continue;

        // Add this pion to the map
        totalWeight += weight;
        auto iter = pionEndStateMap.find(pionEndState);
        if (iter == pionEndStateMap.end())
        {
            pionEndStateMap.emplace(pionEndState, weight);
        }
        else
        {
            iter->second += weight;
        }

        if (pion.isStopping())
        {
            auto stoppingIter = stoppingPionEndStateMap.find(pionEndState);
            if (stoppingIter == stoppingPionEndStateMap.end())
            {
                stoppingPionEndStateMap.emplace(pionEndState, weight);
            }
            else
            {
                stoppingIter->second += weight;
            }
        }

        // Add to the nScatters plots
        pNScattersHist->Fill(static_cast<float>(pion.nElasticScatters()), static_cast<float>(pion.nInelasticScatters()), weight);

        // Add to the momentum fraction lost plots
        for (unsigned int iScatter = 0; iScatter < pion.scatterMomentumFracsLost().size(); ++iScatter)
        {
            const auto isElastic = pion.scatterIsElastic().at(iScatter);
            const auto momentumFracLost = pion.scatterMomentumFracsLost().at(iScatter);

            if (isElastic)
            {
                pMomLostElastic->Fill(momentumFracLost, weight);
            }
            else
            {
                pMomLostInelastic->Fill(momentumFracLost, weight);
            }
        }

        // Set the style based on the end-state
        PlottingHelper::PlotStyle style = PlottingHelper::Default;

        switch (pionEndState)
        {
            // None
            case 0:
                style = PlottingHelper::Primary;
                break;
            // Decays to a muon
            case 1:
                style = PlottingHelper::Secondary;
                break;
            // Inelastic absorption
            case 2:
                style = PlottingHelper::Tertiary;
                break;
            // Charge exchange to pi0
            case 3:
                style = PlottingHelper::Quaternary;
                break;
            // Something else
            default:
                style = PlottingHelper::Other;
                break;
        }

        pionMomentumPlot_endState.Fill(pionMomentum, style, weight);
        pionEndMomentumPlot_endState.Fill(pionEndMomentum, style, weight);
    }

    if (totalWeight <= std::numeric_limits<float>::epsilon())
        throw std::logic_error("SecondaryInteractionsStudy - No pions from CC1Pi events found");

    // Make the momentum plots
    pionMomentumPlot_endState.SaveAsStacked("secondaryInteractions_" + prefix + "pionMomentum_endStates");
    pionEndMomentumPlot_endState.SaveAsStacked("secondaryInteractions_" + prefix + "pionEndMomentum_endStates");

    // Setup the output canvas for the remaining plots
    auto pCanvas = PlottingHelper::GetCanvas();

    // Make the momentum resolution plots
    PlottingHelper::SetLineStyle(pionMomentumResolutionPlots.at("all"), PlottingHelper::Primary);
    PlottingHelper::SetLineStyle(pionMomentumResolutionPlots.at("contained"), PlottingHelper::Secondary);
    PlottingHelper::SetLineStyle(pionMomentumResolutionPlots.at("contained-stopping"), PlottingHelper::Tertiary);
    PlottingHelper::SetLineStyle(pionMomentumResolutionPlots.at("contained-stopping-noScatters"), PlottingHelper::Quaternary);

    pionMomentumResolutionPlots.at("all")->Draw("hist");
    pionMomentumResolutionPlots.at("contained")->Draw("hist same");
    pionMomentumResolutionPlots.at("contained-stopping")->Draw("hist same");
    pionMomentumResolutionPlots.at("contained-stopping-noScatters")->Draw("hist same");
    PlottingHelper::SaveCanvas(pCanvas, "secondaryInteractions_" + prefix + "pionMomentumResolution_range");

    PlottingHelper::SetLineStyle(muonMomentumResolutionPlots.at("all_range"), PlottingHelper::Primary);
    PlottingHelper::SetLineStyle(muonMomentumResolutionPlots.at("contained_range"), PlottingHelper::Secondary);
    muonMomentumResolutionPlots.at("all_range")->Draw("hist");
    muonMomentumResolutionPlots.at("contained_range")->Draw("hist same");
    PlottingHelper::SaveCanvas(pCanvas, "secondaryInteractions_" + prefix + "muonMomentumResolution_range");

    PlottingHelper::SetLineStyle(muonMomentumResolutionPlots.at("all_MCS"), PlottingHelper::Primary);
    PlottingHelper::SetLineStyle(muonMomentumResolutionPlots.at("contained_MCS"), PlottingHelper::Secondary);
    muonMomentumResolutionPlots.at("all_MCS")->Draw("hist");
    muonMomentumResolutionPlots.at("contained_MCS")->Draw("hist same");
    PlottingHelper::SaveCanvas(pCanvas, "secondaryInteractions_" + prefix + "muonMomentumResolution_MCS");

    PlottingHelper::SetLineStyle(protonMomentumResolutionPlots.at("all_range"), PlottingHelper::Primary);
    PlottingHelper::SetLineStyle(protonMomentumResolutionPlots.at("contained_range"), PlottingHelper::Secondary);
    protonMomentumResolutionPlots.at("all_range")->Draw("hist");
    protonMomentumResolutionPlots.at("contained_range")->Draw("hist same");
    PlottingHelper::SaveCanvas(pCanvas, "secondaryInteractions_" + prefix + "protonMomentumResolution_range");

    PlottingHelper::SetLineStyle(protonMomentumResolutionPlots.at("all_range_larsoft"), PlottingHelper::Primary);
    PlottingHelper::SetLineStyle(protonMomentumResolutionPlots.at("contained_range_larsoft"), PlottingHelper::Secondary);
    protonMomentumResolutionPlots.at("all_range_larsoft")->Draw("hist");
    protonMomentumResolutionPlots.at("contained_range_larsoft")->Draw("hist same");
    PlottingHelper::SaveCanvas(pCanvas, "secondaryInteractions_" + prefix + "protonMomentumResolution_range_larsoft");

    // Save the 2D resolution plots
    pPionMomVsRange->Draw("colz");
    auto pFuncPion = AnalysisHelper::GetRangeToMomentumFunctionPion();
    pFuncPion->Draw("same");
    PlottingHelper::SaveCanvas(pCanvas, "secondaryInteractions_" + prefix + "pionMomentumVsRange_golden");
    pPionResVsRange->Draw("colz");
    PlottingHelper::SaveCanvas(pCanvas, "secondaryInteractions_" + prefix + "pionMomentumResolutionVsRange_golden");

    pMuonMomVsRange->Draw("colz");
    auto pFuncMuon = AnalysisHelper::GetRangeToMomentumFunctionMuon();
    pFuncMuon->Draw("same");
    PlottingHelper::SaveCanvas(pCanvas, "secondaryInteractions_" + prefix + "muonMomentumVsRange_contained");
    pMuonResVsRange->Draw("colz");
    PlottingHelper::SaveCanvas(pCanvas, "secondaryInteractions_" + prefix + "muonMomentumResolutionVsRange_contained");

    pProtonMomVsRange->Draw("colz");
    auto pFuncProton = AnalysisHelper::GetRangeToMomentumFunctionProton();
    pFuncProton->Draw("same");
    PlottingHelper::SaveCanvas(pCanvas, "secondaryInteractions_" + prefix + "protonMomentumVsRange_contained");
    pProtonResVsRange->Draw("colz");
    PlottingHelper::SaveCanvas(pCanvas, "secondaryInteractions_" + prefix + "protonMomentumResolutionVsRange_contained");

    // LArsoft momentum function is hard to convert into a single function in momentum, so instead just make a TGraph with the value in each bin
    pProtonMomVsRange->Draw("colz");
    TGraph *pFuncProtonLarsoft = new TGraph(pProtonMomVsRange->GetXaxis()->GetNbins());
    for (auto i_x=0; i_x<pProtonMomVsRange->GetXaxis()->GetNbins(); i_x++){
        auto range_val = pProtonMomVsRange->GetXaxis()->GetBinCenter(i_x+1);
        auto mom_val = AnalysisHelper::GetProtonMomentumFromRangeLarsoft(range_val);
        pFuncProtonLarsoft->SetPoint(i_x,range_val,mom_val);
    }
    pFuncProtonLarsoft->SetLineWidth(2);
    pFuncProtonLarsoft->Draw("same L");
    PlottingHelper::SaveCanvas(pCanvas, "secondaryInteractions_" + prefix + "protonMomentumVsRange_contained_larsoft");
    pProtonResVsRangeLarsoft->Draw("colz");
    PlottingHelper::SaveCanvas(pCanvas, "secondaryInteractions_" + prefix + "protonMomentumResolutionVsRange_contained_larsoft");

    // Make the scatters plot
    pNScattersHist->Scale(100.f / pNScattersHist->GetEntries());

    // Set the bin labels
    for (unsigned int iBin = 1; iBin <= static_cast<unsigned int>(pNScattersHist->GetNbinsX()); ++iBin)
    {
        pNScattersHist->GetXaxis()->SetBinLabel(iBin, std::to_string(iBin - 1).c_str());
        pNScattersHist->GetYaxis()->SetBinLabel(iBin, std::to_string(iBin - 1).c_str());
    }

    gStyle->SetPaintTextFormat("2.2f");
    pNScattersHist->SetMarkerSize(2.2); // Text size
    pNScattersHist->Draw("colz text");
    PlottingHelper::SaveCanvas(pCanvas, "secondaryInteractions_" + prefix + "nScatters");

    // Make the momentum fraction lost plots
    PlottingHelper::SetLineStyle(pMomLostElastic, PlottingHelper::Primary);
    PlottingHelper::SetLineStyle(pMomLostInelastic, PlottingHelper::Secondary);
    pMomLostElastic->Scale(1.f / pMomLostElastic->GetEntries());
    pMomLostInelastic->Scale(1.f / pMomLostInelastic->GetEntries());

    pMomLostElastic->Draw("hist");
    pMomLostInelastic->Draw("hist same");
    PlottingHelper::SaveCanvas(pCanvas, "secondaryInteractions_" + prefix + "scatterMomentumFractionLost");

    // Print the end-state table
    FormattingHelper::Table table({"End state", "Events", "Fraction"});
    for (const auto &entry : pionEndStateMap)
    {
        const auto &endState = entry.first;
        const auto weight = entry.second;
        const auto fraction = weight / totalWeight;

        table.AddEmptyRow();
        table.SetEntry("End state", endState);
        table.SetEntry("Events", weight);
        table.SetEntry("Fraction", fraction);
    }

    FormattingHelper::PrintLine();
    std::cout << "All CC1pi pions" << std::endl;
    FormattingHelper::PrintLine();
    table.WriteToFile("secondaryInteractions_" + prefix + "pionEndStates.md");

    FormattingHelper::Table stoppingTable({"End state", "Events", "Fraction"});
    for (const auto &entry : stoppingPionEndStateMap)
    {
        const auto &endState = entry.first;
        const auto weight = entry.second;
        const auto fraction = weight / totalWeight;

        stoppingTable.AddEmptyRow();
        stoppingTable.SetEntry("End state", endState);
        stoppingTable.SetEntry("Events", weight);
        stoppingTable.SetEntry("Fraction", fraction);
    }

    FormattingHelper::PrintLine();
    std::cout << "All CC1pi stopping pions" << std::endl;
    FormattingHelper::PrintLine();
    stoppingTable.WriteToFile("secondaryInteractions_" + prefix + "stoppingPionEndStates.md");
}

} // namespace ubcc1pi_macros
