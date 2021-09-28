/**
 *  @file  ubcc1pi_standalone/Macros/PlotEBRequests.cxx
 *
 *  @brief The implementation file of the PlotEBRequests macro
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

using namespace ubcc1pi;

namespace ubcc1pi_macros
{

void PlotEBRequests(const Config &config)
{
    // Open the file
    FileReader reader(config.files.overlaysFileName);
    auto pEvent = reader.GetBoundEventAddress();
    const auto nEvents = reader.GetNumberOfEvents();

    //
    // Setup the data structures for the table
    //
    // std::map<int, float> pionEndStateMap;
    // std::map<int, float> stoppingPionEndStateMap;
    float totalWeight = 0.f;

    //
    // Setup the plots
    //

    // // Reco-true momentum resolution plots (index of map is a descriptor, eg. "contained")
    // std::map<std::string, std::shared_ptr<TH1F> > pionMomentumResolutionPlots, muonMomentumResolutionPlots;

    // pionMomentumResolutionPlots.emplace("all", new TH1F("pionRes_all", "", 50, -1.f, 1.5f));
    // pionMomentumResolutionPlots.emplace("contained", new TH1F("pionRes_contained", "", 50, -1.f, 1.5f));
    // pionMomentumResolutionPlots.emplace("contained-stopping", new TH1F("pionRes_contained-stopping", "", 50, -1.f, 1.5f));
    // pionMomentumResolutionPlots.emplace("contained-stopping-noScatters", new TH1F("pionRes_contained-stopping-noScatters", "", 50, -1.f, 1.5f));

    // muonMomentumResolutionPlots.emplace("all_range", new TH1F("muonRes_all_range", "", 50, -1.f, 0.7f));
    // muonMomentumResolutionPlots.emplace("contained_range", new TH1F("muonRes_contained_range", "", 50, -1.f, 0.7f));
    // muonMomentumResolutionPlots.emplace("all_MCS", new TH1F("muonRes_all_MCS", "", 50, -1.f, 0.7f));
    // muonMomentumResolutionPlots.emplace("contained_MCS", new TH1F("muonRes_contained_MCS", "", 50, -1.f, 0.7f));

    // Pion reco-momentum vs wiggliness
    const int bins = 50;
    std::shared_ptr<TH2F> pPionMomVsWiggliness_golden(new TH2F("pionMomVsWiggliness_golden", "", bins, 0.f, 0.28f, bins, 0.f, 1.5f));
    std::shared_ptr<TH2F> pPionMomVsWiggliness_nonGolden(new TH2F("pionMomVsWiggliness_nonGolden", "", bins, 0.f, 0.28f, bins, 0.f, 1.5f));
    std::shared_ptr<TH2F> pPionMomVsWiggliness_all(new TH2F("pionMomVsWiggliness_all", "", bins, 0.f, 0.28f, bins, 0.f, 1.5f));
    std::shared_ptr<TH2F> pMuonMomVsWiggliness(new TH2F("pionMuonVsWiggliness", "", bins, 0.f, 0.28f, bins, 0.f, 2.0f));
    
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
                    throw std::logic_error("PlotEBRequests - Found multiple pions in a signal event!");

                pionIndex = index;
                foundPion = true;
            }

            if (pdgCode == 13)
            {
                if (foundMuon)
                    throw std::logic_error("PlotEBRequests - Found multiple muons in a signal event!");

                muonIndex = index;
                foundMuon = true;
            }
        }

        if (!foundPion)
            throw std::logic_error("PlotEBRequests - Couldn't find the pion in a signal event!");

        if (!foundMuon)
            throw std::logic_error("PlotEBRequests - Couldn't find the muon in a signal event!");

        // Get the relevant information
        const auto &pion = truthParticles.at(pionIndex);
        const auto &muon = truthParticles.at(muonIndex);

        const auto pionMomentum = pion.momentum();
        // const auto pionEndMomentum = pion.endMomentum();
        // const auto pionEndState = pion.endState();

        const auto muonMomentum = muon.momentum();

        // Get the reco particles
        for (const auto &recoParticle : pEvent->reco.particles)
        {
            if (!AnalysisHelper::HasTrackFit(recoParticle))
                continue;

            const auto style = PlottingHelper::GetPlotStyle(recoParticle, AnalysisHelper::Overlay, truthParticles, false, config.global.useAbsPdg);

            if (style == PlottingHelper::GoldenPion || style == PlottingHelper::NonGoldenPion)
            {
                // const auto pionRecoMomentum = AnalysisHelper::GetPionMomentumFromRange(recoParticle.range());
                // const auto pionResolution = (pionRecoMomentum - pionMomentum) / pionMomentum;

                // Fill the pion momentum resolution plots
                // pionMomentumResolutionPlots.at("all")->Fill(pionResolution, weight);

                // if (AnalysisHelper::IsContained(pion))
                // {
                //     // pionMomentumResolutionPlots.at("contained")->Fill(pionResolution, weight);

                //     if (pion.isStopping())
                //     {
                //         // pionMomentumResolutionPlots.at("contained-stopping")->Fill(pionResolution, weight);

                //         if (pion.nElasticScatters() == 0 && pion.nInelasticScatters() == 0)
                //         {
                //             pionMomentumResolutionPlots.at("contained-stopping-noScatters")->Fill(pionResolution, weight);
                //         }
                //     }
                // }

                const auto pionRecoMomentum = AnalysisHelper::GetPionMomentumFromRange(recoParticle.range());
                if (AnalysisHelper::IsGolden(pion))
                {
                    pPionMomVsWiggliness_golden->Fill(recoParticle.wiggliness(), pionRecoMomentum, weight);
                    //pPionMomVsWiggliness_golden->Fill(recoParticle.wiggliness(), pionMomentum, weight);
                    // pPionResVsRange->Fill(recoParticle.range(), pionResolution, weight);
                }
                else
                {
                    pPionMomVsWiggliness_nonGolden->Fill(recoParticle.wiggliness(), pionRecoMomentum, weight);
                    // pPionMomVsWiggliness_nonGolden->Fill(recoParticle.wiggliness(), pionMomentum, weight);
                }

                pPionMomVsWiggliness_all->Fill(recoParticle.wiggliness(), pionRecoMomentum, weight);
            }
            else if (style == PlottingHelper::Muon)
            {
                const auto muonRecoMomentum = AnalysisHelper::GetMuonMomentumFromRange(recoParticle.range());
                pMuonMomVsWiggliness->Fill(recoParticle.wiggliness(), muonRecoMomentum, weight);
            }
        }
    }
    
    auto pCanvas = PlottingHelper::GetCanvas();
    // Save the 2D resolution plots
    pPionMomVsWiggliness_golden->Draw("colz");
    PlottingHelper::SaveCanvas(pCanvas, "EBRequest_recoPionMomentumVsWiggliness_golden");

    pPionMomVsWiggliness_nonGolden->Draw("colz");
    PlottingHelper::SaveCanvas(pCanvas, "EBRequest_recoPionMomentumVsWiggliness_nonGolden");

    pPionMomVsWiggliness_all->Draw("colz");
    PlottingHelper::SaveCanvas(pCanvas, "EBRequest_recoPionMomentumVsWiggliness_all");

    pMuonMomVsWiggliness->Draw("colz");
    PlottingHelper::SaveCanvas(pCanvas, "EBRequest_recoMuonMomentumVsWiggliness");
}

} // namespace ubcc1pi_macros