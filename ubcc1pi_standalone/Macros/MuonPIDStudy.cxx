/**
 *  @file  ubcc1pi_standalone/Macros/MuonPIDStudy.cxx
 *
 *  @brief The implementation file of the MuonPIDStudy macro
 */

#include "ubcc1pi_standalone/Macros/Macros.h"

#include "ubcc1pi_standalone/Objects/FileReader.h"

#include "ubcc1pi_standalone/Helpers/AnalysisHelper.h"
#include "ubcc1pi_standalone/Helpers/FormattingHelper.h"
#include "ubcc1pi_standalone/Helpers/PlottingHelper.h"
#include "ubcc1pi_standalone/Helpers/BDTHelper.h"

#include <string>
#include <vector>
#include <utility>
#include <map>

using namespace ubcc1pi;

namespace ubcc1pi_macros
{

void MuonPIDStudy(const Config &config)
{
    // Open the file
    FileReader reader(config.files.overlaysFileName);
    auto pEvent = reader.GetBoundEventAddress();
    const auto nEvents = reader.GetNumberOfEvents();
    
    // Setup the muon BDT
    const auto featureNames = BDTHelper::ParticleBDTFeatureNames;
    BDTHelper::BDT muonBDT("muon", featureNames); 

    // Setup the counters
    float nEventsTotal = 0.f;
    float nEventsHasTrack = 0.f;
    float nEventsMax1Escaping = 0.f;
    float nEvents0Escaping = 0.f;
    float nEvents0EscapingHasBDT = 0.f;

    float nEventsUsedCCInclusiveCandidate = 0.f;
    float nEventsUsedEscapingCandidate = 0.f;
    float nEventsUsedBDTCandidate = 0.f;
    
    float nEventsUsedCCInclusiveCandidate_ccincCorrect = 0.f;
    float nEventsUsedEscapingCandidate_ccincCorrect = 0.f;
    float nEventsUsedBDTCandidate_ccincCorrect = 0.f;
    float nEventsUsedCCInclusiveCandidate_correct = 0.f;
    float nEventsUsedEscapingCandidate_correct = 0.f;
    float nEventsUsedBDTCandidate_correct = 0.f;
    

    std::map<PlottingHelper::PlotStyle, float> trueParticleTypeToCountMap_ccinc;
    std::map<PlottingHelper::PlotStyle, float> trueParticleTypeToCountMap_ours;

    // Loop over the events
    for (unsigned int i = 0; i < nEvents; ++i)
    {
        AnalysisHelper::PrintLoadingBar(i, nEvents);
        reader.LoadEvent(i);

        // Only consider true CC1pi events
        if (!AnalysisHelper::IsTrueCC1Pi(pEvent, config.global.useAbsPdg))
            continue;

        // Only consider events that pass the CC inclusive filter
        if (!pEvent->reco.passesCCInclusive())
            continue;

        // Flags describing which end-state we are in
        bool usedCCInclusiveCandidate = false;
        bool usedEscapingCandidate    = false;
        bool usedBDTCandidate         = false;

        const auto weight = AnalysisHelper::GetNominalEventWeight(pEvent);
        nEventsTotal += weight;
        
        const auto recoParticles = pEvent->reco.particles;

        // Start of by assuming that the CC inclusive muon candidate is correct
        unsigned int ccInclusiveMuonIndex = std::numeric_limits<unsigned int>::max();
        bool foundCCInclusiveMuon = false;
        for (unsigned int index = 0; index < recoParticles.size(); ++index)
        {
            const auto &particle = recoParticles.at(index);
            if (!particle.isCCInclusiveMuonCandidate())
                continue;
            
            if (foundCCInclusiveMuon)
                throw std::logic_error("MuonPIDStudy - found multiple CC inclusive muon candidates");

            foundCCInclusiveMuon = true;
            ccInclusiveMuonIndex = index;
        }

        if (!foundCCInclusiveMuon)
            throw std::logic_error("MuonPIDStudy - couldn't find the CC inclusive muon candidate");

        unsigned int muonIndex = ccInclusiveMuonIndex;


        // Now decide if we should choose a different particle as the muon
        // Determine if the event has any particles with a fitted track
        std::vector<unsigned int> trackFitParticleIndices;
        for (unsigned int index = 0; index < recoParticles.size(); ++index)
        {
            const auto &particle = recoParticles.at(index);

            if (AnalysisHelper::HasTrackFit(particle))
                trackFitParticleIndices.push_back(index);
        }
        
        // Insist at least one has a fitted track
        if (!trackFitParticleIndices.empty())
        {
            nEventsHasTrack += weight;

            // Now find the particles that escape
            std::vector<unsigned int> containedIndices, uncontainedIndices;
            for (const auto &index : trackFitParticleIndices)
            {
                const auto &particle = recoParticles.at(index);

                auto &vector = AnalysisHelper::IsContained(particle) ? containedIndices : uncontainedIndices;
                vector.push_back(index);
            }

            // Insist that at most one particle escapes
            if (uncontainedIndices.size() <= 1)
            {
                nEventsMax1Escaping += weight;

                // If there is exactly one escaping particle, then call is the muon
                if (uncontainedIndices.size() == 1)
                {
                    muonIndex = uncontainedIndices.front();
                    usedEscapingCandidate = true;
                }
                else
                {
                    // Here there are no escaping particles
                    nEvents0Escaping += weight;

                    // Choose the contained particle with the best muon BDT response
                    bool foundMuonCandidate = false;
                    float maxMuonBDTResponse = -std::numeric_limits<float>::max();
                    for (const auto &index : containedIndices)
                    {
                        const auto &particle = recoParticles.at(index);

                        std::vector<float> features;
                        const auto hasFeatures = BDTHelper::GetBDTFeatures(particle, featureNames, features);

                        if (!hasFeatures)
                            continue;
                    
                        const auto muonBDTResponse = muonBDT.GetResponse(features);
                        if (muonBDTResponse < maxMuonBDTResponse)
                            continue;

                        maxMuonBDTResponse = muonBDTResponse;
                        foundMuonCandidate = true;
                        muonIndex = index;
                    }

                    if (foundMuonCandidate)
                    {
                        usedBDTCandidate = true;
                        nEvents0EscapingHasBDT += weight;
                    }
                    else
                    {
                        usedCCInclusiveCandidate = true;
                    }
                }
            }
            else
            {
                usedCCInclusiveCandidate = true;
            }
        }
        else
        {
            usedCCInclusiveCandidate = true;
        }


        // Sanity check the result
        if ((usedCCInclusiveCandidate && (usedEscapingCandidate || usedBDTCandidate)) ||
            (usedEscapingCandidate && (usedBDTCandidate || usedCCInclusiveCandidate)) ||
            (usedBDTCandidate && (usedCCInclusiveCandidate || usedEscapingCandidate)) ||
            (!usedCCInclusiveCandidate && !usedEscapingCandidate && !usedBDTCandidate))
        {
            throw std::logic_error("MuonPIDStudy - sanity check failed!");
        }

    
        // Now determine if true origin of the muon candidate is correct
        const auto truthParticles = pEvent->truth.particles;

        const auto ccInclusiveMuon = recoParticles.at(ccInclusiveMuonIndex);
        const auto ccInclusiveMuonStyle = PlottingHelper::GetPlotStyle(ccInclusiveMuon, AnalysisHelper::Overlay, truthParticles, false, config.global.useAbsPdg);
        const auto ccInclusiveMuonCorrect = (ccInclusiveMuonStyle == PlottingHelper::Muon);
            
        auto iter_ccinc = trueParticleTypeToCountMap_ccinc.find(ccInclusiveMuonStyle);
        if (iter_ccinc == trueParticleTypeToCountMap_ccinc.end())
        {
            trueParticleTypeToCountMap_ccinc.emplace(ccInclusiveMuonStyle, weight);
        }
        else
        {
            iter_ccinc->second += weight;
        }
        
        const auto ourMuon = recoParticles.at(muonIndex);
        const auto ourMuonStyle = PlottingHelper::GetPlotStyle(ourMuon, AnalysisHelper::Overlay, truthParticles, false, config.global.useAbsPdg);
        const auto ourMuonCorrect = (ourMuonStyle == PlottingHelper::Muon);
        
        auto iter_ours = trueParticleTypeToCountMap_ours.find(ourMuonStyle);
        if (iter_ours == trueParticleTypeToCountMap_ours.end())
        {
            trueParticleTypeToCountMap_ours.emplace(ourMuonStyle, weight);
        }
        else
        {
            iter_ours->second += weight;
        }

    
        if (usedCCInclusiveCandidate)
        {
            nEventsUsedCCInclusiveCandidate += weight;

            nEventsUsedCCInclusiveCandidate_ccincCorrect += (ccInclusiveMuonCorrect ? weight : 0.f);
            nEventsUsedCCInclusiveCandidate_correct += (ourMuonCorrect ? weight : 0.f);
        }

        if (usedEscapingCandidate)
        {
            nEventsUsedEscapingCandidate += weight;
            
            nEventsUsedEscapingCandidate_ccincCorrect += (ccInclusiveMuonCorrect ? weight : 0.f);
            nEventsUsedEscapingCandidate_correct += (ourMuonCorrect ? weight : 0.f);
        }

        if (usedBDTCandidate)
        {
            nEventsUsedBDTCandidate += weight;
            
            nEventsUsedBDTCandidate_ccincCorrect += (ccInclusiveMuonCorrect ? weight : 0.f);
            nEventsUsedBDTCandidate_correct += (ourMuonCorrect ? weight : 0.f);
        }
    }

    FormattingHelper::Table tableSummary({"Method", "Particle type", "Events", "Fraction"});

    std::map<std::string, std::map<PlottingHelper::PlotStyle, float> > methodToCountMap;
    methodToCountMap.emplace("CC inclusive", trueParticleTypeToCountMap_ccinc);
    methodToCountMap.emplace("Ours", trueParticleTypeToCountMap_ours);

    for (const auto &methodMap : methodToCountMap)
    {
        const auto &method = methodMap.first;
        const auto &map = methodMap.second;
        
        for (const auto &entry : map)
        {
            const auto style = entry.first;
            const auto nEvents = entry.second;

            std::string particleType;
            switch (style)
            {
                case PlottingHelper::Muon:
                    particleType = "Muon";
                    break;
                case PlottingHelper::Proton:
                    particleType = "Proton";
                    break;
                case PlottingHelper::GoldenPion:
                    particleType = "Golden pion";
                    break;
                case PlottingHelper::NonGoldenPion:
                    particleType = "Non golden pion";
                    break;
                case PlottingHelper::External:
                    particleType = "External";
                    break;
                default:
                    particleType = "Other (" + std::to_string(style) + ")";
                    break;
            }

            tableSummary.AddEmptyRow();
            tableSummary.SetEntry("Method", method);
            tableSummary.SetEntry("Particle type", particleType);
            tableSummary.SetEntry("Events", nEvents);
            tableSummary.SetEntry("Fraction", nEvents / nEventsTotal);
        }
    }

    tableSummary.WriteToFile("muonPIDStudy_summary.md");



    FormattingHelper::Table table({"Label", "Events"});

    table.AddEmptyRow();
    table.SetEntry("Label", "All CC1Pi events passing CC inclusive");
    table.SetEntry("Events", nEventsTotal);
    
    table.AddEmptyRow();
    table.SetEntry("Label", "... with at least one fitted track");
    table.SetEntry("Events", nEventsHasTrack);
    
    table.AddEmptyRow();
    table.SetEntry("Label", "... and <= 1 escaping particle");
    table.SetEntry("Events", nEventsMax1Escaping);
    
    table.AddEmptyRow();
    table.SetEntry("Label", "... and 0 escaping particles");
    table.SetEntry("Events", nEvents0Escaping);
    
    table.AddEmptyRow();
    table.SetEntry("Label", "... and has BDT response");
    table.SetEntry("Events", nEvents0EscapingHasBDT);
    
    table.AddEmptyRow();
    table.AddEmptyRow();
    table.SetEntry("Label", "Events using CC inclusive candidate");
    table.SetEntry("Events", nEventsUsedCCInclusiveCandidate);
    
    table.AddEmptyRow();
    table.SetEntry("Label", "... CC inclusive candidate is correct");
    table.SetEntry("Events", nEventsUsedCCInclusiveCandidate_ccincCorrect);
    
    table.AddEmptyRow();
    table.SetEntry("Label", "... our candidate is correct");
    table.SetEntry("Events", nEventsUsedCCInclusiveCandidate_correct);
    
    table.AddEmptyRow();
    table.SetEntry("Label", "Events using escaping candidate");
    table.SetEntry("Events", nEventsUsedEscapingCandidate);
    
    table.AddEmptyRow();
    table.SetEntry("Label", "... CC inclusive candidate is correct");
    table.SetEntry("Events", nEventsUsedEscapingCandidate_ccincCorrect);
    
    table.AddEmptyRow();
    table.SetEntry("Label", "... our candidate is correct");
    table.SetEntry("Events", nEventsUsedEscapingCandidate_correct);
    
    table.AddEmptyRow();
    table.SetEntry("Label", "Events using BDT candidate");
    table.SetEntry("Events", nEventsUsedBDTCandidate);
    
    table.AddEmptyRow();
    table.SetEntry("Label", "... CC inclusive candidate is correct");
    table.SetEntry("Events", nEventsUsedBDTCandidate_ccincCorrect);
    
    table.AddEmptyRow();
    table.SetEntry("Label", "... our candidate is correct");
    table.SetEntry("Events", nEventsUsedBDTCandidate_correct);

    table.WriteToFile("muonPIDStudy_breakdown.md");
}

} // namespace ubcc1pi_macros
