#include "ubcc1pi_standalone/Objects/FileReader.h"

#include "ubcc1pi_standalone/Helpers/PlottingHelper.h"
#include "ubcc1pi_standalone/Helpers/AnalysisHelper.h"

using namespace ubcc1pi;

int PlotInputVariables(const std::string &overlayFileName, const float overlayWeight, const std::string &dataEXTFileName, const float extWeight, const std::string &dataBNBFileName, const bool cleanSignalOnly)
{
    // Number of hits
    PlottingHelper::ParticlePlot plot_nHitsU("Number of collection (U) plane hits", 100, 0, 1000);
    PlottingHelper::ParticlePlot plot_nHitsV("Number of collection (V) plane hits", 100, 0, 1000);
    PlottingHelper::ParticlePlot plot_nHitsW("Number of collection (W) plane hits", 100, 0, 1000);
    PlottingHelper::ParticlePlot plot_nDescendents("Number of descendents", 6, 0, 6);
    PlottingHelper::ParticlePlot plot_nDescendentHits("Number of descendent hits", 100, 0, 200);
    PlottingHelper::ParticlePlot plot_nDescendentHitsNonZero("Number of descendent hits", 100, 1, 200);
    PlottingHelper::ParticlePlot plot_nHitsInLargestDescendent("nHitsInLargestDescendent", 100, 1, 200);
    PlottingHelper::ParticlePlot plot_nSpacePointsNearEnd("Number of spacepoints near track end", 60, 0, 120);
    PlottingHelper::ParticlePlot plot_wiggliness("Wiggliness", 60, 0, 0.005);
    PlottingHelper::ParticlePlot plot_trackScore("Track score", 100, 0, 1);
    PlottingHelper::ParticlePlot plot_longitudinalVertexDist("Longitudinal vertex dist", 100, 0, 10);
    PlottingHelper::ParticlePlot plot_transverseVertexDist("Transverse vertex dist", 100, 0, 10);
    PlottingHelper::ParticlePlot plot_truncatedMeandEdx("Truncated Mean dEdx", 100, 0, 10);
    PlottingHelper::ParticlePlot plot_logLikelihoodpMIP("log(L_p / L_MIP)", 100, -8, 8);
    PlottingHelper::ParticlePlot plot_logLikelihoodpiMIP("log(L_pi / L_MIP)", 100, -4, 8);
    PlottingHelper::ParticlePlot plot_forwardProton("Proton forward likelihood", 100, 0.3, 0.7);
    PlottingHelper::ParticlePlot plot_forwardMuon("Muon forward likelihood", 100, 0.3, 0.7);

    for (const auto fileName : {dataEXTFileName, dataBNBFileName, overlayFileName})
    {
        std::cout << "Reading input file: " << fileName << std::endl;

        const bool isBNBData = (fileName == dataBNBFileName);
        const bool isOverlay = (fileName == overlayFileName);
        const bool isEXTData = (fileName == dataEXTFileName);

        const float bodgeFactor = 1.273f; // ATTN this factor is a normalisation added so we can compare the shape of the distributions, can't exist in the final result!
        float weight = 1.f;
        if (isOverlay) weight = overlayWeight * bodgeFactor;
        if (isEXTData) weight = extWeight * bodgeFactor;

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
                
            if (cleanSignalOnly)
            {
                // Only use MC events 
                if (!pEvent->metadata.hasTruthInfo())
                    continue;

                // Only use signal events
                if (!AnalysisHelper::IsTrueCC1Pi(pEvent))
                    continue;
            }


            const auto truthParticles = pEvent->truth.particles;
            const auto recoParticles = pEvent->reco.particles;

            for (const auto &particle : recoParticles)
            {
                bool isContained = false;
                try
                {
                    isContained = AnalysisHelper::IsContained(particle);
                }
                catch (const std::invalid_argument &) {}

                if (!isContained)
                    continue;

                if (cleanSignalOnly)
                {
                    try
                    {
                        const auto truthParticleIndex = AnalysisHelper::GetBestMatchedTruthParticleIndex(particle, truthParticles);
                        const auto completeness = particle.truthMatchCompletenesses().at(truthParticleIndex);

                        if (completeness < 0.5f)
                            continue;
                    }
                    catch (const std::logic_error &)
                    {
                        continue;
                    }
                }

                const auto particleStyle = isBNBData ? PlottingHelper::BNBData : PlottingHelper::GetParticleStyle(particle, truthParticles);

                if (particle.nHitsU.IsSet())
                    plot_nHitsU.Fill(particle.nHitsU(), particleStyle, weight);

                if (particle.nHitsV.IsSet())
                    plot_nHitsV.Fill(particle.nHitsV(), particleStyle, weight);

                if (particle.nHitsW.IsSet())
                    plot_nHitsW.Fill(particle.nHitsW(), particleStyle, weight);

                if (particle.nDescendents.IsSet())
                    plot_nDescendents.Fill(particle.nDescendents(), particleStyle, weight);

                if (particle.nDescendentHitsU.IsSet() && particle.nDescendentHitsV.IsSet() && particle.nDescendentHitsW.IsSet())
                {
                    plot_nDescendentHits.Fill(particle.nDescendentHitsU() + particle.nDescendentHitsV() + particle.nDescendentHitsW(), particleStyle, weight);
                    plot_nDescendentHitsNonZero.Fill(particle.nDescendentHitsU() + particle.nDescendentHitsV() + particle.nDescendentHitsW(), particleStyle, weight);
                }

                if (particle.nHitsInLargestDescendent.IsSet())
                    plot_nHitsInLargestDescendent.Fill(particle.nHitsInLargestDescendent(), particleStyle, weight);

                if (particle.nSpacePointsNearEnd.IsSet())
                    plot_nSpacePointsNearEnd.Fill(particle.nSpacePointsNearEnd(), particleStyle, weight);

                if (particle.wiggliness.IsSet())
                    plot_wiggliness.Fill(particle.wiggliness(), particleStyle, weight);

                if (particle.trackScore.IsSet())
                    plot_trackScore.Fill(particle.trackScore(), particleStyle, weight);

                if (particle.longitudinalVertexDist.IsSet())
                    plot_longitudinalVertexDist.Fill(particle.longitudinalVertexDist(), particleStyle, weight);

                if (particle.transverseVertexDist.IsSet())
                    plot_transverseVertexDist.Fill(particle.transverseVertexDist(), particleStyle, weight);

                if (particle.truncatedMeandEdx.IsSet())
                    plot_truncatedMeandEdx.Fill(particle.truncatedMeandEdx(), particleStyle, weight);

                // TODO refactor with new analysis helper functions
                if (particle.likelihoodForwardProton.IsSet() && particle.likelihoodMIP.IsSet())
                {
                    const float denom = particle.likelihoodMIP();
                    if (denom > std::numeric_limits<float>::epsilon())
                        plot_logLikelihoodpMIP.Fill(std::log(particle.likelihoodForwardProton() / denom), particleStyle, weight);
                }

                if (particle.likelihoodForwardPion.IsSet() && particle.likelihoodMIP.IsSet())
                {
                    const float denom = particle.likelihoodMIP();
                    if (denom > std::numeric_limits<float>::epsilon())
                        plot_logLikelihoodpiMIP.Fill(std::log(particle.likelihoodForwardPion() / denom), particleStyle, weight);
                }

                if (particle.likelihoodForwardProton.IsSet() && particle.likelihoodBackwardProton.IsSet())
                {
                    const float denom = std::exp(particle.likelihoodForwardProton()) + std::exp(particle.likelihoodBackwardProton());
                    if (denom > std::numeric_limits<float>::epsilon())
                        plot_forwardProton.Fill(std::exp(particle.likelihoodForwardProton()) / denom, particleStyle, weight);
                }

                if (particle.likelihoodForwardMuon.IsSet() && particle.likelihoodBackwardMuon.IsSet())
                {
                    const float denom = std::exp(particle.likelihoodForwardMuon()) + std::exp(particle.likelihoodBackwardMuon());
                    if (denom > std::numeric_limits<float>::epsilon())
                        plot_forwardMuon.Fill(std::exp(particle.likelihoodForwardMuon()) / denom, particleStyle, weight);
                }
            }
        }
    }

    plot_nHitsU.SaveAs("nHitsU");
    plot_nHitsV.SaveAs("nHitsV");
    plot_nHitsW.SaveAs("nHitsW");
    plot_nDescendents.SaveAs("nDescendents");
    plot_nDescendentHits.SaveAs("nDescendentHits");
    plot_nDescendentHitsNonZero.SaveAs("nDescendentHitsNonZero");
    plot_nHitsInLargestDescendent.SaveAs("nHitsInLargestDescendent");
    plot_nSpacePointsNearEnd.SaveAs("nSpacePointsNearEnd");
    plot_wiggliness.SaveAs("wiggliness");
    plot_trackScore.SaveAs("trackScore");
    plot_longitudinalVertexDist.SaveAs("longitudinalVertexDist");
    plot_transverseVertexDist.SaveAs("transverseVertexDist");
    plot_truncatedMeandEdx.SaveAs("truncatedMeandEdx");
    plot_logLikelihoodpMIP.SaveAs("logLikelihoodpMIP");
    plot_logLikelihoodpiMIP.SaveAs("logLikelihoodpiMIP");
    plot_forwardProton.SaveAs("forwardProton");
    plot_forwardMuon.SaveAs("forwardMuon");
    
    plot_nHitsU.SaveAsStacked("nHitsU_stacked");
    plot_nHitsV.SaveAsStacked("nHitsV_stacked");
    plot_nHitsW.SaveAsStacked("nHitsW_stacked");
    plot_nDescendents.SaveAsStacked("nDescendents_stacked");
    plot_nDescendentHits.SaveAsStacked("nDescendentHits_stacked");
    plot_nDescendentHitsNonZero.SaveAsStacked("nDescendentHitsNonZero_stacked");
    plot_nHitsInLargestDescendent.SaveAsStacked("nHitsInLargestDescendent_stacked");
    plot_nSpacePointsNearEnd.SaveAsStacked("nSpacePointsNearEnd_stacked");
    plot_wiggliness.SaveAsStacked("wiggliness_stacked");
    plot_trackScore.SaveAsStacked("trackScore_stacked");
    plot_longitudinalVertexDist.SaveAsStacked("longitudinalVertexDist_stacked");
    plot_transverseVertexDist.SaveAsStacked("transverseVertexDist_stacked");
    plot_truncatedMeandEdx.SaveAsStacked("truncatedMeandEdx_stacked");
    plot_logLikelihoodpMIP.SaveAsStacked("logLikelihoodpMIP_stacked");
    plot_logLikelihoodpiMIP.SaveAsStacked("logLikelihoodpiMIP_stacked");
    plot_forwardProton.SaveAsStacked("forwardProton_stacked");
    plot_forwardMuon.SaveAsStacked("forwardMuon_stacked");

    return 0;
}
