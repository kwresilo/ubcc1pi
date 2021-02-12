/**
 *  @file  ubcc1pi_standalone/Macros/PlotPionInputVariables.cxx
 *
 *  @brief The implementation file of the PlotPionInputVariables macro
 */

#include "ubcc1pi_standalone/Macros/Macros.h"

#include "ubcc1pi_standalone/Objects/FileReader.h"

#include "ubcc1pi_standalone/Helpers/PlottingHelper.h"
#include "ubcc1pi_standalone/Helpers/BDTHelper.h"
#include "ubcc1pi_standalone/Helpers/AnalysisHelper.h"
#include "ubcc1pi_standalone/Helpers/NormalisationHelper.h"

#include <TH2F.h>

using namespace ubcc1pi;

namespace ubcc1pi_macros
{

void PlotPionInputVariables(const Config &config)
{
    //
    // Setup the input files
    //
    std::vector< std::tuple<AnalysisHelper::SampleType, std::string, float> > inputData;

    inputData.emplace_back(AnalysisHelper::Overlay, config.files.overlaysFileName, NormalisationHelper::GetOverlaysNormalisation(config));
    inputData.emplace_back(AnalysisHelper::Dirt,    config.files.dirtFileName,     NormalisationHelper::GetDirtNormalisation(config));
    inputData.emplace_back(AnalysisHelper::DataEXT, config.files.dataEXTFileName,  NormalisationHelper::GetDataEXTNormalisation(config));
    inputData.emplace_back(AnalysisHelper::DataBNB, config.files.dataBNBFileName,  1.f);

    //
    // Set up the plots for each BDT feature
    //
    const auto featureNames = BDTHelper::ParticleBDTFeatureNames;
    const std::string yLabel = "Fraction of reco particles";
    std::vector< PlottingHelper::MultiPlot > plotVectorGolden, plotVectorNonGolden;

    for (const auto &featureName : featureNames)
    {
        if (featureName == "logBragg_pToMIP")
        {
            plotVectorGolden.emplace_back("log(L_p / L_MIP)", yLabel, 60, -9, 7, false);
            plotVectorNonGolden.emplace_back("log(L_p / L_MIP)", yLabel, 60, -9, 7, false);
            continue;
        }

        if (featureName == "logBragg_piToMIP")
        {
            plotVectorGolden.emplace_back("log(L_pi / L_MIP)", yLabel, 60, -4, 6, false);
            plotVectorNonGolden.emplace_back("log(L_pi / L_MIP)", yLabel, 60, -4, 6, false);
            continue;
        }

        if (featureName == "truncMeandEdx")
        {
            plotVectorGolden.emplace_back("Truncated Mean dEdx", yLabel, 60, 0, 10, false);
            plotVectorNonGolden.emplace_back("Truncated Mean dEdx", yLabel, 60, 0, 10, false);
            continue;
        }

        if (featureName == "protonForward")
        {
            plotVectorGolden.emplace_back("Proton forward likelihood", yLabel, 60, 0.42, 0.62, false);
            plotVectorNonGolden.emplace_back("Proton forward likelihood", yLabel, 60, 0.42, 0.62, false);
            continue;
        }

        if (featureName == "muonForward")
        {
            plotVectorGolden.emplace_back("Muon forward likelihood", yLabel, 60, 0.35, 0.65, false);
            plotVectorNonGolden.emplace_back("Muon forward likelihood", yLabel, 60, 0.35, 0.65, false);
            continue;
        }

        if (featureName == "nDescendents")
        {
            plotVectorGolden.emplace_back("Number of descendent particles", yLabel, 4, 0, 4, false);
            plotVectorNonGolden.emplace_back("Number of descendent particles", yLabel, 4, 0, 4, false);
            continue;
        }

        if (featureName == "nSpacePointsNearEnd")
        {
            plotVectorGolden.emplace_back("Number of spacepoints near track end", yLabel, 45, 0, 90, false);
            plotVectorNonGolden.emplace_back("Number of spacepoints near track end", yLabel, 45, 0, 90, false);
            continue;
        }

        if (featureName == "wiggliness")
        {
            // ATTN we use log-x for clarity, so use non-linear binning
            const auto binEdges = PlottingHelper::GenerateLogBinEdges(50, 1e-4, 0.08);
            plotVectorGolden.emplace_back("Wiggliness", yLabel, binEdges, false);
            plotVectorNonGolden.emplace_back("Wiggliness", yLabel, binEdges, false);

            continue;
        }

        if (featureName == "trackScore")
        {
            plotVectorGolden.emplace_back("Track score", yLabel, 30, 0, 1, false);
            plotVectorNonGolden.emplace_back("Track score", yLabel, 30, 0, 1, false);
            continue;
        }

        throw std::logic_error("PlotPionInputVariables - unknown feature: \"" + featureName + "\"");
    }

    // Setup the BDT outputs
    const auto goldenPionFeatureNames = BDTHelper::GoldenPionBDTFeatureNames;
    const auto protonFeatureNames = BDTHelper::ProtonBDTFeatureNames;
    const auto muonFeatureNames = BDTHelper::MuonBDTFeatureNames;

    PlottingHelper::MultiPlot muonBDTPlotGolden("Muon BDT response", yLabel, 40, -0.85f, 0.50f, false);
    PlottingHelper::MultiPlot protonBDTPlotGolden("Proton BDT response", yLabel, 40, -0.60f, 0.60f, false);
    PlottingHelper::MultiPlot goldenPionBDTPlotGolden("Golden pion BDT response", yLabel, 40, -0.8f, 0.4f, false);

    PlottingHelper::MultiPlot muonBDTPlotNonGolden("Muon BDT response", yLabel, 40, -0.85f, 0.50f, false);
    PlottingHelper::MultiPlot protonBDTPlotNonGolden("Proton BDT response", yLabel, 40, -0.60f, 0.60f, false);
    PlottingHelper::MultiPlot goldenPionBDTPlotNonGolden("Golden pion BDT response", yLabel, 40, -0.8f, 0.4f, false);

    std::shared_ptr<BDTHelper::BDT> pGoldenPionBDT, pProtonBDT, pMuonBDT;
    if (config.plotInputVariables.plotBDTResponses)
    {
        // Setup the BDTs
        pGoldenPionBDT = std::make_shared<BDTHelper::BDT>("goldenPion", goldenPionFeatureNames);
        pProtonBDT = std::make_shared<BDTHelper::BDT>("proton", protonFeatureNames);
        pMuonBDT = std::make_shared<BDTHelper::BDT>("muon", muonFeatureNames);
    }

    //
    // Fill the plots using overlay events
    //
    FileReader reader(config.files.overlaysFileName);
    auto pEvent = reader.GetBoundEventAddress();

    const auto normalisation = NormalisationHelper::GetOverlaysNormalisation(config);

    const auto nEvents = reader.GetNumberOfEvents();
    for (unsigned int i = 0; i < nEvents; ++i)
    {
        AnalysisHelper::PrintLoadingBar(i, nEvents);
        reader.LoadEvent(i);

        // Only use true CC1Pi events
        if (!AnalysisHelper::IsTrueCC1Pi(pEvent, config.global.useAbsPdg))
            continue;

        // Only use events passing the CC inclusive selection
        if (!pEvent->reco.passesCCInclusive())
            continue;

        // Get the reco and tru particles
        const auto weight = normalisation * AnalysisHelper::GetNominalEventWeight(pEvent);
        const auto recoParticles = pEvent->reco.particles;
        const auto truthParticles = pEvent->truth.particles;

        // Loop over the reco particles
        for (unsigned int index = 0; index < recoParticles.size(); ++index)
        {
            const auto &particle = recoParticles.at(index);

            // Insist the particle has a fitted track
            if (!AnalysisHelper::HasTrackFit(particle))
                continue;

            // Insist the particle is contained
            if (!AnalysisHelper::IsContained(particle))
                continue;

            // Get the best matched truth particle
            unsigned int truthParticleIndex = std::numeric_limits<unsigned int>::max();
            try
            {
                truthParticleIndex = AnalysisHelper::GetBestMatchedTruthParticleIndex(particle, truthParticles);
            }
            catch (const std::logic_error &)
            {
                // If there is no matched truth particle then skip it
                continue;
            }

            const auto truthParticle = truthParticles.at(truthParticleIndex);
            const auto truePdgCode = config.global.useAbsPdg ? std::abs(truthParticle.pdgCode()) : truthParticle.pdgCode();

            // Only use particles that match to pions
            if (truePdgCode != 211)
                continue;

            // Now breakdown the pions by type
            const auto isGolden = AnalysisHelper::IsGolden(truthParticle);
            const auto muonDecayTruth = (truthParticle.endState() == 1); // This value corresponds to pions that decay to a muon

            // Set thresholds to define a "visible" secondary
            const auto thresholdHitWeight = 5.f;
            const auto thresholdViews = 2u;
            const auto thresholdHitWeightTotal = 15.f;

            // Check if the secondary is visible in truth
            const auto nVisibleViewsTruth = ((truthParticle.endStateProductsHitWeightU() > thresholdHitWeight) ? 1u : 0u) +
                                            ((truthParticle.endStateProductsHitWeightV() > thresholdHitWeight) ? 1u : 0u) +
                                            ((truthParticle.endStateProductsHitWeightW() > thresholdHitWeight) ? 1u : 0u);

            const auto endStateProducesHitWeightTotal = truthParticle.endStateProductsHitWeightU() + truthParticle.endStateProductsHitWeightV() + truthParticle.endStateProductsHitWeightW();
            const auto isSecondaryVisibleTruth =  (nVisibleViewsTruth >= thresholdViews) && (endStateProducesHitWeightTotal >= thresholdHitWeightTotal);

            // Check if a secondary has been reconstructed
            const auto hasSecondaryReco = (particle.nDaughters() != 0);

            // Set the plot style - here we define a number for each possible classification that (when experessed as a 3-bit binary number)
            // the bits indicate if a given boolean is true or false. The first bit is muonDecayTruth, the second is isSecondaryVisibleTruth
            // and the third is hasSecondaryReco. Below are the possible combinations, and the corresponding styles.
            // ATTN here we order the styles so that the more common classifications get the lower order (primary, secondary) styles
            const std::vector<PlottingHelper::PlotStyle> plotStyles = {
                // Boolean formula                                                  -> Binary   Decimal
                // !muonDecayTruth && !isSecondaryVisibleTruth && !hasSecondaryReco -> 000    = 0
                PlottingHelper::Quinary,
                // !muonDecayTruth && !isSecondaryVisibleTruth &&  hasSecondaryReco -> 001    = 1
                PlottingHelper::Senary,
                // !muonDecayTruth &&  isSecondaryVisibleTruth && !hasSecondaryReco -> 010    = 2
                PlottingHelper::Tertiary,
                // !muonDecayTruth &&  isSecondaryVisibleTruth &&  hasSecondaryReco -> 011    = 3
                PlottingHelper::Quaternary,
                //  muonDecayTruth && !isSecondaryVisibleTruth && !hasSecondaryReco -> 100    = 4
                PlottingHelper::Octonary,
                //  muonDecayTruth && !isSecondaryVisibleTruth &&  hasSecondaryReco -> 101    = 5
                PlottingHelper::Septenary,
                //  muonDecayTruth &&  isSecondaryVisibleTruth && !hasSecondaryReco -> 110    = 6
                PlottingHelper::Primary,
                //  muonDecayTruth &&  isSecondaryVisibleTruth &&  hasSecondaryReco -> 111    = 7
                PlottingHelper::Secondary
            };

            const unsigned int classification = (muonDecayTruth ? 4u : 0u)
                                              + (isSecondaryVisibleTruth ? 2u : 0u)
                                              + (hasSecondaryReco ? 1u : 0u);

            const auto plotStyle = plotStyles.at(classification);

            // Get the BDT features (and skip particles for which the features aren't calculable)
            std::vector<float> features;
            if (!BDTHelper::GetBDTFeatures(particle, featureNames, features))
                continue;

            // Fill the feature plots
            for (unsigned int iFeature = 0; iFeature < featureNames.size(); ++iFeature)
            {
                if (isGolden)
                {
                    plotVectorGolden.at(iFeature).Fill(features.at(iFeature), plotStyle, weight);
                }
                else
                {
                    plotVectorNonGolden.at(iFeature).Fill(features.at(iFeature), plotStyle, weight);
                }
            }

            // Fill the BDT plots
            if (config.plotInputVariables.plotBDTResponses)
            {
                std::vector<float> goldenPionFeatures;
                const auto areAllFeaturesAvailableGoldenPion = BDTHelper::GetBDTFeatures(particle, goldenPionFeatureNames, goldenPionFeatures);

                std::vector<float> protonFeatures;
                const auto areAllFeaturesAvailableProton = BDTHelper::GetBDTFeatures(particle, protonFeatureNames, protonFeatures);

                std::vector<float> muonFeatures;
                const auto areAllFeaturesAvailableMuon = BDTHelper::GetBDTFeatures(particle, muonFeatureNames, muonFeatures);

                if (areAllFeaturesAvailableGoldenPion)
                {
                    const auto goldenPionBDTResponse = pGoldenPionBDT->GetResponse(goldenPionFeatures);

                    if (isGolden)
                    {
                        goldenPionBDTPlotGolden.Fill(goldenPionBDTResponse, plotStyle, weight);
                    }
                    else
                    {
                        goldenPionBDTPlotNonGolden.Fill(goldenPionBDTResponse, plotStyle, weight);
                    }

                }

                if (areAllFeaturesAvailableProton)
                {
                    const auto protonBDTResponse = pProtonBDT->GetResponse(protonFeatures);

                    if (isGolden)
                    {
                        protonBDTPlotGolden.Fill(protonBDTResponse, plotStyle, weight);
                    }
                    else
                    {
                        protonBDTPlotNonGolden.Fill(protonBDTResponse, plotStyle, weight);
                    }
                }

                if (areAllFeaturesAvailableMuon)
                {
                    const auto muonBDTResponse = pMuonBDT->GetResponse(muonFeatures);

                    if (isGolden)
                    {
                        muonBDTPlotGolden.Fill(muonBDTResponse, plotStyle, weight);
                    }
                    else
                    {
                        muonBDTPlotNonGolden.Fill(muonBDTResponse, plotStyle, weight);
                    }
                }
            }
        }
    }

    // Save the plots
    for (unsigned int iFeature = 0; iFeature < featureNames.size(); ++iFeature)
    {
        const auto &featureName = featureNames.at(iFeature);

        const bool useLogY = (featureName == "trackScore");
        const bool useLogX = (featureName == "wiggliness");

        plotVectorGolden.at(iFeature).SaveAsStacked("inputVariables_pionBreakdown_golden_stack_" + featureName, useLogX, false, useLogY);
        plotVectorNonGolden.at(iFeature).SaveAsStacked("inputVariables_pionBreakdown_nonGolden_stack_" + featureName, useLogX, false, useLogY);

        plotVectorGolden.at(iFeature).SaveAs("inputVariables_pionBreakdown_golden_" + featureName, useLogX, false, 100u, useLogY);
        plotVectorNonGolden.at(iFeature).SaveAs("inputVariables_pionBreakdown_nonGolden_" + featureName, useLogX, false, 100u, useLogY);
    }

    goldenPionBDTPlotGolden.SaveAsStacked("inputVariables_pionBreakdown_golden_stack_goldenPionBDTResponse");
    protonBDTPlotGolden.SaveAsStacked("inputVariables_pionBreakdown_golden_stack_protonBDTResponse");
    muonBDTPlotGolden.SaveAsStacked("inputVariables_pionBreakdown_golden_stack_muonBDTResponse");

    goldenPionBDTPlotNonGolden.SaveAsStacked("inputVariables_pionBreakdown_nonGolden_stack_goldenPionBDTResponse");
    protonBDTPlotNonGolden.SaveAsStacked("inputVariables_pionBreakdown_nonGolden_stack_protonBDTResponse");
    muonBDTPlotNonGolden.SaveAsStacked("inputVariables_pionBreakdown_nonGolden_stack_muonBDTResponse");

    goldenPionBDTPlotGolden.SaveAs("inputVariables_pionBreakdown_golden_goldenPionBDTResponse", false, false, 100u);
    protonBDTPlotGolden.SaveAs("inputVariables_pionBreakdown_golden_protonBDTResponse", false, false, 100u);
    muonBDTPlotGolden.SaveAs("inputVariables_pionBreakdown_golden_muonBDTResponse", false, false, 100u);

    goldenPionBDTPlotNonGolden.SaveAs("inputVariables_pionBreakdown_nonGolden_goldenPionBDTResponse", false, false, 100u);
    protonBDTPlotNonGolden.SaveAs("inputVariables_pionBreakdown_nonGolden_protonBDTResponse", false, false, 100u);
    muonBDTPlotNonGolden.SaveAs("inputVariables_pionBreakdown_nonGolden_muonBDTResponse", false, false, 100u);
}

} // ubcc1pi macros
