/**
 *  @file  ubcc1pi_standalone/Macros/ExtractXSecs.cxx
 *
 *  @brief The implementation file of the ExtractXSecs macro
 */

#include "ubcc1pi_standalone/Macros/Macros.h"
#include "ubcc1pi_standalone/Objects/FileReader.h"

#include "ubcc1pi_standalone/Helpers/AnalysisHelper.h"
#include "ubcc1pi_standalone/Helpers/SelectionHelper.h"
#include "ubcc1pi_standalone/Helpers/CrossSectionHelper.h"
#include "ubcc1pi_standalone/Helpers/NormalisationHelper.h"
#include "ubcc1pi_standalone/Helpers/FormattingHelper.h"

#include <TLine.h>
#include <TArrow.h>

using namespace ubcc1pi;

namespace ubcc1pi_macros
{

void ExtractXSecs(const Config &config)
{
    //
    // Setup the input files
    // 
    std::vector< std::tuple<AnalysisHelper::SampleType, std::string, float> > inputData;
    
    inputData.emplace_back(AnalysisHelper::Overlay, config.files.overlaysFileName, NormalisationHelper::GetOverlaysNormalisation(config)); 
    inputData.emplace_back(AnalysisHelper::Dirt,    config.files.dirtFileName,     NormalisationHelper::GetDirtNormalisation(config)); 
    inputData.emplace_back(AnalysisHelper::DataEXT, config.files.dataEXTFileName,  NormalisationHelper::GetDataEXTNormalisation(config)); 
    inputData.emplace_back(AnalysisHelper::DataBNB, config.files.dataBNBFileName,  1.f);
    
    // Get the selection
    auto selection = SelectionHelper::GetDefaultSelection();

    // Get the common input data and do any required unit conversions
    //
    // - Flux:            [10^-10 cm^-2 POT^-1]
    // - Exposure:        [10^20 POT]                 (stored in config as [POT])
    // - Target density:  [10^31 nucleons/cm^3]       (stored in config as [10^23 nucleons/cm^3])
    // - Fiducial volume: [cm^3]
    //
    // - Cross-section:   [Flux * Exposure * Target density * Fiducial volume]^-1 = [10^-41 cm^2 / nucleon]   // TODO check if this is right
    //
    CrossSectionHelper::CrossSection::InputData inputXSecData;
    inputXSecData.m_flux = CrossSectionHelper::GetFluxHistogram(config.flux);
    inputXSecData.m_exposurePOT = config.norms.dataBNBTor875WCut / (1e20);
    inputXSecData.m_nTargets = config.global.targetDensity * (1e-8) * AnalysisHelper::GetFiducialVolume();
    inputXSecData.m_fluxParameters = config.extractXSecs.fluxParams;

    // Print out the common input data
    std::cout << "flux = " << CrossSectionHelper::GetTotalFlux(inputXSecData.m_flux) << " 10^-10 cm^-2 POT^-1" << std::endl;
    std::cout << "exposure = " << inputXSecData.m_exposurePOT << " 10^20 POT" << std::endl;
    std::cout << "FV = " << AnalysisHelper::GetFiducialVolume() << " cm^3" << std::endl;
    std::cout << "nTargets = " << inputXSecData.m_nTargets << " 10^31 nucleons" << std::endl;

    // Build the systematic universe sizes map - this defines up-front the systematic parameters we want to apply and the number of
    // universes the we expect for each
    auto &systUniverseSizesMap = inputXSecData.m_systUniverseSizesMap;

    // Add the bootstrap universes (for the MC stat uncertainty)
    systUniverseSizesMap.emplace("bootstrap", config.extractXSecs.nBootstrapUniverses);

    // Add the independent systematic parameters
    for (const auto &[param, nUniverses] : config.extractXSecs.systematicParams)
        systUniverseSizesMap.emplace(param, nUniverses);
    
    // Add the mutually exclusive systematic parameters
    for (const auto &[param, nUniverses, exclusiveParams] : config.extractXSecs.mutuallyExclusiveSystematicParams)
        systUniverseSizesMap.emplace(param, nUniverses);

    //
    // Setup the cross-section objects
    //

    // Total cross-section - here we re-use the machinary of the differential cross-section by making a single bin measurement. The bin
    // edges we supply are just dummy values as we aren't measuring a kinematic paramter, here chosen to be [-1 -> 1]. We also request that
    // the cross-section isn't scaled by bin width.
    CrossSectionHelper::CrossSection xsec_total({-1.f, 1.f}, false, false, false, inputXSecData);
    
    // Muon momentum
    const auto &[muonMomentum_extendedBinEdges, muonMomentum_hasUnderflow, muonMomentum_hasOverflow] = CrossSectionHelper::GetExtendedBinEdges(config.global.muonMomentum);
    CrossSectionHelper::CrossSection xsec_muonMomentum(muonMomentum_extendedBinEdges, muonMomentum_hasUnderflow, muonMomentum_hasOverflow, true, inputXSecData);
    
    // Loop over the events
    for (const auto &[sampleType, fileName, normalisation] : inputData)
    {
        std::cout << "Reading input file: " << fileName << std::endl;
        const auto isOverlay = (sampleType == AnalysisHelper::Overlay);
        const auto isDataBNB = (sampleType == AnalysisHelper::DataBNB);

        // Open the input file for reading and enable the branches with systematic event weights (if required)
        FileReader reader(fileName);
        if (isOverlay) reader.EnableSystematicBranches();
        auto pEvent = reader.GetBoundEventAddress();

        const auto nEvents = reader.GetNumberOfEvents();
        for (unsigned int i = 0; i < nEvents; ++i)
        {
            AnalysisHelper::PrintLoadingBar(i, nEvents);

            //// BEGIN DEBUG
            // Reduce the rate of events for speed while debugging
            const unsigned int nEventsTarget = isOverlay ? 10000u : nEvents;
            const unsigned int batchSize = std::round(static_cast<float>(nEvents) / static_cast<float>(std::min(nEvents, nEventsTarget)));
            if (i % batchSize != 0)
                continue;

            // Work out the amount by which we need to up-weight the events we do use to make the scaling work
            const auto nEventsUsed = nEvents / batchSize; // Integer division
            const auto debugScaling = static_cast<float>(nEvents) / static_cast<float>(nEventsUsed);
            //// END DEBUG

            reader.LoadEvent(i);
            
            // Get the nominal event weight, scaled by the sample normalisation
            const auto weight = AnalysisHelper::GetNominalEventWeight(pEvent) * normalisation * debugScaling; //// DEBUG remove debugScaling
//            const auto weight = AnalysisHelper::GetNominalEventWeight(pEvent) * normalisation;
            
            // Get the systematic event weights non BNB data events
            CrossSectionHelper::SystematicWeightsMap systWeightsMap;
            if (!isDataBNB)
            {
                if (isOverlay)
                {
                    // For overlay events get the event weights from the input file
                    CrossSectionHelper::AddSystematicWeights(pEvent->truth, config.extractXSecs.systematicParams, systWeightsMap);
                    CrossSectionHelper::AddMutuallyExclusiveSystematicWeights(pEvent->truth, config.extractXSecs.mutuallyExclusiveSystematicParams, systWeightsMap);

                    // Generate random boostrap weights
                    CrossSectionHelper::AddBootstrapWeights(config.extractXSecs.nBootstrapUniverses, systWeightsMap);

                    // When we get the nominal event weight, this includes a parameter to make our MC look like the MCC9 tune of GENIE
                    // This parameter is stored under event.truth.genieTuneEventWeight. In general we want to apply this weight to every
                    // event. But, the systematic weights relating to GENIE are calculated using the MCC9 tune - so if we apply those
                    // weights AND genieTuneEventWeight we will "double count". To account for this, here we scale down the GENIE systematic
                    // weights by genieTuneEventWeight so they can behave like all other weights.
                    const auto genieTuneEventWeight = pEvent->truth.splineEventWeight.IsSet() ? pEvent->truth.genieTuneEventWeight() : 0.f;
                    if (pEvent->truth.splineEventWeight.IsSet() && std::abs(genieTuneEventWeight) > std::numeric_limits<float>::epsilon())
                    {
                        for (const auto &genieParam : config.extractXSecs.genieParams)
                        {
                            // Get the weights for this genie parameter 
                            auto genieParamIter = systWeightsMap.find(genieParam);
                            if (genieParamIter == systWeightsMap.end())
                                throw std::logic_error("ExtractXSecs - Supplied GENIE parameter: " + genieParam + " not found");

                            // Divide each weight down by the GENIE tune event weight
                            for (auto &weight : genieParamIter->second)
                                weight /= genieTuneEventWeight; 
                        }
                    }

                    // Add all overlay events in the fiducial volume to the cross-section objects so that we can use to total true neutrino
                    // energy spectrum in each universe to estimate the impact of the flux systematic variations on the total flux.
                    if (AnalysisHelper::IsFiducial(pEvent->truth.nuVertex()))
                    {
                        for (auto &xsec : std::vector<CrossSectionHelper::CrossSection>{
                            xsec_total,
                            xsec_muonMomentum
                        }) xsec.AddNeutrinoEvent(pEvent->truth.nuEnergy(), weight, systWeightsMap);
                    }
                }
                else
                {
                    // For non-overlay events just use unit weights
                    CrossSectionHelper::AddUnitWeights(config.extractXSecs.systematicParams, systWeightsMap);
                    CrossSectionHelper::AddUnitWeights(config.extractXSecs.mutuallyExclusiveSystematicParams, systWeightsMap);

                    // Generate random boostrap weights
                    CrossSectionHelper::AddBootstrapWeights(config.extractXSecs.nBootstrapUniverses, systWeightsMap);
                }
            }


            // Run the event selection and store which cuts are passed
            std::vector<std::string> cutsPassed;
            std::vector<int> assignedPdgCodes;
            auto passedGoldenSelection = selection.Execute(pEvent, cutsPassed, assignedPdgCodes);
            auto passedGenericSelection = (std::find(cutsPassed.begin(), cutsPassed.end(), config.global.lastCutGeneric) != cutsPassed.end());
            
            // Determine if the phase-space restrictions are met in reco
            auto recoData = AnalysisHelper::GetDummyAnalysisData();
            bool passesPhaseSpaceReco = true;
            if (passedGenericSelection)
            {
                // Get the reconstructed kinematic parameters
                recoData = AnalysisHelper::GetRecoAnalysisData(pEvent->reco, assignedPdgCodes, passedGoldenSelection);

                if (recoData.muonMomentum < config.global.muonMomentum.min || recoData.muonMomentum > config.global.muonMomentum.max)
                    passesPhaseSpaceReco = false;
            }

            // Apply the phase-space restrictions in reco as an additional event selection cut
            const bool isSelectedGeneric = passedGenericSelection && passesPhaseSpaceReco;
            const bool isSelectedGolden = passedGoldenSelection && passesPhaseSpaceReco;


            // Determine if the event is truly a signal event
            const auto isTrueCC1Pi = isOverlay && AnalysisHelper::IsTrueCC1Pi(pEvent, config.global.useAbsPdg);

            // Determine of the phase-space restrictions are met in truth
            auto truthData = AnalysisHelper::GetDummyAnalysisData();
            bool passesPhaseSpaceTruth = true;
            if (isTrueCC1Pi)
            {
                // Get the true kinematic parameters
                truthData = AnalysisHelper::GetTruthAnalysisData(pEvent->truth, config.global.useAbsPdg, config.global.protonMomentumThreshold);

                if (truthData.muonMomentum < config.global.muonMomentum.min || truthData.muonMomentum > config.global.muonMomentum.max)
                    passesPhaseSpaceTruth = false;
            }

            // A signal event is a true CC1Pi event that passes the phase space restrictions
            const bool isSignal = isTrueCC1Pi && passesPhaseSpaceTruth;

            
            // We only need to count signal events, or any event that passes the generic selection - everything else we can safely skip
            const bool shouldCount = isSignal || isSelectedGeneric;
            if (!shouldCount)
                continue;

            // Count selected BNB data events
            if (isDataBNB)
            {
                if (isSelectedGeneric)
                {
                    xsec_total.AddSelectedBNBDataEvent(0.f); // ATTN here a dummy value is used of 0.f
                    xsec_muonMomentum.AddSelectedBNBDataEvent(recoData.muonMomentum);
                }
            }
            else
            {
                if (isSignal)
                {
                    // Count selected signal events

                    xsec_total.AddSignalEvent(0.f, 0.f, isSelectedGeneric, weight, systWeightsMap); // ATTN here a dummy value is used of 0.f

                    xsec_muonMomentum.AddSignalEvent(truthData.muonMomentum, recoData.muonMomentum, isSelectedGeneric, weight, systWeightsMap);
                }
                else
                {
                    // Count selected background events
                    
                    if (isSelectedGeneric)
                    {
                        xsec_total.AddSelectedBackgroundEvent(0.f, weight, systWeightsMap); // ATTN here a dummy value is used of 0.f
                        
                        xsec_muonMomentum.AddSelectedBackgroundEvent(recoData.muonMomentum, weight, systWeightsMap);
                    }
                }
            }
        }
    }

    // Get the flux variations for each parameter and save them
    auto pCanvas = PlottingHelper::GetCanvas();
    for (const auto &[paramName, pHist] : xsec_total.GetFluxVariations())
    {
        pHist->Draw("colz");
        inputXSecData.m_flux->Draw("hist same");
        PlottingHelper::SaveCanvas(pCanvas, "xsec_fluxVariations_" + paramName);
    }

    /*
    // Define a function to make the scatter plots
    const auto MakeScatterPlots = [](const auto &xsec, const auto &covariances, const auto &scatterPlots, const std::string &namePrefix)
    {
        auto pSquareCanvas = PlottingHelper::GetCanvas(960, 960);

        for (const auto &[param, covarianceBiasPair] : covariances)
        {
            for (const auto &[iBin, jBin, graph] : scatterPlots.at(param))
            {
                // Get the eigenvectors of the covariance matrix
                const auto eigenVectors = CrossSectionHelper::GetUncertaintyEigenVectors(covarianceBiasPair, iBin, jBin);

                // Get the bias vector
                const auto iBias = covarianceBiasPair.second->GetBinContent(iBin);
                const auto jBias = covarianceBiasPair.second->GetBinContent(jBin);

                // Get the cross-section in the nominal universe in these bins
                const auto iValue = xsec->GetBinContent(iBin);
                const auto jValue = xsec->GetBinContent(jBin);

                // Draw the scatter graph
                graph->Draw("ap");

                // Add the bias line
                std::vector< std::shared_ptr<TArrow> > lines;
                lines.emplace_back(std::make_shared<TArrow>(
                            iValue, jValue,                  // X, Y start of line
                            iValue + iBias, jValue + jBias,  // X, Y end of line
                            0.01, "->"                       // Arrow style
                ));

                auto &pBiasLine = lines.back();
                pBiasLine->SetLineColor(kBlue);
                pBiasLine->Draw();

                // Add the eigenvector lines
                for (const auto &vector : {eigenVectors.first, eigenVectors.second})
                {
                    lines.emplace_back(std::make_shared<TArrow>(
                                iValue + iBias - vector.X(), jValue + jBias - vector.Y(),  // X, Y start of line
                                iValue + iBias + vector.X(), jValue + jBias + vector.Y(),  // X, Y end of line
                                0.01, "<>"                                                 // Arrow style
                    ));

                    auto &pLine = lines.back();
                    pLine->SetLineColor(kRed);
                    pLine->Draw();
                }

                PlottingHelper::SaveCanvas(pSquareCanvas, "xsec_scatter_" + namePrefix + "_" + param + "_" + std::to_string(iBin) + "-" + std::to_string(jBin));
            }
        }
    };

    // Define a function to save the cross-section plots
    const auto MakeXSecPlot = [] (const auto &xsec, const auto &plot, const auto &covariances, const std::string &namePrefix)
    {
        auto pCanv = PlottingHelper::GetCanvas();

        // Loop over the normal bins and get the bin edges
        const unsigned int nBinsTotal = plot->GetXaxis()->GetNbins();
        std::vector<float> binEdges;
        for (unsigned int iBin = 1; iBin <= nBinsTotal; ++iBin)
        {
            if (xsec.IsUnderOverflowBin(iBin))
                continue;

            binEdges.push_back(plot->GetBinLowEdge(iBin));
        }
        
        const auto hasOverflow = xsec.IsUnderOverflowBin(nBinsTotal);
        const auto lastBinIndex = hasOverflow ? (nBinsTotal-1) : nBinsTotal;
        const auto upperEdge = plot->GetBinLowEdge(lastBinIndex) + plot->GetBinWidth(lastBinIndex);
        binEdges.push_back(upperEdge);

        // Make new plots skipping the under-overflow bins
        const auto nNewBins = binEdges.size() - 1;
        auto pXSecHist = std::make_shared<TH1F>(("xsec_plot_" + namePrefix).c_str(), "", nNewBins, binEdges.data());
        std::map< std::string, std::shared_ptr<TH1F> > biasMap, fracBiasMap;
        std::map< std::string, std::shared_ptr<TH2F> > covarianceMap, fracCovarianceMap;

        for (const auto &[paramName, covarianceBiasPair] : covariances)
        {
            biasMap.emplace(paramName, std::make_shared<TH1F>(("xsec_bias_" + namePrefix).c_str(), "", nNewBins, 0, nNewBins));
            fracBiasMap.emplace(paramName, std::make_shared<TH1F>(("xsec_fracBias_" + namePrefix).c_str(), "", nNewBins, 0, nNewBins));
            covarianceMap.emplace(paramName, std::make_shared<TH2F>(("xsec_cov_" + namePrefix).c_str(), "", nNewBins, 0, nNewBins, nNewBins, 0, nNewBins));
            fracCovarianceMap.emplace(paramName, std::make_shared<TH2F>(("xsec_fracCov_" + namePrefix).c_str(), "", nNewBins, 0, nNewBins, nNewBins, 0, nNewBins));
        }

        unsigned int iBinNew = 1;
        for (unsigned int iBin = 1; iBin <= nBinsTotal; ++iBin)
        {
            if (xsec.IsUnderOverflowBin(iBin))
                continue;

            pXSecHist->SetBinContent(iBinNew, plot->GetBinContent(iBin));
            
            for (const auto &[paramName, covarianceBiasPair] : covariances)
            {
                biasMap.at(paramName)->SetBinContent(iBinNew, covarianceBiasPair.second->GetBinContent(iBin));
                fracBiasMap.at(paramName)->SetBinContent(iBinNew, covarianceBiasPair.second->GetBinContent(iBin) / plot->GetBinContent(iBin));
                
                unsigned int jBinNew = 1;
                for (unsigned int jBin = 1; jBin <= nBinsTotal; ++jBin)
                {
                    if (xsec.IsUnderOverflowBin(jBin))
                        continue;
                
                    covarianceMap.at(paramName)->SetBinContent(iBinNew, jBinNew, covarianceBiasPair.first->GetBinContent(iBin, jBin));
                    fracCovarianceMap.at(paramName)->SetBinContent(iBinNew, jBinNew, covarianceBiasPair.first->GetBinContent(iBin, jBin) / (plot->GetBinContent(iBin) * plot->GetBinContent(jBin)));

                    jBinNew++;
                }
            }

            iBinNew++;
        }
   
        // Save the plots
        pXSecHist->Draw("hist");
        PlottingHelper::SaveCanvas(pCanv, "xsec_" + namePrefix);

        for (const auto &[paramName, pHist] : covarianceMap)
        {
            pHist->Draw("colz");
            PlottingHelper::SaveCanvas(pCanv, "xsec_covariance_" + namePrefix + "_" + paramName);
        }
        
        for (const auto &[paramName, pHist] : fracCovarianceMap)
        {
            pHist->Draw("colz");
            PlottingHelper::SaveCanvas(pCanv, "xsec_fracCovariance_" + namePrefix + "_" + paramName);
        }
        
        for (const auto &[paramName, pHist] : biasMap)
        {
            pHist->Draw("colz");
            PlottingHelper::SaveCanvas(pCanv, "xsec_bias_" + namePrefix + "_" + paramName);
        }

        for (const auto &[paramName, pHist] : fracBiasMap)
        {
            pHist->Draw("colz");
            PlottingHelper::SaveCanvas(pCanv, "xsec_fracBias_" + namePrefix + "_" + paramName);
        }
    };
*/

    // Calculate the cross-section and print it out!

    for (auto &[namePrefix, xsec] : std::vector<std::pair<std::string, CrossSectionHelper::CrossSection> >({
        {"muonMomentum", xsec_muonMomentum},
        {"total", xsec_total},
    }))
    {
        FormattingHelper::PrintLine();
        std::cout << namePrefix << std::endl;
        FormattingHelper::PrintLine();

        // Get the cross-section itself
        std::cout << "Cross-section" << std::endl;
        const auto &xsec_plot = xsec.GetCrossSection();
        FormattingHelper::SaveHistAsTable(xsec_plot, "xsec_" + namePrefix + ".md");

        // Get the covariance matrices and bias vectors for each systematic parameter
        std::cout << "Getting covariance & biases" << std::endl;
        const auto &xsec_covariances = xsec.GetCrossSectionCovarianceMatricies();
        for (const auto &[paramName, covarianceBiasPair] : xsec_covariances)
        {
            const auto &covarianceMatrix = covarianceBiasPair.first;
            const auto &biasVector = covarianceBiasPair.second;

            std::cout << "Bias vector" << std::endl;
            FormattingHelper::SaveHistAsTable(biasVector, "xsec_" + namePrefix + "_" + paramName + "_bias.md");
            PlottingHelper::SaveBiasVector(biasVector, xsec_plot, xsec.HasUnderflowBin(), xsec.HasOverflowBin(), "xsec_" + namePrefix + "_" + paramName);
            
            std::cout << "Covariance matrix" << std::endl;
            FormattingHelper::SaveHistAsTable(covarianceMatrix, "xsec_" + namePrefix + "_"+ paramName + "_covariance.md");
            PlottingHelper::SaveCovarianceMatrix(covarianceMatrix, xsec_plot, xsec.HasUnderflowBin(), xsec.HasOverflowBin(), "xsec_" + namePrefix + "_" + paramName);
        }

        // Get the total covariance matrix, accounting for the biases
        std::cout << "Getting total covariance matrix" << std::endl;
        const auto &xsec_totalCovariance = CrossSectionHelper::GetTotalCovarianceMatrix(xsec_covariances);
        FormattingHelper::SaveHistAsTable(xsec_totalCovariance, "xsec_" + namePrefix + "_allParams_covariance.md");
        PlottingHelper::SaveCovarianceMatrix(xsec_totalCovariance, xsec_plot, xsec.HasUnderflowBin(), xsec.HasOverflowBin(), "xsec_" + namePrefix + "_allParams");

        // For performance clear the cached cross-sections from memory
        xsec.ClearCache();
    }
}

} // namespace ubcc1pi_macros
