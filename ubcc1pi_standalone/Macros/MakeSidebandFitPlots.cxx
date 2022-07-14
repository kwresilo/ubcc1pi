/**
 *  @file  ubcc1pi_standalone/Macros/MakeSidebandFitPlots.cxx
 *
 *  @brief The implementation file of the MakeSidebandFitPlots macro
 */

#include "ubcc1pi_standalone/Macros/Macros.h"
#include "ubcc1pi_standalone/Objects/FileReader.h"

#include "ubcc1pi_standalone/Helpers/PlottingHelper.h"
#include "ubcc1pi_standalone/Helpers/FormattingHelper.h"
#include "ubcc1pi_standalone/Helpers/CrossSectionHelper.h"

#include "ubsmear.h"

#include <TStyle.h>

// Boost libraries
#include "binary_iarchive.hpp"
// #include "binary_oarchive.hpp"
#include "binary_object.hpp"
#include "map.hpp"
#include "vector.hpp"

using namespace ubcc1pi;

namespace ubcc1pi_macros
{

void MakeSidebandFitPlots(const Config &config)
{

    CrossSectionHelper::CrossSection::SystParams systParams;
    systParams.nBootstrapUniverses = config.extractXSecs.nBootstrapUniverses;
    systParams.fluxDimensions = config.extractXSecs.fluxDimensions;
    systParams.xsecDimensions = config.extractXSecs.xsecDimensions;
    systParams.reintDimensions = config.extractXSecs.reintDimensions;
    systParams.detVarDimensions = config.extractXSecs.detVarDimensions;
    systParams.potFracUncertainty = config.extractXSecs.potFracUncertainty;

    // -------------------------------------------------------------------------------------------------------------------------------------
    // Get the binning of each cross-section
    // -------------------------------------------------------------------------------------------------------------------------------------
    std::map<std::string, ubsmear::UBXSecMeta> metadataMap;

    for (const auto &[name, binning, scaleByBinWidth] : std::vector< std::tuple<std::string, Config::Global::Binning, bool> > {

        // The names of the cross-section kinematic parameters, and their binning information.
        // The third (boolean) parameter indicates if the cross-section bins should be scaled by their width
        { "muonCosTheta",  config.global.muonCosTheta,  true  },
        { "muonPhi",       config.global.muonPhi,       true  },
        { "muonMomentum",  config.global.muonMomentum,  true  },

        { "pionCosTheta",  config.global.pionCosTheta,  true  },
        { "pionPhi",       config.global.pionPhi,       true  },
        { "pionMomentum",  config.global.pionMomentum,  true  },

        { "muonPionAngle", config.global.muonPionAngle, true  },
        { "nProtons",      config.global.nProtons,      false }

    })
    {
        std::cout << "Getting bins for: " << name << std::endl;

        // Get the bin edges from the input configuration
        const auto &[extendedBinEdges, hasUnderflow, hasOverflow] = CrossSectionHelper::GetExtendedBinEdges(binning.min, binning.max, binning.binEdges);
        metadataMap.emplace(name, ubsmear::UBXSecMeta(extendedBinEdges, hasUnderflow, hasOverflow, scaleByBinWidth));
    }

    // Add the dummy metadata for the total cross-section (see ExtractXSecs for more details)
    metadataMap.emplace("total", ubsmear::UBXSecMeta({-1.f, 1.f}, false, false, false));


    std::cout<<"\n\n\n\nUsing NuWro parameters!!!!!!!!!!!!!!"<<std::endl;
    // // -------------------------------------------------------------------------------------------------------------------------------------
    // // Get the sideband weights
    // // -------------------------------------------------------------------------------------------------------------------------------------
    // // Loop over all cross-section objects
    typedef std::pair<std::vector<Double_t>,std::vector<Double_t>> paramAndErrorPair; // Todo: Improve code!
    //Parameters: selectionName name
    std::map<std::string,std::map<std::string, std::vector<float>>> cc0piCovarianceMap; 
    std::map<std::string, std::map<std::string, paramAndErrorPair>> cc0piNominalConstraintMap;
    //Parameters: selectionName name paramName (i.e. golden muonMomentum hadronProduction)
    std::map<std::string, std::map<std::string, std::map<std::string, std::vector<paramAndErrorPair>>>> cc0piUniverseConstraintMap;

    std::cout<<"\n\nDone with CC0pi constraint."<<std::endl;

	std::ifstream ifs1("cc0piCovarianceMapNuWro.bin", std::ios::binary);
    std::ifstream ifs2("cc0piNominalConstraintMapNuWro.bin", std::ios::binary);
    std::ifstream ifs3("cc0piUniverseConstraintMapNuWro.bin", std::ios::binary);
	
    boost::archive::binary_iarchive iarch1(ifs1);
    boost::archive::binary_iarchive iarch2(ifs2);
    boost::archive::binary_iarchive iarch3(ifs3);

    iarch1 >> cc0piCovarianceMap;
    iarch2 >> cc0piNominalConstraintMap;
    iarch3 >> cc0piUniverseConstraintMap;

	ifs1.close();
    ifs2.close();
    ifs3.close();




    // -----------------------------------------------------------------------------------------------------------------------------
    // Make nominal plots
    // -----------------------------------------------------------------------------------------------------------------------------
    std::cout<<"Debug point 0"<<std::endl;
    for (const auto &[selectionName, enabledMap] : config.extractXSecs.crossSectionIsEnabled)
    {
        if(selectionName=="golden") continue;

        for (const auto &[xsecName, isEnabled] : enabledMap)
        {
            // Skip cross-sections that aren't enabled
            if (!isEnabled)
                continue;

            std::cout<<"Debug point 1"<<std::endl;
            // Define a prefix for the names of the plots
            const std::string prefix = "nuWroSidebandFitPlots_" + selectionName + "_" + xsecName;
            const auto &metadata = metadataMap.at(xsecName);

            const auto selectedEventsSignalTruth = ubsmear::UBFileHelper::ReadMatrix("NuWroSidebandFit_" + selectionName + "_" + xsecName + "_selectedEventsSignalTruth.txt");
            const auto selectedEventsData = ubsmear::UBFileHelper::ReadMatrix("NuWroSidebandFit_" + selectionName + "_" + xsecName + "_selectedEventsData.txt");
            const auto selectedEventsBackgroundReco = ubsmear::UBFileHelper::ReadMatrix("NuWroSidebandFit_" + selectionName + "_" + xsecName + "_selectedEventsBackgroundReco.txt");
            const auto smearingMatrix = ubsmear::UBFileHelper::ReadMatrix("NuWroSidebandFit_" + selectionName + "_" + xsecName + "_smearingMatrix.txt");

            const auto signalData = selectedEventsData - selectedEventsBackgroundReco;
            const auto selectedEventsSignalReco = smearingMatrix * selectedEventsSignalTruth;

            auto cc0piNominalConstraintParam = cc0piNominalConstraintMap.at(selectionName).at(xsecName).first;
            auto cc0piNominalConstraintParamError = cc0piNominalConstraintMap.at(selectionName).at(xsecName).second;
            const vector<float> cc0piNominalConstraintParamFloat(cc0piNominalConstraintParam.begin(), cc0piNominalConstraintParam.end());
            const ubsmear::UBMatrix sidebandParamVectorTruth(cc0piNominalConstraintParamFloat, cc0piNominalConstraintParam.size(), 1);
            const auto scaledSelectedEventsSignalTruth = ElementWiseOperation(selectedEventsSignalTruth, sidebandParamVectorTruth, [](const auto &l, const auto& r) { return l * r; });
            const auto scaledSelectedEventsSignalReco = smearingMatrix * scaledSelectedEventsSignalTruth;

            // -----------------------------------------------------------------------------------------------------------------------------
            // Make the comparison plot between data and smeared prediction
            // -----------------------------------------------------------------------------------------------------------------------------
            // Now we have the data and the prediction on the same footing we can compare them on a plot!
            // Get the bin edges (excluding any underflow/overflow bins)
            const auto extendedBinEdges = metadata.GetBinEdges();
            std::vector<float> binEdges;
            for (unsigned int iBin = 0; iBin < metadata.GetNBins(); ++iBin)
            {
                // // Skip underflow/overflow bins
                // if (metadata.IsUnderOverflowBin(iBin))
                //     continue;

                // If this is the first bin then add the lower edge
                if (binEdges.empty())
                {
                    if (metadata.IsScaledByBinWidth())
                    {
                        binEdges.push_back(std::min(5.f, extendedBinEdges.at(iBin)));
                    }
                    else
                    {
                        // If we don't scale by bin width, then just use zero as the first bin edge
                        binEdges.push_back(0.f);
                    }
                }

                // Add the upper bin edge
                if (metadata.IsScaledByBinWidth())
                {
                    binEdges.push_back(std::min(5.f, extendedBinEdges.at(iBin + 1)));
                }
                else
                {
                    // If we don't scale by bin width, then just use unit width bins
                    binEdges.push_back(binEdges.back() + 1.f);
                }
            }



            std::cout<<"Debug point 2"<<std::endl;
            // -----------------------------------------------------------------------------------------------------------------------------
            // Data vs Prediction
            // -----------------------------------------------------------------------------------------------------------------------------
            // Setup the data histogram and the prediction histogram
            auto pDataHist = std::make_shared<TH1F>((prefix + "_raw_data").c_str(), "", binEdges.size() - 1, binEdges.data());
            // auto pDataStatOnlyHist = std::make_shared<TH1F>((prefix + "_dataStatOnly").c_str(), "", binEdges.size() - 1, binEdges.data());
            auto pPredictionHist = std::make_shared<TH1F>((prefix + "_raw_prediction").c_str(), "", binEdges.size() - 1, binEdges.data());

            // Fill the bins
            // ATTN here we only show the diagonals of the error matrices
            float minY = +std::numeric_limits<float>::max();
            float maxY = -std::numeric_limits<float>::max();
            // const auto dataStatErrorMatrix = errorMatrixMap.at("data").at("stat");
            for (unsigned int iBin = 1; iBin <= binEdges.size() - 1; ++iBin)
            {
                // Set the values for the data histogram
                const auto dataValue = signalData.At(iBin - 1, 0);
                // const auto dataError = std::pow(totalDataErrorMatrix.At(iBin - 1, iBin - 1), 0.5f);
                // const auto dataStatOnlyError = std::pow(dataStatErrorMatrix.At(iBin - 1, iBin - 1), 0.5f);

                pDataHist->SetBinContent(iBin, dataValue);
                pDataHist->SetBinError(iBin, 0);

                // pDataStatOnlyHist->SetBinContent(iBin, dataValue);
                // pDataStatOnlyHist->SetBinError(iBin, dataStatOnlyError);

                // Set the values of the prediction
                const auto predictionValue = selectedEventsSignalReco.At(iBin - 1, 0);
                // const auto predictionError = std::pow(smearedPredictionErrorMatrix.At(iBin - 1, iBin - 1), 0.5f);

                pPredictionHist->SetBinContent(iBin, predictionValue);
                pPredictionHist->SetBinError(iBin, 0);

                // For the proton multiplicity plot, use explicit bin labels
                if (xsecName == "nProtons")
                {
                    for (auto &pHist : {pDataHist, pPredictionHist})
                    {
                        pHist->GetXaxis()->SetBinLabel(1, "0");
                        pHist->GetXaxis()->SetBinLabel(2, "1");
                        pHist->GetXaxis()->SetBinLabel(3, ">1");
                    }
                }

                // Get the limiting values
                minY = std::min(minY, dataValue);
                minY = std::min(minY, predictionValue);

                maxY = std::max(maxY, dataValue);
                maxY = std::max(maxY, predictionValue);
            }

            // Set the y-range
            auto padding = (maxY - minY) * 0.05;
            maxY += padding;
            minY -= padding;
            minY = std::max(minY, 0.f);
            minY = 0.f; // Remove this line to get a dynamic lower y-range
            pDataHist->GetYaxis()->SetRangeUser(minY, maxY);
            // pDataStatOnlyHist->GetYaxis()->SetRangeUser(minY, maxY);
            pPredictionHist->GetYaxis()->SetRangeUser(minY, maxY);

            // Set the colours of the histograms
            PlottingHelper::SetLineStyle(pDataHist, PlottingHelper::Primary);
            // PlottingHelper::SetLineStyle(pDataStatOnlyHist, PlottingHelper::Primary);
            PlottingHelper::SetLineStyle(pPredictionHist, PlottingHelper::Secondary);

            // Make the plot!
            auto pCanvas1 = PlottingHelper::GetCanvas();
            gStyle->SetEndErrorSize(4);

            // Draw the smeared prediction
            pPredictionHist->Draw("hist");

            // Draw the prediction uncertainties as a semi-transparent band
            auto pHistClone = static_cast<TH1F *>(pPredictionHist->Clone());
            auto col = pHistClone->GetLineColor();
            pHistClone->SetFillStyle(1001);
            pHistClone->SetLineColorAlpha(col, 0.f);
            pHistClone->SetFillColorAlpha(col, 0.3f);
            pHistClone->Draw("e2 same");

            // Draw the data as points with error bars
            // pDataStatOnlyHist->Draw("e1 same");
            pDataHist->Draw("e1 same");

            PlottingHelper::SaveCanvas(pCanvas1, prefix + "_raw_data-vs-smearedPrediction");



            std::cout<<"Debug point 3"<<std::endl;
            // -----------------------------------------------------------------------------------------------------------------------------
            // Data vs Scaled Prediction
            // -----------------------------------------------------------------------------------------------------------------------------
            // Setup the data histogram and the prediction histogram
            auto pDataHist2 = std::make_shared<TH1F>((prefix + "_raw_data2").c_str(), "", binEdges.size() - 1, binEdges.data());
            // auto pDataStatOnlyHist = std::make_shared<TH1F>((prefix + "_dataStatOnly").c_str(), "", binEdges.size() - 1, binEdges.data());
            auto pPredictionHist2 = std::make_shared<TH1F>((prefix + "_scaled_raw_prediction").c_str(), "", binEdges.size() - 1, binEdges.data());

            // Fill the bins
            // ATTN here we only show the diagonals of the error matrices
            minY = +std::numeric_limits<float>::max();
            maxY = -std::numeric_limits<float>::max();
            // const auto dataStatErrorMatrix = errorMatrixMap.at("data").at("stat");
            for (unsigned int iBin = 1; iBin <= binEdges.size() - 1; ++iBin)
            {
                // Set the values for the data histogram
                const auto dataValue = signalData.At(iBin - 1, 0);
                // const auto dataError = std::pow(totalDataErrorMatrix.At(iBin - 1, iBin - 1), 0.5f);
                // const auto dataStatOnlyError = std::pow(dataStatErrorMatrix.At(iBin - 1, iBin - 1), 0.5f);

                pDataHist2->SetBinContent(iBin, dataValue);
                pDataHist2->SetBinError(iBin, 0);

                // pDataStatOnlyHist->SetBinContent(iBin, dataValue);
                // pDataStatOnlyHist->SetBinError(iBin, dataStatOnlyError);

                // Set the values of the prediction
                const auto predictionValue = scaledSelectedEventsSignalReco.At(iBin - 1, 0);
                // const auto predictionError = std::pow(smearedPredictionErrorMatrix.At(iBin - 1, iBin - 1), 0.5f);

                pPredictionHist2->SetBinContent(iBin, predictionValue);
                pPredictionHist2->SetBinError(iBin, 0);

                // For the proton multiplicity plot, use explicit bin labels
                if (xsecName == "nProtons")
                {
                    for (auto &pHist2 : {pDataHist2, pPredictionHist2})
                    {
                        pHist2->GetXaxis()->SetBinLabel(1, "0");
                        pHist2->GetXaxis()->SetBinLabel(2, "1");
                        pHist2->GetXaxis()->SetBinLabel(3, ">1");
                    }
                }

                // Get the limiting values
                minY = std::min(minY, dataValue);
                minY = std::min(minY, predictionValue);

                maxY = std::max(maxY, dataValue);
                maxY = std::max(maxY, predictionValue);
            }

            // Set the y-range
            padding = (maxY - minY) * 0.05;
            maxY += padding;
            minY -= padding;
            minY = std::max(minY, 0.f);
            minY = 0.f; // Remove this line to get a dynamic lower y-range
            pDataHist2->GetYaxis()->SetRangeUser(minY, maxY);
            // pDataStatOnlyHist->GetYaxis()->SetRangeUser(minY, maxY);
            pPredictionHist2->GetYaxis()->SetRangeUser(minY, maxY);

            // Set the colours of the histograms
            PlottingHelper::SetLineStyle(pDataHist2, PlottingHelper::Primary);
            // PlottingHelper::SetLineStyle(pDataStatOnlyHist, PlottingHelper::Primary);
            PlottingHelper::SetLineStyle(pPredictionHist2, PlottingHelper::Secondary);

            // Make the plot!
            auto pCanvas2 = PlottingHelper::GetCanvas();
            gStyle->SetEndErrorSize(4);

            // Draw the smeared prediction
            pPredictionHist2->Draw("hist");

            // Draw the prediction uncertainties as a semi-transparent band
            auto pHistClone2 = static_cast<TH1F *>(pPredictionHist2->Clone());
            col = pHistClone2->GetLineColor();
            pHistClone2->SetFillStyle(1001);
            pHistClone2->SetLineColorAlpha(col, 0.f);
            pHistClone2->SetFillColorAlpha(col, 0.3f);
            pHistClone2->Draw("e2 same");

            // Draw the data as points with error bars
            // pDataStatOnlyHist->Draw("e1 same");
            pDataHist2->Draw("e1 same");

            PlottingHelper::SaveCanvas(pCanvas2, prefix + "_scaled_raw_data-vs-smearedPrediction");

            std::cout<<"Debug point 4"<<std::endl;
            // -----------------------------------------------------------------------------------------------------------------------------
            // Scaling factors
            // -----------------------------------------------------------------------------------------------------------------------------
            // Setup the data histogram and the prediction histogram
            auto pFactorHist = std::make_shared<TH1F>((prefix + "_factors").c_str(), "", binEdges.size() - 1, binEdges.data());
            // auto pDataStatOnlyHist = std::make_shared<TH1F>((prefix + "_dataStatOnly").c_str(), "", binEdges.size() - 1, binEdges.data());
            // auto pPredictionHist = std::make_shared<TH1F>((prefix + "_prediction").c_str(), "", binEdges.size() - 1, binEdges.data());

            // Fill the bins
            // ATTN here we only show the diagonals of the error matrices
            minY = +std::numeric_limits<float>::max();
            maxY = -std::numeric_limits<float>::max();
            // const auto dataStatErrorMatrix = errorMatrixMap.at("data").at("stat");
            for (unsigned int iBin = 1; iBin <= binEdges.size() - 1; ++iBin)
            {
                // Set the values for the data histogram
                const auto factorValue = (float) cc0piNominalConstraintParam.at(iBin - 1);
                const auto factorError = (float) cc0piNominalConstraintParamError.at(iBin - 1);//std::pow( ... , 0.5f);
                // const auto dataStatOnlyError = std::pow(dataStatErrorMatrix.At(iBin - 1, iBin - 1), 0.5f);

                pFactorHist->SetBinContent(iBin, factorValue);
                pFactorHist->SetBinError(iBin, factorError);

                // pDataStatOnlyHist->SetBinContent(iBin, dataValue);
                // pDataStatOnlyHist->SetBinError(iBin, dataStatOnlyError);

                // Set the values of the prediction
                // const auto predictionValue = smearedPrediction.At(iBin - 1, 0);
                // const auto predictionError = std::pow(smearedPredictionErrorMatrix.At(iBin - 1, iBin - 1), 0.5f);

                // pPredictionHist->SetBinContent(iBin, predictionValue);
                // pPredictionHist->SetBinError(iBin, predictionError);

                // For the proton multiplicity plot, use explicit bin labels
                if (xsecName == "nProtons")
                {
                    for (auto &pHist : {pFactorHist})//, pDataStatOnlyHist, pPredictionHist})
                    {
                        pHist->GetXaxis()->SetBinLabel(1, "0");
                        pHist->GetXaxis()->SetBinLabel(2, "1");
                        pHist->GetXaxis()->SetBinLabel(3, ">1");
                    }
                }

                // Get the limiting values
                minY = std::min(minY, factorValue - factorError);
                // minY = std::min(minY, predictionValue - predictionError);

                maxY = std::max(maxY, factorValue + factorError);
                // maxY = std::max(maxY, predictionValue + predictionError);
            }

            // Set the y-range
            padding = (maxY - minY) * 0.05;
            maxY += padding;
            minY -= padding;
            minY = std::max(minY, 0.f);
            minY = 0.f; // Remove this line to get a dynamic lower y-range
            pFactorHist->GetYaxis()->SetRangeUser(minY, maxY);
            // pDataStatOnlyHist->GetYaxis()->SetRangeUser(minY, maxY);
            // pPredictionHist->GetYaxis()->SetRangeUser(minY, maxY);

            // Set the colours of the histograms
            PlottingHelper::SetLineStyle(pFactorHist, PlottingHelper::Primary);
            // PlottingHelper::SetLineStyle(pDataStatOnlyHist, PlottingHelper::Primary);
            // PlottingHelper::SetLineStyle(pPredictionHist, PlottingHelper::Secondary);

            // Make the plot!
            auto pCanvas3 = PlottingHelper::GetCanvas();
            gStyle->SetEndErrorSize(4);

            // // Draw the smeared prediction
            // pPredictionHist->Draw("hist");

            // // Draw the prediction uncertainties as a semi-transparent band
            // auto pHistClone = static_cast<TH1F *>(pPredictionHist->Clone());
            // const auto col = pHistClone->GetLineColor();
            // pHistClone->SetFillStyle(1001);
            // pHistClone->SetLineColorAlpha(col, 0.f);
            // pHistClone->SetFillColorAlpha(col, 0.3f);
            // pHistClone->Draw("e2 same");

            // Draw the data as points with error bars
            // pDataStatOnlyHist->Draw("e1 same");
            pFactorHist->Draw("e1 same");

            PlottingHelper::SaveCanvas(pCanvas3, prefix + "_factors");
        }
    }


    // -----------------------------------------------------------------------------------------------------------------------------
        // -----------------------------------------------------------------------------------------------------------------------------
            // -----------------------------------------------------------------------------------------------------------------------------
                // -----------------------------------------------------------------------------------------------------------------------------
                    // -----------------------------------------------------------------------------------------------------------------------------
                        // -----------------------------------------------------------------------------------------------------------------------------
                    // -----------------------------------------------------------------------------------------------------------------------------
                // -----------------------------------------------------------------------------------------------------------------------------
            // -----------------------------------------------------------------------------------------------------------------------------
        // -----------------------------------------------------------------------------------------------------------------------------
    // -----------------------------------------------------------------------------------------------------------------------------

    // -----------------------------------------------------------------------------------------------------------------------------
    // Make universe plots
    // -----------------------------------------------------------------------------------------------------------------------------
    std::cout<<"Debug point 5"<<std::endl;
    for (const auto &[selectionName, enabledMap] : config.extractXSecs.crossSectionIsEnabled)
    {
        if(selectionName=="golden") continue;

        for (const auto &[xsecName, isEnabled] : enabledMap)
        {
            // Skip cross-sections that aren't enabled
            if (!isEnabled)
                continue;
            // Define a prefix for the names of the plots
            const auto &metadata = metadataMap.at(xsecName);
            const auto selectedEventsData = ubsmear::UBFileHelper::ReadMatrix("NuWroSidebandFit_" + selectionName + "_" + xsecName + "_selectedEventsData.txt");

            // -----------------------------------------------------------------------------------------------------------------------------
            // Make the comparison plot between data and smeared prediction
            // -----------------------------------------------------------------------------------------------------------------------------
            // Now we have the data and the prediction on the same footing we can compare them on a plot!
            // Get the bin edges (excluding any underflow/overflow bins)
            const auto extendedBinEdges = metadata.GetBinEdges();
            std::vector<float> binEdges;
            for (unsigned int iBin = 0; iBin < metadata.GetNBins(); ++iBin)
            {
                // // Skip underflow/overflow bins
                // if (metadata.IsUnderOverflowBin(iBin))
                //     continue;

                // If this is the first bin then add the lower edge
                if (binEdges.empty())
                {
                    if (metadata.IsScaledByBinWidth())
                    {
                        binEdges.push_back(std::min(5.f, extendedBinEdges.at(iBin)));
                    }
                    else
                    {
                        // If we don't scale by bin width, then just use zero as the first bin edge
                        binEdges.push_back(0.f);
                    }
                }

                // Add the upper bin edge
                if (metadata.IsScaledByBinWidth())
                {
                    binEdges.push_back(std::min(5.f, extendedBinEdges.at(iBin + 1)));
                }
                else
                {
                    // If we don't scale by bin width, then just use unit width bins
                    binEdges.push_back(binEdges.back() + 1.f);
                }
            }

            for (unsigned int iUni = 0; iUni < 20; ++iUni)
            {
                const auto weightDimensions = {std::make_pair("xsec", systParams.xsecDimensions), std::make_pair("reint", systParams.reintDimensions), std::make_pair("flux", systParams.fluxDimensions)};
                for (const auto &[group, dimensions] : weightDimensions)
                {
                    for (const auto &[paramName, nUniverses] : dimensions)
                    {
                        const std::string prefix = "nuWroSidebandFitPlots_" + selectionName + "_" + xsecName + "_" + group + "_" + paramName + "_" + std::to_string(iUni);

                        try //todo improve code here
                        {
                            const auto _selectedEventsSignalTruth = ubsmear::UBFileHelper::ReadMatrix("NuWroSidebandFit_" + selectionName + "_" + xsecName  + "_" + group + "_" + paramName + "_" + std::to_string(iUni) + "_selectedEventsSignalTruth.txt");
                            const auto _selectedEventsBackgroundReco = ubsmear::UBFileHelper::ReadMatrix("NuWroSidebandFit_" + selectionName + "_" + xsecName + "_" + group + "_" + paramName + "_" + std::to_string(iUni) + "_selectedEventsBackgroundReco.txt");
                            const auto _smearingMatrix = ubsmear::UBFileHelper::ReadMatrix("NuWroSidebandFit_" + selectionName + "_" + xsecName + "_" + group + "_" + paramName + "_" + std::to_string(iUni) + "_smearingMatrix.txt");
                        }
                        catch (const std::exception &e)
                        {
                            std::cout << "Skipping " << prefix << std::endl;
                            continue;
                        }

                        const auto selectedEventsSignalTruth = ubsmear::UBFileHelper::ReadMatrix("NuWroSidebandFit_" + selectionName + "_" + xsecName  + "_" + group + "_" + paramName + "_" + std::to_string(iUni) + "_selectedEventsSignalTruth.txt");
                        const auto selectedEventsBackgroundReco = ubsmear::UBFileHelper::ReadMatrix("NuWroSidebandFit_" + selectionName + "_" + xsecName + "_" + group + "_" + paramName + "_" + std::to_string(iUni) + "_selectedEventsBackgroundReco.txt");
                        const auto smearingMatrix = ubsmear::UBFileHelper::ReadMatrix("NuWroSidebandFit_" + selectionName + "_" + xsecName + "_" + group + "_" + paramName + "_" + std::to_string(iUni) + "_smearingMatrix.txt");

                        const auto signalData = selectedEventsData - selectedEventsBackgroundReco;
                        const auto selectedEventsSignalReco = smearingMatrix * selectedEventsSignalTruth;

                        auto cc0piUniverseConstraintParam = cc0piUniverseConstraintMap.at(selectionName).at(xsecName).at(paramName).at(iUni).first;
                        auto cc0piUniverseConstraintParamError = cc0piUniverseConstraintMap.at(selectionName).at(xsecName).at(paramName).at(iUni).second;
                        const vector<float> cc0piUniverseConstraintParamFloat(cc0piUniverseConstraintParam.begin(), cc0piUniverseConstraintParam.end());
                        const ubsmear::UBMatrix sidebandParamVectorTruth(cc0piUniverseConstraintParamFloat, cc0piUniverseConstraintParam.size(), 1);
                        const auto scaledSelectedEventsSignalTruth = ElementWiseOperation(selectedEventsSignalTruth, sidebandParamVectorTruth, [](const auto &l, const auto& r) { return l * r; });
                        const auto scaledSelectedEventsSignalReco = smearingMatrix * scaledSelectedEventsSignalTruth;


                        // -----------------------------------------------------------------------------------------------------------------------------
                        // Data vs Prediction
                        // -----------------------------------------------------------------------------------------------------------------------------
                        // Setup the data histogram and the prediction histogram
                        auto pDataHist = std::make_shared<TH1F>((prefix + "_raw_data").c_str(), "", binEdges.size() - 1, binEdges.data());
                        // auto pDataStatOnlyHist = std::make_shared<TH1F>((prefix + "_dataStatOnly").c_str(), "", binEdges.size() - 1, binEdges.data());
                        auto pPredictionHist = std::make_shared<TH1F>((prefix + "_raw_prediction").c_str(), "", binEdges.size() - 1, binEdges.data());

                        // Fill the bins
                        // ATTN here we only show the diagonals of the error matrices
                        float minY = +std::numeric_limits<float>::max();
                        float maxY = -std::numeric_limits<float>::max();
                        // const auto dataStatErrorMatrix = errorMatrixMap.at("data").at("stat");
                        for (unsigned int iBin = 1; iBin <= binEdges.size() - 1; ++iBin)
                        {
                            // Set the values for the data histogram
                            const auto dataValue = signalData.At(iBin - 1, 0);
                            // const auto dataError = std::pow(totalDataErrorMatrix.At(iBin - 1, iBin - 1), 0.5f);
                            // const auto dataStatOnlyError = std::pow(dataStatErrorMatrix.At(iBin - 1, iBin - 1), 0.5f);

                            pDataHist->SetBinContent(iBin, dataValue);
                            pDataHist->SetBinError(iBin, 0);

                            // pDataStatOnlyHist->SetBinContent(iBin, dataValue);
                            // pDataStatOnlyHist->SetBinError(iBin, dataStatOnlyError);

                            // Set the values of the prediction
                            const auto predictionValue = selectedEventsSignalReco.At(iBin - 1, 0);
                            // const auto predictionError = std::pow(smearedPredictionErrorMatrix.At(iBin - 1, iBin - 1), 0.5f);

                            pPredictionHist->SetBinContent(iBin, predictionValue);
                            pPredictionHist->SetBinError(iBin, 0);

                            // For the proton multiplicity plot, use explicit bin labels
                            if (xsecName == "nProtons")
                            {
                                for (auto &pHist : {pDataHist, pPredictionHist})
                                {
                                    pHist->GetXaxis()->SetBinLabel(1, "0");
                                    pHist->GetXaxis()->SetBinLabel(2, "1");
                                    pHist->GetXaxis()->SetBinLabel(3, ">1");
                                }
                            }

                            // Get the limiting values
                            minY = std::min(minY, dataValue);
                            minY = std::min(minY, predictionValue);

                            maxY = std::max(maxY, dataValue);
                            maxY = std::max(maxY, predictionValue);
                        }

                        // Set the y-range
                        auto padding = (maxY - minY) * 0.05;
                        maxY += padding;
                        minY -= padding;
                        minY = std::max(minY, 0.f);
                        minY = 0.f; // Remove this line to get a dynamic lower y-range
                        pDataHist->GetYaxis()->SetRangeUser(minY, maxY);
                        // pDataStatOnlyHist->GetYaxis()->SetRangeUser(minY, maxY);
                        pPredictionHist->GetYaxis()->SetRangeUser(minY, maxY);

                        // Set the colours of the histograms
                        PlottingHelper::SetLineStyle(pDataHist, PlottingHelper::Primary);
                        // PlottingHelper::SetLineStyle(pDataStatOnlyHist, PlottingHelper::Primary);
                        PlottingHelper::SetLineStyle(pPredictionHist, PlottingHelper::Secondary);

                        // Make the plot!
                        auto pCanvas1 = PlottingHelper::GetCanvas();
                        gStyle->SetEndErrorSize(4);

                        // Draw the smeared prediction
                        pPredictionHist->Draw("hist");

                        // Draw the prediction uncertainties as a semi-transparent band
                        auto pHistClone = static_cast<TH1F *>(pPredictionHist->Clone());
                        auto col = pHistClone->GetLineColor();
                        pHistClone->SetFillStyle(1001);
                        pHistClone->SetLineColorAlpha(col, 0.f);
                        pHistClone->SetFillColorAlpha(col, 0.3f);
                        pHistClone->Draw("e2 same");

                        // Draw the data as points with error bars
                        // pDataStatOnlyHist->Draw("e1 same");
                        pDataHist->Draw("e1 same");

                        PlottingHelper::SaveCanvas(pCanvas1, prefix + "_raw_data-vs-smearedPrediction");




                        // -----------------------------------------------------------------------------------------------------------------------------
                        // Data vs Scaled Prediction
                        // -----------------------------------------------------------------------------------------------------------------------------
                        // Setup the data histogram and the prediction histogram
                        auto pDataHist2 = std::make_shared<TH1F>((prefix + "_raw_data2").c_str(), "", binEdges.size() - 1, binEdges.data());
                        // auto pDataStatOnlyHist = std::make_shared<TH1F>((prefix + "_dataStatOnly").c_str(), "", binEdges.size() - 1, binEdges.data());
                        auto pPredictionHist2 = std::make_shared<TH1F>((prefix + "_scaled_raw_prediction").c_str(), "", binEdges.size() - 1, binEdges.data());

                        // Fill the bins
                        // ATTN here we only show the diagonals of the error matrices
                        minY = +std::numeric_limits<float>::max();
                        maxY = -std::numeric_limits<float>::max();
                        // const auto dataStatErrorMatrix = errorMatrixMap.at("data").at("stat");
                        for (unsigned int iBin = 1; iBin <= binEdges.size() - 1; ++iBin)
                        {
                            // Set the values for the data histogram
                            const auto dataValue = signalData.At(iBin - 1, 0);
                            // const auto dataError = std::pow(totalDataErrorMatrix.At(iBin - 1, iBin - 1), 0.5f);
                            // const auto dataStatOnlyError = std::pow(dataStatErrorMatrix.At(iBin - 1, iBin - 1), 0.5f);

                            pDataHist2->SetBinContent(iBin, dataValue);
                            pDataHist2->SetBinError(iBin, 0);

                            // pDataStatOnlyHist->SetBinContent(iBin, dataValue);
                            // pDataStatOnlyHist->SetBinError(iBin, dataStatOnlyError);

                            // Set the values of the prediction
                            const auto predictionValue = scaledSelectedEventsSignalReco.At(iBin - 1, 0);
                            // const auto predictionError = std::pow(smearedPredictionErrorMatrix.At(iBin - 1, iBin - 1), 0.5f);

                            pPredictionHist2->SetBinContent(iBin, predictionValue);
                            pPredictionHist2->SetBinError(iBin, 0);

                            // For the proton multiplicity plot, use explicit bin labels
                            if (xsecName == "nProtons")
                            {
                                for (auto &pHist2 : {pDataHist2, pPredictionHist2})
                                {
                                    pHist2->GetXaxis()->SetBinLabel(1, "0");
                                    pHist2->GetXaxis()->SetBinLabel(2, "1");
                                    pHist2->GetXaxis()->SetBinLabel(3, ">1");
                                }
                            }

                            // Get the limiting values
                            minY = std::min(minY, dataValue);
                            minY = std::min(minY, predictionValue);

                            maxY = std::max(maxY, dataValue);
                            maxY = std::max(maxY, predictionValue);
                        }

                        // Set the y-range
                        padding = (maxY - minY) * 0.05;
                        maxY += padding;
                        minY -= padding;
                        minY = std::max(minY, 0.f);
                        minY = 0.f; // Remove this line to get a dynamic lower y-range
                        pDataHist2->GetYaxis()->SetRangeUser(minY, maxY);
                        // pDataStatOnlyHist->GetYaxis()->SetRangeUser(minY, maxY);
                        pPredictionHist2->GetYaxis()->SetRangeUser(minY, maxY);

                        // Set the colours of the histograms
                        PlottingHelper::SetLineStyle(pDataHist2, PlottingHelper::Primary);
                        // PlottingHelper::SetLineStyle(pDataStatOnlyHist, PlottingHelper::Primary);
                        PlottingHelper::SetLineStyle(pPredictionHist2, PlottingHelper::Secondary);

                        // Make the plot!
                        auto pCanvas2 = PlottingHelper::GetCanvas();
                        gStyle->SetEndErrorSize(4);

                        // Draw the smeared prediction
                        pPredictionHist2->Draw("hist");

                        // Draw the prediction uncertainties as a semi-transparent band
                        auto pHistClone2 = static_cast<TH1F *>(pPredictionHist2->Clone());
                        col = pHistClone2->GetLineColor();
                        pHistClone2->SetFillStyle(1001);
                        pHistClone2->SetLineColorAlpha(col, 0.f);
                        pHistClone2->SetFillColorAlpha(col, 0.3f);
                        pHistClone2->Draw("e2 same");

                        // Draw the data as points with error bars
                        // pDataStatOnlyHist->Draw("e1 same");
                        pDataHist2->Draw("e1 same");

                        PlottingHelper::SaveCanvas(pCanvas2, prefix + "_scaled_raw_data-vs-smearedPrediction");


                        // -----------------------------------------------------------------------------------------------------------------------------
                        // Scaling factors
                        // -----------------------------------------------------------------------------------------------------------------------------
                        // Setup the data histogram and the prediction histogram
                        auto pFactorHist = std::make_shared<TH1F>((prefix + "_factors").c_str(), "", binEdges.size() - 1, binEdges.data());
                        // auto pDataStatOnlyHist = std::make_shared<TH1F>((prefix + "_dataStatOnly").c_str(), "", binEdges.size() - 1, binEdges.data());
                        // auto pPredictionHist = std::make_shared<TH1F>((prefix + "_prediction").c_str(), "", binEdges.size() - 1, binEdges.data());

                        // Fill the bins
                        // ATTN here we only show the diagonals of the error matrices
                        minY = +std::numeric_limits<float>::max();
                        maxY = -std::numeric_limits<float>::max();
                        // const auto dataStatErrorMatrix = errorMatrixMap.at("data").at("stat");
                        for (unsigned int iBin = 1; iBin <= binEdges.size() - 1; ++iBin)
                        {
                            // Set the values for the data histogram
                            const auto factorValue = (float) cc0piUniverseConstraintParam.at(iBin - 1);
                            const auto factorError = (float) cc0piUniverseConstraintParamError.at(iBin - 1);//std::pow( ... , 0.5f);
                            // const auto dataStatOnlyError = std::pow(dataStatErrorMatrix.At(iBin - 1, iBin - 1), 0.5f);

                            pFactorHist->SetBinContent(iBin, factorValue);
                            pFactorHist->SetBinError(iBin, factorError);

                            // pDataStatOnlyHist->SetBinContent(iBin, dataValue);
                            // pDataStatOnlyHist->SetBinError(iBin, dataStatOnlyError);

                            // Set the values of the prediction
                            // const auto predictionValue = smearedPrediction.At(iBin - 1, 0);
                            // const auto predictionError = std::pow(smearedPredictionErrorMatrix.At(iBin - 1, iBin - 1), 0.5f);

                            // pPredictionHist->SetBinContent(iBin, predictionValue);
                            // pPredictionHist->SetBinError(iBin, predictionError);

                            // For the proton multiplicity plot, use explicit bin labels
                            if (xsecName == "nProtons")
                            {
                                for (auto &pHist : {pFactorHist})//, pDataStatOnlyHist, pPredictionHist})
                                {
                                    pHist->GetXaxis()->SetBinLabel(1, "0");
                                    pHist->GetXaxis()->SetBinLabel(2, "1");
                                    pHist->GetXaxis()->SetBinLabel(3, ">1");
                                }
                            }

                            // Get the limiting values
                            minY = std::min(minY, factorValue - factorError);
                            // minY = std::min(minY, predictionValue - predictionError);

                            maxY = std::max(maxY, factorValue + factorError);
                            // maxY = std::max(maxY, predictionValue + predictionError);
                        }

                        // Set the y-range
                        padding = (maxY - minY) * 0.05;
                        maxY += padding;
                        minY -= padding;
                        minY = std::max(minY, 0.f);
                        minY = 0.f; // Remove this line to get a dynamic lower y-range
                        pFactorHist->GetYaxis()->SetRangeUser(minY, maxY);
                        // pDataStatOnlyHist->GetYaxis()->SetRangeUser(minY, maxY);
                        // pPredictionHist->GetYaxis()->SetRangeUser(minY, maxY);

                        // Set the colours of the histograms
                        PlottingHelper::SetLineStyle(pFactorHist, PlottingHelper::Primary);
                        // PlottingHelper::SetLineStyle(pDataStatOnlyHist, PlottingHelper::Primary);
                        // PlottingHelper::SetLineStyle(pPredictionHist, PlottingHelper::Secondary);

                        // Make the plot!
                        auto pCanvas3 = PlottingHelper::GetCanvas();
                        gStyle->SetEndErrorSize(4);

                        // // Draw the smeared prediction
                        // pPredictionHist->Draw("hist");

                        // // Draw the prediction uncertainties as a semi-transparent band
                        // auto pHistClone = static_cast<TH1F *>(pPredictionHist->Clone());
                        // const auto col = pHistClone->GetLineColor();
                        // pHistClone->SetFillStyle(1001);
                        // pHistClone->SetLineColorAlpha(col, 0.f);
                        // pHistClone->SetFillColorAlpha(col, 0.3f);
                        // pHistClone->Draw("e2 same");

                        // Draw the data as points with error bars
                        // pDataStatOnlyHist->Draw("e1 same");
                        pFactorHist->Draw("e1 same");

                        PlottingHelper::SaveCanvas(pCanvas3, prefix + "_factors");
                    }
                }
            }
        }
    }

}

} // namespace ubcc1pi_macros

