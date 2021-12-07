/**
 *  @file  ubcc1pi_standalone/Macros/MakeXSecPlots.cxx
 *
 *  @brief The implementation file of the MakeXSecPlots macro
 */

#include "ubcc1pi_standalone/Macros/Macros.h"
#include "ubcc1pi_standalone/Objects/FileReader.h"

#include "ubcc1pi_standalone/Helpers/PlottingHelper.h"
#include "ubcc1pi_standalone/Helpers/FormattingHelper.h"
#include "ubcc1pi_standalone/Helpers/CrossSectionHelper.h"

#include "ubsmear.h"

#include <TStyle.h>

using namespace ubcc1pi;

namespace ubcc1pi_macros
{

void MakeXSecPlots(const Config &config)
{
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

    // -------------------------------------------------------------------------------------------------------------------------------------
    // Setup a table to hold the chi2 values between the data and prediction for each cross-section
    // -------------------------------------------------------------------------------------------------------------------------------------
    FormattingHelper::Table table({"Selection", "Cross-section", "nBins", "", "Chi2", "DoF", "Chi2/DoF", "p-value"});

    // -------------------------------------------------------------------------------------------------------------------------------------
    // Loop over all possible cross-sections
    // -------------------------------------------------------------------------------------------------------------------------------------
    for (const auto &[selectionName, enabledMap] : config.extractXSecs.crossSectionIsEnabled)
    {
        for (const auto &[xsecName, isEnabled] : enabledMap)
        {
            // Skip cross-sections that aren't enabled
            if (!isEnabled)
                continue;

            std::cout << selectionName << " - " << xsecName << std::endl;

            // Define a lambda function to get a matrix from a file for brevity
            const auto getMatrix = [&, selectionName = selectionName, xsecName = xsecName](const std::string &identifier) -> ubsmear::UBMatrix {
                return ubsmear::UBFileHelper::ReadMatrix("xsec_" + selectionName + "_" + xsecName + "_" + identifier + ".txt");
            };

            // Define a lambda function to get a matrix from a file an trim any overflow / underflow bins
            const auto &metadata = metadataMap.at(xsecName);
            const auto getTrimmedMatrix = [&] (const std::string &identifier) -> ubsmear::UBMatrix {
                return ubsmear::UBSmearingHelper::TrimUnderOverflowBins(getMatrix(identifier), metadata);
            };

            // Get the data cross-section
            const auto data = getTrimmedMatrix("data");

            // Get the smearing matrix
            const auto smearingMatrix = getMatrix("smearingMatrix");

            // Build the total error matrix for each group of systematic parameters
            // Here errorMatrixMap is indexed by [quantity][group], where quantity = "data" or "smearingMatrix"
            // Here totalErrorMatrixMap is index by [quantity] and contains the sum (over all groups) of the entries in errorMatrixMap
            std::map<std::string, std::map<std::string, ubsmear::UBMatrix> > errorMatrixMap;
            std::map<std::string, ubsmear::UBMatrix> totalErrorMatrixMap;

            // -----------------------------------------------------------------------------------------------------------------------------
            // Read the matrices from disk
            // -----------------------------------------------------------------------------------------------------------------------------
            // Loop over the two quantities that have systematic uncertainties
            for (const std::string &quantity : {"data", "smearingMatrix"})
            {
                // Get the number of bins (for the smearing matrix there are N^2 bins when flattened)
                const auto nBins = (quantity == "data" ? data.GetRows() : std::pow(smearingMatrix.GetRows(), 2));

                // Define the function to use when getting a matrix (i.e. should we trim the under/overflow bins or not?)
                const auto &getMatrixFunction = (quantity == "data"
                    ? std::function<ubsmear::UBMatrix(const std::string&)>(getTrimmedMatrix)
                    : std::function<ubsmear::UBMatrix(const std::string&)>(getMatrix));

                // Setup an empty error matrix for this quantity
                auto errorMatrixTotalSum =  ubsmear::UBMatrixHelper::GetZeroMatrix(nBins, nBins);

                // For the data cross-section, we also have a stat uncertainty
                // For the smearing matrix, just use a zero vector
                const auto statUncertainties = (quantity == "data") ? getMatrixFunction("data_stat") : ubsmear::UBMatrixHelper::GetZeroMatrix(nBins, 1);

                // To convert these stat uncertainties into an error matrix, we produce a diagonal matrix whose diagonal entries contain the
                // variances (i.e. square of the stat uncerainties)
                const auto statVariances = ubsmear::ElementWiseOperation(statUncertainties, statUncertainties, [](const auto &l, const auto &r) { return l * r; });
                const auto statErrorMatrix = ubsmear::UBMatrixHelper::GetDiagonalMatrix(statVariances);

                // Store this in the map
                errorMatrixMap[quantity].emplace("stat", statErrorMatrix);

                // Add this to the grand total error matrix
                errorMatrixTotalSum = errorMatrixTotalSum + statErrorMatrix;

                // Handle the multisim parameters
                for (const auto &[group, dimensions] : std::map<std::string, CrossSectionHelper::SystDimensionsMap>(
                    {
                        {"flux", config.extractXSecs.fluxDimensions},
                        {"xsec", config.extractXSecs.xsecDimensions},
                        {"reint", config.extractXSecs.reintDimensions},
                        {"misc", {
                            {"bootstrap", config.extractXSecs.nBootstrapUniverses},
                            {"dirt", 2},
                            {"POT", 0}
                        }}
                    }))
                {
                    // Setup an empty error matrix for this group
                    auto errorMatrixTotal =  ubsmear::UBMatrixHelper::GetZeroMatrix(nBins, nBins);

                    // Loop over the parameters in this group
                    for (const auto &[paramName, nUniverses] : dimensions)
                    {
                        // Get the bias vector and covariance matrix
                        const auto biasVector = getMatrixFunction(quantity + "_" + group + "_" + paramName + "_bias");
                        const auto covarianceMatrix = getMatrixFunction(quantity + "_" + group + "_" + paramName + "_covariance");

                        // Get the total error matrix from the bias and covariance
                        const auto errorMatrix = CrossSectionHelper::GetErrorMatrix(biasVector, covarianceMatrix);

                        // Add this error matrix to the total
                        errorMatrixTotal = errorMatrixTotal + errorMatrix;
                    }

                    // Store this total in the map
                    errorMatrixMap[quantity].emplace(group, errorMatrixTotal);

                    // Add this total to the grand total error matrix
                    errorMatrixTotalSum = errorMatrixTotalSum + errorMatrixTotal;
                }

                // Handle the unisim parameters
                for (const auto &[group, dimensions] : std::map<std::string, CrossSectionHelper::SystUnisimDimensionsMap>(
                    {
                        {"detector", config.extractXSecs.detVarDimensions}
                    }))
                {
                    // Setup an empty error matrix for this group
                    auto errorMatrixTotal =  ubsmear::UBMatrixHelper::GetZeroMatrix(nBins, nBins);

                    // Loop over the parameters in this group
                    for (const auto &[paramName, cvName] : dimensions)
                    {
                        // Get the bias vector and (dummy) covariance matrix
                        const auto biasVector = getMatrixFunction(quantity + "_" + group + "_" + paramName + "_bias");
                        const auto covarianceMatrix = getMatrixFunction(quantity + "_" + group + "_" + paramName + "_covariance");

                        // Get the total error matrix from the bias and covariance
                        const auto errorMatrix = CrossSectionHelper::GetErrorMatrix(biasVector, covarianceMatrix);

                        // Add this error matrix to the total
                        errorMatrixTotal = errorMatrixTotal + errorMatrix;
                    }

                    // Store this total in the map
                    errorMatrixMap[quantity].emplace(group, errorMatrixTotal);

                    // Add this total to the grand total error matrix
                    errorMatrixTotalSum = errorMatrixTotalSum + errorMatrixTotal;
                }

                // Add the grand total to the map
                totalErrorMatrixMap.emplace(quantity, errorMatrixTotalSum);
            }

            // Define a prefix for the names of the plots
            const std::string prefix = "xsecPlots_" + selectionName + "_" + xsecName;

            // -----------------------------------------------------------------------------------------------------------------------------
            // Plot the error matrices
            // -----------------------------------------------------------------------------------------------------------------------------
            for (const auto &[quantity, groupToMatrixMap] : errorMatrixMap)
            {
                // Get the quantity in question as a vector
                // ATTN here we flatten the smearing matrix to a column vector
                const auto quantityVector = (quantity == "data" ? data : ubsmear::UBSmearingHelper::Flatten(smearingMatrix));

                // Plot the total error matrix (summed over all groups)
                const auto totalErrorMatrix = totalErrorMatrixMap.at(quantity);
                PlottingHelper::PlotErrorMatrix(totalErrorMatrix, prefix + "_" + quantity + "_totalErrorMatrix", metadata);
                PlottingHelper::PlotFractionalErrorMatrix(totalErrorMatrix, quantityVector, prefix + "_" + quantity + "_totalFracErrorMatrix", metadata);

                // Plot the total error matrix for each group individually
                for (const auto &[group, errorMatrix] : groupToMatrixMap)
                {
                    PlottingHelper::PlotErrorMatrix(errorMatrix, prefix + "_" + quantity + "_" + group + "_totalErrorMatrix", metadata);
                    PlottingHelper::PlotFractionalErrorMatrix(errorMatrix, quantityVector, prefix + "_" + quantity + "_" + group + "_totalFracErrorMatrix", metadata);
                }
            }

            // -----------------------------------------------------------------------------------------------------------------------------
            // Forward-fold the prediction and plot it
            // -----------------------------------------------------------------------------------------------------------------------------
            // for now we just re-use the error matrix function (and use the default palette)
            PlottingHelper::PlotErrorMatrix(smearingMatrix, prefix + "_smearingMatrix", metadata, true, false);

            // Now get the predicted cross-section and it's error matrix
            const auto prediction = getMatrix("prediction");
            const auto predictionBiasVector = getMatrix("prediction_stat_bias");
            const auto predictionCovarianceMatrix = getMatrix("prediction_stat_covariance");
            const auto predictionErrorMatrix = CrossSectionHelper::GetErrorMatrix(predictionBiasVector, predictionCovarianceMatrix);

            // Plot the error matrix on the prediction
            PlottingHelper::PlotErrorMatrix(predictionErrorMatrix, prefix + "_prediction_stat_totalErrorMatrix", metadata);
            PlottingHelper::PlotFractionalErrorMatrix(predictionErrorMatrix, prediction, prefix + "_prediction_stat_totalFracErrorMatrix", metadata);

            // Now smear the prediction so it can be compared to the data
            std::cout << "Forward folding prediction" << std::endl;
            const auto &[smearedPrediction, smearedPredictionErrorMatrix] = ubsmear::ForwardFold(
                metadata,                                                  // The metadata that defines the binning
                prediction, predictionErrorMatrix,                         // The prediction and it's error matrix
                smearingMatrix, totalErrorMatrixMap.at("smearingMatrix"),  // The smearing matrix and it's error matrix
                config.makeXSecPlots.nUniverses,                           // The number of universes to use when propagating the uncertainties
                config.makeXSecPlots.precision);                           // The precision to use when finding eigenvalues and eigenvectors

            // Plot the error matrix on the smeared prediction
            PlottingHelper::PlotErrorMatrix(smearedPredictionErrorMatrix, prefix + "_smearedPrediction_totalErrorMatrix", metadata);
            PlottingHelper::PlotFractionalErrorMatrix(smearedPredictionErrorMatrix, smearedPrediction, prefix + "_smearedPrediction_totalFracErrorMatrix", metadata);

            // Save the forward-folded prediction
            FormattingHelper::SaveMatrix(smearedPrediction, "xsec_" + selectionName + "_" + xsecName + "_forwardFoldedPrediction.txt");
            FormattingHelper::SaveMatrix(smearedPredictionErrorMatrix, "xsec_" + selectionName + "_" + xsecName + "_forwardFoldedPrediction_error.txt");

            // -----------------------------------------------------------------------------------------------------------------------------
            // Get the chi2 for the data-prediction comparison and add it to the table
            // -----------------------------------------------------------------------------------------------------------------------------
            const auto totalDataErrorMatrix = totalErrorMatrixMap.at("data");
            std::cout << "Getting chi2" << std::endl;
            const auto &[chi2, degreesOfFreedom] = ubsmear::GetChi2(smearedPrediction, smearedPredictionErrorMatrix, data, totalDataErrorMatrix, config.makeXSecPlots.precision);

            std::cout << "Getting p-value" << std::endl;
            const auto pValue = ubsmear::UBStatisticsHelper::GetPValue(chi2, degreesOfFreedom, config.makeXSecPlots.precision);

            table.AddEmptyRow();
            table.SetEntry("Selection", selectionName);
            table.SetEntry("Cross-section", xsecName);
            table.SetEntry("nBins", data.GetRows());
            table.SetEntry("Chi2", chi2);
            table.SetEntry("DoF", degreesOfFreedom);
            table.SetEntry("Chi2/DoF", chi2 / static_cast<float>(degreesOfFreedom));
            table.SetEntry("p-value", pValue);

            // -----------------------------------------------------------------------------------------------------------------------------
            // Make the comparison plot between data and smeared prediction
            // -----------------------------------------------------------------------------------------------------------------------------
            // Now we have the data and the prediction on the same footing we can compare them on a plot!
            // Get the bin edges (excluding any underflow/overflow bins)
            const auto extendedBinEdges = metadata.GetBinEdges();
            std::vector<float> binEdges;
            for (unsigned int iBin = 0; iBin < metadata.GetNBins(); ++iBin)
            {
                // Skip underflow/overflow bins
                if (metadata.IsUnderOverflowBin(iBin))
                    continue;

                // If this is the first bin then add the lower edge
                if (binEdges.empty())
                {
                    if (metadata.IsScaledByBinWidth())
                    {
                        binEdges.push_back(extendedBinEdges.at(iBin));
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
                    binEdges.push_back(extendedBinEdges.at(iBin + 1));
                }
                else
                {
                    // If we don't scale by bin width, then just use unit width bins
                    binEdges.push_back(binEdges.back() + 1.f);
                }
            }

            // Setup the data histogram and the prediction histogram
            auto pDataHist = std::make_shared<TH1F>((prefix + "_data").c_str(), "", binEdges.size() - 1, binEdges.data());
            auto pDataStatOnlyHist = std::make_shared<TH1F>((prefix + "_dataStatOnly").c_str(), "", binEdges.size() - 1, binEdges.data());
            auto pPredictionHist = std::make_shared<TH1F>((prefix + "_prediction").c_str(), "", binEdges.size() - 1, binEdges.data());

            // Fill the bins
            // ATTN here we only show the diagonals of the error matrices
            float minY = +std::numeric_limits<float>::max();
            float maxY = -std::numeric_limits<float>::max();
            const auto dataStatErrorMatrix = errorMatrixMap.at("data").at("stat");
            for (unsigned int iBin = 1; iBin <= binEdges.size() - 1; ++iBin)
            {
                // Set the values for the data histogram
                const auto dataValue = data.At(iBin - 1, 0);
                const auto dataError = std::pow(totalDataErrorMatrix.At(iBin - 1, iBin - 1), 0.5f);
                const auto dataStatOnlyError = std::pow(dataStatErrorMatrix.At(iBin - 1, iBin - 1), 0.5f);

                pDataHist->SetBinContent(iBin, dataValue);
                pDataHist->SetBinError(iBin, dataError);

                pDataStatOnlyHist->SetBinContent(iBin, dataValue);
                pDataStatOnlyHist->SetBinError(iBin, dataStatOnlyError);

                // Set the values of the prediction
                const auto predictionValue = smearedPrediction.At(iBin - 1, 0);
                const auto predictionError = std::pow(smearedPredictionErrorMatrix.At(iBin - 1, iBin - 1), 0.5f);

                pPredictionHist->SetBinContent(iBin, predictionValue);
                pPredictionHist->SetBinError(iBin, predictionError);

                // For the proton multiplicity plot, use explicit bin labels
                if (xsecName == "nProtons")
                {
                    for (auto &pHist : {pDataHist, pDataStatOnlyHist, pPredictionHist})
                    {
                        pHist->GetXaxis()->SetBinLabel(1, "0");
                        pHist->GetXaxis()->SetBinLabel(2, "1");
                        pHist->GetXaxis()->SetBinLabel(3, ">1");
                    }
                }

                // Get the limiting values
                minY = std::min(minY, dataValue - dataError);
                minY = std::min(minY, predictionValue - predictionError);

                maxY = std::max(maxY, dataValue + dataError);
                maxY = std::max(maxY, predictionValue + predictionError);
            }

            // Set the y-range
            const auto padding = (maxY - minY) * 0.05;
            maxY += padding;
            minY -= padding;
            minY = std::max(minY, 0.f);
            minY = 0.f; // Remove this line to get a dynamic lower y-range
            pDataHist->GetYaxis()->SetRangeUser(minY, maxY);
            pDataStatOnlyHist->GetYaxis()->SetRangeUser(minY, maxY);
            pPredictionHist->GetYaxis()->SetRangeUser(minY, maxY);

            // Set the colours of the histograms
            PlottingHelper::SetLineStyle(pDataHist, PlottingHelper::Primary);
            PlottingHelper::SetLineStyle(pDataStatOnlyHist, PlottingHelper::Primary);
            PlottingHelper::SetLineStyle(pPredictionHist, PlottingHelper::Secondary);

            // Make the plot!
            auto pCanvas = PlottingHelper::GetCanvas();
            gStyle->SetEndErrorSize(4);

            // Draw the smeared prediction
            pPredictionHist->Draw("hist");

            // Draw the prediction uncertainties as a semi-transparent band
            auto pHistClone = static_cast<TH1F *>(pPredictionHist->Clone());
            const auto col = pHistClone->GetLineColor();
            pHistClone->SetFillStyle(1001);
            pHistClone->SetLineColorAlpha(col, 0.f);
            pHistClone->SetFillColorAlpha(col, 0.3f);
            pHistClone->Draw("e2 same");

            // Draw the data as points with error bars
            pDataStatOnlyHist->Draw("e1 same");
            pDataHist->Draw("e1 same");

            PlottingHelper::SaveCanvas(pCanvas, prefix + "_data-vs-smearedPrediction");
        }
    }

    // Save the table
    table.WriteToFile("xsecPlots_goodnessOfFitStatistics.md");
}

} // namespace ubcc1pi_macros
