/**
 *  @file  ubcc1pi_standalone/Macros/MakeNuWroXSecPlots.cxx
 *
 *  @brief The implementation file of the MakeNuWroXSecPlots macro
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

void MakeNuWroXSecPlots(const Config &config)
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
    FormattingHelper::Table tableScaled({"Selection", "Cross-section", "nBins", "", "Chi2", "DoF", "Chi2/DoF", "p-value"});
    FormattingHelper::Table tableUnscaled({"Selection", "Cross-section", "nBins", "", "Chi2", "DoF", "Chi2/DoF", "p-value"});

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
                return ubsmear::UBFileHelper::ReadMatrix("xsecNuWro_" + selectionName + "_" + xsecName + "_" + identifier + ".txt");
            };

            const auto getMatrixTrue = [&, selectionName = selectionName, xsecName = xsecName](const std::string &identifier) -> ubsmear::UBMatrix {
                return ubsmear::UBFileHelper::ReadMatrix("xsecNuWro_true_" + selectionName + "_" + xsecName + "_" + identifier + ".txt");
            };

            // Define a lambda function to get a matrix from a file an trim any overflow / underflow bins
            const auto &metadata = metadataMap.at(xsecName);
            const auto getTrimmedMatrix = [&] (const std::string &identifier) -> ubsmear::UBMatrix {
                return ubsmear::UBSmearingHelper::TrimUnderOverflowBins(getMatrix(identifier), metadata);
            };

            // Define a lambda function to get a matrix from a file an trim any overflow / underflow bins
            // const auto &metadataTrue = metadataMap.at(xsecName);
            const auto getTrimmedMatrixTrue = [&] (const std::string &identifier) -> ubsmear::UBMatrix {
                return ubsmear::UBSmearingHelper::TrimUnderOverflowBins(getMatrixTrue(identifier), metadata);
            };

            // Get the data cross-section
            const auto dataScaled = getTrimmedMatrix("data_scaled");
            const auto dataUnscaled = getTrimmedMatrix("data_unscaled");
            // const auto dataTrue = getTrimmedMatrixTrue("data"); //getTrimmedMatrixTrue("data");

            // Get the smearing matrix
            const auto smearingMatrixScaled = getMatrix("smearingMatrix_scaled");
            const auto smearingMatrixUnscaled = getMatrix("smearingMatrix_unscaled");
            const auto smearingMatrixNuwro = getMatrixTrue("smearingMatrix");

            // Build the total error matrix for each group of systematic parameters
            // Here errorMatrixMap is indexed by [quantity][group], where quantity = "data" or "smearingMatrix"
            // Here totalErrorMatrixMap is index by [quantity] and contains the sum (over all groups) of the entries in errorMatrixMap
            std::map<std::string, std::map<std::string, ubsmear::UBMatrix> > errorMatrixMapScaled;
            std::map<std::string, ubsmear::UBMatrix> totalErrorMatrixMapScaled;
            std::map<std::string, std::map<std::string, ubsmear::UBMatrix> > errorMatrixMapUnscaled;
            std::map<std::string, ubsmear::UBMatrix> totalErrorMatrixMapUnscaled;

            // -----------------------------------------------------------------------------------------------------------------------------
            // Read the matrices from disk
            // -----------------------------------------------------------------------------------------------------------------------------
            // Loop over the two quantities that have systematic uncertainties
            for (const std::string &quantity : {"data", "smearingMatrix"})
            {
                // Get the number of bins (for the smearing matrix there are N^2 bins when flattened)
                const auto nBins = (quantity == "data" ? dataUnscaled.GetRows() : std::pow(smearingMatrixUnscaled.GetRows(), 2));

                // Define the function to use when getting a matrix (i.e. should we trim the under/overflow bins or not?)
                const auto &getMatrixFunction = (quantity == "data"
                    ? std::function<ubsmear::UBMatrix(const std::string&)>(getTrimmedMatrix)
                    : std::function<ubsmear::UBMatrix(const std::string&)>(getMatrix));

                // Setup an empty error matrix for this quantity
                auto errorMatrixTotalSumScaled =  ubsmear::UBMatrixHelper::GetZeroMatrix(nBins, nBins);
                auto errorMatrixTotalSumUnscaled =  ubsmear::UBMatrixHelper::GetZeroMatrix(nBins, nBins);

                // For the data cross-section, we also have a stat uncertainty
                // For the smearing matrix, just use a zero vector
                const auto statUncertaintiesScaled = (quantity == "data") ? getMatrixFunction("data_stat_scaled") : ubsmear::UBMatrixHelper::GetZeroMatrix(nBins, 1);
                const auto statUncertaintiesUnscaled = (quantity == "data") ? getMatrixFunction("data_stat_unscaled") : ubsmear::UBMatrixHelper::GetZeroMatrix(nBins, 1);

                // To convert these stat uncertainties into an error matrix, we produce a diagonal matrix whose diagonal entries contain the
                // variances (i.e. square of the stat uncerainties)
                const auto statVariancesScaled = ubsmear::ElementWiseOperation(statUncertaintiesScaled, statUncertaintiesScaled, [](const auto &l, const auto &r) { return l * r; });
                const auto statErrorMatrixScaled = ubsmear::UBMatrixHelper::GetDiagonalMatrix(statVariancesScaled);
                const auto statVariancesUnscaled = ubsmear::ElementWiseOperation(statUncertaintiesUnscaled, statUncertaintiesUnscaled, [](const auto &l, const auto &r) { return l * r; });
                const auto statErrorMatrixUnscaled = ubsmear::UBMatrixHelper::GetDiagonalMatrix(statVariancesUnscaled);

                // Store this in the map
                errorMatrixMapScaled[quantity].emplace("stat", statErrorMatrixScaled);
                errorMatrixMapUnscaled[quantity].emplace("stat", statErrorMatrixUnscaled);

                // Add this to the grand total error matrix
                errorMatrixTotalSumScaled = errorMatrixTotalSumScaled + statErrorMatrixScaled;
                errorMatrixTotalSumUnscaled = errorMatrixTotalSumUnscaled + statErrorMatrixUnscaled;

                // Handle the multisim parameters
                for (const auto &[group, dimensions] : std::map<std::string, CrossSectionHelper::SystDimensionsMap>(
                    {
                        {"flux", config.extractXSecs.fluxDimensions},
                        {"xsec", config.extractXSecs.xsecDimensions},
                        {"reint", config.extractXSecs.reintDimensions},
                        {"misc", {
                            {"bootstrap", config.extractXSecs.nBootstrapUniverses},
                            {"sidebandWeights", config.extractXSecs.nBootstrapUniverses},
                            {"dirt", 2},
                            {"POT", 0}
                        }}
                    }))
                {
                    // Setup an empty error matrix for this group
                    auto errorMatrixTotalScaled =  ubsmear::UBMatrixHelper::GetZeroMatrix(nBins, nBins);
                    auto errorMatrixTotalUnscaled =  ubsmear::UBMatrixHelper::GetZeroMatrix(nBins, nBins);

                    // Loop over the parameters in this group
                    for (const auto &[paramName, nUniverses] : dimensions)
                    {
                        // Get the bias vector and covariance matrix
                        const auto biasVectorScaled = getMatrixFunction(quantity + "_" + group + "_" + paramName + "_bias_scaled");
                        const auto covarianceMatrixScaled = getMatrixFunction(quantity + "_" + group + "_" + paramName + "_covariance_scaled");

                        const auto biasVectorUnscaled = getMatrixFunction(quantity + "_" + group + "_" + paramName + "_bias_unscaled");
                        const auto covarianceMatrixUnscaled = getMatrixFunction(quantity + "_" + group + "_" + paramName + "_covariance_unscaled");

                        // Get the total error matrix from the bias and covariance
                        const auto errorMatrixScaled = CrossSectionHelper::GetErrorMatrix(biasVectorScaled, covarianceMatrixScaled);
                        const auto errorMatrixUnscaled = CrossSectionHelper::GetErrorMatrix(biasVectorUnscaled, covarianceMatrixUnscaled);

                        // Add this error matrix to the total
                        errorMatrixTotalScaled = errorMatrixTotalScaled + errorMatrixScaled;
                        errorMatrixTotalUnscaled = errorMatrixTotalUnscaled + errorMatrixUnscaled;
                    }

                    // Store this total in the map
                    errorMatrixMapScaled[quantity].emplace(group, errorMatrixTotalScaled);
                    errorMatrixMapUnscaled[quantity].emplace(group, errorMatrixTotalUnscaled);

                    // Add this total to the grand total error matrix
                    errorMatrixTotalSumScaled = errorMatrixTotalSumScaled + errorMatrixTotalScaled;
                    errorMatrixTotalSumUnscaled = errorMatrixTotalSumUnscaled + errorMatrixTotalUnscaled;
                }

                // Handle the unisim parameters
                for (const auto &[group, dimensions] : std::map<std::string, CrossSectionHelper::SystUnisimDimensionsMap>(
                    {
                        {"detector", config.extractXSecs.detVarDimensions}
                    }))
                {
                    // Setup an empty error matrix for this group
                    auto errorMatrixTotalScaled = ubsmear::UBMatrixHelper::GetZeroMatrix(nBins, nBins);
                    auto errorMatrixTotalUnscaled = ubsmear::UBMatrixHelper::GetZeroMatrix(nBins, nBins);

                    // Loop over the parameters in this group
                    for (const auto &[paramName, cvName] : dimensions)
                    {
                        // Get the bias vector and (dummy) covariance matrix
                        const auto biasVectorScaled = getMatrixFunction(quantity + "_" + group + "_" + paramName + "_bias_scaled");
                        const auto covarianceMatrixScaled = getMatrixFunction(quantity + "_" + group + "_" + paramName + "_covariance_scaled");

                        const auto biasVectorUnscaled = getMatrixFunction(quantity + "_" + group + "_" + paramName + "_bias_unscaled");
                        const auto covarianceMatrixUnscaled = getMatrixFunction(quantity + "_" + group + "_" + paramName + "_covariance_unscaled");

                        // Get the total error matrix from the bias and covariance
                        const auto errorMatrixScaled = CrossSectionHelper::GetErrorMatrix(biasVectorScaled, covarianceMatrixScaled);
                        const auto errorMatrixUnscaled = CrossSectionHelper::GetErrorMatrix(biasVectorUnscaled, covarianceMatrixUnscaled);

                        // Add this error matrix to the total
                        errorMatrixTotalScaled = errorMatrixTotalScaled + errorMatrixScaled;
                        errorMatrixTotalUnscaled = errorMatrixTotalUnscaled + errorMatrixUnscaled;
                    }

                    // Store this total in the map
                    errorMatrixMapScaled[quantity].emplace(group, errorMatrixTotalScaled);
                    errorMatrixMapUnscaled[quantity].emplace(group, errorMatrixTotalUnscaled);

                    // Add this total to the grand total error matrix
                    errorMatrixTotalSumScaled = errorMatrixTotalSumScaled + errorMatrixTotalScaled;
                    errorMatrixTotalSumUnscaled = errorMatrixTotalSumUnscaled + errorMatrixTotalUnscaled;
                }

                // Add the grand total to the map
                totalErrorMatrixMapScaled.emplace(quantity, errorMatrixTotalSumScaled);
                totalErrorMatrixMapUnscaled.emplace(quantity, errorMatrixTotalSumUnscaled);
            }

            // Define a prefix for the names of the plots
            const std::string prefix = "NuWroXSecPlots_" + selectionName + "_" + xsecName;

            // -----------------------------------------------------------------------------------------------------------------------------
            // Plot the error matrices
            // -----------------------------------------------------------------------------------------------------------------------------
            for (const auto &[quantity, groupToMatrixMap] : errorMatrixMapScaled)
            {
                // Get the quantity in question as a vector
                // ATTN here we flatten the smearing matrix to a column vector
                const auto quantityVector = (quantity == "data" ? dataScaled : ubsmear::UBSmearingHelper::Flatten(smearingMatrixScaled));

                // Plot the total error matrix (summed over all groups)
                const auto totalErrorMatrix = totalErrorMatrixMapScaled.at(quantity);
                PlottingHelper::PlotErrorMatrix(totalErrorMatrix, prefix + "_" + quantity + "_totalErrorMatrix_scaled", metadata);
                PlottingHelper::PlotFractionalErrorMatrix(totalErrorMatrix, quantityVector, prefix + "_" + quantity + "_totalFracErrorMatrix_scaled", metadata);

                // Plot the total error matrix for each group individually
                for (const auto &[group, errorMatrix] : groupToMatrixMap)
                {
                    PlottingHelper::PlotErrorMatrix(errorMatrix, prefix + "_" + quantity + "_" + group + "_totalErrorMatrix_scaled", metadata);
                    PlottingHelper::PlotFractionalErrorMatrix(errorMatrix, quantityVector, prefix + "_" + quantity + "_" + group + "_totalFracErrorMatrix_scaled", metadata);
                }
            }

            for (const auto &[quantity, groupToMatrixMap] : errorMatrixMapUnscaled)
            {
                // Get the quantity in question as a vector
                // ATTN here we flatten the smearing matrix to a column vector
                const auto quantityVector = (quantity == "data" ? dataUnscaled : ubsmear::UBSmearingHelper::Flatten(smearingMatrixUnscaled));

                // Plot the total error matrix (summed over all groups)
                const auto totalErrorMatrix = totalErrorMatrixMapUnscaled.at(quantity);
                PlottingHelper::PlotErrorMatrix(totalErrorMatrix, prefix + "_" + quantity + "_totalErrorMatrix_unscaled", metadata);
                PlottingHelper::PlotFractionalErrorMatrix(totalErrorMatrix, quantityVector, prefix + "_" + quantity + "_totalFracErrorMatrix_unscaled", metadata);

                // Plot the total error matrix for each group individually
                for (const auto &[group, errorMatrix] : groupToMatrixMap)
                {
                    PlottingHelper::PlotErrorMatrix(errorMatrix, prefix + "_" + quantity + "_" + group + "_totalErrorMatrix_unscaled", metadata);
                    PlottingHelper::PlotFractionalErrorMatrix(errorMatrix, quantityVector, prefix + "_" + quantity + "_" + group + "_totalFracErrorMatrix_unscaled", metadata);
                }
            }

            // -----------------------------------------------------------------------------------------------------------------------------
            // Forward-fold the prediction and plot it
            // -----------------------------------------------------------------------------------------------------------------------------
            // for now we just re-use the error matrix function (and use the default palette)
            PlottingHelper::PlotErrorMatrix(smearingMatrixScaled, prefix + "_smearingMatrix_scaled", metadata, true, false);
            PlottingHelper::PlotErrorMatrix(smearingMatrixUnscaled, prefix + "_smearingMatrix_unscaled", metadata, true, false);

            // Now get the predicted cross-section and it's error matrix
            const auto predictionScaled = getMatrix("prediction_scaled");
            const auto predictionBiasVectorScaled = getMatrix("prediction_stat_bias_scaled");
            const auto predictionCovarianceMatrixScaled = getMatrix("prediction_stat_covariance_scaled");
            const auto predictionErrorMatrixScaled = CrossSectionHelper::GetErrorMatrix(predictionBiasVectorScaled, predictionCovarianceMatrixScaled);

            const auto predictionUnscaled = getMatrix("prediction_unscaled");
            const auto predictionBiasVectorUnscaled = getMatrix("prediction_stat_bias_unscaled");
            const auto predictionCovarianceMatrixUnscaled = getMatrix("prediction_stat_covariance_unscaled");
            const auto predictionErrorMatrixUnscaled = CrossSectionHelper::GetErrorMatrix(predictionBiasVectorUnscaled, predictionCovarianceMatrixUnscaled);

            const auto predictionNuWro = getMatrixTrue("prediction");

            // Plot the error matrix on the prediction
            PlottingHelper::PlotErrorMatrix(predictionErrorMatrixScaled, prefix + "_prediction_stat_totalErrorMatrix_scaled", metadata);
            PlottingHelper::PlotFractionalErrorMatrix(predictionErrorMatrixScaled, predictionScaled, prefix + "_prediction_stat_totalFracErrorMatrix_scaled", metadata);
            
            PlottingHelper::PlotErrorMatrix(predictionErrorMatrixUnscaled, prefix + "_prediction_stat_totalErrorMatrix_unscaled", metadata);
            PlottingHelper::PlotFractionalErrorMatrix(predictionErrorMatrixUnscaled, predictionUnscaled, prefix + "_prediction_stat_totalFracErrorMatrix_unscaled", metadata);

            // Now smear the prediction so it can be compared to the data
            std::cout << "Forward folding prediction" << std::endl;
            const auto &[smearedPredictionScaled, smearedPredictionErrorMatrixScaled] = ubsmear::ForwardFold(
                metadata,                                                  // The metadata that defines the binning
                predictionScaled, predictionErrorMatrixScaled,             // The prediction and it's error matrix
                smearingMatrixScaled, totalErrorMatrixMapScaled.at("smearingMatrix"),  // The smearing matrix and it's error matrix
                config.makeXSecPlots.nUniverses,                           // The number of universes to use when propagating the uncertainties
                config.makeXSecPlots.precision);                           // The precision to use when finding eigenvalues and eigenvectors

            // Now smear the prediction so it can be compared to the data
            std::cout << "Forward folding prediction" << std::endl;
            const auto &[smearedPredictionUnscaled, smearedPredictionErrorMatrixUnscaled] = ubsmear::ForwardFold(
                metadata,                                                  // The metadata that defines the binning
                predictionUnscaled, predictionErrorMatrixUnscaled,                         // The prediction and it's error matrix
                smearingMatrixUnscaled, totalErrorMatrixMapUnscaled.at("smearingMatrix"),  // The smearing matrix and it's error matrix
                config.makeXSecPlots.nUniverses,                           // The number of universes to use when propagating the uncertainties
                config.makeXSecPlots.precision);                           // The precision to use when finding eigenvalues and eigenvectors

            std::cout << "Forward folding NuWro prediction" << std::endl;
                const auto &[smearedPredictionNuWro, smearedPredictionErrorMatrixNuWro] = ubsmear::ForwardFold(
                metadata,                                                       // The metadata that defines the binning
                predictionNuWro, predictionErrorMatrixUnscaled,                 // The prediction and it's error matrix
                smearingMatrixNuwro, totalErrorMatrixMapUnscaled.at("smearingMatrix"),  // The smearing matrix and it's error matrix
                config.makeXSecPlots.nUniverses,                                // The number of universes to use when propagating the uncertainties
                config.makeXSecPlots.precision);                                // The precision to use when finding eigenvalues and eigenvectors

            // Plot the error matrix on the smeared prediction
            PlottingHelper::PlotErrorMatrix(smearedPredictionErrorMatrixScaled, prefix + "_smearedPrediction_totalErrorMatrix_scaled", metadata);
            PlottingHelper::PlotFractionalErrorMatrix(smearedPredictionErrorMatrixScaled, smearedPredictionScaled, prefix + "_smearedPrediction_totalFracErrorMatrix_scaled", metadata);

            PlottingHelper::PlotErrorMatrix(smearedPredictionErrorMatrixUnscaled, prefix + "_smearedPrediction_totalErrorMatrix_unscaled", metadata);
            PlottingHelper::PlotFractionalErrorMatrix(smearedPredictionErrorMatrixUnscaled, smearedPredictionUnscaled, prefix + "_smearedPrediction_totalFracErrorMatrix_unscaled", metadata);

            // Save the forward-folded prediction
            FormattingHelper::SaveMatrix(smearedPredictionNuWro, "xsecNuWro_" + selectionName + "_" + xsecName + "_forwardFoldedPrediction.txt");
            // FormattingHelper::SaveMatrix(smearedPredictionErrorMatrix, "xsecNuWro_" + selectionName + "_" + xsecName + "_forwardFoldedPrediction_error.txt");

            // -----------------------------------------------------------------------------------------------------------------------------
            // Get the chi2 for the data-prediction comparison and add it to the table
            // -----------------------------------------------------------------------------------------------------------------------------
            const auto totalDataErrorMatrixScaled = totalErrorMatrixMapScaled.at("data");
            std::cout << "Getting chi2" << std::endl;
            const auto &[chi2Scaled, degreesOfFreedomScaled] = ubsmear::GetChi2(smearedPredictionScaled, smearedPredictionErrorMatrixScaled, dataScaled, totalDataErrorMatrixScaled, config.makeXSecPlots.precision);

            std::cout << "Getting p-value" << std::endl;
            const auto pValueScaled = ubsmear::UBStatisticsHelper::GetPValue(chi2Scaled, degreesOfFreedomScaled, config.makeXSecPlots.precision);

            tableScaled.AddEmptyRow();
            tableScaled.SetEntry("Selection", selectionName);
            tableScaled.SetEntry("Cross-section", xsecName);
            tableScaled.SetEntry("nBins", dataScaled.GetRows());
            tableScaled.SetEntry("Chi2", chi2Scaled);
            tableScaled.SetEntry("DoF", degreesOfFreedomScaled);
            tableScaled.SetEntry("Chi2/DoF", chi2Scaled / static_cast<float>(degreesOfFreedomScaled));
            tableScaled.SetEntry("p-value", pValueScaled);

            const auto totalDataErrorMatrixUnscaled = totalErrorMatrixMapUnscaled.at("data");
            std::cout << "Getting chi2" << std::endl;
            const auto &[chi2Unscaled, degreesOfFreedomUnscaled] = ubsmear::GetChi2(smearedPredictionUnscaled, smearedPredictionErrorMatrixUnscaled, dataUnscaled, totalDataErrorMatrixUnscaled, config.makeXSecPlots.precision);

            std::cout << "Getting p-value" << std::endl;
            const auto pValueUnscaled = ubsmear::UBStatisticsHelper::GetPValue(chi2Unscaled, degreesOfFreedomUnscaled, config.makeXSecPlots.precision);

            tableUnscaled.AddEmptyRow();
            tableUnscaled.SetEntry("Selection", selectionName);
            tableUnscaled.SetEntry("Cross-section", xsecName);
            tableUnscaled.SetEntry("nBins", dataUnscaled.GetRows());
            tableUnscaled.SetEntry("Chi2", chi2Unscaled);
            tableUnscaled.SetEntry("DoF", degreesOfFreedomUnscaled);
            tableUnscaled.SetEntry("Chi2/DoF", chi2Unscaled / static_cast<float>(degreesOfFreedomUnscaled));
            tableUnscaled.SetEntry("p-value", pValueUnscaled);

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
            auto pDataHistScaled = std::make_shared<TH1F>((prefix + "_data_scaled").c_str(), "", binEdges.size() - 1, binEdges.data());
            auto pDataStatOnlyHistScaled = std::make_shared<TH1F>((prefix + "_dataStatOnly_scaled").c_str(), "", binEdges.size() - 1, binEdges.data());
            auto pPredictionHistScaled = std::make_shared<TH1F>((prefix + "_prediction_scaled").c_str(), "", binEdges.size() - 1, binEdges.data());

            auto pDataHistUnscaled = std::make_shared<TH1F>((prefix + "_data_unscaled").c_str(), "", binEdges.size() - 1, binEdges.data());
            auto pDataStatOnlyHistUnscaled = std::make_shared<TH1F>((prefix + "_dataStatOnly_unscaled").c_str(), "", binEdges.size() - 1, binEdges.data());
            auto pPredictionHistUnscaled = std::make_shared<TH1F>((prefix + "_prediction_unscaled").c_str(), "", binEdges.size() - 1, binEdges.data());

            auto pNuWroPredictionHist = std::make_shared<TH1F>((prefix + "_predictionNuWro").c_str(), "", binEdges.size() - 1, binEdges.data());

            // Fill the bins
            // ATTN here we only show the diagonals of the error matrices
            float minY = +std::numeric_limits<float>::max();
            float maxY = -std::numeric_limits<float>::max();
            const auto dataStatErrorMatrixScaled = errorMatrixMapScaled.at("data").at("stat");
            const auto dataStatErrorMatrixUnscaled = errorMatrixMapUnscaled.at("data").at("stat");
            for (unsigned int iBin = 1; iBin <= binEdges.size() - 1; ++iBin)
            {
                // Set the values for the data histogram
                const auto dataValueScaled = dataScaled.At(iBin - 1, 0);
                const auto dataErrorScaled = std::pow(totalDataErrorMatrixScaled.At(iBin - 1, iBin - 1), 0.5f);
                const auto dataStatOnlyErrorScaled = std::pow(dataStatErrorMatrixScaled.At(iBin - 1, iBin - 1), 0.5f);

                const auto dataValueUnscaled = dataUnscaled.At(iBin - 1, 0);
                const auto dataErrorUnscaled = std::pow(totalDataErrorMatrixUnscaled.At(iBin - 1, iBin - 1), 0.5f);
                const auto dataStatOnlyErrorUnscaled = std::pow(dataStatErrorMatrixUnscaled.At(iBin - 1, iBin - 1), 0.5f);

                pDataHistScaled->SetBinContent(iBin, dataValueScaled);
                pDataHistScaled->SetBinError(iBin, dataErrorScaled);

                pDataHistUnscaled->SetBinContent(iBin, dataValueUnscaled);
                pDataHistUnscaled->SetBinError(iBin, dataErrorUnscaled);

                pDataStatOnlyHistScaled->SetBinContent(iBin, dataValueScaled);
                pDataStatOnlyHistScaled->SetBinError(iBin, dataStatOnlyErrorScaled);

                pDataStatOnlyHistUnscaled->SetBinContent(iBin, dataValueUnscaled);
                pDataStatOnlyHistUnscaled->SetBinError(iBin, dataStatOnlyErrorUnscaled);

                // Set the values of the prediction
                const auto predictionValueScaled = smearedPredictionScaled.At(iBin - 1, 0);
                const auto predictionErrorScaled = std::pow(smearedPredictionErrorMatrixScaled.At(iBin - 1, iBin - 1), 0.5f);

                const auto predictionValueUnscaled = smearedPredictionUnscaled.At(iBin - 1, 0);
                const auto predictionErrorUnscaled = std::pow(smearedPredictionErrorMatrixUnscaled.At(iBin - 1, iBin - 1), 0.5f);

                const auto predictionValueNuWro = smearedPredictionNuWro.At(iBin - 1, 0); //dataTrue.At(iBin - 1, 0);

                pPredictionHistScaled->SetBinContent(iBin, predictionValueScaled);
                pPredictionHistScaled->SetBinError(iBin, predictionErrorScaled);

                pPredictionHistUnscaled->SetBinContent(iBin, predictionValueUnscaled);
                pPredictionHistUnscaled->SetBinError(iBin, predictionErrorUnscaled);

                pNuWroPredictionHist->SetBinContent(iBin, predictionValueNuWro);
                pNuWroPredictionHist->SetBinError(iBin, 0);

                // For the proton multiplicity plot, use explicit bin labels
                if (xsecName == "nProtons")
                {
                    for (auto &pHist : {pDataHistScaled, pDataStatOnlyHistScaled, pPredictionHistScaled, pDataHistUnscaled, pDataStatOnlyHistUnscaled, pPredictionHistUnscaled, pNuWroPredictionHist})
                    {
                        pHist->GetXaxis()->SetBinLabel(1, "0");
                        pHist->GetXaxis()->SetBinLabel(2, "1");
                        pHist->GetXaxis()->SetBinLabel(3, ">1");
                    }
                }

                // Get the limiting values
                minY = std::min(minY, dataValueScaled - dataErrorScaled);
                minY = std::min(minY, predictionValueScaled - predictionErrorScaled);
                minY = std::min(minY, dataValueUnscaled - dataErrorUnscaled);
                minY = std::min(minY, predictionValueUnscaled - predictionErrorUnscaled);
                minY = std::min(minY, predictionValueNuWro);


                maxY = std::max(maxY, dataValueScaled + dataErrorScaled);
                maxY = std::max(maxY, predictionValueScaled + predictionErrorScaled);
                maxY = std::max(maxY, dataValueUnscaled + dataErrorUnscaled);
                maxY = std::max(maxY, predictionValueUnscaled + predictionErrorUnscaled);
                maxY = std::max(maxY, predictionValueNuWro);
            }

            // Set the y-range
            const auto padding = (maxY - minY) * 0.05;
            maxY += padding;
            minY -= padding;
            minY = std::max(minY, 0.f);
            minY = 0.f; // Remove this line to get a dynamic lower y-range
            pDataHistScaled->GetYaxis()->SetRangeUser(minY, maxY);
            pDataStatOnlyHistScaled->GetYaxis()->SetRangeUser(minY, maxY);
            pPredictionHistScaled->GetYaxis()->SetRangeUser(minY, maxY);
            pDataHistUnscaled->GetYaxis()->SetRangeUser(minY, maxY);
            pDataStatOnlyHistUnscaled->GetYaxis()->SetRangeUser(minY, maxY);
            pPredictionHistUnscaled->GetYaxis()->SetRangeUser(minY, maxY);
            pNuWroPredictionHist->GetYaxis()->SetRangeUser(minY, maxY);

            // Set the colours of the histograms
            PlottingHelper::SetLineStyle(pDataHistScaled, PlottingHelper::Primary);
            PlottingHelper::SetLineStyle(pDataStatOnlyHistScaled, PlottingHelper::Primary);
            PlottingHelper::SetLineStyle(pPredictionHistScaled, PlottingHelper::Secondary);
            PlottingHelper::SetLineStyle(pDataHistUnscaled, PlottingHelper::Quaternary);
            PlottingHelper::SetLineStyle(pDataStatOnlyHistUnscaled, PlottingHelper::Quaternary);
            PlottingHelper::SetLineStyle(pPredictionHistUnscaled, PlottingHelper::Quinary);

            PlottingHelper::SetLineStyle(pNuWroPredictionHist, PlottingHelper::Tertiary);
            pNuWroPredictionHist->SetLineStyle(4); // Dashed line

            // Make the plot!
            auto pCanvas = PlottingHelper::GetCanvas();
            gStyle->SetEndErrorSize(4);

            // Draw the smeared prediction
            pPredictionHistScaled->Draw("hist");

            // Draw the prediction uncertainties as a semi-transparent band
            auto pHistCloneScaled = static_cast<TH1F *>(pPredictionHistScaled->Clone());
            const auto colScaled = pHistCloneScaled->GetLineColor();
            pHistCloneScaled->SetFillStyle(1001);
            pHistCloneScaled->SetLineColorAlpha(colScaled, 0.f);
            pHistCloneScaled->SetFillColorAlpha(colScaled, 0.3f);
            pHistCloneScaled->Draw("e2 same");

            // Draw the smeared prediction
            pPredictionHistUnscaled->Draw("hist same");

            // Draw the prediction uncertainties as a semi-transparent band
            auto pHistCloneUnscaled = static_cast<TH1F *>(pPredictionHistUnscaled->Clone());
            const auto colUnscaled = pHistCloneUnscaled->GetLineColor();
            pHistCloneUnscaled->SetFillStyle(1001);
            pHistCloneUnscaled->SetLineColorAlpha(colUnscaled, 0.f);
            pHistCloneUnscaled->SetFillColorAlpha(colUnscaled, 0.3f);
            pHistCloneUnscaled->Draw("e2 same");

            // Draw the data as points with error bars
            pDataStatOnlyHistScaled->Draw("e1 same");
            pDataHistScaled->Draw("e1 same");
            pDataStatOnlyHistUnscaled->Draw("e1 same");
            pDataHistUnscaled->Draw("e1 same");
            pNuWroPredictionHist->Draw("e1 same");

            PlottingHelper::SaveCanvas(pCanvas, prefix + "_data-vs-smearedPrediction");
        }
    }

    // Save the table
    tableScaled.WriteToFile("NuWroXSecPlots_goodnessOfFitStatistics_scaled.md");
    tableUnscaled.WriteToFile("NuWroXSecPlots_goodnessOfFitStatistics_unscaled.md");
}

} // namespace ubcc1pi_macros
