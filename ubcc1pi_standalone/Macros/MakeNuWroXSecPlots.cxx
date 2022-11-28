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

    const std::map<std::string, CrossSectionHelper::SystDimensionsMap> systDimensionMapBNB(
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
        });

    const std::map<std::string, CrossSectionHelper::SystDimensionsMap> systDimensionMapFakeData(
        {
            {"xsec", config.extractXSecs.xsecDimensions},
            {"misc", {
                {"bootstrap", config.extractXSecs.nBootstrapUniverses},
                {"sidebandWeights", config.extractXSecs.nBootstrapUniverses},
            }}
        });

    std::vector<std::string> dataTypeList;
    if(config.global.useBNBAsData) dataTypeList.push_back("NuWro");
    if(config.global.useNuWroAsData) dataTypeList.push_back("BNB");
    if(config.global.useGenieAsData) dataTypeList.push_back("Genie");

    for(const auto dataTypeName: dataTypeList)
    {
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
                std::cout<<"DEBug Point 0"<<std::endl;
                // Skip cross-sections that aren't enabled
                if (!isEnabled)
                    continue;

                std::cout << selectionName << " - " << xsecName << std::endl;

                // Define a lambda function to get a matrix from a file for brevity
                const auto getMatrix = [&, selectionName = selectionName, xsecName = xsecName](const std::string &identifier) -> ubsmear::UBMatrix {
                    return ubsmear::UBFileHelper::ReadMatrix("xsec_" + selectionName + "_" + xsecName + "_" + identifier + ".txt");
                };

                const auto getMatrixTrue = [&, selectionName = selectionName, xsecName = xsecName](const std::string &identifier) -> ubsmear::UBMatrix {
                    return ubsmear::UBFileHelper::ReadMatrix("xsec_" + selectionName + "_" + xsecName + "_" + identifier + "_NuWroTruth.txt");
                };

                // Define a lambda function to get a matrix from a file an trim any overflow / underflow bins
                const auto &metadata = metadataMap.at(xsecName);
                const auto getTrimmedMatrix = [&] (const std::string &identifier) -> ubsmear::UBMatrix {
                    return ubsmear::UBSmearingHelper::TrimUnderOverflowBins(getMatrix(identifier), metadata);
                };

                // Define a lambda function to get a matrix from a file an trim any over- / underflow bins
                // const auto &metadataTrue = metadataMap.at(xsecName);
                const auto getTrimmedMatrixTrue = [&] (const std::string &identifier) -> ubsmear::UBMatrix {
                    return ubsmear::UBSmearingHelper::TrimUnderOverflowBins(getMatrixTrue(identifier), metadata);
                };

                // Get the data cross-section
                const auto dataScaled = getTrimmedMatrix("data_"+dataTypeName+"Scaled");
                const auto dataUnscaled = getTrimmedMatrix("data_"+dataTypeName+"Unscaled");
                // const auto dataTrue = getTrimmedMatrixTrue("data"); //getTrimmedMatrixTrue("data");

                // Get the smearing matrix
                const auto smearingMatrixScaled = getMatrix("smearingMatrix_"+dataTypeName+"Scaled");
                const auto smearingMatrixUnscaled = getMatrix("smearingMatrix_"+dataTypeName+"Unscaled");
                // const auto smearingMatrixNuwro = getMatrixTrue("smearingMatrix");

                // Build the total error matrix for each group of systematic parameters
                // Here errorMatrixMap is indexed by [quantity][group], where quantity = "data" or "smearingMatrix"
                // Here totalErrorMatrixMap is index by [quantity] and contains the sum (over all groups) of the entries in errorMatrixMap
                std::map<std::string, std::map<std::string, ubsmear::UBMatrix> > errorMatrixMapScaled;
                std::map<std::string, ubsmear::UBMatrix> totalErrorMatrixMapScaled;
                std::map<std::string, std::map<std::string, ubsmear::UBMatrix> > errorMatrixMapUnscaled;
                std::map<std::string, ubsmear::UBMatrix> totalErrorMatrixMapUnscaled;

                std::cout<<"DEBug Point 1"<<std::endl;
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
                    const auto statUncertaintiesScaled = (quantity == "data") ? getMatrixFunction("data_stat_"+dataTypeName+"Scaled") : ubsmear::UBMatrixHelper::GetZeroMatrix(nBins, 1);
                    const auto statUncertaintiesUnscaled = (quantity == "data") ? getMatrixFunction("data_stat_"+dataTypeName+"Unscaled") : ubsmear::UBMatrixHelper::GetZeroMatrix(nBins, 1);

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

                    const auto systDimensionMap = (dataTypeName == "BNB") ? systDimensionMapBNB : systDimensionMapFakeData;

                    // Handle the multisim parameters
                    for (const auto &[group, dimensions] : systDimensionMap)
                    {
                        // Setup an empty error matrix for this group
                        auto errorMatrixTotalScaled =  ubsmear::UBMatrixHelper::GetZeroMatrix(nBins, nBins);
                        auto errorMatrixTotalUnscaled =  ubsmear::UBMatrixHelper::GetZeroMatrix(nBins, nBins);

                        // Loop over the parameters in this group
                        for (const auto &[paramName, nUniverses] : dimensions)
                        {
                            // Get the bias vector and covariance matrix
                            const auto biasVectorScaled = getMatrixFunction(quantity + "_" + group + "_" + paramName + "_bias_"+dataTypeName+"Scaled");
                            const auto covarianceMatrixScaled = getMatrixFunction(quantity + "_" + group + "_" + paramName + "_covariance_"+dataTypeName+"Scaled");

                            const auto biasVectorUnscaled = getMatrixFunction(quantity + "_" + group + "_" + paramName + "_bias_"+dataTypeName+"Unscaled");
                            const auto covarianceMatrixUnscaled = getMatrixFunction(quantity + "_" + group + "_" + paramName + "_covariance_"+dataTypeName+"Unscaled");

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

                    std::cout<<"DEBug Point 2"<<std::endl;
                    // Handle the unisim parameters
                    if(dataTypeName == "BNB") // Not needed for fake-data (NuWro/Genie)
                    {
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
                                const auto biasVectorScaled = getMatrixFunction(quantity + "_" + group + "_" + paramName + "_bias_"+dataTypeName+"Scaled");
                                const auto covarianceMatrixScaled = getMatrixFunction(quantity + "_" + group + "_" + paramName + "_covariance_"+dataTypeName+"Scaled");

                                const auto biasVectorUnscaled = getMatrixFunction(quantity + "_" + group + "_" + paramName + "_bias_"+dataTypeName+"Unscaled");
                                const auto covarianceMatrixUnscaled = getMatrixFunction(quantity + "_" + group + "_" + paramName + "_covariance_"+dataTypeName+"Unscaled");

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
                    }
                    // Add the grand total to the map
                    totalErrorMatrixMapScaled.emplace(quantity, errorMatrixTotalSumScaled);
                    totalErrorMatrixMapUnscaled.emplace(quantity, errorMatrixTotalSumUnscaled);
                }

                // Define a prefix for the names of the plots
                const std::string prefix = "xsecPlots_" + selectionName + "_" + xsecName;

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
                    PlottingHelper::PlotErrorMatrix(totalErrorMatrix, prefix + "_" + quantity + "_totalErrorMatrix_"+dataTypeName+"Scaled", metadata);
                    PlottingHelper::PlotFractionalErrorMatrix(totalErrorMatrix, quantityVector, prefix + "_" + quantity + "_totalFracErrorMatrix_"+dataTypeName+"Scaled", metadata);

                    // Plot the total error matrix for each group individually
                    for (const auto &[group, errorMatrix] : groupToMatrixMap)
                    {
                        PlottingHelper::PlotErrorMatrix(errorMatrix, prefix + "_" + quantity + "_" + group + "_totalErrorMatrix_"+dataTypeName+"Scaled", metadata);
                        PlottingHelper::PlotFractionalErrorMatrix(errorMatrix, quantityVector, prefix + "_" + quantity + "_" + group + "_totalFracErrorMatrix_"+dataTypeName+"Scaled", metadata);
                    }
                }

                for (const auto &[quantity, groupToMatrixMap] : errorMatrixMapUnscaled)
                {
                    // Get the quantity in question as a vector
                    // ATTN here we flatten the smearing matrix to a column vector
                    const auto quantityVector = (quantity == "data" ? dataUnscaled : ubsmear::UBSmearingHelper::Flatten(smearingMatrixUnscaled));

                    // Plot the total error matrix (summed over all groups)
                    const auto totalErrorMatrix = totalErrorMatrixMapUnscaled.at(quantity);
                    PlottingHelper::PlotErrorMatrix(totalErrorMatrix, prefix + "_" + quantity + "_totalErrorMatrix_"+dataTypeName+"Unscaled", metadata);
                    PlottingHelper::PlotFractionalErrorMatrix(totalErrorMatrix, quantityVector, prefix + "_" + quantity + "_totalFracErrorMatrix_"+dataTypeName+"Unscaled", metadata);

                    // Plot the total error matrix for each group individually
                    for (const auto &[group, errorMatrix] : groupToMatrixMap)
                    {
                        PlottingHelper::PlotErrorMatrix(errorMatrix, prefix + "_" + quantity + "_" + group + "_totalErrorMatrix_"+dataTypeName+"Unscaled", metadata);
                        PlottingHelper::PlotFractionalErrorMatrix(errorMatrix, quantityVector, prefix + "_" + quantity + "_" + group + "_totalFracErrorMatrix_"+dataTypeName+"Unscaled", metadata);
                    }
                }

                // -----------------------------------------------------------------------------------------------------------------------------
                // Forward-fold the prediction and plot it
                // -----------------------------------------------------------------------------------------------------------------------------
                // for now we just re-use the error matrix function (and use the default palette)
                PlottingHelper::PlotErrorMatrix(smearingMatrixScaled, prefix + "_smearingMatrix_"+dataTypeName+"Scaled", metadata, true, false);
                PlottingHelper::PlotErrorMatrix(smearingMatrixUnscaled, prefix + "_smearingMatrix_"+dataTypeName+"Unscaled", metadata, true, false);

                // Now get the predicted cross-section and it's error matrix
                const auto predictionScaled = getMatrix("prediction_"+dataTypeName+"Scaled");
                const auto predictionBiasVectorScaled = getMatrix("prediction_stat_bias_"+dataTypeName+"Scaled");
                const auto predictionCovarianceMatrixScaled = getMatrix("prediction_stat_covariance_"+dataTypeName+"Scaled");
                const auto predictionErrorMatrixScaled = CrossSectionHelper::GetErrorMatrix(predictionBiasVectorScaled, predictionCovarianceMatrixScaled);

                const auto predictionUnscaled = getMatrix("prediction_"+dataTypeName+"Unscaled");
                const auto predictionBiasVectorUnscaled = getMatrix("prediction_stat_bias_"+dataTypeName+"Unscaled");
                const auto predictionCovarianceMatrixUnscaled = getMatrix("prediction_stat_covariance_"+dataTypeName+"Unscaled");
                const auto predictionErrorMatrixUnscaled = CrossSectionHelper::GetErrorMatrix(predictionBiasVectorUnscaled, predictionCovarianceMatrixUnscaled);

                // Plot the error matrix on the prediction
                PlottingHelper::PlotErrorMatrix(predictionErrorMatrixScaled, prefix + "_prediction_stat_totalErrorMatrix_"+dataTypeName+"Scaled", metadata);
                PlottingHelper::PlotFractionalErrorMatrix(predictionErrorMatrixScaled, predictionScaled, prefix + "_prediction_stat_totalFracErrorMatrix_"+dataTypeName+"Scaled", metadata);
                
                PlottingHelper::PlotErrorMatrix(predictionErrorMatrixUnscaled, prefix + "_prediction_stat_totalErrorMatrix_"+dataTypeName+"Unscaled", metadata);
                PlottingHelper::PlotFractionalErrorMatrix(predictionErrorMatrixUnscaled, predictionUnscaled, prefix + "_prediction_stat_totalFracErrorMatrix_"+dataTypeName+"Unscaled", metadata);

                // Now smear the prediction so it can be compared to the data
                std::cout << "Forward folding scaled prediction" << std::endl;
                const auto &[smearedPredictionScaled, smearedPredictionErrorMatrixScaled] = ubsmear::ForwardFold(
                    metadata,                                                  // The metadata that defines the binning
                    predictionScaled, predictionErrorMatrixScaled,             // The prediction and it's error matrix
                    smearingMatrixScaled, totalErrorMatrixMapScaled.at("smearingMatrix"),  // The smearing matrix and it's error matrix
                    config.makeXSecPlots.nUniverses,                           // The number of universes to use when propagating the uncertainties
                    config.makeXSecPlots.precision);                           // The precision to use when finding eigenvalues and eigenvectors

                // Now smear the prediction so it can be compared to the data
                std::cout << "Forward folding unscaled prediction" << std::endl;
                const auto &[smearedPredictionUnscaled, smearedPredictionErrorMatrixUnscaled] = ubsmear::ForwardFold(
                    metadata,                                                  // The metadata that defines the binning
                    predictionUnscaled, predictionErrorMatrixUnscaled,                         // The prediction and it's error matrix
                    smearingMatrixUnscaled, totalErrorMatrixMapUnscaled.at("smearingMatrix"),  // The smearing matrix and it's error matrix
                    config.makeXSecPlots.nUniverses,                           // The number of universes to use when propagating the uncertainties
                    config.makeXSecPlots.precision);                           // The precision to use when finding eigenvalues and eigenvectors

                // Only used for NuWro fake data study - use simple smear instead of forward fold as we don't need smearedPredictionErrorMatrix
                const auto zeroMatrix = ubsmear::UBMatrixHelper::GetZeroMatrix(1, 1);
                const auto smearingMatrixNuWroTruth = dataTypeName == "NuWro" ? getMatrix("smearingMatrix_NuWroTruth") : zeroMatrix; // todo: improve code
                const auto predictionNuWroTruth = dataTypeName == "NuWro" ? getMatrixTrue("prediction") : zeroMatrix;
                std::cout << "Forward folding NuWro prediction" << std::endl;
                const auto smearedPredictionNuWroTruthUntrimmed = dataTypeName == "NuWro" ? ubsmear::UBSmearingHelper::Smear(metadata, predictionNuWroTruth, smearingMatrixNuWroTruth): zeroMatrix;
                const auto smearedPredictionNuWroTruth = dataTypeName == "NuWro" ? ubsmear::UBSmearingHelper::TrimUnderOverflowBins(smearedPredictionNuWroTruthUntrimmed, metadata) : zeroMatrix;
                FormattingHelper::SaveMatrix(smearedPredictionNuWroTruth, "xsec_" + selectionName + "_" + xsecName + "_forwardFoldedPrediction_NuWroTruth.txt");

                // const auto &[smearedPredictionNuWroTruth, smearedPredictionErrorMatrixNuWro] = ubsmear::ForwardFold(
                //     metadata,                                                       // The metadata that defines the binning
                //     predictionNuWroTruth, predictionErrorMatrixUnscaled,                 // The prediction and it's error matrix
                //     smearingMatrixNuwroTruth, totalErrorMatrixMapUnscaled.at("smearingMatrix"),  // The smearing matrix and it's error matrix
                //     config.makeXSecPlots.nUniverses,                                // The number of universes to use when propagating the uncertainties
                //     config.makeXSecPlots.precision);                                // The precision to use when finding eigenvalues and eigenvectors

                // Save the forward-folded prediction

                // Plot the error matrix on the smeared prediction
                PlottingHelper::PlotErrorMatrix(smearedPredictionErrorMatrixScaled, prefix + "_smearedPrediction_totalErrorMatrix_"+dataTypeName+"Scaled", metadata);
                PlottingHelper::PlotFractionalErrorMatrix(smearedPredictionErrorMatrixScaled, smearedPredictionScaled, prefix + "_smearedPrediction_totalFracErrorMatrix_"+dataTypeName+"Scaled", metadata);

                PlottingHelper::PlotErrorMatrix(smearedPredictionErrorMatrixUnscaled, prefix + "_smearedPrediction_totalErrorMatrix_"+dataTypeName+"Unscaled", metadata);
                PlottingHelper::PlotFractionalErrorMatrix(smearedPredictionErrorMatrixUnscaled, smearedPredictionUnscaled, prefix + "_smearedPrediction_totalFracErrorMatrix_"+dataTypeName+"Unscaled", metadata);

                std::cout<<"DEBug Point 3"<<std::endl;
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
                auto pDataHistScaled = std::make_shared<TH1F>((prefix + "_data_"+dataTypeName+"Scaled").c_str(), "", binEdges.size() - 1, binEdges.data());
                auto pDataStatOnlyHistScaled = std::make_shared<TH1F>((prefix + "_dataStatOnly_"+dataTypeName+"Scaled").c_str(), "", binEdges.size() - 1, binEdges.data());
                auto pPredictionHistScaled = std::make_shared<TH1F>((prefix + "_prediction_"+dataTypeName+"Scaled").c_str(), "", binEdges.size() - 1, binEdges.data());

                auto pDataHistUnscaled = std::make_shared<TH1F>((prefix + "_data_"+dataTypeName+"Unscaled").c_str(), "", binEdges.size() - 1, binEdges.data());
                auto pDataStatOnlyHistUnscaled = std::make_shared<TH1F>((prefix + "_dataStatOnly_"+dataTypeName+"Unscaled").c_str(), "", binEdges.size() - 1, binEdges.data());
                auto pPredictionHistUnscaled = std::make_shared<TH1F>((prefix + "_prediction_"+dataTypeName+"Unscaled").c_str(), "", binEdges.size() - 1, binEdges.data());

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

                    pPredictionHistScaled->SetBinContent(iBin, predictionValueScaled);
                    pPredictionHistScaled->SetBinError(iBin, predictionErrorScaled);

                    pPredictionHistUnscaled->SetBinContent(iBin, predictionValueUnscaled);
                    pPredictionHistUnscaled->SetBinError(iBin, predictionErrorUnscaled);

                    float predictionValueNuWroTruth = 0.f; // placeholder value - stodo clean up;
                    if(dataTypeName == "NuWro")
                    {
                        predictionValueNuWroTruth = smearedPredictionNuWroTruth.At(iBin - 1, 0); //dataTrue.At(iBin - 1, 0);
                        pNuWroPredictionHist->SetBinContent(iBin, predictionValueNuWroTruth);
                        pNuWroPredictionHist->SetBinError(iBin, 0);
                    }

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
                    if(dataTypeName == "NuWro") minY = std::min(minY, predictionValueNuWroTruth);


                    maxY = std::max(maxY, dataValueScaled + dataErrorScaled);
                    maxY = std::max(maxY, predictionValueScaled + predictionErrorScaled);
                    maxY = std::max(maxY, dataValueUnscaled + dataErrorUnscaled);
                    maxY = std::max(maxY, predictionValueUnscaled + predictionErrorUnscaled);
                    if(dataTypeName == "NuWro") maxY = std::max(maxY, predictionValueNuWroTruth);
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
                PlottingHelper::SetLineStyle(pDataHistUnscaled, 2);//PlottingHelper::Tertiary);
                PlottingHelper::SetLineStyle(pDataStatOnlyHistUnscaled, 2);//PlottingHelper::Tertiary);
                PlottingHelper::SetLineStyle(pPredictionHistUnscaled, PlottingHelper::Quinary);


                PlottingHelper::SetLineStyle(pNuWroPredictionHist, PlottingHelper::Quaternary);
                pNuWroPredictionHist->SetLineStyle(4); // Dashed line
                pDataHistUnscaled->SetLineWidth(3);
                pDataStatOnlyHistUnscaled->SetLineWidth(3);
                pDataHistScaled->SetLineWidth(2);
                pDataStatOnlyHistScaled->SetLineWidth(2);

                // Make the plot!
                auto pCanvas = PlottingHelper::GetCanvas();

                gStyle->SetEndErrorSize(4);

                // Draw the smeared prediction
                pPredictionHistScaled->Draw("hist");

                // Draw the prediction uncertainties as a semi-transparent band
                auto pHistCloneScaled = static_cast<TH1F *>(pPredictionHistScaled->Clone((dataTypeName+xsecName).c_str()));
                PlottingHelper::SetLineStyle(pHistCloneScaled, PlottingHelper::Secondary); // needed because otherwise this end up white for some plots some reason
                const auto colScaled = pHistCloneScaled->GetLineColor();
                pHistCloneScaled->SetFillStyle(1001);
                pHistCloneScaled->SetLineColorAlpha(colScaled, 0.f);
                pHistCloneScaled->SetFillColorAlpha(colScaled, 0.3f);
                pHistCloneScaled->Draw("e2 same");

                // Same values as above - only used to confirm
                // // Draw the smeared prediction
                // pPredictionHistUnscaled->Draw("hist same");

                // // Draw the prediction uncertainties as a semi-transparent band
                // auto pHistCloneUnscaled = static_cast<TH1F *>(pPredictionHistUnscaled->Clone());
                // const auto colUnscaled = pHistCloneUnscaled->GetLineColor();
                // pHistCloneUnscaled->SetFillStyle(1001);
                // pHistCloneUnscaled->SetLineColorAlpha(colUnscaled, 0.f);
                // pHistCloneUnscaled->SetFillColorAlpha(colUnscaled, 0.3f);
                // pHistCloneUnscaled->Draw("e2 same");

                // Draw the data as points with error bars
                if(dataTypeName == "NuWro") pNuWroPredictionHist->Draw("e1 same");
                pDataStatOnlyHistUnscaled->Draw("e1 same");
                pDataHistUnscaled->Draw("e1 same");
                pDataStatOnlyHistScaled->Draw("e1 same");
                pDataHistScaled->Draw("e1 same");

                PlottingHelper::SaveCanvas(pCanvas, prefix + "_data-vs-smearedPrediction_" + dataTypeName);
            }
        }

        // Save the table
        tableScaled.WriteToFile("XSecPlots_goodnessOfFitStatistics_"+dataTypeName+"Scaled.md");
        tableUnscaled.WriteToFile("XSecPlots_goodnessOfFitStatistics_"+dataTypeName+"Unscaled.md");
    }
}

} // namespace ubcc1pi_macros
