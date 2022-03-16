/**
 *  @file  ubcc1pi_standalone/Macros/MakeSidebandParameterPlots.cxx
 *
 *  @brief The implementation file of the MakeSidebandParameterPlots macro
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

void MakeSidebandParameterPlots(const Config &config)
{
    std::cout<<"MakeSidebandParameterPlots - Point 0"<<std::endl;
    // -------------------------------------------------------------------------------------------------------------------------------------
    // Setup an object that holds the details of the systematic parameters to apply
    // -------------------------------------------------------------------------------------------------------------------------------------
    // Here we read in the "dimensions" of the systematic parameters to apply. For multisim parameters (flux & xsec), this is a map from
    // each parameter name to the expected number of universes. For unisim parameters (detector variations), this is a map from the
    // identidier (i.e. name) of the detector variation sample, to the identifier of the corresponding central-value sample. Additionally
    // we set the number of bootstrap universes (for the MC stat uncertainty), the corresponding weights are generated for each event.
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
    std::cout<<"MakeSidebandParameterPlots - Point 1"<<std::endl;

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
    std::cout<<"MakeSidebandParameterPlots - Point 2"<<std::endl;
    // // -------------------------------------------------------------------------------------------------------------------------------------
    // // Get the sideband weights
    // // -------------------------------------------------------------------------------------------------------------------------------------
    // // Loop over all cross-section objects
    typedef std::pair<std::vector<Double_t>,std::vector<Double_t>> paramAndErrorPair; // Todo: Improve code!
    //Parameters: selectionName name
    // std::map<std::string,std::map<std::string, std::vector<float>>> cc0piCovarianceMap; 
    std::map<std::string, std::map<std::string, paramAndErrorPair>> cc0piNominalConstraintMap;
    //Parameters: selectionName name paramName (i.e. golden muonMomentum hadronProduction)
    std::map<std::string, std::map<std::string, std::map<std::string, std::vector<paramAndErrorPair>>>> cc0piUniverseConstraintMap;

    std::cout<<"\n\nDone with CC0pi constraint."<<std::endl;

	// std::ifstream ifs1("cc0piCovarianceMap.bin", std::ios::binary);
    std::ifstream ifs2("cc0piNominalConstraintMap.bin", std::ios::binary);
    std::ifstream ifs3("cc0piUniverseConstraintMap.bin", std::ios::binary);
	
    // boost::archive::binary_iarchive iarch1(ifs1);
    boost::archive::binary_iarchive iarch2(ifs2);
    boost::archive::binary_iarchive iarch3(ifs3);

    // iarch1 >> cc0piCovarianceMap;
    iarch2 >> cc0piNominalConstraintMap;
    iarch3 >> cc0piUniverseConstraintMap;

	// ifs1.close();
    ifs2.close();
    ifs3.close();



    // -------------------------------------------------------------------------------------------------------------------------------------
    for (const auto &[selectionName, enabledMap] : config.extractXSecs.crossSectionIsEnabled)
    {
        if(selectionName=="golden") continue;

        for (const auto &[xsecName, isEnabled] : enabledMap)
        {
            // Skip cross-sections that aren't enabled
            if (!isEnabled) continue;

            const auto weightDimensions = {std::make_pair("xsec", systParams.xsecDimensions), std::make_pair("reint", systParams.reintDimensions), std::make_pair("flux", systParams.fluxDimensions)};
            for (const auto &[group, dimensions] : weightDimensions)
            {
                for (const auto &[paramName, nUniverses] : dimensions)
                {

                    // Define a prefix for the names of the plots
                    const std::string prefix = "sidebandFitPlots_" + selectionName + "_" + xsecName + "_" + group + "_" + paramName;
                    const auto &metadata = metadataMap.at(xsecName);

                    auto cc0piNominalConstraintParam = cc0piNominalConstraintMap.at(selectionName).at(xsecName).first;
                    auto cc0piNominalConstraintParamError = cc0piNominalConstraintMap.at(selectionName).at(xsecName).second;
                    const vector<float> cc0piNominalConstraintParamFloat(cc0piNominalConstraintParam.begin(), cc0piNominalConstraintParam.end());
                    const vector<float> cc0piNominalConstraintParamErrorFloat(cc0piNominalConstraintParamError.begin(), cc0piNominalConstraintParamError.end());

                    std::cout<<"&*%$ --- "<<xsecName<<" "<<group<<" "<<paramName<<" "<<nUniverses<<" ---"<<std::endl;
                    // std::cout<<"&*%$ Nominal "<<": ";
                    // for(const auto& param: cc0piNominalConstraintParam) {
                    //     std::cout<<param<<" ";
                    // }
                    // std::cout<<std::endl;

                    // Calculate mean parameter value across all universes
                    vector<float> meanParameters(metadata.GetNBins(), 0.f);
                    unsigned int totalUniverses=0; 
                    // for (const auto &[paramName, nUniverses] : systParams.xsecDimensions)
                    // {
                    if(cc0piUniverseConstraintMap.at(selectionName).at(xsecName).at(paramName).at(0).first.at(0)>=0)
                    {
                        totalUniverses += nUniverses;
                        for (unsigned int iBin = 0; iBin < metadata.GetNBins(); ++iBin)
                        {
                            std::cout<<"&*%$ Bin "<<iBin<<" - Nom.: "<<cc0piNominalConstraintParam.at(iBin)<<" - Uni.: ";
                            for (unsigned int iUni = 0; iUni < nUniverses; ++iUni)
                            {
                                std::cout<<cc0piUniverseConstraintMap.at(selectionName).at(xsecName).at(paramName).at(iUni).first.at(iBin)<<" ";
                                meanParameters.at(iBin) += cc0piUniverseConstraintMap.at(selectionName).at(xsecName).at(paramName).at(iUni).first.at(iBin);
                            }
                            std::cout<<std::endl;
                        }
                    }
                    // }

                    for (unsigned int iBin = 0; iBin < metadata.GetNBins(); ++iBin)
                    {
                        meanParameters.at(iBin) /= totalUniverses;
                    }

                    // Calculate standard deviation across all universes
                    vector<float> stddevParameters(metadata.GetNBins(), 0.f);
                    // for (const auto &[paramName, nUniverses] : systParams.xsecDimensions)
                    // {
                    if(cc0piUniverseConstraintMap.at(selectionName).at(xsecName).at(paramName).at(0).first.at(0)>=0)
                    {
                        for (unsigned int iBin = 0; iBin < metadata.GetNBins(); ++iBin)
                        {
                            for (unsigned int iUni = 0; iUni < nUniverses; ++iUni)
                            {
                                stddevParameters.at(iBin) += std::pow(cc0piUniverseConstraintMap.at(selectionName).at(xsecName).at(paramName).at(iUni).first.at(iBin)-meanParameters.at(iBin),2);
                                // std::cout<<"stddev: "<<stddevParameters.at(iBin)<<" - "<<cc0piUniverseConstraintMap.at(selectionName).at(xsecName).at(paramName).at(iUni).first.at(iBin)<<" - "<<meanParameters.at(iBin)<<std::endl;
                            }
                        }
                    }
                    // }
                    for (unsigned int iBin = 0; iBin < metadata.GetNBins(); ++iBin)
                    {
                        stddevParameters.at(iBin) = std::pow(stddevParameters.at(iBin)/totalUniverses, 0.5);
                    }
                    
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
                        const auto dataValue = cc0piNominalConstraintParamFloat.at(iBin - 1);
                        const auto dataError = cc0piNominalConstraintParamErrorFloat.at(iBin - 1);

                        pDataHist->SetBinContent(iBin, dataValue);
                        pDataHist->SetBinError(iBin, dataError);

                        // Set the values of the prediction
                        const auto predictionValue = meanParameters.at(iBin-1);
                        const auto predictionError = stddevParameters.at(iBin-1);

                        pPredictionHist->SetBinContent(iBin, predictionValue);
                        pPredictionHist->SetBinError(iBin, predictionError);

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
                    auto padding = (maxY - minY) * 0.2;
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
                    pPredictionHist->Draw("e1");

                    // // Draw the prediction uncertainties as a semi-transparent band
                    // auto pHistClone = static_cast<TH1F *>(pPredictionHist->Clone());
                    // auto col = pHistClone->GetLineColor();
                    // pHistClone->SetFillStyle(1001);
                    // pHistClone->SetLineColorAlpha(col, 0.f);
                    // pHistClone->SetFillColorAlpha(col, 0.3f);
                    // pHistClone->Draw("e2 same");

                    // Draw the data as points with error bars
                    // pDataStatOnlyHist->Draw("e1 same");
                    pDataHist->Draw("e1 same");

                    PlottingHelper::SaveCanvas(pCanvas1, prefix + "_SidebandParameters");
                }
            }
        }
    }

}

} // namespace ubcc1pi_macros

