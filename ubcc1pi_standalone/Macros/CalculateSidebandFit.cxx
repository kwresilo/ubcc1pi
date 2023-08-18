/**
 *  @file  ubcc1pi_standalone/Macros/CalculateSidebandFit.cxx
 *
 *  @brief The implementation file of the CalculateSidebandFit macro
 */

#include "ubcc1pi_standalone/Macros/Macros.h"
#include "ubcc1pi_standalone/Objects/FileReader.h"

#include "ubcc1pi_standalone/Helpers/AnalysisHelper.h"
#include "ubcc1pi_standalone/Helpers/NormalisationHelper.h"
#include "ubcc1pi_standalone/Helpers/SelectionHelper.h"
#include "ubcc1pi_standalone/Helpers/CrossSectionHelper.h"
#include "ubcc1pi_standalone/Helpers/FormattingHelper.h"
#include "ubcc1pi_standalone/Helpers/FittingHelper.h"
#include "ubsmear.h"

#include <fstream> // Todo: not use txt files

// Boost libraries
// #include "binary_iarchive.hpp"
#include "binary_oarchive.hpp"
#include "binary_object.hpp"
#include "map.hpp"
#include "vector.hpp"

using namespace ubcc1pi;

std::vector<float> x, y, errorY, S;

//______________________________________________________________________________

void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
    //calculate chi2
    if(x.size()!=y.size() || y.size()!=errorY.size() || x.size()*x.size()!=S.size())
        throw std::logic_error("Fitting function for CalculateSidebandFit - Incompatible input dimenstions.");

    Double_t chisq = 0;
    Double_t delta;
    Int_t nBins = x.size();

    auto xScaled = x;
    for (Int_t i=0; i<nBins; i++)
    {
        xScaled[i]*=par[i];
    }

    std::vector<float> xSmeared(nBins, 0);

    for (Int_t i=0; i<nBins; i++)
    {
        for (Int_t j=0; j<nBins; j++)
        {
            xSmeared[i] +=  S[i+j*nBins]*xScaled[j];
        }
    }
    // const auto xSmeared = S*xScaled;

    for (Int_t i=0; i<nBins; i++)
    {
        delta  = (y[i]-xSmeared[i])/errorY[i];
        chisq += delta*delta;
    }
    f = chisq;
}

namespace ubcc1pi_macros
{

void CalculateSidebandFit(const Config &config)
{
    std::map<std::string, std::map<std::string, CrossSectionHelper::CrossSection> > xsecMapSideband;
	std::ifstream ifs1("xsecMapSideband.bin", std::ios::binary);
    boost::archive::binary_iarchive iarch1(ifs1);
    iarch1 >> xsecMapSideband;
	ifs1.close();
    std::cout<<"Sucess Loading xsecMapSideband"<<std::endl;



    // -------------------------------------------------------------------------------------------------------------------------------------
    // Calculate the sideband weights
    // -------------------------------------------------------------------------------------------------------------------------------------
    // Loop over all cross-section objects
    typedef std::pair<std::vector<Double_t>,std::vector<Double_t>> paramAndErrorPair; // Todo: Improve code!
    //Parameters: selectionName name
    std::map<std::string,std::map<std::string, std::vector<float>>> cc0piCovarianceMap;
    std::map<std::string, std::map<std::string, paramAndErrorPair>> cc0piNominalConstraintMap;
    //Parameters: selectionName name paramName (i.e. golden muonMomentum hadronProduction)
    std::map<std::string, std::map<std::string, std::map<std::string, std::vector<paramAndErrorPair>>>> cc0piUniverseConstraintMap;

    try
    {
        for (const auto &[selectionName, xsecs] : xsecMapSideband)
        {
            for (const auto &[name, xsec] : xsecs)
            {
                // -------------------------------------------------------------------------------------------------------------------------------------
                // Fit nominal
                // -------------------------------------------------------------------------------------------------------------------------------------

                std::cout<<"_______________________Fitting: "<<selectionName<<" - "<<name<<"_______________________"<<std::endl;
                // const auto sidebandCVFit = xsec.GetPredictedCrossSection(scalingData);

                std::cout<<"_______________________Fitting Point 0.1"<<std::endl;
                const auto selectedEventsData = xsec.GetSelectedBNBDataEvents();
                std::cout<<"_______________________Fitting Point 0.2"<<std::endl;
                // Get the smearing matrix of selected events
                const auto smearingMatrixAllSelected = xsec.GetSmearingMatrixAllSelected();
                std::cout<<"_______________________Fitting Point 0.3"<<std::endl;
                const auto selectedEventsBackgroundReco = xsec.GetSelectedBackgroundEvents();
                std::cout<<"_______________________Fitting Point 0.4"<<std::endl;
                const auto selectedEventsSignalTruth = xsec.GetSelectedSignalEvents();
                std::cout<<"_______________________Fitting Point 0.5"<<std::endl;
                auto signalData = selectedEventsData - selectedEventsBackgroundReco;


                for(unsigned int r = 0; r<signalData.GetRows(); r++) // DEBUG - TODO: REMOVE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                {
                    for(unsigned int c = 0; c<signalData.GetColumns(); c++) // DEBUG - TODO: REMOVE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    {
                        if(signalData.At(r, c) < 0)
                        {
                            std::cout<<"signalData.At("<<r<<","<<c<<") = "<<signalData.At(r, c)<<std::endl;
                        }
                        signalData.SetElement(r, c, std::max(signalData.At(r, c), 0.f));
                    }
                }

                std::cout<<"_______________________Fitting Point 1"<<std::endl;
                // const auto &[pPredictionStatBias, pPredictionStatCovariance] = xsec.GetPredictedCrossSectionStatUncertainty(scalingData);
                // const auto predictionErrorMatrix = CrossSectionHelper::GetErrorMatrix(*pPredictionBiasVector, *pPredictionCovarianceMatrix);\

                std::vector<float> elements;
                const auto nBins = signalData.GetRows();
                std::cout<<"_______________________Fitting Point 1.1"<<std::endl;
                for (unsigned int iBin = 0; iBin < nBins; ++iBin)
                {
                    std::cout<<"_______________________Fitting Point 1.2"<<std::endl;
                    const auto value = signalData.At(iBin, 0);
                    if (value<0)
                    {
                        std::cout<<"ERROR: ExtractXSec - Background-removed signal data is negative."<<std::endl;
                        // throw std::logic_error("ERROR: ExtractXSec - Background-removed signal data is negative."); // TODO: Uncomment
                    }
                    elements.push_back(AnalysisHelper::GetCountUncertainty(std::max(value,0.f)));
                    // elements.push_back(AnalysisHelper::GetCountUncertainty(value));
                    std::cout<<"_______________________Fitting Point 1.3"<<std::endl;
                }
                std::cout<<"_______________________Fitting Point 2"<<std::endl;
                // const ubsmear::UBMatrix signalDataUncertainty(elements, nBins, 1);

                std::cout<<"_______________________Fitting Point 3"<<std::endl;
                x = selectedEventsSignalTruth.GetValues();
                y = signalData.GetValues();
                errorY = elements;//signalDataUncertainty.GetValues();
                S = smearingMatrixAllSelected.GetValues();

                std::cout<<"\nx (selectedEventsSignalTruth): \n";
                for (const auto &xValue : x)
                    std::cout<<xValue<<" ";

                std::cout<<"\ny (signalData): \n";
                for (const auto &yValue : y)
                    std::cout<<yValue<<" ";

                std::cout<<"\nerrorY (signalDataUncertainty): \n";
                for (const auto &errorYValue : errorY)
                    std::cout<<errorYValue<<" ";

                std::cout<<"\nSmearing matrix: \n";
                for(unsigned int i = 0; i<S.size(); i++)
                {
                    if(i%nBins==0)
                        std::cout<<"\n";
                    std::cout<<S[i]<<" ";
                }


                std::cout<<"_______________________Fitting Point 4"<<std::endl;
                // auto minimizer = FittingHelper(selectedEventsSignal, signalData, signalDataUncertainty, smearingMatrix);
                auto minimizer = FittingHelper(nBins);
                std::pair<std::vector<Double_t>, std::vector<Double_t>> result;

                std::vector<float> fitCovMatrixVector;
                minimizer.Fit(fcn, result, fitCovMatrixVector, 0);
                cc0piCovarianceMap[selectionName].emplace(name, fitCovMatrixVector);

                std::cout<<"\nFitting covariance matrix: \n";
                for(unsigned int i = 0; i<fitCovMatrixVector.size(); i++)
                {
                    if(i%nBins==0)
                        std::cout<<"\n";
                    std::cout<<fitCovMatrixVector[i]<<" ";
                }
                std::cout<<"_______________________Fitting Point 4.1"<<std::endl;
                // vector<float> fitCovMatrixFloat(fitCovMatrix.begin(), fitCovMatrix.end()); //Todo avoid this
                const ubsmear::UBMatrix sidebandCovMatrix(fitCovMatrixVector, nBins, nBins);
                std::cout<<"_______________________Fitting Point 4.2"<<std::endl;
                FormattingHelper::SaveMatrix(sidebandCovMatrix, "CC0Pi_" + selectionName + "_" + name + "_sideband_stat_covariance.txt");
                std::cout<<"_______________________Fitting Point 4.3"<<std::endl;
                FormattingHelper::SaveMatrix(smearingMatrixAllSelected, "CC0Pi_" + selectionName + "_" + name + "_smearingMatrixAllSelected.txt");
                std::cout<<"_______________________Fitting Point 4.4"<<std::endl;
                vector<float> paramVector(result.first.begin(), result.first.end()); //Todo avoid this
                vector<float> paramErrorVector(result.second.begin(), result.second.end()); //Todo avoid this
                std::cout<<"_______________________Fitting Point 4.5"<<std::endl;
                std::cout<<"nBins :"<<nBins<<std::endl;
                std::cout<<"\nerrorY (paramErrorVector Double_t): \n";
                for (const auto &e : result.second)
                    std::cout<<e<<" ";
                std::cout<<"\nerrorY (paramErrorVector float): \n";
                for (const auto &e : paramErrorVector)
                    std::cout<<e<<" ";
                std::cout<<std::endl;
                const ubsmear::UBMatrix sidebandParamVectorTruth(paramVector, nBins, 1);
                const ubsmear::UBMatrix sidebandErrorVectorTruth(paramErrorVector, nBins, 1);
                std::cout<<"_______________________Fitting Point 4.6"<<std::endl;
                // const auto sidebandErrorVectorReco = smearingMatrixAllSelected*sidebandErrorVectorTruth; // Todo: check this multiplication is correctly computed
                FormattingHelper::SaveMatrix(sidebandParamVectorTruth, "CC0Pi_" + selectionName + "_" + name + "_sideband_parameterVector.txt");
                FormattingHelper::SaveMatrix(sidebandErrorVectorTruth, "CC0Pi_" + selectionName + "_" + name + "_sideband_parameterErrorVector.txt");
                // const ubsmear::UBMatrix sidebandCovMatrix(fitCovMatrixVector, nBins, nBins);
                // FormattingHelper::SaveMatrix(sidebandCovMatrix, "CC0Pi_" + selectionName + "_" + name + "_sideband_stat_errorVector.txt");

                std::cout<<"_______________________Fitting Point 5"<<std::endl;

                std::cout<<"param: \n";
                for (const auto &p : result.first)
                    std::cout<<p<<" ";

                std::cout<<"paramError: \n";
                for (const auto &p : result.second)
                    std::cout<<p<<" ";


                cc0piNominalConstraintMap[selectionName].emplace(name, result);
                // FittingHelper::Fit(selectedEventsSignal, signalData, signalDataUncertainty, smearingMatrix);

                // const auto sidebandWeights = xsec.GetSidebandWeights(scalingData);

                for (const auto &[paramName, nUniverses] : systParams.xsecDimensions)
                {
                    // -------------------------------------------------------------------------------------------------------------------------------------
                    // Fit each universe
                    // -------------------------------------------------------------------------------------------------------------------------------------

                    // const auto nUniverses = config.extractXSecs.nBootstrapUniverses;

                    const auto selectedSignalTruthUniverses = xsec.GetSelectedSignalRecoTruthMap().at("xsec").at(paramName);
                    const auto selectedBackgroundRecoUniverses = xsec.GetSelectedBackgroundRecoMap().at("xsec").at(paramName);

                    std::vector<paramAndErrorPair> resultVector;
                    for (unsigned int iUni = 0; iUni < nUniverses; ++iUni)
                    {
                        // AnalysisHelper::PrintLoadingBar(iUni, nUniverses);
                        const auto selectedSignalTruth = xsec.GetSignalSelectedTrue(selectedSignalTruthUniverses.at(iUni));
                        const auto selectedBackgoundReco = CrossSectionHelper::GetMatrixFromHist(selectedBackgroundRecoUniverses.at(iUni));
                        const auto signalData = selectedEventsData - selectedBackgoundReco;

                        std::vector<float> elements;
                        const auto nBins = signalData.GetRows();
                        for (unsigned int iBin = 0; iBin < nBins; ++iBin)
                        {
                            const auto value = signalData.At(iBin, 0);
                            if (value<0)
                                std::cout<<"Value below zero: "<<value<<std::endl;

                            elements.push_back(AnalysisHelper::GetCountUncertainty(std::max(value,0.f)));
                            // elements.push_back(AnalysisHelper::GetCountUncertainty(value));
                        }

                        x = selectedSignalTruth.GetValues();
                        y = signalData.GetValues();
                        errorY = elements;
                        std::pair<std::vector<Double_t>, std::vector<Double_t>> result;
                        std::vector<float> covMatrixInUniverse;
                        minimizer.Fit(fcn, result, covMatrixInUniverse, 0);

                        // std::cout<<"\nParameter("<<iUni<<"):";
                        // for (const auto &r : result.first)
                        //     std::cout<<" "<<r;

                        // std::cout<<"\nUncertainty:";
                        // for (const auto &r : result.second)
                        //     std::cout<<" "<<r;

                        resultVector.push_back(result);
                    }
                    cc0piUniverseConstraintMap[selectionName][name].emplace(paramName, resultVector);
                }
            }
        }
    }
    catch(exception &e)
    {
        cout << "CC0pi - Caught exception Point 1: "<<e.what();
        return;
    }

	std::ofstream ofs1("cc0piCovarianceMap.bin", std::ios::binary);
    std::ofstream ofs2("cc0piNominalConstraintMap.bin", std::ios::binary);
    std::ofstream ofs3("cc0piUniverseConstraintMap.bin", std::ios::binary);


    boost::archive::binary_oarchive oarch1(ofs1);
    boost::archive::binary_oarchive oarch2(ofs2);
    boost::archive::binary_oarchive oarch3(ofs3);

    oarch1 << cc0piCovarianceMap;
    oarch2 << cc0piNominalConstraintMap;
    oarch3 << cc0piUniverseConstraintMap;

	ofs1.close();
    ofs2.close();
    ofs3.close();

    std::cout<<"------------- All Done -------------"<<std::endl;

    return;
}

} // namespace ubcc1pi_macros
