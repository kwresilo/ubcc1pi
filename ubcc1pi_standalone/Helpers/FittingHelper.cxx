/**
 *  @file  ubcc1pi_standalone/Helpers/FittingHelper.cxx
 *
 *  @brief The implementation file of the fitting helper class
 */

#include "ubcc1pi_standalone/Helpers/FittingHelper.h"
#include <TMinuit.h>


namespace ubcc1pi
{

FittingHelper::FittingHelper(const Int_t binNumber) :
    nBins(binNumber){}

void FittingHelper::Fit(void(*fcn)(Int_t &, Double_t *, Double_t &f, Double_t *, Int_t), std::pair<std::vector<Double_t>, std::vector<Double_t>> &result, bool &successful, std::vector<float> &covMatrix, const int &printlevel, const std::vector<Double_t> &initialParameters)
{
    // if(smearingMatrix.IsSquare() && truthRecoMatrix.GetRows() == data.GetRows() && data.GetRows()==dataStatUncertainty.GetRows() && data.GetColumns()==dataStatUncertainty.GetColumns())
    //     throw std::logic_error("FittingHelper::Fit - Incompatible input dimenstions.");
    // std::cout<<"FittingHelper::Fit ((0))"<<std::endl;
    successful = false;
    // auto nBins = x.GetRows();
    TMinuit minuit(nBins);
    // std::cout<<"FittingHelper::Fit ((1))"<<std::endl;
    minuit.SetPrintLevel(printlevel);
    // minuit.SetMaxIterations(5000);
    minuit.SetFCN(fcn);
    // std::cout<<"FittingHelper::Fit ((2))"<<std::endl;

    Double_t arglist[10];
    Int_t ierflg = 0;

    arglist[0]=1;
    minuit.mnexcm("SET ERR", arglist, 1, ierflg);

    if(initialParameters.empty())
    {
        std::cout<<"FittingHelper::Fit ((1))"<<std::endl;
        // Set starting values and step sizes for parameters
        for (unsigned int iBin=0; iBin<nBins; iBin++)
        {
            minuit.mnparm(iBin, std::to_string(iBin), 1.0, 0.1, -9999.0, 9999.0, ierflg); // random limit: 9999.0 mnparm // parameter number here follows normals convention starting from 0 ...
        }
    }
    else
    {
        std::cout<<"FittingHelper::Fit ((2))"<<std::endl;
        if(initialParameters.size() != nBins)
        {
            std::cout<<"FittingHelper::Fit - initialParameters.size != nBins"<<std::endl;
            throw std::logic_error("FittingHelper::Fit - initialParameters.size != nBins");
        }
        // Set starting values and step sizes for parameters
        for (unsigned int iBin=0; iBin<nBins; iBin++)
        {
            minuit.mnparm(iBin, std::to_string(iBin), initialParameters.at(iBin), 0.05, -9999.0, 9999.0, ierflg); // random limit: 9999.0 mnparm // parameter number here follows normals convention starting from 0 ...
        }
    }
    std::cout<<"FittingHelper::Fit ((3))"<<std::endl;                                                                          // ... specifically: "Parameter number as referenced by user in FCN"

    std::cout<<"FittingHelper::Fit ((4))"<<std::endl;
    std::vector<bool> fixedParameters(nBins, false);
    while (!successful){
        arglist[0]=5000;
        arglist[1]=0.1;
        std::cout<<"FittingHelper::Fit ((4.1))"<<std::endl;
        minuit.mnexcm("MIGRAD", arglist, 2, ierflg);
        std::cout<<"FittingHelper::Fit ((4.2))"<<std::endl;
        if(ierflg!=0)
        {
            std::cout<<"FittingHelper::Fit - Did not converge. Execution failed with ierflg: "<<ierflg<<std::endl;
            return;
        }
        successful = true;
        for (unsigned int iBin=0; iBin<nBins; iBin++)
        {
            std::cout<<"FittingHelper::Fit ((4.3))"<<std::endl;
            Double_t param, paramError;
            minuit.GetParameter(iBin, param, paramError);
            std::cout<<"FittingHelper::Fit ((4.4))"<<std::endl;
            // Double_t paramDebug, paramErrorDebug; // todo remove
            // minuit.GetParameter(nBins, paramDebug, paramErrorDebug); // todo remove
            // std::cout<<"Debug - paramDebug: "<<paramDebug<<" paramErrorDebug: "<<paramErrorDebug<<std::endl; // todo remove
            if(param<0)
            {
                std::cout<<"FittingHelper::Fit - Parameter "<<iBin<<" has negative value: "<<param<<std::endl;
                Int_t ierflgSet;
                fixedParameters.at(iBin) = true;
                successful = false;
                arglist[0] = iBin+1; // Numeration starts from 1 instead of 0
                arglist[1] = 0.0;
                std::cout<<"FittingHelper::Fit ((5))"<<std::endl;
                minuit.mnexcm("SET PARAM", arglist, 2, ierflgSet);
                if(ierflgSet!=0)
                {
                    std::cout<<"FittingHelper::Fit - minuit SET PARAM failed."<<std::endl;
                    throw std::logic_error("FittingHelper::Fit - minuit SET PARAM failed.");
                }
                std::cout<<"FittingHelper::Fit ((6))"<<std::endl;
                minuit.FixParameter(iBin);
            }
            std::cout<<"FittingHelper::Fit ((7))"<<std::endl;
        }
    }
    std::cout<<"FittingHelper::Fit ((10))"<<std::endl;
    for (unsigned int iBin=0; iBin<nBins; iBin++)
    {
       if(fixedParameters.at(iBin))
       {
            Int_t ierflgScan;
            minuit.Release(iBin);
            arglist[0] = iBin+1; // Numeration starts from 1 instead of 0
            minuit.mnexcm("SCAN", arglist, 1, ierflgScan);
            if(ierflgScan!=0)
            {
                std::cout<<"FittingHelper::Fit - minuit SCAN failed."<<std::endl;
                throw std::logic_error("FittingHelper::Fit - minuit SCAN failed.");
            }
       }
    }
    std::cout<<"FittingHelper::Fit ((11))"<<std::endl;

    std::vector<Double_t>paramVector, paramErrorVector;
    for (unsigned int iBin=0; iBin<nBins; iBin++)
    {
        Double_t param, paramError;
        minuit.GetParameter(iBin, param, paramError);
        paramVector.push_back(param);
        paramErrorVector.push_back(paramError);
    }

    std::cout<<"FittingHelper::Fit ((12))"<<std::endl;

    // FMIN: the best function value found so far
    // FEDM: the estimated vertical distance remaining to minimum
    // ERRDEF: the value of UP defining parameter uncertainties
    // NPARI: the number of currently variable parameters
    // NPARX: the highest (external) parameter number defined by user
    // ISTAT: a status integer indicating how good is the covariance matrix:
    //     0= not calculated at all
    //     1= approximation only, not accurate
    //     2= full matrix, but forced positive-definite
    //     3= full accurate covariance matrix

    Double_t amin,edm,errdef;
    Int_t nvpar,nparx,icstat;
    minuit.mnstat(amin,edm,errdef,nvpar,nparx,icstat);
    // minuit.mnprin(3,amin);

    // covMatrix = std::vector<Double_t>(nBins*nBins);
    Double_t covMatrixFit[nBins][nBins];
    minuit.mnemat(&covMatrixFit[0][0], nBins);
    minuit.mnmatu(1);

    covMatrix.clear();
    for (unsigned int i=0; i<nBins; i++)
    {
        for (unsigned int j=0; j<nBins; j++)
        {
            covMatrix.push_back((float)covMatrixFit[i][j]);
        }
    }

    result = std::make_pair(paramVector, paramErrorVector);
    std::cout<<"FittingHelper::Fit ((13))"<<std::endl;
}


// //______________________________________________________________________________
// void FittingHelper::fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
// {
//     //calculate chi2
//     Double_t chisq = 0;
//     Double_t delta;
//     auto nBins = x.GetRows();

//     auto xScaled = x;
//     for (Int_t i=0; i<nBins; i++)
//     {
//         xScaled.SetElement(i,0,xScaled.At(i,0)*par[i]);
//     }
//     const auto xSmeared = S*xScaled;

//     for (Int_t i=0; i<nBins; i++)
//     {
//         delta  = (y.At(i,0)-xSmeared.At(i,0))/errorY.At(i,0);
//         chisq += delta*delta;
//     }
//     f = chisq;
// }

}
