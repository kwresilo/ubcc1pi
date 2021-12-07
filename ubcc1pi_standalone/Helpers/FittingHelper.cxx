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

void FittingHelper::Fit(void(*fcn)(Int_t &, Double_t *, Double_t &f, Double_t *, Int_t), std::pair<std::vector<Double_t>, std::vector<Double_t>> &result, std::vector<float> &covMatrix, const int printlevel)
{
    // if(smearingMatrix.IsSquare() && truthRecoMatrix.GetRows() == data.GetRows() && data.GetRows()==dataStatUncertainty.GetRows() && data.GetColumns()==dataStatUncertainty.GetColumns())
    //     throw std::logic_error("FittingHelper::Fit - Incompatible input dimenstions.");

    // auto nBins = x.GetRows();
    TMinuit minuit(nBins);
    minuit.SetPrintLevel(printlevel);
    minuit.SetFCN(fcn);

    Double_t arglist[10];
    Int_t ierflg = 0;

    arglist[0]=1;
    minuit.mnexcm("SET ERR", arglist ,1, ierflg);

    // Set starting values and step sizes for parameters
    for (Int_t iBin=0; iBin<nBins; iBin++)
    {
        minuit.mnparm(iBin, std::to_string(iBin), 1.f, 0.2, 0.f, 0.f, ierflg);
    }

    arglist[0]=500;
    arglist[1]=0.01;
    minuit.mnexcm("MIGRAD", arglist,2,ierflg);
    if(ierflg!=0)
    {
        std::cout<<"FittingHelper::Fit - MIGRAD not execution failed with ierflg: "<<ierflg<<std::endl;
        throw std::logic_error("FittingHelper::Fit - MIGRAD not execution failed with ierflg: "+std::to_string(ierflg));
    }

    
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
    //minuit.mnprin(3,amin);
    
    std::vector<Double_t>paramVector, paramErrorVector;
    for (Int_t iBin=0; iBin<nBins; iBin++)
    {
        Double_t param, paramError;
        minuit.GetParameter(iBin, param, paramError);
        paramVector.push_back(param);
        paramErrorVector.push_back(paramError);
    }
 
    // covMatrix = std::vector<Double_t>(nBins*nBins); 
    Double_t covMatrixFit[nBins][nBins];
    minuit.mnemat(&covMatrixFit[0][0], nBins);
    std::cout<<"@@@@@@@@ Covariance matrix:"<<std::endl;
    minuit.mnmatu(1);

    covMatrix.clear();
    for (int i=0; i<nBins; i++)
    {
        for (int j=0; j<nBins; j++)
        {
            covMatrix.push_back((float)covMatrixFit[i][j]);       
        }
    }
    result = std::make_pair(paramVector, paramErrorVector);
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
