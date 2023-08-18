/**
 *  @file  ubcc1pi_standalone/Macros/MakeSidebandTemplateFit.cxx
 *
 *  @brief The implementation file of the ExtractXSecs macro
 */

#include "ubcc1pi_standalone/Macros/Macros.h"
#include "ubcc1pi_standalone/Objects/FileReader.h"

#include "ubcc1pi_standalone/Helpers/AnalysisHelper.h"
// #include "ubcc1pi_standalone/Helpers/NormalisationHelper.h"
// #include "ubcc1pi_standalone/Helpers/SelectionHelper.h"
// #include "ubcc1pi_standalone/Helpers/CrossSectionHelper.h"
// #include "ubcc1pi_standalone/Helpers/FormattingHelper.h"

#include <fstream> // Todo: not use txt files
#include <TMinuit.h>

std::vector<Double_t> x, y, errorY;
Int_t nBins;

/// https://root.cern/doc/v608/Ifit_8C_source.html
//______________________________________________________________________________
Double_t func(Double_t x, Double_t par)
{
    return x*par;
}

//______________________________________________________________________________
void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
    //calculate chisquare
    Double_t chisq = 0;
    Double_t delta;
    for (Int_t i=0; i<nBins; i++)
    {
        delta  = (y[i]-func(x[i],par[i]))/errorY[i];
        chisq += delta*delta;
    }
    f = chisq;
}

using namespace ubcc1pi;

namespace ubcc1pi_macros
{

void MakeSidebandTemplateFit(const Config &config)
{
        // -------------------------------------------------------------------------------------------------------------------------------------
    // Setup the cross-section objects
    // -------------------------------------------------------------------------------------------------------------------------------------
    // Here we make a map from a name of the cross-section to the cross-section object itself. In this way, we can iterate through the
    // cross-section objects and reduce code-bloat. The first index is an identifier for the selection that's applied (generic or goldlen),
    // the second index is an identifier for the kinematic quantity that's relevant for the cross-section (e.g. muonMomentum), and the
    // mapped type is the cross-section object.
    std::map<std::string, std::map<std::string, CrossSectionHelper::CrossSection> > xsecMap;

    // We additionally make a map from each cross-section to the limits of the phase-space that we should consider. The key is the
    // identifier for the kinematic quantity and the mapped value is a pair containing the limits [min, max]
    std::map< std::string, std::pair<Double_t, Double_t> > phaseSpaceMap;

    // ATTN the configuration allows the user to enable or disable each cross-section. If a cross-section has been disabled, then it won't
    // be added the xSecMap. However, the phaseSpaceMap always includes all kinematic paramters!

    // Add the differential cross-sections
    for (const auto &name : std::vector <std::string> {"muonCosTheta"})//,"muonPhi","muonMomentum", "pionCosTheta", "pionPhi", "pionMomentum", "muonPionAngle", "nProtons"})
    {
        x.clear();
        y.clear();
        errorY.clear();
        std::vector<Double_t> cc0piData;
        std::vector<Double_t> cc0piDataStatUncertainty;
        std::vector<Double_t> cc0piPrediction;
        std::vector<Double_t> cc0piPredictionStatUncertainty;

        std::ifstream signalSelected("CC0pi_" + name + "_signal_selected_eventRate.txt");
        std::ifstream dataSelected("CC0pi_" + name + "_data_selected_eventRate.txt");
        std::ifstream backgroundSelected("CC0pi_" + name + "_background_selected_eventRate.txt");
        std::string signalBinValue;
        std::string dataBinValue;
        std::string backgroundBinValue;
        // Read the next line from File untill it reaches the end.
        while (std::getline(signalSelected, signalBinValue) && std::getline(dataSelected, dataBinValue) && std::getline(backgroundSelected, backgroundBinValue))
        {
            const auto binValueData = std::stof(dataBinValue);
            const auto binValuePrediction = std::stof(signalBinValue) + std::stof(backgroundBinValue);
            cc0piData.push_back(binValueData);
            cc0piDataStatUncertainty.push_back(AnalysisHelper::GetCountUncertainty(binValueData));
            cc0piPrediction.push_back(binValuePrediction);
            cc0piPredictionStatUncertainty.push_back(AnalysisHelper::GetCountUncertainty(binValuePrediction));
        }
        if(std::getline(signalSelected, signalBinValue) || std::getline(dataSelected, dataBinValue) || std::getline(backgroundSelected, backgroundBinValue))
            throw std::logic_error("ExtractXSecs - CC0pi weight files have different numbers of entries.");

        // cc0piConstraInt_tMap.emplace(name, {cc0piData, cc0piDataStatUncertainty, cc0piPrediction, cc0piPredictionStatUncertainty});

        x = cc0piPrediction;
        y = cc0piData;
        errorY = cc0piDataStatUncertainty;
        nBins = x.size();

        std::cout<<"#############################################################\n";
        std::cout<<"#############################################################\n";
        std::cout<<"cc0piPrediction: ";
        for (const auto &xBin: x)
        {
            std::cout<<xBin<<", ";
        }
        std::cout<<"\n";

        std::cout<<"cc0piData: ";
        for (const auto &yBin: y)
        {
            std::cout<<yBin<<", ";
        }
        std::cout<<"\n";

        std::cout<<"cc0picc0piDataStatUncertaintyData: ";
        for (const auto &errorYBin: errorY)
        {
            std::cout<<errorYBin<<", ";
        }
        std::cout<<"\n";

        std::cout<<"#############################################################\n";
        std::cout<<"#############################################################\n";

        TMinuit minuit(nBins);  //initialize TMinuit with a maximum of 5 params
        minuit.SetFCN(fcn);

        //std::vector<Double_t> arglist(10);
        Double_t arglist[10];
        Int_t ierflg = 0;

        //arglist.at(0) = 1;
        arglist[0]=1;
        minuit.mnexcm("SET ERR", arglist ,1, ierflg);

        // Set starting values and step sizes for parameters
        for (Int_t iBin=0; iBin<nBins; iBin++)
        {
            minuit.mnparm(iBin, std::to_string(iBin), 1.f, 0.2, 0.f, 0.f, ierflg);
        }

        // Now ready for minimization step
        // arglist.at(0) = 500;
        // arglist.at(1) = 1;
        arglist[0]=500;
        arglist[1]=0.01;
        minuit.mnexcm("MIGRAD", arglist ,2,ierflg);

        // PrInt_t results
        Double_t amin,edm,errdef;
        Int_t nvpar,nparx,icstat;
        minuit.mnstat(amin,edm,errdef,nvpar,nparx,icstat);
        //minuit.mnprin(3,amin);
    }
}

} // namespace ubcc1pi_macros