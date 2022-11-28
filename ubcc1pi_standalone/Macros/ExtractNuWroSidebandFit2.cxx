/**
 *  @file  ubcc1pi_standalone/Macros/ExtractNuWroSidebandFit2.cxx
 *
 *  @brief The implementation file of the ExtractNuWroSidebandFit2 macro
 */

#include "ubcc1pi_standalone/Macros/Macros.h"
#include "ubcc1pi_standalone/Objects/FileReader.h"

#include "ubcc1pi_standalone/Helpers/AnalysisHelper.h"
#include "ubcc1pi_standalone/Helpers/NormalisationHelper.h"
#include "ubcc1pi_standalone/Helpers/SelectionHelper.h"
#include "ubcc1pi_standalone/Helpers/CrossSectionHelper.h"
#include "ubcc1pi_standalone/Helpers/FormattingHelper.h"
#include "ubcc1pi_standalone/Helpers/FittingHelper.h"
#include "ubcc1pi_standalone/Helpers/ExtractionHelper.h"
#include "ubsmear.h"
// #include "ubcc1pi_standalone/ubsmear/inc/ubsmear/Helpers/UBSmearingHelper.h"

// Boost libraries
// #include "binary_iarchive.hpp"
#include "binary_oarchive.hpp"
#include "binary_object.hpp"
#include "map.hpp"
#include "vector.hpp"

using namespace ubcc1pi;

const ubsmear::UBMatrix *pPredictionNuWro, *pDataNuWro, *pSmearingMatrixNuWro, *pPredictionErrorMatrixNuWro; //, *pTotalErrorMatrixMapDataNuWro, *pTotalErrorMatrixMapSmearingMatrixNuWro;//, *pPredictionErrorMatrixNuWro;
const std::map<std::string, ubsmear::UBMatrix> *pTotalErrorMatrixMap;
const ubsmear::UBXSecMeta *pMetadataNuWro;

const auto testingFraction = 0.05f;
// void GetChi2NuWro(Int_t &npar, Double_t *gin, Double_t &chi2, Double_t *par, Int_t iflag)
// {
//     std::vector<float> paramVec;
//     const auto nBins = pMetadataNuWro->GetNBins();
//     for (unsigned int i = 0; i<nBins; ++i) paramVec.push_back(par[i]);
//     const ubsmear::UBMatrix parameters(paramVec, paramVec.size(), 1);
    
//     // for(const auto &p: paramVec) std::cout << p << " ";
//     // std::cout<<std::endl;

//     const float precision = std::numeric_limits<float>::epsilon(); // std::numeric_limits<float>::epsilon();  ///< The precision to use when finding eigenvalues & eigenvectors

//     const auto parameterScaledPrediction = ElementWiseOperation(*pPredictionNuWro, parameters, [](const auto &l, const auto& r) { return l * r; });
//     const auto smearedPrediction = *pSmearingMatrixNuWro * parameterScaledPrediction;
//     // std::cout<<"parameterScaledPrediction:"<<std::endl;
//     // parameterScaledPrediction.Print();
//     // std::cout<<"smearedPrediction:"<<std::endl;
//     // smearedPrediction.Print();
//     // const auto smearedPredictionTrimmed = ubsmear::UBSmearingHelper::TrimUnderOverflowBins(smearedPrediction, *pMetadataNuWro);
//     // const auto dataNuWroTrimmed = ubsmear::UBSmearingHelper::TrimUnderOverflowBins(*pDataNuWro, *pMetadataNuWro);
//     // const auto statErrorMatrixTrimmed = ubsmear::UBSmearingHelper::TrimUnderOverflowBins(pTotalErrorMatrixMap->at("data"), *pMetadataNuWro);

//     const auto smearedPredictionTrimmed = smearedPrediction;
//     const auto dataNuWroTrimmed = *pDataNuWro;
//     const auto statErrorMatrixTrimmed = pTotalErrorMatrixMap->at("data");

//     const auto nBinsTrimmed = dataNuWroTrimmed.GetRows();
//     auto zeroMatrix =  ubsmear::UBMatrixHelper::GetZeroMatrix(nBinsTrimmed, nBinsTrimmed);
//     std::tie(chi2, std::ignore) = ubsmear::GetChi2(smearedPredictionTrimmed, zeroMatrix, dataNuWroTrimmed, statErrorMatrixTrimmed, precision);
//     // std::cout<<"GetChi2NuWro: "<<chi2<<std::endl;

//     if(!std::isfinite(chi2))
//     {
//         throw std::logic_error("ERROR: GetChi2NuWro - chi2 value not finite.");
//     }
// }

void GetChi2NuWro(Int_t &npar, Double_t *gin, Double_t &chi2, Double_t *par, Int_t iflag)
{
    // std::cout<<"GetChi2NuWro started"<<std::endl;
    const auto trim = true;
    const auto precision = 1e-5;//std::numeric_limits<float>::epsilon();
    const auto nUniverses = 10000u; // 10000u;

    std::vector<float> paramVec;
    // std::cout<<"GetChi2NuWro - Point 0"<<std::endl;
    const auto nBins = pMetadataNuWro->GetNBins();
    for (unsigned int i = 0; i<nBins; ++i) paramVec.push_back(par[i]);
    const ubsmear::UBMatrix parameters(paramVec, paramVec.size(), 1);
    // std::cout<<"GetChi2NuWro - Point 1"<<std::endl;
    const auto parameterScaledPrediction = ElementWiseOperation(*pPredictionNuWro, parameters, [](const auto &l, const auto& r) { return l * r; });
    std::cout<<"GetChi2NuWro - Point 1.1 - pPredictionNuWro:"<<std::endl;
    pPredictionNuWro->Print();
    std::cout<<"GetChi2NuWro - Point 1.1.1 - parameters:"<<std::endl;
    parameters.Print();
    std::cout<<"GetChi2NuWro - Point 1.2 - parameterScaledPrediction:"<<std::endl;
    parameterScaledPrediction.Print();
    // std::cout<<"GetChi2NuWro - Point 2"<<std::endl;
    // const auto &[smearedPrediction, smearedPredictionErrorMatrix] = ubsmear::ForwardFold(
    //     *pMetadataNuWro,                                                     // The metadata that defines the binning
    //     parameterScaledPrediction, *pPredictionErrorMatrixNuWro,             // The prediction and it's error matrix
    //     *pSmearingMatrixNuWro, pTotalErrorMatrixMap->at("smearingMatrix"),   // The smearing matrix and it's error matrix
    //     nUniverses,                                                          // The number of universes to use when propagating the uncertainties
    //     precision);                                                         // The precision to use when finding eigenvalues and eigenvectors

    const auto smearedPrediction = ubsmear::UBSmearingHelper::Smear(*pMetadataNuWro, parameterScaledPrediction, *pSmearingMatrixNuWro);

    // std::cout<<"GetChi2NuWro - Point 2.1 - smearedPrediction:"<<std::endl;
    // smearedPrediction.Print();
    // std::cout<<"GetChi2NuWro - Point 3"<<std::endl;
    // const auto totalDataErrorMatrixScaled = pTotalErrorMatrixMap->at("data");

    // std::cout<<"GetChi2NuWro - Point 4 - pTotalErrorMatrixMap->at(\"data\"):"<<std::endl;
    // pTotalErrorMatrixMap->at("data").Print();
    // std::cout<<"GetChi2NuWro - Point 4.1 - smearedPredictionErrorMatrix:"<<std::endl;
    // smearedPredictionErrorMatrix.Print();
    // std::cout<<"GetChi2NuWro - Point 4.3 - pDataNuWro:"<<std::endl;
    // pDataNuWro->Print();
    // const auto &[chi2Scaled, degreesOfFreedomScaled] = ubsmear::GetChi2(smearedPrediction, smearedPredictionErrorMatrix, ubsmear::UBSmearingHelper::TrimUnderOverflowBins(*pDataNuWro , *pMetadataNuWro), 
    //                                                                     ubsmear::UBSmearingHelper::TrimUnderOverflowBins(pTotalErrorMatrixMap->at("data") , *pMetadataNuWro), precision);

    // const auto nRows = pTotalErrorMatrixMap->at("data").GetRows();
    // const auto nColumns = pTotalErrorMatrixMap->at("data").GetColumns();
    // auto zeroMatrix =  ubsmear::UBMatrixHelper::GetZeroMatrix(nRows, nColumns);
    // const auto &[chi2Scaled, degreesOfFreedomScaled] = ubsmear::GetChi2(smearedPrediction, zeroMatrix, *pDataNuWro, pTotalErrorMatrixMap->at("data"), precision);

    // const auto &[chi2Scaled, degreesOfFreedomScaled] = ubsmear::GetChi2(smearedPrediction, zeroMatrix, ubsmear::UBSmearingHelper::TrimUnderOverflowBins(*pDataNuWro , *pMetadataNuWro), // for ForwardFold() where over-/underflow bins are trimmed
    //                                                                     ubsmear::UBSmearingHelper::TrimUnderOverflowBins(pTotalErrorMatrixMap->at("data") , *pMetadataNuWro), precision);
    
    Double_t chi2Scaled = 0;
    const auto diff = smearedPrediction - *pDataNuWro; //trimmedDataNuWro;
    std::cout<<"GetChi2NuWro - Point 3.1 - diff:"<<std::endl;
    diff.Print();
    // const auto chi2Scaled = ElementWiseOperation(diff, trimmedDataNuWro, [](const auto &l, const auto& r) { return l * r; });

    const auto trimmedDiff = trim ? ubsmear::UBSmearingHelper::TrimUnderOverflowBins(diff, *pMetadataNuWro) : diff;
    const auto trimmedDataErrorNuWro = trim ? ubsmear::UBSmearingHelper::TrimUnderOverflowBins(pTotalErrorMatrixMap->at("data"), *pMetadataNuWro): pTotalErrorMatrixMap->at("data");
    const auto trimmedParameters = trim ? ubsmear::UBSmearingHelper::TrimUnderOverflowBins(parameters, *pMetadataNuWro) : parameters;

    // const auto trimmedDiff = ubsmear::UBSmearingHelper::TrimUnderOverflowBins(diff, *pMetadataNuWro);
    // const auto trimmedDataErrorNuWro = ubsmear::UBSmearingHelper::TrimUnderOverflowBins(pTotalErrorMatrixMap->at("data"), *pMetadataNuWro);
    // const auto trimmedParameters = ubsmear::UBSmearingHelper::TrimUnderOverflowBins(parameters, *pMetadataNuWro);

    for (unsigned int i=0; i<trimmedDiff.GetRows(); i++)
    {
        chi2Scaled += std::pow(trimmedDiff.At(i, 0), 2)/trimmedDataErrorNuWro.At(i, i);
        // std::cout<<"GetChi2NuWro - chi2Scaled: "<<chi2Scaled<<" - trimmedDiff["<<i<<"]: "<<trimmedDiff.At(i, 0)<<" - trimmedDataErrorNuWro.At("<<i<<", "<<i<<"): "<<trimmedDataErrorNuWro.At(i, i);
        // std::cout<<" - trimmedParameters.At("<<i<<", 0):"<<trimmedParameters.At(i, 0)<<std::endl;
    }

    // for (unsigned int i=0; i<pDataNuWro->GetRows(); i++)
    // {
    //     chi2Scaled += std::pow(diff.At(i, 0), 2)/pTotalErrorMatrixMap->at("data").At(i, i);//trimmedDataErrorNuWro.At(i, i);
    //     std::cout<<"GetChi2NuWro - chi2Scaled: "<<chi2Scaled<<" - diff["<<i<<"]: "<<diff.At(i, 0)<<" - pTotalErrorMatrixMap->at(\"data\").At("<<i<<", "<<i<<"): "<<pTotalErrorMatrixMap->at("data").At(i, i);
    //     std::cout<<" - parameters.At("<<i<<", 0):"<<parameters.At(i, 0)<<std::endl;
    // }

    // std::cout<<"* chi2Scaled: "<<chi2Scaled<<" - chi2ScaledToo: "<<chi2ScaledToo<<std::endl;
    // for(const auto &p: paramVec) std::cout << p << " ";
    // std::cout<<std::endl;
    // std::cout<<"GetChi2NuWro finished"<<std::endl;
    chi2 = chi2Scaled;
}

namespace ubcc1pi_macros
{

void ExtractNuWroSidebandFit2(const Config &config)
{
    // -------------------------------------------------------------------------------------------------------------------------------------
    // Setup an object that holds the details of the systematic parameters to apply
    // -------------------------------------------------------------------------------------------------------------------------------------
    // Here we read in the "dimensions" of the systematic parameters to apply. For multisim parameters (flux & xsec), this is a map from
    // each parameter name to the expected number of universes. For unisim parameters (detector variations), this is a map from the
    // identidier (i.e. name) of the detector variation sample, to the identifier of the corresponding central-value sample. Additionally
    // we set the number of bootstrap universes (for the MC stat uncertainty), the corresponding weights are generated for each event.
    CrossSectionHelper::CrossSection::SystParams systParams;
    systParams.fluxDimensions = config.extractXSecs.fluxDimensions;
    systParams.xsecDimensions = config.extractXSecs.xsecDimensions;
    systParams.reintDimensions = config.extractXSecs.reintDimensions;
    systParams.detVarDimensions = config.extractXSecs.detVarDimensions;
    systParams.potFracUncertainty = config.extractXSecs.potFracUncertainty;
    systParams.nBootstrapUniverses = config.extractXSecs.nBootstrapUniverses;
    systParams.nSidebandFitUniverses = config.extractXSecs.nSidebandFitUniverses;

    CrossSectionHelper::SystDimensionsMap emptySystDimMap = {};
    auto systParamsNuWro = systParams; // shallow copy sufficient since struct contains no pointers
    systParamsNuWro.potFracUncertainty = 0.0;
    systParamsNuWro.fluxDimensions = emptySystDimMap;
    systParamsNuWro.reintDimensions = emptySystDimMap;

    auto systParamsGenie = systParamsNuWro; // shallow copy sufficient since struct contains no pointers
    systParamsGenie.xsecDimensions = emptySystDimMap;

    // -------------------------------------------------------------------------------------------------------------------------------------
    // Setup an object that holds the information about how we should scale an event rate to obtain a cross-section
    // -------------------------------------------------------------------------------------------------------------------------------------
    // Here we specify the:
    // - Flux             [10^-10 cm^-2 POT^-1]
    // - Exposure POT     [10^20 POT]               (stored in config as [POT])
    // - Target density   [10^31 nucleons/cm^3]     (stored in config as [10^23 nucleons/cm^3])
    // - Fiducial volume  [cm^3]
    //
    // Hence the units of the eventual cross-seciton are:
    // - Cross-section    [Flux * Exposure * Target density * Fiducial volume]^-1 = [10^-41 cm^2 / nucleon]
    //
    // Here we use a FluxReweightor to specify the flux. For each universe of each flux systematic paramter, this uses the ratio of the
    // total neutrino event rate (as a function of the true neutrino energy) to the same rate in the nominal simulation to reweight the
    // input neutrino flux distribution. The integrated flux is each universe is used to scale the selected event rate when calculating the
    // cross-section in that universe. For all non-flux parameters, the nominal integrated flux is used.
    CrossSectionHelper::CrossSection::ScalingData scalingData;

    // -------------------------------------------------------------------------------------------------------------------------------------
    // Setup the event selection
    // -------------------------------------------------------------------------------------------------------------------------------------
    // Here we are using the default CC1pi selection. The final cut of this selection (named "likelyGoldenPion") is used to identify events
    // in which the reconstructed & selected pion candidate is "golden". Here a golden pion is one that is: contained within the TPC, 
    // doesn't undergo any scatters on the argon, and comes to rest (before any secondary interaction). If an event passes all the cuts of
    // this selection, it's said to pass the "golden" selection. If an event passes the penultimate cut, it is said to have passed the
    // "generic" selection. Events passing the generic selection are likely to be CC1pi, but may or may not have a golden pion. All events
    // that pass the golden selection also pass the generic selection - but not all events that pass the generic selection also pass the
    // golden selection.
    // auto selection = SelectionHelper::GetDefaultSelection();
    // std::cout<<"..........................................\nUSING CC0pi Selection!\n.........................................."<<std::endl;
    // auto sidebandSelection = SelectionHelper::GetSelection("CC0pi");
    std::cout<<"..........................................\nUSING Modified CC0pi Selection: muonLikeProtonValue=-0.48f, barelyResemblingProtonValue=0.12f\n.........................................."<<std::endl;
    auto sidebandSelection = SelectionHelper::GetCC0piSelectionModified(-0.48f, 0.12f);

    // -------------------------------------------------------------------------------------------------------------------------------------
    // Setup the cross-section objects
    // -------------------------------------------------------------------------------------------------------------------------------------
    // Here we make a map from a name of the cross-section to the cross-section object itself. In this way, we can iterate through the
    // cross-section objects and reduce code-bloat. The first index is an identifier for the selection that's applied (generic or goldlen), 
    // the second index is an identifier for the kinematic quantity that's relevant for the cross-section (e.g. muonMomentum), and the
    // mapped type is the cross-section object.
    typedef std::map<std::string, std::map<std::string, CrossSectionHelper::CrossSection> > CrossSectionMap;
    CrossSectionMap xsecMapSidebandBNB;
    CrossSectionMap xsecMapSidebandNuWro;
    CrossSectionMap xsecMapSidebandGenie;
    CrossSectionMap xsecMapSidebandNuWroTrue;

    // We additionally make a map from each cross-section to the limits of the phase-space that we should consider. The key is the
    // identifier for the kinematic quantity and the mapped value is a pair containing the limits [min, max]
    std::map< std::string, std::pair<float, float> > phaseSpaceMap;

    // ATTN the configuration allows the user to enable or disable each cross-section. If a cross-section has been disabled, then it won't
    // be added the xSecMap. However, the phaseSpaceMap always includes all kinematic parameters!

    // Add the differential cross-sections
    for (const auto &[name, binning, scaleByBinWidth] : std::vector< std::tuple<std::string, Config::Global::Binning, bool> > {

        // The names of the cross-section kinematic parameters, and their binning information.
        // The third (boolean) parameter indicates if the cross-section bins should be scaled by their width
        { "muonCosTheta", config.global.muonCosTheta, true  }, 
        { "muonPhi", config.global.muonPhi, true  }, 
        { "muonMomentum", config.global.muonMomentum, true  }, // todo: check this works

        { "pionCosTheta", config.global.pionCosTheta, true  }, 
        { "pionPhi", config.global.pionPhi, true  }, 
        { "pionMomentum", config.global.pionMomentum, true  }, // todo: check this works 

        { "muonPionAngle", config.global.muonPionAngle, true  }, 
        { "nProtons", config.global.nProtons, false }

    })
    {
        // Add to the phase-space map
        phaseSpaceMap.emplace(name, std::pair<float, float>({binning.min, binning.max}));

        // Here we calculate every cross-section using both the generic and golden selection. In the end we only use the golden selection for
        // the pion momentum, but we additionally apply it to other cross-sections as a cross-check.
        // Add the cross-section object to the map using the binning from the input configuration
        const auto &[extendedBinEdges, hasUnderflow, hasOverflow] = CrossSectionHelper::GetExtendedBinEdges(binning.min, binning.max, binning.binEdges);
        for (const auto &selectionName : {"generic", "golden"})
        {
            // Don't setup a cross-section object if it's been disabled in the configuration
            if (config.extractXSecs.crossSectionIsEnabled.at(selectionName).at(name))
            {
                xsecMapSidebandBNB[selectionName].emplace(name, CrossSectionHelper::CrossSection(systParams, extendedBinEdges, hasUnderflow, hasOverflow, scaleByBinWidth));
                xsecMapSidebandNuWro[selectionName].emplace(name, CrossSectionHelper::CrossSection(systParamsNuWro, extendedBinEdges, hasUnderflow, hasOverflow, scaleByBinWidth));
                xsecMapSidebandGenie[selectionName].emplace(name, CrossSectionHelper::CrossSection(systParamsGenie, extendedBinEdges, hasUnderflow, hasOverflow, scaleByBinWidth));
                xsecMapSidebandNuWroTrue[selectionName].emplace(name, CrossSectionHelper::CrossSection(systParamsNuWro, extendedBinEdges, hasUnderflow, hasOverflow, scaleByBinWidth));
            }
        }
    }
    

    // ATTN here we use the machinary for a differential cross-section, and treat the total cross-section as a single-bin measurement.
    // The "kinematic quantity" in this case is just a dummy parameter. Here we define a single bin with edges arbitrarily chosen to be
    // (-1 -> +1), and we request that the cross-section object does not apply bin-width scaling. When we fill this object, we will use
    // the dummy kinematic quantity with a value of 0 for all events. This is arbitrary, as long as it's within the bin edges we chose.
    // In this way the single bin contains all events. It's just a trick to avoid implementing extra logic for the total cross-section.
    for (const auto &selectionName : {"generic", "golden"})
    {
        // Don't setup a cross-section object if it's been disabled in the configuration
        if (config.extractXSecs.crossSectionIsEnabled.at(selectionName).at("total"))
        {
            xsecMapSidebandBNB[selectionName].emplace("total", CrossSectionHelper::CrossSection(systParams, {-1.f, 1.f}, false, false, false));
            xsecMapSidebandNuWro[selectionName].emplace("total", CrossSectionHelper::CrossSection(systParamsNuWro, {-1.f, 1.f}, false, false, false));
            xsecMapSidebandGenie[selectionName].emplace("total", CrossSectionHelper::CrossSection(systParamsGenie, {-1.f, 1.f}, false, false, false));
            xsecMapSidebandNuWroTrue[selectionName].emplace("total", CrossSectionHelper::CrossSection(systParamsNuWro, {-1.f, 1.f}, false, false, false));
        }
    }

    // The dummy value that will be used as the "kinematic quantity" for the total cross-section
    const auto dummyValue = 0.f;

    // -------------------------------------------------------------------------------------------------------------------------------------
    // Setup the relevent "getters" for each the sideband
    // -------------------------------------------------------------------------------------------------------------------------------------
    // Here we define a map from the name of each cross-section to a function which pulls out the relevant kinematic quanitity from an input
    // analysis data object. Again this is done up-front to reduce code-bloat below.
    ExtractionHelper::AnalysisValueMap getSidebandValue;
    ExtractionHelper::PopulateAnalysisValueMap(getSidebandValue, true); // true: creates sideband getters

    // -------------------------------------------------------------------------------------------------------------------------------------
    // Setup the input files
    // -------------------------------------------------------------------------------------------------------------------------------------
    // Here we define a vector of tuples with 4 entries
    //   - First, the sample type (e.g. overlay)
    //   - Second, a string which is used to identify a given detector variation sample (for other sample type, this is unused)
    //   - Third, the path to the input file
    //   - Fourth, the normalisation factor to apply to all events in that file
    const auto fluxHistNames = CrossSectionHelper::GetNominalFluxHistNames(config.flux.nuPdgsSignal, config.flux.nuPdgToHistName, config.flux.nomHistPattern);
    const auto &[fluxBinEdges, fluxValues] = CrossSectionHelper::ReadNominalFlux(config.flux.fileName, fluxHistNames, config.flux.pot);
    scalingData.pFluxReweightor = std::make_shared<CrossSectionHelper::FluxReweightor>(fluxBinEdges, fluxValues, systParams.fluxDimensions);


    ExtractionHelper::InputFileList inputData;
    float totalExposurePOT;
    ExtractionHelper::PopulateInputFileList(config, inputData, totalExposurePOT);
    
    scalingData.exposurePOT = totalExposurePOT*testingFraction;
    scalingData.nTargets = config.global.targetDensity * (1e-8) * AnalysisHelper::GetFiducialVolume();


    // -------------------------------------------------------------------------------------------------------------------------------------
    // -------------------------------------------------------------------------------------------------------------------------------------
    // -------------------------------------------------------------------------------------------------------------------------------------
    // -------------------------------------------------------------------------------------------------------------------------------------
    // -------------------------------------------------------------------------------------------------------------------------------------
    // -------------------------------------------------------------------------------------------------------------------------------------
    // Count the events for CC0pi
    // -------------------------------------------------------------------------------------------------------------------------------------
    // Loop over the files
    for (const auto &[sampleType, sampleName, fileName, normalisation] : inputData)
    {
        std::cout << "Reading input file: " << fileName << std::endl;

        const auto isOverlay = (sampleType == AnalysisHelper::Overlay);
        const auto isDirt    = (sampleType == AnalysisHelper::Dirt);
        const auto isNuWro   = (sampleType == AnalysisHelper::NuWro);
        const auto isDataBNB = (sampleType == AnalysisHelper::DataBNB);
        const auto isDetVar  = (sampleType == AnalysisHelper::DetectorVariation);
        const auto isDataEXT = (sampleType == AnalysisHelper::DataEXT);

        // if(isDetVar) continue; // Todo: Remove!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        // if(isDataBNB || isDirt || isDataEXT) continue; // Todo: Remove!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        // Open the input file for reading and enable the branches with systematic event weights (if required)
        FileReader reader(fileName);

        if (isOverlay || isNuWro) // Todo: Is isNuWro needed ?
            reader.EnableSystematicBranches();

        auto pEvent = reader.GetBoundEventAddress();

        // Loop over the events in the file
        const auto nEvents = reader.GetNumberOfEvents();

        //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        
        std::cout<<"############################\nOnly using "<<testingFraction*100<<"% of events!\n############################"<<std::endl;
        for (unsigned int i = 0; i < nEvents*testingFraction; ++i) // todo change back to all events!!!!!!!!!!!!!!!!!!!!!!
        {
            AnalysisHelper::PrintLoadingBar(i, int(nEvents*testingFraction));
            reader.LoadEvent(i);

            // -----------------------------------------------------------------------------------------------------------------------------
            // Work out if this event passed the selection and apply any additional phase-space cuts based on the input binning
            // -----------------------------------------------------------------------------------------------------------------------------

            // Run the selection
            const auto &[passedGoldenSidebandSelection, sidebandCutsPassed, sidebandAssignedPdgCodes] = sidebandSelection.Execute(pEvent);

            const auto passedGenericSidebandSelection = SelectionHelper::IsCutPassed(sidebandCutsPassed, config.global.lastCutGeneric);

            // Get the reco analysis data (if available, otherwise set to dummy values)
            const auto recoData = (
                passedGenericSidebandSelection
                    ? AnalysisHelper::GetRecoAnalysisDataCC0Pi(pEvent->reco, sidebandAssignedPdgCodes, passedGoldenSidebandSelection)
                    : AnalysisHelper::GetDummyAnalysisData()
            );

            // Here we apply reco-level phase-space restrictions
            // For any event that passes the generic selection, get the value of the kinematic quantity and check if it is outside of the
            // min/max values supplied in the binning. If so, then reject the event.
            bool passesPhaseSpaceReco = false;
            if (passedGenericSidebandSelection)
            {
                // Start by assuming the event passes the phase-space cuts
                passesPhaseSpaceReco = true;

                // Check the value of the kinematic quantities are within the phase-space limits
                for (const auto &[name, minMax] : phaseSpaceMap)
                {
                    const auto &[min, max] = minMax;
                    const auto value = getSidebandValue.at(name)(recoData);

                    if (value < min || value > max)
                    {
                        passesPhaseSpaceReco = false;
                        break;
                    }
                }
            }

            const auto isSelectedSidebandGeneric = passedGenericSidebandSelection && passesPhaseSpaceReco;
            const auto isSelectedSidebandGolden = passedGoldenSidebandSelection && passesPhaseSpaceReco;

            // std::map<std::string, bool> isSelectedMap = {{"sideband", isSelectedSideband}};
            std::map<std::string, bool> isSelectedMap = {{"generic", isSelectedSidebandGeneric}, {"golden", isSelectedSidebandGolden}};

            // -----------------------------------------------------------------------------------------------------------------------------
            // Handle 'data' (BNB & fake-data)
            // -----------------------------------------------------------------------------------------------------------------------------
            const auto weight = AnalysisHelper::GetNominalEventWeight(pEvent) * normalisation;
            for (auto &[selectionName, xsecs] : xsecMapSidebandBNB)
            {
                // Determine if we passed the relevant selection
                const auto isSelected = isSelectedMap.at(selectionName);
                // Only count events passing the selection
                if (!isSelected) continue;

                for (auto &[name, xsec] : xsecs)
                {
                    if (isDataBNB)
                    {
                        xsec.AddSelectedBNBDataEvent(getSidebandValue.at(name)(recoData));
                    }
                    else if (isNuWro)
                    {
                        xsecMapSidebandNuWro.at(selectionName).at(name).AddWeightedSelectedDataEvent(getSidebandValue.at(name)(recoData), weight);
                    }
                    else if(config.global.useGenieAsData && isOverlay)
                    {
                        std::cout<<"DEBUG Part 1 - event: "<<i<<" - weight: "<<weight<<" - selectionName: "<<selectionName<<" - name: "<<name<<std::endl;
                        xsecMapSidebandGenie.at(selectionName).at(name).AddWeightedSelectedDataEvent(getSidebandValue.at(name)(recoData), weight);
                    }
                }
            }            
            if (isDataBNB) continue; // For BNB data that's all we need to do!

            // if (isDataBNB)
            // {
            //     for (auto &[selectionName, xsecs] : xsecMapSidebandBNB)
            //     {
            //         // Determine if we passed the relevant selection
            //         const auto isSelected = isSelectedMap.at(selectionName);
            //         // Only count events passing the selection
            //         if (!isSelected) continue;

            //         for (auto &[name, xsec] : xsecs)
            //             xsec.AddSelectedBNBDataEvent(getSidebandValue.at(name)(recoData));
            //     }

            //     // For BNB data that's all we need to do!
            //     continue;
            // }

            // // Get the nominal event weight, scaled by the sample normalisation
            // const auto weight = AnalysisHelper::GetNominalEventWeight(pEvent) * normalisation;

            // if (isNuWro)
            // {
            //     for (auto &[selectionName, xsecs] : xsecMapSidebandNuWro)
            //     {
            //         // Determine if we passed the relevant selection
            //         const auto isSelected = isSelectedMap.at(selectionName);
            //         // Only count events passing the selection
            //         if (!isSelected) continue;

            //         for (auto &[name, xsec] : xsecs)
            //             xsec.AddWeightedSelectedDataEvent(getSidebandValue.at(name)(recoData), weight);
            //     }
            // }

            // if (isOverlay)
            // {
            //     for (auto &[selectionName, xsecs] : xsecMapSidebandGenie)
            //     {
            //         // Determine if we passed the relevant selection
            //         const auto isSelected = isSelectedMap.at(selectionName);
            //         // Only count events passing the selection
            //         if (!isSelected) continue;

            //         for (auto &[name, xsec] : xsecs)
            //             xsec.AddWeightedSelectedDataEvent(getSidebandValue.at(name)(recoData), weight);
            //     }
            // }


            // -----------------------------------------------------------------------------------------------------------------------------
            // Work out if this event is signal, and apply any phase-space restrictions based on the input binning
            // -----------------------------------------------------------------------------------------------------------------------------

            // Determine if this is truly a CC0Pi event
            const auto isTrueCC0Pi = (isOverlay || isDetVar || isNuWro) && AnalysisHelper::IsTrueCC0Pi(pEvent, config.global.useAbsPdg, config.global.protonMomentumThreshold);

            // Get the truth analysis data (if available, otherwise set to dummy values)
            const auto truthData = (
                (isTrueCC0Pi)
                    ? AnalysisHelper::GetTruthAnalysisDataCC0Pi(pEvent->truth, config.global.useAbsPdg, config.global.protonMomentumThreshold)
                    : AnalysisHelper::GetDummyAnalysisData()
            );

            // Here we apply truth-level phase-space restrictions
            // For all true CC1Pi events, we check if the values of each kinematic variable are within the supplied limits. If not then the
            // event is not classed as "signal"
            bool passesPhaseSpaceTruth = false;
            if (isTrueCC0Pi)
            {
                // Start by assuming the event passes the phase-space cuts
                passesPhaseSpaceTruth = true;

                // Check the value of the kinematic quantities are within the phase-space limits
                for (const auto &[name, minMax] : phaseSpaceMap)
                {
                    if(name == "pionMomentum") continue; // Not really compatible with the pion momentum in the CC1pi selection
                    const auto &[min, max] = minMax;
                    const auto value = getSidebandValue.at(name)(truthData);

                    if (value < min || value > max)
                    {
                        passesPhaseSpaceTruth = false;
                        break;
                    }
                }
            }

            const auto isCC0PiSignal = isTrueCC0Pi && passesPhaseSpaceTruth;
            // std::map<std::string, bool> isSignalMap = {{"sideband", isCC0PiSignal}};
            std::map<std::string, bool> isSignalMap = {{"generic", isCC0PiSignal}, {"golden", isCC0PiSignal}};

            // -----------------------------------------------------------------------------------------------------------------------------
            // Handle the detector variation samples (unisims)
            // -----------------------------------------------------------------------------------------------------------------------------
            if (isDetVar)
            {
                for (auto &[selectionName, xsecs] : xsecMapSidebandBNB)
                {
                    const auto isSignal = isSignalMap.at(selectionName);
                    // Handle signal events
                    if (isSignal)
                    {
                        // Determine if we passed the relevant selection
                        const auto isSelected = isSelectedMap.at(selectionName);
                        for (auto &[name, xsec] : xsecs)
                        {
                            const auto recoValue = getSidebandValue.at(name)(recoData);
                            const auto trueValue = getSidebandValue.at(name)(truthData);
                            xsec.AddSignalEventDetVar(recoValue, trueValue, isSelected, weight, sampleName);
                        }
                    }
                    // Handle selected background events
                    else
                    {
                        for (auto &[selectionName, xsecs] : xsecMapSidebandBNB)
                        {
                            // Only use selected background events
                            const auto isSelected = isSelectedMap.at(selectionName);
                            if (!isSelected)
                                continue;

                            for (auto &[name, xsec] : xsecs)
                            {
                                const auto recoValue = getSidebandValue.at(name)(recoData);
                                xsec.AddSelectedBackgroundEventDetVar(recoValue, weight, sampleName);
                            }
                        }
                    }
                }

                // For detector variation samples, that's all we need to do!
                continue;
            }

            // -----------------------------------------------------------------------------------------------------------------------------
            // Handle all other events (i.e those from the nominal simulation): Overlays, dirt, EXT data
            // -----------------------------------------------------------------------------------------------------------------------------
            // Get the flux weights
            const auto fluxWeights = (
                isOverlay
                    ? CrossSectionHelper::GetWeightsMap(pEvent->truth, systParams.fluxDimensions, config.extractXSecs.mutuallyExclusiveDimensions)
                    : CrossSectionHelper::GetUnitWeightsMap(systParams.fluxDimensions)
            );

            // Fill the flux-reweightor with all overlay events from a neutrinos of the desired flavour in the active volume
            if (isOverlay && AnalysisHelper::IsInActiveVolume(pEvent->truth.nuVertex()) &&
                std::find(config.flux.nuPdgsSignal.begin(), config.flux.nuPdgsSignal.end(), pEvent->truth.nuPdgCode()) != config.flux.nuPdgsSignal.end())
            {
                scalingData.pFluxReweightor->AddEvent(pEvent->truth.nuEnergy(), weight, fluxWeights);
            }

            // Get the cross-section weights
            // ATTN here we optionally scale the cross-section weights down by the genieTuneEventWeight - this is done so we don't
            // double count this weight (once in the nominal event weight, and once in the xsec systematic event weights)
            // const auto xsecWeightsScaleFactor = 1.f;
            auto xsecWeightsScaleFactor = (isOverlay && config.extractXSecs.scaleXSecWeights) ? pEvent->truth.genieTuneEventWeight() : 1.f;
            // if(xsecWeightsScaleFactor<0.1f && isCC0PiSignal) std::cout<<"xsecWeightsScaleFactor: "<<xsecWeightsScaleFactor<< " - event: " << i <<" - true protonMultiplicity: "<< getSidebandValue.at("nProtons")(truthData) << std::endl;
            // xsecWeightsScaleFactor = std::max(xsecWeightsScaleFactor, 0.0001f);
            if(xsecWeightsScaleFactor<=std::numeric_limits<float>::epsilon())// || std::isinf(xsecWeightsScaleFactor))
            {
                if(weight>std::numeric_limits<float>::epsilon())
                {
                    std::cout<<"xsecWeightsScaleFactor is close to zero but nominal weight is not!"<<std::endl;
                    std::cout<<"Skipped event: "<<i<<" - weight: "<<weight<<" - xsecWeightsScaleFactor: "<<xsecWeightsScaleFactor<<std::endl;
                    throw std::runtime_error("xsecWeightsScaleFactor is close to zero but nominal weight is not!"); 
                }
                std::cout<<"Skipped: xsecWeightsScaleFactor and nominal weight are close to zero!"<<std::endl;
                continue;
            }

            const auto xsecWeights = (
                isOverlay
                    ? CrossSectionHelper::ScaleWeightsMap(CrossSectionHelper::GetWeightsMap(pEvent->truth, systParams.xsecDimensions, config.extractXSecs.mutuallyExclusiveDimensions), xsecWeightsScaleFactor, config.extractXSecs.xsecUBGenieScalingMap)
                    : CrossSectionHelper::GetUnitWeightsMap(systParams.xsecDimensions)
            );

            // const auto weightVector = xsecWeights.at("VecFFCCQEshape_UBGenie");
            // for(const auto &v: weightVector)
            // {
            //     if(v>5.f) std::cout<<"weightVector: "<<v<<" - isCC0PiSignal: "<<isCC0PiSignal<<" - run: "<< fileName << " - event: " << i <<" - true protonMultiplicity: "<< getSidebandValue.at("nProtons")(truthData) << std::endl;
            // }

            // Get the reinteraction weights
            const auto reintWeights = (
                isOverlay
                    ? CrossSectionHelper::GetWeightsMap(pEvent->truth, systParams.reintDimensions, config.extractXSecs.mutuallyExclusiveDimensions)
                    : CrossSectionHelper::GetUnitWeightsMap(systParams.reintDimensions)
            );
            
            // Unit weights to be used for fake-data studies (NuWro/Genie)
            // const auto fluxUnitWeights = CrossSectionHelper::GetUnitWeightsMap(systParams.fluxDimensions);
            // const auto reintUnitWeights = CrossSectionHelper::GetUnitWeightsMap(systParams.reintDimensions);
            // const auto xsecUnitWeights = CrossSectionHelper::GetUnitWeightsMap(systParams.xsecDimensions);
            CrossSectionHelper::SystFloatMap emptySystFloatMap = {};
            std::cout<<"DEBUG Part 1.1 - event: "<<i<<std::endl;
            for (auto &[selectionName, xsecs] : xsecMapSidebandNuWro)
            {
                const auto isSignal = isSignalMap.at(selectionName);
                const auto isSelected = isSelectedMap.at(selectionName);
                // Handle signal events
                if (isSignal)
                {
                    for (auto &[name, xsec] : xsecs)
                    {
                        const auto recoValue = getSidebandValue.at(name)(recoData);
                        const auto trueValue = getSidebandValue.at(name)(truthData);
                        if(isNuWro)
                        {
                            xsecMapSidebandNuWroTrue.at(selectionName).at(name).AddSignalEvent(recoValue, trueValue, isSelected, weight, emptySystFloatMap, xsecWeights, emptySystFloatMap);
                        } else
                        {
                            if(isOverlay) // Genie MC for NuWro (no EXT or dirt events in NuWro)
                            {
                                if(config.global.useNuWroAsData) xsec.AddSignalEvent(recoValue, trueValue, isSelected, weight, emptySystFloatMap, xsecWeights, emptySystFloatMap);
                                std::cout<<"DEBUG Part 2A - event: "<<i<<" - weight: "<<weight<<" - selectionName: "<<selectionName<<" - name: "<<name<<" - isSelected: "<<isSelected<<std::endl;
                                if(config.global.useGenieAsData) xsecMapSidebandGenie.at(selectionName).at(name).AddSignalEvent(recoValue, trueValue, isSelected, weight, emptySystFloatMap, emptySystFloatMap, emptySystFloatMap);
                            }
                            if(config.global.useBNBAsData && (isOverlay || isDataEXT || isDirt)) // Genie MC for Data
                            {
                                xsecMapSidebandBNB.at(selectionName).at(name).AddSignalEvent(recoValue, trueValue, isSelected, weight, fluxWeights, xsecWeights, reintWeights);
                            }
                        }
                    }
                }
                // Handle selected background events
                else if (isSelected)
                {
                    for (auto &[name, xsec] : xsecs)
                    {
                        const auto recoValue = getSidebandValue.at(name)(recoData);
                        std::vector<float> bootstrapWeights; // Parameter only needed for CC1pi
                        if(isNuWro)
                        {
                            xsecMapSidebandNuWroTrue.at(selectionName).at(name).AddSelectedBackgroundEvent(recoValue, false, weight, emptySystFloatMap, xsecWeights, emptySystFloatMap, bootstrapWeights);
                        } 
                        else
                        {
                            if(isOverlay) // Genie MC for NuWro (no EXT or dirt events in NuWro)
                            {
                                if(config.global.useNuWroAsData) xsec.AddSelectedBackgroundEvent(recoValue, false, weight, emptySystFloatMap, xsecWeights, emptySystFloatMap, bootstrapWeights);
                                std::cout<<"DEBUG Part 2B - event: "<<i<<" - weight: "<<weight<<" - selectionName: "<<selectionName<<" - name: "<<name<<std::endl;
                                if(config.global.useGenieAsData) xsecMapSidebandGenie.at(selectionName).at(name).AddSelectedBackgroundEvent(recoValue, false, weight, emptySystFloatMap, emptySystFloatMap, emptySystFloatMap, bootstrapWeights);
                            }
                            if((isOverlay || isDataEXT || isDirt) && config.global.useBNBAsData) // Genie MC for Data
                            {
                                xsecMapSidebandBNB.at(selectionName).at(name).AddSelectedBackgroundEvent(recoValue, isDirt, weight, fluxWeights, xsecWeights, reintWeights, bootstrapWeights);
                            }
                        }
                    }
                }
            }
        }
    }
    std::cout<<"-----------------Finished processing events-----------------"<<std::endl;




    // -------------------------------------------------------------------------------------------------------------------------------------
    // Calculate the sideband weights
    // -------------------------------------------------------------------------------------------------------------------------------------
    // Loop over all cross-section objects
    typedef std::pair<std::vector<Double_t>, std::vector<Double_t>> paramAndErrorPair; // Todo: Improve code!
    //Parameters: dataTypeName selectionName name
    typedef std::map<std::string, std::map<std::string, std::map<std::string, paramAndErrorPair>>> nominalFitMap;
    //Parameters: dataTypeName selectionName name paramName (i.e. golden muonMomentum hadronProduction)
    typedef std::map<std::string, std::map<std::string, std::map<std::string, std::map<std::string, std::vector<paramAndErrorPair>>>>> universeFitMap;

    nominalFitMap cc0piNominalConstraintMap;
    universeFitMap cc0piUniverseConstraintMap;

    std::cout<<"DEBUG Pre-Fit 0"<<std::endl;
    std::vector<std::pair<std::string, CrossSectionMap>> xsecMapVector;
    if(config.global.useBNBAsData) xsecMapVector.push_back(std::make_pair(std::string("BNB"), xsecMapSidebandBNB));
    if(config.global.useNuWroAsData) xsecMapVector.push_back(std::make_pair(std::string("NuWro"), xsecMapSidebandNuWro));
    if(config.global.useGenieAsData) xsecMapVector.push_back(std::make_pair(std::string("Genie"), xsecMapSidebandGenie));

    for (const auto &[dataTypeName, xsecMap] : xsecMapVector) // {std::make_pair(std::string("NuWro"), xsecMapSidebandNuWro), std::make_pair(std::string("BNB"), xsecMapSidebandBNB)})
    {
        for (const auto &[selectionName, xsecs] : xsecMap)
        {
            for (const auto &[name, xsec] : xsecs)
            {
                // -------------------------------------------------------------------------------------------------------------------------------------
                // Fit nominal
                // -------------------------------------------------------------------------------------------------------------------------------------

                std::cout<<"_______________________ExtractNuWroSidebandFit2 Fitting: "<<dataTypeName<<" - "<<selectionName<<" - "<<name<<"_______________________"<<std::endl;

                // Get the smearing matrix of selected events
                const auto smearingMatrix = xsec.GetSmearingMatrixAllSelected();
                const auto smearingMatrixGeneral = xsec.GetSmearingMatrix();

                const auto selectedEventsBackgroundReco = xsec.GetSelectedBackgroundEvents();
                const auto selectedEventsData = xsec.GetSelectedBNBDataEvents();
                const auto selectedEventsSignalTruth = xsec.GetSelectedSignalEvents();
                const auto eventsSignalTruth = xsec.GetSignalEvents();
                const auto metadata = xsec.GetMetadata();
                auto signalData = selectedEventsData - selectedEventsBackgroundReco;

                FormattingHelper::SaveMatrix(eventsSignalTruth, dataTypeName + "_SidebandFit_" + selectionName + "_" + name + "_eventsSignalTruth.txt");
                FormattingHelper::SaveMatrix(selectedEventsSignalTruth, dataTypeName + "_SidebandFit_" + selectionName + "_" + name + "_selectedEventsSignalTruth.txt");
                FormattingHelper::SaveMatrix(selectedEventsData, dataTypeName + "_SidebandFit_" + selectionName + "_" + name + "_selectedEventsData.txt");
                FormattingHelper::SaveMatrix(selectedEventsBackgroundReco, dataTypeName + "_SidebandFit_" + selectionName + "_" + name + "_selectedEventsBackgroundReco.txt");
                FormattingHelper::SaveMatrix(smearingMatrix, dataTypeName + "_SidebandFit_" + selectionName + "_" + name + "_smearingMatrix.txt");
                FormattingHelper::SaveMatrix(smearingMatrixGeneral, dataTypeName + "_SidebandFit_" + selectionName + "_" + name + "_smearingMatrixGeneral.txt");

                // for(unsigned int r = 0; r<signalData.GetRows(); r++) // DEBUG - TODO: REMOVE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                // {
                //     for(unsigned int c = 0; c<signalData.GetColumns(); c++) // DEBUG - TODO: REMOVE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                //     {
                //         if(signalData.At(r, c) < 0)
                //         {
                //             std::cout<<"signalData.At("<<r<<", "<<c<<") = "<<signalData.At(r, c)<<std::endl;
                //             signalData.SetElement(r, c, std::max(signalData.At(r, c), 0.f));
                //         }
                //     }
                // }

                // std::map<std::string, ubsmear::UBMatrix> totalErrorMatrixMap;
                // std::vector<float> elements;
                const auto nBins = selectedEventsData.GetRows();
                std::cout<<"DEBUG - nBins: "<<nBins<<std::endl;
                // for (unsigned int iBin = 0; iBin < nBins; ++iBin)
                // {
                //     const auto value = selectedEventsData.At(iBin, 0);
                //     if (value<0.f)
                //     {
                //         std::cout<<"ERROR: ExtractXSec - Background-removed signal data is negative."<<std::endl;
                //         throw std::logic_error("ERROR: ExtractXSec - Background-removed signal data is negative.");
                //     }
                //     elements.push_back(AnalysisHelper::GetCountUncertainty(value));
                // }
                // const ubsmear::UBMatrix statUncertainties(elements, nBins, 1);
                // const auto statVariances = ubsmear::ElementWiseOperation(statUncertainties, statUncertainties, [](const auto &l, const auto &r) { return l * r; });
                // const auto statErrorMatrix = ubsmear::UBMatrixHelper::GetDiagonalMatrix(statVariances);
                // totalErrorMatrixMap.emplace("data", statErrorMatrix); // todo improve code

                // pMetadataNuWro = &metadata;
                // pSmearingMatrixNuWro = &smearingMatrix;
                // pPredictionNuWro = &selectedEventsSignalTruth;
                // pDataNuWro = &signalData;
                // pTotalErrorMatrixMap = &totalErrorMatrixMap; // todo check whether individual universe uncertainties need to be generated (due to differenced is genie backgrounds in each universe)


                ///////////////////////////////////////////
                /////////////////////////////////////////// Start cross-section values
                auto data = xsec.GetBNBDataCrossSection(scalingData);
                std::cout<<dataTypeName<<": nominal data:"<<std::endl;
                data.Print();
                for(unsigned int r = 0; r<data.GetRows(); r++) // DEBUG - TODO: REMOVE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                {
                    for(unsigned int c = 0; c<data.GetColumns(); c++) // DEBUG - TODO: REMOVE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    {
                        if(data.At(r, c) < 0)
                        {
                            std::cout<<"data.At("<<r<<", "<<c<<") = "<<data.At(r, c)<<std::endl;
                            data.SetElement(r, c, std::max(data.At(r, c), 0.f));
                        }
                    }
                }
                // std::cout<<dataTypeName<<": nominal data (after max()):"<<std::endl;
                // data.Print();

                // FormattingHelper::SaveMatrix(data, dataTypeName + "_SidebandFit_" + selectionName + "_" + name + "_data_"+postfix+".txt");

                // const auto dataStatUncertainties = xsec.GetBNBDataCrossSectionStatUncertainty(scalingData);
                // FormattingHelper::SaveMatrix(dataStatUncertainties, dataTypeName + "_SidebandFit_" + selectionName + "_" + name + "_data_stat_"+postfix+".txt");

                // const auto dataSystBiasCovariances = xsec.GetBNBDataCrossSectionSystUncertainties(scalingData);
                // for (const auto &[group, map] : dataSystBiasCovariances)
                // {
                //     for (const auto &[paramName, biasCovariance] : map)
                //     {
                //         const auto &[pBias, pCovariance] = biasCovariance;
                //         FormattingHelper::SaveMatrix(*pBias, dataTypeName + "_SidebandFit_" + selectionName + "_" + name + "_data_" + group + "_" + paramName + "_bias_"+postfix+".txt");
                //         FormattingHelper::SaveMatrix(*pCovariance, dataTypeName + "_SidebandFit_" + selectionName + "_" + name + "_data_" + group + "_" + paramName + "_covariance_"+postfix+".txt");
                //     }
                // }
                // const auto smearingMatrixSystBiasCovariances = xsec.GetSmearingMatrixSystUncertainties();
                const auto prediction = xsec.GetPredictedCrossSection(scalingData);
                std::cout<<dataTypeName<<": nominal prediction:"<<std::endl;
                prediction.Print();

                std::cout<<dataTypeName<<": nominal smearingMatrixGeneral:"<<std::endl;
                smearingMatrixGeneral.Print();

                const auto smearedPrediction = smearingMatrixGeneral * prediction;
                std::cout<<dataTypeName<<": nominal smearedPrediction:"<<std::endl;
                smearedPrediction.Print();
                // FormattingHelper::SaveMatrix(prediction, dataTypeName + "_SidebandFit_" + selectionName + "_" + name + "_prediction_"+postfix+".txt");

                // const auto &[pPredictionStatBias, pPredictionStatCovariance] = xsec.GetPredictedCrossSectionStatUncertainty(scalingData);
                // const auto predictionErrorMatrix = CrossSectionHelper::GetErrorMatrix(*pPredictionStatBias, *pPredictionStatCovariance);
                
                const auto totalErrorMatrixMap = xsec.GetTotalErrorMatrixMap(scalingData);
                std::cout<<dataTypeName<<": nominal totalErrorMatrixMap.at(\"data\"):"<<std::endl;
                totalErrorMatrixMap.at("data").Print();
                // std::cout<<dataTypeName<<": nominal totalErrorMatrixMap.at(\"smearingMatrix\"):"<<std::endl;
                // totalErrorMatrixMap.at("smearingMatrix").Print();
                /////////////////////////////////////////// End cross-section values
                ///////////////////////////////////////////


                // pPredictionErrorMatrixNuWro = &predictionErrorMatrix;
                pMetadataNuWro = &metadata;
                pSmearingMatrixNuWro = &smearingMatrixGeneral; //todo check which smearing matrix is needed
                pPredictionNuWro = &prediction;
                pDataNuWro = &data;
                pTotalErrorMatrixMap = &totalErrorMatrixMap; // todo check whether individual universe uncertainties need to be generated (due to differenced is genie backgrounds in each universe)

                // std::cout<<"pSmearingMatrixNuWro: "<<std::endl;
                // pSmearingMatrixNuWro->Print();
                // std::cout<<"pDataNuWro: "<<std::endl;
                // pDataNuWro->Print();
                // std::cout<<"pPredictionNuWro: "<<std::endl;
                // pPredictionNuWro->Print();

                // Something like this needed ???
                // const auto dataStatUncertainties = xsec.GetBNBDataCrossSectionStatUncertainty(scalingData);
                // const auto dataSystBiasCovariances = xsec.GetBNBDataCrossSectionSystUncertainties(scalingData);
                // e.g. named:
                // const auto dataStatUncertainties = xsec.GetBNBDataStatUncertainty(scalingData);
                // const auto dataSystBiasCovariances = xsec.GetBNBDataSystUncertainties(scalingData);

                auto minimizer = FittingHelper(nBins);
                std::pair<std::vector<Double_t>, std::vector<Double_t>> result;
                
                std::vector<float> fitCovMatrixVector;
                bool successful = false;
                std::cout<<"DEBUG Nom Fit Point -1"<<std::endl;
                minimizer.Fit(GetChi2NuWro, result, successful, fitCovMatrixVector, 1);
                if(xsec.HasUnderflow())
                {
                    result.first.at(0) = 1.;
                    result.second.at(0) = 0.;
                }
                if(xsec.HasOverflow())
                {
                    result.first.at(nBins-1) = 1.;
                    result.second.at(nBins-1) = 0.;
                }
                std::cout<<"DEBUG Nom Fit Point 0"<<std::endl;
                if(!successful)
                {
                    std::cout<<"ERROR: ExtractNuWroSidebandFit2 - Fit failed."<<std::endl;
                    throw std::logic_error("ERROR: ExtractNuWroSidebandFit2 - Nominal fit failed.");
                }
                std::cout<<"DEBUG Nom Fit Point 1"<<std::endl;

                const ubsmear::UBMatrix sidebandCovMatrix(fitCovMatrixVector, nBins, nBins);
                FormattingHelper::SaveMatrix(sidebandCovMatrix, dataTypeName + "_SidebandFit_" + selectionName + "_" + name + "_sideband_stat_covariance.txt");
                vector<float> paramVector(result.first.begin(), result.first.end()); // Todo avoid this
                vector<float> paramErrorVector(result.second.begin(), result.second.end()); // Todo avoid this
                const ubsmear::UBMatrix sidebandParamVectorTruth(paramVector, nBins, 1);
                const ubsmear::UBMatrix sidebandErrorVectorTruth(paramErrorVector, nBins, 1);
                FormattingHelper::SaveMatrix(sidebandParamVectorTruth, dataTypeName + "_SidebandFit_" + selectionName + "_" + name + "_sideband_parameterVector.txt");
                FormattingHelper::SaveMatrix(sidebandErrorVectorTruth, dataTypeName + "_SidebandFit_" + selectionName + "_" + name + "_sideband_parameterErrorVector.txt");
                std::cout<<"nominal params:"<<std::endl;
                sidebandParamVectorTruth.Print();
                std::cout<<"nominal param errors:"<<std::endl;
                sidebandErrorVectorTruth.Print();
                
                std::cout<<"DEBUG Nom Fit Point 2"<<std::endl;

                cc0piNominalConstraintMap[dataTypeName][selectionName].emplace(name, result);

                std::cout<<"DEBUG Nom Fit Point 3"<<std::endl;
                
                if(!config.global.fitInSystematicUniverses) continue;

                // -------------------------------------------------------------------------------------------------------------------------------------
                // Fit each universe
                // -------------------------------------------------------------------------------------------------------------------------------------
                std::cout<<"DEBUG Nom Fit Point 4"<<std::endl;
                std::map<std::string, CrossSectionHelper::SystDimensionsMap> weightDimensionMap;
                // const auto weightDimensions = {std::make_pair("xsec", systParams.xsecDimensions)};
                if(dataTypeName!="Genie") // Do not fit this for Genie fake-data
                {
                    weightDimensionMap.emplace("xsec", systParams.xsecDimensions);
                }
                if(dataTypeName=="BNB") // Do not fit this for fake-data 
                {
                    weightDimensionMap.emplace("reint", systParams.reintDimensions);
                    weightDimensionMap.emplace("flux", systParams.fluxDimensions);
                }

                for (const auto &[group, dimensions] : weightDimensionMap)
                {
                    for (const auto &[paramName, nUniverses] : dimensions)
                    {
                        std::cout<<"Fitting In Universes: "<<selectionName<<" - "<<name<<" - group: "<<group<<" - paramName: "<<paramName<<"++++++++++++++"<<std::endl;
                        std::vector<paramAndErrorPair> resultVector;
                        for (unsigned int iUni = 0; iUni < nUniverses; ++iUni)
                        {
                            AnalysisHelper::PrintLoadingBar(iUni, nUniverses);
                            // const auto selectedSignalTruthUniverses = xsec.GetSelectedSignalRecoTruthMap().at(group).at(paramName).at(iUni);
                            // const auto selectedBackgroundRecoUniverses = xsec.GetSelectedBackgroundRecoMap().at(group).at(paramName).at(iUni);
                            std::pair<std::vector<Double_t>, std::vector<Double_t>> result;
                            // const auto pSmearingMatrixInUniverse = xsec.GetSmearingMatrixInUniverseAllSelected(group, paramName, iUni);//xsec.GetSmearingMatrixAllSelected();
                            const auto pSmearingMatrixInUniverse = xsec.GetSmearingMatrixInUniverse(group, paramName, iUni);//todo check if this is righ or GetSmearingMatrixInUniverseAllSelected should be used 

                            // std::cout<<"+ExtractNuWroSidebandFit2 Fitting Point 0"<<std::endl;
                            if(pSmearingMatrixNuWro)
                            {
                                const auto smearingMatrixInUniverse = *pSmearingMatrixInUniverse;

                                // // AnalysisHelper::PrintLoadingBar(iUni, nUniverses);
                                // const auto selectedSignalTruthInUniverse = xsec.GetSignalSelectedTrue(selectedSignalTruthUniverses);
                                // const auto selectedBackgoundRecoInUniverse = CrossSectionHelper::GetMatrixFromHist(selectedBackgroundRecoUniverses);
                                // auto signalDataInUniverse = selectedEventsData - selectedBackgoundRecoInUniverse;

                                // std::cout<<"+ExtractNuWroSidebandFit2 Fitting Point 1"<<std::endl;
                                // for(unsigned int r = 0; r<signalDataInUniverse.GetRows(); r++) // DEBUG - TODO: REMOVE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                                // {
                                //     for(unsigned int c = 0; c<signalDataInUniverse.GetColumns(); c++) // DEBUG - TODO: REMOVE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                                //     {
                                //         if(signalDataInUniverse.At(r, c) < 0)
                                //         {
                                //             std::cout<<"universe signalData.At("<<r<<", "<<c<<") = "<<signalDataInUniverse.At(r, c)<<std::endl;
                                //             signalDataInUniverse.SetElement(r, c, std::max(signalDataInUniverse.At(r, c), 0.f));
                                //         }
                                //     }
                                // }

                                // pPredictionNuWro = &selectedSignalTruthInUniverse;
                                // pDataNuWro = &signalDataInUniverse;
                                // pSmearingMatrixNuWro = &smearingMatrixInUniverse;
 
                                const auto dataInUniverse = xsec.GetBNBDataCrossSectionInUniverse(group, paramName, iUni, scalingData);
                                const auto predictionInUniverse = xsec.GetPredictedCrossSectionInUniverse(group, paramName, iUni, scalingData);
                                pDataNuWro = &dataInUniverse; 
                                pSmearingMatrixNuWro = &smearingMatrixInUniverse;
                                pPredictionNuWro = &predictionInUniverse;
                                //todo: check if totalerrormatrixmap also needs to be calaculated in each universe

                                if(paramName == "VecFFCCQEshape_UBGenie")
                                {
                                    std::cout<<"iUni: "<<iUni<<std::endl;
                                    std::cout<<"pSmearingMatrixInUniverse:"<<std::endl;
                                    pSmearingMatrixInUniverse->Print();
                                    std::cout<<"dataInUniverse:"<<std::endl;
                                    dataInUniverse.Print();
                                    std::cout<<"predictionInUniverse:"<<std::endl;
                                    predictionInUniverse.Print();
                                }

                                // std::cout<<"+ExtractNuWroSidebandFit2 Fitting Point 2"<<std::endl;
                                std::vector<float> covMatrixInUniverse;
                                bool successful = false;
                                // std::cout<<"+ExtractNuWroSidebandFit2 Fitting Point 3"<<std::endl;
                                minimizer.Fit(GetChi2NuWro, result, successful, covMatrixInUniverse, -1);

                                if(!successful)
                                {
                                    std::cout<<"Fit not possible - Minimization failed."<<std::endl;
                                    result = std::make_pair(std::vector<Double_t>(nBins, -1), std::vector<Double_t>(nBins, -1));
                                }

                                if(xsec.HasUnderflow()) //todo check if this is correct
                                {
                                    result.first.at(0) = 1.;
                                    result.second.at(0) = 0.;
                                }
                                if(xsec.HasOverflow())
                                {
                                    result.first.at(nBins-1) = 1.;
                                    result.second.at(nBins-1) = 0.;
                                }
                                // std::cout<<"+ExtractNuWroSidebandFit2 Fitting Point 4"<<std::endl;
                                // std::cout<<"+ExtractNuWroSidebandFit2 Fitting Point 5"<<std::endl;
                            }
                            else
                            {
                                std::cout<<"Fit not possible - Smearing matrix not available."<<std::endl;
                                result = std::make_pair(std::vector<Double_t>(nBins, -1), std::vector<Double_t>(nBins, -1));
                            }
                            resultVector.push_back(result);
                        }

                        cc0piUniverseConstraintMap[dataTypeName][selectionName][name].emplace(paramName, resultVector);
                    }
                }
            }
        }
    }
    std::cout<<"-----------------Finished calculating sideband fit-----------------"<<std::endl;

    //
    if(config.global.useNuWroAsData)
    {

        // std::ofstream ofs1("cc0piGenericNominalConstraintMapNuWro.bin", std::ios::binary);
        // std::ofstream ofs2("cc0piGoldenNominalConstraintMapNuWro.bin", std::ios::binary);
        // boost::archive::binary_oarchive oarch1(ofs1);
        // boost::archive::binary_oarchive oarch2(ofs2);
        // oarch1 << cc0piGenericNominalConstraintMapNuWro;
        // oarch2 << cc0piGoldenNominalConstraintMapNuWro;
        // ofs1.close();
        // ofs2.close();
        // if(config.global.fitInSystematicUniverses) 
        // {
        //     std::ofstream ofs3("cc0piGenericUniverseConstraintMapNuWro.bin", std::ios::binary);
        //     std::ofstream ofs4("cc0piGoldenUniverseConstraintMapNuWro.bin", std::ios::binary);
        //     boost::archive::binary_oarchive oarch3(ofs3);
        //     boost::archive::binary_oarchive oarch4(ofs4);
        //     oarch3 << cc0piGenericUniverseConstraintMapNuWro;
        //     oarch4 << cc0piGoldenUniverseConstraintMapNuWro;
        //     ofs3.close();
        //     ofs4.close();
        // }

        // Loop over all cross-section objects
        for (const auto &[selectionName, xsecs] : xsecMapSidebandNuWroTrue)
        {
            for (const auto &[name, xsec] : xsecs)
            {
                std::cout << "True: Processing sideband: "<<selectionName<< " - " << name << std::endl;

                
                // Get the smearing matrix of selected events
                const auto smearingMatrix = xsec.GetSmearingMatrixAllSelected();
                const auto smearingMatrixGeneral = xsec.GetSmearingMatrix();
                const auto selectedEventsSignalTruth = xsec.GetSelectedSignalEvents();
                const auto eventsSignalTruth = xsec.GetSignalEvents();

                FormattingHelper::SaveMatrix(eventsSignalTruth, "NuWro_SidebandFitTruth_" + selectionName + "_" + name + "_eventsSignalTruth.txt");
                FormattingHelper::SaveMatrix(selectedEventsSignalTruth, "NuWro_SidebandFitTruth_" + selectionName + "_" + name + "_selectedEventsSignalTruth.txt");
                FormattingHelper::SaveMatrix(smearingMatrix, "NuWro_SidebandFitTruth_" + selectionName + "_" + name + "_smearingMatrix.txt");
                FormattingHelper::SaveMatrix(smearingMatrixGeneral, "NuWro_SidebandFitTruth_" + selectionName + "_" + name + "_smearingMatrixGeneral.txt");
            }
        }
    }
    
    // if(config.global.useBNBAsData)
    // {
    //     std::ofstream ofs5("cc0piGenericNominalConstraintMapBNB.bin", std::ios::binary);
    //     std::ofstream ofs6("cc0piGoldenNominalConstraintMapBNB.bin", std::ios::binary);
    //     boost::archive::binary_oarchive oarch5(ofs5);
    //     boost::archive::binary_oarchive oarch6(ofs6);
    //     oarch5 << cc0piGenericNominalConstraintMapBNB;
    //     oarch6 << cc0piGoldenNominalConstraintMapBNB;
    //     ofs5.close();
    //     ofs6.close();
    //     if(config.global.fitInSystematicUniverses) 
    //     {
    //         std::ofstream ofs7("cc0piGenericUniverseConstraintMapBNB.bin", std::ios::binary);
    //         std::ofstream ofs8("cc0piGoldenUniverseConstraintMapBNB.bin", std::ios::binary);
    //         boost::archive::binary_oarchive oarch7(ofs7);
    //         boost::archive::binary_oarchive oarch8(ofs8);
    //         oarch7 << cc0piGenericUniverseConstraintMapBNB;
    //         oarch8 << cc0piGoldenUniverseConstraintMapBNB;
    //         ofs7.close();
    //         ofs8.close();
    //     }
    // }

    std::ofstream ofs1("cc0piNominalConstraintMap.bin", std::ios::binary);
    boost::archive::binary_oarchive oarch1(ofs1);
    oarch1 << cc0piNominalConstraintMap;
    ofs1.close();
    if(config.global.fitInSystematicUniverses) 
    {
        std::ofstream ofs2("cc0piUniverseConstraintMap.bin", std::ios::binary);
        boost::archive::binary_oarchive oarch2(ofs2);
        oarch2 << cc0piUniverseConstraintMap;
        ofs2.close();
    }

    std::cout<<"------------- All done -------------"<<std::endl;
    // return;
}

} // namespace ubcc1pi_macros
