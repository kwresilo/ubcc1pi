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
#include "ubsmear.h"

// #include <fstream> // Todo: not use txt files

// Boost libraries
// #include "binary_iarchive.hpp"
#include "binary_oarchive.hpp"
#include "binary_object.hpp"
#include "map.hpp"
#include "vector.hpp"

using namespace ubcc1pi;

const ubsmear::UBMatrix *pPredictionNuWro, *pDataNuWro, *pSmearingMatrixNuWro;//, *pTotalErrorMatrixMapDataNuWro, *pTotalErrorMatrixMapSmearingMatrixNuWro;//, *pPredictionErrorMatrixNuWro;
const std::map<std::string, ubsmear::UBMatrix> *pTotalErrorMatrixMap;
const ubsmear::UBXSecMeta *pMetadataNuWro;

void GetChi2NuWro(Int_t &npar, Double_t *gin, Double_t &chi2, Double_t *par, Int_t iflag)
{
    // std::cout<<"GetChi2NuWro [[0]]"<<std::endl;
    std::vector<float> paramVec;
    // assert(par->size() != pMetadataNuWro->GetNBins());
    for (unsigned int i =0; i<pMetadataNuWro->GetNBins(); ++i) paramVec.push_back(par[i]);

    // std::cout<<"\n\nGetChi2NuWro [[1]]"<<std::endl;
    // for(const auto &v: paramVec) std::cout<<v<<" ";
    // std::cout<<"paramVec.size(): "<<paramVec.size()<<std::endl;

    const ubsmear::UBMatrix parameters(paramVec, paramVec.size(), 1);
    // std::cout<<"GetChi2NuWro [[1.1]]"<<std::endl;
    // const auto binWidths = pMetadataNuWro->GetBinWidths();

    // std::cout<<"\nGetChi2NuWro [[2]] pPredictionNuWro: "<<std::endl;
    pPredictionNuWro->Print();

    // const auto scaledPrediction = ElementWiseOperation(parameters, *pPredictionNuWro, [](const auto &l, const auto& r) { return l * r; });
    // unsigned int nUniverses = 100u; //10000u; ///< The number of universes to use when propagating the smearing matrix uncertainties to the smeared prediction
    float precision = 0.001f; // std::numeric_limits<float>::epsilon();  ///< The precision to use when finding eigenvalues & eigenvectors

    const auto start = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();

    // std::cout<<"GetChi2NuWro [[3]]"<<std::endl;
    // const auto &[smearedPrediction, smearedPredictionErrorMatrix] = ubsmear::ForwardFold(
    //     *pMetadataNuWro,                                            // The metadata that defines the binning
    //     scaledPrediction, zeroMatrix, //*pPredictionErrorMatrixNuWro,              // The prediction and it's error matrix
    //     pSmearingMatrixNuWro->GetTranspose(), pTotalErrorMatrixMap->at("smearingMatrix")->GetTranspose(),  // The smearing matrix and it's error matrix
    //     nUniverses,                                             // The number of universes to use when propagating the uncertainties
    //     precision);                                             // The precision to use when finding eigenvalues and eigenvectors

    // ######## Start smearing prediction
    // Get the bin widths
    const auto binWidths = pMetadataNuWro->GetBinWidths();

    // Scale the prediction by the fitting parameters
    const auto parameterScaledPrediction = ElementWiseOperation(*pPredictionNuWro, parameters, [](const auto &l, const auto& r) { return l * r; });
    // Scale the prediction up by the bin widths (to get a quantity that's proportional to an event rate so it can be smeared)
    const auto binScaledPrediction = ElementWiseOperation(parameterScaledPrediction, binWidths, [](const auto &l, const auto& r) { return l * r; });

    // Smear the scaled prediction using the smearing matrix and scale the smeared prediction back down by the bin widths
    const auto smearedPrediction = ElementWiseOperation(*pSmearingMatrixNuWro * binScaledPrediction, binWidths, [](const auto &l, const auto& r)
    {
        if (r <= std::numeric_limits<float>::epsilon())
            throw std::invalid_argument("UBSmearingHelper::Smear - input metadata contains a bin with a zero or negative width");
        return l / r;
    });
    // ######## Start smearing prediction

    // std::cout<<"\nGetChi2NuWro [[3.1]] scaledPrediction: "<<std::endl;
    // binScaledPrediction.Print();
    // std::cout<<"\nGetChi2NuWro [[3.2]] smearedPrediction: "<<std::endl;
    // smearedPrediction.Print();
    // std::cout<<"\nGetChi2NuWro [[3.3]] pDataNuWro: "<<std::endl;
    // pDataNuWro->Print();
    // std::cout<<"\nGetChi2NuWro [[3.4]]"<<std::endl;
    // // std::cout<<"\nGetChi2NuWro [[3.3.1]] pSmearingMatrixNuWro: "<<std::endl;
    // pSmearingMatrixNuWro->Print();
    // // std::cout<<"\nGetChi2NuWro [[3.4]] smearedPredictionErrorMatrix: "<<std::endl;
    // smearedPredictionErrorMatrix.Print();
    // std::cout<<"\nGetChi2NuWro [[3.4.1]] pTotalErrorMatrixMap->at(\"data\").GetTranspose(): "<<std::endl;
    // pTotalErrorMatrixMap->at("data").GetTranspose().Print();
    // std::cout<<"\nGetChi2NuWro [[3.4.2]]"<<std::endl;
    // const auto covMatrix = smearedPredictionErrorMatrix.GetTranspose()+pTotalErrorMatrixMapDataNuWro->GetTranspose();
    // // std::cout<<"\nGetChi2NuWro [[3.5.1]] Combined covariance matrix: "<<std::endl;
    // covMatrix.Print();

    // const auto &[chi2, degreesOfFreedom] = ubsmear::GetChi2(smearedPrediction, smearedPredictionErrorMatrix.GetTranspose(), *pDataNuWro, pTotalErrorMatrixMapDataNuWro->GetTranspose(), precision);
    // std::tie(chi2, std::ignore) = ubsmear::GetChi2(smearedPrediction, smearedPredictionErrorMatrix.GetTranspose(), *pDataNuWro, pTotalErrorMatrixMapDataNuWro->GetTranspose(), precision);
    const auto nBins = smearedPrediction.GetRows();
    auto zeroMatrix =  ubsmear::UBMatrixHelper::GetZeroMatrix(nBins, nBins);
    // std::cout<<"\nGetChi2NuWro [[3.5]]"<<std::endl;
    std::tie(chi2, std::ignore) = ubsmear::GetChi2(smearedPrediction, zeroMatrix, *pDataNuWro, pTotalErrorMatrixMap->at("data").GetTranspose(), precision);
    // std::cout<<"\nGetChi2NuWro [[3.5.1]]"<<std::endl;

    const auto stop = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
    // std::cout <<"GetChi2NuWro [[3.6]]: "<< stop-start << std::endl;
    // std::cout<<"GetChi2NuWro [[4]]: "<<chi2<<std::endl;
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
    systParams.nBootstrapUniverses = config.extractXSecs.nBootstrapUniverses;
    systParams.fluxDimensions = config.extractXSecs.fluxDimensions;
    systParams.xsecDimensions = config.extractXSecs.xsecDimensions;
    systParams.reintDimensions = config.extractXSecs.reintDimensions;
    systParams.detVarDimensions = config.extractXSecs.detVarDimensions;
    systParams.potFracUncertainty = config.extractXSecs.potFracUncertainty;

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
    // Here we use a FluxReweightor to specify the flux. For each universe of each flux systematic paramter, this is uses the ratio of the
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
    auto sidebandSelection = SelectionHelper::GetSelection("CC0pi");

    // -------------------------------------------------------------------------------------------------------------------------------------
    // Setup the cross-section objects
    // -------------------------------------------------------------------------------------------------------------------------------------
    // Here we make a map from a name of the cross-section to the cross-section object itself. In this way, we can iterate through the
    // cross-section objects and reduce code-bloat. The first index is an identifier for the selection that's applied (generic or goldlen),
    // the second index is an identifier for the kinematic quantity that's relevant for the cross-section (e.g. muonMomentum), and the
    // mapped type is the cross-section object.
    std::map<std::string, std::map<std::string, CrossSectionHelper::CrossSection> > xsecMapSidebandData;
    std::map<std::string, std::map<std::string, CrossSectionHelper::CrossSection> > xsecMapSidebandNuWro;
    std::map<std::string, std::map<std::string, CrossSectionHelper::CrossSection> > xsecMapSidebandNuWroTrue;

    // We additionally make a map from each cross-section to the limits of the phase-space that we should consider. The key is the
    // identifier for the kinematic quantity and the mapped value is a pair containing the limits [min, max]
    std::map< std::string, std::pair<float, float> > phaseSpaceMap;

    // ATTN the configuration allows the user to enable or disable each cross-section. If a cross-section has been disabled, then it won't
    // be added the xSecMap. However, the phaseSpaceMap always includes all kinematic parameters!

    // Add the differential cross-sections
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
        // Add to the phase-space map
        phaseSpaceMap.emplace(name, std::pair<float, float>({binning.min, binning.max}));

        // Here we calculate every cross-section using both the generic and golden selection. In the end we only use the golden selection for
        // the pion momentum, but we additionally apply it to other cross-sections as a cross-check.
        // Add the cross-section object to the map using the binning from the input configuration
        const auto &[extendedBinEdges, hasUnderflow, hasOverflow] = CrossSectionHelper::GetExtendedBinEdges(binning.min, binning.max, binning.binEdges);
        for (const auto &selectionName : {"generic"})
        {
            // Don't setup a cross-section object if it's been disabled in the configuration
            if (config.extractXSecs.crossSectionIsEnabled.at(selectionName).at(name))
            {
                xsecMapSidebandData[selectionName].emplace(name, CrossSectionHelper::CrossSection(systParams, extendedBinEdges, hasUnderflow, hasOverflow, scaleByBinWidth));
                xsecMapSidebandNuWro[selectionName].emplace(name, CrossSectionHelper::CrossSection(systParams, extendedBinEdges, hasUnderflow, hasOverflow, scaleByBinWidth));
                xsecMapSidebandNuWroTrue[selectionName].emplace(name, CrossSectionHelper::CrossSection(systParams, extendedBinEdges, hasUnderflow, hasOverflow, scaleByBinWidth));
            }
        }
    }

    // ATTN here we use the machinary for a differential cross-section, and treat the total cross-section as a single-bin measurement.
    // The "kinematic quantity" in this case is just a dummy parameter. Here we define a single bin with edges arbitrarily chosen to be
    // (-1 -> +1), and we request that the cross-section object does not apply bin-width scaling. When we fill this object, we will use
    // the dummy kinematic quantity with a value of 0 for all events. This is arbitrary, as long as it's within the bin edges we chose.
    // In this way the single bin contains all events. It's just a trick to avoid implementing extra logic for the total cross-section.
    for (const auto &selectionName : {"generic"})
    {
        // Don't setup a cross-section object if it's been disabled in the configuration
        if (config.extractXSecs.crossSectionIsEnabled.at(selectionName).at("total"))
        {
            xsecMapSidebandData[selectionName].emplace("total", CrossSectionHelper::CrossSection(systParams, {-1.f, 1.f}, false, false, false));
            xsecMapSidebandNuWro[selectionName].emplace("total", CrossSectionHelper::CrossSection(systParams, {-1.f, 1.f}, false, false, false));
            xsecMapSidebandNuWroTrue[selectionName].emplace("total", CrossSectionHelper::CrossSection(systParams, {-1.f, 1.f}, false, false, false));
        }
    }

    // The dummy value that will be used as the "kinematic quantity" for the total cross-section
    const auto dummyValue = 0.f;

    // -------------------------------------------------------------------------------------------------------------------------------------
    // Setup the relevent "getters" for each the sideband
    // -------------------------------------------------------------------------------------------------------------------------------------
    // Here we define a map from the name of each cross-section to a function which pulls out the relevant kinematic quanitity from an input
    // analysis data object. Again this is done up-front to reduce code-bloat below.
    std::unordered_map< std::string, std::function<float(const AnalysisHelper::AnalysisData &)> > getSidebandValue;

    // Differential cross-section kinematic parameters
    getSidebandValue.emplace("muonCosTheta",  [](const auto &data) { return data.muonCosTheta;    });
    getSidebandValue.emplace("muonPhi",       [](const auto &data) { return data.muonPhi;         });
    getSidebandValue.emplace("muonMomentum",  [](const auto &data) { return data.muonMomentum;    });
    getSidebandValue.emplace("pionCosTheta",  [](const auto &data) { return data.protonCosTheta;  }); // Getting proton instead of pion values
    getSidebandValue.emplace("pionPhi",       [](const auto &data) { return data.protonPhi;       }); // Getting proton instead of pion values
    getSidebandValue.emplace("pionMomentum",  [](const auto &data) { return data.protonMomentum;  }); // Getting proton instead of pion values
    getSidebandValue.emplace("muonPionAngle", [](const auto &data) { return data.muonProtonAngle; }); // Getting proton instead of pion values
    getSidebandValue.emplace("nProtons",      [](const auto &data) { return data.nProtons-1;      }); // Leading proton treated as pion in CC0pi analysis
    // ATTN as described above, for the total cross-section we don't have an associated kinematic quantity so we just return a dummy value
    getSidebandValue.emplace("total", [=](const auto &) { return dummyValue; });


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
    float totalExposurePOT = 0;
    std::vector< std::tuple<AnalysisHelper::SampleType, std::string, std::string, float> > inputData;

    if(std::binary_search(config.global.runs.begin(), config.global.runs.end(), 1))
    {
        inputData.emplace_back(AnalysisHelper::Overlay, "", config.filesRun1.overlaysFileName, NormalisationHelper::GetOverlaysNormalisation(config, 1));
        inputData.emplace_back(AnalysisHelper::Dirt,    "", config.filesRun1.dirtFileName, NormalisationHelper::GetDirtNormalisation(config, 1));
        inputData.emplace_back(AnalysisHelper::DataEXT, "", config.filesRun1.dataEXTFileName, NormalisationHelper::GetDataEXTNormalisation(config, 1));
        inputData.emplace_back(AnalysisHelper::NuWro,   "", config.filesRun1.nuWroFileName, NormalisationHelper::GetNuWroNormalisation(config, 1));
        inputData.emplace_back(AnalysisHelper::DataBNB, "", config.filesRun1.dataBNBFileName, 1.f);
        // Add the detector variation files
        for (const auto &[name, fileName] : config.filesRun1.detVarFiles)
        {
            inputData.emplace_back(AnalysisHelper::DetectorVariation, name, fileName, NormalisationHelper::GetDetectorVariationNormalisation(config, name, 1));
        }
        totalExposurePOT += config.normsRun1.dataBNBTor875WCut / (1e20);
    }
    if(std::binary_search(config.global.runs.begin(), config.global.runs.end(), 2))
    {
        inputData.emplace_back(AnalysisHelper::Overlay, "", config.filesRun2.overlaysFileName, NormalisationHelper::GetOverlaysNormalisation(config, 2));
        inputData.emplace_back(AnalysisHelper::Dirt,    "", config.filesRun2.dirtFileName, NormalisationHelper::GetDirtNormalisation(config, 2));
        inputData.emplace_back(AnalysisHelper::DataEXT, "", config.filesRun2.dataEXTFileName, NormalisationHelper::GetDataEXTNormalisation(config, 2));
        inputData.emplace_back(AnalysisHelper::NuWro,   "", config.filesRun2.nuWroFileName, NormalisationHelper::GetNuWroNormalisation(config, 2));
        inputData.emplace_back(AnalysisHelper::DataBNB, "", config.filesRun2.dataBNBFileName, 1.f);
        // Add the detector variation files
        for (const auto &[name, fileName] : config.filesRun2.detVarFiles)
        {
            inputData.emplace_back(AnalysisHelper::DetectorVariation, name, fileName, NormalisationHelper::GetDetectorVariationNormalisation(config, name, 2));
        }
        totalExposurePOT += config.normsRun2.dataBNBTor875WCut / (1e20);
    }
    if(std::binary_search(config.global.runs.begin(), config.global.runs.end(), 3))
    {
        inputData.emplace_back(AnalysisHelper::Overlay, "", config.filesRun3.overlaysFileName, NormalisationHelper::GetOverlaysNormalisation(config, 3));
        inputData.emplace_back(AnalysisHelper::Dirt,    "", config.filesRun3.dirtFileName, NormalisationHelper::GetDirtNormalisation(config, 3));
        inputData.emplace_back(AnalysisHelper::DataEXT, "", config.filesRun3.dataEXTFileName, NormalisationHelper::GetDataEXTNormalisation(config, 3));
        inputData.emplace_back(AnalysisHelper::NuWro,   "", config.filesRun3.nuWroFileName, NormalisationHelper::GetNuWroNormalisation(config, 3));
        inputData.emplace_back(AnalysisHelper::DataBNB, "", config.filesRun3.dataBNBFileName, 1.f);

        // Add the detector variation files
        for (const auto &[name, fileName] : config.filesRun3.detVarFiles)
        {
            inputData.emplace_back(AnalysisHelper::DetectorVariation, name, fileName, NormalisationHelper::GetDetectorVariationNormalisation(config, name, 3));
        }
        totalExposurePOT += config.normsRun3.dataBNBTor875WCut / (1e20);
    }
    
    scalingData.exposurePOT = totalExposurePOT;
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

        if(isDetVar) continue; // Todo: Remove!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
        std::cout<<"\n##############\nOnly counting every 20th event!\n##############"<<std::endl;
        for (unsigned int i = 0; i < nEvents/20; ++i) //todo change back to every event!!!!!!!!!!!!!!!!!!!!!!
        {
            AnalysisHelper::PrintLoadingBar(i, nEvents);
            reader.LoadEvent(i);

            // -----------------------------------------------------------------------------------------------------------------------------
            // Work out if this event passed the selection and apply any additional phase-space cuts based on the input binning
            // -----------------------------------------------------------------------------------------------------------------------------

            // Run the selection
            const auto &[passedSidebandSelection, sidebandCutsPassed, sidebandAssignedPdgCodes] = sidebandSelection.Execute(pEvent);

            // Get the reco analysis data (if available, otherwise set to dummy values)
            const auto recoData = (
                passedSidebandSelection
                    ? AnalysisHelper::GetRecoAnalysisDataCC0Pi(pEvent->reco, sidebandAssignedPdgCodes)
                    : AnalysisHelper::GetDummyAnalysisData()
            );

            // Here we apply reco-level phase-space restrictions
            // For any event that passes the generic selection, get the value of the kinematic quantity and check if it is outside of the
            // min/max values supplied in the binning. If so, then reject the event.
            bool passesPhaseSpaceReco = false;
            if (passedSidebandSelection)
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

            const auto isSelectedSideband = passedSidebandSelection && passesPhaseSpaceReco;
            // std::map<std::string, bool> isSelectedMap = {{"sideband",isSelectedSideband}};
            std::map<std::string, bool> isSelectedMap = {{"generic",isSelectedSideband},{"golden",isSelectedSideband}};

            // -----------------------------------------------------------------------------------------------------------------------------
            // Handle BNB data
            // -----------------------------------------------------------------------------------------------------------------------------
            // std::cout<<"Debug {{-1}}"<<std::endl;
            if (isDataBNB)
            {
                for (auto &[selectionName, xsecs] : xsecMapSidebandData)
                {
                    // Determine if we passed the relevant selection
                    const auto isSelected = isSelectedMap.at(selectionName);
                    // Only count events passing the selection
                    if (!isSelected) continue;

                    for (auto &[name, xsec] : xsecs)
                        xsec.AddSelectedBNBDataEvent(getSidebandValue.at(name)(recoData));
                }

                // For BNB data that's all we need to do!
                continue;
            }

            // Get the nominal event weight, scaled by the sample normalisation
            const auto weight = AnalysisHelper::GetNominalEventWeight(pEvent) * normalisation;

            if (isNuWro)
            {
                for (auto &[selectionName, xsecs] : xsecMapSidebandNuWro)
                {
                    // Determine if we passed the relevant selection
                    const auto isSelected = isSelectedMap.at(selectionName);
                    // Only count events passing the selection
                    if (!isSelected) continue;

                    for (auto &[name, xsec] : xsecs)
                        xsec.AddWeightedSelectedBNBDataEvent(getSidebandValue.at(name)(recoData), weight);
                }

                // For BNB data that's all we need to do!
                // continue;
            }


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
            // std::map<std::string, bool> isSignalMap = {{"sideband",isCC0PiSignal}};
            std::map<std::string, bool> isSignalMap = {{"generic",isCC0PiSignal},{"golden",isCC0PiSignal}};

            // std::cout<<"Debug {{-0.9}}"<<std::endl;
            // -----------------------------------------------------------------------------------------------------------------------------
            // Handle the detector variation samples (unisims)
            // -----------------------------------------------------------------------------------------------------------------------------
            if (isDetVar)
            {
                // std::cout<<"Debug {{-0.8}}"<<std::endl;
                for (auto &[selectionName, xsecs] : xsecMapSidebandData)
                {
                    const auto isSignal = isSignalMap.at(selectionName);
                    // Handle signal events
                    if (isSignal)
                    {
                        // std::cout<<"Debug {{-0.7}}"<<std::endl;
                        // Determine if we passed the relevant selection
                        const auto isSelected = isSelectedMap.at(selectionName);
                        // std::cout<<"Debug {{-0.64}}"<<std::endl;
                        for (auto &[name, xsec] : xsecs)
                        {
                            // std::cout<<"Debug {{-0.63}}"<<std::endl;
                            const auto recoValue = getSidebandValue.at(name)(recoData);
                            const auto trueValue = getSidebandValue.at(name)(truthData);
                            // std::cout<<"Debug {{-0.62}}"<<std::endl;
                            xsec.AddSignalEventDetVar(recoValue, trueValue, isSelected, weight, sampleName);
                        }
                        // std::cout<<"Debug {{-0.61}}"<<std::endl;
                    }
                    // Handle selected background events
                    else
                    {
                        // std::cout<<"Debug {{-0.6}}"<<std::endl;
                        for (auto &[selectionName, xsecs] : xsecMapSidebandData)
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
                        // std::cout<<"Debug {{-0.5}}"<<std::endl;
                    }
                }

                // For detector variation samples, that's all we need to do!
                // std::cout<<"Debug {{-0.4}}"<<std::endl;
                continue;
            }
            // -----------------------------------------------------------------------------------------------------------------------------
            // Handle all other events (i.e those from the nominal simulation): Overlays, dirt, EXT data
            // -----------------------------------------------------------------------------------------------------------------------------
            // std::cout<<"Debug {{0}}"<<std::endl;
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
            if(xsecWeightsScaleFactor<0.001f) std::cout<<"xsecWeightsScaleFactor: "<<xsecWeightsScaleFactor<<std::endl;            
            xsecWeightsScaleFactor = std::max(xsecWeightsScaleFactor, 0.0001f); //Todo: use continue instead

            const auto xsecWeights = (
                isOverlay
                    ? CrossSectionHelper::ScaleWeightsMap(CrossSectionHelper::GetWeightsMap(pEvent->truth, systParams.xsecDimensions, config.extractXSecs.mutuallyExclusiveDimensions), xsecWeightsScaleFactor)
                    : CrossSectionHelper::GetUnitWeightsMap(systParams.xsecDimensions)
            );

            // Get the reinteraction weights
            const auto reintWeights = (
                isOverlay
                    ? CrossSectionHelper::GetWeightsMap(pEvent->truth, systParams.reintDimensions, config.extractXSecs.mutuallyExclusiveDimensions)
                    : CrossSectionHelper::GetUnitWeightsMap(systParams.reintDimensions)
            );
            
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
                        const auto seedString =  selectionName + name + std::to_string(i);
                        if(isNuWro)
                        {
                            xsecMapSidebandNuWroTrue.at(selectionName).at(name).AddSignalEvent(recoValue, trueValue, isSelected, weight, fluxWeights, xsecWeights, reintWeights, seedString);
                        } else
                        {
                            if(isOverlay) // Genie MC for NuWro (no EXT or dirt events in NuWro)
                            {
                                xsec.AddSignalEvent(recoValue, trueValue, isSelected, weight, fluxWeights, xsecWeights, reintWeights, seedString);
                            }
                            if(isOverlay || isDataEXT || isDirt) // Genie MC for Data
                            {
                                xsecMapSidebandData.at(selectionName).at(name).AddSignalEvent(recoValue, trueValue, isSelected, weight, fluxWeights, xsecWeights, reintWeights, seedString);
                            }
                        }
                    }
                }
                // HandleisSelectedd background events
                else if (isSelected)
                {
                    for (auto &[name, xsec] : xsecs)
                    {
                        const auto recoValue = getSidebandValue.at(name)(recoData);
                        const auto seedString =  selectionName + name + std::to_string(i);
                        std::vector<float> bootstrapWeights; // Parameter only needed for CC1pi
                        if(isNuWro)
                        {
                            xsecMapSidebandNuWroTrue.at(selectionName).at(name).AddSelectedBackgroundEvent(recoValue, false, weight, fluxWeights, xsecWeights, reintWeights, bootstrapWeights, seedString);
                        }else
                        {
                            if(isOverlay) // Genie MC for NuWro (no EXT or dirt events in NuWro)
                            {
                                xsec.AddSelectedBackgroundEvent(recoValue, false, weight, fluxWeights, xsecWeights, reintWeights, bootstrapWeights, seedString);
                            }
                            if(isOverlay || isDataEXT || isDirt) // Genie MC for Data
                            {
                                xsecMapSidebandData.at(selectionName).at(name).AddSelectedBackgroundEvent(recoValue, isDirt, weight, fluxWeights, xsecWeights, reintWeights, bootstrapWeights, seedString);
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
    //Parameters: selectionName name
    std::map<std::string,std::map<std::string, std::vector<float>>> cc0piCovarianceMapData, cc0piCovarianceMapNuWro; 
    std::map<std::string, std::map<std::string, paramAndErrorPair>> cc0piNominalConstraintMapData, cc0piNominalConstraintMapNuWro;
    //Parameters: selectionName name paramName (i.e. golden muonMomentum hadronProduction)
    std::map<std::string, std::map<std::string, std::map<std::string, std::vector<paramAndErrorPair>>>> cc0piUniverseConstraintMapData, cc0piUniverseConstraintMapNuWro;
    
    for (const auto &[dataTypeName, xsecMap] : {std::make_pair(std::string("NuWro"), xsecMapSidebandNuWro), std::make_pair(std::string("Data"), xsecMapSidebandData)})
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


                for(unsigned int r = 0; r<signalData.GetRows(); r++) // DEBUG - TODO: REMOVE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                {
                    for(unsigned int c = 0; c<signalData.GetColumns(); c++) // DEBUG - TODO: REMOVE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    {
                        if(signalData.At(r, c) < 0)
                        {
                            std::cout<<"signalData.At("<<r<<","<<c<<") = "<<signalData.At(r, c)<<std::endl;
                            signalData.SetElement(r, c, std::max(signalData.At(r, c), 0.f));
                        }
                    }
                }

                // // std::cout<<"Debug >>-2<<"<<std::endl;
                // std::vector<float> elements;
                // const auto nBins = selectedEventsData.GetRows();
                // // std::cout<<"_______________________Fitting Point 1.1"<<std::endl;
                // for (unsigned int iBin = 0; iBin < selectedEventsData.GetRows(); ++iBin)
                // {
                //     const auto value = selectedEventsData.At(iBin, 0);
                //     if (value<0)
                //     {
                //         std::cout<<"ERROR: ExtractNuWroSidebandFit2 - Background-removed signal data is negative."<<std::endl;
                //         // throw std::logic_error("ERROR: ExtractNuWroSidebandFit2 - Background-removed signal data is negative."); // TODO: Uncomment
                //     }
                //     elements.push_back(AnalysisHelper::GetCountUncertainty(std::max(value,0.f))); //TODO: Check if this is correct
                // }
                
                ///////////////////////////////////////////
                //////// Start: Get the totalErrorMatrixMap
                // // std::cout<<"Debug >>-1.2<<"<<std::endl;
                // const auto &[pPredictionStatBias, pPredictionStatCovariance] = xsec.GetPredictedCrossSectionStatUncertainty(scalingData);
                // const auto predictionErrorMatrix = CrossSectionHelper::GetErrorMatrix(*pPredictionStatBias, *pPredictionStatCovariance);
                // // std::cout<<"Debug >>-1.1<<"<<std::endl;

                // Build the total error matrix for each group of systematic parameters
                // Here errorMatrixMap is indexed by [quantity][group], where quantity = "data" or "smearingMatrix"
                // Here totalErrorMatrixMap is index by [quantity] and contains the sum (over all groups) of the entries in errorMatrixMap
                // std::map<std::string, std::map<std::string, ubsmear::UBMatrix> > errorMatrixMap;
                std::map<std::string, ubsmear::UBMatrix> totalErrorMatrixMap;

                // std::cout<<"Debug >>-1<<"<<std::endl;
                const auto dataSystBiasCovariances = xsec.GetBNBDataCrossSectionSystUncertainties(scalingData);//todo remove "cross-section"? 
                // const auto smearingMatrixSystBiasCovariances = xsec.GetSmearingMatrixSystUncertainties();//todo remove "cross-section"?
                for (const auto &[quantity, covariances] :
                {
                    std::make_pair(std::string("data"), dataSystBiasCovariances), 
                    // std::make_pair(std::string("smearingMatrix"), smearingMatrixSystBiasCovariances)
                })
                {
                    // Get the number of bins (for the smearing matrix there are N^2 bins when flattened)
                    const auto nBins = (quantity == "data" ? selectedEventsData.GetRows() : std::pow(smearingMatrix.GetRows(), 2));

                    auto errorMatrixTotalSum =  ubsmear::UBMatrixHelper::GetZeroMatrix(nBins, nBins);
                    const auto statUncertainties = (quantity == "data") ? xsec.GetBNBDataCrossSectionStatUncertainty(scalingData) : ubsmear::UBMatrixHelper::GetZeroMatrix(nBins, 1);
                    const auto statVariances = ubsmear::ElementWiseOperation(statUncertainties, statUncertainties, [](const auto &l, const auto &r) { return l * r; });
                    const auto statErrorMatrix = ubsmear::UBMatrixHelper::GetDiagonalMatrix(statVariances);
                    errorMatrixTotalSum = errorMatrixTotalSum + statErrorMatrix;

                    // Handle the multisim parameters        
                    for (const auto &[group, map] : covariances)
                    {
                        // Setup an empty error matrix for this group
                        auto errorMatrixTotal =  ubsmear::UBMatrixHelper::GetZeroMatrix(nBins, nBins);
                        // Loop over the parameters in this group
                        for (const auto &[paramName, biasCovariance] : map)
                        {
                            std::cout<<"DEBUG --- "<<group<<" - "<<paramName<<std::endl;
                            // Get the bias vector and covariance matrix
                            const auto &[pBias, pCovariance] = biasCovariance;

                            // std::cout<<"Debug >>-0.3<<"<<std::endl;
                            // Get the total error matrix from the bias and covariance
                            const auto errorMatrix = CrossSectionHelper::GetErrorMatrix(*pBias, *pCovariance);

                            // std::cout<<"Debug >>-0.2<<      error matrix total:"<<std::endl;
                            // Add this error matrix to the total

                            std::cout<<errorMatrixTotal.GetRows()<<"   "<<errorMatrixTotal.GetColumns()<<std::endl;
                            // std::cout<<"Debug >>-0.11<<      error matrix:"<<std::endl;
                            
                            std::cout<<errorMatrix.GetRows()<<"   "<<errorMatrix.GetColumns()<<std::endl;
                            errorMatrixTotal = errorMatrixTotal + errorMatrix;
                        }
                        // std::cout<<"Debug >>-0.1<<"<<std::endl;

                        // // Store this total in the map
                        // errorMatrixMap[quantity].emplace(group, errorMatrixTotal);

                        // Add this total to the grand total error matrix
                        errorMatrixTotalSum = errorMatrixTotalSum + errorMatrixTotal;
                    }
                    // std::cout<<"Debug >>0<<"<<std::endl;

                    // Add the grand total to the map
                    totalErrorMatrixMap.emplace(quantity, errorMatrixTotalSum);
                }
                //////// End: Get the totalErrorMatrixMap

                // std::cout<<"Debug >>1<<"<<std::endl;
                pMetadataNuWro = &metadata;
                pSmearingMatrixNuWro = &smearingMatrix;
                pPredictionNuWro = &selectedEventsSignalTruth;
                pDataNuWro = &signalData;
                pTotalErrorMatrixMap = &totalErrorMatrixMap; // todo check whether individual universe uncertainties need to be generated (due to differenced is genie backgrounds in each universe)
                // pTotalErrorMatrixMapDataNuWro = &totalErrorMatrixMap.at("data");
                // pTotalErrorMatrixMapSmearingMatrixNuWro = &totalErrorMatrixMap.at("smearingMatrix");
                // pPredictionErrorMatrixNuWro = &predictionErrorMatrix;

                // std::cout<<"Debug >>2<<"<<std::endl;
                const auto nBins = pDataNuWro->GetRows();
                auto minimizer = FittingHelper(nBins);
                // std::cout<<"Debug >>2.1<<"<<std::endl;
                std::pair<std::vector<Double_t>, std::vector<Double_t>> result;
                
                std::vector<float> fitCovMatrixVector;
                bool successful = false;
                minimizer.Fit(GetChi2NuWro, result, successful, fitCovMatrixVector, 0);
                // std::cout<<"Debug >>3<<"<<std::endl;
                if(!successful)
                {
                    std::cout<<"ERROR: ExtractNuWroSidebandFit2 - Fit failed."<<std::endl;
                    throw std::logic_error("ERROR: ExtractNuWroSidebandFit2 - Nominal fit failed.");
                }

                // if(dataTypeName=="NuWro") cc0piCovarianceMapNuWro[selectionName].emplace(name, fitCovMatrixVector);
                // else if (dataTypeName=="Data") cc0piCovarianceMapData[selectionName].emplace(name, fitCovMatrixVector);
                // else throw std::logic_error("ERROR: Unknown dataTypeName");

                // std::cout<<"Debug >>4<<"<<std::endl;
                const ubsmear::UBMatrix sidebandCovMatrix(fitCovMatrixVector, nBins, nBins);
                FormattingHelper::SaveMatrix(sidebandCovMatrix, dataTypeName + "_SidebandFit_" + selectionName + "_" + name + "_sideband_stat_covariance.txt");
                vector<float> paramVector(result.first.begin(), result.first.end()); // Todo avoid this
                vector<float> paramErrorVector(result.second.begin(), result.second.end()); // Todo avoid this
                const ubsmear::UBMatrix sidebandParamVectorTruth(paramVector, nBins, 1);
                const ubsmear::UBMatrix sidebandErrorVectorTruth(paramErrorVector, nBins, 1);
                FormattingHelper::SaveMatrix(sidebandParamVectorTruth, dataTypeName + "_SidebandFit_" + selectionName + "_" + name + "_sideband_parameterVector.txt");
                FormattingHelper::SaveMatrix(sidebandErrorVectorTruth, dataTypeName + "_SidebandFit_" + selectionName + "_" + name + "_sideband_parameterErrorVector.txt");

                // std::cout<<"Debug >>5<<"<<std::endl;
                if(dataTypeName=="NuWro")
                {
                    cc0piNominalConstraintMapNuWro[selectionName].emplace(name, result);
                    cc0piCovarianceMapNuWro[selectionName].emplace(name, fitCovMatrixVector);
                }
                else if (dataTypeName=="Data")
                {
                    cc0piNominalConstraintMapData[selectionName].emplace(name, result);
                    cc0piCovarianceMapData[selectionName].emplace(name, fitCovMatrixVector);
                }
                else throw std::logic_error("ERROR: Unknown dataTypeName");

                // std::cout<<"Debug >>6<<"<<std::endl;
                continue; // todo REMOVE; Just for debugging!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                // const auto sidebandWeights = xsec.GetSidebandWeights(scalingData);

                // -------------------------------------------------------------------------------------------------------------------------------------
                // Fit each universe
                // -------------------------------------------------------------------------------------------------------------------------------------
                const auto weightDimensions = {std::make_pair("xsec", systParams.xsecDimensions), std::make_pair("reint", systParams.reintDimensions), std::make_pair("flux", systParams.fluxDimensions)};
                for (const auto &[group, dimensions] : weightDimensions)
                {
                    for (const auto &[paramName, nUniverses] : dimensions)
                    {
                        std::vector<paramAndErrorPair> resultVector;
                        for (unsigned int iUni = 0; iUni < nUniverses; ++iUni)
                        {
                            std::cout<<"++++++++++++++ExtractNuWroSidebandFit2 Fitting: "<<selectionName<<" - "<<name<<" - universe: "<<iUni<<" - group: "<<group<<" - paramName: "<<paramName<<"++++++++++++++"<<std::endl;
                            const auto selectedSignalTruthUniverses = xsec.GetSelectedSignalRecoTruthMap().at(group).at(paramName).at(iUni);
                            const auto selectedBackgroundRecoUniverses = xsec.GetSelectedBackgroundRecoMap().at(group).at(paramName).at(iUni);
                            std::pair<std::vector<Double_t>, std::vector<Double_t>> result;
                            const auto pSmearingMatrixInUniverse = xsec.GetSmearingMatrixInUniverseAllSelected(group, paramName, iUni);//xsec.GetSmearingMatrixAllSelected();
                            if(pSmearingMatrixInUniverse)
                            {
                                auto smearingMatrixInUniverse = *pSmearingMatrixInUniverse;

                                // AnalysisHelper::PrintLoadingBar(iUni, nUniverses);
                                const auto selectedSignalTruthInUniverse = xsec.GetSignalSelectedTrue(selectedSignalTruthUniverses);
                                const auto selectedBackgoundRecoInUniverse = CrossSectionHelper::GetMatrixFromHist(selectedBackgroundRecoUniverses);
                                auto signalDataInUniverse = selectedEventsData - selectedBackgoundRecoInUniverse;

                                for(unsigned int r = 0; r<signalDataInUniverse.GetRows(); r++) // DEBUG - TODO: REMOVE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                                {
                                    for(unsigned int c = 0; c<signalDataInUniverse.GetColumns(); c++) // DEBUG - TODO: REMOVE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                                    {
                                        if(signalDataInUniverse.At(r, c) < 0)
                                        {
                                            std::cout<<"universe signalData.At("<<r<<","<<c<<") = "<<signalDataInUniverse.At(r, c)<<std::endl;
                                            signalDataInUniverse.SetElement(r, c, std::max(signalDataInUniverse.At(r, c), 0.f));
                                        }
                                    }
                                }

                                pPredictionNuWro = &selectedSignalTruthInUniverse;
                                pDataNuWro = &signalDataInUniverse;
                                pSmearingMatrixNuWro = &smearingMatrixInUniverse;

                                std::vector<float> covMatrixInUniverse;
                                bool successful = false;
                                minimizer.Fit(GetChi2NuWro, result, successful, covMatrixInUniverse, 0);
                                if(!successful)
                                {
                                    std::cout<<"Fit not possible - Minimization failed."<<std::endl;
                                    result = std::make_pair(std::vector<Double_t>(nBins, -1), std::vector<Double_t>(nBins, -1));
                                }
                            }
                            else
                            {
                                std::cout<<"Fit not possible - Smearing matrix not available."<<std::endl;
                                result = std::make_pair(std::vector<Double_t>(nBins, -1), std::vector<Double_t>(nBins, -1));
                            }
                            resultVector.push_back(result);
                        }

                        if(dataTypeName=="NuWro") cc0piUniverseConstraintMapNuWro[selectionName][name].emplace(paramName, resultVector);
                        else if (dataTypeName=="Data") cc0piUniverseConstraintMapData[selectionName][name].emplace(paramName, resultVector);
                        else throw std::logic_error("ERROR: Unknown dataTypeName");
                    }
                }
            }
        }
    }

    std::ofstream ofs1("cc0piCovarianceMapNuWro.bin", std::ios::binary);
    std::ofstream ofs2("cc0piNominalConstraintMapNuWro.bin", std::ios::binary);
    std::ofstream ofs3("cc0piUniverseConstraintMapNuWro.bin", std::ios::binary);
    std::ofstream ofs4("cc0piCovarianceMapData.bin", std::ios::binary);
    std::ofstream ofs5("cc0piNominalConstraintMapData.bin", std::ios::binary);
    std::ofstream ofs6("cc0piUniverseConstraintMapData.bin", std::ios::binary);
    

    boost::archive::binary_oarchive oarch1(ofs1);
    boost::archive::binary_oarchive oarch2(ofs2);
    boost::archive::binary_oarchive oarch3(ofs3);
    boost::archive::binary_oarchive oarch4(ofs4);
    boost::archive::binary_oarchive oarch5(ofs5);
    boost::archive::binary_oarchive oarch6(ofs6);

    oarch1 << cc0piCovarianceMapNuWro;
    oarch2 << cc0piNominalConstraintMapNuWro;
    oarch3 << cc0piUniverseConstraintMapNuWro;
    oarch4 << cc0piCovarianceMapData;
    oarch5 << cc0piNominalConstraintMapData;
    oarch6 << cc0piUniverseConstraintMapData;

    ofs1.close();
    ofs2.close();
    ofs3.close();
    ofs4.close();
    ofs5.close();
    ofs6.close();

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

            FormattingHelper::SaveMatrix(eventsSignalTruth, "NuWroSidebandFitTruth_" + selectionName + "_" + name + "_eventsSignalTruth.txt");
            FormattingHelper::SaveMatrix(selectedEventsSignalTruth, "NuWroSidebandFitTruth_" + selectionName + "_" + name + "_selectedEventsSignalTruth.txt");
            FormattingHelper::SaveMatrix(smearingMatrix, "NuWroSidebandFitTruth_" + selectionName + "_" + name + "_smearingMatrix.txt");
            FormattingHelper::SaveMatrix(smearingMatrixGeneral, "NuWroSidebandFitTruth_" + selectionName + "_" + name + "_smearingMatrixGeneral.txt");
        }
    }

    std::cout<<"------------- All done -------------"<<std::endl;
    return;
}

} // namespace ubcc1pi_macros
