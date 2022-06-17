/**
 *  @file  ubcc1pi_standalone/Macros/ExtractNuWroSidebandFit.cxx
 *
 *  @brief The implementation file of the ExtractNuWroSidebandFit macro
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

std::vector<float> xNuWro, yNuWro, errorYNuWro, SNuWro;
bool underflowNuWro, overflowNuWro;

//______________________________________________________________________________

void fcnNuWro(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
    //calculate chi2
    if(xNuWro.size()!=yNuWro.size() || yNuWro.size()!=errorYNuWro.size() || xNuWro.size()*xNuWro.size()!=SNuWro.size())
        throw std::logic_error("Fitting function for ExtractNuWroSidebandFit - Incompatible input dimenstions.");

    Double_t chisq = 0;
    Double_t delta;
    Int_t nBins = xNuWro.size();
    
    std::vector<Double_t> xScaled(xNuWro.begin(), xNuWro.end()); 
    for (Int_t i=0; i<nBins; i++)
    {
        xScaled[i]*=par[i];
    }
    
    std::vector<Double_t> xSmeared(nBins, 0);

    for (Int_t i=0; i<nBins; i++)
    {
        for (Int_t j=0; j<nBins; j++)
        {
            xSmeared[i] +=  SNuWro[i*nBins+j]*xScaled[j];
        }
    }
    // const auto xSmeared = S*xScaled;

    Int_t relevantBinsMax = overflowNuWro ? nBins-1 : nBins;
    Int_t relevantBinsMin = underflowNuWro ? 1 : 0;

    for (Int_t i=relevantBinsMin; i<relevantBinsMax; i++) 
    {
        delta  = (yNuWro[i]-xSmeared[i])/errorYNuWro[i];
        chisq += delta*delta;
    }

    f = chisq;
}

namespace ubcc1pi_macros
{

void ExtractNuWroSidebandFit(const Config &config)
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
    // systParams.detVarDimensions = config.extractXSecs.detVarDimensions;
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
            std::cout<<"selectionName: "<<selectionName<<" - name: "<<name<<" - extendedBinEdges: ";
            for(const auto &edge: extendedBinEdges) std::cout<<edge<<" ";
            std::cout<<"\n - hasUnderflow"<<hasUnderflow<<" - hasOverflow"<<hasOverflow<<std::endl;

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
            // xsecMap[selectionName].emplace("total", CrossSectionHelper::CrossSection(systParams, {-1.f, 1.f}, false, false, false));
            xsecMapSidebandData[selectionName].emplace("total", CrossSectionHelper::CrossSection(systParams, {-1.f, 1.f}, false, false, false));
            xsecMapSidebandNuWro[selectionName].emplace("total", CrossSectionHelper::CrossSection(systParams, {-1.f, 1.f}, false, false, false));
            xsecMapSidebandNuWroTrue[selectionName].emplace("total", CrossSectionHelper::CrossSection(systParams, {-1.f, 1.f}, false, false, false));
            // xsecMapSideband2[selectionName].emplace("total", CrossSectionHelper::CrossSection(systParams, {-1.f, 1.f}, false, false, false));
        }
    }
    // Don't setup a cross-section object if it's been disabled in the configuration
    // if (config.extractXSecs.crossSectionIsEnabled.at("sideband").at("total"))
    //     xsecMapSideband["sideband"].emplace("total", CrossSectionHelper::CrossSection(systParams, {-1.f, 1.f}, false, false, false));

    // The dummy value that will be used as the "kinematic quantity" for the total cross-section
    const auto dummyValue = 0.f;

    // // Check to see if we have any cross-sections enabled
    // if (xsecMap.empty())
    // {
    //     std::cout << "All cross-sections have been disabled in the configuration! Nothing more to do" << std::endl;
    //     return;
    // }

    // // Print the names of the cross-sections we are going to extract
    // std::cout << "The following cross-sections are enabled:" << std::endl;
    // for (const auto &[selectionName, xsecs] : xsecMap)
    // {
    //     std::cout << "  - Selection: " << selectionName << std::endl;
    //     for (const auto &entry : xsecs)
    //     {
    //         const auto &name = entry.first;
    //         std::cout << "    - " << name << std::endl;
    //     }
    // }

    // -------------------------------------------------------------------------------------------------------------------------------------
    // Setup the relevent "getters" for each the sideband
    // -------------------------------------------------------------------------------------------------------------------------------------
    // Here we define a map from the name of each cross-section to a function which pulls out the relevant kinematic quanitity from an input
    // analysis data object. Again this is done up-front to reduce code-bloat below.
    std::unordered_map< std::string, std::function<float(const AnalysisHelper::AnalysisData &)> > getSidebandValue;

    // Differential cross-section kinematic parameters
    getSidebandValue.emplace("muonCosTheta",  [](const auto &data) { return data.muonCosTheta;  });
    getSidebandValue.emplace("muonPhi",       [](const auto &data) { return data.muonPhi;       });
    getSidebandValue.emplace("muonMomentum",  [](const auto &data) { return data.muonMomentum;  });

    getSidebandValue.emplace("pionCosTheta",  [](const auto &data) { return data.protonCosTheta;  }); // Getting proton instead of pion values
    getSidebandValue.emplace("pionPhi",       [](const auto &data) { return data.protonPhi;       }); // Getting proton instead of pion values
    getSidebandValue.emplace("pionMomentum",  [](const auto &data) { return data.protonMomentum;  }); // Getting proton instead of pion values

    getSidebandValue.emplace("muonPionAngle", [](const auto &data) { return data.muonProtonAngle; }); // Getting proton instead of pion values
    getSidebandValue.emplace("nProtons",      [](const auto &data) { return data.nProtons-1;      }); //Leading proton treated as pion in CC0pi analysis

    // ATTN as described above, for the total cross-section we don't have an associated kinematic quantity so we just return a dummy value
    getSidebandValue.emplace("total", [=](const auto &) { return dummyValue; });


    const auto fluxHistNames = CrossSectionHelper::GetNominalFluxHistNames(config.flux.nuPdgsSignal, config.flux.nuPdgToHistName, config.flux.nomHistPattern);
    const auto &[fluxBinEdges, fluxValues] = CrossSectionHelper::ReadNominalFlux(config.flux.fileName, fluxHistNames, config.flux.pot);
    scalingData.pFluxReweightor = std::make_shared<CrossSectionHelper::FluxReweightor>(fluxBinEdges, fluxValues, systParams.fluxDimensions);
    float totalExposurePOT = 0;
    std::vector< std::tuple<AnalysisHelper::SampleType, std::string, std::string, float> > inputData;
    for (const auto run: config.global.runs)
    {
        // -------------------------------------------------------------------------------------------------------------------------------------
        // Setup the input files
        // -------------------------------------------------------------------------------------------------------------------------------------
        // Here we define a vector of tuples with 4 entries
        //   - First, the sample type (e.g. overlay)
        //   - Second, a string which is used to identify a given detector variation sample (for other sample type, this is unused)
        //   - Third, the path to the input file
        //   - Fourth, the normalisation factor to apply to all events in that file
        if(run == 1)
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
        else if(run == 2)
        {
            inputData.emplace_back(AnalysisHelper::Overlay, "", config.filesRun2.overlaysFileName, NormalisationHelper::GetOverlaysNormalisation(config, 2));
            inputData.emplace_back(AnalysisHelper::Dirt,    "", config.filesRun2.dirtFileName, NormalisationHelper::GetDirtNormalisation(config, 2));
            inputData.emplace_back(AnalysisHelper::DataEXT, "", config.filesRun2.dataEXTFileName, NormalisationHelper::GetDataEXTNormalisation(config, 2));
            inputData.emplace_back(AnalysisHelper::NuWro,   "", config.filesRun2.nuWroFileName, NormalisationHelper::GetNuWroNormalisation(config, 2));
            inputData.emplace_back(AnalysisHelper::DataBNB, "", config.filesRun2.dataBNBFileName, 2.f);
            // Add the detector variation files
            for (const auto &[name, fileName] : config.filesRun2.detVarFiles)
            {
                inputData.emplace_back(AnalysisHelper::DetectorVariation, name, fileName, NormalisationHelper::GetDetectorVariationNormalisation(config, name, 2));
            }
            totalExposurePOT += config.normsRun2.dataBNBTor875WCut / (1e20);
        }
        else if(run == 3)
        {
            inputData.emplace_back(AnalysisHelper::Overlay, "", config.filesRun3.overlaysFileName, NormalisationHelper::GetOverlaysNormalisation(config, 3));
            inputData.emplace_back(AnalysisHelper::Dirt,    "", config.filesRun3.dirtFileName, NormalisationHelper::GetDirtNormalisation(config, 3));
            inputData.emplace_back(AnalysisHelper::DataEXT, "", config.filesRun3.dataEXTFileName, NormalisationHelper::GetDataEXTNormalisation(config, 3));
            inputData.emplace_back(AnalysisHelper::NuWro,   "", config.filesRun3.nuWroFileName, NormalisationHelper::GetNuWroNormalisation(config, 1));
            inputData.emplace_back(AnalysisHelper::DataBNB, "", config.filesRun3.dataBNBFileName, 3.f);

            // Add the detector variation files
            for (const auto &[name, fileName] : config.filesRun3.detVarFiles)
            {
                inputData.emplace_back(AnalysisHelper::DetectorVariation, name, fileName, NormalisationHelper::GetDetectorVariationNormalisation(config, name, 3));
            }
            totalExposurePOT += config.normsRun3.dataBNBTor875WCut / (1e20);
        }
        else throw std::logic_error("ExtractNuWroSidebandFit - Invalid run number");
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
        const auto isDetVar  = (sampleType == AnalysisHelper::DetectorVariation);
        const auto isDataBNB = (sampleType == AnalysisHelper::DataBNB);
        const auto isDataEXT = (sampleType == AnalysisHelper::DataEXT);

        if(isDetVar) continue; //Todo Remove!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        // Open the input file for reading and enable the branches with systematic event weights (if required)
        FileReader reader(fileName);

        if (isOverlay || isNuWro)
            reader.EnableSystematicBranches();

        auto pEvent = reader.GetBoundEventAddress();

        // Loop over the events in the file
        const auto nEvents = reader.GetNumberOfEvents();
        for (unsigned int i = 0; i < nEvents; ++i)
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
                        // std::cout << "Event failed reco phase-space cuts" << std::endl;
                        break;
                    }
                }
            }

            const auto isSelectedSideband = passedSidebandSelection && passesPhaseSpaceReco;
            // std::map<std::string, bool> isSelectedMap = {{"sideband",isSelectedSideband}};
            std::map<std::string, bool> isSelectedMap = {{"generic",isSelectedSideband},{"golden",isSelectedSideband}};

            // -----------------------------------------------------------------------------------------------------------------------------
            // Work out if this event is signal, and apply any phase-space restrictions based on the input binning
            // -----------------------------------------------------------------------------------------------------------------------------

            // -----------------------------------------------------------------------------------------------------------------------------
            // Handle BNB data
            // -----------------------------------------------------------------------------------------------------------------------------
            
            if (isDataBNB)
            {
                for (auto &[selectionName, xsecs] : xsecMapSidebandData)
                {
                    // Determine if we passed the relevant selection
                    const auto isSelected = isSelectedMap.at(selectionName);
                    // Only count events passing the selection
                    if (!isSelected)
                        continue;

                    for (auto &[name, xsec] : xsecs)
                    {
                        xsec.AddSelectedBNBDataEvent(getSidebandValue.at(name)(recoData));
                    }
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
                    if (!isSelected)
                        continue;

                    for (auto &[name, xsec] : xsecs)
                    {
                        xsec.AddWeightedSelectedBNBDataEvent(getSidebandValue.at(name)(recoData), weight);
                    }
                }

                // For BNB data that's all we need to do!
                // continue;
            }


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
                        // std::cout << "Event failed truth phase-space cuts - name: "<<name<<" - min: "<<min<<" - value: "<<value<<" - max: "<<max << std::endl;                        
                        passesPhaseSpaceTruth = false;
                        break;
                    }
                }
            }

            const auto isCC0PiSignal = isTrueCC0Pi && passesPhaseSpaceTruth;
            // std::map<std::string, bool> isSignalMap = {{"sideband",isCC0PiSignal}};
            std::map<std::string, bool> isSignalMap = {{"generic",isCC0PiSignal},{"golden",isCC0PiSignal}};

            // -----------------------------------------------------------------------------------------------------------------------------
            // Handle the detector variation samples (unisims)
            // -----------------------------------------------------------------------------------------------------------------------------
            if (isDetVar)
            {
                for (auto &[selectionName, xsecs] : xsecMapSidebandData)
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
                            // std::cout<<"!isOverlay - name"<<name<<" - weight: "<<weight<<" - isSelected: "<< isSelected<<" - recoValue: "<<recoValue<<" - trueValue: "<<trueValue<<std::endl;
                            xsecMapSidebandNuWroTrue.at(selectionName).at(name).AddSignalEvent(recoValue, trueValue, isSelected, weight, fluxWeights, xsecWeights, reintWeights, seedString);
                        } else
                        {
                            if(isOverlay) // Genie MC for NuWro (no EXT or dirt events in NuWro)
                            {
                                // std::cout<<"!isNuWro - name"<<name<<" - weight: "<<weight<<" - isSelected: "<< isSelected<<" - recoValue: "<<recoValue<<" - trueValue: "<<trueValue<<std::endl;
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




    // -------------------------------------------------------------------------------------------------------------------------------------
    // Calculate the sideband weights
    // -------------------------------------------------------------------------------------------------------------------------------------
    // Loop over all cross-section objects
    typedef std::pair<std::vector<Double_t>,std::vector<Double_t>> paramAndErrorPair; // Todo: Improve code!
    //Parameters: selectionName name
    std::map<std::string,std::map<std::string, std::vector<float>>> cc0piCovarianceMapData, cc0piCovarianceMapNuWro; 
    std::map<std::string, std::map<std::string, paramAndErrorPair>> cc0piNominalConstraintMapData, cc0piNominalConstraintMapNuWro;
    //Parameters: selectionName name paramName (i.e. golden muonMomentum hadronProduction)
    std::map<std::string, std::map<std::string, std::map<std::string, std::vector<paramAndErrorPair>>>> cc0piUniverseConstraintMapData, cc0piUniverseConstraintMapNuWro;
    
    for (const auto &[dataTypeName, xsecMap] : {std::make_pair(std::string("NuWro"),xsecMapSidebandNuWro), std::make_pair(std::string("Data"),xsecMapSidebandData)})
    {
        for (const auto &[selectionName, xsecs] : xsecMap)
        {
            for (const auto &[name, xsec] : xsecs)
            {
                // -------------------------------------------------------------------------------------------------------------------------------------
                // Fit nominal
                // -------------------------------------------------------------------------------------------------------------------------------------

                std::cout<<"_______________________ExtractNuWroSidebandFit Fitting: "<<dataTypeName<<" - "<<selectionName<<" - "<<name<<"_______________________"<<std::endl;

                // Get the smearing matrix of selected events
                const auto smearingMatrix = xsec.GetSmearingMatrixAllSelected();
                const auto smearingMatrixGeneral = xsec.GetSmearingMatrix();

                const auto selectedEventsBackgroundReco = xsec.GetSelectedBackgroundEvents();
                const auto selectedEventsData = xsec.GetSelectedBNBDataEvents();
                const auto selectedEventsSignalTruth = xsec.GetSelectedSignalEvents();
                const auto eventsSignalTruth = xsec.GetSignalEvents();
                auto signalData = selectedEventsData - selectedEventsBackgroundReco;

                FormattingHelper::SaveMatrix(eventsSignalTruth, "NuWroSidebandFit_" + selectionName + "_" + name + "_eventsSignalTruth.txt");
                FormattingHelper::SaveMatrix(selectedEventsSignalTruth, "NuWroSidebandFit_" + selectionName + "_" + name + "_selectedEventsSignalTruth.txt");
                FormattingHelper::SaveMatrix(selectedEventsData, "NuWroSidebandFit_" + selectionName + "_" + name + "_selectedEventsData.txt");
                FormattingHelper::SaveMatrix(selectedEventsBackgroundReco, "NuWroSidebandFit_" + selectionName + "_" + name + "_selectedEventsBackgroundReco.txt");
                FormattingHelper::SaveMatrix(smearingMatrix, "NuWroSidebandFit_" + selectionName + "_" + name + "_smearingMatrix.txt");
                FormattingHelper::SaveMatrix(smearingMatrixGeneral, "NuWroSidebandFit_" + selectionName + "_" + name + "_smearingMatrixGeneral.txt");


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

                std::vector<float> elements;
                const auto nBins = selectedEventsData.GetRows();
                // std::cout<<"_______________________Fitting Point 1.1"<<std::endl;
                for (unsigned int iBin = 0; iBin < nBins; ++iBin)
                {
                    // std::cout<<"_______________________Fitting Point 1.2"<<std::endl;
                    const auto value = selectedEventsData.At(iBin, 0);
                    if (value<0)
                    {
                        std::cout<<"ERROR: ExtractXSec - Background-removed signal data is negative."<<std::endl;
                        // throw std::logic_error("ERROR: ExtractXSec - Background-removed signal data is negative."); // TODO: Uncomment
                    }
                    elements.push_back(AnalysisHelper::GetCountUncertainty(std::max(value,0.f)));
                }
                
                xNuWro = selectedEventsSignalTruth.GetValues();
                yNuWro = signalData.GetValues();
                errorYNuWro = elements;//signalDataUncertainty.GetValues();
                underflowNuWro = false; //xsec.HasUnderflow(); //Also fit under and overflow bins
                overflowNuWro = false; //xsec.HasOverflow();
                SNuWro = smearingMatrix.GetValues();

                // std::cout<<"\nxNuWro (selectedEventsSignalTruth): \n";
                // for (const auto &xValue : xNuWro)
                //     std::cout<<xValue<<" ";

                // std::cout<<"\nyNuWro (signalData): \n";
                // for (const auto &yValue : yNuWro)
                //     std::cout<<yValue<<" ";

                // std::cout<<"\nerrorYNuWro (signalDataUncertainty): \n";
                // for (const auto &errorYValue : errorYNuWro)
                //     std::cout<<errorYValue<<" ";

                // std::cout<<"\nSmearing matrix: \n";
                // for(unsigned int i = 0; i<SNuWro.size(); i++)
                // {
                //     if(i%nBins==0)
                //         std::cout<<"\n";
                //     std::cout<<SNuWro[i]<<" ";
                // }
                

                auto minimizer = FittingHelper(nBins);
                std::pair<std::vector<Double_t>, std::vector<Double_t>> result;
                
                std::vector<float> fitCovMatrixVector;
                bool successful = false;
                minimizer.Fit(fcnNuWro, result, successful, fitCovMatrixVector, 0);
                if(!successful)
                {
                    std::cout<<"ERROR: ExtractXSec - Fit failed."<<std::endl;
                    throw std::logic_error("ERROR: ExtractXSec - Nominal fit failed.");
                }

                if(dataTypeName=="NuWro")
                {
                    cc0piCovarianceMapNuWro[selectionName].emplace(name, fitCovMatrixVector);
                }
                else if (dataTypeName=="Data")
                {
                    cc0piCovarianceMapData[selectionName].emplace(name, fitCovMatrixVector);
                }
                else
                {
                    throw std::logic_error("ERROR: Unknown dataTypeName");
                }

                const ubsmear::UBMatrix sidebandCovMatrix(fitCovMatrixVector, nBins, nBins);
                FormattingHelper::SaveMatrix(sidebandCovMatrix, dataTypeName + "SidebandFit_" + selectionName + "_" + name + "_sideband_stat_covariance.txt");
                vector<float> paramVector(result.first.begin(), result.first.end()); //Todo avoid this
                vector<float> paramErrorVector(result.second.begin(), result.second.end()); //Todo avoid this
                const ubsmear::UBMatrix sidebandParamVectorTruth(paramVector, nBins, 1);
                const ubsmear::UBMatrix sidebandErrorVectorTruth(paramErrorVector, nBins, 1); 
                FormattingHelper::SaveMatrix(sidebandParamVectorTruth, dataTypeName + "SidebandFit_" + selectionName + "_" + name + "_sideband_parameterVector.txt");
                FormattingHelper::SaveMatrix(sidebandErrorVectorTruth, dataTypeName + "SidebandFit_" + selectionName + "_" + name + "_sideband_parameterErrorVector.txt");

                if(dataTypeName=="NuWro")
                {
                    cc0piNominalConstraintMapNuWro[selectionName].emplace(name, result);
                }
                else if (dataTypeName=="Data")
                {
                    cc0piNominalConstraintMapData[selectionName].emplace(name, result);
                }
                else
                {
                    throw std::logic_error("ERROR: Unknown dataTypeName");
                }

                // const auto sidebandWeights = xsec.GetSidebandWeights(scalingData);
                const auto weightDimensions = {std::make_pair("xsec", systParams.xsecDimensions), std::make_pair("reint", systParams.reintDimensions), std::make_pair("flux", systParams.fluxDimensions)};
                for (const auto &[group, dimensions] : weightDimensions)
                {
                    for (const auto &[paramName, nUniverses] : dimensions)
                    {
                        // -------------------------------------------------------------------------------------------------------------------------------------
                        // Fit each universe
                        // -------------------------------------------------------------------------------------------------------------------------------------

                        // const auto nUniverses = config.extractXSecs.nBootstrapUniverses;

                        std::vector<paramAndErrorPair> resultVector;
                        for (unsigned int iUni = 0; iUni < nUniverses; ++iUni)
                        {
                            std::cout<<"++++++++++++++ExtractNuWroSidebandFit Fitting: "<<selectionName<<" - "<<name<<" - universe: "<<iUni<<" - group: "<<group<<" - paramName: "<<paramName<<"++++++++++++++"<<std::endl;
                            const auto selectedSignalTruthUniverses = xsec.GetSelectedSignalRecoTruthMap().at(group).at(paramName).at(iUni);
                            const auto selectedBackgroundRecoUniverses = xsec.GetSelectedBackgroundRecoMap().at(group).at(paramName).at(iUni);
                            std::pair<std::vector<Double_t>, std::vector<Double_t>> result;
                            const auto p_smearingMatrixInUniverse = xsec.GetSmearingMatrixInUniverseAllSelected(group, paramName, iUni);//xsec.GetSmearingMatrixAllSelected();
                            if(p_smearingMatrixInUniverse)
                            {
                                auto smearingMatrixInUniverse = *p_smearingMatrixInUniverse;

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

                                xNuWro = selectedSignalTruthInUniverse.GetValues();
                                yNuWro = signalDataInUniverse.GetValues();
                                SNuWro = smearingMatrixInUniverse.GetValues();
                                //overflow variable same as nominal
                                // errorY = elements;
                                std::vector<float> covMatrixInUniverse;

                                bool successful = false;
                                minimizer.Fit(fcnNuWro, result, successful, covMatrixInUniverse, 0);
                            }
                            else
                            {
                                std::cout<<"Fit not possible - Smearing matrix not available."<<std::endl;
                                result = std::make_pair(std::vector<Double_t>(nBins, -1), std::vector<Double_t>(nBins, -1));
                            }
                            
                            resultVector.push_back(result);
                        }
                        if(dataTypeName=="NuWro")
                        {
                            cc0piUniverseConstraintMapNuWro[selectionName][name].emplace(paramName, resultVector);
                        }
                        else if (dataTypeName=="Data")
                        {
                            cc0piUniverseConstraintMapData[selectionName][name].emplace(paramName, resultVector);
                        }
                        else
                        {
                            throw std::logic_error("ERROR: Unknown dataTypeName");
                        }
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
