/**
 *  @file  ubcc1pi_standalone/Macros/ExtractSideband.cxx
 *
 *  @brief The implementation file of the ExtractSideband macro
 */

#include "ubcc1pi_standalone/Macros/Macros.h"
#include "ubcc1pi_standalone/Objects/FileReader.h"

#include "ubcc1pi_standalone/Helpers/AnalysisHelper.h"
#include "ubcc1pi_standalone/Helpers/NormalisationHelper.h"
#include "ubcc1pi_standalone/Helpers/SelectionHelper.h"
#include "ubcc1pi_standalone/Helpers/CrossSectionHelper.h"
#include "ubcc1pi_standalone/Helpers/FormattingHelper.h"
#include "ubsmear.h"

using namespace ubcc1pi;

namespace ubcc1pi_macros
{

void ExtractSideband(const Config &config)
{
    ROOT::EnableImplicitMT(2);

    std::vector<PlottingHelper::MultiPlot> momentumPlots; //Debug
    const std::string yLabel = "Number of particles (norm. to bin width)";
    momentumPlots.emplace_back("Proton momentum / GeV", yLabel, 10, 0, 1, true, config.global.axisTitles);
    momentumPlots.emplace_back("Muon momentum / GeV", yLabel, 10, 0, 3, true, config.global.axisTitles);
    momentumPlots.emplace_back("Muon momentum xSec-binning/ GeV", yLabel, config.global.muonMomentum.binEdges, true, config.global.axisTitles);

    std::vector<PlottingHelper::MultiPlot> momentumPlotsControll; //Debug
    momentumPlotsControll.emplace_back("Scaled Proton momentum / GeV", yLabel, 10, 0, 1, true, config.global.axisTitles);
    momentumPlotsControll.emplace_back("Scaled Muon momentum / GeV", yLabel, 10, 0, 3, true, config.global.axisTitles);
    momentumPlotsControll.emplace_back("Proton momentum / GeV", yLabel, 10, 0, 1, true, config.global.axisTitles);
    momentumPlotsControll.emplace_back("Muon momentum / GeV", yLabel, 10, 0, 3, true, config.global.axisTitles);
    momentumPlotsControll.emplace_back("Muon momentum xSec-binning/ GeV", yLabel, config.global.muonMomentum.binEdges, true, config.global.axisTitles);

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
    // std::cout << "- Flux:            " << scalingData.pFluxReweightor->GetIntegratedNominalFlux() << " * 10^-10 cm^-2 POT^-1" << std::endl;
    // std::cout << "- Exposure:        " << scalingData.exposurePOT << " * 10^20 POT" << std::endl;
    // std::cout << "- Target density:  " << config.global.targetDensity << " * 10^31 nucleons/cm^3" << std::endl;
    // std::cout << "- Fiducial volume: " << AnalysisHelper::GetFiducialVolume() << " * cm^3" << std::endl;
    // std::cout << "- nTargets:        " << scalingData.nTargets << " * 10^31 nucleons" << std::endl;

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
    // std::map<std::string, std::map<std::string, CrossSectionHelper::CrossSection> > xsecMap;
    std::map<std::string, std::map<std::string, CrossSectionHelper::CrossSection> > xsecMapSideband;
    // std::map<std::string, std::map<std::string, CrossSectionHelper::CrossSection> > xsecMapSideband2;

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
                // xsecMap[selectionName].emplace(name, CrossSectionHelper::CrossSection(systParams, extendedBinEdges, hasUnderflow, hasOverflow, scaleByBinWidth));
                xsecMapSideband[selectionName].emplace(name, CrossSectionHelper::CrossSection(systParams, extendedBinEdges, hasUnderflow, hasOverflow, scaleByBinWidth));
                // xsecMapSideband2[selectionName].emplace(name, CrossSectionHelper::CrossSection(systParams, extendedBinEdges, hasUnderflow, hasOverflow, scaleByBinWidth));
            }
        }
        // if (config.extractXSecs.crossSectionIsEnabled.at("sideband").at(name))
        //     xsecMapSideband["sideband"].emplace(name, CrossSectionHelper::CrossSection(systParams, extendedBinEdges, hasUnderflow, hasOverflow, scaleByBinWidth));

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
            xsecMapSideband[selectionName].emplace("total", CrossSectionHelper::CrossSection(systParams, {-1.f, 1.f}, false, false, false));
            // xsecMapSideband2[selectionName].emplace("total", CrossSectionHelper::CrossSection(systParams, {-1.f, 1.f}, false, false, false));
        }

    }
    // Don't setup a cross-section object if it's been disabled in the configuration
    // if (config.extractXSecs.crossSectionIsEnabled.at("sideband").at("total"))
    //     xsecMapSideband["sideband"].emplace("total", CrossSectionHelper::CrossSection(systParams, {-1.f, 1.f}, false, false, false));

    // The dummy value that will be used as the "kinematic quantity" for the total cross-section
    const auto dummyValue = 0.f;

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

    scalingData.exposurePOT = config.normsRun1.dataBNBTor875WCut / (1e20);
    scalingData.nTargets = config.global.targetDensity * (1e-8) * AnalysisHelper::GetFiducialVolume();

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
        std::vector< std::tuple<AnalysisHelper::SampleType, std::string, std::string, float> > inputData;

        if(run == 1)
        {
            inputData.emplace_back(AnalysisHelper::Overlay, "", config.filesRun1.overlaysFileName, NormalisationHelper::GetOverlaysNormalisation(config, 1));
            inputData.emplace_back(AnalysisHelper::Dirt,    "", config.filesRun1.dirtFileName, NormalisationHelper::GetDirtNormalisation(config, 1));
            inputData.emplace_back(AnalysisHelper::DataEXT, "", config.filesRun1.dataEXTFileName, NormalisationHelper::GetDataEXTNormalisation(config, 1));
            inputData.emplace_back(AnalysisHelper::DataBNB, "", config.filesRun1.dataBNBFileName, 1.f);

            // // Add the detector variation files
            // for (const auto &[name, fileName] : config.filesRun1.detVarFiles)
            //     inputData.emplace_back(AnalysisHelper::DetectorVariation, name, fileName, NormalisationHelper::GetDetectorVariationNormalisation(config, name, 1));

            // // Read in the flux
            // const auto fluxHistNames = CrossSectionHelper::GetNominalFluxHistNames(config.fluxRun1.nuPdgsSignal, config.fluxRun1.nuPdgToHistName, config.fluxRun1.nomHistPattern);
            // const auto &[fluxBinEdges, fluxValues] = CrossSectionHelper::ReadNominalFlux(config.fluxRun1.fileName, fluxHistNames, config.fluxRun1.pot);
            // scalingData.pFluxReweightor = std::make_shared<CrossSectionHelper::FluxReweightor>(fluxBinEdges, fluxValues, systParams.fluxDimensions);

            // scalingData.exposurePOT = config.normsRun1.dataBNBTor875WCut / (1e20);
            // scalingData.nTargets = config.global.targetDensity * (1e-8) * AnalysisHelper::GetFiducialVolume();

        }
        else if(run == 2)
        {
            inputData.emplace_back(AnalysisHelper::Overlay, "", config.filesRun2.overlaysFileName, NormalisationHelper::GetOverlaysNormalisation(config, 2));
            inputData.emplace_back(AnalysisHelper::Dirt,    "", config.filesRun2.dirtFileName, NormalisationHelper::GetDirtNormalisation(config, 2));
            inputData.emplace_back(AnalysisHelper::DataEXT, "", config.filesRun2.dataEXTFileName, NormalisationHelper::GetDataEXTNormalisation(config, 2));
            inputData.emplace_back(AnalysisHelper::DataBNB, "", config.filesRun2.dataBNBFileName, 1.f);

            // // Add the detector variation files
            // for (const auto &[name, fileName] : config.filesRun2.detVarFiles)
            //     inputData.emplace_back(AnalysisHelper::DetectorVariation, name, fileName, NormalisationHelper::GetDetectorVariationNormalisation(config, name, 2));

            // // Read in the flux
            // const auto fluxHistNames = CrossSectionHelper::GetNominalFluxHistNames(config.fluxRun2.nuPdgsSignal, config.fluxRun2.nuPdgToHistName, config.fluxRun2.nomHistPattern);
            // const auto &[fluxBinEdges, fluxValues] = CrossSectionHelper::ReadNominalFlux(config.fluxRun2.fileName, fluxHistNames, config.fluxRun2.pot);
            // scalingData.pFluxReweightor = std::make_shared<CrossSectionHelper::FluxReweightor>(fluxBinEdges, fluxValues, systParams.fluxDimensions);

            // scalingData.exposurePOT = config.normsRun2.dataBNBTor875WCut / (1e20);
            // scalingData.nTargets = config.global.targetDensity * (1e-8) * AnalysisHelper::GetFiducialVolume();
        }
        else if(run == 3)
        {
            inputData.emplace_back(AnalysisHelper::Overlay, "", config.filesRun3.overlaysFileName, NormalisationHelper::GetOverlaysNormalisation(config, 3));
            inputData.emplace_back(AnalysisHelper::Dirt,    "", config.filesRun3.dirtFileName, NormalisationHelper::GetDirtNormalisation(config, 3));
            inputData.emplace_back(AnalysisHelper::DataEXT, "", config.filesRun3.dataEXTFileName, NormalisationHelper::GetDataEXTNormalisation(config, 3));
            inputData.emplace_back(AnalysisHelper::DataBNB, "", config.filesRun3.dataBNBFileName, 1.f);

            // // Add the detector variation files
            // for (const auto &[name, fileName] : config.filesRun3.detVarFiles)
            //     inputData.emplace_back(AnalysisHelper::DetectorVariation, name, fileName, NormalisationHelper::GetDetectorVariationNormalisation(config, name, 3));

            // // Read in the flux
            // const auto fluxHistNames = CrossSectionHelper::GetNominalFluxHistNames(config.fluxRun3.nuPdgsSignal, config.fluxRun3.nuPdgToHistName, config.fluxRun3.nomHistPattern);
            // const auto &[fluxBinEdges, fluxValues] = CrossSectionHelper::ReadNominalFlux(config.fluxRun3.fileName, fluxHistNames, config.fluxRun3.pot);
            // scalingData.pFluxReweightor = std::make_shared<CrossSectionHelper::FluxReweightor>(fluxBinEdges, fluxValues, systParams.fluxDimensions);

            // scalingData.exposurePOT = config.normsRun3.dataBNBTor875WCut / (1e20);
            // scalingData.nTargets = config.global.targetDensity * (1e-8) * AnalysisHelper::GetFiducialVolume();
        }
        else throw std::logic_error("ExtractSidebandFit - Invalid run number");

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
            const auto isDataBNB = (sampleType == AnalysisHelper::DataBNB);
            const auto isDetVar  = (sampleType == AnalysisHelper::DetectorVariation);

            // Open the input file for reading and enable the branches with systematic event weights (if required)
            FileReader reader(fileName);

            if (isOverlay)
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
                            std::cout << "Event failed reco phase-space cuts" << std::endl;
                            break;
                        }
                    }
                }

                float protonMomentum=0.f;
                float muonMomentum=0.f;
                float weight2=0.f;
                if(passedSidebandSelection)
                {
                    const auto muonFeatureNames = BDTHelper::MuonBDTFeatureNames;
                    BDTHelper::BDT muonBDT("muon", muonFeatureNames);
                    const auto recoParticles = pEvent->reco.particles;
                    const auto protonIndex = SelectionHelper::GetLeadingProtonCandidateIndex(recoParticles, sidebandAssignedPdgCodes);
                    const auto muonIndex = SelectionHelper::GetMuonCandidateIndex(recoParticles, muonFeatureNames, muonBDT);
                    if (protonIndex!=std::numeric_limits<unsigned int>::max())
                    {
                        const auto proton = recoParticles.at(protonIndex);
                        const auto muon = recoParticles.at(muonIndex);
                        const auto plotStyle = PlottingHelper::GetPlotStyle(sampleType, pEvent, config.global.useAbsPdg);
                        protonMomentum = AnalysisHelper::GetProtonMomentumFromRange(proton.range());
                        muonMomentum = AnalysisHelper::GetMuonMomentum(muon);
                        weight2 = AnalysisHelper::GetNominalEventWeight(pEvent) * normalisation;
                        momentumPlots.at(0).Fill(protonMomentum, plotStyle, weight2);
                        momentumPlots.at(1).Fill(muonMomentum, plotStyle, weight2);
                        momentumPlots.at(2).Fill(muonMomentum, plotStyle, weight2);
                    }
                    else
                    {
                        std::cout << "Stacked Histogram - No proton candidate found" << std::endl;
                    }
                }

                const auto isSelectedSideband = passedSidebandSelection  && passesPhaseSpaceReco;
                // std::map<std::string, bool> isSelectedMap = {{"sideband",isSelectedSideband}};
                std::map<std::string, bool> isSelectedMap = {{"generic",isSelectedSideband},{"golden",isSelectedSideband}};

                // -----------------------------------------------------------------------------------------------------------------------------
                // Handle BNB data
                // -----------------------------------------------------------------------------------------------------------------------------
                if (isDataBNB)
                {
                    for (auto &[selectionName, xsecs] : xsecMapSideband)
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
                        const auto plotStyle = PlottingHelper::GetPlotStyle2(0);
                        std::cout<<"\n\n ---- weight2: "<<weight2<<std::endl;
                        momentumPlotsControll.at(0).Fill(protonMomentum, plotStyle, weight2);
                        momentumPlotsControll.at(1).Fill(muonMomentum, plotStyle, weight2);
                        momentumPlotsControll.at(2).Fill(protonMomentum, plotStyle, 1.f);
                        momentumPlotsControll.at(3).Fill(muonMomentum, plotStyle, 1.f);
                        momentumPlotsControll.at(4).Fill(muonMomentum, plotStyle, 1.f);
                    }

                    // For BNB data that's all we need to do!
                    continue;
                }

                // -----------------------------------------------------------------------------------------------------------------------------
                // Work out if this event is signal, and apply any phase-space restrictions based on the input binning
                // -----------------------------------------------------------------------------------------------------------------------------

                // Get the nominal event weight, scaled by the sample normalisation
                const auto weight = AnalysisHelper::GetNominalEventWeight(pEvent) * normalisation;

                // Determine if this is truly a CC0Pi event
                const auto isTrueCC0Pi = (isOverlay || isDetVar) && AnalysisHelper::IsTrueCC0Pi(pEvent, config.global.useAbsPdg, config.global.protonMomentumThreshold);            

                // Get the truth analysis data (if available, otherwise set to dummy values)
                const auto truthData = (
                    (isTrueCC0Pi)
                        ? AnalysisHelper::GetTruthAnalysisDataCC0Pi(pEvent->truth, config.global.useAbsPdg, config.global.protonMomentumThreshold)
                        : AnalysisHelper::GetDummyAnalysisData()
                );

                // Here we apply truth-level phase-space restrictions
                // For all true CC1Pi events, we check if the values of each kinematic variable are within the supplied limits. If not then the
                // event is not classed as "signal"

                // bool passesPhaseSpaceTruth = true;

                bool passesPhaseSpaceTruth = false;
                if (isTrueCC0Pi)
                {
                    // Start by assuming the event passes the phase-space cuts
                    passesPhaseSpaceTruth = true;

                    // Check the value of the kinematic quantities are within the phase-space limits
                    for (const auto &[name, minMax] : phaseSpaceMap)
                    {
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

                // // -----------------------------------------------------------------------------------------------------------------------------
                // // Handle the detector variation samples (unisims)
                // // -----------------------------------------------------------------------------------------------------------------------------
                // if (isDetVar)
                // {
                //     for (auto &[selectionName, xsecs] : xsecMapSideband)
                //     {
                //         const auto isSignal = isSignalMap.at(selectionName);
                //         // Handle signal events
                //         if (isSignal)
                //         {
                //             // Determine if we passed the relevant selection
                //             const auto isSelected = isSelectedMap.at(selectionName);

                //             for (auto &[name, xsec] : xsecs)
                //             {
                //                 const auto recoValue = getSidebandValue.at(name)(recoData);
                //                 const auto trueValue = getSidebandValue.at(name)(truthData);
                //                 xsec.AddSignalEventDetVar(recoValue, trueValue, isSelected, weight, sampleName);
                //             }
                //         }
                //         // Handle selected background events
                //         else
                //         {
                //             for (auto &[selectionName, xsecs] : xsecMapSideband)
                //             {
                //                 // Only use selected background events
                //                 const auto isSelected = isSelectedMap.at(selectionName);
                //                 if (!isSelected)
                //                     continue;

                //                 for (auto &[name, xsec] : xsecs)
                //                 {
                //                     const auto recoValue = getSidebandValue.at(name)(recoData);
                //                     xsec.AddSelectedBackgroundEventDetVar(recoValue, weight, sampleName);
                //                 }
                //             }
                //         }
                //     }

                //     // For detector variation samples, that's all we need to do!
                //     continue;
                // }


                // -----------------------------------------------------------------------------------------------------------------------------
                // Handle all other events (i.e those from the nominal simulation): Overlays, dirt, EXT data
                // -----------------------------------------------------------------------------------------------------------------------------

                // Get the flux weights
                const auto fluxWeights = (
                    isOverlay
                        ? CrossSectionHelper::GetWeightsMap(pEvent->truth, systParams.fluxDimensions, config.extractXSecs.mutuallyExclusiveDimensions)
                        : CrossSectionHelper::GetUnitWeightsMap(systParams.fluxDimensions)
                );

                // Fill the flux-reweightor with all overlay events from neutrinos of the desired flavour in the active volume
                // if (isOverlay && AnalysisHelper::IsInActiveVolume(pEvent->truth.nuVertex()))
                // {
                //     bool foundFluxPdg = false;
                //     if(run == 1) foundFluxPdg = std::find(config.fluxRun1.nuPdgsSignal.begin(), config.fluxRun1.nuPdgsSignal.end(), pEvent->truth.nuPdgCode()) != config.fluxRun1.nuPdgsSignal.end();
                //     else if(run == 2) foundFluxPdg = std::find(config.fluxRun2.nuPdgsSignal.begin(), config.fluxRun2.nuPdgsSignal.end(), pEvent->truth.nuPdgCode()) != config.fluxRun2.nuPdgsSignal.end();
                //     else if(run == 3) foundFluxPdg = std::find(config.fluxRun3.nuPdgsSignal.begin(), config.fluxRun3.nuPdgsSignal.end(), pEvent->truth.nuPdgCode()) != config.fluxRun3.nuPdgsSignal.end();
                //     else throw std::logic_error("ExtractSidebandFit - Invalid run number");

                //     if(foundFluxPdg) scalingData.pFluxReweightor->AddEvent(pEvent->truth.nuEnergy(), weight, fluxWeights);
                // }

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
                xsecWeightsScaleFactor = std::max(xsecWeightsScaleFactor, 0.0001f);

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
                for (auto &[selectionName, xsecs] : xsecMapSideband)
                {
                    const auto isSignal = isSignalMap.at(selectionName);
                    // Handle signal events
                    if (isSignal)
                    {
                        for (auto &[selectionName, xsecs] : xsecMapSideband)
                        {
                            // Determine if we passed the relevant selection
                            const auto isSelected = isSelectedMap.at(selectionName);

                            for (auto &[name, xsec] : xsecs)
                            {
                                const auto recoValue = getSidebandValue.at(name)(recoData);
                                const auto trueValue = getSidebandValue.at(name)(truthData);
                                const auto seedString =  selectionName + name + std::to_string(i);
                                // std::cout<<"Added (selected? "<< isSelected <<") Signal Event: "<<std::endl;
                                xsec.AddSignalEvent(recoValue, trueValue, isSelected, weight, fluxWeights, xsecWeights, reintWeights, seedString);
                            }
                            if(isSelected)
                            {
                                const auto plotStyle = PlottingHelper::GetPlotStyle2(1);
                                momentumPlotsControll.at(0).Fill(protonMomentum, plotStyle, weight);
                                momentumPlotsControll.at(1).Fill(muonMomentum, plotStyle, weight);
                                momentumPlotsControll.at(2).Fill(protonMomentum, plotStyle, weight);
                                momentumPlotsControll.at(3).Fill(muonMomentum, plotStyle, weight);
                                momentumPlotsControll.at(4).Fill(muonMomentum, plotStyle, weight);
                            }
                        }
                    }
                    // Handle selected background events
                    else
                    {
                        for (auto &[selectionName, xsecs] : xsecMapSideband)
                        {
                            // Only use selected background events
                            const auto isSelected = isSelectedMap.at(selectionName);
                            if (!isSelected)
                                continue;

                            for (auto &[name, xsec] : xsecs)
                            {
                                const auto recoValue = getSidebandValue.at(name)(recoData);
                                const auto seedString =  selectionName + name + std::to_string(i);
                                std::vector<float> bootstrapWeights; // Parameter only needed for CC1pi
                                // std::cout<<"Added selected background event."<<std::endl;
                                xsec.AddSelectedBackgroundEvent(recoValue, isDirt, weight, fluxWeights, xsecWeights, reintWeights, bootstrapWeights, seedString);
                            }
                            const auto plotStyle = PlottingHelper::GetPlotStyle2(2);
                            momentumPlotsControll.at(0).Fill(protonMomentum, plotStyle, weight);
                            momentumPlotsControll.at(1).Fill(muonMomentum, plotStyle, weight);
                            momentumPlotsControll.at(2).Fill(protonMomentum, plotStyle, weight);
                            momentumPlotsControll.at(3).Fill(muonMomentum, plotStyle, weight);
                            momentumPlotsControll.at(4).Fill(muonMomentum, plotStyle, weight);
                        }
                    }
                }
            }
        }
        std::cout<<"------------- Run completed -------------"<<std::endl;
    }

    // -------------------------------------------------------------------------------------------------------------------------------------
    // Save the sideband values
    // -------------------------------------------------------------------------------------------------------------------------------------
    for (const auto &[selectionName, xsecs] : xsecMapSideband)
    {
        for (const auto &[name, xsec] : xsecs)
        {
            // -------------------------------------------------------------------------------------------------------------------------------------
            // Fit nominal
            // -------------------------------------------------------------------------------------------------------------------------------------

            // std::cout<<"_______________________Fitting: "<<selectionName<<" - "<<name<<"_______________________"<<std::endl;
            // const auto sidebandCVFit = xsec.GetPredictedCrossSection(scalingData);

            // std::cout<<"_______________________Fitting Point 0.1"<<std::endl;
            const auto selectedEventsData = xsec.GetSelectedBNBDataEvents();
            // std::cout<<"_______________________Fitting Point 0.2"<<std::endl;
            // Get the smearing matrix of selected events
            const auto smearingMatrixAllSelected = xsec.GetSmearingMatrixAllSelected();
            // std::cout<<"_______________________Fitting Point 0.3"<<std::endl;
            const auto selectedEventsBackgroundReco = xsec.GetSelectedBackgroundEvents();
            // std::cout<<"_______________________Fitting Point 0.4"<<std::endl;
            const auto selectedEventsSignalTruth = xsec.GetSelectedSignalEvents();
            // std::cout<<"_______________________Fitting Point 0.5"<<std::endl;
            auto signalData = selectedEventsData - selectedEventsBackgroundReco;

            FormattingHelper::SaveMatrix(selectedEventsSignalTruth, "SidebandFit_" + selectionName + "_" + name + "_selectedEventsSignalTruth.txt");
            FormattingHelper::SaveMatrix(selectedEventsData, "SidebandFit_" + selectionName + "_" + name + "_selectedEventsData.txt");
            FormattingHelper::SaveMatrix(selectedEventsBackgroundReco, "SidebandFit_" + selectionName + "_" + name + "_selectedEventsBackgroundReco.txt");
            FormattingHelper::SaveMatrix(smearingMatrixAllSelected, "SidebandFit_" + selectionName + "_" + name + "_smearingMatrixAllSelected.txt");
        }
    }

    momentumPlots.at(0).SaveAsStacked("reco_cc1pi_protonMomentum_extraction_", false, config.global.scaleByBinWidth, false, config.global.axisTitles);
    momentumPlots.at(1).SaveAsStacked("reco_cc1pi_muonMomentum_extraction_", false, config.global.scaleByBinWidth, false, config.global.axisTitles);
    momentumPlots.at(2).SaveAsStacked("reco_cc1pi_muonMomentum_extraction_xSecBinning_", false, config.global.scaleByBinWidth, false, config.global.axisTitles);
    momentumPlotsControll.at(0).SaveAsStacked("reco_cc1pi_protonMomentum_extraction_controll_scaled_", false, config.global.scaleByBinWidth, false, config.global.axisTitles);
    momentumPlotsControll.at(1).SaveAsStacked("reco_cc1pi_muonMomentum_extraction_controll_scaled_", false, config.global.scaleByBinWidth, false, config.global.axisTitles);
    momentumPlotsControll.at(2).SaveAsStacked("reco_cc1pi_protonMomentum_extraction_controll_", false, config.global.scaleByBinWidth, false, config.global.axisTitles);
    momentumPlotsControll.at(3).SaveAsStacked("reco_cc1pi_muonMomentum_extraction_controll_", false, config.global.scaleByBinWidth, false, config.global.axisTitles);
    momentumPlotsControll.at(4).SaveAsStacked("reco_cc1pi_muonMomentum_extraction_controll_xSecBinning_", false, true, false, config.global.axisTitles);
    momentumPlotsControll.at(4).SaveAsStacked("reco_cc1pi_muonMomentum_extraction_controll_xSecBinning_notScaledByBinWidth_", false, false, false, config.global.axisTitles);

    std::cout<<"------------- All done -------------"<<std::endl;
    return;
}

} // namespace ubcc1pi_macros
