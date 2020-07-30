/**
 *  @file  ubcc1pi_standalone/Macros/ExtractXSecs.cxx
 *
 *  @brief The implementation file of the ExtractXSecs macro
 */

#include "ubcc1pi_standalone/Macros/Macros.h"
#include "ubcc1pi_standalone/Objects/FileReader.h"

#include "ubcc1pi_standalone/Helpers/AnalysisHelper.h"
#include "ubcc1pi_standalone/Helpers/SelectionHelper.h"
#include "ubcc1pi_standalone/Helpers/CrossSectionHelper.h"
#include "ubcc1pi_standalone/Helpers/NormalisationHelper.h"

using namespace ubcc1pi;

namespace ubcc1pi_macros
{

void ExtractXSecs(const Config &config)
{
    //
    // Setup the input files
    // 
    std::vector< std::tuple<AnalysisHelper::SampleType, std::string, float> > inputData;
    
    inputData.emplace_back(AnalysisHelper::Overlay, config.files.overlaysFileName, NormalisationHelper::GetOverlaysNormalisation(config)); 
    inputData.emplace_back(AnalysisHelper::Dirt,    config.files.dirtFileName,     NormalisationHelper::GetDirtNormalisation(config)); 
    inputData.emplace_back(AnalysisHelper::DataEXT, config.files.dataEXTFileName,  NormalisationHelper::GetDataEXTNormalisation(config)); 
    inputData.emplace_back(AnalysisHelper::DataBNB, config.files.dataBNBFileName,  1.f);


    // Get the selection
    auto selection = SelectionHelper::GetDefaultSelection();
   
    const auto muonMomentumCutoff = config.global.muonMomentum.min;
    const auto pionMomentumCutoff = config.global.pionMomentum.min;
    const auto muonPionAngleCutoff = config.global.muonPionAngle.max;

    // Set up the cross-section objects
    auto xSec_muonCosTheta = CrossSectionHelper::XSec("xsec_muonCosTheta.root", config.global.muonCosTheta.min, config.global.muonCosTheta.max, config.global.useAbsPdg, config.global.countProtonsInclusively);
    auto xSec_muonPhi = CrossSectionHelper::XSec("xsec_muonPhi.root", config.global.muonPhi.min, config.global.muonPhi.max, config.global.useAbsPdg, config.global.countProtonsInclusively);
    auto xSec_muonMomentum = CrossSectionHelper::XSec("xsec_muonMomentum.root", config.global.muonMomentum.min, config.global.muonMomentum.max, config.global.useAbsPdg, config.global.countProtonsInclusively);
    auto xSec_pionCosTheta = CrossSectionHelper::XSec("xsec_pionCosTheta.root", config.global.pionCosTheta.min, config.global.pionCosTheta.max, config.global.useAbsPdg, config.global.countProtonsInclusively);
    auto xSec_pionPhi = CrossSectionHelper::XSec("xsec_pionPhi.root", config.global.pionPhi.min, config.global.pionPhi.max, config.global.useAbsPdg, config.global.countProtonsInclusively);
    auto xSec_pionMomentum = CrossSectionHelper::XSec("xsec_pionMomentum.root", config.global.pionMomentum.min, config.global.pionMomentum.max, config.global.useAbsPdg, config.global.countProtonsInclusively);
    auto xSec_muonPionAngle = CrossSectionHelper::XSec("xsec_muonPionAngle.root", config.global.muonPionAngle.min, config.global.muonPionAngle.max, config.global.useAbsPdg, config.global.countProtonsInclusively);
    auto xSec_nProtons = CrossSectionHelper::XSec("xsec_nProtons.root", config.global.nProtons.min, config.global.nProtons.max, config.global.useAbsPdg, config.global.countProtonsInclusively);
    
    //
    // Fill the cross-section objects
    //
    for (const auto [sampleType, fileName, normalisation] : inputData)
    {
        std::cout << "Reading input file: " << fileName << std::endl;

        FileReader reader(fileName);
        auto pEvent = reader.GetBoundEventAddress();
        const auto isOverlay = (sampleType == AnalysisHelper::Overlay);
        const auto isDataBNB = (sampleType == AnalysisHelper::DataBNB);

        const auto nEvents = reader.GetNumberOfEvents();
        for (unsigned int i = 0; i < nEvents; ++i)
        {
            AnalysisHelper::PrintLoadingBar(i, nEvents);
            reader.LoadEvent(i);

            const auto weight = AnalysisHelper::GetNominalEventWeight(pEvent);

            // Run the event selection and store which cuts are passed
            std::vector<std::string> cutsPassed;
            std::vector<int> assignedPdgCodes;
            auto passedGoldenSelection = selection.Execute(pEvent, cutsPassed, assignedPdgCodes);
            auto passedGenericSelection = (std::find(cutsPassed.begin(), cutsPassed.end(), config.global.lastCutGeneric) != cutsPassed.end());

            // Start with dummy reco data
            AnalysisHelper::AnalysisData recoData;
            recoData.muonMomentum = -std::numeric_limits<float>::max();
            recoData.muonCosTheta = -std::numeric_limits<float>::max();
            recoData.muonPhi = -std::numeric_limits<float>::max();
            recoData.pionMomentum = -std::numeric_limits<float>::max();
            recoData.pionCosTheta = -std::numeric_limits<float>::max();
            recoData.pionPhi = -std::numeric_limits<float>::max();
            recoData.muonPionAngle = -std::numeric_limits<float>::max();
            recoData.nProtons = std::numeric_limits<unsigned int>::max();
            recoData.hasGoldenPion = false;

            // If we passed the generic selection we should be able to access the reco information
            bool passedThresholds = true;
            if (passedGenericSelection)
            {
                recoData = AnalysisHelper::GetRecoAnalysisData(pEvent->reco, assignedPdgCodes, passedGoldenSelection);

                // Check that we are within the defined thresholds and warn if not.
                if (recoData.muonMomentum <= muonMomentumCutoff)
                {
                    std::cerr << "Warning - found event with reco muon momentum: " << recoData.muonMomentum << " <= " << muonMomentumCutoff << std::endl;
                    passedThresholds = false;
                }
                
                if (recoData.pionMomentum <= pionMomentumCutoff)
                {
                    std::cerr << std::endl << std::endl;
                    std::cerr << "Warning - found event with reco pion momentum: " << recoData.pionMomentum << " <= " << pionMomentumCutoff << std::endl;

                    passedThresholds = false;
                }
                
                if (recoData.muonPionAngle >= muonPionAngleCutoff)
                {
                    std::cerr << "Warning - found event with reco muon-pion angle: " << recoData.muonPionAngle << " >= " << muonPionAngleCutoff << std::endl;
                    passedThresholds = false;
                }
            }
    

            // Here we explicity reject an event if it has reconstructed information outside of the thresholds
            passedGenericSelection = passedThresholds && passedGenericSelection;
            passedGoldenSelection = passedThresholds && passedGoldenSelection;

            // Overlay signal events
            if (isOverlay && AnalysisHelper::IsTrueCC1Pi(pEvent, config.global.useAbsPdg))
            {
                const auto truthData = AnalysisHelper::GetTruthAnalysisData(pEvent->truth, config.global.useAbsPdg, config.global.protonMomentumThreshold);

                // Apply the phase-space restriction to the signal definition
                if (truthData.muonMomentum > muonMomentumCutoff &&
                    truthData.pionMomentum > pionMomentumCutoff &&
                    truthData.muonPionAngle < muonPionAngleCutoff)
                {
                    xSec_muonCosTheta.AddSignalEvent(pEvent, passedGenericSelection, truthData.muonCosTheta, recoData.muonCosTheta, weight, normalisation);
                    xSec_muonPhi.AddSignalEvent(pEvent, passedGenericSelection, truthData.muonPhi, recoData.muonPhi, weight, normalisation);
                    xSec_muonMomentum.AddSignalEvent(pEvent, passedGenericSelection, truthData.muonMomentum, recoData.muonMomentum, weight, normalisation);
                    xSec_pionCosTheta.AddSignalEvent(pEvent, passedGenericSelection, truthData.pionCosTheta, recoData.pionCosTheta, weight, normalisation);
                    xSec_pionPhi.AddSignalEvent(pEvent, passedGenericSelection, truthData.pionPhi, recoData.pionPhi, weight, normalisation);
                    xSec_muonPionAngle.AddSignalEvent(pEvent, passedGenericSelection, truthData.muonPionAngle, recoData.muonPionAngle, weight, normalisation);
                    xSec_nProtons.AddSignalEvent(pEvent, passedGenericSelection, truthData.nProtons, recoData.nProtons, weight, normalisation);
                    
                    xSec_pionMomentum.AddSignalEvent(pEvent, passedGoldenSelection, truthData.pionMomentum, recoData.pionMomentum, weight, normalisation);
                
                    continue;
                }
            }

            if (!passedGenericSelection)
                continue;

            // BNB data
            if (isDataBNB)
            {
                xSec_muonCosTheta.AddSelectedBNBDataEvent(pEvent, recoData.muonCosTheta);
                xSec_muonPhi.AddSelectedBNBDataEvent(pEvent, recoData.muonPhi);
                xSec_pionCosTheta.AddSelectedBNBDataEvent(pEvent, recoData.pionCosTheta);
                xSec_pionPhi.AddSelectedBNBDataEvent(pEvent, recoData.pionPhi);
                xSec_muonMomentum.AddSelectedBNBDataEvent(pEvent, recoData.muonMomentum);
                xSec_muonPionAngle.AddSelectedBNBDataEvent(pEvent, recoData.muonPionAngle);
                xSec_nProtons.AddSelectedBNBDataEvent(pEvent, recoData.nProtons);

                if (passedGoldenSelection)
                    xSec_pionMomentum.AddSelectedBNBDataEvent(pEvent, recoData.pionMomentum);
            }
            // EXT, Dirt, Overlay backgrounds
            else 
            {
                xSec_muonCosTheta.AddSelectedBackgroundEvent(pEvent, sampleType, recoData.muonCosTheta, weight, normalisation);
                xSec_muonPhi.AddSelectedBackgroundEvent(pEvent, sampleType, recoData.muonPhi, weight, normalisation);
                xSec_pionCosTheta.AddSelectedBackgroundEvent(pEvent, sampleType, recoData.pionCosTheta, weight, normalisation);
                xSec_pionPhi.AddSelectedBackgroundEvent(pEvent, sampleType, recoData.pionPhi, weight, normalisation);
                xSec_muonMomentum.AddSelectedBackgroundEvent(pEvent, sampleType, recoData.muonMomentum, weight, normalisation);
                xSec_muonPionAngle.AddSelectedBackgroundEvent(pEvent, sampleType, recoData.muonPionAngle, weight, normalisation);
                xSec_nProtons.AddSelectedBackgroundEvent(pEvent, sampleType, recoData.nProtons, weight, normalisation);

                if (passedGoldenSelection)
                    xSec_pionMomentum.AddSelectedBackgroundEvent(pEvent, sampleType, recoData.pionMomentum, weight, normalisation);
            }
        }
    }
    
    // Set the binning
    std::cout << "Setting bins" << std::endl;
    xSec_muonCosTheta.SetBins(config.global.muonCosTheta.binEdges);
    xSec_muonPhi.SetBins(config.global.muonPhi.binEdges);
    xSec_muonMomentum.SetBins(config.global.muonMomentum.binEdges);
    xSec_pionCosTheta.SetBins(config.global.pionCosTheta.binEdges);
    xSec_pionPhi.SetBins(config.global.pionPhi.binEdges);
    xSec_pionMomentum.SetBins(config.global.pionMomentum.binEdges);
    xSec_muonPionAngle.SetBins(config.global.muonPionAngle.binEdges);
    xSec_nProtons.SetBins(config.global.nProtons.binEdges);

    // Print the results
    std::cout << "Muon cos(theta)" << std::endl;
    xSec_muonCosTheta.PrintBinContents("xsec_muonCosTheta");
    xSec_muonCosTheta.MakePlots("xsec_muonCosTheta");
    
    std::cout << "Muon phi" << std::endl;
    xSec_muonPhi.PrintBinContents("xsec_muonPhi");
    xSec_muonPhi.MakePlots("xsec_muonPhi");
    
    std::cout << "Muon momentum" << std::endl;
    xSec_muonMomentum.PrintBinContents("xsec_muonMomentum");
    xSec_muonMomentum.MakePlots("xsec_muonMomentum");
    
    std::cout << "Pion cos(theta)" << std::endl;
    xSec_pionCosTheta.PrintBinContents("xsec_pionCosTheta");
    xSec_pionCosTheta.MakePlots("xsec_pionCosTheta");
    
    std::cout << "Pion phi" << std::endl;
    xSec_pionPhi.PrintBinContents("xsec_pionPhi");
    xSec_pionPhi.MakePlots("xsec_pionPhi");
    
    std::cout << "Pion momentum" << std::endl;
    xSec_pionMomentum.PrintBinContents("xsec_pionMomentum");
    xSec_pionMomentum.MakePlots("xsec_pionMomentum");
    
    std::cout << "Muon-pion angle" << std::endl;
    xSec_muonPionAngle.PrintBinContents("xsec_muonPionAngle");
    xSec_muonPionAngle.MakePlots("xsec_muonPionAngle");
    
    std::cout << "nProtons" << std::endl;
    xSec_nProtons.PrintBinContents("xsec_nProtons");
    xSec_nProtons.MakePlots("xsec_nProtons");
}

} // namespace ubcc1pi_macros
