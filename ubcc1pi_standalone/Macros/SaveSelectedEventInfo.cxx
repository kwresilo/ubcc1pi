/**
 *  @file  ubcc1pi_standalone/Macros/SaveSelectedEventInfo.cxx
 *
 *  @brief The implementation file of the SaveSelectedEventInfo macro
 */

#include "ubcc1pi_standalone/Macros/Macros.h"

#include "ubcc1pi_standalone/Objects/FileReader.h"

#include "ubcc1pi_standalone/Helpers/PlottingHelper.h"
#include "ubcc1pi_standalone/Helpers/AnalysisHelper.h"
#include "ubcc1pi_standalone/Helpers/BDTHelper.h"
#include "ubcc1pi_standalone/Helpers/SelectionHelper.h"
#include "ubcc1pi_standalone/Helpers/NormalisationHelper.h"


using namespace ubcc1pi;

namespace ubcc1pi_macros
{

void SaveSelectedEventInfo(const Config &config)
{
    //
    // Setup the input files
    // 
    std::vector< std::tuple<AnalysisHelper::SampleType, std::string, float> > inputData;
    
    inputData.emplace_back(AnalysisHelper::Overlay, config.files.overlaysFileName, NormalisationHelper::GetOverlaysNormalisation(config)); 
    inputData.emplace_back(AnalysisHelper::Dirt,    config.files.dirtFileName,     NormalisationHelper::GetDirtNormalisation(config)); 
    inputData.emplace_back(AnalysisHelper::DataEXT, config.files.dataEXTFileName,  NormalisationHelper::GetDataEXTNormalisation(config)); 
    inputData.emplace_back(AnalysisHelper::DataBNB, config.files.dataBNBFileName,  1.f);

    // Setup the output file
    TFile *pFile = new TFile("selectedEventData.root", "RECREATE");
    TTree *pTree = new TTree("events", "events");

    // Setup the variables to bind to the output branches
    int o_event, o_run, o_subRun, o_sampleType;
    float o_normalisation, o_weight;
    bool o_isCC1Pi;
    std::string o_classification;
    AnalysisHelper::AnalysisData o_recoData, o_truthData;
    int o_recoMuonHitsU, o_recoMuonHitsV, o_recoMuonHitsW, o_recoMuonClusters;
    int o_recoPionHitsU, o_recoPionHitsV, o_recoPionHitsW, o_recoPionClusters;

    // Bind the variables to the output branches
    pTree->Branch("event", &o_event);
    pTree->Branch("run", &o_run);
    pTree->Branch("subRun", &o_subRun);
    pTree->Branch("sampleType", &o_sampleType);
    pTree->Branch("normalisation", &o_normalisation);
    pTree->Branch("weight", &o_weight);
    pTree->Branch("classification", &o_classification);
    
    pTree->Branch("reco_muonMomentum", &o_recoData.muonMomentum);
    pTree->Branch("reco_muonCosTheta", &o_recoData.muonCosTheta);
    pTree->Branch("reco_muonPhi", &o_recoData.muonPhi);
    pTree->Branch("reco_pionMomentum", &o_recoData.pionMomentum);
    pTree->Branch("reco_pionCosTheta", &o_recoData.pionCosTheta);
    pTree->Branch("reco_pionPhi", &o_recoData.pionPhi);
    pTree->Branch("reco_muonPionAngle", &o_recoData.muonPionAngle);
    pTree->Branch("reco_nProtons", &o_recoData.nProtons);
    pTree->Branch("reco_hasGoldenPion", &o_recoData.hasGoldenPion);

    pTree->Branch("reco_muonHitsU", &o_recoMuonHitsU);
    pTree->Branch("reco_muonHitsV", &o_recoMuonHitsV);
    pTree->Branch("reco_muonHitsW", &o_recoMuonHitsW);
    pTree->Branch("reco_muonClusters", &o_recoMuonClusters);
    pTree->Branch("reco_pionHitsU", &o_recoPionHitsU);
    pTree->Branch("reco_pionHitsV", &o_recoPionHitsV);
    pTree->Branch("reco_pionHitsW", &o_recoPionHitsW);
    pTree->Branch("reco_pionClusters", &o_recoPionClusters);
    
    pTree->Branch("truth_muonMomentum", &o_truthData.muonMomentum);
    pTree->Branch("truth_muonCosTheta", &o_truthData.muonCosTheta);
    pTree->Branch("truth_muonPhi", &o_truthData.muonPhi);
    pTree->Branch("truth_pionMomentum", &o_truthData.pionMomentum);
    pTree->Branch("truth_pionCosTheta", &o_truthData.pionCosTheta);
    pTree->Branch("truth_pionPhi", &o_truthData.pionPhi);
    pTree->Branch("truth_muonPionAngle", &o_truthData.muonPionAngle);
    pTree->Branch("truth_nProtons", &o_truthData.nProtons);
    pTree->Branch("truth_hasGoldenPion", &o_truthData.hasGoldenPion);
    
    pTree->Branch("truth_isCC1Pi", &o_isCC1Pi);

    //
    // Get the selection
    //
    auto selection = SelectionHelper::GetDefaultSelection();
    const auto allCuts = selection.GetCuts();
                
    // Loop over the events
    for (const auto [sampleType, fileName, normalisation] : inputData)
    {
        std::cout << "Reading input file: " << fileName << std::endl;

        FileReader reader(fileName);
        auto pEvent = reader.GetBoundEventAddress();

        const auto nEvents = reader.GetNumberOfEvents();
        for (unsigned int i = 0; i < nEvents; ++i)
        {
            AnalysisHelper::PrintLoadingBar(i, nEvents);
            reader.LoadEvent(i);

            // Run the event selection and store which cuts are passed
            std::vector<std::string> cutsPassed;
            std::vector<int> assignedPdgCodes;
            const auto passesGoldenPionSelection = selection.Execute(pEvent, cutsPassed, assignedPdgCodes);
            const auto passesGenericSelection = (std::find(cutsPassed.begin(), cutsPassed.end(), config.global.lastCutGeneric) != cutsPassed.end());
         
            // Only care about selected events
            if (!passesGenericSelection)
                continue;

            // Set the output variables
            o_event = pEvent->metadata.event();
            o_run = pEvent->metadata.run();
            o_subRun = pEvent->metadata.subRun();
            o_sampleType = static_cast<int>(sampleType);
            o_normalisation = normalisation;
            o_weight = AnalysisHelper::GetNominalEventWeight(pEvent);
            o_classification = AnalysisHelper::GetClassificationString(pEvent, config.global.useAbsPdg, config.global.countProtonsInclusively);
            o_recoData = AnalysisHelper::GetRecoAnalysisData(pEvent->reco, assignedPdgCodes, passesGoldenPionSelection);
            o_isCC1Pi = (sampleType == AnalysisHelper::Overlay && AnalysisHelper::IsTrueCC1Pi(pEvent, config.global.useAbsPdg));

            //// BEGIN TEST
            const auto &recoParticles = pEvent->reco.particles;
            const auto &muon = recoParticles.at(AnalysisHelper::GetParticleIndexWithPdg(assignedPdgCodes, 13));
            const auto &pion = recoParticles.at(AnalysisHelper::GetParticleIndexWithPdg(assignedPdgCodes, 211));

            o_recoMuonHitsU = muon.nHitsU();
            o_recoMuonHitsV = muon.nHitsV();
            o_recoMuonHitsW = muon.nHitsW();
            o_recoMuonClusters = (muon.nHitsU() > 0 ? 1 : 0) + (muon.nHitsV() > 0 ? 1 : 0) + (muon.nHitsW() > 0 ? 1 : 0);
            
            o_recoPionHitsU = pion.nHitsU();
            o_recoPionHitsV = pion.nHitsV();
            o_recoPionHitsW = pion.nHitsW();
            o_recoPionClusters = (pion.nHitsU() > 0 ? 1 : 0) + (pion.nHitsV() > 0 ? 1 : 0) + (pion.nHitsW() > 0 ? 1 : 0);
            //// END TEST

            if (o_isCC1Pi)
            {
                o_truthData = AnalysisHelper::GetTruthAnalysisData(pEvent->truth, config.global.useAbsPdg, config.global.protonMomentumThreshold);
            }
            else
            {
                // Set dummy values
                o_truthData.muonMomentum = -std::numeric_limits<float>::max();
                o_truthData.muonCosTheta = -std::numeric_limits<float>::max();
                o_truthData.muonPhi = -std::numeric_limits<float>::max();
                o_truthData.pionMomentum = -std::numeric_limits<float>::max();
                o_truthData.pionCosTheta = -std::numeric_limits<float>::max();
                o_truthData.pionPhi = -std::numeric_limits<float>::max();
                o_truthData.muonPionAngle = -std::numeric_limits<float>::max();
                o_truthData.nProtons = std::numeric_limits<unsigned int>::max();
                o_truthData.hasGoldenPion = false;
            }

            pTree->Fill();
        }
    }

    pFile->Write();
    pFile->Close();
}

} // namespace ubcc1pi_macros
