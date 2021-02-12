/**
 *  @file  ubcc1pi_standalone/Macros/DumpSelectedEventInfo.cxx
 *
 *  @brief The implementation file of the DumpSelectedEventInfo macro
 */

#include "ubcc1pi_standalone/Macros/Macros.h"

#include "ubcc1pi_standalone/Objects/FileReader.h"

#include "ubcc1pi_standalone/Helpers/PlottingHelper.h"
#include "ubcc1pi_standalone/Helpers/AnalysisHelper.h"
#include "ubcc1pi_standalone/Helpers/SelectionHelper.h"
#include "ubcc1pi_standalone/Helpers/NormalisationHelper.h"


using namespace ubcc1pi;

namespace ubcc1pi_macros
{

void DumpSelectedEventInfo(const Config &config)
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
    auto pFile = std::make_shared<TFile>("selectedEventInfo.root", "RECREATE");
    auto pTree = std::make_shared<TTree>("events", "events");

    int o_run, o_subRun, o_event;
    pTree->Branch("run", &o_run);
    pTree->Branch("subRun", &o_subRun);
    pTree->Branch("event", &o_event);

    float o_weight;
    pTree->Branch("weight", &o_weight);

    int o_sampleType, o_plotStyle;
    pTree->Branch("sampleType", &o_sampleType);
    pTree->Branch("plotStyle", &o_plotStyle);

    std::string o_classification;
    pTree->Branch("classification", &o_classification);

    AnalysisHelper::AnalysisData o_truthData;
    pTree->Branch("truth_muonMomentum", &o_truthData.muonMomentum);
    pTree->Branch("truth_muonCosTheta", &o_truthData.muonCosTheta);
    pTree->Branch("truth_muonPhi", &o_truthData.muonPhi);
    pTree->Branch("truth_pionMomentum", &o_truthData.pionMomentum);
    pTree->Branch("truth_pionCosTheta", &o_truthData.pionCosTheta);
    pTree->Branch("truth_pionPhi", &o_truthData.pionPhi);
    pTree->Branch("truth_muonPionAngle", &o_truthData.muonPionAngle);
    pTree->Branch("truth_nProtons", &o_truthData.nProtons);
    pTree->Branch("truth_hasGoldenPion", &o_truthData.hasGoldenPion);

    AnalysisHelper::AnalysisData o_recoData;
    pTree->Branch("reco_muonMomentum", &o_recoData.muonMomentum);
    pTree->Branch("reco_muonCosTheta", &o_recoData.muonCosTheta);
    pTree->Branch("reco_muonPhi", &o_recoData.muonPhi);
    pTree->Branch("reco_pionMomentum", &o_recoData.pionMomentum);
    pTree->Branch("reco_pionCosTheta", &o_recoData.pionCosTheta);
    pTree->Branch("reco_pionPhi", &o_recoData.pionPhi);
    pTree->Branch("reco_muonPionAngle", &o_recoData.muonPionAngle);
    pTree->Branch("reco_nProtons", &o_recoData.nProtons);
    pTree->Branch("reco_hasGoldenPion", &o_recoData.hasGoldenPion);

    // BEGIN TEST
    float o_recoPionTrackScore, o_recoPionWiggliness;
    pTree->Branch("reco_pionTrackScore", &o_recoPionTrackScore);
    pTree->Branch("reco_pionWiggliness", &o_recoPionWiggliness);
    // END TEST

    //
    // Get the selection
    //
    auto selection = SelectionHelper::GetDefaultSelection();

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
            const auto &[isSelectedGolden, cutsPassed, assignedPdgCodes] = selection.Execute(pEvent);
            const auto isSelectedGeneric = SelectionHelper::IsCutPassed(cutsPassed, config.global.lastCutGeneric);

            // Only use events that at least pass the generic selection
            if (!isSelectedGeneric)
                continue;

            o_run = pEvent->metadata.run();
            o_subRun = pEvent->metadata.subRun();
            o_event = pEvent->metadata.event();

            o_weight = AnalysisHelper::GetNominalEventWeight(pEvent) * normalisation;
            o_sampleType = sampleType;
            o_plotStyle = PlottingHelper::GetPlotStyle(sampleType, pEvent, config.global.useAbsPdg);
            o_classification = AnalysisHelper::GetClassificationString(pEvent, config.global.useAbsPdg, config.global.countProtonsInclusively);

            // Get the truth and reco analysis data
            o_recoData = AnalysisHelper::GetRecoAnalysisData(pEvent->reco, assignedPdgCodes, isSelectedGolden);
            o_truthData = AnalysisHelper::GetDummyAnalysisData();
            if (sampleType == AnalysisHelper::Overlay && AnalysisHelper::IsTrueCC1Pi(pEvent, config.global.useAbsPdg))
            {
                o_truthData = AnalysisHelper::GetTruthAnalysisData(pEvent->truth, config.global.useAbsPdg, config.global.protonMomentumThreshold);
            }

            // BEGIN TEST
            const auto pion = pEvent->reco.particles.at(AnalysisHelper::GetParticleIndexWithPdg(assignedPdgCodes, 211));
            o_recoPionTrackScore = pion.trackScore();
            o_recoPionWiggliness = pion.wiggliness();
            // END TEST

            pTree->Fill();
        }
    }

    pFile->Write();
    pFile->Close();
}

} // namespace ubcc1pi_macros
