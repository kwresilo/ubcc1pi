// /**
//  *  @file  ubcc1pi_standalone/Macros/PlotProtonVariables.cxx
//  *
//  *  @brief The implementation file of the PlotProtonVariables macro
//  */

// #include "ubcc1pi_standalone/Macros/Macros.h"

// #include "ubcc1pi_standalone/Objects/FileReader.h"

// #include "ubcc1pi_standalone/Helpers/PlottingHelper.h"
// #include "ubcc1pi_standalone/Helpers/AnalysisHelper.h"
// #include "ubcc1pi_standalone/Helpers/SelectionHelper.h"
// #include "ubcc1pi_standalone/Helpers/NormalisationHelper.h"


// using namespace ubcc1pi;

// namespace ubcc1pi_macros
// {

// void PlotProtonVariables(const Config &config)
// {
//     //
//     // Setup the input files
//     //
//     // std::cout<<"PlotProtonVariables Point -1"<<std::endl;
//     std::vector< std::tuple<AnalysisHelper::SampleType, std::string, float> > inputData;

//     std::cout<<"##########################################\nUSING NUWRO AS DATA & Only CC0pi!\n##########################################"<<std::endl;
//     for (const auto run: config.global.runs)
//     {
//         if(run == 1)
//         {
//             inputData.emplace_back(AnalysisHelper::Overlay, config.filesRun1.overlaysFileName, NormalisationHelper::GetOverlaysNormalisationToNuWro(config, 1));
//             inputData.emplace_back(AnalysisHelper::DataBNB, config.filesRun1.nuWroFileName, 1.f);
//         }
//         else if(run == 2)
//         {
//             inputData.emplace_back(AnalysisHelper::Overlay, config.filesRun2.overlaysFileName, NormalisationHelper::GetOverlaysNormalisationToNuWro(config, 2));
//             inputData.emplace_back(AnalysisHelper::DataBNB, config.filesRun2.nuWroFileName, 1.f);
//         }
//         else if(run == 3)
//         {
//             inputData.emplace_back(AnalysisHelper::Overlay, config.filesRun3.overlaysFileName, NormalisationHelper::GetOverlaysNormalisationToNuWro(config, 3));
//             inputData.emplace_back(AnalysisHelper::DataBNB, config.filesRun3.nuWroFileName, 1.f);
//         }
//         else throw std::logic_error("PlotEventSelectionCuts - Invalid run number");
//     }


//     //
//     // Setup the plots
//     //
//     // std::cout<<"PlotProtonVariables Point 0"<<std::endl;
//     const std::string yLabel = "Number of particles";

//     // PlottingHelper::MultiPlot muonMomentumPlot("Muon momentum / GeV", yLabel, config.global.muonMomentum.binEdges);
//     // PlottingHelper::MultiPlot muonCosThetaPlot("Muon cos(theta)", yLabel, config.global.muonCosTheta.binEdges);
//     // PlottingHelper::MultiPlot muonPhiPlot("Muon phi / rad", yLabel, config.global.muonPhi.binEdges);

//     // PlottingHelper::MultiPlot muonMomentumParticlePlot("Muon momentum / GeV", yLabel, config.global.muonMomentum.binEdges);
//     // PlottingHelper::MultiPlot muonCosThetaParticlePlot("Muon cos(theta)", yLabel, config.global.muonCosTheta.binEdges);
//     // PlottingHelper::MultiPlot muonPhiParticlePlot("Muon phi / rad", yLabel, config.global.muonPhi.binEdges);

//     // PlottingHelper::MultiPlot pionMomentumPlot("Pion momentum / GeV", yLabel, config.global.pionMomentum.binEdges);
//     // PlottingHelper::MultiPlot pionCosThetaPlot("Pion cos(theta)", yLabel, config.global.pionCosTheta.binEdges);
//     // PlottingHelper::MultiPlot pionPhiPlot("Pion phi / rad", yLabel, config.global.pionPhi.binEdges);

//     // PlottingHelper::MultiPlot pionMomentumParticlePlot("Pion momentum / GeV", yLabel, config.global.pionMomentum.binEdges);
//     // PlottingHelper::MultiPlot pionCosThetaParticlePlot("Pion cos(theta)", yLabel, config.global.pionCosTheta.binEdges);
//     // PlottingHelper::MultiPlot pionPhiParticlePlot("Pion phi / rad", yLabel, config.global.pionPhi.binEdges);

//     // PlottingHelper::MultiPlot muonPionAnglePlot("Muon-pion opening angle / rad", yLabel, config.global.muonPionAngle.binEdges);
//     // PlottingHelper::MultiPlot nProtonsTruePlot("True proton multiplicity", yLabel, config.global.nProtons.binEdges);

//     auto pNProtonsTruePlotNuWroAll            = std::make_shared<TH1F>("True proton multiplicity nuwro all", "", 6,0,6);
//     auto pNProtonsTruePlotGenieAll            = std::make_shared<TH1F>("True proton multiplicity genie all", "", 6,0,6);
//     auto pNProtonsRecoPlotNuWroGeneric        = std::make_shared<TH1F>("Reco proton multiplicity nuwro generic", "", 6,0,6);
//     auto pNProtonsRecoPlotGenieGeneric        = std::make_shared<TH1F>("Reco proton multiplicity genie generic", "", 6,0,6);
//     auto pNProtonsRecoPlotNuWroSideband       = std::make_shared<TH1F>("Reco proton multiplicity nuwro sideband", "", 6,0,6);
//     auto pNProtonsRecoPlotGenieSideband       = std::make_shared<TH1F>("Reco proton multiplicity genie sideband", "", 6,0,6);
//     auto pNProtonsTruePlotNuWroGeneric        = std::make_shared<TH1F>("True proton multiplicity nuwro generic", "", 6,0,6);
//     auto pNProtonsTruePlotGenieGeneric        = std::make_shared<TH1F>("True proton multiplicity genie generic", "", 6,0,6);
//     auto pNProtonsTruePlotNuWroSideband       = std::make_shared<TH1F>("True proton multiplicity nuwro sideband", "", 6,0,6);
//     auto pNProtonsTruePlotGenieSideband       = std::make_shared<TH1F>("True proton multiplicity genie sideband", "", 6,0,6);


//     auto pProtonMomentumTruePlotNuWroAll      = std::make_shared<TH1F>("True proton momentum nuwro all",      "", 10,0.f,2.5f);
//     auto pProtonMomentumTruePlotGenieAll      = std::make_shared<TH1F>("True proton momentum genie all",      "", 10,0.f,2.5f);
//     auto pProtonMomentumRecoPlotNuWroGeneric  = std::make_shared<TH1F>("Reco proton momentum nuwro generic",  "", 10,0.f,2.5f);
//     auto pProtonMomentumRecoPlotGenieGeneric  = std::make_shared<TH1F>("Reco proton momentum genie generic",  "", 10,0.f,2.5f);
//     auto pProtonMomentumRecoPlotNuWroSideband = std::make_shared<TH1F>("Reco proton momentum nuwro sideband", "", 10,0.f,2.5f);
//     auto pProtonMomentumRecoPlotGenieSideband = std::make_shared<TH1F>("Reco proton momentum genie sideband", "", 10,0.f,2.5f);
//     auto pProtonMomentumTruePlotNuWroGeneric  = std::make_shared<TH1F>("True proton momentum nuwro generic",  "", 10,0.f,2.5f);
//     auto pProtonMomentumTruePlotGenieGeneric  = std::make_shared<TH1F>("True proton momentum genie generic",  "", 10,0.f,2.5f);
//     auto pProtonMomentumTruePlotNuWroSideband = std::make_shared<TH1F>("True proton momentum nuwro sideband", "", 10,0.f,2.5f);
//     auto pProtonMomentumTruePlotGenieSideband = std::make_shared<TH1F>("True proton momentum genie sideband", "", 10,0.f,2.5f);

//     auto pProtonCosThetaTruePlotNuWroAll      = std::make_shared<TH1F>("True proton cos(theta) nuwro all",      "", 10,-1.f,1.f);
//     auto pProtonCosThetaTruePlotGenieAll      = std::make_shared<TH1F>("True proton cos(theta) genie all",      "", 10,-1.f,1.f);
//     auto pProtonCosThetaRecoPlotNuWroGeneric  = std::make_shared<TH1F>("Reco proton cos(theta) nuwro generic",  "", 10,-1.f,1.f);
//     auto pProtonCosThetaRecoPlotGenieGeneric  = std::make_shared<TH1F>("Reco proton cos(theta) genie generic",  "", 10,-1.f,1.f);
//     auto pProtonCosThetaRecoPlotNuWroSideband = std::make_shared<TH1F>("Reco proton cos(theta) nuwro sideband", "", 10,-1.f,1.f);
//     auto pProtonCosThetaRecoPlotGenieSideband = std::make_shared<TH1F>("Reco proton cos(theta) genie sideband", "", 10,-1.f,1.f);
//     auto pProtonCosThetaTruePlotNuWroGeneric  = std::make_shared<TH1F>("True proton cos(theta) nuwro generic",  "", 10,-1.f,1.f);
//     auto pProtonCosThetaTruePlotGenieGeneric  = std::make_shared<TH1F>("True proton cos(theta) genie generic",  "", 10,-1.f,1.f);
//     auto pProtonCosThetaTruePlotNuWroSideband = std::make_shared<TH1F>("True proton cos(theta) nuwro sideband", "", 10,-1.f,1.f);
//     auto pProtonCosThetaTruePlotGenieSideband = std::make_shared<TH1F>("True proton cos(theta) genie sideband", "", 10,-1.f,1.f);

//     auto pProtonPhiTruePlotNuWroAll           = std::make_shared<TH1F>("True proton phi nuwro all",      "", 10,-3.142f,3.142f);
//     auto pProtonPhiTruePlotGenieAll           = std::make_shared<TH1F>("True proton phi genie all",      "", 10,-3.142f,3.142f);
//     auto pProtonPhiRecoPlotNuWroGeneric       = std::make_shared<TH1F>("Reco proton phi nuwro generic",  "", 10,-3.142f,3.142f);
//     auto pProtonPhiRecoPlotGenieGeneric       = std::make_shared<TH1F>("Reco proton phi genie generic",  "", 10,-3.142f,3.142f);
//     auto pProtonPhiRecoPlotNuWroSideband      = std::make_shared<TH1F>("Reco proton phi nuwro sideband", "", 10,-3.142f,3.142f);
//     auto pProtonPhiRecoPlotGenieSideband      = std::make_shared<TH1F>("Reco proton phi genie sideband", "", 10,-3.142f,3.142f);
//     auto pProtonPhiTruePlotNuWroGeneric       = std::make_shared<TH1F>("True proton phi nuwro generic",  "", 10,-3.142f,3.142f);
//     auto pProtonPhiTruePlotGenieGeneric       = std::make_shared<TH1F>("True proton phi genie generic",  "", 10,-3.142f,3.142f);
//     auto pProtonPhiTruePlotNuWroSideband      = std::make_shared<TH1F>("True proton phi nuwro sideband", "", 10,-3.142f,3.142f);
//     auto pProtonPhiTruePlotGenieSideband      = std::make_shared<TH1F>("True proton phi genie sideband", "", 10,-3.142f,3.142f);

//     auto pMuonMomentumTruePlotNuWroAll        = std::make_shared<TH1F>("True muon momentum nuwro all",      "", 10,0.f,2.f);
//     auto pMuonMomentumTruePlotGenieAll        = std::make_shared<TH1F>("True muon momentum genie all",      "", 10,0.f,2.f);
//     auto pMuonMomentumRecoPlotNuWroGeneric    = std::make_shared<TH1F>("Reco muon momentum nuwro generic",  "", 10,0.f,2.f);
//     auto pMuonMomentumRecoPlotGenieGeneric    = std::make_shared<TH1F>("Reco muon momentum genie generic",  "", 10,0.f,2.f);
//     auto pMuonMomentumRecoPlotNuWroSideband   = std::make_shared<TH1F>("Reco muon momentum nuwro sideband", "", 10,0.f,2.f);
//     auto pMuonMomentumRecoPlotGenieSideband   = std::make_shared<TH1F>("Reco muon momentum genie sideband", "", 10,0.f,2.f);
//     auto pMuonMomentumTruePlotNuWroGeneric    = std::make_shared<TH1F>("True muon momentum nuwro generic",  "", 10,0.f,2.f);
//     auto pMuonMomentumTruePlotGenieGeneric    = std::make_shared<TH1F>("True muon momentum genie generic",  "", 10,0.f,2.f);
//     auto pMuonMomentumTruePlotNuWroSideband   = std::make_shared<TH1F>("True muon momentum nuwro sideband", "", 10,0.f,2.f);
//     auto pMuonMomentumTruePlotGenieSideband   = std::make_shared<TH1F>("True muon momentum genie sideband", "", 10,0.f,2.f);

//     auto pMuonCosThetaTruePlotNuWroAll        = std::make_shared<TH1F>("True muon cos(theta) nuwro all",      "", 10,-1.f,1.f);
//     auto pMuonCosThetaTruePlotGenieAll        = std::make_shared<TH1F>("True muon cos(theta) genie all",      "", 10,-1.f,1.f);
//     auto pMuonCosThetaRecoPlotNuWroGeneric    = std::make_shared<TH1F>("Reco muon cos(theta) nuwro generic",  "", 10,-1.f,1.f);
//     auto pMuonCosThetaRecoPlotGenieGeneric    = std::make_shared<TH1F>("Reco muon cos(theta) genie generic",  "", 10,-1.f,1.f);
//     auto pMuonCosThetaRecoPlotNuWroSideband   = std::make_shared<TH1F>("Reco muon cos(theta) nuwro sideband", "", 10,-1.f,1.f);
//     auto pMuonCosThetaRecoPlotGenieSideband   = std::make_shared<TH1F>("Reco muon cos(theta) genie sideband", "", 10,-1.f,1.f);
//     auto pMuonCosThetaTruePlotNuWroGeneric    = std::make_shared<TH1F>("True muon cos(theta) nuwro generic",  "", 10,-1.f,1.f);
//     auto pMuonCosThetaTruePlotGenieGeneric    = std::make_shared<TH1F>("True muon cos(theta) genie generic",  "", 10,-1.f,1.f);
//     auto pMuonCosThetaTruePlotNuWroSideband   = std::make_shared<TH1F>("True muon cos(theta) nuwro sideband", "", 10,-1.f,1.f);
//     auto pMuonCosThetaTruePlotGenieSideband   = std::make_shared<TH1F>("True muon cos(theta) genie sideband", "", 10,-1.f,1.f);

//     auto pMuonPhiTruePlotNuWroAll             = std::make_shared<TH1F>("True muon phi nuwro all",      "", 10,-3.142f,3.142f);
//     auto pMuonPhiTruePlotGenieAll             = std::make_shared<TH1F>("True muon phi genie all",      "", 10,-3.142f,3.142f);
//     auto pMuonPhiRecoPlotNuWroGeneric         = std::make_shared<TH1F>("Reco muon phi nuwro generic",  "", 10,-3.142f,3.142f);
//     auto pMuonPhiRecoPlotGenieGeneric         = std::make_shared<TH1F>("Reco muon phi genie generic",  "", 10,-3.142f,3.142f);
//     auto pMuonPhiRecoPlotNuWroSideband        = std::make_shared<TH1F>("Reco muon phi nuwro sideband", "", 10,-3.142f,3.142f);
//     auto pMuonPhiRecoPlotGenieSideband        = std::make_shared<TH1F>("Reco muon phi genie sideband", "", 10,-3.142f,3.142f);
//     auto pMuonPhiTruePlotNuWroGeneric         = std::make_shared<TH1F>("True muon phi nuwro generic",  "", 10,-3.142f,3.142f);
//     auto pMuonPhiTruePlotGenieGeneric         = std::make_shared<TH1F>("True muon phi genie generic",  "", 10,-3.142f,3.142f);
//     auto pMuonPhiTruePlotNuWroSideband        = std::make_shared<TH1F>("True muon phi nuwro sideband", "", 10,-3.142f,3.142f);
//     auto pMuonPhiTruePlotGenieSideband        = std::make_shared<TH1F>("True muon phi genie sideband", "", 10,-3.142f,3.142f);

//     // PlottingHelper::MultiPlot muonMomentumPlot("Muon momentum / GeV", yLabel, 50u, 0.f, 2.f, true, config.global.axisTitles);
//     // PlottingHelper::MultiPlot muonCosThetaPlot("Muon cos(theta)", yLabel, 50u, config.global.muonCosTheta.min, config.global.muonCosTheta.max, true, config.global.axisTitles);
//     // PlottingHelper::MultiPlot muonPhiPlot("Muon phi / rad", yLabel, 50u, config.global.muonPhi.min, config.global.muonPhi.max, true, config.global.axisTitles);

//     // PlottingHelper::MultiPlot muonMomentumParticlePlot("Muon momentum / GeV", yLabel, 50u, 0.f, 2.f, true, config.global.axisTitles);
//     // PlottingHelper::MultiPlot muonCosThetaParticlePlot("Muon cos(theta)", yLabel, 50u, config.global.muonCosTheta.min, config.global.muonCosTheta.max, true, config.global.axisTitles);
//     // PlottingHelper::MultiPlot muonPhiParticlePlot("Muon phi / rad", yLabel, 50u, config.global.muonPhi.min, config.global.muonPhi.max, true, config.global.axisTitles);

//     // PlottingHelper::MultiPlot pionMomentumPlot("Pion momentum / GeV", yLabel, 50u, 0.f, 0.8f, true, config.global.axisTitles);
//     // PlottingHelper::MultiPlot pionCosThetaPlot("Pion cos(theta)", yLabel, 50u, config.global.pionCosTheta.min, config.global.pionCosTheta.max, true, config.global.axisTitles);
//     // PlottingHelper::MultiPlot pionPhiPlot("Pion phi / rad", yLabel, 50u, config.global.pionPhi.min, config.global.pionPhi.max, true, config.global.axisTitles);

//     // PlottingHelper::MultiPlot pionMomentumParticlePlot("Pion momentum / GeV", yLabel, 50u, 0.f, 0.8f, true, config.global.axisTitles);
//     // PlottingHelper::MultiPlot pionCosThetaParticlePlot("Pion cos(theta)", yLabel, 50u, config.global.pionCosTheta.min, config.global.pionCosTheta.max, true, config.global.axisTitles);
//     // PlottingHelper::MultiPlot pionPhiParticlePlot("Pion phi / rad", yLabel, 50u, config.global.pionPhi.min, config.global.pionPhi.max, true, config.global.axisTitles);

//     // PlottingHelper::MultiPlot muonPionAnglePlot("Muon-pion opening angle / rad", yLabel, 50u, config.global.muonPionAngle.min, config.global.muonPionAngle.max, true, config.global.axisTitles);
//     // PlottingHelper::MultiPlot nProtonsPlot("Proton multiplicity", yLabel, 5u, 0, 5, true, config.global.axisTitles);

//     // The highest energy proton variables are plotted with three different selections
//     const std::vector<std::string> protonPlotNames{"nProtons>0", "nProtons==1", "nProtons>=2"};
//     std::vector<PlottingHelper::MultiPlot> protonMomentumPlots, protonCosThetaPlots, protonPhiPlots;
//     std::vector<PlottingHelper::MultiPlot> protonMomentumParticlePlots, protonCosThetaParticlePlots, protonPhiParticlePlots;
//     std::vector<PlottingHelper::MultiPlot> protonPionAnglePlots, protonMuonAnglePlots;

//     std::vector<PlottingHelper::MultiPlot> protonMomentumPlotsTruth, protonCosThetaPlotsTruth, protonPhiPlotsTruth;
//     std::vector<PlottingHelper::MultiPlot> protonMomentumParticlePlotsTruth, protonCosThetaParticlePlotsTruth, protonPhiParticlePlotsTruth;
//     std::vector<PlottingHelper::MultiPlot> protonMuonAnglePlotsTruth;// protonPionAnglePlotsTruth;


//     std::vector<PlottingHelper::MultiPlot> protonMomentumPlotsSideband, protonCosThetaPlotsSideband, protonPhiPlotsSideband;
//     std::vector<PlottingHelper::MultiPlot> protonMomentumParticlePlotsSideband, protonCosThetaParticlePlotsSideband, protonPhiParticlePlotsSideband;
//     std::vector<PlottingHelper::MultiPlot> protonPionAnglePlotsSideband, protonMuonAnglePlotsSideband;

//     std::vector<PlottingHelper::MultiPlot> protonMomentumPlotsSidebandTruth, protonCosThetaPlotsSidebandTruth, protonPhiPlotsSidebandTruth;
//     std::vector<PlottingHelper::MultiPlot> protonMomentumParticlePlotsSidebandTruth, protonCosThetaParticlePlotsSidebandTruth, protonPhiParticlePlotsSidebandTruth;
//     std::vector<PlottingHelper::MultiPlot> protonPionAnglePlotsSidebandTruth, protonMuonAnglePlotsSidebandTruth;

//     std::shared_ptr<TH2F> pProtonMomvsPhiAllNuWro(new TH2F("protonMomVsPhiAllNuWro", ";Proton phi / rad;Proton momentum / GeV", 20, -3.142f, 3.142f, 20, 0.1f, 2.2f));
//     std::shared_ptr<TH2F> pProtonMomvsPhiGenericNuWro(new TH2F("protonMomVsPhiGenericNuWro", ";Proton phi / rad;Proton momentum / GeV", 20, -3.142f, 3.142f, 20, 0.1f, 2.2f));
//     std::shared_ptr<TH2F> pProtonMomvsPhiSidebandNuWro(new TH2F("protonMomVsPhiSidebandNuWro", ";Proton phi / rad;Proton momentum / GeV", 20, -3.142f, 3.142f, 20, 0.1f, 2.2f));
//     std::shared_ptr<TH2F> pProtonMomvsPhiAllGenie(new TH2F("protonMomVsPhiAllGenie", ";Proton phi / rad;Proton momentum / GeV", 20, -3.142f, 3.142f, 20, 0.1f, 2.2f));
//     std::shared_ptr<TH2F> pProtonMomvsPhiGenericGenie(new TH2F("protonMomVsPhiGenericGenie", ";Proton phi / rad;Proton momentum / GeV", 20, -3.142f, 3.142f, 20, 0.1f, 2.2f));
//     std::shared_ptr<TH2F> pProtonMomvsPhiSidebandGenie(new TH2F("protonMomVsPhiSidebandGenie", ";Proton phi / rad;Proton momentum / GeV", 20, -3.142f, 3.142f, 20, 0.1f, 2.2f));

//     std::shared_ptr<TH2F> pProtonMomvsCosThetaAllNuWro(new TH2F("protonMomVsCosThetaAllNuWro", ";Cos(theta);Proton momentum / GeV", 20, -1.f, 1.f, 20, 0.1f, 2.2f));
//     std::shared_ptr<TH2F> pProtonMomvsCosThetaGenericNuWro(new TH2F("protonMomVsCosThetaGenericNuWro", ";Cos(theta);Proton momentum / GeV", 20, -1.f, 1.f, 20, 0.1f, 2.2f));
//     std::shared_ptr<TH2F> pProtonMomvsCosThetaSidebandNuWro(new TH2F("protonMomVsCosThetaSidebandNuWro", ";Cos(theta);Proton momentum / GeV", 20, -1.f, 1.f, 20, 0.1f, 2.2f));
//     std::shared_ptr<TH2F> pProtonMomvsCosThetaAllGenie(new TH2F("protonMomVsCosThetaAllGenie", ";Cos(theta);Proton momentum / GeV", 20, -1.f, 1.f, 20, 0.1f, 2.2f));
//     std::shared_ptr<TH2F> pProtonMomvsCosThetaGenericGenie(new TH2F("protonMomVsCosThetaGenericGenie", ";Cos(theta);Proton momentum / GeV", 20, -1.f, 1.f, 20, 0.1f, 2.2f));
//     std::shared_ptr<TH2F> pProtonMomvsCosThetaSidebandGenie(new TH2F("protonMomVsCosThetaSidebandGenie", ";Cos(theta);Proton momentum / GeV", 20, -1.f, 1.f, 20, 0.1f, 2.2f));

//     std::shared_ptr<TH2F> pProtonPhivsCosThetaAllNuWro(new TH2F("protonPhiVsCosThetaAllNuWro", ";Cos(theta);Proton phi / rad", 20, -1.f, 1.f, 20, -3.142f, 3.142f));
//     std::shared_ptr<TH2F> pProtonPhivsCosThetaGenericNuWro(new TH2F("protonPhiVsCosThetaGenericNuWro", ";Cos(theta);Proton phi / rad", 20, -1.f, 1.f, 20, -3.142f, 3.142f));
//     std::shared_ptr<TH2F> pProtonPhivsCosThetaSidebandNuWro(new TH2F("protonPhiVsCosThetaSidebandNuWro", ";Cos(theta);Proton phi / rad", 20, -1.f, 1.f, 20, -3.142f, 3.142f));
//     std::shared_ptr<TH2F> pProtonPhivsCosThetaAllGenie(new TH2F("protonPhiVsCosThetaAllGenie", ";Cos(theta);Proton phi / rad", 20, -1.f, 1.f, 20, -3.142f, 3.142f));
//     std::shared_ptr<TH2F> pProtonPhivsCosThetaGenericGenie(new TH2F("protonPhiVsCosThetaGenericGenie", ";Cos(theta);Proton phi / rad", 20, -1.f, 1.f, 20, -3.142f, 3.142f));
//     std::shared_ptr<TH2F> pProtonPhivsCosThetaSidebandGenie(new TH2F("protonPhiVsCosThetaSidebandGenie", ";Cos(theta);Proton phi / rad", 20, -1.f, 1.f, 20, -3.142f, 3.142f));

//     std::shared_ptr<TH2F> pTrueVsRecoProtonMultiplicityAllGenericNuWro(new TH2F("trueVsRecoProtonMultiplicityAllGenericNuWro", ";Reco proton multiplicity;True proton multiplicity", 6, 0, 6, 6, 0, 6));
//     std::shared_ptr<TH2F> pTrueVsRecoProtonMultiplicityAllSidebandNuWro(new TH2F("trueVsRecoProtonMultiplicityAllSidebandNuWro", ";Reco proton multiplicity;True proton multiplicity", 6, 0, 6, 6, 0, 6));
//     std::shared_ptr<TH2F> pTrueVsRecoProtonMultiplicityAllGenericGenie(new TH2F("trueVsRecoProtonMultiplicityAllGenericGenie", ";Reco proton multiplicity;True proton multiplicity", 6, 0, 6, 6, 0, 6));
//     std::shared_ptr<TH2F> pTrueVsRecoProtonMultiplicityAllSidebandGenie(new TH2F("trueVsRecoProtonMultiplicityAllSidebandGenie", ";Reco proton multiplicity;True proton multiplicity", 6, 0, 6, 6, 0, 6));
//     std::shared_ptr<TH2F> pTrueVsRecoProtonMultiplicityGenericNuWro(new TH2F("trueVsRecoProtonMultiplicityGenericNuWro", ";Reco proton multiplicity;True proton multiplicity", 6, 0, 6, 6, 0, 6));
//     std::shared_ptr<TH2F> pTrueVsRecoProtonMultiplicitySidebandNuWro(new TH2F("trueVsRecoProtonMultiplicitySidebandNuWro", ";Reco proton multiplicity;True proton multiplicity", 6, 0, 6, 6, 0, 6));
//     std::shared_ptr<TH2F> pTrueVsRecoProtonMultiplicityGenericGenie(new TH2F("trueVsRecoProtonMultiplicityGenericGenie", ";Reco proton multiplicity;True proton multiplicity", 6, 0, 6, 6, 0, 6));
//     std::shared_ptr<TH2F> pTrueVsRecoProtonMultiplicitySidebandGenie(new TH2F("trueVsRecoProtonMultiplicitySidebandGenie", ";Reco proton multiplicity;True proton multiplicity", 6, 0, 6, 6, 0, 6));

//     std::shared_ptr<TH2F> pTrueProtonMomvsTrueProtonMultiplicityAllNuWro(new TH2F("trueProtonMomvsTrueProtonMultiplicityAllNuWro", ";True proton multiplicity;True proton momentum / GeV", 10, 0, 10, 20, 0.1f, 2.2f));
//     std::shared_ptr<TH2F> pTrueProtonMomvsTrueProtonMultiplicityGenericNuWro(new TH2F("trueProtonMomvsTrueProtonMultiplicityGenericNuWro", ";True proton multiplicity;True proton momentum / GeV", 10, 0, 10, 20, 0.1f, 2.2f));
//     std::shared_ptr<TH2F> pTrueProtonMomvsTrueProtonMultiplicitySidebandNuWro(new TH2F("trueProtonMomvsTrueProtonMultiplicitySidebandNuWro", ";True proton multiplicity;True proton momentum / GeV", 10, 0, 10, 20, 0.1f, 2.2f));
//     std::shared_ptr<TH2F> pTrueProtonMomvsTrueProtonMultiplicityNotGenericNuWro(new TH2F("trueProtonMomvsTrueProtonMultiplicityNotGenericNuWro", ";True proton multiplicity;True proton momentum / GeV", 10, 0, 10, 20, 0.1f, 2.2f));
//     std::shared_ptr<TH2F> pTrueProtonMomvsTrueProtonMultiplicityNotSidebandNuWro(new TH2F("trueProtonMomvsTrueProtonMultiplicityNotSidebandNuWro", ";True proton multiplicity;True proton momentum / GeV", 10, 0, 10, 20, 0.1f, 2.2f));
//     std::shared_ptr<TH2F> pTrueProtonMomvsTrueProtonMultiplicityAllGenie(new TH2F("trueProtonMomvsTrueProtonMultiplicityAllGenie", ";True proton multiplicity;True proton momentum / GeV", 10, 0, 10, 20, 0.1f, 2.2f));
//     std::shared_ptr<TH2F> pTrueProtonMomvsTrueProtonMultiplicityGenericGenie(new TH2F("trueProtonMomvsTrueProtonMultiplicityGenericGenie", ";True proton multiplicity;True proton momentum / GeV", 10, 0, 10, 20, 0.1f, 2.2f));
//     std::shared_ptr<TH2F> pTrueProtonMomvsTrueProtonMultiplicitySidebandGenie(new TH2F("trueProtonMomvsTrueProtonMultiplicitySidebandGenie", ";True proton multiplicity;True proton momentum / GeV", 10, 0, 10, 20, 0.1f, 2.2f));
//     std::shared_ptr<TH2F> pTrueProtonMomvsTrueProtonMultiplicityNotGenericGenie(new TH2F("trueProtonMomvsTrueProtonMultiplicityNotGenericGenie", ";True proton multiplicity;True proton momentum / GeV", 10, 0, 10, 20, 0.1f, 2.2f));
//     std::shared_ptr<TH2F> pTrueProtonMomvsTrueProtonMultiplicityNotSidebandGenie(new TH2F("trueProtonMomvsTrueProtonMultiplicityNotSidebandGenie", ";True proton multiplicity;True proton momentum / GeV", 10, 0, 10, 20, 0.1f, 2.2f));

//     std::shared_ptr<TH2F> pTrueProtonMomvsTrueProtonMultiplicityGenericNuWroRatio(new TH2F("trueProtonMomvsTrueProtonMultiplicityGenericNuWroRatio", ";True proton multiplicity;True proton momentum / GeV", 10, 0, 10, 20, 0.1f, 2.2f));
//     std::shared_ptr<TH2F> pTrueProtonMomvsTrueProtonMultiplicitySidebandNuWroRatio(new TH2F("trueProtonMomvsTrueProtonMultiplicitySidebandNuWroRatio", ";True proton multiplicity;True proton momentum / GeV", 10, 0, 10, 20, 0.1f, 2.2f));
//     std::shared_ptr<TH2F> pTrueProtonMomvsTrueProtonMultiplicityGenericGenieRatio(new TH2F("trueProtonMomvsTrueProtonMultiplicityGenericGenieRatio", ";True proton multiplicity;True proton momentum / GeV", 10, 0, 10, 20, 0.1f, 2.2f));
//     std::shared_ptr<TH2F> pTrueProtonMomvsTrueProtonMultiplicitySidebandGenieRatio(new TH2F("trueProtonMomvsTrueProtonMultiplicitySidebandGenieRatio", ";True proton multiplicity;True proton momentum / GeV", 10, 0, 10, 20, 0.1f, 2.2f));


//     // Proton plots are generated for three different muliplicity ranges
//     // const std::vector<std::pair<std::string, std::function<bool(unsigned int)>> {
//     //     {"protons>0", [](unsigned int n){return n>0;}},
//     //     {"protons=1", [](unsigned int n){return n=1;}},
//     //     {"protons>=2", [](unsigned int n){return n>=2;}}
//     // };

//     // std::cout<<"PlotProtonVariables Point 1"<<std::endl;

//     for(auto const& plotName : protonPlotNames)
//     {
//         //TODO: Some of the values here are placeholders - make sure they are sensible
//         protonMomentumPlots.emplace_back("Proton momentum / GeV (" + plotName + ")", yLabel, 20u, config.global.protonMomentum.min, config.global.protonMomentum.max, true, config.global.axisTitles);
//         protonCosThetaPlots.emplace_back("Proton cos(theta) (" + plotName + ")", yLabel, 10u, config.global.protonCosTheta.min, config.global.protonCosTheta.max, true, config.global.axisTitles);
//         protonPhiPlots.emplace_back("Proton phi / rad (" + plotName + ")", yLabel, 10u, config.global.protonPhi.min, config.global.protonPhi.max, true, config.global.axisTitles);

//         protonMomentumParticlePlots.emplace_back("Proton momentum / GeV (" + plotName + ")", yLabel, 20u, config.global.protonMomentum.min, config.global.protonMomentum.max, true, config.global.axisTitles);
//         protonCosThetaParticlePlots.emplace_back("Proton cos(theta) (" + plotName + ")", yLabel, 10u, config.global.protonCosTheta.min, config.global.protonCosTheta.max, true, config.global.axisTitles);
//         protonPhiParticlePlots.emplace_back("Proton phi / rad (" + plotName + ")", yLabel, 10u, config.global.protonPhi.min, config.global.protonPhi.max, true, config.global.axisTitles);

//         protonPionAnglePlots.emplace_back("Proton-pion opening angle / rad (" + plotName + ")", yLabel, 10u, config.global.protonPionAngle.min, config.global.protonPionAngle.max, true, config.global.axisTitles);
//         protonMuonAnglePlots.emplace_back("Proton-muon opening angle / rad (" + plotName + ")", yLabel, 10u, config.global.protonMuonAngle.min, config.global.protonMuonAngle.max, true, config.global.axisTitles);


//         protonMomentumPlotsSideband.emplace_back("Sideband Proton momentum / GeV (" + plotName + ")", yLabel, 20u, config.global.protonMomentum.min, config.global.protonMomentum.max, true, config.global.axisTitles);
//         protonCosThetaPlotsSideband.emplace_back("Sideband Proton cos(theta) (" + plotName + ")", yLabel, 10u, config.global.protonCosTheta.min, config.global.protonCosTheta.max, true, config.global.axisTitles);
//         protonPhiPlotsSideband.emplace_back("Sideband Proton phi / rad (" + plotName + ")", yLabel, 10u, config.global.protonPhi.min, config.global.protonPhi.max, true, config.global.axisTitles);

//         protonMomentumParticlePlotsSideband.emplace_back("Sideband Proton momentum / GeV (" + plotName + ")", yLabel, 20u, config.global.protonMomentum.min, config.global.protonMomentum.max, true, config.global.axisTitles);
//         protonCosThetaParticlePlotsSideband.emplace_back("Sideband Proton cos(theta) (" + plotName + ")", yLabel, 10u, config.global.protonCosTheta.min, config.global.protonCosTheta.max, true, config.global.axisTitles);
//         protonPhiParticlePlotsSideband.emplace_back("Sideband Proton phi / rad (" + plotName + ")", yLabel, 10u, config.global.protonPhi.min, config.global.protonPhi.max, true, config.global.axisTitles);

//         protonPionAnglePlotsSideband.emplace_back("Sideband Proton-pion opening angle / rad (" + plotName + ")", yLabel, 10u, config.global.protonPionAngle.min, config.global.protonPionAngle.max, true, config.global.axisTitles);
//         protonMuonAnglePlotsSideband.emplace_back("Sideband Proton-muon opening angle / rad (" + plotName + ")", yLabel, 10u, config.global.protonMuonAngle.min, config.global.protonMuonAngle.max, true, config.global.axisTitles);



//         protonMomentumPlotsTruth.emplace_back("True proton momentum / GeV (" + plotName + ")", yLabel, 20u, config.global.protonMomentum.min, config.global.protonMomentum.max, true, config.global.axisTitles);
//         protonCosThetaPlotsTruth.emplace_back("True proton cos(theta) (" + plotName + ")", yLabel, 10u, config.global.protonCosTheta.min, config.global.protonCosTheta.max, true, config.global.axisTitles);
//         protonPhiPlotsTruth.emplace_back("True proton phi / rad (" + plotName + ")", yLabel, 10u, config.global.protonPhi.min, config.global.protonPhi.max, true, config.global.axisTitles);

//         protonMomentumParticlePlotsTruth.emplace_back("True proton momentum / GeV (" + plotName + ")", yLabel, 20u, config.global.protonMomentum.min, config.global.protonMomentum.max, true, config.global.axisTitles);
//         protonCosThetaParticlePlotsTruth.emplace_back("True proton cos(theta) (" + plotName + ")", yLabel, 10u, config.global.protonCosTheta.min, config.global.protonCosTheta.max, true, config.global.axisTitles);
//         protonPhiParticlePlotsTruth.emplace_back("True proton phi / rad (" + plotName + ")", yLabel, 10u, config.global.protonPhi.min, config.global.protonPhi.max, true, config.global.axisTitles);

//         // protonPionAnglePlotsTruth.emplace_back("True proton-pion opening angle / rad (" + plotName + ")", yLabel, 10u, config.global.protonPionAngle.min, config.global.protonPionAngle.max, true, config.global.axisTitles);
//         protonMuonAnglePlotsTruth.emplace_back("True proton-muon opening angle / rad (" + plotName + ")", yLabel, 10u, config.global.protonMuonAngle.min, config.global.protonMuonAngle.max, true, config.global.axisTitles);


//         protonMomentumPlotsSidebandTruth.emplace_back("True sideband Proton momentum / GeV (" + plotName + ")", yLabel, 20u, config.global.protonMomentum.min, config.global.protonMomentum.max, true, config.global.axisTitles);
//         protonCosThetaPlotsSidebandTruth.emplace_back("True sideband Proton cos(theta) (" + plotName + ")", yLabel, 10u, config.global.protonCosTheta.min, config.global.protonCosTheta.max, true, config.global.axisTitles);
//         protonPhiPlotsSidebandTruth.emplace_back("True sideband Proton phi / rad (" + plotName + ")", yLabel, 10u, config.global.protonPhi.min, config.global.protonPhi.max, true, config.global.axisTitles);

//         protonMomentumParticlePlotsSidebandTruth.emplace_back("True sideband Proton momentum / GeV (" + plotName + ")", yLabel, 20u, config.global.protonMomentum.min, config.global.protonMomentum.max, true, config.global.axisTitles);
//         protonCosThetaParticlePlotsSidebandTruth.emplace_back("True sideband Proton cos(theta) (" + plotName + ")", yLabel, 10u, config.global.protonCosTheta.min, config.global.protonCosTheta.max, true, config.global.axisTitles);
//         protonPhiParticlePlotsSidebandTruth.emplace_back("True sideband Proton phi / rad (" + plotName + ")", yLabel, 10u, config.global.protonPhi.min, config.global.protonPhi.max, true, config.global.axisTitles);

//         protonPionAnglePlotsSidebandTruth.emplace_back("True sideband Proton-pion opening angle / rad (" + plotName + ")", yLabel, 10u, config.global.protonPionAngle.min, config.global.protonPionAngle.max, true, config.global.axisTitles);
//         protonMuonAnglePlotsSidebandTruth.emplace_back("True sideband Proton-muon opening angle / rad (" + plotName + ")", yLabel, 10u, config.global.protonMuonAngle.min, config.global.protonMuonAngle.max, true, config.global.axisTitles);
//     }

//     //
//     // Get the selection
//     //
//     auto selection = SelectionHelper::GetDefaultSelection();

//     std::cout<<"..........................................\nUSING Modified CC0pi Selection: muonLikeProtonValue=-0.35f, barelyResemblingProtonValue=0.45f\n.........................................."<<std::endl;
//     auto sidebandSelection = SelectionHelper::GetCC0piSelectionModified(-0.35f, 0.45f);
//     // auto sidebandSelection = SelectionHelper::GetCC0piSelection();

//     // Loop over the events
//     for (const auto [sampleType, fileName, normalisation] : inputData)
//     {
//         std::cout << "Reading input file: " << fileName << std::endl;

//         FileReader reader(fileName);
//         auto pEvent = reader.GetBoundEventAddress();

//         const auto nEvents = reader.GetNumberOfEvents();
//         //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//         //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//         //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//         // std::cout<<"\n##############\nOnly counting every event!\n##############"<<std::endl;
//         for (unsigned int i = 0;i < nEvents;++i)//todo remove !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//         {
//             AnalysisHelper::PrintLoadingBar(i, nEvents);
//             reader.LoadEvent(i);

//             const auto isTrueCC0Pi = AnalysisHelper::IsTrueCC0Pi(pEvent, config.global.useAbsPdg, config.global.protonMomentumThreshold);// todo remove this
//             if(!isTrueCC0Pi) continue;

//             // Run the event selection and store which cuts are passed
//             const auto &[isSelectedGolden, cutsPassed, assignedPdgCodes] = selection.Execute(pEvent);
//             const auto isSelectedGeneric = (std::find(cutsPassed.begin(), cutsPassed.end(), config.global.lastCutGeneric) != cutsPassed.end());

//             const auto &[passedSidebandSelection, sidebandCutsPassed, sidebandAssignedPdgCodes] = sidebandSelection.Execute(pEvent);

//             const auto recoData = (
//                 isSelectedGeneric
//                     ? AnalysisHelper::GetRecoAnalysisData(pEvent->reco, assignedPdgCodes, isSelectedGolden)
//                     : AnalysisHelper::GetDummyAnalysisData()
//             );

//             const auto sidebandRecoData = (
//                 passedSidebandSelection
//                     ? AnalysisHelper::GetRecoAnalysisDataCC0Pi(pEvent->reco, sidebandAssignedPdgCodes)
//                     : AnalysisHelper::GetDummyAnalysisData()
//             );

//             const auto truthData = (
//                 (isTrueCC0Pi)
//                     ? AnalysisHelper::GetTruthAnalysisDataCC0Pi(pEvent->truth, config.global.useAbsPdg, config.global.protonMomentumThreshold)
//                     : AnalysisHelper::GetDummyAnalysisData()
//             );

//             // std::cout<<"assignedPdgCodes:"<<std::endl;
//             // for(const auto &v: assignedPdgCodes) std::cout<<v<<" ";
//             // std::cout<<std::endl;
//             // std::cout<<"sidebandAssignedPdgCodes:"<<std::endl;
//             // for(const auto &v: sidebandAssignedPdgCodes) std::cout<<v<<" ";
//             // std::cout<<std::endl;

//             // Get the truth and reco analysis data
//             const auto weight = AnalysisHelper::GetNominalEventWeight(pEvent) * normalisation;

//             // const auto truthData = AnalysisHelper::GetTruthAnalysisData(pEvent->truth, config.global.useAbsPdg, config.global.protonMomentumThreshold);
//             const auto plotStyle = PlottingHelper::GetPlotStyle(sampleType, pEvent, config.global.useAbsPdg);

//             // Get the true origin of the selected muon and pion candidates
//             const auto &recoParticles = pEvent->reco.particles;
//             const auto &truthParticles = pEvent->truth.particles;

//             const auto &sidebandRecoParticles = pEvent->reco.particles;
//             const auto &sidebandTruthParticles = pEvent->truth.particles;


//             // auto nRecoProtonsGeneric = 0u;
//             // for(const auto &pdg: {2212, 211, -211})
//             //     nRecoProtonsGeneric += std::count(assignedPdgCodes.begin(), assignedPdgCodes.end(), pdg);
//             // auto nRecoProtonsSideband = 0u;
//             // for(const auto &pdg: {2212, 211, -211})
//             //     nRecoProtonsSideband += std::count(sidebandAssignedPdgCodes.begin(), sidebandAssignedPdgCodes.end(), pdg);

//             auto nRecoProtonsGeneric = assignedPdgCodes.size()-1;
//             auto nRecoProtonsSideband = sidebandAssignedPdgCodes.size()-1;

//             if(sampleType == AnalysisHelper::Overlay)
//             {
//                 pTrueVsRecoProtonMultiplicityAllGenericGenie->Fill(std::min(5UL, nRecoProtonsGeneric), std::min(5u, truthData.nProtons), weight);
//                 pTrueVsRecoProtonMultiplicityAllSidebandGenie->Fill(std::min(5UL, nRecoProtonsSideband), std::min(5u, truthData.nProtons), weight);
//             }
//             else if (sampleType == AnalysisHelper::DataBNB)
//             {
//                 pTrueVsRecoProtonMultiplicityAllGenericNuWro->Fill(std::min(5UL, nRecoProtonsGeneric), std::min(5u, truthData.nProtons), weight);
//                 pTrueVsRecoProtonMultiplicityAllSidebandNuWro->Fill(std::min(5UL, nRecoProtonsSideband), std::min(5u, truthData.nProtons), weight);
//             }

//             // const auto muonPlotStyle = PlottingHelper::GetPlotStyle(muon, sampleType, truthParticles, false, config.global.useAbsPdg);
//             // const auto pionPlotStyle = PlottingHelper::GetPlotStyle(pion, sampleType, truthParticles, false, config.global.useAbsPdg);

//             // Fill the plots
//             // muonMomentumPlot.Fill(recoData.muonMomentum, plotStyle, weight);
//             // muonCosThetaPlot.Fill(recoData.muonCosTheta, plotStyle, weight);
//             // muonPhiPlot.Fill(recoData.muonPhi, plotStyle, weight);

//             // muonMomentumParticlePlot.Fill(recoData.muonMomentum, muonPlotStyle, weight);
//             // muonCosThetaParticlePlot.Fill(recoData.muonCosTheta, muonPlotStyle, weight);
//             // muonPhiParticlePlot.Fill(recoData.muonPhi, muonPlotStyle, weight);

//             // if (recoData.hasGoldenPion)
//             // {
//             //     pionMomentumPlot.Fill(recoData.pionMomentum, plotStyle, weight);
//             //     pionMomentumParticlePlot.Fill(recoData.pionMomentum, pionPlotStyle, weight);
//             // }

//             // pionCosThetaPlot.Fill(recoData.pionCosTheta, plotStyle, weight);
//             // pionPhiPlot.Fill(recoData.pionPhi, plotStyle, weight);

//             // pionCosThetaParticlePlot.Fill(recoData.pionCosTheta, pionPlotStyle, weight);
//             // pionPhiParticlePlot.Fill(recoData.pionPhi, pionPlotStyle, weight);

//             // muonPionAnglePlot.Fill(recoData.muonPionAngle, plotStyle, weight);
//             // nProtonsTruePlot.Fill(recoData.nProtons, plotStyle, weight);
//             if(sampleType == AnalysisHelper::Overlay)
//             {
//                 pNProtonsTruePlotGenieAll->Fill(              std::min(6u, truthData.nProtons),                 weight);
//                 pNProtonsRecoPlotGenieGeneric->Fill(          std::min(6u, recoData.nProtons),                  weight);
//                 pNProtonsRecoPlotGenieSideband->Fill(         std::min(6u, sidebandRecoData.nProtons),          weight);

//                 pProtonMomentumTruePlotGenieAll->Fill(        truthData.protonMomentum,                         weight);
//                 pProtonMomentumRecoPlotGenieGeneric->Fill(    recoData.protonMomentum,                          weight);
//                 pProtonMomentumRecoPlotGenieSideband->Fill(   sidebandRecoData.protonMomentum,                  weight);
//                 pProtonCosThetaTruePlotGenieAll->Fill(        truthData.protonCosTheta,                         weight);
//                 pProtonCosThetaRecoPlotGenieGeneric->Fill(    recoData.protonCosTheta,                          weight);
//                 pProtonCosThetaRecoPlotGenieSideband->Fill(   sidebandRecoData.protonCosTheta,                  weight);
//                 pProtonPhiTruePlotGenieAll->Fill(             truthData.protonPhi,                              weight);
//                 pProtonPhiRecoPlotGenieGeneric->Fill(         recoData.protonPhi,                               weight);
//                 pProtonPhiRecoPlotGenieSideband->Fill(        sidebandRecoData.protonPhi,                       weight);

//                 pMuonMomentumTruePlotGenieAll->Fill(          truthData.muonMomentum,                           weight);
//                 pMuonMomentumRecoPlotGenieGeneric->Fill(      recoData.muonMomentum,                            weight);
//                 pMuonMomentumRecoPlotGenieSideband->Fill(     sidebandRecoData.muonMomentum,                    weight);
//                 pMuonCosThetaTruePlotGenieAll->Fill(          truthData.muonCosTheta,                           weight);
//                 pMuonCosThetaRecoPlotGenieGeneric->Fill(      recoData.muonCosTheta,                            weight);
//                 pMuonCosThetaRecoPlotGenieSideband->Fill(     sidebandRecoData.muonCosTheta,                    weight);
//                 pMuonPhiTruePlotGenieAll->Fill(               truthData.muonPhi,                                weight);
//                 pMuonPhiRecoPlotGenieGeneric->Fill(           recoData.muonPhi,                                 weight);
//                 pMuonPhiRecoPlotGenieSideband->Fill(          sidebandRecoData.muonPhi,                         weight);

//             }
//             else if (sampleType == AnalysisHelper::DataBNB)
//             {
//                 pNProtonsTruePlotNuWroAll->Fill(              truthData.nProtons,                               weight);
//                 pNProtonsRecoPlotNuWroGeneric->Fill(          recoData.nProtons,                                weight);
//                 pNProtonsRecoPlotNuWroSideband->Fill(         sidebandRecoData.nProtons,                        weight);

//                 pProtonMomentumTruePlotNuWroAll->Fill(        truthData.protonMomentum,                         weight);
//                 pProtonMomentumRecoPlotNuWroGeneric->Fill(    recoData.protonMomentum,                          weight);
//                 pProtonMomentumRecoPlotNuWroSideband->Fill(   sidebandRecoData.protonMomentum,                  weight);
//                 pProtonCosThetaTruePlotNuWroAll->Fill(        truthData.protonCosTheta,                         weight);
//                 pProtonCosThetaRecoPlotNuWroGeneric->Fill(    recoData.protonCosTheta,                          weight);
//                 pProtonCosThetaRecoPlotNuWroSideband->Fill(   sidebandRecoData.protonCosTheta,                  weight);
//                 pProtonPhiTruePlotNuWroAll->Fill(             truthData.protonPhi,                              weight);
//                 pProtonPhiRecoPlotNuWroGeneric->Fill(         recoData.protonPhi,                               weight);
//                 pProtonPhiRecoPlotNuWroSideband->Fill(        sidebandRecoData.protonPhi,                       weight);

//                 pMuonMomentumTruePlotNuWroAll->Fill(          truthData.muonMomentum,                           weight);
//                 pMuonMomentumRecoPlotNuWroGeneric->Fill(      recoData.muonMomentum,                            weight);
//                 pMuonMomentumRecoPlotNuWroSideband->Fill(     sidebandRecoData.muonMomentum,                    weight);
//                 pMuonCosThetaTruePlotNuWroAll->Fill(          truthData.muonCosTheta,                           weight);
//                 pMuonCosThetaRecoPlotNuWroGeneric->Fill(      recoData.muonCosTheta,                            weight);
//                 pMuonCosThetaRecoPlotNuWroSideband->Fill(     sidebandRecoData.muonCosTheta,                    weight);
//                 pMuonPhiTruePlotNuWroAll->Fill(               truthData.muonPhi,                                weight);
//                 pMuonPhiRecoPlotNuWroGeneric->Fill(           recoData.muonPhi,                                 weight);
//                 pMuonPhiRecoPlotNuWroSideband->Fill(          sidebandRecoData.muonPhi,                         weight);
//             }


//             // Not all selected events have protons. Use only those that do to plot reconstructed proton variables.

//             // Only get the reconstructed proton variables when a suitable candidate is present
//             // std::cout<<"PlotProtonVariables Point 2"<<std::endl;

//             if (isSelectedGeneric)
//             {
//                 if(sampleType == AnalysisHelper::Overlay) pTrueVsRecoProtonMultiplicityGenericGenie->Fill(std::min(5UL,recoParticles.size()-2),std::min(5u,truthData.nProtons),weight);
//                 else if (sampleType == AnalysisHelper::DataBNB) pTrueVsRecoProtonMultiplicityGenericNuWro->Fill(std::min(5UL,recoParticles.size()-2),std::min(5u,truthData.nProtons),weight);

//                 if(sampleType == AnalysisHelper::Overlay)
//                 {
//                     pNProtonsTruePlotGenieGeneric->Fill(        std::min(6u, truthData.nProtons),       weight);

//                     pProtonMomentumTruePlotGenieGeneric->Fill(  truthData.protonMomentum, weight);
//                     pProtonCosThetaTruePlotGenieGeneric->Fill(  truthData.protonCosTheta, weight);
//                     pProtonPhiTruePlotGenieGeneric->Fill(       truthData.protonPhi,      weight);
//                     pMuonMomentumTruePlotGenieGeneric->Fill(    truthData.muonMomentum,   weight);
//                     pMuonCosThetaTruePlotGenieGeneric->Fill(    truthData.muonCosTheta,   weight);
//                     pMuonPhiTruePlotGenieGeneric->Fill(         truthData.muonPhi,        weight);
//                 }
//                 else if (sampleType == AnalysisHelper::DataBNB)
//                 {
//                     pNProtonsTruePlotNuWroGeneric->Fill(std::min(6u, truthData.nProtons), weight);

//                     pProtonMomentumTruePlotNuWroGeneric->Fill(  truthData.protonMomentum, weight);
//                     pProtonCosThetaTruePlotNuWroGeneric->Fill(  truthData.protonCosTheta, weight);
//                     pProtonPhiTruePlotNuWroGeneric->Fill(       truthData.protonPhi,      weight);
//                     pMuonMomentumTruePlotNuWroGeneric->Fill(    truthData.muonMomentum,   weight);
//                     pMuonCosThetaTruePlotNuWroGeneric->Fill(    truthData.muonCosTheta,   weight);
//                     pMuonPhiTruePlotNuWroGeneric->Fill(         truthData.muonPhi,        weight);
//                 }

//                 // std::cout<<"PlotProtonVariables Point 2.1"<<std::endl;
//                 const auto protonIndex = SelectionHelper::GetLeadingProtonCandidateIndex(recoParticles, assignedPdgCodes);
//                 // std::cout<<"PlotProtonVariables Point 2.2"<<std::endl;
//                 if (protonIndex!=std::numeric_limits<unsigned int>::max())
//                 {
//                     // std::cout<<"PlotProtonVariables Point 2.3"<<std::endl;
//                     const auto &muon = recoParticles.at(AnalysisHelper::GetParticleIndexWithPdg(assignedPdgCodes, 13));
//                     const auto &pion = recoParticles.at(AnalysisHelper::GetParticleIndexWithPdg(assignedPdgCodes, 211));
//                     const auto &proton = recoParticles.at(protonIndex);
//                     const auto muonDir = TVector3(muon.directionX(), muon.directionY(), muon.directionZ()).Unit();
//                     const auto pionDir = TVector3(pion.directionX(), pion.directionY(), pion.directionZ()).Unit();
//                     const auto protonDir = TVector3(proton.directionX(), proton.directionY(), proton.directionZ()).Unit();
//                     const auto protonMomentum = AnalysisHelper::GetProtonMomentumFromRange(proton.range());
//                     const auto protonCosTheta = protonDir.Z();
//                     const auto protonPhi = std::atan2(protonDir.Y(), protonDir.X());
//                     const auto protonPionAngle = std::acos(protonDir.Dot(pionDir));
//                     const auto protonMuonAngle = std::acos(protonDir.Dot(muonDir));

//                     // std::cout<<"PlotProtonVariables Point 2.4"<<std::endl;
//                     const auto protonPlotStyle = PlottingHelper::GetPlotStyle(proton, sampleType, truthParticles, false, config.global.useAbsPdg);
//                     for (unsigned int plot = 0;plot < protonPlotNames.size();++plot)
//                     {
//                         // Check proton multiplicity to match conditions described in protonPlotNames
//                         // No need to to test for plot==0 due to: if (protonIndex!=std::numeric_limits<unsigned int>::max())
//                         // std::cout<<"PlotProtonVariables Point 3"<<std::endl;
//                         if(plot==1 && recoData.nProtons!=1) continue;
//                         if(plot==2 && recoData.nProtons<2) continue;
//                         // std::cout<<"PlotProtonVariables Point 4"<<std::endl;

//                         protonMomentumPlots.at(plot).Fill(protonMomentum, plotStyle, weight);
//                         protonCosThetaPlots.at(plot).Fill(protonCosTheta, plotStyle, weight);
//                         protonPhiPlots.at(plot).Fill(protonPhi, plotStyle, weight);
//                         // std::cout<<"PlotProtonVariables Point 4.2"<<std::endl;
//                         protonMomentumParticlePlots.at(plot).Fill(protonMomentum, protonPlotStyle, weight);
//                         protonCosThetaParticlePlots.at(plot).Fill(protonCosTheta, protonPlotStyle, weight);
//                         protonPhiParticlePlots.at(plot).Fill(protonPhi, protonPlotStyle, weight);

//                         protonPionAnglePlots.at(plot).Fill(protonPionAngle, plotStyle, weight);
//                         protonMuonAnglePlots.at(plot).Fill(protonMuonAngle, plotStyle, weight);
//                         // std::cout<<"PlotProtonVariables Point 4.31"<<std::endl;
//                     }
//                 }
//             }


//             auto protonIndexTruth = AnalysisHelper::GetTrueLeadingProtonIndex(pEvent->truth, true, config.global.protonMomentumThreshold);
//             if (protonIndexTruth!=std::numeric_limits<unsigned int>::max())
//             {
//                 // std::cout<<"PlotProtonVariables Point 4.32"<<std::endl;

//                 const auto muonIndex = AnalysisHelper::GetTrueMuonIndex(pEvent->truth, true);
//                 const auto &muonTruth = truthParticles.at(muonIndex);
//                 // std::cout<<"PlotProtonVariables Point 4.32..1"<<std::endl;
//                 const auto &protonTruth = truthParticles.at(protonIndexTruth);
//                 // std::cout<<"PlotProtonVariables Point 4.32.2"<<std::endl;
//                 const auto muonDirTruth = TVector3(muonTruth.momentumX(), muonTruth.momentumY(), muonTruth.momentumZ()).Unit();
//                 const auto protonDirTruth = TVector3(protonTruth.momentumX(), protonTruth.momentumY(), protonTruth.momentumZ()).Unit();
//                 // std::cout<<"PlotProtonVariables Point 4.33"<<std::endl;
//                 const auto protonMomentumTruth = protonTruth.momentum();
//                 const auto protonCosThetaTruth = protonDirTruth.Z();
//                 // std::cout<<"PlotProtonVariables Point 4.34"<<std::endl;
//                 const auto protonPhiTruth = std::atan2(protonDirTruth.Y(), protonDirTruth.X());
//                 // const auto protonPionAngleTruth = std::acos(protonDirTruth.Dot(pionDirTruth));
//                 const auto protonMuonAngleTruth = std::acos(protonDirTruth.Dot(muonDirTruth));

//                 if(sampleType == AnalysisHelper::Overlay) pTrueProtonMomvsTrueProtonMultiplicityAllGenie->Fill(std::min(9u,truthData.nProtons), protonMomentumTruth,weight);
//                 else if (sampleType == AnalysisHelper::DataBNB) pTrueProtonMomvsTrueProtonMultiplicityAllNuWro->Fill(std::min(9u,truthData.nProtons), protonMomentumTruth,weight);
//                 if (isSelectedGeneric)
//                 {
//                     if(sampleType == AnalysisHelper::Overlay) pTrueProtonMomvsTrueProtonMultiplicityGenericGenie->Fill(std::min(9u,truthData.nProtons), protonMomentumTruth,weight);
//                     else if (sampleType == AnalysisHelper::DataBNB) pTrueProtonMomvsTrueProtonMultiplicityGenericNuWro->Fill(std::min(9u,truthData.nProtons), protonMomentumTruth,weight);
//                 }
//                 else
//                 {
//                     if(sampleType == AnalysisHelper::Overlay) pTrueProtonMomvsTrueProtonMultiplicityNotGenericGenie->Fill(std::min(9u,truthData.nProtons), protonMomentumTruth,weight);
//                     else if (sampleType == AnalysisHelper::DataBNB) pTrueProtonMomvsTrueProtonMultiplicityNotGenericNuWro->Fill(std::min(9u,truthData.nProtons), protonMomentumTruth,weight);
//                 }

//                 if (passedSidebandSelection)
//                 {
//                     if(sampleType == AnalysisHelper::Overlay) pTrueProtonMomvsTrueProtonMultiplicitySidebandGenie->Fill(std::min(9u,truthData.nProtons), protonMomentumTruth,weight);
//                     else if (sampleType == AnalysisHelper::DataBNB) pTrueProtonMomvsTrueProtonMultiplicitySidebandNuWro->Fill(std::min(9u,truthData.nProtons), protonMomentumTruth,weight);
//                 }
//                 else
//                 {
//                     if(sampleType == AnalysisHelper::Overlay) pTrueProtonMomvsTrueProtonMultiplicityNotSidebandGenie->Fill(std::min(9u,truthData.nProtons), protonMomentumTruth,weight);
//                     else if (sampleType == AnalysisHelper::DataBNB) pTrueProtonMomvsTrueProtonMultiplicityNotSidebandNuWro->Fill(std::min(9u,truthData.nProtons), protonMomentumTruth,weight);
//                 }

//                 if(sampleType == AnalysisHelper::Overlay)
//                 {
//                     pProtonMomvsPhiAllGenie->Fill(protonPhiTruth, protonMomentumTruth, weight);
//                     pProtonMomvsCosThetaAllGenie->Fill(protonCosThetaTruth, protonMomentumTruth, weight);
//                     pProtonPhivsCosThetaAllGenie->Fill(protonCosThetaTruth, protonPhiTruth, weight);
//                     if (isSelectedGeneric)
//                     {
//                         pProtonMomvsPhiGenericGenie->Fill(protonPhiTruth, protonMomentumTruth, weight);
//                         pProtonMomvsCosThetaGenericGenie->Fill(protonCosThetaTruth, protonMomentumTruth, weight);
//                         pProtonPhivsCosThetaGenericGenie->Fill(protonCosThetaTruth, protonPhiTruth, weight);
//                     }
//                     if (passedSidebandSelection)
//                     {
//                         pProtonMomvsPhiSidebandGenie->Fill(protonPhiTruth, protonMomentumTruth, weight);
//                         pProtonMomvsCosThetaSidebandGenie->Fill(protonCosThetaTruth, protonMomentumTruth, weight);
//                         pProtonPhivsCosThetaSidebandGenie->Fill(protonCosThetaTruth, protonPhiTruth, weight);
//                     }
//                 }
//                 else if (sampleType == AnalysisHelper::DataBNB)
//                 {
//                     pProtonMomvsPhiAllNuWro->Fill(protonPhiTruth, protonMomentumTruth, weight);
//                     pProtonMomvsCosThetaAllNuWro->Fill(protonCosThetaTruth, protonMomentumTruth, weight);
//                     pProtonPhivsCosThetaAllNuWro->Fill(protonCosThetaTruth, protonPhiTruth, weight);
//                     if (isSelectedGeneric) {
//                         pProtonMomvsPhiGenericNuWro->Fill(protonPhiTruth, protonMomentumTruth, weight);
//                         pProtonMomvsCosThetaGenericNuWro->Fill(protonCosThetaTruth, protonMomentumTruth, weight);
//                         pProtonPhivsCosThetaGenericNuWro->Fill(protonCosThetaTruth, protonPhiTruth, weight);
//                     }
//                     if (passedSidebandSelection)
//                     {
//                         pProtonMomvsPhiSidebandNuWro->Fill(protonPhiTruth, protonMomentumTruth, weight);
//                         pProtonMomvsCosThetaSidebandNuWro->Fill(protonCosThetaTruth, protonMomentumTruth, weight);
//                         pProtonPhivsCosThetaSidebandNuWro->Fill(protonCosThetaTruth, protonPhiTruth, weight);
//                     }
//                 }


//                 // const auto protonPlotStyle = PlottingHelper::GetPlotStyle(protonTruth, sampleType, truthParticles, false, config.global.useAbsPdg);
//                 for (unsigned int plot = 0;plot < protonPlotNames.size();++plot)
//                 {
//                     // Check proton multiplicity to match conditions described in protonPlotNames
//                     // No need to to test for plot==0 due to: if (protonIndex!=std::numeric_limits<unsigned int>::max())
//                     if(plot==1 && truthData.nProtons!=1) continue;
//                     if(plot==2 && truthData.nProtons<2) continue;

//                     if (isSelectedGeneric)
//                     {
//                         protonMomentumPlotsTruth.at(plot).Fill(protonMomentumTruth, plotStyle, weight);
//                         protonCosThetaPlotsTruth.at(plot).Fill(protonCosThetaTruth, plotStyle, weight);
//                         protonPhiPlotsTruth.at(plot).Fill(protonPhiTruth, plotStyle, weight);

//                         protonMomentumParticlePlotsTruth.at(plot).Fill(protonMomentumTruth, plotStyle, weight);
//                         protonCosThetaParticlePlotsTruth.at(plot).Fill(protonCosThetaTruth, plotStyle, weight);
//                         protonPhiParticlePlotsTruth.at(plot).Fill(protonPhiTruth, plotStyle, weight);

//                         // protonPionAnglePlotsTruth.at(plot).Fill(protonPionAngleTruth, plotStyle, weight);
//                         protonMuonAnglePlotsTruth.at(plot).Fill(protonMuonAngleTruth, plotStyle, weight);
//                     }

//                     if (passedSidebandSelection)
//                     {
//                         protonMomentumPlotsSidebandTruth.at(plot).Fill(protonMomentumTruth, plotStyle, weight);
//                         protonCosThetaPlotsSidebandTruth.at(plot).Fill(protonCosThetaTruth, plotStyle, weight);
//                         protonPhiPlotsSidebandTruth.at(plot).Fill(protonPhiTruth, plotStyle, weight);

//                         protonMomentumParticlePlotsSidebandTruth.at(plot).Fill(protonMomentumTruth, plotStyle, weight);
//                         protonCosThetaParticlePlotsSidebandTruth.at(plot).Fill(protonCosThetaTruth, plotStyle, weight);
//                         protonPhiParticlePlotsSidebandTruth.at(plot).Fill(protonPhiTruth, plotStyle, weight);

//                         // protonPionAnglePlotsTruth.at(plot).Fill(protonPionAngleTruth, plotStyle, weight);
//                         protonMuonAnglePlotsSidebandTruth.at(plot).Fill(protonMuonAngleTruth, plotStyle, weight);
//                     }
//                 }
//             }

//             // std::cout<<"PlotProtonVariables Point 5"<<std::endl;
//             if (passedSidebandSelection)
//             {
//                 if(sampleType == AnalysisHelper::Overlay) pTrueVsRecoProtonMultiplicitySidebandGenie->Fill(std::min(5UL,recoParticles.size()-1),std::min(5u,truthData.nProtons),weight);
//                 else if (sampleType == AnalysisHelper::DataBNB) pTrueVsRecoProtonMultiplicitySidebandNuWro->Fill(std::min(5UL,recoParticles.size()-1),std::min(5u,truthData.nProtons),weight);

//                 if(sampleType == AnalysisHelper::Overlay)
//                 {
//                     pNProtonsTruePlotGenieSideband->Fill(        std::min(6u, truthData.nProtons),       weight);

//                     pProtonMomentumTruePlotGenieSideband->Fill(  truthData.protonMomentum, weight);
//                     pProtonCosThetaTruePlotGenieSideband->Fill(  truthData.protonCosTheta, weight);
//                     pProtonPhiTruePlotGenieSideband->Fill(       truthData.protonPhi,      weight);
//                     pMuonMomentumTruePlotGenieSideband->Fill(    truthData.muonMomentum,   weight);
//                     pMuonCosThetaTruePlotGenieSideband->Fill(    truthData.muonCosTheta,   weight);
//                     pMuonPhiTruePlotGenieSideband->Fill(         truthData.muonPhi,        weight);
//                 }
//                 else if (sampleType == AnalysisHelper::DataBNB)
//                 {
//                     pNProtonsTruePlotNuWroSideband->Fill(std::min(6u, truthData.nProtons), weight);

//                     pProtonMomentumTruePlotNuWroSideband->Fill(  truthData.protonMomentum, weight);
//                     pProtonCosThetaTruePlotNuWroSideband->Fill(  truthData.protonCosTheta, weight);
//                     pProtonPhiTruePlotNuWroSideband->Fill(       truthData.protonPhi,      weight);
//                     pMuonMomentumTruePlotNuWroSideband->Fill(    truthData.muonMomentum,   weight);
//                     pMuonCosThetaTruePlotNuWroSideband->Fill(    truthData.muonCosTheta,   weight);
//                     pMuonPhiTruePlotNuWroSideband->Fill(         truthData.muonPhi,        weight);
//                 }

//                 const auto protonIndex = SelectionHelper::GetLeadingProtonCandidateIndex(sidebandRecoParticles, sidebandAssignedPdgCodes);
//                 if (protonIndex!=std::numeric_limits<unsigned int>::max())
//                 {
//                     // std::cout<<"PlotProtonVariables Point 6"<<std::endl;
//                     const auto &muon = sidebandRecoParticles.at(AnalysisHelper::GetParticleIndexWithPdg(sidebandAssignedPdgCodes, 13));
//                     // const auto &pion = sidebandRecoParticles.at(AnalysisHelper::GetParticleIndexWithPdg(sidebandAssignedPdgCodes, 211));
//                     const auto &proton = sidebandRecoParticles.at(protonIndex);
//                     const auto muonDir = TVector3(muon.directionX(), muon.directionY(), muon.directionZ()).Unit();
//                     // const auto pionDir = TVector3(pion.directionX(), pion.directionY(), pion.directionZ()).Unit();
//                     const auto protonDir = TVector3(proton.directionX(), proton.directionY(), proton.directionZ()).Unit();

//                     const auto protonMomentum = AnalysisHelper::GetProtonMomentumFromRange(proton.range());
//                     const auto protonCosTheta = protonDir.Z();
//                     const auto protonPhi  = std::atan2(protonDir.Y(), protonDir.X());
//                     // const auto protonPionAngle = std::acos(protonDir.Dot(pionDir));
//                     const auto protonMuonAngle = std::acos(protonDir.Dot(muonDir));

//                     const auto protonPlotStyle = PlottingHelper::GetPlotStyle(proton, sampleType, truthParticles, false, config.global.useAbsPdg);
//                     for (unsigned int plot = 0;plot < protonPlotNames.size();++plot)
//                     {
//                         // std::cout<<"PlotProtonVariables Point 7"<<std::endl;
//                         // Check proton multiplicity to match conditions described in protonPlotNames
//                         // No need to to test for plot==0 due to: if (protonIndex!=std::numeric_limits<unsigned int>::max())
//                         if(plot==1 && sidebandRecoData.nProtons!=1) continue;
//                         if(plot==2 && sidebandRecoData.nProtons<2) continue;
//                         // std::cout<<"PlotProtonVariables Point 8"<<std::endl;

//                         protonMomentumPlotsSideband.at(plot).Fill(protonMomentum, plotStyle, weight);
//                         protonCosThetaPlotsSideband.at(plot).Fill(protonCosTheta, plotStyle, weight);
//                         protonPhiPlotsSideband.at(plot).Fill(protonPhi, plotStyle, weight);

//                         protonMomentumParticlePlotsSideband.at(plot).Fill(protonMomentum, protonPlotStyle, weight);
//                         protonCosThetaParticlePlotsSideband.at(plot).Fill(protonCosTheta, protonPlotStyle, weight);
//                         protonPhiParticlePlotsSideband.at(plot).Fill(protonPhi, protonPlotStyle, weight);

//                         // protonPionAnglePlots.at(plot).Fill(protonPionAngle, plotStyle, weight);
//                         protonMuonAnglePlotsSideband.at(plot).Fill(protonMuonAngle, plotStyle, weight);
//                     }
//                 }
//             }
//         }
//     }

//     // muonMomentumPlot.SaveAsStacked("reco_muonMomentum",false,false,config.global.axisTitles);
//     // muonCosThetaPlot.SaveAsStacked("reco_muonCosTheta",false,false,config.global.axisTitles);
//     // muonPhiPlot.SaveAsStacked("reco_muonPhi",false,false,config.global.axisTitles);

//     // muonMomentumParticlePlot.SaveAsStacked("reco_muonMomentum_particle",false,false,config.global.axisTitles);
//     // muonCosThetaParticlePlot.SaveAsStacked("reco_muonCosTheta_particle",false,false,config.global.axisTitles);
//     // muonPhiParticlePlot.SaveAsStacked("reco_muonPhi_particle",false,false,config.global.axisTitles);

//     // pionMomentumPlot.SaveAsStacked("reco_pionMomentum",false,false,config.global.axisTitles);
//     // pionCosThetaPlot.SaveAsStacked("reco_pionCosTheta",false,false,config.global.axisTitles);
//     // pionPhiPlot.SaveAsStacked("reco_pionPhi",false,false,config.global.axisTitles);

//     // pionMomentumParticlePlot.SaveAsStacked("reco_pionMomentum_particle",false,false,config.global.axisTitles);
//     // pionCosThetaParticlePlot.SaveAsStacked("reco_pionCosTheta_particle",false,false,config.global.axisTitles);
//     // pionPhiParticlePlot.SaveAsStacked("reco_pionPhi_particle",false,false,config.global.axisTitles);

//     // muonPionAnglePlot.SaveAsStacked("reco_muonPionAngle",false,false,config.global.axisTitles);
//     // nProtonsTruePlot.SaveAsStacked("true_nProtons",false,false,config.global.axisTitles);


//     const std::string prefix = "CC0pi";
//     for (unsigned int plot = 0;plot < protonPlotNames.size();++plot)
//     {
//         const auto plotName = protonPlotNames.at(plot);

//         protonMomentumPlots.at(plot).SaveAsStacked("reco_protonMomentum_" + prefix + "_" + plotName,false,false,false,config.global.axisTitles);
//         protonCosThetaPlots.at(plot).SaveAsStacked("reco_protonCosTheta_" + prefix + "_" + plotName,false,false,false,config.global.axisTitles);
//         protonPhiPlots.at(plot).SaveAsStacked("reco_protonPhi_" + prefix + "_" + plotName,false,false,false,config.global.axisTitles);

//         protonMomentumParticlePlots.at(plot).SaveAsStacked("reco_protonMomentum_particle_" + prefix + "_" + plotName,false,false,false,config.global.axisTitles);
//         protonCosThetaParticlePlots.at(plot).SaveAsStacked("reco_protonCosTheta_particle_" + prefix + "_" + plotName,false,false,false,config.global.axisTitles);
//         protonPhiParticlePlots.at(plot).SaveAsStacked("reco_protonPhi_particle_" + prefix + "_" + plotName,false,false,false,config.global.axisTitles);

//         protonPionAnglePlots.at(plot).SaveAsStacked("reco_protonPionAngle_" + prefix + "_" + plotName,false,false,false,config.global.axisTitles);
//         protonMuonAnglePlots.at(plot).SaveAsStacked("reco_protonMuonAngle_" + prefix + "_" + plotName,false,false,false,config.global.axisTitles);


//         protonMomentumPlotsTruth.at(plot).SaveAsStacked("truth_protonMomentum_" + prefix + "_" + plotName,false,false,false,config.global.axisTitles);
//         protonCosThetaPlotsTruth.at(plot).SaveAsStacked("truth_protonCosTheta_" + prefix + "_" + plotName,false,false,false,config.global.axisTitles);
//         protonPhiPlotsTruth.at(plot).SaveAsStacked("truth_protonPhi_" + prefix + "_" + plotName,false,false,false,config.global.axisTitles);

//         protonMomentumParticlePlotsTruth.at(plot).SaveAsStacked("truth_protonMomentum_particle_" + prefix + "_" + plotName,false,false,false,config.global.axisTitles);
//         protonCosThetaParticlePlotsTruth.at(plot).SaveAsStacked("truth_protonCosTheta_particle_" + prefix + "_" + plotName,false,false,false,config.global.axisTitles);
//         protonPhiParticlePlotsTruth.at(plot).SaveAsStacked("truth_protonPhi_particle_" + prefix + "_" + plotName,false,false,false,config.global.axisTitles);

//         // protonPionAnglePlotsTruth.at(plot).SaveAsStacked("truth_protonPionAngle_" + prefix + "_" + plotName,false,false,false,config.global.axisTitles);
//         protonMuonAnglePlotsTruth.at(plot).SaveAsStacked("truth_protonMuonAngle_" + prefix + "_" + plotName,false,false,false,config.global.axisTitles);


//         protonMomentumPlotsSideband.at(plot).SaveAsStacked("reco_sideband_protonMomentum_" + prefix + "_" + plotName,false,false,false,config.global.axisTitles);
//         protonCosThetaPlotsSideband.at(plot).SaveAsStacked("reco_sideband_protonCosTheta_" + prefix + "_" + plotName,false,false,false,config.global.axisTitles);
//         protonPhiPlotsSideband.at(plot).SaveAsStacked("reco_sideband_protonPhi_" + prefix + "_" + plotName,false,false,false,config.global.axisTitles);

//         protonMomentumParticlePlotsSideband.at(plot).SaveAsStacked("reco_sideband_protonMomentum_particle_" + prefix + "_" + plotName,false,false,false,config.global.axisTitles);
//         protonCosThetaParticlePlotsSideband.at(plot).SaveAsStacked("reco_sideband_protonCosTheta_particle_" + prefix + "_" + plotName,false,false,false,config.global.axisTitles);
//         protonPhiParticlePlotsSideband.at(plot).SaveAsStacked("reco_sideband_protonPhi_particle_" + prefix + "_" + plotName,false,false,false,config.global.axisTitles);

//         // protonPionAnglePlotsSideband.at(plot).SaveAsStacked("reco_sideband_protonPionAngle_" + prefix + "_" + plotName,false,false,false,config.global.axisTitles);
//         protonMuonAnglePlotsSideband.at(plot).SaveAsStacked("reco_sideband_protonMuonAngle_" + prefix + "_" + plotName,false,false,false,config.global.axisTitles);


//         protonMomentumPlotsSidebandTruth.at(plot).SaveAsStacked("truth_sideband_protonMomentum_" + prefix + "_" + plotName,false,false,false,config.global.axisTitles);
//         protonCosThetaPlotsSidebandTruth.at(plot).SaveAsStacked("truth_sideband_protonCosTheta_" + prefix + "_" + plotName,false,false,false,config.global.axisTitles);
//         protonPhiPlotsSidebandTruth.at(plot).SaveAsStacked("truth_sideband_protonPhi_" + prefix + "_" + plotName,false,false,false,config.global.axisTitles);

//         protonMomentumParticlePlotsSidebandTruth.at(plot).SaveAsStacked("truth_sideband_protonMomentum_particle_" + prefix + "_" + plotName,false,false,false,config.global.axisTitles);
//         protonCosThetaParticlePlotsSidebandTruth.at(plot).SaveAsStacked("truth_sideband_protonCosTheta_particle_" + prefix + "_" + plotName,false,false,false,config.global.axisTitles);
//         protonPhiParticlePlotsSidebandTruth.at(plot).SaveAsStacked("truth_sideband_protonPhi_particle_" + prefix + "_" + plotName,false,false,false,config.global.axisTitles);

//         // protonPionAnglePlotsTruth.at(plot).SaveAsStacked("truth_protonPionAngle_" + prefix + "_" + plotName,false,false,false,config.global.axisTitles);
//         protonMuonAnglePlotsSidebandTruth.at(plot).SaveAsStacked("truth_sideband_protonMuonAngle_" + prefix + "_" + plotName,false,false,false,config.global.axisTitles);
//     }

//         // Setup the output canvas for the remaining plots
//     auto pCanvas = PlottingHelper::GetCanvas();
//     pProtonMomvsPhiAllGenie->Draw("colz text");
//     PlottingHelper::SaveCanvas(pCanvas, "truth_all_genie_protonMomentumVsPhi");
//     pProtonMomvsPhiGenericGenie->Draw("colz text");
//     PlottingHelper::SaveCanvas(pCanvas, "truth_generic_genie_protonMomentumVsPhi");
//     pProtonMomvsPhiSidebandGenie->Draw("colz text");
//     PlottingHelper::SaveCanvas(pCanvas, "truth_sideband_genie_protonMomentumVsPhi");
//     pProtonMomvsPhiAllNuWro->Draw("colz text");
//     PlottingHelper::SaveCanvas(pCanvas, "truth_all_nuwro_protonMomentumVsPhi");
//     pProtonMomvsPhiGenericNuWro->Draw("colz text");
//     PlottingHelper::SaveCanvas(pCanvas, "truth_generic_nuwro_protonMomentumVsPhi");
//     pProtonMomvsPhiSidebandNuWro->Draw("colz text");
//     PlottingHelper::SaveCanvas(pCanvas, "truth_sideband_nuwro_protonMomentumVsPhi");

//     pProtonMomvsCosThetaAllGenie->Draw("colz text");
//     PlottingHelper::SaveCanvas(pCanvas, "truth_all_genie_protonMomentumVsCosTheta");
//     pProtonMomvsCosThetaGenericGenie->Draw("colz text");
//     PlottingHelper::SaveCanvas(pCanvas, "truth_generic_genie_protonMomentumVsCosTheta");
//     pProtonMomvsCosThetaSidebandGenie->Draw("colz text");
//     PlottingHelper::SaveCanvas(pCanvas, "truth_sideband_genie_protonMomentumVsCosTheta");
//     pProtonMomvsCosThetaAllNuWro->Draw("colz text");
//     PlottingHelper::SaveCanvas(pCanvas, "truth_all_nuwro_protonMomentumVsCosTheta");
//     pProtonMomvsCosThetaGenericNuWro->Draw("colz text");
//     PlottingHelper::SaveCanvas(pCanvas, "truth_generic_nuwro_protonMomentumVsCosTheta");
//     pProtonMomvsCosThetaSidebandNuWro->Draw("colz text");
//     PlottingHelper::SaveCanvas(pCanvas, "truth_sideband_nuwro_protonMomentumVsCosTheta");

//     pProtonPhivsCosThetaAllGenie->Draw("colz text");
//     PlottingHelper::SaveCanvas(pCanvas, "truth_all_genie_protonPhiVsCosTheta");
//     pProtonPhivsCosThetaGenericGenie->Draw("colz text");
//     PlottingHelper::SaveCanvas(pCanvas, "truth_generic_genie_protonPhiVsCosTheta");
//     pProtonPhivsCosThetaSidebandGenie->Draw("colz text");
//     PlottingHelper::SaveCanvas(pCanvas, "truth_sideband_genie_protonPhiVsCosTheta");
//     pProtonPhivsCosThetaAllNuWro->Draw("colz text");
//     PlottingHelper::SaveCanvas(pCanvas, "truth_all_nuwro_protonPhiVsCosTheta");
//     pProtonPhivsCosThetaGenericNuWro->Draw("colz text");
//     PlottingHelper::SaveCanvas(pCanvas, "truth_generic_nuwro_protonPhiVsCosTheta");
//     pProtonPhivsCosThetaSidebandNuWro->Draw("colz text");
//     PlottingHelper::SaveCanvas(pCanvas, "truth_sideband_nuwro_protonPhiVsCosTheta");

//     pTrueVsRecoProtonMultiplicityAllGenericGenie->Draw("colz text");
//     PlottingHelper::SaveCanvas(pCanvas, "all_generic_genie_trueVsRecoProtonMultiplicity");
//     pTrueVsRecoProtonMultiplicityAllGenericNuWro->Draw("colz text");
//     PlottingHelper::SaveCanvas(pCanvas, "all_generic_nuwro_trueVsRecoProtonMultiplicity");
//     pTrueVsRecoProtonMultiplicityAllSidebandGenie->Draw("colz text");
//     PlottingHelper::SaveCanvas(pCanvas, "all_sideband_genie_trueVsRecoProtonMultiplicity");
//     pTrueVsRecoProtonMultiplicityAllSidebandNuWro->Draw("colz text");
//     PlottingHelper::SaveCanvas(pCanvas, "all_sideband_nuwro_trueVsRecoProtonMultiplicity");

//     pTrueVsRecoProtonMultiplicityGenericNuWro->Draw("colz text");
//     PlottingHelper::SaveCanvas(pCanvas, "generic_nuwro_trueVsRecoProtonMultiplicity");
//     pTrueVsRecoProtonMultiplicitySidebandNuWro->Draw("colz text");
//     PlottingHelper::SaveCanvas(pCanvas, "sideband_nuwro_trueVsRecoProtonMultiplicity");
//     pTrueVsRecoProtonMultiplicityGenericGenie->Draw("colz text");
//     PlottingHelper::SaveCanvas(pCanvas, "generic_genie_trueVsRecoProtonMultiplicity");
//     pTrueVsRecoProtonMultiplicitySidebandGenie->Draw("colz text");
//     PlottingHelper::SaveCanvas(pCanvas, "sideband_genie_trueVsRecoProtonMultiplicity");

//     pTrueProtonMomvsTrueProtonMultiplicityAllNuWro->Draw("colz text");
//     PlottingHelper::SaveCanvas(pCanvas, "all_nuwro_trueProtonMomVsTrueProtonMultiplicity");
//     pTrueProtonMomvsTrueProtonMultiplicityGenericNuWro->Draw("colz text");
//     PlottingHelper::SaveCanvas(pCanvas, "generic_nuwro_trueProtonMomVsTrueProtonMultiplicity");
//     pTrueProtonMomvsTrueProtonMultiplicitySidebandNuWro->Draw("colz text");
//     PlottingHelper::SaveCanvas(pCanvas, "sideband_nuwro_trueProtonMomVsTrueProtonMultiplicity");
//     pTrueProtonMomvsTrueProtonMultiplicityNotGenericNuWro->Draw("colz text");
//     PlottingHelper::SaveCanvas(pCanvas, "notGeneric_nuwro_trueProtonMomVsTrueProtonMultiplicity");
//     pTrueProtonMomvsTrueProtonMultiplicityNotSidebandNuWro->Draw("colz text");
//     PlottingHelper::SaveCanvas(pCanvas, "notSideband_nuwro_trueProtonMomVsTrueProtonMultiplicity");

//     pTrueProtonMomvsTrueProtonMultiplicityAllGenie->Draw("colz text");
//     PlottingHelper::SaveCanvas(pCanvas, "all_genie_trueProtonMomVsTrueProtonMultiplicity");
//     pTrueProtonMomvsTrueProtonMultiplicityGenericGenie->Draw("colz text");
//     PlottingHelper::SaveCanvas(pCanvas, "generic_genie_trueProtonMomVsTrueProtonMultiplicity");
//     pTrueProtonMomvsTrueProtonMultiplicitySidebandGenie->Draw("colz text");
//     PlottingHelper::SaveCanvas(pCanvas, "sideband_genie_trueProtonMomVsTrueProtonMultiplicity");
//     pTrueProtonMomvsTrueProtonMultiplicityNotGenericGenie->Draw("colz text");
//     PlottingHelper::SaveCanvas(pCanvas, "notGeneric_genie_trueProtonMomVsTrueProtonMultiplicity");
//     pTrueProtonMomvsTrueProtonMultiplicityNotSidebandGenie->Draw("colz text");
//     PlottingHelper::SaveCanvas(pCanvas, "notSideband_genie_trueProtonMomVsTrueProtonMultiplicity");


//     *pTrueProtonMomvsTrueProtonMultiplicityGenericNuWroRatio = (*pTrueProtonMomvsTrueProtonMultiplicityGenericNuWro)/(*pTrueProtonMomvsTrueProtonMultiplicityAllNuWro);
//     *pTrueProtonMomvsTrueProtonMultiplicitySidebandNuWroRatio = (*pTrueProtonMomvsTrueProtonMultiplicitySidebandNuWro)/(*pTrueProtonMomvsTrueProtonMultiplicityAllNuWro);
//     *pTrueProtonMomvsTrueProtonMultiplicityGenericGenieRatio = (*pTrueProtonMomvsTrueProtonMultiplicityGenericGenie)/(*pTrueProtonMomvsTrueProtonMultiplicityAllGenie);
//     *pTrueProtonMomvsTrueProtonMultiplicitySidebandGenieRatio = (*pTrueProtonMomvsTrueProtonMultiplicitySidebandGenie)/(*pTrueProtonMomvsTrueProtonMultiplicityAllGenie);

//     const std::vector<Double_t> cont = {0.0, 0.001, 0.002, 0.004, 0.006, 0.009, 0.012, 0.015, 0.019, 0.024, 0.029, 0.034, 0.04, 0.047, 0.054, 0.061, 0.069, 0.077, 0.086, 0.095, 0.105, 0.115, 0.126, 0.137, 0.149, 0.161, 0.173, 0.186, 0.2, 1};
//     pTrueProtonMomvsTrueProtonMultiplicityGenericNuWroRatio->SetContour(cont.size(), cont.data());
//     pTrueProtonMomvsTrueProtonMultiplicityGenericNuWroRatio->Draw("colz text");
//     PlottingHelper::SaveCanvas(pCanvas, "generic_nuwro_trueProtonMomVsTrueProtonMultiplicityRatio");
//     pTrueProtonMomvsTrueProtonMultiplicitySidebandNuWroRatio->SetContour(cont.size(), cont.data());
//     pTrueProtonMomvsTrueProtonMultiplicitySidebandNuWroRatio->Draw("colz text");
//     PlottingHelper::SaveCanvas(pCanvas, "sideband_nuwro_trueProtonMomVsTrueProtonMultiplicityRatio");
//     pTrueProtonMomvsTrueProtonMultiplicityGenericGenieRatio->SetContour(cont.size(), cont.data());
//     pTrueProtonMomvsTrueProtonMultiplicityGenericGenieRatio->Draw("colz text");
//     PlottingHelper::SaveCanvas(pCanvas, "generic_genie_trueProtonMomVsTrueProtonMultiplicityRatio");
//     pTrueProtonMomvsTrueProtonMultiplicitySidebandGenieRatio->SetContour(cont.size(), cont.data());
//     pTrueProtonMomvsTrueProtonMultiplicitySidebandGenieRatio->Draw("colz text");
//     PlottingHelper::SaveCanvas(pCanvas, "sideband_genie_trueProtonMomVsTrueProtonMultiplicityRatio");

//     for (auto &pHist : {pNProtonsTruePlotGenieAll, pNProtonsTruePlotNuWroAll, pNProtonsRecoPlotGenieGeneric, pNProtonsRecoPlotNuWroGeneric, pNProtonsRecoPlotGenieSideband, pNProtonsRecoPlotNuWroSideband,
//     pNProtonsTruePlotNuWroGeneric, pNProtonsTruePlotGenieGeneric, pNProtonsTruePlotNuWroSideband, pNProtonsTruePlotGenieSideband})
//     {
//         pHist->GetXaxis()->SetBinLabel(1, "0");
//         pHist->GetXaxis()->SetBinLabel(2, "1");
//         pHist->GetXaxis()->SetBinLabel(3, "2");
//         pHist->GetXaxis()->SetBinLabel(4, "3");
//         pHist->GetXaxis()->SetBinLabel(5, "4");
//         pHist->GetXaxis()->SetBinLabel(6, ">4");
//     }
//     // Draw the prediction uncertainties as a semi-transparent band
//     // auto pNProtonsTruePlotGenieClone = static_cast<TH1F *>(pNProtonsTruePlotGenie->Clone());
//     // const auto col = pNProtonsTruePlotGenieClone->GetLineColor();
//     // pNProtonsTruePlotGenieClone->SetFillStyle(1001);
//     // pNProtonsTruePlotGenieClone->SetLineColorAlpha(col, 0.f);
//     // pNProtonsTruePlotGenieClone->SetFillColorAlpha(col, 0.3f);


//     //################################# nProtons
//     PlottingHelper::SetLineStyle(pNProtonsTruePlotNuWroAll, PlottingHelper::Primary);
//     PlottingHelper::SetLineStyle(pNProtonsTruePlotGenieAll, PlottingHelper::Secondary);
//     auto maxY = 1.05 * std::max(pNProtonsTruePlotGenieAll->GetMaximum(), pNProtonsTruePlotGenieAll->GetMaximum());
//     pNProtonsTruePlotGenieAll->GetYaxis()->SetRangeUser(0, maxY);
//     pNProtonsTruePlotNuWroAll->GetYaxis()->SetRangeUser(0, maxY);
//     pNProtonsTruePlotGenieAll->Draw("hist");
//     pNProtonsTruePlotNuWroAll->Draw("hist same");
//     PlottingHelper::SaveCanvas(pCanvas, prefix+"_nProtonsTrueAll");

//     PlottingHelper::SetLineStyle(pNProtonsRecoPlotNuWroGeneric, PlottingHelper::Primary);
//     PlottingHelper::SetLineStyle(pNProtonsRecoPlotGenieGeneric, PlottingHelper::Secondary);
//     maxY = 1.05 * std::max(pNProtonsRecoPlotGenieGeneric->GetMaximum(), pNProtonsRecoPlotNuWroGeneric->GetMaximum());
//     pNProtonsRecoPlotGenieGeneric->GetYaxis()->SetRangeUser(0, maxY);
//     pNProtonsRecoPlotNuWroGeneric->GetYaxis()->SetRangeUser(0, maxY);
//     pNProtonsRecoPlotGenieGeneric->Draw("hist");
//     pNProtonsRecoPlotNuWroGeneric->Draw("hist same");
//     PlottingHelper::SaveCanvas(pCanvas, prefix+"_nProtonsRecoGeneric");

//     PlottingHelper::SetLineStyle(pNProtonsRecoPlotNuWroSideband, PlottingHelper::Primary);
//     PlottingHelper::SetLineStyle(pNProtonsRecoPlotGenieSideband, PlottingHelper::Secondary);
//     maxY = 1.05 * std::max(pNProtonsRecoPlotGenieSideband->GetMaximum(), pNProtonsRecoPlotNuWroSideband->GetMaximum());
//     pNProtonsRecoPlotGenieSideband->GetYaxis()->SetRangeUser(0, maxY);
//     pNProtonsRecoPlotNuWroSideband->GetYaxis()->SetRangeUser(0, maxY);
//     pNProtonsRecoPlotGenieSideband->Draw("hist");
//     pNProtonsRecoPlotNuWroSideband->Draw("hist same");
//     PlottingHelper::SaveCanvas(pCanvas, prefix+"_nProtonsRecoSideband");

//     PlottingHelper::SetLineStyle(pNProtonsTruePlotNuWroGeneric, PlottingHelper::Primary);
//     PlottingHelper::SetLineStyle(pNProtonsTruePlotGenieGeneric, PlottingHelper::Secondary);
//     maxY = 1.05 * std::max(pNProtonsTruePlotGenieGeneric->GetMaximum(), pNProtonsTruePlotNuWroGeneric->GetMaximum());
//     pNProtonsTruePlotGenieGeneric->GetYaxis()->SetRangeUser(0, maxY);
//     pNProtonsTruePlotNuWroGeneric->GetYaxis()->SetRangeUser(0, maxY);
//     pNProtonsTruePlotGenieGeneric->Draw("hist");
//     pNProtonsTruePlotNuWroGeneric->Draw("hist same");
//     PlottingHelper::SaveCanvas(pCanvas, prefix+"_nProtonsTrueGeneric");

//     PlottingHelper::SetLineStyle(pNProtonsTruePlotNuWroSideband, PlottingHelper::Primary);
//     PlottingHelper::SetLineStyle(pNProtonsTruePlotGenieSideband, PlottingHelper::Secondary);
//     maxY = 1.05 * std::max(pNProtonsTruePlotGenieSideband->GetMaximum(), pNProtonsTruePlotNuWroSideband->GetMaximum());
//     pNProtonsTruePlotGenieSideband->GetYaxis()->SetRangeUser(0, maxY);
//     pNProtonsTruePlotNuWroSideband->GetYaxis()->SetRangeUser(0, maxY);
//     pNProtonsTruePlotGenieSideband->Draw("hist");
//     pNProtonsTruePlotNuWroSideband->Draw("hist same");
//     PlottingHelper::SaveCanvas(pCanvas, prefix+"_nProtonsTrueSideband");


//     //################################# protonMomentum
//     PlottingHelper::SetLineStyle(pProtonMomentumTruePlotNuWroAll, PlottingHelper::Primary);
//     PlottingHelper::SetLineStyle(pProtonMomentumTruePlotGenieAll, PlottingHelper::Secondary);
//     maxY = 1.05 * std::max(pProtonMomentumTruePlotGenieAll->GetMaximum(), pProtonMomentumTruePlotNuWroAll->GetMaximum());
//     pProtonMomentumTruePlotGenieAll->GetYaxis()->SetRangeUser(0, maxY);
//     pProtonMomentumTruePlotNuWroAll->GetYaxis()->SetRangeUser(0, maxY);
//     pProtonMomentumTruePlotGenieAll->Draw("hist");
//     pProtonMomentumTruePlotNuWroAll->Draw("hist same");
//     PlottingHelper::SaveCanvas(pCanvas, prefix+"_protonMomentumTrueAll");

//     PlottingHelper::SetLineStyle(pProtonMomentumRecoPlotNuWroGeneric, PlottingHelper::Primary);
//     PlottingHelper::SetLineStyle(pProtonMomentumRecoPlotGenieGeneric, PlottingHelper::Secondary);
//     maxY = 1.05 * std::max(pProtonMomentumRecoPlotGenieGeneric->GetMaximum(), pProtonMomentumRecoPlotNuWroGeneric->GetMaximum());
//     pProtonMomentumRecoPlotGenieGeneric->GetYaxis()->SetRangeUser(0, maxY);
//     pProtonMomentumRecoPlotNuWroGeneric->GetYaxis()->SetRangeUser(0, maxY);
//     pProtonMomentumRecoPlotGenieGeneric->Draw("hist");
//     pProtonMomentumRecoPlotNuWroGeneric->Draw("hist same");
//     PlottingHelper::SaveCanvas(pCanvas, prefix+"_protonMomentumRecoGeneric");

//     PlottingHelper::SetLineStyle(pProtonMomentumRecoPlotNuWroSideband, PlottingHelper::Primary);
//     PlottingHelper::SetLineStyle(pProtonMomentumRecoPlotGenieSideband, PlottingHelper::Secondary);
//     maxY = 1.05 * std::max(pProtonMomentumRecoPlotGenieSideband->GetMaximum(), pProtonMomentumRecoPlotNuWroSideband->GetMaximum());
//     pProtonMomentumRecoPlotGenieSideband->GetYaxis()->SetRangeUser(0, maxY);
//     pProtonMomentumRecoPlotNuWroSideband->GetYaxis()->SetRangeUser(0, maxY);
//     pProtonMomentumRecoPlotGenieSideband->Draw("hist");
//     pProtonMomentumRecoPlotNuWroSideband->Draw("hist same");
//     PlottingHelper::SaveCanvas(pCanvas, prefix+"_protonMomentumRecoSideband");

//     PlottingHelper::SetLineStyle(pProtonMomentumTruePlotNuWroGeneric, PlottingHelper::Primary);
//     PlottingHelper::SetLineStyle(pProtonMomentumTruePlotGenieGeneric, PlottingHelper::Secondary);
//     maxY = 1.05 * std::max(pProtonMomentumTruePlotGenieGeneric->GetMaximum(), pProtonMomentumTruePlotNuWroGeneric->GetMaximum());
//     pProtonMomentumTruePlotGenieGeneric->GetYaxis()->SetRangeUser(0, maxY);
//     pProtonMomentumTruePlotNuWroGeneric->GetYaxis()->SetRangeUser(0, maxY);
//     pProtonMomentumTruePlotGenieGeneric->Draw("hist");
//     pProtonMomentumTruePlotNuWroGeneric->Draw("hist same");
//     PlottingHelper::SaveCanvas(pCanvas, prefix+"_protonMomentumTrueGeneric");

//     PlottingHelper::SetLineStyle(pProtonMomentumTruePlotNuWroSideband, PlottingHelper::Primary);
//     PlottingHelper::SetLineStyle(pProtonMomentumTruePlotGenieSideband, PlottingHelper::Secondary);
//     maxY = 1.05 * std::max(pProtonMomentumTruePlotGenieSideband->GetMaximum(), pProtonMomentumTruePlotNuWroSideband->GetMaximum());
//     pProtonMomentumTruePlotGenieSideband->GetYaxis()->SetRangeUser(0, maxY);
//     pProtonMomentumTruePlotNuWroSideband->GetYaxis()->SetRangeUser(0, maxY);
//     pProtonMomentumTruePlotGenieSideband->Draw("hist");
//     pProtonMomentumTruePlotNuWroSideband->Draw("hist same");
//     PlottingHelper::SaveCanvas(pCanvas, prefix+"_protonMomentumTrueSideband");


//     //################################# muonMomentum
//     PlottingHelper::SetLineStyle(pMuonMomentumTruePlotNuWroAll, PlottingHelper::Primary);
//     PlottingHelper::SetLineStyle(pMuonMomentumTruePlotGenieAll, PlottingHelper::Secondary);
//     maxY = 1.05 * std::max(pMuonMomentumTruePlotGenieAll->GetMaximum(), pMuonMomentumTruePlotNuWroAll->GetMaximum());
//     pMuonMomentumTruePlotGenieAll->GetYaxis()->SetRangeUser(0, maxY);
//     pMuonMomentumTruePlotNuWroAll->GetYaxis()->SetRangeUser(0, maxY);
//     pMuonMomentumTruePlotGenieAll->Draw("hist");
//     pMuonMomentumTruePlotNuWroAll->Draw("hist same");
//     PlottingHelper::SaveCanvas(pCanvas, prefix+"_muonMomentumTrueAll");

//     PlottingHelper::SetLineStyle(pMuonMomentumRecoPlotNuWroGeneric, PlottingHelper::Primary);
//     PlottingHelper::SetLineStyle(pMuonMomentumRecoPlotGenieGeneric, PlottingHelper::Secondary);
//     maxY = 1.05 * std::max(pMuonMomentumRecoPlotGenieGeneric->GetMaximum(), pMuonMomentumRecoPlotNuWroGeneric->GetMaximum());
//     pMuonMomentumRecoPlotGenieGeneric->GetYaxis()->SetRangeUser(0, maxY);
//     pMuonMomentumRecoPlotNuWroGeneric->GetYaxis()->SetRangeUser(0, maxY);
//     pMuonMomentumRecoPlotGenieGeneric->Draw("hist");
//     pMuonMomentumRecoPlotNuWroGeneric->Draw("hist same");
//     PlottingHelper::SaveCanvas(pCanvas, prefix+"_muonMomentumRecoGeneric");

//     PlottingHelper::SetLineStyle(pMuonMomentumRecoPlotNuWroSideband, PlottingHelper::Primary);
//     PlottingHelper::SetLineStyle(pMuonMomentumRecoPlotGenieSideband, PlottingHelper::Secondary);
//     maxY = 1.05 * std::max(pMuonMomentumRecoPlotGenieSideband->GetMaximum(), pMuonMomentumRecoPlotNuWroSideband->GetMaximum());
//     pMuonMomentumRecoPlotGenieSideband->GetYaxis()->SetRangeUser(0, maxY);
//     pMuonMomentumRecoPlotNuWroSideband->GetYaxis()->SetRangeUser(0, maxY);
//     pMuonMomentumRecoPlotGenieSideband->Draw("hist");
//     pMuonMomentumRecoPlotNuWroSideband->Draw("hist same");
//     PlottingHelper::SaveCanvas(pCanvas, prefix+"_muonMomentumRecoSideband");

//     PlottingHelper::SetLineStyle(pMuonMomentumTruePlotNuWroGeneric, PlottingHelper::Primary);
//     PlottingHelper::SetLineStyle(pMuonMomentumTruePlotGenieGeneric, PlottingHelper::Secondary);
//     maxY = 1.05 * std::max(pMuonMomentumTruePlotGenieGeneric->GetMaximum(), pMuonMomentumTruePlotNuWroGeneric->GetMaximum());
//     pMuonMomentumTruePlotGenieGeneric->GetYaxis()->SetRangeUser(0, maxY);
//     pMuonMomentumTruePlotNuWroGeneric->GetYaxis()->SetRangeUser(0, maxY);
//     pMuonMomentumTruePlotGenieGeneric->Draw("hist");
//     pMuonMomentumTruePlotNuWroGeneric->Draw("hist same");
//     PlottingHelper::SaveCanvas(pCanvas, prefix+"_muonMomentumTrueGeneric");

//     PlottingHelper::SetLineStyle(pMuonMomentumTruePlotNuWroSideband, PlottingHelper::Primary);
//     PlottingHelper::SetLineStyle(pMuonMomentumTruePlotGenieSideband, PlottingHelper::Secondary);
//     maxY = 1.05 * std::max(pMuonMomentumTruePlotGenieSideband->GetMaximum(), pMuonMomentumTruePlotNuWroSideband->GetMaximum());
//     pMuonMomentumTruePlotGenieSideband->GetYaxis()->SetRangeUser(0, maxY);
//     pMuonMomentumTruePlotNuWroSideband->GetYaxis()->SetRangeUser(0, maxY);
//     pMuonMomentumTruePlotGenieSideband->Draw("hist");
//     pMuonMomentumTruePlotNuWroSideband->Draw("hist same");
//     PlottingHelper::SaveCanvas(pCanvas, prefix+"_muonMomentumTrueSideband");


//     //################################# protonCosTheta
//     PlottingHelper::SetLineStyle(pProtonCosThetaTruePlotNuWroAll, PlottingHelper::Primary);
//     PlottingHelper::SetLineStyle(pProtonCosThetaTruePlotGenieAll, PlottingHelper::Secondary);
//     maxY = 1.05 * std::max(pProtonCosThetaTruePlotGenieAll->GetMaximum(), pProtonCosThetaTruePlotNuWroAll->GetMaximum());
//     pProtonCosThetaTruePlotGenieAll->GetYaxis()->SetRangeUser(0, maxY);
//     pProtonCosThetaTruePlotNuWroAll->GetYaxis()->SetRangeUser(0, maxY);
//     pProtonCosThetaTruePlotGenieAll->Draw("hist");
//     pProtonCosThetaTruePlotNuWroAll->Draw("hist same");
//     PlottingHelper::SaveCanvas(pCanvas, prefix+"_protonCosThetaTrueAll");

//     PlottingHelper::SetLineStyle(pProtonCosThetaRecoPlotNuWroGeneric, PlottingHelper::Primary);
//     PlottingHelper::SetLineStyle(pProtonCosThetaRecoPlotGenieGeneric, PlottingHelper::Secondary);
//     maxY = 1.05 * std::max(pProtonCosThetaRecoPlotGenieGeneric->GetMaximum(), pProtonCosThetaRecoPlotNuWroGeneric->GetMaximum());
//     pProtonCosThetaRecoPlotGenieGeneric->GetYaxis()->SetRangeUser(0, maxY);
//     pProtonCosThetaRecoPlotNuWroGeneric->GetYaxis()->SetRangeUser(0, maxY);
//     pProtonCosThetaRecoPlotGenieGeneric->Draw("hist");
//     pProtonCosThetaRecoPlotNuWroGeneric->Draw("hist same");
//     PlottingHelper::SaveCanvas(pCanvas, prefix+"_protonCosThetaRecoGeneric");

//     PlottingHelper::SetLineStyle(pProtonCosThetaRecoPlotNuWroSideband, PlottingHelper::Primary);
//     PlottingHelper::SetLineStyle(pProtonCosThetaRecoPlotGenieSideband, PlottingHelper::Secondary);
//     maxY = 1.05 * std::max(pProtonCosThetaRecoPlotGenieSideband->GetMaximum(), pProtonCosThetaRecoPlotNuWroSideband->GetMaximum());
//     pProtonCosThetaRecoPlotGenieSideband->GetYaxis()->SetRangeUser(0, maxY);
//     pProtonCosThetaRecoPlotNuWroSideband->GetYaxis()->SetRangeUser(0, maxY);
//     pProtonCosThetaRecoPlotGenieSideband->Draw("hist");
//     pProtonCosThetaRecoPlotNuWroSideband->Draw("hist same");
//     PlottingHelper::SaveCanvas(pCanvas, prefix+"_protonCosThetaRecoSideband");

//     PlottingHelper::SetLineStyle(pProtonCosThetaTruePlotNuWroGeneric, PlottingHelper::Primary);
//     PlottingHelper::SetLineStyle(pProtonCosThetaTruePlotGenieGeneric, PlottingHelper::Secondary);
//     maxY = 1.05 * std::max(pProtonCosThetaTruePlotGenieGeneric->GetMaximum(), pProtonCosThetaTruePlotNuWroGeneric->GetMaximum());
//     pProtonCosThetaTruePlotGenieGeneric->GetYaxis()->SetRangeUser(0, maxY);
//     pProtonCosThetaTruePlotNuWroGeneric->GetYaxis()->SetRangeUser(0, maxY);
//     pProtonCosThetaTruePlotGenieGeneric->Draw("hist");
//     pProtonCosThetaTruePlotNuWroGeneric->Draw("hist same");
//     PlottingHelper::SaveCanvas(pCanvas, prefix+"_protonCosThetaTrueGeneric");

//     PlottingHelper::SetLineStyle(pProtonCosThetaTruePlotNuWroSideband, PlottingHelper::Primary);
//     PlottingHelper::SetLineStyle(pProtonCosThetaTruePlotGenieSideband, PlottingHelper::Secondary);
//     maxY = 1.05 * std::max(pProtonCosThetaTruePlotGenieSideband->GetMaximum(), pProtonCosThetaTruePlotNuWroSideband->GetMaximum());
//     pProtonCosThetaTruePlotGenieSideband->GetYaxis()->SetRangeUser(0, maxY);
//     pProtonCosThetaTruePlotNuWroSideband->GetYaxis()->SetRangeUser(0, maxY);
//     pProtonCosThetaTruePlotGenieSideband->Draw("hist");
//     pProtonCosThetaTruePlotNuWroSideband->Draw("hist same");
//     PlottingHelper::SaveCanvas(pCanvas, prefix+"_protonCosThetaTrueSideband");


//     //################################# muonCosTheta
//     PlottingHelper::SetLineStyle(pMuonCosThetaTruePlotNuWroAll, PlottingHelper::Primary);
//     PlottingHelper::SetLineStyle(pMuonCosThetaTruePlotGenieAll, PlottingHelper::Secondary);
//     maxY = 1.05 * std::max(pMuonCosThetaTruePlotGenieAll->GetMaximum(), pMuonCosThetaTruePlotNuWroAll->GetMaximum());
//     pMuonCosThetaTruePlotGenieAll->GetYaxis()->SetRangeUser(0, maxY);
//     pMuonCosThetaTruePlotNuWroAll->GetYaxis()->SetRangeUser(0, maxY);
//     pMuonCosThetaTruePlotGenieAll->Draw("hist");
//     pMuonCosThetaTruePlotNuWroAll->Draw("hist same");
//     PlottingHelper::SaveCanvas(pCanvas, prefix+"_muonCosThetaTrueAll");

//     PlottingHelper::SetLineStyle(pMuonCosThetaRecoPlotNuWroGeneric, PlottingHelper::Primary);
//     PlottingHelper::SetLineStyle(pMuonCosThetaRecoPlotGenieGeneric, PlottingHelper::Secondary);
//     maxY = 1.05 * std::max(pMuonCosThetaRecoPlotGenieGeneric->GetMaximum(), pMuonCosThetaRecoPlotNuWroGeneric->GetMaximum());
//     pMuonCosThetaRecoPlotGenieGeneric->GetYaxis()->SetRangeUser(0, maxY);
//     pMuonCosThetaRecoPlotNuWroGeneric->GetYaxis()->SetRangeUser(0, maxY);
//     pMuonCosThetaRecoPlotGenieGeneric->Draw("hist");
//     pMuonCosThetaRecoPlotNuWroGeneric->Draw("hist same");
//     PlottingHelper::SaveCanvas(pCanvas, prefix+"_muonCosThetaRecoGeneric");

//     PlottingHelper::SetLineStyle(pMuonCosThetaRecoPlotNuWroSideband, PlottingHelper::Primary);
//     PlottingHelper::SetLineStyle(pMuonCosThetaRecoPlotGenieSideband, PlottingHelper::Secondary);
//     maxY = 1.05 * std::max(pMuonCosThetaRecoPlotGenieSideband->GetMaximum(), pMuonCosThetaRecoPlotNuWroSideband->GetMaximum());
//     pMuonCosThetaRecoPlotGenieSideband->GetYaxis()->SetRangeUser(0, maxY);
//     pMuonCosThetaRecoPlotNuWroSideband->GetYaxis()->SetRangeUser(0, maxY);
//     pMuonCosThetaRecoPlotGenieSideband->Draw("hist");
//     pMuonCosThetaRecoPlotNuWroSideband->Draw("hist same");
//     PlottingHelper::SaveCanvas(pCanvas, prefix+"_muonCosThetaRecoSideband");

//     PlottingHelper::SetLineStyle(pMuonCosThetaTruePlotNuWroGeneric, PlottingHelper::Primary);
//     PlottingHelper::SetLineStyle(pMuonCosThetaTruePlotGenieGeneric, PlottingHelper::Secondary);
//     maxY = 1.05 * std::max(pMuonCosThetaTruePlotGenieGeneric->GetMaximum(), pMuonCosThetaTruePlotNuWroGeneric->GetMaximum());
//     pMuonCosThetaTruePlotGenieGeneric->GetYaxis()->SetRangeUser(0, maxY);
//     pMuonCosThetaTruePlotNuWroGeneric->GetYaxis()->SetRangeUser(0, maxY);
//     pMuonCosThetaTruePlotGenieGeneric->Draw("hist");
//     pMuonCosThetaTruePlotNuWroGeneric->Draw("hist same");
//     PlottingHelper::SaveCanvas(pCanvas, prefix+"_muonCosThetaTrueGeneric");

//     PlottingHelper::SetLineStyle(pMuonCosThetaTruePlotNuWroSideband, PlottingHelper::Primary);
//     PlottingHelper::SetLineStyle(pMuonCosThetaTruePlotGenieSideband, PlottingHelper::Secondary);
//     maxY = 1.05 * std::max(pMuonCosThetaTruePlotGenieSideband->GetMaximum(), pMuonCosThetaTruePlotNuWroSideband->GetMaximum());
//     pMuonCosThetaTruePlotGenieSideband->GetYaxis()->SetRangeUser(0, maxY);
//     pMuonCosThetaTruePlotNuWroSideband->GetYaxis()->SetRangeUser(0, maxY);
//     pMuonCosThetaTruePlotGenieSideband->Draw("hist");
//     pMuonCosThetaTruePlotNuWroSideband->Draw("hist same");
//     PlottingHelper::SaveCanvas(pCanvas, prefix+"_muonCosThetaTrueSideband");


//     //################################# protonPhi
//     PlottingHelper::SetLineStyle(pProtonPhiTruePlotNuWroAll, PlottingHelper::Primary);
//     PlottingHelper::SetLineStyle(pProtonPhiTruePlotGenieAll, PlottingHelper::Secondary);
//     maxY = 1.05 * std::max(pProtonPhiTruePlotGenieAll->GetMaximum(), pProtonPhiTruePlotNuWroAll->GetMaximum());
//     pProtonPhiTruePlotGenieAll->GetYaxis()->SetRangeUser(0, maxY);
//     pProtonPhiTruePlotNuWroAll->GetYaxis()->SetRangeUser(0, maxY);
//     pProtonPhiTruePlotGenieAll->Draw("hist");
//     pProtonPhiTruePlotNuWroAll->Draw("hist same");
//     PlottingHelper::SaveCanvas(pCanvas, prefix+"_protonPhiTrueAll");

//     PlottingHelper::SetLineStyle(pProtonPhiRecoPlotNuWroGeneric, PlottingHelper::Primary);
//     PlottingHelper::SetLineStyle(pProtonPhiRecoPlotGenieGeneric, PlottingHelper::Secondary);
//     maxY = 1.05 * std::max(pProtonPhiRecoPlotGenieGeneric->GetMaximum(), pProtonPhiRecoPlotNuWroGeneric->GetMaximum());
//     pProtonPhiRecoPlotGenieGeneric->GetYaxis()->SetRangeUser(0, maxY);
//     pProtonPhiRecoPlotNuWroGeneric->GetYaxis()->SetRangeUser(0, maxY);
//     pProtonPhiRecoPlotGenieGeneric->Draw("hist");
//     pProtonPhiRecoPlotNuWroGeneric->Draw("hist same");
//     PlottingHelper::SaveCanvas(pCanvas, prefix+"_protonPhiRecoGeneric");

//     PlottingHelper::SetLineStyle(pProtonPhiRecoPlotNuWroSideband, PlottingHelper::Primary);
//     PlottingHelper::SetLineStyle(pProtonPhiRecoPlotGenieSideband, PlottingHelper::Secondary);
//     maxY = 1.05 * std::max(pProtonPhiRecoPlotGenieSideband->GetMaximum(), pProtonPhiRecoPlotNuWroSideband->GetMaximum());
//     pProtonPhiRecoPlotGenieSideband->GetYaxis()->SetRangeUser(0, maxY);
//     pProtonPhiRecoPlotNuWroSideband->GetYaxis()->SetRangeUser(0, maxY);
//     pProtonPhiRecoPlotGenieSideband->Draw("hist");
//     pProtonPhiRecoPlotNuWroSideband->Draw("hist same");
//     PlottingHelper::SaveCanvas(pCanvas, prefix+"_protonPhiRecoSideband");

//     PlottingHelper::SetLineStyle(pProtonPhiTruePlotNuWroGeneric, PlottingHelper::Primary);
//     PlottingHelper::SetLineStyle(pProtonPhiTruePlotGenieGeneric, PlottingHelper::Secondary);
//     maxY = 1.05 * std::max(pProtonPhiTruePlotGenieGeneric->GetMaximum(), pProtonPhiTruePlotNuWroGeneric->GetMaximum());
//     pProtonPhiTruePlotGenieGeneric->GetYaxis()->SetRangeUser(0, maxY);
//     pProtonPhiTruePlotNuWroGeneric->GetYaxis()->SetRangeUser(0, maxY);
//     pProtonPhiTruePlotGenieGeneric->Draw("hist");
//     pProtonPhiTruePlotNuWroGeneric->Draw("hist same");
//     PlottingHelper::SaveCanvas(pCanvas, prefix+"_protonPhiTrueGeneric");

//     PlottingHelper::SetLineStyle(pProtonPhiTruePlotNuWroSideband, PlottingHelper::Primary);
//     PlottingHelper::SetLineStyle(pProtonPhiTruePlotGenieSideband, PlottingHelper::Secondary);
//     maxY = 1.05 * std::max(pProtonPhiTruePlotGenieSideband->GetMaximum(), pProtonPhiTruePlotNuWroSideband->GetMaximum());
//     pProtonPhiTruePlotGenieSideband->GetYaxis()->SetRangeUser(0, maxY);
//     pProtonPhiTruePlotNuWroSideband->GetYaxis()->SetRangeUser(0, maxY);
//     pProtonPhiTruePlotGenieSideband->Draw("hist");
//     pProtonPhiTruePlotNuWroSideband->Draw("hist same");
//     PlottingHelper::SaveCanvas(pCanvas, prefix+"_protonPhiTrueSideband");


//     //################################# muonPhi
//     PlottingHelper::SetLineStyle(pMuonPhiTruePlotNuWroAll, PlottingHelper::Primary);
//     PlottingHelper::SetLineStyle(pMuonPhiTruePlotGenieAll, PlottingHelper::Secondary);
//     maxY = 1.05 * std::max(pMuonPhiTruePlotGenieAll->GetMaximum(), pMuonPhiTruePlotNuWroAll->GetMaximum());
//     pMuonPhiTruePlotGenieAll->GetYaxis()->SetRangeUser(0, maxY);
//     pMuonPhiTruePlotNuWroAll->GetYaxis()->SetRangeUser(0, maxY);
//     pMuonPhiTruePlotGenieAll->Draw("hist");
//     pMuonPhiTruePlotNuWroAll->Draw("hist same");
//     PlottingHelper::SaveCanvas(pCanvas, prefix+"_muonPhiTrueAll");

//     PlottingHelper::SetLineStyle(pMuonPhiRecoPlotNuWroGeneric, PlottingHelper::Primary);
//     PlottingHelper::SetLineStyle(pMuonPhiRecoPlotGenieGeneric, PlottingHelper::Secondary);
//     maxY = 1.05 * std::max(pMuonPhiRecoPlotGenieGeneric->GetMaximum(), pMuonPhiRecoPlotNuWroGeneric->GetMaximum());
//     pMuonPhiRecoPlotGenieGeneric->GetYaxis()->SetRangeUser(0, maxY);
//     pMuonPhiRecoPlotNuWroGeneric->GetYaxis()->SetRangeUser(0, maxY);
//     pMuonPhiRecoPlotGenieGeneric->Draw("hist");
//     pMuonPhiRecoPlotNuWroGeneric->Draw("hist same");
//     PlottingHelper::SaveCanvas(pCanvas, prefix+"_muonPhiRecoGeneric");

//     PlottingHelper::SetLineStyle(pMuonPhiRecoPlotNuWroSideband, PlottingHelper::Primary);
//     PlottingHelper::SetLineStyle(pMuonPhiRecoPlotGenieSideband, PlottingHelper::Secondary);
//     maxY = 1.05 * std::max(pMuonPhiRecoPlotGenieSideband->GetMaximum(), pMuonPhiRecoPlotNuWroSideband->GetMaximum());
//     pMuonPhiRecoPlotGenieSideband->GetYaxis()->SetRangeUser(0, maxY);
//     pMuonPhiRecoPlotNuWroSideband->GetYaxis()->SetRangeUser(0, maxY);
//     pMuonPhiRecoPlotGenieSideband->Draw("hist");
//     pMuonPhiRecoPlotNuWroSideband->Draw("hist same");
//     PlottingHelper::SaveCanvas(pCanvas, prefix+"_muonPhiRecoSideband");

//     PlottingHelper::SetLineStyle(pMuonPhiTruePlotNuWroGeneric, PlottingHelper::Primary);
//     PlottingHelper::SetLineStyle(pMuonPhiTruePlotGenieGeneric, PlottingHelper::Secondary);
//     maxY = 1.05 * std::max(pMuonPhiTruePlotGenieGeneric->GetMaximum(), pMuonPhiTruePlotNuWroGeneric->GetMaximum());
//     pMuonPhiTruePlotGenieGeneric->GetYaxis()->SetRangeUser(0, maxY);
//     pMuonPhiTruePlotNuWroGeneric->GetYaxis()->SetRangeUser(0, maxY);
//     pMuonPhiTruePlotGenieGeneric->Draw("hist");
//     pMuonPhiTruePlotNuWroGeneric->Draw("hist same");
//     PlottingHelper::SaveCanvas(pCanvas, prefix+"_muonPhiTrueGeneric");

//     PlottingHelper::SetLineStyle(pMuonPhiTruePlotNuWroSideband, PlottingHelper::Primary);
//     PlottingHelper::SetLineStyle(pMuonPhiTruePlotGenieSideband, PlottingHelper::Secondary);
//     maxY = 1.05 * std::max(pMuonPhiTruePlotGenieSideband->GetMaximum(), pMuonPhiTruePlotNuWroSideband->GetMaximum());
//     pMuonPhiTruePlotGenieSideband->GetYaxis()->SetRangeUser(0, maxY);
//     pMuonPhiTruePlotNuWroSideband->GetYaxis()->SetRangeUser(0, maxY);
//     pMuonPhiTruePlotGenieSideband->Draw("hist");
//     pMuonPhiTruePlotNuWroSideband->Draw("hist same");
//     PlottingHelper::SaveCanvas(pCanvas, prefix+"_muonPhiTrueSideband");
// }

// } // namespace ubcc1pi_macros
