/**
 *  @file  ubcc1pi_standalone/Macros/PlotMuonRecoVariables.cxx
 *
 *  @brief The implementation file of the PlotMuonRecoVariables macro
 */

#include "ubcc1pi_standalone/Macros/Macros.h"

#include "ubcc1pi_standalone/Objects/FileReader.h"

#include "ubcc1pi_standalone/Helpers/PlottingHelper.h"
#include "ubcc1pi_standalone/Helpers/AnalysisHelper.h"
#include "ubcc1pi_standalone/Helpers/SelectionHelper.h"
#include "ubcc1pi_standalone/Helpers/NormalisationHelper.h"

#include <TH2F.h>
#include <TMath.h>

using namespace ubcc1pi;

namespace ubcc1pi_macros
{

void MakeSidebandSamplePlots(const Config &config)
{
    // Setup the input files
    std::vector< std::tuple<AnalysisHelper::SampleType, std::string, float> > inputData;

    inputData.emplace_back(AnalysisHelper::Overlay, config.files.overlaysFileName, NormalisationHelper::GetOverlaysNormalisation(config));
    inputData.emplace_back(AnalysisHelper::Dirt,    config.files.dirtFileName,     NormalisationHelper::GetDirtNormalisation(config));
    inputData.emplace_back(AnalysisHelper::DataEXT, config.files.dataEXTFileName,  NormalisationHelper::GetDataEXTNormalisation(config));
    inputData.emplace_back(AnalysisHelper::DataBNB, config.files.dataBNBFileName,  1.f);

    // Get the selections
    auto selection = SelectionHelper::GetSelection("CC0pi");
    const auto cuts = selection.GetCuts();
    auto CC1piSelection = SelectionHelper::GetSelection("Default");
    const auto CC1piCuts = CC1piSelection.GetCuts();

    const auto featureNames = BDTHelper::ParticleBDTFeatureNames;
    // Get the muon BDT
    const auto muonFeatureNames = BDTHelper::MuonBDTFeatureNames;
    BDTHelper::BDT muonBDT("muon", muonFeatureNames);

    // Get the proton BDT
    const auto protonFeatureNames = BDTHelper::ProtonBDTFeatureNames;
    BDTHelper::BDT protonBDT("proton", protonFeatureNames);

    // Get the golden pion BDT
    const auto goldenPionFeatureNames = BDTHelper::GoldenPionBDTFeatureNames;
    BDTHelper::BDT goldenpionBDT("goldenPion",goldenPionFeatureNames);

    //
    // Setup the plots
    //
    const std::string yLabel = "Number of particles (norm. to bin width)";
    std::vector<PlottingHelper::MultiPlot> muonMomentumPlots, muonCosThetaPlots, muonPhiPlots;
    std::vector<PlottingHelper::MultiPlot> protonMomentumPlots, protonCosThetaPlots, protonPhiPlots;
    for (const auto &cut : cuts)
    {
        // muonMomentumPlots.emplace_back("Muon momentum / GeV", yLabel, config.global.muonMomentum.binEdges, true, config.global.axisTitles);
        // muonCosThetaPlots.emplace_back("Muon cos(theta)", yLabel, config.global.muonCosTheta.binEdges, true, config.global.axisTitles);
        // muonPhiPlots.emplace_back("Muon phi / rad", yLabel, config.global.muonPhi.binEdges, true, config.global.axisTitles);
        //
        // protonMomentumPlots.emplace_back("Proton momentum / GeV", yLabel, config.global.pionMomentum.binEdges, true, config.global.axisTitles);
        // protonCosThetaPlots.emplace_back("Proton cos(theta)", yLabel, config.global.pionCosTheta.binEdges, true, config.global.axisTitles);
        // protonPhiPlots.emplace_back("Proton phi / rad", yLabel, config.global.pionPhi.binEdges, true, config.global.axisTitles);

        muonMomentumPlots.emplace_back("Muon momentum / GeV", yLabel, 10, 0, 1.5, true, config.global.axisTitles);
        muonCosThetaPlots.emplace_back("Muon cos(theta)", yLabel, 10, -1, 1, true, config.global.axisTitles);
        muonPhiPlots.emplace_back("Muon phi / rad", yLabel, 15, -TMath::Pi(), TMath::Pi(), true, config.global.axisTitles);

        protonMomentumPlots.emplace_back("Proton momentum / GeV", yLabel, 10, 0, 1, true, config.global.axisTitles);
        protonCosThetaPlots.emplace_back("Proton cos(theta)", yLabel, 10, -1, 1, true, config.global.axisTitles);
        protonPhiPlots.emplace_back("Proton phi / rad", yLabel, 15, -TMath::Pi(), TMath::Pi(), true, config.global.axisTitles);
    }

    // Kinematic plots for true CC0pi1p events only, in CC1pi and CC0pi selections
    std::map<std::string, std::shared_ptr<TH1F>> TrueCC0pi_MuonMomentum, TrueCC0pi_MuonCosTheta, TrueCC0pi_MuonPhi, TrueCC0pi_ProtonMomentum, TrueCC0pi_ProtonCosTheta, TrueCC0pi_ProtonPhi;

    // TrueCC0pi_MuonMomentum.emplace("SelCC1pi", new TH1F("TrueCC0pi_SelCC1pi_MuonMomentum",(string(";Muon Momentum / GeV;")+yLabel).c_str(), config.global.muonMomentum.binEdges.size()-1,config.global.muonMomentum.binEdges.data()));
    // TrueCC0pi_MuonCosTheta.emplace("SelCC1pi", new TH1F("TrueCC0pi_SelCC1pi_MuonCosTheta",(string(";Muon cos(theta);")+yLabel).c_str(), config.global.muonCosTheta.binEdges.size()-1,config.global.muonCosTheta.binEdges.data()));
    // TrueCC0pi_MuonPhi.emplace("SelCC1pi", new TH1F("TrueCC0pi_SelCC1pi_MuonPhi",(string(";Muon phi / rad;")+yLabel).c_str(), config.global.muonPhi.binEdges.size()-1,config.global.muonPhi.binEdges.data()));
    // TrueCC0pi_ProtonMomentum.emplace("SelCC1pi", new TH1F("TrueCC0pi_SelCC1pi_ProtonMomentum",(string(";Proton Momentum / GeV;")+yLabel).c_str(), config.global.pionMomentum.binEdges.size()-1,config.global.pionMomentum.binEdges.data()));
    // TrueCC0pi_ProtonCosTheta.emplace("SelCC1pi", new TH1F("TrueCC0pi_SelCC1pi_ProtonCosTheta",(string(";Proton cos(theta);")+yLabel).c_str(), config.global.pionCosTheta.binEdges.size()-1,config.global.pionCosTheta.binEdges.data()));
    // TrueCC0pi_ProtonPhi.emplace("SelCC1pi", new TH1F("TrueCC0pi_SelCC1pi_ProtonPhi",(string(";Proton phi / rad;")+yLabel).c_str(), config.global.pionPhi.binEdges.size()-1,config.global.pionPhi.binEdges.data()));
    //
    // TrueCC0pi_MuonMomentum.emplace("SelCC0pi", new TH1F("TrueCC0pi_SelCC0pi_MuonMomentum",(string(";Muon Momentum / GeV;")+yLabel).c_str(), config.global.muonMomentum.binEdges.size()-1,config.global.muonMomentum.binEdges.data()));
    // TrueCC0pi_MuonCosTheta.emplace("SelCC0pi", new TH1F("TrueCC0pi_SelCC0pi_MuonCosTheta",(string(";Muon cos(theta);")+yLabel).c_str(), config.global.muonCosTheta.binEdges.size()-1,config.global.muonCosTheta.binEdges.data()));
    // TrueCC0pi_MuonPhi.emplace("SelCC0pi", new TH1F("TrueCC0pi_SelCC0pi_MuonPhi",(string(";Muon phi / rad;")+yLabel).c_str(), config.global.muonPhi.binEdges.size()-1,config.global.muonPhi.binEdges.data()));
    // TrueCC0pi_ProtonMomentum.emplace("SelCC0pi", new TH1F("TrueCC0pi_SelCC0pi_ProtonMomentum",(string(";Proton Momentum / GeV;")+yLabel).c_str(), config.global.pionMomentum.binEdges.size()-1,config.global.pionMomentum.binEdges.data()));
    // TrueCC0pi_ProtonCosTheta.emplace("SelCC0pi", new TH1F("TrueCC0pi_SelCC0pi_ProtonCosTheta",(string(";Proton cos(theta);")+yLabel).c_str(), config.global.pionCosTheta.binEdges.size()-1,config.global.pionCosTheta.binEdges.data()));
    // TrueCC0pi_ProtonPhi.emplace("SelCC0pi", new TH1F("TrueCC0pi_SelCC0pi_ProtonPhi",(string(";Proton phi / rad;")+yLabel).c_str(), config.global.pionPhi.binEdges.size()-1,config.global.pionPhi.binEdges.data()));

    TrueCC0pi_MuonMomentum.emplace("SelCC1pi", new TH1F("TrueCC0pi_SelCC1pi_MuonMomentum",(string(";Muon Momentum / GeV;")+yLabel).c_str(), 10, 0, 1.5));
    TrueCC0pi_MuonCosTheta.emplace("SelCC1pi", new TH1F("TrueCC0pi_SelCC1pi_MuonCosTheta",(string(";Muon cos(theta);")+yLabel).c_str(), 10, -1, 1));
    TrueCC0pi_MuonPhi.emplace("SelCC1pi", new TH1F("TrueCC0pi_SelCC1pi_MuonPhi",(string(";Muon phi / rad;")+yLabel).c_str(), 15, -TMath::Pi(), TMath::Pi()));
    TrueCC0pi_ProtonMomentum.emplace("SelCC1pi", new TH1F("TrueCC0pi_SelCC1pi_ProtonMomentum",(string(";Proton Momentum / GeV;")+yLabel).c_str(), 15, 0, 1.5));
    TrueCC0pi_ProtonCosTheta.emplace("SelCC1pi", new TH1F("TrueCC0pi_SelCC1pi_ProtonCosTheta",(string(";Proton cos(theta);")+yLabel).c_str(), 10, -1, 1));
    TrueCC0pi_ProtonPhi.emplace("SelCC1pi", new TH1F("TrueCC0pi_SelCC1pi_ProtonPhi",(string(";Proton phi / rad;")+yLabel).c_str(), 15, -TMath::Pi(), TMath::Pi()));

    TrueCC0pi_MuonMomentum.emplace("SelCC0pi", new TH1F("TrueCC0pi_SelCC0pi_MuonMomentum",(string(";Muon Momentum / GeV;")+yLabel).c_str(), 10, 0, 1.5));
    TrueCC0pi_MuonCosTheta.emplace("SelCC0pi", new TH1F("TrueCC0pi_SelCC0pi_MuonCosTheta",(string(";Muon cos(theta);")+yLabel).c_str(), 10, -1, 1));
    TrueCC0pi_MuonPhi.emplace("SelCC0pi", new TH1F("TrueCC0pi_SelCC0pi_MuonPhi",(string(";Muon phi / rad;")+yLabel).c_str(), 15, -TMath::Pi(), TMath::Pi()));
    TrueCC0pi_ProtonMomentum.emplace("SelCC0pi", new TH1F("TrueCC0pi_SelCC0pi_ProtonMomentum",(string(";Proton Momentum / GeV;")+yLabel).c_str(), 15, 0, 1.5));
    TrueCC0pi_ProtonCosTheta.emplace("SelCC0pi", new TH1F("TrueCC0pi_SelCC0pi_ProtonCosTheta",(string(";Proton cos(theta);")+yLabel).c_str(), 10, -1, 1));
    TrueCC0pi_ProtonPhi.emplace("SelCC0pi", new TH1F("TrueCC0pi_SelCC0pi_ProtonPhi",(string(";Proton phi / rad;")+yLabel).c_str(), 15, -TMath::Pi(), TMath::Pi()));


    std::map<std::string, std::shared_ptr<TH2F>>  TrueCC0pi_ProtonMomentumCosTheta, TrueCC0pi_ProtonCosThetaPhi, TrueCC0pi_ProtonMomentumPhi;

    TrueCC0pi_ProtonMomentumCosTheta.emplace("SelCC1pi", new TH2F("TrueCC0pi_SelCC1pi_ProtonMomentumCosTheta",";True Momentum / GeV;True cos(theta) / rad", 30,0,1.5,20, -1, 1));
    TrueCC0pi_ProtonCosThetaPhi.emplace("SelCC1pi", new TH2F("TrueCC0pi_SelCC1pi_ProtonCosThetaPhi",";True cos(theta);True phi / rad", 20,-1,1.5,15, -TMath::Pi(), TMath::Pi()));
    TrueCC0pi_ProtonMomentumPhi.emplace("SelCC1pi", new TH2F("TrueCC0pi_SelCC1pi_ProtonMomentumPhi",";True Momentum / GeV;True phi / rad", 20,0,1.5,15, -TMath::Pi(), TMath::Pi()));

    TrueCC0pi_ProtonMomentumCosTheta.emplace("SelCC0pi", new TH2F("TrueCC0pi_SelCC0pi_ProtonMomentumCosTheta",";Proton Momentum / GeV;Proton cos(theta) / rad", 30,0,1.5,20, -1, 1));
    TrueCC0pi_ProtonCosThetaPhi.emplace("SelCC0pi", new TH2F("TrueCC0pi_SelCC0pi_ProtonCosThetaPhi",";Proton cos(theta);Proton phi / rad", 20,-1,1.5,15, -TMath::Pi(), TMath::Pi()));
    TrueCC0pi_ProtonMomentumPhi.emplace("SelCC0pi", new TH2F("TrueCC0pi_SelCC0pi_ProtonMomentumPhi",";Proton Momentum / GeV;Proton phi / rad", 20,0,1.5,15, -TMath::Pi(), TMath::Pi()));

    // Plot 2D proton kinematics vs PID variables
    // Proton BDT score
    // Muon BDT score
    // Golden pion BDT score
    // logBragg_pToMIP
    // logBragg_piToMIP
    // truncMeandEdx
    // protonForwrd likelihood
    // muonForward likelihood
    // nDescendents
    // nSpacePointsNearEnd
    // wiggliness
    // trackScore
    std::vector<TH2F*> TrueCC0pi_SelCC0pi_ProtonPIDMomentum, TrueCC0pi_SelCC0pi_ProtonPIDCosTheta, TrueCC0pi_SelCC0pi_ProtonPIDPhi;
    std::vector<TH2F*> TrueCC0pi_SelCC1pi_ProtonPIDMomentum, TrueCC0pi_SelCC1pi_ProtonPIDCosTheta, TrueCC0pi_SelCC1pi_ProtonPIDPhi;
    for (size_t i_feat=0; i_feat<featureNames.size()+3; i_feat++){
        std::string xlabel;
        std::vector<float> binning;

        if (i_feat < featureNames.size()){
            if (featureNames.at(i_feat) == "logBragg_pToMIP"){
                xlabel = "log(L_p / L_MIP)";
                binning = {30,-9,7};
            }
            if (featureNames.at(i_feat) == "logBragg_piToMIP"){
                xlabel = "log(L_pi / L_MIP)";
                binning = {30,-4,6};
            }
            if (featureNames.at(i_feat) == "truncMeandEdx"){
                xlabel = "Truncated Mean dEdx";
                binning = {20,0,10};
            }
            if (featureNames.at(i_feat) == "protonForward"){
                xlabel = "Proton forward likelihood";
                binning = {30,0.42,0.62};
            }
            if (featureNames.at(i_feat) == "muonForward"){
                xlabel = "Muon forward likelihood";
                binning = {30,0.35,0.65};
            }
            if (featureNames.at(i_feat) == "nDescendents"){
                xlabel = "Number of descendents";
                binning = {4,0,4};
            }
            if (featureNames.at(i_feat) == "nSpacePointsNearEnd"){
                xlabel = "Number of spacepoints near track end";
                binning = {45,0,90};
            }
            if (featureNames.at(i_feat) == "wiggliness"){
                xlabel = "Wiggliness";
                binning = {30,0,0.1};
            }
            if (featureNames.at(i_feat) == "trackScore"){
                xlabel = "Track Score";
                binning = {20,0,1};
            }
        }
        else if (i_feat == featureNames.size()){ // Proton BDT score
            xlabel = "Proton BDT Score";
            binning = {20,-1,1};
        }
        else if (i_feat == featureNames.size()+1){ // Muon BDT score
            xlabel = "Muon BDT Score";
            binning = {20,-1,1};
        }
        else if (i_feat == featureNames.size()+2){ // Golden pion BDT score
            xlabel = "Golden pion BDT Score";
            binning = {20,-1,1};
        }


        TrueCC0pi_SelCC0pi_ProtonPIDMomentum.emplace_back(new TH2F(std::string(std::string("TrueCC0pi_SelCC0pi_ProtonMomentumPID")+std::to_string(i_feat)).c_str(),std::string(std::string(";Leading Proton Candidate ")+xlabel+std::string(";True Proton Momentum / GeV")).c_str(),binning.at(0),binning.at(1),binning.at(2),30,0,1.5));
        TrueCC0pi_SelCC0pi_ProtonPIDCosTheta.emplace_back(new TH2F(std::string(std::string("TrueCC0pi_SelCC0pi_ProtonCosThetaPID")+std::to_string(i_feat)).c_str(),std::string(std::string(";Leading Proton Candidate ")+xlabel+std::string(";True Proton cos(theta)")).c_str(),binning.at(0),binning.at(1),binning.at(2),20,-1,1));
        TrueCC0pi_SelCC0pi_ProtonPIDPhi.emplace_back(new TH2F(std::string(std::string("TrueCC0pi_SelCC0pi_ProtonPhiPID")+std::to_string(i_feat)).c_str(),std::string(std::string(";Leading Proton Candidate ")+xlabel+std::string(";True Proton Phi / rad.")).c_str(),binning.at(0),binning.at(1),binning.at(2),30,-TMath::Pi(),TMath::Pi()));

        TrueCC0pi_SelCC1pi_ProtonPIDMomentum.emplace_back(new TH2F(std::string(std::string("TrueCC0pi_SelCC1pi_ProtonMomentumPID")+std::to_string(i_feat)).c_str(),std::string(std::string(";Pion Candidate ")+xlabel+std::string(";True Momentum / GeV")).c_str(),binning.at(0),binning.at(1),binning.at(2),30,0,1.5));
        TrueCC0pi_SelCC1pi_ProtonPIDCosTheta.emplace_back(new TH2F(std::string(std::string("TrueCC0pi_SelCC1pi_ProtonCosThetaPID")+std::to_string(i_feat)).c_str(),std::string(std::string(";Pion Candidate ")+xlabel+std::string(";True cos(theta)")).c_str(),binning.at(0),binning.at(1),binning.at(2),20,-1,1));
        TrueCC0pi_SelCC1pi_ProtonPIDPhi.emplace_back(new TH2F(std::string(std::string("TrueCC0pi_SelCC1pi_ProtonPhiPID")+std::to_string(i_feat)).c_str(),std::string(std::string(";Pion Candidate ")+xlabel+std::string(";True Phi / rad.")).c_str(),binning.at(0),binning.at(1),binning.at(2),30,-TMath::Pi(),TMath::Pi()));
    }

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

            if (!pEvent->reco.passesCCInclusive())
                continue;

            const auto recoParticles = pEvent->reco.particles;

            // Run the event selection and store which cuts are passed
            const auto &[isSelected, cutsPassed, assignedPdgCodes] = selection.Execute(pEvent);

            // Get the truth and reco analysis data
            const auto weight = AnalysisHelper::GetNominalEventWeight(pEvent) * normalisation;
            const auto plotStyle = PlottingHelper::GetPlotStyle(sampleType, pEvent, config.global.useAbsPdg);

            // Get the muon reco data
            const auto muonIndex = SelectionHelper::GetMuonCandidateIndex(recoParticles, muonFeatureNames, muonBDT);
            const auto muon = recoParticles.at(muonIndex);
            const auto muonDir = TVector3(muon.directionX(), muon.directionY(), muon.directionZ()).Unit();
            const float muonCosTheta = muonDir.Z();
            const float muonPhi = std::atan2(muonDir.Y(), muonDir.X());
            const float muonMomentum = AnalysisHelper::GetMuonMomentum(muon);

            // Get the proton reco data
            float protonMomentum = -std::numeric_limits<float>::max();
            float protonCosTheta = -std::numeric_limits<float>::max();
            float protonPhi = -std::numeric_limits<float>::max();
            float pbdtResponse = -std::numeric_limits<float>::max();
            float mbdtResponse = -std::numeric_limits<float>::max();
            float gbdtResponse = -std::numeric_limits<float>::max();
            std::vector<float> allfeatures;
            std::vector<float> pfeatures;
            std::vector<float> mfeatures;
            std::vector<float> gfeatures;

            const auto protonIndex = SelectionHelper::GetLeadingProtonCandidateIndex(recoParticles, assignedPdgCodes);

            if (protonIndex!=std::numeric_limits<unsigned int>::max()){
                const auto proton = recoParticles.at(protonIndex);
                const auto protonDir = TVector3(proton.directionX(), proton.directionY(), proton.directionZ()).Unit();
                protonCosTheta = protonDir.Z();
                protonPhi = std::atan2(protonDir.Y(), protonDir.X());
                protonMomentum = AnalysisHelper::GetProtonMomentumFromRange(proton.range());

                const auto hasallFeatures = BDTHelper::GetBDTFeatures(proton, featureNames, allfeatures);

                const auto haspFeatures = BDTHelper::GetBDTFeatures(proton, protonFeatureNames, pfeatures);

                if (haspFeatures){
                    pbdtResponse = protonBDT.GetResponse(pfeatures);
                }

                const auto hasmFeatures = BDTHelper::GetBDTFeatures(proton, muonFeatureNames, mfeatures);
                if (hasmFeatures){
                    mbdtResponse = muonBDT.GetResponse(mfeatures);
                }

                const auto hasgFeatures = BDTHelper::GetBDTFeatures(proton, goldenPionFeatureNames, gfeatures);
                if (hasgFeatures){
                    gbdtResponse = goldenpionBDT.GetResponse(gfeatures);
                }

            }

            for (unsigned int iCut = 0; iCut < cuts.size(); ++iCut)
            {
                const auto passesCut = (std::find(cutsPassed.begin(), cutsPassed.end(), cuts.at(iCut)) != cutsPassed.end());

                if (!passesCut)
                    continue;

                muonMomentumPlots.at(iCut).Fill(muonMomentum, plotStyle, weight);
                muonCosThetaPlots.at(iCut).Fill(muonCosTheta, plotStyle, weight);
                muonPhiPlots.at(iCut).Fill(muonPhi, plotStyle, weight);

                if (protonIndex!=std::numeric_limits<unsigned int>::max()){
                    protonMomentumPlots.at(iCut).Fill(protonMomentum, plotStyle, weight);
                    protonCosThetaPlots.at(iCut).Fill(protonCosTheta, plotStyle, weight);
                    protonPhiPlots.at(iCut).Fill(protonPhi, plotStyle, weight);
                }
            }

            // Fill kinematic comparison and efficiency plots only for true CC0pi events that pass the full selection
            // First, check if event is true CC0pi
            if (plotStyle == PlottingHelper::NumuCC0Pi){
            // if (plotStyle == PlottingHelper::NumuCC1PiChargedGolden || plotStyle == PlottingHelper::NumuCC1PiChargedNonGolden){
                // Get the features we want to plot (truth, not reco)
                const auto truemuonidx = AnalysisHelper::GetTrueMuonIndex(pEvent->truth, config.global.useAbsPdg);
                const auto truemu = pEvent->truth.particles.at(truemuonidx);
                const auto TrueMuDir = TVector3(truemu.momentumX(), truemu.momentumY(), truemu.momentumZ()).Unit();
                const auto TrueMuCosTheta = TrueMuDir.Z();
                const auto TrueMuPhi = std::atan2(TrueMuDir.Y(),TrueMuDir.X());

                // If event passes CC0pi selection, fill plots
                if (isSelected){

                    // For proton kinematics, use true leading proton. If no true leading proton is found, don't fill any kinematics
                    auto leadpidx = AnalysisHelper::GetTrueLeadingProtonIndex(pEvent->truth, config.global.useAbsPdg, config.global.protonMomentumThreshold);

                    if (leadpidx!=std::numeric_limits<unsigned int>::max()){
                        const auto trueleadp = pEvent->truth.particles.at(leadpidx);
                        const auto dir = TVector3(trueleadp.momentumX(), trueleadp.momentumY(), trueleadp.momentumZ()).Unit();
                        const auto cosTheta = dir.Z();
                        const auto phi = std::atan2(dir.Y(),dir.X());

                        TrueCC0pi_MuonMomentum.at("SelCC0pi")->Fill(truemu.momentum());
                        TrueCC0pi_MuonCosTheta.at("SelCC0pi")->Fill(TrueMuCosTheta);
                        TrueCC0pi_MuonPhi.at("SelCC0pi")->Fill(TrueMuPhi);

                        TrueCC0pi_ProtonMomentum.at("SelCC0pi")->Fill(trueleadp.momentum());
                        TrueCC0pi_ProtonCosTheta.at("SelCC0pi")->Fill(cosTheta);
                        TrueCC0pi_ProtonPhi.at("SelCC0pi")->Fill(phi);

                        TrueCC0pi_ProtonMomentumCosTheta.at("SelCC0pi")->Fill(trueleadp.momentum(),cosTheta);
                        TrueCC0pi_ProtonCosThetaPhi.at("SelCC0pi")->Fill(cosTheta,phi);
                        TrueCC0pi_ProtonMomentumPhi.at("SelCC0pi")->Fill(trueleadp.momentum(),phi);

                        for (size_t i_feat=0; i_feat<allfeatures.size()+3;i_feat++){
                            float fillx = -std::numeric_limits<float>::max();

                            if (i_feat < allfeatures.size())
                                 fillx = allfeatures.at(i_feat);
                            else if (i_feat == allfeatures.size()) // Proton BDT score
                                fillx = pbdtResponse;
                            else if (i_feat == allfeatures.size()+1) // Muon BDT score
                                fillx = mbdtResponse;
                            else if (i_feat == allfeatures.size()+2) // Golden pion BDT score
                                fillx = gbdtResponse;

                            TrueCC0pi_SelCC0pi_ProtonPIDMomentum.at(i_feat)->Fill(fillx,trueleadp.momentum());
                            TrueCC0pi_SelCC0pi_ProtonPIDCosTheta.at(i_feat)->Fill(fillx,cosTheta);
                            TrueCC0pi_SelCC0pi_ProtonPIDPhi.at(i_feat)->Fill(fillx,phi);
                        }
                    }
                }

                // Find out if the event passes the generic CC1pi selection
                const auto &[isSelectedCC1pi, cutsPassedCC1pi, assignedPdgCodesCC1pi] = CC1piSelection.Execute(pEvent);
                const auto &isSelectedCC1piGeneric = (std::find(cutsPassedCC1pi.begin(),cutsPassedCC1pi.end(), config.global.lastCutGeneric) != cutsPassedCC1pi.end());

                // If event passes CC1pi+ selection, fill plots
                if (isSelectedCC1piGeneric){
                    // For proton kinematics, use proton that is selected as the pi+ candidate. If the pi+ candidate is a true proton, fill muon and proton kinematic plots
                    const auto recoPion = recoParticles.at(AnalysisHelper::GetParticleIndexWithPdg(assignedPdgCodesCC1pi, 211));

                    Event::Truth::Particle pionmatch;

                    try{
                        pionmatch = AnalysisHelper::GetBestMatchedTruthParticle(recoPion,pEvent->truth.particles);
                    }
                    catch(const std::logic_error &){
                        continue;
                    }


                    if (pionmatch.pdgCode()==2212){
                        TrueCC0pi_MuonMomentum.at("SelCC1pi")->Fill(truemu.momentum());
                        TrueCC0pi_MuonCosTheta.at("SelCC1pi")->Fill(TrueMuCosTheta);
                        TrueCC0pi_MuonPhi.at("SelCC1pi")->Fill(TrueMuPhi);

                        const auto dir = TVector3(pionmatch.momentumX(), pionmatch.momentumY(), pionmatch.momentumZ()).Unit();
                        const auto cosTheta = dir.Z();
                        const auto phi = std::atan2(dir.Y(),dir.X());

                        TrueCC0pi_ProtonMomentum.at("SelCC1pi")->Fill(pionmatch.momentum());
                        TrueCC0pi_ProtonCosTheta.at("SelCC1pi")->Fill(cosTheta);
                        TrueCC0pi_ProtonPhi.at("SelCC1pi")->Fill(phi);

                        TrueCC0pi_ProtonMomentumCosTheta.at("SelCC1pi")->Fill(pionmatch.momentum(),cosTheta);
                        TrueCC0pi_ProtonCosThetaPhi.at("SelCC1pi")->Fill(cosTheta,phi);
                        TrueCC0pi_ProtonMomentumPhi.at("SelCC1pi")->Fill(pionmatch.momentum(),phi);

                        std::vector<float> allfeatures;
                        std::vector<float> pfeatures;
                        std::vector<float> mfeatures;
                        std::vector<float> gfeatures;

                        const auto hasallFeatures = BDTHelper::GetBDTFeatures(recoPion, featureNames, allfeatures);

                        const auto haspFeatures = BDTHelper::GetBDTFeatures(recoPion, protonFeatureNames, pfeatures);
                        if (haspFeatures){
                            pbdtResponse = protonBDT.GetResponse(pfeatures);
                        }

                        const auto hasmFeatures = BDTHelper::GetBDTFeatures(recoPion, muonFeatureNames, mfeatures);
                        if (hasmFeatures){
                            mbdtResponse = muonBDT.GetResponse(mfeatures);
                        }

                        const auto hasgFeatures = BDTHelper::GetBDTFeatures(recoPion, goldenPionFeatureNames, gfeatures);
                        if (hasgFeatures){
                            gbdtResponse = goldenpionBDT.GetResponse(gfeatures);
                        }

                        for (size_t i_feat=0; i_feat<allfeatures.size()+3;i_feat++){
                            float fillx = -std::numeric_limits<float>::max();

                            if (i_feat < allfeatures.size())
                                 fillx = allfeatures.at(i_feat);
                            else if (i_feat == allfeatures.size()) // Proton BDT score
                                fillx = pbdtResponse;
                            else if (i_feat == allfeatures.size()+1) // Muon BDT score
                                fillx = mbdtResponse;
                            else if (i_feat == allfeatures.size()+2) // Golden pion BDT score
                                fillx = gbdtResponse;

                            TrueCC0pi_SelCC1pi_ProtonPIDMomentum.at(i_feat)->Fill(fillx,pionmatch.momentum());
                            TrueCC0pi_SelCC1pi_ProtonPIDCosTheta.at(i_feat)->Fill(fillx,cosTheta);
                            TrueCC0pi_SelCC1pi_ProtonPIDPhi.at(i_feat)->Fill(fillx,phi);
                        }
                    }
                }
            } // end if true numu CC0pi

        }
    }

    for (unsigned int iCut = 0; iCut < cuts.size(); ++iCut)
    {
        const auto &cut = cuts.at(iCut);

        const std::string suffix = std::to_string(iCut) + "_" + cut;
        muonMomentumPlots.at(iCut).SaveAsStacked("reco_cc1pi_muonMomentum_" + suffix, false, config.global.scaleByBinWidth, config.global.axisTitles);
        muonCosThetaPlots.at(iCut).SaveAsStacked("reco_cc1pi_muonCosTheta_" + suffix, false, config.global.scaleByBinWidth, config.global.axisTitles);
        muonPhiPlots.at(iCut).SaveAsStacked("reco_cc1pi_muonPhi_" + suffix, false, config.global.scaleByBinWidth, config.global.axisTitles);

        protonMomentumPlots.at(iCut).SaveAsStacked("reco_cc1pi_protonMomentum_" + suffix, false, config.global.scaleByBinWidth, config.global.axisTitles);
        protonCosThetaPlots.at(iCut).SaveAsStacked("reco_cc1pi_protonCosTheta_" + suffix, false, config.global.scaleByBinWidth, config.global.axisTitles);
        protonPhiPlots.at(iCut).SaveAsStacked("reco_cc1pi_protonPhi_" + suffix, false, config.global.scaleByBinWidth, config.global.axisTitles);
    }

    // Save kinematic comparison plots for true CC0pi
    auto pCanvas = PlottingHelper::GetCanvas();

    PlottingHelper::SetLineStyle(TrueCC0pi_MuonMomentum.at("SelCC1pi"), PlottingHelper::Primary);
    PlottingHelper::SetLineStyle(TrueCC0pi_MuonMomentum.at("SelCC0pi"), PlottingHelper::Secondary);
    TrueCC0pi_MuonMomentum.at("SelCC0pi")->Draw("hist");
    TrueCC0pi_MuonMomentum.at("SelCC1pi")->Draw("hist same");
    PlottingHelper::SaveCanvas(pCanvas,"SidebandComparisons_MuonMomentum");

    PlottingHelper::SetLineStyle(TrueCC0pi_MuonCosTheta.at("SelCC1pi"), PlottingHelper::Primary);
    PlottingHelper::SetLineStyle(TrueCC0pi_MuonCosTheta.at("SelCC0pi"), PlottingHelper::Secondary);
    TrueCC0pi_MuonCosTheta.at("SelCC0pi")->Draw("hist");
    TrueCC0pi_MuonCosTheta.at("SelCC1pi")->Draw("hist same");
    PlottingHelper::SaveCanvas(pCanvas,"SidebandComparisons_MuonCosTheta");

    PlottingHelper::SetLineStyle(TrueCC0pi_MuonPhi.at("SelCC1pi"), PlottingHelper::Primary);
    PlottingHelper::SetLineStyle(TrueCC0pi_MuonPhi.at("SelCC0pi"), PlottingHelper::Secondary);
    TrueCC0pi_MuonPhi.at("SelCC0pi")->Draw("hist");
    TrueCC0pi_MuonPhi.at("SelCC1pi")->Draw("hist same");
    PlottingHelper::SaveCanvas(pCanvas,"SidebandComparisons_MuonPhi");

    PlottingHelper::SetLineStyle(TrueCC0pi_ProtonMomentum.at("SelCC1pi"), PlottingHelper::Primary);
    PlottingHelper::SetLineStyle(TrueCC0pi_ProtonMomentum.at("SelCC0pi"), PlottingHelper::Secondary);
    TrueCC0pi_ProtonMomentum.at("SelCC0pi")->Draw("hist");
    TrueCC0pi_ProtonMomentum.at("SelCC1pi")->Draw("hist same");
    PlottingHelper::SaveCanvas(pCanvas,"SidebandComparisons_ProtonMomentum");

    PlottingHelper::SetLineStyle(TrueCC0pi_ProtonCosTheta.at("SelCC1pi"), PlottingHelper::Primary);
    PlottingHelper::SetLineStyle(TrueCC0pi_ProtonCosTheta.at("SelCC0pi"), PlottingHelper::Secondary);
    TrueCC0pi_ProtonCosTheta.at("SelCC0pi")->Draw("hist");
    TrueCC0pi_ProtonCosTheta.at("SelCC1pi")->Draw("hist same");
    PlottingHelper::SaveCanvas(pCanvas,"SidebandComparisons_ProtonCosTheta");

    PlottingHelper::SetLineStyle(TrueCC0pi_ProtonPhi.at("SelCC1pi"), PlottingHelper::Primary);
    PlottingHelper::SetLineStyle(TrueCC0pi_ProtonPhi.at("SelCC0pi"), PlottingHelper::Secondary);
    TrueCC0pi_ProtonPhi.at("SelCC0pi")->Draw("hist");
    TrueCC0pi_ProtonPhi.at("SelCC1pi")->Draw("hist same");
    PlottingHelper::SaveCanvas(pCanvas,"SidebandComparisons_ProtonPhi");

    // Now normalise plots so we can compare shapes
    TrueCC0pi_MuonMomentum.at("SelCC1pi")->Sumw2();
    TrueCC0pi_MuonMomentum.at("SelCC1pi")->Scale(1.0/TrueCC0pi_MuonMomentum.at("SelCC1pi")->Integral());
    TrueCC0pi_MuonMomentum.at("SelCC0pi")->Sumw2();
    TrueCC0pi_MuonMomentum.at("SelCC0pi")->Scale(1.0/TrueCC0pi_MuonMomentum.at("SelCC0pi")->Integral());
    TrueCC0pi_MuonMomentum.at("SelCC0pi")->GetYaxis()->SetRangeUser(0,0.2);
    TrueCC0pi_MuonMomentum.at("SelCC0pi")->Draw("hist");
    TrueCC0pi_MuonMomentum.at("SelCC1pi")->Draw("hist same");
    PlottingHelper::SaveCanvas(pCanvas,"SidebandComparisons_MuonMomentum_areanorm");

    TrueCC0pi_MuonCosTheta.at("SelCC1pi")->Sumw2();
    TrueCC0pi_MuonCosTheta.at("SelCC1pi")->Scale(1.0/TrueCC0pi_MuonCosTheta.at("SelCC1pi")->Integral());
    TrueCC0pi_MuonCosTheta.at("SelCC0pi")->Sumw2();
    TrueCC0pi_MuonCosTheta.at("SelCC0pi")->Scale(1.0/TrueCC0pi_MuonCosTheta.at("SelCC0pi")->Integral());
    TrueCC0pi_MuonCosTheta.at("SelCC0pi")->GetYaxis()->SetRangeUser(0,0.65);
    TrueCC0pi_MuonCosTheta.at("SelCC0pi")->Draw("hist");
    TrueCC0pi_MuonCosTheta.at("SelCC1pi")->Draw("hist same");
    PlottingHelper::SaveCanvas(pCanvas,"SidebandComparisons_MuonCosTheta_areanorm");

    TrueCC0pi_MuonPhi.at("SelCC1pi")->Sumw2();
    TrueCC0pi_MuonPhi.at("SelCC1pi")->Scale(1.0/TrueCC0pi_MuonPhi.at("SelCC1pi")->Integral());
    TrueCC0pi_MuonPhi.at("SelCC0pi")->Sumw2();
    TrueCC0pi_MuonPhi.at("SelCC0pi")->Scale(1.0/TrueCC0pi_MuonPhi.at("SelCC0pi")->Integral());
    TrueCC0pi_MuonPhi.at("SelCC0pi")->GetYaxis()->SetRangeUser(0,0.11);
    TrueCC0pi_MuonPhi.at("SelCC0pi")->Draw("hist");
    TrueCC0pi_MuonPhi.at("SelCC1pi")->Draw("hist same");
    PlottingHelper::SaveCanvas(pCanvas,"SidebandComparisons_MuonPhi_areanorm");

    TrueCC0pi_ProtonMomentum.at("SelCC1pi")->Sumw2();
    TrueCC0pi_ProtonMomentum.at("SelCC1pi")->Scale(1.0/TrueCC0pi_ProtonMomentum.at("SelCC1pi")->Integral());
    TrueCC0pi_ProtonMomentum.at("SelCC0pi")->Sumw2();
    TrueCC0pi_ProtonMomentum.at("SelCC0pi")->Scale(1.0/TrueCC0pi_ProtonMomentum.at("SelCC0pi")->Integral());
    TrueCC0pi_ProtonMomentum.at("SelCC0pi")->GetYaxis()->SetRangeUser(0,0.21);
    TrueCC0pi_ProtonMomentum.at("SelCC0pi")->Draw("hist");
    TrueCC0pi_ProtonMomentum.at("SelCC1pi")->Draw("hist same");
    PlottingHelper::SaveCanvas(pCanvas,"SidebandComparisons_ProtonMomentum_areanorm");

    TrueCC0pi_ProtonCosTheta.at("SelCC1pi")->Sumw2();
    TrueCC0pi_ProtonCosTheta.at("SelCC1pi")->Scale(1.0/TrueCC0pi_ProtonCosTheta.at("SelCC1pi")->Integral());
    TrueCC0pi_ProtonCosTheta.at("SelCC0pi")->Sumw2();
    TrueCC0pi_ProtonCosTheta.at("SelCC0pi")->Scale(1.0/TrueCC0pi_ProtonCosTheta.at("SelCC0pi")->Integral());
    TrueCC0pi_ProtonCosTheta.at("SelCC0pi")->GetYaxis()->SetRangeUser(0,0.5);
    TrueCC0pi_ProtonCosTheta.at("SelCC0pi")->Draw("hist");
    TrueCC0pi_ProtonCosTheta.at("SelCC1pi")->Draw("hist same");
    PlottingHelper::SaveCanvas(pCanvas,"SidebandComparisons_ProtonCosTheta_areanorm");

    TrueCC0pi_ProtonPhi.at("SelCC1pi")->Sumw2();
    TrueCC0pi_ProtonPhi.at("SelCC1pi")->Scale(1.0/TrueCC0pi_ProtonPhi.at("SelCC1pi")->Integral());
    TrueCC0pi_ProtonPhi.at("SelCC0pi")->Sumw2();
    TrueCC0pi_ProtonPhi.at("SelCC0pi")->Scale(1.0/TrueCC0pi_ProtonPhi.at("SelCC0pi")->Integral());
    TrueCC0pi_ProtonPhi.at("SelCC0pi")->GetYaxis()->SetRangeUser(0,0.13);
    TrueCC0pi_ProtonPhi.at("SelCC0pi")->Draw("hist");
    TrueCC0pi_ProtonPhi.at("SelCC1pi")->Draw("hist same");
    PlottingHelper::SaveCanvas(pCanvas,"SidebandComparisons_ProtonPhi_areanorm");

    TrueCC0pi_ProtonMomentumCosTheta.at("SelCC0pi")->Draw("colz");
    PlottingHelper::SaveCanvas(pCanvas,"SidebandComparisons_ProtonMomentumCosTheta_CC0pi");
    TrueCC0pi_ProtonCosThetaPhi.at("SelCC0pi")->Draw("colz");
    PlottingHelper::SaveCanvas(pCanvas,"SidebandComparisons_ProtonCosThetaPhi_CC0pi");
    TrueCC0pi_ProtonMomentumPhi.at("SelCC0pi")->Draw("colz");
    PlottingHelper::SaveCanvas(pCanvas,"SidebandComparisons_ProtonMomentumPhi_CC0pi");

    TrueCC0pi_ProtonMomentumCosTheta.at("SelCC1pi")->Draw("colz");
    PlottingHelper::SaveCanvas(pCanvas,"SidebandComparisons_ProtonMomentumCosTheta_CC1pi");
    TrueCC0pi_ProtonCosThetaPhi.at("SelCC1pi")->Draw("colz");
    PlottingHelper::SaveCanvas(pCanvas,"SidebandComparisons_ProtonCosThetaPhi_CC1pi");
    TrueCC0pi_ProtonMomentumPhi.at("SelCC1pi")->Draw("colz");
    PlottingHelper::SaveCanvas(pCanvas,"SidebandComparisons_ProtonMomentumPhi_CC1pi");

    pCanvas->SetLogz();
    // TrueCC0pi_ProtonBDTscoreMomentum.at("SelCC0pi")->Draw("colz");
    // PlottingHelper::SaveCanvas(pCanvas,"SidebandComparisons_BDTScoreProtonMomentum_CC0pi");
    // TrueCC0pi_ProtonBDTscoreCosTheta.at("SelCC0pi")->Draw("colz");
    // PlottingHelper::SaveCanvas(pCanvas,"SidebandComparisons_BDTScoreProtonCosTheta_CC0pi");
    // TrueCC0pi_ProtonBDTscorePhi.at("SelCC0pi")->Draw("colz");
    // PlottingHelper::SaveCanvas(pCanvas,"SidebandComparisons_BDTScoreProtonPhi_CC0pi");
    //
    // TrueCC0pi_ProtonBDTscoreMomentum.at("SelCC1pi")->Draw("colz");
    // PlottingHelper::SaveCanvas(pCanvas,"SidebandComparisons_BDTScoreProtonMomentum_CC1pi");
    // TrueCC0pi_ProtonBDTscoreCosTheta.at("SelCC1pi")->Draw("colz");
    // PlottingHelper::SaveCanvas(pCanvas,"SidebandComparisons_BDTScoreProtonCosTheta_CC1pi");
    // TrueCC0pi_ProtonBDTscorePhi.at("SelCC1pi")->Draw("colz");
    // PlottingHelper::SaveCanvas(pCanvas,"SidebandComparisons_BDTScoreProtonPhi_CC1pi");
    for (unsigned int i_hist=0; i_hist<TrueCC0pi_SelCC0pi_ProtonPIDMomentum.size(); i_hist++){
        TrueCC0pi_SelCC0pi_ProtonPIDMomentum.at(i_hist)->Draw("colz");
        PlottingHelper::SaveCanvas(pCanvas,"SidebandComparisons_ProtonPIDMomentum_"+std::to_string(i_hist)+"_CC0pi");
        TrueCC0pi_SelCC0pi_ProtonPIDCosTheta.at(i_hist)->Draw("colz");
        PlottingHelper::SaveCanvas(pCanvas,"SidebandComparisons_ProtonPIDCosTheta_"+std::to_string(i_hist)+"_CC0pi");
        TrueCC0pi_SelCC0pi_ProtonPIDPhi.at(i_hist)->Draw("colz");
        PlottingHelper::SaveCanvas(pCanvas,"SidebandComparisons_ProtonPIDPhi_"+std::to_string(i_hist)+"_CC0pi");

        TrueCC0pi_SelCC1pi_ProtonPIDMomentum.at(i_hist)->Draw("colz");
        PlottingHelper::SaveCanvas(pCanvas,"SidebandComparisons_ProtonPIDMomentum_"+std::to_string(i_hist)+"_CC1pi");
        TrueCC0pi_SelCC1pi_ProtonPIDCosTheta.at(i_hist)->Draw("colz");
        PlottingHelper::SaveCanvas(pCanvas,"SidebandComparisons_ProtonPIDCosTheta_"+std::to_string(i_hist)+"_CC1pi");
        TrueCC0pi_SelCC1pi_ProtonPIDPhi.at(i_hist)->Draw("colz");
        PlottingHelper::SaveCanvas(pCanvas,"SidebandComparisons_ProtonPIDPhi_"+std::to_string(i_hist)+"_CC1pi");
    }


}

} // namespace ubcc1pi_macros
