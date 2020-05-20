#include "ubcc1pi_standalone/Objects/FileReader.h"

#include "ubcc1pi_standalone/Helpers/PlottingHelper.h"
#include "ubcc1pi_standalone/Helpers/AnalysisHelper.h"
#include "ubcc1pi_standalone/Helpers/SelectionHelper.h"
#include "ubcc1pi_standalone/Helpers/BDTHelper.h"

#include <TH1F.h>
#include <TH2F.h>

using namespace ubcc1pi;

int PlotReconstructedVariables(const std::string &overlayFileName, const bool useAbsPdg = true, const float protonMomentumThreshold = 0.2f)
{
    // Set up the plots
    TH2F *pMuonCosThetaHist = new TH2F("muonCosTheta", "", 10, -1.f, 1.f, 10, -1.f, 1.f);
    TH2F *pPionCosThetaHist = new TH2F("pionCosTheta", "", 10, -1.f, 1.f, 10, -1.f, 1.f);
    TH2F *pPionMomentumHist = new TH2F("pionMomentum", "", 10, 0.f, 1.f, 10, 0.f, 1.f);
    TH2F *pMuonPionAngleHist = new TH2F("muonPionAngle", "", 10, 0.f, 3.14f, 10, 0.f, 3.14f);
    TH2F *pNProtonsHist = new TH2F("nProtons", "", 5, 0, 5, 5, 0, 5);
    
    TH1F *pMuonCosThetaResolutionHist = new TH1F("muonCosThetaResolution", "", 60, -0.5f, 0.5f);
    TH1F *pPionCosThetaResolutionHist = new TH1F("pionCosThetaResolution", "", 60, -0.5f, 0.5f);
    TH1F *pPionMomentumResolutionHist = new TH1F("pionMomentumResolution", "", 60, -1.f, 1.f);
    TH1F *pMuonPionAngleResolutionHist = new TH1F("muonPionAngleResolution", "", 60, -0.5f, 0.5f);
    TH1F *pNProtonsResolutionHist = new TH1F("nProtonsResolution", "", 8, -4, 4);

    const float bodgeFactor = 1.273f; // ATTN this factor is a normalisation added so we can compare the shape of the distributions, can't exist in the final result!
    
    // Get the selection
    auto selection = SelectionHelper::GetDefaultSelection();
    const auto allCuts = selection.GetCuts();

    FileReader reader(overlayFileName);
    auto pEvent = reader.GetBoundEventAddress();

    const auto nEvents = reader.GetNumberOfEvents();
    for (unsigned int i = 0; i < nEvents; ++i)
    {
        AnalysisHelper::PrintLoadingBar(i, nEvents);
        reader.LoadEvent(i);

        // Only use true signal events
        if (!AnalysisHelper::IsTrueCC1Pi(pEvent, useAbsPdg))
            continue;

        // Run the event selection and store which cuts are passed
        std::vector<std::string> cutsPassed;
        std::vector<int> assignedPdgCodes;
        const auto isSelected = selection.Execute(pEvent, cutsPassed, assignedPdgCodes);

        // Only use events which pass the generic selection
        const auto lastCut = "noShowers";
        if (std::find(cutsPassed.begin(), cutsPassed.end(), lastCut) == cutsPassed.end())
            continue;

        // Sanity check
        const auto recoParticles = pEvent->reco.particles;
        if (assignedPdgCodes.size() != recoParticles.size())
            throw std::logic_error("PlotReconstructedVariables - The output particle PDG codes is the wrong size");

        // Extract the selected particles
        auto muonIndex = std::numeric_limits<unsigned int>::max();
        auto pionIndex = std::numeric_limits<unsigned int>::max();
        unsigned int nMuons = 0u;
        unsigned int nPions = 0u;
        unsigned int nProtons = 0u;
        unsigned int nOther = 0u;

        for (unsigned int index = 0; index < assignedPdgCodes.size(); ++index)
        {
            const auto recoPdg = assignedPdgCodes.at(index);

            switch (recoPdg)
            {
                case 13:
                    muonIndex = index;
                    nMuons++;
                    break;
                case 211:
                    pionIndex = index;
                    nPions++;
                    break;
                case 2212:
                    nProtons++;
                    break;
                default:
                    nOther++;
                    break;
            }
        }

        // Make sure the reco PDGs make sense
        if (nMuons != 1)
            throw std::logic_error("PlotReconstructedVariables - Reconstructed " + std::to_string(nMuons) + " muons");

        if (nPions != 1)
            throw std::logic_error("PlotReconstructedVariables - Reconstructed " + std::to_string(nMuons) + " pions");

        // Get the muon and pion reconstructed particles
        const auto muon = recoParticles.at(muonIndex);
        const auto pion = recoParticles.at(pionIndex);

        // Find the reconstructed variables
        const auto muonDir = TVector3(muon.directionX(), muon.directionY(), muon.directionZ()).Unit();
        const auto muonCosTheta = muonDir.Z();

        const auto pionDir = TVector3(pion.directionX(), pion.directionY(), pion.directionZ()).Unit();
        const auto pionCosTheta = pionDir.Z();
        
        const auto pionMomentum = AnalysisHelper::GetPionMomentumFromRange(pion.range());

        const auto muonPionAngle = std::acos(muonDir.Dot(pionDir));



        // Get the truth particles
        const auto truthParticles = pEvent->truth.particles;
        auto trueMuonIndex = std::numeric_limits<unsigned int>::max();
        auto truePionIndex = std::numeric_limits<unsigned int>::max();
        unsigned int nTrueMuons = 0u;
        unsigned int nTruePions = 0u;
        unsigned int nTrueProtons = 0u;
        unsigned int nTrueOther = 0u;

        for (unsigned int index = 0; index < truthParticles.size(); ++index)
        {
            const auto truePdg = truthParticles.at(index).pdgCode();
            const auto truePdgAbs = useAbsPdg ? std::abs(truePdg) : truePdg;

            switch (truePdgAbs)
            {
                case 13:
                    trueMuonIndex = index;
                    nTrueMuons++;
                    break;
                case 211:
                    truePionIndex = index;
                    nTruePions++;
                    break;
                case 2212:
                    if (truthParticles.at(index).momentum() > protonMomentumThreshold)
                        nTrueProtons++;

                    break;
                default:
                    nTrueOther++;
                    break;
            }
        }

        // Make sure the true PDGs make sense
        if (nTrueMuons != 1)
            throw std::logic_error("PlotReconstructedVariables - In truth: " + std::to_string(nMuons) + " muons");

        if (nTruePions != 1)
            throw std::logic_error("PlotReconstructedVariables - In truth: " + std::to_string(nMuons) + " pions");

        // Get the muon and pion truth particles
        const auto trueMuon = truthParticles.at(trueMuonIndex);
        const auto truePion = truthParticles.at(truePionIndex);

        // Find the reconstructed variables
        const auto trueMuonDir = TVector3(trueMuon.momentumX(), trueMuon.momentumY(), trueMuon.momentumZ()).Unit();
        const auto trueMuonCosTheta = trueMuonDir.Z();

        const auto truePionDir = TVector3(truePion.momentumX(), truePion.momentumY(), truePion.momentumZ()).Unit();
        const auto truePionCosTheta = truePionDir.Z();
        const auto truePionMomentum = truePion.momentum();

        const auto trueMuonPionAngle = std::acos(trueMuonDir.Dot(truePionDir));

        // Fill the histograms
        pMuonCosThetaHist->Fill(trueMuonCosTheta, muonCosTheta);
        pPionCosThetaHist->Fill(truePionCosTheta, pionCosTheta);

        if (isSelected) // Only fill golden pion selected events
            pPionMomentumHist->Fill(truePionMomentum, pionMomentum);

        pMuonPionAngleHist->Fill(trueMuonPionAngle, muonPionAngle);
        pNProtonsHist->Fill(nTrueProtons, nProtons);

        if (trueMuonCosTheta > std::numeric_limits<float>::epsilon())
            pMuonCosThetaResolutionHist->Fill((muonCosTheta - trueMuonCosTheta));
        
        if (truePionCosTheta > std::numeric_limits<float>::epsilon())
            pPionCosThetaResolutionHist->Fill((pionCosTheta - truePionCosTheta));
        
        if (truePionMomentum > std::numeric_limits<float>::epsilon() && isSelected) // Only fill golden pion selected events
            pPionMomentumResolutionHist->Fill((pionMomentum - truePionMomentum));

        if (trueMuonPionAngle > std::numeric_limits<float>::epsilon())
            pMuonPionAngleResolutionHist->Fill((muonPionAngle - trueMuonPionAngle));

        pNProtonsResolutionHist->Fill(static_cast<int>(nProtons) - static_cast<int>(nTrueProtons));
    }

    auto pCanvas = PlottingHelper::GetCanvas();

    pMuonCosThetaHist->Draw("colz");
    PlottingHelper::SaveCanvas(pCanvas, "trueToReco_muonCosTheta");
    
    pPionCosThetaHist->Draw("colz");
    PlottingHelper::SaveCanvas(pCanvas, "trueToReco_pionCosTheta");
    
    pPionMomentumHist->Draw("colz");
    PlottingHelper::SaveCanvas(pCanvas, "trueToReco_pionMomentum");
    
    pMuonPionAngleHist->Draw("colz");
    PlottingHelper::SaveCanvas(pCanvas, "trueToReco_muonPionAngle");

    pNProtonsHist->Draw("colz");
    PlottingHelper::SaveCanvas(pCanvas, "trueToReco_nProtons");
    
    pMuonCosThetaResolutionHist->Draw("colz");
    PlottingHelper::SaveCanvas(pCanvas, "trueToReco_resolution_muonCosTheta");
    
    pPionCosThetaResolutionHist->Draw("colz");
    PlottingHelper::SaveCanvas(pCanvas, "trueToReco_resolution_pionCosTheta");
    
    pPionMomentumResolutionHist->Draw("colz");
    PlottingHelper::SaveCanvas(pCanvas, "trueToReco_resolution_pionMomentum");
    
    pMuonPionAngleResolutionHist->Draw("colz");
    PlottingHelper::SaveCanvas(pCanvas, "trueToReco_resolution_muonPionAngle");

    pNProtonsResolutionHist->Draw("colz");
    PlottingHelper::SaveCanvas(pCanvas, "trueToReco_resolution_nProtons");

    return 0;
}
