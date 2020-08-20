/**
 *  @file  ubcc1pi_standalone/Macros/MultiPlanePIDDemo.cxx
 *
 *  @brief The implementation file of the MultiPlanePIDDemo macro
 */

#include "ubcc1pi_standalone/Macros/Macros.h"

#include "ubcc1pi_standalone/Objects/FileReader.h"

#include "ubcc1pi_standalone/Helpers/PlottingHelper.h"
#include "ubcc1pi_standalone/Helpers/AnalysisHelper.h"
#include "ubcc1pi_standalone/Helpers/FormattingHelper.h"
#include "ubcc1pi_standalone/Helpers/NormalisationHelper.h"

#include <TH2F.h>
#include <TLine.h>
#include <TArrow.h>

using namespace ubcc1pi;

namespace ubcc1pi_macros
{

void MultiPlanePIDDemo(const Config &config)
{
    //
    // Setup the input files
    // 
    std::vector< std::tuple<AnalysisHelper::SampleType, std::string, float> > inputData;
    
    inputData.emplace_back(AnalysisHelper::Overlay, config.files.overlaysFileName, NormalisationHelper::GetOverlaysNormalisation(config)); 
    inputData.emplace_back(AnalysisHelper::Dirt,    config.files.dirtFileName,     NormalisationHelper::GetDirtNormalisation(config)); 
    inputData.emplace_back(AnalysisHelper::DataEXT, config.files.dataEXTFileName,  NormalisationHelper::GetDataEXTNormalisation(config)); 
    
    const auto angleThreshold = std::asin(std::pow(config.multiPlanePIDDemo.sin2AngleThreshold, 0.5f));
    std::cout << "Angle threshold = " << angleThreshold << std::endl;

    // Choose the number of bins such that the angle threhsold lies on a bin edge
    const auto pi = 2 * std::acos(0);
    const auto nBinsX = 40u;
    const auto nBinsY = 40u;

    const auto minVal = -9.f;
    const auto maxVal = 7.f;
    TH2F *hU = new TH2F("hU", "", nBinsX, -pi, pi, nBinsY, minVal, maxVal);
    TH2F *hV = new TH2F("hV", "", nBinsX, -pi, pi, nBinsY, minVal, maxVal);
    TH2F *hW = new TH2F("hW", "", nBinsX, -pi, pi, nBinsY, minVal, maxVal);
    TH2F *hUV = new TH2F("hUV", "", nBinsX, -pi, pi, nBinsY, minVal, maxVal);
    TH2F *hUVW = new TH2F("hUVW", "", nBinsX, -pi, pi, nBinsY, minVal, maxVal);

    float nTotal = 0.f;
    float nU = 0.f;
    float nV = 0.f;
    float nW = 0.f;
    float nUAngleCut = 0.f;
    float nVAngleCut = 0.f;
    float nWAngleCut = 0.f;
    float nUV = 0.f;
    float nUVW = 0.f;
    
    // Define the special angles
    const auto uAngle = pi / 3.f;
    const auto vAngle = -pi / 3.f;
    const auto wAngle = 0.f;

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

            const auto weight = AnalysisHelper::GetNominalEventWeight(pEvent) * normalisation;
            const auto recoParticles = pEvent->reco.particles;
            const auto truthParticles = pEvent->truth.particles;

            for (const auto &particle : recoParticles)
            {
                if (!AnalysisHelper::HasTrackFit(particle))
                    continue;
                
                // Get the plot style 
                const auto particleStyle = PlottingHelper::GetPlotStyle(particle, sampleType, truthParticles, false, config.global.useAbsPdg);

                // Fold the angle so that antiparallel angles have the same value, and the range is 0-pi
                //const auto yzAngle = particle.yzAngle() - std::floor(particle.yzAngle() / pi) * pi; 
                const auto yzAngle = particle.yzAngle();

                float ratioU = std::numeric_limits<float>::max();
                float ratioV = std::numeric_limits<float>::max();
                float ratioW = std::numeric_limits<float>::max();
                float ratioUVW = std::numeric_limits<float>::max();

                const bool hasU = AnalysisHelper::GetLogLikelihoodRatio(particle.likelihoodForwardProtonU, particle.likelihoodMIPU, ratioU);
                const bool hasV = AnalysisHelper::GetLogLikelihoodRatio(particle.likelihoodForwardProtonV, particle.likelihoodMIPV, ratioV);
                const bool hasW = AnalysisHelper::GetLogLikelihoodRatio(particle.likelihoodForwardProtonW, particle.likelihoodMIPW, ratioW);
                const bool hasUVW = AnalysisHelper::GetLogLikelihoodRatio(particle.likelihoodForwardProton, particle.likelihoodMIP, ratioUVW);

                const bool isTrackAlongUWire = (std::pow(std::sin(yzAngle - uAngle), 2) < config.multiPlanePIDDemo.sin2AngleThreshold);
                const bool isTrackAlongVWire = (std::pow(std::sin(yzAngle - vAngle), 2) < config.multiPlanePIDDemo.sin2AngleThreshold);
                const bool isTrackAlongWWire = (std::pow(std::sin(yzAngle - wAngle), 2) < config.multiPlanePIDDemo.sin2AngleThreshold);

                if (hasU)
                {
                    hU->Fill(yzAngle, ratioU, weight);

                    nU += weight;

                    if (!isTrackAlongUWire)
                        nUAngleCut += weight;
                }

                if (hasV)
                {
                    hV->Fill(yzAngle, ratioV, weight);

                    nV += weight;

                    if (!isTrackAlongVWire)
                        nVAngleCut += weight;
                }

                if (hasW)
                {
                    hW->Fill(yzAngle, ratioW, weight);

                    nW += weight;

                    if (!isTrackAlongWWire)
                        nWAngleCut += weight;
                }

                if (hasUVW)
                {
                    hUVW->Fill(yzAngle, ratioUVW, weight);
                    nUVW += weight;
                }

                // Manually make the UV combination (this has already been done in LArSoft, and we just use the final UVW value)
                const bool useU = !isTrackAlongUWire && hasU;
                const bool useV = !isTrackAlongVWire && hasV;

                if (useU || useV)
                {
                    const auto weightU = useU ? static_cast<float>(particle.nHitsU()) : 0.f;
                    const auto weightV = useV ? static_cast<float>(particle.nHitsV()) : 0.f;
                    const auto weightUV = weightU + weightV;

                    const auto ratioUWeighted = useU ? ((ratioU * weightU) / weightUV) : 0.f;
                    const auto ratioVWeighted = useV ? ((ratioV * weightV) / weightUV) : 0.f;

                    const auto ratioUV = ratioUWeighted + ratioVWeighted;

                    hUV->Fill(yzAngle, ratioUV, weight);
                    nUV += weight;
                }

                nTotal += weight;
            }
        }
    }

    // Normalize the histograms
    const float zScaleFactor = 10000.f;
    hU->Scale(zScaleFactor / hU->GetEntries());
    hV->Scale(zScaleFactor / hV->GetEntries());
    hW->Scale(zScaleFactor / hW->GetEntries());
    hUV->Scale(zScaleFactor / hUV->GetEntries());
    hUVW->Scale(zScaleFactor / hUVW->GetEntries());

    auto pCanvas = PlottingHelper::GetCanvas();

    // U lines
    const auto yArrow = 0.f;
    const auto arrowSize = 0.02f;
    TLine *lU1 = new TLine(uAngle, minVal, uAngle, maxVal);
    TLine *lU1Plus = new TLine(uAngle + angleThreshold, minVal, uAngle + angleThreshold, maxVal);
    TLine *lU1Minus = new TLine(uAngle - angleThreshold, minVal, uAngle - angleThreshold, maxVal);
    TLine *lU2 = new TLine(uAngle - pi, minVal, uAngle - pi, maxVal);
    TLine *lU2Plus = new TLine(uAngle - pi + angleThreshold, minVal, uAngle - pi + angleThreshold, maxVal);
    TLine *lU2Minus = new TLine(uAngle - pi - angleThreshold, minVal, uAngle - pi - angleThreshold, maxVal);

    TArrow *aU1 = new TArrow(-pi, yArrow, uAngle - pi - angleThreshold, yArrow, arrowSize, ">");
    TArrow *aU2 = new TArrow(uAngle - pi + angleThreshold, yArrow, uAngle - angleThreshold, yArrow, arrowSize, "<>");
    TArrow *aU3 = new TArrow(uAngle + angleThreshold, yArrow, pi, yArrow, arrowSize, "<");

    lU1->SetLineColor(kWhite);
    lU1Plus->SetLineColor(kWhite);
    lU1Minus->SetLineColor(kWhite);
    lU2->SetLineColor(kWhite);
    lU2Plus->SetLineColor(kWhite);
    lU2Minus->SetLineColor(kWhite);
    lU1->SetLineWidth(2);
    lU1Plus->SetLineWidth(2);
    lU1Minus->SetLineWidth(2);
    lU1Plus->SetLineStyle(2);
    lU1Minus->SetLineStyle(2);
    lU2->SetLineWidth(2);
    lU2Plus->SetLineWidth(2);
    lU2Minus->SetLineWidth(2);
    lU2Plus->SetLineStyle(2);
    lU2Minus->SetLineStyle(2);

    aU1->SetLineColor(kWhite);
    aU1->SetLineWidth(2);
    aU2->SetLineColor(kWhite);
    aU2->SetLineWidth(2);
    aU3->SetLineColor(kWhite);
    aU3->SetLineWidth(2);
    
    // V lines
    TLine *lV1 = new TLine(vAngle + pi, minVal, vAngle + pi, maxVal);
    TLine *lV1Plus = new TLine(vAngle + pi + angleThreshold, minVal, vAngle + pi + angleThreshold, maxVal);
    TLine *lV1Minus = new TLine(vAngle + pi - angleThreshold, minVal, vAngle + pi - angleThreshold, maxVal);
    TLine *lV2 = new TLine(vAngle, minVal, vAngle, maxVal);
    TLine *lV2Plus = new TLine(vAngle + angleThreshold, minVal, vAngle + angleThreshold, maxVal);
    TLine *lV2Minus = new TLine(vAngle - angleThreshold, minVal, vAngle - angleThreshold, maxVal);
    
    TArrow *aV1 = new TArrow(-pi, yArrow, vAngle - angleThreshold, yArrow, arrowSize, ">");
    TArrow *aV2 = new TArrow(vAngle + angleThreshold, yArrow, vAngle + pi - angleThreshold, yArrow, arrowSize, "<>");
    TArrow *aV3 = new TArrow(vAngle + pi + angleThreshold, yArrow, pi, yArrow, arrowSize, "<");

    lV1->SetLineColor(kWhite);
    lV1Plus->SetLineColor(kWhite);
    lV1Minus->SetLineColor(kWhite);
    lV2->SetLineColor(kWhite);
    lV2Plus->SetLineColor(kWhite);
    lV2Minus->SetLineColor(kWhite);
    lV1->SetLineWidth(2);
    lV1Plus->SetLineWidth(2);
    lV1Minus->SetLineWidth(2);
    lV1Plus->SetLineStyle(2);
    lV1Minus->SetLineStyle(2);
    lV2->SetLineWidth(2);
    lV2Plus->SetLineWidth(2);
    lV2Minus->SetLineWidth(2);
    lV2Plus->SetLineStyle(2);
    lV2Minus->SetLineStyle(2);
    
    aV1->SetLineColor(kWhite);
    aV1->SetLineWidth(2);
    aV2->SetLineColor(kWhite);
    aV2->SetLineWidth(2);
    aV3->SetLineColor(kWhite);
    aV3->SetLineWidth(2);
    
    // W lines
    TLine *lW1 = new TLine(wAngle, minVal, wAngle, maxVal);
    TLine *lW1Plus = new TLine(wAngle + angleThreshold, minVal, wAngle + angleThreshold, maxVal);
    TLine *lW1Minus = new TLine(wAngle - angleThreshold, minVal, wAngle - angleThreshold, maxVal);
    TLine *lW2Plus = new TLine(wAngle - pi + angleThreshold, minVal, wAngle - pi + angleThreshold, maxVal);
    TLine *lW2Minus = new TLine(wAngle + pi - angleThreshold, minVal, wAngle + pi - angleThreshold, maxVal);
    
    TArrow *aW1 = new TArrow(-pi + angleThreshold, yArrow, wAngle - angleThreshold, yArrow, arrowSize, "<>");
    TArrow *aW2 = new TArrow(wAngle + angleThreshold, yArrow, pi - angleThreshold, yArrow, arrowSize, "<>");
    TArrow *aUVW = new TArrow(-pi, yArrow, pi, yArrow, arrowSize, "<>");

    lW1->SetLineColor(kWhite);
    lW1Plus->SetLineColor(kWhite);
    lW1Minus->SetLineColor(kWhite);
    lW2Plus->SetLineColor(kWhite);
    lW2Minus->SetLineColor(kWhite);
    lW1->SetLineWidth(2);
    lW1Plus->SetLineWidth(2);
    lW1Minus->SetLineWidth(2);
    lW1Plus->SetLineStyle(2);
    lW1Minus->SetLineStyle(2);
    lW2Plus->SetLineWidth(2);
    lW2Minus->SetLineWidth(2);
    lW2Plus->SetLineStyle(2);
    lW2Minus->SetLineStyle(2);
    
    aW1->SetLineColor(kWhite);
    aW1->SetLineWidth(2);
    aW2->SetLineColor(kWhite);
    aW2->SetLineWidth(2);
    aUVW->SetLineColor(kWhite);
    aUVW->SetLineWidth(2);


    // Draw the histograms

    hU->Draw("colz");
//    lU1->Draw();
    lU1Plus->Draw();
    lU1Minus->Draw();
//    lU2->Draw();
    lU2Plus->Draw();
    lU2Minus->Draw();
    aU1->Draw();
    aU2->Draw();
    aU3->Draw();
    PlottingHelper::SaveCanvas(pCanvas, "multiPlanePIDDemo_logLikelihood_p_MIP_U");
    
    hV->Draw("colz");
//    lV1->Draw();
    lV1Plus->Draw();
    lV1Minus->Draw();
//    lV2->Draw();
    lV2Plus->Draw();
    lV2Minus->Draw();
    aV1->Draw();
    aV2->Draw();
    aV3->Draw();
    PlottingHelper::SaveCanvas(pCanvas, "multiPlanePIDDemo_logLikelihood_p_MIP_V");
    
    hW->Draw("colz");
//    lW1->Draw();
    lW1Plus->Draw();
    lW1Minus->Draw();
    lW2Plus->Draw();
    lW2Minus->Draw();
    aW1->Draw();
    aW2->Draw();
    PlottingHelper::SaveCanvas(pCanvas, "multiPlanePIDDemo_logLikelihood_p_MIP_W");
    
    hUV->Draw("colz");
    aUVW->Draw();
    PlottingHelper::SaveCanvas(pCanvas, "multiPlanePIDDemo_logLikelihood_p_MIP_UV");
    
    hUVW->Draw("colz");
//    lW1->Draw();
//    lW1Plus->Draw();
//    lW1Minus->Draw();
//    lW2Plus->Draw();
//    lW2Minus->Draw();
    aUVW->Draw();
    PlottingHelper::SaveCanvas(pCanvas, "multiPlanePIDDemo_logLikelihood_p_MIP_UVW");

    FormattingHelper::Table table({"Plane", "Acceptance"});

    table.AddEmptyRow();
    table.SetEntry("Plane", "U");
    table.SetEntry("Acceptance", nU / nTotal);
    
    table.AddEmptyRow();
    table.SetEntry("Plane", "U (good angle)");
    table.SetEntry("Acceptance", nUAngleCut / nTotal);
    
    table.AddEmptyRow();
    table.SetEntry("Plane", "V");
    table.SetEntry("Acceptance", nV / nTotal);
    
    table.AddEmptyRow();
    table.SetEntry("Plane", "V (good angle)");
    table.SetEntry("Acceptance", nVAngleCut / nTotal);
    
    table.AddEmptyRow();
    table.SetEntry("Plane", "W");
    table.SetEntry("Acceptance", nW / nTotal);
    
    table.AddEmptyRow();
    table.SetEntry("Plane", "W (good angle)");
    table.SetEntry("Acceptance", nWAngleCut / nTotal);
    
    table.AddEmptyRow();
    table.SetEntry("Plane", "UV");
    table.SetEntry("Acceptance", nUV / nTotal);
    
    table.AddEmptyRow();
    table.SetEntry("Plane", "UVW");
    table.SetEntry("Acceptance", nUVW / nTotal);

    table.WriteToFile("multiPlanePIDDemo_acceptances.md");
}

} // namespace ubcc1pi_macros
