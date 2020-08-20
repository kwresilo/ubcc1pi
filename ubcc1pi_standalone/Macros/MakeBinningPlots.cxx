/**
 *  @file  ubcc1pi_standalone/Macros/MakeBinningPlots.cxx
 *
 *  @brief The implementation file of the MakeBinningPlots macro
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

void MakeBinningPlots(const Config &config)
{
    //
    // Setup the input files
    // 
    std::vector< std::tuple<AnalysisHelper::SampleType, std::string, float> > inputData;
    
    inputData.emplace_back(AnalysisHelper::Overlay, config.files.overlaysFileName, NormalisationHelper::GetOverlaysNormalisation(config)); 
    inputData.emplace_back(AnalysisHelper::Dirt,    config.files.dirtFileName,     NormalisationHelper::GetDirtNormalisation(config)); 
    inputData.emplace_back(AnalysisHelper::DataEXT, config.files.dataEXTFileName,  NormalisationHelper::GetDataEXTNormalisation(config)); 
    //inputData.emplace_back(AnalysisHelper::DataBNB, config.files.dataBNBFileName,  1.f);

    // Get the selection
    auto selection = SelectionHelper::GetDefaultSelection();

    // Get the fine bin edges
    auto getFineBinEdges = [](const unsigned int nBinsTarget, const float min, const float max, const std::vector<float> &coarseBinEdges) -> std::vector<float> {

        // If required extend the bin edges to the min and max
        std::vector<float> extendedCoarseBinEdges;

        if (coarseBinEdges.front() - min > std::numeric_limits<float>::epsilon())
            extendedCoarseBinEdges.push_back(min);

        extendedCoarseBinEdges.insert(extendedCoarseBinEdges.end(), coarseBinEdges.begin(), coarseBinEdges.end());

        if (max - coarseBinEdges.back() > std::numeric_limits<float>::epsilon())
            extendedCoarseBinEdges.push_back(max);
        
        // Now break up the coarse bins into smaller bins
        // Use the nBinsTarget to determine the ideal width of these smaller bins
        const auto targetBinWidth = (max - min) / static_cast<float>(nBinsTarget);

        std::vector<float> fineBinEdges(1, extendedCoarseBinEdges.front());
        for (unsigned int i = 1; i < extendedCoarseBinEdges.size(); ++i)
        {
            const auto lowerEdge = extendedCoarseBinEdges.at(i - 1);
            const auto upperEdge = extendedCoarseBinEdges.at(i);
            const auto binWidth = upperEdge - lowerEdge;

            const auto nFineBins = std::max(1u, static_cast<unsigned int>(std::round(binWidth / targetBinWidth)));
            const auto fineBinWidth = binWidth / static_cast<float>(nFineBins);
        
            // Add the upper edges of the fine bins
            for (unsigned int j = 1; j <= nFineBins; ++j)
            { 
                const auto fineUpperEdge = lowerEdge + j * fineBinWidth;
                fineBinEdges.push_back(fineUpperEdge);
            }
        }
        
        return fineBinEdges;
    };
    
    const auto edges_muonCosTheta = getFineBinEdges(50u, -1.f, 1.f, config.global.muonCosTheta.binEdges);
    const auto nFineBins_muonCosTheta = edges_muonCosTheta.size() - 1;
    
    const auto edges_muonPhi = getFineBinEdges(50u, -3.142f, 3.142f, config.global.muonPhi.binEdges);
    const auto nFineBins_muonPhi = edges_muonPhi.size() - 1;
    
    const auto edges_muonMomentum = getFineBinEdges(50u, 0.f, 2.f, config.global.muonMomentum.binEdges);
    const auto nFineBins_muonMomentum = edges_muonMomentum.size() - 1;
    
    const auto edges_pionCosTheta = getFineBinEdges(50u, -1.f, 1.f, config.global.pionCosTheta.binEdges);
    const auto nFineBins_pionCosTheta = edges_pionCosTheta.size() - 1;
    
    const auto edges_pionPhi = getFineBinEdges(50u, -3.142f, 3.142f, config.global.pionPhi.binEdges);
    const auto nFineBins_pionPhi = edges_pionPhi.size() - 1;
    
    const auto edges_pionMomentum = getFineBinEdges(50u, 0.f, 0.6f, config.global.pionMomentum.binEdges);
    const auto nFineBins_pionMomentum = edges_pionMomentum.size() - 1;
    
    const auto edges_muonPionAngle = getFineBinEdges(50u, 0.f, 3.142f, config.global.muonPionAngle.binEdges);
    const auto nFineBins_muonPionAngle = edges_muonPionAngle.size() - 1;

    // Setup the efficiency plots
    const auto cutLabels = std::vector<std::string>({"generic", "golden"});
    PlottingHelper::EfficiencyPlot eff_muonCosTheta("True muon cos(theta)", edges_muonCosTheta, cutLabels);
    PlottingHelper::EfficiencyPlot eff_muonPhi("True muon phi / rad", edges_muonPhi, cutLabels);
    PlottingHelper::EfficiencyPlot eff_muonMomentum("True muon momentum / GeV", edges_muonMomentum, cutLabels);
    PlottingHelper::EfficiencyPlot eff_pionCosTheta("True pion cos(theta)", edges_pionCosTheta, cutLabels);
    PlottingHelper::EfficiencyPlot eff_pionPhi("True pion phi / rad", edges_pionPhi, cutLabels);
    PlottingHelper::EfficiencyPlot eff_pionMomentum("True pion momentum / GeV", edges_pionMomentum, cutLabels);
    PlottingHelper::EfficiencyPlot eff_muonPionAngle("True muon-pion opening angle / rad", edges_muonPionAngle, cutLabels);

    // Setup the multi-plots
    const std::string yLabel = "Number of events";

    PlottingHelper::MultiPlot multi_muonCosTheta_generic("Reco muon cos(theta)", yLabel, edges_muonCosTheta);
    PlottingHelper::MultiPlot multi_muonCosTheta_golden("Reco muon cos(theta)", yLabel, edges_muonCosTheta);
    
    PlottingHelper::MultiPlot multi_muonPhi_generic("Reco muon phi / rad", yLabel, edges_muonPhi);
    PlottingHelper::MultiPlot multi_muonPhi_golden("Reco muon phi / rad", yLabel, edges_muonPhi);
    
    PlottingHelper::MultiPlot multi_muonMomentum_generic("Reco muon momentum / GeV", yLabel, edges_muonMomentum);
    PlottingHelper::MultiPlot multi_muonMomentum_golden("Reco muon momentum / GeV", yLabel, edges_muonMomentum);
    
    PlottingHelper::MultiPlot multi_pionCosTheta_generic("Reco pion cos(theta)", yLabel, edges_pionCosTheta);
    PlottingHelper::MultiPlot multi_pionCosTheta_golden("Reco pion cos(theta)", yLabel, edges_pionCosTheta);
    
    PlottingHelper::MultiPlot multi_pionPhi_generic("Reco pion phi / rad", yLabel, edges_pionPhi);
    PlottingHelper::MultiPlot multi_pionPhi_golden("Reco pion phi / rad", yLabel, edges_pionPhi);
    
    PlottingHelper::MultiPlot multi_pionMomentum_generic("Reco pion momentum / GeV", yLabel, edges_pionMomentum);
    PlottingHelper::MultiPlot multi_pionMomentum_golden("Reco pion momentum / GeV", yLabel, edges_pionMomentum);
    
    PlottingHelper::MultiPlot multi_muonPionAngle_generic("Reco muon-pion opening angle / rad", yLabel, edges_muonPionAngle);
    PlottingHelper::MultiPlot multi_muonPionAngle_golden("Reco muon-pion opening angle / rad", yLabel, edges_muonPionAngle);
    
    // Add the lines at the bin edges
    for (const auto &edge : config.global.muonCosTheta.binEdges)
    {
        multi_muonCosTheta_generic.AddCutLine(edge);
        multi_muonCosTheta_golden.AddCutLine(edge);
        eff_muonCosTheta.AddCutLine(edge);
    }
    
    for (const auto &edge : config.global.muonPhi.binEdges)
    {
        multi_muonPhi_generic.AddCutLine(edge);
        multi_muonPhi_golden.AddCutLine(edge);
        eff_muonPhi.AddCutLine(edge);
    }
    
    for (const auto &edge : config.global.muonMomentum.binEdges)
    {
        multi_muonMomentum_generic.AddCutLine(edge);
        multi_muonMomentum_golden.AddCutLine(edge);
        eff_muonMomentum.AddCutLine(edge);
    }
    
    for (const auto &edge : config.global.pionCosTheta.binEdges)
    {
        multi_pionCosTheta_generic.AddCutLine(edge);
        multi_pionCosTheta_golden.AddCutLine(edge);
        eff_pionCosTheta.AddCutLine(edge);
    }
    
    for (const auto &edge : config.global.pionPhi.binEdges)
    {
        multi_pionPhi_generic.AddCutLine(edge);
        multi_pionPhi_golden.AddCutLine(edge);
        eff_pionPhi.AddCutLine(edge);
    }
    
    for (const auto &edge : config.global.pionMomentum.binEdges)
    {
        multi_pionMomentum_generic.AddCutLine(edge);
        multi_pionMomentum_golden.AddCutLine(edge);
        eff_pionMomentum.AddCutLine(edge);
    }

    for (const auto &edge : config.global.muonPionAngle.binEdges)
    {
        multi_muonPionAngle_generic.AddCutLine(edge);
        multi_muonPionAngle_golden.AddCutLine(edge);
        eff_muonPionAngle.AddCutLine(edge);
    }

    // Setup the resolution plots
    TH2F *hRes_muonCosTheta_generic = new TH2F("hRes_muonCosTheta_generic", "", nFineBins_muonCosTheta, edges_muonCosTheta.data(), nFineBins_muonCosTheta, edges_muonCosTheta.data());
    TH2F *hRes_muonCosTheta_golden = new TH2F("hRes_muonCosTheta_golden", "", nFineBins_muonCosTheta, edges_muonCosTheta.data(), nFineBins_muonCosTheta, edges_muonCosTheta.data());
    
    TH2F *hRes_muonPhi_generic = new TH2F("hRes_muonPhi_generic", "", nFineBins_muonPhi, edges_muonPhi.data(), nFineBins_muonPhi, edges_muonPhi.data());
    TH2F *hRes_muonPhi_golden = new TH2F("hRes_muonPhi_golden", "", nFineBins_muonPhi, edges_muonPhi.data(), nFineBins_muonPhi, edges_muonPhi.data());
    
    TH2F *hRes_muonMomentum_generic = new TH2F("hRes_muonMomentum_generic", "", nFineBins_muonMomentum, edges_muonMomentum.data(), nFineBins_muonMomentum, edges_muonMomentum.data());
    TH2F *hRes_muonMomentum_golden = new TH2F("hRes_muonMomentum_golden", "", nFineBins_muonMomentum, edges_muonMomentum.data(), nFineBins_muonMomentum, edges_muonMomentum.data());
    
    TH2F *hRes_pionCosTheta_generic = new TH2F("hRes_pionCosTheta_generic", "", nFineBins_pionCosTheta, edges_pionCosTheta.data(), nFineBins_pionCosTheta, edges_pionCosTheta.data());
    TH2F *hRes_pionCosTheta_golden = new TH2F("hRes_pionCosTheta_golden", "", nFineBins_pionCosTheta, edges_pionCosTheta.data(), nFineBins_pionCosTheta, edges_pionCosTheta.data());
    
    TH2F *hRes_pionPhi_generic = new TH2F("hRes_pionPhi_generic", "", nFineBins_pionPhi, edges_pionPhi.data(), nFineBins_pionPhi, edges_pionPhi.data());
    TH2F *hRes_pionPhi_golden = new TH2F("hRes_pionPhi_golden", "", nFineBins_pionPhi, edges_pionPhi.data(), nFineBins_pionPhi, edges_pionPhi.data());
    
    TH2F *hRes_pionMomentum_generic = new TH2F("hRes_pionMomentum_generic", "", nFineBins_pionMomentum, edges_pionMomentum.data(), nFineBins_pionMomentum, edges_pionMomentum.data());
    TH2F *hRes_pionMomentum_golden = new TH2F("hRes_pionMomentum_golden", "", nFineBins_pionMomentum, edges_pionMomentum.data(), nFineBins_pionMomentum, edges_pionMomentum.data());
    
    TH2F *hRes_muonPionAngle_generic = new TH2F("hRes_muonPionAngle_generic", "", nFineBins_muonPionAngle, edges_muonPionAngle.data(), nFineBins_muonPionAngle, edges_muonPionAngle.data());
    TH2F *hRes_muonPionAngle_golden = new TH2F("hRes_muonPionAngle_golden", "", nFineBins_muonPionAngle, edges_muonPionAngle.data(), nFineBins_muonPionAngle, edges_muonPionAngle.data());
    
    // Loop over the events
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

            const auto isSignal = isOverlay && AnalysisHelper::IsTrueCC1Pi(pEvent, config.global.useAbsPdg);
            const auto weight = AnalysisHelper::GetNominalEventWeight(pEvent);

            // Run the event selection and store which cuts are passed
            std::vector<std::string> cutsPassed;
            std::vector<int> assignedPdgCodes;
            auto passedGoldenSelection = selection.Execute(pEvent, cutsPassed, assignedPdgCodes);
            auto passedGenericSelection = (std::find(cutsPassed.begin(), cutsPassed.end(), config.global.lastCutGeneric) != cutsPassed.end());

            // Get the reco analysis data (if available)
            auto recoData = AnalysisHelper::GetDummyAnalysisData();
            if (passedGenericSelection)
            {
                recoData = AnalysisHelper::GetRecoAnalysisData(pEvent->reco, assignedPdgCodes, passedGoldenSelection);
            }

            // Fill the efficiency plots
            if (isSignal)
            {
                const auto truthData = AnalysisHelper::GetTruthAnalysisData(pEvent->truth, config.global.useAbsPdg, config.global.protonMomentumThreshold);

                eff_muonCosTheta.AddEvent(truthData.muonCosTheta, weight, "generic", passedGenericSelection);
                eff_muonCosTheta.AddEvent(truthData.muonCosTheta, weight, "golden", passedGoldenSelection);
                
                eff_muonPhi.AddEvent(truthData.muonPhi, weight, "generic", passedGenericSelection);
                eff_muonPhi.AddEvent(truthData.muonPhi, weight, "golden", passedGoldenSelection);
                
                eff_muonMomentum.AddEvent(truthData.muonMomentum, weight, "generic", passedGenericSelection);
                eff_muonMomentum.AddEvent(truthData.muonMomentum, weight, "golden", passedGoldenSelection);
                
                eff_pionCosTheta.AddEvent(truthData.pionCosTheta, weight, "generic", passedGenericSelection);
                eff_pionCosTheta.AddEvent(truthData.pionCosTheta, weight, "golden", passedGoldenSelection);
                
                eff_pionPhi.AddEvent(truthData.pionPhi, weight, "generic", passedGenericSelection);
                eff_pionPhi.AddEvent(truthData.pionPhi, weight, "golden", passedGoldenSelection);
                
                eff_pionMomentum.AddEvent(truthData.pionMomentum, weight, "generic", passedGenericSelection);
                eff_pionMomentum.AddEvent(truthData.pionMomentum, weight, "golden", passedGoldenSelection);
                
                eff_muonPionAngle.AddEvent(truthData.muonPionAngle, weight, "generic", passedGenericSelection);
                eff_muonPionAngle.AddEvent(truthData.muonPionAngle, weight, "golden", passedGoldenSelection);

                // Fill the resolution plots
                if (passedGenericSelection)
                {
                    hRes_muonCosTheta_generic->Fill(truthData.muonCosTheta, recoData.muonCosTheta, weight);
                    hRes_muonPhi_generic->Fill(truthData.muonPhi, recoData.muonPhi, weight);
                    hRes_muonMomentum_generic->Fill(truthData.muonMomentum, recoData.muonMomentum, weight);
                    hRes_pionCosTheta_generic->Fill(truthData.pionCosTheta, recoData.pionCosTheta, weight);
                    hRes_pionPhi_generic->Fill(truthData.pionPhi, recoData.pionPhi, weight);
                    hRes_pionMomentum_generic->Fill(truthData.pionMomentum, recoData.pionMomentum, weight);
                    hRes_muonPionAngle_generic->Fill(truthData.muonPionAngle, recoData.muonPionAngle, weight);
                }

                if (passedGoldenSelection)
                {
                    hRes_muonCosTheta_golden->Fill(truthData.muonCosTheta, recoData.muonCosTheta, weight);
                    hRes_muonPhi_golden->Fill(truthData.muonPhi, recoData.muonPhi, weight);
                    hRes_muonMomentum_golden->Fill(truthData.muonMomentum, recoData.muonMomentum, weight);
                    hRes_pionCosTheta_golden->Fill(truthData.pionCosTheta, recoData.pionCosTheta, weight);
                    hRes_pionPhi_golden->Fill(truthData.pionPhi, recoData.pionPhi, weight);
                    hRes_pionMomentum_golden->Fill(truthData.pionMomentum, recoData.pionMomentum, weight);
                    hRes_muonPionAngle_golden->Fill(truthData.muonPionAngle, recoData.muonPionAngle, weight);
                }
            }

            // Fill the multi plots
            const auto style = PlottingHelper::GetPlotStyle(sampleType, pEvent, config.global.useAbsPdg);
            const auto factor = weight*normalisation;

            if (passedGenericSelection)
            {
                multi_muonCosTheta_generic.Fill(recoData.muonCosTheta, style, factor);
                multi_muonPhi_generic.Fill(recoData.muonPhi, style, factor);
                multi_muonMomentum_generic.Fill(recoData.muonMomentum, style, factor);
                multi_pionCosTheta_generic.Fill(recoData.pionCosTheta, style, factor);
                multi_pionPhi_generic.Fill(recoData.pionPhi, style, factor);
                multi_pionMomentum_generic.Fill(recoData.pionMomentum, style, factor);
                multi_muonPionAngle_generic.Fill(recoData.muonPionAngle, style, factor);
            }
            
            if (passedGoldenSelection)
            {
                multi_muonCosTheta_golden.Fill(recoData.muonCosTheta, style, factor);
                multi_muonPhi_golden.Fill(recoData.muonPhi, style, factor);
                multi_muonMomentum_golden.Fill(recoData.muonMomentum, style, factor);
                multi_pionCosTheta_golden.Fill(recoData.pionCosTheta, style, factor);
                multi_pionPhi_golden.Fill(recoData.pionPhi, style, factor);
                multi_pionMomentum_golden.Fill(recoData.pionMomentum, style, factor);
                multi_muonPionAngle_golden.Fill(recoData.muonPionAngle, style, factor);
            }
        }
    }

    // Scale the 2D histograms by bin width
    for (auto &pHist : {
        hRes_muonCosTheta_generic,
        hRes_muonCosTheta_golden,
        
        hRes_muonPhi_generic,
        hRes_muonPhi_golden,
        
        hRes_muonMomentum_generic,
        hRes_muonMomentum_golden,
        
        hRes_pionCosTheta_generic,
        hRes_pionCosTheta_golden,
        
        hRes_pionPhi_generic,
        hRes_pionPhi_golden,
        
        hRes_pionMomentum_generic,
        hRes_pionMomentum_golden,

        hRes_muonPionAngle_generic,
        hRes_muonPionAngle_golden
    })
    {
        pHist->Scale(1.f, "width");
    }

    // Save the efficiency plots
    const auto cutStyles = std::vector<PlottingHelper::PlotStyle>({PlottingHelper::Primary, PlottingHelper::Secondary});
    eff_muonCosTheta.SaveAs(cutLabels, cutStyles, "binningPlots_efficiency_muonCosTheta");
    eff_muonPhi.SaveAs(cutLabels, cutStyles, "binningPlots_efficiency_muonPhi");
    eff_muonMomentum.SaveAs(cutLabels, cutStyles, "binningPlots_efficiency_muonMomentum");
    eff_pionCosTheta.SaveAs(cutLabels, cutStyles, "binningPlots_efficiency_pionCosTheta");
    eff_pionPhi.SaveAs(cutLabels, cutStyles, "binningPlots_efficiency_pionPhi");
    eff_pionMomentum.SaveAs(cutLabels, cutStyles, "binningPlots_efficiency_pionMomentum");
    eff_muonPionAngle.SaveAs(cutLabels, cutStyles, "binningPlots_efficiency_muonPionAngle");
    
    // Save the multi-plots
    multi_muonCosTheta_generic.SaveAsStacked("binningPlots_stack_generic_muonCosTheta", false, true);
    multi_muonCosTheta_golden.SaveAsStacked("binningPlots_stack_golden_muonCosTheta", false, true);
    
    multi_muonPhi_generic.SaveAsStacked("binningPlots_stack_generic_muonPhi", false, true);
    multi_muonPhi_golden.SaveAsStacked("binningPlots_stack_golden_muonPhi", false, true);
    
    multi_muonMomentum_generic.SaveAsStacked("binningPlots_stack_generic_muonMomentum", false, true);
    multi_muonMomentum_golden.SaveAsStacked("binningPlots_stack_golden_muonMomentum", false, true);
    
    multi_pionCosTheta_generic.SaveAsStacked("binningPlots_stack_generic_pionCosTheta", false, true);
    multi_pionCosTheta_golden.SaveAsStacked("binningPlots_stack_golden_pionCosTheta", false, true);
    
    multi_pionPhi_generic.SaveAsStacked("binningPlots_stack_generic_pionPhi", false, true);
    multi_pionPhi_golden.SaveAsStacked("binningPlots_stack_golden_pionPhi", false, true);
    
    multi_pionMomentum_generic.SaveAsStacked("binningPlots_stack_generic_pionMomentum", false, true);
    multi_pionMomentum_golden.SaveAsStacked("binningPlots_stack_golden_pionMomentum", false, true);

    multi_muonPionAngle_generic.SaveAsStacked("binningPlots_stack_generic_muonPionAngle", false, true);
    multi_muonPionAngle_golden.SaveAsStacked("binningPlots_stack_golden_muonPionAngle", false, true);


    // Save the resolution plots
    auto saveWithLines = [](TH2F *pHist, const std::vector<float> &binEdges, const std::string &name) {

        const auto xMin = pHist->GetXaxis()->GetXmin();
        const auto xMax = pHist->GetXaxis()->GetXmax();
        const auto yMin = pHist->GetYaxis()->GetXmin();
        const auto yMax = pHist->GetYaxis()->GetXmax();

        auto pCanvas = PlottingHelper::GetCanvas();
        pHist->Draw("colz");

        std::vector< std::shared_ptr<TLine> > lines;
        for (const auto &edge : binEdges)
        {
            lines.emplace_back(new TLine(edge, yMin, edge, yMax));
            auto &lx = lines.back();
            lx->SetLineWidth(1);
            lx->Draw();

            lines.emplace_back(new TLine(xMin, edge, xMax, edge));
            auto &ly = lines.back();
            ly->SetLineWidth(1);
            ly->Draw();
        }

        PlottingHelper::SaveCanvas(pCanvas, name);
    };

    saveWithLines(hRes_muonCosTheta_generic, config.global.muonCosTheta.binEdges, "binningPlots_resolution_generic_muonCosTheta");
    saveWithLines(hRes_muonCosTheta_golden, config.global.muonCosTheta.binEdges, "binningPlots_resolution_golden_muonCosTheta");
    
    saveWithLines(hRes_muonPhi_generic, config.global.muonPhi.binEdges, "binningPlots_resolution_generic_muonPhi");
    saveWithLines(hRes_muonPhi_golden, config.global.muonPhi.binEdges, "binningPlots_resolution_golden_muonPhi");
    
    saveWithLines(hRes_muonMomentum_generic, config.global.muonMomentum.binEdges, "binningPlots_resolution_generic_muonMomentum");
    saveWithLines(hRes_muonMomentum_golden, config.global.muonMomentum.binEdges, "binningPlots_resolution_golden_muonMomentum");
    
    saveWithLines(hRes_pionCosTheta_generic, config.global.pionCosTheta.binEdges, "binningPlots_resolution_generic_pionCosTheta");
    saveWithLines(hRes_pionCosTheta_golden, config.global.pionCosTheta.binEdges, "binningPlots_resolution_golden_pionCosTheta");
    
    saveWithLines(hRes_pionPhi_generic, config.global.pionPhi.binEdges, "binningPlots_resolution_generic_pionPhi");
    saveWithLines(hRes_pionPhi_golden, config.global.pionPhi.binEdges, "binningPlots_resolution_golden_pionPhi");
    
    saveWithLines(hRes_pionMomentum_generic, config.global.pionMomentum.binEdges, "binningPlots_resolution_generic_pionMomentum");
    saveWithLines(hRes_pionMomentum_golden, config.global.pionMomentum.binEdges, "binningPlots_resolution_golden_pionMomentum");
    
    saveWithLines(hRes_muonPionAngle_generic, config.global.muonPionAngle.binEdges, "binningPlots_resolution_generic_muonPionAngle");
    saveWithLines(hRes_muonPionAngle_golden, config.global.muonPionAngle.binEdges, "binningPlots_resolution_golden_muonPionAngle");
}

} // namespace ubcc1pi_macros
