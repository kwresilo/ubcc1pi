/**
 *  @file  ubcc1pi_standalone/Helpers/PlottingHelper.cxx
 *
 *  @brief The implementation file of the plotting helper class 
 */

#include "ubcc1pi_standalone/Helpers/PlottingHelper.h"

#include "ubcc1pi_standalone/Helpers/AnalysisHelper.h"

#include <stdexcept>

#include <THStack.h>
#include <TLine.h>
#include <TStyle.h>

namespace ubcc1pi
{

PlottingHelper::MultiPlot::MultiPlot(const std::string &xLabel, const std::string &yLabel, unsigned int nBins, float min, float max, bool drawErrors) :
    MultiPlot(xLabel, yLabel, PlottingHelper::GenerateUniformBinEdges(nBins, min, max), drawErrors)
{
}

// -----------------------------------------------------------------------------------------------------------------------------------------

PlottingHelper::MultiPlot::MultiPlot(const std::string &xLabel, const std::string &yLabel, const std::vector<float> &binEdges, bool drawErrors) :
    /// @cond Doxygen can't handle this initilizer list
    m_xLabel(xLabel),
    m_nBins(binEdges.size() - 1),
    m_min(binEdges.front()),
    m_max(binEdges.back()),
    m_binEdges(binEdges),
    m_id(++m_lastId),
    m_cloneCount(0),
    m_drawErrors(drawErrors)
    /// @endcond
{

    for (const auto &style : PlottingHelper::AllPlotStyles)
    {
        // Make a unique name for this plot to avoid collisionsi
        const auto nameStr = "ubcc1pi_multiPlot_" + std::to_string(m_id) + "_" + std::to_string(static_cast<int>(style));
        const auto name = nameStr.c_str();

        auto pHist = std::make_shared<TH1F>(name, "", m_nBins, m_binEdges.data());
        
        pHist->Sumw2();
        PlottingHelper::SetLineStyle(pHist.get(), style);
        //pHist->GetXaxis()->SetTitle(xLabel.c_str());
        //pHist->GetYaxis()->SetTitle(yLabel.c_str());
        
        m_plotToHistMap.emplace(style, pHist);
    }
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void PlottingHelper::MultiPlot::SetIntegerBinLabels()
{
    // Check that the bin edges lie on an integer and get those integers
    std::vector<int> intBinEdges;
    for (const auto &edge : m_binEdges)
    {
        if (std::abs(edge - std::round(edge)) > std::numeric_limits<float>::epsilon())
            throw std::logic_error("MultiPlot::SetIntegerBinLabels - found bin edge of: " + std::to_string(edge) + ", which isn't an integer");

        intBinEdges.push_back(static_cast<int>(std::round(edge)));
    }

    // Get the bin labels
    std::vector<std::string> labels;
    for (unsigned int iBin = 0; iBin < m_nBins; ++iBin)
    {
        const auto lower = intBinEdges.at(iBin);
        const auto upper = intBinEdges.at(iBin + 1);
        const auto width = upper - lower;

        // If this bin contains one integer, then just label it with that value
        if (width == 1)
        {
            labels.push_back(std::to_string(lower));
            continue;
        }

        // If the bin contained multiple integers then use the range
        labels.push_back(std::to_string(lower) + " -> " + std::to_string(upper - 1));
    }

    // Set the bin labels
    this->SetBinLabels(labels);
}

// -----------------------------------------------------------------------------------------------------------------------------------------
                
void PlottingHelper::MultiPlot::SetBinLabels(const std::vector<std::string> &labels)
{
    if (labels.size() != m_nBins)
        throw std::invalid_argument("MultiPlot::SetBinLabels - there are " + std::to_string(m_nBins) + " bins, but only " + std::to_string(labels.size()) + " labels were provided");

    for (const auto &style : PlottingHelper::AllPlotStyles)
    {
        auto &pHist = m_plotToHistMap.at(style);
        auto pXAxis = pHist->GetXaxis();

        for (unsigned int iBin = 0; iBin < m_nBins; ++iBin)
        {
            // ATTN root indexes it's bins from 1
            pXAxis->SetBinLabel(iBin + 1, labels.at(iBin).c_str());
        }
    }
}

// -----------------------------------------------------------------------------------------------------------------------------------------
                
void PlottingHelper::MultiPlot::AddCutLine(const float value)
{
    if (value < m_min || value > m_max)
        throw std::invalid_argument("MultiPlot::AddCutLine - supplied value: " + std::to_string(value) + " is out of range: " + std::to_string(m_min) + " -> " + std::to_string(m_max));

    m_cutValues.push_back(value);
}

// -----------------------------------------------------------------------------------------------------------------------------------------
                
void PlottingHelper::MultiPlot::Fill(const float value, const PlotStyle &plotStyle, const float weight)
{
    auto iter = m_plotToHistMap.find(plotStyle);
    if (iter == m_plotToHistMap.end())
        throw std::invalid_argument("Input plot style is unknown");

    if (value < m_min || value > m_max)
        return;

    // Fill the histogram
    iter->second->Fill(value, weight);
}
                
// -----------------------------------------------------------------------------------------------------------------------------------------

void PlottingHelper::MultiPlot::GetHistogramClones(std::unordered_map<PlotStyle, TH1F*> &plotToHistCloneMap)
{
    for (const auto &style : PlottingHelper::AllPlotStyles)
    {
        const auto pHist = m_plotToHistMap.at(style);
        const auto pHistClone = static_cast<TH1F *>(pHist->Clone());
        const auto cloneNameStr = std::string(pHist->GetName()) + "_clone_" + std::to_string(m_cloneCount);
        const auto cloneName = cloneNameStr.c_str();
        pHistClone->SetName(cloneName);

        plotToHistCloneMap.emplace(style, pHistClone);
    }

    m_cloneCount++;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void PlottingHelper::MultiPlot::ScaleHistograms(std::unordered_map<PlotStyle, TH1F*> &plotToHistCloneMap, const bool scaleByBinWidth) const
{
    for (const auto &style : PlottingHelper::AllPlotStyles)
    {
        auto pHist = plotToHistCloneMap.at(style);
        const auto option = scaleByBinWidth ? "width" : "";
        pHist->Scale(1.f, option);
        pHist->Scale(1.f / pHist->Integral());
    }
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void PlottingHelper::MultiPlot::SetHistogramYRanges(const unsigned int minEntriesToDraw, std::unordered_map<PlotStyle, TH1F*> &plotToHistCloneMap) const
{
    float yMin = std::numeric_limits<float>::max();
    float yMax = -std::numeric_limits<float>::max();

    bool shouldScale = false;
    for (const auto &style : PlottingHelper::AllPlotStyles)
    {
        const auto pHist = plotToHistCloneMap.at(style);

        if (pHist->GetEntries() < minEntriesToDraw || pHist->GetEntries() == 0u)
            continue;

        yMin = std::min(yMin, static_cast<float>(pHist->GetMinimum()));
        yMax = std::max(yMax, static_cast<float>(pHist->GetMaximum()));
        shouldScale = true;
    }
    
    if (!shouldScale)
        return;

    // Add some padding to the top of the histogram
    yMax += (yMax - yMin) * 0.05;
    
    for (const auto &style : PlottingHelper::AllPlotStyles)
    {
        auto pHist = plotToHistCloneMap.at(style);
        pHist->GetYaxis()->SetRangeUser(yMin, yMax);
    }
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void PlottingHelper::MultiPlot::SaveAs(const std::string &fileName, const bool useLogX, const bool scaleByBinWidth, const unsigned int minEntriesToDraw)
{
    auto pCanvas = PlottingHelper::GetCanvas();

    if (useLogX)
    {
        pCanvas->SetLogx();
    }

    // Clone the histogtams and scale them (we clone to the original hisograms can be subsequently filled)
    std::unordered_map<PlotStyle, TH1F* > plotToHistCloneMap;
    this->GetHistogramClones(plotToHistCloneMap);
    this->ScaleHistograms(plotToHistCloneMap, scaleByBinWidth);
    this->SetHistogramYRanges(minEntriesToDraw, plotToHistCloneMap);
    
    // Draw the error bands if required
    bool isFirst = true;
    if (m_drawErrors)
    {
        for (const auto &style : PlottingHelper::AllPlotStyles)
        {
            // The points already have error bars
            if (PlottingHelper::ShouldUsePoints(style))
                continue;

            // Don't draw BNB data
            if (style == BNBData)
                continue;

            auto pHist = plotToHistCloneMap.at(style);
    
            if (pHist->GetEntries() < minEntriesToDraw || pHist->GetEntries() == 0u)
                continue;
        
            // Draw a clone of the histogram so we can safely change it's style
            auto pHistClone = static_cast<TH1F *>(pHist->Clone());
            const auto col = pHistClone->GetLineColor();
            pHistClone->SetFillStyle(1001);
            pHistClone->SetLineColorAlpha(col, 0.f);
            pHistClone->SetFillColorAlpha(col, 0.3f);

            pHistClone->Draw(isFirst ? "e2" : "e2 same");
            isFirst = false;
        }
    }
        
    // Draw the histogram itself
    for (const auto &style : PlottingHelper::AllPlotStyles)
    {
        // Don't draw BNB data
        if (style == BNBData)
            continue;

        auto pHist = plotToHistCloneMap.at(style);

        if (pHist->GetEntries() < minEntriesToDraw || pHist->GetEntries() == 0u)
            continue;


        const bool usePoints = PlottingHelper::ShouldUsePoints(style);
        const TString tstyle = usePoints ? "e1" : "hist";
        pHist->Draw(isFirst ? tstyle : tstyle + "same");
        isFirst = false;
    }

    PlottingHelper::SaveCanvas(pCanvas, fileName);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void PlottingHelper::MultiPlot::SaveAsStacked(const std::string &fileName, const bool useLogX, const bool scaleByBinWidth)
{
    auto pCanvas = PlottingHelper::GetCanvas();

    if (useLogX)
    {
        pCanvas->SetLogx();
    }

    // Clone the histogtams and scale them (we clone to the original hisograms can be subsequently filled)
    std::unordered_map<PlotStyle, TH1F* > plotToHistCloneMap;
    this->GetHistogramClones(plotToHistCloneMap);

    // Scale by bin width
    if (scaleByBinWidth)
    {
        for (auto &entry : plotToHistCloneMap)
            entry.second->Scale(1.f, "width");
    }

    // Work out if we have BNB data to plot, and if we do then get the minimum and maximum Y coordinates
    auto bnbDataHistIter = plotToHistCloneMap.find(BNBData);
    const bool hasBNBData = bnbDataHistIter != plotToHistCloneMap.end();

    auto yMin = hasBNBData ? static_cast<float>(bnbDataHistIter->second->GetMinimum()) : std::numeric_limits<float>::max();
    auto yMax = hasBNBData ? static_cast<float>(bnbDataHistIter->second->GetMaximum()) : -std::numeric_limits<float>::max();

    // Sum the non BNB data histograms to get the "MC" total
    const auto nameTotalStr = "ubcc1pi_multiPlot_" + std::to_string(m_id) + "_total";
    const auto nameTotal = nameTotalStr.c_str();
    auto pHistTotal = std::make_shared<TH1F>(nameTotal, "", m_nBins, m_binEdges.data());
    pHistTotal->Sumw2();

    bool isFirst = true;
    for (const auto &style : PlottingHelper::AllPlotStyles)
    {
        // Don't add BNB data to the total, this will be drawn seperately
        if (style == BNBData)
            continue;

        auto pHist = plotToHistCloneMap.at(style);

        if (pHist->GetEntries() == 0)
            continue;

        if (isFirst)
        {
            // If we are using alphanumeric bin labels, then we need to copy those labels to the total histogram or root will complain
            const auto pAxis = pHist->GetXaxis();
            if (pAxis->IsAlphanumeric())
            {
                for (unsigned int iBin = 1; iBin <= static_cast<unsigned int>(pAxis->GetNbins()); ++iBin)
                    pHistTotal->GetXaxis()->SetBinLabel(iBin, pAxis->GetBinLabel(iBin));
            }

            isFirst = false;
        }

        pHistTotal->Add(pHist);
    }
    
    // Get the maximum and minimum Y coordinates
    yMin = std::min(yMin, static_cast<float>(pHistTotal->GetMinimum()));
    yMax = std::max(yMax, static_cast<float>(pHistTotal->GetMaximum()));
    yMax += (yMax - yMin) * 0.05;
        
    // Draw the stacked histogram
    const auto nameStackStr = "ubcc1pi_plotPlot_" + std::to_string(m_id) + "_stack";
    const auto nameStack = nameStackStr.c_str();
    auto pHistStack = std::make_shared<THStack>(nameStack, "");
    
    for (const auto &style : PlottingHelper::AllPlotStyles)
    {
        // Don't draw BNB data
        if (style == BNBData)
            continue;

        auto pHist = plotToHistCloneMap.at(style);
        
        if (pHist->GetEntries() == 0)
            continue;
        
        const auto col = pHist->GetLineColor();
        pHist->SetLineColorAlpha(col, 0.f);
        pHist->SetFillStyle(1001);
        pHist->SetFillColor(col);
       
        pHistStack->Add(pHist);
    }

    pHistStack->SetMinimum(0.f);
    pHistStack->SetMaximum(yMax / (1+gStyle->GetHistTopMargin()));
    pHistStack->Draw("hist");
           
    // Draw the error bands on the stack if required
    if (m_drawErrors)
    {
        pHistTotal->SetFillStyle(1001);
        pHistTotal->SetLineColorAlpha(kBlack, 0.f);
        pHistTotal->SetFillColorAlpha(kBlack, 0.3f);
        pHistTotal->Draw("e2 same");
    }

    // Draw the BNB data
    if (hasBNBData)
    {
        bnbDataHistIter->second->Draw("e1 same");
    }

    // Draw the cut lines
    for (const auto &cutValue : m_cutValues)
    {
        TLine *pLine = new TLine(cutValue, 0.f, cutValue, yMax);
        pLine->SetLineWidth(2);
        pLine->Draw();
    }

    PlottingHelper::SaveCanvas(pCanvas, fileName);
    
    // Now make the ratio plot
    auto pCanvasRatio = PlottingHelper::GetCanvas(960, 270);
    
    if (useLogX)
    {
        pCanvasRatio->SetLogx();
    }

    auto pHistBNBData = static_cast<TH1F *>(bnbDataHistIter->second->Clone());
    //pHistBNBData->GetYaxis()->SetTitle("Beam on / (Overlay + Beam off)");
    pHistBNBData->Divide(pHistTotal.get());

    // Set the y-axis range
    float minRatio = 0.75f;
    float maxRatio = 1.25f;
    for (unsigned int i = 1; i <= static_cast<unsigned int>(pHistBNBData->GetNbinsX()); ++i)
    {
        minRatio = std::min(minRatio, static_cast<float>(pHistBNBData->GetBinContent(i) - pHistBNBData->GetBinError(i)));
        maxRatio = std::max(maxRatio, static_cast<float>(pHistBNBData->GetBinContent(i) + pHistBNBData->GetBinError(i)));
    }
    const auto padding = (maxRatio - minRatio) * 0.05f;
    const auto ratioYMin = std::max(0.f, minRatio - padding);
    const auto ratioYMax = maxRatio + padding;
    pHistBNBData->GetYaxis()->SetRangeUser(ratioYMin, ratioYMax);

    // Draw
    pHistBNBData->Draw("e1");
    pCanvasRatio->Update();

    // Add the lines at 0.8, 1.0 and 1.2
    const auto ratioXMin = pHistBNBData->GetXaxis()->GetXmin();
    const auto ratioXMax = pHistBNBData->GetXaxis()->GetXmax();
    TLine *l=new TLine(ratioXMin,1.0,ratioXMax,1.0);
    TLine *lPlus=new TLine(ratioXMin,1.2,ratioXMax,1.2);
    TLine *lMinus=new TLine(ratioXMin,0.8,ratioXMax,0.8);
    l->SetLineColor(kBlue);
    lPlus->SetLineColor(kBlue);
    lMinus->SetLineColor(kBlue);
    lPlus->SetLineStyle(2);
    lMinus->SetLineStyle(2);
    l->Draw();
    lPlus->Draw();
    lMinus->Draw();
    
    // Draw the cut lines
    for (const auto &cutValue : m_cutValues)
    {
        TLine *pLine = new TLine(cutValue, ratioYMin, cutValue, ratioYMax);
        pLine->SetLineWidth(2);
        pLine->Draw();
    }

    PlottingHelper::SaveCanvas(pCanvasRatio, fileName + "_ratio");
}

// -----------------------------------------------------------------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------------------------------------------------------------

PlottingHelper::EfficiencyPlot::EfficiencyPlot(const std::string &xLabel, unsigned int nBins, float min, float max, const std::vector<string> &cuts, bool drawErrors) :
    EfficiencyPlot(xLabel, PlottingHelper::GenerateUniformBinEdges(nBins, min, max), cuts, drawErrors)
{
}

// -----------------------------------------------------------------------------------------------------------------------------------------

PlottingHelper::EfficiencyPlot::EfficiencyPlot(const std::string &xLabel, const std::vector<float> &binEdges, const std::vector<string> &cuts, bool drawErrors) :
    /// @cond Doxygen can't handle this initilizer list
    m_xLabel(xLabel),
    m_binEdges(binEdges),
    m_nBins(binEdges.size() - 1),
    m_min(binEdges.front()),
    m_max(binEdges.back()),
    m_cuts(cuts),
    m_drawErrors(drawErrors),
    m_id(++m_lastId)
    /// @endcond
{
    if (m_cuts.empty())
        throw std::logic_error("EfficiencyPlot::EfficiencyPlot - There are no cuts set!");

    // Make histograms for the cut
    for (unsigned int iCut = 0; iCut < m_cuts.size(); ++iCut)
    {
        const auto &cut = m_cuts.at(iCut);

        if (std::count(m_cuts.begin(), m_cuts.end(), cut) != 1)
            throw std::invalid_argument("EfficiencyPlot::EfficiencyPlot - Repeated cut: \"" + cut + "\"");

        // Define the plot names for root
        const std::string nameStr = "ubcc1pi_efficiencyPlot_" + std::to_string(m_id) + "_" + cut;
        const auto nameNumerator = nameStr + "_numerator";
        const auto nameDenominator = nameStr + "_denominator";

        // Make the pair of plots, one for the numerator one for denominator
        auto plotPair = std::pair< std::shared_ptr<TH1F>, std::shared_ptr<TH1F> >(
            std::make_shared<TH1F>(nameNumerator.c_str(), "", m_nBins, m_binEdges.data()),
            std::make_shared<TH1F>(nameDenominator.c_str(), "", m_nBins, m_binEdges.data()));
        
        // Setup the plots
        for (auto &pHist : {plotPair.first, plotPair.second})
        {
            pHist->Sumw2();
            //pHist->GetXaxis()->SetTitle(xLabel.c_str());
            //pHist->GetYaxis()->SetTitle("Efficiency");
        }
        
        m_cutToPlotsMap.emplace(cut, plotPair);
    }
}

// -----------------------------------------------------------------------------------------------------------------------------------------
                
void PlottingHelper::EfficiencyPlot::AddCutLine(const float value)
{
    if (value < m_min || value > m_max)
        throw std::invalid_argument("EfficiencyPlot::AddCutLine - supplied value: " + std::to_string(value) + " is out of range: " + std::to_string(m_min) + " -> " + std::to_string(m_max));

    m_cutValues.push_back(value);
}

// -----------------------------------------------------------------------------------------------------------------------------------------
       
void PlottingHelper::EfficiencyPlot::AddEvent(const float value, const float weight, const std::string &cut, const bool passedCut)
{
    auto plotPairIter = m_cutToPlotsMap.find(cut);

    if (plotPairIter == m_cutToPlotsMap.end())
        throw std::invalid_argument("EfficiencyPlot::AddEvent - Unknown cut: \"" + cut + "\"");

    auto &pHistNumerator = plotPairIter->second.first;
    auto &pHistDenominator = plotPairIter->second.second;

    if (passedCut)
        pHistNumerator->Fill(value, weight);

    pHistDenominator->Fill(value, weight);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void PlottingHelper::EfficiencyPlot::SaveAs(const std::vector<std::string> &cuts, const std::vector<int> &colors, const std::string &fileName)
{
    if (cuts.size() != colors.size())
        throw std::invalid_argument("EfficiencyPlot::SaveAs - number of cuts and colors doesn't match");
        
    for (const auto &cut : cuts)
    {
        if (std::find(m_cuts.begin(), m_cuts.end(), cut) == m_cuts.end())
            throw std::invalid_argument("EfficiencyPlot::SaveAs - Unknown cut: \"" + cut + "\"");
    }

    auto pCanvas = PlottingHelper::GetCanvas();

    // Save the raw histogram - just use the denominator of the first cut
    auto pHistRaw = static_cast<TH1F *>(m_cutToPlotsMap.at(m_cuts.front()).second->Clone());
    //pHistRaw->GetYaxis()->SetTitle("Number of events");
    PlottingHelper::SetLineStyle(pHistRaw, Default);

    if (m_drawErrors)
    {
        auto pHistClone = static_cast<TH1F *>(pHistRaw->Clone());
        const auto deafultCol = PlottingHelper::GetColor(Default);
        pHistClone->SetFillStyle(1001);
        pHistClone->SetLineColorAlpha(deafultCol, 0.f);
        pHistClone->SetFillColorAlpha(deafultCol, 0.3f);

        pHistClone->Draw("e2");
    }
    
    pHistRaw->Draw(m_drawErrors ? "hist same" : "hist"); 
    PlottingHelper::SaveCanvas(pCanvas, fileName + "_raw");

    // Get the efficiency histograms
    std::vector<TH1F *> efficiencyHists;
    std::vector<TH1F *> efficiencyErrorHists;
    for (unsigned int i = 0; i < cuts.size(); ++i)
    {   
        const auto &cut = cuts.at(i);
        const auto col = colors.at(i);

        const auto &plotPair = m_cutToPlotsMap.at(cut);
        const auto &pHistNumerator = plotPair.first;
        const auto &pHistDenominator = plotPair.second;

        // Clone the numerator and divide by the denominator to get the efficiency
        auto pHistEfficiency = static_cast<TH1F *>(pHistNumerator->Clone());
        pHistEfficiency->Divide(pHistDenominator.get());
        PlottingHelper::SetLineStyle(pHistEfficiency, col);
        efficiencyHists.push_back(pHistEfficiency);
    
        if (m_drawErrors)
        {
            auto pHistClone = static_cast<TH1F *>(pHistEfficiency->Clone());
            pHistClone->SetFillStyle(1001);
            pHistClone->SetLineColorAlpha(col, 0.f);
            pHistClone->SetFillColorAlpha(col, 0.3f);

            efficiencyErrorHists.push_back(pHistClone);
        }
    }

    // Get the y-range
    float yMin = std::numeric_limits<float>::max();
    float yMax = -std::numeric_limits<float>::max();

    for (const auto &pHist : efficiencyHists)
    {
        yMin = std::min(yMin, static_cast<float>(pHist->GetMinimum()));
        yMax = std::max(yMax, static_cast<float>(pHist->GetMaximum()));
    }
    
    // Add some padding to the top of the histogram
    yMax += (yMax - yMin) * 0.05;

    //// TEST
    yMin = 0.f;
    yMax = 1.05f;
    //// END TEST
    
    // Make the cut lines
    std::vector< shared_ptr<TLine> > lines;
    for (const auto &cutValue : m_cutValues)
    {
        lines.emplace_back(new TLine(cutValue, yMin, cutValue, yMax));
        auto &pLine = lines.back();
        pLine->SetLineWidth(2);
    }

    // Draw the individual efficiency plots for each cut
    for (unsigned int i = 0; i < cuts.size(); ++i)
    {
        const auto &cut = cuts.at(i);

        if (m_drawErrors)
        {
            auto pHistErr = efficiencyErrorHists.at(i);
            pHistErr->GetYaxis()->SetRangeUser(yMin, yMax); 
            pHistErr->Draw("e2");
        }

        auto pHist = efficiencyHists.at(i);
        pHist->GetYaxis()->SetRangeUser(yMin, yMax);
        pHist->Draw(!m_drawErrors ? "hist" : "hist same");

        // Draw the cut value lines
        for (const auto &pLine : lines)
            pLine->Draw();
    
        PlottingHelper::SaveCanvas(pCanvas, fileName + "_" + std::to_string(i) + "_" + cut);
    }

    // Draw the error histograms all together
    bool isFirst = true;
    for (const auto &pHist : efficiencyErrorHists)
    {
        pHist->GetYaxis()->SetRangeUser(yMin, yMax);
        pHist->Draw(isFirst ? "e2" : "e2 same");
        isFirst = false;
    }
    
    // Draw the histogram lines all together
    for (const auto &pHist : efficiencyHists)
    {
        pHist->GetYaxis()->SetRangeUser(yMin, yMax);
        pHist->Draw(isFirst ? "hist" : "hist same");
        isFirst = false;
    }
        
    // Draw the cut value lines
    for (const auto &pLine : lines)
        pLine->Draw();
    
    PlottingHelper::SaveCanvas(pCanvas, fileName);
}

// -----------------------------------------------------------------------------------------------------------------------------------------
                
void PlottingHelper::EfficiencyPlot::SaveAs(const std::vector<std::string> &cuts, const std::vector<PlotStyle> &styles, const std::string &fileName)
{
    // Convert the styles to colors
    std::vector<int> colors;
    for (const auto &style : styles)
        colors.push_back(PlottingHelper::GetColor(style));
    
    this->SaveAs(cuts, colors, fileName);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void PlottingHelper::EfficiencyPlot::SaveAs(const std::string &fileName)
{
    // Get the colors - one for each plot
    const auto palette = PlottingHelper::GetColorVector();
    std::vector<int> colors;
    for (unsigned int i = 0; i < m_cuts.size(); ++i)
        colors.push_back(palette.at(i % palette.size()));

    this->SaveAs(m_cuts, colors, fileName);
}

// -----------------------------------------------------------------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------------------------------------------------------------
        
PlottingHelper::PlotStyle PlottingHelper::GetPlotStyle(const Event::Reco::Particle &particle, const AnalysisHelper::SampleType &sampleType, const std::vector<Event::Truth::Particle> &truthParticles, const bool usePoints, const bool useAbsPdg)
{
    const auto externalType = usePoints ? ExternalPoints : External;
    const auto dirtType = usePoints ? DirtPoints : Dirt;

    if (sampleType == AnalysisHelper::DataBNB)
        return BNBData;

    if (sampleType == AnalysisHelper::Dirt)
        return dirtType;

    if (sampleType == AnalysisHelper::DataEXT)
        return externalType;

    if (!particle.hasMatchedMCParticle.IsSet())
        return externalType;

    if (!particle.hasMatchedMCParticle())
        return externalType;

    if (truthParticles.empty())
        return externalType;

    // Get the true pdg code applying the visibility conditions
    auto truePdgCode = -std::numeric_limits<int>::max();
    bool isGolden = false;
    try
    {
        const auto truthParticle = AnalysisHelper::GetBestMatchedTruthParticle(particle, truthParticles);
        truePdgCode = useAbsPdg ? std::abs(truthParticle.pdgCode()) : truthParticle.pdgCode();
        isGolden = AnalysisHelper::IsGolden(truthParticle);
    }
    catch (const std::logic_error &)
    {
        return externalType;
    }

    switch (truePdgCode)
    {
        case 13:
            return usePoints ? MuonPoints : Muon;
        case 2212:
            return usePoints ? ProtonPoints : Proton;
        case 211:
            return (isGolden ? (usePoints ? GoldenPionPoints : GoldenPion) : (usePoints ? NonGoldenPionPoints : NonGoldenPion));
        case -211:
            return usePoints ? PiMinusPoints : PiMinus;
        case 11:
            return usePoints ? ElectronPoints : Electron;
        case 111: // TODO Deal with pi0s properly when re-producing files
        case 22:
            return usePoints ? PhotonPoints : Photon;
        default:
            return usePoints ? OtherPoints : Other;
    }
}

// -----------------------------------------------------------------------------------------------------------------------------------------
        
PlottingHelper::PlotStyle PlottingHelper::GetPlotStyle(const AnalysisHelper::SampleType &sampleType, const std::shared_ptr<Event> &pEvent, const bool useAbsPdg)
{
    if (sampleType == AnalysisHelper::Dirt)
        return Dirt;
    
    if (sampleType == AnalysisHelper::DataEXT)
        return External;
    
    if (sampleType == AnalysisHelper::DataBNB)
        return BNBData;

    if (!AnalysisHelper::IsFiducial(pEvent->truth.nuVertex()))
        return NonFiducial;
    
    // Count the particles by PDG code
    std::vector<int> foundPdgs;
    std::unordered_map<int, unsigned int> pdgCodeCountMap;
    AnalysisHelper::GetPdgCodeCountMap(pEvent->truth.particles, useAbsPdg, foundPdgs, pdgCodeCountMap);

    // Quick utility lambda
    const auto CountByPdg = [&](const int pdg) {
        if (std::find(foundPdgs.begin(), foundPdgs.end(), pdg) == foundPdgs.end())
            return 0u;

        return pdgCodeCountMap.at(pdg);
    };

    const auto hasMuon = CountByPdg(13) > 0;
    const auto hasElectron = CountByPdg(11) > 0;

    if (!hasMuon && hasElectron)
        return Nue;
    
    if (!hasMuon && !hasElectron)
        return NC;

    const auto nPiPlus = CountByPdg(211);
    const auto nPiMinus = CountByPdg(-211);
    const auto nPiZero = CountByPdg(111);

    const auto nPion = nPiPlus + nPiMinus + nPiZero;
    const auto nCharedPion = nPiPlus + nPiMinus;

    if (nPion == 0)
        return NumuCC0Pi;

    if (nPiZero == 1 && nCharedPion == 0)
        return NumuCC1PiZero;

    const auto hasGoldenPion = (AnalysisHelper::CountGoldenParticlesWithPdgCode(pEvent->truth.particles, 211, useAbsPdg) != 0);

    if (nCharedPion == 1 && nPiZero == 0)
    {
        if (hasGoldenPion)
            return NumuCC1PiChargedGolden;

        return NumuCC1PiChargedNonGolden;
    }
    
    return NumuCCOther;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

std::shared_ptr<TCanvas> PlottingHelper::GetCanvas(const unsigned int width, const unsigned int height)
{
    const auto nameStr = "canvas_" + std::to_string(++m_lastCanvasId);
    const auto name = nameStr.c_str();

    std::shared_ptr<TCanvas> pCanvas = std::make_shared<TCanvas>(name, "", width, height);
    pCanvas->cd();

    return pCanvas;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void PlottingHelper::SaveCanvas(std::shared_ptr<TCanvas> &pCanvas, const std::string &fileName)
{
    for (const auto &ext : {"png", "pdf", "C"})
    {
        pCanvas->SaveAs((fileName + "." + ext).c_str());
    }
}

// -----------------------------------------------------------------------------------------------------------------------------------------

std::vector<float> PlottingHelper::GenerateUniformBinEdges(const unsigned int nBins, const float min, const float max)
{
    if (min >= max)
        throw std::invalid_argument("PlottingHelper::GenerateUniformBinEdges - The range: " + std::to_string(min) + " -> " + std::to_string(max) + ", is invalid");

    std::vector<float> outVect;

    for (unsigned int i = 0; i <= nBins; ++i)
        outVect.emplace_back(min + (max - min) * static_cast<float>(i) / static_cast<float>(nBins));

    return outVect;
}

// -----------------------------------------------------------------------------------------------------------------------------------------
        
std::vector<float> PlottingHelper::GenerateLogBinEdges(const unsigned int nBins, const float min, const float max)
{
    if (std::min(min, max) <= 0.f)
        throw std::invalid_argument("PlottingHelper::GenerateUniformBinEdges - The range: " + std::to_string(min) + " -> " + std::to_string(max) + ", includes region <= 0");

    // Get uniformly spaced bins in log-space
    const auto logMin = std::log10(min);
    const auto logMax = std::log10(max);
    auto logBinEdges = PlottingHelper::GenerateUniformBinEdges(nBins, logMin, logMax);

    // Transform these edges into non-log space
    std::vector<float> outVect;
    for (const auto &logBinEdge : logBinEdges)
    {
        outVect.push_back(std::pow(10, logBinEdge));
    }

    return outVect;
}


} // namespace ubcc1pi
