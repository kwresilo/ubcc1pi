#include "ubcc1pi_standalone/Helpers/PlottingHelper.h"

#include "ubcc1pi_standalone/Helpers/AnalysisHelper.h"

#include <stdexcept>

#include <THStack.h>

namespace ubcc1pi
{

PlottingHelper::ParticlePlot::ParticlePlot(const std::string &xLabel, unsigned int nBins, float min, float max, bool drawErrors) :
    m_xLabel(xLabel),
    m_nBins(nBins),
    m_min(min),
    m_max(max),
    m_id(++m_lastId),
    m_cloneCount(0),
    m_drawErrors(drawErrors)
{
    for (const auto &particle : PlottingHelper::AllParticleStyles)
    {
        // Make a unique name for this plot to avoid collisionsi
        const auto nameStr = "ubcc1pi_particlePlot_" + std::to_string(m_id) + "_" + std::to_string(static_cast<int>(particle));
        const auto name = nameStr.c_str();

        auto pHist = std::make_shared<TH1F>(name, "", nBins, min, max);
        
        pHist->Sumw2();
        PlottingHelper::SetLineStyle(pHist.get(), particle);
        pHist->GetXaxis()->SetTitle(xLabel.c_str());
        pHist->GetYaxis()->SetTitle("Fraction of reco particles");
        
        m_particleToHistMap.emplace(particle, pHist);
    }
}

// -----------------------------------------------------------------------------------------------------------------------------------------
                
void PlottingHelper::ParticlePlot::Fill(const float value, const ParticleStyle &particleStyle, const float weight)
{
    auto iter = m_particleToHistMap.find(particleStyle);
    if (iter == m_particleToHistMap.end())
        throw std::invalid_argument("Input particle style is unknown");

    if (value < m_min || value > m_max)
        return;

    // Fill the histogram
    iter->second->Fill(value, weight);
}
                
// -----------------------------------------------------------------------------------------------------------------------------------------

void PlottingHelper::ParticlePlot::GetHistogramClones(std::unordered_map<ParticleStyle, TH1F*> &particleToHistCloneMap)
{
    for (const auto &particle : PlottingHelper::AllParticleStyles)
    {
        const auto pHist = m_particleToHistMap.at(particle);
        const auto pHistClone = static_cast<TH1F *>(pHist->Clone());
        const auto cloneNameStr = std::string(pHist->GetName()) + "_clone_" + std::to_string(m_cloneCount);
        const auto cloneName = cloneNameStr.c_str();
        pHistClone->SetName(cloneName);

        particleToHistCloneMap.emplace(particle, pHistClone);
    }

    m_cloneCount++;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void PlottingHelper::ParticlePlot::ScaleHistograms(std::unordered_map<ParticleStyle, TH1F*> &particleToHistCloneMap) const
{
    for (const auto &particle : PlottingHelper::AllParticleStyles)
    {
        auto pHist = particleToHistCloneMap.at(particle);
        pHist->Scale(1.f / pHist->Integral());
    }
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void PlottingHelper::ParticlePlot::SetHistogramYRanges(std::unordered_map<ParticleStyle, TH1F*> &particleToHistCloneMap) const
{
    float yMin = std::numeric_limits<float>::max();
    float yMax = -std::numeric_limits<float>::max();
    
    for (const auto &particle : PlottingHelper::AllParticleStyles)
    {
        const auto pHist = particleToHistCloneMap.at(particle);
        yMin = std::min(yMin, static_cast<float>(pHist->GetMinimum()));
        yMax = std::max(yMax, static_cast<float>(pHist->GetMaximum()));
    }
    
    // Add some padding to the top of the histogram
    yMax += (yMax - yMin) * 0.05;
    
    for (const auto &particle : PlottingHelper::AllParticleStyles)
    {
        auto pHist = particleToHistCloneMap.at(particle);
        pHist->GetYaxis()->SetRangeUser(yMin, yMax);
    }
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void PlottingHelper::ParticlePlot::SaveAs(const std::string &fileName)
{
    auto pCanvas = PlottingHelper::GetCanvas();

    // Clone the histogtams and scale them (we clone to the original hisograms can be subsequently filled)
    std::unordered_map<ParticleStyle, TH1F* > particleToHistCloneMap;
    this->GetHistogramClones(particleToHistCloneMap);
    this->ScaleHistograms(particleToHistCloneMap);
    this->SetHistogramYRanges(particleToHistCloneMap);

    // Draw the error bands if required
    bool isFirst = true;
    if (m_drawErrors)
    {
        for (const auto &particle : PlottingHelper::AllParticleStyles)
        {
            // Don't draw BNB data
            if (particle == BNBData)
                continue;

            auto pHist = particleToHistCloneMap.at(particle);
    
            if (pHist->GetEntries() == 0)
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
    for (const auto &particle : PlottingHelper::AllParticleStyles)
    {
        // Don't draw BNB data
        if (particle == BNBData)
            continue;

        auto pHist = particleToHistCloneMap.at(particle);

        if (pHist->GetEntries() == 0)
            continue;
        
        pHist->Draw(isFirst ? "hist" : "hist same");
        isFirst = false;
    }

    PlottingHelper::SaveCanvas(pCanvas, fileName);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void PlottingHelper::ParticlePlot::SaveAsStacked(const std::string &fileName)
{
    auto pCanvas = PlottingHelper::GetCanvas();

    // Clone the histogtams and scale them (we clone to the original hisograms can be subsequently filled)
    std::unordered_map<ParticleStyle, TH1F* > particleToHistCloneMap;
    this->GetHistogramClones(particleToHistCloneMap);

    // Work out if we have BNB data to plot, and if we do then get the minimum and maximum Y coordinates
    auto bnbDataHistIter = particleToHistCloneMap.find(BNBData);
    const bool hasBNBData = bnbDataHistIter != particleToHistCloneMap.end();
    auto yMin = hasBNBData ? static_cast<float>(bnbDataHistIter->second->GetMinimum()) : std::numeric_limits<float>::max();
    auto yMax = hasBNBData ? static_cast<float>(bnbDataHistIter->second->GetMaximum()) : -std::numeric_limits<float>::max();

    // Sum the non BNB data histograms to get the "MC" total
    const auto nameTotalStr = "ubcc1pi_particlePlot_" + std::to_string(m_id) + "_total";
    const auto nameTotal = nameTotalStr.c_str();
    auto pHistTotal = std::make_shared<TH1F>(nameTotal, "", m_nBins, m_min, m_max);
    pHistTotal->Sumw2();

    for (const auto &particle : PlottingHelper::AllParticleStyles)
    {
        // Don't add BNB data to the total, this will be drawn seperately
        if (particle == BNBData)
            continue;

        auto pHist = particleToHistCloneMap.at(particle);
        
        std::cout << "Histogram: " << pHist->GetName() << " - nEntries = " << pHist->GetEntries() << std::endl;

        if (pHist->GetEntries() == 0)
            continue;

        pHistTotal->Add(pHist);
    }

    // Get the maximum and minimum Y coordinates
    yMin = std::min(yMin, static_cast<float>(pHistTotal->GetMinimum()));
    yMax = std::min(yMax, static_cast<float>(pHistTotal->GetMaximum()));
    yMax += (yMax - yMin) * 0.05;

    // Draw the stacked histogram
    const auto nameStackStr = "ubcc1pi_particlePlot_" + std::to_string(m_id) + "_stack";
    const auto nameStack = nameStackStr.c_str();
    auto pHistStack = std::make_shared<THStack>(nameStack, "");
    
    for (const auto &particle : PlottingHelper::AllParticleStyles)
    {
        // Don't draw BNB data
        if (particle == BNBData)
            continue;

        auto pHist = particleToHistCloneMap.at(particle);

        if (pHist->GetEntries() == 0)
            continue;
        
        const auto col = pHist->GetLineColor();
        pHist->SetLineColorAlpha(col, 0.f);
        pHist->SetFillStyle(1001);
        pHist->SetFillColor(col);
       
        pHistStack->Add(pHist);
    }

    pHistStack->Draw("hist");
    pHistStack->GetYaxis()->SetRangeUser(yMin, yMax);
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

    PlottingHelper::SaveCanvas(pCanvas, fileName);

    // Now make the ratio plot
    auto pCanvasRatio = PlottingHelper::GetCanvas(960, 270);
    auto pHistBNBData = bnbDataHistIter->second;
    pHistBNBData->GetYaxis()->SetTitle("Beam on / (Overlay + Beam off)");
    pHistBNBData->Divide(pHistTotal.get());
    pHistBNBData->Draw("e1");
    PlottingHelper::SaveCanvas(pCanvasRatio, fileName + "_ratio");
}

// -----------------------------------------------------------------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------------------------------------------------------------
        
PlottingHelper::ParticleStyle PlottingHelper::GetParticleStyle(const Event::Reco::Particle &particle, const std::vector<Event::Truth::Particle> &truthParticles)
{
    if (!particle.hasMatchedMCParticle.IsSet())
        return External;

    if (!particle.hasMatchedMCParticle())
        return External;

    if (truthParticles.empty())
        return External;

    // Get the true pdg code applying the visibility conditions
    auto truePdgCode = -std::numeric_limits<int>::max();
    bool isGolden = false;
    try
    {
        const auto truthParticle = AnalysisHelper::GetBestMatchedTruthParticle(particle, truthParticles);
        truePdgCode = truthParticle.pdgCode();
        isGolden = AnalysisHelper::IsGolden(truthParticle);
    }
    catch (const std::logic_error &)
    {
        return External;
    }

    switch (truePdgCode)
    {
        case 13:
            return Muon;
        case 2212:
            return Proton;
        case 211:
            return (isGolden ? GoldenPion : NonGoldenPion);
        case -211:
            return PiMinus;
        case 11:
            return Electron;
        case 22:
            return Photon;
        default:
            return Other;
    }
}

// -----------------------------------------------------------------------------------------------------------------------------------------

std::shared_ptr<TCanvas> PlottingHelper::GetCanvas(const unsigned int width, const unsigned int height)
{
    const auto nameStr = "canvas_" + std::to_string(++m_lastCanvasId);
    const auto name = nameStr.c_str();

    std::shared_ptr<TCanvas> pCanvas = std::make_shared<TCanvas>(name, "", width, height);

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

} // namespace ubcc1pi
