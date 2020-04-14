#include "ubcc1pi_standalone/Objects/FileReader.h"

#include "ubcc1pi_standalone/Helpers/PlottingHelper.h"
#include "ubcc1pi_standalone/Helpers/AnalysisHelper.h"

using namespace ubcc1pi;

int CCInclusiveTruthStudy(const std::string &overlayFileName)
{
    FileReader reader(overlayFileName);
    auto pEvent = reader.GetBoundEventAddress();

    auto hMuonMomentum = std::make_shared<TH1F>("hMuonMomentum", "", 100, 0, 1.2);
    auto hProtonMomentum = std::make_shared<TH1F>("hProtonMomentum", "", 100, 0, 0.6);
    auto hGoldenPionMomentum = std::make_shared<TH1F>("hGoldenPionMomentum", "", 100, 0, 0.8);
    auto hNonGoldenPionMomentum = std::make_shared<TH1F>("hNonGoldenPionMomentum", "", 100, 0, 0.8);

    auto hMuonMomentumSel = std::make_shared<TH1F>("hMuonMomentumSel", "", 100, 0, 1.2);
    auto hProtonMomentumSel = std::make_shared<TH1F>("hProtonMomentumSel", "", 100, 0, 0.6);
    auto hGoldenPionMomentumSel = std::make_shared<TH1F>("hGoldenPionMomentumSel", "", 100, 0, 0.8);
    auto hNonGoldenPionMomentumSel = std::make_shared<TH1F>("hNonGoldenPionMomentumSel", "", 100, 0, 0.8);

    hMuonMomentum->Sumw2();
    hProtonMomentum->Sumw2();
    hGoldenPionMomentum->Sumw2();
    hNonGoldenPionMomentum->Sumw2();
    
    hMuonMomentumSel->Sumw2();
    hProtonMomentumSel->Sumw2();
    hGoldenPionMomentumSel->Sumw2();
    hNonGoldenPionMomentumSel->Sumw2();

    unsigned int nSignal = 0;
    unsigned int nSignalSelected = 0;

    const auto nEvents = reader.GetNumberOfEvents();
    for (unsigned int i = 0; i < nEvents; ++i)
    {
        AnalysisHelper::PrintLoadingBar(i, nEvents);

        reader.LoadEvent(i);
            
        // Only use signal events
        if (!AnalysisHelper::IsTrueCC1Pi(pEvent))
            continue;

        nSignal++;

        // Check if event passes the CC inclusive selection
        const bool passesCCInclusive = pEvent->reco.passesCCInclusive();

        if (passesCCInclusive)
            nSignalSelected++;

        const auto truthParticles = pEvent->truth.particles;
        for (const auto &particle : truthParticles)
        {
            if (!AnalysisHelper::PassesVisibilityThreshold(particle))
                continue;

            switch (particle.pdgCode())
            {
                case 13:
                    hMuonMomentum->Fill(particle.momentum());

                    if (passesCCInclusive)
                        hMuonMomentumSel->Fill(particle.momentum());
                    break;
                case 2212:
                    hProtonMomentum->Fill(particle.momentum());
                    
                    if (passesCCInclusive)
                        hProtonMomentumSel->Fill(particle.momentum());
                    break;
                case 211:
                    if (AnalysisHelper::IsGolden(particle))
                    {
                        hGoldenPionMomentum->Fill(particle.momentum());
                        
                        if (passesCCInclusive)
                            hGoldenPionMomentumSel->Fill(particle.momentum());
                    }
                    else
                    {
                        hNonGoldenPionMomentum->Fill(particle.momentum());
                        
                        if (passesCCInclusive)
                            hNonGoldenPionMomentumSel->Fill(particle.momentum());
                    }
                    break;
                default:
                    break;
            }
        }
    }

    PlottingHelper::SetLineStyle(hMuonMomentum, PlottingHelper::Muon);
    PlottingHelper::SetLineStyle(hMuonMomentumSel, PlottingHelper::Muon);
    hMuonMomentumSel->SetLineStyle(2);
    
    PlottingHelper::SetLineStyle(hProtonMomentum, PlottingHelper::Proton);
    PlottingHelper::SetLineStyle(hProtonMomentumSel, PlottingHelper::Proton);
    hProtonMomentumSel->SetLineStyle(2);
    
    PlottingHelper::SetLineStyle(hGoldenPionMomentum, PlottingHelper::GoldenPion);
    PlottingHelper::SetLineStyle(hGoldenPionMomentumSel, PlottingHelper::GoldenPion);
    hGoldenPionMomentumSel->SetLineStyle(2);
    
    PlottingHelper::SetLineStyle(hNonGoldenPionMomentum, PlottingHelper::NonGoldenPion);
    PlottingHelper::SetLineStyle(hNonGoldenPionMomentumSel, PlottingHelper::NonGoldenPion);
    hNonGoldenPionMomentumSel->SetLineStyle(2);


    auto pCanvas = PlottingHelper::GetCanvas();

    hMuonMomentum->Draw("hist");
    hMuonMomentumSel->Draw("hist same");
    PlottingHelper::SaveCanvas(pCanvas, "ccinc_muonMomentum");
    
    hProtonMomentum->Draw("hist");
    hProtonMomentumSel->Draw("hist same");
    PlottingHelper::SaveCanvas(pCanvas, "ccinc_protonMomentum");
    
    hGoldenPionMomentum->Draw("hist");
    hGoldenPionMomentumSel->Draw("hist same");
    PlottingHelper::SaveCanvas(pCanvas, "ccinc_goldenPionMomentum");
    
    hNonGoldenPionMomentum->Draw("hist");
    hNonGoldenPionMomentumSel->Draw("hist same");
    PlottingHelper::SaveCanvas(pCanvas, "ccinc_nonGoldenPionMomentum");

    // Ratio plots
    auto pCanvasRatio = PlottingHelper::GetCanvas(960, 270);

    hMuonMomentumSel->SetLineStyle(1);
    hProtonMomentumSel->SetLineStyle(1);
    hGoldenPionMomentumSel->SetLineStyle(1);
    hNonGoldenPionMomentumSel->SetLineStyle(1);

    hMuonMomentumSel->Divide(hMuonMomentum.get());
    hProtonMomentumSel->Divide(hProtonMomentum.get());
    hGoldenPionMomentumSel->Divide(hGoldenPionMomentum.get());
    hNonGoldenPionMomentumSel->Divide(hNonGoldenPionMomentum.get());

    hMuonMomentumSel->Draw("e1");
    PlottingHelper::SaveCanvas(pCanvasRatio, "ccinc_muonMomentum_ratio");
    
    hProtonMomentumSel->Draw("e1");
    PlottingHelper::SaveCanvas(pCanvasRatio, "ccinc_protonMomentum_ratio");
    
    hGoldenPionMomentumSel->Draw("e1");
    PlottingHelper::SaveCanvas(pCanvasRatio, "ccinc_goldenPionMomentum_ratio");
    
    hNonGoldenPionMomentumSel->Draw("e1");
    PlottingHelper::SaveCanvas(pCanvasRatio, "ccinc_nonGoldenPionMomentum_ratio");

    std::cout << "N signal events: " << nSignal << std::endl;
    std::cout << "N signal events selected: " << nSignalSelected << std::endl;

    return 0;
}
