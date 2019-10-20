TTree *GetTree(TFile *file, const std::string &dirName)
{
    TDirectoryFile *dir;
    TTree *tree;
    file->GetObject(dirName.c_str(), dir);
    dir->GetObject("events", tree);

    return tree;
}

void makePlots()
{
    TFile *file = TFile::Open("/uboone/data/users/asmith/ubcc1pi/09102019/patrecBenchmarkStudy.root");

    const auto treeDefault = GetTree(file, "patrecBenchmarkStudyDefault");
    const auto treeCut1 = GetTree(file, "patrecBenchmarkStudy1");
    const auto treeCut2 = GetTree(file, "patrecBenchmarkStudy2");
    const auto treeCut5 = GetTree(file, "patrecBenchmarkStudy5");
    const auto treeCut10 = GetTree(file, "patrecBenchmarkStudy10");
    const auto treeCut20 = GetTree(file, "patrecBenchmarkStudy20");

    TH1F *hDefaultMatches = new TH1F("hDefaultMatches", "", 6, 0, 6);
    TH1F *hCutMatches = new TH1F("hCutMatches", "", 6, 0, 6);
    
    hDefaultMatches->SetLineWidth(2);
    hDefaultMatches->SetLineColor(kBlue);
    hCutMatches->SetLineWidth(2);
    hCutMatches->SetLineColor(kOrange);
    
    TH1F *hDefaultPurity = new TH1F("hDefaultPurity", "", 100, -0.05, 1.05);
    TH1F *hCutPurity = new TH1F("hCutPurity", "", 100, -0.05, 1.05);
    
    hDefaultPurity->SetLineWidth(2);
    hDefaultPurity->SetLineColor(kBlue);
    hCutPurity->SetLineWidth(2);
    hCutPurity->SetLineColor(kOrange);
    
    TH1F *hDefaultCompleteness = new TH1F("hDefaultCompleteness", "", 100, -0.05, 1.05);
    TH1F *hCutCompleteness = new TH1F("hCutCompleteness", "", 100, -0.05, 1.05);
    
    hDefaultCompleteness->SetLineWidth(2);
    hDefaultCompleteness->SetLineColor(kBlue);
    hCutCompleteness->SetLineWidth(2);
    hCutCompleteness->SetLineColor(kOrange);

    const auto trees = std::vector<TTree *>({treeCut1, treeCut2, treeCut5, treeCut10, treeCut20});
    const auto cuts = std::vector<float>({1, 2, 5, 10, 20});

    const auto pdgs = std::vector<int>({13, 211, 2212});
    const auto maxs = std::vector<float>({20000, 17000, 45000});

    TCanvas *c = new TCanvas();
    for (unsigned int j = 0; j < pdgs.size(); ++j)
    {
        const auto pdg = pdgs.at(j);
        const auto max = maxs.at(j);

        for (unsigned int i = 0; i < trees.size(); ++i)
        {
            const auto tree = trees.at(i);
            const auto cut = cuts.at(i);

            // Number of matches
            treeDefault->Draw("mcp_nPFPMatches >> hDefaultMatches", ("isMCSSelected && mcp_pdg == " + std::to_string(pdg)).c_str());
            tree->Draw("mcp_nPFPMatches >> hCutMatches", ("isMCSSelected && mcp_pdg == " + std::to_string(pdg)).c_str());

            hDefaultMatches->GetYaxis()->SetRangeUser(0, max);
            hCutMatches->GetYaxis()->SetRangeUser(0, max);

            hDefaultMatches->Draw();
            hCutMatches->Draw("same");

            c->SaveAs(("nPFPMatches_default-vs-cut" + std::to_string(cut) + "_" + std::to_string(pdg) + ".png").c_str());
            
            // Purity
            treeDefault->Draw("mcp_bestMatchPurity >> hDefaultPurity", ("isMCSSelected && mcp_hasMatchedPFP && mcp_hitWeight > 0 && mcp_pdg == " + std::to_string(pdg)).c_str());
            tree->Draw("mcp_bestMatchPurity >> hCutPurity", ("isMCSSelected && mcp_hasMatchedPFP && mcp_hitWeight > 0 && mcp_pdg == " + std::to_string(pdg)).c_str());

            hDefaultPurity->Scale(1.f / hDefaultPurity->GetEntries());
            hCutPurity->Scale(1.f / hCutPurity->GetEntries());
            hDefaultPurity->GetYaxis()->SetRangeUser(0, 0.35);
            hCutPurity->GetYaxis()->SetRangeUser(0, 0.35);

            hDefaultPurity->Draw("hist");
            hCutPurity->Draw("hist same");

            c->SaveAs(("purity_default-vs-cut" + std::to_string(cut) + "_" + std::to_string(pdg) + ".png").c_str());
            
            // Completeness
            treeDefault->Draw("mcp_bestMatchCompleteness >> hDefaultCompleteness", ("isMCSSelected && mcp_hasMatchedPFP && mcp_hitWeight > 0 && mcp_pdg == " + std::to_string(pdg)).c_str());
            tree->Draw("mcp_bestMatchCompleteness >> hCutCompleteness", ("isMCSSelected && mcp_hasMatchedPFP && mcp_hitWeight > 0 && mcp_pdg == " + std::to_string(pdg)).c_str());

            hDefaultCompleteness->Scale(1.f / hDefaultCompleteness->GetEntries());
            hCutCompleteness->Scale(1.f / hCutCompleteness->GetEntries());
            hDefaultCompleteness->GetYaxis()->SetRangeUser(0, 0.35);
            hCutCompleteness->GetYaxis()->SetRangeUser(0, 0.35);

            hDefaultCompleteness->Draw("hist");
            hCutCompleteness->Draw("hist same");

            c->SaveAs(("completeness_default-vs-cut" + std::to_string(cut) + "_" + std::to_string(pdg) + ".png").c_str());
        }
    }
}
