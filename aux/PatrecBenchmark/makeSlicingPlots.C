TTree *GetTree(TFile *file, const std::string &dirName)
{
    TDirectoryFile *dir;
    TTree *tree;
    file->GetObject(dirName.c_str(), dir);
    dir->GetObject("events", tree);

    return tree;
}

void MakePlot(TTree *tree, const std::string &variable, const unsigned int nBins, const float min, const float max, const std::string &fileName)
{
    TCanvas *c = new TCanvas();
    TPad *pad1 = new TPad("pad1","",0,0,1,1);
    TPad *pad2 = new TPad("pad2","",0,0,1,1);

    // Make Pad2 transparent
    pad2->SetFillStyle(4000);
    pad2->SetFrameFillStyle(0);
    gStyle->SetOptStat(0);

    TH1F *h = new TH1F("h", "", nBins, min, max);
    TH1F *hSel = new TH1F("hSel", "", nBins, min, max);
    TH1F *hRat = new TH1F("hRat", "", nBins, min, max);

    h->SetLineColor(kRed);
    h->SetLineWidth(2);
    h->GetXaxis()->SetTitle(variable.c_str());
    h->GetYaxis()->SetTitle("Fraction of events");
    h->GetYaxis()->SetTitleOffset(1.3);

    hRat->SetLineColor(kBlue);
    hRat->SetLineWidth(2);
    hRat->GetYaxis()->SetTitle("Selection efficiency");
    hRat->GetYaxis()->SetTitleOffset(1.3);

    tree->Draw((variable + " >> h").c_str(), "hasSlice && hasNuHits");
    tree->Draw((variable + " >> hSel").c_str(), "hasSlice && hasNuHits && isMCSSelected");
    hRat->Add(hSel);
    hRat->Sumw2();
    hRat->Divide(h);
    h->Scale(1.f / h->GetEntries());

    pad1->Draw();
    pad1->cd();
    hRat->GetYaxis()->SetRangeUser(0, 1);
    hRat->DrawCopy("hist");
    hRat->SetFillColor(kBlue);
    hRat->SetFillStyle(3002);
    hRat->Draw("E2 same");

    pad2->Draw();
    pad2->cd();
    h->DrawCopy("hist Y+");
    h->SetFillColor(kRed);
    h->SetFillStyle(3002);
    h->Draw("E2 same");

    c->SaveAs(fileName.c_str());

    h->Delete();
    hSel->Delete();
    hRat->Delete();
}

void makeSlicingPlots()
{
    TFile *file = TFile::Open("/uboone/data/users/asmith/ubcc1pi/09102019/patrecBenchmarkStudy.root");
    const auto treeDefault = GetTree(file, "patrecBenchmarkStudyDefault");

    MakePlot(treeDefault, "mcsPurity", 50, 0, 1, "slicePurity.png");
    MakePlot(treeDefault, "mcsCompleteness", 50, 0, 1, "sliceCompleteness.png");
    MakePlot(treeDefault, "nuEnergy", 50, 0.3, 2.5, "nuEnergy.png");
    MakePlot(treeDefault, "nNuHitsTotal", 50, 0, 3000, "nNuHitsTotal.png");
}
