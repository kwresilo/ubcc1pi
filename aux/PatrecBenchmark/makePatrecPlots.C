TTree *GetTree(TFile *file, const std::string &dirName)
{
    TDirectoryFile *dir;
    TTree *tree;
    file->GetObject(dirName.c_str(), dir);
    dir->GetObject("events", tree);

    return tree;
}

void MakeNMatchesPlot(TTree *tree, const std::string &variable, const unsigned int nBins, const float min, const float max, const int pdg, const std::string &fileName)
{
    TCanvas *c = new TCanvas();
    gStyle->SetOptStat(0);

    THStack *hs = new THStack("hs", "");
    TH1F *h0 = new TH1F("h0", "", nBins, min, max);
    TH1F *h1 = new TH1F("h1", "", nBins, min, max);
    TH1F *h2 = new TH1F("h2", "", nBins, min, max);
    TH1F *h3plus = new TH1F("h3plus", "", nBins, min, max);

    h0->SetFillColor(kRed);
    h1->SetFillColor(kBlue);
    h2->SetFillColor(kOrange);
    h3plus->SetFillColor(kGreen);

    const auto n0 = tree->Draw((variable + " >> h0").c_str(), ("isMCSSelected && mcp_nPFPMatches == 0 && mcp_pdg == " + std::to_string(pdg)).c_str());
    const auto n1 = tree->Draw((variable + " >> h1").c_str(), ("isMCSSelected && mcp_nPFPMatches == 1 && mcp_pdg == " + std::to_string(pdg)).c_str());
    const auto n2 = tree->Draw((variable + " >> h2").c_str(), ("isMCSSelected && mcp_nPFPMatches == 2 && mcp_pdg == " + std::to_string(pdg)).c_str());
    const auto n3plus = tree->Draw((variable + " >> h3plus").c_str(), ("isMCSSelected && mcp_nPFPMatches >= 3 && mcp_pdg == " + std::to_string(pdg)).c_str());

    std::cout << "For PDG = " << pdg << std::endl;
    std::cout << "  - 0 match    : " << n0 << std::endl;
    std::cout << "  - 1 match    : " << n1 << std::endl;
    std::cout << "  - 2 matches  : " << n2 << std::endl;
    std::cout << "  - 3+ matches : " << n3plus << std::endl;

    hs->Add(h0);
    hs->Add(h1);
    hs->Add(h2);
    hs->Add(h3plus);

    hs->Draw();
    c->SaveAs(fileName.c_str());

    hs->Delete();
    h0->Delete();
    h1->Delete();
    h2->Delete();
    h3plus->Delete();
}

void makePatrecPlots()
{
    TFile *file = TFile::Open("/uboone/data/users/asmith/ubcc1pi/09102019/patrecBenchmarkStudy.root");
    const auto treeDefault = GetTree(file, "patrecBenchmarkStudyDefault");

    MakeNMatchesPlot(treeDefault, "mcp_momentum", 60, 0, 2, 13, "nPFPMatches-vs-momentum_13.png");
    MakeNMatchesPlot(treeDefault, "mcp_momentum", 60, 0, 2, 211, "nPFPMatches-vs-momentum_211.png");
    MakeNMatchesPlot(treeDefault, "mcp_momentum", 60, 0, 2, 2212, "nPFPMatches-vs-momentum_2212.png");
    MakeNMatchesPlot(treeDefault, "mcp_hitWeight", 60, 0, 3000, 13, "nPFPMatches-vs-hitWeight_13.png");
    MakeNMatchesPlot(treeDefault, "mcp_hitWeight", 60, 0, 2000, 211, "nPFPMatches-vs-hitWeight_211.png");
    MakeNMatchesPlot(treeDefault, "mcp_hitWeight", 60, 0, 600, 2212, "nPFPMatches-vs-hitWeight_2212.png");
}
