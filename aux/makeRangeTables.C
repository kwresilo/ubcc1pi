// Instructions
//
// Compile this code in root
// >> root -l
// >> .L makeRangeTables.C
//
// Make a default parameters object and print it's details
// >> rtable::Parameters p
// >> rtable::PrintParameters(p)
//
// Modify any of the parameters you wish, e.g.
// >> p.m_inputFile = "/path/to/my/input/file.root"
//
// Actually make the range tables
// >> rtable::MakeRangeTables(p)

namespace rtable
{

struct Parameters
{
    std::string m_inputFile = "energyEstimator.root"; ///< The input file name
    std::string m_directory = "energyEstimator";      ///< The TDirectoryFile into which we have to `cd` to access the tree
    std::string m_treeName  = "particles";            ///< The name of the tree containing the particles

    float       m_muonRangeBinWidth = 1.f;            ///< The bin width to use when sampling muon ranges (cm)
    float       m_pionRangeBinWidth = 1.f;            ///< The bin width to use when sampling pion ranges (cm)
    float       m_fitRangeCut = 200.f;                ///< The range below which pions will be considered when making the fit
    float       m_maxMuonRange = 400.f;               ///< The maximum muon range to plot
    float       m_maxPionRange = 200.f;               ///< The maximum pion range to plot
};

// -----------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
void PrintParameter(const std::string &type, const std::string &name, const T &value)
{
    std::cout << std::setw(16) << type << std::setw(24) << name << " = " << value << std::endl;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void PrintParameters(const Parameters &p)
{
    std::cout << "Range table parameters" << std::endl;
    PrintParameter("string", "m_inputFile", p.m_inputFile);
    PrintParameter("string", "m_directory", p.m_directory);
    PrintParameter("string", "m_treeName", p.m_treeName);
    PrintParameter("float", "m_muonRangeBinWidth", p.m_muonRangeBinWidth);
    PrintParameter("float", "m_pionRangeBinWidth", p.m_pionRangeBinWidth);
    PrintParameter("float", "m_fitRangeCut", p.m_fitRangeCut);
    PrintParameter("float", "m_maxMuonRange", p.m_maxMuonRange);
    PrintParameter("float", "m_maxPionRange", p.m_maxPionRange);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void FitData(const unsigned int pdg, const float fitRangeCut, const float maxRange, const std::vector<float> &ranges, const std::vector<float> &KEs, TF1 *&pFit, TGraph *&pGraph)
{
    // Apply the cut to the fitted point
    std::vector<float> cutRanges, cutKEs;
    for (unsigned int i = 0; i < ranges.size(); ++i)
    {
        if (ranges.at(i) > fitRangeCut)
            continue;

        cutRanges.push_back(ranges.at(i));
        cutKEs.push_back(KEs.at(i));
    }

    // Convert the vectors to arrays
    const unsigned int N = ranges.size();
    float rangesArray[N];
    float KEsArray[N];
    for (unsigned int i = 0; i < N; ++i)
    {
        const auto range = ranges.at(i);
        rangesArray[i] = range;
        KEsArray[i] = KEs.at(i);
    }
    
    const unsigned int cutN = cutRanges.size();
    float cutRangesArray[cutN];
    float cutKEsArray[cutN];
    for (unsigned int i = 0; i < cutN; ++i)
    {
        cutRangesArray[i] = cutRanges.at(i);
        cutKEsArray[i] = cutKEs.at(i);
    }

    pGraph = new TGraph(N, rangesArray, KEsArray);
    TGraph *pCutGraph = new TGraph(cutN, cutRangesArray, cutKEsArray);

    pFit = new TF1(("func" + std::to_string(pdg)).c_str(), "[0]*pow(x, [1]) + [2]*pow(x, [3])", 0.f, maxRange);
    pFit->SetParameter(0, 0.0003);
    pFit->SetParameter(1, 1.2);
    pFit->SetParameter(2, 0.01);
    pFit->SetParameter(3, 0.5);
    pCutGraph->Fit(pFit);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void GetKEsFromFit(const std::vector<float> &ranges, TF1 *pFit, std::vector<float> &KEsFromFit)
{
    for (unsigned int i = 0; i < ranges.size(); ++i)
    {
        const float range = ranges.at(i);
        const float KE = pFit->Eval(range);
        KEsFromFit.push_back(KE);
    }
    
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void MakePlot(TCanvas *c, const std::vector<float> &x, const float xMin, const float xMax, const float dx, const std::string &name)
{
    const unsigned int nBinsX = std::floor((xMax - xMin) / dx);

    TH1F *h = new TH1F(name.c_str(), name.c_str(), nBinsX, xMin, xMax);

    for (unsigned int i = 0; i < x.size(); ++i)
    {
        h->Fill(x.at(i));
    }

    h->Draw();
    c->SaveAs((name + ".png").c_str());
}

// -----------------------------------------------------------------------------------------------------------------------------------------

    
void MakePlot(TCanvas *c, const std::vector<float> &x, const float xMin, const float xMax, const float dx, const std::vector<float> &y, const float yMin, const float yMax, const float dy, const std::string &name)
{
    const unsigned int nBinsX = std::floor((xMax - xMin) / dx);
    const unsigned int nBinsY = std::floor((yMax - yMin) / dy);

    TH2F *h = new TH2F(name.c_str(), name.c_str(), nBinsX, xMin, xMax, nBinsY, yMin, yMax);

    if (x.size() != y.size())
    {
        std::cerr << "Make plot error. x & y input data sizes don't match" << std::endl;
        exit(2);
    }

    for (unsigned int i = 0; i < x.size(); ++i)
    {
        h->Fill(x.at(i), y.at(i));
    }

    h->Draw("colz");
    c->SaveAs((name + ".png").c_str());
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void GetRatios(const std::vector<float> &a, const std::vector<float> &b, std::vector<float> &ratio)
{
    if (a.size() != b.size())
    {
        std::cerr << "Get ratios error. a & b input data sizes don't match" << std::endl;
        exit(3);
    }

    for (unsigned int i = 0; i < a.size(); ++i)
    {
        ratio.push_back((a.at(i) - b.at(i)) / b.at(i));
    }
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void MakeRangeTables(const Parameters &p)
{
    // Access the tree
    TFile *pFile = TFile::Open(p.m_inputFile.c_str());
    TDirectoryFile *pDir = static_cast<TDirectoryFile *>(pFile->Get(p.m_directory.c_str()));
    TTree *pTree = static_cast<TTree *>(pDir->Get(p.m_treeName.c_str()));

    // Set up the branches
    bool hasMatchedMCParticle;
    pTree->SetBranchAddress("hasMatchedMCParticle", &hasMatchedMCParticle);
    
    float trueMatchCompleteness;
    pTree->SetBranchAddress("trueMatchCompleteness", &trueMatchCompleteness);
    
    int truePdgCode;
    pTree->SetBranchAddress("truePdgCode", &truePdgCode);

    bool isGolden;
    pTree->SetBranchAddress("isGolden", &isGolden);

    float trueKE;
    pTree->SetBranchAddress("trueKE", &trueKE);
    
    float trueRange;
    pTree->SetBranchAddress("trueRange", &trueRange);
    
    float range;
    pTree->SetBranchAddress("range", &range);


    // Extract the KEs and ranges
    std::vector<float> muonRanges, pionRanges, muonKEs, pionKEs, muonRecoRanges, pionRecoRanges;
    for (unsigned int i = 0; i < pTree->GetEntries(); ++i)
    {
        pTree->GetEntry(i);

        // Ensure we only have PFParticles that match to golden MCParticles
        if (!hasMatchedMCParticle || trueMatchCompleteness <= 0.5f || !isGolden)
            continue;

        if (truePdgCode == 13)
        {
            if (trueRange > p.m_maxMuonRange)
                continue;

            muonRanges.push_back(trueRange);
            muonKEs.push_back(trueKE);
            muonRecoRanges.push_back(range);
        }
        
        if (truePdgCode == 211)
        {
            if (trueRange > p.m_maxPionRange)
                continue;

            pionRanges.push_back(trueRange);
            pionKEs.push_back(trueKE);
            pionRecoRanges.push_back(range);
        }
    }


    // Do the fit and make the plots
    TF1 *pMuonFit = nullptr;
    TGraph *pMuonGraph = nullptr;
    FitData(13, p.m_fitRangeCut, p.m_maxMuonRange, muonRanges, muonKEs, pMuonFit, pMuonGraph);
    
    TF1 *pPionFit = nullptr;
    TGraph *pPionGraph = nullptr;
    FitData(211, p.m_fitRangeCut, p.m_maxPionRange, pionRanges, pionKEs, pPionFit, pPionGraph);

    TCanvas *c = new TCanvas();
    gStyle->SetOptStat(0);

    pMuonGraph->Draw("AP");
    pMuonFit->SetLineColor(kRed);
    pMuonFit->Draw("same");
    c->SaveAs("muonRangeTable.png");
    
    pPionGraph->Draw("AP");
    pPionFit->SetLineColor(kRed);
    pPionFit->Draw("same");
    c->SaveAs("pionRangeTable.png");

    // Use these functions to go from reco range -> KE
    std::vector<float> muonKEsFromFit;
    GetKEsFromFit(muonRanges, pMuonFit, muonKEsFromFit);
    
    std::vector<float> muonRecoKEsFromFit;
    GetKEsFromFit(muonRecoRanges, pMuonFit, muonRecoKEsFromFit);
    
    std::vector<float> pionKEsFromFit;
    GetKEsFromFit(pionRanges, pPionFit, pionKEsFromFit);
    
    std::vector<float> pionRecoKEsFromFit;
    GetKEsFromFit(pionRecoRanges, pPionFit, pionRecoKEsFromFit);

    // Get the ratios of some of these variables: (a-b) / b
    std::vector<float> muonRangeRatio;
    GetRatios(muonRecoRanges, muonRanges, muonRangeRatio);
    
    std::vector<float> muonKEsFromFitRatio;
    GetRatios(muonKEsFromFit, muonKEs, muonKEsFromFitRatio);
    
    std::vector<float> muonRecoKEsFromFitRatio;
    GetRatios(muonRecoKEsFromFit, muonKEs, muonRecoKEsFromFitRatio);
    
    std::vector<float> pionRangeRatio;
    GetRatios(pionRecoRanges, pionRanges, pionRangeRatio);
    
    std::vector<float> pionKEsFromFitRatio;
    GetRatios(pionKEsFromFit, pionKEs, pionKEsFromFitRatio);
    
    std::vector<float> pionRecoKEsFromFitRatio;
    GetRatios(pionRecoKEsFromFit, pionKEs, pionRecoKEsFromFitRatio);

    // Make the relevant plots
    MakePlot(c, muonRanges, 0.f, p.m_maxMuonRange, p.m_muonRangeBinWidth, muonKEs, 0.f, 0.9f, 0.005f, "muon_trueKE-vs-trueRange");
    MakePlot(c, muonRanges, 0.f, p.m_maxMuonRange, p.m_muonRangeBinWidth, muonKEsFromFitRatio, -1.f, 1.f, 0.005f, "muon_trueKE-to-trueFitKE-vs-trueRange");
    MakePlot(c, muonRanges, 0.f, p.m_maxMuonRange, p.m_muonRangeBinWidth, muonRecoRanges, 0.f, p.m_maxMuonRange, p.m_muonRangeBinWidth, "muon_recoRange-vs-trueRange");
    MakePlot(c, muonRanges, 0.f, p.m_maxMuonRange, p.m_muonRangeBinWidth, muonRangeRatio, -1.f, 1.f, 0.005f, "muon_rangeRatio-vs-trueRange");
    MakePlot(c, muonRanges, 0.f, p.m_maxMuonRange, p.m_muonRangeBinWidth, muonRecoKEsFromFitRatio, -2.f, 2.f, 0.005f, "muon_trueKE-to-recoFitKE-vs-trueRange");
    
    MakePlot(c, pionRanges, 0.f, p.m_maxPionRange, p.m_pionRangeBinWidth, pionKEs, 0.f, 0.5f, 0.005f, "pion_trueKE-vs-trueRange");
    MakePlot(c, pionRanges, 0.f, p.m_maxPionRange, p.m_pionRangeBinWidth, pionKEsFromFitRatio, -1.f, 1.f, 0.005f, "pion_trueKE-to-trueFitKE-vs-trueRange");
    MakePlot(c, pionRanges, 0.f, p.m_maxPionRange, p.m_pionRangeBinWidth, pionRecoRanges, 0.f, p.m_maxPionRange, p.m_pionRangeBinWidth, "pion_recoRange-vs-trueRange");
    MakePlot(c, pionRanges, 0.f, p.m_maxPionRange, p.m_pionRangeBinWidth, pionRangeRatio, -1.f, 1.f, 0.01f, "pion_rangeRatio-vs-trueRange");
    MakePlot(c, pionRanges, 0.f, p.m_maxPionRange, p.m_pionRangeBinWidth, pionRecoKEsFromFitRatio, -2.f, 2.f, 0.01f, "pion_trueKE-to-recoFitKE-vs-trueRange");
    
    gStyle->SetOptStat(1);
    MakePlot(c, muonKEsFromFitRatio, -0.2f, 0.2f, 0.005f, "muon_trueKE-to-trueFitKE");
    MakePlot(c, pionKEsFromFitRatio, -0.2f, 0.2f, 0.005f, "pion_trueKE-to-trueFitKE");
}


}
