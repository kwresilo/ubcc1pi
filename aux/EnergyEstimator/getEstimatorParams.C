// Instructions
//
// Compile this code in root
// >> root -l
// >> .L getEstimatorParams.C
//
// Make a default parameters object and print it's details
// >> ubcc1pi::Parameters p
// >> ubcc1pi::PrintParameters(p)
//
// Modify any of the parameters you wish, e.g.
// >> p.m_inputFile = "/path/to/my/input/file.root"
//
// Actually make the range tables
// >> ubcc1pi::GetEstimatorParams(p)

namespace ubcc1pi
{

struct Parameters
{
    std::string m_inputFile = "energyEstimator.root"; ///< The input file name
    std::string m_directory = "energyEstimator";      ///< The TDirectoryFile into which we have to `cd` to access the tree
    std::string m_treeName  = "particles";            ///< The name of the tree containing the particles

    std::vector<int> m_pdgCodes = {13, 211};    ///< The pdg codes for which to do the fit

    std::map<int, float> m_pdgToRangeBinWidth = {
        {13,   1.f},
        {211,  1.f}
    };                                                ///< The bin width to use when sampling ranges (cm)

    std::map<int, float> m_pdgToMinRange = {
        {13,   0.f},
        {211,  0.f}
    };                                                ///< The minimum range to consider in the fit (cm)
    
    std::map<int, float> m_pdgToMaxRange = {
        {13,   400.f},
        {211,  200.f}
    };                                                ///< The maximum range to consider in the fit (cm)
};

// -----------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
void PrintParameter(const std::string &type, const std::string &name, const T &value)
{
    std::cout << std::setw(16) << type << std::setw(48) << name << " = " << value << std::endl;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void ValidateParameters(const Parameters &p, const bool shouldPrint)
{
    if (shouldPrint)
    {
        std::cout << "Input parameters:" << std::endl;
        PrintParameter("string", "m_inputFile", p.m_inputFile);
        PrintParameter("string", "m_directory", p.m_directory);
        PrintParameter("string", "m_treeName", p.m_treeName);

        // Print the PDG codes
        std::string pdgCodeString;
        for (const auto &pdgCode : p.m_pdgCodes)
            pdgCodeString += std::to_string(pdgCode) + " ";

        PrintParameter("vector<int>", "m_pdgCodes", pdgCodeString); 
    }

    // Print the map values for each pdg code
    for (const auto &pdgCode : p.m_pdgCodes)
    {
        // Range bin width
        const auto rangeBinWidthIter = p.m_pdgToRangeBinWidth.find(pdgCode);
        const auto rangeBinWidthString = "m_pdgToRangeBinWidth.at(" + std::to_string(pdgCode) + ")";

        if (rangeBinWidthIter == p.m_pdgToRangeBinWidth.end())
            throw std::invalid_argument(rangeBinWidthString + " not set");

        if (shouldPrint)
            PrintParameter("float", rangeBinWidthString, rangeBinWidthIter->second);
        
        // Min range
        const auto minRangeIter = p.m_pdgToMinRange.find(pdgCode);
        const auto minRangeString = "m_pdgToMinRange.at(" + std::to_string(pdgCode) + ")";

        if (minRangeIter == p.m_pdgToMinRange.end())
            throw std::invalid_argument(minRangeString + " not set");

        if (shouldPrint)
            PrintParameter("float", minRangeString, minRangeIter->second);
        
        // Max range
        const auto maxRangeIter = p.m_pdgToMaxRange.find(pdgCode);
        const auto maxRangeString = "m_pdgToMaxRange.at(" + std::to_string(pdgCode) + ")";

        if (maxRangeIter == p.m_pdgToMaxRange.end())
            throw std::invalid_argument(maxRangeString + " not set");

        if (shouldPrint)
            PrintParameter("float", maxRangeString, maxRangeIter->second);
    }
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void PrintParameters(const Parameters &p)
{
    const auto shouldPrint = true;
    ValidateParameters(p, shouldPrint);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void FitData(const Parameters &p, const unsigned int pdg, const std::vector<float> &trueRanges, const std::vector<float> &trueMomenta, TF1 *&pFit, TGraph *&pGraph)
{
    const auto N = trueRanges.size();
    if (trueMomenta.size() != N)
        throw std::invalid_argument("Input vectors of true ranges and momenta don't have the same size");

    // Convert the vectors to arrays
    const auto &rangesArray = trueRanges.data();
    const auto &momentaArray = trueMomenta.data();

    // Set up the fit
    const auto functionName = ("func" + std::to_string(pdg)).c_str();
    const auto minRange = p.m_pdgToMinRange.at(pdg);
    const auto maxRange = p.m_pdgToMaxRange.at(pdg);

    pGraph = new TGraph(N, rangesArray, momentaArray);
    pFit = new TF1(functionName, "[0]*pow(x, [1]) + [2]*pow(x, [3])", minRange, maxRange);

    // Choose some initial parameters to help the fit converge quickly
    pFit->SetParameter(0, 0.0003);
    pFit->SetParameter(1, 1.2);
    pFit->SetParameter(2, 0.01);
    pFit->SetParameter(3, 0.5);

    // Do the fit!
    pGraph->Fit(pFit);

    // Print the fit to screen
    const auto fit0 = pFit->GetParameter(0);
    const auto fit1 = pFit->GetParameter(1);
    const auto fit2 = pFit->GetParameter(2);
    const auto fit3 = pFit->GetParameter(3);
    std::cout << "Final fit parameters for PDG: " << pdg << std::endl;
    std::cout << "momentum = " << fit0 << " * std::pow(range, " << fit1 << ") + " << fit2 << " * std::pow(range, " << fit3 << ")" << std::endl;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void GetEstimatorParams(const Parameters &p)
{
    // Print the parameters and raise an exeption if they are invalid
    const bool shouldPrintInputParams = false;
    ValidateParameters(p, shouldPrintInputParams);

    // Access the tree
    TFile *pFile = TFile::Open(p.m_inputFile.c_str());
    TDirectoryFile *pDir = static_cast<TDirectoryFile *>(pFile->Get(p.m_directory.c_str()));
    TTree *pTree = static_cast<TTree *>(pDir->Get(p.m_treeName.c_str()));

    // Set up the branches
    bool isSignal;
    pTree->SetBranchAddress("isSignal", &isSignal);

    bool hasMatchedMCParticle;
    pTree->SetBranchAddress("hasMatchedMCParticle", &hasMatchedMCParticle);
    
    float trueMatchCompleteness;
    pTree->SetBranchAddress("trueMatchCompleteness", &trueMatchCompleteness);
    
    int truePdgCode;
    pTree->SetBranchAddress("truePdgCode", &truePdgCode);

    bool isGolden;
    pTree->SetBranchAddress("isGolden", &isGolden);
    
    float trueMomentum;
    pTree->SetBranchAddress("trueMomentum", &trueMomentum);
    
    float trueRange;
    pTree->SetBranchAddress("trueRange", &trueRange);
    
    bool isPrimary;
    pTree->SetBranchAddress("isPrimary", &isPrimary);
    
    // Get the true range and momenta of the particles we want to fit
    std::map<int, std::vector<float> > pdgToTrueRangesMap, pdgToTrueMomentaMap;
    for (unsigned int i = 0; i < pTree->GetEntries(); ++i)
    {
        pTree->GetEntry(i);

        // Ensure we are looking at CC1Pi events
        if (!isSignal)
            continue;

        // Ensure we only have primary PFParticles that match to golden MCParticles
        if (!hasMatchedMCParticle || trueMatchCompleteness <= 0.5f || !isGolden || !isPrimary)
            continue;
        
        // Ensure that this is a PDG code we care about
        if (std::find(p.m_pdgCodes.begin(), p.m_pdgCodes.end(), truePdgCode) == p.m_pdgCodes.end())
            continue;

        // Apply the range thresholds to make the plot
        if (trueRange < p.m_pdgToMinRange.at(truePdgCode) || trueRange > p.m_pdgToMaxRange.at(truePdgCode))
            continue;

        pdgToTrueRangesMap[truePdgCode].push_back(trueRange);
        pdgToTrueMomentaMap[truePdgCode].push_back(trueMomentum);
    }

    // Do the fit for each PDG code and make the plots
    TCanvas *c = new TCanvas();
    gStyle->SetOptStat(0);

    for (const auto &pdgCode : p.m_pdgCodes)
    {
        TF1 *pFit = nullptr;
        TGraph *pGraph = nullptr;

        const auto &trueRanges = pdgToTrueRangesMap.at(pdgCode);
        const auto &trueMomenta = pdgToTrueMomentaMap.at(pdgCode);

        FitData(p, pdgCode, trueRanges, trueMomenta, pFit, pGraph);
    
        pGraph->Draw("AP");
        pFit->SetLineColor(kOrange - 3);
        pFit->Draw("same");

        c->SaveAs(("rangeTable_" + std::to_string(pdgCode) + ".png").c_str());
        c->SaveAs(("rangeTable_" + std::to_string(pdgCode) + ".pdf").c_str());
        c->SaveAs(("rangeTable_" + std::to_string(pdgCode) + ".C").c_str());
    }
}

}
