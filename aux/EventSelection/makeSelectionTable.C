std::string GetInteractionType(const int nMuMinus,
                               const int nMuPlus,
                               const int nPiPlus,
                               const int nPiMinus,
                               const int nKPlus,
                               const int nKMinus,
                               const int nProton,
                               const int nNeutron,
                               const int nPhoton,
                               const int nElectron,
                               const int nPositron,
                               const int nTotal,
                               const bool isTrueNuFiducial,
                               const bool isSignal)
{
    const int nOther = nTotal - nMuMinus - nMuPlus - nPiPlus - nPiMinus - nKPlus - nKMinus - nProton - nNeutron - nPhoton - nElectron - nPositron;
    std::string interaction;
    
    if (!isTrueNuFiducial)
        return "non-fiducial  ";
 
    if (isSignal)
        interaction += "signal  ";

    if (nMuMinus != 0)
        interaction += std::to_string(nMuMinus) + " Mu-  ";

    if (nMuPlus != 0)
        interaction += std::to_string(nMuPlus) + " Mu+  ";

    if (nPiPlus != 0)
        interaction += std::to_string(nPiPlus) + " Pi+  ";

    if (nPiMinus != 0)
        interaction += std::to_string(nPiMinus) + " Pi-  ";

    if (nKPlus != 0)
        interaction += std::to_string(nKPlus) + " K+  ";

    if (nKMinus != 0)
        interaction += std::to_string(nKMinus) + " K-  ";

    if (nProton != 0)
        interaction += std::to_string(nProton) + " p  ";

    if (nNeutron != 0)
        interaction += std::to_string(nNeutron) + " n  ";

    if (nPhoton != 0)
        interaction += std::to_string(nPhoton) + " gamma  ";

    if (nElectron != 0)
        interaction += std::to_string(nElectron) + " e-  ";

    if (nPositron != 0)
        interaction += std::to_string(nPositron) + " e+  ";
    
    if (nOther != 0)
        interaction += std::to_string(nOther) + " ?  ";

    return interaction;
}
    

bool IsNear(const TVector3 &posA, const TVector3 &posB, const float cut)
{
    return ((posA - posB).Mag() < cut);
}

void makeSelectionTable()
{
    // Get the input tree
    TFile *file = TFile::Open("/uboone/data/users/asmith/ubcc1pi/19062019/eventSelection.root");    
    TDirectoryFile *dir;
    TTree *tree;
    file->GetObject("eventSelection", dir);
    dir->GetObject("events", tree);
    std::cout << "Opened input tree" << std::endl;

    // True event level info
    int run = -std::numeric_limits<int>::max();
    int subRun = -std::numeric_limits<int>::max();
    int event = -std::numeric_limits<int>::max();
    bool isTrueNuFiducial = false;
    bool isSignal = false;
    float trueNuE = -std::numeric_limits<float>::max();
    int nMuMinus = -std::numeric_limits<int>::max();
    int nMuPlus = -std::numeric_limits<int>::max();
    int nPiPlus = -std::numeric_limits<int>::max();
    int nPiMinus = -std::numeric_limits<int>::max();
    int nKPlus = -std::numeric_limits<int>::max();
    int nKMinus = -std::numeric_limits<int>::max();
    int nProton = -std::numeric_limits<int>::max();
    int nNeutron = -std::numeric_limits<int>::max();
    int nPhoton = -std::numeric_limits<int>::max();
    int nElectron = -std::numeric_limits<int>::max();
    int nPositron = -std::numeric_limits<int>::max();
    int nTotal = -std::numeric_limits<int>::max();

    // Reco event level info
    float topologicalScore = -std::numeric_limits<float>::max();
    bool isRecoNuFiducial = false;
    int nFinalStatePFPs = -std::numeric_limits<int>::max();
    TVector3 *recoNuVtx = nullptr;

    // Track info
    std::vector<bool> *hasTrackInfoVect = nullptr;
    std::vector<bool> *isContainedVect = nullptr;
    std::vector<float> *startXVect = nullptr;
    std::vector<float> *startYVect = nullptr;
    std::vector<float> *startZVect = nullptr;
    std::vector<float> *endXVect = nullptr;
    std::vector<float> *endYVect = nullptr;
    std::vector<float> *endZVect = nullptr;
    std::vector<float> *directionXVect = nullptr;
    std::vector<float> *directionYVect = nullptr;
    std::vector<float> *directionZVect = nullptr;
    std::vector<float> *trackShowerVect = nullptr;
    std::vector<float> *chi2pWVect = nullptr;
    std::vector<float> *chi2pUVVect = nullptr;

    // Set up the branches
    tree->SetBranchAddress("run", &run);
    tree->SetBranchAddress("subRun", &subRun);
    tree->SetBranchAddress("event", &event);
    tree->SetBranchAddress("isTrueNuFiducial", &isTrueNuFiducial);
    tree->SetBranchAddress("isSignal", &isSignal);
    tree->SetBranchAddress("trueNuE", &trueNuE);
    tree->SetBranchAddress("nMuMinus", &nMuMinus);
    tree->SetBranchAddress("nMuPlus", &nMuPlus);
    tree->SetBranchAddress("nPiPlus", &nPiPlus);
    tree->SetBranchAddress("nPiMinus", &nPiMinus);
    tree->SetBranchAddress("nKPlus", &nKPlus);
    tree->SetBranchAddress("nKMinus", &nKMinus);
    tree->SetBranchAddress("nProton", &nProton);
    tree->SetBranchAddress("nNeutron", &nNeutron);
    tree->SetBranchAddress("nPhoton", &nPhoton);
    tree->SetBranchAddress("nElectron", &nElectron);
    tree->SetBranchAddress("nPositron", &nPositron);
    tree->SetBranchAddress("nTotal", &nTotal);
    tree->SetBranchAddress("topologicalScore", &topologicalScore);
    tree->SetBranchAddress("isRecoNuFiducial", &isRecoNuFiducial);
    tree->SetBranchAddress("nFinalStatePFPs", &nFinalStatePFPs);
    tree->SetBranchAddress("recoNuVtx", &recoNuVtx);
    tree->SetBranchAddress("hasTrackInfoVect", &hasTrackInfoVect);
    tree->SetBranchAddress("isContainedVect", &isContainedVect);
    tree->SetBranchAddress("startXVect", &startXVect);
    tree->SetBranchAddress("startYVect", &startYVect);
    tree->SetBranchAddress("startZVect", &startZVect);
    tree->SetBranchAddress("endXVect", &endXVect);
    tree->SetBranchAddress("endYVect", &endYVect);
    tree->SetBranchAddress("endZVect", &endZVect);
    tree->SetBranchAddress("directionXVect", &directionXVect);
    tree->SetBranchAddress("directionYVect", &directionYVect);
    tree->SetBranchAddress("directionZVect", &directionZVect);
    tree->SetBranchAddress("trackShowerVect", &trackShowerVect);
    tree->SetBranchAddress("chi2pWVect", &chi2pWVect);
    tree->SetBranchAddress("chi2pUVVect", &chi2pUVVect);
    
    std::cout << "Set up branches" << std::endl;

    // Event loop
    std::map<std::string, std::vector<unsigned int> > totalMap;
    std::map<std::string, std::vector<unsigned int> > selectedMap;

    for (unsigned int i = 0; i < std::min(9999999u, static_cast<unsigned int>(tree->GetEntries())); ++i)
    {
        tree->GetEntry(i);

        const auto interaction = GetInteractionType(nMuMinus, nMuPlus, nPiPlus, nPiMinus, nKPlus, nKMinus, nProton, nNeutron, nPhoton, nElectron, nPositron, nTotal, isTrueNuFiducial, isSignal);
        totalMap[interaction].push_back(i);

        if (!isRecoNuFiducial)
            continue;

        // Collect the PFParticles that start near the vertex - these are the candidate primaries
        std::vector<unsigned int> pfpsNearVertex;
        for (unsigned int j = 0; j < nFinalStatePFPs; ++j)
        {
            const auto hasTrackInfo = hasTrackInfoVect->at(j);
            if (!hasTrackInfo)
                continue;

            const auto startX = startXVect->at(j);
            const auto startY = startYVect->at(j);
            const auto startZ = startZVect->at(j);
            const auto start = TVector3(startX, startY, startZ);

            if (IsNear(*recoNuVtx, start, 4.f))
                pfpsNearVertex.push_back(j);
        }

        // Insist that we have at least 2 PFParticles starting near to the vertex (muon & pion)
        if (pfpsNearVertex.size() < 2)
            continue;

        // Insist that all other PFParticles start near the end point of one of the candidate primaries
        // or are far away and track like (CR in the slice)
        bool areSecondariesViable = true;
        for (unsigned int j = 0; j < nFinalStatePFPs; ++j)
        {
            // Don't consider the primary candidates in this loop
            if (std::find(pfpsNearVertex.begin(), pfpsNearVertex.end(), j) != pfpsNearVertex.end())
                continue;
            
            const auto startX = startXVect->at(j);
            const auto startY = startYVect->at(j);
            const auto startZ = startZVect->at(j);
            const auto start = TVector3(startX, startY, startZ);
            
            const auto endX = endXVect->at(j);
            const auto endY = endYVect->at(j);
            const auto endZ = endZVect->at(j);
            const auto end = TVector3(endX, endY, endZ);
            
            const auto trackShower = trackShowerVect->at(j);

            // Check if it's near an existing primary (the ones near the vertex)
            bool foundNearPrimary = false;
            for (unsigned int k = 0; k < pfpsNearVertex.size(); ++k)
            {
                const auto primaryIndex = pfpsNearVertex.at(k);
                const auto primaryEndX = endXVect->at(primaryIndex);
                const auto primaryEndY = endYVect->at(primaryIndex);
                const auto primaryEndZ = endZVect->at(primaryIndex);
                const auto primaryEnd = TVector3(primaryEndX, primaryEndY, primaryEndZ);

                if (IsNear(start, primaryEnd, 4.f) || IsNear(end, primaryEnd, 4.f))
                {
                    foundNearPrimary = true;
                    break;
                }
            }
        
            // Check if it's far from all other existing primaries
            bool foundFarAway = true;
            for (unsigned int k = 0; k < pfpsNearVertex.size(); ++k)
            {
                const auto primaryIndex = pfpsNearVertex.at(k);
                const auto primaryStartX = endXVect->at(primaryIndex);
                const auto primaryStartY = endYVect->at(primaryIndex);
                const auto primaryStartZ = endZVect->at(primaryIndex);
                const auto primaryStart = TVector3(primaryStartX, primaryStartY, primaryStartZ);
                const auto primaryEndX = endXVect->at(primaryIndex);
                const auto primaryEndY = endYVect->at(primaryIndex);
                const auto primaryEndZ = endZVect->at(primaryIndex);
                const auto primaryEnd = TVector3(primaryEndX, primaryEndY, primaryEndZ);

                if (IsNear(start, primaryStart, 42.f) || IsNear(start, primaryEnd, 42.f) ||
                    IsNear(end, primaryStart, 42.f) || IsNear(end, primaryEnd, 42.f))
                {
                    foundFarAway = false;
                    break;
                }
            }

            // Must be nearby an existing primary particle end point, or far away and track like (CR)
            if (!foundNearPrimary && !(foundFarAway && trackShower > 0.5))
            {
                areSecondariesViable = false;
                break;
            }
        }
        
        if (!areSecondariesViable)
            continue;

        // Find the number of uncontained particles, and flag the index of the last one as the muon (used later)
        unsigned int nUncontained = 0;
        unsigned int muonIndex = std::numeric_limits<unsigned int>::max();
        for (unsigned int k = 0; k < pfpsNearVertex.size(); ++k)
        {
            const auto primaryIndex = pfpsNearVertex.at(k);
            const auto isContained = isContainedVect->at(primaryIndex);

            if (!isContained)
            {
                nUncontained++;
                muonIndex = primaryIndex;
            }
        }

        // Insist that only one primary PFParticle can be uncontained, and if so call it the muon
        if (nUncontained > 1)
            continue;

        // Find the proton candidates
        std::vector<unsigned int> protons;
        for (unsigned int k = 0; k < pfpsNearVertex.size(); ++k)
        {
            const auto primaryIndex = pfpsNearVertex.at(k);

            const auto trackShower = trackShowerVect->at(primaryIndex);
            const auto chi2pW = chi2pWVect->at(primaryIndex);
            const auto chi2pUV = chi2pUVVect->at(primaryIndex);

            const auto hasTrackShower = (trackShower > -0.5);
            const auto hasChi2pW = (chi2pW > -0.5);
            const auto hasChi2pUV = (chi2pUV > -0.5);

            // The proton can't be the muon
            if (nUncontained == 1 && primaryIndex == muonIndex)
                continue;
        
            // Must have the PID info
            if (!hasTrackShower || (!hasChi2pW && !hasChi2pUV))
                continue;

            // Proton must be vaguely track-like
            if (trackShower < 0.2)
                continue;

            // Proton must have a decent chi2 under proton hypothesis
            if (((hasChi2pW && chi2pW > 60) || !hasChi2pW) && ((hasChi2pUV && chi2pUV > 30) || !hasChi2pUV))
                continue;

            protons.push_back(primaryIndex);
        }

        // Insist that we have 2 non-proton primaries (the muon & the pion)
        if (protons.size() + 2 != pfpsNearVertex.size())
            continue;

        // Make sure that the muon & pion candidates aren't back-to-back
        std::vector<unsigned int> muonPions;
        for (unsigned int k = 0; k < pfpsNearVertex.size(); ++k)
        {
            const auto primaryIndex = pfpsNearVertex.at(k);

            if (std::find(protons.begin(), protons.end(), primaryIndex) != protons.end())
                continue;

            muonPions.push_back(primaryIndex);
        }

        // Check to see if we ballzed up
        if (muonPions.size() != 2)
        {
            std::cerr << "Logic error, got " << muonPions.size() << " muon/pion candidates. We should have 2"  << std::endl;
            return 1;
        }

        const auto muonPionIndex0 = muonPions.at(0);
        const auto muonPionIndex1 = muonPions.at(1);

        /*
        std::cout << "----------------------" << std::endl;
        std::cout << run << " " << subRun << " " << event << std::endl;
        std::cout << interaction << std::endl;
        std::cout << "nPFPs     : " << nFinalStatePFPs << std::endl;
        std::cout << "primaries : " << pfpsNearVertex.size() << std::endl;
        std::cin.get();
        */

        // No back to backsies
        /*
        const auto muonPionDirection0 = TVector3(directionXVect->at(muonPionIndex0), directionYVect->at(muonPionIndex0), directionZVect->at(muonPionIndex0));
        const auto muonPionDirection1 = TVector3(directionXVect->at(muonPionIndex1), directionYVect->at(muonPionIndex1), directionZVect->at(muonPionIndex1));

        if (muonPionDirection0.Dot(muonPionDirection1 * (-1)) > 0.99)
            continue;
            */

        if (topologicalScore < 0.1)
            continue;

        selectedMap[interaction].push_back(i);
    }


    // Sort the interactions, signal first, then backgrounds ordered by the number selected
    std::vector<std::pair<std::string, unsigned int> > selectedVector;
    unsigned int totalSelected = 0;
    for (const auto &entry : totalMap)
    {
        const auto interaction = entry.first;
        const auto selectedIter = selectedMap.find(entry.first);
        const auto nSelected = (selectedIter == selectedMap.end()) ? 0 : selectedIter->second.size();
        totalSelected += nSelected;
        selectedVector.emplace_back(interaction, nSelected);
    }
   
    // Do the sort
    std::sort(selectedVector.begin(), selectedVector.end(), [](const std::pair<std::string, unsigned int> &a, const std::pair<std::string, unsigned int> &b) -> bool {
        const bool isASignal = (a.first.substr(0, std::min(6u, static_cast<unsigned int>(a.first.length()))) == "signal");
        const bool isBSignal = (b.first.substr(0, std::min(6u, static_cast<unsigned int>(b.first.length()))) == "signal");

        if (isASignal && !isBSignal)
            return true;
        
        if (!isASignal && isBSignal)
            return false;

        if (a.second == b.second || (isASignal && isASignal))
            return a.first < b.first;

        return a.second > b.second;
    });

    // Print the outcome
    std::cout << std::setw(60) << "interaction  " << "| " << std::setw(8) << "selected" << " / " << std::setw(6) << "total" << " | " << std::setw(10) << "efficiency" << " | " << std::setw(10) << "purity" << std::endl;
    for (unsigned int i = 0; i < selectedVector.size(); ++i)
    {
        const auto entry = selectedVector.at(i);
        const auto interaction = entry.first;
        const auto nTotal = totalMap.at(interaction).size();
        const auto nSelected = entry.second;
        const auto efficiency = static_cast<float>(nSelected) / static_cast<float>(nTotal);
        const auto purity = static_cast<float>(nSelected) / static_cast<float>(totalSelected);

        if (i < 20)
            std::cout << std::setw(60) << interaction << "| " << std::setw(8) << nSelected << " / " << std::setw(6) << nTotal << " | " << std::setw(10) << efficiency << " | " << std::setw(10) << purity << std::endl;
    }
}
