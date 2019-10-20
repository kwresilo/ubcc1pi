namespace cc1pievsel
{

/**
 *  @brief  Particle class holding all the details of a PFParticle
 */
class Particle
{
    public:
        unsigned int index;

        float trackShower;
        bool hasTrackInfo;
        bool isContained;
        TVector3 start;
        TVector3 end;
        TVector3 direction;
        float chi2pW;
        float chi2pUV;
};

// =========================================================================================================================================

/**
 *  @brief  EventManager class holding all the details of the current event
 */
class EventManager
{
    public:
        /**
         *  @brief  Constructor
         *
         *  @param  fileName the file name to load
         *  @param  dirName the directory name that holds the tree we want
         *  @param  treeName the tree name
         */
        EventManager(const std::string &fileName, const std::string &dirName = "eventSelection", const std::string &treeName = "events");

        /**
         *  @brief  Get the number of events
         *
         *  @return number of events
         */
        unsigned int GetNEvents() const;

        /**
         *  @brief  Load the event with the given index
         *
         *  @param  eventIndex the event index to load
         */
        void LoadEvent(const unsigned int eventIndex);

        /**
         *  @brief  Get the index of the event currently loaded
         *
         *  @return the current event index
         */
        unsigned int GetCurrentEventIndex() const;

        /**
         *  @brief  Get the vector of PFParticles
         *
         *  @return the particle vector
         */
        std::vector<Particle> GetParticles() const;

        /**
         *  @brief  Get a string uniquely describing the topology of the current event
         *
         *  @return the topology string
         */
        std::string GetTopologyString() const;

        // ATTN Here we haven't used the m_ prefix for readability in the analysis code. Beware these variables are non-const so they can be
        // changed by the analysis code, poor design... but I'm too lazy to write the getter functions.

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

    private:
        /**
         *  @brief  Load the input tree from the file
         */
        void LoadInputTree();

        /**
         *  @brief  Set the branch addresses
         */
        void SetBranchAddresses();

        // The input file details
        std::string m_fileName;
        std::string m_dirName;
        std::string m_treeName;

        // The tree from which we read the data
        TTree *m_tree;

        // Information about the current event
        unsigned int m_currentEventIndex;
        std::vector<Particle> m_particles;
        
        // Particle info
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
        
};

// -----------------------------------------------------------------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------------------------------------------------------------

EventManager::EventManager(const std::string &fileName, const std::string &dirName, const std::string &treeName) : 
    m_fileName(fileName),
    m_dirName(dirName),
    m_treeName(treeName)
{
    this->LoadInputTree();
    this->SetBranchAddresses();

    if (this->GetNEvents() == 0)
        throw std::invalid_argument("Input file has no events to read!");
    
    this->LoadEvent(0);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  Load the input tree from the file specified
 */
void EventManager::LoadInputTree()
{
    TFile *file = TFile::Open(m_fileName.c_str());
    TDirectoryFile *dir;
    file->GetObject(m_dirName.c_str(), dir);
    dir->GetObject(m_treeName.c_str(), m_tree);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  Set the branch addresses for the input tree
 */
void EventManager::SetBranchAddresses()
{
    m_tree->SetBranchAddress("run", &run);
    m_tree->SetBranchAddress("subRun", &subRun);
    m_tree->SetBranchAddress("event", &event);
    m_tree->SetBranchAddress("isTrueNuFiducial", &isTrueNuFiducial);
    m_tree->SetBranchAddress("isSignal", &isSignal);
    m_tree->SetBranchAddress("trueNuE", &trueNuE);
    m_tree->SetBranchAddress("nMuMinus", &nMuMinus);
    m_tree->SetBranchAddress("nMuPlus", &nMuPlus);
    m_tree->SetBranchAddress("nPiPlus", &nPiPlus);
    m_tree->SetBranchAddress("nPiMinus", &nPiMinus);
    m_tree->SetBranchAddress("nKPlus", &nKPlus);
    m_tree->SetBranchAddress("nKMinus", &nKMinus);
    m_tree->SetBranchAddress("nProton", &nProton);
    m_tree->SetBranchAddress("nNeutron", &nNeutron);
    m_tree->SetBranchAddress("nPhoton", &nPhoton);
    m_tree->SetBranchAddress("nElectron", &nElectron);
    m_tree->SetBranchAddress("nPositron", &nPositron);
    m_tree->SetBranchAddress("nTotal", &nTotal);
    m_tree->SetBranchAddress("topologicalScore", &topologicalScore);
    m_tree->SetBranchAddress("isRecoNuFiducial", &isRecoNuFiducial);
    m_tree->SetBranchAddress("nFinalStatePFPs", &nFinalStatePFPs);
    m_tree->SetBranchAddress("recoNuVtx", &recoNuVtx);
    m_tree->SetBranchAddress("hasTrackInfoVect", &hasTrackInfoVect);
    m_tree->SetBranchAddress("isContainedVect", &isContainedVect);
    m_tree->SetBranchAddress("startXVect", &startXVect);
    m_tree->SetBranchAddress("startYVect", &startYVect);
    m_tree->SetBranchAddress("startZVect", &startZVect);
    m_tree->SetBranchAddress("endXVect", &endXVect);
    m_tree->SetBranchAddress("endYVect", &endYVect);
    m_tree->SetBranchAddress("endZVect", &endZVect);
    m_tree->SetBranchAddress("directionXVect", &directionXVect);
    m_tree->SetBranchAddress("directionYVect", &directionYVect);
    m_tree->SetBranchAddress("directionZVect", &directionZVect);
    m_tree->SetBranchAddress("trackShowerVect", &trackShowerVect);
    m_tree->SetBranchAddress("chi2pWVect", &chi2pWVect);
    m_tree->SetBranchAddress("chi2pUVVect", &chi2pUVVect);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void EventManager::LoadEvent(const unsigned int eventIndex)
{
    m_tree->GetEntry(eventIndex);
    m_currentEventIndex = eventIndex;
    
    // Read the particle-level vectors
    m_particles.clear();
    for (unsigned int index = 0; index < nFinalStatePFPs; ++index)
    {
        Particle p;
        p.index = index;
        p.trackShower = trackShowerVect->at(index);
        p.hasTrackInfo = hasTrackInfoVect->at(index);
        p.isContained = isContainedVect->at(index);
        p.start = TVector3(startXVect->at(index), startYVect->at(index), startZVect->at(index));
        p.end = TVector3(endXVect->at(index), endYVect->at(index), endZVect->at(index));
        p.direction = TVector3(directionXVect->at(index), directionYVect->at(index), directionZVect->at(index));
        p.chi2pW = chi2pWVect->at(index);
        p.chi2pUV = chi2pUVVect->at(index);

        m_particles.push_back(p);
    }
}

// -----------------------------------------------------------------------------------------------------------------------------------------

unsigned int EventManager::GetNEvents() const
{
    return m_tree->GetEntries();
}

// -----------------------------------------------------------------------------------------------------------------------------------------
        
unsigned int EventManager::GetCurrentEventIndex() const
{
    return m_currentEventIndex;
}

// -----------------------------------------------------------------------------------------------------------------------------------------
        
std::vector<Particle> EventManager::GetParticles() const
{
    return m_particles;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

std::string EventManager::GetTopologyString() const
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

// =========================================================================================================================================

/**
 *  @brief  A selection counter keeps track of the events selected after different cuts have been applied
 */
class SelectionCounter
{
    public:
        typedef std::unordered_map<std::string, std::vector<unsigned int> > StringToIndicesMap;

        /**
         *  @brief  Add an event passing the supplied cut with the supplied topology
         *
         *  @param  cut the cut that the event passes
         *  @param  topology the of the event
         *  @param  eventIndex the index of the event
         */
        void AddEventPassingCut(const std::string &cut, const std::string &topology, const unsigned int eventIndex);

        /**
         *  @brief  Print the details of the cut for all collected events
         *
         *  @param  nTopologies the number of topologies to print per entry (including signal & background)
         */
        void PrintBreakdown(const unsigned int nTopologies = 20);
        
        /**
         *  @brief  Print the performance of a given cut
         *
         *  @param  cut the cut to show
         *  @param  nTopologies the number of topologies to print per entry (including signal & background)
         */
        void PrintPerformance(const std::string &cut, const unsigned int nTopologies = 20);

    private:
        /**
         *  @brief  Get the vector of topologies ordered by the number that pass a given cut
         *
         *  @param  cut the cut to order by
         *
         *  @return the topologies
         */
        std::vector<std::string> GetOrderedTopologies(const std::string &cut);

        std::vector<std::string>                            m_cuts; ///< The vector of cuts in the order they were seen
        std::vector<std::string>                            m_topologies; ///< The vector of topologies seen
        std::unordered_map<std::string, StringToIndicesMap> m_cutMap; ///< The map from the name of a cut to the breakdown of the events passing the cut organised by their topology
};

// -----------------------------------------------------------------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------------------------------------------------------------

void SelectionCounter::AddEventPassingCut(const std::string &cut, const std::string &topology, const unsigned int eventIndex)
{
    // Add this cut to the vector if we haven't seen it before
    if (std::find(m_cuts.begin(), m_cuts.end(), cut) == m_cuts.end())
        m_cuts.push_back(cut);
    
    // Add this topology to the vector if we haven't seen it before
    if (std::find(m_topologies.begin(), m_topologies.end(), topology) == m_topologies.end())
        m_topologies.push_back(topology);

    m_cutMap[cut][topology].push_back(eventIndex);
}

// -----------------------------------------------------------------------------------------------------------------------------------------
        
std::vector<std::string> SelectionCounter::GetOrderedTopologies(const std::string &cut)
{
    // ATTN non const as we want to add empty entries for topologiesn not seen by this cut
    auto iter = m_cutMap.find(cut);
    if (iter == m_cutMap.end())
        throw std::invalid_argument("Can't get topologies ordered by events passing unknown cut: " + cut);

    // Convert the map to a vector so we can sort it
    auto topologyToEventsMap = iter->second;
    std::vector<std::pair<std::string, unsigned int> > topologyNEventsVector;
    for (const auto &topology : m_topologies)
        topologyNEventsVector.emplace_back(topology, topologyToEventsMap[topology].size());

    // Do the sort
    std::sort(topologyNEventsVector.begin(), topologyNEventsVector.end(), [](const std::pair<std::string, unsigned int> &a, const std::pair<std::string, unsigned int> &b) -> bool {

        // Work out if either are signal topologies
        const bool isASignal = (a.first.substr(0, std::min(6u, static_cast<unsigned int>(a.first.length()))) == "signal");
        const bool isBSignal = (b.first.substr(0, std::min(6u, static_cast<unsigned int>(b.first.length()))) == "signal");

        // Always order a signal topology before a background
        if (isASignal && !isBSignal)
            return true;
        
        if (!isASignal && isBSignal)
            return false;

        // If both are signal, or both have the same number of events passing the cut - then just order alphabetically by the topology name
        if (a.second == b.second || (isASignal && isASignal))
            return a.first < b.first;

        // Order by the number of events passing the cut
        return a.second > b.second;
    });

    // Extract the topologies in the correct order
    std::vector<std::string> orderedTopologies;
    for (const auto &entry : topologyNEventsVector)
        orderedTopologies.push_back(entry.first);

    return orderedTopologies;
}

// -----------------------------------------------------------------------------------------------------------------------------------------
        
void SelectionCounter::PrintBreakdown(const unsigned int nTopologies = 20)
{
    if (m_cuts.empty())
    {
        std::cout << "No events added" << std::endl;
        return;
    }

    for (const auto &cut : m_cuts)
        this->PrintPerformance(cut, nTopologies);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void SelectionCounter::PrintPerformance(const std::string &cut, const unsigned int nTopologies = 20)
{
    auto iter = m_cutMap.find(cut);
    if (iter == m_cutMap.end())
        throw std::invalid_argument("Unknown cut: " + cut);
    
    const auto firstCut = m_cuts.front();
        
    // Get the topologies ordered accoring to the cut
    const auto topologies = this->GetOrderedTopologies(cut);

    std::cout << std::string(140, '-') << std::endl;
    std::cout << "Comparing performance between : " << firstCut << " & " << cut << std::endl;

    // Get the total number of events passing the cut, and count the integrated signal
    unsigned int nTotal = 0;
    unsigned int nTotalInitial = 0;
    unsigned int nSignal = 0;
    unsigned int nSignalInitial = 0;
    for (const auto &topology : topologies)
    {
        nTotal += m_cutMap[cut][topology].size();
        nTotalInitial += m_cutMap[firstCut][topology].size();
        
        const bool isSignal = (topology.substr(0, std::min(6u, static_cast<unsigned int>(topology.length()))) == "signal");
        if (!isSignal)
            continue;

        nSignal += m_cutMap[cut][topology].size();
        nSignalInitial += m_cutMap[firstCut][topology].size();
    }
    
    std::cout << " - N events initial     : " << nTotalInitial << std::endl;
    std::cout << " - N events passing cut : " << nTotal << std::endl;
    std::cout << " - N signal initial     : " << nSignalInitial << std::endl;
    std::cout << " - N signal passing cut : " << nSignal << std::endl;

    // Get the purity and efficiency of the selection at this cut
    float purity = -1;
    if (nTotal != 0)
        purity = static_cast<float>(nSignal) / static_cast<float>(nTotal);
    
    float efficiency = -1;
    if (nSignalInitial != 0)
        efficiency = static_cast<float>(nSignal) / static_cast<float>(nSignalInitial);
    
    std::cout << " - Purity               : " << purity << std::endl;
    std::cout << " - Efficiency           : " << efficiency << std::endl;

    // Print the detailed breakdown per topology
    std::cout << std::string(140, '-') << std::endl;
    std::cout << std::setw(60) << "topology" << " | ";
    std::cout << std::setw(8) << "nPassed" << " / ";
    std::cout << std::setw(8) << "nInitial" << " | ";
    std::cout << std::setw(12) << "fracPassed" << " | ";
    std::cout << std::setw(12) << "fracOfTotal" << std::endl;

    // Get the total number of events passing the cut for each topology
    for (unsigned int i = 0; i < std::min(nTopologies, static_cast<unsigned int>(topologies.size())); ++i)
    {
        const auto &topology(topologies.at(i));
        const auto nFirst = m_cutMap[firstCut][topology].size();
        const auto nCut = m_cutMap[cut][topology].size();

        // Get the fraction selected (selection efficiency)
        float fracSelected = -1;
        if (nFirst != 0)
            fracSelected = static_cast<float>(nCut) / static_cast<float>(nFirst);

        // Get the fraction of all events after the cut that are of this topology (selectionPurity)
        float fracTotal = -1;
        if (nTotal != 0)
            fracTotal = static_cast<float>(nCut) / static_cast<float>(nTotal);

        std::cout << std::setw(60) << topology << " | ";
        std::cout << std::setw(8) << nCut << " / ";
        std::cout << std::setw(8) << nFirst << " | ";
        std::cout << std::setw(12) << fracSelected << " | ";
        std::cout << std::setw(12) << fracTotal << std::endl;
    }
}

// =========================================================================================================================================

/**
 *  @brief  Determine if two 3D positions are with the specified separation
 *
 *  @param  posA the first position
 *  @param  posB the second position
 *  @param  cut the maximum separation threshold
 *
 *  @return bool - true if the separation between posA and posB is within the cut
 */
/*
bool IsNear(const TVector3 &posA, const TVector3 &posB, const float cut)
{
    return ((posA - posB).Mag2() < cut * cut);
}
*/

// -----------------------------------------------------------------------------------------------------------------------------------------


} // namespace cc1pievsel

// -----------------------------------------------------------------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------------------------------------------------------------

using namespace cc1pievsel;

/**
 *  @brief  The main function called from the command line: root -l makeSelectionTableNew.C
 */
void makeSelectionTableNew()
{
    const unsigned int nEvents = std::numeric_limits<unsigned int>::max();

    EventManager em("/uboone/data/users/asmith/ubcc1pi/19062019/eventSelection.root");
    SelectionCounter counter;

    // Event loop
    for (unsigned int eventIndex = 0; eventIndex < std::min(nEvents, em.GetNEvents()); ++eventIndex)
    {
        em.LoadEvent(eventIndex);
        const auto topology = em.GetTopologyString();

        counter.AddEventPassingCut("all", topology, eventIndex);

        const auto allParticles = em.GetParticles();
        if (allParticles.size() < 2)
            continue;
        
        counter.AddEventPassingCut("min2Particles", topology, eventIndex);
    }

    counter.PrintBreakdown();
    

    /*
    std::map<std::string, std::vector<unsigned int> > totalMap;
    std::map<std::string, std::vector<unsigned int> > selectedMap;

    for (unsigned int i = 0; i < std::min(9999999u, static_cast<unsigned int>(tree->GetEntries())); ++i)
    {
        tree->GetEntry(i);

        const auto interaction = GetInteractionType();
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

        std::cout << "----------------------" << std::endl;
        std::cout << run << " " << subRun << " " << event << std::endl;
        std::cout << interaction << std::endl;
        std::cout << "nPFPs     : " << nFinalStatePFPs << std::endl;
        std::cout << "primaries : " << pfpsNearVertex.size() << std::endl;
        std::cin.get();

        // No back to backsies
        const auto muonPionDirection0 = TVector3(directionXVect->at(muonPionIndex0), directionYVect->at(muonPionIndex0), directionZVect->at(muonPionIndex0));
        const auto muonPionDirection1 = TVector3(directionXVect->at(muonPionIndex1), directionYVect->at(muonPionIndex1), directionZVect->at(muonPionIndex1));

        if (muonPionDirection0.Dot(muonPionDirection1 * (-1)) > 0.99)
            continue;

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
    
    */
}
