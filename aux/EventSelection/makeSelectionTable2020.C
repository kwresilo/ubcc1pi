namespace cc1pievsel
{

/**
 *  @brief  Particle class holding all the details of a PFParticle
 */
class Particle
{
    public:
        unsigned int index;

        int nHitsU;
        int nHitsV;
        int nHitsW;
        float trackShower;
        bool hasTrackInfo;
        bool isContained;
        TVector3 start;
        TVector3 end;
        TVector3 direction;
        float length;
        float braggpW;
        float braggMIPW;
        float braggpUV;
        float braggMIPUV;
        bool hasMatchedMCParticle;
        int truePdgCode;
        float truthMatchPurity;
        float truthMatchCompleteness;
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
        float trueMuEnergy = -std::numeric_limits<float>::max();
        float trueMuTheta = -std::numeric_limits<float>::max();
        float trueMuPhi = -std::numeric_limits<float>::max();
        float truePiEnergy = -std::numeric_limits<float>::max();
        float truePiTheta = -std::numeric_limits<float>::max();
        float truePiPhi = -std::numeric_limits<float>::max();
        float trueMuPiAngle = -std::numeric_limits<float>::max();
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
        std::vector<int> *nHitsUVect = nullptr;
        std::vector<int> *nHitsVVect = nullptr;
        std::vector<int> *nHitsWVect = nullptr;
        std::vector<bool> *hasTrackInfoVect = nullptr;
        std::vector<bool> *isContainedVect = nullptr;
        std::vector<float> *lengthVect = nullptr;
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
        std::vector<float> *braggpWVect = nullptr;
        std::vector<float> *braggMIPWVect = nullptr;
        std::vector<float> *braggpUVVect = nullptr;
        std::vector<float> *braggMIPUVVect = nullptr;
        std::vector<bool> *hasMatchedMCParticleVect = nullptr;
        std::vector<int> *truePdgCodeVect = nullptr;
        std::vector<float> *truthMatchPurityVect = nullptr;
        std::vector<float> *truthMatchCompletenessVect = nullptr;
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
    m_tree->SetBranchAddress("trueMuEnergy", &trueMuEnergy);
    m_tree->SetBranchAddress("trueMuTheta", &trueMuTheta);
    m_tree->SetBranchAddress("trueMuPhi", &trueMuPhi);
    m_tree->SetBranchAddress("truePiEnergy", &truePiEnergy);
    m_tree->SetBranchAddress("truePiTheta", &truePiTheta);
    m_tree->SetBranchAddress("truePiPhi", &truePiPhi);
    m_tree->SetBranchAddress("trueMuPiAngle", &trueMuPiAngle);
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
    m_tree->SetBranchAddress("nHitsUVect", &nHitsUVect);
    m_tree->SetBranchAddress("nHitsVVect", &nHitsVVect);
    m_tree->SetBranchAddress("nHitsWVect", &nHitsWVect);
    m_tree->SetBranchAddress("hasTrackInfoVect", &hasTrackInfoVect);
    m_tree->SetBranchAddress("lengthVect", &lengthVect);
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
    m_tree->SetBranchAddress("braggpWVect", &braggpWVect);
    m_tree->SetBranchAddress("braggMIPWVect", &braggMIPWVect);
    m_tree->SetBranchAddress("braggpUVVect", &braggpUVVect);
    m_tree->SetBranchAddress("braggMIPUVVect", &braggMIPUVVect);
    m_tree->SetBranchAddress("hasMatchedMCParticleVect", &hasMatchedMCParticleVect);
    m_tree->SetBranchAddress("truePdgCodeVect", &truePdgCodeVect);
    m_tree->SetBranchAddress("truthMatchPurityVect", &truthMatchPurityVect);
    m_tree->SetBranchAddress("truthMatchCompletenessVect", &truthMatchCompletenessVect);
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

        p.nHitsU = nHitsUVect->at(index);
        p.nHitsV = nHitsVVect->at(index);
        p.nHitsW = nHitsWVect->at(index);
        p.trackShower = trackShowerVect->at(index);
        p.hasTrackInfo = hasTrackInfoVect->at(index);
        p.length = lengthVect->at(index);
        p.isContained = isContainedVect->at(index);
        p.start = TVector3(startXVect->at(index), startYVect->at(index), startZVect->at(index));
        p.end = TVector3(endXVect->at(index), endYVect->at(index), endZVect->at(index));
        p.direction = TVector3(directionXVect->at(index), directionYVect->at(index), directionZVect->at(index));
        p.braggpW = braggpWVect->at(index);
        p.braggMIPW = braggMIPWVect->at(index);
        p.braggpUV = braggpUVVect->at(index);
        p.braggMIPUV = braggMIPUVVect->at(index);
        p.hasMatchedMCParticle = hasMatchedMCParticleVect->at(index);
        p.truePdgCode = truePdgCodeVect->at(index);
        p.truthMatchPurity = truthMatchPurityVect->at(index);
        p.truthMatchCompleteness = truthMatchCompletenessVect->at(index);

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

bool IsSignal(const std::string &topology)
{
    return (topology.substr(0, std::min(6u, static_cast<unsigned int>(topology.length()))) == "signal");
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
    
    /*
    if (nProton != 0)
        interaction += std::to_string(nProton) + " p  ";
    */
    interaction += "X p  ";

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
        void PrintBreakdown(const unsigned int nTopologies = 5);
        
        /**
         *  @brief  Print the performance of a given cut
         *
         *  @param  cut the cut to show
         *  @param  nTopologies the number of topologies to print per entry (including signal & background)
         */
        void PrintPerformance(const std::string &cut, const unsigned int nTopologies = 20);

        /**
         *  @brief  Get the list of cuts applied in the order they were first seen
         *
         *  @return the names of the cuts
         */
        std::vector<std::string> GetCuts() const;

        /**
         *  @brief  Get the indices of the events passing the cut
         *
         *  @param  cut the cut
         *
         *  @return the event indices
         */
        std::vector<unsigned int> GetEventIndices(const std::string &cut) const;
        
        /**
         *  @brief  Get the indices of the signal events passing the cut
         *
         *  @param  cut the cut
         *
         *  @return the event indices
         */
        std::vector<unsigned int> GetSignalEventIndices(const std::string &cut) const;

    private:
        /**
         *  @brief  Get the vector of topologies ordered by the number that pass a given cut
         *
         *  @param  cut the cut to order by
         *
         *  @return the topologies
         */
        std::vector<std::string> GetOrderedTopologies(const std::string &cut) const;

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
        
std::vector<std::string> SelectionCounter::GetCuts() const
{
    return m_cuts;
}

// -----------------------------------------------------------------------------------------------------------------------------------------
       
std::vector<unsigned int> SelectionCounter::GetEventIndices(const std::string &cut) const
{
    const auto iter = m_cutMap.find(cut);
    if (iter == m_cutMap.end())
        throw std::invalid_argument("Can't get event indices for unknown cut: " + cut);

    std::vector<unsigned int> eventIndices;
    for (const auto &entry : iter->second)
        eventIndices.insert(eventIndices.end(), entry.second.begin(), entry.second.end());
    
    std::sort(eventIndices.begin(), eventIndices.end());

    return eventIndices;
}

// -----------------------------------------------------------------------------------------------------------------------------------------
       
std::vector<unsigned int> SelectionCounter::GetSignalEventIndices(const std::string &cut) const
{
    const auto iter = m_cutMap.find(cut);
    if (iter == m_cutMap.end())
        throw std::invalid_argument("Can't get event indices for unknown cut: " + cut);

    std::vector<unsigned int> eventIndices;
    for (const auto &entry : iter->second)
    {
        if (IsSignal(entry.first))
            eventIndices.insert(eventIndices.end(), entry.second.begin(), entry.second.end());
    }

    std::sort(eventIndices.begin(), eventIndices.end());

    return eventIndices;
}

// -----------------------------------------------------------------------------------------------------------------------------------------
        
std::vector<std::string> SelectionCounter::GetOrderedTopologies(const std::string &cut) const
{
    const auto iter = m_cutMap.find(cut);
    if (iter == m_cutMap.end())
        throw std::invalid_argument("Can't get topologies ordered by events passing unknown cut: " + cut);

    // Convert the map to a vector so we can sort it
    const auto &topologyToEventsMap = iter->second;
    std::vector<std::pair<std::string, unsigned int> > topologyNEventsVector;
    for (const auto &topology : m_topologies)
    {
        const auto iterTop = topologyToEventsMap.find(topology);
        topologyNEventsVector.emplace_back(topology, (iterTop == topologyToEventsMap.end()) ? 0 : iterTop->second.size());
    }

    // Do the sort
    std::sort(topologyNEventsVector.begin(), topologyNEventsVector.end(), [](const std::pair<std::string, unsigned int> &a, const std::pair<std::string, unsigned int> &b) -> bool {

        // Work out if either are signal topologies
        const bool isASignal = IsSignal(a.first); 
        const bool isBSignal = IsSignal(b.first);

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
    std::cout << "Comparing performance between : " << firstCut << " -> " << cut << std::endl;

    // Get the total number of events passing the cut, and count the integrated signal
    unsigned int nTotal = 0;
    unsigned int nTotalInitial = 0;
    unsigned int nSignal = 0;
    unsigned int nSignalInitial = 0;
    for (const auto &topology : topologies)
    {
        nTotal += m_cutMap[cut][topology].size();
        nTotalInitial += m_cutMap[firstCut][topology].size();
        
        const bool isSignal = IsSignal(topology); 
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
 *  @brief  The event level plot class. Holds event-level plots after various cuts
 */
class EventPlot
{
    public:
        /**
         *  @brief  Constructor
         *
         *  @param  title the title of the plot
         *  @param  xAxisLabel the x-axis label
         *  @param  nBins the number of bins
         *  @param  min the minimum value
         *  @param  max the maximum value
         */
        EventPlot(const std::string &title, const std::string &xAxisLabel, const unsigned int nBins, const float min, const float max);

        /**
         *  @brief  Destructor
         */
        ~EventPlot();

        /**
         *  @brief  Fill the histogram related to the supplied cut at the given value
         *
         *  @param  value value to fill
         *  @param  cut the cut to use
         */
        void Fill(const float &value, const std::string &cut);

        /**
         *  @brief  Draw the histograms and save them to a file
         *
         *  @param  cuts the cuts to draw, also make comparison plots to the first cut
         */
        void Draw(const std::vector<std::string> &cuts) const;

    private:
        typedef std::unordered_map<std::string, TH1F*> StringToHistMap;

        std::string   m_title;
        std::string   m_xAxisLabel;
        unsigned int  m_nBins;
        float         m_min;
        float         m_max;

        StringToHistMap m_cutToHistMap;
};

// -----------------------------------------------------------------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------------------------------------------------------------

EventPlot::EventPlot(const std::string &title, const std::string &xAxisLabel, const unsigned int nBins, const float min, const float max) :
    m_title(title),
    m_xAxisLabel(xAxisLabel),
    m_nBins(nBins),
    m_min(min),
    m_max(max)
{
}

// -----------------------------------------------------------------------------------------------------------------------------------------

EventPlot::~EventPlot()
{
    for (auto &entry : m_cutToHistMap)
        entry.second->Delete();
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void EventPlot::Fill(const float &value, const std::string &cut)
{
    auto iter = m_cutToHistMap.find(cut);
    if (iter == m_cutToHistMap.end())
    {
        const auto name = (m_title + "_" + cut).c_str();
        m_cutToHistMap.emplace(cut, new TH1F(name, m_title.c_str(), m_nBins, m_min, m_max));
    }

    m_cutToHistMap.at(cut)->Fill(value);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void EventPlot::Draw(const std::vector<std::string> &cuts) const
{
    if (cuts.empty())
        throw std::invalid_argument("Can't draw, input cut vector is empty");
    
    // Remove nasty characters from the title
    auto cleanTitle = m_title;
    cleanTitle.erase(std::remove_if(cleanTitle.begin(), cleanTitle.end(), [](char c) { return !isalpha(c); } ), cleanTitle.end());
    
    // Make a scaled copy of the histograms
    StringToHistMap cutToHistMapCopy;
    for (unsigned int cutNumber = 0; cutNumber < cuts.size(); ++cutNumber)
    {
        const auto cut = cuts.at(cutNumber);
        auto cleanCut = cut;
        cleanCut.erase(std::remove_if(cleanCut.begin(), cleanCut.end(), [](char c) { return !isalpha(c); } ), cleanCut.end());

        // Find the histogram for this cut
        const auto iter = m_cutToHistMap.find(cut);
        if (iter == m_cutToHistMap.end())
            throw std::invalid_argument("Can't draw, unknown cut: " + cut);

        auto pHist = static_cast<TH1F*>(iter->second->Clone((cleanTitle + "_" + cleanCut + "_clone").c_str()));
        pHist->Scale(1.f / pHist->GetEntries());
        pHist->Sumw2();
        cutToHistMapCopy.emplace(cut, pHist);
    }

    // Get the first histogram (to compare to)
    const auto firstCut = cuts.front();
    const auto iterFirst = cutToHistMapCopy.find(firstCut);
    if (iterFirst == cutToHistMapCopy.end())
        throw std::invalid_argument("Can't draw, unknown first cut: " + firstCut);

    const auto pHistFirst = static_cast<TH1F*>(iterFirst->second->Clone((cleanTitle + "_clone").c_str()));
    pHistFirst->SetLineWidth(2);
    pHistFirst->SetLineColor(kAzure - 2);
    pHistFirst->SetFillColor(kAzure - 2);
    pHistFirst->SetFillStyle(0);
    
    pHistFirst->GetXaxis()->SetTitle(m_xAxisLabel.c_str());
    pHistFirst->GetYaxis()->SetTitle("Fraction of events");
    const auto firstMax = pHistFirst->GetMaximum();

    // Set up the canvas
    TCanvas *c = new TCanvas();
    gStyle->SetOptStat(false);

    // Draw the histograms
    for (unsigned int cutNumber = 0; cutNumber < cuts.size(); ++cutNumber)
    {
        const auto cut = cuts.at(cutNumber);

        // Find the histogram for this cut
        const auto iter = cutToHistMapCopy.find(cut);
        if (iter == cutToHistMapCopy.end())
            throw std::invalid_argument("Can't draw, unknown cut: " + cut);

        auto cleanCut = iter->first;
        cleanCut.erase(std::remove_if(cleanCut.begin(), cleanCut.end(), [](char c) { return !isalpha(c); } ), cleanCut.end());
        
        // Make the plot
        auto pHist = iter->second;
        pHist->SetLineWidth(2);
        pHist->SetLineColor(kOrange - 2);
        pHist->SetFillColor(kOrange - 2);
        pHist->SetFillStyle(0);
        
        const auto yMax = 1.1f * std::max(pHist->GetMaximum(), firstMax);
        pHistFirst->GetYaxis()->SetRangeUser(0, yMax);
        
        pHistFirst->DrawCopy("hist");
        pHistFirst->SetFillStyle(3001);
        pHistFirst->Draw("e2 same");

        pHist->DrawCopy("hist same");
        pHist->SetFillStyle(3002);
        pHist->Draw("e2 same");

        c->SaveAs((cleanTitle + "_" + std::to_string(cutNumber) + "_" + cleanCut + ".png").c_str());

        pHistFirst->SetFillStyle(0);
        pHist->SetFillStyle(0);
    }

    delete c;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

// =========================================================================================================================================
// =========================================================================================================================================

/**
 *  @brief  Get the particles in the input vector that start within the the supplied threshold of the neutrino vertex
 *
 *  @param  allParticles the complete list of primary PFParticles from Pandora
 *  @param  recoNuVtx the reconstructed neutrino vertex position
 *  @param  primaryDist the threshold separation from the vertex
 *
 *  @return the particles passing the threshold
 */
std::vector<Particle> GetParticlesNearVertex(const std::vector<Particle> &allParticles, const TVector3 &recoNuVtx, const float primaryDist)
{
    std::vector<Particle> outputParticles;
    const float cut2 = primaryDist * primaryDist;

    for (const auto &particle : allParticles)
    {       
        if ((particle.start - recoNuVtx).Mag2() < cut2)
            outputParticles.push_back(particle);
    }

    return outputParticles;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  Check if the input vector of particles contains the test particle
 *
 *  @param  particles the input vector of particles
 *  @param  testParticle the test particle
 *
 *  @return boolean, true if particles contains testParticle
 */
bool ContainsParticle(const std::vector<Particle> &particles, const Particle &testParticle)
{
    for (const auto &particle : particles)
    {
        if (particle.index == testParticle.index)
            return true;
    }

    return false;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  Check if two collections of particles are identical
 *
 *  @param  collectionA vector of vector of particles
 *  @param  collectionB vector of vector of particles
 *
 *  @return boolean, true if all particles in collectionA are the same as those in collectionB
 */
bool AreParticlesIdentical(const std::vector< std::vector<Particle> > &collectionA, const std::vector< std::vector<Particle> > &collectionB)
{
    // Get the particle indices from the A collection
    std::vector<unsigned int> indicesA;
    for (const auto &particlesA : collectionA)
    {
        for (const auto &particleA : particlesA)
        {
            indicesA.push_back(particleA.index);
        }
    }
    
    // Get the particle indices from the B collection
    std::vector<unsigned int> indicesB;
    for (const auto &particlesB : collectionB)
    {
        for (const auto &particleB : particlesB)
        {
            indicesB.push_back(particleB.index);
        }
    }

    if (indicesA.size() != indicesB.size())
        return false;

    // Sort the vectors and check if they match
    std::sort(indicesA.begin(), indicesA.end());
    std::sort(indicesB.begin(), indicesB.end());

    return (indicesA == indicesB);
}

// -----------------------------------------------------------------------------------------------------------------------------------------
        
/**
 *  @brief  Get the cosmic ray candidates
 *
 *  @param  allParticles the complete list of particles from Pandora
 *  @param  primaryCandidates the primary candidate particles
 *  @param  recoNuVtx the reconstructed neutrino vertex
 *  @param  crDist the threshold distance beyond which a particle is considered a cosmic-ray
 *
 *  @return the CR candidates
 */
std::vector<Particle> GetCRCandidates(const std::vector<Particle> &allParticles, const std::vector<Particle> &primaryCandidates, const TVector3 &recoNuVtx, const float crDist)
{
    std::vector<Particle> outputParticles;
    const float cut2 = crDist * crDist;

    for (const auto &particle : allParticles)
    {
        // Skip the primary candidates
        if (ContainsParticle(primaryCandidates, particle))
            continue;

        // The particle must be well separated from the neutrino vertex
        if ((particle.start - recoNuVtx).Mag2() < cut2)
            continue;

        // The particle must be well separated from the primary candidates
        bool isWellSeparated = true;
        for (const auto &primary : primaryCandidates)
        {
            if ((particle.start - primary.start).Mag2() < cut2 || (particle.start - primary.end).Mag2() < cut2)
            {
                isWellSeparated = false;
                break;
            }
        }

        if (!isWellSeparated)
            continue;

        // The particle is well separated from the rest of the neutrino interaction, so call it ca cosmic-ray
        outputParticles.push_back(particle);
    }

    return outputParticles;
}

// -----------------------------------------------------------------------------------------------------------------------------------------
        
/**
 *  @brief  Get the candidate secondary particles (daughters of primaries)
 *
 *  @param  allParticles the full list of particles from Pandora
 *  @param  primaryCandidates the list of primary candidate particles
 *  @param  secondaryDist the threshold distance
 *
 *  @return the secondary candidates
 */
std::vector<Particle> GetSecondaryCandidates(const std::vector<Particle> &allParticles, const std::vector<Particle> &primaryCandidates, const float secondaryDist)
{
    std::vector<Particle> outputParticles;
    const auto dist2 = secondaryDist * secondaryDist;

    for (const auto &particle : allParticles)
    {
        // Skip the primary candidates
        if (ContainsParticle(primaryCandidates, particle))
            continue;
        
        // Check if this particle starts nearer to the end of a primary than to the vertex
        bool isNearPrimaryEnd = false;
        for (const auto &primary : primaryCandidates)
        {
            if ((particle.start - primary.end).Mag2() < dist2)
            {
                isNearPrimaryEnd = true;
                break;
            }
        }

        if (!isNearPrimaryEnd)
            continue;

        outputParticles.push_back(particle);
    }
    
    return outputParticles;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  Get the particles from the input vector that are uncontained
 *
 *  @param  particles the input vector of particles
 *
 *  @return the uncontained particles
 */
std::vector<Particle> GetUncontainedParticles(const std::vector<Particle> &particles)
{
    std::vector<Particle> outputParticles;

    for (const auto &particle : particles)
    {
        if (!particle.isContained)
            outputParticles.push_back(particle);
    }

    return outputParticles;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  Check if the input particle has PID info available from the W plane
 *
 *  @param  particle the input particle
 *
 *  @return bool
 */
bool HasPIDWAvailable(const Particle &particle)
{
    const auto hasBraggpW = (particle.braggpW >= -std::numeric_limits<float>::epsilon());
    const auto hasBraggMIPW = (particle.braggMIPW >= -std::numeric_limits<float>::epsilon());
    
    return (hasBraggpW && hasBraggMIPW);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  Check if the input particle has PID info available from the UV planes
 *
 *  @param  particle the input particle
 *
 *  @return bool
 */
bool HasPIDUVAvailable(const Particle &particle)
{
    const auto hasBraggpUV = (particle.braggpUV >= -std::numeric_limits<float>::epsilon());
    const auto hasBraggMIPUV = (particle.braggMIPUV >= -std::numeric_limits<float>::epsilon());

    return (hasBraggpUV && hasBraggMIPUV);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  Check if the input particle has all required PID info available
 *
 *  @param  particle the input particle
 *
 *  @return bool
 */
bool HasPIDAvailable(const Particle &particle)
{
    const auto hasTrackShower = (particle.trackShower >= -std::numeric_limits<float>::epsilon());
    const auto hasLength = (particle.length >= -std::numeric_limits<float>::epsilon());

    const auto hasPIDW = HasPIDWAvailable(particle);
    const auto hasPIDUV = HasPIDUVAvailable(particle);

    return (hasTrackShower && hasLength && (hasPIDW || hasPIDUV));
}

// -----------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  Check if all of the input particles have all of the required PID information available
 *
 *  @param  particles the input particles
 *
 *  @return bool
 */
bool HasPIDAvailable(const std::vector<Particle> &particles)
{
    for (const auto &particle : particles)
    {
        if (!HasPIDAvailable(particle))
            return false;
    }

    return true;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  Get the muon candidates
 *
 *  @param  particles the input vector of particles
 *
 *  @return the muon candidates
 */
std::vector<Particle> GetMuonCandidates(const std::vector<Particle> &particles)
{
    std::vector<Particle> outputParticles;

    for (const auto &particle : particles)
    {
        if (!particle.isContained)
        {
            outputParticles.push_back(particle);
            continue;
        }

        if (!HasPIDAvailable(particle))
            continue;

        if (particle.trackShower < 0.8)
            continue;
        
        const auto hasPIDW = HasPIDWAvailable(particle);
        if (hasPIDW && std::log(particle.braggpW / particle.braggMIPW) > 0.f && particle.length < 250.f)
            continue;

        if (!hasPIDW)
        {
            const auto hasPIDUV = HasPIDUVAvailable(particle);
            if (!hasPIDUV)
                throw std::logic_error("Particle has no PID for W or UV, these should have been filtered out already!");

            if (std::log(particle.braggpUV / particle.braggMIPUV) > 0.f && particle.length < 250.f)
                continue;
        }

        outputParticles.push_back(particle);
    }

    return outputParticles;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  Get the proton candidates
 *
 *  @param  particles the input vector of particles
 *  @param  muonCandidates the muon candidates
 *
 *  @return the proton candidates
 */
std::vector<Particle> GetProtonCandidates(const std::vector<Particle> &particles, const std::vector<Particle> &muonCandidates)
{
    std::vector<Particle> outputParticles;

    for (const auto &particle : particles)
    {
        if (ContainsParticle(muonCandidates, particle))
            continue;

        if (!particle.isContained)
            continue;

        if (!HasPIDAvailable(particle))
            continue;

        if (particle.trackShower < 0.5)
            continue;
        
        const auto hasPIDW = HasPIDWAvailable(particle);
        if (hasPIDW && std::log(particle.braggpW / particle.braggMIPW) < 2.f)
            continue;

        if (!hasPIDW)
        {
            const auto hasPIDUV = HasPIDUVAvailable(particle);
            if (!hasPIDUV)
                throw std::logic_error("Particle has no PID for W or UV, these should have been filtered out already!");

            if (std::log(particle.braggpUV / particle.braggMIPUV) < 2.f)
                continue;
        }

        outputParticles.push_back(particle);
    }

    return outputParticles;
}

} // namespace cc1pievsel

// -----------------------------------------------------------------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------------------------------------------------------------

using namespace cc1pievsel;

/**
 *  @brief  The main function called from the command line: root -l makeSelectionTableNew.C
 */
void makeSelectionTableNew()
{
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Configuration
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    const unsigned int nEvents = std::numeric_limits<unsigned int>::max();

    const float minTopologicalScore = 0.05f;     // The minimum topological score for an event to be selected
    const float primaryDist = 5.f;               // The maximum distance from which a particle can start from the vertex to be called a primary
    const float secondaryDist = 5.f;             // The maximum distance from which a particle can start from a primary to be called a secondary
    const float crDist = 14.f * 3;               // The distance beyond which a particle is considered a cosmic in the slice
    const float maxOpeningAngle = 2.6f;          // The maximum reconstructed muon-pion opening angle

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    EventManager em("/uboone/data/users/asmith/ubcc1pi/19062019/eventSelection.root");
    SelectionCounter counter;

    // Event loop
    for (unsigned int eventIndex = 0; eventIndex < std::min(nEvents, em.GetNEvents()); ++eventIndex)
    {
        em.LoadEvent(eventIndex);
        const auto topology = em.GetTopologyString();
        counter.AddEventPassingCut("all", topology, eventIndex);
        
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        // Insist that the reconstructed neutrino is fiducial
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if (!em.isRecoNuFiducial)
            continue;

        counter.AddEventPassingCut("fiducial", topology, eventIndex);
        
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        // Insist as most one primary particle is uncontained
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        const auto allParticles = em.GetParticles();
        const auto recoNuVtx = *em.recoNuVtx;
        const auto primaryCandidates = GetParticlesNearVertex(allParticles, recoNuVtx, primaryDist);

        const auto uncontainedParticles = GetUncontainedParticles(primaryCandidates);
        if (uncontainedParticles.size() > 1)
            continue;
        
        counter.AddEventPassingCut("maxOneUncontained", topology, eventIndex);
        
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        // Insist that all primary candidates have the required PID info
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        const auto isPIDAvailable = HasPIDAvailable(primaryCandidates);
        if (!isPIDAvailable)
            continue;
        
        counter.AddEventPassingCut("pidAvailable", topology, eventIndex);
        
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        // Insist that there are at least 2 PFParticles that start near the vertex
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if (primaryCandidates.size() < 2)
            continue;
        
        counter.AddEventPassingCut("minTwoPrimaries", topology, eventIndex);
        
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        // Insist that there is at least one muon candidate
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        const auto muonCandidates = GetMuonCandidates(primaryCandidates);
        if (muonCandidates.empty())
            continue;
        
        counter.AddEventPassingCut("minOneMuon", topology, eventIndex);

        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        // Insist that there are exactly 2 non-proton candidates
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        const auto protonCandidates = GetProtonCandidates(primaryCandidates, muonCandidates);
        if (protonCandidates.size() + 2 != primaryCandidates.size())
            continue;
        
        counter.AddEventPassingCut("twoNonProtons", topology, eventIndex);
        
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        // Insist that muon and pion opening angle isn't too large
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        std::vector<TVector3> directions;
        for (const auto &particle : primaryCandidates)
        {
            if (ContainsParticle(protonCandidates, particle))
                continue;

            directions.push_back(particle.direction);
        }
        
        if (directions.size() != 2)
            throw std::logic_error("More than 2 mu/pi candidates selected. Something has gone wrong");

        const auto theta = std::acos(directions.front().Dot(directions.back()));
        if (theta > maxOpeningAngle)
            continue;

        counter.AddEventPassingCut("notBackToBack", topology, eventIndex);
        
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        // Insist that there aren't shower like primary looking particles we haven't considered as primaries
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        const auto secondaryCandidates = GetSecondaryCandidates(allParticles, primaryCandidates, secondaryDist);
        const auto crCandidates = GetCRCandidates(allParticles, primaryCandidates, recoNuVtx, crDist);

        float nShowers = 0;
        for (const auto &particle : allParticles)
        {
            if (ContainsParticle(primaryCandidates, particle))
                continue;
            
            if (ContainsParticle(secondaryCandidates, particle))
                continue;
            
            if (ContainsParticle(crCandidates, particle))
                continue;

            if (particle.trackShower > 0.5)
                continue;

            nShowers++;
        }

        if (nShowers > 0)
            continue;
        
        counter.AddEventPassingCut("showerCut", topology, eventIndex);

        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        // Insist that the pandora topological neutrino score isn't very low
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        /*
        if (em.topologicalScore < minTopologicalScore)
            continue;
        
        counter.AddEventPassingCut("topological", topology, eventIndex);
        */
    }

    counter.PrintBreakdown();
    
    // Now go through the events again and make the plots. This could be done more efficiently, by making the plots at the same time as
    // doing the selection, but this is good enough for me.
    /* 
    const auto pi = std::acos(0) * 2;
    EventPlot trueNuEPlot("True Neutrino Energy", "Neutrino Energy [GeV]", 100, 0, 4);
    EventPlot trueMuEnergyPlot("True Muon Energy", "Muon Energy [GeV]", 100, 0, 3);
    EventPlot trueMuThetaPlot("True Muon Theta", "Muon Theta [rad]", 100, 0, pi);
    EventPlot trueMuPhiPlot("True Muon Phi", "Muon Phi [rad]", 100, -pi, pi);
    EventPlot truePiEnergyPlot("True Pion Energy", "Pion Energy [GeV]", 100, 0, 3);
    EventPlot truePiThetaPlot("True Pion Theta", "Pion Theta [rad]", 100, 0, pi);
    EventPlot truePiPhiPlot("True Pion Phi", "Pion Phi [rad]", 100, -pi, pi);
    EventPlot trueMuPiAnglePlot("True Muon-Pion Opening Angle", "Opening angle [rad]", 100, 0, pi);
    EventPlot trueNProtonPlot("Proton Multiplicity", "Number of protons", 6, 0, 6);

    EventPlot pionCompletenessPlot("Pion Completeness", "Pion completeness", 100, 0, 1);
    EventPlot muonCompletenessPlot("Muon Completeness", "Muon completeness", 100, 0, 1);

    // For speed make the mapping from event index to the cuts passed for the signal events
    // TODO make this cleaner
    std::map<unsigned int, std::vector<std::string> > eventIndexToCutMap;
    for (const auto &cut : counter.GetCuts())
    {
        const auto selectedIndices = counter.GetSignalEventIndices(cut);
        for (const auto index : selectedIndices)
            eventIndexToCutMap[index].push_back(cut);
    }
   
    // Now fill the plots
    for (const auto entry : eventIndexToCutMap)
    {
        const auto eventIndex = entry.first;
        const auto cutsPassed = entry.second;
        em.LoadEvent(eventIndex);

        for (const auto &cut : cutsPassed)
        {
            trueNuEPlot.Fill(em.trueNuE, cut);
            trueMuEnergyPlot.Fill(em.trueMuEnergy, cut);
            trueMuThetaPlot.Fill(em.trueMuTheta, cut);
            trueMuPhiPlot.Fill(em.trueMuPhi, cut);
            truePiEnergyPlot.Fill(em.truePiEnergy, cut);
            truePiThetaPlot.Fill(em.truePiTheta, cut);
            truePiPhiPlot.Fill(em.truePiPhi, cut);
            trueMuPiAnglePlot.Fill(em.trueMuPiAngle, cut);
            trueNProtonPlot.Fill(em.nProton, cut);
    
            // Sum the completenesses of the best matched particles
            float pionCompleteness = 0.f;
            float muonCompleteness = 0.f;

            for (const auto &particle : em.GetParticles())
            {
                // Add up the pion completeness
                if (particle.hasMatchedMCParticle && particle.truePdgCode == 211)
                    pionCompleteness += particle.truthMatchCompleteness;
                
                // Add up the muon completeness
                if (particle.hasMatchedMCParticle && particle.truePdgCode == 13)
                    muonCompleteness += particle.truthMatchCompleteness;
            }

            pionCompletenessPlot.Fill(pionCompleteness, cut);
            muonCompletenessPlot.Fill(muonCompleteness, cut);
        }
    }
   
    const auto allCuts = counter.GetCuts();
    trueNuEPlot.Draw(allCuts);
    trueMuEnergyPlot.Draw(allCuts);
    trueMuThetaPlot.Draw(allCuts);
    trueMuPhiPlot.Draw(allCuts);
    truePiEnergyPlot.Draw(allCuts);
    truePiThetaPlot.Draw(allCuts);
    truePiPhiPlot.Draw(allCuts);
    trueMuPiAnglePlot.Draw(allCuts);
    trueNProtonPlot.Draw(allCuts);
    pionCompletenessPlot.Draw(allCuts);
    muonCompletenessPlot.Draw(allCuts);
    */
}
