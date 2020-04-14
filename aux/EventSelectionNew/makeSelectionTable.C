namespace cc1pievsel
{

/**
 *  @brief  Particle class holding all the details of a PFParticle
 */
class Particle
{
    public:
        unsigned int index;

        // True variables
        bool     t_hasMatchedMCParticle;
        int      t_pdgCode;
        bool     t_isGolden;
        float    t_truthMatchPurity;
        float    t_truthMatchCompleteness;
        float    t_range;
        TVector3 t_momentum;

        // Reco spectator variables (aren't used in the selection, but we want to look at them afterward)
        float    r_range;
    
        // Reco feature variables
        TVector3 r_start;
        TVector3 r_direction;
        bool     r_isContained;

        float    r_forwardLikelihood;
        float    r_logBragg_pToMIP;
        float    r_logBragg_piToMIP;
        int      r_nDescendents;
        int      r_nDownstreamHits;
        int      r_nSpacePointsInSphere5;
        float    r_rmsTrackDeviation;
        float    r_trackShower;
        float    r_dEdxTruncMeanStart;

        bool     r_areFeaturesAvailable;
};

// =========================================================================================================================================

/**
 *  @brief  Wrapper class for at root TMVA BDT
 */
class BDT
{
    public:
        /**
         *  @brief  Constructor
         *
         *  @param  weightsFile the path to the trained BDT weights file
         */
        BDT(const std::string &weightsFile);

        /**
         *  @brief  Destructor
         */
        ~BDT();

        /**
         *  @brief  Get the BDT response for the input particle
         *
         *  @param  particle the input particle
         *
         *  @return the BDT response
         */
        float GetResponse(const Particle &particle);

    private:
        TMVA::Reader *m_reader = nullptr; ///< The TMVA reader
        std::string   m_weightsFile;      ///< The weights file
        
        // The parameters to feed to the BDT, must all be float
        float    r_logBragg_pToMIP;
        float    r_logBragg_piToMIP;
        float    r_nDescendents;
        float    r_nDownstreamHits;
        float    r_nSpacePointsInSphere5;
        float    r_rmsTrackDeviation;
        float    r_trackShower;

        // The dummy spectator variables, again must all be float
        float    t_hasMatchedMCParticleDummy;
        float    t_pdgCodeDummy;
        float    t_isGoldenDummy;
        float    t_momentumDummy;
};

// ------------------------------------------------------------------------------------------------------------------------------------------
// ------------------------------------------------------------------------------------------------------------------------------------------

BDT::BDT(const std::string &weightsFile) :
    m_weightsFile(weightsFile),
    t_hasMatchedMCParticleDummy(-std::numeric_limits<float>::max()),
    t_pdgCodeDummy(-std::numeric_limits<float>::max()),
    t_isGoldenDummy(-std::numeric_limits<float>::max()),
    t_momentumDummy(-std::numeric_limits<float>::max())
{
    TMVA::Tools::Instance();
    m_reader = new TMVA::Reader("Silent");

    // Set up the variables that were used for training
    m_reader->AddVariable("r_logBragg_pToMIP", &r_logBragg_pToMIP);
    m_reader->AddVariable("r_logBragg_piToMIP", &r_logBragg_piToMIP);
    m_reader->AddVariable("r_nDescendents", &r_nDescendents);
    m_reader->AddVariable("r_nDownstreamHits", &r_nDownstreamHits);
    m_reader->AddVariable("r_nSpacePointsInSphere5", &r_nSpacePointsInSphere5);
    m_reader->AddVariable("r_rmsTrackDeviation", &r_rmsTrackDeviation);
    m_reader->AddVariable("r_trackShower", &r_trackShower);
    
    // Add the spectators (not used but have to be there for TMVA to run)
    m_reader->AddSpectator("t_hasMatchedMCParticle", &t_hasMatchedMCParticleDummy);
    m_reader->AddSpectator("t_pdgCode", &t_pdgCodeDummy);
    m_reader->AddSpectator("t_isGolden", &t_isGoldenDummy);
    m_reader->AddSpectator("t_momentum", &t_momentumDummy);

    // Hook up the trained weights
    m_reader->BookMVA("BDT", weightsFile.c_str());
}

// ------------------------------------------------------------------------------------------------------------------------------------------

BDT::~BDT()
{
    delete m_reader;
}

// ------------------------------------------------------------------------------------------------------------------------------------------

float BDT::GetResponse(const Particle &particle)
{
    r_logBragg_pToMIP = static_cast<float>(particle.r_logBragg_pToMIP);
    r_logBragg_piToMIP = static_cast<float>(particle.r_logBragg_piToMIP);
    r_nDescendents = static_cast<float>(particle.r_nDescendents);
    r_nDownstreamHits = static_cast<float>(particle.r_nDownstreamHits);
    r_nSpacePointsInSphere5 = static_cast<float>(particle.r_nSpacePointsInSphere5);
    r_rmsTrackDeviation = static_cast<float>(particle.r_rmsTrackDeviation);
    r_trackShower = static_cast<float>(particle.r_trackShower);

    return m_reader->EvaluateMVA("BDT");
}

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
         *  @brief  Destructor
         */
        ~EventManager();

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
         *  @brief  Calculate the truncated mean dEdx feature
         *
         *  @param  dedxPerHit the input dEdx per hit
         *  @param  residualRangePerHit the input residual ranges per hit
         *
         *  @return the truncated mean dEdx
         */
        float GetTruncatedMeandEdx(const std::vector<float> dedxPerHit, const std::vector<float> &residualRangePerHit) const;
        
        /**
         *  @brief  Fill the event tree
         *
         *  @param  goldenPionBDT the golden pion BDR
         *  @param  protonBDT the proton BDT
         *  @param  muonBDT the muon BDT
         *  @param  isSelected if the event has been selected
         *  @param  lastStage the name of the last stage obtained by the event
         *  @param  lastStageIndex the index of the last stage obtained by the event
         *  @param  goldenPionId the reco golden pion ID
         *  @param  muonId the reco muon ID
         *  @param  protonIds the reco proton IDs
         */
        void FillEventTree(BDT &goldenPionBDT, BDT &protonBDT, BDT &muonBDT, const bool isSelected, const std::string &lastStage, const unsigned int lastStageIndex, const unsigned int goldenPionId, const unsigned int muonId, const std::vector<unsigned int> &protonIds);
    
        /**
         *  @brief  Fill the particle tree with the particles in the current event
         *
         *  @param  goldenPionBDT the golden pion BDR
         *  @param  protonBDT the proton BDT
         *  @param  muonBDT the muon BDT
         *  @param  isSelected if the event has been selected
         *  @param  lastStage the name of the last stage obtained by the event
         *  @param  lastStageIndex the index of the last stage obtained by the event
         *  @param  goldenPionId the reco golden pion ID
         *  @param  muonId the reco muon ID
         *  @param  protonIds the reco proton IDs
         */
        void FillParticleTree(BDT &goldenPionBDT, BDT &protonBDT, BDT &muonBDT, const bool isSelected, const std::string &lastStage, const unsigned int lastStageIndex, const unsigned int goldenPionId, const unsigned int muonId, const std::vector<unsigned int> &protonIds);

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
//        float topologicalScore = -std::numeric_limits<float>::max();
        bool isRecoNuFiducial = false;
        int nFinalStatePFPs = -std::numeric_limits<int>::max();
        TVector3 *recoNuVtx = nullptr;

        // The output event information
        bool        outputEvent_isSignal = false;
        std::string outputEvent_topology;

        int         outputEvent_lastStageIndex = -std::numeric_limits<int>::max();
        std::string outputEvent_lastStage;

        float       outputEvent_muonCandidate_goldenPionBDTResponse = -std::numeric_limits<float>::max();
        float       outputEvent_muonCandidate_muonBDTResponse = -std::numeric_limits<float>::max();
        float       outputEvent_muonCandidate_protonBDTResponse = -std::numeric_limits<float>::max();
        
        float       outputEvent_goldenPionCandidate_goldenPionBDTResponse = -std::numeric_limits<float>::max();
        float       outputEvent_goldenPionCandidate_muonBDTResponse = -std::numeric_limits<float>::max();
        float       outputEvent_goldenPionCandidate_protonBDTResponse = -std::numeric_limits<float>::max();

    private:
        /**
         *  @brief  Load the input tree from the file
         */
        void LoadInputTree();

        /**
         *  @brief  Set the branch addresses
         */
        void SetBranchAddresses();
        
        /**
         *  @brief  Setup the output particle tree
         */
        void SetupOutputParticleTree();
        
        /**
         *  @brief  Setup the output event tree
         */
        void SetupOutputEventTree();

        // The input file details
        std::string m_fileName;
        std::string m_dirName;
        std::string m_treeName;

        // The tree from which we read the data
        TTree *m_tree;

        // Information about the current event
        unsigned int m_currentEventIndex;
        std::vector<Particle> m_particles;
        
        // True particle info
        std::vector<bool>  *hasMatchedMCParticleVect = nullptr;
        std::vector<int>   *truePdgCodeVect = nullptr;
        std::vector<bool>  *trueIsGoldenVect = nullptr;
        std::vector<float> *truthMatchPurityVect = nullptr;
        std::vector<float> *truthMatchCompletenessVect = nullptr;
        std::vector<float> *trueRangeVect = nullptr;
        std::vector<float> *trueMomentumXVect = nullptr;
        std::vector<float> *trueMomentumYVect = nullptr;
        std::vector<float> *trueMomentumZVect = nullptr;

        // Reco particle info
        std::vector<bool>  *isContainedVect = nullptr;
        std::vector<float> *rangeVect = nullptr;
        std::vector<float> *startXVect = nullptr;
        std::vector<float> *startYVect = nullptr;
        std::vector<float> *startZVect = nullptr;
        std::vector<float> *directionXVect = nullptr;
        std::vector<float> *directionYVect = nullptr;
        std::vector<float> *directionZVect = nullptr;
        std::vector<float> *yzAngleVect = nullptr;
        std::vector<float> *trackShowerVect = nullptr;
        std::vector<float> *braggpWVect = nullptr;
        std::vector<float> *braggpBackwardWVect = nullptr;
        std::vector<float> *braggpiWVect = nullptr;
        std::vector<float> *braggMIPWVect = nullptr;
        std::vector<float> *braggpUVVect = nullptr;
        std::vector<float> *braggpBackwardUVVect = nullptr;
        std::vector<float> *braggpiUVVect = nullptr;
        std::vector<float> *braggMIPUVVect = nullptr; 
        std::vector<int>   *nDescendentsVect = nullptr;
        std::vector<int>   *nHitsUVect = nullptr;
        std::vector<int>   *nHitsVVect = nullptr;
        std::vector<int>   *nHitsWVect = nullptr;
        std::vector<int>   *nDescendentHitsUVect = nullptr;
        std::vector<int>   *nDescendentHitsVVect = nullptr;
        std::vector<int>   *nDescendentHitsWVect = nullptr;
        std::vector<int>   *nSpacePointsInSphere5Vect = nullptr;
        std::vector<float> *rmsTrackDeviationVect = nullptr;

        std::vector<std::vector<float> > *dedxPerHitUVect = nullptr;
        std::vector<std::vector<float> > *dedxPerHitVVect = nullptr;
        std::vector<std::vector<float> > *dedxPerHitWVect = nullptr;
        std::vector<std::vector<float> > *residualRangePerHitUVect = nullptr;
        std::vector<std::vector<float> > *residualRangePerHitVVect = nullptr;
        std::vector<std::vector<float> > *residualRangePerHitWVect = nullptr;

        // The output particle for filling the output tree
        TFile      *m_outputFile = nullptr;
        TTree      *m_outputParticleTree = nullptr;;
        TTree      *m_outputEventTree = nullptr;;
        Particle    m_outputParticle;
        std::string m_outputParticle_lastStage;
        int         m_outputParticle_lastStageIndex;
        float       m_outputParticle_goldenPionBDTResponse;
        float       m_outputParticle_muonBDTResponse;
        float       m_outputParticle_protonBDTResponse;
        bool        m_outputParticle_isMaxGoldenPionBDTResponse;
        float       m_outputParticle_maxGoldenPionBDTResponse;
        bool        m_outputParticle_isEventSelected;
        bool        m_outputParticle_isGoldenPionCandidate;
        bool        m_outputParticle_isMuonCandidate;
        bool        m_outputParticle_isProtonCandidate;
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
   
    this->SetupOutputParticleTree();
    this->SetupOutputEventTree();

    this->LoadEvent(0);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

EventManager::~EventManager()
{
    m_outputFile->Write();

    delete m_outputParticleTree;
    delete m_outputEventTree;
    delete m_outputFile;
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
    // Event level info
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
//    m_tree->SetBranchAddress("topologicalScore", &topologicalScore);
    m_tree->SetBranchAddress("isRecoNuFiducial", &isRecoNuFiducial);
    m_tree->SetBranchAddress("nFinalStatePFPs", &nFinalStatePFPs);
    m_tree->SetBranchAddress("recoNuVtx", &recoNuVtx);

    // Particle level info
    m_tree->SetBranchAddress("hasMatchedMCParticleVect", &hasMatchedMCParticleVect);
    m_tree->SetBranchAddress("truePdgCodeVect", &truePdgCodeVect);
    m_tree->SetBranchAddress("trueIsGoldenVect", &trueIsGoldenVect);
    m_tree->SetBranchAddress("truthMatchPurityVect", &truthMatchPurityVect);
    m_tree->SetBranchAddress("truthMatchCompletenessVect", &truthMatchCompletenessVect);
    m_tree->SetBranchAddress("trueRangeVect", &trueRangeVect);
    m_tree->SetBranchAddress("trueMomentumXVect", &trueMomentumXVect);
    m_tree->SetBranchAddress("trueMomentumYVect", &trueMomentumYVect);
    m_tree->SetBranchAddress("trueMomentumZVect", &trueMomentumZVect);

    m_tree->SetBranchAddress("rangeVect", &rangeVect);
    m_tree->SetBranchAddress("startXVect", &startXVect);
    m_tree->SetBranchAddress("startYVect", &startYVect);
    m_tree->SetBranchAddress("startZVect", &startZVect);
    m_tree->SetBranchAddress("directionXVect", &directionXVect);
    m_tree->SetBranchAddress("directionYVect", &directionYVect);
    m_tree->SetBranchAddress("directionZVect", &directionZVect);
    m_tree->SetBranchAddress("isContainedVect", &isContainedVect);
    m_tree->SetBranchAddress("yzAngleVect", &yzAngleVect);
    m_tree->SetBranchAddress("trackShowerVect", &trackShowerVect);
    m_tree->SetBranchAddress("braggpWVect", &braggpWVect);
    m_tree->SetBranchAddress("braggpBackwardWVect", &braggpBackwardWVect);
    m_tree->SetBranchAddress("braggpiWVect", &braggpiWVect);
    m_tree->SetBranchAddress("braggMIPWVect", &braggMIPWVect);
    m_tree->SetBranchAddress("braggpUVVect", &braggpUVVect);
    m_tree->SetBranchAddress("braggpBackwardUVVect", &braggpBackwardUVVect);
    m_tree->SetBranchAddress("braggpiUVVect", &braggpiUVVect);
    m_tree->SetBranchAddress("braggMIPUVVect", &braggMIPUVVect);
    m_tree->SetBranchAddress("nDescendentsVect", &nDescendentsVect);
    m_tree->SetBranchAddress("nHitsUVect", &nHitsUVect);
    m_tree->SetBranchAddress("nHitsVVect", &nHitsVVect);
    m_tree->SetBranchAddress("nHitsWVect", &nHitsWVect);
    m_tree->SetBranchAddress("nDescendentHitsUVect", &nDescendentHitsUVect);
    m_tree->SetBranchAddress("nDescendentHitsVVect", &nDescendentHitsVVect);
    m_tree->SetBranchAddress("nDescendentHitsWVect", &nDescendentHitsWVect);
    m_tree->SetBranchAddress("nSpacePointsInSphere5Vect", &nSpacePointsInSphere5Vect);
    m_tree->SetBranchAddress("rmsTrackDeviationVect", &rmsTrackDeviationVect);

    m_tree->SetBranchAddress("dedxPerHitUVect", &dedxPerHitUVect);
    m_tree->SetBranchAddress("dedxPerHitVVect", &dedxPerHitVVect);
    m_tree->SetBranchAddress("dedxPerHitWVect", &dedxPerHitWVect);
    m_tree->SetBranchAddress("residualRangePerHitUVect", &residualRangePerHitUVect);
    m_tree->SetBranchAddress("residualRangePerHitVVect", &residualRangePerHitVVect);
    m_tree->SetBranchAddress("residualRangePerHitWVect", &residualRangePerHitWVect);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void EventManager::LoadEvent(const unsigned int eventIndex)
{
    m_tree->GetEntry(eventIndex);
    m_currentEventIndex = eventIndex;
    
    // Read the particle-level vectors
    bool hasRecoGoldenPion = false;
    m_particles.clear();
    for (unsigned int index = 0; index < nFinalStatePFPs; ++index)
    {
        Particle p;
        p.index = index;
        
        p.t_hasMatchedMCParticle = hasMatchedMCParticleVect->at(index);
        p.t_pdgCode = truePdgCodeVect->at(index);
        p.t_isGolden = trueIsGoldenVect->at(index);
        p.t_truthMatchPurity = truthMatchPurityVect->at(index);
        p.t_truthMatchCompleteness = truthMatchCompletenessVect->at(index);
        p.t_range = trueRangeVect->at(index);
        p.t_momentum = TVector3(trueMomentumXVect->at(index), trueMomentumYVect->at(index), trueMomentumZVect->at(index));
        p.r_range = rangeVect->at(index); 
        p.r_start = TVector3(startXVect->at(index), startYVect->at(index), startZVect->at(index));
        p.r_direction = TVector3(directionXVect->at(index), directionYVect->at(index), directionZVect->at(index));
        p.r_isContained = isContainedVect->at(index);

        // Where required calculate the features from the input tree
        const auto yzAngle = yzAngleVect->at(index);
        const auto isTrackAlongWWire = (std::pow(std::sin(yzAngle), 2) < 0.175);

        const auto braggpW = braggpWVect->at(index);
        const auto braggpBackwardW = braggpBackwardWVect->at(index);
        const auto braggpiW = braggpiWVect->at(index);
        const auto braggMIPW = braggMIPWVect->at(index);
        
        const auto braggpUV = braggpUVVect->at(index);
        const auto braggpBackwardUV = braggpBackwardUVVect->at(index);
        const auto braggpiUV = braggpiUVVect->at(index);
        const auto braggMIPUV = braggMIPUVVect->at(index);

        const auto braggp = (!isTrackAlongWWire && braggpW > -1) ? braggpW : braggpUV;
        const auto braggpBackward = (!isTrackAlongWWire && braggpBackwardW > -1) ? braggpBackwardW : braggpBackwardUV;
        const auto braggpi = (!isTrackAlongWWire && braggpiW > -1) ? braggpiW : braggpiUV;
        const auto braggMIP = (!isTrackAlongWWire && braggMIPW > -1) ? braggMIPW : braggMIPUV;

        const auto forwardLikelihoodDenom = std::exp(braggp) + std::exp(braggpBackward);
        const auto isForwardLikelihoodAvailable = braggp > -1 && braggpBackward > -1 && std::abs(forwardLikelihoodDenom) > std::numeric_limits<float>::epsilon();
        p.r_forwardLikelihood = isForwardLikelihoodAvailable ? (std::exp(braggp) / forwardLikelihoodDenom) : -std::numeric_limits<float>::max();

        const auto isLogBragg_pToMIPAvailable = braggp > -1 && braggMIP > -1 && std::abs(braggMIP) > std::numeric_limits<float>::epsilon();
        p.r_logBragg_pToMIP = isLogBragg_pToMIPAvailable ? (std::log(braggp / braggMIP)) : -std::numeric_limits<float>::max();
        
        const auto isLogBragg_piToMIPAvailable = braggpi > -1 && braggMIP > -1 && std::abs(braggMIP) > std::numeric_limits<float>::epsilon();
        p.r_logBragg_piToMIP = isLogBragg_piToMIPAvailable ? (std::log(braggpi / braggMIP)) : -std::numeric_limits<float>::max();

        p.r_nDownstreamHits = (nDescendentHitsUVect->at(index) + nDescendentHitsVVect->at(index) + nDescendentHitsWVect->at(index))
                              - (nHitsUVect->at(index) + nHitsVVect->at(index) + nHitsWVect->at(index));

        p.r_nDescendents = nDescendentsVect->at(index);
        p.r_nSpacePointsInSphere5 = nSpacePointsInSphere5Vect->at(index);
        p.r_rmsTrackDeviation = rmsTrackDeviationVect->at(index);
        p.r_trackShower = trackShowerVect->at(index);

        // This should really be moved upstream to the analyser module, doing it here for now!
        // ATTN this is commented out for speed while it's not in use
        p.r_dEdxTruncMeanStart = 0.f;
        /*
        const auto dEdxTruncMeanStartW = this->GetTruncatedMeandEdx(dedxPerHitWVect->at(index), residualRangePerHitWVect->at(index));
        p.r_dEdxTruncMeanStart = (!isTrackAlongWWire && dEdxTruncMeanStartW > -1) ? dEdxTruncMeanStartW : -std::numeric_limits<float>::max();

        if (p.r_dEdxTruncMeanStart < -1)
        {
            const auto dEdxTruncMeanStartU = this->GetTruncatedMeandEdx(dedxPerHitUVect->at(index), residualRangePerHitUVect->at(index));
        
            const auto dEdxTruncMeanStartV = this->GetTruncatedMeandEdx(dedxPerHitVVect->at(index), residualRangePerHitVVect->at(index));
    
            const auto dEdxTruncMeanStartUWeight = (dEdxTruncMeanStartU > -1) ? dedxPerHitUVect->at(index).size() : 0;
            const auto dEdxTruncMeanStartVWeight = (dEdxTruncMeanStartV > -1) ? dedxPerHitUVect->at(index).size() : 0;
            const auto dEdxTruncMeanStartUVWeight = dEdxTruncMeanStartUWeight + dEdxTruncMeanStartVWeight;

            if (dEdxTruncMeanStartUVWeight > 0)
            {
                p.r_dEdxTruncMeanStart = (dEdxTruncMeanStartUWeight * dEdxTruncMeanStartU + dEdxTruncMeanStartVWeight * dEdxTruncMeanStartV) / static_cast<float>(dEdxTruncMeanStartUVWeight);
            }
        }
        */

        // For all features other than track-shower, invalid entries are at -floatMax. For track-shower the invalid entry is at -1
        p.r_areFeaturesAvailable = isLogBragg_pToMIPAvailable &&
                                   isLogBragg_piToMIPAvailable &&
                                   /*isForwardLikelihoodAvailable &&*/
                                   p.r_nDownstreamHits > -1 &&
                                   p.r_nDescendents > -1 &&
                                   p.r_nSpacePointsInSphere5 > -1 &&
                                   p.r_rmsTrackDeviation > -1 &&
                                   p.r_trackShower > -0.5; // &&
                                   /*p.r_dEdxTruncMeanStart > -1;*/

        hasRecoGoldenPion = hasRecoGoldenPion || (p.t_pdgCode == 211 && p.t_isGolden);

        m_particles.push_back(p);
    }
  
    // TODO fix below
    //outputEvent_isSignal = isSignal && hasRecoGoldenPion; // ATTN this should be done usign truth in the grid jobs, this is a hack for now
    outputEvent_isSignal = isSignal && hasRecoGoldenPion && nProton > 0; // ATTN this should be done usign truth in the grid jobs, this is a hack for now
}

// -----------------------------------------------------------------------------------------------------------------------------------------

float EventManager::GetTruncatedMeandEdx(const std::vector<float> dedxPerHit, const std::vector<float> &residualRangePerHit) const
{
    const unsigned int nHitsToSkip = 3;
    const float lengthFraction = 1.f / 3.f; // The fraction of the track length to use

    if (dedxPerHit.size() != residualRangePerHit.size())
        throw std::invalid_argument("dEdx and residual range vectors have different sizes");

    // Check if the variable is calculable
    if (residualRangePerHit.size() <= nHitsToSkip)
        return -std::numeric_limits<float>::max();

    // Make make the vector of pairs to keep track of the incides
    std::vector<std::pair<float, unsigned int> > residualRangeIndices;
    float maxResidualRange = -std::numeric_limits<float>::max();
    for (unsigned int i = 0; i < residualRangePerHit.size(); ++i)
    {
        const auto residualRange = residualRangePerHit.at(i);
        maxResidualRange = std::max(maxResidualRange, residualRange);

        residualRangeIndices.emplace_back(residualRange, i);
    }
    const auto residualRangeCutoff = maxResidualRange * lengthFraction;

    // Sort the residual ranges such that the largest residual range (closest to the start of the track) is first
    std::sort(residualRangeIndices.begin(), residualRangeIndices.end(), [](auto &a, auto &b) {
        return a.first > b.first;
    });

    // Get the dEdx of the hits at the start of the track
    std::vector<float> dedxPerHitAtStart;
    for (unsigned int i = nHitsToSkip; i < residualRangeIndices.size(); ++i)
    {
        const auto entry = residualRangeIndices.at(i);
        const auto residualRange = entry.first;
        const auto hitIndex = entry.second;

        // ATTN small residual ranges are at the start of the track
        if (residualRange < residualRangeCutoff)
            continue;

        dedxPerHitAtStart.push_back(dedxPerHit.at(hitIndex));
    }

    const auto nHits = dedxPerHitAtStart.size();
    if (nHits == 0)
        return -std::numeric_limits<float>::max();

    // Sort the dEdx so we can find the median
    std::sort(dedxPerHitAtStart.begin(), dedxPerHitAtStart.end());
    const auto median = dedxPerHitAtStart.at(nHits / 2);
   
    // Now find the mean
    float total = 0.f;
    for (const auto &dEdx : dedxPerHitAtStart)
        total += dEdx;

    const auto mean = total / static_cast<float>(nHits);
    
    // Now find the variance
    float squareSum = 0.f;
    for (const auto &dEdx : dedxPerHitAtStart)
        squareSum += std::pow(dEdx - mean, 2);

    const auto variance = squareSum / static_cast<float>(nHits);
    
    // Get the mean dEdx of the hits within one standard deviation of the median
    float truncatedTotal = 0.f;
    unsigned int nTruncatedHits = 0;
    for (const auto &dEdx : dedxPerHitAtStart)
    {
        if (std::pow(dEdx - median, 2) > variance)
            continue;
        
        truncatedTotal += dEdx;
        nTruncatedHits++;
    }
    
    if (nTruncatedHits == 0)
        return -std::numeric_limits<float>::max();
    
    const auto truncatedMean = truncatedTotal / static_cast<float>(nTruncatedHits);

    return truncatedMean;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void EventManager::SetupOutputEventTree()
{
    if (!m_outputFile)
        throw std::logic_error("Can't set up output event tree as the output file has yet to be opened");

    m_outputEventTree = new TTree("events", "");

    m_outputEventTree->Branch("run", &run);
    m_outputEventTree->Branch("subRun", &subRun);
    m_outputEventTree->Branch("event", &event);
    m_outputEventTree->Branch("t_topology", &outputEvent_topology);
    m_outputEventTree->Branch("t_isCC1Pi", &isSignal);
    m_outputEventTree->Branch("t_isSignal", &outputEvent_isSignal);
    m_outputEventTree->Branch("t_isNuFiducial", &isTrueNuFiducial);
    m_outputEventTree->Branch("t_nuE", &trueNuE);
    m_outputEventTree->Branch("t_muEnergy", &trueMuEnergy);
    m_outputEventTree->Branch("t_muTheta", &trueMuTheta);
    m_outputEventTree->Branch("t_muPhi", &trueMuPhi);
    m_outputEventTree->Branch("t_piEnergy", &truePiEnergy);
    m_outputEventTree->Branch("t_piTheta", &truePiTheta);
    m_outputEventTree->Branch("t_piPhi", &truePiPhi);
    m_outputEventTree->Branch("t_muPiAngle", &trueMuPiAngle);
    m_outputEventTree->Branch("t_nMuMinus", &nMuMinus);
    m_outputEventTree->Branch("t_nMuPlus", &nMuPlus);
    m_outputEventTree->Branch("t_nPiPlus", &nPiPlus);
    m_outputEventTree->Branch("t_nPiMinus", &nPiMinus);
    m_outputEventTree->Branch("t_nKPlus", &nKPlus);
    m_outputEventTree->Branch("t_nKMinus", &nKMinus);
    m_outputEventTree->Branch("t_nProton", &nProton);
    m_outputEventTree->Branch("t_nNeutron", &nNeutron);
    m_outputEventTree->Branch("t_nPhoton", &nPhoton);
    m_outputEventTree->Branch("t_nElectron", &nElectron);
    m_outputEventTree->Branch("t_nPositron", &nPositron);
    m_outputEventTree->Branch("t_nTotal", &nTotal);

    m_outputEventTree->Branch("r_isNuFiducial", &isRecoNuFiducial);
    m_outputEventTree->Branch("r_nFinalStatePFPs", &nFinalStatePFPs);
    m_outputEventTree->Branch("r_recoNuVtx", &recoNuVtx);

    // Selection info
    m_outputEventTree->Branch("r_lastStage", &outputEvent_lastStage);
    m_outputEventTree->Branch("r_lastStageIndex", &outputEvent_lastStageIndex);

    m_outputEventTree->Branch("r_muonCandidate_goldenPionBDTResponse", &outputEvent_muonCandidate_goldenPionBDTResponse);
    m_outputEventTree->Branch("r_muonCandidate_muonBDTResponse", &outputEvent_muonCandidate_muonBDTResponse);
    m_outputEventTree->Branch("r_muonCandidate_protonBDTResponse", &outputEvent_muonCandidate_protonBDTResponse);
    m_outputEventTree->Branch("r_goldenPionCandidate_goldenPionBDTResponse", &outputEvent_goldenPionCandidate_goldenPionBDTResponse);
    m_outputEventTree->Branch("r_goldenPionCandidate_muonBDTResponse", &outputEvent_goldenPionCandidate_muonBDTResponse);
    m_outputEventTree->Branch("r_goldenPionCandidate_protonBDTResponse", &outputEvent_goldenPionCandidate_protonBDTResponse);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void EventManager::SetupOutputParticleTree()
{
    m_outputFile = new TFile("eventSelectionOutput.root", "RECREATE");

    m_outputParticleTree = new TTree("particles", "");

    m_outputParticleTree->Branch("run", &run);
    m_outputParticleTree->Branch("subRun", &subRun);
    m_outputParticleTree->Branch("event", &event);
    m_outputParticleTree->Branch("isSignalEvent", &outputEvent_isSignal);
    m_outputParticleTree->Branch("isCC1PiEvent", &isSignal);
    m_outputParticleTree->Branch("eventTopology", &outputEvent_topology);
    m_outputParticleTree->Branch("t_hasMatchedMCParticle", &m_outputParticle.t_hasMatchedMCParticle);
    m_outputParticleTree->Branch("t_pdgCode", &m_outputParticle.t_pdgCode);
    m_outputParticleTree->Branch("t_isGolden", &m_outputParticle.t_isGolden);
    m_outputParticleTree->Branch("t_truthMatchPurity", &m_outputParticle.t_truthMatchPurity);
    m_outputParticleTree->Branch("t_truthMatchCompleteness", &m_outputParticle.t_truthMatchCompleteness);
    m_outputParticleTree->Branch("t_range", &m_outputParticle.t_range);
    m_outputParticleTree->Branch("t_momentum", &m_outputParticle.t_momentum);
    m_outputParticleTree->Branch("r_range", &m_outputParticle.r_range);
    m_outputParticleTree->Branch("r_start", &m_outputParticle.r_start);
    m_outputParticleTree->Branch("r_direction", &m_outputParticle.r_direction);
    m_outputParticleTree->Branch("r_isContained", &m_outputParticle.r_isContained);
    m_outputParticleTree->Branch("r_forwardLikelihood", &m_outputParticle.r_forwardLikelihood);
    m_outputParticleTree->Branch("r_logBragg_pToMIP", &m_outputParticle.r_logBragg_pToMIP);
    m_outputParticleTree->Branch("r_logBragg_piToMIP", &m_outputParticle.r_logBragg_piToMIP);
    m_outputParticleTree->Branch("r_nDescendents", &m_outputParticle.r_nDescendents);
    m_outputParticleTree->Branch("r_nDownstreamHits", &m_outputParticle.r_nDownstreamHits);
    m_outputParticleTree->Branch("r_nSpacePointsInSphere5", &m_outputParticle.r_nSpacePointsInSphere5);
    m_outputParticleTree->Branch("r_rmsTrackDeviation", &m_outputParticle.r_rmsTrackDeviation);
    m_outputParticleTree->Branch("r_trackShower", &m_outputParticle.r_trackShower);
    m_outputParticleTree->Branch("r_dEdxTruncMeanStart", &m_outputParticle.r_dEdxTruncMeanStart);
    m_outputParticleTree->Branch("r_areFeaturesAvailable", &m_outputParticle.r_areFeaturesAvailable);
    m_outputParticleTree->Branch("r_lastStage", &m_outputParticle_lastStage);
    m_outputParticleTree->Branch("r_lastStageIndex", &m_outputParticle_lastStageIndex);
    m_outputParticleTree->Branch("r_goldenPionBDTResponse", &m_outputParticle_goldenPionBDTResponse);
    m_outputParticleTree->Branch("r_muonBDTResponse", &m_outputParticle_muonBDTResponse);
    m_outputParticleTree->Branch("r_protonBDTResponse", &m_outputParticle_protonBDTResponse);
    m_outputParticleTree->Branch("r_isMaxGoldenPionBDTResponse", &m_outputParticle_isMaxGoldenPionBDTResponse);
    m_outputParticleTree->Branch("r_maxGoldenPionBDTResponse", &m_outputParticle_maxGoldenPionBDTResponse);
    m_outputParticleTree->Branch("r_isEventSelected", &m_outputParticle_isEventSelected);
    m_outputParticleTree->Branch("r_isGoldenPionCandidate", &m_outputParticle_isGoldenPionCandidate);
    m_outputParticleTree->Branch("r_isMuonCandidate", &m_outputParticle_isMuonCandidate);
    m_outputParticleTree->Branch("r_isProtonCandidate", &m_outputParticle_isProtonCandidate);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void EventManager::FillEventTree(BDT &goldenPionBDT, BDT &protonBDT, BDT &muonBDT, const bool isSelected, const std::string &lastStage, const unsigned int lastStageIndex, const unsigned int goldenPionId, const unsigned int muonId, const std::vector<unsigned int> &protonIds)
{
    outputEvent_lastStage = lastStage;
    outputEvent_lastStageIndex = lastStageIndex;
    outputEvent_topology = this->GetTopologyString();
    
    outputEvent_muonCandidate_goldenPionBDTResponse = -std::numeric_limits<float>::max();
    outputEvent_muonCandidate_muonBDTResponse = -std::numeric_limits<float>::max();
    outputEvent_muonCandidate_protonBDTResponse = -std::numeric_limits<float>::max();

    outputEvent_goldenPionCandidate_goldenPionBDTResponse = -std::numeric_limits<float>::max();
    outputEvent_goldenPionCandidate_muonBDTResponse = -std::numeric_limits<float>::max();
    outputEvent_goldenPionCandidate_protonBDTResponse = -std::numeric_limits<float>::max();

    const auto particles = this->GetParticles();
    for (unsigned int i = 0; i < particles.size(); ++i)
    {
        const auto particle = particles.at(i);
        
        if (particle.index != goldenPionId && particle.index != muonId)
            continue;

        const auto goldenPionBDTResponse = goldenPionBDT.GetResponse(particle);
        const auto muonBDTResponse = muonBDT.GetResponse(particle);
        const auto protonBDTResponse = protonBDT.GetResponse(particle);

        if (particle.index == goldenPionId)
        {
            outputEvent_goldenPionCandidate_goldenPionBDTResponse = goldenPionBDTResponse;
            outputEvent_goldenPionCandidate_muonBDTResponse = muonBDTResponse;
            outputEvent_goldenPionCandidate_protonBDTResponse = protonBDTResponse;
        }
        else if (particle.index == muonId)
        {
            outputEvent_muonCandidate_goldenPionBDTResponse = goldenPionBDTResponse;
            outputEvent_muonCandidate_muonBDTResponse = muonBDTResponse;
            outputEvent_muonCandidate_protonBDTResponse = protonBDTResponse;
        }
    }

    m_outputEventTree->Fill();
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void EventManager::FillParticleTree(BDT &goldenPionBDT, BDT &protonBDT, BDT &muonBDT, const bool isSelected, const std::string &lastStage, const unsigned int lastStageIndex, const unsigned int goldenPionId, const unsigned int muonId, const std::vector<unsigned int> &protonIds)
{
    outputEvent_topology = this->GetTopologyString();


    // Find the index of the particle with the maximum BDT response that is contained
    const auto particles = this->GetParticles();

    // Get the particle with the maximum golden pion BDT response
    unsigned int maxIndex = std::numeric_limits<unsigned int>::max();
    float maxResponse = -std::numeric_limits<float>::max();
    for (unsigned int i = 0; i < particles.size(); ++i)
    {
        if (!particles.at(i).r_isContained)
            continue;

        const auto response = goldenPionBDT.GetResponse(particles.at(i));
        if (response > maxResponse)
        {
            maxResponse = response;
            maxIndex = i;       
        }
    }

    m_outputParticle_isEventSelected = isSelected;
    for (unsigned int i = 0; i < particles.size(); ++i)
    {
        const auto particle = particles.at(i);

        m_outputParticle = particle;
        m_outputParticle_lastStage = lastStage;
        m_outputParticle_lastStageIndex = lastStageIndex;

        m_outputParticle_goldenPionBDTResponse = goldenPionBDT.GetResponse(m_outputParticle);
        m_outputParticle_muonBDTResponse = muonBDT.GetResponse(m_outputParticle);
        m_outputParticle_protonBDTResponse = protonBDT.GetResponse(m_outputParticle);
        m_outputParticle_isMaxGoldenPionBDTResponse = (i == maxIndex);
        m_outputParticle_maxGoldenPionBDTResponse = maxResponse;

        m_outputParticle_isGoldenPionCandidate = (particle.index == goldenPionId);
        m_outputParticle_isMuonCandidate = (particle.index == muonId);
        m_outputParticle_isProtonCandidate = false;
        for (const auto protonId : protonIds)
        {
            if (particle.index == protonId)
            {
                m_outputParticle_isProtonCandidate = true;
                break;
            }
        }

        m_outputParticleTree->Fill();
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

    if (outputEvent_isSignal)
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
    //interaction += "X p  ";
   
    if (nProton != 0)
        interaction += "N p  ";
    else
        interaction += "0 p  ";

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
        
void SelectionCounter::PrintBreakdown(const unsigned int nTopologies)
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

void SelectionCounter::PrintPerformance(const std::string &cut, const unsigned int nTopologies)
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
    std::cout << " - Purity * Efficiency  : " << purity * efficiency << std::endl;

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
// SELECTION STARTS HERE
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
        if ((particle.r_start - recoNuVtx).Mag2() < cut2)
            outputParticles.push_back(particle);
    }

    return outputParticles;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  Get the particles from the input vector that are contained
 *
 *  @param  particles the input vector of particles
 *  @param  containedParticles the contained particles to populate
 *  @param  uncontainedParticles uncontained particles to populate
 */
void GetContainedParticles(const std::vector<Particle> &particles, std::vector<Particle> &containedParticles, std::vector<Particle> &uncontainedParticles)
{
    for (const auto &particle : particles)
    {
        if (particle.r_isContained)
        {
            containedParticles.push_back(particle);
        }
        else
        {
            uncontainedParticles.push_back(particle);
        }
    }
}

} // namespace cc1pievsel

// -----------------------------------------------------------------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------------------------------------------------------------

using namespace cc1pievsel;

bool SelectEvent(const EventManager &em, const unsigned int eventIndex, BDT &goldenPionBDT, BDT &protonBDT, BDT &muonBDT, SelectionCounter &counter, std::string &lastStage, unsigned int &lastStageIndex, unsigned int &goldenPionId, unsigned int &muonId, std::vector<unsigned int> &protonIds)
{
    // TODO make configurable
    const float trackStartToVertexDist = 5.f; // The maximum distance from which a primary particle can start from the vertex

    // Keep track of the stages
    std::vector<std::string> stages;

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    const auto topology = em.GetTopologyString();

    // ...........................................................
    stages.push_back("all");

    lastStageIndex = stages.size() - 1;
    lastStage = std::to_string(lastStageIndex) + "_" + stages.back();
    counter.AddEventPassingCut(lastStage, topology, eventIndex);
    // ...........................................................

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Insist that the reconstructed neutrino is fiducial
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if (!em.isRecoNuFiducial)
        return false;

    // ...........................................................
    stages.push_back("fiducial");

    lastStageIndex = stages.size() - 1;
    lastStage = std::to_string(lastStageIndex) + "_" + stages.back();
    counter.AddEventPassingCut(lastStage, topology, eventIndex);
    // ...........................................................

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Insist that all particles start near the vertex
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    const auto allParticles = em.GetParticles();
    const auto recoNuVtx = *em.recoNuVtx;
    const auto particlesNearVertex = GetParticlesNearVertex(allParticles, recoNuVtx, trackStartToVertexDist);

    if (particlesNearVertex.size() != allParticles.size())
        return false;

    // ...........................................................
    stages.push_back("allPrimariesNearVertex");

    lastStageIndex = stages.size() - 1;
    lastStage = std::to_string(lastStageIndex) + "_" + stages.back();
    counter.AddEventPassingCut(lastStage, topology, eventIndex);
    // ...........................................................

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Insist that there are at least 2 primary PFParticles
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if (allParticles.size() < 2)
        return false;

    // ...........................................................
    stages.push_back("minTwoPrimaries");

    lastStageIndex = stages.size() - 1;
    lastStage = std::to_string(lastStageIndex) + "_" + stages.back();
    counter.AddEventPassingCut(lastStage, topology, eventIndex);
    // ...........................................................

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Insist as most one primary particle is uncontained
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    std::vector<Particle> containedParticles, uncontainedParticles;
    GetContainedParticles(allParticles, containedParticles, uncontainedParticles);
    if (uncontainedParticles.size() > 1)
        return false;

    // ...........................................................
    stages.push_back("maxOneUncontained");

    lastStageIndex = stages.size() - 1;
    lastStage = std::to_string(lastStageIndex) + "_" + stages.back();
    counter.AddEventPassingCut(lastStage, topology, eventIndex);
    // ...........................................................

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Insist that at all contained primary PFParticles have the features available to run the BDT
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    for (const auto &particle : containedParticles)
    {
        if (!particle.r_areFeaturesAvailable)
            return false;
    }
    
    // ...........................................................
    stages.push_back("allFeaturesAvailable");

    lastStageIndex = stages.size() - 1;
    lastStage = std::to_string(lastStageIndex) + "_" + stages.back();
    counter.AddEventPassingCut(lastStage, topology, eventIndex);
    // ...........................................................
    
   
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Insist that there are exactly 2 non-proton-like particles
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    // Get the golden pion candidate
    unsigned int goldenPionIndex = std::numeric_limits<unsigned int>::max();
    float maxResponse = -std::numeric_limits<float>::max();

    for (unsigned int i = 0; i < containedParticles.size(); ++i)
    {
        const auto particle = containedParticles.at(i);
        const auto response = goldenPionBDT.GetResponse(particle);

        if (response > maxResponse)
        {
            goldenPionIndex = i;
            maxResponse = response;
        }
    }

    // Separate out the contained particles that are not the golden pion candidate
    std::vector<Particle> containedNonGoldenPionParticles;
    for (unsigned int i = 0; i < containedParticles.size(); ++i)
    {
        if (i != goldenPionIndex)
            containedNonGoldenPionParticles.push_back(containedParticles.at(i));
    }

    // Separate out the proton candidates and the remaining particles (assumed to be muon candidates)
    std::vector<Particle> protonCandidates;
    std::vector<Particle> muonCandidates = uncontainedParticles;

    for (const auto &particle : containedNonGoldenPionParticles)
    {
        const auto response = protonBDT.GetResponse(particle);

        if (response > 0.2f)
        {
            protonCandidates.push_back(particle);
        }
        else
        {
            muonCandidates.push_back(particle);
        }
    }

    // Check that we have the right number of proton candidates
    if (protonCandidates.size() != allParticles.size() - 2)
        return false;

    // ...........................................................
    stages.push_back("correctProtonCount");

    lastStageIndex = stages.size() - 1;
    lastStage = std::to_string(lastStageIndex) + "_" + stages.back();
    counter.AddEventPassingCut(lastStage, topology, eventIndex);
    // ...........................................................
    
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Identify the particles, and insist that the golden pion has a large enough score
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if (muonCandidates.size() != 1)
        throw std::logic_error("Sanity check failed: There are " + std::to_string(muonCandidates.size()) + " muon candidates!");

    const auto goldenPion = containedParticles.at(goldenPionIndex);
    const auto muon = muonCandidates.front();
    const auto protons = protonCandidates;

    // Do some sanity checks before we move on
    if (protons.size() + 2 != allParticles.size())
        throw std::logic_error("Sanity check failted: Number of proton candidates + 2 isn't the number of input particles!");
    
    if (goldenPion.index == muon.index)
        throw std::logic_error("Sanity check failed: Muon and golden pion candidate are the same particle!");

    for (const auto &proton : protons)
    {
        if (proton.index == goldenPion.index)
            throw std::logic_error("Sanity check failed: One of the proton candidates is the golden pion candidate!");
        
        if (proton.index == muon.index)
            throw std::logic_error("Sanity check failed: One of the muon candidates is the golden pion candidate!");
    }

    // Save the output particle indices
    goldenPionId = goldenPion.index;
    muonId = muon.index;
    for (const auto &proton : protons)
        protonIds.push_back(proton.index);

    // Apply the golden pion cut
    if (goldenPionBDT.GetResponse(goldenPion) < 0.1f)
        return false;

    // ...........................................................
    stages.push_back("likelyGoldenPion");

    lastStageIndex = stages.size() - 1;
    lastStage = std::to_string(lastStageIndex) + "_" + stages.back();
    counter.AddEventPassingCut(lastStage, topology, eventIndex);
    // ...........................................................
    
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Cut on the muon-pion opening angle
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    const auto openingAngle = std::acos(muon.r_direction.Dot(goldenPion.r_direction));
    if (openingAngle > 2.6f)
        return false;
    
    // ...........................................................
    stages.push_back("noBackToBack");

    lastStageIndex = stages.size() - 1;
    lastStage = std::to_string(lastStageIndex) + "_" + stages.back();
    counter.AddEventPassingCut(lastStage, topology, eventIndex);
    // ...........................................................
    
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Veto showers
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
    unsigned int nShowers = 0;
    for (const auto &particle : allParticles)
    {
        if (particle.r_trackShower < 0.5f)
            nShowers++;
    }

    if (nShowers > 1)
        return false;
    
    // ...........................................................
    stages.push_back("showerVeto");

    lastStageIndex = stages.size() - 1;
    lastStage = std::to_string(lastStageIndex) + "_" + stages.back();
    counter.AddEventPassingCut(lastStage, topology, eventIndex);
    // ...........................................................


    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Insist on one proton candidate
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    if (protons.size() == 0)
        return false;
    
    // ...........................................................
    stages.push_back("hasProton");

    lastStageIndex = stages.size() - 1;
    lastStage = std::to_string(lastStageIndex) + "_" + stages.back();
    counter.AddEventPassingCut(lastStage, topology, eventIndex);
    // ...........................................................
    
    return true;
}

/**
 *  @brief  The main function called from the command line: root -l makeSelectionTableNew.C
 */
void makeSelectionTable()
{
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Configuration
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    const unsigned int nEvents = std::numeric_limits<unsigned int>::max();
    const unsigned int nSkip = 0;
    const bool shouldRunSelection = true;
    const bool shouldWriteOutputTrees = true;
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    EventManager em("/uboone/data/users/asmith/ubcc1pi/17022020/eventSelection.root");
    
    BDT goldenPionBDT("/uboone/data/users/asmith/ubcc1pi/24022020/goldenPionBDT_weights.xml");
    BDT protonBDT("/uboone/data/users/asmith/ubcc1pi/24022020/protonBDT_weights.xml");
    BDT muonBDT("/uboone/data/users/asmith/ubcc1pi/24022020/muonBDT_weights.xml");

    SelectionCounter counter;

    // Event loop
    const auto maxEventToProcess = std::min(nSkip + nEvents, em.GetNEvents());
    const auto firstEventIndex = std::min(nSkip, maxEventToProcess);
    for (unsigned int eventIndex = firstEventIndex; eventIndex < maxEventToProcess; ++eventIndex)
    {
        em.LoadEvent(eventIndex);

        if (eventIndex % 10000 == 0)
            std::cout << eventIndex << " / " << (maxEventToProcess - firstEventIndex) << std::endl;

        // Run the event selection
        std::string lastStage;
        auto lastStageIndex = std::numeric_limits<unsigned int>::max();
        auto goldenPionId = std::numeric_limits<unsigned int>::max();
        auto muonId = std::numeric_limits<unsigned int>::max();
        std::vector<unsigned int> protonIds;

        bool isSelected = false;
        if (shouldRunSelection)
            isSelected = SelectEvent(em, eventIndex, goldenPionBDT, protonBDT, muonBDT, counter, lastStage, lastStageIndex, goldenPionId, muonId, protonIds);

        // Write the flat particle tree for use when training the particle BDTs
        if (shouldWriteOutputTrees)
        {
            em.FillParticleTree(goldenPionBDT, protonBDT, muonBDT, isSelected, lastStage, lastStageIndex, goldenPionId, muonId, protonIds);
            em.FillEventTree(goldenPionBDT, protonBDT, muonBDT, isSelected, lastStage, lastStageIndex, goldenPionId, muonId, protonIds);
        }
    }

    counter.PrintBreakdown();
    
}
