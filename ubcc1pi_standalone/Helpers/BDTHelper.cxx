/**
 *  @file  ubcc1pi_standalone/Helpers/BDTHelper.cxx
 *
 *  @brief Implementation of the BDT helper class
 */

#include "ubcc1pi_standalone/Helpers/BDTHelper.h"
#include "ubcc1pi_standalone/Helpers/AnalysisHelper.h"

#include <TString.h>
#include <random>
#include <algorithm>

namespace ubcc1pi
{

BDTHelper::EventShuffler::EventShuffler(const unsigned int nEvents, const unsigned int nTraining)
{
    if (nTraining > nEvents)
        throw std::invalid_argument("BDTHelper::EventShuffler::EventShuffler - Number of training events greater than total number of events");
    
    // Setup the training and testing events
    const auto nTesting = nEvents - nTraining;
    m_isTrainingEvent.insert(m_isTrainingEvent.end(), nTraining, true);
    m_isTrainingEvent.insert(m_isTrainingEvent.end(), nTesting, false);

    // Randomize
    std::shuffle(m_isTrainingEvent.begin(), m_isTrainingEvent.end(), std::default_random_engine());
}
// -----------------------------------------------------------------------------------------------------------------------------------------

bool BDTHelper::EventShuffler::IsTrainingEvent(const unsigned int eventIndex)
{
    if (eventIndex >= m_isTrainingEvent.size())
        throw std::invalid_argument("BDTHelper::EventShuffler::IsTrainingEvent - Input eventIndex out of range");

    return m_isTrainingEvent.at(eventIndex);
}


// -----------------------------------------------------------------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------------------------------------------------------------
                
BDTHelper::BDTBase::BDTBase(const std::string &bdtName, const std::vector<std::string> &featureNames) :
    m_bdtName(bdtName),
    m_featureNames(featureNames)
{
}

// -----------------------------------------------------------------------------------------------------------------------------------------

std::string BDTHelper::BDTBase::GetName() const
{
    return m_bdtName;
}

// -----------------------------------------------------------------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------------------------------------------------------------
                
BDTHelper::BDTFactory::BDTFactory(const std::string &bdtName, const std::vector<std::string> &featureNames) :
    BDTBase(bdtName, featureNames),
    m_isBooked(false),
    m_nSignalTraining(0u),
    m_nSignalTest(0u),
    m_nBackgroundTraining(0u),
    m_nBackgroundTest(0)
{
    // Setup the output file and make the TMVA factory
    m_pOutputFile = TFile::Open((bdtName + "_BDT.root").c_str(), "RECREATE" );
    m_pFactory = new TMVA::Factory(TString(bdtName), m_pOutputFile, "!V:!Silent:Color:DrawProgressBar:AnalysisType=Classification");
    m_pDataLoader = new TMVA::DataLoader((bdtName + "_dataset").c_str());

    // Add the features to the data loader
    for (const auto &feature : featureNames)
        m_pDataLoader->AddVariable(TString(feature));
}

// -----------------------------------------------------------------------------------------------------------------------------------------

BDTHelper::BDTFactory::~BDTFactory()
{
    m_pOutputFile->Close();
}

// -----------------------------------------------------------------------------------------------------------------------------------------
                
void BDTHelper::BDTFactory::AddEntry(const std::vector<float> &features, const bool isSignal, const bool isTraining, const float weight)
{
    if (isSignal && isTraining)
        this->AddSignalTrainingEntry(features, weight);
    
    if (!isSignal && isTraining)
        this->AddBackgroundTrainingEntry(features, weight);
    
    if (isSignal && !isTraining)
        this->AddSignalTestEntry(features, weight);
    
    if (!isSignal && !isTraining)
        this->AddBackgroundTestEntry(features, weight);
}

// -----------------------------------------------------------------------------------------------------------------------------------------
                
void BDTHelper::BDTFactory::AddSignalTrainingEntry(const std::vector<float> &features, const float weight)
{
    if (features.size() != m_featureNames.size())
        throw std::invalid_argument("BDTHelper::BDTFactory::AddSignalTrainingEntry - Wrong number of features");

    m_pDataLoader->AddSignalTrainingEvent(std::vector<double>(features.begin(), features.end()), weight);
    m_nSignalTraining++;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void BDTHelper::BDTFactory::AddBackgroundTrainingEntry(const std::vector<float> &features, const float weight)
{
    if (features.size() != m_featureNames.size())
        throw std::invalid_argument("BDTHelper::BDTFactory::AddBackgroundTrainingEntry - Wrong number of features");

    m_pDataLoader->AddBackgroundTrainingEvent(std::vector<double>(features.begin(), features.end()), weight);
    m_nBackgroundTraining++;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void BDTHelper::BDTFactory::AddSignalTestEntry(const std::vector<float> &features, const float weight)
{
    if (features.size() != m_featureNames.size())
        throw std::invalid_argument("BDTHelper::BDTFactory::AddSignalTestEntry - Wrong number of features");

    m_pDataLoader->AddSignalTestEvent(std::vector<double>(features.begin(), features.end()), weight);
    m_nSignalTest++;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void BDTHelper::BDTFactory::AddBackgroundTestEntry(const std::vector<float> &features, const float weight)
{
    if (features.size() != m_featureNames.size())
        throw std::invalid_argument("BDTHelper::BDTFactory::AddBackgroundTestEntry - Wrong number of features");

    m_pDataLoader->AddBackgroundTestEvent(std::vector<double>(features.begin(), features.end()), weight);
    m_nBackgroundTest++;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void BDTHelper::BDTFactory::BookMethod()
{
    if (!m_isBooked)
    {
        std::cout << "DEBUG: Preparing" << std::endl;
        m_pDataLoader->PrepareTrainingAndTestTree("0==0", m_nSignalTraining, m_nBackgroundTraining, m_nSignalTest, m_nBackgroundTest);
        std::cout << "DEBUG: Booking method" << std::endl;
        m_pFactory->BookMethod(m_pDataLoader, TMVA::Types::kBDT, "BDT", "!H:!V:NTrees=505:MinNodeSize=1.26436%:MaxDepth=4:BoostType=AdaBoost:AdaBoostBeta=0.2:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=20");
    
        m_isBooked = true;
        std::cout << "DEBUG: Method booked" << std::endl;
    }
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void BDTHelper::BDTFactory::OptimizeParameters()
{
    this->BookMethod();
    m_pFactory->OptimizeAllMethodsForClassification();
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void BDTHelper::BDTFactory::TrainAndTest()
{
    this->BookMethod();
    
    std::cout << "DEBUG: Training" << std::endl;
    m_pFactory->TrainAllMethods();
    std::cout << "DEBUG: Testing" << std::endl;
    m_pFactory->TestAllMethods();
    std::cout << "DEBUG: Evaluating" << std::endl;
    m_pFactory->EvaluateAllMethods();
}

// -----------------------------------------------------------------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------------------------------------------------------------

BDTHelper::BDT::BDT(const std::string &bdtName, const std::vector<std::string> &featureNames) :
    BDTBase(bdtName, featureNames),
    m_pReader(new TMVA::Reader("Silent")),
    m_weightsFile(bdtName + "_dataset/weights/" + bdtName + "_BDT.weights.xml"),
    m_features(featureNames.size(), -std::numeric_limits<float>::max())
{
    // Bind the feature vector to the reader
    for (unsigned int iFeature = 0; iFeature < featureNames.size(); ++iFeature)
        m_pReader->AddVariable(featureNames.at(iFeature).c_str(), &m_features.at(iFeature));

    // Book the MVA using the weights file
    m_pReader->BookMVA("BDT", m_weightsFile.c_str());
}

// -----------------------------------------------------------------------------------------------------------------------------------------
                
float BDTHelper::BDT::GetResponse(const std::vector<float> &features)
{
    if (features.size() != m_features.size())
        throw std::invalid_argument("BDTHelper::BDT::GetResponse - Wrong number of features, expected " + std::to_string(m_features.size()) + " but passed " + std::to_string(features.size()));

    // Shallow copy the features
    for (unsigned int iFeature = 0; iFeature < m_features.size(); ++iFeature)
        m_features.at(iFeature) = features.at(iFeature);

    // Run the BDT
    return m_pReader->EvaluateMVA("BDT");
}

// -----------------------------------------------------------------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------------------------------------------------------------
        
bool BDTHelper::GetBDTFeatures(const Event::Reco::Particle &recoParticle, const std::vector<std::string> &featureNames, std::vector<float> &features, const bool shouldDebug)
{
    if (!features.empty())
        throw std::logic_error("BDTHelper::GetBDTFeatures - Input feature vector isn't empty");

    for (const auto &name : featureNames)
    {
        if (shouldDebug)
            std::cout << "DEBUG - Calculating: " << name << std::endl;

        if (name == "logBragg_pToMIP")
        {
            float feature = -std::numeric_limits<float>::max();
            if (!AnalysisHelper::GetLikelihoodRatio(recoParticle.likelihoodForwardProtonW, recoParticle.likelihoodMIPW, feature))
            {
                if (shouldDebug)
                    std::cout << "DEBUG - Can't calculate: " << name << std::endl;
                
                return false;
            }

            features.push_back(feature);
            continue;
        }

        if (name == "logBragg_piToMIP")
        {
            float feature = -std::numeric_limits<float>::max();
            if (!AnalysisHelper::GetLikelihoodRatio(recoParticle.likelihoodForwardPionW, recoParticle.likelihoodMIPW, feature))
            {
                if (shouldDebug)
                    std::cout << "DEBUG - Can't calculate: " << name << std::endl;
                
                return false;
            }

            features.push_back(feature);
            continue;
        }

        if (name == "truncMeandEdx")
        {
            if (!recoParticle.truncatedMeandEdxW.IsSet())
            {
                if (shouldDebug)
                    std::cout << "DEBUG - Can't calculate: " << name << std::endl;
                
                return false;
            }

            features.push_back(recoParticle.truncatedMeandEdxW());
            continue;
        }

        if (name == "protonForward")
        {
            float feature = -std::numeric_limits<float>::max();
            if (!AnalysisHelper::GetSoftmax(recoParticle.likelihoodForwardPionW, recoParticle.likelihoodBackwardProtonW, feature))
            {
                if (shouldDebug)
                    std::cout << "DEBUG - Can't calculate: " << name << std::endl;
                
                return false;
            }

            features.push_back(feature);
            continue;
        }

        if (name == "muonForward")
        {
            float feature = -std::numeric_limits<float>::max();
            if (!AnalysisHelper::GetSoftmax(recoParticle.likelihoodForwardMuonW, recoParticle.likelihoodBackwardMuonW, feature))
            {
                if (shouldDebug)
                    std::cout << "DEBUG - Can't calculate: " << name << std::endl;
                
                return false;
            }

            features.push_back(feature);
            continue;
        }

        if (name == "nDescendents")
        {
            if (!recoParticle.nDescendents.IsSet())
            {
                if (shouldDebug)
                    std::cout << "DEBUG - Can't calculate: " << name << std::endl;
                
                return false;
            }

            features.push_back(recoParticle.nDescendents());
            continue;
        }

        if (name == "nSpacePointsNearEnd")
        {
            if (!recoParticle.nSpacePointsNearEnd.IsSet())
            {
                if (shouldDebug)
                    std::cout << "DEBUG - Can't calculate: " << name << std::endl;
                
                return false;
            }

            features.push_back(recoParticle.nSpacePointsNearEnd());
            continue;
        }

        if (name == "wiggliness")
        {
            if (!recoParticle.wiggliness.IsSet())
            {
                if (shouldDebug)
                    std::cout << "DEBUG - Can't calculate: " << name << std::endl;
                
                return false;
            }

            features.push_back(recoParticle.wiggliness());
            continue;
        }

        if (name == "trackScore")
        {
            if (!recoParticle.trackScore.IsSet())
            {
                if (shouldDebug)
                    std::cout << "DEBUG - Can't calculate: " << name << std::endl;
                
                return false;
            }

            features.push_back(recoParticle.trackScore());
            continue;
        }

        throw std::invalid_argument("BDTHelper::GetBDTFeatures - Unknown feature: \"" + name + "\"");
    }

    // This should never happen - just a sanity check
    if (features.size() != featureNames.size())
        throw std::logic_error("BDTHelper::GetBDTFeatures - Number of features calculated is incorrect!");

    return true;
}


} // namespace ubcc1pi
