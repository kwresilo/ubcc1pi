/**
 *  @file  ubcc1pi_standalone/Helpers/BDTHelper.h
 *
 *  @brief The header file for the bdt helper class
 */

#ifndef UBCC1PI_STANDALONE_HELPERS_BDT_HELPER
#define UBCC1PI_STANDALONE_HELPERS_BDT_HELPER

#include <string>
#include <vector>
#include <unordered_map>

#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Reader.h"
#include <TFile.h>

#include "ubcc1pi_standalone/Interface/Event.h"

namespace ubcc1pi
{

/**
 *  @brief  The BDT helper class
 */
class BDTHelper
{
    public:
        /**
         *  @brief  Base wrapper class for at root TMVA BDT
         */
        class BDTBase
        {
            public:
                /**
                 *  @brief  Constructor
                 *
                 *  @param  bdtName the name of the BDT
                 *  @param  featureNames the vector of feature variable names
                 */
                BDTBase(const std::string &bdtName, const std::vector<std::string> &featureNames);
                
                /**
                 *  @brief  Get the name
                 *
                 *  @return BDT name
                 */
                std::string GetName() const;
                
            protected:

                std::string               m_bdtName;      ///< The BDT name
                std::vector<std::string>  m_featureNames; ///< The feature variable names
        };

        /**
         *  @brief  Event shuffler class selects a random set of event indices for training and for testing
         */
        class EventShuffler
        {
            public:
                /**
                 *  @brief  Constructor
                 *
                 *  @param  nEvents the total number of events to shuffle
                 *  @param  nTraining the number of events that should be used for training
                 */
                EventShuffler(const unsigned int nEvents, const unsigned int nTraining);

                /**
                 *  @brief  Determine if the input event should be used for training or testing
                 *
                 *  @param  eventIndex the input eventIndex
                 *
                 *  @return boolean, true if event should be used for training
                 */
                bool IsTrainingEvent(const unsigned int eventIndex);

            private:
                std::vector<bool> m_isTrainingEvent; ///< If the event with a given index should be used for training
        };

        /**
         *  @brief  BDT factory class used to train BDTs
         */
        class BDTFactory : public BDTBase
        {
            public:
                /**
                 *  @brief  Constructor
                 *
                 *  @param  bdtName the name of the BDT
                 *  @param  featureNames the vector of feature variable names
                 */
                BDTFactory(const std::string &bdtName, const std::vector<std::string> &featureNames);

                /**
                 *  @brief  Destructor
                 */
                ~BDTFactory();

                /**
                 *  @brief  Add an entry
                 *
                 *  @param  features the features of the entry
                 *  @param  isSignal if the entry is signal or background
                 *  @param  isTraining if the entry is for training or testing
                 *  @param  weight the weight of the entry
                 */
                void AddEntry(const std::vector<float> &features, const bool isSignal, const bool isTraining, const float weight = 1.f);

                /**
                 *  @brief  Add a signal training entry
                 *
                 *  @param  features the features of the entry
                 *  @param  weight the weight of the entry
                 */
                void AddSignalTrainingEntry(const std::vector<float> &features, const float weight = 1.f);
                
                /**
                 *  @brief  Add a background training entry
                 *
                 *  @param  features the features of the entry
                 *  @param  weight the weight of the entry
                 */
                void AddBackgroundTrainingEntry(const std::vector<float> &features, const float weight = 1.f);
                
                /**
                 *  @brief  Add a signal test entry
                 *
                 *  @param  features the features of the entry
                 *  @param  weight the weight of the entry
                 */
                void AddSignalTestEntry(const std::vector<float> &features, const float weight = 1.f);
                
                /**
                 *  @brief  Add a background test entry
                 *
                 *  @param  features the features of the entry
                 *  @param  weight the weight of the entry
                 */
                void AddBackgroundTestEntry(const std::vector<float> &features, const float weight = 1.f);

                /**
                 *  @brief  Optimize the BDT parameters
                 */
                void OptimizeParameters();

                /**
                 *  @brief  Train and test the BDT using the entries added so far
                 */
                void TrainAndTest();

            private:
                
                /**
                 *  @brief  Book the BDT method
                 */
                void BookMethod();

                TFile             *m_pOutputFile; ///< The TMVA output file
                TMVA::Factory     *m_pFactory;    ///< The TMVA factory
                TMVA::DataLoader  *m_pDataLoader; ///< The data loader
                bool               m_isBooked;    ///< If the BDT has been booked

                unsigned int       m_nSignalTraining;     ///< The number of signal training events
                unsigned int       m_nSignalTest;         ///< The number of signal testing events
                unsigned int       m_nBackgroundTraining; ///< The number of background training events
                unsigned int       m_nBackgroundTest;     ///< The number of background testing events
        };

        /**
         *  @brief  The class for using a trained BDT
         */
        class BDT : public BDTBase
        {
            public:
                /**
                 *  @brief  Constructor
                 *
                 *  @param  bdtName the name of the BDT
                 *  @param  featureNames the vector of feature variable names
                 */
                BDT(const std::string &bdtName, const std::vector<std::string> &featureNames);

                /**
                 *  @brief  Get the BDT response for a given set of features
                 *
                 *  @param  features the input vector of features
                 *
                 *  @return the BDT response
                 */
                float GetResponse(const std::vector<float> &features);

            private:
                TMVA::Reader      *m_pReader;     ///< The TMVA reader
                std::string        m_weightsFile; ///< The weights file
                std::vector<float> m_features;    ///< The BDT features
        };

        /**
         *  @brief  The particle BDT feature names 
         */
        static const std::vector<std::string> ParticleBDTFeatureNames;

        /**
         *  @brief  Get the features corresponding to the given names for the input reco particle
         *
         *  @param  recoParticle the input reco particle 
         *  @param  featureNames the feature names to calculate
         *  @param  features the output features
         *
         *  @return boolean, true if all features are available
         */
        static bool GetBDTFeatures(const Event::Reco::Particle &recoParticle, const std::vector<std::string> &featureNames, std::vector<float> &features, const bool shouldDebug = false);
};

// -----------------------------------------------------------------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------------------------------------------------------------

const std::vector<std::string> BDTHelper::ParticleBDTFeatureNames = {
    "logBragg_pToMIP",
    "logBragg_piToMIP",
    "truncMeandEdx",
    "protonForward",
    "muonForward",
    "nDescendents",
    "nSpacePointsNearEnd",
    "wiggliness",
    "trackScore"
};

} // namespace ubcc1pi

#endif
