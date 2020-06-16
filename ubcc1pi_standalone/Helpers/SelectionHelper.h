/**
 *  @file  ubcc1pi_standalone/Helpers/SelectionHelper.h
 *
 *  @brief The header file for the selection helper class
 */

#ifndef UBCC1PI_STANDALONE_HELPERS_SELECTION_HELPER
#define UBCC1PI_STANDALONE_HELPERS_SELECTION_HELPER

#include "ubcc1pi_standalone/Interface/Event.h"

#include "ubcc1pi_standalone/Helpers/AnalysisHelper.h"
#include "ubcc1pi_standalone/Helpers/BDTHelper.h"

#include <memory>
#include <string>
#include <functional>

namespace ubcc1pi
{

/**
 *  @brief  The selection helper class
 */
class SelectionHelper
{
    public:
        /**
         *  @brief  The event selection class
         *          This class isn't very well designed and is rather "stateful" which isn't nice, work is needed to make it more friendly
         */
        class EventSelection
        {
            public:
                class CutManager
                {
                    public:
                        /**
                         *  @brief  Constructor
                         */
                        CutManager();

                        /**
                         *  @brief  Get the result of a cut using the supplied method and possibly store the result if opimizing
                         *
                         *  @param  name the name of the cut
                         *  @param  method the cut method supplied as a lambda function returun true if the cut passes
                         *
                         *  @return boolean, the result of the cut - true if passes
                         */
                        bool GetCutResult(const std::string &name, const std::function<bool()> &method);
                        
                        /**
                         *  @brief  Get the result of a cut using the supplied method and possibly store the result if opimizing
                         *
                         *  @param  name the name of the cut
                         *  @param  method the cut method supplied as a lambda function returun true if the cut passes
                         *
                         *  @return boolean, the result of the cut - true if passes
                         */
                        bool GetCutResult(const std::string &name, const std::function<bool(const float &)> &method);

                        /**
                         *  @brief  Set the determined PDG code of a given reco particle
                         *
                         *  @param  recoParticleIndex the index of the reco particle
                         *  @param  pdgCode the PDG code to set
                         */
                        void SetParticlePdg(const unsigned int recoParticleIndex, const int pdgCode);

                        /**
                         *  @brief  The cut structure
                         */
                        struct Cut
                        {
                            public:
                                std::string m_name;            ///< The name of the cut
                                bool        m_canDisable;      ///< If the cut can be switched off

                                bool        m_hasValue;        ///< If the cut has a floating point value
                                float       m_nominal;         ///< The nominal value of the cut
                                float       m_min;             ///< The minimum cut value
                                float       m_max;             ///< The maximum cut value

                                bool        m_shouldOptimize;  ///< If we should optimize the parameters of the cut
                                std::string m_searchQuery;     ///< Insist that the signal classification also contains this query when optimizing
                                bool        m_isEnabled;       ///< If the cut is enabled
                                float       m_value;           ///< The value if the cut
                        };
                    
                    private:

                        friend EventSelection;

                        /**
                         *  @brief  Determine if the input cut name exists
                         *
                         *  @param  name the cut name
                         *
                         *  @return boolean, true if cut exists
                         */
                        bool HasCut(const std::string &name) const;

                        /**
                         *  @brief  Get the cut with the specified name
                         *
                         *  @param  name the name of the cut
                         *
                         *  @return reference to the cut
                         */
                        Cut& GetCut(const std::string &name);
                        
                        /**
                         *  @brief  Get the result of the cut using the supplied method
                         *
                         *  @param  cut the cut name
                         *  @param  method the method
                         *
                         *  @return booean, true if cut passes
                         */
                        bool GetCutResult(const Cut &cut, const std::function<bool()> &method) const;
    
                        /**
                         *  @brief  Get the result of the cut using the supplied method
                         *
                         *  @param  cut the cut name
                         *  @param  method the method
                         *
                         *  @return booean, true if cut passes
                         */
                        bool GetCutResult(const Cut &cut, const std::function<bool(const float&)> &method) const;
    
                        std::vector<Cut>                          m_cuts;                  ///< The cuts declared by the user
                        unsigned int                              m_nScanPoints;           ///< The number of points to scan while optimizing
                
                        std::shared_ptr<Event>                    m_pEvent;                ///< The current event
                        AnalysisHelper::SampleType                m_sampleType;            ///< The sample type of the current event
                        float                                     m_weight;                ///< The weight of the current event


                        AnalysisHelper::EventCounter              m_defaultEventCounter;   ///< The event counter when not optimizing

                        std::string                               m_cutOptimizing;         ///< The name of the cut we are currently optimizing
                        AnalysisHelper::EventCounter              m_disabledEventCounter;  ///< The event counter with the current cut disabled
                        std::vector<AnalysisHelper::EventCounter> m_enabledEventCounters;  ///< The event counters with the current cut enabled for various cut values
                        std::vector<float>                        m_values;                ///< The cut values samples
                        std::vector<int>                          m_assignedPdgCodes;      ///< The assigned PDG codes of the reco particles
                };

                /**
                 *  @brief  The BDT manager class
                 */
                class BDTManager
                {
                    public:
                        /**
                         *  @brief  Add a BDT with the given name and features
                         *
                         *  @param  name the BDT name
                         *  @param  featureNames the features
                         */
                        void Add(const std::string &name, const std::vector<std::string> &featureNames);
    
                        /**
                         *  @brief  Get the BDT with supplied name
                         *
                         *  @param  name the name
                         *
                         *  @return the BDT
                         */
                        BDTHelper::BDT & Get(const std::string &name);

                    private:

                        std::vector< std::shared_ptr<BDTHelper::BDT> > m_bdts; ///< The BDTs
                };

                /**
                 *  @brief  An event selection method
                 */
                typedef std::function<bool(const std::shared_ptr<Event> &, BDTManager &, CutManager &cuts)> SelectionMethod;

                /**
                 *  @brief  Declare a cut which can be optionally be disabled during optimisation
                 *
                 *  @param  name the name of the cut
                 */
                void DeclareCut(const std::string &name);
                
                /**
                 *  @brief  Declare a cut which can be optionally be disabled, and have the cut value varied during optimisation
                 *
                 *  @param  name the name of the cut
                 *  @param  nominal the nominal cut value used when we aren't optimizing the parameter
                 */
                void DeclareCut(const std::string &name, const float &nominal);
               
                /**
                 *  @brief  Enable the optimization of a cut
                 *
                 *  @param  name
                 *  @param  searchQuery additional requirement on classification to use as signal while optimizing
                 */
                void EnableOptimization(const std::string &name, const std::string &searchQuery = "");

                /**
                 *  @brief  Enable the optimization of a cut
                 *
                 *  @param  name
                 *  @param  canDisable if the cut can be disabled
                 *  @param  min the minimum cut value to test
                 *  @param  max the maximum cut value to test
                 *  @param  searchQuery additional requirement on classification to use as signal while optimizing
                 */
                void EnableOptimization(const std::string &name, const bool canDisable, const float &min, const float &max, const std::string &searchQuery = "");

                /**
                 *  @brief  Assign a BDT to this selection
                 *
                 *  @param  bdtName the name of the BDT
                 *  @param  featureNames the features
                 */
                void AssignBDT(const std::string &bdtName, const std::vector<std::string> &featureNames);

                /**
                 *  @brief  Get the cut names
                 *
                 *  @return the cuts
                 */
                std::vector<std::string> GetCuts() const;

                /**
                 *  @brief  Define the event selection method
                 *
                 *  @param  method the actual event selection code, supplied as a lambda function returning true if the selection passes
                 */
                void DefineSelectionMethod(const SelectionMethod &method);

                /**
                 *  @brief  Run the event selection method and optimize the free parameters storing the result
                 *
                 *  @param  dataBNBFileName the BNB data file name
                 *  @param  overlayFileName the overlay file name
                 *  @param  overlayWeight the overlay weight
                 *  @param  dataEXTFileName the EXT data file name
                 *  @param  dataEXTWeight the EXT data weight
                 *  @param  dirtFileName the dirt file name
                 *  @param  dirtWeight the dirt weight
                 *  @param  nScanPoints the number of points to scan when optimizing cut values
                 *  @param  processFraction the fraction of events to use for optimization
                 */
                void Optimize(const std::string &dataBNBFileName, const std::string &overlayFileName,
                        const float overlayWeight, const std::string &dataEXTFileName, const float dataEXTWeight,
                        const std::string &dirtFileName, const float dirtWeight, const unsigned int nScanPoints = 20u, const float processFraction = 0.2f);

                /**
                 *  @brief  Run the event selection method with the stored parameters on a single event
                 *
                 *  @param  pEvent the event to run on
                 *  @param  cutsPassed the output vector of cuts that were passes
                 *  @param  assignedPdgCodes the output vector of PDG codes assigned to each reco particle in the event
                 *
                 *  @return boolean, true if all cuts were passed
                 */
                bool Execute(const std::shared_ptr<Event> &pEvent, std::vector<std::string> &cutsPassed, std::vector<int> &assignedPdgCodes);

                /**
                 *  @brief  Run the event selection method with the stored parameters (nominal is used if not optimized) and print the performance
                 *
                 *  @param  dataBNBFileName the BNB data file name
                 *  @param  overlayFileName the overlay file name
                 *  @param  overlayWeight the overlay weight
                 *  @param  dataEXTFileName the EXT data file name
                 *  @param  dataEXTWeight the EXT data weight
                 *  @param  dirtFileName the dirt file name
                 *  @param  dirtWeight the dirt weight
                 *  @param  shouldPrint if we should print the result
                 *  @param  processFraction the fraction of events to use for optimization
                 *  @param  nEntriesToPrint the number of entries to print in the output table
                 */
                void Execute(const std::string &dataBNBFileName, const std::string &overlayFileName, const float overlayWeight,
                        const std::string &dataEXTFileName, const float dataEXTWeight, const std::string &dirtFileName,
                        const float dirtWeight, const bool shouldPrint = true, const float processFraction = 1.f,
                        const unsigned int nEntriesToPrint = 20u);

            private:

                CutManager       m_cutManager;      ///< The cut manager
                SelectionMethod  m_selectionMethod; ///< The event selection method
                BDTManager       m_bdtManager;      ///< The BDTs manager
        };

        /**
         *  @brief  Get the CC inclusive selection
         *
         *  @return the event selection
         */
        static EventSelection GetCCInclusiveSelection();

        /**
         *  @brief  Get the default event selection
         *
         *  @return the event selection
         */
        static EventSelection GetDefaultSelection();
};

} // namespace ubcc1pi

#endif
