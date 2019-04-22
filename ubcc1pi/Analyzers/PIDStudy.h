/**
 *  @file  ubcc1pi/Analyzers/PIDStudy.h
 *
 *  @brief The header file for the PID study analyzer.
 */

#ifndef UBCC1PI_ANALYZERS_PID_STUDY
#define UBCC1PI_ANALYZERS_PID_STUDY

#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDAnalyzer.h"

#include "art/Framework/Services/Optional/TFileService.h"
#include <TTree.h>

#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Table.h"
#include "canvas/Utilities/InputTag.h"

namespace ubcc1pi
{

/**
 *  @brief  The PID study class
 */
class PIDStudy : public art::EDAnalyzer
{
    public:
        /**
         *  @brief  The configuration structure
         */
        struct Config
        {
            fhicl::Atom<art::InputTag> MCTruthLabel
            {
                fhicl::Name("MCTruthLabel"),
                fhicl::Comment("The label for the neutrino MCTruth producer")
            };
            
            fhicl::Atom<art::InputTag> MCParticleLabel
            {
                fhicl::Name("MCParticleLabel"),
                fhicl::Comment("The label for the MCParticle producer")
            };
            
            fhicl::Atom<art::InputTag> BacktrackerLabel
            {
                fhicl::Name("BacktrackerLabel"),
                fhicl::Comment("The label for the MCParticle to Hit backtracker producer")
            };
            
            fhicl::Atom<art::InputTag> PFParticleLabel
            {
                fhicl::Name("PFParticleLabel"),
                fhicl::Comment("The label for the PFParticle producer (Pandora)")
            };
            
            fhicl::Atom<art::InputTag> TrackLabel
            {
                fhicl::Name("TrackLabel"),
                fhicl::Comment("The label for the Track producer")
            };
            
            fhicl::Atom<art::InputTag> PIDLabel
            {
                fhicl::Name("PIDLabel"),
                fhicl::Comment("The label for the PID producer")
            };
            
            fhicl::Atom<art::InputTag> SliceLabel
            {
                fhicl::Name("SliceLabel"),
                fhicl::Comment("The label for the Slice producer (Pandora)")
            };
            
            fhicl::Atom<art::InputTag> HitLabel
            {
                fhicl::Name("HitLabel"),
                fhicl::Comment("The label for the Hit producer")
            };
        };

        /**
         *  @brief  The output PID algorithm level structure
         */
        struct OutputAlgorithm
        {
            // Event metadata
            int          m_run;                   ///< The run number
            int          m_subRun;                ///< The subrun number
            int          m_event;                 ///< The event number
            bool         m_isSignal;              ///< If the event is a true CC1Pi signal
        
            // PFParticle details
            int          m_nPFPHits;              ///< The number of hits associated to the PFParticle

            // Matched True Particle details
            int          m_truePdgCode;           ///< The particle PDG code
            float        m_trueMomentum;          ///< The particle momentum
            float        m_trueMatchPurity;       ///< Match purity to the true particle
            float        m_trueMatchCompleteness; ///< Match completeness to the true particle

            // PID algorithm details
            std::string  m_name;                  ///< The algorithm name
            int          m_variableType;          ///< The variable type from the PID object
            int          m_trackDir;              ///< The track direction (forward, backward, no direction)
            int          m_nDOF;                  ///< The number of degrees of freedom
            int          m_assumedPdg;            ///< The PDG assumed by the algorithm
            bool         m_planeWUsed;            ///< If the W plane is used
            bool         m_planeUUsed;            ///< If the U plane is used
            bool         m_planeVUsed;            ///< If the V plane is used

            // The PID algorithm output
            float        m_value;                 ///< The output value of the PID algorithm
        };

        /**
         *  @brief  Constructor
         *
         *  @param  config the set of input fhicl parameters
         */
        PIDStudy(const art::EDAnalyzer::Table<Config> &config);

        /**
         *  @brief  Analyze an event
         *
         *  @param  event the event to analyze
         */
        void analyze(const art::Event &event);

    private:
            
        art::EDAnalyzer::Table<Config>  m_config;          ///< The FHiCL configuration options
        
        TTree                          *m_pAlgorithmTree;  ///< The output tree for all PID algorithm level data
        OutputAlgorithm                 m_outputAlgorithm; ///< The output algorithm-level object
};

} // namespace ubcc1pi

DEFINE_ART_MODULE(ubcc1pi::PIDStudy)

#endif
