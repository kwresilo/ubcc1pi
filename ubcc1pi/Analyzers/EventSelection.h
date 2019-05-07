/**
 *  @file  ubcc1pi/Analyzers/EventSelection.h
 *
 *  @brief The header file for the event selection analyzer.
 */

#ifndef UBCC1PI_ANALYZERS_EVENT_SELECTION_STUDY
#define UBCC1PI_ANALYZERS_EVENT_SELECTION_STUDY

#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDAnalyzer.h"

#include "art/Framework/Services/Optional/TFileService.h"
#include <TTree.h>

#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Table.h"
#include "canvas/Utilities/InputTag.h"

#include "ubcc1pi/Objects/EventSelector.h"

namespace ubcc1pi
{

/**
 *  @brief  The PID study class
 */
class EventSelection : public art::EDAnalyzer
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
            
            fhicl::Atom<art::InputTag> ShowerLabel
            {
                fhicl::Name("ShowerLabel"),
                fhicl::Comment("The label for the Shower producer")
            };
            
            fhicl::Atom<art::InputTag> PIDLabel
            {
                fhicl::Name("PIDLabel"),
                fhicl::Comment("The label for the PID producer")
            };
        };
        
        /**
         *  @brief  The output particle level structure
         */
        struct OutputParticle
        {
            // Event metadata
            int          m_run;                   ///< The run number
            int          m_subRun;                ///< The subrun number
            int          m_event;                 ///< The event number
            bool         m_isSignal;              ///< If the event is a true CC1Pi signal
        
            // Matched True Particle details
            bool         m_hasMatchedMCParticle;  ///< If the particle has a matched MCParticle
            int          m_truePdgCode;           ///< The particle PDG code
            float        m_trueMomentum;          ///< The particle momentum
            float        m_trueMatchPurity;       ///< Match purity to the true particle
            float        m_trueMatchCompleteness; ///< Match completeness to the true particle

            // The PID algorithm outputs
            float        m_trackShowerScore;      ///< The track vs. shower score
            float        m_protonMIPScore;        ///< The proton vs. MIP score
            float        m_muonPionScore;         ///< The muon vs. pion score
        };

        /**
         *  @brief  Constructor
         *
         *  @param  config the set of input fhicl parameters
         */
        EventSelection(const art::EDAnalyzer::Table<Config> &config);

        /**
         *  @brief  Analyze an event
         *
         *  @param  event the event to analyze
         */
        void analyze(const art::Event &event);

    private:
        
        art::EDAnalyzer::Table<Config>  m_config;          ///< The FHiCL configuration options
        
        TTree                          *m_pParticleTree;   ///< The output tree for all particle level data
        OutputParticle                  m_outputParticle;  ///< The output particle-level object
};

} // namespace ubcc1pi

DEFINE_ART_MODULE(ubcc1pi::EventSelection)

#endif
