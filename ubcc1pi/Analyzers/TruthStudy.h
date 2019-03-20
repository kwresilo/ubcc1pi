/**
 *  @file  ubcc1pi/Analyzers/TruthStudy.h
 *
 *  @brief The header file for the truth study analyzer.
 */

#ifndef UBCC1PI_ANALYZERS_TRUTH_STUDY
#define UBCC1PI_ANALYZERS_TRUTH_STUDY

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
 *  @brief  The truth study class
 */
class TruthStudy : public art::EDAnalyzer
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
        };

        /**
         *  @brief  The interaction output structure
         */
        struct InteractionOutput
        {
            // Event metadata
            int          m_run;          ///< The run number
            int          m_subRun;       ///< The subrun number
            int          m_event;        ///< The event number
            // Interaction information
            std::string  m_interaction;  ///< The truth interaction type
            bool         m_isNuFiducial; ///< If the neutrino is fiducial
            bool         m_isSignal;     ///< If the event is a CC1Pi signal event
            int          m_nMuMinus;     ///< The number of mu-
            int          m_nMuPlus;      ///< The number of mu+
            int          m_nPiPlus;      ///< The number of pi+
            int          m_nPiMinus;     ///< The number of pi-
            int          m_nKPlus;       ///< The number of K+
            int          m_nKMinus;      ///< The number of K-
            int          m_nProton;      ///< The number of p
            int          m_nNeutron;     ///< The number of n
            int          m_nPhoton;      ///< The number of gamma
            int          m_nTotal;       ///< The total number of particles
            // Kinematics
            float        m_nuE;          ///< The neutrino energy
        };
        
        /**
         *  @brief  The signal output structure
         */
        struct SignalOutput
        {
            // Event metadata
            int          m_run;          ///< The run number
            int          m_subRun;       ///< The subrun number
            int          m_event;        ///< The event number
            // Interaction information
            std::string  m_interaction;  ///< The truth interaction type
            int          m_nProton;      ///< The number of p
            // Kinematics
            float        m_nuE;          ///< The neutrino energy
            float        m_muMomX;       ///< The muon momentum in x
            float        m_muMomY;       ///< The muon momentum in y
            float        m_muMomZ;       ///< The muon momentum in z
            float        m_piMomX;       ///< The pion momentum in x
            float        m_piMomY;       ///< The pion momentum in y
            float        m_piMomZ;       ///< The pion momentum in z
        };

        /**
         *  @brief  Constructor
         *
         *  @param  config the set of input fhicl parameters
         */
        TruthStudy(const art::EDAnalyzer::Table<Config> &config);

        /**
         *  @brief  Analyze an event
         *
         *  @param  event the event to analyze
         */
        void analyze(const art::Event &event);

    private:
            
        art::EDAnalyzer::Table<Config>  m_config;            ///< The FHiCL configuration options
        TTree                          *m_pInteractionTree;  ///< The output tree for all interactions
        TTree                          *m_pSignalTree;       ///< The output tree for signal interactions

        InteractionOutput               m_interactionOutput; ///< The output objects for all interactions
        SignalOutput                    m_signalOutput;      ///< The output objects for signal interactions
};

} // namespace ubcc1pi

DEFINE_ART_MODULE(ubcc1pi::TruthStudy)

#endif
