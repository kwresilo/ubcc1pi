/**
 *  @file  ubcc1pi/Analyzers/PatrecBenchmarkStudy.h
 *
 *  @brief The header file for the patrec benchmark study analyzer.
 */

#ifndef UBCC1PI_ANALYZERS_PATREC_BENCHMARK_STUDY
#define UBCC1PI_ANALYZERS_PATREC_BENCHMARK_STUDY

#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDAnalyzer.h"

#include "art/Framework/Services/Optional/TFileService.h"
#include <TTree.h>

#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Table.h"
#include "canvas/Utilities/InputTag.h"

#include "ubcc1pi/Helpers/TruthHelper.h"
#include "ubcc1pi/Helpers/CollectionHelper.h"

namespace ubcc1pi
{

/**
 *  @brief  The patrec benchmark study class
 */
class PatrecBenchmarkStudy : public art::EDAnalyzer
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
                fhicl::Comment("The label for the MCParticle to hit backtracker producer")
            };
            
            fhicl::Atom<art::InputTag> HitLabel
            {
                fhicl::Name("HitLabel"),
                fhicl::Comment("The label for the hit producer")
            };
            
            fhicl::Atom<art::InputTag> SliceLabel
            {
                fhicl::Name("SliceLabel"),
                fhicl::Comment("The label for the slice producer")
            };
            
            fhicl::Atom<art::InputTag> PFParticleLabel
            {
                fhicl::Name("PFParticleLabel"),
                fhicl::Comment("The label for the PFParticle producer")
            };
            
            fhicl::OptionalAtom<art::InputTag> AlternatePFParticleLabel
            {
                fhicl::Name("AlternatePFParticleLabel"),
                fhicl::Comment("The label for the new PFParticle producer for a re-interpreted hierarchy")
            };
        };
        
        /**
         *  @brief  The output event level structure
         */
        struct OutputEvent
        {
            // Event metadata
            int                m_run;             ///< The run number
            int                m_subRun;          ///< The subrun number
            int                m_event;           ///< The event number

            // Kinematics
            float              m_nuEnergy;        ///< The neutrino energy

            // Slicing details
            bool               m_hasSlice;        ///< If the event has a slice
            int                m_nNuHitsTotal;    ///< The total number of neutrino induced hits in the whole event
            bool               m_hasNuHits;       ///< If there are any neutrino induced hits in the whole event
                           
            float              m_mcsPurity;       ///< The purity of the most complete slice (mcs) in the event
            float              m_mcsCompleteness; ///< The completeness of the most complete slice (mcs) in the event
            bool               m_isMCSSelected;   ///< If the most complete slice is selected
            
            // Particle details
            std::vector<int>   m_mcp_pdg;                   ///< The pdg codes of the MCParticles
            std::vector<float> m_mcp_momentum;              ///< The momentum of the MCParticles
            std::vector<float> m_mcp_hitWeight;             ///< The weighted number of hits from each of the MCParticles
            std::vector<float> m_mcp_fracHitsInMCS;         ///< The fraction of hits from the MCParticle that are in the most complete slice
            std::vector<int>   m_mcp_nPFPMatches;           ///< The number of primary PFParticles that match to this MCParticle
            std::vector<bool>  m_mcp_hasMatchedPFP;         ///< If the MCParticle has any matched PFParticles
            std::vector<float> m_mcp_bestMatchPurity;       ///< The purity of the best match
            std::vector<float> m_mcp_bestMatchCompleteness; ///< The completeness of the best match
        };

        /**
         *  @brief  Constructor
         *
         *  @param  config the set of input fhicl parameters
         */
        PatrecBenchmarkStudy(const art::EDAnalyzer::Table<Config> &config);

        /**
         *  @brief  Analyze an event
         *
         *  @param  event the event to analyze
         */
        void analyze(const art::Event &event);

    private:

        /**
         *  @brief  Reset the output event struture, and set the run, subrun and event numbers for the input event
         *
         *  @param  event the art event
         */
        void ResetOutputEvent(const art::Event &event);

        art::EDAnalyzer::Table<Config>  m_config;      ///< The FHiCL configuration options

        TTree                          *m_pEventTree;  ///< The output tree
        OutputEvent                     m_outputEvent; ///< The output event structure
};

} // namespace ubcc1pi

DEFINE_ART_MODULE(ubcc1pi::PatrecBenchmarkStudy)

#endif
