/**
 *  @file  ubcc1pi/Analyzers/RecoStudy.h
 *
 *  @brief The header file for the reco study analyzer.
 */

#ifndef UBCC1PI_ANALYZERS_RECO_STUDY
#define UBCC1PI_ANALYZERS_RECO_STUDY

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
 *  @brief  The reco study class
 */
class RecoStudy : public art::EDAnalyzer
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
         *  @brief  The output event level structure
         */
        struct OutputEvent
        {
            // Event metadata
            int          m_run;            ///< The run number
            int          m_subRun;         ///< The subrun number
            int          m_event;          ///< The event number

            // Kinematics
            float        m_nuEnergy;       ///< The neutrino energy

            // Multiplicities
            int          m_nProtons;       ///< The number of target proton MCParticle
            int          m_nPFPs;          ///< The number of PFParticles

            // Slice details
            int          m_nNuHitsTotal;              ///< The total number of neutrino induced hits in the whole event
            bool         m_hasChosenSlice;            ///< If a slice was identified as the neutrino
            float        m_chosenSlicePurity;         ///< The purity of the chosen neutrino slice
            float        m_chosenSliceCompleteness;   ///< The completeness of the chosen neutrino slice
            bool         m_isChosenSliceMostComplete; ///< If the chosen slice is the most complete in the event

            // Match details
            int          m_nMuMatches;          ///< The number of PFParticles matched to the muon MCParticle
            int          m_nPiMatches;          ///< The number of PFParticles matched to the pion MCParticle
            int          m_nProtonsMatched;     ///< The number of proton MCParticles with at least one PFParticle matched
            int          m_nProtonsMatchedOnce; ///< The number of proton MCParticles with exactly one PFParticle matched
            int          m_nUnmatchedPFPs;      ///< The number of PFParticles that don't share hits with any neutrino primary MCParticle ==> cosmic or induced
        };
        
        /**
         *  @brief  The output particle level structure
         */
        struct OutputParticle
        {
            // Event metadata
            int          m_run;                       ///< The run number
            int          m_subRun;                    ///< The subrun number
            int          m_event;                     ///< The event number
            int          m_nUnmatchedPFPsInEvent;     ///< The number of unmatched neutrino induced PFParticles in the event
            int          m_nPFPsInEvent;              ///< The number of neutrino induced PFParticles in the event
            bool         m_eventHasChosenSlice;       ///< If a slice was identified as the neutrino
            float        m_chosenSlicePurity;         ///< The purity of the chosen neutrino slice
            float        m_chosenSliceCompleteness;   ///< The completeness of the chosen neutrino slice
            bool         m_isChosenSliceMostComplete; ///< If the chosen slice is the most complete in the event
                                            
            // True Particle details             
            int          m_pdgCode;               ///< The particle PDG code
            float        m_momentum;              ///< The particle momentum
            float        m_mcHitWeight;           ///< The total weighted number of hits associated to the MCParticle
            float        m_mcHitWeightFromSlice;  ///< The total weighted number of hits associated to the MCParticle that are in the chosen slice

            // Match details
            int          m_nMatches;              ///< The number of PFParticles matches
            bool         m_isBestMatchTrack;      ///< If the best matched PFParticle has been identified as track-like (false => shower-like)
            float        m_bestMatchPurity;       ///< The purity of the best match (NB. best match has highest purity by definition)
            float        m_bestMatchCompleteness; ///< The completeness of the best match
        };

        /**
         *  @brief  Constructor
         *
         *  @param  config the set of input fhicl parameters
         */
        RecoStudy(const art::EDAnalyzer::Table<Config> &config);

        /**
         *  @brief  Analyze an event
         *
         *  @param  event the event to analyze
         */
        void analyze(const art::Event &event);

    private:
            
        art::EDAnalyzer::Table<Config>  m_config;         ///< The FHiCL configuration options
        
        TTree                          *m_pEventTree;     ///< The output tree for all event level data
        TTree                          *m_pParticleTree;  ///< The output tree for all particle level data
        
        OutputEvent                     m_outputEvent;    ///< The output event-level object
        OutputParticle                  m_outputParticle; ///< The output particle-level object
};

} // namespace ubcc1pi

DEFINE_ART_MODULE(ubcc1pi::RecoStudy)

#endif
