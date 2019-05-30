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

//#include "ubcc1pi/Objects/EventSelector.h"
#include "ubcc1pi/Helpers/CollectionHelper.h"

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
            
            fhicl::Atom<float> Chi2ProtonMIPCut
            {
                fhicl::Name("Chi2ProtonMIPCut"),
                fhicl::Comment("The chi2 (minimum over all views) under the proton hypothesis above which a track is considered a MIP")
            };
        };
        
        /**
         *  @brief  The output event level structure
         */
        struct OutputEvent
        {
            // Event metadata
            int          m_run;                                    ///< The run number
            int          m_subRun;                                 ///< The subrun number
            int          m_event;                                  ///< The event number
                                                                   
            // The truth interaction information                   
            bool         m_isSignal;                               ///< If the event is a true CC1Pi signal
            std::string  m_interaction;                            ///< The interaction type string
            bool         m_isNuFiducial;                           ///< If the neutrino interaction is fiducial
            float        m_nuE;                                    ///< The true neutrino energy
            float        m_nuX;                                    ///< The true neutrino vertex x-position
            float        m_nuY;                                    ///< The true neutrino vertex y-position
            float        m_nuZ;                                    ///< The true neutrino vertex z-position
                                                                   
            // The true particle multiplicities                    
            int          m_nMuMinus;                               ///< The number of mu-
            int          m_nMuPlus;                                ///< The number of mu+
            int          m_nPiPlus;                                ///< The number of pi+
            int          m_nPiMinus;                               ///< The number of pi-
            int          m_nKPlus;                                 ///< The number of K+
            int          m_nKMinus;                                ///< The number of K-
            int          m_nProton;                                ///< The number of p
            int          m_nNeutron;                               ///< The number of n
            int          m_nPhoton;                                ///< The number of gamma
            int          m_nTotal;                                 ///< The total number of particles

            // The event selection info
            float        m_trueMuEnergy;                           ///< The true muon energy
            float        m_trueMuTheta;                            ///< The true muon theta (angle to beam)
            float        m_trueMuPhi;                              ///< The true muon phi (angle about beam)
    
            float        m_truePiEnergy;                           ///< The true pion energy
            float        m_truePiTheta;                            ///< The true pion theta
            float        m_truePiPhi;                              ///< The true pion phi
            float        m_trueMuPiAngle;                          ///< The muon-pion opening angle

            float        m_recoNuX;                                ///< The reco neutrino vertex x-position
            float        m_recoNuY;                                ///< The reco neutrino vertex y-position
            float        m_recoNuZ;                                ///< The reco neutrino vertex z-position

            int          m_nFinalStatePFPs;                        ///< The number of reconstructed final state PFParticles
            bool         m_isRecoNuFiducial;                       ///< If the reconstructed vertex is fiducial
            int          m_nRecoMIPs;                              ///< The number of PFPs passing the MIP selection
            bool         m_isSelected;                             ///< If the event is selected

            bool         m_isSelectedMuonTruthMatched;             ///< If the selected muon PFParticle is matched to any MCParticle
            int          m_selectedMuonTruePdg;                    ///< The true pdg code of the PFParticle selected as the muon
            float        m_selectedMuonCompleteness;               ///< The match completeness for the selected muon
            float        m_selectedMuonPurity;                     ///< The match purity for the selected muon
    
            bool         m_isSelectedPionTruthMatched;             ///< If the selected pion PFParticle is matched to any MCParticle
            int          m_selectedPionTruePdg;                    ///< The true pdg code of the PFParticle selected as the pion
            float        m_selectedPionCompleteness;               ///< The match completeness for the selected pion
            float        m_selectedPionPurity;                     ///< The match purity for the selected pion
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
        
        // TODO Doxygen comments
        void PerformPID(const art::Event &event, const PFParticleVector &finalStates, const float chi2ProtonCut, PFParticleVector &muons, PFParticleVector &pions, PFParticleVector &protons, PFParticleVector &showerLikes) const;
        void SelectTracksAndShowers(const PFParticleVector &finalStates, PFParticleVector &trackLikes, PFParticleVector &showerLikes) const;
        PFParticleVector SelectMIPs(const art::Event &event, const PFParticleVector &finalStates, const float chi2ProtonCut) const;
        bool IsTrackLike(const art::Ptr<recob::PFParticle> &pfParticle) const;
        bool IsShowerLike(const art::Ptr<recob::PFParticle> &pfParticle) const;
        bool PassesMIPSelection(const art::Ptr<recob::PFParticle> &pfParticle, const Association<recob::PFParticle, recob::Track> &pfpToTracks, const Association<recob::Track, anab::ParticleID> &trackToPIDs, const float chi2ProtonCut) const;
        float GetMinChi2Proton(const art::Ptr<anab::ParticleID> &pid) const;
        art::Ptr<recob::PFParticle> SelectMuon(const art::Event &event, const PFParticleVector &mips) const;
        TVector3 GetRecoNeutrinoVertex(const art::Event &event, const PFParticleVector &allPFParticles) const;
        
        art::EDAnalyzer::Table<Config>  m_config;          ///< The FHiCL configuration options
        
        TTree                          *m_pEventTree;      ///< The output tree for all event level data

        OutputEvent                     m_outputEvent;     ///< The output event-level object
};

} // namespace ubcc1pi

DEFINE_ART_MODULE(ubcc1pi::EventSelection)

#endif
