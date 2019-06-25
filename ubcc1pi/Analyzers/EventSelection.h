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

#include "ubcc1pi/Helpers/CollectionHelper.h"
#include "ubcc1pi/Helpers/BacktrackHelper.h"
#include "ubcc1pi/Helpers/RecoHelper.h"

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
            
            fhicl::Atom<art::InputTag> PIDLabel
            {
                fhicl::Name("PIDLabel"),
                fhicl::Comment("The label for the PID producer")
            };
            
            fhicl::Atom<float> Sin2YZAngleCut
            {
                fhicl::Name("Sin2YZAngleCut"),
                fhicl::Comment("Cut on the squared sine of the angle between a track and a wire in the YZ plane for the PID information to be used on that plane")
            };
            
            fhicl::Atom<float> ContainmentBorder
            {
                fhicl::Name("ContainmentBorder"),
                fhicl::Comment("The maximum distance a position can be from a detector face and be considered uncontained")
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
            float        m_trueNuE;                                ///< The true neutrino energy
            TVector3     m_trueNuVtx;                              ///< The true neutrino vertex position
            bool         m_isTrueNuFiducial;                       ///< If the neutrino interaction is fiducial
                                                                   
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
            int          m_nElectron;                              ///< The number of e
            int          m_nPositron;                              ///< The number of e+
            int          m_nTotal;                                 ///< The total number of particles

            // True particle kinematics
            float        m_trueMuEnergy;                           ///< The true muon energy
            float        m_trueMuTheta;                            ///< The true muon theta (angle to beam)
            float        m_trueMuPhi;                              ///< The true muon phi (angle about beam)
    
            float        m_truePiEnergy;                           ///< The true pion energy
            float        m_truePiTheta;                            ///< The true pion theta
            float        m_truePiPhi;                              ///< The true pion phi

            float        m_trueMuPiAngle;                          ///< The muon-pion opening angle

            // Reconstructed information
            bool                   m_hasRecoNeutrino;              ///< If the event has a reconstructed neutrino
            TVector3               m_recoNuVtx;                    ///< The SCE corrected reco neutrino vertex position
            bool                   m_isRecoNuFiducial;             ///< If the reconstructed vertex is fiducial
            float                  m_topologicalScore;             ///< The pandora neutrino ID topological score for this slice

            // Particle level information
            int                    m_nFinalStatePFPs;              ///< The number of reconstructed final state PFParticles
            
            // Truth matching
            std::vector<bool>      m_hasMatchedMCParticleVect;     ///< If the particles have a matched MCParticle 
            std::vector<int>       m_truePdgCodeVect;              ///< The best matched MCParticles PDG code
            std::vector<float>     m_truthMatchCompletenessVect;   ///< The completeness of the match
            std::vector<float>     m_truthMatchPurityVect;         ///< The purity of the match
            std::vector<float>     m_trueMomentumXVect;            ///< The true momentum of the particle - X
            std::vector<float>     m_trueMomentumYVect;            ///< The true momentum of the particle - Y
            std::vector<float>     m_trueMomentumZVect;            ///< The true momentum of the particle - Z
            std::vector<float>     m_trueStartXVect;                ///< The true start position of the particle - X
            std::vector<float>     m_trueStartYVect;                ///< The true start position of the particle - Y
            std::vector<float>     m_trueStartZVect;                ///< The true start position of the particle - Z
            
            // Pandora info
            std::vector<int>       m_nHitsUVect;                   ///< The number of hits in the U view
            std::vector<int>       m_nHitsVVect;                   ///< The number of hits in the V view
            std::vector<int>       m_nHitsWVect;                   ///< The number of hits in the W view
            std::vector<float>     m_trackShowerVect;              ///< The pandora track vs. shower score

            // Track info
            std::vector<bool>      m_hasTrackInfoVect;             ///< If the PFParticle has an associated track
            std::vector<float>     m_startXVect;                   ///< The SCE corrected reconstructed start position of the particle - X
            std::vector<float>     m_startYVect;                   ///< The SCE corrected reconstructed start position of the particle - Y
            std::vector<float>     m_startZVect;                   ///< The SCE corrected reconstructed start position of the particle - Z
            std::vector<float>     m_endXVect;                     ///< The SCE corrected reconstructed end position of the particle - X
            std::vector<float>     m_endYVect;                     ///< The SCE corrected reconstructed end position of the particle - Y
            std::vector<float>     m_endZVect;                     ///< The SCE corrected reconstructed end position of the particle - Z
            std::vector<float>     m_directionXVect;               ///< The reconstructed direction of the particle - X
            std::vector<float>     m_directionYVect;               ///< The reconstructed direction of the particle - Y
            std::vector<float>     m_directionZVect;               ///< The reconstructed direction of the particle - Z
            std::vector<float>     m_yzAngleVect;                  ///< The reconstructed angle in the YZ plane to the vertical
            std::vector<float>     m_lengthVect;                   ///< The reconstructed length
            std::vector<bool>      m_isContainedVect;              ///< If the particle is contained within the detector

            std::vector<bool>      m_hasPIDInfoVect;             ///< If the PFParticle has an associated PID object

            // Chi2 proton
            std::vector<bool>      m_isWChi2pAvailableVect;        ///< If the chi2 proton is available for the W plane
            std::vector<float>     m_chi2pWVect;                   ///< The chi2 proton for the W plane

            std::vector<bool>      m_isUChi2pAvailableVect;        ///< If the chi2 proton is available for the U plane
            std::vector<float>     m_chi2pUVect;                   ///< The chi2 proton for the U plane
            
            std::vector<bool>      m_isVChi2pAvailableVect;        ///< If the chi2 proton is available for the V plane
            std::vector<float>     m_chi2pVVect;                   ///< The chi2 proton for the V plane

            std::vector<bool>      m_isUVChi2pAvailableVect;       ///< If the weighted average chi2 proton over the U & V planes is available
            std::vector<float>     m_chi2pUVVect;                  ///< The weighted average chi2 proton over the U & V planes

            // Bragg likelihood proton
            std::vector<bool>      m_isWBraggpAvailableVect;       ///< If the Bragg proton is available for the W plane                        
            std::vector<float>     m_braggpWVect;                  ///< The Bragg proton for the W plane
                                                                                                                                               
            std::vector<bool>      m_isUBraggpAvailableVect;       ///< If the Bragg proton is available for the U plane
            std::vector<float>     m_braggpUVect;                  ///< The Bragg proton for the U plane
                                                                                                                                               
            std::vector<bool>      m_isVBraggpAvailableVect;       ///< If the Bragg proton is available for the V plane
            std::vector<float>     m_braggpVVect;                  ///< The Bragg proton for the V plane
                                                                                                                                               
            std::vector<bool>      m_isUVBraggpAvailableVect;      ///< If the weighted average Bragg proton over the U & V planes is available
            std::vector<float>     m_braggpUVVect;                 ///< The weighted average Bragg proton over the U & V planes
            
            // Bragg likelihood MIP
            std::vector<bool>      m_isWBraggMIPAvailableVect;     ///< If the Bragg MIP is available for the W plane                        
            std::vector<float>     m_braggMIPWVect;                ///< The Bragg MIP for the W plane
                                                                                                                                               
            std::vector<bool>      m_isUBraggMIPAvailableVect;     ///< If the Bragg MIP is available for the U plane
            std::vector<float>     m_braggMIPUVect;                ///< The Bragg MIP for the U plane
                                                                                                                                               
            std::vector<bool>      m_isVBraggMIPAvailableVect;     ///< If the Bragg MIP is available for the V plane
            std::vector<float>     m_braggMIPVVect;                ///< The Bragg MIP for the V plane
                                                                                                                                               
            std::vector<bool>      m_isUVBraggMIPAvailableVect;    ///< If the weighted average Bragg MIP over the U & V planes is available
            std::vector<float>     m_braggMIPUVVect;               ///< The weighted average Bragg MIP over the U & V planes
            
            // Bragg likelihood MIP backward fit
            std::vector<bool>      m_isWBraggMIPBackwardAvailableVect;     ///< If the Bragg MIP is available for the W plane fitted backwards
            std::vector<float>     m_braggMIPBackwardWVect;                ///< The Bragg MIP for the W plane fitted backwards
                                                                                                                                               
            std::vector<bool>      m_isUBraggMIPBackwardAvailableVect;     ///< If the Bragg MIP is available for the U plane fitted backwards
            std::vector<float>     m_braggMIPBackwardUVect;                ///< The Bragg MIP for the U plane fitted backwards
                                                                                                                                               
            std::vector<bool>      m_isVBraggMIPBackwardAvailableVect;     ///< If the Bragg MIP is available for the V plane fitted backwards
            std::vector<float>     m_braggMIPBackwardVVect;                ///< The Bragg MIP for the V plane fitted backwards
                                                                                                                                               
            std::vector<bool>      m_isUVBraggMIPBackwardAvailableVect;    ///< If the weighted average Bragg MIP over the U & V planes is available fitted backwards
            std::vector<float>     m_braggMIPBackwardUVVect;               ///< The weighted average Bragg MIP over the U & V planes fitted backwards

            // Bragg likelihood ratio
            std::vector<bool>      m_isWBraggRatioAvailableVect;   ///< If the Bragg ratio is available for the W plane
            std::vector<float>     m_braggRatioWVect;              ///< The Bragg ratio for the W plane
            
            std::vector<bool>      m_isUVBraggRatioAvailableVect;  ///< If the weighted average Bragg ratio over the U & V planes is available
            std::vector<float>     m_braggRatioUVVect;             ///< The weighted average Bragg MIP over the U & V planes
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
      
        // TODO doxygen comments
        void ResetEventTree();
        void SetEventMetadata(const art::Event &event);
        void SetEventTruthInfo(const art::Event &event);
        art::Ptr<simb::MCParticle> SelectMCParticleWithPdgCode(const MCParticleVector &mcParticles, const int pdgCode) const;
        float GetTheta(const art::Ptr<simb::MCParticle> &mcParticle) const;
        float GetPhi(const art::Ptr<simb::MCParticle> &mcParticle) const;
        float GetOpeningAngle(const art::Ptr<simb::MCParticle> &muon, const art::Ptr<simb::MCParticle> &pion) const;
        void SetRecoInfo(const art::Event &event);
        void SetEventRecoInfo(const art::Event &event, const PFParticleVector &allPFParticles, const PFParticleVector &finalStates, const Association<recob::PFParticle, larpandoraobj::PFParticleMetadata> &pfpToMetadata, const spacecharge::SpaceChargeService::provider_type *const pSpaceChargeService);
        void SetPFParticleInfo(const unsigned int index, const art::Ptr<recob::PFParticle> &finalState, const BacktrackHelper::BacktrackerData &backtrackerData, const Association<recob::PFParticle, recob::Track> &pfpToTracks, const Association<recob::Track, anab::ParticleID> &trackToPIDs, const Association<recob::PFParticle, larpandoraobj::PFParticleMetadata> &pfpToMetadata, const spacecharge::SpaceChargeService::provider_type *const pSpaceChargeService);
        void SetPFParticleMCParticleMatchInfo(const art::Ptr<recob::PFParticle> &finalState, const BacktrackHelper::BacktrackerData &backtrackerData);
        void SetPFParticlePandoraInfo(const art::Ptr<recob::PFParticle> &finalState, const BacktrackHelper::BacktrackerData &backtrackerData, const Association<recob::PFParticle, larpandoraobj::PFParticleMetadata> &pfpToMetadata);
        void SetDummyTrackInfo();
        void SetTrackInfo(const art::Ptr<recob::Track> &track, const spacecharge::SpaceChargeService::provider_type *const pSpaceChargeService);
        float GetYZAngle(const TVector3 &dir);
        void SetDummyPIDInfo();
        void SetPIDInfo(const art::Ptr<anab::ParticleID> &pid, const float yzAngle, const int nHitsU, const int nHitsV);
        geo::View_t GetView(const std::bitset<8> &planeMask) const;
        void SetPIDVariables(const geo::View_t &view, const float algoValue, float &wValue, float &uValue, float &vValue, bool &wSet, bool &uSet, bool &vSet) const;
        void CombineInductionPlanes(const float yzAngle, const int nHitsU, const int nHitsV, const float uValue, const float vValue, const bool uSet, const bool vSet, float &uvValue, bool &uvSet) const;
        void ValidateOutputVectorSizes(const unsigned int index) const;

        
        art::EDAnalyzer::Table<Config>  m_config;          ///< The FHiCL configuration options
        
        TTree                          *m_pEventTree;      ///< The output tree for all event level data
        OutputEvent                     m_outputEvent;     ///< The output event-level object
};

} // namespace ubcc1pi

DEFINE_ART_MODULE(ubcc1pi::EventSelection)

#endif
