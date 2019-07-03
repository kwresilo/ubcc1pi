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
            
            fhicl::Atom<art::InputTag> CalorimetryLabel
            {
                fhicl::Name("CalorimetryLabel"),
                fhicl::Comment("The label for the calorimetry producer")
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
            int          m_run;                                            ///< The run number
            int          m_subRun;                                         ///< The subrun number
            int          m_event;                                          ///< The event number
                                                                           
            // The truth interaction information                           
            bool         m_isSignal;                                       ///< If the event is a true CC1Pi signal
            std::string  m_interaction;                                    ///< The interaction type string
            float        m_trueNuE;                                        ///< The true neutrino energy
            TVector3     m_trueNuVtx;                                      ///< The true neutrino vertex position
            bool         m_isTrueNuFiducial;                               ///< If the neutrino interaction is fiducial
                                                                           
            // The true particle multiplicities                            
            int          m_nMuMinus;                                       ///< The number of mu-
            int          m_nMuPlus;                                        ///< The number of mu+
            int          m_nPiPlus;                                        ///< The number of pi+
            int          m_nPiMinus;                                       ///< The number of pi-
            int          m_nKPlus;                                         ///< The number of K+
            int          m_nKMinus;                                        ///< The number of K-
            int          m_nProton;                                        ///< The number of p
            int          m_nNeutron;                                       ///< The number of n
            int          m_nPhoton;                                        ///< The number of gamma
            int          m_nElectron;                                      ///< The number of e
            int          m_nPositron;                                      ///< The number of e+
            int          m_nTotal;                                         ///< The total number of particles
                                                                           
            // True particle kinematics                                    
            float        m_trueMuEnergy;                                   ///< The true muon energy
            float        m_trueMuTheta;                                    ///< The true muon theta (angle to beam)
            float        m_trueMuPhi;                                      ///< The true muon phi (angle about beam)
                                                                           
            float        m_truePiEnergy;                                   ///< The true pion energy
            float        m_truePiTheta;                                    ///< The true pion theta
            float        m_truePiPhi;                                      ///< The true pion phi
                                                                           
            float        m_trueMuPiAngle;                                  ///< The muon-pion opening angle
            
            // MCParticles
            int                            m_nMCParticles;                 ///< The number of MCParticles
            std::vector<int>               m_mcpIdVect;                    ///< The unique ID of the MCParticle
            std::vector<int>               m_mcpPDGVect;                   ///< The PDG code of the MCParticle
            std::vector<bool>              m_mcpIsTargetFinalStateVect;    ///< If the MCParticle is a neutrino final state passing the reconstructability momentum thresholds
            std::vector<std::string>       m_mcpProcessVect;               ///< The process name
            std::vector<int>               m_mcpMotherIdVect;              ///< The unique ID of the mother of the MCParticle
            std::vector<int>               m_mcpNDaughtersVect;            ///< The number of daughter MCParticles
            std::vector<std::vector<int> > m_mcpDaughterIdsVect;           ///< The unique IDs of the daughter MCParticles
            std::vector<float>             m_mcpEnergyVect;                ///< The energy of the MCParticle
            std::vector<float>             m_mcpMomentumVect;              ///< The total momentum of the MCParticle
            std::vector<float>             m_mcpMomentumXVect;             ///< The momentum of the MCParticle - X
            std::vector<float>             m_mcpMomentumYVect;             ///< The momentum of the MCParticle - Y
            std::vector<float>             m_mcpMomentumZVect;             ///< The momentum of the MCParticle - Z

            // TODO Fill in the MCParticle -> Hit variables by updating the backtracker helper
            /*
            std::vector<int>               m_mcpNHitsUVect;                ///< The number of reconstructed hits that this MCParticle contributed to - U view
            std::vector<int>               m_mcpNHitsVVect;                ///< The number of reconstructed hits that this MCParticle contributed to - V view
            std::vector<int>               m_mcpNHitsWVect;                ///< The number of reconstructed hits that this MCParticle contributed to - W view
            std::vector<float>             m_mcpNHitsUWeightVect;          ///< Same as NHits but weighted by the fraction of the hit's charge that came from this MCParticle - U view
            std::vector<float>             m_mcpNHitsVWeightVect;          ///< Same as NHits but weighted by the fraction of the hit's charge that came from this MCParticle - V view
            std::vector<float>             m_mcpNHitsWWeightVect;          ///< Same as NHits but weighted by the fraction of the hit's charge that came from this MCParticle - W view
            std::vector<bool>              m_mcpIsReconstructable;         ///< If the MCParticle is deemed reconstructable based on how many hits it contributes to
            */

            // Reconstructed information                                   
            bool                   m_hasRecoNeutrino;                      ///< If the event has a reconstructed neutrino
            TVector3               m_recoNuVtx;                            ///< The SCE corrected reco neutrino vertex position
            bool                   m_isRecoNuFiducial;                     ///< If the reconstructed vertex is fiducial
            float                  m_topologicalScore;                     ///< The pandora neutrino ID topological score for this slice
                                                                           
            // Particle level information                                  
            int                    m_nFinalStatePFPs;                      ///< The number of reconstructed final state PFParticles
            
            // Truth matching
            std::vector<bool>      m_hasMatchedMCParticleVect;             ///< If the particles have a matched MCParticle
            std::vector<int>       m_matchedMCParticleIdVect;              ///< The ID of the matched MCParticle
            std::vector<int>       m_truePdgCodeVect;                      ///< The best matched MCParticles PDG code
            std::vector<float>     m_truthMatchCompletenessVect;           ///< The completeness of the match
            std::vector<float>     m_truthMatchPurityVect;                 ///< The purity of the match
            std::vector<float>     m_trueEnergyVect;                       ///< The true energy of the particle
            std::vector<float>     m_trueKEVect;                           ///< The true kinetic energy of the particle
            std::vector<float>     m_trueMomentumXVect;                    ///< The true momentum of the particle - X
            std::vector<float>     m_trueMomentumYVect;                    ///< The true momentum of the particle - Y
            std::vector<float>     m_trueMomentumZVect;                    ///< The true momentum of the particle - Z
            std::vector<float>     m_trueStartXVect;                       ///< The true start position of the particle - X
            std::vector<float>     m_trueStartYVect;                       ///< The true start position of the particle - Y
            std::vector<float>     m_trueStartZVect;                       ///< The true start position of the particle - Z
            
            // Pandora info
            std::vector<int>       m_nHitsUVect;                           ///< The number of hits in the U view
            std::vector<int>       m_nHitsVVect;                           ///< The number of hits in the V view
            std::vector<int>       m_nHitsWVect;                           ///< The number of hits in the W view
            std::vector<float>     m_trackShowerVect;                      ///< The pandora track vs. shower score
                                                                           
            // Track info                                                  
            std::vector<bool>      m_hasTrackInfoVect;                     ///< If the PFParticle has an associated track
            std::vector<float>     m_startXVect;                           ///< The SCE corrected reconstructed start position of the particle - X
            std::vector<float>     m_startYVect;                           ///< The SCE corrected reconstructed start position of the particle - Y
            std::vector<float>     m_startZVect;                           ///< The SCE corrected reconstructed start position of the particle - Z
            std::vector<float>     m_endXVect;                             ///< The SCE corrected reconstructed end position of the particle - X
            std::vector<float>     m_endYVect;                             ///< The SCE corrected reconstructed end position of the particle - Y
            std::vector<float>     m_endZVect;                             ///< The SCE corrected reconstructed end position of the particle - Z
            std::vector<float>     m_directionXVect;                       ///< The reconstructed direction of the particle - X
            std::vector<float>     m_directionYVect;                       ///< The reconstructed direction of the particle - Y
            std::vector<float>     m_directionZVect;                       ///< The reconstructed direction of the particle - Z
            std::vector<float>     m_thetaVect;                            ///< The reconstructed angle to the beam direction
            std::vector<float>     m_phiVect;                              ///< The reconstructed angle aroun the beam direction
            std::vector<float>     m_yzAngleVect;                          ///< The reconstructed angle in the YZ plane to the vertical
            std::vector<float>     m_lengthVect;                           ///< The reconstructed length
            std::vector<bool>      m_isContainedVect;                      ///< If the particle is contained within the detector

            // Calorimetry info
            std::vector<bool>                m_hasCalorimetryInfoVect;     ///< If the PFParticle has an associated Calorimetry object
            std::vector<std::vector<float> > m_dedxPerHitUVect;            ///< The dEdx at each hit point along the tracjectory in the U plane
            std::vector<std::vector<float> > m_dedxPerHitVVect;            ///< The dEdx at each hit point along the tracjectory in the V plane
            std::vector<std::vector<float> > m_dedxPerHitWVect;            ///< The dEdx at each hit point along the tracjectory in the W plane
            std::vector<std::vector<float> > m_residualRangePerHitUVect;   ///< The residual range at each hit point along the tracjectory in the U plane
            std::vector<std::vector<float> > m_residualRangePerHitVVect;   ///< The residual range at each hit point along the tracjectory in the V plane
            std::vector<std::vector<float> > m_residualRangePerHitWVect;   ///< The residual range at each hit point along the tracjectory in the W plane

            // PID info
            std::vector<bool>      m_hasPIDInfoVect;                       ///< If the PFParticle has an associated PID object

            // Chi2 proton
            std::vector<bool>      m_isWChi2pAvailableVect;                ///< If the chi2 proton is available for the W plane
            std::vector<float>     m_chi2pWVect;                           ///< The chi2 proton for the W plane
                                                                           
            std::vector<bool>      m_isUChi2pAvailableVect;                ///< If the chi2 proton is available for the U plane
            std::vector<float>     m_chi2pUVect;                           ///< The chi2 proton for the U plane
                                                                           
            std::vector<bool>      m_isVChi2pAvailableVect;                ///< If the chi2 proton is available for the V plane
            std::vector<float>     m_chi2pVVect;                           ///< The chi2 proton for the V plane
                                                                           
            std::vector<bool>      m_isUVChi2pAvailableVect;               ///< If the weighted average chi2 proton over the U & V planes is available
            std::vector<float>     m_chi2pUVVect;                          ///< The weighted average chi2 proton over the U & V planes

            // Bragg likelihood proton
            std::vector<bool>      m_isWBraggpAvailableVect;               ///< If the Bragg proton is available for the W plane                        
            std::vector<float>     m_braggpWVect;                          ///< The Bragg proton for the W plane
                                                                                                                                                       
            std::vector<bool>      m_isUBraggpAvailableVect;               ///< If the Bragg proton is available for the U plane
            std::vector<float>     m_braggpUVect;                          ///< The Bragg proton for the U plane
                                                                                                                                                       
            std::vector<bool>      m_isVBraggpAvailableVect;               ///< If the Bragg proton is available for the V plane
            std::vector<float>     m_braggpVVect;                          ///< The Bragg proton for the V plane
                                                                                                                                                       
            std::vector<bool>      m_isUVBraggpAvailableVect;              ///< If the weighted average Bragg proton over the U & V planes is available
            std::vector<float>     m_braggpUVVect;                         ///< The weighted average Bragg proton over the U & V planes
            
            // Bragg likelihood proton backward fit
            std::vector<bool>      m_isWBraggpBackwardAvailableVect;       ///< If the backward Bragg proton is available for the W plane                        
            std::vector<float>     m_braggpBackwardWVect;                  ///< The backward Bragg proton for the W plane
                                                                                                                                               
            std::vector<bool>      m_isUBraggpBackwardAvailableVect;       ///< If the backward Bragg proton is available for the U plane
            std::vector<float>     m_braggpBackwardUVect;                  ///< The backward Bragg proton for the U plane
                                                                                                                                               
            std::vector<bool>      m_isVBraggpBackwardAvailableVect;       ///< If the backward Bragg proton is available for the V plane
            std::vector<float>     m_braggpBackwardVVect;                  ///< The backward Bragg proton for the V plane
                                                                                                                                               
            std::vector<bool>      m_isUVBraggpBackwardAvailableVect;      ///< If the weighted average backward Bragg proton over the U & V planes is available
            std::vector<float>     m_braggpBackwardUVVect;                 ///< The weighted average backward Bragg proton over the U & V planes
            
            // Bragg likelihood MIP
            std::vector<bool>      m_isWBraggMIPAvailableVect;             ///< If the Bragg MIP is available for the W plane                        
            std::vector<float>     m_braggMIPWVect;                        ///< The Bragg MIP for the W plane
                                                                                                                                                       
            std::vector<bool>      m_isUBraggMIPAvailableVect;             ///< If the Bragg MIP is available for the U plane
            std::vector<float>     m_braggMIPUVect;                        ///< The Bragg MIP for the U plane
                                                                                                                                                       
            std::vector<bool>      m_isVBraggMIPAvailableVect;             ///< If the Bragg MIP is available for the V plane
            std::vector<float>     m_braggMIPVVect;                        ///< The Bragg MIP for the V plane
                                                                                                                                                       
            std::vector<bool>      m_isUVBraggMIPAvailableVect;            ///< If the weighted average Bragg MIP over the U & V planes is available
            std::vector<float>     m_braggMIPUVVect;                       ///< The weighted average Bragg MIP over the U & V planes
            
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
            std::vector<bool>      m_isWBraggRatioAvailableVect;           ///< If the Bragg ratio is available for the W plane
            std::vector<float>     m_braggRatioWVect;                      ///< The Bragg ratio for the W plane
                                                                          
            std::vector<bool>      m_isUVBraggRatioAvailableVect;          ///< If the weighted average Bragg ratio over the U & V planes is available
            std::vector<float>     m_braggRatioUVVect;                     ///< The weighted average Bragg MIP over the U & V planes
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

        /**
         *  @brief  Reset the output branches for the next event
         */
        void ResetEventTree();

        /**
         *  @brief  Set the event metadata in the output branches (event number, run number, etc...)_
         *
         *  @param  event the art event
         */
        void SetEventMetadata(const art::Event &event);

        /**
         *  @brief  Set the event level truth info about the neutrino interaction
         *
         *  @param  event the art event
         */
        void SetEventTruthInfo(const art::Event &event);

        /**
         *  @brief  Set the information about all MCParticles in the event
         *
         *  @param  allMCParticles all of the MCParticles in the event
         *  @param  reconstrutableFinalStates the final state MCParticles that pass the momentum thresholds
         */
        void SetMCParticleInfo(const MCParticleVector &allMCParticles, const MCParticleVector &reconstrutableFinalStates);

        /**
         *  @brief  Select the MCParticle with the given PDG code (assumes each PDG is in the input vector at most once)
         *
         *  @param  mcParticles the input vector of MCParticles from which to selected the chosen particle
         *  @param  pdgCode the PDG code to select
         *
         *  @return the MCParticle in the input vector with the desired PDG code
         *  @throws if there wasn't exactly one MCParticle with the desired PDG code
         */
        art::Ptr<simb::MCParticle> SelectMCParticleWithPdgCode(const MCParticleVector &mcParticles, const int pdgCode) const;

        /**
         *  @brief  Get the angle to the beam direction of the momentum vector of the input MCParticle
         *
         *  @param  mcParticle the input MCParticle
         *
         *  @return the theta angle
         */
        float GetTheta(const art::Ptr<simb::MCParticle> &mcParticle) const;

        /**
         *  @brief  Get the angle around the beam direction of the momentum vector of the input MCParticle
         *
         *  @param  mcParticle the input MCParticle
         *
         *  @return the phi angle
         */
        float GetPhi(const art::Ptr<simb::MCParticle> &mcParticle) const;

        /**
         *  @brief  Get the opening angle between two input MCParticles (muon & pion)
         *
         *  @param  muon the muon MCParticle
         *  @param  pion the pion MCParticle
         *
         *  @return the opening angle between the particles (radians)
         */
        float GetOpeningAngle(const art::Ptr<simb::MCParticle> &muon, const art::Ptr<simb::MCParticle> &pion) const;

        /**
         *  @brief  Set the reconstructed information about the event in the output branches
         *
         *  @param  event the art event
         */
        void SetRecoInfo(const art::Event &event);

        /**
         *  @brief  Set the event-level reconstructed information in the output branches
         *
         *  @param  event the art event
         *  @param  allPFParticles the complete vector of all PFParticles
         *  @param  finalStates the neutrino final state PFParticles
         *  @param  pfpToMetadata the mapping from PFParticles to metadata
         *  @param  pSpaceChargeService the space charge service
         */
        void SetEventRecoInfo(const art::Event &event, const PFParticleVector &allPFParticles, const PFParticleVector &finalStates, const Association<recob::PFParticle, larpandoraobj::PFParticleMetadata> &pfpToMetadata, const spacecharge::SpaceChargeService::provider_type *const pSpaceChargeService);

        /**
         *  @brief  Set the particle-level reconstructed information in the output branches
         *
         *  @param  index the index of the final state PFParticle 
         *  @param  finalState the final state PFParticle to output
         *  @param  backtrackerData the backtracking data
         *  @param  pfpToTracks the mapping from PFParticles to Tracks
         *  @param  trackToPIDs the mapping from Tracks to PID
         *  @param  trackToCalorimetries the mapping from Tracks to Calorimetry
         *  @param  pfpToMetadata the mapping from PFParticles to metadata
         *  @param  pSpaceChargeService the space charge service
         */
        void SetPFParticleInfo(const unsigned int index, const art::Ptr<recob::PFParticle> &finalState, const BacktrackHelper::BacktrackerData &backtrackerData, const Association<recob::PFParticle, recob::Track> &pfpToTracks, const Association<recob::Track, anab::ParticleID> &trackToPIDs, const Association<recob::Track, anab::Calorimetry> &trackToCalorimetries, const Association<recob::PFParticle, larpandoraobj::PFParticleMetadata> &pfpToMetadata, const spacecharge::SpaceChargeService::provider_type *const pSpaceChargeService);

        /**
         *  @brief  Set the MCParticle matching info for the current PFParticle to the output branches
         *
         *  @param  finalState the final state PFParticle to output
         *  @param  backtrackerData the backtracking data
         */
        void SetPFParticleMCParticleMatchInfo(const art::Ptr<recob::PFParticle> &finalState, const BacktrackHelper::BacktrackerData &backtrackerData);

        /**
         *  @brief  Set the PFParticle info from Pandora in the output branches
         *
         *  @param  finalState the final state PFParticle to output
         *  @param  backtrackerData the backtracking data
         *  @param  pfpToMetadata the mapping from PFParticles to metadata
         */
        void SetPFParticlePandoraInfo(const art::Ptr<recob::PFParticle> &finalState, const BacktrackHelper::BacktrackerData &backtrackerData, const Association<recob::PFParticle, larpandoraobj::PFParticleMetadata> &pfpToMetadata);

        /**
         *  @brief  Set dummy track info for the current PFParticle 
         */
        void SetDummyTrackInfo();

        /**
         *  @brief  Set the track info for the current PFParticle
         *
         *  @param  track the track to output
         *  @param  pSpaceChargeService the space charge service
         */
        void SetTrackInfo(const art::Ptr<recob::Track> &track, const spacecharge::SpaceChargeService::provider_type *const pSpaceChargeService);

        /**
         *  @brief  Get the angle of in the YZ plane from an input direction
         *
         *  @param  dir the input direction
         *
         *  @return the YZ angle
         */
        float GetYZAngle(const TVector3 &dir);

        /**
         *  @brief  Set dummy calorimetry info for the current PFParticle
         */
        void SetDummyCalorimetryInfo();

        /**
         *  @brief  Set the calorimetry info for the current PFParticle
         *
         *  @param  calos the input vector of calorimetry object
         */
        void SetCalorimetryInfo(const CalorimetryVector &calos);

        /**
         *  @brief  Set dummy PID info for the current PFParticle
         */
        void SetDummyPIDInfo();

        /**
         *  @brief  Set the PID info for the current PFParticle
         *
         *  @param  pid the PID to output
         *  @param  yzAngle the yz angle of the track
         *  @param  nHitsU the number of hits in the U plane
         *  @param  nHitsV the number of hits in the V plane
         */
        void SetPIDInfo(const art::Ptr<anab::ParticleID> &pid, const float yzAngle, const int nHitsU, const int nHitsV);

        /**
         *  @brief  Get the view from the PID plane mask
         *
         *  @param  planeMask the plane mask from the PID
         *
         *  @return the view represented by the plane mask
         */
        geo::View_t GetView(const std::bitset<8> &planeMask) const;

        /**
         *  @brief  Shorthand to set the PID variables for the given view in the output branches
         *
         *  @param  view the view to set (only the fields for this view will be populated)
         *  @param  algoValue the value to set
         *  @param  wValue the output value in the W plane
         *  @param  uValue the output value in the U plane
         *  @param  vValue the output value in the V plane
         *  @param  wSet if the value in the W plane has been set
         *  @param  uSet if the value in the U plane has been set
         *  @param  vSet if the value in the V plane has been set
         */
        void SetPIDVariables(const geo::View_t &view, const float algoValue, float &wValue, float &uValue, float &vValue, bool &wSet, bool &uSet, bool &vSet) const;

        /**
         *  @brief  Combine information from the induction planes using a YZ angle dependent weighted sum
         *
         *  @param  yzAngle the YZ angle of the track
         *  @param  nHitsU the number of hits in the U plane
         *  @param  nHitsV the number of hits in the V plane
         *  @param  uValue the value of the variable in the U plane
         *  @param  vValue the value of the variable in the V plane
         *  @param  uSet if the U plane variable is available
         *  @param  vSet if the V plane variable is available
         *  @param  uvValue the ouput combined UV variable
         *  @param  uvSet if the output combined UV variable is available
         */
        void CombineInductionPlanes(const float yzAngle, const int nHitsU, const int nHitsV, const float uValue, const float vValue, const bool uSet, const bool vSet, float &uvValue, bool &uvSet) const;

        /**
         *  @brief  Sanity check that all of the output vectors have the expected size at this point
         *
         *  @param  index the index of the PFParticle we have just filled
         */
        void ValidateOutputVectorSizes(const unsigned int index) const;

        
        art::EDAnalyzer::Table<Config>  m_config;          ///< The FHiCL configuration options
        
        TTree                          *m_pEventTree;      ///< The output tree for all event level data
        OutputEvent                     m_outputEvent;     ///< The output event-level object
};

} // namespace ubcc1pi

DEFINE_ART_MODULE(ubcc1pi::EventSelection)

#endif
