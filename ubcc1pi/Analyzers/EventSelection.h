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
#include "ubcc1pi/Helpers/TruthHelper.h"

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

            fhicl::Sequence<float> SphereSizes
            {
                fhicl::Name("SphereSizes"),
                fhicl::Comment("The sizes of spheres to use at the end of tracks in which to identify secondary interactions")
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
                                                                           
            // True particle kinematics (for signal events) 
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

            std::vector<int>               m_mcpNHitsUVect;                ///< The number of reconstructed hits that this MCParticle contributed to - U view
            std::vector<int>               m_mcpNHitsVVect;                ///< The number of reconstructed hits that this MCParticle contributed to - V view
            std::vector<int>               m_mcpNHitsWVect;                ///< The number of reconstructed hits that this MCParticle contributed to - W view
            std::vector<int>               m_mcpNGoodHitsUVect;            ///< Only cont hits for which the MCParticle contribute more than 1/2 of the charge - U view
            std::vector<int>               m_mcpNGoodHitsVVect;            ///< Only cont hits for which the MCParticle contribute more than 1/2 of the charge - V view
            std::vector<int>               m_mcpNGoodHitsWVect;            ///< Only cont hits for which the MCParticle contribute more than 1/2 of the charge - W view
            std::vector<float>             m_mcpHitWeightUVect;            ///< The total hit weight that the MCParticle contributed - U view
            std::vector<float>             m_mcpHitWeightVVect;            ///< The total hit weight that the MCParticle contributed - V view
            std::vector<float>             m_mcpHitWeightWVect;            ///< The total hit weight that the MCParticle contributed - W view

            // Reconstructed information                                   
            bool                   m_hasRecoNeutrino;                      ///< If the event has a reconstructed neutrino
            TVector3               m_recoNuVtx;                            ///< The SCE corrected reco neutrino vertex position
            bool                   m_isRecoNuFiducial;                     ///< If the reconstructed vertex is fiducial
//            float                  m_topologicalScore;                     ///< The pandora neutrino ID topological score for this slice
                                                                           
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
            std::vector<float>     m_trueRangeVect;                           ///< The true range of the particle
            std::vector<float>     m_trueMomentumXVect;                    ///< The true momentum of the particle - X
            std::vector<float>     m_trueMomentumYVect;                    ///< The true momentum of the particle - Y
            std::vector<float>     m_trueMomentumZVect;                    ///< The true momentum of the particle - Z
            std::vector<float>     m_trueStartXVect;                       ///< The true start position of the particle - X
            std::vector<float>     m_trueStartYVect;                       ///< The true start position of the particle - Y
            std::vector<float>     m_trueStartZVect;                       ///< The true start position of the particle - Z
            std::vector<bool>      m_trueIsContainedVect;                  ///< If the true particle is contained
            std::vector<bool>      m_trueIsStoppingVect;                   ///< If the true particle stops
            std::vector<int>       m_trueNScattersVect;                    ///< The true number of elastic scatters
            std::vector<bool>      m_trueIsGoldenVect;                     ///< If the true particle is golden
            
            // Pandora info
            std::vector<int>       m_nHitsUVect;                           ///< The number of hits in the U view
            std::vector<int>       m_nHitsVVect;                           ///< The number of hits in the V view
            std::vector<int>       m_nHitsWVect;                           ///< The number of hits in the W view
            std::vector<int>       m_nDescendentHitsUVect;                 ///< The number of descendent hits in the U view
            std::vector<int>       m_nDescendentHitsVVect;                 ///< The number of descendent hits in the V view
            std::vector<int>       m_nDescendentHitsWVect;                 ///< The number of descendent hits in the W view
            std::vector<float>     m_trackShowerVect;                      ///< The pandora track vs. shower score
            std::vector<int>       m_nDaughtersVect;                       ///< The number of daughter PFParticles
            std::vector<int>       m_nDescendentsVect;                     ///< The number of descendents PFParticles follow all downstream links
                                                                           
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
            std::vector<float>     m_xyAngleVect;                          ///< The reconstructed angle in the XY plane
            std::vector<float>     m_xzAngleVect;                          ///< The reconstructed angle in the XZ plane
            std::vector<float>     m_lengthVect;                           ///< The reconstructed length
            std::vector<float>     m_rangeVect;                            ///< The reconstructed range
            std::vector<bool>      m_isContainedVect;                      ///< If the particle is contained within the detector
            std::vector<float>     m_rmsSequentialTrackDeviationVect;      ///< The RMS sin angle of each track segment with the previous segment
            std::vector<float>     m_minSequentialTrackDeviationVect;      ///< The min sin angle of each track segment with the previous segment
            std::vector<float>     m_maxSequentialTrackDeviationVect;      ///< The max sin angle of each track segment with the previous segment
            std::vector<float>     m_rmsTrackDeviationVect;                ///< The RMS sin angle of each track segment with the mean direction
            std::vector<float>     m_minTrackDeviationVect;                ///< The min sin angle of all track segments with the mean direction
            std::vector<float>     m_maxTrackDeviationVect;                ///< The max sin angle of all track segments with the mean direction
    
            std::vector<std::vector<int>*> m_nSpacePointsInSphereAll;      ///< Number of points in a sphere around the track end, each entry of the outer vector is a different radius
            std::vector<std::vector<int>*> m_nOtherSpacePointsInSphereAll; ///< Number of points in sphere not including those from this particle
            std::vector<std::vector<float>*> m_offAxisSpacePointsRMSInSphereAll; ///< The RMS of the distances of the space points in the sphere to the track direction axis
            
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
            
            // Chi2 muon
            std::vector<bool>      m_isWChi2muAvailableVect;                ///< If the chi2 muon is available for the W plane
            std::vector<float>     m_chi2muWVect;                           ///< The chi2 muon for the W plane
                                                                           
            std::vector<bool>      m_isUChi2muAvailableVect;                ///< If the chi2 muon is available for the U plane
            std::vector<float>     m_chi2muUVect;                           ///< The chi2 muon for the U plane
                                                                           
            std::vector<bool>      m_isVChi2muAvailableVect;                ///< If the chi2 muon is available for the V plane
            std::vector<float>     m_chi2muVVect;                           ///< The chi2 muon for the V plane
                                                                           
            std::vector<bool>      m_isUVChi2muAvailableVect;               ///< If the weighted average chi2 muon over the U & V planes is available
            std::vector<float>     m_chi2muUVVect;                          ///< The weighted average chi2 muon over the U & V planes
            
            // Chi2 pion
            std::vector<bool>      m_isWChi2piAvailableVect;                ///< If the chi2 pion is available for the W plane
            std::vector<float>     m_chi2piWVect;                           ///< The chi2 pion for the W plane
                                                                           
            std::vector<bool>      m_isUChi2piAvailableVect;                ///< If the chi2 pion is available for the U plane
            std::vector<float>     m_chi2piUVect;                           ///< The chi2 pion for the U plane
                                                                           
            std::vector<bool>      m_isVChi2piAvailableVect;                ///< If the chi2 pion is available for the V plane
            std::vector<float>     m_chi2piVVect;                           ///< The chi2 pion for the V plane
                                                                           
            std::vector<bool>      m_isUVChi2piAvailableVect;               ///< If the weighted average chi2 pion over the U & V planes is available
            std::vector<float>     m_chi2piUVVect;                          ///< The weighted average chi2 pion over the U & V planes
            
            // Chi2 MIP
            std::vector<bool>      m_isWChi2MIPAvailableVect;                ///< If the chi2 MIP is available for the W plane
            std::vector<float>     m_chi2MIPWVect;                           ///< The chi2 MIP for the W plane
                                                                           
            std::vector<bool>      m_isUChi2MIPAvailableVect;                ///< If the chi2 MIP is available for the U plane
            std::vector<float>     m_chi2MIPUVect;                           ///< The chi2 MIP for the U plane
                                                                           
            std::vector<bool>      m_isVChi2MIPAvailableVect;                ///< If the chi2 MIP is available for the V plane
            std::vector<float>     m_chi2MIPVVect;                           ///< The chi2 MIP for the V plane
                                                                           
            std::vector<bool>      m_isUVChi2MIPAvailableVect;               ///< If the weighted average chi2 MIP over the U & V planes is available
            std::vector<float>     m_chi2MIPUVVect;                          ///< The weighted average chi2 MIP over the U & V planes

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
            
            // Bragg likelihood muon
            std::vector<bool>      m_isWBraggmuAvailableVect;              ///< If the Bragg muon is available for the W plane                        
            std::vector<float>     m_braggmuWVect;                         ///< The Bragg muon for the W plane
                                                                                                                                                      
            std::vector<bool>      m_isUBraggmuAvailableVect;              ///< If the Bragg muon is available for the U plane
            std::vector<float>     m_braggmuUVect;                         ///< The Bragg muon for the U plane
                                                                                                                                                      
            std::vector<bool>      m_isVBraggmuAvailableVect;              ///< If the Bragg muon is available for the V plane
            std::vector<float>     m_braggmuVVect;                         ///< The Bragg muon for the V plane
                                                                                                                                                      
            std::vector<bool>      m_isUVBraggmuAvailableVect;             ///< If the weighted average Bragg muon over the U & V planes is available
            std::vector<float>     m_braggmuUVVect;                        ///< The weighted average Bragg muon over the U & V planes
            
            // Bragg likelihood pion
            std::vector<bool>      m_isWBraggpiAvailableVect;              ///< If the Bragg pion is available for the W plane                        
            std::vector<float>     m_braggpiWVect;                         ///< The Bragg pion for the W plane
                                                                                                                                                      
            std::vector<bool>      m_isUBraggpiAvailableVect;              ///< If the Bragg pion is available for the U plane
            std::vector<float>     m_braggpiUVect;                         ///< The Bragg pion for the U plane
                                                                                                                                                      
            std::vector<bool>      m_isVBraggpiAvailableVect;              ///< If the Bragg pion is available for the V plane
            std::vector<float>     m_braggpiVVect;                         ///< The Bragg pion for the V plane
                                                                                                                                                      
            std::vector<bool>      m_isUVBraggpiAvailableVect;             ///< If the weighted average Bragg pion over the U & V planes is available
            std::vector<float>     m_braggpiUVVect;                        ///< The weighted average Bragg pion over the U & V planes
        };

        /**
         *  @brief  Constructor
         *
         *  @param  config the set of input fhicl parameters
         */
        EventSelection(const art::EDAnalyzer::Table<Config> &config);
        
        /**
         *  @brief  Destructor
         */
        ~EventSelection();

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
         *  @param  mcParticleToHits the mapping from MCParticle to hits along with the backtracker matching data
         */
        void SetMCParticleInfo(const MCParticleVector &allMCParticles, const MCParticleVector &reconstrutableFinalStates, const AssociationData<simb::MCParticle, recob::Hit, anab::BackTrackerHitMatchingData> &mcParticleToHits);

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
         *  @param  mcParticleMap the MCParticle hierarchy map
         *  @param  pfParticleMap the PFParticle hierarchy map
         *  @param  spacePoints all of the spacepoints in the event
         *  @param  pfpToHits the mapping from PFParticle to hits
         *  @param  pfpToSpacePoints the mapping from PFParticles to spacepoints
         *  @param  pfpToTracks the mapping from PFParticles to Tracks
         *  @param  trackToPIDs the mapping from Tracks to PID
         *  @param  trackToCalorimetries the mapping from Tracks to Calorimetry
         *  @param  pfpToMetadata the mapping from PFParticles to metadata
         *  @param  pSpaceChargeService the space charge service
         */
        void SetPFParticleInfo(const unsigned int index, const art::Ptr<recob::PFParticle> &finalState, const BacktrackHelper::BacktrackerData &backtrackerData, const MCParticleMap &mcParticleMap, const PFParticleMap &pfParticleMap, const Collection<recob::SpacePoint> &spacePoints, const Association<recob::PFParticle, recob::Hit> &pfpToHits, const Association<recob::PFParticle, recob::SpacePoint> &pfpToSpacePoints, const Association<recob::PFParticle, recob::Track> &pfpToTracks, const Association<recob::Track, anab::ParticleID> &trackToPIDs, const Association<recob::Track, anab::Calorimetry> &trackToCalorimetries, const Association<recob::PFParticle, larpandoraobj::PFParticleMetadata> &pfpToMetadata, const spacecharge::SpaceChargeService::provider_type *const pSpaceChargeService);

        /**
         *  @brief  Set the MCParticle matching info for the current PFParticle to the output branches
         *
         *  @param  finalState the final state PFParticle to output
         *  @param  backtrackerData the backtracking data
         *  @param  mcParticleMap
         */
        void SetPFParticleMCParticleMatchInfo(const art::Ptr<recob::PFParticle> &finalState, const BacktrackHelper::BacktrackerData &backtrackerData, const MCParticleMap &mcParticleMap);

        /**
         *  @brief  Set the PFParticle info from Pandora in the output branches
         *
         *  @param  finalState the final state PFParticle to output
         *  @param  pfpToHits the mapping from PFParticles to hits
         *  @param  pfParticleMap the PFParticle hierarchy map
         *  @param  pfpToMetadata the mapping from PFParticles to metadata
         */
        void SetPFParticlePandoraInfo(const art::Ptr<recob::PFParticle> &finalState, const Association<recob::PFParticle, recob::Hit> &pfpToHits, const PFParticleMap &pfParticleMap, const Association<recob::PFParticle, larpandoraobj::PFParticleMetadata> &pfpToMetadata);

        /**
         *  @brief  Set dummy track info for the current PFParticle 
         */
        void SetDummyTrackInfo();

        /**
         *  @brief  Set the track info for the current PFParticle
         *
         *  @param  track the track to output
         *  @param  spacePoints all of the spacepoints in the event
         *  @param  spacePointsInParticle the spacepoints that went toward making the track
         *  @param  pSpaceChargeService the space charge service
         */
        void SetTrackInfo(const art::Ptr<recob::Track> &track, const Collection<recob::SpacePoint> &spacePoints, const Collection<recob::SpacePoint> &spacePointsInParticle, const spacecharge::SpaceChargeService::provider_type *const pSpaceChargeService);

        /**
         *  @brief  Get the spacepoints in the input list that are within a given distance of a given point
         *
         *  @param  spacePoints the input spacepoints
         *  @param  point the point which the spacepoints must be near
         *  @param  dist the threshold distance
         *  @param  pSpaceChargeService the spacecharge service
         *
         *  @return the nearby spacepoints
         */
        SpacePointVector GetSpacePointsNearPoint(const SpacePointVector &spacePoints, const TVector3 &point, const float dist, const spacecharge::SpaceChargeService::provider_type *const pSpaceChargeService) const;

        /**
         *  @brief  Get the angle of in the YZ plane from an input direction
         *
         *  @param  dir the input direction
         *
         *  @return the YZ angle
         */
        float GetYZAngle(const TVector3 &dir);
        
        /**
         *  @brief  Get the angle of in the XY plane from an input direction
         *
         *  @param  dir the input direction
         *
         *  @return the XY angle
         */
        float GetXYAngle(const TVector3 &dir);
        
        /**
         *  @brief  Get the angle of in the XZ plane from an input direction
         *
         *  @param  dir the input direction
         *
         *  @return the XZ angle
         */
        float GetXZAngle(const TVector3 &dir);

        /**
         *  @brief  Get the RMS of the sin angle between sequential track segments directions
         *
         *  @param  track the track
         *  @param  validPoints the indices of the valid points
         *
         *  @return the RMS sequential deviation
         */
        float GetRMSSequentialTrackDeviation(const art::Ptr<recob::Track> &track, const std::vector<size_t> &validPoints);

        /**
         *  @brief  Get the minimum of the absolute sin angle between sequential track direction
         *
         *  @param  track the track
         *  @param  validPoints the indices of the valid points
         *
         *  @return the min sequential deviation
         */
        float GetMinSequentialTrackDeviation(const art::Ptr<recob::Track> &track, const std::vector<size_t> &validPoints);
        
        /**
         *  @brief  Get the maximum of the absolute sin angle between sequential track direction
         *
         *  @param  track the track
         *  @param  validPoints the indices of the valid points
         *
         *  @return the max sequential deviation
         */
        float GetMaxSequentialTrackDeviation(const art::Ptr<recob::Track> &track, const std::vector<size_t> &validPoints);

        /**
         *  @brief  Get the mean direction of the track along all valid points
         *
         *  @param  track the track
         *  @param  validPoints the indices of the valid points
         *
         *  @return the mean track direction
         */
        TVector3 GetMeanTrackDirection(const art::Ptr<recob::Track> &track, const std::vector<size_t> &validPoints);

        /**
         *  @brief  Get the RMS of the sin angle between the mean track direction and each valid track segment
         *
         *  @param  track the track
         *  @param  validPoints the indices of the valid points
         *  @param  meanDir the mean track direction
         *
         *  @return the RMS devation
         */
        float GetRMSTrackDeviation(const art::Ptr<recob::Track> &track, const std::vector<size_t> &validPoints, const TVector3 &meanDir);

        /**
         *  @brief  Get the minimum absolute sin squared angle to the mean track direction
         *
         *  @param  track the track
         *  @param  validPoints the indices of the valid points
         *  @param  meanDir the mean track direction
         *
         *  @return the minimum deviation
         */
        float GetMinTrackDeviation(const art::Ptr<recob::Track> &track, const std::vector<size_t> &validPoints, const TVector3 &meanDir);

        /**
         *  @brief  Get the maximum absolute sin squared angle to the mean track direction
         *
         *  @param  track the track
         *  @param  validPoints the indices of the valid points
         *  @param  meanDir the mean track direction
         *
         *  @return the maximum deviation
         */
        float GetMaxTrackDeviation(const art::Ptr<recob::Track> &track, const std::vector<size_t> &validPoints, const TVector3 &meanDir);

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
         *  @brief  Get the range of the input MCParticle
         *
         *  @param  mcParticle the MCParticle
         *
         *  @return the range
         */
        float GetRange(const art::Ptr<simb::MCParticle> &mcParticle) const;

        /**
         *  @brief  Get the range of the input track
         *
         *  @param  track the track
         *  @param  pSpaceChargeService the space charge service
         *
         *  @return the range
         */
        float GetRange(const art::Ptr<recob::Track> &track, const spacecharge::SpaceChargeService::provider_type *const pSpaceChargeService) const;

        /**
         *  @brief  Get the valid points in the track
         *
         *  @param  track the track
         *
         *  @return the indices of the valid points in the track
         */
        std::vector<size_t> GetValidPoints(const art::Ptr<recob::Track> &track) const;

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
