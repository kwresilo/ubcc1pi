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

#include "ubcc1pi/Helpers/CollectionHelper.h"
#include "ubcc1pi/Helpers/BacktrackHelper.h"

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
            
            fhicl::Atom<art::InputTag> CalorimetryLabel
            {
                fhicl::Name("CalorimetryLabel"),
                fhicl::Comment("The label for the calorimetry producer")
            };
        };

        /**
         *  @brief  The output PID algorithm level structure
         */
        struct OutputParticle
        {
            // Event metadata
            int          m_run;                   ///< The run number
            int          m_subRun;                ///< The subrun number
            int          m_event;                 ///< The event number
            bool         m_isSignal;              ///< If the event is a true CC1Pi signal

            TVector3     m_nuVertex;              ///< The reconstructed neutrino vertex
            TVector3     m_nuVertexCorrected;     ///< The SCE corrected reconstructed neutrino vertex
        
            // Matched True Particle details
            bool         m_hasMatchedMCParticle;  ///< If the PFParticle matches to any neutrino MCParticles
            int          m_truePdgCode;           ///< The particle PDG code
            float        m_trueMomentum;          ///< The particle momentum
            float        m_trueMatchPurity;       ///< Match purity to the true particle
            float        m_trueMatchCompleteness; ///< Match completeness to the true particle

            // Geometric variables
            TVector3     m_start;                 ///< The start position of the particle
            TVector3     m_end;                   ///< The end position of the particle
            TVector3     m_startCorrected;        ///< The SCE corrected start position of the particle
            TVector3     m_endCorrected;          ///< The SCE corrected end position of the particle

            // PID variables
            int          m_nHitsU;                ///< The number of hits in the U view
            int          m_nHitsV;                ///< The number of hits in the U view
            int          m_nHitsW;                ///< The number of hits in the U view
            float        m_length;                ///< The track-length
            float        m_trackShower;           ///< The track-shower score
            float        m_mipFraction;           ///< The fraction of dEdx sample points that are in the MIP region
            float        m_primaryFraction;       ///< The fraction of hits in the PFParticle's hierarchy that are in the primary itself
            float        m_chi2_mu_U;             ///< The chi2 under the muon hypothesis in the U view
            float        m_chi2_pi_U;             ///< The chi2 under the pion hypothesis in the U view
            float        m_chi2_p_U;              ///< The chi2 under the proton hypothesis in the U view
            float        m_chi2_mu_V;             ///< The chi2 under the muon hypothesis in the V view
            float        m_chi2_pi_V;             ///< The chi2 under the pion hypothesis in the V view
            float        m_chi2_p_V;              ///< The chi2 under the proton hypothesis in the V view
            float        m_chi2_mu_W;             ///< The chi2 under the muon hypothesis in the W view
            float        m_chi2_pi_W;             ///< The chi2 under the pion hypothesis in the W view
            float        m_chi2_p_W;              ///< The chi2 under the proton hypothesis in the W view
            float        m_braggPeakLLH_mu_U;     ///< The forwards bragg peak likelihood under the muon hypothesis in the U view
            float        m_braggPeakLLH_pi_U;     ///< The forwards bragg peak likelihood under the pion hypothesis in the U view
            float        m_braggPeakLLH_p_U;      ///< The forwards bragg peak likelihood under the proton hypothesis in the U view
            float        m_braggPeakLLH_MIP_U;    ///< The forwards bragg peak likelihood under the MIP hypothesis in the U view 
            float        m_braggPeakLLH_mu_V;     ///< The forwards bragg peak likelihood under the muon hypothesis in the V view
            float        m_braggPeakLLH_pi_V;     ///< The forwards bragg peak likelihood under the pion hypothesis in the V view
            float        m_braggPeakLLH_p_V;      ///< The forwards bragg peak likelihood under the proton hypothesis in the V view
            float        m_braggPeakLLH_MIP_V;    ///< The forwards bragg peak likelihood under the MIP hypothesis in the V view 
            float        m_braggPeakLLH_mu_W;     ///< The forwards bragg peak likelihood under the muon hypothesis in the W view
            float        m_braggPeakLLH_pi_W;     ///< The forwards bragg peak likelihood under the pion hypothesis in the W view
            float        m_braggPeakLLH_p_W;      ///< The forwards bragg peak likelihood under the proton hypothesis in the W view
            float        m_braggPeakLLH_MIP_W;    ///< The forwards bragg peak likelihood under the MIP hypothesis in the W view 
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

        void SetEventInfo(const art::Event &event, const PFParticleVector &allPFParticles);
        void ResetParticleInfo();
        void SetMatchedMCParticleInfo(const art::Ptr<recob::PFParticle> &pfParticle, const BacktrackHelper::BacktrackerData &backtrackerData);
        bool GetTrack(const art::Ptr<recob::PFParticle> &pfParticle, const Association<recob::PFParticle, recob::Track> &pfpToTrack, art::Ptr<recob::Track> &track);
        geo::View_t GetView(const std::bitset<8> &planeMask) const;
            
        art::EDAnalyzer::Table<Config>  m_config;          ///< The FHiCL configuration options
        
        TTree                          *m_pParticleTree;  ///< The output tree for all particle level data
        OutputParticle                  m_outputParticle; ///< The output particle-level object
};

} // namespace ubcc1pi

DEFINE_ART_MODULE(ubcc1pi::PIDStudy)

#endif
