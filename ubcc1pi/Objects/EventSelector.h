/**
 *  @file  ubcc1pi/Objects/EventSelector.h
 *
 *  @brief The header file for the event selector class
 */

#ifndef UBCC1PI_OBJECTS_EVENT_SELECTOR
#define UBCC1PI_OBJECTS_EVENT_SELECTOR

#include "ubcc1pi/Helpers/CollectionHelper.h"

#include <TFile.h>
#include <TH1.h>

namespace ubcc1pi
{

/**
 *  @brief  The event selector class
 *
 *          By working under the assumption that the event is truly CC1Pi+ this object determines the most likely PID of each particle type.
 *          It also gives the likelihood that the particle is truly CC1Pi+ using the most likely PID.
 */
class EventSelector
{
    public:
        /**
         *  @brief  Constructor
         *
         *  @param  finalState the final state PFParticles
         *  @param  event the art event
         *  @param  pfParticleLabel the PFParticle producer label
         *  @param  trackLabel the Track producer label
         *  @param  pidLabel the PID producer label
         */
        EventSelector(const PFParticleVector &finalStates, const art::Event &event, const art::InputTag &pfParticleLabel, const art::InputTag &trackLabel, const art::InputTag &pidLabel);

        /**
         *  @breif  Get the track vs. shower score for the given PFParticle
         *
         *  @param  pfParticle the input PFParticle
         */
        float GetTrackShowerScore(const art::Ptr<recob::PFParticle> &pfParticle) const;
        
        /**
         *  @breif  Get the proton vs. MIP score for the given PFParticle
         *
         *  @param  pfParticle the input PFParticle
         */
        float GetProtonMIPScore(const art::Ptr<recob::PFParticle> &pfParticle) const;
        
        /**
         *  @breif  Get the muon vs. pion score for the given PFParticle
         *
         *  @param  pfParticle the input PFParticle
         */
        float GetMuonPionScore(const art::Ptr<recob::PFParticle> &pfParticle) const;

        /**
         *  @breif  Get the likelihood that the pfParticle would have the relevant scores, given that it's truly a muon
         *
         *  @param  pfParticle the input PFParticle
         */
        float GetMuonLikelihood(const art::Ptr<recob::PFParticle> &pfParticle) const;
        
        /**
         *  @breif  Get the likelihood that the pfParticle would have the relevant scores, given that it's truly a pion
         *
         *  @param  pfParticle the input PFParticle
         */
        float GetPionLikelihood(const art::Ptr<recob::PFParticle> &pfParticle) const;
        
        /**
         *  @breif  Get the likelihood that the pfParticle would have the relevant scores, given that it's truly a proton
         *
         *  @param  pfParticle the input PFParticle
         */
        float GetProtonLikelihood(const art::Ptr<recob::PFParticle> &pfParticle) const;

        /**
         *  @breif  Get the likelihood of the event being CC1Pi with the most likely PID assignment (assuming the event is truly CC1Pi)
         *
         *  @param  bestMuon the pfparticle most likely to be the muon assuming the event is cc1pi
         *  @param  bestPion the pfparticle most likely to be the pion assuming the event is cc1pi
         *  @param  bestProtons the pfparticles most likely to be protons assuming the event is cc1pi
         */
        float GetBestCC1PiScore(art::Ptr<recob::PFParticle> &bestMuon, art::Ptr<recob::PFParticle> &bestPion, PFParticleVector &bestProtons) const;

    private:

        /**
         *  @brief  Get the total likelihood for a given PFParticle using the supplied histograms
         *
         *  @param  pfParticle the input PFParticle
         *  @param  pTrackShower the track-shower likelihood histogram
         *  @param  pProtonMIP the proton-MIP likelihood histogram
         *  @param  pMuonPion the muon-pion likelihood histogram
         */
        float GetLikelihood(const art::Ptr<recob::PFParticle> &pfParticle, const TH1F *pTrackShower, const TH1F *pProtonMIP, const TH1F *pMuonPion) const;
        
        /**
         *  @brief  Get the single likelihood from the supplied histogram
         *
         *  @param  value the value to look up in the histogram
         *  @param  pHistogram the input likelihood histogram
         */
        float GetLikelihood(const float value, const TH1F *pHistogram) const;

        /**
         *  @brief  Get the full file path for the PID histogram root file
         */
        std::string GetPIDHistFileName() const;

        PFParticleVector                                                   m_finalStates;           ///< The input final state PFParticles
        const art::Event                                                  *m_pEvent;                ///< The art event

        Association<recob::PFParticle, larpandoraobj::PFParticleMetadata>  m_pfParticleToMetadata;  ///< The mapping from PFParticles to metadata
        Association<recob::PFParticle, recob::Track>                       m_pfParticleToTrack;     ///< The mapping from PFParticles to tracks
        Association<recob::Track, anab::ParticleID>                        m_trackToPID;            ///< The mapping from Tracks to PID

        TH1F                                                              *m_hTrackShower_muon;     ///< The track-shower likelihood histogram for muons
        TH1F                                                              *m_hTrackShower_pion;     ///< The track-shower likelihood histogram for pions
        TH1F                                                              *m_hTrackShower_proton;   ///< The track-shower likelihood histogram for protons
        
        TH1F                                                              *m_hProtonMIP_muon;       ///< The proton-MIP likelihood histogram for muons
        TH1F                                                              *m_hProtonMIP_pion;       ///< The proton-MIP likelihood histogram for pions
        TH1F                                                              *m_hProtonMIP_proton;     ///< The proton-MIP likelihood histogram for protons
        
        TH1F                                                              *m_hMuonPion_muon;        ///< The muon-pion likelihood histogram for muons
        TH1F                                                              *m_hMuonPion_pion;        ///< The muon-pion likelihood histogram for pions
        TH1F                                                              *m_hMuonPion_proton;      ///< The muon-pion likelihood histogram for protons
};

} // namespace ubcc1pi

#endif
