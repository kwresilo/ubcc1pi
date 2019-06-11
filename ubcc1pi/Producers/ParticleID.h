/**
 *  @file  ubcc1pi/Producers/ParticleID.h
 *
 *  @brief The header file for the particle ID producer.
 */

#ifndef UBCC1PI_PRODUCERS_PARTICLE_ID
#define UBCC1PI_PRODUCERS_PARTICLE_ID

#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDProducer.h"
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
 *  @brief  The particle ID producer class
 */
class ParticleID : public art::EDProducer
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
                fhicl::Comment("The label for the Calorimetry producer")
            };
        };
        
        /**
         *  @brief  The output calorimetry point level structure
         */
        struct OutputPoint
        {
            // Event metadata
            int          m_run;                   ///< The run number
            int          m_subRun;                ///< The subrun number
            int          m_event;                 ///< The event number

            // Particle metadata
            int          m_pfParticleId;          ///< The PFParticle ID
            bool         m_hasMatchedMCParticle;  ///< If the PFParticle has a matched MCParticle
            int          m_truePdgCode;           ///< The particle PDG code
            float        m_trueMomentum;          ///< The particle momentum
            float        m_trueMatchPurity;       ///< Match purity to the true particle
            float        m_trueMatchCompleteness; ///< Match completeness to the true particle

            // Track metadata
            float        m_trackTheta;            ///< Angle of track to z-axis
            float        m_trackPhi;              ///< Azimuthal angle of track around z-axis

            // Calorimetry point information
            int          m_view;                  ///< The plane from which this point was created
            float        m_residualRange;         ///< The residual range of the point
            float        m_dEdx;                  ///< The dEdx of the point
        };

        /**
         *  @brief  Constructor
         *
         *  @param  config the set of input fhicl parameters
         */
        ParticleID(const art::EDProducer::Table<Config> &config);

        /**
         *  @brief  Produce collections in the event record
         *
         *  @param  event the event to add to
         */
        void produce(art::Event &event);

    private:
        std::vector<std::pair<float, float> > GetCalorimetryPoints(const CalorimetryVector &calos) const;
        std::vector<std::pair<float, float> > GetCalorimetryPoints(const art::Ptr<anab::Calorimetry> &calo) const;
        void SortPoints(std::vector< std::pair<float, float> > &points) const;
        std::vector<std::pair<float, float> > GetPointsInRange(const std::vector< std::pair<float, float> > &points, const float min, const float max) const;
        float GetMeandEdx(const std::vector<std::pair<float, float> > &points) const;
        float GetMostProbabledEdx(const std::vector<std::pair<float, float> > &points, const float kernalWidth, const float maxStepSize) const;
        void SetMatchedMCParticleInfo(const art::Ptr<recob::PFParticle> &pfParticle, const BacktrackHelper::BacktrackerData &backtrackerData);
        void OutputPoints(const CalorimetryVector &calos);
            
        art::EDProducer::Table<Config>  m_config;          ///< The FHiCL configuration options
        
        TTree                          *m_pPointTree;      ///< The output tree for all calorimetry-point level data
        OutputPoint                     m_outputPoint;     ///< The output calorimetry-point level object
};

} // namespace ubcc1pi

DEFINE_ART_MODULE(ubcc1pi::ParticleID)

#endif
