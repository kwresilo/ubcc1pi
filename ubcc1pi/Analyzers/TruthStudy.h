/**
 *  @file  ubcc1pi/Analyzers/TruthStudy.h
 *
 *  @brief The header file for the truth study analyzer.
 */

#ifndef UBCC1PI_ANALYZERS_TRUTH_STUDY
#define UBCC1PI_ANALYZERS_TRUTH_STUDY

#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDAnalyzer.h"

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

        art::EDAnalyzer::Table<Config>  m_config; ///< The FHiCL configuration options
};

} // namespace ubcc1pi

DEFINE_ART_MODULE(ubcc1pi::TruthStudy)

#endif
