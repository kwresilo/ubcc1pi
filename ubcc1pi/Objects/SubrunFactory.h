/**
 *  @file  ubcc1pi/Objects/SubrunFactory.h
 *
 *  @brief The header file for the subrun factory class
 */

#ifndef UBCC1PI_OBJECTS_SUBRUN_FACTORY
#define UBCC1PI_OBJECTS_SUBRUN_FACTORY

#include "ubcc1pi_standalone/Interface/Subrun.h"
#include "art/Framework/Principal/SubRun.h"

#include "fhiclcpp/types/Atom.h"
#include "canvas/Utilities/InputTag.h"

namespace ubcc1pi
{

/**
 *  @brief  The subrun factory class is used to populate a ubcc1pi::Subrun from an art::Subrun 
 *          This implementation is separated from the Subrun class itself so that one can use the Subrun without any coupling to LArSoft
 */
class SubrunFactory
{
    public:

        /**
         *  @brief  The configuration required to build a ubcc1pi::Subrun
         */
        struct Config
        {
            /**
             *  @brief  If the file has truth info
             */
            fhicl::Atom<bool> HasTruthInfo
            {
                fhicl::Name("HasTruthInfo"),
                fhicl::Comment("If the input events have truth info that we should look for")
            };

            /**
             *  @brief  The POT summary label
             */
            fhicl::Atom<art::InputTag> POTSummaryLabel
            {
                fhicl::Name("POTSummaryLabel"),
                fhicl::Comment("The label for the POT summary producer")
            };
        };

        /**
         *  @brief  Populate the output subrun with information from the input subrun
         *
         *  @param  subrun the input subrun
         *  @param  config the configuration options
         *  @param  pOutputSubrun the output subrun to populate
         */
        static void PopulateSubrun(const art::SubRun &subrun, const Config &config, Subrun *pOutputSubrun);
};

} // namespace ubcc1pi

#endif
