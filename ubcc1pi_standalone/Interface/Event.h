/**
 *  @file  ubcc1pi_standalone/Interface/Event.h
 *
 *  @brief The header file for the event class
 */

#ifndef UBCC1PI_STANDALONE_INTERFACE_EVENT
#define UBCC1PI_STANDALONE_INTERFACE_EVENT

#include "ubcc1pi_standalone/Interface/EventMembers.h"
#include "ubcc1pi_standalone/Interface/Member.h"

#include <iostream>
#include <stdexcept>
#include <limits>
#include <string>
#include <vector>
#include <memory>

#include <TTree.h>

namespace ubcc1pi
{

class FileWriter;
class FileReader;
class EventFactory;

class Event
{
    public:

        /**
         *  @brief  Constructor
         */
        Event();

        /**
         *  @brief  Print the member variables to the terminal
         */
        void Print() const;

        /**
         *  @brief  The metadata information structure
         */
        struct Metadata
        {
            UBCC1PI_MACRO_EVENT_METADATA_MEMBERS("", "", UBCC1PI_MACRO_DECLARE_MEMBER)
        };

        /**
         *  @brief  The truth information structure
         */
        struct Truth
        {
            UBCC1PI_MACRO_EVENT_TRUTH_MEMBERS("", "", UBCC1PI_MACRO_DECLARE_MEMBER)

            /**
             *  @brief  The truth particle information structure
             */
            struct Particle
            {
                UBCC1PI_MACRO_EVENT_TRUTH_PARTICLE_MEMBERS("", "", UBCC1PI_MACRO_DECLARE_MEMBER)
            };

            std::vector<Particle> particles; ///< The truth particles
        };
        
        /**
         *  @brief  The reco information structure
         */
        struct Reco
        {
            UBCC1PI_MACRO_EVENT_RECO_MEMBERS("", "", UBCC1PI_MACRO_DECLARE_MEMBER)

            /**
             *  @brief  The reco particle information structure
             */
            struct Particle
            {
                UBCC1PI_MACRO_EVENT_RECO_PARTICLE_MEMBERS("", "", UBCC1PI_MACRO_DECLARE_MEMBER)
            };

            std::vector<Particle> particles; ///< The reco particles
        };

        Metadata metadata;  ///< The metadata
        Truth truth;        ///< The truth information
        Reco  reco;         ///< The reco information

    private:

        friend FileWriter;
        friend FileReader;
        friend EventFactory;

        /**
         *  @brief  Bind this event to an output tree
         *
         *  @param  pTree the tree with which to bind
         */
        void BindToOutputTree(TTree * pTree);
        
        /**
         *  @brief  Bind this event to an input tree
         *
         *  @param  pTree the tree with which to bind
         */
        void BindToInputTree(TTree * pTree);

        /**
         *  @brief  Reset the member variables as if the event were new
         */
        void Reset();

        /**
         *  @brief  Fill the output vectors with the data from the particles ready for a fill
         */
        void PrepareForTreeFill();

        /**
         *  @brief  Fill the input particles with the info from the vectors in the tree
         */
        void PrepareAfterTreeRead();

        // Here we define private member variables for each of the particle parameters as a vector so they can be read from the root file
        UBCC1PI_MACRO_EVENT_TRUTH_PARTICLE_MEMBERS(truth_particle, "", UBCC1PI_MACRO_DECLARE_MEMBER_VECTOR)
        UBCC1PI_MACRO_EVENT_RECO_PARTICLE_MEMBERS(reco_particle, "", UBCC1PI_MACRO_DECLARE_MEMBER_VECTOR)
};

} // namespace ubcc1pi

#endif
