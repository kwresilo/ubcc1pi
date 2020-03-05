/**
 *  @file  ubcc1pi/Objects/Event.h
 *
 *  @brief The header file for the event class
 */

#ifndef UBCC1PI_OBJECTS_EVENT
#define UBCC1PI_OBJECTS_EVENT

#include "ubcc1pi/Objects/EventMembers.h"

#include <iostream>
#include <stdexcept>
#include <limits>
#include <string>
#include <vector>

#include <TTree.h>

namespace ubcc1pi
{

class FileWriter;
class EventFactory;

class Event
{
    public:

        /**
         *  @brief  Print the member variables to the terminal
         */
        void Print() const;

        /**
         *  @brief  The class for member variable which will throw if not set
         */
        template <typename T>
        class Member
        {
            public:
                /**
                 *  @brief  Default constructor
                 */
                Member();

                /**
                 *  @brief  Returns if this member has been set
                 *
                 *  @return true if set, false otherwise
                 */
                bool IsSet() const;

                /**
                 *  @brief  Set the value of the member
                 *
                 *  @param  value the value to set
                 */
                void Set(const T &value);

                /**
                 *  @brief  Gets the value of the member if set and throws otherwise
                 *
                 *  @return the value
                 */
                const T &Get() const;

                /**
                 *  @brief  Overload the () operator as a shorthand Get()
                 */
                const T &operator()() const;

                /**
                 *  @brief  Reset the member value to default
                 */
                void Reset();

            private:

                // Allow the event class to access private members of this class, all other classes will have to go through the public
                // accessor functions Get() and Set()
                friend class Event;

                /**
                 *  @brief  Set the input value to default
                 *
                 *  @param  value the value to set to default
                 */
                void SetDefault(T &value);

                /**
                 *  @brief  Set the input string value to default
                 *
                 *  @param  value the value to set to default
                 */
                void SetDefault(std::string &value);

                T    m_value; ///< The value of the member
                bool m_isSet; ///< If the value has been set
        };

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

        Metadata metadata;  ///< The metadata
        Truth truth;        ///< The truth information

    private:

        friend FileWriter;
        friend EventFactory;

        /**
         *  @brief  Bind this event to an output tree
         *
         *  @param  pTree the tree with which to bind
         */
        void BindToOutputTree(TTree * pTree);

        /**
         *  @brief  Reset the member variables as if the event were new
         */
        void Reset();

        /**
         *  @brief  Fill the output vectors with the data from the particles ready for a fill
         */
        void PrepareForTreeFill();

        // Here we define private member variables for each of the particle parameters as a vector so they can be read from the root file
        UBCC1PI_MACRO_EVENT_TRUTH_PARTICLE_MEMBERS(truth_particle, "", UBCC1PI_MACRO_DECLARE_MEMBER_VECTOR)
};

// -----------------------------------------------------------------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
inline Event::Member<T>::Member()
{
    this->Reset();
}   

// -----------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
inline bool Event::Member<T>::IsSet() const
{
    return m_isSet;
}

// -----------------------------------------------------------------------------------------------------------------------------------------
        
template <typename T>
inline void Event::Member<T>::Set(const T &value)
{
    m_value = value;
    m_isSet = true;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
inline const T & Event::Member<T>::Get() const
{
    if (this->IsSet())
        return m_value;

    throw std::logic_error("You are trying to access a member value of ubcc1pi::Event that hasn't been set");
}

// -----------------------------------------------------------------------------------------------------------------------------------------

template <typename T>        
inline const T & Event::Member<T>::operator()() const
{
    return this->Get();
}
        
// -----------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
inline void Event::Member<T>::Reset()
{
    this->SetDefault(m_value);
    m_isSet = false;
}

// -----------------------------------------------------------------------------------------------------------------------------------------
        
template <typename T>
inline void Event::Member<T>::SetDefault(T &value)
{
    value = std::numeric_limits<T>::max();
}

// -----------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
inline void Event::Member<T>::SetDefault(std::string &value)
{
    value = "";
}

} // namespace ubcc1pi

#endif
