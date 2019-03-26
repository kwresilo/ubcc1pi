/**
 *  @file  ubcc1pi/Helpers/DebugHelper.h
 *
 *  @brief The header file for the debug helper class
 */

#ifndef UBCC1PI_HELPERS_DEBUG_HELPER
#define UBCC1PI_HELPERS_DEBUG_HELPER

#include "art/Framework/Principal/Event.h"
#include "nusimdata/SimulationBase/MCParticle.h"

#include "ubcc1pi/Helpers/TruthHelper.h"
#include "ubcc1pi/Helpers/BacktrackHelper.h"

#include <iomanip>

namespace ubcc1pi
{

/**
 *  @brief  The debug helper class
 */
class DebugHelper
{
    public:
        /**
         *  @brief  Print an MCParticle
         *
         *  @param  particle the particle to print
         *  @param  depth the number of indents to use, for hierarchical information
         */
        static void Print(const art::Ptr<simb::MCParticle> &particle, const unsigned int depth = 0);

        // TODO doxygen comments
        static void Print(const MCParticleVector &particles, const unsigned int depth = 0);
        static void PrintSummary(const MCParticleVector &particles, const unsigned int depth = 0);
        static void Print(const TruthHelper::Interaction &interaction, const unsigned int depth = 0);
        static void Print(const art::Ptr<recob::PFParticle> &particle, const unsigned int depth = 0);
        static void Print(const BacktrackHelper::BacktrackerData &data, const unsigned int depth = 0);
        static std::string GetInteractionString(const TruthHelper::Interaction &interaction, const bool brief = false);
    private:

        static void PrintHeader(const std::string &title, const unsigned int depth = 0);
        static void PrintProperty(const std::string &key, const unsigned int depth = 0);

        template <typename T>
        static void PrintProperty(const std::string &key, const T &value, const unsigned int depth = 0);
        
        static unsigned int m_indent;     ///< The number of spaces to an indent
        static unsigned int m_lineWidth;  ///< The number of characters in a horizontal rule
        static unsigned int m_keyWidth;   ///< The number of characters to a key
        static char         m_lineChar;   ///< The character to use for a horizontal rule
};

// Static member variables
unsigned int DebugHelper::m_indent = 4;
unsigned int DebugHelper::m_lineWidth = 90;
unsigned int DebugHelper::m_keyWidth = 16;
char DebugHelper::m_lineChar = '-';


// -----------------------------------------------------------------------------------------------------------------------------------------

inline void DebugHelper::PrintHeader(const std::string &title, const unsigned int depth)
{
    std::cout << std::string(DebugHelper::m_lineWidth, DebugHelper::m_lineChar) << '\r' << std::string(depth * DebugHelper::m_indent, ' ') << title << " " << std::endl;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

inline void DebugHelper::PrintProperty(const std::string &key, const unsigned int depth)
{
    DebugHelper::PrintProperty(key, "", depth);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
inline void DebugHelper::PrintProperty(const std::string &key, const T &value, const unsigned int depth)
{
    std::cout << std::string(depth * DebugHelper::m_indent, ' ') << "- " << std::setw(DebugHelper::m_keyWidth) << std::left << key << " : " << value << std::endl;
}

} // namespace ubcc1pi

#endif
