/**
 *  @file  ubcc1pi_standalone/Helpers/PlottingHelper.h
 *
 *  @brief The header file for the plotting helper class
 */

#ifndef UBCC1PI_STANDALONE_HELPERS_PLOTTING_HELPER
#define UBCC1PI_STANDALONE_HELPERS_PLOTTING_HELPER

#include <string>
#include <memory>
#include <unordered_map>

#include <TH1F.h>
#include <TCanvas.h>

#include "ubcc1pi_standalone/Interface/Event.h"

namespace ubcc1pi
{

/**
 *  @brief  The plotting helper class
 */
class PlottingHelper
{
    public:
        /**
         *  @brief  The particle style enumeration
         */
        enum ParticleStyle
        {
            External,
            ExternalPoints,
            Muon,
            MuonPoints,
            Proton,
            ProtonPoints,
            GoldenPion,
            GoldenPionPoints,
            NonGoldenPion,
            NonGoldenPionPoints,
            PiMinus,
            PiMinusPoints,
            Electron,
            ElectronPoints,
            Photon,
            PhotonPoints,
            Other,
            OtherPoints,
            BNBData
        };

        /**
         *  @brief  The vector of all particle styles
         */
        static const std::vector<ParticleStyle> AllParticleStyles;

        /**
         *  @brief  If we should draw the input style with points (or a line)
         *
         *  @param  style the input style
         *
         *  @return boolean, true if points should be used
         */
        static bool ShouldUsePoints(const ParticleStyle &style);

        /**
         *  @brief  The plot class that manages histograms of particle characteristics
         */
        class ParticlePlot
        {
            public:
                /**
                 *  @brief  Constructor
                 *
                 *  @param  xLabel the x-label of the histogram
                 *  @param  nBins the number of bins
                 *  @param  min the minimum value 
                 *  @param  max the maximum value
                 *  @param  drawErrors whether to draw the error bands
                 */
                ParticlePlot(const std::string &xLabel, unsigned int nBins, float min, float max, bool drawErrors = true);

                /**
                 *  @brief  Fill the histogram with the given value and particle style
                 *
                 *  @param  value the value to fill
                 *  @param  particleStyle the style of the hisogram to fill
                 *  @param  weight the weight
                 */
                void Fill(const float value, const ParticleStyle &particleStyle, const float weight = 1.f);

                /**
                 *  @brief  Draw and save the plot
                 *
                 *  @param  fileName the output file name (don't include an extension)
                 */
                void SaveAs(const std::string &fileName);
                
                /**
                 *  @brief  Draw and save the plot as a stacked histogram
                 *
                 *  @param  fileName the output file name (don't include an extension)
                 */
                void SaveAsStacked(const std::string &fileName);

            private:
                /**
                 *  @brief  Get clones of the histograms so they can be safely manipulated 
                 *
                 *  @param  particleToHistCloneMap the output map of cloned histograms
                 */
                void GetHistogramClones(std::unordered_map<ParticleStyle, TH1F* > &particleToHistCloneMap);

                /**
                 *  @brief  Normalise the histograms by their number of entries
                 *
                 *  @param  particleToHistCloneMap the histograms to normalise
                 */
                void ScaleHistograms(std::unordered_map<ParticleStyle, TH1F*> &particleToHistCloneMap) const;

                /**
                 *  @brief  Set the Y-range of the histograms so all fit on the canvas
                 *
                 *  @param  particleToHistCloneMap the histograms to modify
                 */
                void SetHistogramYRanges(std::unordered_map<ParticleStyle, TH1F*> &particleToHistCloneMap) const;

                std::string  m_xLabel;     ///< The x-label of the histogram
                unsigned int m_nBins;      ///< The number of bins
                float        m_min;        ///< The minimum value
                float        m_max;        ///< The maximum value
                unsigned int m_id;         ///< The ID of this particle plot
                unsigned int m_cloneCount; ///< A count of the number of clones to avoid name collisions
                bool         m_drawErrors; ///< Whether to draw the error bands

                std::unordered_map<ParticleStyle, std::shared_ptr<TH1F> > m_particleToHistMap; ///< The mapping from particle style to hist

                static unsigned int m_lastId; ///< The last ParticlePlot ID that was set
        };

        /**
         *  @brief  Get the particle style type for the input reco particle using it's truth match information
         *
         *  @param  particle the input reco particle
         *  @param  truthParticles the input list of all truth particles
         *  @param  usePoints if we should use datapoints instead of a line
         *
         *  @return the particle style
         */
        static ParticleStyle GetParticleStyle(const Event::Reco::Particle &particle, const std::vector<Event::Truth::Particle> &truthParticles, const bool usePoints = false);

        /**
         *  @brief  Set the line style for a given particle type
         *
         *  @tparam T the ROOT object class
         *  @param  pObject the ROOT object
         *  @param  particleStyle the style
         */
        template<typename T>
        static void SetLineStyle(T *pObject, const ParticleStyle particleStyle);
        
        /**
         *  @brief  Set the line style for a given particle type
         *
         *  @tparam T the ROOT object class
         *  @param  pObject the ROOT object
         *  @param  particleStyle the style
         */
        template<typename T>
        static void SetLineStyle(std::shared_ptr<T> &pObject, const ParticleStyle particleStyle);

        /**
         *  @brief  Get a canvas to draw on
         *
         *  @param  width the canvas width (px)
         *  @param  height the canvas height (px)
         *
         *  @return the canvas
         */
        static std::shared_ptr<TCanvas> GetCanvas(const unsigned int width = 960, const unsigned int height = 540);

        /**
         *  @brief  Save the canvas
         *
         *  @param  canvas the canvas to save
         *  @param  fileName the output file name (no extension)
         */
        static void SaveCanvas(std::shared_ptr<TCanvas> &canvas, const std::string &fileName);

    private:

        static unsigned int m_lastCanvasId;  ///< The last canvas ID
};

// -----------------------------------------------------------------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------------------------------------------------------------

const std::vector<PlottingHelper::ParticleStyle> PlottingHelper::AllParticleStyles = {
    External,
    ExternalPoints,
    Muon,
    MuonPoints,
    Proton,
    ProtonPoints,
    GoldenPion,
    GoldenPionPoints,
    NonGoldenPion,
    NonGoldenPionPoints,
    PiMinus,
    PiMinusPoints,
    Electron,
    ElectronPoints,
    Photon,
    PhotonPoints,
    Other,
    OtherPoints,
    BNBData
};

// -----------------------------------------------------------------------------------------------------------------------------------------

unsigned int PlottingHelper::ParticlePlot::m_lastId = 0;

// -----------------------------------------------------------------------------------------------------------------------------------------
        
unsigned int PlottingHelper::m_lastCanvasId = 0;

// -----------------------------------------------------------------------------------------------------------------------------------------
        
template<typename T>
inline void PlottingHelper::SetLineStyle(std::shared_ptr<T> &pObject, const ParticleStyle particleStyle)
{
    PlottingHelper::SetLineStyle(pObject.get(), particleStyle);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

template<typename T>
inline void PlottingHelper::SetLineStyle(T *pObject, const ParticleStyle particleStyle)
{
    auto col = static_cast<int>(kBlack);

    switch (particleStyle)
    {
        case Muon:
        case MuonPoints:
            col = kAzure - 2;
            break;

        case Proton:
        case ProtonPoints:
            col = kOrange - 3;
            break;

        case GoldenPion:
        case GoldenPionPoints:
            col = kGreen + 1;
            break;

        case NonGoldenPion:
        case NonGoldenPionPoints:
            col = kMagenta + 1;
            break;
        
        case PiMinus:
        case PiMinusPoints:
            col = kRed - 4;
            break;
        
        case Electron:
        case ElectronPoints:
            col = kCyan + 1;
            break;
        
        case Photon:
        case PhotonPoints:
            col = kYellow + 1;
            break;
        
        case External:
        case ExternalPoints:
            col = kGray + 3;
            break;

        case Other:
        case OtherPoints:
            col = kGray;
            break;
        
        case BNBData:
            col = kBlack;
            break;

        default: break;
    }
    
    pObject->SetLineWidth(2);
    pObject->SetLineColor(col);

    if (PlottingHelper::ShouldUsePoints(particleStyle))
    {
        pObject->SetMarkerColor(col);
    }
}

// -----------------------------------------------------------------------------------------------------------------------------------------
        
bool PlottingHelper::ShouldUsePoints(const ParticleStyle &style)
{
    if (style == ExternalPoints ||
        style == MuonPoints ||
        style == ProtonPoints ||
        style == GoldenPionPoints ||
        style == NonGoldenPionPoints ||
        style == PiMinusPoints ||
        style == ElectronPoints ||
        style == PhotonPoints ||
        style == OtherPoints)
        return true;

    return false;
}

} // namespace ubcc1pi

#endif
