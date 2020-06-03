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
#include "ubcc1pi_standalone/Helpers/AnalysisHelper.h"

namespace ubcc1pi
{

/**
 *  @brief  The plotting helper class
 */
class PlottingHelper
{
    public:
        /**
         *  @brief  The plot style enumeration
         */
        enum PlotStyle
        {
            // Particle types
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

            // Event types
            NonFiducial,
            NumuCC0Pi,
            NumuCC1PiChargedGolden,
            NumuCC1PiChargedNonGolden,
            NumuCC1PiZero,
            NumuCCOther,
            Nue,
            NC,
            
            // Common types
            External,
            ExternalPoints,
            Dirt,
            DirtPoints,
            BNBData,
            Other,
            OtherPoints
        };
        
        /**
         *  @brief  The vector of all plot styles
         */
        static const std::vector<PlotStyle> AllPlotStyles;
        
        /**
         *  @brief  If we should draw the input style with points (or a line)
         *
         *  @param  style the input style
         *
         *  @return boolean, true if points should be used
         */
        static bool ShouldUsePoints(const PlotStyle &style);

        /**
         *  @brief  The plot class that manages multiple related histograms 
         */
        class MultiPlot
        {
            public:
                /**
                 *  @brief  Constructor
                 *
                 *  @param  xLabel the x-label of the histogram
                 *  @param  yLabel the y-label of the histogram
                 *  @param  nBins the number of bins
                 *  @param  min the minimum value 
                 *  @param  max the maximum value
                 *  @param  drawErrors whether to draw the error bands
                 */
                MultiPlot(const std::string &xLabel, const std::string &yLabel, unsigned int nBins, float min, float max, bool drawErrors = true);
                
                /**
                 *  @brief  Constructor
                 *
                 *  @param  xLabel the x-label of the histogram
                 *  @param  yLabel the y-label of the histogram
                 *  @param  binEdges the bin edges (for variable binning)
                 *  @param  drawErrors whether to draw the error bands
                 */
                MultiPlot(const std::string &xLabel, const std::string &yLabel, const std::vector<float> &binEdges, bool drawErrors = true);

                /**
                 *  @brief  Fill the histogram with the given value and plot style
                 *
                 *  @param  value the value to fill
                 *  @param  plotStyle the style of the hisogram to fill
                 *  @param  weight the weight
                 */
                void Fill(const float value, const PlotStyle &plotStyle, const float weight = 1.f);

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
                 *  @param  plotToHistCloneMap the output map of cloned histograms
                 */
                void GetHistogramClones(std::unordered_map<PlotStyle, TH1F* > &plotToHistCloneMap);

                /**
                 *  @brief  Normalise the histograms by their number of entries
                 *
                 *  @param  plotToHistCloneMap the histograms to normalise
                 */
                void ScaleHistograms(std::unordered_map<PlotStyle, TH1F*> &plotToHistCloneMap) const;

                /**
                 *  @brief  Set the Y-range of the histograms so all fit on the canvas
                 *
                 *  @param  plotToHistCloneMap the histograms to modify
                 */
                void SetHistogramYRanges(std::unordered_map<PlotStyle, TH1F*> &plotToHistCloneMap) const;

                std::string        m_xLabel;     ///< The x-label of the histogram
                unsigned int       m_nBins;      ///< The number of bins
                float              m_min;        ///< The minimum value
                float              m_max;        ///< The maximum value
                std::vector<float> m_binEdges;   ///< The edges of the bins (for variable binning)
                unsigned int       m_id;         ///< The ID of this plot plot
                unsigned int       m_cloneCount; ///< A count of the number of clones to avoid name collisions
                bool               m_drawErrors; ///< Whether to draw the error bands

                std::unordered_map<PlotStyle, std::shared_ptr<TH1F> > m_plotToHistMap; ///< The mapping from plot style to hist

                static unsigned int m_lastId; ///< The last plotPlot ID that was set
        };

        /**
         *  @brief  Efficiency plot class
         */
        class EfficiencyPlot
        {
            public:
                /**
                 *  @brief  Constructor
                 *
                 *  @param  xLabel the x-label of the histogram
                 *  @param  nBins the number of bins
                 *  @param  min the minimum value 
                 *  @param  max the maximum value
                 *  @param  cuts the names of all possible cuts in order
                 *  @param  drawErrors whether to draw the error bands
                 */
                EfficiencyPlot(const std::string &xLabel, unsigned int nBins, float min, float max, const std::vector<string> &cuts, bool drawErrors = true);
       
                /**
                 *  @brief  Add an event to the plot for the given cut
                 *
                 *  @param  value the value at which to add the event
                 *  @param  cut the cut 
                 *  @param  passedCut if the cut was passed
                 */
                void AddEvent(const float value, const std::string &cut, const bool passedCut);

                /**
                 *  @brief  Draw and save the plot
                 *
                 *  @param  fileName the file name without an extension
                 */
                void SaveAs(const std::string &fileName);

            private:


                /**
                 *  @brief  Mapping from a cut name to a pair of histograms, first being the numerator, second being the denominator
                 */
                typedef std::unordered_map<std::string, std::pair< std::shared_ptr<TH1F>, std::shared_ptr<TH1F> > > CutToPlotsMap;

                std::string              m_xLabel;        ///< The x label
                unsigned int             m_nBins;         ///< The number of bins
                float                    m_min;           ///< The minimum histogram value
                float                    m_max;           ///< The maximum histogram value
                std::vector<std::string> m_cuts;          ///< The cuts
                bool                     m_drawErrors;    ///< If we should draw error bands

                CutToPlotsMap            m_cutToPlotsMap; ///< The mapping from cut name to the numerator and denominator histograms
                unsigned int             m_id;            ///< The ID of the plot
                static unsigned int      m_lastId;        ///< The last plot ID that was set - avoids ROOT name collisions
        };

        /**
         *  @brief  Get the particle style type for the input reco particle using it's truth match information
         *
         *  @param  particle the input reco particle
         *  @param  sampleType the input sample type
         *  @param  truthParticles the input list of all truth particles
         *  @param  usePoints if we should use datapoints instead of a line
         *  @param  useAbsPdg if we should use the absolute PDG codes
         *
         *  @return the particle style
         */
        static PlotStyle GetPlotStyle(const Event::Reco::Particle &particle, const AnalysisHelper::SampleType &sampleType, const std::vector<Event::Truth::Particle> &truthParticles, const bool usePoints = false, const bool useAbsPdg = false);

        /**
         *  @brief  Get a color for a given style
         *
         *  @param  plotStyle the input style
         *
         *  @return the color
         */
        static int GetColor(const PlotStyle plotStyle);

        /**
         *  @brief  Get the complete list of unique colors used
         *
         *  @return the colors
         */
        static std::vector<int> GetColorVector();
        
        /**
         *  @brief  Get the plot style of an event
         *
         *  @param  sampleType the input sample type
         *  @param  pEvent the input event
         *  @param  useAbsPdg if we should use the absolute PDG codes
         *
         *  @return the plot style
         */
        static PlotStyle GetPlotStyle(const AnalysisHelper::SampleType &sampleType, const std::shared_ptr<Event> &pEvent, const bool useAbsPdg);
        
        /**
         *  @brief  Set the line style for a given plot type
         *
         *  @tparam T the ROOT object class
         *  @param  pObject the ROOT object
         *  @param  plotStyle the plot style
         */
        template<typename T>
        static void SetLineStyle(T *pObject, const PlotStyle plotStyle);

        /**
         *  @brief  Set the line style for a given color
         *
         *  @tparam T the ROOT object class
         *  @param  pObject the ROOT object
         *  @param  col the color
         */
        template<typename T>
        static void SetLineStyle(T *pObject, const int col);
        
        /**
         *  @brief  Set the line style for a given plot type
         *
         *  @tparam T the ROOT object class
         *  @param  pObject the ROOT object
         *  @param  plotStyle the style
         */
        template<typename T>
        static void SetLineStyle(std::shared_ptr<T> &pObject, const PlotStyle plotStyle);
        
        /**
         *  @brief  Set the line style for a given plot type
         *
         *  @tparam T the ROOT object class
         *  @param  pObject the ROOT object
         *  @param  col the color
         */
        template<typename T>
        static void SetLineStyle(std::shared_ptr<T> &pObject, const int col);

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

const std::vector<PlottingHelper::PlotStyle> PlottingHelper::AllPlotStyles = {
    External,
    ExternalPoints,
    Dirt,
    DirtPoints,
    NonFiducial,
    Other,
    OtherPoints,
    PiMinus,
    PiMinusPoints,
    Electron,
    ElectronPoints,
    Photon,
    PhotonPoints,
    Proton,
    ProtonPoints,
    Muon,
    MuonPoints,
    NonGoldenPion,
    NonGoldenPionPoints,
    GoldenPion,
    GoldenPionPoints,
    Nue,
    NC,
    NumuCCOther,
    NumuCC0Pi,
    NumuCC1PiZero,
    NumuCC1PiChargedGolden,
    NumuCC1PiChargedNonGolden,
    BNBData
};

// -----------------------------------------------------------------------------------------------------------------------------------------

unsigned int PlottingHelper::MultiPlot::m_lastId = 0;

// -----------------------------------------------------------------------------------------------------------------------------------------

unsigned int PlottingHelper::EfficiencyPlot::m_lastId = 0;

// -----------------------------------------------------------------------------------------------------------------------------------------
        
unsigned int PlottingHelper::m_lastCanvasId = 0;

// -----------------------------------------------------------------------------------------------------------------------------------------
        
template<typename T>
inline void PlottingHelper::SetLineStyle(std::shared_ptr<T> &pObject, const PlotStyle plotStyle)
{
    PlottingHelper::SetLineStyle(pObject.get(), plotStyle);
}

// -----------------------------------------------------------------------------------------------------------------------------------------
        
template<typename T>
inline void PlottingHelper::SetLineStyle(std::shared_ptr<T> &pObject, const int col)
{
    PlottingHelper::SetLineStyle(pObject.get(), col);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

int PlottingHelper::GetColor(const PlotStyle plotStyle)
{
    auto col = static_cast<int>(kBlack);
    
    switch (plotStyle)
    {
        case Muon:
        case MuonPoints:
        case NumuCC0Pi:
            col = kAzure - 2;
            break;

        case Proton:
        case ProtonPoints:
        case NumuCC1PiZero:
            col = kOrange - 3;
            break;

        case GoldenPion:
        case GoldenPionPoints:
        case NumuCC1PiChargedGolden:
            col = kGreen + 1;
            break;

        case NonGoldenPion:
        case NonGoldenPionPoints:
        case NumuCC1PiChargedNonGolden:
            col = kMagenta + 1;
            break;
        
        case PiMinus:
        case PiMinusPoints:
        case NonFiducial:
            col = kRed - 4;
            break;
        
        case Electron:
        case ElectronPoints:
        case NumuCCOther:
            col = kCyan + 1;
            break;
        
        case Photon:
        case PhotonPoints:
        case Nue:
            col = kYellow + 1;
            break;

        case NC:
            col = kGreen + 3;
            break;
        
        case External:
        case ExternalPoints:
            col = kGray + 3;
            break;
        
        case Dirt:
        case DirtPoints:
            col = kMagenta - 5;
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

    return col;
}
        
// -----------------------------------------------------------------------------------------------------------------------------------------

std::vector<int> PlottingHelper::GetColorVector()
{
    std::vector<int> colors;

    for (const auto &style : PlottingHelper::AllPlotStyles)
    {
        const auto col = PlottingHelper::GetColor(style);

        if (std::find(colors.begin(), colors.end(), col) == colors.end())
            colors.push_back(col);
    }

    return colors;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

template<typename T>
inline void PlottingHelper::SetLineStyle(T *pObject, const int col)
{
    pObject->SetLineWidth(2);
    pObject->SetLineColor(col);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

template<typename T>
inline void PlottingHelper::SetLineStyle(T *pObject, const PlotStyle plotStyle)
{
    const auto col = PlottingHelper::GetColor(plotStyle);

    pObject->SetLineWidth(2);
    pObject->SetLineColor(col);

    if (PlottingHelper::ShouldUsePoints(plotStyle))
    {
        pObject->SetMarkerColor(col);
    }
}

// -----------------------------------------------------------------------------------------------------------------------------------------
        
bool PlottingHelper::ShouldUsePoints(const PlotStyle &style)
{
    if (style == ExternalPoints ||
        style == DirtPoints ||
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

// -----------------------------------------------------------------------------------------------------------------------------------------

} // namespace ubcc1pi

#endif
