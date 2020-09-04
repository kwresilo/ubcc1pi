/**
 *  @file  ubcc1pi_standalone/Objects/Config.h
 *
 *  @brief The header file for the config structure
 */

#ifndef UBCC1PI_STANDALONE_OBJECTS_CONFIG
#define UBCC1PI_STANDALONE_OBJECTS_CONFIG

#include <string>
#include <unordered_map>

#include "ubcc1pi_standalone/Helpers/PlottingHelper.h"

namespace ubcc1pi
{

/**
 *  @brief  The config struture
 */
struct Config
{
    /**
     *  @brief  Input files structure
     *
     * Specify the location of special ubcc1pi analysis root files that the standalone analysis code uses as input
     * You can produce them from art root files by running the ubcc1pi::AnalysisFileWriter module
     *
     */
    struct Files
    {
        std::string  overlaysFileName = "/uboone/data/users/asmith/ubcc1pi/samples/may2020/samples/ubcc1piAnalysis_overlays.root"; ///< Overlays file name input
        std::string  dirtFileName     = "/uboone/data/users/asmith/ubcc1pi/samples/may2020/samples/ubcc1piAnalysis_dirt.root";     ///< Dirt file name input
        std::string  dataEXTFileName  = "/uboone/data/users/asmith/ubcc1pi/samples/may2020/samples/ubcc1piAnalysis_dataEXT.root";  ///< EXT data file name input
        std::string  dataBNBFileName  = "/uboone/data/users/asmith/ubcc1pi/samples/may2020/samples/ubcc1piAnalysis_dataBNB.root";  ///< BNB data file name input
    };
    Files files; ///< The input files

    // -------------------------------------------------------------------------------------------------------------------------------------

    /**
     *  @brief  The sample normalisations structure
     *
     * These values can be found by:
     *      - Running the ubcc1pi::CountPOT macro (for overlays and dirt)
     *      - Running the ubcc1pi::GetRunSubrunList macro (for data), and then running Zarko's tool:
     *          /uboone/app/users/zarko/getDataInfo.py -v2 --run-subrun-list runSubrunList.txt
     */
    struct Norms
    {
        float  overlaysPOT        = 1.22447e+21;   ///< The total POT for the overlays MC
        float  dirtPOT            = 2.85049e+20;   ///< The total POT for the dirt MC
        float  dataEXTTriggers    = 62540367.0;    ///< The EXT triggers for the EXT data
        float  dataBNBTor875WCut  = 1.455e+20;     ///< The POT measured by the 875m toroid (with quality cuts)
        float  dataBNBE1DCNTWCut  = 32339256.0;    ///< The BNB spills sent by the accelerator division (with quality cuts)
    };
    Norms norms; ///< The sample normalisations

    // -------------------------------------------------------------------------------------------------------------------------------------

    /**
     *  @brief  Global configuration options structure
     */
    struct Global
    {
        bool        useAbsPdg               = true;               ///< If we should use absolute PDG codes (this makes pi+ == pi- in the signal definition)
        bool        countProtonsInclusively = true;               ///< If we should count protons inclusively (as Xp), or exclusively as (0p, 1p, 2p, ...)
        std::string lastCutGeneric          = "startNearVertex";  ///< The last cut of the generic selection (remaining cuts are part of the golden selection)
        float       protonMomentumThreshold = 0.3f;               ///< The minimum proton momentum to be counted [GeV]

        /**
         *  @brief  The muonCosTheta plot limits structure
         */
        struct MuonCosTheta
        {
            float               min = -1.f;                                                                     ///< Minimum possible value
            float               max =  1.f;                                                                     ///< Maximum possible value
//            std::vector<float>  binEdges = {-1.f, 0.5f, 0.64f, 0.75f, 0.83f, 0.88f, 0.93f, 0.96f, 0.98f, 1.f};  ///< The bin edges
            std::vector<float>  binEdges = {-1.f, -0.27f, 0.29f, 0.46f, 0.58f, 0.67f, 0.77f, 0.82f, 0.88f, 0.93f, 0.97f, 1.f};  ///< The bin edges
        };
        MuonCosTheta muonCosTheta; ///< The muonCosTheta plot limits

        /**
         *  @brief  The muonPhi plot limits structure
         */
        struct MuonPhi
        {
            float               min = -3.142f;                                                           ///< Minimum possible value
            float               max =  3.142f;                                                           ///< Maximum possible value
            //std::vector<float>  binEdges = PlottingHelper::GenerateUniformBinEdges(10, -3.142f, 3.142f); ///< The bin edges
            std::vector<float>  binEdges = PlottingHelper::GenerateUniformBinEdges(15, -3.142f, 3.142f); ///< The bin edges
        };
        MuonPhi muonPhi; ///< The muonPhi plot limits

        /**
         *  @brief  The muonMomentum plot limits structure
         */
        struct MuonMomentum
        {
            /*
            float               min = 0.090f;                                         ///< Minimum possible value
            float               max = 100.f;                                          ///< Maximum possible value
            std::vector<float>  binEdges = {0.090f, 0.29f, 0.37f, 0.56f, 0.8f, 1.5f}; ///< The bin edges
            */
            float               min = 0.15f;                                          ///< Minimum possible value
            float               max = 100.f;                                          ///< Maximum possible value
//            std::vector<float>  binEdges = {0.15f, 0.29f, 0.37f, 0.56f, 0.8f, 1.5f};  ///< The bin edges
            std::vector<float>  binEdges = {0.15f, 0.23f, 0.32f, 0.45f, 0.66f, 1.5f};  ///< The bin edges
        };
        MuonMomentum muonMomentum; ///< The muonMomentum plot limits
       
        /**
         *  @brief  The pionCosTheta plot limits structure
         */
        struct PionCosTheta
        {
            float               min = -1.f;                                                       ///< Minimum possible value
            float               max =  1.f;                                                       ///< Maximum possible value
//            std::vector<float>  binEdges = {-1.f, -0.2f, 0.23f, 0.37f, 0.6f, 0.84f, 0.92f, 1.f};  ///< The bin edges
            std::vector<float>  binEdges = {-1.f, -0.47f, 0.f, 0.39f, 0.65f, 0.84f, 0.93f, 1.f};  ///< The bin edges
        };
        PionCosTheta pionCosTheta; ///< The pionCosTheta plot limits

        /**
         *  @brief  The pionPhi plot limits structure
         */
        struct PionPhi
        {
            float               min = -3.142f;                                                           ///< Minimum possible value
            float               max =  3.142f;                                                           ///< Maximum possible value
            std::vector<float>  binEdges = PlottingHelper::GenerateUniformBinEdges(10, -3.142f, 3.142f); ///< The bin edges
        };
        PionPhi pionPhi; ///< The pionPhi plot limits

        /**
         *  @brief  The pionMomentum plot limits structure
         */
        struct PionMomentum
        {
            /*
            float               min = 0.114f;                                    ///< Minimum possible value
            float               max = 10.f;                                      ///< Maximum possible value
            std::vector<float>  binEdges = {0.114f, 0.16f, 0.19f, 0.21f, 0.5f};  ///< The bin edges
            */
            float               min = 0.f;                                            ///< Minimum possible value
            float               max = 10.f;                                           ///< Maximum possible value
            //std::vector<float>  binEdges = {0.f, 0.114f, 0.16f, 0.19f, 0.21f, 0.5f};  ///< The bin edges
            //std::vector<float>  binEdges = {0.1f, 0.15f, 0.2f, 0.25f, 0.6f};  ///< The bin edges
            std::vector<float>  binEdges = {0.1f, 0.16f, 0.19f, 0.22f, 0.6f};  ///< The bin edges
        };
        PionMomentum pionMomentum; ///< The pionMomentum plot limits
   
        /**
         *  @brief  The muonPionAngle plot limits structure
         */
        struct MuonPionAngle
        {
            float               min = 0.f;                                                         ///< Minimum possible value
            float               max = 2.65f;                                                      ///< Maximum possible value
            //std::vector<float>  binEdges = {0.f, 0.7f, 0.9f, 1.1f, 1.4f, 1.5f, 1.7f, 2.f, 2.65f};  ///< The bin edges
            std::vector<float>  binEdges = {0.f, 0.49f, 0.93f, 1.26f, 1.57f, 1.88f, 2.21f, 2.65f};  ///< The bin edges
        };
        MuonPionAngle muonPionAngle; ///< The muonPionAngle plot limits

        /**
         *  @brief  The nProtons plot limits structure
         */
        struct NProtons
        {
            float               min = 0;                  ///< Minimum possible value
            float               max = 12;                 ///< Maximum possible value
            std::vector<float>  binEdges = {0, 1, 2, 3};  ///< The bin edges
        };
        NProtons nProtons; ///< The nProtons plot limits

    };
    Global global; ///< The global configuration options
        
    // -------------------------------------------------------------------------------------------------------------------------------------

    /**
     *  @brief  Configuration for the CountPOT macro
     */
    struct CountPOT
    {
        bool  useOverlays = true; ///< If we should count the POT for the overlays
        bool  useDirt     = true; ///< If we should count the POT for the dirt
    };
    CountPOT countPOT; ///< The configuration options for the CountPOT macro
    
    // -------------------------------------------------------------------------------------------------------------------------------------

    /**
     *  @brief  Configuration for the GetRunSubrunList macro
     */
    struct GetRunSubrunList
    {
        bool  useDataEXT = true; ///< If we should run on the EXT data
        bool  useDataBNB = true; ///< If we should run on the BNB data
    };
    GetRunSubrunList getRunSubrunList; ///< The configuration options for the GetRunSubrunList macro
    
    // -------------------------------------------------------------------------------------------------------------------------------------
    
    /**
     *  @brief  Configuration for the MultiPlanePIDDemo macro
     */
    struct MultiPlanePIDDemo
    {
        float sin2AngleThreshold = 0.175; ///< The squared sin angular threshold between a particle and a wire in the YZ plane to use dEdx information
    };
    MultiPlanePIDDemo multiPlanePIDDemo; ///< The configuration options for the MultiPlanePIDDemo macro
    
    // -------------------------------------------------------------------------------------------------------------------------------------

    /**
     *  @brief  Configuration for the PlotInputVariables macro
     */
    struct PlotInputVariables
    {
        bool plotBDTResponses = true; ///< If we should plot the responses of the trained BDTs, set to true if you haven't already trained the BDTs
    };
    PlotInputVariables plotInputVariables; ///< The configuration options for the PlotInputVariables macro
    
    // -------------------------------------------------------------------------------------------------------------------------------------

    /**
     *  @brief  Configuration for the NMinusOneBDTStudy macro
     */
    struct NMinusOneBDTStudy
    {
        bool shouldTrainBDTs = true;                                          ///< If we should run the BDT training (if false then look for trained BDTs)
        PlottingHelper::PlotStyle signalType = PlottingHelper::GoldenPion;    ///< The type of particle considered signal by the BDT in question
        std::vector< std::string > featureNames = {
            "logBragg_pToMIP",
            "logBragg_piToMIP",
            "truncMeandEdx",
            "protonForward",
            "muonForward",
            "nDescendents",
            "nSpacePointsNearEnd",
            "wiggliness",
            "trackScore"
        };                                                                    ///< The features to consider by the BDT
        unsigned int nSamplePoints = 300u;                                    ///< The number of sampling points to use when finding the ROC curves
    };
    NMinusOneBDTStudy nMinusOneBDTStudy; ///< The configuration options for the NMinusOneBDTStudy macro
    
    // -------------------------------------------------------------------------------------------------------------------------------------

    /**
     *  @brief  Configuration for TrainBDTs macro
     */
    struct TrainBDTs
    {
        float trainingFraction     = 0.5f;  ///< Fraction of the sample on which we should train the BDTs
        bool  onlyGoodTruthMatches = false; ///< If we should only train on reco particles with a completeness >50%
        bool  weightByCompleteness = true;  ///< If we should weight the training examples by the reco-truth match completeness
        bool  shouldOptimize       = false; ///< If we should optimize the BDT parameters (this can take a while)
        bool  shouldMakePlots      = true;  ///< If we should make plots of the BDT responses
    };
    TrainBDTs trainBDTs; ///< The configuration options for the TrainBDTs macro

    // -------------------------------------------------------------------------------------------------------------------------------------

    /**
     *  @brief  Configuration for the MakeEventSelectionTable macro
     */
    struct MakeEventSelectionTable
    {
        bool         shouldOptimize  = false; ///< If we should optimize the cuts
        unsigned int nScanPoints     = 20u;   ///< The number of scan points to use while optimizing
        float        processFraction = 0.2f;  ///< The fraction of events to process while optimizing
    };
    MakeEventSelectionTable makeEventSelectionTable; ///< The configuration options for the MakeEventSelectionTable macro
    
    // -------------------------------------------------------------------------------------------------------------------------------------

    /**
     *  @brief  Configuration for the MakeSelectedPIDTable macro
     */
    struct MakeSelectedPIDTable
    {
        bool  useGenericSelection = false;  ///< If we should use the generic selection (if false, we use full golden selection)
        bool  goldenPionIsSignal = false;   ///< If we should only treat events containing golden pions as signal
        bool  onlyLowMomentumPions = false; ///< If we should only treat events with low momentum pions as "signal" - to check the performance in a restricted region of phase-space
        float pionMomentumThreshold = 0.1f; ///< The threshold pion momentum below which we consider "signal" if onlyLowMomentumPions == true
    };
    MakeSelectedPIDTable makeSelectedPIDTable; ///< The configuration options for the MakeSelectedPIDTable macro
    
    // -------------------------------------------------------------------------------------------------------------------------------------
    
    /**
     *  @brief  Configuration for the EfficiencyPlots macro
     */
    struct EfficiencyPlots
    {
        bool drawErrors = true; ///< If we should draw errors
    };
    EfficiencyPlots efficiencyPlots; ///< The configuration options for the EfficiencyPlots macro
    
    // -------------------------------------------------------------------------------------------------------------------------------------

    /**
     *  @brief  Configuration for the MakeBinningPlots macro 
     */
    struct MakeBinningPlots
    {
        bool useFineBinEdges = true; ///< If we break up the analsis bins into finer sub-bins
    };
    MakeBinningPlots makeBinningPlots; ///< The configuration options for the MakeBinningPlots macro
};

} // namespace ubcc1pi

#endif
