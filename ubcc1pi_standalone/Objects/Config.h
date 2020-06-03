/**
 *  @file  ubcc1pi_standalone/Objects/Config.h
 *
 *  @brief The header file for the config structure
 */

#ifndef UBCC1PI_STANDALONE_OBJECTS_CONFIG
#define UBCC1PI_STANDALONE_OBJECTS_CONFIG

#include <string>
#include <unordered_map>

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
        std::string  overlaysFileName = "/uboone/data/users/asmith/ubcc1pi/samples/may2020/samples/ubcc1piAnalysis_overlays.root";
        std::string  dirtFileName     = "/uboone/data/users/asmith/ubcc1pi/samples/may2020/samples/ubcc1piAnalysis_dirt.root";
        std::string  dataEXTFileName  = "/uboone/data/users/asmith/ubcc1pi/samples/may2020/samples/ubcc1piAnalysis_dataEXT.root";
        std::string  dataBNBFileName  = "/uboone/data/users/asmith/ubcc1pi/samples/may2020/samples/ubcc1piAnalysis_dataBNB.root";
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
        bool        useAbsPdg               = true;        ///< If we should use absolute PDG codes (this makes pi+ == pi- in the signal definition)
        bool        countProtonsInclusively = true;        ///< If we should count protons inclusively (as Xp), or exclusively as (0p, 1p, 2p, ...)
        std::string lastCutGeneric          = "noShowers"; ///< The last cut of the generic selection (remaining cuts are part of the golden selection)
        float       protonMomentumThreshold = 0.3f;        ///< The minimum proton momentum to be counted [GeV]

        /**
         *  @brief  The muonCosTheta plot limits structure
         */
        struct MuonCosTheta
        {
            float               min = -1.f;
            float               max =  1.f;
            std::vector<float>  binEdges = {-1.f, 0.5f, 0.64f, 0.75f, 0.83f, 0.88f, 0.93f, 0.96f, 0.98f, 1.f};
        };
        MuonCosTheta muonCosTheta; ///< The muonCosTheta plot limits

        /**
         *  @brief  The muonPhi plot limits structure
         */
        struct MuonPhi
        {
            float               min = -3.142f;
            float               max =  3.142f;
            std::vector<float>  binEdges = {-3.142f, -2.8f, -2.2f, -1.5f, -1.f, -0.4f, 0.1f, 0.7f, 2.f, 2.7f, 3.142f};
        };
        MuonPhi muonPhi; ///< The muonPhi plot limits

        /**
         *  @brief  The muonMomentum plot limits structure
         */
        struct MuonMomentum
        {
            float               min = 0.f;
            float               max =  100.f;
            std::vector<float>  binEdges = {0.f, 0.29f, 0.37f, 0.56f, 0.8f, 1.5f};
        };
        MuonMomentum muonMomentum; ///< The muonMomentum plot limits
       
        /**
         *  @brief  The pionCosTheta plot limits structure
         */
        struct PionCosTheta
        {
            float               min = -1.f;
            float               max =  1.f;
            std::vector<float>  binEdges = {-1.f, -0.2f, 0.23f, 0.37f, 0.6f, 0.84f, 0.92f, 1.f};
        };
        PionCosTheta pionCosTheta; ///< The pionCosTheta plot limits

        /**
         *  @brief  The pionPhi plot limits structure
         */
        struct PionPhi
        {
            float               min = -3.142f;
            float               max =  3.142f;
            std::vector<float>  binEdges = {-3.142f, -2.4f, -1.5f, -0.7f, -0.3f, 0.2f, 0.8f, 1.4f, 2.2f, 2.6f, 3.142f};
        };
        PionPhi pionPhi; ///< The pionPhi plot limits

        /**
         *  @brief  The pionMomentum plot limits structure
         */
        struct PionMomentum
        {
            float               min = 0.1f;
            float               max =  10.f;
            std::vector<float>  binEdges = {0.1f, 0.16f, 0.19f, 0.21f, 0.5f};
        };
        PionMomentum pionMomentum; ///< The pionMomentum plot limits
   
        /**
         *  @brief  The muonPionAngle plot limits structure
         */
        struct MuonPionAngle
        {
            float               min = 0.f;
            float               max =  2.65f;
            std::vector<float>  binEdges = {0.f, 0.7f, 0.9f, 1.1f, 1.4f, 1.5f, 1.7f, 2.f, 2.65f};
        };
        MuonPionAngle muonPionAngle; ///< The muonPionAngle plot limits

        /**
         *  @brief  The nProtons plot limits structure
         */
        struct NProtons
        {
            float               min = 0;
            float               max =  12;
            std::vector<float>  binEdges = {0, 1, 12};
        };
        NProtons nProtons; ///< The nProtons plot limits

    };
    Global global; ///< The global configuration options
        
    // -------------------------------------------------------------------------------------------------------------------------------------

    /**
     *  @brief  Configuration fo the CountPOT macro
     */
    struct CountPOT
    {
        bool  useOverlays = true; ///< If we should count the POT for the overlays
        bool  useDirt     = true; ///< If we should count the POT for the dirt
    };
    CountPOT countPOT; ///< The configuration options for the CountPOT macro
    
    // -------------------------------------------------------------------------------------------------------------------------------------

    /**
     *  @brief  Configuration fo the GetRunSubrunList macro
     */
    struct GetRunSubrunList
    {
        bool  useDataEXT = true; ///< If we should run on the EXT data
        bool  useDataBNB = true; ///< If we should run on the BNB data
    };
    GetRunSubrunList getRunSubrunList; ///< The configuration options for the GetRunSubrunList macro
    
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
        bool useGenericSelection = false;  ///< If we should use the generic selection (if false, we use full golden selection)
        bool goldenPionIsSignal = false;   ///< If we should only treat events containing golden pions as signal
    };
    MakeSelectedPIDTable makeSelectedPIDTable;
    
    // -------------------------------------------------------------------------------------------------------------------------------------
    
    /**
     *  @brief  Configuration for the EfficiencyPlots macro
     */
    struct EfficiencyPlots
    {
        bool drawErrors = false; ///< If we should draw errors
    };
    EfficiencyPlots efficiencyPlots; ///< The configuration options for the EfficiencyPlots macro
};

} // namespace ubcc1pi

#endif
