/**
 *  @file  ubcc1pi_standalone/Objects/Config.h
 *
 *  @brief The header file for the config structure
 */

#ifndef UBCC1PI_STANDALONE_OBJECTS_CONFIG
#define UBCC1PI_STANDALONE_OBJECTS_CONFIG

#include <string>

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
        std::string     overlaysFileName = "/uboone/data/users/asmith/ubcc1pi/samples/may2020/ubcc1piAnalysis_overlays.root";
        std::string     dirtFileName     = "/uboone/data/users/asmith/ubcc1pi/samples/may2020/ubcc1piAnalysis_dirt.root";
        std::string     dataEXTFileName  = "/uboone/data/users/asmith/ubcc1pi/samples/may2020/ubcc1piAnalysis_dataEXT.root";
        std::string     dataBNBFileName  = "/uboone/data/users/asmith/ubcc1pi/samples/may2020/ubcc1piAnalysis_dataBNB.root";
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
        float           overlaysPOT        = 1.22447e+21;   ///< The total POT for the overlays MC
        float           dirtPOT            = 2.85049e+20;   ///< The total POT for the dirt MC
        float           dataEXTTriggers    = 62540367.0;    ///< The EXT triggers for the EXT data
        float           dataBNBTor875WCut  = 1.455e+20;     ///< The POT measured by the 875m toroid (with quality cuts)
        float           dataBNBE1DCNTWCut  = 32339256.0;    ///< The BNB spills sent by the accelerator division (with quality cuts)
    };
    Norms norms; ///< The sample normalisations

    // -------------------------------------------------------------------------------------------------------------------------------------

    /**
     *  @brief  Global configuration options structure
     */
    struct Global
    {
        bool            useAbsPdg               = true;     ///< If we should use absolute PDG codes (this makes pi+ == pi- in the signal definition)
        bool            countProtonsInclusively = true;     ///< If we should count protons inclusively (as Xp), or exclusively as (0p, 1p, 2p, ...)

        // TODO momentum thresholds
    };
    Global global; ///< The global configuration options
    
    // -------------------------------------------------------------------------------------------------------------------------------------

};

} // namespace ubcc1pi
