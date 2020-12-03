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
#include "ubcc1pi_standalone/Helpers/CrossSectionHelper.h"

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
        // September overlays files contain systematic weights
        std::string  overlaysFileName = "/uboone/data/users/asmith/ubcc1pi/samples/sep2020/samples/ubcc1piAnalysis_overlays.root"; ///< Overlays file name input
        std::string  dirtFileName     = "/uboone/data/users/asmith/ubcc1pi/samples/may2020/samples/ubcc1piAnalysis_dirt.root";     ///< Dirt file name input
        std::string  dataEXTFileName  = "/uboone/data/users/asmith/ubcc1pi/samples/may2020/samples/ubcc1piAnalysis_dataEXT.root";  ///< EXT data file name input
        std::string  dataBNBFileName  = "/uboone/data/users/asmith/ubcc1pi/samples/may2020/samples/ubcc1piAnalysis_dataBNB.root";  ///< BNB data file name input

        /**
         *  @brief  The detector variation files
         *
         *          The events in each detector variation sample (in the same run) are identical apart from the detector parameter that's
         *          been changed - this way we can limit the effects of statistical variations. As a result we also need to have a
         *          central-value (CV) sample in which nothing is changed that we can compare to. Here, only some of the variations are
         *          available for run 1. Instead we use run3b, for which the varaitions are available. In every case the variation is
         *          considered wrt to the relevant "CV" file, and the fractional difference is quantity we care about.
         */
        std::vector< std::pair<std::string, std::string> > detVarFiles = {
            // Run-1 files
            {"CVRun1",         "/uboone/data/users/asmith/ubcc1pi/samples/oct2020/samples/ubcc1piAnalysis_overlays_DetVar_CV_run1.root"},
            {"LYDown",         "/uboone/data/users/asmith/ubcc1pi/samples/oct2020/samples/ubcc1piAnalysis_overlays_DetVar_LYDown_run1.root"},
            {"LYRayleigh",     "/uboone/data/users/asmith/ubcc1pi/samples/oct2020/samples/ubcc1piAnalysis_overlays_DetVar_LYRayleigh_run1.root"},

            // Run-3b files
            {"CVRun3b",        "/uboone/data/users/asmith/ubcc1pi/samples/oct2020/samples/ubcc1piAnalysis_overlays_DetVar_CV_run3b.root"},
            {"SCE",            "/uboone/data/users/asmith/ubcc1pi/samples/oct2020/samples/ubcc1piAnalysis_overlays_DetVar_SCE_run3b.root"},
            {"Recomb2",        "/uboone/data/users/asmith/ubcc1pi/samples/oct2020/samples/ubcc1piAnalysis_overlays_DetVar_Recomb2_run3b.root"},
            {"WireModX",       "/uboone/data/users/asmith/ubcc1pi/samples/oct2020/samples/ubcc1piAnalysis_overlays_DetVar_WireModX_run3b.root"},
            {"WireModYZ",      "/uboone/data/users/asmith/ubcc1pi/samples/oct2020/samples/ubcc1piAnalysis_overlays_DetVar_WireModYZ_run3b.root"},
            {"WireModThetaXZ", "/uboone/data/users/asmith/ubcc1pi/samples/oct2020/samples/ubcc1piAnalysis_overlays_DetVar_WireModThetaXZ_run3b.root"},
            {"WireModThetaYZ", "/uboone/data/users/asmith/ubcc1pi/samples/oct2020/samples/ubcc1piAnalysis_overlays_DetVar_WireModThetaYZ_run3b.root"}
        };
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

        /**
         *  @brief  The detector variation POTs
         */
        std::unordered_map<std::string, float> detVarPOTs = {
            {"CVRun1",         1.14339e+20},
            {"LYDown",         1.05031e+20},
            {"LYRayleigh",     1.06661e+20},
            {"CVRun3b",        9.82298e+19},
            {"SCE",            1.02517e+20},
            {"Recomb2",        1.00832e+20},
            {"WireModX",       1.09739e+20},
            {"WireModYZ",      1.10877e+20},
            {"WireModThetaXZ", 1.12906e+20},
            {"WireModThetaYZ", 1.09244e+20}
        };
    };
    Norms norms; ///< The sample normalisations

    // -------------------------------------------------------------------------------------------------------------------------------------

    /**
     *  @brief  The flux structure
     *
     *  The values stored here are taken from MICROBOONE-NOTE-1031-PUB
     */
    struct Flux
    {
        /**
         *  @brief  The neutrino energy bin edges for the flux [GeV]
         */
        std::vector<float> binEdges = {
            0.00f, 0.05f, 0.10f, 0.15f, 0.20f, 0.25f, 0.30f, 0.35f, 0.40f, 0.45f, 0.50f, 0.55f, 0.60f, 0.65f, 0.70f, 0.75f, 0.80f, 0.85f,
            0.90f, 0.95f, 1.00f, 1.05f, 1.10f, 1.15f, 1.20f, 1.25f, 1.30f, 1.35f, 1.40f, 1.45f, 1.50f, 1.55f, 1.60f, 1.65f, 1.70f, 1.75f,
            1.80f, 1.85f, 1.90f, 1.95f, 2.00f, 2.05f, 2.10f, 2.15f, 2.20f, 2.25f, 2.30f, 2.35f, 2.40f, 2.45f, 2.50f, 2.55f, 2.60f, 2.65f,
            2.70f, 2.75f, 2.80f, 2.85f, 2.90f, 2.95f, 3.00f
        };

        /**
         *  @brief  The value of the flux in each energy bin [e-10 cm^-2 POT^-1]
         */
        std::vector<float> energyBins = {
            3.09e-2, 1.19e-1, 1.53e-1, 1.83e-1, 2.27e-1, 2.50e-1, 2.67e-1, 2.80e-1, 2.96e-1, 3.09e-1, 3.16e-1, 3.16e-1,
            3.12e-1, 3.09e-1, 3.06e-1, 2.99e-1, 2.87e-1, 2.75e-1, 2.63e-1, 2.49e-1, 2.36e-1, 2.22e-1, 2.06e-1, 1.92e-1,
            1.78e-1, 1.62e-1, 1.47e-1, 1.32e-1, 1.17e-1, 1.02e-1, 8.85e-2, 7.65e-2, 6.50e-2, 5.48e-2, 4.63e-2, 3.83e-2,
            3.18e-2, 2.57e-2, 2.10e-2, 1.70e-2, 1.35e-2, 1.11e-2, 9.11e-3, 7.23e-3, 6.21e-3, 5.35e-3, 4.63e-3, 4.04e-3,
            3.67e-3, 3.33e-3, 3.08e-3, 2.89e-3, 2.78e-3, 2.69e-3, 2.58e-3, 2.40e-3, 2.37e-3, 2.29e-3, 2.25e-3, 2.09e-3
        };
    };
    Flux flux; ///< The flux

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
        float       targetDensity           = 8.44191f;           ///< The number of target nuclei per unit volume - units [10^23 nucleons / cm^3]

        /**
         *  @brief  The Binning structure
         *          The supplied binEdges define the "analysis bins" in which the data will be presented. The min/max values give the limits
         *          of the phase-space that's includede in the measurment. If min/max extend beyond the bin edges, then an
         *          underflow/overflow bin will be added.
         */
        struct Binning
        {
            float               min;      ///< Minimum possible value
            float               max;      ///< Maximum possible value
            std::vector<float>  binEdges; ///< The analysis bin edges
        };

        // Here we define the binning for each kinematic variable. We're initializing each Binning object using member initialization lists.
        // I.e. we supply the member variable of the Binning struct, in the order they are deined (min, max, binEdges).

        Binning muonCosTheta {
            -1.f,                                                                              // min
             1.f,                                                                              // max
            {-1.f, -0.27f, 0.29f, 0.46f, 0.58f, 0.67f, 0.77f, 0.82f, 0.88f, 0.93f, 0.97f, 1.f} // binEdges
        }; ///< The muonCosTheta binning

        Binning muonPhi {
            -3.142f,                                                     // min
             3.142f,                                                     // max
            PlottingHelper::GenerateUniformBinEdges(15, -3.142f, 3.142f) // binEdges
        }; ///< The muonPhi binning

        Binning muonMomentum {
            0.15f,                                    // min
            std::numeric_limits<float>::max(),        // max
            {0.15f, 0.23f, 0.32f, 0.45f, 0.66f, 1.5f} // binEdges
        }; ///< The muonMomentum binning

        Binning pionCosTheta {
            -1.f,                                                // min
             1.f,                                                // max
            {-1.f, -0.47f, 0.f, 0.39f, 0.65f, 0.84f, 0.93f, 1.f} // binEdges
        }; ///< The pionCosTheta binning

        Binning pionPhi {
            -3.142f,                                                     // min
             3.142f,                                                     // max
            PlottingHelper::GenerateUniformBinEdges(10, -3.142f, 3.142f) // binEdges
        }; ///< The pionPhi binning

        Binning pionMomentum {
            0.f,                               // min
            std::numeric_limits<float>::max(), // max
            {0.1f, 0.16f, 0.19f, 0.22f, 0.6f}  // binEdges
        }; ///< The pionMomentum binning

        Binning muonPionAngle {
            0.f,                                                   // min
            2.65f,                                                 // max
            {0.f, 0.49f, 0.93f, 1.26f, 1.57f, 1.88f, 2.21f, 2.65f} // binEdges
        }; ///< The muonPionAngle binning

        Binning nProtons {
            0,                                           // min
            std::numeric_limits<float>::max(),           // max
            {0, 1, 2, std::numeric_limits<float>::max()} // binEdges
        }; ///< The nProtons binning

    };
    Global global; ///< The global configuration options

    // -------------------------------------------------------------------------------------------------------------------------------------

    /**
     *  @brief  Configuration for the CountPOT macro
     */
    struct CountPOT
    {
        bool  useOverlays           = true; ///< If we should count the POT for the overlays
        bool  useDirt               = true; ///< If we should count the POT for the dirt
        bool  useDetectorVariations = true; ///< If we should count the POT for the detector variations
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

    // -------------------------------------------------------------------------------------------------------------------------------------

    /**
     *  @brief  Configuration for the ExtractXSecs macro
     */
    struct ExtractXSecs
    {
        /**
        *  @brief  If we should scale the GENIE weights so not to double count the genieTuneEventWeight
        *          By default we apply (splineEventWeight * genieTuneEventWeight) as the "nominal weight" to all events. Then we apply the
        *          multisim universe weights on top of the nominal weight. For the GENIE systematic parameters, the universe weights already
        *          include the genieTuneEventWeight, and so require special treatment to avoid double counting this weight. If this option
        *          is set to true, then we scale the GENIE universe weights down by genieTuneEventWeight to put them on the same footing as
        *          all other parameters. If this option is false then the GENIE weights recieve no special treatment.
        */
        bool scaleXSecWeights = true;

        unsigned int nBootstrapUniverses = 1000u; ///< The number of bootrap universes to generate for the MC stat uncertainty

        CrossSectionHelper::SystDimensionsMap fluxDimensions = {
            {"hadronProduction",          1000u},
            {"expskin_FluxUnisim",        1000u},
            {"horncurrent_FluxUnisim",    1000u},
            {"nucleoninexsec_FluxUnisim", 1000u},
            {"nucleonqexsec_FluxUnisim",  1000u},
            {"nucleontotxsec_FluxUnisim", 1000u},
            {"pioninexsec_FluxUnisim",    1000u},
            {"pionqexsec_FluxUnisim",     1000u},
            {"piontotxsec_FluxUnisim",    1000u}
        }; ///< A mapping from the flux parameter names to the number of universes

        CrossSectionHelper::SystDimensionsMap xsecDimensions = {
            {"All_Genie",             100u},
            {"AxFFCCQEshape_Genie",   2u},
            {"DecayAngMEC_Genie",     2u},
            {"MaNCRES_Genie",         2u},
            {"Theta_Delta2Npi_Genie", 2u},
            {"VecFFCCQEshape_Genie",  2u}
        }; ///< A mapping from the cross-section parameter names to the number of universes

        CrossSectionHelper::SystUnisimDimensionsMap detVarDimensions = {
            {"LYDown",         "CVRun1"},
            {"LYRayleigh",     "CVRun1"},
            {"SCE",            "CVRun3b"},
            {"Recomb2",        "CVRun3b"},
            {"WireModX",       "CVRun3b"},
            {"WireModYZ",      "CVRun3b"},
            {"WireModThetaXZ", "CVRun3b"},
            {"WireModThetaYZ", "CVRun3b"},
        }; ///< A mapping from the detector variation sample identifiers, to the identifiers for their relevant central-value sample

        /**
        *  @brief  A mapping from a user-defined parameter name, to corresponding set of mutually exclusive parameter names and the number of universes
        *
        *          A set of parameters is mutually exclusive if for any given universe at most one parameter has a weight that's not equal
        *          to one. In this case, the set of parameters can be "combined" into one parameter (with the user-defined name). The value
        *          of the combined parameter is taken from whichever parameter in the set has been varied(or equivalently just the product
        *          of the weights in the input parameter set).
        *
        *          For example, the flux hadron production channel can be any one of kminus, kplus, kzero, piminus, piplus, but never more
        *          than one channel at once. As a result, a given event will have a (possibly) non-unit weight for the relevant channel, and
        *          a weight of one for all other channels. In this sense, there is a "correlation" between the weights for each channel.
        *          Instead of applying each channel as a separate systematic parameter, we apply their product as a combined parameter which
        *          accounts for their "correlations"
        */
        CrossSectionHelper::SystMutuallyExclusiveDimensionsMap mutuallyExclusiveDimensions = {
            {
                "hadronProduction",
                {
                    {
                        "kminus_PrimaryHadronNormalization",
                        "kplus_PrimaryHadronFeynmanScaling",
                        "kzero_PrimaryHadronSanfordWang",
                        "piminus_PrimaryHadronSWCentralSplineVariation",
                        "piplus_PrimaryHadronSWCentralSplineVariation"
                    }, 1000u
                }
            }
        };

    };
    ExtractXSecs extractXSecs; ///< The configuration options for the ExtractXSecs macro
};

} // namespace ubcc1pi

#endif
