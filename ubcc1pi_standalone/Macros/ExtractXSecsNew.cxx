/**
 *  @file  ubcc1pi_standalone/Macros/ExtractXSecsNew.cxx
 *
 *  @brief The implementation file of the ExtractXSecs macro
 */

#include "ubcc1pi_standalone/Macros/Macros.h"
#include "ubcc1pi_standalone/Objects/FileReader.h"

#include "ubcc1pi_standalone/Helpers/AnalysisHelper.h"
#include "ubcc1pi_standalone/Helpers/NormalisationHelper.h"
#include "ubcc1pi_standalone/Helpers/CrossSectionHelperNew.h"

using namespace ubcc1pi;

namespace ubcc1pi_macros
{

void ExtractXSecsNew(const Config &config)
{
    //
    // Setup the input files
    //
    std::vector< std::tuple<AnalysisHelper::SampleType, std::string, float> > inputData;
    inputData.emplace_back(AnalysisHelper::Overlay, config.files.overlaysFileName, NormalisationHelper::GetOverlaysNormalisation(config));


    // TODO make these configurable
    // Setup an object that holds the details of the systematic parameters to apply
    CrossSectionHelperNew::CrossSection::SystParams systParams;

    // The number of "bootstrap" universes, used to find the MC statistical uncertainty
    systParams.nBootstrapUniverses = 1000u;

    // The flux parameters
    systParams.fluxDimensions = {
        {"hadronProduction",          1000u},
        {"expskin_FluxUnisim",        1000u},
        {"horncurrent_FluxUnisim",    1000u},
        {"nucleoninexsec_FluxUnisim", 1000u},
        {"nucleonqexsec_FluxUnisim",  1000u},
        {"nucleontotxsec_FluxUnisim", 1000u},
        {"pioninexsec_FluxUnisim",    1000u},
        {"pionqexsec_FluxUnisim",     1000u},
        {"piontotxsec_FluxUnisim",    1000u}
    };

    // The cross-section parameters
    systParams.xsecDimensions = {
        {"All_Genie",             100u},
        {"AxFFCCQEshape_Genie",   2u},
        {"DecayAngMEC_Genie",     2u},
        {"MaNCRES_Genie",         2u},
        {"Theta_Delta2Npi_Genie", 2u},
        {"VecFFCCQEshape_Genie",  2u}
    };

    // The detector variation parameters
    systParams.detVarDimensions = {
        {"LYDown",         "CVRun1"},
        {"LYRayleigh",     "CVRun1"},
        {"SCE",            "CVRun3b"},
        {"Recomb2",        "CVRun3b"},
        {"WireModX",       "CVRun3b"},
        {"WireModYZ",      "CVRun3b"},
        {"WireModThetaXZ", "CVRun3b"},
        {"WireModThetaYZ", "CVRun3b"},
    };

    // Defined the dimensions of any mutually exclusive parameters
    CrossSectionHelperNew::SystMutuallyExclusiveDimensionsMap mutuallyExclusiveDimensions = {
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

    // Setup the flux reweightor
    CrossSectionHelperNew::FluxReweightor fluxReweightor(config.flux.binEdges, config.flux.energyBins, fluxDimensions);
    std::cout << "Integrated flux: " << fluxReweightor.GetIntegratedNominalFlux() << std::endl;

    // Loop over the events
    for (const auto &[sampleType, fileName, normalisation] : inputData)
    {
        std::cout << "Reading input file: " << fileName << std::endl;
        const auto isOverlay = (sampleType == AnalysisHelper::Overlay);

        // Open the input file for reading and enable the branches with systematic event weights (if required)
        FileReader reader(fileName);
        if (isOverlay) reader.EnableSystematicBranches();
        auto pEvent = reader.GetBoundEventAddress();

        const auto nEvents = reader.GetNumberOfEvents();
        for (unsigned int i = 0; i < nEvents; ++i)
        {
            AnalysisHelper::PrintLoadingBar(i, nEvents);
            reader.LoadEvent(i);

            // Get the nominal event weight, scaled by the sample normalisation
            const auto weight = AnalysisHelper::GetNominalEventWeight(pEvent) * normalisation;

            // Get the flux weights
            const auto fluxWeights = CrossSectionHelperNew::GetWeightsMap(pEvent->truth, fluxDimensions);

            // Add the true neutrino energy to the flux reweightor
            fluxReweightor.AddEvent(pEvent->truth.nuEnergy(), weight, fluxWeights);
        }
    }
}

} // namespace ubcc1pi_macros
