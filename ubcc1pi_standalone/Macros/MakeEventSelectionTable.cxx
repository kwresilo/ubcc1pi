#include "ubcc1pi_standalone/Objects/FileReader.h"

#include "ubcc1pi_standalone/Helpers/AnalysisHelper.h"
#include "ubcc1pi_standalone/Helpers/BDTHelper.h"

using namespace ubcc1pi;

int MakeEventSelectionTable(const std::string &overlayFileName, const float overlayWeight, const std::string &dataEXTFileName, const float extWeight, const std::string &dataBNBFileName)
{
    // Set up the BDTs
    const auto featureNames = BDTHelper::ParticleBDTFeatureNames;
    BDTHelper::BDT goldenPionBDT("goldenPion", featureNames); 
    BDTHelper::BDT protonBDT("proton", featureNames); 

    AnalysisHelper::EventCounter counter; 
    for (const auto fileName : {dataEXTFileName, dataBNBFileName, overlayFileName})
    {
        std::cout << "Reading input file: " << fileName << std::endl;

        const bool isBNBData = (fileName == dataBNBFileName);
        const bool isOverlay = (fileName == overlayFileName);
        const bool isEXTData = (fileName == dataEXTFileName);

        const auto sampleType = isBNBData ? AnalysisHelper::DataBNB : (isEXTData ? AnalysisHelper::DataEXT : (isOverlay ? AnalysisHelper::Overlay : AnalysisHelper::Dirt));

        const float bodgeFactor = 1.273f; // ATTN this factor is a normalisation added so we can compare the shape of the distributions, can't exist in the final result!
        float weight = 1.f;

        if (isOverlay) weight = overlayWeight * bodgeFactor;
        if (isEXTData) weight = extWeight * bodgeFactor;

        FileReader reader(fileName);
        auto pEvent = reader.GetBoundEventAddress();

        const auto nEvents = reader.GetNumberOfEvents();
        for (unsigned int i = 0; i < nEvents; ++i)
        {
            AnalysisHelper::PrintLoadingBar(i, nEvents);

            reader.LoadEvent(i);
            
            // For overlays count the events before any selection
            if (isOverlay)
            {
                counter.CountEvent("all", sampleType, pEvent, weight);
            }


            // Only use events passing the CC inclusive selection
            if (!pEvent->reco.passesCCInclusive())
                continue;

            counter.CountEvent("passesCCInclusive", sampleType, pEvent, weight);


            // Insist there are at least 2 reco particles
            const auto recoParticles = pEvent->reco.particles;
            if (recoParticles.size() < 2)
                continue;
            
            counter.CountEvent("min2Particles", sampleType, pEvent, weight);

            // Insist that the BDT features are calculable for all particles
            std::vector< std::vector<float> > particleBDTFeatures;
            bool areFeaturesAvailable = true;
            for (const auto &particlce : recoParticles)
            {
                std::vector<float> features;
                if (!BDTHelper::GetBDTFeatures(particle, featureNames, features))
                {
                    areFeaturesAvailable = false;
                    break;
                }
            
                particleBDTFeatures.push_back(features);
            }

            if (!areFeaturesAvailable)
                continue;

            counter.CountEvent("allBDTFeaturesAvailable", sampleType, pEvent, weight);
        }
    }

    std::cout << "BREAKDOWN" << std::endl;
    counter.PrintBreakdownWithWeights(8);

    return 0;
}
