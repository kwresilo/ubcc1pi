#include "ubcc1pi_standalone/Objects/FileReader.h"

#include "ubcc1pi_standalone/Helpers/AnalysisHelper.h"
#include "ubcc1pi_standalone/Helpers/SelectionHelper.h"
#include "ubcc1pi_standalone/Helpers/CrossSectionHelper.h"

using namespace ubcc1pi;

int ExtractXSecs(const std::string &overlayFileName, const float overlayWeight, const std::string &dataEXTFileName, const float extWeight, const std::string &dataBNBFileName, const bool useAbsPdg = true, const bool countProtonsInclusively = true)
{
    const float bodgeFactor = 1.273f; // ATTN this factor is a normalisation added so we can compare the shape of the distributions, can't exist in the final result!
    
    // Get the selection
    auto selection = SelectionHelper::GetDefaultSelection();

    // Set up the cross-section objects
    auto xSec_muCosTheta = CrossSectionHelper::XSec<float>("xsec_muonCosTheta.root", -1.f, 1.f, useAbsPdg, countProtonsInclusively);
    xSec_muCosTheta.SetBins({-1.f, 0.f, 0.27f, 0.45f, 0.62f, 0.76f, 0.86f, 0.94f, 1.f});
    
    for (const auto fileName : {dataEXTFileName, overlayFileName, dataBNBFileName})
    {
        std::cout << "Reading input file: " << fileName << std::endl;

        const bool isBNBData = (fileName == dataBNBFileName);
        const bool isOverlay = (fileName == overlayFileName);
        const bool isEXTData = (fileName == dataEXTFileName);

        const auto normalisation = isOverlay ? overlayWeight : (isEXTData ? extWeight : 1.f);

        FileReader reader(fileName);
        auto pEvent = reader.GetBoundEventAddress();

        const auto nEvents = reader.GetNumberOfEvents();
        for (unsigned int i = 0; i < nEvents; ++i)
        {
            AnalysisHelper::PrintLoadingBar(i, nEvents);
            reader.LoadEvent(i);

            // Run the event selection and store which cuts are passed
            std::vector<std::string> cutsPassed;
            std::vector<int> assignedPdgCodes;
            const auto isSelected = selection.Execute(pEvent, cutsPassed, assignedPdgCodes);
        
            // Determine if this event passed the generic selection
            const auto lastCut = "noShowers";
            const auto passedGenericSelection = (std::find(cutsPassed.begin(), cutsPassed.end(), lastCut) != cutsPassed.end());

            // Get the variables we care about
            float muCosTheta_reco = -std::numeric_limits<float>::max();
            if (passedGenericSelection)
            {
                for (unsigned int particleIndex = 0; particleIndex < pEvent->reco.particles.size(); ++particleIndex)
                {
                    if (assignedPdgCodes.at(particleIndex) == 13)
                    {
                        const auto muon = pEvent->reco.particles.at(particleIndex);
                        muCosTheta_reco = muon.directionZ();
                    }
                }
            }

            // Add the event to the outputs
            const float weight = 1.f; // TODO assign event weights
            
            if (isOverlay && AnalysisHelper::IsTrueCC1Pi(pEvent, useAbsPdg))
            {
                float muCosTheta_true = std::numeric_limits<float>::max();

                for (const auto particle : pEvent->truth.particles)
                {
                    const auto truePdg = particle.pdgCode();
                    const auto truePdgAbs = useAbsPdg ? std::abs(truePdg) : truePdg;
    
                    if (truePdgAbs == 13)
                    {
                        const auto muDir = TVector3(particle.momentumX(), particle.momentumY(), particle.momentumZ()).Unit();
                        muCosTheta_true = muDir.Z();
                    }
                }

                xSec_muCosTheta.AddSignalEvent(pEvent, passedGenericSelection, muCosTheta_true, muCosTheta_reco, weight, normalisation);
            }
            else if (passedGenericSelection)
            {
                if (isBNBData)
                {
                    xSec_muCosTheta.AddSelectedBNBDataEvent(pEvent, muCosTheta_reco);
                }
                else
                {
                    const auto sampleType = isOverlay ? AnalysisHelper::Overlay : (isEXTData ? AnalysisHelper::DataEXT : AnalysisHelper::Dirt);
                    xSec_muCosTheta.AddSelectedBackgroundEvent(pEvent, sampleType, muCosTheta_reco, weight, normalisation);
                }
            }
        }
    }

    xSec_muCosTheta.PrintBinContents();

    return 0;
}

