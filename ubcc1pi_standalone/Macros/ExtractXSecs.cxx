#include "ubcc1pi_standalone/Objects/FileReader.h"

#include "ubcc1pi_standalone/Helpers/AnalysisHelper.h"
#include "ubcc1pi_standalone/Helpers/SelectionHelper.h"
#include "ubcc1pi_standalone/Helpers/CrossSectionHelper.h"

using namespace ubcc1pi;

int ExtractXSecs(const std::string &overlayFileName, const float overlayWeight, const std::string &dataEXTFileName, const float extWeight, const std::string &dataBNBFileName, const bool useAbsPdg = true, const bool countProtonsInclusively = true)
{
    const float bodgeFactor = 28678.f / 22534.f; // ATTN this factor is a normalisation added so we can compare the shape of the distributions, can't exist in the final result!
    
    // Get the selection
    auto selection = SelectionHelper::GetDefaultSelection();

    // Set up the cross-section objects
    auto xSec_muCosTheta = CrossSectionHelper::XSec("xsec_muonCosTheta.root", -1.f, 1.f, useAbsPdg, countProtonsInclusively);
    auto xSec_piCosTheta = CrossSectionHelper::XSec("xsec_pionCosTheta.root", -1.f, 1.f, useAbsPdg, countProtonsInclusively);
    auto xSec_piMomentum = CrossSectionHelper::XSec("xsec_pionMomentum.root", 0.f, 10.f, useAbsPdg, countProtonsInclusively);
    
    for (const auto fileName : {dataEXTFileName, overlayFileName, dataBNBFileName})
    {
        std::cout << "Reading input file: " << fileName << std::endl;

        const bool isBNBData = (fileName == dataBNBFileName);
        const bool isOverlay = (fileName == overlayFileName);
        const bool isEXTData = (fileName == dataEXTFileName);

        const auto normalisation = isOverlay ? (overlayWeight * bodgeFactor): (isEXTData ? (extWeight * bodgeFactor): 1.f);

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
            const auto passedGoldenSelection = selection.Execute(pEvent, cutsPassed, assignedPdgCodes);
        
            // Determine if this event passed the generic selection
            const auto lastCut = "noShowers";
            const auto passedGenericSelection = (std::find(cutsPassed.begin(), cutsPassed.end(), lastCut) != cutsPassed.end());

            // Get the variables we care about
            float muCosTheta_reco = -std::numeric_limits<float>::max();
            float piCosTheta_reco = -std::numeric_limits<float>::max();
            float piMomentum_reco = -std::numeric_limits<float>::max();

            if (passedGenericSelection)
            {
                for (unsigned int particleIndex = 0; particleIndex < pEvent->reco.particles.size(); ++particleIndex)
                {
                    if (assignedPdgCodes.at(particleIndex) == 13)
                    {
                        const auto muon = pEvent->reco.particles.at(particleIndex);
                        muCosTheta_reco = muon.directionZ();
                    }
                    
                    if (assignedPdgCodes.at(particleIndex) == 211)
                    {
                        const auto pion = pEvent->reco.particles.at(particleIndex);
                        piCosTheta_reco = pion.directionZ();
                        piMomentum_reco = AnalysisHelper::GetPionMomentumFromRange(pion.range());
                    }
                }
            }

            // Add the event to the outputs
            const float weight = 1.f; // TODO assign event weights
            
            if (isOverlay && AnalysisHelper::IsTrueCC1Pi(pEvent, useAbsPdg))
            {
                float muCosTheta_true = std::numeric_limits<float>::max();
                float piCosTheta_true = std::numeric_limits<float>::max();
                float piMomentum_true = std::numeric_limits<float>::max();

                for (const auto particle : pEvent->truth.particles)
                {
                    const auto truePdg = particle.pdgCode();
                    const auto truePdgAbs = useAbsPdg ? std::abs(truePdg) : truePdg;
    
                    if (truePdgAbs == 13)
                    {
                        const auto muDir = TVector3(particle.momentumX(), particle.momentumY(), particle.momentumZ()).Unit();
                        muCosTheta_true = muDir.Z();
                    }
                    
                    if (truePdgAbs == 211)
                    {
                        const auto piDir = TVector3(particle.momentumX(), particle.momentumY(), particle.momentumZ()).Unit();
                        piCosTheta_true = piDir.Z();
                        piMomentum_true = particle.momentum();
                    }
                }

                xSec_muCosTheta.AddSignalEvent(pEvent, passedGenericSelection, muCosTheta_true, muCosTheta_reco, weight, normalisation);
                xSec_piCosTheta.AddSignalEvent(pEvent, passedGenericSelection, piCosTheta_true, piCosTheta_reco, weight, normalisation);

                xSec_piMomentum.AddSignalEvent(pEvent, passedGoldenSelection, piMomentum_true, piMomentum_reco, weight, normalisation);
            }
            else if (passedGenericSelection)
            {
                if (isBNBData)
                {
                    xSec_muCosTheta.AddSelectedBNBDataEvent(pEvent, muCosTheta_reco);
                    xSec_piCosTheta.AddSelectedBNBDataEvent(pEvent, piCosTheta_reco);

                    if (passedGoldenSelection)
                    {
                        xSec_piMomentum.AddSelectedBNBDataEvent(pEvent, piMomentum_reco);
                    }
                }
                else
                {
                    const auto sampleType = isOverlay ? AnalysisHelper::Overlay : (isEXTData ? AnalysisHelper::DataEXT : AnalysisHelper::Dirt);
                    xSec_muCosTheta.AddSelectedBackgroundEvent(pEvent, sampleType, muCosTheta_reco, weight, normalisation);
                    xSec_piCosTheta.AddSelectedBackgroundEvent(pEvent, sampleType, piCosTheta_reco, weight, normalisation);
                    
                    if (passedGoldenSelection)
                    {
                        xSec_piMomentum.AddSelectedBackgroundEvent(pEvent, sampleType, piMomentum_reco, weight, normalisation);
                    }
                }
            }
        }
    }
    
    // Set the binning
    const auto piMomentumCutoff = AnalysisHelper::GetPionMomentumFromRange(5.f); // TODO make this configurable

    /*
    xSec_muCosTheta.SetBins({-1.f, 0.f, 0.27f, 0.45f, 0.62f, 0.76f, 0.86f, 0.94f, 1.f});
    xSec_piCosTheta.SetBins({-1.f, 0.f, 0.27f, 0.45f, 0.62f, 0.76f, 0.86f, 0.94f, 1.f});
    xSec_piMomentum.SetBins({piMomentumCutoff, 0.135f, 0.17f, 0.2f, 0.25f, 0.5f});
    */
    xSec_muCosTheta.SetBinsAuto(-1.f, 1.f, 0.2f, 0.5f);
    xSec_piCosTheta.SetBinsAuto(-1.f, 1.f, 0.2f, 0.5f);
    xSec_piMomentum.SetBinsAuto(piMomentumCutoff, 0.5f, 0.2f, 0.5f);

    // Print the results
    std::cout << "Muon cos(theta)" << std::endl;
    xSec_muCosTheta.PrintBinContents();
    xSec_muCosTheta.MakePlots("muCosTheta");
    
    std::cout << "Pion cos(theta)" << std::endl;
    xSec_piCosTheta.PrintBinContents();
    xSec_piCosTheta.MakePlots("piCosTheta");
    
    std::cout << "Pion momentum" << std::endl;
    xSec_piMomentum.PrintBinContents();
    xSec_piMomentum.MakePlots("piMomentum");

    return 0;
}

