#include "ubcc1pi_standalone/Objects/FileReader.h"

#include "ubcc1pi_standalone/Helpers/AnalysisHelper.h"
#include "ubcc1pi_standalone/Helpers/BDTHelper.h"

using namespace ubcc1pi;

int MakeEventSelectionTable(const std::string &overlayFileName, const float overlayWeight, const std::string &dataEXTFileName, const float extWeight, const std::string &dataBNBFileName, const float eventFraction = 1.f)
{
    // Set up the BDTs
    const auto featureNames = BDTHelper::ParticleBDTFeatureNames;
    BDTHelper::BDT goldenPionBDT("goldenPion", featureNames); 
    BDTHelper::BDT protonBDT("proton", featureNames); 
    BDTHelper::BDT muonBDT("muon", featureNames); 

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

        const auto nEvents = static_cast<unsigned int>(std::floor(static_cast<float>(reader.GetNumberOfEvents()) * std::max(0.f, std::min(1.f, eventFraction))));
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


            // Get the particles that have a track fit, and insist that there are at least 3
            const auto recoParticles = pEvent->reco.particles;

            std::vector<Event::Reco::Particle> fittedParticles;
            for (const auto &particle : recoParticles)
            {
                if (AnalysisHelper::HasTrackFit(particle))
                    fittedParticles.push_back(particle);
            }
            
            if (fittedParticles.size() < 3)
                continue;
            
            counter.CountEvent("min3Tracks", sampleType, pEvent, weight);
            
            
            // Insist that no more than one particle is uncontained
            std::vector<Event::Reco::Particle> uncontainedParticles, containedParticles;
            for (const auto &particle : fittedParticles)
            {
                auto &particleVector = AnalysisHelper::IsContained(particle) ? containedParticles : uncontainedParticles;
                particleVector.push_back(particle);
            }

            if (uncontainedParticles.size() > 1)
                continue;
            
            counter.CountEvent("max1Uncontained", sampleType, pEvent, weight);
           

            // Insist we have at least 1 longish contained track
            std::vector<Event::Reco::Particle> containedLongParticles, containedShortParticles;
            for (const auto &particle : containedParticles)
            {
                const auto isLong = (particle.length() > 5.f);
                auto &particleVector = isLong ? containedLongParticles : containedShortParticles;
                particleVector.push_back(particle);
            }

            if (containedLongParticles.empty())
                continue;

            counter.CountEvent("min1LongUncontained", sampleType, pEvent, weight);


            // Identify the clear protons
            std::vector<Event::Reco::Particle> muons = uncontainedParticles;
            std::vector<Event::Reco::Particle> protons = containedShortParticles;
            std::vector<Event::Reco::Particle> others;
            for (const auto &particle : containedLongParticles)
            {
                bool isProton = true;

                std::vector<float> features;
                if (BDTHelper::GetBDTFeatures(particle, featureNames, features))
                {
                    const auto protonBDTResponse = protonBDT.GetResponse(features);
                    isProton = (protonBDTResponse > 0.2f);
                }

                auto &particleVector = isProton ? protons : others;
                particleVector.push_back(particle);
            }

            // Insist there's the right number of non-protons
            if ((others.size() + muons.size()) < 2)
                continue;

            if (others.empty())
                continue;
            
            counter.CountEvent("nonProtonMultiplicity", sampleType, pEvent, weight);
            
            
            // Identify the golden pion and muon candidate
            auto goldenPionIndex = std::numeric_limits<unsigned int>::max();
            auto goldenPionBDTResponse = -std::numeric_limits<float>::max();
            for (unsigned int i = 0; i < others.size(); ++i)
            {
                const auto &particle = others.at(i);

                std::vector<float> features;
                if (!BDTHelper::GetBDTFeatures(particle, featureNames, features))
                    throw std::logic_error("MakeEventSeletionTable - can't run BDT on non-proton candidate");

                const auto response = goldenPionBDT.GetResponse(features);
                if (response < goldenPionBDTResponse)
                    continue;

                goldenPionBDTResponse = response;
                goldenPionIndex = i;
            }

            const auto goldenPion = others.at(goldenPionIndex);

            // Add the remaining others to the muon candidates
            for (unsigned int i = 0; i < others.size(); ++i)
            {
                if (i == goldenPionIndex)
                    continue;

                muons.push_back(others.at(i));
            }
            others.clear();

            if (muons.empty())
                throw std::logic_error("MakeEventSeletionTable - no muon candidates!");

            // Select the longest remaining particle as the muon
            auto muonIndex = std::numeric_limits<unsigned int>::max();
            auto muonLength = -std::numeric_limits<float>::max();
            for (unsigned int i = 0; i < muons.size(); ++i)
            {
                const auto &particle = muons.at(i);

                if (particle.length() < muonLength)
                    continue;

                muonLength = particle.length();
                muonIndex = i;
            }

            const auto muon = muons.at(muonIndex);

            std::vector<Event::Reco::Particle> ambiguousParticles;
            for (unsigned int i = 0; i < muons.size(); ++i)
            {
                if (i == muonIndex)
                    continue;

                ambiguousParticles.push_back(muons.at(i));
            }

            // We should now have identified the golden pion candidate, the muon candidate, the clear protons and the remaining
            if (ambiguousParticles.size() + protons.size() + 2 != fittedParticles.size())
                throw std::logic_error("MakeEventSeletionTable - something doesn't add up!");

            // Insist the muon and pion aren't back to back
            const auto muonDir = TVector3(muon.directionX(), muon.directionY(), muon.directionZ()).Unit();
            const auto pionDir = TVector3(goldenPion.directionX(), goldenPion.directionY(), goldenPion.directionZ()).Unit();

            const auto openingAngle = muonDir.Angle(pionDir);
            const auto backToBack = (openingAngle > 2.6f);

            if (backToBack)
                continue;

            counter.CountEvent("notBackToBack", sampleType, pEvent, weight);


            // Insist that the golden pion candidate is likely
            if (goldenPionBDTResponse < -0.1f)
                continue;

            counter.CountEvent("likelyGoldenPion", sampleType, pEvent, weight);

        }
    }

    std::cout << "Breakdown summary" << std::endl;
    counter.PrintBreakdownSummary();
    
    std::cout << "Breakdown details" << std::endl;
    counter.PrintBreakdownDetails();

    return 0;
}
