#include "ubcc1pi_standalone/Macros/Macros.h"

#include "ubcc1pi_standalone/Objects/FileReader.h"
#include "ubcc1pi_standalone/Helpers/AnalysisHelper.h"
#include "ubcc1pi_standalone/Helpers/NormalisationHelper.h"
#include "ubcc1pi_standalone/Helpers/SelectionHelper.h"

using namespace ubcc1pi;

namespace ubcc1pi_macros
{

void PrintRunSubruns(const Config &config)
{
    //
    // Setup the input files
    // 
    std::vector< std::tuple<AnalysisHelper::SampleType, std::string, float> > inputData;
    
    inputData.emplace_back(AnalysisHelper::Overlay, config.files.overlaysFileName, NormalisationHelper::GetOverlaysNormalisation(config)); 
    inputData.emplace_back(AnalysisHelper::Dirt,    config.files.dirtFileName,     NormalisationHelper::GetDirtNormalisation(config)); 
    inputData.emplace_back(AnalysisHelper::DataEXT, config.files.dataEXTFileName,  NormalisationHelper::GetDataEXTNormalisation(config)); 
    inputData.emplace_back(AnalysisHelper::DataBNB, config.files.dataBNBFileName,  1.f);

    //
    // Get the selection
    //
    const std::string chosenCut = "2NonProtons"; // TODO make configurable
    auto selection = SelectionHelper::GetDefaultSelection();
    const auto allCuts = selection.GetCuts();

    if (std::find(allCuts.begin(), allCuts.end(), chosenCut) == allCuts.end())
        throw std::invalid_argument("PrintRunSubruns - chosen cut \"" + chosenCut + "\" isn't known to the selection");



    //
    // Fill the cross-section objects
    //
    for (const auto [sampleType, fileName, normalisation] : inputData)
    {
        std::cout << "Reading input file: " << fileName << std::endl;

        FileReader reader(fileName);
        auto pEvent = reader.GetBoundEventAddress();
        
        const auto nEvents = reader.GetNumberOfEvents();
        for (unsigned int i = 0; i < nEvents; ++i)
        {
            reader.LoadEvent(i);

            const auto weight = AnalysisHelper::GetNominalEventWeight(pEvent);

            // Run the event selection and store which cuts are passed
            std::vector<std::string> cutsPassed;
            std::vector<int> assignedPdgCodes;
            selection.Execute(pEvent, cutsPassed, assignedPdgCodes);
            const auto passedChosenCut = (std::find(cutsPassed.begin(), cutsPassed.end(), chosenCut) != cutsPassed.end());

            // Insist we pass the selection up to this point
            if (!passedChosenCut)
                continue;

            // Get the muon reco data
            bool foundRecoMuon = false;
            unsigned int muonIndex = std::numeric_limits<unsigned int>::max();
            const auto recoParticles = pEvent->reco.particles;
            for (unsigned int index = 0; index < recoParticles.size(); ++index)
            {
                if (recoParticles.at(index).isCCInclusiveMuonCandidate())
                {
                    if (foundRecoMuon)
                        throw std::logic_error("Found multiple muon candidates");

                    muonIndex = index;
                    foundRecoMuon = true;
                }
            }
            
            if (!foundRecoMuon)
                continue;

            const auto muon = recoParticles.at(muonIndex);
            const auto muonDir = TVector3(muon.directionX(), muon.directionY(), muon.directionZ()).Unit();
            //const float muonCosTheta = muonDir.Z();
            const float muonPhi = std::atan2(muonDir.Y(), muonDir.X());
            //const float muonMomentum = AnalysisHelper::GetMuonMomentum(muon);


            // Now insist the muon has the properties we are looking for
            // TODO make configurable
            const auto piBy2 = std::acos(0);
            if (std::abs(muonPhi + piBy2) > 0.25f)
                continue;

            const auto metadata = pEvent->metadata;
            std::cout << muonPhi << " | " <<  metadata.run() << ", " << metadata.subRun() << ", " << metadata.event() << std::endl;

        }
    }
}

} // namespace ubcc1pi_macros
