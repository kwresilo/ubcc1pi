/**
 *  @file  ubcc1pi_standalone/Macros/PlotProtonMomentumByInteraction.cxx
 *
 *  @brief The implementation file of the PlotProtonMomentumByInteraction macro
 */

#include "ubcc1pi_standalone/Macros/Macros.h"

#include "ubcc1pi_standalone/Objects/FileReader.h"

#include "ubcc1pi_standalone/Helpers/PlottingHelper.h"
#include "ubcc1pi_standalone/Helpers/AnalysisHelper.h"
#include "ubcc1pi_standalone/Helpers/SelectionHelper.h"
#include "ubcc1pi_standalone/Helpers/NormalisationHelper.h"


using namespace ubcc1pi;

namespace ubcc1pi_macros
{

void PlotProtonMomentumByInteraction(const Config &config)
{
    //
    // Setup the input files
    //
    std::vector< std::tuple<AnalysisHelper::SampleType, std::string, float> > inputData;

    inputData.emplace_back(AnalysisHelper::Overlay, config.files.overlaysFileName, NormalisationHelper::GetOverlaysNormalisation(config));
    // inputData.emplace_back(AnalysisHelper::Dirt,    config.files.dirtFileName,     NormalisationHelper::GetDirtNormalisation(config));
    // inputData.emplace_back(AnalysisHelper::DataEXT, config.files.dataEXTFileName,  NormalisationHelper::GetDataEXTNormalisation(config));
    // inputData.emplace_back(AnalysisHelper::DataBNB, config.files.dataBNBFileName,  1.f);

    //
    // Setup the plots
    //
    const std::string yLabel = "Number of particles";
    
    // The highest energy proton variables are plotted with three different selections  
    const std::vector<std::string> protonSelectionNames{"protons>0", "protons=1", "protons>=2"};
    std::vector<PlottingHelper::MultiPlot> protonMomentumPlots

    for(auto const& selectionName : protonSelectionNames)
    {
        //TODO: Some of the values here are placeholders - make sure they are sensible
        protonMomentumPlots.emplace_back("Proton momentum / GeV (#" + selectionName + ")", yLabel, 50u, 0.f, 0.8f, true, config.global.axisTitles);
    }

    //
    // Get the selection
    //
    auto selection = SelectionHelper::GetDefaultSelection();

    // Loop over the events
    for (const auto [sampleType, fileName, normalisation] : inputData)
    {
        std::cout << "Reading input file: " << fileName << std::endl;

        FileReader reader(fileName);
        auto pEvent = reader.GetBoundEventAddress();

        const auto nEvents = reader.GetNumberOfEvents();
        for (unsigned int i = 0; i < nEvents; ++i)
        {
            AnalysisHelper::PrintLoadingBar(i, nEvents);
            reader.LoadEvent(i);

            // Run the event selection and store which cuts are passed
            const auto &[isSelectedGolden, cutsPassed, assignedPdgCodes] = selection.Execute(pEvent);
            const auto isSelectedGeneric = (std::find(cutsPassed.begin(), cutsPassed.end(), config.global.lastCutGeneric) != cutsPassed.end());

            // Only use events that at least pass the generic selection
            if (!isSelectedGeneric)
                continue;

            // Get the truth and reco analysis data
            const auto weight = AnalysisHelper::GetNominalEventWeight(pEvent) * normalisation;
            const auto recoData = AnalysisHelper::GetRecoAnalysisData(pEvent->reco, assignedPdgCodes, isSelectedGolden);
            const auto truthData = AnalysisHelper::GetTruthAnalysisData(pEvent->truth, config.global.useAbsPdg, config.global.protonMomentumThreshold);

            const auto plotStyle = PlottingHelper::GetPlotStyle(sampleType, pEvent, config.global.useAbsPdg);

            // Get the true origin of the selected muon and pion candidates
            const auto &recoParticles = pEvent->reco.particles;
            const auto &truthParticles = pEvent->truth.particles;

            // Not all selected events have protons. Use only those that do to plot reconstructed proton variables. 
            if(recoData.nProtons>0)
            {
                const auto proton = recoParticles.at(AnalysisHelper::GetHighestMomentumParticleIndexWithPdg(pEvent->reco, assignedPdgCodes, 2212));
                const auto protonPlotStyle = PlottingHelper::GetPlotStyle(proton, sampleType, truthParticles, false, config.global.useAbsPdg);

                for (unsigned int selection = 0; selection < protonSelectionNames.size(); ++selection)
                {
                    // Check proton multiplicity to match conditions described in protonSelectionNames
                    // No need to to test for selection==0 since if(recoData.nProtons>0) was used previously
                    if(selection==1 && recoData.nProtons!=1) continue;
                    if(selection==2 && recoData.nProtons<2) continue;

                    protonMomentumPlots.at(selection).Fill(recoData.protonMomentum, truthData.interactionMode, weight);
                }
            }
        }
    }

    for (unsigned int selection = 0; selection < protonSelectionNames.size(); ++selection)
    {
        const auto selectionName = protonSelectionNames.at(selection);
        protonMomentumPlots.at(selection).SaveAsStacked("reco_protonMomentum_" + selectionName,false,false,config.global.axisTitles);
    }
}

} // namespace ubcc1pi_macros
