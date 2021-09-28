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
// #include "ubcc1pi_standalone/Helpers/NormalisationHelper.h"


using namespace ubcc1pi;

namespace ubcc1pi_macros
{

void PlotProtonMomentumByInteraction(const Config &config)
{
    //
    // Setup the input files
    //
    // std::vector< std::tuple<AnalysisHelper::SampleType, std::string, float> > inputData;
    // inputData.emplace_back(AnalysisHelper::Overlay, config.files.overlaysFileName, NormalisationHelper::GetOverlaysNormalisation(config));

    FileReader reader(config.files.overlaysFileName);
    auto pEvent = reader.GetBoundEventAddress();
    const auto nEvents = reader.GetNumberOfEvents();

    //
    // Setup the plots
    //
    const std::string yLabel = "Number of particles";
    
    // The highest energy proton variables are plotted with three different selections  
    const std::vector<std::string> protonPlotNames{"protons>0", "protons==1", "protons>=2"};
    std::vector<PlottingHelper::MultiPlot> protonMomentumPlots;

    for(auto const& plotName : protonPlotNames)
    {
        //TODO: Some of the values here are placeholders - make sure they are sensible
        protonMomentumPlots.emplace_back("Proton momentum / GeV (#" + plotName + ")", yLabel, 50u, config.global.protonMomentum.min, config.global.protonMomentum.max, true, config.global.axisTitles);
    }

    //
    // Get the selection
    //
    auto selection = SelectionHelper::GetDefaultSelection();


    // // Loop over the events
    // for (const auto [sampleType, fileName, normalisation] : inputData)
    // {
    //     std::cout << "Reading input file: " << fileName << std::endl;

    //     FileReader reader(fileName);
    //     auto pEvent = reader.GetBoundEventAddress();

        // const auto nEvents = reader.GetNumberOfEvents();
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

        // Check that at least one proton is present 
        const auto &recoParticles = pEvent->reco.particles;
        const auto protonIndex = SelectionHelper::GetLeadingProtonCandidateIndex(recoParticles, assignedPdgCodes);
        if (protonIndex!=std::numeric_limits<unsigned int>::max())
        {
            const auto proton = recoParticles.at(protonIndex);
            //const auto truthData = AnalysisHelper::GetTruthAnalysisData(pEvent->truth, config.global.useAbsPdg, config.global.protonMomentumThreshold);
            const auto weight = AnalysisHelper::GetNominalEventWeight(pEvent);
            const auto recoData = AnalysisHelper::GetRecoAnalysisData(pEvent->reco, assignedPdgCodes, isSelectedGolden);
            // const auto &truthParticles = pEvent->truth.particles;
            // const auto protonPlotStyle = PlottingHelper::GetPlotStyle(proton, sampleType, truthParticles, false, config.global.useAbsPdg);
            // const auto plotStyle = PlottingHelper::GetPlotStyle(sampleType, pEvent, config.global.useAbsPdg);
            const auto protonMomentum = AnalysisHelper::GetProtonMomentumFromRange(proton.range());

            for (unsigned int plot = 0; plot < protonPlotNames.size(); ++plot)
            {
                // Check proton multiplicity to match conditions described in protonPlotNames
                // No need to to test for plot==0 since if(recoData.nProtons>0) was used previously
                if(plot==1 && recoData.nProtons!=1) continue;
                if(plot==2 && recoData.nProtons<2) continue;

                // const auto intMode = pEvent->truth.interactionMode();
                // const auto interactionMode = static_cast<PlottingHelper::PlotStyle>((intMode!=10) ? intMode:4);
                const auto style = PlottingHelper::GetPlotStyle(pEvent->truth.interactionMode());
                std::cout<<"Interaction mode "<<pEvent->truth.interactionMode()<<": "<<pEvent->truth.interactionString()<<std::endl;
                protonMomentumPlots.at(plot).Fill(protonMomentum, style, weight);
            }
        }
    }
    // }

    for (unsigned int plot = 0; plot < protonPlotNames.size(); ++plot)
    {
        const auto plotName = protonPlotNames.at(plot);
        protonMomentumPlots.at(plot).SaveAsStacked("reco_protonMomentumByInteraction_" + plotName,false,false,false,config.global.axisTitles);
    }
}

} // namespace ubcc1pi_macros
