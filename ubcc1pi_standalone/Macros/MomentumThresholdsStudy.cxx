/**
 *  @file  ubcc1pi_standalone/Macros/MomentumThresholdsStudy.cxx
 *
 *  @brief The implementation file of the MomentumThresholdsStudy macro
 */

#include "ubcc1pi_standalone/Macros/Macros.h"

#include "ubcc1pi_standalone/Objects/FileReader.h"

#include "ubcc1pi_standalone/Helpers/AnalysisHelper.h"
#include "ubcc1pi_standalone/Helpers/PlottingHelper.h"
#include "ubcc1pi_standalone/Helpers/SelectionHelper.h"
#include "ubcc1pi_standalone/Helpers/FormattingHelper.h"

#include <map>

#include <TH2F.h>

using namespace ubcc1pi;

namespace ubcc1pi_macros
{

void MomentumThresholdsStudy(const Config &config)
{
    // Open the file
    FileReader reader(config.files.overlaysFileName);
    auto pEvent = reader.GetBoundEventAddress();
    const auto nEvents = reader.GetNumberOfEvents();

    // Get the selection and switch off the momentum threshold cuts
    auto selection = SelectionHelper::GetDefaultSelection();

    // Setup the plots
    const std::vector<string> cuts({"withHit", "withRecoParticle"});
    PlottingHelper::EfficiencyPlot protonMomentumPlot("Proton momentum / GeV", 50u, 0.f, 1.f, cuts, false);

    // Loop over the events
    for (unsigned int i = 0; i < nEvents; ++i)
    {
        AnalysisHelper::PrintLoadingBar(i, nEvents);
        reader.LoadEvent(i);

        // Only consider true CC1Pi events
        if (!AnalysisHelper::IsTrueCC1Pi(pEvent, config.global.useAbsPdg))
            continue;

        // Insist we pass the CC inclusive
        if (!pEvent->reco.passesCCInclusive())
            continue;

        // Insist we pass the generic selection
        const auto &[passedGoldenSelection, cutsPassed, assignedPdgCodes] = selection.Execute(pEvent);
        const auto passedGenericSelection = SelectionHelper::IsCutPassed(cutsPassed, config.global.lastCutGeneric);

        if (!passedGenericSelection)
            continue;

        // Make the mapping from truth -> reco particles
        std::map<unsigned int, std::vector<unsigned int> > truthToRecoParticleIndices;
        const auto truthParticles = pEvent->truth.particles;
        const auto recoParticles = pEvent->reco.particles;
        for (unsigned int iReco = 0u; iReco < recoParticles.size(); ++iReco)
        {
            const auto &recoParticle = recoParticles.at(iReco);

            try
            {
                const auto iTrue = AnalysisHelper::GetBestMatchedTruthParticleIndex(recoParticle, truthParticles);
                truthToRecoParticleIndices[iTrue].push_back(iReco);
            }
            catch (const std::exception &)
            {
            }
        }

        // Loop over the true particles
        const auto weight = AnalysisHelper::GetNominalEventWeight(pEvent);
        for (unsigned int iTrue = 0u; iTrue < truthParticles.size(); ++iTrue)
        {
            const auto &particle = truthParticles.at(iTrue);

            // Insist the particle is a proton
            const auto pdg = config.global.useAbsPdg ? std::abs(particle.pdgCode()) : particle.pdgCode();
            if (pdg != 2212)
                continue;

            // Get the true momentum
            const auto momentum = particle.momentum();

            // Determine if the particle has a hit
            const auto hasHit = (particle.hitWeightU() > 1.f) || (particle.hitWeightV() > 1.f) || (particle.hitWeightW() > 1.f);
            protonMomentumPlot.AddEvent(momentum, weight, "withHit", hasHit);

            // Determine if the particle has a matched reco particle
            const auto &recoParticleIndices = truthToRecoParticleIndices[iTrue];
            const auto hasRecoParticle = !recoParticleIndices.empty();
            protonMomentumPlot.AddEvent(momentum, weight, "withRecoParticle", hasRecoParticle);
        }
    }

    // Setup the styles
    const std::vector<PlottingHelper::PlotStyle> styles({
        PlottingHelper::Primary,                                // withHit
        PlottingHelper::Secondary                               // withRecoParticle
    });
    protonMomentumPlot.SaveAs(cuts, styles, "momentumThresholds_proton");
}

} // namespace ubcc1pi_macros
