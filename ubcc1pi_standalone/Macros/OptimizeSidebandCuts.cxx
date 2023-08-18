/**
 *  @file  ubcc1pi_standalone/Macros/OptimizeSidebandCuts.cxx
 *
 *  @brief The implementation file of the OptimizeSidebandCuts macro
 */

#include "ubcc1pi_standalone/Macros/Macros.h"

#include "ubcc1pi_standalone/Objects/FileReader.h"

#include "ubcc1pi_standalone/Helpers/PlottingHelper.h"
#include "ubcc1pi_standalone/Helpers/AnalysisHelper.h"
#include "ubcc1pi_standalone/Helpers/SelectionHelper.h"
#include "ubcc1pi_standalone/Helpers/NormalisationHelper.h"
// #include <fstream> // Todo: not use txt files
#include <iomanip>

#include <TH2F.h>
#include <TMath.h>

float muonCutValue(const int k) { return  -0.3f-k*0.02f; }
float protonCutValue(const int k) { return 0.0f+k*0.02f; }
// float muonCutValue(const int k) { return  -0.3f-k*0.01f; }
// float protonCutValue(const int k) { return 0.4f+k*0.01f; }

using namespace ubcc1pi;

namespace ubcc1pi_macros
{

void OptimizeSidebandCuts(const Config &config)
{
    ROOT::EnableImplicitMT(2);
    // Setup the input files
    std::vector< std::tuple<AnalysisHelper::SampleType, std::string, float> > inputData;

    std::cout<<"##########################################\nUSING NUWRO AS DATA & Only CC0pi!\n##########################################"<<std::endl;
    for (const auto run: config.global.runs)
    {
        if(run == 1)
        {
            inputData.emplace_back(AnalysisHelper::Overlay, config.filesRun1.overlaysFileName, NormalisationHelper::GetOverlaysNormalisationToNuWro(config, 1));
        }
        else if(run == 2)
        {
            inputData.emplace_back(AnalysisHelper::Overlay, config.filesRun2.overlaysFileName, NormalisationHelper::GetOverlaysNormalisationToNuWro(config, 2));
        }
        else if(run == 3)
        {
            inputData.emplace_back(AnalysisHelper::Overlay, config.filesRun3.overlaysFileName, NormalisationHelper::GetOverlaysNormalisationToNuWro(config, 3));
        }
        else throw std::logic_error("PlotEventSelectionCuts - Invalid run number");
    }

    // Get the selections
    const std::string yLabel = "Number of particles (norm. to bin width)";
    auto selectionPart1 = SelectionHelper::GetCC0piSelectionModifiedPart1();
    auto selectionPart2 = SelectionHelper::GetCC0piSelectionModifiedPart2(0.f, 0.f);

    std::map<std::pair<int, int>, TH1F> TrueCC0pi_ProtonMomentum_SelCC0pi, TrueCC0pi_nProtons_SelCC0pi;
    for(int i=0; i<15; i++)
    {
        for(int j=0; j<15; j++)
        {
            const auto coords = std::make_pair(i, j);
            auto histProtonMomentum = TH1F(("TrueCC0pi_SelCC0pi_ProtonMomentum"+std::to_string(i)+"+"+std::to_string(j)).c_str(),(string(";Proton Momentum / GeV;")+yLabel).c_str(), 15, 0, 1.5);
            auto histnProtons = TH1F(("TrueCC0pi_SelCC0pi_nProtons"+std::to_string(i)+"+"+std::to_string(j)).c_str(),(string(";Proton Multiplicity (-1);")+yLabel).c_str(), 3, 0, 3);
            TrueCC0pi_ProtonMomentum_SelCC0pi.emplace(coords, histProtonMomentum);
            TrueCC0pi_nProtons_SelCC0pi.emplace(coords, histnProtons);
        }
    }

    auto TrueCC0pi_ProtonMomentum_SelCC1pi = std::make_shared<TH1F>("TrueCC0pi_ProtonMomentum_SelCC1pi",(string(";Proton Momentum / GeV;")+yLabel).c_str(), 15, 0, 1.5);
    auto TrueCC0pi_nProtons_SelCC1pi = std::make_shared<TH1F>("TrueCC0pi_nProtons_SelCC1pi",(string(";Proton Multiplicity (-1);")+yLabel).c_str(), 3, 0, 3);

    auto CC1piSelection = SelectionHelper::GetSelection("Default");

    const auto featureNames = BDTHelper::ParticleBDTFeatureNames;
    // Get the muon BDT
    const auto muonFeatureNames = BDTHelper::MuonBDTFeatureNames;
    BDTHelper::BDT muonBDT("muon", muonFeatureNames);

    // Get the proton BDT
    const auto protonFeatureNames = BDTHelper::ProtonBDTFeatureNames;
    BDTHelper::BDT protonBDT("proton", protonFeatureNames);

    // Get the golden pion BDT
    const auto goldenPionFeatureNames = BDTHelper::GoldenPionBDTFeatureNames;
    BDTHelper::BDT goldenpionBDT("goldenPion", goldenPionFeatureNames);

    auto summedWeights = 0.f;
    // Loop over the events
    for (const auto [sampleType, fileName, normalisation] : inputData)
    {
        std::cout << "Reading input file: " << fileName << std::endl;

        FileReader reader(fileName);
        auto pEvent = reader.GetBoundEventAddress();

        const auto nEvents = reader.GetNumberOfEvents();
        std::cout<<"\n##############\nOnly counting every 6th event!\n##############"<<std::endl;
        for (unsigned int x = 0; x < nEvents/6; ++x)
        {
            AnalysisHelper::PrintLoadingBar(x, nEvents);
            reader.LoadEvent(x);

            if (!pEvent->reco.passesCCInclusive())
                continue;

            // Get the truth and reco analysis data
            const auto plotStyle = PlottingHelper::GetPlotStyle(sampleType, pEvent, config.global.useAbsPdg);

            // First, check if event is true CC0pi
            if (plotStyle == PlottingHelper::NumuCC0Pi)
            {
                const auto recoParticles = pEvent->reco.particles;
                const auto weight = AnalysisHelper::GetNominalEventWeight(pEvent) * normalisation;
                // Get the muon reco data
                const auto muonIndex = SelectionHelper::GetMuonCandidateIndex(recoParticles, muonFeatureNames, muonBDT);
                const auto muon = recoParticles.at(muonIndex);
                const auto muonDir = TVector3(muon.directionX(), muon.directionY(), muon.directionZ()).Unit();
                const float muonCosTheta = muonDir.Z();
                const float muonPhi = std::atan2(muonDir.Y(), muonDir.X());
                const float muonMomentum = AnalysisHelper::GetMuonMomentum(muon);

                // Get the proton reco data
                float protonMomentum = -std::numeric_limits<float>::max();
                float protonCosTheta = -std::numeric_limits<float>::max();
                float protonPhi = -std::numeric_limits<float>::max();

                const auto visibleParticles = AnalysisHelper::SelectVisibleParticles(pEvent->truth.particles);
                const auto nProtons = std::min(AnalysisHelper::CountParticlesAboveMomentumThreshold(visibleParticles, 2212, config.global.useAbsPdg, config.global.protonMomentumThreshold)-1, 2u); //-1 as one proton is treated as the pion

                const auto &[isSelectedPart1, cutsPassedPart1, assignedPdgCodesPart1] = selectionPart1.Execute(pEvent);
                if(isSelectedPart1)
                {
                    auto leadpidx = AnalysisHelper::GetTrueLeadingProtonIndex(pEvent->truth, config.global.useAbsPdg, config.global.protonMomentumThreshold);
                    if (leadpidx!=std::numeric_limits<unsigned int>::max())
                    {
                        summedWeights += weight;
                        for(int i=0; i<15; i++)
                        {
                            for(int j=0; j<15; j++)
                            {
                                // Run the event selection and store which cuts are passed
                                const auto coords = std::make_pair(i, j);
                                selectionPart2.SetCutValue("muonLikeProton", muonCutValue(i));
                                selectionPart2.SetCutValue("barelyResemblingProton", protonCutValue(j));
                                const auto &[isSelectedPart2, cutsPassed, assignedPdgCodes] = selectionPart2.Execute(pEvent);

                                // If event passes CC0pi selection, fill plots
                                if (isSelectedPart2)
                                {
                                    // For proton kinematics, use true leading proton. If no true leading proton is found, don't fill any kinematics
                                    const auto trueleadp = pEvent->truth.particles.at(leadpidx);
                                    TrueCC0pi_ProtonMomentum_SelCC0pi[coords].Fill(trueleadp.momentum(), weight);
                                    TrueCC0pi_nProtons_SelCC0pi[coords].Fill(nProtons, weight);
                                }
                            }
                        }
                    }
                }

                // Find out if the event passes the generic CC1pi selection
                const auto &[isSelectedCC1pi, cutsPassedCC1pi, assignedPdgCodesCC1pi] = CC1piSelection.Execute(pEvent);
                const auto &isSelectedCC1piGeneric = (std::find(cutsPassedCC1pi.begin(),cutsPassedCC1pi.end(), config.global.lastCutGeneric) != cutsPassedCC1pi.end());

                // If event passes CC1pi+ selection, fill plots
                if (isSelectedCC1piGeneric){
                    // For proton kinematics, use proton that is selected as the pi+ candidate. If the pi+ candidate is a true proton, fill muon and proton kinematic plots
                    const auto recoPion = recoParticles.at(AnalysisHelper::GetParticleIndexWithPdg(assignedPdgCodesCC1pi, 211));
                    const auto recoMuon = recoParticles.at(AnalysisHelper::GetParticleIndexWithPdg(assignedPdgCodesCC1pi, 13));
                    Event::Truth::Particle pionMatch;

                    try{
                        pionMatch = AnalysisHelper::GetBestMatchedTruthParticle(recoPion,pEvent->truth.particles);
                    }
                    catch(const std::logic_error &){
                        continue;
                    }

                    if (pionMatch.pdgCode()==2212){
                        TrueCC0pi_ProtonMomentum_SelCC1pi->Fill(pionMatch.momentum(), weight);
                        TrueCC0pi_nProtons_SelCC1pi->Fill(nProtons, weight);
                    }
                }
            } // end if true numu CC0pi

        }
    }

    std::cout<<"-----summedWeights: "<<summedWeights<<std::endl;

    // TrueCC0pi_ProtonMomentum_SelCC1pi->Sumw2();
    const auto selCC1piIntegralProtonMomentum = TrueCC0pi_ProtonMomentum_SelCC1pi->Integral();
    const auto selCC1piIntegralnProtons = TrueCC0pi_nProtons_SelCC1pi->Integral();
    TrueCC0pi_ProtonMomentum_SelCC1pi->Scale(1.0/selCC1piIntegralProtonMomentum);
    TrueCC0pi_nProtons_SelCC1pi->Scale(1.0/selCC1piIntegralnProtons);
    std::cout<<"cc0piIntegralProtonMomentum: "<<selCC1piIntegralProtonMomentum<<std::endl;
    std::cout<<"cc0piIntegralnProtons: "<<selCC1piIntegralnProtons<<std::endl;
    std::map<std::pair<int, int>, float> selcc0piIntegralProtonMomentum, selcc0piIntegralnProtons;

    std::cout<<"Results (proton momentum):"<<std::endl;
    std::vector<float> values(6, std::numeric_limits<float>::max());
    std::vector<float> lowestValues(6, std::numeric_limits<float>::max());
    std::vector<std::pair<int, int>> bestCoords(6, std::make_pair(-1,-1));
    // std::pair<int, int> bestCoordsProtonMomentumOnly, bestCoordsEfficiencyWeightedProtonMomentumOnly, bestCoords, bestCoordsEfficiencyWeighted;
    std::cout<<"\t"<<std::scientific<<std::setprecision(2);
    for(int j=0; j<15; j++)
        std::cout<<"j:"<<protonCutValue(j)<<"\t";
    std::cout<<std::endl;
    for(int i=0; i<15; i++)
    {
        std::cout<<"i:"<<muonCutValue(i)<<"\t";
        for(int j=0; j<15; j++)
        {
            const auto coords = std::make_pair(i, j);

            const auto cc0piIntegralProtonMomentum = TrueCC0pi_ProtonMomentum_SelCC0pi[coords].Integral();
            const auto cc0piIntegralnProtons = TrueCC0pi_nProtons_SelCC0pi[coords].Integral();
            selcc0piIntegralProtonMomentum.emplace(coords, cc0piIntegralProtonMomentum);
            selcc0piIntegralnProtons.emplace(coords, cc0piIntegralnProtons);
            // std::cout<<"cc0piIntegralProtonMomentum - i: "<<i<<" - j: "<<j<<" "<<cc0piIntegralProtonMomentum<<std::endl;
            TrueCC0pi_ProtonMomentum_SelCC0pi[coords].Scale(1.0/cc0piIntegralProtonMomentum);
            TrueCC0pi_ProtonMomentum_SelCC0pi[coords].GetYaxis()->SetRangeUser(0,0.23);
            TrueCC0pi_nProtons_SelCC0pi[coords].Scale(1.0/cc0piIntegralnProtons);
            TrueCC0pi_nProtons_SelCC0pi[coords].GetYaxis()->SetRangeUser(0,0.8);

            const auto efficiencyWeightProtonMomentum = selCC1piIntegralProtonMomentum/cc0piIntegralProtonMomentum;
            const auto efficiencyWeightnProtons = selCC1piIntegralnProtons/cc0piIntegralnProtons;

            const auto subtractedProtonMomentum = *TrueCC0pi_ProtonMomentum_SelCC1pi - TrueCC0pi_ProtonMomentum_SelCC0pi[coords];
            const auto subtractednProtons = *TrueCC0pi_nProtons_SelCC1pi - TrueCC0pi_nProtons_SelCC0pi[coords];
            const auto normalisedDifferenceProtonMomentum = subtractedProtonMomentum/TrueCC0pi_ProtonMomentum_SelCC0pi[coords];
            const auto normalisedDifferencenProtons = subtractednProtons/TrueCC0pi_nProtons_SelCC0pi[coords];
            const auto squaredProtonMomentum = normalisedDifferenceProtonMomentum*normalisedDifferenceProtonMomentum;
            const auto squarednProtons = normalisedDifferencenProtons*normalisedDifferencenProtons;
            const auto squaredProtonMomentumIntegral = squaredProtonMomentum.Integral()/squaredProtonMomentum.GetNbinsX();
            const auto squarednProtonsIntegral = squarednProtons.Integral()/squarednProtons.GetNbinsX();
            // const auto value = squared.Integral()/squared.GetNbinsX();

            std::cout<<"i: "<<i<<"j: "<<j<<" eff.W.nP: "<<efficiencyWeightnProtons<<" sqrdPMomentum: "<<squaredProtonMomentumIntegral<<" sqrdnP: "<<squarednProtonsIntegral<<std::endl;
            // std::cout<<"subtractedProtonMomentum: "<<subtractedProtonMomentum.Integral()<<"  subtractednProtons: "<<subtractednProtons.Integral()<<std::endl;
            // std::cout<<"normalisedDifferenceProtonMomentum: "<<normalisedDifferenceProtonMomentum.Integral()<<"  normalisedDifferencenProtons: "<<normalisedDifferencenProtons.Integral()<<std::endl;
            // std::cout<<"squaredProtonMomentum: "<<squaredProtonMomentumIntegral<<" squarednProtons: "<<squarednProtonsIntegral<<std::endl;

            values[0] = squaredProtonMomentumIntegral;
            values[1] = efficiencyWeightProtonMomentum*squaredProtonMomentumIntegral*squaredProtonMomentumIntegral;
            values[2] = squaredProtonMomentumIntegral*squaredProtonMomentumIntegral*squarednProtonsIntegral;
            values[3] = efficiencyWeightProtonMomentum*squaredProtonMomentumIntegral*squarednProtonsIntegral*squaredProtonMomentumIntegral;
            values[4] = efficiencyWeightProtonMomentum*std::sqrt(squarednProtonsIntegral)*squaredProtonMomentumIntegral;
            values[5] = std::sqrt(efficiencyWeightProtonMomentum*squarednProtonsIntegral)*squaredProtonMomentumIntegral;
            for(int m=0; m<6; m++)
            {
                if(values[m] < lowestValues[m])
                {
                    lowestValues[m] = values[m];
                    bestCoords[m] = coords;
                }
            }
        }
        std::cout<<std::endl;
    }
    auto i = bestCoords[0].first;
    auto j = bestCoords[0].second;
    std::cout<<"----------------------------------------------------"<<std::endl;
    std::cout<<"Lowest value (proton momentum considered): "<<lowestValues[0]<<" at "<<i<<"(muonCutValue: "<<muonCutValue(i)<<"),"<<j<<"(protonCutValue: "<<protonCutValue(j)<<")"<<std::endl;
    std::cout<<"----------------------------------------------------"<<std::endl;

    i = bestCoords[1].first;
    j = bestCoords[1].second;
    std::cout<<"----------------------------------------------------"<<std::endl;
    std::cout<<"Lowest value (proton momentum + number of selected events considered): "<<lowestValues[1]<<" at "<<i<<"(muonCutValue: "<<muonCutValue(i)<<"),"<<j<<"(protonCutValue: "<<protonCutValue(j)<<")"<<std::endl;
    std::cout<<"----------------------------------------------------"<<std::endl;

    i = bestCoords[2].first;
    j = bestCoords[2].second;
    std::cout<<"----------------------------------------------------"<<std::endl;
    std::cout<<"Lowest value (proton momentum + proton multiplicity considered): "<<lowestValues[2]<<" at "<<i<<"(muonCutValue: "<<muonCutValue(i)<<"),"<<j<<"(protonCutValue: "<<protonCutValue(j)<<")"<<std::endl;
    std::cout<<"----------------------------------------------------"<<std::endl;

    i = bestCoords[3].first;
    j = bestCoords[3].second;
    std::cout<<"----------------------------------------------------"<<std::endl;
    std::cout<<"Lowest value (proton momentum + proton multiplicity + number of selected events considered): "<<lowestValues[3]<<" at "<<i<<"(muonCutValue: "<<muonCutValue(i)<<"),"<<j<<"(protonCutValue: "<<protonCutValue(j)<<")"<<std::endl;
    std::cout<<"----------------------------------------------------"<<std::endl;

    i = bestCoords[4].first;
    j = bestCoords[4].second;
    std::cout<<"----------------------------------------------------"<<std::endl;
    std::cout<<"Lowest value (proton momentum + sqrt(proton multiplicity) + number of selected events considered): "<<lowestValues[4]<<" at "<<i<<"(muonCutValue: "<<muonCutValue(i)<<"),"<<j<<"(protonCutValue: "<<protonCutValue(j)<<")"<<std::endl;
    std::cout<<"----------------------------------------------------"<<std::endl;

    i = bestCoords[5].first;
    j = bestCoords[5].second;
    std::cout<<"----------------------------------------------------"<<std::endl;
    std::cout<<"Lowest value (proton momentum * sqrt(proton multiplicity * number of selected events considered)): "<<lowestValues[5]<<" at "<<i<<"(muonCutValue: "<<muonCutValue(i)<<"),"<<j<<"(protonCutValue: "<<protonCutValue(j)<<")"<<std::endl;
    std::cout<<"----------------------------------------------------"<<std::endl;


    std::cout<<"# events:"<<std::endl;
    auto highestValue = -std::numeric_limits<float>::max();
    std::pair<int, int> bestCoordsNumEvents;
    std::cout<<"           \t";
    for(int j=0; j<15; j++)
        std::cout<<"j:"<<protonCutValue(j)<<"\t";
    std::cout<<std::endl;
    for(int i=0; i<15; i++)
    {
        std::cout<<"i:"<<muonCutValue(i)<<"\t";
        for(int j=0; j<15; j++)
        {
            const auto coords = std::make_pair(i, j);
            const auto value = selcc0piIntegralProtonMomentum[coords];
            if(value>highestValue)
            {
                highestValue = value;
                bestCoordsNumEvents = coords;
            }
            std::cout<<value<<"\t";
        }
        std::cout<<std::endl;
    }
    i = bestCoordsNumEvents.first;
    j = bestCoordsNumEvents.second;
    std::cout<<"----------------------------------------------------"<<std::endl;
    std::cout<<"Largest number of events (CC0pi selection): "<<highestValue<<" at "<<i<<"(muonCutValue: "<<muonCutValue(i)<<"),"<<j<<"(protonCutValue: "<<protonCutValue(j)<<")"<<std::endl;
    std::cout<<"----------------------------------------------------"<<std::endl;



    auto pCanvas = PlottingHelper::GetCanvas();

    PlottingHelper::SetLineStyle(TrueCC0pi_ProtonMomentum_SelCC1pi, PlottingHelper::Primary);
    PlottingHelper::SetLineStyle(&TrueCC0pi_ProtonMomentum_SelCC0pi[bestCoords[0]], PlottingHelper::Secondary);
    PlottingHelper::SetLineStyle(&TrueCC0pi_ProtonMomentum_SelCC0pi[bestCoords[1]], PlottingHelper::Secondary);
    PlottingHelper::SetLineStyle(&TrueCC0pi_ProtonMomentum_SelCC0pi[bestCoords[2]], PlottingHelper::Secondary);
    PlottingHelper::SetLineStyle(&TrueCC0pi_ProtonMomentum_SelCC0pi[bestCoords[3]], PlottingHelper::Secondary);
    PlottingHelper::SetLineStyle(&TrueCC0pi_ProtonMomentum_SelCC0pi[bestCoords[4]], PlottingHelper::Secondary);
    PlottingHelper::SetLineStyle(&TrueCC0pi_ProtonMomentum_SelCC0pi[bestCoords[5]], PlottingHelper::Secondary);
    PlottingHelper::SetLineStyle(TrueCC0pi_nProtons_SelCC1pi, PlottingHelper::Primary);
    PlottingHelper::SetLineStyle(&TrueCC0pi_nProtons_SelCC0pi[bestCoords[0]], PlottingHelper::Secondary);
    PlottingHelper::SetLineStyle(&TrueCC0pi_nProtons_SelCC0pi[bestCoords[1]], PlottingHelper::Secondary);
    PlottingHelper::SetLineStyle(&TrueCC0pi_nProtons_SelCC0pi[bestCoords[2]], PlottingHelper::Secondary);
    PlottingHelper::SetLineStyle(&TrueCC0pi_nProtons_SelCC0pi[bestCoords[3]], PlottingHelper::Secondary);
    PlottingHelper::SetLineStyle(&TrueCC0pi_nProtons_SelCC0pi[bestCoords[4]], PlottingHelper::Secondary);
    PlottingHelper::SetLineStyle(&TrueCC0pi_nProtons_SelCC0pi[bestCoords[5]], PlottingHelper::Secondary);

    TrueCC0pi_ProtonMomentum_SelCC0pi[bestCoords[0]].Draw("hist");
    TrueCC0pi_ProtonMomentum_SelCC1pi->Draw("hist same");
    PlottingHelper::SaveCanvas(pCanvas,"OptimizeSidebandCuts_ProtonMomentum_areanorm_raw_ProtonMomentumOnly");

    TrueCC0pi_nProtons_SelCC0pi[bestCoords[0]].Draw("hist");
    TrueCC0pi_nProtons_SelCC1pi->Draw("hist same");
    PlottingHelper::SaveCanvas(pCanvas,"OptimizeSidebandCuts_nProtons_areanorm_raw_ProtonMomentumOnly");

    TrueCC0pi_ProtonMomentum_SelCC0pi[bestCoords[1]].Draw("hist");
    TrueCC0pi_ProtonMomentum_SelCC1pi->Draw("hist same");
    PlottingHelper::SaveCanvas(pCanvas,"OptimizeSidebandCuts_ProtonMomentum_areanorm_weighted_ProtonMomentumOnly");

    TrueCC0pi_nProtons_SelCC0pi[bestCoords[1]].Draw("hist");
    TrueCC0pi_nProtons_SelCC1pi->Draw("hist same");
    PlottingHelper::SaveCanvas(pCanvas,"OptimizeSidebandCuts_nProtons_areanorm_weighted_ProtonMomentumOnly");

    TrueCC0pi_ProtonMomentum_SelCC0pi[bestCoords[2]].Draw("hist");
    TrueCC0pi_ProtonMomentum_SelCC1pi->Draw("hist same");
    PlottingHelper::SaveCanvas(pCanvas,"OptimizeSidebandCuts_ProtonMomentum_areanorm_raw");

    TrueCC0pi_nProtons_SelCC0pi[bestCoords[2]].Draw("hist");
    TrueCC0pi_nProtons_SelCC1pi->Draw("hist same");
    PlottingHelper::SaveCanvas(pCanvas,"OptimizeSidebandCuts_nProtons_areanorm_raw");

    TrueCC0pi_ProtonMomentum_SelCC0pi[bestCoords[3]].Draw("hist");
    TrueCC0pi_ProtonMomentum_SelCC1pi->Draw("hist same");
    PlottingHelper::SaveCanvas(pCanvas,"OptimizeSidebandCuts_ProtonMomentum_areanorm_weighted");

    TrueCC0pi_nProtons_SelCC0pi[bestCoords[3]].Draw("hist");
    TrueCC0pi_nProtons_SelCC1pi->Draw("hist same");
    PlottingHelper::SaveCanvas(pCanvas,"OptimizeSidebandCuts_nProtons_areanorm_weighted");

    TrueCC0pi_ProtonMomentum_SelCC0pi[bestCoords[4]].Draw("hist");
    TrueCC0pi_ProtonMomentum_SelCC1pi->Draw("hist same");
    PlottingHelper::SaveCanvas(pCanvas,"OptimizeSidebandCuts_ProtonMomentum_areanorm_weighted_sqrt1");

    TrueCC0pi_nProtons_SelCC0pi[bestCoords[4]].Draw("hist");
    TrueCC0pi_nProtons_SelCC1pi->Draw("hist same");
    PlottingHelper::SaveCanvas(pCanvas,"OptimizeSidebandCuts_nProtons_areanorm_weighted_sqrt1");

    TrueCC0pi_ProtonMomentum_SelCC0pi[bestCoords[5]].Draw("hist");
    TrueCC0pi_ProtonMomentum_SelCC1pi->Draw("hist same");
    PlottingHelper::SaveCanvas(pCanvas,"OptimizeSidebandCuts_ProtonMomentum_areanorm_weighted_sqrt2");

    TrueCC0pi_nProtons_SelCC0pi[bestCoords[5]].Draw("hist");
    TrueCC0pi_nProtons_SelCC1pi->Draw("hist same");
    PlottingHelper::SaveCanvas(pCanvas,"OptimizeSidebandCuts_nProtons_areanorm_weighted_sqrt2");
}

} // namespace ubcc1pi_macros
