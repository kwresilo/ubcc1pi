/**
 *  @file  ubcc1pi/Analyzers/EventSelection_module.cc
 *
 *  @brief The implementation file for the event selection analyzer.
 */

#include "ubcc1pi/Analyzers/EventSelection.h"

#include "ubcc1pi/Helpers/TruthHelper.h"
#include "ubcc1pi/Helpers/AnalysisHelper.h"
#include "ubcc1pi/Helpers/DebugHelper.h"

namespace ubcc1pi
{

EventSelection::EventSelection(const art::EDAnalyzer::Table<Config> &config) :
    art::EDAnalyzer(config),
    m_config(config)
{
    // Setup the output trees
    art::ServiceHandle<art::TFileService> fileService;

    // Event tree
    m_pEventTree = fileService->make<TTree>("events", "");

    // Metadata
    m_pEventTree->Branch("run", &m_outputEvent.m_run);
    m_pEventTree->Branch("subRun", &m_outputEvent.m_subRun);
    m_pEventTree->Branch("event", &m_outputEvent.m_event);

    // Interaction info
    m_pEventTree->Branch("isSignal", &m_outputEvent.m_isSignal);
    m_pEventTree->Branch("interaction", &m_outputEvent.m_interaction);
    m_pEventTree->Branch("isNuFiducial", &m_outputEvent.m_isNuFiducial);
    m_pEventTree->Branch("nuE", &m_outputEvent.m_nuE);
    m_pEventTree->Branch("nuX", &m_outputEvent.m_nuX);
    m_pEventTree->Branch("nuY", &m_outputEvent.m_nuY);
    m_pEventTree->Branch("nuZ", &m_outputEvent.m_nuZ);

    // True particle multiplicities
    m_pEventTree->Branch("nMuMinus", &m_outputEvent.m_nMuMinus);
    m_pEventTree->Branch("nMuPlus", &m_outputEvent.m_nMuPlus);
    m_pEventTree->Branch("nPiPlus", &m_outputEvent.m_nPiPlus);
    m_pEventTree->Branch("nPiMinus", &m_outputEvent.m_nPiMinus);
    m_pEventTree->Branch("nKPlus", &m_outputEvent.m_nKPlus);
    m_pEventTree->Branch("nKMinus", &m_outputEvent.m_nKMinus);
    m_pEventTree->Branch("nProton", &m_outputEvent.m_nProton);
    m_pEventTree->Branch("nNeutron", &m_outputEvent.m_nNeutron);
    m_pEventTree->Branch("nPhoton", &m_outputEvent.m_nPhoton);
    m_pEventTree->Branch("nTotal", &m_outputEvent.m_nTotal);

    // True muon kinematics
    m_pEventTree->Branch("trueMuEnergy", &m_outputEvent.m_trueMuEnergy);
    m_pEventTree->Branch("trueMuTheta", &m_outputEvent.m_trueMuTheta);
    m_pEventTree->Branch("trueMuPhi", &m_outputEvent.m_trueMuPhi);
    
    // True pion kinematics
    m_pEventTree->Branch("truePiEnergy", &m_outputEvent.m_truePiEnergy);
    m_pEventTree->Branch("truePiTheta", &m_outputEvent.m_truePiTheta);
    m_pEventTree->Branch("truePiPhi", &m_outputEvent.m_truePiPhi);
    m_pEventTree->Branch("trueMuPiAngle", &m_outputEvent.m_trueMuPiAngle);
    
    // Reconstruction details
    m_pEventTree->Branch("recoNuX", &m_outputEvent.m_recoNuX);
    m_pEventTree->Branch("recoNuY", &m_outputEvent.m_recoNuY);
    m_pEventTree->Branch("recoNuZ", &m_outputEvent.m_recoNuZ);
    m_pEventTree->Branch("isRecoNuFiducial", &m_outputEvent.m_isRecoNuFiducial);
    
    m_pEventTree->Branch("nFinalStatePFPs", &m_outputEvent.m_nFinalStatePFPs);
    m_pEventTree->Branch("nRecoMIPs", &m_outputEvent.m_nRecoMIPs);

    // Selection details
    m_pEventTree->Branch("isSelected", &m_outputEvent.m_isSelected);

    m_pEventTree->Branch("isSelectedMuonTruthMatched", &m_outputEvent.m_isSelectedMuonTruthMatched);
    m_pEventTree->Branch("selectedMuonTruePdg", &m_outputEvent.m_selectedMuonTruePdg);
    m_pEventTree->Branch("selectedMuonCompleteness", &m_outputEvent.m_selectedMuonCompleteness);
    m_pEventTree->Branch("selectedMuonPurity", &m_outputEvent.m_selectedMuonPurity);
    
    m_pEventTree->Branch("isSelectedPionTruthMatched", &m_outputEvent.m_isSelectedPionTruthMatched);
    m_pEventTree->Branch("selectedPionTruePdg", &m_outputEvent.m_selectedPionTruePdg);
    m_pEventTree->Branch("selectedPionCompleteness", &m_outputEvent.m_selectedPionCompleteness);
    m_pEventTree->Branch("selectedPionPurity", &m_outputEvent.m_selectedPionPurity);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void EventSelection::analyze(const art::Event &event)
{
    // Extract the labels
    const auto mcTruthLabel = m_config().MCTruthLabel();
    const auto mcParticleLabel = m_config().MCParticleLabel();
    const auto backtrackerLabel = m_config().BacktrackerLabel();
    const auto pfParticleLabel = m_config().PFParticleLabel();
    const auto trackLabel = m_config().TrackLabel();
    const auto pidLabel = m_config().PIDLabel();
    const auto chi2ProtonCut = m_config().Chi2ProtonMIPCut();
    
    // Get the final state PFParticles
    const auto allPFParticles = CollectionHelper::GetCollection<recob::PFParticle>(event, pfParticleLabel);
    const auto finalStates = RecoHelper::GetNeutrinoFinalStates(allPFParticles);
    
    // Get the true neutrino interaction from the event
    const TruthHelper::Interaction interaction(event, mcTruthLabel, mcParticleLabel);
    const auto isTrueSignal = AnalysisHelper::IsCC1PiSignal(interaction);
   
    // Get the truth matching information for the PFParticles
    const auto backtrackerData = AnalysisHelper::GetBacktrackerData(event, mcTruthLabel, mcParticleLabel, backtrackerLabel, pfParticleLabel);

    // Do the event selection
    PFParticleVector muons, pions, protons, showerLikes;
    this->PerformPID(event, finalStates, chi2ProtonCut, muons, pions, protons, showerLikes);

    const auto recoVertex = this->GetRecoNeutrinoVertex(event, allPFParticles);
    const bool isRecoNuFiducial = AnalysisHelper::IsFiducial(recoVertex);
    const bool isSelected = (muons.size() == 1 && pions.size() == 1 && showerLikes.size() == 0 && isRecoNuFiducial);

    // Output event level information
    m_outputEvent.m_run = event.run();
    m_outputEvent.m_subRun = event.subRun();
    m_outputEvent.m_event = event.event();
    m_outputEvent.m_isSignal = isTrueSignal;
    m_outputEvent.m_interaction = DebugHelper::GetInteractionString(interaction, true);
    m_outputEvent.m_isNuFiducial = AnalysisHelper::IsNeutrinoVertexFiducial(interaction);

    const auto mcNeutrino = interaction.GetNeutrino();
    m_outputEvent.m_nuE = mcNeutrino.E();
    m_outputEvent.m_nuX = mcNeutrino.Vx(); 
    m_outputEvent.m_nuY = mcNeutrino.Vy(); 
    m_outputEvent.m_nuZ = mcNeutrino.Vz(); 

    const auto mcParticles = AnalysisHelper::GetReconstructableFinalStates(interaction);
    m_outputEvent.m_nMuMinus = AnalysisHelper::CountParticlesWithPDG(mcParticles, 13);
    m_outputEvent.m_nMuPlus = AnalysisHelper::CountParticlesWithPDG(mcParticles, -13);
    m_outputEvent.m_nPiPlus = AnalysisHelper::CountParticlesWithPDG(mcParticles, 211);
    m_outputEvent.m_nPiMinus = AnalysisHelper::CountParticlesWithPDG(mcParticles, -211);
    m_outputEvent.m_nKPlus = AnalysisHelper::CountParticlesWithPDG(mcParticles, 321);
    m_outputEvent.m_nKMinus = AnalysisHelper::CountParticlesWithPDG(mcParticles, -321);  
    m_outputEvent.m_nProton = AnalysisHelper::CountParticlesWithPDG(mcParticles, 2212); 
    m_outputEvent.m_nNeutron = AnalysisHelper::CountParticlesWithPDG(mcParticles, 2112);
    m_outputEvent.m_nPhoton = AnalysisHelper::CountParticlesWithPDG(mcParticles, 22);
    m_outputEvent.m_nTotal = mcParticles.size();
    
    // Set the MCParticle kinematic information
    m_outputEvent.m_trueMuEnergy = -std::numeric_limits<float>::max();
    m_outputEvent.m_trueMuTheta = -std::numeric_limits<float>::max();
    m_outputEvent.m_trueMuPhi = -std::numeric_limits<float>::max();
    
    m_outputEvent.m_truePiEnergy = -std::numeric_limits<float>::max();
    m_outputEvent.m_truePiTheta = -std::numeric_limits<float>::max();
    m_outputEvent.m_truePiPhi = -std::numeric_limits<float>::max();
    m_outputEvent.m_trueMuPiAngle = -std::numeric_limits<float>::max(); 

    if (isTrueSignal)
    {
        TVector3 muAngle, piAngle;

        for (const auto &mcParticle : mcParticles)
        {
            const auto mcpE = mcParticle->E();
            const auto mcpTheta = (mcParticle->P() < std::numeric_limits<float>::epsilon()) ? -std::numeric_limits<float>::max() : std::acos(mcParticle->Pz() / mcParticle->P());
            const auto mcpPhi = std::atan2(mcParticle->Px(), mcParticle->Py());

            switch (mcParticle->PdgCode())
            {
                case 13:
                    m_outputEvent.m_trueMuEnergy = mcpE;
                    m_outputEvent.m_trueMuTheta = mcpTheta; 
                    m_outputEvent.m_trueMuPhi = mcpPhi;
                    muAngle = TVector3(mcParticle->Px(), mcParticle->Py(), mcParticle->Pz()).Unit();
                    break;
                case 211:
                    m_outputEvent.m_truePiEnergy = mcpE; 
                    m_outputEvent.m_truePiTheta = mcpTheta;
                    m_outputEvent.m_truePiPhi = mcpPhi;
                    piAngle = TVector3(mcParticle->Px(), mcParticle->Py(), mcParticle->Pz()).Unit();
                    break;
                default: break;
            }
        }

        m_outputEvent.m_trueMuPiAngle = std::acos(muAngle.Dot(piAngle));
    }

    // Set the reconstructed vertex position
    m_outputEvent.m_recoNuX = recoVertex.X();
    m_outputEvent.m_recoNuY = recoVertex.Y();
    m_outputEvent.m_recoNuZ = recoVertex.Z();

    // Set the selection information
    m_outputEvent.m_nFinalStatePFPs = finalStates.size();
    m_outputEvent.m_isRecoNuFiducial = isRecoNuFiducial;
    m_outputEvent.m_nRecoMIPs = muons.size() + pions.size();
    m_outputEvent.m_isSelected = isSelected;

    // The selected particle info
    m_outputEvent.m_isSelectedMuonTruthMatched = false;
    m_outputEvent.m_selectedMuonTruePdg = -std::numeric_limits<int>::max();
    m_outputEvent.m_selectedMuonCompleteness = -std::numeric_limits<float>::max();
    m_outputEvent.m_selectedMuonPurity = -std::numeric_limits<float>::max();
    
    m_outputEvent.m_isSelectedPionTruthMatched = false;
    m_outputEvent.m_selectedPionTruePdg = -std::numeric_limits<int>::max();
    m_outputEvent.m_selectedPionCompleteness = -std::numeric_limits<float>::max();
    m_outputEvent.m_selectedPionPurity = -std::numeric_limits<float>::max();
    
    if (isSelected)
    {
        // Save the muon information
        try
        {
            const auto matchedMCP = backtrackerData.GetBestMatchedMCParticle(muons.front());
            m_outputEvent.m_selectedMuonTruePdg = matchedMCP->PdgCode();
            m_outputEvent.m_selectedMuonCompleteness = backtrackerData.GetMatchCompleteness(muons.front(), matchedMCP);
            m_outputEvent.m_selectedMuonPurity = backtrackerData.GetMatchPurity(muons.front(), matchedMCP);
            m_outputEvent.m_isSelectedMuonTruthMatched = true;
        }
        catch (const cet::exception &)
        {
        }
        
        // Save the pion information
        try
        {
            const auto matchedMCP = backtrackerData.GetBestMatchedMCParticle(pions.front());
            m_outputEvent.m_selectedPionTruePdg = matchedMCP->PdgCode();
            m_outputEvent.m_selectedPionCompleteness = backtrackerData.GetMatchCompleteness(pions.front(), matchedMCP);
            m_outputEvent.m_selectedPionPurity = backtrackerData.GetMatchPurity(pions.front(), matchedMCP);
            m_outputEvent.m_isSelectedPionTruthMatched = true;
        }
        catch (const cet::exception &)
        {
        }
    }

    m_pEventTree->Fill();
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void EventSelection::PerformPID(const art::Event &event, const PFParticleVector &finalStates, const float chi2ProtonCut, PFParticleVector &muons, PFParticleVector &pions, PFParticleVector &protons, PFParticleVector &showerLikes) const
{
    if (!muons.empty() || !pions.empty() || !protons.empty())
        throw cet::exception("EventSelection::PerformPID") << " - Output vectors for muons, pions and protons were not all empty" << std::endl;
    
    // No need to do anything if there are no final states supplied
    if (finalStates.empty())
        return;

    // Split the PFParticles into track and shower like based on the decision made in Pandora
    PFParticleVector trackLikes;
    this->SelectTracksAndShowers(finalStates, trackLikes, showerLikes);

    // Select the PFParticles that pass the MIP selection
    const auto mips = this->SelectMIPs(event, trackLikes, chi2ProtonCut);

    // If there are no mips, then just call everything a proton
    if (mips.empty())
    {
        protons = finalStates;
        return;
    }
    
    // Get the muon candidate
    const auto muon = this->SelectMuon(event, mips);

    // Populate the output vectors
    for (const auto &pfParticle : finalStates)
    {
        if (pfParticle == muon)
        {
            muons.push_back(pfParticle);
            continue;
        }

        // Any non-muon MIP is a pion
        if (std::find(mips.begin(), mips.end(), pfParticle) != mips.end())
        {
            pions.push_back(pfParticle);
            continue;
        }

        // Any non-MIP track is a proton
        if (this->IsTrackLike(pfParticle))
        {
            protons.push_back(pfParticle);
            continue;
        }
    }
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void EventSelection::SelectTracksAndShowers(const PFParticleVector &finalStates, PFParticleVector &trackLikes, PFParticleVector &showerLikes) const
{
    for (const auto &pfParticle : finalStates)
    {
        if (this->IsTrackLike(pfParticle))
            trackLikes.push_back(pfParticle);
        
        if (this->IsShowerLike(pfParticle))
            showerLikes.push_back(pfParticle);
    }
}

// -----------------------------------------------------------------------------------------------------------------------------------------

PFParticleVector EventSelection::SelectMIPs(const art::Event &event, const PFParticleVector &finalStates, const float chi2ProtonCut) const
{
    // Get the required associations
    const auto pfpToTracks = CollectionHelper::GetAssociation<recob::PFParticle, recob::Track>(event, m_config().PFParticleLabel(), m_config().TrackLabel());
    const auto trackToPIDs = CollectionHelper::GetAssociation<recob::Track, anab::ParticleID>(event, m_config().TrackLabel(), m_config().PIDLabel());

    // Select MIPs
    PFParticleVector mipCandidates;
    for (const auto &pfParticle : finalStates)
    {
        if (this->PassesMIPSelection(pfParticle, pfpToTracks, trackToPIDs, chi2ProtonCut))
            mipCandidates.push_back(pfParticle);
    }

    return mipCandidates;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

bool EventSelection::IsTrackLike(const art::Ptr<recob::PFParticle> &pfParticle) const
{
    const auto absPdg = std::abs(pfParticle->PdgCode());
    return (absPdg == 13 || absPdg == 211 || absPdg == 2212 || absPdg == 321);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

bool EventSelection::IsShowerLike(const art::Ptr<recob::PFParticle> &pfParticle) const
{
    const auto absPdg = std::abs(pfParticle->PdgCode());
    return (absPdg == 11 || absPdg == 22);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

bool EventSelection::PassesMIPSelection(const art::Ptr<recob::PFParticle> &pfParticle, const Association<recob::PFParticle, recob::Track> &pfpToTracks, const Association<recob::Track, anab::ParticleID> &trackToPIDs, const float chi2ProtonCut) const
{
    // Insist that Pandora has selected this PFParticle as track-like
    if (!this->IsTrackLike(pfParticle))
        return false;

    // Insist that the PFParticle has an associated Track, and this track has a PID object
    try
    {
        const auto track = CollectionHelper::GetSingleAssociated(pfParticle, pfpToTracks);
        const auto pid = CollectionHelper::GetSingleAssociated(track, trackToPIDs);

        // Get the minimum Chi2 under the proton hypothesis using all available views
        const auto minChi2p = this->GetMinChi2Proton(pid);

        return (minChi2p > chi2ProtonCut);
    }
    catch (const cet::exception &)
    {
        return false;
    }
}

// -----------------------------------------------------------------------------------------------------------------------------------------

float EventSelection::GetMinChi2Proton(const art::Ptr<anab::ParticleID> &pid) const
{
    float minChi2Proton = std::numeric_limits<float>::max();
    bool isAvailable = false;
    for (const auto &algo : pid->ParticleIDAlgScores())
    {
        if (algo.fAlgName == "Chi2" && algo.fAssumedPdg == 2212)
        {
            if (algo.fValue >= 0.f && algo.fValue < minChi2Proton)
            {
                minChi2Proton = algo.fValue;
                isAvailable = true;
            }
        }
    }

    if (!isAvailable)
        throw cet::exception("EventSelection::GetMinChi2Proton") << " - Couldn't find any Chi2 under the proton hypothesis for the given PID" << std::endl;

    return minChi2Proton;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

art::Ptr<recob::PFParticle> EventSelection::SelectMuon(const art::Event &event, const PFParticleVector &mips) const
{
    if (mips.empty())
        throw cet::exception("EventSelection::SelectMuon") << " - Input vector of MIP PFParticles is empty" << std::endl;

    const auto pfpToTracks = CollectionHelper::GetAssociation<recob::PFParticle, recob::Track>(event, m_config().PFParticleLabel(), m_config().TrackLabel());
    float maxTrackLength = -std::numeric_limits<float>::max();
    art::Ptr<recob::PFParticle> muon;
    bool foundMuon = false;

    // Find the PFParticle with the longest track
    for (const auto &pfParticle : mips)
    {
        try
        {
            const auto track = CollectionHelper::GetSingleAssociated(pfParticle, pfpToTracks);
            if (track->Length() > maxTrackLength)
            {
                maxTrackLength = track->Length();
                muon = pfParticle;
                foundMuon = true;
            }
        }
        catch (const cet::exception &)
        {
        }
    }

    if (!foundMuon)
        throw cet::exception("EventSelection::SelectMuon") << " - None of the input PFParticles had an associated track" << std::endl;
    
    return muon;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

TVector3 EventSelection::GetRecoNeutrinoVertex(const art::Event &event, const PFParticleVector &allPFParticles) const
{
    const auto pfpToVertex = CollectionHelper::GetAssociation<recob::PFParticle, recob::Vertex>(event, m_config().PFParticleLabel());

    try
    {
        const auto neutrino = RecoHelper::GetNeutrino(allPFParticles);
        const auto vertex = CollectionHelper::GetSingleAssociated(neutrino, pfpToVertex);

        return TVector3(vertex->position().X(), vertex->position().Y(), vertex->position().Z());
    }
    catch (const cet::exception &)
    {
        return TVector3(-std::numeric_limits<float>::max(), -std::numeric_limits<float>::max(), -std::numeric_limits<float>::max());
    }
}

} // namespace ubcc1pi
