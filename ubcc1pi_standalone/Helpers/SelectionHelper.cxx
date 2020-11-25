/**
 *  @file  ubcc1pi_standalone/Helpers/SelectionHelper.cxx
 *
 *  @brief The implementation file for the selection helper class
 */

#include "ubcc1pi_standalone/Helpers/SelectionHelper.h"

namespace ubcc1pi
{
        /**
        class EventSelection
        {
                EventSelection(const std::vector<Cut> &cuts, std::vector<BDTHelper::BDT> &bdts, const SelectionLogic &logic);
                float GetCutValue(const std::string &name) const;
                void SetCutValue(const std::string &name, const float value);
                SelectionResult Execute(const std::shared_ptr<Event> &pEvent);

        };
*/

SelectionHelper::EventSelection::Cut::Cut() :
    m_name(""),
    m_hasValue(false),
    m_value(-std::numeric_limits<float>::max())
{
}

// -----------------------------------------------------------------------------------------------------------------------------------------

SelectionHelper::EventSelection::Cut::Cut(const std::string &name) :
    m_name(name),
    m_hasValue(false),
    m_value(-std::numeric_limits<float>::max())
{
}

// -----------------------------------------------------------------------------------------------------------------------------------------

SelectionHelper::EventSelection::Cut::Cut(const std::string &name, const float value) :
    m_name(name),
    m_hasValue(true),
    m_value(value)
{
}

// -----------------------------------------------------------------------------------------------------------------------------------------

std::string SelectionHelper::EventSelection::Cut::GetName() const
{
    if (m_name.empty())
        throw std::logic_error("Cut::GetName - Name has not been set");

    return m_name;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

float SelectionHelper::EventSelection::Cut::GetValue() const
{
    if (!m_hasValue)
        throw std::logic_error("Cut::GetValue - Can't get value of cut: \"" + m_name + "\" - it has no value");

    return m_value;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void SelectionHelper::EventSelection::Cut::SetValue(const float value)
{
    if (!m_hasValue)
        throw std::logic_error("Cut::SetValue - Can't set value of cut: \"" + m_name + "\" - it has no value");

    m_value = value;
}

// -----------------------------------------------------------------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------------------------------------------------------------

SelectionHelper::EventSelection::CutTracker::CutTracker(const std::vector<Cut> &cuts, const unsigned int nRecoParticles) :
    m_cuts(cuts),
    m_assignedPdgs(nRecoParticles, 0) // ATTN here we use a PDG code of zero to mean "unassigned"
{
}

// -----------------------------------------------------------------------------------------------------------------------------------------

float SelectionHelper::EventSelection::CutTracker::GetCutValue(const std::string &name) const
{
    // Find the cut by name
    const auto iter = std::find_if(m_cuts.begin(), m_cuts.end(), [&](const auto &x) {
        return x.GetName() == name;
    });

    if (iter == m_cuts.end())
        throw std::invalid_argument("CutTracker::GetCutValue - Unknown cut: \"" + name + "\"");

    return iter->GetValue();
}

// -----------------------------------------------------------------------------------------------------------------------------------------

std::vector<std::string> SelectionHelper::EventSelection::CutTracker::GetCutsPassed() const
{
    return m_cutsPassed;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

std::vector<int> SelectionHelper::EventSelection::CutTracker::GetAssignedPDGCodes() const
{
    return m_assignedPdgs;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void SelectionHelper::EventSelection::CutTracker::MarkCutAsPassed(const std::string &name)
{
    // Count the number of cuts that have already passed
    const auto nCutsPassed = m_cutsPassed.size();

    // Check we haven't passed all cuts already
    if (nCutsPassed == m_cuts.size())
        throw std::logic_error("CutTracker::MarkCutAsPassed - All cuts are already marked as passed");

    // Check that the cut we are marking as "passed" matches up with the next cut in the vector
    if (name != m_cuts.at(nCutsPassed).GetName())
        throw std::logic_error("CutTracker::MarkCutAsPassed - The cut \"" + name + "\" is not the next cut in the sequence");

    m_cutsPassed.push_back(name);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void SelectionHelper::EventSelection::CutTracker::AssignPDGCode(const unsigned int recoParticleIndex, const int pdg)
{
    if (recoParticleIndex >= m_assignedPdgs.size())
        throw std::out_of_range("CutTracker::AssignPDGCode - The input recoParticleIndex is out of range");


    // Set the PDG code
    m_assignedPdgs.at(recoParticleIndex) = pdg;
}

// -----------------------------------------------------------------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------------------------------------------------------------

SelectionHelper::EventSelection::EventSelection(const std::vector<Cut> &cuts, BDTMap &bdtMap, const SelectionLogic &logic) :
    m_cuts(cuts),
    m_bdtMap(bdtMap),
    m_logic(logic)
{
    // Make sure none of the input cut names are not repeated
    for (const auto &cut : m_cuts)
    {
        const auto cutName = cut.GetName();

        if (std::count_if(m_cuts.begin(), m_cuts.end(), [&](const auto &x) { return x.GetName() == cutName; }) != 1)
            throw std::invalid_argument("EventSelection::EventSelection - Repeated cut name: \"" + cutName + "\"");
    }
}

// -----------------------------------------------------------------------------------------------------------------------------------------

float SelectionHelper::EventSelection::GetCutValue(const std::string &name) const
{
    // Find the cut by name
    const auto iter = std::find_if(m_cuts.begin(), m_cuts.end(), [&](const auto &x) {
        return x.GetName() == name;
    });

    if (iter == m_cuts.end())
        throw std::invalid_argument("EventSelection::GetCutValue - Unknown cut: \"" + name + "\"");

    return iter->GetValue();
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void SelectionHelper::EventSelection::SetCutValue(const std::string &name, const float value)
{
    // Find the cut by name
    const auto iter = std::find_if(m_cuts.begin(), m_cuts.end(), [&](const auto &x) {
        return x.GetName() == name;
    });

    if (iter == m_cuts.end())
        throw std::invalid_argument("EventSelection::GetCutValue - Unknown cut: \"" + name + "\"");

    return iter->SetValue(value);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

SelectionHelper::EventSelection::SelectionResult SelectionHelper::EventSelection::Execute(const std::shared_ptr<Event> &pEvent)
{
    // Setup the cut tracker
    CutTracker cutTracker(m_cuts, pEvent->reco.particles.size());

    // Pass the event, BDTs and cut tracker through to the selection logic and run the selection
    const auto passed = m_logic(pEvent, m_bdtMap, cutTracker);

    // Return the selection result
    return {
        passed,
        cutTracker.GetCutsPassed(),
        cutTracker.GetAssignedPDGCodes()
    };
}

// -----------------------------------------------------------------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------------------------------------------------------------

unsigned int SelectionHelper::GetMuonCandidateIndex(const std::vector<Event::Reco::Particle> &particles, const std::vector<std::string> &featureNames, BDTHelper::BDT &muonBDT)
{
    bool foundCCInclusiveMuon = false;

    unsigned int ccInclusiveMuonIndex = std::numeric_limits<unsigned int>::max();
    unsigned int muonIndex = std::numeric_limits<unsigned int>::max();

    unsigned int nUncontainedParticles = 0;

    for (unsigned int index = 0; index < particles.size(); ++index)
    {
        const auto &particle = particles.at(index);

        // Check if this particle is the CC inclusive muon candidate
        if (particle.isCCInclusiveMuonCandidate())
        {
            if (foundCCInclusiveMuon)
                throw std::logic_error("SelectionHelper::GetMuonCandidateIndex - found multiple CC inclusive muon candidates");

            foundCCInclusiveMuon = true;
            ccInclusiveMuonIndex = index;
        }

        if (!AnalysisHelper::HasTrackFit(particle))
            continue;

        // For now make the last escaping particle the muon candidate, if there are multiple we will fall back on the CC inclusive candidate
        if (!AnalysisHelper::IsContained(particle))
        {
            muonIndex = index;
            nUncontainedParticles++;
        }
    }

    if (!foundCCInclusiveMuon)
        throw std::logic_error("SelectionHelper::GetMuonCandidateIndex - found no CC inclusive muon candidate");

    if (nUncontainedParticles > 1)
        return ccInclusiveMuonIndex;

    if (nUncontainedParticles == 1)
        return muonIndex;

    // If we are here all particles are contained, choose the muon using the BDT
    float maxMuonBDTResponse = -std::numeric_limits<float>::max();
    bool foundMuon = false;

    for (unsigned int index = 0; index < particles.size(); ++index)
    {
        const auto &particle = particles.at(index);

        if (!AnalysisHelper::HasTrackFit(particle))
            continue;

        if (!AnalysisHelper::IsContained(particle))
            throw std::logic_error("SelectionHelper::GetMuonCandidateIndex - found escaping particle when not expecting to!");

        std::vector<float> features;
        const auto hasFeatures = BDTHelper::GetBDTFeatures(particle, featureNames, features);

        if (!hasFeatures)
            continue;

        const auto muonBDTResponse = muonBDT.GetResponse(features);

        if (muonBDTResponse < maxMuonBDTResponse)
            continue;

        maxMuonBDTResponse = muonBDTResponse;
        muonIndex = index;
        foundMuon = true;
    }

    // If no muon can be found, then default to the CC inclusive candidate
    if (!foundMuon)
        return ccInclusiveMuonIndex;

    return muonIndex;
}

} // namespace ubcc1pi
