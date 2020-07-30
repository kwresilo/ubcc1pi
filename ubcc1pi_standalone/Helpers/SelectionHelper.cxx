/**
 *  @file  ubcc1pi_standalone/Helpers/SelectionHelper.cxx
 *
 *  @brief The implementation file for the selection helper class
 */

#include "ubcc1pi_standalone/Helpers/SelectionHelper.h"
#include "ubcc1pi_standalone/Helpers/FormattingHelper.h"
#include "ubcc1pi_standalone/Objects/FileReader.h"

#include <stdexcept>
#include <algorithm>

namespace ubcc1pi
{
                        
SelectionHelper::EventSelection::CutManager::CutManager() :
    m_nScanPoints(2u)
{
}

// -----------------------------------------------------------------------------------------------------------------------------------------

bool SelectionHelper::EventSelection::CutManager::HasCut(const std::string &name) const
{
    return (std::find_if(m_cuts.begin(), m_cuts.end(), [&](const Cut &cut) {return cut.m_name == name;}) != m_cuts.end());
}

// -----------------------------------------------------------------------------------------------------------------------------------------
                        
SelectionHelper::EventSelection::CutManager::Cut& SelectionHelper::EventSelection::CutManager::GetCut(const std::string &name)
{
    auto iter = std::find_if(m_cuts.begin(), m_cuts.end(), [&](const Cut &cut) {return cut.m_name == name;});
    if (iter == m_cuts.end())
        throw std::invalid_argument("CutManager::GetCut - Unknown cut: \"" + name + "\"");

    return *iter;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

const SelectionHelper::EventSelection::CutManager::Cut& SelectionHelper::EventSelection::CutManager::GetCut(const std::string &name) const
{
    const auto iter = std::find_if(m_cuts.begin(), m_cuts.end(), [&](const Cut &cut) {return cut.m_name == name;});
    if (iter == m_cuts.end())
        throw std::invalid_argument("CutManager::GetCut - Unknown cut: \"" + name + "\"");

    return *iter;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

bool SelectionHelper::EventSelection::CutManager::GetCutResult(const Cut &cut, const std::function<bool()> &method) const
{
    if (cut.m_hasValue)
        throw std::invalid_argument("CutManager::GetCutResult - The cut: \"" + cut.m_name + "\" takes a value, non supplied.");

    // Automatically pass disabled cuts
    if (!cut.m_isEnabled)
        return true;

    // Get the result of the cut
    return method();
}

// -----------------------------------------------------------------------------------------------------------------------------------------

bool SelectionHelper::EventSelection::CutManager::GetCutResult(const Cut &cut, const std::function<bool(const float &)> &method) const
{
    if (!cut.m_hasValue)
        throw std::invalid_argument("CutManager::GetCutResult - The cut: \"" + cut.m_name + "\" doesn't take a value");

    // Automatically pass disabled cuts
    if (!cut.m_isEnabled)
        return true;

    // Get the result of the cut
    return method(cut.m_value);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

bool SelectionHelper::EventSelection::CutManager::GetCutResult(const std::string &name, const std::function<bool()> &method)
{
    // Make a copy of the cut to be modified
    const auto cutRef = this->GetCut(name);
    auto cut = cutRef;

    // If we aren't optimizing, then just run the cut
    if (cut.m_name != m_cutOptimizing)
    {
        const auto passes = this->GetCutResult(cut, method);

        if (passes)
            m_defaultEventCounter.CountEvent(name, m_sampleType, m_pEvent, m_weight);
            
        return passes;
    }
    
    // Store the event counts with the cut disabled
    m_disabledEventCounter.CountEvent(name, m_sampleType, m_pEvent, m_weight);
    
    // Store the event counts with cut enabled
    cut.m_isEnabled = true;
    if (this->GetCutResult(cut, method))
        m_enabledEventCounters.front().CountEvent(name, m_sampleType, m_pEvent, m_weight);

    // Always stop at the cut we are currently optimizing
    return false;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

bool SelectionHelper::EventSelection::CutManager::GetCutResult(const std::string &name, const std::function<bool(const float &)> &method)
{
    // Make a copy of the cut to be modified
    const auto cutRef = this->GetCut(name);
    auto cut = cutRef;
    
    // If we aren't optimizing, then just run the cut
    if (name != m_cutOptimizing)
    {
        const auto passes = this->GetCutResult(cut, method);

        if (passes)
            m_defaultEventCounter.CountEvent(name, m_sampleType, m_pEvent, m_weight);

        return passes;
    }

    // Store the event counts with the cut disabled
    m_disabledEventCounter.CountEvent(name, m_sampleType, m_pEvent, m_weight);

    // Store the event counts with the cut enabled for various cut values
    cut.m_isEnabled = true;
    for (unsigned int i = 0; i < m_enabledEventCounters.size(); ++i)
    {
        cut.m_value = m_values.at(i);

        if (this->GetCutResult(cut, method))
            m_enabledEventCounters.at(i).CountEvent(name, m_sampleType, m_pEvent, m_weight);
    }

    // Always stop at the cut we are currently optimizing
    // This is important to ensure that any changes made to variables captured by `method` are always at the correct cut value
    return false;
}

// -----------------------------------------------------------------------------------------------------------------------------------------
                        
void SelectionHelper::EventSelection::CutManager::SetParticlePdg(const unsigned int recoParticleIndex, const int pdgCode)
{
    if (!m_pEvent)
        throw std::logic_error("CutManager::SetParticlePdg - No event currently assigned");

    const auto nRecoParticles = m_pEvent->reco.particles.size();
    if (m_assignedPdgCodes.size() != nRecoParticles)
        throw std::invalid_argument("CutManager::SetParticlePdg - The stored reco particle PDGs are out of sync with the loaded event");

    if (recoParticleIndex > nRecoParticles)
        throw std::invalid_argument("CutManager::SetParticlePdg - reco particle index is out of bounds");

    m_assignedPdgCodes.at(recoParticleIndex) = pdgCode;
}

// -----------------------------------------------------------------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------------------------------------------------------------
                        
void SelectionHelper::EventSelection::BDTManager::Add(const std::string &name, const std::vector<std::string> &featureNames)
{
    for (const auto &pBdt : m_bdts)
    {
        if (pBdt->GetName() == name)
            throw std::invalid_argument("BDTManager::Add - repeated BDT name: \"" + name + "\"");
    }

    m_bdts.emplace_back(std::make_shared<BDTHelper::BDT>(name, featureNames));
}

// -----------------------------------------------------------------------------------------------------------------------------------------

BDTHelper::BDT & SelectionHelper::EventSelection::BDTManager::Get(const std::string &name)
{
    for (auto &pBdt : m_bdts)
    {
        if (pBdt->GetName() == name)
            return *pBdt;
    }

    throw std::invalid_argument("BDTManager::Get - no BDT found with name: \"" + name + "\"");
}

// -----------------------------------------------------------------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------------------------------------------------------------
                
void SelectionHelper::EventSelection::DeclareCut(const std::string &name)
{
    if (m_cutManager.HasCut(name))
        throw std::invalid_argument("EventSelection::DeclareCut - repeated cut name: \"" + name + "\"");

    if (name.empty())
        throw std::invalid_argument("EventSelection::DeclareCut - empty cut name");

    CutManager::Cut cut;
    cut.m_name = name;
    cut.m_canDisable = false;
    cut.m_hasValue = false;
    cut.m_nominal = -std::numeric_limits<float>::max();
    cut.m_min = -std::numeric_limits<float>::max();
    cut.m_max = -std::numeric_limits<float>::max();
    cut.m_shouldOptimize = false;
    cut.m_searchQuery = "";
    cut.m_isEnabled = true;
    cut.m_value = -std::numeric_limits<float>::max();

    m_cutManager.m_cuts.push_back(cut);
}

// -----------------------------------------------------------------------------------------------------------------------------------------
                
void SelectionHelper::EventSelection::DeclareCut(const std::string &name, const float &nominal)
{
    if (m_cutManager.HasCut(name))
        throw std::invalid_argument("EventSelection::DeclareCut - repeated cut name: \"" + name + "\"");
    
    if (name.empty())
        throw std::invalid_argument("EventSelection::DeclareCut - empty cut name");

    CutManager::Cut cut;
    cut.m_name = name;
    cut.m_canDisable = false;
    cut.m_hasValue = true;
    cut.m_nominal = nominal;
    cut.m_min = -std::numeric_limits<float>::max();
    cut.m_max = -std::numeric_limits<float>::max();
    cut.m_shouldOptimize = false;
    cut.m_searchQuery = "";
    cut.m_isEnabled = true;
    cut.m_value = nominal;

    m_cutManager.m_cuts.push_back(cut);
}

// -----------------------------------------------------------------------------------------------------------------------------------------
                
void SelectionHelper::EventSelection::SetCutNominalValue(const std::string &name, const float &nominal)
{
    if (!m_cutManager.HasCut(name))
        throw std::invalid_argument("EventSelection::SetCutNominalValue - unknown cut name: \"" + name + "\"");
    
    auto &cut = m_cutManager.GetCut(name);
    cut.m_nominal = nominal;
    cut.m_value = nominal;
}

// -----------------------------------------------------------------------------------------------------------------------------------------
                
void SelectionHelper::EventSelection::EnableOptimization(const std::string &name, const std::string &searchQuery)
{
    auto &cut = m_cutManager.GetCut(name);

    if (cut.m_hasValue)
        throw std::invalid_argument("EventSelection::EnableOptimization - no min or max supplied for optimization");
    
    cut.m_canDisable = true;
    cut.m_shouldOptimize = true;
    cut.m_searchQuery = searchQuery;
}

// -----------------------------------------------------------------------------------------------------------------------------------------
                
void SelectionHelper::EventSelection::EnableOptimization(const std::string &name, const bool canDisable, const float &min, const float &max, const std::string &searchQuery)
{
    auto &cut = m_cutManager.GetCut(name);

    if (!cut.m_hasValue)
        throw std::invalid_argument("EventSelection::EnableOptimization - can set min and max of cut without no value");
    
    cut.m_canDisable = canDisable;
    cut.m_min = min;
    cut.m_max = max;
    cut.m_shouldOptimize = true;
    cut.m_searchQuery = searchQuery;
}

// -----------------------------------------------------------------------------------------------------------------------------------------
                
void SelectionHelper::EventSelection::AssignBDT(const std::string &bdtName, const std::vector<std::string> &featureNames)
{
    m_bdtManager.Add(bdtName, featureNames);
}

// -----------------------------------------------------------------------------------------------------------------------------------------
                
std::vector<std::string> SelectionHelper::EventSelection::GetCuts() const
{
    std::vector<std::string> cutNames;
    for (const auto &cut : m_cutManager.m_cuts)
        cutNames.push_back(cut.m_name);

    return cutNames;
}

// -----------------------------------------------------------------------------------------------------------------------------------------
                
float SelectionHelper::EventSelection::GetCutNominalValue(const std::string &name) const
{
    if (!m_cutManager.HasCut(name))
        throw std::invalid_argument("EventSelection::GetCutNominalValue - unknown cut name: \"" + name + "\"");
    
    const auto cut = m_cutManager.GetCut(name);
    return cut.m_nominal;
}

// -----------------------------------------------------------------------------------------------------------------------------------------
                
void SelectionHelper::EventSelection::DefineSelectionMethod(const SelectionMethod &method)
{
    m_selectionMethod = method;
}

// -----------------------------------------------------------------------------------------------------------------------------------------
                
void SelectionHelper::EventSelection::Optimize(const std::string &dataBNBFileName, const std::string &overlayFileName,
    const float overlayWeight, const std::string &dataEXTFileName, const float dataEXTWeight, const std::string &dirtFileName,
    const float dirtWeight, const unsigned int nScanPoints, const float processFraction)
{
    if (nScanPoints <= 1)
        throw std::invalid_argument("EventSelection::Optimize - Need more than 1 scanning point");

    m_cutManager.m_nScanPoints = nScanPoints;

    // Optimize each cut separately
    for (auto &cut : m_cutManager.m_cuts)
    {
        if (!cut.m_shouldOptimize)
            continue;
        
        m_cutManager.m_cutOptimizing = cut.m_name;
        FormattingHelper::PrintLine();
        std::cout << "Optimizing cut: " << m_cutManager.m_cutOptimizing << " for signal events with query: " << cut.m_searchQuery << std::endl;
        FormattingHelper::PrintLine();

        // Reset the event counters and set the values to scan if required
        cut.m_isEnabled = true;
        m_cutManager.m_disabledEventCounter = AnalysisHelper::EventCounter();
        m_cutManager.m_enabledEventCounters.clear();
        m_cutManager.m_values.clear();

        if (cut.m_hasValue)
        {
            // Add the nominal value
            m_cutManager.m_enabledEventCounters.push_back(AnalysisHelper::EventCounter());
            m_cutManager.m_values.push_back(cut.m_nominal);

            for (unsigned int i = 0; i < m_cutManager.m_nScanPoints; ++i)
            {
                m_cutManager.m_enabledEventCounters.push_back(AnalysisHelper::EventCounter());
    
                const float value = cut.m_min + (cut.m_max - cut.m_min) * (static_cast<float>(i) / static_cast<float>(m_cutManager.m_nScanPoints - 1));
                m_cutManager.m_values.push_back(value);
            }
        }
        else
        {
            // For simple on-off cuts we only need one counter
            m_cutManager.m_enabledEventCounters.push_back(AnalysisHelper::EventCounter()); 
        }

        // Run the event selection and store the results
        this->Execute(dataBNBFileName, overlayFileName, overlayWeight, dataEXTFileName, dataEXTWeight, dirtFileName, dirtWeight, false, processFraction);

        // Print the classifications we are optimizing for
        const auto signalClassifications = m_cutManager.m_disabledEventCounter.GetSignalClassifications(cut.m_searchQuery);
        std::cout << "Optimizing for signal classifications:";
        for (const auto &classification : signalClassifications)
            std::cout << "  " << classification;
        std::cout << std::endl;

        // Find the best cut value
        FormattingHelper::Table table({"Is Enabled", "Value", "", "Efficiency", "Purity", "E*P"});

        float bestScore = -std::numeric_limits<float>::max();
        for (unsigned int iCounter = 0; iCounter < m_cutManager.m_enabledEventCounters.size(); ++ iCounter)
        {
            const auto &counter = m_cutManager.m_enabledEventCounters.at(iCounter);

            const auto efficiency = counter.GetSignalEfficiency(cut.m_name, cut.m_searchQuery);
            const auto purity = counter.GetSignalPurity(cut.m_name, cut.m_searchQuery);
            const auto score = (purity < 0 || efficiency < 0) ? -std::numeric_limits<float>::max() : (efficiency * purity);

            table.AddEmptyRow();
            table.SetEntry("Is Enabled", true);

            if (cut.m_hasValue)
                table.SetEntry("Value", m_cutManager.m_values.at(iCounter));

            table.SetEntry("Efficiency", efficiency);
            table.SetEntry("Purity", purity);
            table.SetEntry("E*P", score);
            
            // If this is the best yet then save it
            if (score > bestScore)
            {
                cut.m_isEnabled = true;
                bestScore = score;

                if (cut.m_hasValue)
                    cut.m_value = m_cutManager.m_values.at(iCounter);
            }
        }
        
        // If we can disable, then check if that's better
        if (cut.m_canDisable)
        {
            const auto disabledEfficiency = m_cutManager.m_disabledEventCounter.GetSignalEfficiency(cut.m_name, cut.m_searchQuery);
            const auto disabledPurity = m_cutManager.m_disabledEventCounter.GetSignalPurity(cut.m_name, cut.m_searchQuery);
            const auto disabledScore = disabledEfficiency * disabledPurity;

            table.AddEmptyRow();
            table.SetEntry("Is Enabled", false);
            table.SetEntry("Efficiency", disabledEfficiency);
            table.SetEntry("Purity", disabledPurity);
            table.SetEntry("E*P", disabledScore);

            if (disabledScore > bestScore)
            {
                cut.m_isEnabled = false;
                cut.m_value = -std::numeric_limits<float>::max();
                bestScore = disabledScore;
            }
        }

        table.Print();

        FormattingHelper::PrintLine();
        std::cout << "Optimization result:" << std::endl;
        std::cout << "  - Best E*P:   " << bestScore << std::endl;
        std::cout << "  - Enable cut: " << cut.m_isEnabled << std::endl;
        if (cut.m_isEnabled && cut.m_hasValue)
            std::cout << "  - Cut value:  " << cut.m_value << std::endl;
        FormattingHelper::PrintLine();
    }

    // Done optimizing
    m_cutManager.m_cutOptimizing.clear();
    m_cutManager.m_enabledEventCounters.clear();
    m_cutManager.m_values.clear();
}

// -----------------------------------------------------------------------------------------------------------------------------------------

bool SelectionHelper::EventSelection::Execute(const std::shared_ptr<Event> &pEvent, std::vector<std::string> &cutsPassed, std::vector<int> &assignedPdgCodes)
{
    if (!m_selectionMethod)
        throw std::logic_error("EventSelection::Execute - Selection method hasn't been set");
   
    if (!cutsPassed.empty())
        throw std::invalid_argument("EventSelection::Execute - Input vector of cut names isn't empty");
    
    if (!assignedPdgCodes.empty())
        throw std::invalid_argument("EventSelection::Execute - Input vector of assigned PDG codes isn't empty");

    m_cutManager.m_pEvent = pEvent;

    // Here we use the event counter machinary to see if the event passed, it's not very pretty but it works. The sample type and weight
    // don't actually matter as we check to see if *any* weight passes each cut... Could be done better!
    m_cutManager.m_defaultEventCounter = AnalysisHelper::EventCounter(); // This is probably super slow... consider a better way
    m_cutManager.m_sampleType = AnalysisHelper::Overlay;
    m_cutManager.m_weight = 1.f;

    m_cutManager.m_assignedPdgCodes = std::vector<int>(m_cutManager.m_pEvent->reco.particles.size(), -std::numeric_limits<int>::max());

    const auto passed = m_selectionMethod(m_cutManager.m_pEvent, m_bdtManager, m_cutManager);
    
    for (const auto &tag : m_cutManager.m_defaultEventCounter.GetTags())
    {
        const auto mcWeight = m_cutManager.m_defaultEventCounter.GetTotalMCWeight(tag);
        const auto bnbDataWeight = m_cutManager.m_defaultEventCounter.GetBNBDataWeight(tag);
        
        // ATTN this is not a nice way of doing this!!
        // Here we check if any weight has been added to the counter for this tag
        const auto hasMCWeight = std::abs(m_cutManager.m_defaultEventCounter.GetTotalMCWeight(tag)) > std::numeric_limits<float>::epsilon();
        const auto hasBNBDataWeight = std::abs(m_cutManager.m_defaultEventCounter.GetBNBDataWeight(tag)) > std::numeric_limits<float>::epsilon();
        
        if (hasMCWeight && hasBNBDataWeight)
            throw std::logic_error("EventSelection::Execute - Input event was counted as BNB data and MC!");

        if (hasMCWeight || hasBNBDataWeight)
            cutsPassed.push_back(tag);
    }

    assignedPdgCodes = m_cutManager.m_assignedPdgCodes;

    return passed;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void SelectionHelper::EventSelection::Execute(const std::string &dataBNBFileName, const std::string &overlayFileName,
    const float overlayWeight, const std::string &dataEXTFileName, const float dataEXTWeight, const std::string &dirtFileName,
    const float dirtWeight, const bool shouldPrint, const float processFraction, const unsigned int nEntriesToPrint)
{
    if (!m_selectionMethod)
        throw std::logic_error("EventSelection::Execute - Selection method hasn't been set");

    // Reset the default event counter
    m_cutManager.m_defaultEventCounter = AnalysisHelper::EventCounter();

    for (const auto sampleType : AnalysisHelper::AllSampleTypes)
    {
        // Store the current sample type
        m_cutManager.m_sampleType = sampleType;

        // Get the appropriate file name and normalisation
        std::string fileName = "";
        float normalisation = 0.f;
        switch (sampleType)
        {
            case AnalysisHelper::DataBNB:
                fileName = dataBNBFileName;
                normalisation = 1.f;
                break;
            case AnalysisHelper::Overlay:
                fileName = overlayFileName;
                normalisation = overlayWeight;
                break;
            case AnalysisHelper::DataEXT:
                fileName = dataEXTFileName;
                normalisation = dataEXTWeight;
                break;
            case AnalysisHelper::Dirt:
                fileName = dirtFileName;
                normalisation = dirtWeight;
                break;
            default:
                throw std::logic_error("EventSelection::Execute - Unknown sample type");
        }

        // Don't attempt to process files that aren't set
        if (fileName.empty())
            continue;

        // Read the file
        std::cout << "Processing file - " << fileName << std::endl;
        FileReader reader(fileName);

        // Store the current event
        m_cutManager.m_pEvent = reader.GetBoundEventAddress();

        // Run the selection
        const auto nEvents = reader.GetNumberOfEvents();
        const auto nEventsToProcess = std::min(nEvents, static_cast<unsigned int>(std::ceil(static_cast<float>(nEvents) * std::max(0.f, processFraction))));
        for (unsigned int i = 0; i < nEventsToProcess; ++i)
        {
            AnalysisHelper::PrintLoadingBar(i, nEventsToProcess);

            reader.LoadEvent(i);

            // Set the event weight
            m_cutManager.m_weight = normalisation * AnalysisHelper::GetNominalEventWeight(m_cutManager.m_pEvent);

            // Mark the special "all" cut for overlays
            if (sampleType == AnalysisHelper::Overlay)
            {
                m_cutManager.m_defaultEventCounter.CountEvent("all", m_cutManager.m_sampleType, m_cutManager.m_pEvent, m_cutManager.m_weight);
                m_cutManager.m_disabledEventCounter.CountEvent("all", m_cutManager.m_sampleType, m_cutManager.m_pEvent, m_cutManager.m_weight);

                for (auto &counter : m_cutManager.m_enabledEventCounters)
                    counter.CountEvent("all", m_cutManager.m_sampleType, m_cutManager.m_pEvent, m_cutManager.m_weight);
            }

            m_cutManager.m_assignedPdgCodes = std::vector<int>(m_cutManager.m_pEvent->reco.particles.size(), -std::numeric_limits<int>::max());
            m_selectionMethod(m_cutManager.m_pEvent, m_bdtManager, m_cutManager);
        }
    }

    // Print the result
    if (!shouldPrint)
        return;

    FormattingHelper::PrintLine();
    std::cout << "Cuts" << std::endl;
    FormattingHelper::PrintLine();
    FormattingHelper::Table table({"Cut", "", "Is Enabled", "Optimized value", "Nominal value"});
    for (const auto &cut : m_cutManager.m_cuts)
    {
        table.AddEmptyRow();
        table.SetEntry("Cut", cut.m_name);
        table.SetEntry("Is Enabled", cut.m_isEnabled);

        if (cut.m_hasValue)
        {
            table.SetEntry("Optimized value", cut.m_value);
            table.SetEntry("Nominal value", cut.m_nominal);
        }
    }
    table.WriteToFile("eventSelection_cuts.md");

    FormattingHelper::PrintLine();
    std::cout << "Summary" << std::endl;
    FormattingHelper::PrintLine();
    m_cutManager.m_defaultEventCounter.PrintBreakdownSummary("eventSelection_summary.md");
    
    FormattingHelper::PrintLine();
    std::cout << "Details" << std::endl;
    FormattingHelper::PrintLine();
    m_cutManager.m_defaultEventCounter.PrintBreakdownDetails("eventSelection_details.md", nEntriesToPrint);
}

// -----------------------------------------------------------------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------------------------------------------------------------

SelectionHelper::EventSelection SelectionHelper::GetCCInclusiveSelection()
{
    EventSelection selection;
    
    // Set up the cuts
    selection.DeclareCut("passesCCInclusive");
    
    // Define the selection
    selection.DefineSelectionMethod([](const std::shared_ptr<Event> &pEvent, EventSelection::BDTManager &bdtManager, EventSelection::CutManager &cuts) {
        
        // Insist the event passes the CC inclusive selection
        if (!cuts.GetCutResult("passesCCInclusive", [&](){
            return pEvent->reco.passesCCInclusive();

        })) return false;

        // Identify the muon
        const auto &recoParticles = pEvent->reco.particles;
        for (unsigned int index = 0; index < recoParticles.size(); ++index)
        {
            const auto &particle = recoParticles.at(index);

            if (particle.isCCInclusiveMuonCandidate())
            {
                cuts.SetParticlePdg(index, 13);
            }
        }
        
        return true;
    });

    return selection;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

SelectionHelper::EventSelection SelectionHelper::GetDefaultSelection()
{
    // Set up the selection
    EventSelection selection;
    
    // Set up the cuts
    selection.DeclareCut("passesCCInclusive");
    selection.DeclareCut("min2Tracks");
    selection.DeclareCut("max1Uncontained");
    selection.DeclareCut("2NonProtons", -0.06f);
//    selection.DeclareCut("fakePionRejection", -0.5f);
    selection.DeclareCut("openingAngle", 2.65f);
    selection.DeclareCut("topologicalScore", 0.67f);
    selection.DeclareCut("startNearVertex", 9.5f);
//    selection.DeclareCut("pionRange", 5.f);
//    selection.DeclareCut("muonRange", 20.f);
    selection.DeclareCut("likelyGoldenPion", -0.03f);
    
    // Get the BDT feature names
    const auto muonFeatureNames = BDTHelper::MuonBDTFeatureNames;
    const auto protonFeatureNames = BDTHelper::ProtonBDTFeatureNames;
    const auto goldenPionFeatureNames = BDTHelper::GoldenPionBDTFeatureNames;
    selection.AssignBDT("muon", muonFeatureNames);
    selection.AssignBDT("proton", protonFeatureNames);
    selection.AssignBDT("goldenPion", goldenPionFeatureNames);
    
    // Define the selection
    selection.DefineSelectionMethod([](const std::shared_ptr<Event> &pEvent, EventSelection::BDTManager &bdtManager, EventSelection::CutManager &cuts) {

        // Get the BDTs owned by the event selection object
        const auto muonFeatureNames = BDTHelper::MuonBDTFeatureNames;
        const auto protonFeatureNames = BDTHelper::ProtonBDTFeatureNames;
        const auto goldenPionFeatureNames = BDTHelper::GoldenPionBDTFeatureNames;

        auto &muonBDT = bdtManager.Get("muon");
        auto &protonBDT = bdtManager.Get("proton");
        auto &goldenPionBDT = bdtManager.Get("goldenPion");

        // Make a map from reconstructed particle index to the assigned PDG code
        std::unordered_map<unsigned int, int> pdgCodeMap;

        // -------------------------------------------------------------
        // Insist the event passes the CC inclusive selection
        // -------------------------------------------------------------
        if (!cuts.GetCutResult("passesCCInclusive", [&](){
            return pEvent->reco.passesCCInclusive();

        })) return false;
        
        // -------------------------------------------------------------
        // Find the particles with a track fit and check if they are contained
        // -------------------------------------------------------------
        const auto &recoParticles = pEvent->reco.particles;
        unsigned int nTrackParticles = 0u;
        unsigned int nUncontainedParticles = 0u;
        for (unsigned int index = 0; index < recoParticles.size(); ++index)
        {
            const auto &particle = recoParticles.at(index);

            if (AnalysisHelper::HasTrackFit(particle))
            {
                // Here we use a PDG code of 0 to mean "as yet unknown"
                pdgCodeMap[index] = 0;

                nTrackParticles++;

                if (!AnalysisHelper::IsContained(particle))
                    nUncontainedParticles++;
            }
            else
            {
                // Assume particles without tracks are just small protons
                pdgCodeMap[index] = 2212;
            }
        }
        
        // -------------------------------------------------------------
        // Insist at least 2 particles have a track fit
        // -------------------------------------------------------------
        if (!cuts.GetCutResult("min2Tracks", [&](){
            return (nTrackParticles >= 2);

        })) return false;
        
        // -------------------------------------------------------------
        // Insist at most one particle is uncontained
        // -------------------------------------------------------------
        if (!cuts.GetCutResult("max1Uncontained", [&](){
            return (nUncontainedParticles <= 1);
        })) return false;

        // -------------------------------------------------------------
        // Find the muon candidate 
        // -------------------------------------------------------------
        const auto muonIndex = SelectionHelper::GetMuonCandidateIndex(recoParticles, muonFeatureNames, muonBDT);
        pdgCodeMap.at(muonIndex) = 13;
      
        // -------------------------------------------------------------
        // Cache the proton BDT responses of the remaining particles
        // -------------------------------------------------------------
        std::unordered_map<unsigned int, bool> protonBDTHasResponseCache;
        std::unordered_map<unsigned int, float> protonBDTResponseCache;
        for (unsigned int index = 0; index < recoParticles.size(); ++index)
        {
            // Skip particles we have already identified
            if (pdgCodeMap.at(index) != 0)
                continue;

            const auto &particle = recoParticles.at(index);
            
            if (!AnalysisHelper::IsContained(particle))
                throw std::logic_error("MakeSelectionTable - Found an uncontained particle that hasn't been identified already!");

            std::vector<float> features;
            const auto hasFeatures = BDTHelper::GetBDTFeatures(particle, protonFeatureNames, features);
            protonBDTHasResponseCache[index] = hasFeatures;

            if (!hasFeatures)
                continue;

            const auto protonBDTResponse = protonBDT.GetResponse(features);
            protonBDTResponseCache[index] = protonBDTResponse;
        }
        
        // -------------------------------------------------------------
        // Identify the protons, and insist there are exactly 2 non-protons
        // -------------------------------------------------------------
        std::vector<unsigned int> protonIndices;

        if (!cuts.GetCutResult("2NonProtons", [&](const float &cut) {
            protonIndices.clear(); // Important! This lambda can be evaluated multiple times during optimization, and variable is captured by reference - need to reset

            unsigned int nonProtons = 0u;

            for (unsigned int index = 0; index < recoParticles.size(); ++index)
            {
                // Count the protons we have found already
                if (pdgCodeMap.at(index) == 2212)
                    protonIndices.push_back(index);

                // Skip particles we have already identified
                if (pdgCodeMap.at(index) != 0)
                    continue;
                
                const auto &particle = recoParticles.at(index);

                // Assume the particles without a BDT response are protons
                if (!protonBDTHasResponseCache.at(index))
                {
                    protonIndices.push_back(index);
                    continue;
                }

                // Apply the cut on the proton BDT response
                if (protonBDTResponseCache.at(index) >= cut)
                    protonIndices.push_back(index);
            }

            // Insist we have exactly 2 non-protons
            const auto nNonProtons = static_cast<int>(recoParticles.size()) - static_cast<int>(protonIndices.size());
            return (nNonProtons == 2);

        })) return false;
        
        
        // -------------------------------------------------------------
        // Set the remaining particle types
        // -------------------------------------------------------------
        bool foundPion = false;
        unsigned int pionIndex = -std::numeric_limits<unsigned int>::max();
            
        for (unsigned int index = 0; index < recoParticles.size(); ++index)
        {
            // Skip particles we have already identified
            if (pdgCodeMap.at(index) != 0)
                continue;

            const auto isProton = (std::find(protonIndices.begin(), protonIndices.end(), index) != protonIndices.end());
            if (isProton)
            {
                pdgCodeMap.at(index) = 2212;
                continue;
            }

            // Anything left is our pion
            if (foundPion)
                throw std::logic_error("MakeSelectionTable - Found multiple pion candidates!");

            pdgCodeMap.at(index) = 211;
            pionIndex = index;
            foundPion = true;
        }

        if (!foundPion)
            throw std::logic_error("MakeSelectionTable - Couldn't find pion candidate!");

        // Set the particle types so they are available outside of this function
        for (unsigned int index = 0; index < recoParticles.size(); ++index)
            cuts.SetParticlePdg(index, pdgCodeMap.at(index));
        
        // -------------------------------------------------------------
        // Sanity check the PID
        // -------------------------------------------------------------
        unsigned int nMuonsIdentified = 0u;
        unsigned int nPionsIdentified = 0u;
        unsigned int nProtonsIdentified = 0u;

        for (unsigned int index = 0; index < recoParticles.size(); ++index)
        {
            switch (pdgCodeMap.at(index))
            {
                case 13:
                    nMuonsIdentified++;
                    break;
                case 211:
                    nPionsIdentified++;
                    break;
                case 2212:
                    nProtonsIdentified++;
                    break;
                default:
                    throw std::logic_error("MakeSelectionTable - Sanity check failed, found reco particle assigned a PDG of: " + std::to_string(pdgCodeMap.at(index)));
            }
        }

        if (nMuonsIdentified != 1)
            throw std::logic_error("MakeSelectionTable - Found " + std::to_string(nMuonsIdentified) + " muon candidates!");

        if (nPionsIdentified != 1)
            throw std::logic_error("MakeSelectionTable - Found " + std::to_string(nPionsIdentified) + " pion candidates!");

        if (nMuonsIdentified + nPionsIdentified + nProtonsIdentified != recoParticles.size())
            throw std::logic_error("MakeSelectionTable - Not all particles have been identified");
        
        const auto &muon = recoParticles.at(muonIndex);
        const auto &pion = recoParticles.at(pionIndex);

        // -------------------------------------------------------------
        // Insist the muon-pion opening angle isn't too wide
        // -------------------------------------------------------------
        const auto muonDir = TVector3(muon.directionX(), muon.directionY(), muon.directionZ()).Unit();
        const auto pionDir = TVector3(pion.directionX(), pion.directionY(), pion.directionZ()).Unit();
        const auto openingAngle = muonDir.Angle(pionDir);
        
        if (!cuts.GetCutResult("openingAngle", [&](const float &cut){
            return (openingAngle < cut);

        })) return false;
        
        // -------------------------------------------------------------
        // Cut on the topological score to reduce cosmic backgrounds
        // -------------------------------------------------------------
        const auto topologicalScore = pEvent->reco.selectedTopologicalScore();
        if (!cuts.GetCutResult("topologicalScore", [&](const float &cut){
            return (topologicalScore > cut);

        })) return false;

        // -------------------------------------------------------------
        // Insist all particles start near the vertex
        // -------------------------------------------------------------
        const auto recoVertex = pEvent->reco.nuVertex();
        if (!cuts.GetCutResult("startNearVertex", [&](const float &cut){

            const auto cut2 = cut * cut;
            for (unsigned int index = 0; index < recoParticles.size(); ++index)
            {
                const auto &particle = recoParticles.at(index);

                if (!AnalysisHelper::HasTrackFit(particle))
                    continue;

                const TVector3 start(particle.startX(), particle.startY(), particle.startZ());
                const auto vertexDist2 = (start - recoVertex).Mag2();

                if (vertexDist2 > cut2)
                    return false;
            }

            return true;

        })) return false;
       
        // -------------------------------------------------------------
        // Insist the pion is sufficiently long
        // -------------------------------------------------------------
        /*
        if (!cuts.GetCutResult("pionRange", [&](const float &cut){
            return pion.range() >= cut;

        })) return false;
        */

        // -------------------------------------------------------------
        // Insist the muon is sufficiently long
        // -------------------------------------------------------------
        /*
        if (!cuts.GetCutResult("muonRange", [&](const float &cut){
            return muon.range() >= cut;

        })) return false;
        */

        // -------------------------------------------------------------
        // Get the golden pion BDT response of the pion candidate
        // -------------------------------------------------------------
        std::vector<float> features;
        if (!BDTHelper::GetBDTFeatures(pion, goldenPionFeatureNames, features))
            throw std::logic_error("MakeSelectionTable - can't get golden pion BDT features for pion candidate");

        const auto goldenPionBDTResponse = goldenPionBDT.GetResponse(features);
        
        // -------------------------------------------------------------
        // Insist this is a likely golden pion
        // -------------------------------------------------------------
        if (!cuts.GetCutResult("likelyGoldenPion", [&](const float &cut){
            return (goldenPionBDTResponse > cut);

        })) return false;

        return true;
    });

    return selection;
}

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
        
    if (!foundMuon)
        return ccInclusiveMuonIndex;
//        throw std::logic_error("SelectionHelper::GetMuonCandidateIndex - all particles are contained, but couldn't find a muon candidate using the BDT");

    return muonIndex;
}

} // namespace ubcc1pi
