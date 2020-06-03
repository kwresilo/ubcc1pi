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

SelectionHelper::EventSelection SelectionHelper::GetDefaultSelection()
{
    // Set up the selection
    EventSelection selection;
    
    // Set up the cuts
    selection.DeclareCut("passesCCInclusive");
    selection.DeclareCut("min2Tracks");
    selection.DeclareCut("max1Uncontained");
    selection.DeclareCut("startNearVertex", 5.f);
    selection.DeclareCut("2NonProtons", -0.06f);
    selection.DeclareCut("nonProtonLength", 5.f);
    selection.DeclareCut("openingAngle", 2.65f);
    selection.DeclareCut("topologicalScore", 0.67f);
    selection.DeclareCut("noShowers", 0.04f);
    selection.DeclareCut("likelyGoldenPion", -0.03f);
    
    // Get the BDT feature names
    const auto featureNames = BDTHelper::ParticleBDTFeatureNames;
    selection.AssignBDT("proton", featureNames);
    selection.AssignBDT("goldenPion", featureNames);
    
    // Define the selection
    selection.DefineSelectionMethod([](const std::shared_ptr<Event> &pEvent, EventSelection::BDTManager &bdtManager, EventSelection::CutManager &cuts) {

        // Get the BDTs owned by the event selection object
        const auto featureNames = BDTHelper::ParticleBDTFeatureNames;
        auto &protonBDT = bdtManager.Get("proton");
        auto &goldenPionBDT = bdtManager.Get("goldenPion");

        // Insist the event passes the CC inclusive selection
        if (!cuts.GetCutResult("passesCCInclusive", [&](){
            return pEvent->reco.passesCCInclusive();

        })) return false;
        
        
        // Find the particles with a track
        const auto &recoParticles = pEvent->reco.particles;
        std::vector<unsigned int> trackParticleIndices;
        for (unsigned int index = 0; index < recoParticles.size(); ++index)
        {
            const auto &particle = recoParticles.at(index);

            if (AnalysisHelper::HasTrackFit(particle))
            {
                trackParticleIndices.push_back(index);
            }
            else
            {
                // Assume particles without tracks are just small protons
                cuts.SetParticlePdg(index, 2212);
            }
        }
        
        // Insist at least 2 track-like particles
        if (!cuts.GetCutResult("min2Tracks", [&](){
            return (trackParticleIndices.size() >= 2);

        })) return false;
        
        
        // Get the contained and uncontained particles
        std::vector<unsigned int> containedParticleIndices, uncontainedParticleIndices;
        for (const auto &index : trackParticleIndices)
        {
            const auto &particle = recoParticles.at(index);

            auto &indexVector = (AnalysisHelper::IsContained(particle) ? containedParticleIndices : uncontainedParticleIndices);
            indexVector.push_back(index);
        }

        // Insist at most one particle is uncontained
        if (!cuts.GetCutResult("max1Uncontained", [&](){
            return (uncontainedParticleIndices.size() <= 1);

        })) return false;
        
        // Just a sanity check in case someone tries to disable an above cut 
        if (uncontainedParticleIndices.size() > 1)
            throw std::logic_error("MakeSelectionTable - more than 1 uncontained particle");


        // If there is an uncontained then identify it as the muon
        if (!uncontainedParticleIndices.empty())
            cuts.SetParticlePdg(uncontainedParticleIndices.front(), 13);

        // If there's only one contained particles, identify it as the pion
        if (containedParticleIndices.size() == 1)
            cuts.SetParticlePdg(containedParticleIndices.front(), 211);
       
        // Insist all particles start near the vertex
        const auto recoVertex = pEvent->reco.nuVertex();
        if (!cuts.GetCutResult("startNearVertex", [&](const float &cut){

            const auto cut2 = cut * cut;
            for (const auto &index : trackParticleIndices)
            {
                const auto &particle = recoParticles.at(index);

                const TVector3 start(particle.startX(), particle.startY(), particle.startZ());
                const auto vertexDist2 = (start - recoVertex).Mag2();

                if (vertexDist2 > cut2)
                    return false;
            }

            return true;

        })) return false;

        // Get the BDT responses of the contained particles
        std::vector<bool> hasFeaturesVect;
        std::vector<float> protonBDTResponseVect;
            
        for (unsigned int index = 0; index < recoParticles.size(); ++index)
        {
            const auto &particle = recoParticles.at(index);

            std::vector<float> features;
            bool hasFeatures = false;
            float protonBDTResponse = -std::numeric_limits<float>::max();
           
            // Only use contained particles
            if (std::find(containedParticleIndices.begin(), containedParticleIndices.end(), index) != containedParticleIndices.end())
            {
                hasFeatures = BDTHelper::GetBDTFeatures(particle, featureNames, features);

                if (hasFeatures)
                    protonBDTResponse = protonBDT.GetResponse(features);
            }

            hasFeaturesVect.push_back(hasFeatures);
            protonBDTResponseVect.push_back(protonBDTResponse);
        }

        // Insist that there are at least two particles that are non-protons
        std::vector<unsigned int> containedNonProtonIndices;
        if (!cuts.GetCutResult("2NonProtons", [&](const float &cut){
    
            containedNonProtonIndices.clear(); // Important! This lambda gets evaluated multiple times, and variable is captured by reference - need to reset

            for (const auto &index : containedParticleIndices)
            {
                if (!hasFeaturesVect.at(index))
                    continue;

                const auto protonBDTResponse = protonBDTResponseVect.at(index);
                const auto &particle = recoParticles.at(index);

                if (protonBDTResponse < cut)
                    containedNonProtonIndices.push_back(index);
            }

            return (containedNonProtonIndices.size() + uncontainedParticleIndices.size() == 2);

        })) return false;

        // Assign the PDG codes of the protons
        for (const auto &index : containedParticleIndices)
        {
            if (std::find(containedNonProtonIndices.begin(), containedNonProtonIndices.end(), index) == containedNonProtonIndices.end())
                cuts.SetParticlePdg(index, 2212);
        }

        // If possible assign the PDG code of the pion
        if (containedNonProtonIndices.size() == 1)
            cuts.SetParticlePdg(containedNonProtonIndices.front(), 211);
        
        
        // Combine the non-protons
        auto nonProtonIndices = containedNonProtonIndices;
        nonProtonIndices.insert(nonProtonIndices.end(), uncontainedParticleIndices.begin(), uncontainedParticleIndices.end());

        // Just a sanity check in case someone tries to disable an above cut 
        if (nonProtonIndices.size() != 2)
            throw std::logic_error("MakeSelectionTable - more than 2 non-protons");
        
        // Another sanity check
        if (containedNonProtonIndices.empty())
            throw std::logic_error("MakeSelectionTable - no contained non-protons");
        
        if (containedNonProtonIndices.size() > 2)
            throw std::logic_error("MakeSelectionTable - more than 2 contained non-protons");
    
       
        // Insist the non protons are sufficiency long
        if (!cuts.GetCutResult("nonProtonLength", [&](const float &cut){
            for (const auto &index : nonProtonIndices)
            {
                const auto &particle = recoParticles.at(index);

                if (particle.length() < cut)
                    return false;
            }
            return true;

        })) return false;
        
        // For reproducibility order by length
        const auto frontParticle = recoParticles.at(nonProtonIndices.front());
        const auto backParticle = recoParticles.at(nonProtonIndices.back());

        const auto frontLength = frontParticle.length();
        const auto backLength = backParticle.length();
        const bool isLongerFront = frontLength > backLength;

        const auto longParticle = isLongerFront ? frontParticle : backParticle;
        const auto shortParticle = isLongerFront ? backParticle : frontParticle;

        const auto longDir = TVector3(longParticle.directionX(), longParticle.directionY(), longParticle.directionZ()).Unit();
        const auto shortDir = TVector3(shortParticle.directionX(), shortParticle.directionY(), shortParticle.directionZ()).Unit();
        const auto openingAngle = longDir.Angle(shortDir);
        
        if (!cuts.GetCutResult("openingAngle", [&](const float &cut){
            return (openingAngle < cut);

        })) return false;
        
        
        // Cut on the topological score to reduce cosmic backgrounds
        const auto topologicalScore = pEvent->reco.selectedTopologicalScore();
        if (!cuts.GetCutResult("topologicalScore", [&](const float &cut){
            return (topologicalScore > cut);

        })) return false;
        
        
        // Insist there are no shower like particles
        if (!cuts.GetCutResult("noShowers", [&](const float &cut){
            for (unsigned int index = 0; index < recoParticles.size(); ++index)
            {
                const auto &particle = recoParticles.at(index);

                if (!particle.trackScore.IsSet())
                    continue;

                if (particle.trackScore() < cut)
                    return false;
            }

            return true;

        })) return false;
        
        
        // Get the most likely golden pion
        float goldenPionBDTResponse = -std::numeric_limits<float>::max();
        unsigned int goldenPionIndex = std::numeric_limits<unsigned int>::max();

        for (const auto &index : containedNonProtonIndices)
        {
            const auto &particle = recoParticles.at(index);

            std::vector<float> features;
            if (!BDTHelper::GetBDTFeatures(particle, featureNames, features))
                throw std::logic_error("MakeSelectionTable - can't get features for a contained non-proton");

            const auto response = goldenPionBDT.GetResponse(features);
            if (response > goldenPionBDTResponse)
            {
                goldenPionBDTResponse = response;
                goldenPionIndex = index;
            }
        }

        // Set the PDG codes
        cuts.SetParticlePdg(goldenPionIndex, 211);
        for (const auto &index : containedNonProtonIndices)
        {
            if (index == goldenPionIndex)
                continue;
            
            cuts.SetParticlePdg(index, 13);
        }

        // Insist this is a likely golden pion
        if (!cuts.GetCutResult("likelyGoldenPion", [&](const float &cut){
            return (goldenPionBDTResponse > cut);

        })) return false;

        return true;
    });

    return selection;
}

} // namespace ubcc1pi
