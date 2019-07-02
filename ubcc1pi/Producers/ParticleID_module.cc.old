/**
 *  @file  ubcc1pi/Producers/ParticleID_module.cc
 *
 *  @brief The implementation file for the particle id producer.
 */

#include "ubcc1pi/Producers/ParticleID.h"

#include "ubcc1pi/Helpers/CollectionHelper.h"
#include "ubcc1pi/Helpers/RecoHelper.h"
#include "ubcc1pi/Helpers/AnalysisHelper.h"

namespace ubcc1pi
{

ParticleID::ParticleID(const art::EDProducer::Table<Config> &config) :
    art::EDProducer(config),
    m_config(config)
{
    // Setup the output trees
    art::ServiceHandle<art::TFileService> fileService;
    
    m_pPointTree = fileService->make<TTree>("particles", "");
    m_pPointTree->Branch("run", &m_outputPoint.m_run);
    m_pPointTree->Branch("subRun", &m_outputPoint.m_subRun);
    m_pPointTree->Branch("event", &m_outputPoint.m_event);
    m_pPointTree->Branch("pfParticleId", &m_outputPoint.m_pfParticleId);
    m_pPointTree->Branch("hasMatchedMCParticle", &m_outputPoint.m_hasMatchedMCParticle);
    m_pPointTree->Branch("truePdgCode", &m_outputPoint.m_truePdgCode);
    m_pPointTree->Branch("trueMomentum", &m_outputPoint.m_trueMomentum);
    m_pPointTree->Branch("trueMatchPurity", &m_outputPoint.m_trueMatchPurity);
    m_pPointTree->Branch("trueMatchCompleteness", &m_outputPoint.m_trueMatchCompleteness);
    m_pPointTree->Branch("trackTheta", &m_outputPoint.m_trackTheta);
    m_pPointTree->Branch("trackPhi", &m_outputPoint.m_trackPhi);
    m_pPointTree->Branch("view", &m_outputPoint.m_view);
    m_pPointTree->Branch("residualRange", &m_outputPoint.m_residualRange);
    m_pPointTree->Branch("dEdx", &m_outputPoint.m_dEdx);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void ParticleID::produce(art::Event &event)
{
    const auto mcTruthLabel = m_config().MCTruthLabel();
    const auto mcParticleLabel = m_config().MCParticleLabel();
    const auto backtrackerLabel = m_config().BacktrackerLabel();

    const auto pfParticleLabel = m_config().PFParticleLabel();
    const auto trackLabel = m_config().TrackLabel();
    const auto calorimetryLabel = m_config().CalorimetryLabel();

    const float kernalWidth = 0.1;
    const float maxStepSize = 0.05;
    const float braggPeakLength = 30.f;

    m_outputPoint.m_run = event.run();
    m_outputPoint.m_subRun = event.subRun();
    m_outputPoint.m_event = event.event();

    std::cout << "PFParticleLabel = " << pfParticleLabel << std::endl;
    std::cout << "TrackLabel = " << trackLabel << std::endl;
    std::cout << "CaloLabel = " << calorimetryLabel << std::endl;

    const auto allPFParticles = CollectionHelper::GetCollection<recob::PFParticle>(event, pfParticleLabel);
    const auto finalStatePFParticles = RecoHelper::GetNeutrinoFinalStates(allPFParticles);

    const auto pfpToTrack = CollectionHelper::GetAssociation<recob::PFParticle, recob::Track>(event, pfParticleLabel, trackLabel);
    const auto trackToCalorimetry = CollectionHelper::GetAssociation<recob::Track, anab::Calorimetry>(event, trackLabel, calorimetryLabel);
    
    const auto backtrackerData = AnalysisHelper::GetBacktrackerData(event, mcTruthLabel, mcParticleLabel, backtrackerLabel, pfParticleLabel);

    for (const auto &pfParticle : finalStatePFParticles)
    {
        m_outputPoint.m_pfParticleId = pfParticle->Self();
        this->SetMatchedMCParticleInfo(pfParticle, backtrackerData);

        std::cout << "-------------------------------------------------------" << std::endl;
        std::cout << "PFParticle " << pfParticle->Self() << std::endl;
        std::cout << "  - pdg code = " << (m_outputPoint.m_hasMatchedMCParticle ? m_outputPoint.m_truePdgCode : -999999) << std::endl;
        std::cout << "-------------------------------------------------------" << std::endl;

        art::Ptr<recob::Track> track;
        CalorimetryVector calos;

        try
        {
            track = CollectionHelper::GetSingleAssociated(pfParticle, pfpToTrack);    
            calos = CollectionHelper::GetManyAssociated(track, trackToCalorimetry);
        }
        catch (const cet::exception &ex)
        {
            continue;
        }

        // Output the info
        m_outputPoint.m_trackTheta = track->Theta();
        m_outputPoint.m_trackPhi = track->Phi();
        this->OutputPoints(calos);

        // Test out some stuff
        std::cout << "Track length: " << track->Length() << std::endl;
        const auto points = this->GetCalorimetryPoints(calos);
            
        if (points.empty())
            continue;

        const auto maxResidualRange = points.back().first;
        const auto braggPeakRange = std::min(maxResidualRange, braggPeakLength);

        std::cout << "BraggPeakRange = " << braggPeakRange << std::endl;
        const auto pointsA = this->GetPointsInRange(points, 0.f, braggPeakRange * 0.5f);
        const auto pointsB = this->GetPointsInRange(points, braggPeakRange * 0.5f, braggPeakRange);

        auto dEdxRatio = -std::numeric_limits<float>::max();
        if (!pointsA.empty() && !pointsB.empty())
            dEdxRatio = this->GetMeandEdx(pointsA) / this->GetMeandEdx(pointsB);

        auto truncatedMeandEdx = -std::numeric_limits<float>::max();
        auto truncatedMPVdEdx = -std::numeric_limits<float>::max();

        const auto pointsC = this->GetPointsInRange(points, braggPeakRange * 0.5f, maxResidualRange);
        if (!pointsC.empty())
        {
            truncatedMeandEdx = this->GetMeandEdx(pointsC);
            truncatedMPVdEdx = this->GetMostProbabledEdx(pointsC, kernalWidth, maxStepSize);
        }

        const auto meandEdx = this->GetMeandEdx(points);

        // Count points with dEdx > 2.5
        unsigned int nPointsAboveCut = 0;
        for (const auto &point : points)
        {
            if (point.second > 2.5)
                nPointsAboveCut++;
        }
        const float protonFrac = static_cast<float>(nPointsAboveCut) / static_cast<float>(points.size());

        std::cout << "Total points        = " << points.size() << std::endl;
        std::cout << " - In tip           = " << pointsA.size() << std::endl;
        std::cout << " - In end, not tip  = " << pointsB.size() << std::endl;
        std::cout << " - Not in tip       = " << pointsC.size() << std::endl;

        std::cout << "Mean dEdx           = " << meandEdx << std::endl;
        std::cout << "dEdx ratio          = " << dEdxRatio << std::endl;
        std::cout << "Truncated mean dEdx = " << truncatedMeandEdx << std::endl;
        std::cout << "Truncated MPV dEdx  = " << truncatedMPVdEdx << std::endl;
        std::cout << "protonFrac          = " << protonFrac << std::endl;
    }
}

// -----------------------------------------------------------------------------------------------------------------------------------------

std::vector<std::pair<float, float> > ParticleID::GetCalorimetryPoints(const CalorimetryVector &calos) const
{
    std::vector<std::pair<float, float> > points;

    for (const auto &calo : calos)
    {
        auto pointsInView = this->GetCalorimetryPoints(calo);
        points.insert(points.end(), pointsInView.begin(), pointsInView.end());
    }

    this->SortPoints(points);
    return points;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

std::vector<std::pair<float, float> > ParticleID::GetCalorimetryPoints(const art::Ptr<anab::Calorimetry> &calo) const
{
    std::vector<std::pair<float, float> > points;

    const auto nPoints = calo->dEdx().size();
    if (nPoints != calo->ResidualRange().size())
        throw cet::exception("ParticleID::GetCalorimetryPoints") << " - Number of dEdx points doesn't match number of residual ranges!" << std::endl;

    // Get the dEdx and residual ranges
    for (unsigned int i = 0; i < nPoints; ++i)
    {
        // ATTN this cut may be a bad idea
        if (calo->dEdx().at(i) > 0.5)
            points.emplace_back(calo->ResidualRange().at(i), calo->dEdx().at(i));
    }

    // Sort the points by residual range
    this->SortPoints(points);

    return points;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void ParticleID::SortPoints(std::vector< std::pair<float, float> > &points) const
{
    std::sort(points.begin(), points.end(), [](const std::pair<float, float> &pointA, const std::pair<float, float> &pointB){
        // Sort by residual range
        if (std::abs(pointA.first - pointB.first) > std::numeric_limits<float>::epsilon())
            return pointA.first < pointB.first;

        // If we have two points at the same range, then sort by dEdx
        return pointA.second < pointB.second;
    });
}

// -----------------------------------------------------------------------------------------------------------------------------------------

std::vector<std::pair<float, float> > ParticleID::GetPointsInRange(const std::vector< std::pair<float, float> > &points, const float min, const float max) const
{
    std::vector<std::pair<float, float> > endPoints;

    for (const auto &point : points)
    {
        if (point.first > min && point.first <= max)
            endPoints.push_back(point);
    }

    return endPoints;
} 

// -----------------------------------------------------------------------------------------------------------------------------------------

float ParticleID::GetMeandEdx(const std::vector<std::pair<float, float> > &points) const
{
    if (points.empty())
        throw cet::exception("ParticleID::GetMeandEdX") << " - The input vector of calorimetry points is empty!" << std::endl;

    float total = 0.f;
    for (const auto &point : points)
        total += point.second;

    return total / static_cast<float>(points.size());
}

// -----------------------------------------------------------------------------------------------------------------------------------------

float ParticleID::GetMostProbabledEdx(const std::vector<std::pair<float, float> > &points, const float kernalWidth, const float maxStepSize) const
{
    if (kernalWidth <= std::numeric_limits<float>::epsilon())
        throw cet::exception("ParticleID::GetMostProbabledEdx") << " - Invalid kernal width: " << kernalWidth << std::endl;

    if (maxStepSize <= std::numeric_limits<float>::epsilon())
        throw cet::exception("ParticleID::GetMostProbabledEdx") << " - Invalid maximum step size: " << maxStepSize << std::endl;

    if (points.empty())
        throw cet::exception("ParticleID::GetMostProbabledEdx") << " - The input vector of calorimetry points is empty!" << std::endl;

    // If there is only one point, then this is the most probable
    if (points.size() == 1)
        return points.front().second;

    // Get the maxiumum and the minimum dEdxs
    float max = -std::numeric_limits<float>::max();
    float min = std::numeric_limits<float>::max();

    for (const auto &point : points)
    {
        const auto dEdx = point.second;

        if (dEdx > max)
            max = dEdx;

        if (dEdx < min)
            min = dEdx;
    } 

    const auto range = max - min;

    // Determine the number of steps we will to use
    const auto nSteps = std::floor(range / maxStepSize);

    // If the range of the data is smaller than the step size, then just return something in the middle
    if (static_cast<int>(nSteps) == 0)
        return (max + min) * 0.5f;
    
    // Use Gaussian kernal density estimation to find the most probable value
    const float normalization = 1.f / (std::pow(4 * std::asin(1), 2) * kernalWidth * points.size());
    float maxDensity = -std::numeric_limits<float>::max();
    float mostProbableValue = -std::numeric_limits<float>::max();
    
    const auto stepSize = range / nSteps;
    for (unsigned int i = 0; i <= nSteps; ++i)
    {
        const auto dEdxSample = min + i * stepSize;
        float density = 0.f;

        for (const auto &point : points)
        {
            // Get the contribution from this point using a gaussian kernal centered at the sample dEdx
            const auto dEdx = point.second; 
            density += std::exp(std::pow(dEdx - dEdxSample, 2) / (-2.f * kernalWidth * kernalWidth)) * normalization;
        }

        if (density < maxDensity)
            continue;

        maxDensity = density;
        mostProbableValue = dEdxSample;
    }
    
    return mostProbableValue;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void ParticleID::SetMatchedMCParticleInfo(const art::Ptr<recob::PFParticle> &pfParticle, const BacktrackHelper::BacktrackerData &backtrackerData)
{
    try
    {
        const auto mcParticle = backtrackerData.GetBestMatchedMCParticle(pfParticle);

        m_outputPoint.m_truePdgCode = mcParticle->PdgCode();
        m_outputPoint.m_trueMomentum = mcParticle->P();
        m_outputPoint.m_trueMatchPurity = backtrackerData.GetMatchPurity(pfParticle, mcParticle);
        m_outputPoint.m_trueMatchCompleteness = backtrackerData.GetMatchCompleteness(pfParticle, mcParticle);
        m_outputPoint.m_hasMatchedMCParticle = true;
    }
    catch (const cet::exception &)
    {
        m_outputPoint.m_truePdgCode = -std::numeric_limits<int>::max();
        m_outputPoint.m_trueMomentum = -std::numeric_limits<float>::max();
        m_outputPoint.m_trueMatchPurity = -std::numeric_limits<float>::max();
        m_outputPoint.m_trueMatchCompleteness = -std::numeric_limits<float>::max();
        m_outputPoint.m_hasMatchedMCParticle = false;
    }
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void ParticleID::OutputPoints(const CalorimetryVector &calos)
{
    for (const auto &calo : calos)
    {
        m_outputPoint.m_view = calo->PlaneID().deepestIndex();

        for (const auto &point : this->GetCalorimetryPoints(calo))
        {
            m_outputPoint.m_residualRange = point.first;
            m_outputPoint.m_dEdx = point.second;

            m_pPointTree->Fill();
        }
    }
}

} // namespace ubcc1pi
