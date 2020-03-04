void parseFile()
{
    // -------------------------------------------------------------------------------------------------------------------------------------
    // Setup the input file
    // -------------------------------------------------------------------------------------------------------------------------------------
    TFile *f = new TFile("/uboone/data/users/asmith/ubcc1pi/17022020/eventSelection.root");
    TDirectory *d = static_cast<TDirectory *>(f->Get("eventSelection"));
    TTree *t = static_cast<TTree *>(d->Get("events"));

    const size_t maxEntries = std::numeric_limits<size_t>::max();

    const auto nEntries = static_cast<size_t>(t->GetEntries());
    const auto nEntriesToProcess = std::min(maxEntries, nEntries);
    
    // -------------------------------------------------------------------------------------------------------------------------------------
    // Set up the output root file
    // -------------------------------------------------------------------------------------------------------------------------------------
    TFile *fOut = new TFile("parsedData.root", "RECREATE");
    TTree *tOut = new TTree("particles", "");

    bool isSignalOut;
    tOut->Branch("isSignal", &isSignalOut);
    
    bool isAvailableOut;
    tOut->Branch("isAvailable", &isAvailableOut);
            
    bool isBraggpAvailableOut;
    tOut->Branch("isBraggpAvailable", &isBraggpAvailableOut);

    bool isBraggpiAvailableOut;
    tOut->Branch("isBraggpiAvailable", &isBraggpiAvailableOut);

    bool isBraggMIPAvailableOut;
    tOut->Branch("isBraggMIPAvailable", &isBraggMIPAvailableOut);

    bool isBraggMIPNonZeroOut;
    tOut->Branch("isBraggMIPNonZero", &isBraggMIPNonZeroOut);

    bool isTrackShowerAvailableOut;
    tOut->Branch("isTrackShowerAvailable", &isTrackShowerAvailableOut);

    bool isRMSTrackDeviationAvailableOut;
    tOut->Branch("isRMSTrackDeviationAvailable", &isRMSTrackDeviationAvailableOut);

    bool isRMSSequentialTrackDeviationAvailableOut;
    tOut->Branch("isRMSSequentialTrackDeviationAvailable", &isRMSSequentialTrackDeviationAvailableOut);

    bool isNOtherSpacePointsInSphere5AvailableOut;
    tOut->Branch("isNOtherSpacePointsInSphere5Available", &isNOtherSpacePointsInSphere5AvailableOut);

    bool isNSpacePointsInSphere5AvailableOut;
    tOut->Branch("isNSpacePointsInSphere5Available", &isNSpacePointsInSphere5AvailableOut);
    
    bool isContainedOut;
    tOut->Branch("isContained", &isContainedOut);
    
    bool shouldTrainOut;
    tOut->Branch("shouldTrain", &shouldTrainOut);
    
    bool trueIsGoldenOut;
    tOut->Branch("trueIsGolden", &trueIsGoldenOut);
    
    int truePdgCodeOut;
    tOut->Branch("truePdgCode", &truePdgCodeOut);
    
    float truthMatchCompletenessOut;
    tOut->Branch("truthMatchCompleteness", &truthMatchCompletenessOut);

    // Prefix all features with `f_` for readability
    float logBragg_pToMIPOut;
    tOut->Branch("f_logBragg_pToMIP", &logBragg_pToMIPOut);

    float logBragg_piToMIPOut;
    tOut->Branch("f_logBragg_piToMIP", &logBragg_piToMIPOut);

    float trackShowerOut;
    tOut->Branch("f_trackShower", &trackShowerOut);
    
    float rmsTrackDeviationOut;
    tOut->Branch("f_rmsTrackDeviation", &rmsTrackDeviationOut);
    
    float rmsSequentialTrackDeviationOut;
    tOut->Branch("f_rmsSequentialTrackDeviation", &rmsSequentialTrackDeviationOut);

    float offAxisSpacePointsRMSInSphere5Out;
    tOut->Branch("f_offAxisSpacePointsRMSInSphere5", &offAxisSpacePointsRMSInSphere5Out);
    
    float offAxisSpacePointsRMSInSphere20Out;
    tOut->Branch("f_offAxisSpacePointsRMSInSphere20", &offAxisSpacePointsRMSInSphere20Out);
    
    float offAxisSpacePointsRMSInSphere30Out;
    tOut->Branch("f_offAxisSpacePointsRMSInSphere30", &offAxisSpacePointsRMSInSphere30Out);
    
    int nOtherSpacePointsInSphere5Out;
    tOut->Branch("f_nOtherSpacePointsInSphere5", &nOtherSpacePointsInSphere5Out);
    
    int nSpacePointsInSphere5Out;
    tOut->Branch("f_nSpacePointsInSphere5", &nSpacePointsInSphere5Out);

    int nDownstreamHitsOut;
    tOut->Branch("f_nDownstreamHits", &nDownstreamHitsOut);

    int nDescendentsOut;
    tOut->Branch("f_nDescendents", &nDescendentsOut);

    float dEdxMeanAtStartOut;
    //tOut->Branch("f_dEdxMeanAtStart", &dEdxMeanAtStartOut);

    // -------------------------------------------------------------------------------------------------------------------------------------
    // Set up the input variables
    // -------------------------------------------------------------------------------------------------------------------------------------
    bool isSignal;
    t->SetBranchAddress("isSignal", &isSignal);
    
    int nFinalStatePFPs;
    t->SetBranchAddress("nFinalStatePFPs", &nFinalStatePFPs);

    std::vector<int> *truePdgCodeVect = nullptr;
    t->SetBranchAddress("truePdgCodeVect", &truePdgCodeVect);

    std::vector<float> *truthMatchCompletenessVect = nullptr;
    t->SetBranchAddress("truthMatchCompletenessVect", &truthMatchCompletenessVect);

    std::vector<bool> *trueIsGoldenVect = nullptr;
    t->SetBranchAddress("trueIsGoldenVect", &trueIsGoldenVect);

    std::vector<bool> *isContainedVect = nullptr;
    t->SetBranchAddress("isContainedVect", &isContainedVect);

    std::vector<float> *yzAngleVect = nullptr;
    t->SetBranchAddress("yzAngleVect", &yzAngleVect);

    std::vector<float> *braggpWVect = nullptr;
    t->SetBranchAddress("braggpWVect", &braggpWVect);

    std::vector<float> *braggpUVVect = nullptr;
    t->SetBranchAddress("braggpUVVect", &braggpUVVect);

    std::vector<float> *braggpiWVect = nullptr;
    t->SetBranchAddress("braggpiWVect", &braggpiWVect);

    std::vector<float> *braggpiUVVect = nullptr;
    t->SetBranchAddress("braggpiUVVect", &braggpiUVVect);

    std::vector<float> *braggMIPWVect = nullptr;
    t->SetBranchAddress("braggMIPWVect", &braggMIPWVect);

    std::vector<float> *braggMIPUVVect = nullptr;
    t->SetBranchAddress("braggMIPUVVect", &braggMIPUVVect);

    std::vector<float> *trackShowerVect = nullptr;
    t->SetBranchAddress("trackShowerVect", &trackShowerVect);

    std::vector<float> *rmsTrackDeviationVect = nullptr;
    t->SetBranchAddress("rmsTrackDeviationVect", &rmsTrackDeviationVect);
    
    std::vector<float> *rmsSequentialTrackDeviationVect = nullptr;
    t->SetBranchAddress("rmsSequentialTrackDeviationVect", &rmsSequentialTrackDeviationVect);

    std::vector<float> *offAxisSpacePointsRMSInSphere5Vect = nullptr;
    t->SetBranchAddress("offAxisSpacePointsRMSInSphere5Vect", &offAxisSpacePointsRMSInSphere5Vect);
    
    std::vector<float> *offAxisSpacePointsRMSInSphere20Vect = nullptr;
    t->SetBranchAddress("offAxisSpacePointsRMSInSphere20Vect", &offAxisSpacePointsRMSInSphere20Vect);
    
    std::vector<float> *offAxisSpacePointsRMSInSphere30Vect = nullptr;
    t->SetBranchAddress("offAxisSpacePointsRMSInSphere30Vect", &offAxisSpacePointsRMSInSphere30Vect);
    
    std::vector<int> *nSpacePointsInSphere5Vect = nullptr;
    t->SetBranchAddress("nSpacePointsInSphere5Vect", &nSpacePointsInSphere5Vect);
    
    std::vector<int> *nOtherSpacePointsInSphere5Vect = nullptr;
    t->SetBranchAddress("nOtherSpacePointsInSphere5Vect", &nOtherSpacePointsInSphere5Vect);

    std::vector<int> *nHitsUVect = nullptr;
    t->SetBranchAddress("nHitsUVect", &nHitsUVect);

    std::vector<int> *nHitsVVect = nullptr;
    t->SetBranchAddress("nHitsVVect", &nHitsVVect);

    std::vector<int> *nHitsWVect = nullptr;
    t->SetBranchAddress("nHitsWVect", &nHitsWVect);

    std::vector<int> *nDescendentHitsUVect = nullptr;
    t->SetBranchAddress("nDescendentHitsUVect", &nDescendentHitsUVect);

    std::vector<int> *nDescendentHitsVVect = nullptr;
    t->SetBranchAddress("nDescendentHitsVVect", &nDescendentHitsVVect);

    std::vector<int> *nDescendentHitsWVect = nullptr;
    t->SetBranchAddress("nDescendentHitsWVect", &nDescendentHitsWVect);

    std::vector<int> *nDescendentsVect = nullptr;
    t->SetBranchAddress("nDescendentsVect", &nDescendentsVect);

    std::vector<std::vector<float> > *dedxPerHitWVect = nullptr;
    t->SetBranchAddress("dedxPerHitWVect", &dedxPerHitWVect);
    
    std::vector<std::vector<float> > *residualRangePerHitWVect = nullptr;
    t->SetBranchAddress("residualRangePerHitWVect", &residualRangePerHitWVect);
    
    std::vector<std::vector<float> > *dedxPerHitUVect = nullptr;
    t->SetBranchAddress("dedxPerHitUVect", &dedxPerHitUVect);
    
    std::vector<std::vector<float> > *residualRangePerHitUVect = nullptr;
    t->SetBranchAddress("residualRangePerHitUVect", &residualRangePerHitUVect);
    
    std::vector<std::vector<float> > *dedxPerHitVVect = nullptr;
    t->SetBranchAddress("dedxPerHitVVect", &dedxPerHitVVect);
    
    std::vector<std::vector<float> > *residualRangePerHitVVect = nullptr;
    t->SetBranchAddress("residualRangePerHitVVect", &residualRangePerHitVVect);

    // -------------------------------------------------------------------------------------------------------------------------------------
    // Extract the relevant info
    // -------------------------------------------------------------------------------------------------------------------------------------
    const auto printInterval = static_cast<int>(std::ceil(static_cast<float>(nEntriesToProcess) / 100.f));
    const auto startTime = std::time(nullptr);

    for (size_t iEvent = 0; iEvent < nEntriesToProcess; ++iEvent)
    {
        t->GetEntry(iEvent);

        if (iEvent % printInterval == 0 && iEvent != 0)
        {
            const auto timeNow = std::time(nullptr);
            const auto fracComplete = static_cast<float>(iEvent) / static_cast<float>(nEntriesToProcess);

            const auto secondsElapsed = timeNow - startTime;
            const auto secondsPerEvent = static_cast<float>(secondsElapsed) / static_cast<float>(iEvent);
            const auto eventsRemaining = nEntriesToProcess - iEvent;
            const auto secondsTotalRemaining = static_cast<float>(eventsRemaining * secondsPerEvent);

            const auto minutesRemaining = static_cast<int>(std::floor(secondsTotalRemaining / 60.f));
            const auto secondsRemaining = static_cast<int>(std::floor(secondsTotalRemaining - 60 * minutesRemaining));

            std::cout << "Event " << iEvent << " of " << nEntriesToProcess << " (" << (100.f * fracComplete) << "%). ";
            std::cout << minutesRemaining << "m " << secondsRemaining << "s remaining" << std::endl;
        }

        // Only look for signal events
        if (!isSignal)
            continue;
        
        // Loop over the final state PFParticles
        for (size_t iPFP = 0; iPFP < nFinalStatePFPs; ++iPFP)
        {
            // Get the PFP from the vectors
            const auto truePdgCode = truePdgCodeVect->at(iPFP);
            const auto truthMatchCompleteness = truthMatchCompletenessVect->at(iPFP);
            const auto trueIsGolden = trueIsGoldenVect->at(iPFP);
            const auto isContained = isContainedVect->at(iPFP);
            const auto yzAngle = yzAngleVect->at(iPFP);
            const auto braggpW = braggpWVect->at(iPFP);
            const auto braggpUV = braggpUVVect->at(iPFP);
            const auto braggpiW = braggpiWVect->at(iPFP);
            const auto braggpiUV = braggpiUVVect->at(iPFP);
            const auto braggMIPW = braggMIPWVect->at(iPFP);
            const auto braggMIPUV = braggMIPUVVect->at(iPFP);
            const auto trackShower = trackShowerVect->at(iPFP);
            const auto rmsTrackDeviation = rmsTrackDeviationVect->at(iPFP);
            const auto rmsSequentialTrackDeviation = rmsSequentialTrackDeviationVect->at(iPFP);
            auto offAxisSpacePointsRMSInSphere5 = offAxisSpacePointsRMSInSphere5Vect->at(iPFP);
            auto offAxisSpacePointsRMSInSphere20 = offAxisSpacePointsRMSInSphere20Vect->at(iPFP);
            auto offAxisSpacePointsRMSInSphere30 = offAxisSpacePointsRMSInSphere30Vect->at(iPFP);
            const auto nHitsU = nHitsUVect->at(iPFP);
            const auto nHitsV = nHitsVVect->at(iPFP);
            const auto nHitsW = nHitsWVect->at(iPFP);
            const auto nDescendentHitsU = nDescendentHitsUVect->at(iPFP);
            const auto nDescendentHitsV = nDescendentHitsVVect->at(iPFP);
            const auto nDescendentHitsW = nDescendentHitsWVect->at(iPFP);
            const auto nOtherSpacePointsInSphere5 = nOtherSpacePointsInSphere5Vect->at(iPFP);
            const auto nSpacePointsInSphere5 = nSpacePointsInSphere5Vect->at(iPFP);
            const auto nDescendents = nDescendentsVect->at(iPFP);
            const auto dedxPerHitW = dedxPerHitWVect->at(iPFP);
            const auto residualRangePerHitW = residualRangePerHitWVect->at(iPFP);
            const auto dedxPerHitU = dedxPerHitUVect->at(iPFP);
            const auto residualRangePerHitU = residualRangePerHitUVect->at(iPFP);
            const auto dedxPerHitV = dedxPerHitVVect->at(iPFP);
            const auto residualRangePerHitV = residualRangePerHitVVect->at(iPFP);
            // END_OF_VARIABLES

            // Golden pions are our signal particle, everything else is background
            const auto isSignalPFP = (truePdgCode == 211 && trueIsGolden);

            const auto isTrackAlongWWire = (std::pow(std::sin(yzAngle), 2) < 0.175);

            /* BEGIN TEST */
            const auto nHitsUsed = dedxPerHitW.size();

            auto maxResidualRange = -std::numeric_limits<float>::max();
            for (size_t iHit = 0; iHit < nHitsUsed; ++iHit)
            {
                const auto residualRange = residualRangePerHitW.at(iHit);
    
                if (residualRange > maxResidualRange)
                    maxResidualRange = residualRange;   
            }
            
            const auto residualRangeThreshold = 0.5f * maxResidualRange;
            float dEdxSum = 0.f;
            int nHitsAboveThreshold = 0;
            for (size_t iHit = 0; iHit < nHitsUsed; ++iHit)
            {
                const auto dEdx = dedxPerHitW.at(iHit);
                const auto residualRange = residualRangePerHitW.at(iHit);
    
                if (residualRange > residualRangeThreshold)
                {
                    dEdxSum += dEdx;
                    nHitsAboveThreshold++;
                }
            }

            auto dEdxMeanAtStart = -std::numeric_limits<float>::max();
            if (nHitsAboveThreshold > 0)
            {
                dEdxMeanAtStart = dEdxSum / static_cast<float>(nHitsAboveThreshold);
            }

            /* END TEST */

            // Check availablity of braggp variable
            const auto isBraggpWAvailable = (braggpW > -1 && !isTrackAlongWWire);
            const auto isBraggpUVAvailable = (braggpUV > -1);
            const auto isBraggpAvailable = isBraggpWAvailable || isBraggpUVAvailable;
            const auto braggp = isBraggpAvailable ? (isBraggpWAvailable ? braggpW : braggpUV) : -std::numeric_limits<float>::max();
            
            // Check availablity of braggpi variable
            const auto isBraggpiWAvailable = (braggpiW > -1 && !isTrackAlongWWire);
            const auto isBraggpiUVAvailable = (braggpiUV > -1);
            const auto isBraggpiAvailable = isBraggpiWAvailable || isBraggpiUVAvailable;
            const auto braggpi = isBraggpiAvailable ? (isBraggpiWAvailable ? braggpiW : braggpiUV) : -std::numeric_limits<float>::max();
            
            // Check availablity of braggMIP variable
            const auto isBraggMIPWAvailable = (braggMIPW > -1 && !isTrackAlongWWire);
            const auto isBraggMIPUVAvailable = (braggMIPUV > -1);
            const auto isBraggMIPAvailable = isBraggMIPWAvailable || isBraggMIPUVAvailable;
            const auto braggMIP = isBraggMIPAvailable ? (isBraggMIPWAvailable ? braggMIPW : braggMIPUV) : -std::numeric_limits<float>::max();

            // Calculate the ratio variables
            const auto isBraggMIPZero = std::abs(braggMIP) < std::numeric_limits<float>::epsilon();

            const auto logBragg_pToMIP = (isBraggpAvailable && isBraggMIPAvailable && !isBraggMIPZero) ? std::log(braggp / braggMIP) : (-std::numeric_limits<float>::max());
            const auto logBragg_piToMIP = (isBraggpiAvailable && isBraggMIPAvailable && !isBraggMIPZero) ? std::log(braggpi / braggMIP) : (-std::numeric_limits<float>::max());
        
            const auto nDownstreamHits = (nDescendentHitsU + nDescendentHitsV + nDescendentHitsW) - (nHitsU + nHitsV + nHitsW);

            // Check the availability of the variables
            const auto isTrackShowerAvailable = (trackShower > -0.5);
            const auto isRMSTrackDeviationAvailable = (rmsTrackDeviation > -1);
            const auto isRMSSequentialTrackDeviationAvailable = (rmsSequentialTrackDeviation > -1);
            const auto isNOtherSpacePointsInSphere5Available = (nOtherSpacePointsInSphere5 > -1);
            const auto isNSpacePointsInSphere5Available = (nSpacePointsInSphere5 > -1);

            if (!isNOtherSpacePointsInSphere5Available || nOtherSpacePointsInSphere5 == 0)
            {
                offAxisSpacePointsRMSInSphere5 = 0.f;
                offAxisSpacePointsRMSInSphere20 = 0.f;
                offAxisSpacePointsRMSInSphere20 = 0.f;
            }

            auto isAvailable = true;
            isAvailable = isAvailable && isBraggpAvailable;
            isAvailable = isAvailable && isBraggpiAvailable;
            isAvailable = isAvailable && isBraggMIPAvailable;
            isAvailable = isAvailable && !isBraggMIPZero;
            isAvailable = isAvailable && isTrackShowerAvailable;
            isAvailable = isAvailable && isRMSTrackDeviationAvailable;
            isAvailable = isAvailable && isRMSSequentialTrackDeviationAvailable;
            isAvailable = isAvailable && isNOtherSpacePointsInSphere5Available;
            isAvailable = isAvailable && isNSpacePointsInSphere5Available;

            // Only train on mostly complete, contained particles that have all of their features available
            const auto shouldTrain = (truthMatchCompleteness > 0.5f && isAvailable && isContained);

            // Set the output variables (NB. I copy with the "Out" suffix here just to help with readability)
            isSignalOut                        = isSignalPFP;
            
            isBraggpAvailableOut                          = isBraggpAvailable;
            isBraggpiAvailableOut                         = isBraggpiAvailable;
            isBraggMIPAvailableOut                        = isBraggMIPAvailable;
            isBraggMIPNonZeroOut                          = !isBraggMIPZero;
            isTrackShowerAvailableOut                     = isTrackShowerAvailable;
            isRMSTrackDeviationAvailableOut               = isRMSTrackDeviationAvailable;
            isRMSSequentialTrackDeviationAvailableOut     = isRMSSequentialTrackDeviationAvailable;
            isNOtherSpacePointsInSphere5AvailableOut      = isNOtherSpacePointsInSphere5Available;
            isNSpacePointsInSphere5AvailableOut           = isNSpacePointsInSphere5Available;

            isAvailableOut                     = isAvailable;
            shouldTrainOut                     = shouldTrain;
            isContainedOut                     = isContained;
            trueIsGoldenOut                    = trueIsGolden;
            truePdgCodeOut                     = truePdgCode;
            truthMatchCompletenessOut          = truthMatchCompleteness;
            logBragg_pToMIPOut                 = logBragg_pToMIP;
            logBragg_piToMIPOut                = logBragg_piToMIP;
            trackShowerOut                     = trackShower;
            rmsTrackDeviationOut               = rmsTrackDeviation;
            rmsSequentialTrackDeviationOut     = rmsSequentialTrackDeviation;
            offAxisSpacePointsRMSInSphere5Out  = offAxisSpacePointsRMSInSphere5;
            offAxisSpacePointsRMSInSphere20Out = offAxisSpacePointsRMSInSphere20;
            offAxisSpacePointsRMSInSphere30Out = offAxisSpacePointsRMSInSphere30;
            nOtherSpacePointsInSphere5Out      = nOtherSpacePointsInSphere5;
            nSpacePointsInSphere5Out           = nSpacePointsInSphere5;
            nDownstreamHitsOut                 = nDownstreamHits;
            nDescendentsOut                    = nDescendents;

            dEdxMeanAtStartOut                 = dEdxMeanAtStart; //// TEST

            tOut->Fill();
        }
    }

    // Save the output file
    fOut->Write();
    fOut->Close();
}
