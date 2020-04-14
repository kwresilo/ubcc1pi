 #include <cstdlib>
 #include <iostream>
 #include <map>
 #include <string>
 
 #include "TChain.h"
 #include "TFile.h"
 #include "TTree.h"
 #include "TString.h"
 #include "TObjString.h"
 #include "TSystem.h"
 #include "TROOT.h"
 
 #include "TMVA/Factory.h"
 #include "TMVA/DataLoader.h"
 #include "TMVA/Tools.h"
 #include "TMVA/TMVAGui.h"
 
 int TMVATrain()
 {
    // This loads the library
    TMVA::Tools::Instance();
 
    // Here the preparation phase begins
 
    // Read training and test data
    // (it is also possible to use ASCII format as input -> see TMVA Users Guide)
    TFile *input(0);
    TString fname = "/uboone/data/users/asmith/ubcc1pi/21022020/eventSelectionOutput.root";
    if (!gSystem->AccessPathName( fname )) {
       input = TFile::Open( fname ); // check if file in local directory exists
    }
    if (!input) {
       std::cout << "ERROR: could not open data file" << std::endl;
       exit(1);
    }

    std::cout << "--- TMVAClassification       : Using input file: " << input->GetName() << std::endl;
 
    // Register the training and test trees
    TTree *signalTree     = (TTree*)input->Get("particles");
    TTree *background     = (TTree*)input->Get("particles");
 
    // Apply additional cuts on the signal and background samples (can be different)
    //TCut isSignal = "t_pdgCode == 211 && t_isGolden";
    TCut isSignal = "t_pdgCode == 13";
    TCut quality = "isSignalEvent && r_areFeaturesAvailable && r_isContained && t_truthMatchCompleteness > 0.5";

    TCut mycuts = quality && isSignal;
    TCut mycutb = quality && !isSignal;

    // Create a ROOT output file where TMVA will store ntuples, histograms, etc.
    TString outfileName( "TMVA.root" );
    TFile* outputFile = TFile::Open( outfileName, "RECREATE" );
 
    // Create the factory object. Later you can choose the methods
    // whose performance you'd like to investigate. The factory is
    // the only TMVA object you have to interact with
    //
    // The first argument is the base of the name of all the
    // weightfiles in the directory weight/
    //
    // The second argument is the output file for the training results
    // All TMVA output can be suppressed by removing the "!" (not) in
    // front of the "Silent" argument in the option string
    TMVA::Factory *factory = new TMVA::Factory( "TMVAClassification", outputFile,
                                                "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification" );
 
    TMVA::DataLoader *dataloader=new TMVA::DataLoader("dataset");
 
    // Define the input variables that shall be used for the MVA training
    // note that you may also use variable expressions, such as: "3*var1/var2*abs(var3)"
    // [all types of expressions that can also be parsed by TTree::Draw( "expression" )]
    dataloader->AddVariable("r_logBragg_pToMIP");
    dataloader->AddVariable("r_logBragg_piToMIP");
    dataloader->AddVariable("r_nDescendents");
    dataloader->AddVariable("r_nDownstreamHits");
    dataloader->AddVariable("r_nSpacePointsInSphere5");
    dataloader->AddVariable("r_rmsTrackDeviation");
    dataloader->AddVariable("r_trackShower");
 
    // You can add so-called "Spectator variables", which are not used in the MVA training,
    // but will appear in the final "TestTree" produced by TMVA. This TestTree will contain the
    // input variables, the response values of all trained MVAs, and the spectator variables
    dataloader->AddSpectator("t_hasMatchedMCParticle");
    dataloader->AddSpectator("t_pdgCode");
    dataloader->AddSpectator("t_isGolden");
    dataloader->AddSpectator("t_momentum");
 
    // global event weights per tree (see below for setting event-wise weights)
    Double_t signalWeight     = 1.0;
    Double_t backgroundWeight = 1.0;
 
    // You can add an arbitrary number of signal or background trees
    dataloader->AddSignalTree    ( signalTree,     signalWeight );
    dataloader->AddBackgroundTree( background, backgroundWeight );
 
    // Tell the dataloader how to use the training and testing events
    dataloader->PrepareTrainingAndTestTree( mycuts, mycutb, "nTrain_Signal=2000:nTrain_Background=2000:SplitMode=Random:NormMode=NumEvents:!V");
 
    // Book the BDT
    //factory->BookMethod( dataloader, TMVA::Types::kBDT, "BDT", "!H:!V:NTrees=850:MinNodeSize=2.5%:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=20");
    
    factory->BookMethod( dataloader, TMVA::Types::kBDT, "BDT",
        "!H:!V:NTrees=505:MinNodeSize=1.26436%:MaxDepth=4:BoostType=AdaBoost:AdaBoostBeta=0.2:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=20");

    // Optimize for classification
    //factory->OptimizeAllMethodsForClassification("ROCIntegral","Scan");
 
    // Now you can tell the factory to train, test, and evaluate the MVAs
    //
    // Train MVAs using the set of training events
    factory->TrainAllMethods();
 
    // Evaluate all MVAs using the set of test events
    factory->TestAllMethods();
 
    // Evaluate and compare performance of all configured MVAs
    factory->EvaluateAllMethods();
 
    // --------------------------------------------------------------
 
    // Save the output
    outputFile->Close();
 
    std::cout << "==> Wrote root file: " << outputFile->GetName() << std::endl;
    std::cout << "==> TMVAClassification is done!" << std::endl;
 
    delete factory;
    delete dataloader;

    // Launch the GUI for the root macros
    if (!gROOT->IsBatch()) TMVA::TMVAGui( outfileName );
    
    return 0;
 }
