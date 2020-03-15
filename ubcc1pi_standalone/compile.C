{
    gSystem->AddIncludePath("  -I../");

    gROOT->ProcessLine("#include <vector>");

    gROOT->ProcessLine(".L Interface/Event.cxx+");
    gROOT->ProcessLine(".L Objects/FileReader.cxx+");

    // Compile the macros
    gROOT->ProcessLine(".L Macros/Test.cxx+");
}
