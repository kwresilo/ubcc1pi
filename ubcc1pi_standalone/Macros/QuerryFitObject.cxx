/**
 *  @file  ubcc1pi_standalone/Macros/ExtractNuWroXSecs2.cxx
 *
 *  @brief The implementation file of the ExtractNuWroXSecs2 macro
 */

#include "ubcc1pi_standalone/Macros/Macros.h"
#include "ubcc1pi_standalone/Objects/FileReader.h"

// #include "ubcc1pi_standalone/Helpers/AnalysisHelper.h"
// #include "ubcc1pi_standalone/Helpers/NormalisationHelper.h"
// #include "ubcc1pi_standalone/Helpers/SelectionHelper.h"
// #include "ubcc1pi_standalone/Helpers/CrossSectionHelper.h"
// #include "ubcc1pi_standalone/Helpers/FormattingHelper.h"
// #include "ubcc1pi_standalone/Helpers/FittingHelper.h"
// #include "ubsmear.h"

#include <fstream> // Todo: not use txt files
// Boost libraries
#include "binary_iarchive.hpp"
// #include "binary_oarchive.hpp"
#include "binary_object.hpp"
#include "map.hpp"
#include "vector.hpp"

using namespace ubcc1pi;

namespace ubcc1pi_macros
{

void QuerryFitObject(const Config &config)
{
    // // -------------------------------------------------------------------------------------------------------------------------------------
    // // Load the sideband weights from file
    // // -------------------------------------------------------------------------------------------------------------------------------------
    // // Loop over all cross-section objects

    // Loop over all cross-section objects
    typedef std::pair<std::vector<Double_t>, std::vector<Double_t>> paramAndErrorPair; // Todo: Improve code!
    //Parameters: dataTypeName selectionName name
    typedef std::map<std::string, std::map<std::string, std::map<std::string, paramAndErrorPair>>> nominalFitMap;
    //Parameters: dataTypeName selectionName name paramName (i.e. golden muonMomentum hadronProduction)
    typedef std::map<std::string, std::map<std::string, std::map<std::string, std::map<std::string, std::vector<paramAndErrorPair>>>>> universeFitMap;

    nominalFitMap cc0piNominalConstraintMap;
    universeFitMap cc0piUniverseConstraintMap;
    // std::ifstream ifs1("cc0piCovarianceMap.bin", std::ios::binary);
    std::ifstream ifs2("cc0piNominalConstraintMap.bin", std::ios::binary);
    std::ifstream ifs3("cc0piUniverseConstraintMap.bin", std::ios::binary);
    // std::ifstream ifs4("cc0piCovarianceMapData.bin", std::ios::binary);
    // std::ifstream ifs5("cc0piNominalConstraintMap.bin", std::ios::binary);
    // std::ifstream ifs6("cc0piUniverseConstraintMap.bin", std::ios::binary);


    // boost::archive::binary_iarchive iarch1(ifs1);
    boost::archive::binary_iarchive iarch2(ifs2);
    boost::archive::binary_iarchive iarch3(ifs3);
    // boost::archive::binary_iarchive iarch4(ifs4);
    // boost::archive::binary_iarchive iarch5(ifs5);
    // boost::archive::binary_iarchive iarch6(ifs6);

    // iarch1 >> cc0piCovarianceMap;
    iarch2 >> cc0piNominalConstraintMap;
    iarch3 >> cc0piUniverseConstraintMap;
    // iarch4 >> cc0piCovarianceMapBNB;
    // iarch5 >> cc0piNominalConstraintMapBNB;
    // iarch6 >> cc0piUniverseConstraintMapBNB;

    const auto name = "muonPhi";
    const auto paramName = "DecayAngMEC_UBGenie"; // "VecFFCCQEshape_UBGenie";

    for(const auto &dataTypeName: {"Genie"})//{"BNB", "NuWro"})
    {
        for(const auto &selectionName: {"generic"})//, "golden"})
        {
            const auto cc0piNominalConstraintVector = cc0piNominalConstraintMap.at(dataTypeName).at(selectionName).at(name);

            std::cout<<dataTypeName<<" nominal "<<selectionName<<": "<<std::endl;
            std::cout<<"Parameters:"<<"\n\t";
            for(const auto &v: cc0piNominalConstraintVector.first)
            {
                std::cout << v << " ";
            }
            std::cout<<std::endl;

            std::cout<<"Parameter errors:"<<"\n\t";
            for(const auto &v: cc0piNominalConstraintVector.second)
            {
                std::cout << v << " ";
            }
            std::cout<<std::endl;

            const auto cc0piUniverseConstraintVector = cc0piUniverseConstraintMap.at(dataTypeName).at(selectionName).at(name).at(paramName);
            std::cout<<dataTypeName<<" universes "<<selectionName<<": "<<std::endl;
            for(const auto &u: cc0piUniverseConstraintVector)
            {
                std::cout<<"Universe values:"<<"\t";
                for(unsigned int i=0; i<u.first.size(); i++)
                {
                    std::cout << u.first.at(i) << "(" << u.second.at(i) << ") ";
                }
                std::cout<<std::endl;
            }
        }
    }


    // ifs1.close();
    ifs2.close();
    ifs3.close();
    // ifs4.close();
    // ifs5.close();
    // ifs6.close();
}

} // namespace ubcc1pi_macros
