/**
 *  @file  ubcc1pi_standalone/Macros/ExtractXSecs.cxx
 *
 *  @brief The implementation file of the ExtractXSecs macro
 */

#include "ubcc1pi_standalone/Macros/Macros.h"
#include "ubcc1pi_standalone/Objects/FileReader.h"

#include "ubcc1pi_standalone/Helpers/AnalysisHelper.h"
#include "ubcc1pi_standalone/Helpers/NormalisationHelper.h"
#include "ubcc1pi_standalone/Helpers/SelectionHelper.h"
#include "ubcc1pi_standalone/Helpers/CrossSectionHelper.h"
#include "ubcc1pi_standalone/Helpers/FormattingHelper.h"
// #include "ubcc1pi_standalone/Helpers/FittingHelper.h"
#include "ubsmear.h"

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

void ParamsToTxt(const Config &config)
{

    // // -------------------------------------------------------------------------------------------------------------------------------------
    // // Get the sideband weights
    // // -------------------------------------------------------------------------------------------------------------------------------------
    // // Loop over all cross-section objects
    typedef std::pair<std::vector<Double_t>,std::vector<Double_t>> paramAndErrorPair; // Todo: Improve code!
    //Parameters: selectionName name
    // std::map<std::string,std::map<std::string, std::vector<float>>> cc0piCovarianceMap; 
    std::map<std::string, std::map<std::string, paramAndErrorPair>> cc0piNominalConstraintMap;
    //Parameters: selectionName name paramName (i.e. golden muonMomentum hadronProduction)
    // std::map<std::string, std::map<std::string, std::map<std::string, std::vector<paramAndErrorPair>>>> cc0piUniverseConstraintMap;

    // std::ifstream ifs1("cc0piCovarianceMap.bin", std::ios::binary);
    std::ifstream ifs2("cc0piNominalConstraintMap.bin", std::ios::binary);
    // std::ifstream ifs3("cc0piUniverseConstraintMap.bin", std::ios::binary);
    
    // boost::archive::binary_iarchive iarch1(ifs1);
    boost::archive::binary_iarchive iarch2(ifs2);
    // boost::archive::binary_iarchive iarch3(ifs3);

    // iarch1 >> cc0piCovarianceMap;
    iarch2 >> cc0piNominalConstraintMap;
    // iarch3 >> cc0piUniverseConstraintMap;

    // ifs1.close();
    ifs2.close();
    // ifs3.close();

    for (const std::string name : {"muonCosTheta", "muonPhi", "muonMomentum", "pionCosTheta", "pionPhi", "pionMomentum", "muonPionAngle", "nProtons"})
    {
        if(!config.extractXSecs.crossSectionIsEnabled.at("generic").at(name)) continue;
        std::cout << "Processing cross-section: "<< name << std::endl;

        const auto cc0piNominalConstraintParam = cc0piNominalConstraintMap.at("generic").at(name).first;
        const vector<float> cc0piNominalConstraintParamFloat(cc0piNominalConstraintParam.begin(), cc0piNominalConstraintParam.end());
        const ubsmear::UBMatrix cc0piNominalConstraintParamMatrix(cc0piNominalConstraintParamFloat, cc0piNominalConstraintParamFloat.size(), 1);
        const auto cc0piNominalConstraintParamError = cc0piNominalConstraintMap.at("generic").at(name).second;
        const vector<float> cc0piNominalConstraintParamErrorFloat(cc0piNominalConstraintParamError.begin(), cc0piNominalConstraintParamError.end());
        const ubsmear::UBMatrix cc0piNominalConstraintParamErrorMatrix(cc0piNominalConstraintParamErrorFloat, cc0piNominalConstraintParamErrorFloat.size(), 1);

        FormattingHelper::SaveMatrix(cc0piNominalConstraintParamMatrix, "Sideband_" + name + "_nominal_parameters.txt");
        FormattingHelper::SaveMatrix(cc0piNominalConstraintParamErrorMatrix, "Sideband_" + name + "_nominal_parameter_errors.txt");
    }

    std::cout<<"------------- All Done -------------"<<std::endl;
    return;
}

} // namespace ubcc1pi_macros
