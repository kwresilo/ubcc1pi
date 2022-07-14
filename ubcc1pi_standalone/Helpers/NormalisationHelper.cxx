/**
 *  @file  ubcc1pi_standalone/Helpers/NormalisationHelper.cxx
 *
 *  @brief The implementation file of the normalisation helper class
 */

#include "ubcc1pi_standalone/Helpers/NormalisationHelper.h"

#include <stdexcept>

namespace ubcc1pi
{

float NormalisationHelper::GetOverlaysNormalisation(const Config &config, const unsigned int run)
{
    float overlaysPOT;
    float dataBNBTor875WCut;
    switch (run)
    {
        case 1:
            dataBNBTor875WCut = config.normsRun1.dataBNBTor875WCut;
            overlaysPOT = config.normsRun1.overlaysPOT;
            break;
        case 2:
            dataBNBTor875WCut = config.normsRun2.dataBNBTor875WCut;
            overlaysPOT = config.normsRun2.overlaysPOT;
            break;
        case 3:
            dataBNBTor875WCut = config.normsRun3.dataBNBTor875WCut;
            overlaysPOT = config.normsRun3.overlaysPOT;
            break;
        default:
            std::cout<<"NormalisationHelper::GetOverlaysNormalisation - Invalid run number"<<std::endl;
            throw std::invalid_argument("NormalisationHelper::GetOverlaysNormalisation - Invalid run number");
    }

    if (overlaysPOT <= std::numeric_limits<float>::epsilon())
        throw std::invalid_argument("NormalisationHelper::GetOverlaysNormalisation - Overlay POT is invalid");

    return dataBNBTor875WCut / overlaysPOT;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

float NormalisationHelper::GetNuWroNormalisation(const Config &config, const unsigned int run)
{
    float nuWroPOT;
    float dataBNBTor875WCut;
    switch (run)
    {
        case 1:
            dataBNBTor875WCut = config.normsRun1.dataBNBTor875WCut;
            nuWroPOT = config.normsRun1.nuWroPOT;
            break;
        case 2:
            dataBNBTor875WCut = config.normsRun2.dataBNBTor875WCut;
            nuWroPOT = config.normsRun2.nuWroPOT;
            break;
        case 3:
            dataBNBTor875WCut = config.normsRun3.dataBNBTor875WCut;
            nuWroPOT = config.normsRun3.nuWroPOT;
            break;
        default:
            std::cout<<"NormalisationHelper::GetNuWroNormalisation - Invalid run number"<<std::endl;
            throw std::invalid_argument("NormalisationHelper::GetNuWroNormalisation - Invalid run number");
    }

    if (nuWroPOT <= std::numeric_limits<float>::epsilon())
        throw std::invalid_argument("NormalisationHelper::GetNuWroNormalisation - Overlay POT is invalid");

    return dataBNBTor875WCut / nuWroPOT;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

float NormalisationHelper::GetOverlaysNormalisationToNuWro(const Config &config, const unsigned int run)
{
    float overlaysPOT;
    float nuWroPOT;
    switch (run)
    {
        case 1:
            nuWroPOT = config.normsRun1.nuWroPOT;
            overlaysPOT = config.normsRun1.overlaysPOT;
            break;
        case 2:
            nuWroPOT = config.normsRun2.nuWroPOT;
            overlaysPOT = config.normsRun2.overlaysPOT;
            break;
        case 3:
            nuWroPOT = config.normsRun3.nuWroPOT;
            overlaysPOT = config.normsRun3.overlaysPOT;
            break;
        default:
            std::cout<<"NormalisationHelper::GetOverlaysNormalisation - Invalid run number"<<std::endl;
            throw std::invalid_argument("NormalisationHelper::GetOverlaysNormalisation - Invalid run number");
    }

    if (overlaysPOT <= std::numeric_limits<float>::epsilon())
        throw std::invalid_argument("NormalisationHelper::GetOverlaysNormalisation - Overlay POT is invalid");

    return nuWroPOT / overlaysPOT;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

float NormalisationHelper::GetDirtNormalisation(const Config &config, const unsigned int run)
{
    float dirtPOT;
    float dataBNBTor875WCut;
    switch (run)
    {
        case 1:
            dataBNBTor875WCut = config.normsRun1.dataBNBTor875WCut;
            dirtPOT = config.normsRun1.dirtPOT;
            break;
        case 2:
            dataBNBTor875WCut = config.normsRun2.dataBNBTor875WCut;
            dirtPOT = config.normsRun2.dirtPOT;
            break;
        case 3:
            dataBNBTor875WCut = config.normsRun3.dataBNBTor875WCut;
            dirtPOT = config.normsRun3.dirtPOT;
            break;
        default:
            std::cout<<"NormalisationHelper::GetOverlaysNormalisation - Invalid run number"<<std::endl;
            throw std::invalid_argument("NormalisationHelper::GetOverlaysNormalisation - Invalid run number");
    }

    if (dirtPOT <= std::numeric_limits<float>::epsilon())
        throw std::invalid_argument("NormalisationHelper::GetDirtNormalisation - Dirt POT is invalid");

    return dataBNBTor875WCut / dirtPOT;
}


// -----------------------------------------------------------------------------------------------------------------------------------------

float NormalisationHelper::GetDataEXTNormalisation(const Config &config, const unsigned int run)
{
    float dataEXTTriggers;
    float dataBNBE1DCNTWCut;
    switch (run)
    {
        case 1:
            dataBNBE1DCNTWCut = config.normsRun1.dataBNBE1DCNTWCut;
            dataEXTTriggers = config.normsRun1.dataEXTTriggers;
            break;
        case 2:
            dataBNBE1DCNTWCut = config.normsRun2.dataBNBE1DCNTWCut;
            dataEXTTriggers = config.normsRun2.dataEXTTriggers;
            break;
        case 3:
            dataBNBE1DCNTWCut = config.normsRun3.dataBNBE1DCNTWCut;
            dataEXTTriggers = config.normsRun3.dataEXTTriggers;
            break;
        default:
            std::cout<<"NormalisationHelper::GetOverlaysNormalisation - Invalid run number"<<std::endl;
            throw std::invalid_argument("NormalisationHelper::GetOverlaysNormalisation - Invalid run number");
    }

    if (dataEXTTriggers <= std::numeric_limits<float>::epsilon())
        throw std::invalid_argument("NormalisationHelper::GetDataEXTNormalisation - Data EXT trigger count is invalid");

    return dataBNBE1DCNTWCut / dataEXTTriggers;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

float NormalisationHelper::GetDetectorVariationNormalisation(const Config &config, const std::string &paramName, const unsigned int run)
{
    float dataBNBTor875WCut;
    // std::unordered_map<std::string, float>::iterator iter;
    auto iter = config.normsRun1.detVarPOTs.find(paramName); //TODO: Fix this; dont use find
    switch (run)
    {
        case 1:
            iter = config.normsRun1.detVarPOTs.find(paramName);
            if (config.normsRun1.detVarPOTs.find(paramName) == config.normsRun1.detVarPOTs.end())
                throw std::invalid_argument("NormalisationHelper::GetDetectorVariationNormalisation - Unknown paramter: " + paramName);
            dataBNBTor875WCut = config.normsRun1.dataBNBTor875WCut;
            break;
        case 2:
            iter = config.normsRun2.detVarPOTs.find(paramName);
            if (config.normsRun2.detVarPOTs.find(paramName) == config.normsRun2.detVarPOTs.end())
                throw std::invalid_argument("NormalisationHelper::GetDetectorVariationNormalisation - Unknown paramter: " + paramName);
            dataBNBTor875WCut = config.normsRun2.dataBNBTor875WCut;
            break;
        case 3:
            iter = config.normsRun3.detVarPOTs.find(paramName);
            if (config.normsRun3.detVarPOTs.find(paramName) == config.normsRun3.detVarPOTs.end())
                throw std::invalid_argument("NormalisationHelper::GetDetectorVariationNormalisation - Unknown paramter: " + paramName);
            dataBNBTor875WCut = config.normsRun3.dataBNBTor875WCut;
            break;
        default:
            std::cout<<"NormalisationHelper::GetOverlaysNormalisation - Invalid run number"<<std::endl;
            throw std::invalid_argument("NormalisationHelper::GetOverlaysNormalisation - Invalid run number");
    }

    // // Check we have an entry for this param name
    // if (iter == config.norms.detVarPOTs.end())
    //     throw std::invalid_argument("NormalisationHelper::GetDetectorVariationNormalisation - Unknown paramter: " + paramName);

    // Get the POT from the map
    const auto pot = iter->second;
    if (pot <= std::numeric_limits<float>::epsilon())
        throw std::invalid_argument("NormalisationHelper::GetDetectorVariationNormalisation - POT is invalid for detector parameter: " + paramName);

    return dataBNBTor875WCut / pot;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

float NormalisationHelper::GetDetectorVariationNormalisationToNuWro(const Config &config, const std::string &paramName, const unsigned int run)
{
    float nuWroPOT;
    // std::unordered_map<std::string, float>::iterator iter;
    auto iter = config.normsRun1.detVarPOTs.find(paramName); //TODO: Fix this; dont use find
    switch (run)
    {
        case 1:
            iter = config.normsRun1.detVarPOTs.find(paramName);
            if (config.normsRun1.detVarPOTs.find(paramName) == config.normsRun1.detVarPOTs.end())
                throw std::invalid_argument("NormalisationHelper::GetDetectorVariationNormalisation - Unknown paramter: " + paramName);
            nuWroPOT = config.normsRun1.nuWroPOT;
            break;
        case 2:
            iter = config.normsRun2.detVarPOTs.find(paramName);
            if (config.normsRun2.detVarPOTs.find(paramName) == config.normsRun2.detVarPOTs.end())
                throw std::invalid_argument("NormalisationHelper::GetDetectorVariationNormalisation - Unknown paramter: " + paramName);
            nuWroPOT = config.normsRun2.nuWroPOT;
            break;
        case 3:
            iter = config.normsRun3.detVarPOTs.find(paramName);
            if (config.normsRun3.detVarPOTs.find(paramName) == config.normsRun3.detVarPOTs.end())
                throw std::invalid_argument("NormalisationHelper::GetDetectorVariationNormalisation - Unknown paramter: " + paramName);
            nuWroPOT = config.normsRun3.nuWroPOT;
            break;
        default:
            std::cout<<"NormalisationHelper::GetOverlaysNormalisation - Invalid run number"<<std::endl;
            throw std::invalid_argument("NormalisationHelper::GetOverlaysNormalisation - Invalid run number");
    }

    // // Check we have an entry for this param name
    // if (iter == config.norms.detVarPOTs.end())
    //     throw std::invalid_argument("NormalisationHelper::GetDetectorVariationNormalisation - Unknown paramter: " + paramName);

    // Get the POT from the map
    const auto pot = iter->second;
    if (pot <= std::numeric_limits<float>::epsilon())
        throw std::invalid_argument("NormalisationHelper::GetDetectorVariationNormalisation - POT is invalid for detector parameter: " + paramName);

    return nuWroPOT / pot;
}


}
