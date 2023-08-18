/**
 *  @file  ubcc1pi_standalone/Helpers/ExtractionHelper.cxx
 *
 *  @brief The implementation file of the Extraction helper class
 */

#include "ubcc1pi_standalone/Helpers/ExtractionHelper.h"
#include "ubcc1pi_standalone/Helpers/NormalisationHelper.h"
#include "ubcc1pi_standalone/Helpers/AnalysisHelper.h"
#include "ubcc1pi_standalone/Helpers/CrossSectionHelper.h"
#include "ubcc1pi_standalone/Helpers/FormattingHelper.h"

#include <stdexcept>

namespace ubcc1pi
{

void ExtractionHelper::ExtractionHelper::PopulateInputFileList(const Config &config, InputFileList &inputData, float &totalExposurePOT)
{
    totalExposurePOT = 0.f;
    if(std::binary_search(config.global.runs.begin(), config.global.runs.end(), 1))
    {
        inputData.emplace_back(AnalysisHelper::Overlay, "", config.filesRun1.overlaysFileName, NormalisationHelper::GetOverlaysNormalisation(config, 1));
        if(config.global.useBNBAsData)
        {
            inputData.emplace_back(AnalysisHelper::DataBNB, "", config.filesRun1.dataBNBFileName, 1.f);
            inputData.emplace_back(AnalysisHelper::Dirt,    "", config.filesRun1.dirtFileName, NormalisationHelper::GetDirtNormalisation(config, 1));
            inputData.emplace_back(AnalysisHelper::DataEXT, "", config.filesRun1.dataEXTFileName, NormalisationHelper::GetDataEXTNormalisation(config, 1));

            // Add the detector variation files
            if(config.global.useDetVar) // useBNBAsData: Detvar not needed when doing only NuWro/Genie fake-data studies
            {
                for (const auto &[name, fileName] : config.filesRun1.detVarFiles)
                {
                    inputData.emplace_back(AnalysisHelper::DetectorVariation, name, fileName, NormalisationHelper::GetDetectorVariationNormalisation(config, name, 1));
                }
            }
        }
        if(config.global.useNuWroAsData)
        {
            inputData.emplace_back(AnalysisHelper::NuWro,   "", config.filesRun1.nuWroFileName, NormalisationHelper::GetNuWroNormalisation(config, 1));
        }

        totalExposurePOT += config.normsRun1.dataBNBTor875WCut / (1e20);
    }

    if(std::binary_search(config.global.runs.begin(), config.global.runs.end(), 2))
    {
        inputData.emplace_back(AnalysisHelper::Overlay, "", config.filesRun2.overlaysFileName, NormalisationHelper::GetOverlaysNormalisation(config, 2));
        if(config.global.useBNBAsData)
        {
            inputData.emplace_back(AnalysisHelper::Dirt,    "", config.filesRun2.dirtFileName, NormalisationHelper::GetDirtNormalisation(config, 2));
            inputData.emplace_back(AnalysisHelper::DataEXT, "", config.filesRun2.dataEXTFileName, NormalisationHelper::GetDataEXTNormalisation(config, 2));
            inputData.emplace_back(AnalysisHelper::DataBNB, "", config.filesRun2.dataBNBFileName, 1.f);

            // Add the detector variation files
            if(config.global.useDetVar) // useBNBAsData: Detvar not needed when doing only NuWro/Genie fake-data studies
            {
                for (const auto &[name, fileName] : config.filesRun2.detVarFiles)
                {
                    inputData.emplace_back(AnalysisHelper::DetectorVariation, name, fileName, NormalisationHelper::GetDetectorVariationNormalisation(config, name, 2));
                }
            }
        }
        if(config.global.useNuWroAsData)
        {
            inputData.emplace_back(AnalysisHelper::NuWro,   "", config.filesRun2.nuWroFileName, NormalisationHelper::GetNuWroNormalisation(config, 2));
        }
        totalExposurePOT += config.normsRun2.dataBNBTor875WCut / (1e20);
    }

    if(std::binary_search(config.global.runs.begin(), config.global.runs.end(), 3))
    {
        inputData.emplace_back(AnalysisHelper::Overlay, "", config.filesRun3.overlaysFileName, NormalisationHelper::GetOverlaysNormalisation(config, 3));
        if(config.global.useBNBAsData)
        {
            inputData.emplace_back(AnalysisHelper::Dirt,    "", config.filesRun3.dirtFileName, NormalisationHelper::GetDirtNormalisation(config, 3));
            inputData.emplace_back(AnalysisHelper::DataEXT, "", config.filesRun3.dataEXTFileName, NormalisationHelper::GetDataEXTNormalisation(config, 3));
            inputData.emplace_back(AnalysisHelper::DataBNB, "", config.filesRun3.dataBNBFileName, 1.f);

            // Add the detector variation files
            if(config.global.useDetVar) // useBNBAsData: Detvar not needed when doing only NuWro/Genie fake-data studies
            {
                for (const auto &[name, fileName] : config.filesRun3.detVarFiles)
                {
                    inputData.emplace_back(AnalysisHelper::DetectorVariation, name, fileName, NormalisationHelper::GetDetectorVariationNormalisation(config, name, 3));
                }
            }
        }
        if(config.global.useNuWroAsData)
        {
            inputData.emplace_back(AnalysisHelper::NuWro,   "", config.filesRun3.nuWroFileName, NormalisationHelper::GetNuWroNormalisation(config, 3));
        }
        totalExposurePOT += config.normsRun3.dataBNBTor875WCut / (1e20);
    }
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void ExtractionHelper::ExtractionHelper::PopulateAnalysisValueMap(AnalysisValueMap &analysisValueMap, const bool createSidebandVersion)
{
    // Differential cross-section kinematic parameters
    analysisValueMap.emplace("muonCosTheta",  [](const auto &data) { return data.muonCosTheta;  });
    analysisValueMap.emplace("muonPhi",       [](const auto &data) { return data.muonPhi;       });
    analysisValueMap.emplace("muonMomentum",  [](const auto &data) { return data.muonMomentum;  });

    if(createSidebandVersion)
    {
        analysisValueMap.emplace("pionCosTheta",  [](const auto &data) { return data.protonCosTheta;  }); // Getting proton instead of pion values
        analysisValueMap.emplace("pionPhi",       [](const auto &data) { return data.protonPhi;       }); // Getting proton instead of pion values
        analysisValueMap.emplace("pionMomentum",  [](const auto &data) { return data.protonMomentum;  }); // Getting proton instead of pion values
        analysisValueMap.emplace("muonPionAngle", [](const auto &data) { return data.muonProtonAngle; }); // Getting proton instead of pion values
        analysisValueMap.emplace("nProtons",      [](const auto &data) { return data.nProtons-1;      }); // Leading proton treated as pion in CC0pi analysis
    }
    else
    {
        analysisValueMap.emplace("pionCosTheta",  [](const auto &data) { return data.pionCosTheta;  });
        analysisValueMap.emplace("pionPhi",       [](const auto &data) { return data.pionPhi;       });
        analysisValueMap.emplace("pionMomentum",  [](const auto &data) { return data.pionMomentum;  });
        analysisValueMap.emplace("muonPionAngle", [](const auto &data) { return data.muonPionAngle; });
        analysisValueMap.emplace("nProtons",      [](const auto &data) { return data.nProtons;      });
    }

    // ATTN for the total cross-section we don't have an associated kinematic quantity so we just return a dummy value
    const auto dummyValue = 0.f;
    analysisValueMap.emplace("total", [=](const auto &) { return dummyValue; });
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void ExtractionHelper::ExtractionHelper::SaveCrossSectionMatricies(const CrossSectionHelper::CrossSection &xsec, const CrossSectionHelper::CrossSection::ScalingData &scalingData,
    const std::string &selectionName, const std::string &name, const std::string &postfix, const bool disableUncertainties)
{
    // -----------------------------------------------------------------------------------------------------------------------------
    // Get the event rates for data, backgrounds, and signal
    // -----------------------------------------------------------------------------------------------------------------------------
    const auto selectedEventsData = xsec.GetSelectedBNBDataEvents();
    std::cout << "Selected data events" << std::endl;
    FormattingHelper::SaveMatrix(selectedEventsData, "xsec_" + selectionName + "_" + name + "_data_selected_eventRate_" + postfix + ".txt");

    const auto selectedEventsBackground = xsec.GetSelectedBackgroundEvents();
    std::cout << "Selected background events" << std::endl;
    FormattingHelper::SaveMatrix(selectedEventsBackground, "xsec_" + selectionName + "_" + name + "_background_selected_eventRate_" + postfix + ".txt");

    const auto selectedEventsSignal = xsec.GetSelectedSignalEvents();
    std::cout << "Selected signal events" << std::endl;
    FormattingHelper::SaveMatrix(selectedEventsSignal, "xsec_" + selectionName + "_" + name + "_signal_selected_eventRate_" + postfix + ".txt");

    const auto allEventsSignal = xsec.GetSignalEvents();
    std::cout << "All signal events" << std::endl;
    FormattingHelper::SaveMatrix(allEventsSignal, "xsec_" + selectionName + "_" + name + "_signal_all_eventRate_" + postfix + ".txt");
    std::cout << "After all signal events" << std::endl;
    // -----------------------------------------------------------------------------------------------------------------------------
    // Get the cross-section as measured with data along with it's uncertainties
    // -----------------------------------------------------------------------------------------------------------------------------
    const auto data = xsec.GetBNBDataCrossSection(scalingData);
    std::cout << "Data cross-section (reco-space)" << std::endl;
    FormattingHelper::SaveMatrix(data, "xsec_" + selectionName + "_" + name + "_data_" + postfix + ".txt");

    if(!disableUncertainties)
    {
        const auto dataStatUncertainties = xsec.GetBNBDataCrossSectionStatUncertainty(scalingData);
        std::cout << "Data stat uncertainty" << std::endl;
        FormattingHelper::SaveMatrix(dataStatUncertainties, "xsec_" + selectionName + "_" + name + "_data_stat_" + postfix + ".txt");

        const auto dataSystBiasCovariances = xsec.GetBNBDataCrossSectionSystUncertainties(scalingData);
        for (const auto &[group, map] : dataSystBiasCovariances)
        {
            for (const auto &[paramName, biasCovariance] : map)
            {
                const auto &[pBias, pCovariance] = biasCovariance;

                std::cout << "Data syst uncertainty: " << group << " " << paramName << std::endl;
                std::cout << "Bias vector" << std::endl;
                FormattingHelper::SaveMatrix(*pBias, "xsec_" + selectionName + "_" + name + "_data_" + group + "_" + paramName + "_bias_" + postfix + ".txt");
                std::cout << "Covariance matrix" << std::endl;
                FormattingHelper::SaveMatrix(*pCovariance, "xsec_" + selectionName + "_" + name + "_data_" + group + "_" + paramName + "_covariance_" + postfix + ".txt");
            }
        }
    }

    // -----------------------------------------------------------------------------------------------------------------------------
    // Get the predicted cross-section along with its MC stat uncertainty
    // -----------------------------------------------------------------------------------------------------------------------------
    const auto prediction = xsec.GetPredictedCrossSection(scalingData);
    std::cout << "Predicted cross-section (truth-space)" << std::endl;
    FormattingHelper::SaveMatrix(prediction, "xsec_" + selectionName + "_" + name + "_prediction_" + postfix + ".txt");

    if(!disableUncertainties)
    {
        const auto &[pPredictionStatBias, pPredictionStatCovariance] = xsec.GetPredictedCrossSectionStatUncertainty(scalingData);
        std::cout << "Predicted cross-section stat uncertainty" << std::endl;
        std::cout << "Bias vector" << std::endl;
        FormattingHelper::SaveMatrix(*pPredictionStatBias, "xsec_" + selectionName + "_" + name + "_prediction_stat_bias_" + postfix + ".txt");
        std::cout << "Covariance matrix" << std::endl;
        FormattingHelper::SaveMatrix(*pPredictionStatCovariance, "xsec_" + selectionName + "_" + name + "_prediction_stat_covariance_" + postfix + ".txt");


        const auto &[pPredictionSidebandStatBias, pPredictionSidebandStatCovariance] = xsec.GetPredictedSidebandCrossSectionStatUncertainty(scalingData);
        std::cout << "Predicted cross-section stat uncertainty" << std::endl;
        std::cout << "Bias vector" << std::endl;
        FormattingHelper::SaveMatrix(*pPredictionSidebandStatBias, "xsec_" + selectionName + "_" + name + "_prediction_sideband_stat_bias_" + postfix + ".txt");
        std::cout << "Covariance matrix" << std::endl;
        FormattingHelper::SaveMatrix(*pPredictionSidebandStatCovariance, "xsec_" + selectionName + "_" + name + "_prediction_sideband_stat_covariance_" + postfix + ".txt");
    }

    // -----------------------------------------------------------------------------------------------------------------------------
    // Get the smearing matrix along with its uncertainties
    // -----------------------------------------------------------------------------------------------------------------------------
    std::cout << "Smearing Matrix (reco-space rows, truth-space columns)" << std::endl;
    const auto smearingMatrix = xsec.GetSmearingMatrix();
    FormattingHelper::SaveMatrix(smearingMatrix, "xsec_" + selectionName + "_" + name + "_smearingMatrix_" + postfix + ".txt");

    std::cout << "Smearing Matrix AllSelected" << std::endl;
    const auto smearingMatrixAllSelected = xsec.GetSmearingMatrixAllSelected();
    FormattingHelper::SaveMatrix(smearingMatrixAllSelected, "xsec_" + selectionName + "_" + name + "_smearingMatrixAllSelected_" + postfix + ".txt");

    if(!disableUncertainties)
    {
        std::cout << "Smearing Matrix SystBiasCovariances" << std::endl;
        const auto smearingMatrixSystBiasCovariances = xsec.GetSmearingMatrixSystUncertainties();
        std::cout << "Smearing Matrix SystBiasCovariances - After" << std::endl;
        for (const auto &[group, map] : smearingMatrixSystBiasCovariances)
        {
            for (const auto &[paramName, biasCovariance] : map)
            {
                const auto &[pBias, pCovariance] = biasCovariance;

                std::cout << "Smearing matrix syst uncertainty: " << group << " " << paramName << std::endl;
                std::cout << "Bias vector" << std::endl;
                FormattingHelper::SaveMatrix(*pBias, "xsec_" + selectionName + "_" + name + "_smearingMatrix_" + group + "_" + paramName + "_bias_" + postfix + ".txt");
                std::cout << "Covariance matrix" << std::endl;
                FormattingHelper::SaveMatrix(*pCovariance, "xsec_" + selectionName + "_" + name + "_smearingMatrix_" + group + "_" + paramName + "_covariance_" + postfix + ".txt");
            }
        }
    }
}


}
