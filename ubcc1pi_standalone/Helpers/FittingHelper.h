/**
 *  @file  ubcc1pi_standalone/Helpers/FittingHelper.h
 *
 *  @brief The header file for the fitting helper class
 */

#ifndef UBCC1PI_STANDALONE_HELPERS_FITTING_HELPER
#define UBCC1PI_STANDALONE_HELPERS_FITTING_HELPER

#include "ubcc1pi_standalone/Objects/Config.h"
#include "ubsmear.h"
#include <TVector.h>

namespace ubcc1pi
{
    /**
     *  @brief  The fitting helper class
     */
    class FittingHelper
    {
        public:

            /**
             *  @brief  Constructor
             *
             *  @param  xxxxxxxxxxxxx
             *
             *  @return xxxxxxxxxxxxx
             */
            FittingHelper(const Int_t binNumber);

            /**
             *  @brief  xxxxxxxxxxxxx
             *
             *  @param  xxxxxxxxxxxxx
             *
             *  @return xxxxxxxxxxxxx
             */
            void Fit(void(*fcn)(Int_t &, Double_t *, Double_t &f, Double_t *, Int_t), std::pair<std::vector<Double_t>, std::vector<Double_t>> &results, bool &successful, std::vector<float> &covMatrix, const int printlevel);
        
        private:

        //     /**
        //      *  @brief  xxxxxxxxxxxxx
        //      *
        //      *  @param  xxxxxxxxxxxxx
        //      *
        //      *  @return xxxxxxxxxxxxx
        //      */
        //     static void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag) const;
            

            // TVector<ubsmear::UBMatrix>;
            // std::vector<ubsmear::UBMatrix> fitData;
            // auto pFitData = ;
            // ubsmear::UBMatrix x, y, errorY, S;
            Int_t nBins;
    };

} // namespace ubcc1pi

#endif
