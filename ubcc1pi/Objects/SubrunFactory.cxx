/**
 *  @file  ubcc1pi/Objects/SubrunFactory.cxx
 *
 *  @brief The implementation of the subrun factory class
 */

#include "ubcc1pi/Objects/SubrunFactory.h"

#include "ubcc1pi/Helpers/CollectionHelper.h"

namespace ubcc1pi
{

void SubrunFactory::PopulateSubrun(const art::SubRun &subrun, const Config &config, Subrun *pOutputSubrun)
{
    pOutputSubrun->Reset();

    pOutputSubrun->run.Set(subrun.run());
    pOutputSubrun->subRun.Set(subrun.subRun());
    pOutputSubrun->hasTruthInfo.Set(config.HasTruthInfo());

    if (!config.HasTruthInfo())
        return;

    const auto potSummary = CollectionHelper::GetObject<sumdata::POTSummary>(subrun, config.POTSummaryLabel());
    pOutputSubrun->totalPOT.Set(potSummary.totpot);
    pOutputSubrun->totalGoodPOT.Set(potSummary.totgoodpot);
    pOutputSubrun->totalSpills.Set(potSummary.totspills);
    pOutputSubrun->goodSpills.Set(potSummary.goodspills);
}

} // namespace ubcc1pi
