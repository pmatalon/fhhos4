#pragma once
#include "ConjugateGradient.h"
#include "FlexibleConjugateGradient.h"
#include "MultigridForHHO.h"
#include "CondensedAMG.h"
#include "AggregAMG.h"
#include "BlockJacobi.h"
#include "EigenCG.h"
#ifdef AGMG_ENABLED
#include "AGMG.h"
#endif