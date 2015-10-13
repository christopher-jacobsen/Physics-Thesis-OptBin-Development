//
//  OptBin.h
//  OptBin
//
//  Created by Christopher Jacobsen on 10/8/2015.
//  Copyright (c) 2015 Christopher Jacobsen. All rights reserved.
//

#ifndef OPT_BIN_H
#define OPT_BIN_H

#include "common.h"
#include "ModelCompare.h"

////////////////////////////////////////////////////////////////////////////////

namespace OptBin
{

////////////////////////////////////////////////////////////////////////////////

void OptBin( const ModelCompare::ObservableVector & observables,
             const ModelCompare::ModelFileVector &  models,
             size_t startBins,
             const char * cacheFileName );

void OptBinUnbinned( const ModelCompare::ObservableVector & observables,
                     const ModelCompare::ModelFileVector &  models,
                     const char * cacheFileName );

////////////////////////////////////////////////////////////////////////////////

} // namespace OptBin

#endif // OPT_BIN_H
