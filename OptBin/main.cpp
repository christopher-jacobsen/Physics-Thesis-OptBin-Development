//
//  main.cpp
//  OptBin
//
//  Created by Christopher Jacobsen on 10/8/2015.
//  Copyright (c) 2015 Christopher Jacobsen. All rights reserved.
//

#include "OptBin.h"

#include "common.h"
#include "RootUtil.h"
#include "ModelCompare.h"
#include "ReweightEFT.h"

using namespace RootUtil;
using namespace ModelCompare;
using namespace ReweightEFT;
using namespace OptBin;

////////////////////////////////////////////////////////////////////////////////

static const ModelCompare::ObservableVector Observables1 =
{
    // phase-space observables

    { "PTZ", "P_{T}(Z)",  150,     0,  750, "P_{T}(Z) [GeV/c]", "Events per 5 GeV/c",       GETOBS{ GetObs(s,v,c, GetObsPT,   24);     } },
    { "MWZ", "M(WZ)",     150,     0, 3000, "M(WZ) [GeV/c^2]",  "Events per 20 GeV/c^2",    GETOBS{ GetObs(s,v,c, GetObsMass, 24, 23); } },
    { "RAZ", "Y(Z)",      100,    -5,    5, "Y(Z)",             "Events per bin",           GETOBS{ GetObs(s,v,c, GetObsRap,  24);     } },

    // optimal observables

    { "cWWW_O1",    "O_{1}(cWWW)",  1000, -6E4,  7E3,   "O_{1}(cWWW) [GeV^{2}]", "Events per bin",  GETOBS{ GetObs(s,v,c, GetObsOpt, "F_0_1_ocWWW"); } },
    { "cWWW_O2",    "O_{2}(cWWW)",  1000,  0,   2E12,   "O_{2}(cWWW) [GeV^{4}]", "Events per bin",  GETOBS{ GetObs(s,v,c, GetObsOpt, "F_1_1_ocWWW"); } },

    { "cW_O1",      "O_{1}(cW)",    1000, -2E6,  2E4,   "O_{1}(cW) [GeV^{2}]",   "Events per bin",  GETOBS{ GetObs(s,v,c, GetObsOpt, "F_0_2_ocW");   } },
    { "cW_O2",      "O_{2}(cW)",    1000,  0,   6E11,   "O_{2}(cW) [GeV^{4}]",   "Events per bin",  GETOBS{ GetObs(s,v,c, GetObsOpt, "F_2_2_ocW");   } },

    { "cB_O1",      "O_{1}(cB)",    1000, -8E2, 5E3,    "O_{1}(cB) [GeV^{2}]",   "Events per bin",  GETOBS{ GetObs(s,v,c, GetObsOpt, "F_0_3_ocB");   } },
    { "cB_O2",      "O_{2}(cB)",    1000,  0,   2E8,    "O_{2}(cB) [GeV^{4}]",   "Events per bin",  GETOBS{ GetObs(s,v,c, GetObsOpt, "F_3_3_ocB");   } },
};

////////////////////////////////////////////////////////////////////////////////

static const ModelCompare::ModelFileVector Models_1E4 =
{
    { "weight/SM_220_weight_1E4.hepmc2g.gz",    "SM_220",   "SM",   18.2613, 0.154079, 10000 },
    { "weight/EFT_all_weight_1E4.hepmc2g.gz",   "EFT_all",  "EFT",  42.9091, 0.379962, 10000 },
};

static const ModelCompare::ModelFileVector Models_1E6 =
{
    { "weight/SM_220_weight_1E6.hepmc2g.gz",    "SM_220",   "SM",   18.5537, 0.0156025, 1000000 },
    { "weight/EFT_all_weight_1E6.hepmc2g.gz",   "EFT_all",  "EFT",  42.8051, 0.0379205, 1000000 },
};

////////////////////////////////////////////////////////////////////////////////
int main(int argc, const char * argv[])
{
    //OptBin( Observables1, Models_1E6, 16384, "optbin/Cache_16384_1E6.root" );

    //OptBin( Observables1, Models_1E6, 10000, "optbin/Cache_10000_1E6.root" );

    //OptBin( Observables1, Models_1E6, 1000000, "optbin/Cache_1000000_1E6.root" );

  //OptBinUnbinned( Observables1, Models_1E4, "optbin/Cache_Unbinned_1E4.root" );
    OptBinUnbinned( Observables1, Models_1E6, "optbin/Cache_Unbinned_1E6.root" );

    LogMsgInfo( "\nDone." );
    return 0;
}
