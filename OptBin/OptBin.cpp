//
//  OptBin.cpp
//  OptBin
//
//  Created by Christopher Jacobsen on 10/8/2015.
//  Copyright (c) 2015 Christopher Jacobsen. All rights reserved.
//

#include "OptBin.h"

#include <TSystem.h>
#include <TDirectory.h>
#include <TFile.h>

#include <TF1.h>
#include <TH1.h>

#include <TNtupleD.h>
#include <TObjArray.h>
#include <TBranch.h>
#include <TLeaf.h>

#include <HFitInterface.h>
#include <TFitResult.h>
#include <Fit/UnBinData.h>
#include <Foption.h>
#include <Math/MinimizerOptions.h>

#include <Fit/LogLikelihoodFCN.h>
#include <Math/WrappedMultiTF1.h>

#include <TVirtualFitter.h>


using namespace RootUtil;

////////////////////////////////////////////////////////////////////////////////

namespace OptBin
{

////////////////////////////////////////////////////////////////////////////////
void GetBinStats( const TH1D & hist, double & meanEntries, double & varEntries )
{
    meanEntries = 0;
    varEntries  = 0;

    double sum  = 0;    // sum of entries
    double sum2 = 0;    // sum of entries^2

    Int_t nBins = hist.GetSize() - 2;
    for (Int_t bin = 1; bin <= nBins; ++bin)
    {
        double binEntries = hist.GetBinContent(bin);
        sum  += binEntries;
        sum2 += binEntries * binEntries;
    }

    meanEntries =  sum  / nBins;
    varEntries  = (sum2 / nBins) - (meanEntries * meanEntries);
}

////////////////////////////////////////////////////////////////////////////////
double GetErrorMeasure( double meanEntries, double varEntries, double binWidth )
{
    return (2 * meanEntries - varEntries) / (binWidth * binWidth);
}

////////////////////////////////////////////////////////////////////////////////
Int_t OptBinPass1( const TH1D & hist )
{
    LogMsgInfo( "\nPass 1: Optimizing bin size for %hs", FMT_HS(hist.GetName()) );
    LogMsgInfo( "nBins  range  width        mean     stddev   measure");
    LogMsgInfo( "------------------------------------------------------------");

    std::unique_ptr<TH1D> upRebin( (TH1D *)hist.Clone() );

    Int_t  nOptBins    = 0;
    double optWidth    = 0;
    double bestMeasure = std::numeric_limits<double>::max();

    Int_t nBins = upRebin->GetSize() - 2;
    while (nBins >= 4)
    {
        double range    = upRebin->GetXaxis()->GetXmax() - upRebin->GetXaxis()->GetXmin();
        double binWidth = range / nBins;

        double meanEntries, varEntries;
        GetBinStats( *upRebin, meanEntries, varEntries );

        double measure = GetErrorMeasure( meanEntries, varEntries, binWidth );

        LogMsgInfo( "%-6i %-6g %-12g %-8g %-8g %g",
            FMT_I(nBins), FMT_F(range), FMT_F(binWidth),
            FMT_F(meanEntries), FMT_F(std::sqrt(varEntries)),
            FMT_F(measure) );

        if (measure < bestMeasure)
        {
            bestMeasure = measure;
            nOptBins    = nBins;
            optWidth    = binWidth;
        }

        upRebin->Rebin(2);  // half the number of bins

        nBins = upRebin->GetSize() - 2;
    }

    LogMsgInfo( "Pass1: Optimal bins=%i width=%g\n", FMT_I(nOptBins), FMT_F(optWidth) );

    return nOptBins;
}

////////////////////////////////////////////////////////////////////////////////
Int_t OptBinPass2( const TH1D & hist, Int_t binPass1 )
{
    LogMsgInfo( "\nPass 2: Optimizing bin size for %hs", FMT_HS(hist.GetName()) );
    LogMsgInfo( "nBins  range  width        mean     stddev   measure");
    LogMsgInfo( "------------------------------------------------------------");

    std::unique_ptr<TH1D> upOrig( (TH1D *)hist.Clone() );

    Int_t  nOptBins    = 0;
    double optWidth    = 0;
    double bestMeasure = std::numeric_limits<double>::max();

    Int_t orgBins  = hist.GetSize() - 2;
    Int_t ctrGroup = orgBins  / binPass1;
    Int_t minGroup = ctrGroup / 2;
    Int_t maxGroup = ctrGroup * 2;

    gErrorIgnoreLevel = kError;

    for (Int_t group = minGroup; group <= maxGroup; ++group)
    {
        std::unique_ptr<TH1D> upRebin( (TH1D *)upOrig->Rebin(group, "rebin") );

        Int_t nBins = upRebin->GetSize() - 2;

        double range    = upRebin->GetXaxis()->GetXmax() - upRebin->GetXaxis()->GetXmin();
        double binWidth = range / nBins;

        double meanEntries, varEntries;
        GetBinStats( *upRebin, meanEntries, varEntries );

        double measure = GetErrorMeasure( meanEntries, varEntries, binWidth );

        LogMsgInfo( "%-6i %-6g %-12g %-8g %-8g %g",
            FMT_I(nBins), FMT_F(range), FMT_F(binWidth),
            FMT_F(meanEntries), FMT_F(std::sqrt(varEntries)),
            FMT_F(measure) );

        if (measure < bestMeasure)
        {
            bestMeasure = measure;
            nOptBins    = nBins;
            optWidth    = binWidth;
        }
    }

    LogMsgInfo( "Pass2: Optimal bins=%i width=%g\n", FMT_I(nOptBins), FMT_F(optWidth) );

    return nOptBins;
}


////////////////////////////////////////////////////////////////////////////////
Int_t OptBin3( const TH1D & hist )
{
    LogMsgInfo( "\nOptimizing bin size for %hs", FMT_HS(hist.GetName()) );
    LogMsgInfo( "nBins  range  width        mean     stddev   measure");
    LogMsgInfo( "------------------------------------------------------------");

    std::unique_ptr<TH1D> upOrig( (TH1D *)hist.Clone() );

    Int_t  nOptBins    = 0;
    double optWidth    = 0;
    double bestMeasure = std::numeric_limits<double>::max();

    const double orgMin   = upOrig->GetXaxis()->GetXmin();
    const double orgMax   = upOrig->GetXaxis()->GetXmax();
    const double orgRange = orgMax - orgMin;
    const double minWidth = orgRange / 10000;
    const double maxWidth = orgRange / 10;

    const int steps = 100;
    for (int i = 0; i <= steps; ++i)
    {
        double binWidth = minWidth * pow( maxWidth/minWidth, double(i)/steps );

        double fBins = orgRange / binWidth;
        Int_t  nBins = (Int_t)std::round(fBins);

        binWidth = orgRange / nBins;

        std::vector<double> edges(nBins+1);

        for (Int_t bin = 0; bin <= nBins; ++bin)
            edges[bin] = orgMin + binWidth * bin;
        edges[nBins] = orgMax;

        std::unique_ptr<TH1D> upRebin( (TH1D *)upOrig->Rebin( nBins, "rebin", edges.data() ) );

        Int_t xBins = upRebin->GetSize() - 2;
        if (xBins != nBins)
            ThrowError( "Something is wrong" );

        double range = upRebin->GetXaxis()->GetXmax() - upRebin->GetXaxis()->GetXmin();

        double meanEntries, varEntries;
        GetBinStats( *upRebin, meanEntries, varEntries );

        double measure = GetErrorMeasure( meanEntries, varEntries, binWidth );

        LogMsgInfo( "%-6i %-6g %-12g %-8g %-8g %g",
            FMT_I(nBins), FMT_F(range), FMT_F(binWidth),
            FMT_F(meanEntries), FMT_F(std::sqrt(varEntries)),
            FMT_F(measure) );

        if (measure < bestMeasure)
        {
            bestMeasure = measure;
            nOptBins    = nBins;
            optWidth    = binWidth;
        }
    }


    LogMsgInfo( "Optimal bins=%i width=%g\n", FMT_I(nOptBins), FMT_F(optWidth) );

    return nOptBins;
}

////////////////////////////////////////////////////////////////////////////////
void OptBin( const ModelCompare::ObservableVector & observables,
             const ModelCompare::ModelFileVector &  models,
             size_t startBins,
             const char * cacheFileName )
{
    // replace the nBins in the observables with startBins
    ModelCompare::ObservableVector setObs(observables);
    {
        for (ModelCompare::Observable & obs: setObs)
            obs.nBins = (Int_t)startBins;
    }

    std::vector<TH1DVector> modelData;  // modelData[model][observable]

    // load the model data for each model and observable
    LoadHistData( models, setObs, modelData, cacheFileName );

    for (const TH1DVector & hists : modelData)
    {
        for (const TH1D * pHist : hists)
        {
            Int_t nOpt1 = OptBinPass1( *pHist );
            Int_t nOpt2 = OptBinPass2( *pHist, nOpt1 );
            //OptBin3( *pHist );
        }
    }
}

////////////////////////////////////////////////////////////////////////////////
TNtupleD * LoadTuple( const char * fileName, const char * tupleName )
{
    if (gSystem->AccessPathName( fileName ))
        return nullptr;

    struct Cleanup
    {
        TDirectory * oldDir = gDirectory;

        ~Cleanup()
        {
            if (oldDir)
                oldDir->cd();
        }

    } cleanup;

    std::unique_ptr<TFile> upFile( new TFile(fileName, "READ") );
    if (upFile->IsZombie() || !upFile->IsOpen())    // IsZombie is true if constructor failed
        return nullptr;

    TNtupleD * pTuple = nullptr;
    upFile->GetObject( tupleName, pTuple );
    if (!pTuple)
        return nullptr;

    TNtupleD * pClone = (TNtupleD *)pTuple->Clone();
    pClone->SetDirectory( nullptr );
    pClone->Reset();
    pClone->ResetBranchAddresses();

    Long64_t nEntries = pTuple->GetEntries();
    for (Long64_t entry = 0; entry < nEntries; ++entry)
    {
        pTuple->GetEntry(entry);
        pClone->Fill( pTuple->GetArgs() );
    }

    return pClone;
}

////////////////////////////////////////////////////////////////////////////////
bool LoadCacheTuple( const char * cacheFileName, TNtupleD * & pTuple )
{
    if (!cacheFileName || !cacheFileName[0] || !pTuple)
        return false;

    std::unique_ptr<TNtupleD> pTupleCache( LoadTuple( cacheFileName, pTuple->GetName() ) );
    if (!pTupleCache)
        return false;

    // compare cached object to current object
    {
        if (pTupleCache->Class() != pTuple->Class())
            return false;

        if (pTupleCache->GetNvar() != pTuple->GetNvar())
            return false;

        TObjArray * pBranches      = pTuple     ->GetListOfBranches();
        TObjArray * pBranchesCache = pTupleCache->GetListOfBranches();

        if (pBranchesCache->GetEntries() != pBranches->GetEntries())
            return false;

        Int_t nBranches = pBranches->GetEntries();
        for (Int_t branchIndex = 0; branchIndex < nBranches; ++branchIndex)
        {
            TBranch * pBranch      = (TBranch *)pBranches     ->At(branchIndex);
            TBranch * pBranchCache = (TBranch *)pBranchesCache->At(branchIndex);

            if (strcmp( pBranchCache->GetName(),  pBranch->GetName()  ) != 0)
                return false;
            if (strcmp( pBranchCache->GetTitle(), pBranch->GetTitle() ) != 0)
                return false;

            TObjArray * pLeaves      = pBranch     ->GetListOfLeaves();
            TObjArray * pLeavesCache = pBranchCache->GetListOfLeaves();

            if (pLeavesCache->GetEntries() != pLeaves->GetEntries())
                return false;

            Int_t nLeaves = pLeaves->GetEntries();
            for (Int_t leafIndex = 0; leafIndex < nLeaves; ++leafIndex)
            {
                TLeaf * pLeaf      = (TLeaf *)pLeaves     ->At(leafIndex);
                TLeaf * pLeafCache = (TLeaf *)pLeavesCache->At(leafIndex);

                if (strcmp( pLeafCache->GetName(),  pLeaf->GetName()  ) != 0)
                    return false;
                if (strcmp( pLeafCache->GetTitle(), pLeaf->GetTitle() ) != 0)
                    return false;
            }
        }
    }

    // replace current tuple object with cached one
    delete pTuple;
    pTuple = pTupleCache.release();
    return true;
}

////////////////////////////////////////////////////////////////////////////////
void SaveTuples( const char * fileName, const ConstTupleVector & tuples, const char * option = "UPDATE" )
{
    struct Cleanup
    {
        TDirectory * oldDir = gDirectory;

        ~Cleanup()
        {
            if (oldDir)
                oldDir->cd();
        }

    } cleanup;

    TFile file( fileName, option );
    if (file.IsZombie() || !file.IsOpen())    // IsZombie is true if constructor failed
    {
        LogMsgError( "Failed to create file (%hs).", FMT_HS(fileName) );
        ThrowError( std::invalid_argument( fileName ) );
    }

    for ( const TNtupleD * pTuple : tuples )
    {
        TNtupleD * pClone = (TNtupleD *)pTuple->Clone();
        pClone->SetDirectory( &file );  // owned by output file, which will call delete
        pClone->Write( 0, TObject::kOverwrite );
    }

    file.Close();
}

////////////////////////////////////////////////////////////////////////////////
void LoadTupleData( const ModelCompare::ObservableVector & observables, const ModelCompare::ModelFile & model,
                   TupleVector & modelData,
                   const char * cacheFileName )
{
    modelData.clear();

    bool        bLoadEvents = false;
    TupleVector loadData;

    for (const ModelCompare::Observable & obs : observables)
    {
        // make a n-tuple
        TNtupleD * pTuple = nullptr;
        {
            std::string sName  = obs.BuildHistName(  model.modelName  );
            std::string sTitle = obs.BuildHistTitle( model.modelTitle );

            pTuple = new TNtupleD( sName.c_str(), sTitle.c_str(), obs.name );
            pTuple->SetDirectory( nullptr );  // decouple from any directory
        }

        if (LoadCacheTuple( cacheFileName, pTuple ))
        {
            LogMsgInfo( "Loaded %hs from cache", FMT_HS(pTuple->GetName()) );
            loadData.push_back( nullptr );  // skip this histogram
        }
        else
        {
            loadData.push_back( pTuple );
            bLoadEvents = true;
        }

        modelData.push_back(pTuple);
    }

    // fill the n-tuple

    auto FillFunc = [&](const HepMC::GenVertex & signal)
    {
        size_t obsIndex = 0;
        for (const ModelCompare::Observable & obs : observables)
        {
            TNtupleD * pTuple = loadData[obsIndex++];
            if (pTuple)
            {
                double value(0);
                obs.getFunction( signal, &value, 1 );
                pTuple->Fill(&value);
            }
        }
    };

    if (bLoadEvents)
    {
        LoadEvents( model.fileName, FillFunc );

        SaveTuples( cacheFileName, ToConstTupleVector(modelData) );
    }

    SaveTuples( "optbin/Cross_Unbinned.root", ToConstTupleVector(modelData) );
}

////////////////////////////////////////////////////////////////////////////////
void LoadTupleData( const ModelCompare::ObservableVector & observables, const ModelCompare::ModelFileVector & models,
                   std::vector<TupleVector> & allData,
                   const char * cacheFileName )
{
    allData.clear();

    for (const ModelCompare::ModelFile & model : models)
    {
        TupleVector modelData;
        LoadTupleData( observables, model, modelData, cacheFileName );
        allData.push_back( modelData );
    }
}

////////////////////////////////////////////////////////////////////////////////
TH1D * MakePDFHist( TNtupleD & data, Double_t binWidth, Double_t xMin, Double_t xMax )
{
    // calculate bin edges

    if (binWidth <= 0)
        ThrowError( "Invalid bin width" );

    Double_t fBins = (xMax - xMin) / binWidth;
    Int_t  nBins = (Int_t)std::floor(fBins);
    if (fBins - nBins > 0)  // if remainder
        ++nBins;

    if (nBins <= 0)
        nBins = 1;

    //LogMsgInfo( "%hs: binWidth = %f nBins=%i", FMT_HS(data.GetName()), FMT_F(binWidth), FMT_I(nBins) );

    std::vector<double> edges(nBins+1);

    for (Int_t bin = 0; bin < nBins; ++bin)
        edges[bin] = xMin + binWidth * bin;
    edges[nBins] = xMax;

    // create histogram

    TH1D * pHist = new TH1D( "PDF", "PDF", nBins, edges.data() );
    pHist->SetDirectory( nullptr );
    pHist->Sumw2(kFALSE);  // save performance, we don't need errors

    // fill histogram

    Long64_t nEntries = data.GetEntries();
    for (Long64_t entry = 0; entry < nEntries; ++entry)
    {
        //data.LoadTree(entry);
        data.GetEntry(entry);
        Double_t value = *data.GetArgs();
        pHist->Fill( value );
    }

    if (pHist->GetEntries() != (Double_t)nEntries)
        ThrowError( "Inconsistent number of entries in histogram" );

    // convert to PDF

    //LogMsgHistStats( *pHist );
    //LogMsgHistBinCounts(*pHist);
    //LogMsgHistEffectiveEntries(*pHist);
    //LogMsgHistUnderOverflow(*pHist);

    double integral = pHist->Integral();

    //pHist->Scale( 1.0 / integral );
    //pHist->Scale( 1.0 / integral, "width" );
    //pHist->Scale( 1.0, "width" );
    //pHist->Scale( 1.0, "width" );

    //LogMsgHistDump( *pHist );

    //double integral2 = pHist->Integral("width");

    return pHist;
}

////////////////////////////////////////////////////////////////////////////////
void OptBinUnbinned1( const ModelCompare::Observable & observable, const TNtupleD & data )
{
    LogMsgInfo( "\nOptimizing bin size for %hs", FMT_HS(data.GetName()) );
    LogMsgInfo( "nBins  range  width        mean     stddev   measure");
    LogMsgInfo( "------------------------------------------------------------");

    const Double_t xMin   = observable.xMin;
    const Double_t xMax   = observable.xMax;
    const Double_t xRange = xMax - xMin;

    // clone data so non-const methods can be used
    std::unique_ptr<TNtupleD> upFitData( (TNtupleD *)data.Clone() );

    // define and setup fit lambda function

    size_t              rejectCount(0);
    double_t            lastBinWidth(0);
    ConstTH1DUniquePtr  upLastHist;

    auto PDFFunc = [&](const Double_t * x, const Double_t * par) -> Double_t
    {
        Double_t binWidth = par[0];

        if (!upLastHist || (lastBinWidth != binWidth))
        {
            upLastHist.reset( MakePDFHist( *upFitData, binWidth, xMin, xMax ) );
            lastBinWidth = binWidth;
        }

        const TH1D * pHist = upLastHist.get();

        Double_t xVal = x[0];
        if ((xVal < xMin) || !(xVal < xMax))
        {
            LogMsgInfo( "Rejecting x=%f", FMT_F(xVal) );
            ++rejectCount;
            TF1::RejectPoint();
            return 0;
        }

        Int_t bin = pHist->FindFixBin( xVal );
        if ((bin < 1) || (bin > pHist->GetSize() - 2))
        {
            ++rejectCount;
            TF1::RejectPoint();
            return 0;
        }

        Double_t content = pHist->GetBinContent( bin );

        //LogMsgInfo( "x=%f b=%i n=%f c=%f", FMT_F(xVal), FMT_I(bin), FMT_F(nEff), FMT_F(content) );

        /*
        // reduce NDF
        if (content <= 0)
        {
            // LogMsgInfo( ">>>> Rejecting bin %i <<<<", bin );
            ++rejectCount;
            TF1::RejectPoint();
            return 0;
        }
        */

        return content;
    };

    // define a fit objective
    // definition: TVirtualFitter::FCNFunc_t
    // typedef void (* FCNFunc_t )(Int_t &npar, Double_t *gin, Double_t &f, Double_t *u, Int_t flag);

    //auto FitObjective = [](Int_t & npar, Double_t * gin, Double_t & f, Double_t * u, Int_t flag) -> void
    //{
    //};

    //TVirtualFitter * pFitter = TVirtualFitter::GetFitter();
    //pFitter->SetFCN( FitObjective );

    //////

    // setup fit function

    TF1 fitFunc( "PDF", PDFFunc, xMin, xMax, 1 );
    {
        fitFunc.SetParName(   0, "BinWidth" );
        fitFunc.SetParameter( 0, xRange / 100 );
        //fitFunc.SetParLimits( 0, xRange / 1E4, xRange / 2 );
        //fitFunc.SetParError(  0, xRange / 50 - xRange / 100 );    // set initial step (helps fit converge)
        //fitFunc.SetParError( 0, std::numeric_limits<double>::min() );
    }

    //upFitData->UnbinnedFit( "PDF", observable.name );

    ROOT::Fit::UnBinData unBinData( (unsigned int)upFitData->GetEntries() );
    {
        //ROOT::Fit::DataRange(xMin, xMax),
        //unBinData.Opt().fUseRange = true;

        Long64_t nEntries = upFitData->GetEntries();
        for (Long64_t entry = 0; entry < nEntries; ++entry)
        {
            //upFitData->LoadTree(entry);       // non-const
            upFitData->GetEntry(entry);         // non-const
            Double_t value = *upFitData->GetArgs();
            if ((value < xMin) || !(value < xMax))
                continue;
            unBinData.Add(value);
        }
    }

    /*
    Foption_t fitOption;
    {
        fitOption.StoreResult   = 1;
        fitOption.Nograph       = 1;
        fitOption.Verbose       = 1;
      //fitOption.Like          = 1;  // 1 = extended. Does not function. Calls PDFFunc with x=NaN
    }

    ROOT::Math::MinimizerOptions minOption;
    {
        minOption.SetPrecision( std::numeric_limits<double>::min() );
    }

    TFitResultPtr fitResult = ROOT::Fit::UnBinFit( &unBinData, &fitFunc, fitOption, minOption );

    fitResult->Print();
    */

    ROOT::Math::WrappedMultiTF1 wrapFitFunc(fitFunc);

    ROOT::Fit::LogLikelihoodFCN<ROOT::Math::IMultiGenFunction> logl( unBinData, wrapFitFunc );

    {
        std::map<double,size_t> myMap;

        double binWidth = xRange;
        for (unsigned int index = 0; index < unBinData.NPoints(); ++index)
        {
            double v = logl.DataElement( &binWidth, index, nullptr );
            myMap[v]++;
        }
        for (const auto & entry : myMap)
        {
            LogMsgInfo( "%u of %g", FMT_U(entry.second), FMT_F(entry.first) );
        }
    }

    Int_t nBinInc = 1;
    for (Int_t nBins = 1; nBins <= 10000; nBins += nBinInc)
    {
        if (nBins % (10 * nBinInc) == 0)
            nBinInc *= 10;

        double binWidth = xRange / nBins;

        double llValue = logl( &binWidth );

        LogMsgInfo( "nBins=%-6i  binWidth=%-8g logl=%e", FMT_I(nBins), FMT_F(binWidth), FMT_F(llValue) );
    }

    exit(1);
}

////////////////////////////////////////////////////////////////////////////////
void OptBinUnbinned2( const ModelCompare::Observable & observable, const TNtupleD & data )
{
    LogMsgInfo( "\nOptimizing bin size for %hs", FMT_HS(data.GetName()) );
    LogMsgInfo( "nBins  range  width        mean     stddev   measure");
    LogMsgInfo( "------------------------------------------------------------");

    const Double_t xMin   = observable.xMin;
    const Double_t xMax   = observable.xMax;
    const Double_t xRange = xMax - xMin;

    // clone data so non-const methods can be used
    std::unique_ptr<TNtupleD> upFitData( (TNtupleD *)data.Clone() );

    ConstTH1DUniquePtr upPDF;

    Int_t  nOptBins    = 0;
    double optWidth    = 0;
    double bestMeasure = std::numeric_limits<double>::max();

    Int_t nBinInc = 1;
    Int_t nBinAdj = 10;
    for (Int_t nBins = 1; nBins <= 10000; nBins += nBinInc)
    {
        if (nBins % nBinAdj == 0)
        {
            nBinInc *= 5;
            nBinInc = std::max( nBinInc, nBinAdj / 4 );
            nBinAdj *= 10;
        }

        double binWidth = xRange / nBins;

        upPDF.reset( MakePDFHist( *upFitData, binWidth, xMin, xMax ) );

        double rangePDF = upPDF->GetXaxis()->GetXmax() - upPDF->GetXaxis()->GetXmin();

        double meanEntries, varEntries;
        GetBinStats( *upPDF, meanEntries, varEntries );

        double measure = GetErrorMeasure( meanEntries, varEntries, binWidth );

        LogMsgInfo( "%-6i %-6g %-12g %-8g %-8g %g",
            FMT_I(nBins), FMT_F(rangePDF), FMT_F(binWidth),
            FMT_F(meanEntries), FMT_F(std::sqrt(varEntries)),
            FMT_F(measure) );

        if (measure < bestMeasure)
        {
            bestMeasure = measure;
            nOptBins    = nBins;
            optWidth    = binWidth;
        }
    }

    LogMsgInfo( "Optimal bins=%i width=%g\n", FMT_I(nOptBins), FMT_F(optWidth) );

    //exit(1);
}

////////////////////////////////////////////////////////////////////////////////
void OptBinUnbinned( const ModelCompare::ObservableVector & observables,
                     const ModelCompare::ModelFileVector &  models,
                     const char * cacheFileName )
{
    std::vector<TupleVector> allData;   // allData[model][observable]

    LoadTupleData( observables, models, allData, cacheFileName );

    for (const TupleVector & tuples : allData)
    {
        auto obsItr = observables.cbegin();
        for (const TNtupleD * pTuple : tuples)
        {
            OptBinUnbinned2( *obsItr++, *pTuple );
        }
    }
}

////////////////////////////////////////////////////////////////////////////////

}  // namespace OptBin
