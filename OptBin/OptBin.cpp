//
//  OptBin.cpp
//  OptBin
//
//  Created by Christopher Jacobsen on 10/8/2015.
//  Copyright (c) 2015 Christopher Jacobsen. All rights reserved.
//

#include "OptBin.h"
#include "ModelCompare.h"
#include "FitEFT.h"

#include <TF1.h>
#include <TH1.h>
#include <TNtupleD.h>

#include <HFitInterface.h>
#include <TFitResult.h>
#include <Fit/UnBinData.h>
#include <Foption.h>
#include <Math/MinimizerOptions.h>

#include <Fit/LogLikelihoodFCN.h>
#include <Math/WrappedMultiTF1.h>

#include <TMinuitMinimizer.h>
#include <TVirtualFitter.h>

#include <Fit/Fitter.h>
#include <Fit/FitConfig.h>
#include <Fit/FitResult.h>

#include <TGraph.h>
#include <Math/Minimizer.h>
#include <TDirectory.h>
#include <TFile.h>

using namespace RootUtil;

////////////////////////////////////////////////////////////////////////////////

namespace OptBin
{

////////////////////////////////////////////////////////////////////////////////
static void GetBinStats( const TH1D & hist, double & meanEntries, double & varEntries, double binWidth )
{
    meanEntries = 0;
    varEntries  = 0;

    double sum  = 0;    // sum of entries
    double sum2 = 0;    // sum of entries^2

    Int_t  nCount = 0;

    Int_t nBins = hist.GetSize() - 2;
    for (Int_t bin = 1; bin <= nBins; ++bin)
    {
        double binEntries = hist.GetBinContent(bin);

        if (binEntries != 0)
        {
            sum  += binEntries;
            sum2 += binEntries * binEntries;
            ++nCount;
        }
    }

    double fBins = nBins;
    if (binWidth > 0.0)
        fBins = (hist.GetXaxis()->GetXmax() - hist.GetXaxis()->GetXmin()) / binWidth;

    meanEntries =  sum  / fBins;
    varEntries  = (sum2 / fBins) - (meanEntries * meanEntries);
}

////////////////////////////////////////////////////////////////////////////////
static double GetErrorMeasure( double meanEntries, double varEntries, double binWidth )
{
    return (2 * meanEntries - varEntries) / (binWidth * binWidth);
}

////////////////////////////////////////////////////////////////////////////////
static double GetCrossValidation( const TH1D & hist, double binWidth )
{
    double sum  = 0;    // sum of entries
    double sum2 = 0;    // sum of entries^2

    Int_t nBins = hist.GetSize() - 2;
    for (Int_t bin = 1; bin <= nBins; ++bin)
    {
        double binEntries = hist.GetBinContent(bin);
        sum  += binEntries;
        sum2 += binEntries * binEntries;
    }

    double n = sum;
    double A = 2 / (n - 1);
    double B = (n + 1) / (n * n) / (n - 1);

    double result = (A - B * sum2) / binWidth;

    return result;
}

////////////////////////////////////////////////////////////////////////////////
static void SaveGraph( const char * fileName, const TGraph & graph )
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

    TFile file( fileName, "UPDATE" );
    if (file.IsZombie() || !file.IsOpen())    // IsZombie is true if constructor failed
    {
        LogMsgError( "Failed to create file (%hs).", FMT_HS(fileName) );
        ThrowError( std::invalid_argument( fileName ) );
    }

    graph.Write( nullptr, TObject::kOverwrite );

    file.Close();
}

////////////////////////////////////////////////////////////////////////////////
static void MergeZeroBins( TH1D & hist )
{
    // get sorted list of emtpy bins

    std::list<Int_t> emptyBins;

    const Int_t nBins = hist.GetSize() - 2;
    for (Int_t bin = 1; bin <= nBins; ++bin)
    {
        double binEntries = hist.GetBinContent(bin);
        if (binEntries == 0.0)
            emptyBins.push_back(bin);
    }

    if (emptyBins.empty())
        return;

    // construct sorted list of zero regions

    typedef std::pair<Int_t,Int_t> BinPair;

    std::list< BinPair > mergeBins;
    {
        BinPair * pCurrent = nullptr;;

        for (Int_t bin : emptyBins)
        {
            if (pCurrent && (bin == pCurrent->second + 1))
            {
                // extend current region
                pCurrent->second = bin;
                continue;
            }

            mergeBins.push_back( { bin, bin } );
            pCurrent = &mergeBins.back();
        }
    }

    // extend zero regions one bin to each side, respecting [1,nBins] limit
    // now regions will have at least on non-zero bin
    for (BinPair & region : mergeBins)
    {
        region.first  = std::max( region.first  - 1, 1     );
        region.second = std::min( region.second + 1, nBins );
    }

    // merge adjacent regions sharing the same border bin
    if (mergeBins.size() > 1)
    {
        auto itr1 = mergeBins.begin();
        auto itr2 = itr1;
        ++itr2;
        while (itr2 != mergeBins.end())
        {
            if (itr1->second == itr2->first)
            {
                itr1->second = itr2->second;
                itr2 = mergeBins.erase( itr2 );
                continue;
            }

            ++itr1;
            ++itr2;
        }
    }

    // now merge the designated bins
    for (const BinPair & merge : mergeBins)
    {
        double sum   = 0;
        Int_t  count = 0;
        for (Int_t bin = merge.first; bin <= merge.second; ++bin)
        {
            sum += hist.GetBinContent(bin);
            ++count;
        }

        double avg = sum / count;

        for (Int_t bin = merge.first; bin <= merge.second; ++bin)
        {
            hist.SetBinContent( bin, avg );
        }
    }
}

////////////////////////////////////////////////////////////////////////////////
static Int_t OptBinPass1( const TH1D & hist )
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
        GetBinStats( *upRebin, meanEntries, varEntries, binWidth );

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
static Int_t OptBinPass2( const TH1D & hist, Int_t binPass1 )
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
        GetBinStats( *upRebin, meanEntries, varEntries, binWidth );

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
static Int_t OptBin3( const TH1D & hist )
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
        GetBinStats( *upRebin, meanEntries, varEntries, binWidth );

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
static TH1D * MakePDFHistExact( TNtupleD & data, Int_t nBins, Double_t xMin, Double_t xMax )
{
    TH1D * pHist = new TH1D( "PDF", "PDF", nBins, xMin, xMax );
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

    return pHist;
}

////////////////////////////////////////////////////////////////////////////////
static TH1D * MakePDFHist( TNtupleD & data, Double_t binWidth, Double_t xMin, Double_t xMax )
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

    //double integral = pHist->Integral();

    //pHist->Scale( 1.0 / integral );
    //pHist->Scale( 1.0 / integral, "width" );
    //pHist->Scale( 1.0, "width" );
    //pHist->Scale( 1.0, "width" );

    //LogMsgHistDump( *pHist );

    //double integral2 = pHist->Integral("width");

    return pHist;
}

////////////////////////////////////////////////////////////////////////////////
static void OptBinUnbinned1( const ModelCompare::Observable & observable, const TNtupleD & data )
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
static void OptBinUnbinned2( const ModelCompare::Observable & observable, const TNtupleD & data )
{
    LogMsgInfo( "\nOptimizing bin size for %hs", FMT_HS(data.GetName()) );
    LogMsgInfo( "------------------------------------------------------------");

    const Double_t xMin   = observable.xMin;
    const Double_t xMax   = observable.xMax;
    const Double_t xRange = xMax - xMin;

    // clone data so non-const methods can be used
    std::unique_ptr<TNtupleD> upFitData( (TNtupleD *)data.Clone() );

    TH1DUniquePtr upPDF;

    Int_t  nOptBins1    = 0;
    Int_t  nOptBins2    = 0;
    double bestMeasure1 = std::numeric_limits<double>::max();
    double bestMeasure2 = std::numeric_limits<double>::max();

    Int_t nBinInc = 5;
    Int_t nBinAdj = 100;

    for (Int_t nBins = 10; nBins <= 1000000; nBins += nBinInc)
    {
        if (nBins % nBinAdj == 0)
        {
            nBinInc *= 5;
            nBinInc = std::max( nBinInc, nBinAdj / 4 );
            nBinAdj *= 10;
        }

        double binWidth = xRange / nBins;

        upPDF.reset( MakePDFHistExact( *upFitData, nBins, xMin, xMax ) );

        //MergeZeroBins( *upPDF );

        double meanEntries, varEntries;
        GetBinStats( *upPDF, meanEntries, varEntries, 0.0 );

        double measure1 = GetErrorMeasure( meanEntries, varEntries, binWidth );
        double measure2 = GetCrossValidation( *upPDF, binWidth );

        if (measure1 < bestMeasure1)
        {
            bestMeasure1 = measure1;
            nOptBins1    = nBins;
        }

        if (measure2 < bestMeasure2)
        {
            bestMeasure2 = measure2;
            nOptBins2    = nBins;
        }

        LogMsgInfo( "n=%-6i  w=%-12g  m1=%-16g  m2=%-16g",
                    FMT_I(nBins), FMT_F(binWidth),
                    FMT_F(measure1), FMT_F(measure2) );
    }

    LogMsgInfo( "\nOptimal bins for %hs:", FMT_HS(data.GetName()) );
    LogMsgInfo( "MISE:  n=%-6i  w=%-12g  m1=%-16g", FMT_I(nOptBins1), FMT_F(xRange / nOptBins1), FMT_F(bestMeasure1) );
    LogMsgInfo( "Cross: n=%-6i  w=%-12g  m2=%-16g", FMT_I(nOptBins2), FMT_F(xRange / nOptBins2), FMT_F(bestMeasure2) );

    //exit(1);
}

////////////////////////////////////////////////////////////////////////////////
static void OptBinUnbinned3( const ModelCompare::ModelFile & model, const ModelCompare::Observable & observable, const TNtupleD & data )
{
    LogMsgInfo( "\nOptimizing bin size for %hs", FMT_HS(data.GetName()) );
    LogMsgInfo( "------------------------------------------------------------");

    const Double_t xMin   = observable.xMin;
    const Double_t xMax   = observable.xMax;
    const Double_t xRange = xMax - xMin;

    const Int_t    nBinMin     = 100;
    const Int_t    nBinMax     = 10000;
    const Int_t    nBinInit    = 1000;
    const Int_t    nBinStep    = 100;

#define MISE_OBJECTIVE   2
#define MISE_OBJ_NAME   "2"

#define MISE_NBINS
#define MISE_CROSS
//#define MISE_MERGE

#ifndef MISE_NBINS
    const Double_t binWidthMin  = xRange / nBinMax;
    const Double_t binWidthMax  = xRange / nBinMin;
    const Double_t binWidthInit = xRange / nBinInit;
    const Double_t binWidthStep = xRange / nBinInit - xRange / (nBinInit + nBinStep);
#endif

    // clone data so non-const methods can be used
    std::unique_ptr<TNtupleD> upFitData( (TNtupleD *)data.Clone() );

    ////////////////////////////
    // OBJECTIVE 1: Use a single PDF, with all bins even except the last

#if MISE_OBJECTIVE == 1
    auto ObjectiveFunc = [&](const double * par) -> double
    {
#ifndef MISE_NBINS
        Double_t fBinWidth = par[0];
        Double_t fBins     = xRange / fBinWidth;
#else
        Double_t fBins     = par[0];
        Double_t fBinWidth = xRange / fBins;
#endif

        TH1DUniquePtr upPDF( MakePDFHist( *upFitData, fBinWidth, xMin, xMax ) );

#ifdef MISE_MERGE
        MergeZeroBins( *upPDF );
#endif

        //double integral = upPDF->Integral();
        //upPDF->Scale( 1.0, "width" );

#ifndef MISE_CROSS
        double meanEntries, varEntries;
        GetBinStats( *upPDF, meanEntries, varEntries, fBinWidth );

        double measure = GetErrorMeasure( meanEntries, varEntries, fBinWidth );
        //double measure = 2 * meanEntries - varEntries;

        LogMsgInfo( "w=%-12g n=%-8g u=%-8g v=%-8g c=%g",
            FMT_F(fBinWidth), FMT_F(fBins),
            FMT_F(meanEntries), FMT_F(std::sqrt(varEntries)),
            FMT_F(measure) );
#else
        double measure = GetCrossValidation( *upPDF, fBinWidth );

        LogMsgInfo( "w=%-12g n=%-8g m=%g", FMT_F(fBinWidth), FMT_F(fBins), FMT_F(measure) );
#endif

        return measure;
    };
#endif

    ////////////////////////////
    // OBJECTIVE 2: Use two PDFs with bin numbers one apart, and
    // interpolate between the PDFs

#if MISE_OBJECTIVE == 2
    auto ObjectiveFunc = [&](const double * par) -> double
    {
#ifndef MISE_NBINS
        Double_t fBinWidth = par[0];
        Double_t fBins     = xRange / fBinWidth;
#else
        Double_t fBins     = par[0];
        Double_t fBinWidth = xRange / fBins;
#endif

        Int_t nBinsLow  = std::ceil( fBins);    // low bin-width
        Int_t nBinsHigh = std::floor(fBins);    // high bin-width

        Double_t binWidthLow  = xRange / nBinsLow;
        Double_t binWidthHigh = xRange / nBinsHigh;

        TH1DUniquePtr upPDFLow(  MakePDFHistExact( *upFitData, nBinsLow,  xMin, xMax ) );
        TH1DUniquePtr upPDFHigh( MakePDFHistExact( *upFitData, nBinsHigh, xMin, xMax ) );

#ifdef MISE_MERGE
        MergeZeroBins( *upPDFLow  );
        MergeZeroBins( *upPDFHigh );
#endif

#ifndef MISE_CROSS
        double meanEntries, varEntries;

        GetBinStats( *upPDFLow, meanEntries, varEntries, 0.0 );
        double measureLow  = GetErrorMeasure( meanEntries, varEntries, binWidthLow  );

        GetBinStats( *upPDFHigh, meanEntries, varEntries, 0.0 );
        double measureHigh = GetErrorMeasure( meanEntries, varEntries, binWidthHigh );
#else
        double measureLow  = GetCrossValidation( *upPDFLow,  binWidthLow  );
        double measureHigh = GetCrossValidation( *upPDFHigh, binWidthHigh );
#endif

#ifndef MISE_NBINS
        double measure;
        if (nBinsLow == nBinsHigh)
            measure = measureLow;
        else
        {
            double slope = (measureHigh - measureLow) / (binWidthHigh - binWidthLow);
            measure = (fBinWidth - binWidthLow) * slope + measureLow;
        }
#else
        double fractionLow   = fBins - nBinsHigh;
        double fractionHigh  = 1.0 - fractionLow;

        double measure = fractionLow * measureLow + fractionHigh * measureHigh;
#endif

        LogMsgInfo( "w=%-12g n=%-8g mL=%-8g mH=%-8g m=%g",
                    FMT_F(fBinWidth), FMT_F(fBins),
                    FMT_F(measureLow), FMT_F(measureHigh),
                    FMT_F(measure) );

        return measure;
    };
#endif

    //////

    ROOT::Fit::Fitter fitter;

    ROOT::Fit::FitConfig & fitConfig = fitter.Config();
    {
        std::vector<ROOT::Fit::ParameterSettings> params;

#ifndef MISE_NBINS
        ROOT::Fit::ParameterSettings par0( "BinWidth_" MISE_OBJ_NAME, binWidthInit, binWidthStep, binWidthMin, binWidthMax );
#else
        ROOT::Fit::ParameterSettings par0( "nBins_" MISE_OBJ_NAME,    nBinInit, nBinStep, nBinMin, nBinMax );
#endif
        params.push_back(par0);

        fitConfig.SetParamsSettings( params );

        //fitConfig.SetMinimizer("Minuit","MigradImproved");

        // run Hesse and Minos
        fitConfig.SetParabErrors(true);
        //fitConfig.SetMinosErrors(true);

        fitConfig.MinimizerOptions().SetStrategy(2);
        fitConfig.MinimizerOptions().SetPrintLevel(3);

        fitConfig.MinimizerOptions().SetTolerance(1);
    }

    bool fitOK = fitter.FitFCN( 1, ObjectiveFunc );
    if (!fitOK)
        LogMsgError( "Fit failed" );

    const ROOT::Fit::FitResult & fitResult = fitter.Result();

    fitResult.Print(std::cout);

    ROOT::Math::Minimizer * pMinimizer = fitter.GetMinimizer();

    // create a minimization scan
    if (pMinimizer)
    {
        std::unique_ptr<TGraph> upGraph( new TGraph( (Int_t)1000 ) );
        TGraph * pGraph = upGraph.get();  // alias

        std::string parName = fitResult.ParName(0);

        // do +/- 2x error scan (default range)
        double scanMin = fitResult.Parameter(0) - 2 * fitResult.ParError(0);
        double scanMax = fitResult.Parameter(0) + 2 * fitResult.ParError(0);

        // extend the range to +/- 1E-6
        //scanMin = std::min( scanMin, -1E-6 );
        //scanMax = std::max( scanMax,  1E-6 );

        // scan symmetric around 0.0 point
        //scanMin = std::min( scanMin, -scanMax );
        //scanMax = std::max( scanMax, -scanMin );

#ifndef MISE_NBINS
        scanMin = binWidthMin * 10; scanMax = binWidthMax;
#else
        scanMin = nBinMin; scanMax = nBinMax;
#endif

        unsigned int nStep = pGraph->GetN();
        bool scanResult = pMinimizer->Scan( 0, nStep, pGraph->GetX(), pGraph->GetY(), scanMin, scanMax );
        if (scanResult && (nStep != 0))
        {
            // success
            pGraph->Set( (Int_t)nStep );  // resize to actual steps

            std::string sName  = "min_" + std::string(model.modelName)  + "_" + std::string(observable.name)  + "_"               + parName;
            std::string sTitle =          std::string(model.modelTitle) + " " + std::string(observable.title) + ": fit min. for " + parName;

            pGraph->SetName(  sName .c_str() );
            pGraph->SetTitle( sTitle.c_str() );

            pGraph->GetXaxis()->SetTitle( parName.c_str() );
            pGraph->GetYaxis()->SetTitle( "MISE Cost"     );


            std::string graphFile( "optbin/graphs" );
#ifdef MISE_NBINS
            graphFile += "_nbins";
#endif
#ifdef MISE_CROSS
            graphFile += "_cross";
#endif
#ifdef MISE_MERGE
            graphFile += "_merge";
#endif
            graphFile += ".root";

            SaveGraph( graphFile.c_str(), *pGraph );
        }
    }

    //exit(1);
}

////////////////////////////////////////////////////////////////////////////////
void PrintStatistics( const ModelCompare::ModelFile & model, const ModelCompare::Observable & obs,
                      const TNtupleD & tuple )
{
    struct Statistics
    {
        double min;
        double max;
        double mean;
        double stddev;

        double P005;
        double P05;
        double P10;
        double P25;
        double P50;
        double P75;
        double P90;
        double P95;
        double P995;

        double U10;     // underflow max 10
        double U100;    // underflow max 100
        double O100;    // overflow  max 100
        double O10;     // overflow  max 10

        double n_sqrtn;
        double n_Rice;
        double h_Scott;
        double h_FD;

        Statistics( const std::vector<double> & v, double eventScale )
        {
            GetStats( v, min, max, mean, stddev );

            P005 = GetValueAtPos( v, v.size() * 0.005 );
            P05  = GetValueAtPos( v, v.size() * 0.05 );
            P10  = GetValueAtPos( v, v.size() * 0.10 );
            P25  = GetValueAtPos( v, v.size() * 0.25 );
            P50  = GetValueAtPos( v, v.size() * 0.50 );
            P75  = GetValueAtPos( v, v.size() * 0.75 );
            P90  = GetValueAtPos( v, v.size() * 0.90 );
            P95  = GetValueAtPos( v, v.size() * 0.95 );
            P995 = GetValueAtPos( v, v.size() * 0.995 );

            int event10  = int(10  * eventScale);
            int event100 = int(100 * eventScale);

            U10 = GetValueAtPos( v, event10 );
            O10 = GetValueAtPos( v, v.size() - 1 - event10 );

            U100 = GetValueAtPos( v, event100 );
            O100 = GetValueAtPos( v, v.size() - 1 - event100 );

            double n1_3 = pow( v.size(), 1.0 / 3.0 );

            n_sqrtn = std::sqrt(v.size());
            n_Rice  = 2.0 * n1_3;
            h_Scott = 3.5 * stddev / n1_3;
            h_FD    = 2.0 * (P75 - P25) / n1_3;
        }

        static double GetValueAtPos( const std::vector<double> & v, double pos )
        {
            int p1 = (int)std::floor(pos);
            int p2 = (int)std::ceil(pos);

            if (p1 == p2)
                return v[p1];

            double v1 = v[p1];
            double v2 = v[p2];

            double result = (pos - p1) * (v2 - v1) + v1;
            return result;
        }

        static void GetStats( const std::vector<double> & v, double & min, double & max, double & mean, double & stddev )
        {
            min  = std::numeric_limits<double>::max();
            max  = -min;

            double sum  = 0;
            double sum2 = 0;

            for (double x : v)
            {
                sum  += x;
                sum2 += x * x;
                min = std::min( min, x );
                max = std::max( max, x );
            }

            mean        = sum  / v.size();
            double var  = sum2 / v.size() - mean * mean;
            stddev      = std::sqrt( var );
        }
    };

    // clone data so non-const methods can be used
    std::unique_ptr<TNtupleD> upTuple( (TNtupleD *)tuple.Clone() );

    std::vector<double> dataAll;
    {
        Long64_t nEntries = upTuple->GetEntries();

        dataAll.reserve((size_t)nEntries);

        for (Long64_t entry = 0; entry < nEntries; ++entry)
        {
            //upTuple->LoadTree(entry);
            upTuple->GetEntry(entry);
            Double_t value = *upTuple->GetArgs();
            dataAll.push_back(value);
        }
    }

    std::sort( dataAll.begin(), dataAll.end() );

    std::vector<double> dataCut(dataAll);
    {
        auto itr = dataCut.begin();
        while (itr != dataCut.end())
        {
            double v = *itr;
            if ((v < obs.xMin) || !(v < obs.xMax))
            {
                itr = dataCut.erase(itr);
                continue;
            }
            ++itr;
        }
    }

    const double luminosity = 10;
    const double eventScale = model.crossSectionEvents / (model.crossSection * 1000 * luminosity);  // file events per luminosity event

    Statistics statAll( dataAll, eventScale );
    Statistics statCut( dataCut, eventScale );

    double rangeAll = statAll.max - statAll.min;
    double rangeCut = obs.xMax - obs.xMin;

    LogMsgInfo( "------------------------------------------------------------");
    LogMsgInfo( "Statistics for %hs %hs  [%g to %g]", FMT_HS(model.modelName), FMT_HS(obs.name),
                                                      FMT_F(obs.xMin), FMT_F(obs.xMax) );
    LogMsgInfo( "------------------------------------------------------------");

    LogMsgInfo( "Range: %g to %g  [%g to %g]",          FMT_F(statAll.min), FMT_F(statAll.max),
                                                        FMT_F(statCut.min), FMT_F(statCut.max) );

    LogMsgInfo( "U/O 10  range: %g to %g",              FMT_F(statAll.U10),  FMT_F(statAll.O10)  );
    LogMsgInfo( "U/O 100 range: %g to %g",              FMT_F(statAll.U100), FMT_F(statAll.O100) );

    LogMsgInfo( "99%% range: %g to %g  [%g to %g]",     FMT_F(statAll.P005), FMT_F(statAll.P995),
                                                        FMT_F(statCut.P005), FMT_F(statCut.P995) );

    LogMsgInfo( "90%% range: %g to %g  [%g to %g]",     FMT_F(statAll.P05), FMT_F(statAll.P95),
                                                        FMT_F(statCut.P05), FMT_F(statCut.P95) );

    LogMsgInfo( "Quarters: %g, %g, %g  [%g, %g, %g]",   FMT_F(statAll.P25), FMT_F(statAll.P50), FMT_F(statAll.P75),
                                                        FMT_F(statCut.P25), FMT_F(statCut.P50), FMT_F(statCut.P75) );

    LogMsgInfo( "IQR:    %g  [%g]",     FMT_F(statAll.P75 - statAll.P25),   FMT_F(statCut.P75 - statCut.P25) );

    LogMsgInfo( "Mean:   %g  [%g]",     FMT_F(statAll.mean),   FMT_F(statCut.mean) );
    LogMsgInfo( "Stddev: %g  [%g]",     FMT_F(statAll.stddev), FMT_F(statCut.stddev)  );

    LogMsgInfo( "--- optimal bins ---" );

    LogMsgInfo( "sqrt(n): n = %g  [%g]    w = %g  [%g]",    FMT_F(statAll.n_sqrtn),             FMT_F(statCut.n_sqrtn),
                                                            FMT_F(rangeAll/statAll.n_sqrtn),    FMT_F(rangeCut/statCut.n_sqrtn) );

    LogMsgInfo( "Rice:    n = %g  [%g]    w = %g  [%g]",    FMT_F(statAll.n_Rice),              FMT_F(statCut.n_Rice),
                                                            FMT_F(rangeAll/statAll.n_Rice),     FMT_F(rangeCut/statCut.n_Rice) );

    LogMsgInfo( "Scott:   n = %g  [%g]    w = %g  [%g]",    FMT_F(rangeAll/statAll.h_Scott),    FMT_F(rangeCut/statCut.h_Scott),
                                                            FMT_F(statAll.h_Scott),             FMT_F(statCut.h_Scott) );

    LogMsgInfo( "F-D:     n = %g  [%g]    w = %g  [%g]",    FMT_F(rangeAll/statAll.h_FD),       FMT_F(rangeCut/statCut.h_FD),
                                                            FMT_F(statAll.h_FD),                FMT_F(statCut.h_FD) );

}

////////////////////////////////////////////////////////////////////////////////
void OptBinUnbinned( const ModelCompare::ObservableVector & observables,
                     const ModelCompare::ModelFileVector &  models,
                     const char * cacheFileName )
{
    // disable reusing the same minuit object
    // this fixes the observed behavior/bug where a fit affected subsequent fits
    TMinuitMinimizer::UseStaticMinuit(false);

    std::vector<TupleVector> allData;   // allData[model][observable]

    FitEFT::LoadTupleData( observables, models, allData, cacheFileName );

    auto modelItr = models.cbegin();
    for (const TupleVector & tuples : allData)
    {
        auto model  = *modelItr++;
        auto obsItr = observables.cbegin();
        for (const TNtupleD * pTuple : tuples)
        {
            OptBinUnbinned2( *obsItr++, *pTuple );
            //OptBinUnbinned3( model, *obsItr++, *pTuple );
            //PrintStatistics( model, *obsItr++, *pTuple );
        }
    }
}

////////////////////////////////////////////////////////////////////////////////

}  // namespace OptBin
