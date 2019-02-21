// ALN Library
// Copyright (C) 2018 William W. Armstrong.
// 
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// Version 3 of the License, or (at your option) any later version.
// 
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
// 
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
// 
// For further information contact 
// William W. Armstrong
// 3624 - 108 Street NW
// Edmonton, Alberta, Canada  T6J 1B4

// alnlfnanalysis.cpp


#ifdef ALNDLL
#define ALNIMP __declspec(dllexport)
#endif

#include <aln.h>
#include "alnpriv.h"
#include <boost/math/special_functions/beta.hpp>

using namespace boost::math; 

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

static int ALNAPI ValidateALNLFNAnalysis(const ALN* pALN,
                                      const ALNDATAINFO* pDataInfo,
                                      const ALNCALLBACKINFO* pCallbackInfo,
                                      void*& pvAnalysis,
                                      int& nLFNStats);

#ifdef _DEBUG
static void DebugValidateALNLFNAnalysis(const ALN* pALN,
                                     const ALNDATAINFO* pDataInfo,
                                     const ALNCALLBACKINFO* pCallbackInfo,
                                     void*& pvAnalysis,
                                     int& nLFNStats);
#endif


// helper to sort data set by active LFN
struct LFNSORT
{
  ALNNODE* pActiveLFN;
  double* adblInputRow;
  double* adblOutputRow;
	double* adblResultRow; //Added August 28, 1999.
};
static int __cdecl CompareLFNs(const void* pElem1, const void* pElem2);

// helpers to calc T and F probability functions
static double CalcProbF(double dblF, double dblV1, double dblV2);
static double CalcProbT(double dblT, double dblV);

// helpers to calc error, regression, and total sum of squares
static void CalcSquares(double* adblDesired, double* adblResult,
                        int nStart, int nEnd,
                        double& dblESS, double& dblRSS, double& dblTSS);

// helpers to allocate / deallocate analysis storage
struct LFNINFO
{
  ALNNODE* pLFN;
  LFNSTATS LFNStats;
  int nWStats;
  LFNWEIGHTSTATS aWStats[1];   // actually, nWStats elements
};

struct LFNANALYSIS
{
  int nLFNInfoSize;     // size of each LFNINFO element
  int nLFNInfo;         // nuumber of LFNINFO elements
  LFNINFO aLFNInfo[1];  // actually, nLFNInfo elements
};

static LFNANALYSIS* AllocLFNAnalysis(const ALN* pALN);
static void FreeLFNAnalysis(LFNANALYSIS* pLFNAnalysis);

inline
LFNINFO* GetLFNInfo(LFNANALYSIS* pLFNAnalysis, int nLFN)
{
  return (LFNINFO*)(((char*)(pLFNAnalysis->aLFNInfo)) + (nLFN * pLFNAnalysis->nLFNInfoSize));
}


// analyze LFNs
ALNIMP int ALNAPI ALNLFNAnalysis(const ALN* pALN,
                                 const ALNDATAINFO* pDataInfo,
                                 const ALNCALLBACKINFO* pCallbackInfo,
                                 void*& pvAnalysis,
                                 int& nLFNStats)
{
	
  int nReturn = ValidateALNLFNAnalysis(pALN, pDataInfo, pCallbackInfo,
                                       pvAnalysis, nLFNStats);
  if (nReturn != ALN_NOERROR)
    return nReturn;
  
  #ifdef _DEBUG
    DebugValidateALNLFNAnalysis(pALN, pDataInfo, pCallbackInfo,
                                pvAnalysis, nLFNStats);
  #endif

  // arrays

  // evaltree results
  double* adblResult = NULL;
  ALNNODE** apActiveLFNs = NULL;
  double* adblInput = NULL;
  double* adblOutput = NULL;
  LFNSORT* aLFNSort = NULL;

  // covariance calculation (on single block of inputs and outputs for one LFN)
  double* adblX = NULL;      // RMA based ALN input vectors (nRows * nCols); nCols = nDim
  double* adblY = NULL;      // RMA based desired column vector (nRows); nRows = #points this LFN
  double* adblC = NULL;      // RMA covariance matrix (nCols * nCols)
  double* adblA = NULL;      // RMA fitted parameter vector (nCols); adblA can be computed by SVD
  double* adblS = NULL;      // RMA std dev vector (nRows); output of SVD using adblX
  double* adblU = NULL;      // RMA U matrix (nRows * nCols); output of SVD using adblX
  double* adblV = NULL;      // RMA V matrix (nCols * nCols); output of SVD using adblX
  double* adblW = NULL;      // RMA W matrix (nCols); singular values from the SVD of matrix adblX 
  double* adblALN = NULL;    // RMA based ALN result column vector (nRows); the given fit to the data

	// idea: Use  std::vector< std::vector<Indices> > indices(layers, std::vector<Indices>(corrections));

  // LFN analysis
  LFNANALYSIS* pLFNAnalysis = NULL;

  try
  {
    // allocate analysis struct
    pLFNAnalysis = AllocLFNAnalysis(pALN);
    if (pLFNAnalysis == NULL)
      ThrowALNMemoryException();

    // see how many points there are
    int nPoints = pDataInfo->nPoints;
    
    // get dim... we use all the input variables of the ALN,
    // plus an explicit bias variable always equal to one,
    int nDim = pALN->nDim;
		int nCols = nDim - 1;
		int nOutput = pALN->nOutput;

    // allocate arrays
    adblResult = new double[nPoints];
    apActiveLFNs = new ALNNODE*[nPoints];
    adblInput = new double[nPoints * nDim];
    adblOutput = new double[nPoints];
    aLFNSort = new LFNSORT[nPoints];
    adblX = new double[nPoints * nCols];
		adblY = new double[nPoints];
    adblC = new double[nCols * nCols];
    adblA = new double[nCols];
    adblS = new double[nPoints];
    adblU = new double[nPoints * nPoints];
    adblV = new double[nCols * nCols];
    adblW = new double[nCols];
    adblALN = new double[nPoints];


    if (!(adblResult && apActiveLFNs && adblInput && adblOutput && 
          aLFNSort && adblX && adblALN && adblC && adblA && adblS &&
          adblU && adblV && adblW && adblY))
    {
      ThrowALNMemoryException();
    }

    // init S
    for (int i = 0; i < nPoints; i++)
    {
      adblS[i] = 1.0; // start off with std dev 1.0 in each axis
    }
    
    // evaluate on data
    int nStart, nEnd;
    nReturn = EvalTree(pALN->pTree, pALN, pDataInfo, pCallbackInfo, 
                       adblResult, &nStart, &nEnd, FALSE,
                       apActiveLFNs, adblInput, adblOutput);
    if (nReturn != ALN_NOERROR)
    {
      ThrowALNException();
    }

    // sort input/output arrays by responsible LFN
		// pointers in LFNSORT entries refer to  *unsorted* arrays
		// In the adblInput array, the nOutput components have been set to 1.0 by EvalTree.
    for (int i = nStart; i <= nEnd; i++)
    {
      LFNSORT& lfnsort = aLFNSort[i];
      lfnsort.pActiveLFN = apActiveLFNs[i];
      lfnsort.adblInputRow = adblInput + (i * nDim);
			// WWA I think it would be better to use adblOutputCol and adblResultCol as names:
			// The way the SVD is presented as
			// Matrix of input rows x parameter column = Result column.
			// The Output column is originally in the matrix of input rows
			// as the column of desired values (later replaced during evaluation of the DTREE)
      lfnsort.adblOutputRow = adblOutput+i; // desired
			adblOutput[i];
			lfnsort.adblResultRow = adblResult+i; // ALN output
    }
    // sort the LFN index structure
    qsort(aLFNSort + nStart, nEnd - nStart + 1, sizeof(LFNSORT), CompareLFNs);
  
    // calc covariance matrix for each block of data points
    int nBlockStart = nStart; // Monroe had 0, changed by wwa Aug 21 1999
    ALNNODE* pBlockLFN = aLFNSort[nStart].pActiveLFN;
    nLFNStats = 0;
    for (int i = nStart; i <= nEnd; i++)
    {
      if (i == nEnd || aLFNSort[i + 1].pActiveLFN != pBlockLFN)
      {
        // get pointer to current LFN info
        LFNINFO* pLFNInfo = GetLFNInfo(pLFNAnalysis, nLFNStats);

        // increment LFN stat count
        nLFNStats += 1;

        // end of LFN block... copy data to scratch input and output
        int nBlockPoints = i - nBlockStart + 1;
        for (int j = 0; j < nBlockPoints; j++)
        {
          // input row (with 1.0 in nOutput component)  WWA we don't copy the 1.0 into X
          memcpy(adblX + j * nCols,                             
                 aLFNSort[nBlockStart + j].adblInputRow, 
                 nCols * sizeof(double)); 
          // ALN output
          adblALN[j] = *(aLFNSort[nBlockStart + j].adblResultRow); //contains ALN output

				  // desired output
					adblY[j] = *(aLFNSort[nBlockStart + j].adblOutputRow); // contains desired
        }
				// The result of t4 above is correct
        // calc covariance

				// NB The actual ALN outputs play no role in the covariance calculation
				// Y below is the desired output
				// The ALN outputs could be compared to the outputs using parameter vector A found below
        double dblChiSq = 0;
        BOOL bSuccess = CalcCovariance(nCols,
                                       nBlockPoints, 
                                       adblX,
                                       adblY,
                                       adblC,
                                       adblA,
                                       adblS,
                                       adblU,
                                       adblV,
                                       adblW,
                                       dblChiSq);
        // redo the Y copy,just in case it changed WWA because of a memcpy that had nDim instead of nCols in
				// some places, the Y values for j = 0 to j = 4 were wrong, both inside and leaving the calc routine.
				// After correction, this part is no longer necessary.
        // for (int j = 0; j < nBlockPoints; j++)
				//           
				  // desired 
					// 	adblY[j] = *(aLFNSort[nBlockStart + j].adblOutputRow); 

          // ALN output
					//    adblALN[j] = *(aLFNSort[nBlockStart + j].adblResultRow); 
					// }

        // calc RSS, ESS, TSS
        double dblRSS, dblESS, dblTSS;
        CalcSquares(adblY, adblALN, 0, nBlockPoints -1, dblESS, dblRSS, dblTSS);
        
        pLFNInfo->LFNStats.dblRSS = dblRSS;
        pLFNInfo->LFNStats.dblESS = dblESS;
        
        // calc R2 for the LFN
        pLFNInfo->LFNStats.dblR2 = dblRSS / dblTSS;

        // calc degrees of freedom
        pLFNInfo->LFNStats.dblDF = nBlockPoints - nDim; 

				pLFNInfo->LFNStats.dblSEE = NAN; // set to NAN if we don't have dblDF>0
        pLFNInfo->LFNStats.dblF   = NAN;
        pLFNInfo->LFNStats.dblFp  = NAN;
        if (pLFNInfo->LFNStats.dblDF > 0)
        {
          ASSERT(nDim >= 1);
          double dblV1 = nDim - 1;
          double dblV2 = pLFNInfo->LFNStats.dblDF;

          // calc SEE "standard error of estimate" for the LFN
					pLFNInfo->LFNStats.dblSEE = sqrt(dblESS / dblV2);

					// calc F and corresponding p value for the LFN
          pLFNInfo->LFNStats.dblF = (dblRSS / dblV1) / (dblESS / dblV2);
          pLFNInfo->LFNStats.dblFp = CalcProbF(pLFNInfo->LFNStats.dblF, dblV1, dblV2);
        }
              
        // set the node
        pLFNInfo->pLFN = pBlockLFN;

        // set number of weight stats
        pLFNInfo->nWStats = nDim;

        // calc T, and corresponding p values for each LFN weight
        const double* adblLFNW = LFN_W(pBlockLFN);
        double dblBias = *adblLFNW++;  // skip past bias
        double dblV = pLFNInfo->LFNStats.dblDF;
				double dblLFNweight; // used for ALN weights to be compared to SVD result
        for (int k = 0; k < nDim; k++)
				{
          // get LFNWEIGHTSTATS pointer
          LFNWEIGHTSTATS* pWStat = &(pLFNInfo->aWStats[k]);

          // set weights from LFN!!!
          if (k == pALN->nOutput)
          dblLFNweight = dblBias;		// bias replaces output weight for stats
          else
          dblLFNweight = adblLFNW[k];	// kth input weight

				 // the weights are now obtained from SVD and are in adblA
				 // (as Monroe suggested doing). 
				 // dblBias should approximate adblA[k]if k = nOutput 
				 pWStat->dblW = adblA[k];

          // set T stat
          pWStat->dblSEw = 0;
          pWStat->dblT = NAN;
          pWStat->dblTp = 1.0;
          if (pLFNInfo->LFNStats.dblDF > 0)
          {
            // calc standard error
            double dblCkk = adblC[k*nDim + k];  // diagonal element Ckk 

            pWStat->dblSEw = sqrt(dblCkk) * pLFNInfo->LFNStats.dblSEE;

            if (pWStat->dblSEw != 0.0)
            {
              // calc T, Tp
							// WWA the following line didn't include the term adblLFNW[k],
							// thus making it refer to the k-th weight, not to the
							// difference of ALN computed and SVD computed weights.
							// In the original interpretation, that would decide whether
							// the weight could just as well be zero, ie not significant.
							// Now we want to use it to say if the difference of ALN
							// weight from the SVD computed one is significant.  If ALN
							// does well, then all of the T's should be less than 2.0 or 3.0 
							// and the probabilities Tp should be not to far from 1.0.
							// Is this OK?
              pWStat->dblT = (dblLFNweight - pWStat->dblW) / pWStat->dblSEw;
              pWStat->dblTp = CalcProbT(fabs(pWStat->dblT), dblV);
            }
          }
        }
        
        // start of new LFN block
        if (i < nEnd)
        {
          nBlockStart = i + 1;  
          pBlockLFN = aLFNSort[nBlockStart].pActiveLFN;
        }
      }   
    }          
	}

  catch (CALNUserException* e)	  // user abort exception
	{
		nReturn = ALN_USERABORT;
    e->Delete();
	}
	catch (CALNMemoryException* e)	// memory specific exceptions
	{
		nReturn = ALN_OUTOFMEM;
    e->Delete();
	}
	catch (CALNException* e)	      // anything other exception we recognize
	{
		nReturn = ALN_GENERIC;
    e->Delete();
	}
	catch (...)		                  // anything else, including FP errs
	{
		nReturn = ALN_GENERIC;
	}
	// deallocate mem
	delete[] adblALN;
	delete[] adblW;
	delete[] adblV;
	delete[] adblU; // this is where the heap corruption happens
	delete[] adblS;
	delete[] adblA;
	delete[] adblC;
	delete[] adblY;
	delete[] adblX;
	delete[] aLFNSort;
	delete[] adblOutput;
	delete[] adblInput;
	delete[] apActiveLFNs;
	delete[] adblResult;

  // set the analysis value
  if (nReturn != ALN_NOERROR && pLFNAnalysis != NULL)
  {
    FreeLFNAnalysis(pLFNAnalysis);
  }
  else
  {
    pvAnalysis = (void*)pLFNAnalysis;
  }
  
  return nReturn;
}

// get LFN stats
ALNIMP int ALNAPI ALNLFNStats(void* pvAnalysis, 
                              int nLFNStat,
                              LFNSTATS* pLFNStats,
                              int* pnWeightStats,
                              ALNNODE** ppLFN)
{
  if (pvAnalysis == NULL)
    return ALN_GENERIC;

  // make sure LFN index in range
  LFNANALYSIS* pLFNAnalysis = (LFNANALYSIS*)pvAnalysis;
  if (nLFNStat < 0 || nLFNStat >= pLFNAnalysis->nLFNInfo)
    return ALN_GENERIC;

  // get LFNINFO
  LFNINFO* pLFNInfo = GetLFNInfo(pLFNAnalysis, nLFNStat);
  
  // copy stat values
  if (pLFNStats != NULL)
    memcpy(pLFNStats, &(pLFNInfo->LFNStats), sizeof(LFNSTATS));
  
  if (pnWeightStats != NULL)
    *pnWeightStats = pLFNInfo->nWStats;

  if (ppLFN != NULL)
    *ppLFN = pLFNInfo->pLFN;

  return ALN_NOERROR;
}

// get LFN weight stats
ALNIMP int ALNAPI ALNLFNWeightStats(void* pvAnalysis, 
                                    int nLFNStat,
                                    int nWeightStat,
                                    LFNWEIGHTSTATS* pLFNWeightStats)
{
  if (pvAnalysis == NULL || pLFNWeightStats == NULL)
    return ALN_GENERIC;

  // make sure LFN index in range
  LFNANALYSIS* pLFNAnalysis = (LFNANALYSIS*)pvAnalysis;
  if (nLFNStat < 0 || nLFNStat >= pLFNAnalysis->nLFNInfo)
    return ALN_GENERIC;

  // get LFNINFO
  LFNINFO* pLFNInfo = GetLFNInfo(pLFNAnalysis, nLFNStat);

  // make sure W index is in range
  if (nWeightStat < 0 || nWeightStat >= pLFNInfo->nWStats)
    return ALN_GENERIC;

  // copy stat values
  memcpy(pLFNWeightStats, 
         &(pLFNInfo->aWStats[nWeightStat]), 
         sizeof(LFNWEIGHTSTATS));

  return ALN_NOERROR;
}

// free analysis memory
ALNIMP int ALNAPI ALNLFNFreeAnalysis(void* pvAnalysis)
{
  if (pvAnalysis != NULL)
    FreeLFNAnalysis((LFNANALYSIS*)pvAnalysis);
    
  return ALN_NOERROR;
}


// validate params
static int ALNAPI ValidateALNLFNAnalysis(const ALN* pALN,
                                      const ALNDATAINFO* pDataInfo,
                                      const ALNCALLBACKINFO* pCallbackInfo,
                                      void*& pvAnalysis,
                                      int& nLFNStats)
{
  int nReturn = ValidateALNDataInfo(pALN, pDataInfo, pCallbackInfo);
  if (nReturn != ALN_NOERROR)
    return nReturn;
  
  //if (pvAnalysis == NULL || nLFNStats != 0)
  //  return ALN_GENERIC;

  return ALN_NOERROR;
}

// debug version ASSERTS if bad params
#ifdef _DEBUG
static void DebugValidateALNLFNAnalysis(const ALN* pALN,
                                     const ALNDATAINFO* pDataInfo,
                                     const ALNCALLBACKINFO* pCallbackInfo,
                                     void*& pvAnalysis,
                                     int& nLFNStats)
{
  DebugValidateALNDataInfo(pALN, pDataInfo, pCallbackInfo);
}
#endif

static int __cdecl CompareLFNs(const void* pElem1, const void* pElem2)
{
  LFNSORT& lfnsort1 = *(LFNSORT*)pElem1;
  LFNSORT& lfnsort2 = *(LFNSORT*)pElem2;

  if (lfnsort1.pActiveLFN < lfnsort2.pActiveLFN)
    return -1;
  else if (lfnsort1.pActiveLFN > lfnsort2.pActiveLFN)
    return 1;

  return 0;
}

// helpers to calc T and F probability functions
static double CalcProbF(double dblF, double dblV1, double dblV2)
{
  // F-Distribution probability function - not two tailed as in
  // Press et al p229, p619.  Here it is used in regression and is 1-tailed.
	// It is testing H0: all weights = 0 (except for the bias weight).
  
	return ibeta(dblV2 / 2.0, dblV1 / 2.0, dblV2 / (dblV2 + (dblV1 * dblF))); //ibeta ??

//
}

static double CalcProbT(double dblT, double dblV)
{
  // T-Distribution probability function - single tailed
  // see Press et al p228, p616
  return ibeta(dblV / 2.0, 0.5, dblV / (dblV + (dblT * dblT)));  //ibetac??
}

// helpers to calc error, regression, and total sum of squares
static void CalcSquares(double* adblDesired, double* adblResult,
                        int nStart, int nEnd,
                        double& dblESS, double& dblRSS, double& dblTSS)
{
  ASSERT(adblDesired && adblResult);
  ASSERT(nStart <= nEnd);

  dblESS = 0;
  dblTSS = 0;
  dblRSS = 0;

  double dblMeanDes = 0;
  int nElem = nEnd - nStart + 1;
  ASSERT(nElem >= 1);

  // calc error sum of squares and average of desired
  for (int i = nStart; i <= nEnd ;i++)
  {
    dblMeanDes += adblDesired[i];
    double dblError = adblResult[i] - adblDesired[i];
    dblESS += dblError * dblError;
  }
  dblMeanDes /= nElem;

  // calc total sum of squares
  for (int i = nStart; i <= nEnd; i++)
  {
    double dblError = adblDesired[i] - dblMeanDes; // Monroe had adblResult[i] - dblMeanDes;
    dblTSS += dblError * dblError;
  }

  // calc regression sum of squares
  dblRSS = dblTSS - dblESS;
}

static LFNANALYSIS* AllocLFNAnalysis(const ALN* pALN)
{
  // calc number of LFNs
  int nTotalLFNs = 0;
  int nAdaptedLFNs = 0;
  CountLFNs(pALN->pTree, nTotalLFNs, nAdaptedLFNs);

  if (nTotalLFNs == 0)
    return NULL;

  // calc number of bytes needed for entire struct
  int nWStatBytes = pALN->nDim * sizeof(LFNWEIGHTSTATS);
  int nLFNInfoBytes = sizeof(LFNINFO) - sizeof(LFNWEIGHTSTATS) + nWStatBytes;
  int nTotalBytes = sizeof(LFNANALYSIS) - sizeof(LFNINFO) + (nTotalLFNs * nLFNInfoBytes);

  // alloc block
  LFNANALYSIS* pLFNAnalysis = (LFNANALYSIS*)calloc(1, nTotalBytes);               
  if (pLFNAnalysis == NULL)
    return NULL;

  // init block
  pLFNAnalysis->nLFNInfoSize = nLFNInfoBytes;
  pLFNAnalysis->nLFNInfo = nTotalLFNs;
  
  return pLFNAnalysis;
}

static void FreeLFNAnalysis(LFNANALYSIS* pLFNAnalysis)
{
  if (pLFNAnalysis != NULL)
    free(pLFNAnalysis);
}

