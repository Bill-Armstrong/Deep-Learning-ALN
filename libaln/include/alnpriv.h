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
// alnpriv.h

#ifndef __ALNPRIV_H__
#define __ALNPRIV_H__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <malloc.h>
#include <limits>
#define ALNAPI __stdcall


// SunOS math.h has conflict with C++ rtti and exception handling 
// when using GNU C
#ifdef sun
extern "C" {
#endif

#include <math.h>

#ifdef sun
}
#endif

#ifdef __GNUC__
#include <typeinfo>
#endif

#ifdef sun
#include <memory.h>
#endif

// id string
static char szIDString[] = "ALN Library "
                           "Copyright (C) 2018 William W. Armstrong "
                           "See source code for LGPL license";

///////////////////////////////////////////////////////////////////////////////
// diagnostic functions

#ifndef TRACE
#define TRACE ALNTRACE
#endif

#ifndef ASSERT
#define ASSERT ALNASSERT
#endif

#ifndef VERIFY
#define VERIFY ALNVERIFY
#endif

/////////////////////////////////////////////////////////////////////////////
// extra types

#ifndef BOOL
typedef int BOOL;
#endif

#undef TRUE
#undef FALSE
#define TRUE 1
#define FALSE 0


/////////////////////////////////////////////////////////////////////////////
// Minimum and maximum templates

template<class T> inline
const T& max(const T& x, const T& y)
  { return ((x > y) ? x : y); }

template<class T> inline
const T& min(const T& x, const T& y)
  { return ((x < y) ? x : y); }


///////////////////////////////////////////////////////////////////////////////
// C++ exception classes for internal ALN library use
// preferred usage:
//          throw new CALNException(TRUE, "exception cause");
//    or    throw new CALNException();
//    or    throw new CALNMemoryException();
//    or    throw new CALNUserException();
// make sure to call Delete() on exception object in handler 
// if you do not intend to re-throw the exception

class CALNException
{
private:
  BOOL m_bAutoDelete;
	char m_szReason[256];

public:
	CALNException(BOOL bAutoDelete = TRUE, const char* pszReason = NULL);

	virtual ~CALNException() {};
  virtual void Delete();

	const char* GetReason() const
		{ return m_szReason; }
};
void ALNAPI ThrowALNException();
  // throws CALNException*
  // make sure to call Delete() on exception object in handler

class CALNMemoryException : public CALNException
{
public:
	CALNMemoryException(BOOL bAutoDelete = TRUE);
	virtual ~CALNMemoryException() {};
};
void ALNAPI ThrowALNMemoryException();
  // throws CALNMemoryException*
  // make sure to call Delete() on exception object in handler

class CALNUserException : public CALNException
{
public:
	CALNUserException(BOOL bAutoDelete = TRUE);
	virtual ~CALNUserException() {};
};
void ALNAPI ThrowALNUserException();
  // throws CALNUserException*
  // make sure to call Delete() on exception object in handler


///////////////////////////////////////////////////////////////////////////////
// data handling routines

int ALNAPI ValidateALNDataInfo(const ALN* pALN,
                               const ALNDATAINFO* pDataInfo,
                               const ALNCALLBACKINFO* pCallbackInfo);

// debug version ASSERTS if bad params
#ifdef _DEBUG
void ALNAPI DebugValidateALNDataInfo(const ALN* pALN,
                                     const ALNDATAINFO* pDataInfo,
                                     const ALNCALLBACKINFO* pCallbackInfo);
#endif


// calculate start and end points of data set given varinfo deltas
void ALNAPI CalcDataEndPoints(long& nStart, long& nEnd, 
                              const ALN* pALN,
                              const ALNDATAINFO* pDataInfo);

// allocate and init data column base pointers
const double** ALNAPI AllocColumnBase(int nStart, 
                                      const ALN* pALN,
                                      const ALNDATAINFO* pDataInfo);

// free data column base pointers
void ALNAPI FreeColumnBase(const double**);

// fill input vector
void ALNAPI FillInputVector(const ALN* pALN,
                            double* adblX, 
                            int nPoint,
                            int nStart,
                            const double** apdblBase,
                            const ALNDATAINFO* pDataInfo,
                            const ALNCALLBACKINFO* pCallbackInfo);

///////////////////////////////////////////////////////////////////////////////
// eval routines

// evaluation cutoff context info
struct CEvalCutoff
{
  int bMin;
  double dblMin;
  int bMax;
  double dblMax;

  CEvalCutoff()
    { bMin = bMax = FALSE; dblMin = dblMax = 0; }
};

// build cutoff route up tree
void ALNAPI BuildCutoffRoute(ALNNODE* pNode);

// check if value meets cutoff criteria for pNode...
// assumes that cutoff bounds have already been loosened for child evaluation
BOOL ALNAPI Cutoff(double dbl, const ALNNODE* pNode, CEvalCutoff& cutoff, 
                   double dbl4SE);

// CCutoffInfo struct to store last known LFN and value for a pattern
struct CCutoffInfo
{
  ALNNODE* pLFN;
  double dblValue;
};

// LFN specific eval - returns distance to surface
//  - non-destructive, ie, does not change ALN structure
double ALNAPI CutoffEvalLFN(const ALNNODE* pNode, const ALN* pALN, 
                            const double* adblX, ALNNODE** ppActiveLFN);

// minmax node specific eval - returns distance to surface
//  - non-destructive, ie, does not change ALN structure
// NOTE: cutoff always passed on stack!
double ALNAPI CutoffEvalMinMax(const ALNNODE* pNode, const ALN* pALN, 
                             const double* adblX, CEvalCutoff cutoff, 
                             ALNNODE** ppActiveLFN);

// generic cutoff eval
inline double CutoffEval(const ALNNODE* pNode, const ALN* pALN, 
                         const double* adblX, CEvalCutoff cutoff, 
                         ALNNODE** ppActiveLFN)
{
  return (pNode->fNode & NF_LFN) ? 
           CutoffEvalLFN(pNode, pALN, adblX, ppActiveLFN) : 
           CutoffEvalMinMax(pNode, pALN, adblX, cutoff, ppActiveLFN);
}

// cutoff eval with cutoff info
double ALNAPI CutoffEval(const ALNNODE* pNode, const ALN* pALN, 
                         const double* adblX, CCutoffInfo* pCutoffInfo, 
                         ALNNODE** ppActiveLFN);

// LFN specific eval - returns distance to surface
// - adaptive, ie, will change ALN structure
double ALNAPI AdaptEvalLFN(ALNNODE* pNode, ALN* pALN, const double* adblX,
                           ALNNODE** ppActiveLFN);

// minmax node specific adapt eval - returns distance to surface
// - adaptive, ie, will change ALN structure
// NOTE: cutoff always passed on stack!
double ALNAPI AdaptEvalMinMax(ALNNODE* pNode, ALN* pALN, const double* adblX,
	CEvalCutoff cutoff, ALNNODE** ppActiveLFN);
                          
// generic adapt eval
inline double AdaptEval(ALNNODE* pNode, ALN* pALN, 
                        const double* adblX, CEvalCutoff cutoff, 
                        ALNNODE** ppActiveLFN)
{
  return (pNode->fNode & NF_LFN) ? 
           AdaptEvalLFN(pNode, pALN, adblX, ppActiveLFN) : 
           AdaptEvalMinMax(pNode, pALN, adblX, cutoff, ppActiveLFN);
}

// adapt eval with cutoff info
double ALNAPI AdaptEval(ALNNODE* pNode, ALN* pALN, const double* adblX, 
                        CCutoffInfo* pCutoffInfo, ALNNODE** ppActiveLFN);



#ifdef _DEBUG
double ALNAPI DebugEvalMinMax(const ALNNODE* pNode, const ALN* pALN, 
                            const double* adblX, ALNNODE** ppActiveLFN);

// generic debug eval
inline double DebugEval(const ALNNODE* pNode, const ALN* pALN, 
                        const double* adblX, ALNNODE** ppActiveLFN)
{
  return (pNode->fNode & NF_LFN) ? 
           CutoffEvalLFN(pNode, pALN, adblX, ppActiveLFN) : 
           DebugEvalMinMax(pNode, pALN, adblX, ppActiveLFN);
}
#endif

// calculate active child, response of active child, and distance
int ALNAPI CalcActiveChild(double& dblRespActive, double& dblDistance, 
                           double dbl0, double dbl1, const ALNNODE* pNode, 
                           double dblSE, double dbl4SE, double dblOV16SE);

// evaluate a tree on a dataset 
// if bErrorResults is true, then errors are returned in adblResults 
//   instead of actual values
// if adblInput is non-NULL, then it must have enough room to store
//   all the input vectors passed to the ALN... the vectors have nDim
//   elements... the output variable element is replaced by the bias
//   element with value == 1.0
// if adblOutput is non-NULL, then it must have enough room to store
//   all the desired output values... the application is responsible
//   for providing correct output values through pDataInfo or during
//   the callback
int ALNAPI EvalTree(const ALNNODE* pNode, const ALN* pALN,
                    const ALNDATAINFO* pDataInfo,
                    const ALNCALLBACKINFO* pCallbackInfo,
                    double* adblResult, int* pnStart, int* pnEnd,
                    BOOL bErrorResults = FALSE,
                    ALNNODE** apActiveLFNs = NULL,
                    double* adblInput = NULL,
                    double* adblOutput = NULL);


///////////////////////////////////////////////////////////////////////////////
// monotonicity checking and support routines

// check monotonicity of a subtree
int ALNAPI CheckMonotonicity(const ALNNODE* pNode, const ALN* pALN, int nVar);



///////////////////////////////////////////////////////////////////////////////
// Training and support routines

// used to prepare ALN for training
BOOL ALNAPI PrepALN(ALN* pALN);

// set smoothing epsilon for a region
inline void ALNAPI SetSmoothingEpsilon(ALNREGION* pRegion)
{
  ASSERT(pRegion);
  double dblSE = pRegion->dblSmoothEpsilon;
  if (dblSE > 0.0)
  {
    pRegion->dbl4SE = 4 * dblSE;
    pRegion->dblOV16SE = 1.0 / (16.0 * dblSE);
  }
  else
  {
    pRegion->dbl4SE = 0.0;
    pRegion->dblOV16SE = 0.0;
  }
} 

// used to count number of LFNs in an ALN
void ALNAPI CountLFNs(const ALNNODE* pNode, int& nTotal, int& nAdapted);

// init any uninitialized LFN's
void ALNAPI InitLFNs(ALNNODE* pNode, ALN* pALN, const double* adblX);

// used to reset resp counters, and other stats (alntrain.cpp)
void ALNAPI ResetCounters(ALNNODE* pNode, ALN* pALN, 
                          BOOL bMarkAsUseful = FALSE);

// calcs RMS err on data set
double ALNAPI DoCalcRMSError(const ALN* pALN,
                             const ALNDATAINFO* pDataInfo,
                             const ALNCALLBACKINFO* pCallbackInfo);

// callback - throws CALNUserException if callback returns 0
inline BOOL CanCallback(int nCode, ALNNOTIFYPROC pfnNotifyProc,
                        int nNotifyMask)
{
  return (pfnNotifyProc && (nNotifyMask & nCode));
}

inline void Callback(const ALN* pALN, int nCode, void* pParam,
                     ALNNOTIFYPROC pfnNotifyProc, void* pvData)
{
  if (!(*(pfnNotifyProc))(pALN, nCode, pParam, pvData))
  {
    ThrowALNUserException();
  }
}

// jitter
void ALNAPI Jitter(ALN* pALN, double* adblX);

// shuffle
void ALNAPI Shuffle(int nStart, int nEnd, int* anShuffle);

// training context info
typedef struct tagTRAINDATA
{
  int nNotifyMask;
  ALNNOTIFYPROC pfnNotifyProc;
  void* pvData;
  double dblLearnRate;
  double dblGlobalError;
  // Potentially add Dmitri's flywheel
} TRAINDATA;


///////////////////////////////////////////////////////////////////////////////
// node adaptation routines

// constraint retrieval 
ALNCONSTRAINT* ALNAPI GetVarConstraint(int nRegion, const ALN* pALN,
                                       int nVar);

// minmax specific adapt
void ALNAPI AdaptMinMax(ALNNODE* pNode, ALN* pALN, const double* adblX, 
                      double dblResponse, BOOL bUsefulAdapt, 
                      const TRAINDATA* ptdata);

// LFN specific adapt
void ALNAPI AdaptLFN(ALNNODE* pNode, ALN* pALN, const double* adblX, 
                     double dblResponse, BOOL bUsefulAdapt, 
                     const TRAINDATA* ptdata);

// generic adapt routine
inline void Adapt(ALNNODE* pNode, ALN* pALN, const double* adblX, 
                  double dblResponse, BOOL bUsefulAdapt, 
                  const TRAINDATA* ptdata)
{
  (pNode->fNode & NF_LFN) ? 
	  AdaptLFN(pNode, pALN, adblX, dblResponse, bUsefulAdapt, ptdata) : 
    AdaptMinMax(pNode, pALN, adblX, dblResponse, bUsefulAdapt, ptdata);
}

// split routines
int ALNAPI SplitLFN(ALN* pALN, ALNNODE* pNode);

// find LFNs to split according to number of points on them and the F-test
void ALNAPI FindSplitLFN(ALN* pALN);


///////////////////////////////////////////////////////////////////////////////
// DTREE conversion routines

#define DTREE_MINDEPTH    1
#define DTREE_MAXDEPTH    30

DTREE* ALNAPI BuildDtree(const ALN* pALN, int nMaxDepth);


///////////////////////////////////////////////////////////////////////////////
// numerical routines 
using namespace std;
#ifndef NAN
#define NAN numeric_limits<double>::quiet_NaN( )
#endif  // NAN


// calculate probability p of an event occuring, such that
// the probablity of m or less such events occuring in n trials
// is x; currently limited to accuracy of 1.e-7
double ALNAPI PLimit(int n, int m, double dblX);
  // returns indefinite (quiet Nan) if x < 0 or x > 1 or n < 0 
  // returns 0 if m < 0
  // returns 1 if m >= n

// calculate covariance matrix C and fitted parameters A to a dataset
// with independent variables X and dependent variable Y; std dev of each 
// data point is in S
// arrays U,V,W are scratch space used by SVD
BOOL ALNAPI CalcCovariance(int nCols,   // number of input vars
                    int nRows,          // number of rows
                    double* adblX,      // RMA based input vectors (nCols * nRows)
                    double* adblY,      // RMA based result vector (nRows)
                    double* adblC,      // RMA covariance matrix (nCols * nCols)
                    double* adblA,      // RMA fitted parameter vector (nCols)
                    double* adblS,      // RMA std dev vector (nRows)
                    double* adblU,      // RMA U matrix (nCols * nRows)
                    double* adblV,      // RMA V matrix (nCols * nCols)
                    double* adblW,      // RMA W matrix (nCols)
                    double& dblChiSq);  // chi square of fit



#endif  // __ALNPRIV_H__
