// ALNfit Learning Engine for approximation of functions defined by samples.
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


// aln.h

#ifndef __ALN_H__
#define __ALN_H__

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

	/* aln version information */
#ifndef __ALNVER_H__
#include <alnver.h>
#endif

/* aln compiler / platform configuration */
#ifndef __ALNCFG_H__
#include <alncfg.h>
#endif

/* aln diagnostic services */
#ifndef __ALNDBG_H__
#include <alndbg.h>
#endif

/* dtree support */
#ifndef __DTREE_H__
#ifdef ALNDLL
#define DTREEDLL
#define DTRIMP ALNIMP
#endif  /* ALNDLL */
#include <dtree.h>
#else
#error Please do not include dtree.h before including aln.h
#endif

	typedef int BOOL;

#ifndef TRUE
#define TRUE 1
#endif

#ifndef FALSE
#define FALSE 0
#endif

	/*
	/////////////////////////////////////////////////////////////////////////////
	// MANIFEST CONSTANT DEFINITIONS
	*/

	/* node flags ---------------------------------------------------------- */
#define NF_LFN      0x00010000      /* LFN node														 */
#define NF_MINMAX     0x00020000    /* minmax node                         */
#define NF_CONSTANT 0x00040000      /* node constant                       */
#define NF_EVAL     0x00080000      /* node evaluated                      */

/* LFN flags ------------------------------------------------------------- */
#define LF_INIT     0x00000010      /* LFN initialized                     */
#define LF_SPLIT    0x00000020      /* LFN splittable                      */

/* minmax flags ---------------------------------------------------------- */
#define GF_MIN      0x00000100      /* AND minmax                          */
#define GF_MAX       0x00000200      /* OR minmax                          */

/* multi layer flags ----------------------------------------------------- */
#define MULTILAYER_FULL 0
#define MULTILAYER_RAGGED 1

/* variable monotonicity types ------------------------------------------- */
#define MONO_CONSTANT     0   /* weights == 0                              */
#define MONO_FREE         1   /* weights are both +ve and -ve              */
#define MONO_WEAKINC      2   /* weights are >= 0                          */
#define MONO_STRONGINC    3   /* weights are > 0                           */
#define MONO_WEAKDEC      4   /* weights are <= 0                          */
#define MONO_STRONGDEC    5   /* weights are < 0                           */

/* error codes ----------------------------------------------------------- */
#define ALN_NOERROR       0   /* no errors occured                         */
#define ALN_OUTOFMEM      10  /* out of memory                             */
#define ALN_USERABORT     20  /* user aborted operation                    */
#define ALN_ERRFILE       30  /* file error, check errno                   */
#define ALN_BADFILEFORMAT 31  /* invalid ALN file format                   */
#define ALN_GENERIC       100 /* unspecified error                         */


/*
/////////////////////////////////////////////////////////////////////////////
// MACROS
*/

/*
// helper macro to calc byte count of bitmap fields for n elements
// ... assumes char width is 8 bits
*/
#define MAPBYTECOUNT(n) (((n) + 7) >> 3)

/*
// helper macros test, set, and reset bits in map
// use following speedups: n / 8 == n >> 3   and   n % 8 == n & 0x07
*/
#define TESTMAP(pMap, n) (((char*)pMap)[(n)>>3] & (0x01<<((n)&0x07)))
#define SETMAP(pMap, n) (((char*)pMap)[(n)>>3] |= (0x01<<((n)&0x07)))
#define RESETMAP(pMap, n) (((char*)pMap)[(n)>>3] &= ~(0x01<<((n)&0x07)))

/* ALNREGION constraint helper macros */
#define CONSTR(pALN, nVar, nRegion) (pALN)->aRegions[nRegion].aConstr[nVar]

/* ALNNODE data access helper macros */
#define NODE_PARENT(pNode) ((pNode)->pParent)
#define NODE_FLAGS(pNode) ((pNode)->fNode)
#define NODE_MINMAXTYPE(pNode) ((pNode)->fNode & (NF_MINMAX | NF_LFN))
#define NODE_ISLFN(pNode) ((pNode)->fNode & NF_LFN)
#define NODE_ISMINMAX(pNode) ((pNode)->fNode & NF_MINMAX)
#define NODE_ISCONSTANT(pNode) ((pNode)->fNode & NF_CONSTANT)
#define NODE_REGION(pNode) ((pNode)->nParentRegion)
#define LFN_FLAGS(pNode) ((pNode)->fNode)
#define LFN_VDIM(pNode) ((pNode)->DATA.LFN.nVDim)
#define LFN_WDIM(pNode) ((pNode)->DATA.LFN.nVDim + 1)
#define LFN_VARMAP(pNode) ((pNode)->DATA.LFN.afVarMap)
#define LFN_W(pNode) ((pNode)->DATA.LFN.adblW)
#define LFN_C(pNode) ((pNode)->DATA.LFN.adblC)
#define LFN_D(pNode) ((pNode)->DATA.LFN.adblD)
#define MINMAX_FLAGS(pNode) ((pNode)->fNode)
#define MINMAX_TYPE(pNode) ((pNode)->fNode & (GF_MIN | GF_MAX))
#define MINMAX_ISMAX(pNode) ((pNode)->fNode & GF_MAX)
#define MINMAX_ISMIN(pNode) ((pNode)->fNode & GF_MIN)
#define MINMAX_LEFT(pNode) ((pNode)->DATA.MINMAX.CHILDREN.CHILDSEPARATE.pLeftChild)
#define MINMAX_RIGHT(pNode) ((pNode)->DATA.MINMAX.CHILDREN.CHILDSEPARATE.pRightChild)
#define MINMAX_NUMCHILDREN(pNode) (2)
#define MINMAX_CHILDREN(pNode) ((pNode)->DATA.MINMAX.CHILDREN.CHILDARRAY.apChildren)


#define NODE_ISEVAL(pNode) ((pNode)->fNode & NF_EVAL)
#define NODE_DISTANCE(pNode) ((pNode)->dblDistance)
#define NODE_RESPCOUNTLASTEPOCH(pNode) ((pNode)->nRespCountLastEpoch)
#define NODE_RESPCOUNT(pNode) ((pNode)->nRespCount)
#define NODE_ISUSELESS(pNode, nDim) ((pNode)->nRespCountLastEpoch + (pNode)->nRespCount < (nDim))

#define LFN_ISINIT(pNode) ((pNode)->fNode & LF_INIT)
#define LFN_CANSPLIT(pNode) ((pNode)->fNode & LF_SPLIT)
#define LFN_SPLIT(pNode) ((pNode)->DATA.LFN.pSplit)
#define LFN_SPLIT_COUNT(pNode) ((pNode)->DATA.LFN.pSplit->nCount)
#define LFN_SPLIT_SQERR(pNode) ((pNode)->DATA.LFN.pSplit->dblSqError)
#define LFN_SPLIT_RESPTOTAL(pNode) ((pNode)->DATA.LFN.pSplit->dblRespTotal)
#define LFN_SPLIT_T(pNode) ((pNode)->DATA.LFN.pSplit->dblT)

#define MINMAX_ACTIVE(pNode) ((pNode)->DATA.MINMAX.pActiveChild)
#define MINMAX_RESPACTIVE(pNode) ((pNode)->DATA.MINMAX.dblRespActive)
#define MINMAX_GOAL(pNode) ((pNode)->DATA.MINMAX.pGoalChild)
#define MINMAX_EVAL(pNode) ((pNode)->DATA.MINMAX.pEvalChild)

/*
/////////////////////////////////////////////////////////////////////////////
// STRUCTURES
*/

/* constraint structure -------------------------------------------------- */
	typedef struct tagALNCONSTRAINT
	{
		int    nVarIndex;                 /* variable index                      */
		double dblEpsilon;  							/* epsilon constraint                  */
		double dblMin, dblMax;            /* range constraints                   */
		double dblWMin, dblWMax;          /* weight constraints                  */

		/* auto calculated quantities */
		double dblSqEpsilon;  						/* store for epsilon squared value      */
	} ALNCONSTRAINT;

	/* LFN split structure - dblRespTotal is used in two different ways:
		 1. during training of linear pieces and
		 2. during tree growth by pieces splitting */
	typedef struct tagALNLFNSPLIT
	{
		int nCount;
		double dblSqError;                /* squared error                       */
		double dblRespTotal;              /* total response                      */
		double dblT;                      /* convexity criterion								 */
	} ALNLFNSPLIT;

	/* node structure -------------------------------------------------------- */
	typedef struct tagALNNODE
	{
		struct tagALNNODE* pParent;       /* pointer to parent node              */
		int nParentRegion;                /* index of parent region in ALN,      */
																			/*   currently unused, must be 0       */
		int fNode;                        /* node flags (NF_*, LF_*, GF_*)       */
		int nRespCount;                   /* responsibility count current epoch  */
		int nRespCountLastEpoch;          /* responsibility count previous epoch */
		double dblDistance;               /* distance to current input point     */

		union tagDATA
		{
			struct tagLFN
			{
				char* afVarMap;               /* var index bitmap,                   */
																			/*   currently unused, must be NULL    */
				int nVDim;                    /* vector dimensions, except adblW     */
				double* adblW;                /* weight vector, nVDim + 1 elements   */
				double* adblC;                /* centroid vector                     */
				double* adblD;                /* ave sq dist from centroid vector    */
				ALNLFNSPLIT* pSplit;          /* split structure                     */
			} LFN;
			struct tagMINMAX
			{
				double dblRespActive;           /* response on active child          */
				struct tagALNNODE* pActiveChild;/* active child on current input     */
				struct tagALNNODE* pGoalChild;  /* goal child on current input       */
				struct tagALNNODE* pEvalChild;  /* first eval child on current input */
				union tagCHILDREN
				{
					struct tagCHILDSEPARATE
					{
						struct tagALNNODE* pLeftChild;  /* ptr to left child             */
						struct tagALNNODE* pRightChild; /* ptr to right child            */
					} CHILDSEPARATE;
					struct tagCHILDARRAY
					{
						struct tagALNNODE* apChildren[2]; /* array of children           */
					} CHILDARRAY;
				} CHILDREN;
			} MINMAX;
		} DATA;
	} ALNNODE;

	/* ALNREGION structure --------------------------------------------------- */
	typedef struct tagALNREGION
	{
		int nParentRegion;                /* index of parent region in ALN       */
																			/*   currently unused, must be -1      */
		double dblLearnFactor;	          /* learning rate factor                */
		char* afVarMap;                   /* var index bitmap,                   */
																			/*   currently unused, must be NULL    */
		int nConstr;                      /* number of constraints in array      */
		ALNCONSTRAINT* aConstr;           /* array of constraints                */

		double dblSmoothEpsilon;          /* smoothing epsilon                   */

		/* auto calculated quantities */
		double dbl4SE;                    /* 4 * smoothing epsilon               */
		double dblOV16SE;                 /* 1 / 16 * smoothing epsilon          */
	} ALNREGION;

	/* confidence interval structure ----------------------------------------- */
	typedef struct tagALNCONFIDENCE
	{
		double dblP;          /* one-sided tail probability, the actual          */
													/* symmetric confidence interval is 1 - 2p         */
		double dblLowerBound; /* upper bound on error                            */
		double dblUpperBound; /* lower bound on error                            */
		int    nSamples;      /* number of samples used when calculating bounds  */
	} ALNCONFIDENCE;

	/* LFN analysis structures ----------------------------------------------- */
	typedef struct tagLFNSTATS
	{
		double dblRSS;        /* regression sum of squares                       */
		double dblESS;        /* Error (residual) sum of squares                 */
		double dblDF;         /* degrees of freedom                              */
		double dblR2;         /* doefficient of determination                    */
		double dblSEE;        /* standard error of estimation                    */
		double dblF;          /* F-statistic                                     */
		double dblFp;         /* probability of F                                */
	} LFNSTATS;

	typedef struct tagLFNWEIGHTSTATS
	{
		double dblW;          /* value of weight                                 */
		double dblSEw;        /* standard error of the weight                    */
		double dblT;          /* T-statistic                                     */
		double dblTp;         /* probability of T                                */
	} LFNWEIGHTSTATS;


	/* ALN struct ------------------------------------------------------------ */
	typedef struct tagALN
	{
		int nVersion;                     /* version of ALN                      */
		int nDim;                         /* number of ALN inputs + 1 for output */
		int nOutput;                      /* index of output var, usually nDim-1 */
		int nRegions;                     /* number of regions only partly       */
																			/* implemented;  currently must be 1   */
		ALNREGION* aRegions;              /* array of regions, nRegions elements */
		ALNNODE* pTree;                   /* pointer to root node of tree        */
	} ALN;

	/*
	// structure used in training and evaluation for indicating
	// column index and time shift for each variable
	*/
	typedef struct tagVARINFO
	{
		int nCol;   /* associated data array column, -1 if same as var index     */
		int nDelta; /* "time" shift                                              */
	} VARINFO;

	/* structure used for passing training and eval data                       */
	typedef struct tagVECTORINFO
	{
		int nPoint;               /* training or eval sequence number            */
		int bNeedData;            /* TRUE if callback must supply data           */
		const VARINFO* aVarInfo;  /* VARINFO array, may be NULL                  */
		double* adblX;	          /* input vector, can be modified               */
	} VECTORINFO;

	/* structures used for passing info to training notification procedure     */
	typedef struct tagEPOCHINFO
	{
		int nEpoch;					      /* epoch number                                */
		int nLFNs;					      /* number of LFNs in ALN                       */
		int nActiveLFNs;			    /* number of active LFNs in ALN                */
		double dblEstRMSErr;	    /* estimated RMS error                         */
	} EPOCHINFO;

	typedef struct tagTRAININFO
	{
		int nEpochs;					    /* total number of epochs required,            */
		int nLFNs;						    /* number of LFNs in ALN                       */
		int nActiveLFNs;			    /* number of resp. LFNs in ALN, 0 at start     */
		double dblRMSErr;			    /* RMS error, 0 at start                       */
	} TRAININFO;

	typedef struct tagADAPTINFO
	{
		int nAdapt;						    /* adaptation sequence number                  */
		const double* adblX;	    /* training vector                             */
		double dblErr;				    /* signed distance from point to surface       */
	} ADAPTINFO;

	typedef struct tagLFNADAPTINFO
	{
		const double* adblX;      /* training vector                             */
		ALNNODE* pLFN;            /* pointer to adapting LFN                     */
		double dblResponse;       /* the weight on the adapting LFN              */
		double dblError;          /* global error we are adapting to             */
	} LFNADAPTINFO;

	typedef struct tagALNDATAINFO
	{
		int nPoints;              /* number of data points                       */
		const VARINFO* aVarInfo;  /* variable info... may be NULL                */
		const double* adblData;   /* data array... may be NULL                   */
		int nCols;                /* number of columns in non-NULL data array    */
	} ALNDATAINFO;

	/*
	/////////////////////////////////////////////////////////////////////////////
	// ALN notification callback prototype
	*/

	/* return 0 if caller should halt, non-zero if caller should continue      */
	typedef int (ALNAPI* ALNNOTIFYPROC)(const ALN* pALN, int nCode, void* pParam,
		void* pvData);

	typedef struct tagALNCALLBACKINFO
	{
		int nNotifyMask;              /* notifications to send                   */
		ALNNOTIFYPROC pfnNotifyProc;  /* notification callback                   */
		void* pvData;                 /* user data                               */
	} ALNCALLBACKINFO;


	/* ALN notification callback codes and corresponding pParam meanings       */
#define AN_TRAINSTART    0x0001
	/* TRAININFO* pTrainInfo = (TRAININFO*)pParam                            */

#define AN_TRAINEND			 0x0002
	/* TRAININFO* pTrainInfo = (TRAININFO*)pParam                            */

#define AN_EPOCHSTART		 0x0004
	/* EPOCHINFO* pEpochInfo = (EPOCHINFO*)pParam                            */

#define AN_EPOCHEND			 0x0008
	/* EPOCHINFO* pEpochInfo = (EPOCHINFO*)pParam                            */

#define AN_ADAPTSTART		 0x0010
	/* ADAPTINFO* pAdaptInfo = (ADAPTINFO*)pParam                            */

#define AN_ADAPTEND			 0x0020
	/* ADAPTINFO* pAdaptInfo = (ADAPTINFO*)pParam                            */

#define AN_LFNADAPTSTART 0x0040
	/* LFNADAPTINFO* pLFNAdaptInfo = (LFNADAPTINFO*)pParam                   */

#define AN_LFNADAPTEND	 0x0080
	/* LFNADAPTINFO* pLFNAdaptInfo = (LFNADAPTINFO*)pParam                   */

#define AN_VECTORINFO    0x0100
	/* VECTORINFO* pVectorInfo = (VECTORINFO*)pParam                         */

#define AN_NONE     0
#define AN_TRAIN    (AN_TRAINSTART|AN_TRAINEND)
#define AN_EPOCH    (AN_EPOCHSTART|AN_EPOCHEND)
#define AN_ADAPT    (AN_ADAPTSTART|AN_ADAPTEND)
#define AN_LFNADAPT (AN_LFNADAPTSTART|AN_LFNADAPTEND)
#define AN_ALL      (AN_TRAIN|AN_EPOCH|AN_ADAPT|AN_LFNADAPT|AN_VECTORINFO)

/*
/////////////////////////////////////////////////////////////////////////////
// FUNCTIONS
*/

/*
// test ALN structure validity
*/
	ALNIMP int ALNAPI ALNTestValid(const ALN* pALN);

	/*
	// seeding for ALN internal pseudo-random number generator
	*/
	ALNIMP void ALNAPI ALNSRand(unsigned int nSeed);

	/*
	// return next value from ALN internal pseudo-random number generator
	*/
	ALNIMP unsigned long ALNAPI ALNRand(void);

	/*
	// max value returned by ALNRand()
	*/
#define ALNRAND_MAX 0xffffffff

	/*
	// quick way to get random float value in [0, 1]
	*/
	ALNIMP float ALNAPI ALNRandFloat(void);

	/*
	// TrainALN
	*/

	ALNIMP int ALNAPI ALNTrain(ALN* pALN,
		const ALNDATAINFO* pDataInfo,
		const ALNCALLBACKINFO* pCallbackInfo,
		int nMaxEpochs,
		double dblMinRMSErr,
		double dblLearnRate,
		BOOL bJitter);

	/*
	// ALNCalcRMSError
	*/

	ALNIMP int ALNAPI ALNCalcRMSError(const ALN* pALN,
		const ALNDATAINFO* pDataInfo,
		const ALNCALLBACKINFO* pCallbackInfo,
		double* pdblRMSErr);

	/*
	// ALN variable monotonicity type
	*/
	ALNIMP int ALNAPI ALNVarMono(const ALN* pALN, int nVar, int* pnMono);

	/*
	// ALN inversion
	*/
	ALNIMP int ALNAPI ALNInvert(ALN* pALN, int nVar);

	/*
	// evaluation of ALN on data
	*/

	ALNIMP int ALNAPI ALNEval(const ALN* pALN,
		const ALNDATAINFO* pDataInfo,
		const ALNCALLBACKINFO* pCallbackInfo,
		double* adblResult,
		int* pnStart, int* pnEnd);

	/*
	// quick evaluation of ALN on single vector
	*/
	ALNIMP double ALNAPI ALNQuickEval(const ALN* pALN, const double* adblX,
		ALNNODE** ppActiveLFN);


	/*
	/////////////////////////////////////////////////////////////////////////////
	// confidence intervals
	*/

	/*
	// calc the confidence intervals for a data set
	// the ALNCONFIDENCE structure must have its dblP member intialized
	// to contain the desired area in each tail outside the confidence
	// interval - dblP must be greater than 0 and less than 0.5
	// returns ALN_* error code, (ALN_NOERROR on success)
	*/
	ALNIMP int ALNAPI ALNCalcConfidence(const ALN* pALN,
		const ALNDATAINFO* pDataInfo,
		const ALNCALLBACKINFO* pCallbackInfo,
		ALNCONFIDENCE* pConfidence);

	/*
	// confidence in the confidence intervals:
	// calculates a value PLimit such that there is only a chance
	// dblSignificance (usually 0.05 or 0.01) that the interval bounds
	// in pConfidence would have been chosen to bound tails of area
	// dblPLimit
	*/
	ALNIMP int ALNAPI ALNConfidencePLimit(const ALNCONFIDENCE* pConfidence,
		double dblSignificance,
		double* pdblPLimit);


	/*
	// confidence in the confidence intervals:
	// calculates probability that the confidence interval in pConfidence
	// covers a specified minimum fraction (tolerance interval) of the
	// symmetric distribution having area dblInterval (usually
	// 1 - 4 * pConfidence->dblP)
	*/

	ALNIMP int ALNAPI ALNConfidenceTLimit(const ALNCONFIDENCE* pConfidence,
		double dblInterval,
		double* pdblTLimit);


	/*
	/////////////////////////////////////////////////////////////////////////////
	// LFN analysis
	*/

	ALNIMP int ALNAPI ALNLFNAnalysis(const ALN* pALN,
		const ALNDATAINFO* pDataInfo,
		const ALNCALLBACKINFO* pCallbackInfo,
		void*& ppvAnalysis,
		int& pnLFNStats);

	ALNIMP int ALNAPI ALNLFNStats(void* pvAnalysis,
		int nLFNStat,
		LFNSTATS* pLFNStats,
		int* pnWeightStats,
		ALNNODE** ppLFN);

	// if nWeight == the output index of the ALN, 
	// then the weight is the bias weight of the LFN
	ALNIMP int ALNAPI ALNLFNWeightStats(void* pvAnalysis,
		int nLFNStat,
		int nWeightStat,
		LFNWEIGHTSTATS* pLFNWeightStats);

	ALNIMP int ALNAPI ALNLFNFreeAnalysis(void* pvAnalysis);


	/*
	/////////////////////////////////////////////////////////////////////////////
	// I/O
	*/

	/*
	// saving ALN to disk file
	*/
	ALNIMP int ALNAPI ALNWrite(const ALN* pALN, const char* pszFileName);

	/*
	// loading ALN from disk file
	*/
	ALNIMP int ALNAPI ALNRead(const char* pszFileName, ALN** ppALN);

	/*
	// conversion to dtree
	*/
	ALNIMP int ALNAPI ALNConvertDtree(const ALN* pALN, int nMaxDepth,
		DTREE** ppDtree);

	/*
	/////////////////////////////////////////////////////////////////////////////
	// memory management
	*/

	/*
	// creates a new ALN
	*/
	ALNIMP ALN* ALNAPI ALNCreateALN(int nDim, int nOutput);

	/*
	// destroys an ALN
	// returns 0 on failure, non-zero on success
	*/
	ALNIMP int ALNAPI ALNDestroyALN(ALN* pALN);

#ifdef ENABLE_REGIONS

	/*
	// adding a new region to the ALN
	// returns index of new region in ALN, -1 on failure
	*/
	ALNIMP int ALNAPI ALNAddRegion(ALN* pALN, int nParentRegion,
		double dblLearnFactor,
		int nConstr, int* anConstr);

#endif /* ENABLE_REGIONS */

	/*
	// adding LFNs to a tree
	*/
	ALNIMP int ALNAPI ALNAddLFNs(ALN* pALN, ALNNODE* pParent,
		int nParentMinMaxType, int nLFNs,
		ALNNODE** apLFNs);

	/*
	// adding multiple layers to a tree
	*/
	ALNIMP int ALNAPI ALNAddMultiLayer(ALN* pALN, ALNNODE* pParent,
		int nParentMinMaxType, int nLayers,
		int nMinMaxFanin, int nLFNFanin,
		int nFlags);

	/*
	// make a tree from a string
	*/

	ALNIMP int ALNAPI ALNAddTreeString(ALN* pALN, ALNNODE* pParent,
		const char* pszTreeString,
		int* pnParsed);

	/*
	// make a growable subtree
	*/

	ALNIMP int ALNAPI ALNSetGrowable(ALN* pALN, ALNNODE* pParent);


	/*
	///////////////////////////////////////////////////////////////////////////////
	// Abort handling
	*/

	/*
	// set ALN abort procedure: if ALN needs to abort, it will call this
	//   function; if not set, ALN will simply call abort()
	// you can use this to chain an abort proc for another library to the ALN
	//   library, and simply call ALNAbort
	*/
	typedef void(__stdcall* ALNABORTPROC)(void);
	ALNIMP void ALNAPI ALNSetAbortProc(ALNABORTPROC pfnAbortProc);

	/*
	// ALN library abort... cleans up internal ALN library, then calls
	// abort procedure defined by call to ALNSetAbortProc()
	*/
	ALNIMP void __stdcall ALNAbort(void);


#ifdef __cplusplus
}
#endif /* __cplusplus */

/*
///////////////////////////////////////////////////////////////////////////////
*/

#endif  /* __ALN_H__ */