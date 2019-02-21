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

// builddtree.cpp
// dtree generation routines

#ifdef ALNDLL
#define ALNIMP __declspec(dllexport)
#endif

#include <aln.h>
#include "alnpriv.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

///////////////////////////////////////////////////////////////////////////////
// Building DTREE from ALN

// helpers for building dtree from ALN
VARDEF* ALNAPI BuildVarDefs(const ALN* pALN);
MINMAXNODE* ALNAPI CopyALN(const ALNNODE* pNode, const ALN* pALN, 
                           LINEARFORM* aLF, int nOutput, int nMono, 
                           int* pnLFN);

// pSrc points to a dtree with only one node and one block
// *ppDest will contain the new, optimized, dtree 
// return DTE_* error code, DTE_NOERR on success
int ALNAPI SplitDtree(DTREE** ppDest, DTREE* pSrc, int nMaxDepth);

// building the dtree from an ALN
DTREE* ALNAPI BuildDtree(const ALN* pALN, int nMaxDepth)
{
  ASSERT(pALN);
  ASSERT(pALN->pTree);
  ASSERT(nMaxDepth >= DTREE_MINDEPTH && nMaxDepth <= DTREE_MAXDEPTH);

  // reset DTREE error code here
  dtree_errno = DTR_NOERROR;

  // this version of library, output is always ALN output
  int nOutput = pALN->nOutput;

  // check monotonicity... we can cope with anything but free!
  int nMono = CheckMonotonicity(pALN->pTree, pALN, nOutput);
  if (nMono == MONO_FREE)
    return NULL;

  // this version of library, since output is always ALN output, we expect 
  // monotonicity to be strongly decreasing
  ASSERT(nMono == MONO_STRONGDEC);
  
  DTREE* pSrc = CreateDtree();
  if (pSrc == NULL)
    return NULL;
  
  // dimension, output and vars
  pSrc->nDim = pALN->nDim;
  pSrc->nOutputIndex = nOutput;
  pSrc->aVarDefs = BuildVarDefs(pALN);
  if (pSrc->aVarDefs == NULL)
  {
    DestroyDtree(pSrc);
    return NULL;
  }
           
  // linear forms           
  int nLFNs = 0, nAdaptedLFNs = 0;
  CountLFNs(pALN->pTree, nLFNs, nAdaptedLFNs);
  pSrc->aLinearForms = CreateLinearFormArray(nLFNs, pSrc->nDim);
  if (pSrc->aLinearForms == NULL)                   
  {
    DestroyDtree(pSrc);
    return NULL;
  }
  pSrc->nLinearForms = nLFNs;
  
  // blocks
  int nOneBlock = 1;
  pSrc->aBlocks = CreateBlockArray(nOneBlock);
  if (pSrc->aBlocks == NULL)
  {
    DestroyDtree(pSrc);
    return NULL;
  }
  pSrc->nBlocks = 1;
  
  // nodes
  int nOneNode = 1;
  pSrc->aNodes = CreateDtreeNodeArray(nOneNode);
  if (pSrc->aNodes == NULL)
  {        
    DestroyDtree(pSrc);
    return NULL;
  }
  pSrc->nNodes = 1;
  
  // single minmax tree copied from the ALN
  int nLFNsCopied = 0;
  pSrc->aBlocks->pMinMaxTree = CopyALN(pALN->pTree, pALN, pSrc->aLinearForms,
                                       nOutput, nMono, &nLFNsCopied);
  if (pSrc->aBlocks->pMinMaxTree == NULL)
  {             
    DestroyDtree(pSrc);
    return NULL;
  }
  ASSERT(nLFNsCopied == nLFNs);
  pSrc->aBlocks->nDtreeIndex = 0;
                  
  // 1st node                  
  pSrc->aNodes->nLeaf = 1;    // this is a leaf!
  DNODE_BLOCKINDEX(pSrc->aNodes) = 0;
  
  // now optimize when nMaxDepth is >1... set dtree_errno if any problems
  DTREE* pDest = NULL;
	if(nMaxDepth > 1)
	{
		// this does a lot of work
		dtree_errno = SplitDtree(&pDest, pSrc, nMaxDepth);
		if (dtree_errno != DTR_NOERROR) 
		{
			DestroyDtree(pDest);
			pDest = NULL;
		}
	  DestroyDtree(pSrc);
	}
	else
	{
		// this avoids work, and also the question of
		// whether the DTREE constructed directly from
		// the ALN is optimal.
		pDest = pSrc;
	}
  return pDest;  
}

// vardef array builder
static 
VARDEF* ALNAPI BuildVarDefs(const ALN* pALN)
{
  ASSERT(pALN);
  ASSERT(pALN->aRegions);
  ASSERT(pALN->aRegions->nConstr == pALN->nDim);

  // allocate space for nDim vars
  int nDim = pALN->nDim;
  VARDEF* aVD = CreateVarDefArray(nDim);
  if (aVD == NULL)
    return NULL;
    
  // copy var bounds, set name
  for (int i = 0; i < nDim; i++)
  {
    // get topmost region var constraint
    ALNCONSTRAINT* pConstr = &(pALN->aRegions[0].aConstr[i]);
    ASSERT(pConstr->nVarIndex == i);

    VARDEF_MIN(aVD + i) = pConstr->dblMin;
    VARDEF_MAX(aVD + i) = pConstr->dblMax;
    
    char szNameBuf[128];
    sprintf(szNameBuf, "x%d", i);
    SetVarDefName(aVD + i, szNameBuf);
  }

  return aVD;
}

// builds min max tree and fills linear form array
static 
MINMAXNODE* ALNAPI CopyALN(const ALNNODE* pNode, const ALN* pALN, 
                           LINEARFORM* aLF, int nOutput, int nMono, 
                           int* pnLFN)
{
  ASSERT(pALN);
  ASSERT(pNode);
  ASSERT(nOutput >= 0 && nOutput < pALN->nDim);
  ASSERT(nMono == MONO_CONSTANT || nMono == MONO_STRONGINC || 
         nMono == MONO_WEAKINC || nMono == MONO_STRONGDEC ||
         nMono == MONO_WEAKDEC);
  ASSERT(pnLFN && *pnLFN >= 0);

  // allocate dtree node  
  MINMAXNODE* pMMN = CreateMinMaxNode();
  if (pMMN == NULL)
    return NULL;

  // copy ALN node  
  if (NODE_ISMINMAX(pNode))
  {   
    int nMinMaxType = MINMAX_TYPE(pNode);
    
    // this node is a MIN iff AND and decreasing mono, or OR and increasing mono
    // else node is a MAX iff AND and increasing mono, or OR and decreasing mono
    // NOTE: MONO_CONSTANT will be treated as decreasing!

    BOOL bInc = nMono == MONO_STRONGINC || nMono == MONO_WEAKINC;
    if (((nMinMaxType & GF_MIN) && !bInc) || ((nMinMaxType & GF_MAX) && bInc))
      pMMN->nType = DTREE_MIN;
    else
      pMMN->nType = DTREE_MAX;
             
    // set children              
    int nChildren = MINMAX_NUMCHILDREN(pNode);
    ASSERT(nChildren > 0);
    ALNNODE* const* apChildren = MINMAX_CHILDREN(pNode);
    ASSERT(apChildren);

    MINMAXNODE* pList = MMN_CHILDLIST(pMMN);              
    for (int i = 0; i < nChildren; i++)
    { 
      MINMAXNODE* pChild = CopyALN(apChildren[i], pALN, aLF, nOutput, nMono, 
                                   pnLFN);
      if (pList == NULL)
      {
        MMN_CHILDLIST(pMMN) = pChild;
        pList = pChild;
      }
      else
      {
        pList->pNext = pChild;
        pList = pChild;
      }
    }

    return pMMN;
  }

  ASSERT(NODE_ISLFN(pNode));
  
  // set op
  pMMN->nType = DTREE_LINEAR;
    
  // set linear form index
  MMN_LFINDEX(pMMN) = *pnLFN;   
  
  // get associated linear form
  LINEARFORM* pLF = &(aLF[*pnLFN]);  
  ASSERT(pLF);

  // inc LF count
  (*pnLFN)++;
  ASSERT(*pnLFN - MMN_LFINDEX(pMMN) == 1);
  
  // get output var weight
  double dblWOutput = LFN_W(pNode)[nOutput + 1];  // account for bias weight
  if (dblWOutput == 0)
  {                   
    // adjust weight to be non-zero by dividing current output var epsilon
    // with range of desired output var at the topmost region... this is the
    // smallest meaningful value that a weight can take (the largest is the
    // range of current output divided by epsilon of desired output)

    // get current output var epsilon
    ASSERT(nOutput != pALN->nOutput); 
      // otherwise weight would be non-zero (-1) on def output var!
    
    ALNCONSTRAINT* pConstrOut = GetVarConstraint(NODE_REGION(pNode), pALN, 
                                                 pALN->nOutput);
    ASSERT(pConstrOut);
    double dblEpsilon = pConstrOut->dblEpsilon;
      
    // get min and max on new output var in topmost region
    ASSERT(pALN->aRegions[0].aConstr[nOutput].nVarIndex == nOutput);
    double dblMin = pALN->aRegions[0].aConstr[nOutput].dblMin;
    double dblMax = pALN->aRegions[0].aConstr[nOutput].dblMax;
    
    dblWOutput = dblEpsilon / (dblMax - dblMin);
    ASSERT(dblWOutput > 0);
    
    // determine sign of weight
    if (nMono != MONO_WEAKINC && nMono != MONO_STRONGINC)
      dblWOutput *= -1; // supposed to be decreasing

    // otherwise leave it positive
  } 
  ASSERT(dblWOutput != 0);
  
  int nDim = pALN->nDim;
  double* adblW = LFN_W(pNode) + 1; // skip bias
  double* adblC = LFN_C(pNode);

  for (int i = 0; i < nDim; i++)
  {
    if (i == nOutput)
      pLF->adblW[i] = dblWOutput;
    else
      pLF->adblW[i] = adblW[i];
    
    if (dblWOutput > 0) 
      pLF->adblW[i] *= -1; // invert so output weight is effectively negative
    
    pLF->adblC[i] = adblC[i];
  }
  pLF->dblBias = LFN_W(pNode)[0];
  if (dblWOutput > 0) 
    pLF->dblBias *= -1;    // invert so output weight is effectively negative
    
  return pMMN;
}
///////////////////////////////////////////////////////////////////////////////
// optimizing DTREE

#define CHILDABOVEBESTBOUND  1
#define BOUNDSOVERLAP        0
#define CHILDBELOWBESTBOUND -1

// splitting helper prototypes 
int Split(int nNode, int nMaxNodes, DTREENODE* aNodes, int* pnNodes, 
          BLOCK* aBlocks, int* pnBlocks, LINEARFORM* aLF, int nLF, 
          double* adblMin, double* adblMax, double* adblX,
          double* adblRespMin, double* adblRespMax, char* aRespLF,
          double& dblBiasBound, double* adblWBound, double& dblHalfWidth,
          int nDim, int nOutput, int nDepth, int nMaxDepth,  int* aNewIndex,
          int* pnCount, BOOL &bSmaller);

void FindBestSplit(int* pnVarIndex, double* pdblT, double* adblMin, 
                   double* adblMax, double* adblRespMin, double* adblRespMax,
                   char* aRespLF, int nLF, int nDim, int nOutput, int nLines);

void FindBestT(int nVarIndex, double* pdblT, double* adblMin, double* adblMax, 
               double* adblRespMin, double* adblRespMax, char* aRespLF, 
               int nLF, int nDim, int nOutput, int nLines, int* pnLeft, 
               int* pnRight);

void SetResp(MINMAXNODE* pMMN, LINEARFORM* aLF, int nLF, double* adblMin, 
             double* adblMax, double* adblX, double* adblRespMin, 
             double* adblRespMax, char* aRespLF, int nDim, int nOutput, 
             int nLines);

void CountLeftRightResp(double dblT, int* pnLeft, int* pnRight, char* aRespLF, 
                        double* adblRespMin, double* adblRespMax, 
                        int nVarIndex, int nDim, int nLF);

void CountMinMaxTreeLines(MINMAXNODE* pMMN, int& nCount);
   
void MapNewLinearForms(MINMAXNODE* pMMN, int* pnCount,int* aNewIndex);

void Reindex(MINMAXNODE* pMMN, int* aNewIndex);

void Amalgamate(MINMAXNODE* pMMN, LINEARFORM* aLF, 
                   int nDim, int nOutput, BOOL &bSmaller);

void MinMaxNodeBoundCylinder(MINMAXNODE* pMMN, LINEARFORM* aLF, 
                             double* adblMin, double* adblMax,
                             double& dblBiasBound, double* adblWBound, double& dblHalfWidth,
                             int nDim, int nOutput, BOOL &bSmaller);

void LinearFormBoundCylinder(LINEARFORM* pLF, double* adblMin, double* adblMax, 
                             double& dblBiasBound, double* adblWBound,
                             int nDim, int nOutput);

int aboveORbelow(double* adblMin, double* adblMax,
                 double& dblMin1, double& dblMax1, double dblMin2, double dblMax2,
                 double& dblBias1, double* adblW1, double& dlbHalfWidth1,
                 double dblBias2, double* adblW2, double dlbHalfWidth2,
                 int nDim, int nOutput, int nMinMaxType);

static
int ALNAPI SplitDtree(DTREE** ppDest, DTREE* pSrc, int nMaxDepth)
{
  // This is the main routine for creating an optimized, multilevel DTREE.
  // pSrc points to a source DTREE with only one block and one DTREE node.
	// A DTREE node indicates a split of the input space where one DTREE on a block 
	// has been replaced by two or more DTREEs on several blocks partitioning the original one.
	//(i.e. a DTREE node is not a MINMAX node, of which the DTREE on a block can have many!).
  // The returned ppDest will point to a pointer to the new, optimized, DTREE. 
  int i;
  int nDim;
  int nOutput;
  int nBlocks;
  int nMaxBlocks;
  int nNodes;
  int nMaxNodes;
  int nErr;
  BOOL bSmaller;
  BLOCK* aBlocks = NULL;
  DTREENODE* aNodes = NULL;
  double* adblMin = NULL;
  double* adblMax = NULL;
  char* aRespLF = NULL;
  double* adblRespMin = NULL;
  double* adblRespMax = NULL;
  double* adblX = NULL;
  double dblBiasBound = 0;
  double dblHalfWidth = 0;
  double* adblWBound = NULL;
  // is the nMaxDepth of the DTREE to be created within bounds?
  if (nMaxDepth < DTREE_MINDEPTH || nMaxDepth > DTREE_MAXDEPTH)
    return DTR_GENERIC;
  // does the source DTREE check out?         
  if (pSrc == NULL)
    return DTR_GENERIC;
    
  if (pSrc->nNodes != 1 || pSrc->nBlocks != 1)
    return DTR_GENERIC;// we optimize only single-node, single-block DTREEs
    
  if (pSrc->aBlocks == NULL || pSrc->aBlocks[0].pMinMaxTree == NULL)
    return DTR_GENERIC;

  // create the destination DTREE        
  if ((*ppDest = CreateDtree()) == NULL)
    return DTR_MALLOCFAILED;
  
  // dimension, output and vars
  nDim = pSrc->nDim;
  nOutput = pSrc->nOutputIndex;
  (*ppDest)->nDim = nDim;
  (*ppDest)->nOutputIndex = nOutput;
  
  // variable defs
  if (((*ppDest)->aVarDefs = CreateVarDefArray(nDim)) == NULL)
  {
    DestroyDtree(*ppDest);
    return DTR_MALLOCFAILED;
  }
  
  // copy var bounds, names
  for (i = 0; i < nDim; i++)
  {
    (*ppDest)->aVarDefs[i].bound.dblMin = pSrc->aVarDefs[i].bound.dblMin;
    (*ppDest)->aVarDefs[i].bound.dblMax = pSrc->aVarDefs[i].bound.dblMax;
    if (SetVarDefName((*ppDest)->aVarDefs + i, pSrc->aVarDefs[i].pszName) == NULL)
    {
      DestroyDtree(*ppDest);
      return DTR_MALLOCFAILED;
    }
  }
           
  // copy linear forms           
  if (((*ppDest)->aLinearForms = CreateLinearFormArray(pSrc->nLinearForms, nDim)) == NULL)
  {
    DestroyDtree(*ppDest);
    return DTR_MALLOCFAILED;
  } 
  (*ppDest)->nLinearForms = pSrc->nLinearForms;
  // as the DTREE splits, we track the linear forms that end up in the final
  // blocks.  Some linear forms that do cutoffs simply vanish in the smaller blocks
  // and some linear forms are optimized away
  int* aNewIndex;
  if((aNewIndex = (int*) malloc(pSrc->nLinearForms * sizeof(int))) == NULL)
  {
    DestroyDtree(*ppDest);
    return DTR_MALLOCFAILED;
  }
  // as final blocks of the DTREE are created, the linear forms used are renumbered wrt
  // the order of creation of these final blocks.  New indexes are initialized to -1.
  for (i = 0; i < pSrc->nLinearForms; i++)
  {
    aNewIndex[i] = -1; // the NewIndex of the linear form originally at i is invalid if it is -1.
  }
  int nNewLFCount = 0; // counts the linear forms that are finally used (see Split(0,..)
  for (i = 0; i < pSrc->nLinearForms; i++)
  {
    int j;
    // check for zero output weight
    if (pSrc->aLinearForms[i].adblW[nOutput] == 0)
    {
      DestroyDtree(*ppDest);
      return DTR_ZEROOUTPUTWEIGHT;
    }
    for (j = 0; j < nDim; j++)
    {
      if(j != nOutput) // adjust the normalization to output weight -1
      {
        (*ppDest)->aLinearForms[i].adblW[j] = 
          -pSrc->aLinearForms[i].adblW[j]/pSrc->aLinearForms[i].adblW[nOutput];
      }
      else
      {
        (*ppDest)->aLinearForms[i].adblW[j] = -1;
      }
      (*ppDest)->aLinearForms[i].adblC[j] = pSrc->aLinearForms[i].adblC[j];
    }
    (*ppDest)->aLinearForms[i].dblBias = pSrc->aLinearForms[i].dblBias;
  }
  
  // blocks -- the greatest number of blocks for a DTREE of depth nMaxDepth is
  nMaxBlocks = (int)pow((double)2.0, (double)(nMaxDepth - 1)); // DTREE depth 1 gives one block.
	// However, since the splitting of the input space goes to depths that
	// vary with the complexity of the function on a block being split,
	// this balanced tree model is poor.  
  aBlocks = CreateBlockArray(nMaxBlocks);
  if (aBlocks == NULL)
  {
    DestroyDtree(*ppDest);
    return DTR_MALLOCFAILED;
  }
  
  // DTREE nodes -- the number of DTREE nodes is limited by the number of blocks.
	// Here is the reasoning: we allow nMaxBlocks from splitting the original
	// single DTREE block.  Each time we split a block (a leaf node), we remove it and add
	// one DTREE non-leaf node and two blocks (leaf nodes). In sum we add two new nodes.
	// Starting at one node and one block and splitting k times we get
	// k+1 blocks (leaf nodes)and k non-leaf nodes for a total of 2*k +1 nodes.
	// Since we don't know how many times optimization will split blocks, we have
	// to allocate the maximum number of blocks at the start of splitting.
	nMaxNodes = 2* nMaxBlocks - 1; 
  aNodes = CreateDtreeNodeArray(nMaxNodes);
  if (aNodes == NULL)
  {        
    DestroyBlockArray(aBlocks, nMaxBlocks);
    DestroyDtree(*ppDest);
    return DTR_MALLOCFAILED;
  }
  
  // 1st minmax tree
  nBlocks = 1;
  aBlocks[0].pMinMaxTree = CopyMinMaxNode(pSrc->aBlocks[0].pMinMaxTree);
  if (aBlocks[0].pMinMaxTree == NULL)
  {             
    DestroyDtreeNodeArray(aNodes);
    DestroyBlockArray(aBlocks, nMaxBlocks);
    DestroyDtree(*ppDest);
    return DTR_MALLOCFAILED;
  }
  aBlocks[0].nDtreeIndex = 0;
                  
  // 1st node                  
  nNodes = 1;
  aNodes[0].nLeaf = 1;
  DNODE_BLOCKINDEX(aNodes + 0) = 0;
  
  // variable bounds
  adblMin = (double*)malloc(nDim * sizeof(double));
  adblMax = (double*)malloc(nDim * sizeof(double));
  // responsibilty flag map
  aRespLF = (char*)malloc(MAPBYTECOUNT(pSrc->nLinearForms));
  // responsibility bounds    
  adblRespMin = (double*)malloc(pSrc->nLinearForms * nDim * sizeof(double));
  adblRespMax = (double*)malloc(pSrc->nLinearForms * nDim * sizeof(double));
  // input vector
  adblX = (double*)malloc(nDim * sizeof(double));

  // variables and array for storage of bounds in the form of linear functions
  adblWBound = (double*) malloc(nDim * sizeof(double));
  // The centroid used for all linear function bounds is the average of the current max and min
  // bounds of the block in each axis.  This is done to reduce the effect of numerical errors.
  // There are two types of bounds used for the values in a box.  The max and the min on the box
  // are the simplest to evaluate.  Two parallel hyperplanes bounding the function values are
  // also used.  They are defined by giving the bias and weights of a linear function,
  // and a "half width" which is an offset to be added to the above linear function to
  // get an upper bound, or subtracted from the linear function to get a lower bound.
  
  if (adblMin == NULL || adblMax == NULL || aRespLF == NULL ||
      adblRespMin == NULL || adblRespMax == NULL || adblX == NULL|| adblWBound == NULL)
  {
    if (adblMin != NULL) free(adblMin);
    if (adblMax != NULL) free(adblMax);
    if (aRespLF != NULL) free(aRespLF);
    if (adblRespMin != NULL) free(adblRespMin);
    if (adblRespMax != NULL) free(adblRespMax);
    if (adblX != NULL) free(adblX);
    if (adblWBound != NULL) free(adblWBound);
    DestroyDtreeNodeArray(aNodes);
    DestroyBlockArray(aBlocks, nMaxBlocks);
    DestroyDtree(*ppDest);
    return DTR_MALLOCFAILED;
  }
  
  // set var bounds
  for (i = 0; i < nDim; i++)
  {
    adblMin[i] = pSrc->aVarDefs[i].bound.dblMin;
    adblMax[i] = pSrc->aVarDefs[i].bound.dblMax;
  } 
  
  // split at node 0
	// note that this gives us the new indexing of linear forms in
	// aNewIndex and the count nNewLFCount of the linear forms needed after
	// optimization
  if ((nErr = Split(0, nMaxNodes, aNodes, &nNodes, aBlocks, &nBlocks, 
                    (*ppDest)->aLinearForms, (*ppDest)->nLinearForms,
                    adblMin, adblMax, adblX, 
                    adblRespMin, adblRespMax, aRespLF,
                    dblBiasBound, adblWBound, dblHalfWidth,
                    nDim, pSrc->nOutputIndex, 0, nMaxDepth, aNewIndex,
                    &nNewLFCount, bSmaller)) != DTR_NOERROR)
  {
    // case where DTREE formation failed
    free(adblMin);
    free(adblMax);
    free(aRespLF);
    free(adblRespMin);
    free(adblRespMax);
    free(adblX);
    free(adblWBound);
    free(aNewIndex);
    DestroyDtreeNodeArray(aNodes);
    DestroyBlockArray(aBlocks, nMaxBlocks);
    DestroyDtree(*ppDest);
    return DTR_MALLOCFAILED;
  }
  // success! clean up some arrays
  free(adblMin);
  free(adblMax);
  free(aRespLF);
  free(adblRespMin);
  free(adblRespMax);
  free(adblX);
  free(adblWBound);
	ASSERT(nNodes == 2 * nBlocks - 1);
    // allocate new space for actual number of blocks
  (*ppDest)->aBlocks = CreateBlockArray(nBlocks);
  if ((*ppDest)->aBlocks == NULL)
  {
    DestroyDtreeNodeArray(aNodes);
    DestroyBlockArray(aBlocks, nBlocks);
    DestroyDtree(*ppDest);
    return DTR_MALLOCFAILED;
  }
  (*ppDest)->nBlocks = nBlocks;
  
  // allocate new space for actual number of nodes
  (*ppDest)->aNodes = CreateDtreeNodeArray(nNodes);
  if ((*ppDest)->aNodes == NULL)
  {
    DestroyDtreeNodeArray(aNodes);
    DestroyBlockArray(aBlocks, nMaxBlocks);
    DestroyDtree(*ppDest);
    return DTR_MALLOCFAILED;
  }
  (*ppDest)->nNodes = nNodes;
  
  // copy blocks
  for (i = 0; i < nBlocks; i++)
  {
    // convert the linear forms to the new index values
    Reindex(aBlocks[i].pMinMaxTree, aNewIndex);
    (*ppDest)->aBlocks[i] = aBlocks[i];
    aBlocks[i].pMinMaxTree = NULL; // mark as empty
  }
	DestroyBlockArray(aBlocks, nMaxBlocks);// aBlocks from Split(0,...) could be large, so get rid of it early
  // copy nodes
  for (i = 0; i < nNodes; i++)
    (*ppDest)->aNodes[i] = aNodes[i];
	DestroyDtreeNodeArray(aNodes); // aNodes could be large, so get rid of it early
        
  // create a downsized array of linear forms that will be needed
  // on the blocks of the optimized multilevel DTREE produced by Split(0,...)

  LINEARFORM* aNewLFArray;
  if((aNewLFArray = CreateLinearFormArray(nNewLFCount, nDim)) == NULL)
  {
    DestroyDtree(*ppDest);
    return DTR_MALLOCFAILED;
  }
  for(i = 0; i < (*ppDest)->nLinearForms; i++)
  {
    // copy the old *ppDest linear form to the array if needed
    // which is the case if aNewIndex[i] is not -1
    if(aNewIndex[i] >= 0)
    {
      for (int j = 0; j < nDim; j++)
      {
        if(j == nOutput) ASSERT( (*ppDest)->aLinearForms[i].adblW[j] == -1);
        (aNewLFArray[aNewIndex[i]]).adblW[j] = (*ppDest)->aLinearForms[i].adblW[j];
        (aNewLFArray[aNewIndex[i]]).adblC[j] = (*ppDest)->aLinearForms[i].adblC[j];
      }
      (aNewLFArray[aNewIndex[i]]).dblBias = (*ppDest)->aLinearForms[i].dblBias;
    }
  }
  // get rid of the old array
	DestroyLinearFormArray((*ppDest)->aLinearForms,(*ppDest)->nLinearForms);
	// make the old pointer point to the new array
	(*ppDest)->aLinearForms = aNewLFArray;
	(*ppDest)->nLinearForms = nNewLFCount; // above changes made on Aug 27, 2007
  // clean up the remaining arrays
  free(aNewIndex);
  return DTR_NOERROR;
}

/////////////////////////////////////////////////////////////////////
// splitting routines

static 
int Split(int nNode, int nMaxNodes, DTREENODE* aNodes, int* pnNodes, 
          BLOCK* aBlocks, int* pnBlocks,
          LINEARFORM* aLF, int nLF, 
          double* adblMin, double* adblMax, double* adblX,
          double* adblRespMin, double* adblRespMax, char* aRespLF,
          double& dblBiasBound, double* adblWBound, double& dblHalfWidth,
          int nDim, int nOutput, int nDepth, int nMaxDepth, int* aNewIndex, int* pnCount, BOOL &bSmaller)
{                      
  double dblMin;
  double dblMax;
  double dblT;
  int nLines;
  int nBlockIndex;
  int nVarIndex;
  // int nReduce;
  int nErr;
  
  if (nDepth == nMaxDepth)
    return DTR_NOERROR;
  
  nBlockIndex = DNODE_BLOCKINDEX((aNodes + nNode)); // cache block index
  
  adblMax[nOutput] = 1e38;
  adblMin[nOutput] = -1e38;
  // optimize tree                     
  do
  {
    bSmaller = FALSE;
    MinMaxNodeBoundCylinder(aBlocks[nBlockIndex].pMinMaxTree, aLF, 
                            adblMin, adblMax,
                            dblBiasBound, adblWBound, dblHalfWidth,
                            nDim, nOutput, bSmaller);
    Amalgamate(aBlocks[nBlockIndex].pMinMaxTree, aLF, 
                             nDim, nOutput, bSmaller);
  } while(bSmaller == TRUE); // keep doing it if the tree gets smaller

  MapNewLinearForms(aBlocks[nBlockIndex].pMinMaxTree, pnCount, aNewIndex);
  if (*pnNodes >= nMaxNodes) // should actually never be greater
  {
    return DTR_NOERROR;                 
    // stop splitting; maximum number of nodes reached
  }
  nLines = 0;
  CountMinMaxTreeLines(aBlocks[nBlockIndex].pMinMaxTree, nLines);
  if (nLines <= nDim)
  {
    return DTR_NOERROR;                 
    // stop splitting; block has no more than nDim linear pieces
  }
                        
  // set responsibility
  SetResp(aBlocks[nBlockIndex].pMinMaxTree,
          aLF, nLF, adblMin, adblMax, adblX, adblRespMin, adblRespMax,
          aRespLF, nDim, nOutput, nLines);
  
  // find best split variable and threshold
  FindBestSplit(&nVarIndex, &dblT, adblMin, adblMax, adblRespMin, adblRespMax, aRespLF,
  nLF, nDim, nOutput, nLines);
  if (nVarIndex == -1) 
  {
    // no good split found
    return DTR_NOERROR;               // could call for  a better split sampling
  }
  // we split further
  dblMin = adblMin[nVarIndex];
  dblMax = adblMax[nVarIndex];
  aNodes[nNode].nLeaf = 0;
  DNODE_VARINDEX(aNodes + nNode) = nVarIndex;
  DNODE_THRESHOLD(aNodes + nNode) = dblT;
  DNODE_LEFTINDEX(aNodes + nNode) = *pnNodes;
  DNODE_RIGHTINDEX(aNodes + nNode) = *pnNodes + 1;
  // left child
  aNodes[*pnNodes].nLeaf = 1;
  aNodes[*pnNodes].nParentIndex = nNode;
  DNODE_BLOCKINDEX(aNodes + *pnNodes) = nBlockIndex;     // left child keeps same min max tree
  aBlocks[nBlockIndex].nDtreeIndex = *pnNodes;
  // right child
  aNodes[*pnNodes + 1].nLeaf = 1;              
  aNodes[*pnNodes + 1].nParentIndex = nNode;
  DNODE_BLOCKINDEX(aNodes + *pnNodes + 1) = *pnBlocks;
  aBlocks[*pnBlocks].nDtreeIndex = *pnNodes + 1;
  aBlocks[*pnBlocks].pMinMaxTree = CopyMinMaxNode(aBlocks[nBlockIndex].pMinMaxTree); 
  if (aBlocks[*pnBlocks].pMinMaxTree == NULL)
    return DTR_MALLOCFAILED;

  // inc number of blocks and nodes  
  *pnBlocks += 1;
  *pnNodes += 2;
  
  // split left child
  adblMax[nVarIndex] = dblT;   // decrease max
  nErr = Split(DNODE_LEFTINDEX(aNodes + nNode), nMaxNodes, aNodes, pnNodes,
               aBlocks, pnBlocks, aLF, nLF,
               adblMin, adblMax, adblX, 
               adblRespMin, adblRespMax, aRespLF,
               dblBiasBound, adblWBound, dblHalfWidth,
               nDim, nOutput, nDepth + 1, nMaxDepth, aNewIndex, pnCount, bSmaller);
  adblMax[nVarIndex] = dblMax; // restore max                               
  if (nErr != DTR_NOERROR)
    return nErr;
  
  // split right child 
  adblMin[nVarIndex] = dblT;   // increase min
  nErr = Split(DNODE_RIGHTINDEX(aNodes + nNode), nMaxNodes, aNodes, pnNodes,
               aBlocks, pnBlocks, aLF, nLF,
               adblMin, adblMax, adblX, 
               adblRespMin, adblRespMax, aRespLF,
               dblBiasBound, adblWBound, dblHalfWidth,
               nDim, nOutput, nDepth + 1, nMaxDepth, aNewIndex, pnCount, bSmaller);
  adblMin[nVarIndex] = dblMin; // restore min                              
  if (nErr != DTR_NOERROR)
    return nErr;
    
  return DTR_NOERROR;
}

static 
void FindBestSplit(int* pnVarIndex, double* pdblT, 
                   double* adblMin, double* adblMax, 
                   double* adblRespMin, 
                   double* adblRespMax, char* aRespLF,
                   int nLF, int nDim, int nOutput, int nLines)
{ 
  // best threshold over all vars is when we have
  // min(max(Nl + No, Nr + No)
  // Nl is number of linear pieces whose resp is entirely to the left of the threshold
  // Nr is number of linear pieces whose resp is entirely to the right of the threshold
  // No is number of linear pieces whose resp overlaps the threshold
  int nLeft, nRight, nOverlap;
  int nMinMax;
  int i;
  double dblSplit;       
  *pnVarIndex = -1;
  for (i = 0; i < nDim; i++)
  {    
    int nMax;
    if (i == nOutput) continue;
    
    FindBestT(i, &dblSplit, adblMin, adblMax, adblRespMin, adblRespMax,
              aRespLF, nLF, nDim, nOutput, nLines, &nLeft, &nRight);
    nOverlap = nLines - nLeft - nRight;
    nMax = __max(nLeft + nOverlap, nRight + nOverlap);
    if (*pnVarIndex == -1)      // first time through
    {
      nMinMax = nMax;
      *pnVarIndex = i;
      *pdblT = dblSplit;
    }
    else if (nMinMax > nMax)
    {
      nMinMax = nMax;
      *pnVarIndex = i;
      *pdblT = dblSplit;
    }
  }  
  if(nMinMax == nLines)
  {
    // no good split found above (ie, nLeft == nRight == 0 for all dimensions)
    if(nLines >= 2 * nDim)
    {
      // failing a good split, we only split blocks in half if they have
      // at least 2 * nDim linear pieces
      i= nOutput; // previously, we just gave up and returned *pnVarIndex = -1
      // pick a random axis, not the output
      while(i == nOutput)
      {
        i = (int) (nDim * ALNRandFloat());
      }
      *pnVarIndex = i;
      ASSERT((0 <= i) && (i <= nDim -1) && (i != nOutput));
      // split in the middle
      *pdblT = 0.5 * (adblMin[i] + adblMax[i]);
    }
    else
    {
      *pnVarIndex = -1; // stops splitting this block
    }
  }
}

static 
void FindBestT(int nVarIndex, double* pdblT, 
               double* adblMin, double* adblMax, 
               double* adblRespMin, 
               double* adblRespMax, char* aRespLF, 
               int nLF, int nDim, int nOutput, int nLines, 
               int* pnLeft, int* pnRight)
{
  // best threshold for this var, over all thresholds
  // min(max(Nl + No, Nr + No)
  // Nl is number of lines whose resp is entirely to the left of the threshold
  // Nr is number of lines whose resp is entirely to the right of the threshold
  // No is number of lines whose resp overlaps the threshold
  
  int nLO, nRO;
  double dblSplitMin = adblMin[nVarIndex];
  double dblSplitMax = adblMax[nVarIndex];
  const int nMaxSteps = 32;     // max steps of binary search
  int nStep = 0;                // current step
  do
  {
    *pdblT = 0.5 * (dblSplitMax + dblSplitMin);
    // count lines to left and right
    CountLeftRightResp(*pdblT, pnLeft, pnRight, aRespLF, adblRespMin, adblRespMax, 
                       nVarIndex, nDim, nLF);
    // minimize the max    
    nLO = nLines - *pnRight;
    nRO = nLines - *pnLeft;
    if (nLO > nRO)      
      dblSplitMax = *pdblT; // new right boundary to the search interval
    else if (nLO < nRO) 
      dblSplitMin = *pdblT; // new left boundary to the search interval
    
    nStep++;
  } while((nLO != nRO) && (nStep < nMaxSteps));
  *pdblT = 0.5 *(dblSplitMax + dblSplitMin);
  // count lines to left and right
  CountLeftRightResp(*pdblT, pnLeft, pnRight, aRespLF, 
                     adblRespMin, adblRespMax, nVarIndex, nDim, nLF);
}

static
void SetResp(MINMAXNODE* pMMN, LINEARFORM* aLF, int nLF,
             double* adblMin, double* adblMax, double* adblX,
             double* adblRespMin, double* adblRespMax,
             char* aRespLF, int nDim, int nOutput, int nLines)
{   
  int i;
  int nPoints;
  // clear resp flags
  memset(aRespLF, 0, MAPBYTECOUNT(nLF));
    
  // sample space and set up responsibility
  // nPoints = __min(5 * nDim * nLines, 500); changed Nov. 4, 2004 (see notes)
  double dblPoints = pow((double)4, (double) nDim) * (double) nLines;
  nPoints = (dblPoints > 8192)?8192:(int) dblPoints;
  for (i = 0; i < nPoints; i++)
  {
    int j;
    int nRespIndex;
    double dblResult;
    for (j = 0; j < nDim; j++)
    {    
      double dblFactor = (double)ALNRandFloat();
      adblX[j] = adblMin[j] + dblFactor * (adblMax[j] - adblMin[j]);
    }
    EvalMinMaxTree(pMMN, aLF, nDim, nOutput,
                   adblX, &dblResult, &nRespIndex);
    for (j = 0; j < nDim; j++)
    { 
      int nIndex;           
      if (j == nOutput) continue;
      nIndex = nRespIndex * nDim + j;
      if (!TESTMAP(aRespLF, nRespIndex))
      {
        adblRespMin[nIndex] = adblX[j];
        adblRespMax[nIndex] = adblX[j];
      }
      else 
      {
        if (adblRespMin[nIndex] > adblX[j])
          adblRespMin[nIndex] = adblX[j];
        if (adblRespMax[nIndex] < adblX[j])
          adblRespMax[nIndex] = adblX[j];
      }
    }
    SETMAP(aRespLF, nRespIndex);
  }
}   

static
void CountLeftRightResp(double dblT, int* pnLeft, int* pnRight, char* aRespLF, 
                        double* adblRespMin, double* adblRespMax, int nVarIndex,
                        int nDim, int nLF)
{     
  int i;
  *pnLeft = 0;
  *pnRight = 0;
  double dblMargin;
  for(i = 0; i < nLF; i++)
  {   
    int nIndex;
    if (!TESTMAP(aRespLF, i)) continue;   // find next resp line
    nIndex = i * nDim + nVarIndex;
    // we enlarge the interval covered by the samples
    // hoping to cover the entire set where the linear piece is active.
    // This is too small for nDim large, but not enlarging is surely
    // an error, since the limits of the piece are not at random sample points
    dblMargin = 0.1 * (adblRespMax[nIndex] - adblRespMin[nIndex]);
    ASSERT(dblMargin >=0);
    if (adblRespMax[nIndex] + dblMargin < dblT)
      *pnLeft += 1;
    else if (adblRespMin[nIndex] - dblMargin > dblT)
      *pnRight += 1;
  }
}

static    
void CountMinMaxTreeLines(MINMAXNODE* pMMN, int& nCount)
{
  if (pMMN->nType == DTREE_LINEAR)
  {
    nCount ++;
  }
  else
  {
    MINMAXNODE* pList = MMN_CHILDLIST(pMMN);
    while(pList != NULL)
    {
      CountMinMaxTreeLines(pList, nCount);
      pList = pList->pNext;
    }
  }
}
 
static    
void MapNewLinearForms(MINMAXNODE* pMMN, int* pnCount,int* aNewIndex)
{
  // this accesses all the linear forms in one final DTREE block's minmax tree
  // and, if a new index has not already been assigned to that linear form,
  // it assigns a new index to be used in aNewLFArray
  if (pMMN->nType == DTREE_LINEAR)
  {
    if(aNewIndex[pMMN->info.nLFIndex] < 0) // -1 is the initial value
    {
      // assign a new index to the linear form
      aNewIndex[pMMN->info.nLFIndex] = *pnCount;
      (*pnCount)++;
    }
  }
  else
  {
    MINMAXNODE* pList = MMN_CHILDLIST(pMMN);
    while(pList != NULL)
    {
      MapNewLinearForms(pList, pnCount, aNewIndex);
      pList = pList->pNext;
    }
  }
}

static
void Reindex(MINMAXNODE* pMMN, int* aNewIndex)
{
  // this recursively accesses all the leaf nodes in the minmax tree of a block
  // and assigns the new index value in aNewIndex to the linear forms
  // This does not alter the DTREE tree structure
  if (pMMN->nType == DTREE_LINEAR)
  {
    if(aNewIndex[pMMN->info.nLFIndex] >= 0)
    {
      pMMN->info.nLFIndex = aNewIndex[pMMN->info.nLFIndex];
    }
  }
  else
  {
    MINMAXNODE* pList = MMN_CHILDLIST(pMMN);
    while(pList != NULL)
    {
      Reindex(pList, aNewIndex);
      pList = pList->pNext;
    }
  }
}

/////////////////////////////////////////////////////////////////////
// setting cylinder bounds

static
void MinMaxNodeBoundCylinder(MINMAXNODE* pMMN, LINEARFORM* aLF, 
                             double* adblMin, double* adblMax,
                             double& dblBiasBound, double* adblWBound, double& dblHalfWidth,
                             int nDim, int nOutput, BOOL &bSmaller)
{
  // This routine updates the max and min and the linear function bounds of a node
  // and passes those bounds back to the caller.
  // In the case of a linear function node, it returns the bias and weights of a linear function
  // and a zero value for the dblHalfWidth in the third line of the parameter declarations.
  // In the case of a min max node it does some elimination of useless children
  // and returns the max and min and linear bounds on the node's function values.
  // There must be at least one child of the node left when it has done all eliminations,
  // and if there is only one, the node is useless and the child takes its place in the tree.
  if(pMMN->nType == DTREE_LINEAR)
  {
    //fprintf(2,"Leaf node reached\n");
    LinearFormBoundCylinder(aLF + MMN_LFINDEX(pMMN), adblMin, adblMax,
      dblBiasBound, adblWBound,
      nDim, nOutput);
      dblHalfWidth = 0;
    return;
  }
  // Accumulators for best bounds so far determined in the sequence of
  // children of this node are declared below. The accumulators contain
  // global bounds  only after the contributions
  // from all children have been received in a previous pass.
  // The first pass gets the min and max global bounds,
  // the second pass gets global linear bounds (maybe improving the
  // global max and min in the process).
  // If a child is certain to be eliminated, its bounds are reset to 1e38 or -1e38 so
  // elimination is faster in the third pass.
  // The third pass may actually eliminate some children.  Changes in the tree structure
  // may follow.
  double dblMin;
  double dblMax;
  double dblBestBiasBound;
  double* adblBestWBound;
  double dblBestHalfWidth;
  // storage of children's bounds
  double* adblChildMin;
  double* adblChildMax;
  double* adblChildBiasBound;
  double* adblChildWBound;
  double* adblChildHalfWidth;
	// count the children of node pMMN
  MINMAXNODE* pList = MMN_CHILDLIST(pMMN);
	int numChild = 0;
	while(pList != NULL)
	{
		pList = pList->pNext;
		numChild++;
	}
  // fprintf(2,"Min/max node reached with %d children\n",numChild);
  // initialize the best bounds achieved so far for this min/max node
  dblMin = -1e38;
  dblMax = 1e38;
  dblBestBiasBound = 0;
  dblBestHalfWidth = 0;
  // allocate array which accumulates W for the best linear bound
  adblBestWBound = (double*) malloc(nDim * sizeof(double));
  // allocate two arrays which will store min and max bounds on child values
	adblChildMin = (double*) malloc(numChild * sizeof(double));
	adblChildMax = (double*) malloc(numChild * sizeof(double));
  // allocate three arrays which store linear function bounds on chilren
  adblChildBiasBound = (double*) malloc(numChild * sizeof(double));
  adblChildWBound = (double*) malloc(numChild * nDim * sizeof(double));
  adblChildHalfWidth = (double*) malloc(numChild * sizeof(double));
  if(adblBestWBound == NULL ||
      adblChildMin == NULL ||
      adblChildMax == NULL ||
      adblChildBiasBound == NULL ||
      adblChildHalfWidth == NULL)
  {
    free(adblBestWBound);
    free(adblChildMin);
	  free(adblChildMax);
    free(adblChildWBound);
    free(adblChildBiasBound);
    free(adblChildWBound);
    free(adblChildHalfWidth);
  }
  // LOOP 1 of 2
  // fprintf(2,"Starting loop 1 of 2 with %d children\n", numChild);
  // fflush(2);
	// set up a counter for traversing the children
	int counter = 0;
	// reinitialize the pointer to the first child
  pList = MMN_CHILDLIST(pMMN);
  while(pList != NULL)
  {
    int nIgnore; // we update linear bounds during loop 1, but we don't eliminate nodes
    // get the bounds from the child at pList
    double dblTempBiasBound = 0;
    double dblTempHalfWidth = 0;
    MinMaxNodeBoundCylinder(pList, aLF, adblMin, adblMax,
      dblTempBiasBound, adblChildWBound + counter * nDim, dblTempHalfWidth,
      nDim, nOutput, bSmaller);
    adblChildBiasBound[counter] = dblTempBiasBound;
    adblChildHalfWidth[counter] = dblTempHalfWidth;
    // at the start, initialize dblMin and dblMax and the linear bounds
    if (pList == MMN_CHILDLIST(pMMN))
    {
      ASSERT(counter == 0);
      dblMin = adblChildMin[counter] = adblMin[nOutput];  // get the values from MinMaxNodeBoundCylinder
      dblMax = adblChildMax[counter] = adblMax[nOutput];  // the counter should be 0 here
      // we can initialize the best bound accumulator using this child
      dblBestBiasBound = adblChildBiasBound[counter];
      adblChildWBound[counter * nDim + nOutput] = -1;
      for(int ii=0; ii < nDim; ii++)
      {
        // if(ii == nOutput)continue; not needed now
        adblBestWBound[ii] = adblChildWBound[counter * nDim + ii];
      }
      // adblBestWBound[nOutput] = -1; not needed now
      dblBestHalfWidth = adblChildHalfWidth[counter];
    }
    else // this is not the first child on the list, work towards global bounds
    {
      ASSERT(counter < numChild);
      // update the min and max arrays for the children's values
      adblChildMin[counter] = adblMin[nOutput];
      adblChildMax[counter] = adblMax[nOutput];
      // dblMin and dblMax are tracking the best min and max found so far for the children
      if (pMMN->nType == DTREE_MIN) 
      {   
        // update dblMin and dblMax if they are less for the child
        if (adblMax[nOutput] < dblMax) dblMax = adblMax[nOutput];
        if (adblMin[nOutput] < dblMin) dblMin = adblMin[nOutput]; 
        // update the linear bounds
        // subtle point: we can't use a new max and min bound to eliminate any
        // node contributing to it, so we use dblTempMin, dblTempMax to absorb any changes
        nIgnore = aboveORbelow(adblMin, adblMax,
          dblMin, dblMax, adblChildMin[counter], adblChildMax[counter],
          dblBestBiasBound, adblBestWBound, dblBestHalfWidth,
          adblChildBiasBound[counter], adblChildWBound + counter * nDim, adblChildHalfWidth[counter],
          nDim,nOutput, DTREE_MIN);
      }
      else // it's a maximum node
      {
        // update dblMin and dblMax if they are greater for the child
        if (adblMin[nOutput] > dblMin) dblMin = adblMin[nOutput];
        if (adblMax[nOutput] > dblMax) dblMax = adblMax[nOutput];
        nIgnore = aboveORbelow(adblMin, adblMax,
          dblMin, dblMax, adblChildMin[counter], adblChildMax[counter],
          dblBestBiasBound, adblBestWBound, dblBestHalfWidth,
          adblChildBiasBound[counter], adblChildWBound + counter * nDim, adblChildHalfWidth[counter],
          nDim,nOutput, DTREE_MAX);
      }
    }
    pList = pList->pNext;
		counter++;
  }
  // Up to now we have found global max and min bounds which must bound
  // any values produced by this node, and we have a good starting point
  // to accumulate a best global linear bound.
  // LOOP 2 of 2
  // fprintf(2,"Starting loop 2 of 2 with %d children\n", numChild);
  // fflush(2);
  // now we go through the children again, eliminating some using
  // both the global min and max bounds and the global linear bounds obtained in loop 1
	// reinitialize the pointer to the first child and the counter to 0
  // now we need a prior pointer so we can manipulate the tree structure
  // we keep track of whether the tree got smaller so we can then repeat MinMaxNodeBoundCylinder
  pList = MMN_CHILDLIST(pMMN);
  MINMAXNODE* pPrior = pList;
	counter = 0;
  while(pList != NULL)
  {
		if(pList == MMN_CHILDLIST(pMMN)) // if this is the first on the list
		{
      if(((pMMN->nType == DTREE_MIN) && (adblChildMin[counter] > dblMax))
			|| ((pMMN->nType == DTREE_MAX) && (adblChildMax[counter] < dblMin)))
			{
				// delete this child
			  MMN_CHILDLIST(pMMN) = pList->pNext;
				DestroyMinMaxNode(pList);
				pPrior = pList = MMN_CHILDLIST(pMMN);
        bSmaller = TRUE;
			}
			else // check the linear bounds
      {
        if(pMMN->nType == DTREE_MIN)
        {
          if(aboveORbelow(adblMin, adblMax,
               dblMin, dblMax, adblChildMin[counter], adblChildMax[counter],
               dblBestBiasBound, adblBestWBound, dblBestHalfWidth,
               adblChildBiasBound[counter], adblChildWBound + counter * nDim, adblChildHalfWidth[counter],
               nDim,nOutput, DTREE_MIN) == CHILDABOVEBESTBOUND) // the child is too great to be active
          {
	          // delete this child
	          MMN_CHILDLIST(pMMN) = pList->pNext;
				    DestroyMinMaxNode(pList);
				    pPrior = pList = MMN_CHILDLIST(pMMN);
            bSmaller = TRUE;
          }
          else
          {
            // move on
            pPrior = pList;
            pList = pList->pNext;     
          }
        }
        else  // pMMN->nType == DTREE_MAX
        {
          if(aboveORbelow(adblMin, adblMax,
               dblMin, dblMax, adblChildMin[counter], adblChildMax[counter],
               dblBestBiasBound, adblBestWBound, dblBestHalfWidth,
               adblChildBiasBound[counter], adblChildWBound + counter * nDim, adblChildHalfWidth[counter],
               nDim,nOutput,DTREE_MAX) == CHILDBELOWBESTBOUND)
          {
	          // delete this child
	          MMN_CHILDLIST(pMMN) = pList->pNext;
				    DestroyMinMaxNode(pList);
				    pPrior = pList = MMN_CHILDLIST(pMMN);         
            bSmaller = TRUE;
          }
          else
          {
	          //  move on
	          pPrior = pList;
	          pList = pList->pNext;
          }
        }
      } // end of checking linear bounds
    } // end of processing for first on the list
    else // this is not the first on the list
    {
      if(((pMMN->nType == DTREE_MIN) && (adblChildMin[counter] > dblMax))
			|| ((pMMN->nType == DTREE_MAX) && (adblChildMax[counter] < dblMin)))
			{
				// delete this child
				pPrior->pNext = pList->pNext;
				DestroyMinMaxNode(pList);
				pList = pPrior->pNext;        
        bSmaller = TRUE;
			}
			else // check the linear bounds
      {
        if(pMMN->nType == DTREE_MIN)
        {
          if(aboveORbelow(adblMin, adblMax,
               dblMin, dblMax, adblChildMin[counter], adblChildMax[counter],
               dblBestBiasBound, adblBestWBound, dblBestHalfWidth,
               adblChildBiasBound[counter], adblChildWBound + counter * nDim, adblChildHalfWidth[counter],
               nDim,nOutput, DTREE_MIN) == CHILDABOVEBESTBOUND) // the child is too great to be active
          {
	          // delete this child
	          pPrior->pNext = pList->pNext;
	          DestroyMinMaxNode(pList);
	          pList = pPrior->pNext;     
            bSmaller = TRUE;
          }
          else
          {
            // move on
            pPrior = pList;
            pList = pList->pNext;
          }
        }
        else  // pMMN->nType == DTREE_MAX
        {
          if(aboveORbelow(adblMin, adblMax,
               dblMin, dblMax, adblChildMin[counter], adblChildMax[counter],
               dblBestBiasBound, adblBestWBound, dblBestHalfWidth,
               adblChildBiasBound[counter], adblChildWBound + counter * nDim, adblChildHalfWidth[counter],
               nDim,nOutput,DTREE_MAX) == CHILDBELOWBESTBOUND)
          {
            // delete this child
	          pPrior->pNext = pList->pNext;
	          DestroyMinMaxNode(pList);
	          pList = pPrior->pNext;              
            bSmaller = TRUE;
          }
          else
          {
	          //  move on
	          pPrior = pList;
	          pList = pList->pNext;
          }
        }
      } // end of checking linear bounds
    } // end of processing when this is not the first on the list
  	counter++;
  } 
	// if pMMN has only one child left, it becomes the child
	while(pMMN->nType != DTREE_LINEAR && MMN_CHILDLIST(pMMN)->pNext == NULL)
	{   
    ASSERT(bSmaller == TRUE); // the only way to get here is to have eliminated a child
		MINMAXNODE* pChild = MMN_CHILDLIST(pMMN);     // save pointer to (singular) child list
		MINMAXNODE* pNext = pMMN->pNext;              // save pointer to sibling
		memcpy(pMMN, pChild, sizeof(MINMAXNODE));     // copy child info into this
		pMMN->pNext = pNext;                          // restore sibling
		MMN_CHILDLIST(pChild) = NULL;                 // mark child's child list as empty before destroy
    ASSERT(pChild != NULL);
		DestroyMinMaxNode(pChild);
	}
  // return the best minimum and maximum bounds on this node
  adblMin[nOutput] = dblMin;
  adblMax[nOutput] = dblMax;
  // return the best linear bound
  dblBiasBound = dblBestBiasBound;
  adblBestWBound[nOutput] = -1;
  for(int ii=0; ii < nDim; ii++)
  {
    // if(ii == nOutput)continue; not needed
    adblWBound[ii] = adblBestWBound[ii];
  }
  dblHalfWidth = dblBestHalfWidth;

  //fprintf(2,"Freeing arrays for MinMaxNodeBoundCylinder\n");
	free(adblChildMin);
	free(adblChildMax);
  free(adblChildBiasBound);
  free(adblChildWBound);
  free(adblChildHalfWidth);
  free(adblBestWBound);
}

static 
void LinearFormBoundCylinder(LINEARFORM* pLF, 
                             double* adblMin, double* adblMax,
                             double& dblBiasBound, double* adblWBound,
                             int nDim, int nOutput)
{ 
  // This routine passes up the bounds on a linear piece, including
  // the minimum and maximum values on the piece as well as data giving
  // the upper and lower linear function bounds.
  // The min and the max are recorded in the nOutput component
  // of the adblMin and adblMax arrays, the other components being the box bounds.
  // Linear function bounds are expressed as a linear piece with a centroid
  // equal to the centre of the current box. The box is given by the values in
  // adblMin and adblMax except for the nOutput index. The linear piece is shifted up by the half width
  // to get an upper bound and down by the half width to get a parallel lower bound.
  // The half width is 0 in the case of a linear form, since the upper and lower linear
  // bounds are equal.  Hence, there is no need to pass this back as a parameter.
  double dblMin, dblMax, dblCentroidThisDim, dblChange;
  // calc min and max of the intersection of the linear piece with the cylinder
  // formed by the max and min array values in the input directions
  dblMin = dblMax = dblChange = 0;
  for (int i = 0; i < nDim; i++)
  {
    if (i == nOutput) continue;
    adblWBound[i] = pLF->adblW[i];
    dblCentroidThisDim = 0.5 * (adblMin[i] + adblMax[i]);
    dblChange += adblWBound[i] * (dblCentroidThisDim - pLF->adblC[i]); // change due to changed centroid
    dblMin += adblWBound[i] * (((adblWBound[i] > 0) ? adblMin[i] : adblMax[i]) - dblCentroidThisDim);
    dblMax += adblWBound[i] * (((adblWBound[i] > 0) ? adblMax[i] : adblMin[i]) - dblCentroidThisDim);
  }
  adblWBound[nOutput] = -1;
  dblBiasBound = pLF->adblC[nOutput] + dblChange; 
  //dblMin /= -pLF->adblW[nOutput]; // normally this should be a division by 1
  //dblMax /= -pLF->adblW[nOutput]; // since adblW[nOutput] == -1; but let's leave this in anyway
  // later, we may be able to use this to extract inverses without forming a new ALN.
  adblMin[nOutput] = dblMin + dblBiasBound;
  adblMax[nOutput] = dblMax + dblBiasBound;
}

static
int aboveORbelow(double* adblMin, double* adblMax,
                 double& dblMin1, double& dblMax1, double dblMin2, double dblMax2,
                 double& dblBias1, double* adblW1, double& dblHalfWidth1,
                 double dblBias2, double* adblW2, double dblHalfWidth2,
                 int nDim, int nOutput, int nMinMaxType)
{
  // The minmax types are DTREE_MAX (= 1) and DTREE_MIN (= 0) and DTREE_LINEAR (= 2, not used here).
  // This code is mainly written for the case DTREE_MIN, but the beginning and end parts
  // turn the functions upside-down to handle the case of a DTREE_MAX.
  // Bound 1 is the best bound found so far.  Bound 2 is the child bound.
  // The bounds are all computed relative to the box centre (which changes in each DTREE block).
  // If the bound 2 is better than bound 1, it replaces bound 1. Otherwise a combination
  // of the bounds replaces bound 1. Routine aboveORbelow doesn't allow equality
  //  elimination, so we don't delete a child possibly contributing to any best bound.
  // The variables with X in them are used to handle both max and min nodes,
  // the max nodes being inverted. Only the nOutput component is affected
  // by inversions (sign change, min becomes max and vice-versa).
  double dblXMin1, dblXMax1, dblXMin2, dblXMax2;
  double dblXBias1, dblXHalfWidth1, dblXBias2, dblXHalfWidth2;
  double dblCentroidThisDim;
  double* adblXW1 = (double*) malloc(nDim * sizeof(double));
  double* adblXW2 = (double*) malloc(nDim * sizeof(double));
  double* adblC1 = (double*) malloc(nDim * sizeof(double)); // corner where B1-B2 is greatest
  double* adblC2 = (double*) malloc(nDim * sizeof(double)); //  or least
  int nReturn; // this is the value that will be returned (after inversion for a max node)
  nReturn = BOUNDSOVERLAP; // set this as default
  // fprintf(2,"Node type %s *************************  \nInputs:\n",
  //  (nMinMaxType == 0?"DTREE_MIN":"DTREE_MAX"));
  if(nMinMaxType == DTREE_MIN)
  {
    dblXMin1 = dblMin1; // the current min and max bounds on input
    dblXMax1 = dblMax1;
    dblXMin2 = dblMin2;   // the invariable current child min and max
    dblXMax2 = dblMax2;
    dblXBias1 = dblBias1; // the current bias for the best linear bound on input
    dblXBias2 = dblBias2;   // the invariable child bias
    dblXHalfWidth1 = dblHalfWidth1; // the current half-width for the best linear bound on input
    dblXHalfWidth2 = dblHalfWidth2;   // the invariable child half-width
    adblXW1[nOutput] = -1;
    adblXW2[nOutput] = -1;
    for(int ii = 0; ii < nDim; ii++)
    {
     if(ii == nOutput) continue;
     adblXW1[ii] = adblW1[ii];   // the initial weights for the best linear bound
     adblXW2[ii] = adblW2[ii];   // the invariable weights on the child
    }
  }
  else // the node is type DTREE_MAX, and we invert the output dimension
  {
    dblXMin1 = -dblMax1;
    dblXMax1 = -dblMin1;
    dblXMin2 = -dblMax2;
    dblXMax2 = -dblMin2;
    dblXBias1 = -dblBias1;
    dblXBias2 = -dblBias2;
    dblXHalfWidth1 = dblHalfWidth1; // the signs don't change
    dblXHalfWidth2 = dblHalfWidth2; // for the half-widths
    adblXW1[nOutput] = -1;
    adblXW2[nOutput] = -1;
    for(int ii = 0; ii < nDim; ii++)
    {
     if(ii == nOutput) continue;
     adblXW1[ii] = -adblW1[ii];
     adblXW2[ii] = -adblW2[ii];
    }
  }
  // store the input values for checking
  double dblYMin1, dblYMax1, dblYBias1, dblYHalfWidth1;
  double* adblYW1 = (double*) malloc(nDim * sizeof(double));
  dblYMin1 = dblXMin1;
  dblYMax1 = dblXMax1;
  dblYBias1 = dblXBias1;
  dblYHalfWidth1 = dblXHalfWidth1;
  for(int ii = 0; ii < nDim; ii++)
  {
   adblYW1[ii] = adblXW1[ii];
  }
  /*int nSection; indicator of active part for debugging
  // print out all values *****************************
  fprintf(2,"Child bounds before & after for index 0 dblXMin2= %f, dblXMax2 = %f\n",dblXMin2,dblXMax2);
  fprintf(2,"dblXBias1 = %f, dblXBias2 = %f, dblXHalfWidth1 = %f, dblXHalfWidth2 = %f\n",
                      dblXBias1, dblXBias2, dblXHalfWidth1, dblXHalfWidth2); 
  if(nDim == 2)
  {
    dblCentroidThisDim = 0.5 * (adblMin[0] + adblMax[0]);
    fprintf(2,"Child upper linear bound before left value for index 0 %f,          right value = %f\n",
              dblXBias2 + dblXHalfWidth2 + adblXW2[0] * (adblMin[0] - dblCentroidThisDim),
              dblXBias2 + dblXHalfWidth2 + adblXW2[0] * (adblMax[0] - dblCentroidThisDim));
    fprintf(2,"Child lower linear bound before left value for index 0 %f,          right value = %f\n",
              dblXBias2 - dblXHalfWidth2 + adblXW2[0] * (adblMin[0] - dblCentroidThisDim),
              dblXBias2 - dblXHalfWidth2 + adblXW2[0] * (adblMax[0] - dblCentroidThisDim));
    fprintf(2,"Best min/max bounds before for index 0 dblXMin1= %f, dblXMax1= %f\n",dblXMin1,dblXMax1);
    fprintf(2,"Best upper linear bound before for index 0:     left value %f,        right value = %f\n",
              dblXBias1 + dblXHalfWidth1 + adblXW1[0] * (adblMin[0] - dblCentroidThisDim),
              dblXBias1 + dblXHalfWidth1 + adblXW1[0] * (adblMax[0] - dblCentroidThisDim));
    fprintf(2,"Best lower linear bound before for index 0:     left value %f,                          right value = %f\n",
              dblXBias1 - dblXHalfWidth1 + adblXW1[0] * (adblMin[0] - dblCentroidThisDim),
              dblXBias1 - dblXHalfWidth1 + adblXW1[0] * (adblMax[0] - dblCentroidThisDim));
  }*/


  // set up the intermediate variables
  double dblB1C1, dblB1C2, dblB2C1, dblB2C2;
  double dblTempBias, dblTempBiasTilde, dblTempHalfWidth, dblAlpha,
    dblAlphaTilde, dblQ1, dblQ2, dblV1, dblV2;
  // we compare the best and the child upper bounds
  // adblC1[ii] is where the best bound so far (B1)
  // is greatest with respect to the child bound (B2)
  // ie the corner of the input box where B1 - B2 is greatest
  // and adblC2[ii] is where B1 - B2 is least
  dblB1C1 = dblB1C2 = dblXBias1 + dblXHalfWidth1;
  dblB2C1 = dblB2C2 = dblXBias2 + dblXHalfWidth2;
  //dblW1dotCdiff = dblW2dotCdiff = 0;
  for(int ii = 0; ii < nDim; ii++)
  {
    if(ii == nOutput)continue;
    dblCentroidThisDim = 0.5 * (adblMin[ii] + adblMax[ii]);
    if(adblXW1[ii] - adblXW2[ii] > 0)
    {
      adblC1[ii] = adblMax[ii] - dblCentroidThisDim;
      adblC2[ii] = adblMin[ii] - dblCentroidThisDim;
    }
    else
    {
      adblC1[ii] = adblMin[ii] - dblCentroidThisDim;
      adblC2[ii] = adblMax[ii] - dblCentroidThisDim;
    }
    dblB1C1 += adblXW1[ii] * adblC1[ii];
    dblB1C2 += adblXW1[ii] * adblC2[ii];
    dblB2C1 += adblXW2[ii] * adblC1[ii];
    dblB2C2 += adblXW2[ii] * adblC2[ii];
  }
  adblC1[nOutput] = adblC2[nOutput] = 1e38;  // just for safety; value, if used, causes havoc
  // now we analyse which linear bound is bound is better: the current best (1) or the child bound (2)
  if((dblB1C2 - 2 * (1 + 1e-12) * dblXHalfWidth1) > dblB2C2 )  // lower the best bound minimium for this test
  {
    // nSection = 0;
    // lower bound of B1 is always > upper bound of B2 with
    // a certain excess to allow for numerical error.
    // For a Min node, which we are assuming, B2 is better
    // so we copy bound 2 into bound 1 
    dblXBias1 = dblXBias2;
    for(int ii = 0; ii < nDim; ii++)
    {
      if (ii == nOutput) continue;
      adblXW1[ii] = adblXW2[ii];
    }
    dblXHalfWidth1 = dblXHalfWidth2;
    // fprintf(2,"Child upper linear bound is below best lower linear bound; check new best bound\n"); 
    nReturn = CHILDBELOWBESTBOUND;  // means B2 is < B1 (with strict inequality), will be inverted for max
  }
  else
  {
    if(dblB2C1 - 2 * (1 + 1e-12) * dblXHalfWidth2 > dblB1C1 ) // minimum of child bound lowered a bit to make
    // sure of separation even when there are small numerical errors in the values
    {
      // nSection = 1;
      // fprintf(2,"Child lower linear bound is above best linear upper bound (for a MIN); best unchanged\n"); 
      // bound 2 always > bound 1
      nReturn = CHILDABOVEBESTBOUND;  // child 2 can be deleted; return value will be inverted for max
    }
    else // this finishes the cases where the child range is entirely above or below the best range
    {
      nReturn = BOUNDSOVERLAP;
      // in the above cases, the *lower* bound of one is > the upper bound of the other by a small margin 1e -7
      // in what follows, this is not the case and we always return 0 = BOUNDSOVERLAP
      // we next deal with the cases where the *upper* bound of one is >= the upper bound of the other
      // the clear choice for upper bound is then the lower of the two for a minimum node
      if(dblB1C2 >= dblB2C2)
      {
        // nSection = 3;
        // fprintf(2,"Bounds overlap; child upper linear bound is below best upper linear bound (for a MIN)\n"); 
        // in this case bound 2 is everywhere the lesser upper bound and should be used as the best
        dblTempHalfWidth = 0.5 * (dblB2C2 - dblB1C2) + dblXHalfWidth1;
        // return the new best bound after gathering info for check
        dblXHalfWidth1 = __max(dblTempHalfWidth, dblXHalfWidth2);
        dblXBias1 = dblXBias2 + dblXHalfWidth2 - dblXHalfWidth1;
        for(int ii = 0; ii < nDim; ii++)
        {
          if(ii == nOutput)continue;
          adblXW1[ii] = adblXW2[ii];
        }
        // expand the best bound just to make sure numerical errors don't lead to a bad bound
        dblXHalfWidth1 *= (1 + 1e-12);
      }
      else
      {
        if(dblB1C1 <= dblB2C1)
        {
          // nSection = 4;
          // fprintf(2,"Bounds overlap, child upper linear bound is above best linear upper bound (for a MIN)\n"); 
          // in this case bound 1 is the lesser upper bound and should be used as the best
          dblTempHalfWidth = 0.5 * (dblB1C1 - dblB2C1) + dblXHalfWidth2;
          // return the new best bound
          dblXBias1 += dblXHalfWidth1; // old dblXHalfWidth1
          dblXHalfWidth1 = __max(dblTempHalfWidth, dblXHalfWidth1);
          dblXBias1 -= dblXHalfWidth1; // new dblXHalfWidth1
          // expand the best bound just to make sure numerical errors don't lead to a bad bound
          dblXHalfWidth1 *= (1 + 1e-12);
          // the weights returned for W1 are unchanged
          // however the half-width changes possibly and so do the min and max
        }
        else 
        {
          ASSERT((dblB1C1 - dblB2C1) * (dblB1C2 - dblB2C2) <= 0);
          // fprintf(2,"Bounds overlap, upper linear bounds intersect");
          // now we know the upper bounds intersect inside the box
          // we compute the lowest points of the lower bounds relative to the other
          dblV1 = __min(dblB1C1 - 2*dblXHalfWidth1, dblB2C1 - 2*dblXHalfWidth2);
          dblV2 = __min(dblB1C2 - 2*dblXHalfWidth1, dblB2C2 - 2*dblXHalfWidth2);
          if((dblB1C1 - dblB1C2) * ( dblB2C1 - dblB2C2) >= 0)
          {
            // nSection = 5;
            // there is no peak since the slopes have the same sign or one is 0
            // and we can choose one of the two bounds for the best
            // we take the one that has the lesser centroid value (B2 in case of equality)
            if((dblB1C1 + dblB1C2) >= (dblB2C1 + dblB2C2))
            {
              // choose 2 since it has the lesser 
              dblTempBias = dblXBias2 + dblXHalfWidth2;
              dblXHalfWidth1 = 0.5 * __max(dblB2C2 - dblV2, dblB2C1 - dblV1);
              dblXBias1 = dblTempBias - dblXHalfWidth1;
              for(int ii = 0; ii < nDim; ii++)
              {
                if(ii == nOutput)continue;
                adblXW1[ii] = adblXW2[ii];
              }
            }
            else
            {
              // we keep the weights of B1 but change the bias and half-width
              dblTempBias = dblXBias1 + dblXHalfWidth1;
              dblXHalfWidth1 = 0.5 * __max(dblB1C2 - dblV2, dblB1C1 - dblV1);
              dblXBias1 = dblTempBias - dblXHalfWidth1;
            }
            // expand the best bound just to make sure numerical errors don't lead to a bad bound
            dblXHalfWidth1 *= (1 + 1e-12);
          }
          else
          {
            ASSERT((dblB1C1 - dblB1C2) * ( dblB2C1 - dblB2C2) < 0 ); // there is a peak
            // Compute a new maximum using the linear bounds
            //  nSection = 6;
            double dblMaxVal, dblMinVal, dblMaxValTemp, dblMinValTemp, dblAvgWt;
            dblMaxVal = 1e38;
            dblMinVal = -1e38;
            // try tilting to equalize heights of upper bound along various axes
            for(int ii = 0; ii < nDim; ii++)
            {
              if(ii == nOutput) continue;
              // for other dimensions, we solve the following to get
              // several tilts for which we can compute the minimum maximum value of bound
              if((fabs(adblXW1[ii] - adblXW2[ii])) < 3)
              {
                dblAlphaTilde = 0.5; // the slopes of the bounds are so close, just use 50-50
              }
              else
              {
                dblAlphaTilde = -adblXW2[ii] / (adblXW1[ii] - adblXW2[ii]);
              }
              if(dblAlphaTilde > 1.0) dblAlphaTilde = 1.0;
              if(dblAlphaTilde < 0) dblAlphaTilde = 0;
              // but we might at least get an improvement
              dblTempBiasTilde = dblAlphaTilde * (dblXBias1 + dblXHalfWidth1) +
                                 (1.0 - dblAlphaTilde) * (dblXBias2 + dblXHalfWidth2);
              dblMaxValTemp = dblMinValTemp = dblTempBiasTilde;
              
              for(int jj = 0; jj < nDim; jj++)
              {
                if(jj == nOutput) continue; // in direction ii the slope of the bound is 0, so use the centroid
                dblCentroidThisDim = 0.5 * (adblMax[jj] + adblMin[jj]);
                dblAvgWt = dblAlphaTilde * adblXW1[jj] + (1.0 - dblAlphaTilde) * adblXW2[jj];
                if(dblAvgWt > 0)
                {
                  dblMaxValTemp += dblAvgWt * (adblMax[jj] - dblCentroidThisDim); 
                  dblMinValTemp += dblAvgWt * (adblMin[jj] - dblCentroidThisDim);
                }
                else
                {
                  dblMaxValTemp += dblAvgWt * (adblMin[jj] - dblCentroidThisDim); 
                  dblMinValTemp += dblAvgWt * (adblMax[jj] - dblCentroidThisDim);
                }
              }
              // now dblMaxVal for the axis ii is an upper bound
              // if it improves it, we can reduce the dblXMax1, maybe significantly
              // However, because it is not used until the next child, we can't
              // use the new bound to eliminate any child contributing to it
              // Does the new tilt give a lower maximum value?
              // This operates at least once, when dblMaxVal = 1e38
              if(dblMaxValTemp < dblMaxVal)
              {
                dblMaxVal = dblMaxValTemp;
                if(dblMaxVal < dblXMax1) dblXMax1 = dblMaxVal; // to be corrected by later values which are lower
                // fprintf(2,"The maximum of a bound obtained at a peak along axis %d = %f\n", ii, dblMaxVal);
                // fprintf(2,"Better Alpha for tilt = %f\n", dblAlphaTilde);
                dblAlpha = dblAlphaTilde; // this will give us the linear bound with the lowest maximum
                dblTempBias = dblTempBiasTilde; // this gives us the right alpha value
                dblMinVal = dblMinValTemp; // to be corrected by subtracting 2 * final half-width
              }
            }
            // now compute the new linear upper bound using the Alpha from above
            for(int ii = 0; ii < nDim; ii++)
            {
               if(ii == nOutput)continue;
               adblXW1[ii] = dblAlpha * adblXW1[ii] + (1.0 - dblAlpha) * adblXW2[ii]; 
            }
            // this is the only case where we can possibly improve the upper bound on the DTREE_MIN node
            // now we look at all the axes, and get upper bounds using dblAlphaTilde which depends on ii
            // make sure the equation is solved
            dblQ1 = dblAlpha * dblB1C1 + (1.0 - dblAlpha) * dblB2C1;
            dblQ2 = dblAlpha * dblB1C2 + (1.0 - dblAlpha) * dblB2C2;
            dblXHalfWidth1 = 0.5 * __max(dblQ1 - dblV1,dblQ2 - dblV2);
            dblXBias1 = dblTempBias - dblXHalfWidth1;
            // expand the best bound just to make sure numerical errors don't lead to a bad bound
            dblXHalfWidth1 *= (1 + 1e-12);
            if(dblMinVal - 2 * dblXHalfWidth1 > dblXMin1) dblXMin1 = dblMinVal - 2 * dblXHalfWidth1;
          } // ends the block where there is a peak
        } // ends the else where the upper bounds intersect in the box
      } // ends the block where bounds overlap
    } // ends else where the bounds overlap
  } // ends the else where the child is not below the best bound

/*
  fprintf(2, "Outputs of node\n");
  fprintf(2,"Best min/max bounds after (for a MIN): dblXMin1= %f, dblXMax1= %f\n",dblXMin1,dblXMax1);
  fprintf(2,"dblXBias1 = %f, dblXHalfWidth1 = %f\n", dblXBias1, dblXHalfWidth1); 
  if(nDim == 2)
  {
    dblCentroidThisDim = 0.5 * (adblMin[0] + adblMax[0]);
    fprintf(2,"Best upper linear bound (for a MIN in 1D) after                  left value %f,               right value = %f\n",
            dblXBias1 + dblXHalfWidth1 + adblXW1[0] * (adblMin[0] - dblCentroidThisDim),
            dblXBias1 + dblXHalfWidth1 + adblXW1[0] * (adblMax[0] - dblCentroidThisDim));
    fprintf(2,"Best lower linear bound (for a MIN in 1D) after for index 0:     left value %f,                          right value = %f\n",
            dblXBias1 - dblXHalfWidth1 + adblXW1[0] * (adblMin[0] - dblCentroidThisDim),
            dblXBias1 - dblXHalfWidth1 + adblXW1[0] * (adblMax[0] - dblCentroidThisDim));
  }
  fflush(2);

  // check all upper and lower bounds in the corners C1 and C2
  // Y indicates best before values
  // B indicates best values after the routine finishes
  // 2 indicates the child values, presumably unchanged
  double  dblYUC1, dblYLC1, dblYUC2, dblYLC2,
          dbl2UC1, dbl2LC1, dbl2UC2, dbl2LC2,
          dblBUC1, dblBLC1, dblBUC2, dblBLC2;
  dblYUC1 = dblYBias1 + dblYHalfWidth1;
  dblYLC1 = dblYBias1 - dblYHalfWidth1;
  dblYUC2 = dblYBias1 + dblYHalfWidth1;
  dblYLC2 = dblYBias1 - dblYHalfWidth1;
  dblBUC1 = dblXBias1 + dblXHalfWidth1;
  dblBLC1 = dblXBias1 - dblXHalfWidth1;
  dblBUC2 = dblXBias1 + dblXHalfWidth1;
  dblBLC2 = dblXBias1 - dblXHalfWidth1;
  dbl2UC1 = dblXBias2 + dblXHalfWidth2;
  dbl2LC1 = dblXBias2 - dblXHalfWidth2;
  dbl2UC2 = dblXBias2 + dblXHalfWidth2;
  dbl2LC2 = dblXBias2 - dblXHalfWidth2;
  for(ii = 0; ii < nDim; ii++)
  {
    if(ii == nOutput) continue;
    ASSERT(adblC1[ii] * adblC2[ii] < 0);
    dblYUC1 += adblYW1[ii] * adblC1[ii];
    dblYLC1 += adblYW1[ii] * adblC1[ii];
    dblYUC2 += adblYW1[ii] * adblC2[ii];
    dblYLC2 += adblYW1[ii] * adblC2[ii];
    dblBUC1 += adblXW1[ii] * adblC1[ii];
    dblBLC1 += adblXW1[ii] * adblC1[ii];
    dblBUC2 += adblXW1[ii] * adblC2[ii];
    dblBLC2 += adblXW1[ii] * adblC2[ii];
    dbl2UC1 += adblXW2[ii] * adblC1[ii]; 
    dbl2LC1 += adblXW2[ii] * adblC1[ii];
    dbl2UC2 += adblXW2[ii] * adblC2[ii];
    dbl2LC2 += adblXW2[ii] * adblC2[ii];
  }

  if(nSection != 6) // the upper bounds don't work the same way for section 6
  {
    ASSERT(__min(dblYUC1, dbl2UC1) <= dblBUC1 + 0.0001);
    ASSERT(__min(dblYUC2, dbl2UC2) <= dblBUC2 + 0.0001);
  }
  ASSERT(__min(dblYLC1, dbl2LC1) >= dblBLC1 - 0.0001); // section 6 failure
  ASSERT(__min(dblYLC2, dbl2LC2) >= dblBLC2 - 0.0001);
*/


  // now possibly invert the results before the return
  if(nMinMaxType == DTREE_MIN)
  {
    // improvements of the current best bounds (min/max and linear)
    dblMin1 = dblXMin1; 
    dblMax1 = dblXMax1; 
    dblBias1 = dblXBias1;
    dblHalfWidth1 = dblXHalfWidth1;
    for(int ii = 0; ii < nDim; ii++)
    {
      adblW1[ii] = adblXW1[ii];
    }
  }
  else // the node is type DTREE_MAX: invert the output dimension
  {
    dblMin1 = -dblXMax1; 
    dblMax1 = -dblXMin1;
    dblBias1 = -dblXBias1;
    dblHalfWidth1 = dblXHalfWidth1; // the sign doesn't change
    for(int ii = 0; ii < nDim; ii++)
    {
     if(ii == nOutput) continue;
     adblW1[ii] = -adblXW1[ii];
    }
    adblW1[nOutput] = -1;
    nReturn = - nReturn; // this inverts whether the child is above/below the current best bound
  }
  free(adblC1);
  free(adblC2);
  free(adblXW1);
  free(adblXW2);
  free(adblYW1);
  return nReturn;
}


static
void Amalgamate(MINMAXNODE* pMMN, LINEARFORM* aLF, 
                   int nDim, int nOutput, BOOL &bSmaller)
{                                     
  MINMAXNODE* pChild;
  MINMAXNODE* pChildPrev;
  // When children of a min or max node are removed, it could be that there is only one child left.
  // In that case, the min/max node is unnecessary, and is eliminated by moving the child up
  // into the place of the node which "becomes the child".  A while-loop makes sure this is done
  // enough times to get rid of the single-child anomaly to any number of levels.
  // Following the above action, it can be that a child is of the same type, max or min,
  // as the node, and so the child's list can be put in place of the child (amalgamation).
  while(pMMN->nType != DTREE_LINEAR && MMN_CHILDLIST(pMMN)->pNext == NULL)
  {   
    MINMAXNODE* pChild = MMN_CHILDLIST(pMMN);     // save pointer to (singular) child list
    MINMAXNODE* pNext = pMMN->pNext;              // save pointer to sibling
    memcpy(pMMN, pChild, sizeof(MINMAXNODE));     // copy child info into this
    pMMN->pNext = pNext;                          // restore sibling
    MMN_CHILDLIST(pChild) = NULL;                 // mark child's child list as empty before destroy
    ASSERT(pChild != NULL);
    DestroyMinMaxNode(pChild);
    bSmaller = TRUE;
  }
  if(pMMN->nType != DTREE_LINEAR)
  {
    // it's a min-max node; check the children  
    pChild = MMN_CHILDLIST(pMMN);
    pChildPrev = MMN_CHILDLIST(pMMN);
    while(pChild != NULL)
    {
      if (pChild->nType != DTREE_LINEAR)         // min/max child
      {
        // do amalgamations on the child first
        Amalgamate(pChild, aLF,
                   nDim, nOutput, bSmaller);
        // can we amalgamate child's child list with ours?
        if (pChild->nType == pMMN->nType ||           // child is same op as us, or
            MMN_CHILDLIST(pChild)->pNext == NULL)     // child has only one child (can't happen here!)
        {
          MINMAXNODE* pChildNext;
          MINMAXNODE* pListTemp = MMN_CHILDLIST(pMMN);  // save our child list head
          MINMAXNODE* pInsertList = MMN_CHILDLIST(pChild); // get child's child list 
          MMN_CHILDLIST(pMMN) = pInsertList;        // insert at head (child list has already been reduced)
          while (pInsertList->pNext != NULL)
            pInsertList = pInsertList->pNext;       // find end of inserted list
          pInsertList->pNext = pListTemp;           // append our child list
          if (pChild == pListTemp)                  // if child was at old child list head,
            pChildPrev = pInsertList;               //   fixup previous pointer
          // remove child
          pChildNext = pChild->pNext;               // save child's sibling pointer
          pChildPrev->pNext = pChildNext;           // fixup sibling list
          MMN_CHILDLIST(pChild) = NULL;             // mark child's child list as empty before destroy
          bSmaller = TRUE;
          DestroyMinMaxNode(pChild);                // destroy unlinked child
          pChild = pChildNext;  // next child to examine
        }
        else  // child's child list cannot be amalgamated
        {
          pChildPrev = pChild;                      // advance pointers
          pChild = pChild->pNext;                   // next child to examine
        }
      }
      else //child is a linear piece, move on
      {
        pChildPrev = pChild;                      // advance pointers
        pChild = pChild->pNext;                   // next child to examine
      }
    } // loop if we haven't reached the end of the (possibly amalgamated) child list
  } // end what is done if node is not a linear form
  
  // if we have only one child, then we become the child!
  while(pMMN->nType != DTREE_LINEAR && MMN_CHILDLIST(pMMN)->pNext == NULL)
  {   
    MINMAXNODE* pChild = MMN_CHILDLIST(pMMN);     // save pointer to (singular) child list
    MINMAXNODE* pNext = pMMN->pNext;              // save pointer to sibling
    memcpy(pMMN, pChild, sizeof(MINMAXNODE));     // copy child info into this
    pMMN->pNext = pNext;                          // restore sibling
    MMN_CHILDLIST(pChild) = NULL;                 // mark child's child list as empty before destroy
    DestroyMinMaxNode(pChild);
    bSmaller = TRUE; // can't happen as a result of amalgamation
  } 
}
