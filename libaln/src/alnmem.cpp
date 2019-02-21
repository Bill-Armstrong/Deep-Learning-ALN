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

// alnmem.cpp

#ifdef ALNDLL
#define ALNIMP __declspec(dllexport)
#endif

#include <aln.h>
#include "alnpriv.h"
#include <float.h>

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

// creates a new ALN
// nDim variables, 
// nOutput is default output variable
// ALN contains one region (index 0) with nDim constraints and one tree node 
// that is an LFN, constraints have min,max=[-DBL_MAX, DBL_MAX], epsilon=0.0001, 
// wmin,wmax=[-1000000, 1000000]
// ... returns pointer to ALN or NULL on failure
ALNIMP ALN* ALNAPI ALNCreateALN(int nDim, int nOutput)
{
  if (nDim < 2) return NULL;
  if (nOutput < 0 || nOutput >= nDim) return NULL;

  ALN* pALN = (ALN*)malloc(sizeof(ALN));
  if (pALN == NULL) return NULL;

  memset(pALN, 0, sizeof(ALN));
  pALN->nVersion = ALNVER;
  pALN->nDim = nDim;
  pALN->nOutput = nOutput;
  // allocate first region
  pALN->aRegions = (ALNREGION*)malloc(sizeof(ALNREGION));
  if (pALN->aRegions == NULL)
  {
    ALNDestroyALN(pALN);
    return NULL;
  }
  memset(pALN->aRegions, 0, sizeof(ALNREGION));
  pALN->nRegions = 1;

  // region parent
  pALN->aRegions->nParentRegion = -1;  // no parent

  // default region learn factor
  pALN->aRegions->dblLearnFactor = 1.0;
  
  // region smoothing epsilon - disabled
  pALN->aRegions->dblSmoothEpsilon = 0.0;
  
  // allocate and init region var constraints
  pALN->aRegions->aConstr = (ALNCONSTRAINT*)malloc(nDim * sizeof(ALNCONSTRAINT));
  if (pALN->aRegions->aConstr == NULL)
  {
    ALNDestroyALN(pALN);
    return NULL;
  }
  pALN->aRegions->nConstr = nDim;
  for (int i = 0; i < nDim; i++)
  {
    pALN->aRegions->aConstr[i].nVarIndex = i;
    pALN->aRegions->aConstr[i].dblMin = -DBL_MAX;
    pALN->aRegions->aConstr[i].dblMax = DBL_MAX;
    pALN->aRegions->aConstr[i].dblEpsilon = 0.0001;
    
    if (i != nOutput)
    {
      pALN->aRegions->aConstr[i].dblWMin = -1000000.0;
      pALN->aRegions->aConstr[i].dblWMax = 1000000.0;
    }
    else
    {
      pALN->aRegions->aConstr[i].dblWMin = -1.0;
      pALN->aRegions->aConstr[i].dblWMax = -1.0;
    }
  }

  // allocate and init tree node
  pALN->pTree = (ALNNODE*)malloc(sizeof(ALNNODE));
  if (pALN->pTree == NULL)
  {
    ALNDestroyALN(pALN);
    return NULL;
  }  
  memset(pALN->pTree, 0, sizeof(ALNNODE));
  pALN->pTree->pParent = NULL;
  pALN->pTree->fNode |= NF_LFN;
  pALN->pTree->nParentRegion = 0;
  LFN_SPLIT(pALN->pTree) = NULL;
  LFN_VDIM(pALN->pTree) = nDim;
  LFN_W(pALN->pTree) = (double*)malloc((nDim + 1) * sizeof(double));
  LFN_C(pALN->pTree) = (double*)malloc(nDim * sizeof(double));
  LFN_D(pALN->pTree) = (double*)malloc(nDim * sizeof(double));
  if (LFN_W(pALN->pTree) == NULL || LFN_C(pALN->pTree) == NULL || LFN_D(pALN->pTree) == NULL)
  {
    ALNDestroyALN(pALN);
    return NULL;
  }
  memset(LFN_W(pALN->pTree), 0, (nDim + 1) * sizeof(double));
  memset(LFN_C(pALN->pTree), 0, nDim * sizeof(double));
  memset(LFN_D(pALN->pTree), 0, nDim * sizeof(double));
  return pALN;
}

// helper: destroys tree node
int ALNAPI DestroyTree(ALNNODE* pTree)
{
  if (pTree == NULL)
    return 0;

  if (pTree->fNode & NF_LFN)
  {
    if (LFN_VARMAP(pTree) != NULL)
      free(LFN_VARMAP(pTree));

    if (LFN_SPLIT(pTree) != NULL)
      free(LFN_SPLIT(pTree));

    if (LFN_W(pTree) != NULL)
      free(LFN_W(pTree));

    if (LFN_C(pTree) != NULL)
      free(LFN_C(pTree));

    if (LFN_D(pTree) != NULL)
      free(LFN_D(pTree));

	}
  else
  {
    // destroy children 
    ASSERT(pTree->fNode & NF_MINMAX);

    int nChildren = MINMAX_NUMCHILDREN(pTree);
    for (int i = 0; i < nChildren; i++)
    {
      DestroyTree(MINMAX_CHILDREN(pTree)[i]);
      MINMAX_CHILDREN(pTree)[i] = NULL;
    }
  }

  // free node memory
  free(pTree);

  return 1;
}

// destroys an ALN
//   ... returns 0 on failure, non-zero on success
ALNIMP int ALNAPI ALNDestroyALN(ALN* pALN)
{
  if (pALN == NULL)
    return 0;

  // regions and constraints
  for (int i = 0; i < pALN->nRegions; i++)
  {
    if (pALN->aRegions[i].nConstr > 0)
    {
      free(pALN->aRegions[i].aConstr);
      if (pALN->aRegions[i].afVarMap)
        free(pALN->aRegions[i].afVarMap);
    }
  }
  if (pALN->nRegions > 0)
    free(pALN->aRegions);

  // tree
  if (pALN->pTree)
    DestroyTree(pALN->pTree);

  // ALN
  free(pALN);

  return 1;
}

#ifdef ENABLE_REGIONS

// adding a new region to the ALN
//   - pALN
//   - nParentRegion >= 0
//   - nConstr = number of constraints, or -1 if all vars constrained
//   - anConstr = array of constrained var indexes, ignored if nConstr == -1
//   constraints inherit their attributes from parent region
//   ... returns index of new region in ALN, -1 on failure
//   
ALNIMP int ALNAPI ALNAddRegion(ALN* pALN, int nParentRegion, 
                               double dblLearnFactor, 
                               int nConstr, int* anConstr)
{
  if (pALN == NULL)
    return -1;

  if (nParentRegion < 0 || nParentRegion >= pALN->nRegions)
    return -1;

  if (nConstr < 0 || nConstr > pALN->nDim)
    return -1;

  if ((nConstr > 0 && nConstr < pALN->nDim) && anConstr == NULL)
    return -1;

  if (dblLearnFactor < 0)
    return -1;
  
  // allocate var map?
  char* afVarMap = (char*)((nConstr == 0 || nConstr == pALN->nDim) ? NULL :
                      malloc(MAPBYTECOUNT(pALN->nDim)));
		// failure tolerable
	if (afVarMap != NULL)
  {
    memset(afVarMap, 0, MAPBYTECOUNT(pALN->nDim));
    for (int i = 0; i < nConstr; i++)
    {
      if (TESTMAP(afVarMap, anConstr[i]))
      {
        // multiple occurences of same var index!
        free(afVarMap);
        return -1;
      }

      SETMAP(afVarMap, anConstr[i]);
    }
  }

  // allocate constraints
  CONSTRAINT* aConstr = NULL;
  if (nConstr > 0)
  {
    aConstr = (CONSTRAINT*)malloc(nConstr * sizeof(CONSTRAINT));
    if (aConstr == NULL)
    {
      if (afVarMap != NULL) free(afVarMap);
      return -1;
    }

    // inherit parent constraints
    for (int i = 0, j = 0; i < pALN->nDim; i++)
    {
      if (afVarMap && !TESTMAP(afVarMap, i))
        continue;

			if (afVarMap == NULL && nConstr < pALN->nDim)
			{
				// check constraint array
				for (int k = 0; k < nConstr; k++)
				{
					if (anConstr[k] == i)
						break;
				}
				if (k == nConstr)
					continue;	// constraint for var i not found
			}

			// var i is constrained in this region, copy from parent
      CONSTRAINT* pConstr = GetVarConstraint(nParentRegion, pALN, i);
      ASSERT(pConstr && pConstr->nVarIndex == i);
      ASSERT(nConstr == 0 || j < nConstr);
      memcpy(&(aConstr[j++]), pConstr, sizeof(CONSTRAINT));
    }
  }

  // allocate new space for region
  ALNREGION* aRegions = (ALNREGION*)realloc(pALN->aRegions, 
                           (pALN->nRegions + 1) * sizeof(ALNREGION));
  if (aRegions == NULL)
  {
    if (aConstr != NULL) free(aConstr);
    if (afVarMap != NULL) free(afVarMap);
    return -1;
  }
  pALN->aRegions = aRegions;
  int nNewRegion = pALN->nRegions++;  // inc region count

  // assign region members
  aRegions[nNewRegion].nParentRegion = nParentRegion;
  aRegions[nNewRegion].dblLearnFactor = dblLearnFactor;
  aRegions[nNewRegion].nConstr = nConstr;
  aRegions[nNewRegion].aConstr = aConstr;
  aRegions[nNewRegion].afVarMap = afVarMap;

  // smoothing memebrs initially uncalculated (will be calc'ed in PrepALN)
  aRegions[nNewRegion].dblSmoothEpsilon = 0;   // TODO: this should be a param!
  aRegions[nNewRegion].dbl4SE = 0;                    
  aRegions[nNewRegion].dblOV16SE = 0;                 

  return nNewRegion;
}

#endif /* ENABLE_REGIONS */

// adding LFNs to a tree
//   - pALN
//   - pParent = parent node, must be an LFN 
//   - pParentMinMaxType = one of GF_MIN, GF_MAX,
//   - nLFNs = number of LFNs to add
//   - apLFNs = array of pointers to all LFN children created
// LFN parent automatically converted to minmax (and vectors freed)
// all new children are LFNs... vectors are automatically allocated, child
// parent regions are the same as the parent node
ALNIMP int ALNAPI ALNAddLFNs(ALN* pALN, ALNNODE* pParent, 
                             int nParentMinMaxType, int nLFNs,
                             ALNNODE** apLFNs)
{
  if (nLFNs < 2)
    return ALN_GENERIC;

  if (!pALN)
    return ALN_GENERIC;

  if (!pParent)
    return ALN_GENERIC;

  if (pParent->nParentRegion < 0)
    return ALN_GENERIC;  // invalid parent region index?

  if (NODE_ISMINMAX(pParent))
    return ALN_GENERIC;  // can't add children to existing minmax

  if (nParentMinMaxType != GF_MIN && nParentMinMaxType != GF_MAX)
    return ALN_GENERIC;  // must specify minmax type 

  if (nParentMinMaxType == (GF_MIN | GF_MAX))
    return ALN_GENERIC;  // can't specify both!

  // are we splitting?
  ASSERT(NODE_ISLFN(pParent));
  BOOL bSplit = LFN_CANSPLIT(pParent);

  // add children
  ALNNODE* apChildren[2];
  apChildren[0] = apChildren[1] = NULL;

  // allocate new children for this node and convert to LFNs
  for (int i = 0; i < 2; i++)
  {
    // alloc new child
    ALNNODE*& pChild = apChildren[i];
    try
    {
      pChild = (ALNNODE*)malloc(sizeof(ALNNODE));
      if (pChild == NULL) ThrowALNMemoryException();
    
      memset(pChild, 0, sizeof(ALNNODE));    
      pChild->pParent = pParent;
      pChild->fNode |= NF_LFN;
    
      pChild->nParentRegion = pParent->nParentRegion;
      pChild->nRespCount = 0;
      pChild->nRespCountLastEpoch = 0;
      pChild->dblDistance = 0;

      // allocate vectors
      LFN_SPLIT(pChild) = NULL;
      LFN_VARMAP(pChild) = NULL;
      LFN_VDIM(pChild) = pALN->nDim;
      LFN_W(pChild) = (double*)malloc((pALN->nDim + 1) * sizeof(double));
      LFN_C(pChild) = (double*)malloc(pALN->nDim * sizeof(double));
      LFN_D(pChild) = (double*)malloc(pALN->nDim * sizeof(double));
      if (LFN_W(pChild) == NULL || LFN_C(pChild) == NULL ||  LFN_D(pChild) == NULL )
      {
        ThrowALNMemoryException();
      }

      // set split
      if (bSplit) 
      {
        LFN_SPLIT(pChild) = (ALNLFNSPLIT*)malloc(sizeof(ALNLFNSPLIT));
        if (LFN_SPLIT(pChild) == NULL)
          ThrowALNMemoryException();
		
        pChild->fNode |= LF_SPLIT;
        LFN_SPLIT_COUNT(pChild) = 0;
        LFN_SPLIT_SQERR(pChild) = 0.0;
        LFN_SPLIT_RESPTOTAL(pChild) = 0.0;

        // copy parent vectors
        memcpy(LFN_W(pChild), LFN_W(pParent), (pALN->nDim + 1) * sizeof(double));
        memcpy(LFN_C(pChild), LFN_C(pParent), pALN->nDim * sizeof(double));
        memcpy(LFN_D(pChild), LFN_D(pParent), pALN->nDim * sizeof(double));
        // shift LFN up or down depending on parent minmax type
				// the shifts are different but close so the two children differentiate
				// and the combined effect is not to change the value of the single LFN
        double dblSE = pALN->aRegions[pChild->nParentRegion].dblSmoothEpsilon;
				int nOutput = pALN->nOutput;
				double dblChange;
				// 2009.11.19  This change has to be tiny because it affects the fillets!
				//dblChange = (0.9343727 + i * 0.1334818) * dblSE; old values where did they come from?
				// When a piece splits into two equal leaf nodes, there is a fillet inserted so both pieces have to move
				// in the output direction to leave the function unchanged
				dblChange = (1.0 + i * 0.0001) * dblSE;
        if (nParentMinMaxType == GF_MIN)
				{
          *LFN_W(pChild) += dblChange;
					LFN_C(pChild)[nOutput] += dblChange;
				}
        else
				{
          *LFN_W(pChild) -= dblChange;
					LFN_C(pChild)[nOutput] -= dblChange;
				}

        // set LFN initialized
        NODE_FLAGS(pChild) |= LF_INIT;
      }
      else
      {
        LFN_SPLIT(pChild) = NULL;

        // zero vectors
        memset(LFN_W(pChild), 0, (pALN->nDim + 1) * sizeof(double));
        memset(LFN_C(pChild), 0, pALN->nDim * sizeof(double));
        memset(LFN_D(pChild), 0, pALN->nDim * sizeof(double));
      }

      
    }
    catch(CALNMemoryException* e)
    {
      e->Delete();

      // free any successful new child allocations
      for (i = 0; i < 2; i++)
      {
        pChild = apChildren[i];

        if (pChild)
        {
          ASSERT(pChild->fNode & NF_LFN);
          if (LFN_VARMAP(pChild)) free(LFN_VARMAP(pChild));
          if (LFN_SPLIT(pChild)) free(LFN_SPLIT(pChild));
          if (LFN_W(pChild)) free(LFN_W(pChild));
          if (LFN_C(pChild)) free(LFN_C(pChild));
          if (LFN_D(pChild)) free(LFN_D(pChild));
          free(pChild);
          pChild = NULL;
        }
      }
      
      return NULL;
    }
  }

  // pParent is unmodified at this point
  // no further memory allocations are required, so it is safe to convert
  // LFN to a minmax without any errors or exceptions

  ASSERT(NODE_ISLFN(pParent));
  
  // free existing vectors
  if (LFN_VARMAP(pParent)) free(LFN_VARMAP(pParent));
  if (LFN_SPLIT(pParent)) free(LFN_SPLIT(pParent));
  if (LFN_W(pParent)) free(LFN_W(pParent));
  if (LFN_C(pParent)) free(LFN_C(pParent));
  if (LFN_D(pParent)) free(LFN_D(pParent));
 
  // convert node type
  pParent->fNode &= ~NF_LFN;
  pParent->fNode |= NF_MINMAX | nParentMinMaxType;

  // NOTE: we leave any split flags present in converted LFN so that we may
  // trace the effects of any split algorithms

  // assign children
  memset(MINMAX_CHILDREN(pParent), 0, 2 * sizeof(ALNNODE*));
  MINMAX_LEFT(pParent) = apChildren[0];
  MINMAX_RIGHT(pParent) = apChildren[1];
  
  MINMAX_EVAL(pParent) = NULL;
  MINMAX_ACTIVE(pParent) = NULL;
  MINMAX_GOAL(pParent) = NULL;
  
  ASSERT(NODE_ISMINMAX(pParent) && MINMAX_TYPE(pParent) == nParentMinMaxType);

  int nResult = ALN_NOERROR;
  
  if (nLFNs > 2)
  {
    // add left child to array of LFNs 
    if (apLFNs != NULL)
    {
      // add at the end of the array, so recursive calls with fewer children
      // have empty slots near the beginning of the array
      apLFNs[nLFNs - 1] = apChildren[0];
    }

    // recurse, adding more LFNs to right child
    nResult = ALNAddLFNs(pALN, apChildren[1], nParentMinMaxType, 
                         nLFNs - 1, apLFNs);
  }
  else
  {
    // add both children to array of LFNs
    if (apLFNs != NULL)
    {
      apLFNs[0] = apChildren[0];
      apLFNs[1] = apChildren[1];
    }
  }

  return nResult;
}

// adding multiple layers tree to a tree
//   - pALN
//   - pParent = parent node, must be an LFN 
//   - pParentMinMaxType = one of GF_MIN, GF_MAX
//   - nLayers = maximum number of layers to add, the existing node pParent
//               counts as layer 1
//   - nMinMaxFanin = number of logic minmaxs that fanin to a parent, must be >= 2
//   - nLFNFanin = number of LFNs that fanin to a parent, must be >= 2 
//   - nFlags = MULTILAYER_FULL or MULTILAYER_RAGGED
//              MULTILAYER_RAGGED has nLFNFanin + 1 children at (nLayer - 1)
//              
//   ... returns pointer to child array on success, NULL if node already has
//       children or other failure
// LFN nodes automatically converted to minmax if necessary (and vectors freed)
// child parent regions are the same as the parent node

ALNIMP int ALNAPI ALNAddMultiLayer(ALN* pALN, ALNNODE* pParent, 
                                   int nParentMinMaxType, int nLayers,
                                   int nMinMaxFanin, int nLFNFanin,
                                   int nFlags)
{
  if (nMinMaxFanin < 2 || nLFNFanin < 2)
    return ALN_GENERIC;

  if (!pALN)
    return ALN_GENERIC;

  if (!pParent)
    return ALN_GENERIC;

  if (pParent->nParentRegion < 0)
    return ALN_GENERIC;  // invalid parent region index?

  if (NODE_ISMINMAX(pParent))
    return ALN_GENERIC;  // can't add children to existing minmax

  if (nParentMinMaxType != GF_MIN && nParentMinMaxType != GF_MAX)
    return ALN_GENERIC;  // must specify minmax type 

  if (nParentMinMaxType == (GF_MIN | GF_MAX))
    return ALN_GENERIC;  // can't specify both!

  if (nFlags != MULTILAYER_FULL && nFlags != MULTILAYER_RAGGED)
    return ALN_GENERIC;  // unknown flags

  if (nLayers <= 1)
    return ALN_NOERROR;  // nothing to do!
  
  // determine fanin at this level
  int nFanin = nMinMaxFanin;  // by default
  if (nLayers == 2) // next layer will be LFNs
  {
    nFanin = nLFNFanin;
  }
  else if (nLayers == 3 && nFlags == MULTILAYER_RAGGED)
  {
    // ragged layer... we'll have LFNs fanning in, plus one minmax
    nFanin = nLFNFanin + 1;
  }

  // space for LFNs
  ALNNODE** apLFNs = NULL;

  // allocate space for child LFN array
  apLFNs = (ALNNODE**)malloc(nFanin * sizeof(ALNNODE*));
  if (apLFNs == NULL)
    return ALN_OUTOFMEM;

  // current result
  int nResult = ALN_NOERROR;

  // add the children
  nResult = ALNAddLFNs(pALN, pParent, nParentMinMaxType, nFanin, apLFNs);
  if (nResult != ALN_NOERROR)
  {
    free(apLFNs);  // clean up
    return nResult;
  }

  // set child minmax type
  int nChildMinMaxType = GF_MAX;   // assume child is OR...
  if (nParentMinMaxType & GF_MAX)  // unless parent is an OR...
    nChildMinMaxType = GF_MIN;    // in which case child is an AND

  // recurse over child LFNs...
  for (int i = 0; i < nFanin; i++)
  {
    // calc child depth
    int nChildLayers = nLayers - 1;
    
    // check for special ragged tree conditions
    if ((nChildLayers == 2) && (i > 0) && (nFlags == MULTILAYER_RAGGED))
    {
      // only first child of ragged tree goes completely to two layers... 
      // all others are LFNs, and stay at one layer
      break;
    }

    if (nChildLayers > 1)
    { 
      ALNNODE* pChild = apLFNs[i];

      nResult = ALNAddMultiLayer(pALN, pChild, nChildMinMaxType,
                                 nLayers - 1, nMinMaxFanin, nLFNFanin, 
                                 nFlags);

      if (nResult != ALN_NOERROR)
      {
        // unsuccessful... destroy allocated children
        for (int j = 0; j < nFanin; j++)
        {
          pChild = apLFNs[j];

          DestroyTree(pChild);
          pChild = NULL;
        }

        break;
        // for (i...) loop will terminate, 
        // function will exit normally at bottom with desired return code
      }
    }
  }
  
  // clean up
  free(apLFNs);

  return nResult;
}


// make a growable subtree
//   - pALN
//   - pParent = parent node 
//   ... returns 0 on failure, non-zero on success

ALNIMP int ALNAPI ALNSetGrowable(ALN* pALN, ALNNODE* pParent)
{
  if (!pALN)
    return 0;

  if (!pParent)
    return 0;

  if (pParent->nParentRegion < 0)
    return 0;

  if ((NODE_MINMAXTYPE(pParent) & NF_LFN) == 0)
    return 0;  // needs to be an LFN!

  if (LFN_CANSPLIT(pParent))
  {
    ASSERT(LFN_SPLIT(pParent) != NULL);
    return 1;
  }

  LFN_SPLIT(pParent) = (ALNLFNSPLIT*)malloc(sizeof(ALNLFNSPLIT));
  if (LFN_SPLIT(pParent) == NULL)
  {
    return 0;
  }

  pParent->fNode |= LF_SPLIT;
  LFN_SPLIT_COUNT(pParent) = 0;
  LFN_SPLIT_SQERR(pParent) = 0.0;
  LFN_SPLIT_RESPTOTAL(pParent) = 0.0;

  return 1;
}

