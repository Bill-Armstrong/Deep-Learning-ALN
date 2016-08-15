// ALN Library
// Copyright (C) 1995 - 2010 William W. Armstrong.
// 
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
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
// findsplitlfn.cpp

///////////////////////////////////////////////////////////////////////////////
//  File version info:
// 
//  $Archive: /ALN Development/libaln/src/findsplitlfn.cpp $
//  $Workfile: findsplitlfn.cpp $
//  $Revision: 13 $
//  $Date: 8/18/07 4:27p $
//  $Author: Arms $
//
///////////////////////////////////////////////////////////////////////////////

#ifdef ALNDLL
#define ALNIMP __declspec(dllexport)
#endif

#include <aln.h>
#include "alnpriv.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

static BOOL CanSplitLFN(ALN* pALN, ALNNODE* pNode);
static ALNNODE* ALNAPI CompareSplitLFNs(ALNNODE* pNode1, ALNNODE* pNode2);
static void ALNAPI DoFindSplitLFN(ALN* pALN, ALNNODE* pNode, 
                                  ALNNODE*& pSplitLFN);

ALNNODE* ALNAPI FindSplitLFN(ALN* pALN)
{
  ASSERT(pALN);
  ASSERT(pALN->pTree);

  ALNNODE* pSplitLFN = NULL;
  DoFindSplitLFN(pALN, pALN->pTree, pSplitLFN);

  ASSERT(pSplitLFN == NULL || LFN_CANSPLIT(pSplitLFN));
  
  return pSplitLFN;
}

inline 
static BOOL CanSplitLFN(ALN* pALN, ALNNODE* pNode)
{
  // Code reviewed by WWA March 10, 2015.  The definition of LFN_SPLIT_T has completely changed in the last few days.
  // In order to split, the leaf must have at least certain number of hits during the epoch.
  // The number depends on the learning rate, so splitting should occur about every
	// 1/(learning rate) epochs. Since the error of a piece is not changed by a split
	// because the two pieces are duplicated and adjusted for fillets, the error should decrease after a split.
  // There have to be enough total adaptations ( taking account of responsibility)
	// to equal the dimension of the problem multiplied by 4 to have pieces be determined after the split.
	// A higher factor, e.g. 4 allows the split to be unequal and still have both parts determined,
	if (LFN_CANSPLIT(pNode)&& (LFN_SPLIT_COUNT(pNode) > 4.0 * pALN->nDim)
		&& (LFN_SPLIT_RESPTOTAL(pNode) > 4.0 * pALN->nDim ))  // more experimentation would be good
  {
    ALNCONSTRAINT* pConstr = GetVarConstraint(NODE_REGION(pNode), pALN, pALN->nOutput);
		ASSERT(pConstr->dblSqEpsilon == pConstr->dblEpsilon * pConstr->dblEpsilon);
		// Now, if the square error is sufficiently large, the leaf is split (in splitlfn.cpp.
		if(LFN_SPLIT_SQERR(pNode) > LFN_SPLIT_RESPTOTAL(pNode)* pConstr->dblSqEpsilon) 
		{
			// LFN_SPLIT_T measures how much the ALN surface(smoothed)is currently above the data points far from the centroid.
			// Its units are output units. If LFN_SPLIT_T is positive, the split results in a MIN node, otherwise a MAX.
			// If the criterion LFN_SPLIT_T is below the level of noise, it is set to zero
			// which will signal splitlfn.cpp to choose the same type, MAX or MIN, 
			// as the parent node (if the present node is not the root).
			if(fabs(LFN_SPLIT_T(pNode)) < pConstr->dblEpsilon ) LFN_SPLIT_T(pNode) = 0; 
			return TRUE;
		}
	} 
  return FALSE;
}

static void ALNAPI DoFindSplitLFN(ALN* pALN, ALNNODE* pNode,
                                  ALNNODE*& pSplitLFN)
{
  ASSERT(pNode);
  if (NODE_ISMINMAX(pNode))
  {
    DoFindSplitLFN(pALN, MINMAX_LEFT(pNode), pSplitLFN);
    DoFindSplitLFN(pALN, MINMAX_RIGHT(pNode), pSplitLFN);
  }
  else
  {
    ASSERT(NODE_ISLFN(pNode));
    // check if splittable
    if (CanSplitLFN(pALN, pNode))
      SplitLFN(pALN, pNode);
  }
}
