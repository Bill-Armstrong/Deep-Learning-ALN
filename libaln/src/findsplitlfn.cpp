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
static BOOL CanSplitLFN(ALN* pALN, ALNNODE* pNode) //routine
{
  // In order to split, the leaf must have at least certain number of adaptation hits during the epoch
	// such that the training error of the piece is made almost as small as possible without splitting.
  // That number of hits depends on the learning rate, so the earliest splitting should occur is about every
	// 1/(learning rate) epochs. The error remains the same after a split,
	// because the two pieces are duplicated and adjusted for fillets. After the slightest change,
	// the error should be divided between the two pieces resulting from the split
	// The total training error should then decrease after a split.
	if (LFN_CANSPLIT(pNode)&& (LFN_SPLIT_COUNT(pNode) > pALN->nDim)	&& (LFN_SPLIT_RESPTOTAL(pNode) > pALN->nDim ))
  {
		return TRUE;
	}
	return FALSE;
}

// The following is probably of no use.  It was removed above
//ALNCONSTRAINT* pConstr = GetVarConstraint(NODE_REGION(pNode), pALN, pALN->nOutput);
//ASSERT(pConstr->dblSqEpsilon == pConstr->dblEpsilon * pConstr->dblEpsilon);
// Now, if the square error is sufficiently large, the leaf is split (in splitlfn.cpp.
// Review Sept 30 2018  This is where the noise variance estimate should enter instead of the dblEpsilon.
//if(LFN_SPLIT_SQERR(pNode) > LFN_SPLIT_RESPTOTAL(pNode)* pConstr->dblSqEpsilon) 
//{
// LFN_SPLIT_T measures how much the ALN surface(smoothed)is currently above the data points far from the centroid.
// Its units are output units. If LFN_SPLIT_T is positive, the split results in a MIN node, otherwise a MAX.
// If the criterion LFN_SPLIT_T is below the level of noise, it is set to zero
// which will signal splitlfn.cpp to choose the same type, MAX or MIN, 
// as the parent node (if the present node is not the root).
//if(fabs(LFN_SPLIT_T(pNode)) < pConstr->dblEpsilon ) LFN_SPLIT_T(pNode) = 0; Forget this
//return TRUE;
//}


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
