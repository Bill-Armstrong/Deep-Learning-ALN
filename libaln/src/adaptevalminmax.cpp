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

// adaptevalminmax.cpp


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
// minmax node specific eval - evaluation and adaptation setup
//  - sets the distance, active and goal surfaces, and eval flag
// NOTE: cutoff always passed on stack!

double ALNAPI AdaptEvalMinMax(ALNNODE* pNode, ALN* pALN, const double* adblX, CEvalCutoff cutoff, ALNNODE** ppActiveLFN)
{
	ASSERT(NODE_ISMINMAX(pNode));

	// set node eval flags
 	NODE_FLAGS(pNode) |= NF_EVAL;  //NODE_FLAGS(pNode) ((pNode)->fNode)
	NODE_FLAGS(MINMAX_LEFT(pNode)) &= ~NF_EVAL;
	NODE_FLAGS(MINMAX_RIGHT(pNode)) &= ~NF_EVAL;//((pNode)->DATA.MINMAX.CHILDREN.CHILDSEPARATE.pRightChild)
	// set first child
	ALNNODE* pChild0;
	if (MINMAX_EVAL(pNode))    // ((pNode)->DATA.MINMAX.pEvalChild)
		pChild0 = MINMAX_EVAL(pNode);
	else
		pChild0 = MINMAX_LEFT(pNode);

	// set next child
	ALNNODE* pChild1;
	if (pChild0 == MINMAX_LEFT(pNode))
		pChild1 = MINMAX_RIGHT(pNode);
	else
		pChild1 = MINMAX_LEFT(pNode);

	// get reference to region for this node
	ALNREGION& region = pALN->aRegions[NODE_REGION(pNode)];
	/*
	if (region.dbl4SE > 0.0) // if smoothing is used
	{
		// loosen cutoff constraint for children
		if (MINMAX_ISMAX(pNode) && cutoff.bMax)
			cutoff.dblMax -= region.dbl4SE;
		else if (MINMAX_ISMIN(pNode) && cutoff.bMin)
			cutoff.dblMin += region.dbl4SE;
	} REMOVED THIS TO SEE WHAT HAPPENs
	*/
	// eval first child
	ALNNODE* pActiveLFN0;
	double dbl0 = AdaptEval(pChild0, pALN, adblX, cutoff, &pActiveLFN0);
	/*
	// see if we can cutoff...
	if (Cutoff(dbl0, pNode, cutoff, region.dbl4SE))
	{
		*ppActiveLFN = pActiveLFN0;
		MINMAX_ACTIVE(pNode) = pChild0;
		NODE_DISTANCE(pNode) = dbl0;
		MINMAX_RESPACTIVE(pNode) = 1.0;	 // we can't have < 1 without additional evaluation
		return dbl0;  
	}   Removed the cutoff to see what happens
	*/
	// eval second child
	ALNNODE* pActiveLFN1;
	double dbl1 = AdaptEval(pChild1, pALN, adblX, cutoff, &pActiveLFN1);

	if (region.dblSmoothEpsilon > 0.0) // if this is true, we are using smoothing
	{
		// calc active child, active child response, and distance
		int nActive = CalcActiveChild(MINMAX_RESPACTIVE(pNode),
			NODE_DISTANCE(pNode),
			dbl0, dbl1, pNode,
			region.dblSmoothEpsilon,
			region.dbl4SE, region.dblOV16SE);

		if (nActive == 0)
		{
			*ppActiveLFN = pActiveLFN0;
			MINMAX_ACTIVE(pNode) = pChild0;
			// distance and respactive already set in call to CalcActiveChild
		}
		else
		{
			*ppActiveLFN = pActiveLFN1;
			MINMAX_ACTIVE(pNode) = pChild1;
			// distance and respactive already set in call to CalcActiveChild
		}
	}
	else  // there is no smoothing and we can simplify without calling CalcActiveChild()
	{
		// Recall that dbl0 == dbl1 is not a rare event!  It always happens after a split,
		// however it happens then only once as the first adapt will likely destroy equality.
		MINMAX_RESPACTIVE(pNode) = 1.0;
		if ((MINMAX_ISMAX(pNode) > 0) == (dbl1 > dbl0)) // int MINMAX_ISMAX is used as a bit-vector!
		{
			NODE_DISTANCE(pNode) = dbl1;
			*ppActiveLFN = pActiveLFN1;
			MINMAX_ACTIVE(pNode) = pChild1;
		}
		else
		{
			NODE_DISTANCE(pNode) = dbl0;
			*ppActiveLFN = pActiveLFN0;
			MINMAX_ACTIVE(pNode) = pChild0;
		}
	}
	return NODE_DISTANCE(pNode);
}



