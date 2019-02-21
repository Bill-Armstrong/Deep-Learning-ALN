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

// cutoffevalminmax.cpp

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
// minmax node specific eval - returns distance to surface
//  - non-destructive, ie, does not change ALN structure
// NOTE: cutoff always passed on stack!

double ALNAPI CutoffEvalMinMax(const ALNNODE* pNode, const ALN* pALN,
	const double* adblX, CEvalCutoff cutoff,
	ALNNODE** ppActiveLFN)
{
	ASSERT(NODE_ISMINMAX(pNode));

	// set first child
	const ALNNODE* pChild0;
	if (MINMAX_EVAL(pNode))
		pChild0 = MINMAX_EVAL(pNode);
	else
		pChild0 = MINMAX_LEFT(pNode);

	// set next child
	const ALNNODE* pChild1;
	if (pChild0 == MINMAX_LEFT(pNode))
		pChild1 = MINMAX_RIGHT(pNode);
	else
		pChild1 = MINMAX_LEFT(pNode);

	// get reference to region for this node
	const ALNREGION& region = pALN->aRegions[NODE_REGION(pNode)];
	double dblDist, dblRespActive;
	if (region.dbl4SE > 0.0) // smoothing is used
	{
		// I think this step is wrong.  Smoothing makes a max function greater and a min function smaller.
		// So I am commenting it out to see what happens.  WWA
		/*
		// loosen cutoff constraint for children
		if (MINMAX_ISMAX(pNode) && cutoff.bMax)
			cutoff.dblMax -= region.dbl4SE;
		else if (MINMAX_ISMIN(pNode) && cutoff.bMin)
			cutoff.dblMin += region.dbl4SE;
		*/
		// eval first child
		ALNNODE* pActiveLFN0;
		double dbl0 = CutoffEval(pChild0, pALN, adblX, cutoff, &pActiveLFN0);

		// see if we can cutoff...
		if (Cutoff(dbl0, pNode, cutoff, region.dbl4SE))
		{
			*ppActiveLFN = pActiveLFN0;
			return dbl0;
		}

		// eval second child
		ALNNODE* pActiveLFN1;
		double dbl1 = CutoffEval(pChild1, pALN, adblX, cutoff, &pActiveLFN1);

		// calc active child, active child response, and distance
		int nActive = CalcActiveChild(dblRespActive,
			dblDist,
			dbl0, dbl1, pNode,
			region.dblSmoothEpsilon,
			region.dbl4SE, region.dblOV16SE);

		if (nActive == 0)
		{
			*ppActiveLFN = pActiveLFN0;
		}
		else
		{
			*ppActiveLFN = pActiveLFN1;
		}
	}
	else  //  dbl4SE == 0, i.e. there is zero smoothing
	{
		// eval first child
		ALNNODE* pActiveLFN0;
		double dbl0 = CutoffEval(pChild0, pALN, adblX, cutoff, &pActiveLFN0);

		// see if we can cutoff...
		if (Cutoff(dbl0, pNode, cutoff, 0.0))
		{
			*ppActiveLFN = pActiveLFN0;
			return dbl0;
		}

		// eval second child
		ALNNODE* pActiveLFN1;
		double dbl1 = CutoffEval(pChild1, pALN, adblX, cutoff, &pActiveLFN1);

		// calc active child and distance without using CalcActiveChild()
		if((MINMAX_ISMAX(pNode) > 0) == (dbl1 > dbl0)) // int MINMAX_ISMAX is used as a bit-vector!
		{
			*ppActiveLFN = pActiveLFN1;
			dblDist = dbl1;
		}
		else
		{
			*ppActiveLFN = pActiveLFN0;
			dblDist = dbl0;
		}

	}
  
  return dblDist;
}

