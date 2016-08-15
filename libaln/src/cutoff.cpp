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

// cutoff.cpp

///////////////////////////////////////////////////////////////////////////////
//  File version info:
// 
//  $Archive: /ALN Development/libaln/src/cutoff.cpp $
//  $Workfile: cutoff.cpp $
//  $Revision: 6 $
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

BOOL ALNAPI Cutoff(double dbl, const ALNNODE* pNode, CEvalCutoff& cutoff, 
                   double dbl4SE)
{
  if(MINMAX_ISMAX(pNode))  // if pNode is a MAX
	{   
    // cutoff if we're greater than or equal to existing min
		if (cutoff.bMin && (dbl >= cutoff.dblMin))
		{ 
			return TRUE;  // cutoff!
		}

    // no cutoff... set new max

    // assume max adjusted by node prior to evaluation:
		// if (cutoff.bMax) 
		//   cutoff.dblMax -= dbl4SE;

		// I think all of the dbl4SE entries are wrong, even when there is smoothing. WWA
		/*
		if (!cutoff.bMax)
		{
		cutoff.bMax = TRUE;
		cutoff.dblMax = dbl - dbl4SE;
		}
		else if (dbl > cutoff.dblMax + dbl4SE)
		{
		cutoff.dblMax = dbl - dbl4SE;
		}
		else if(dbl > cutoff.dblMax + 0.75 * dbl4SE)
		{
		cutoff.dblMax = -dbl4SE +
		4 * (cutoff.dblMax + dbl4SE) -
		3 * dbl;
		}
		else if(dbl > cutoff.dblMax)
		{
		cutoff.dblMax = -0.3333333 * dbl +
		1.3333333 * cutoff.dblMax +
		dbl4SE;
		}
		else
		{
		cutoff.dblMax += dbl4SE;
		}
		}
		else  // pNode is a MIN
		{
		ASSERT(MINMAX_ISMIN(pNode));

		// cutoff if we're less than or equal to existing max
		if (cutoff.bMax && (dbl <= cutoff.dblMax))
		{
		return TRUE;
		}

		// no cutoff... set new min

		// assume min adjusted by node prior to evaluation:
		// if(cutoff.bMin)
		//   cutoff.dblMin += dbl4SE;

		if(!cutoff.bMin)
		{
		cutoff.bMin = TRUE;
		cutoff.dblMin = dbl + dbl4SE;
		}
		else if(dbl < cutoff.dblMin - dbl4SE)
		{
		cutoff.dblMin = dbl + dbl4SE;
		}
		else if(dbl < cutoff.dblMin - 0.75 * dbl4SE)
		{
		cutoff.dblMin = dbl4SE +
		4 * (cutoff.dblMin - dbl4SE) -
		3 * dbl;
		}
		else if (dbl < cutoff.dblMin)
		{
		cutoff.dblMin = -0.3333333 * dbl +
		1.3333333 * cutoff.dblMin -
		dbl4SE;
		}
		else
		{
		cutoff.dblMin -= dbl4SE;
		}
		}
	*/
		if (!cutoff.bMax)
		{
			cutoff.bMax = TRUE;
			cutoff.dblMax = dbl;
		}
		else if (dbl > cutoff.dblMax)  // we set a higher lower bound on the value of the current MAX node
		{
			cutoff.dblMax = dbl;
		}
	}
	else  // pNode is a MIN
	{
		ASSERT(MINMAX_ISMIN(pNode));

		// cutoff if we're less than or equal to existing max
		if (cutoff.bMax && (dbl <= cutoff.dblMax))
		{
			return TRUE;
		}

		// no cutoff... set new min

		// assume min adjusted by node prior to evaluation:
		if (!cutoff.bMin)
		{
			cutoff.bMin = TRUE;
			cutoff.dblMin = dbl;
		}
		else if (dbl < cutoff.dblMin)  // we set a lower upper bound on the value of the current MIN node
		{
			cutoff.dblMin = dbl;
		}
	}
  return FALSE; // no cutoff
}
