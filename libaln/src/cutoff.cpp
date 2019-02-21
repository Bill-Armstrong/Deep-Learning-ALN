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

// cutoff.cpp


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
