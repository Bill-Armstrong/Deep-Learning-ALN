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
// splitlfn.cpp

///////////////////////////////////////////////////////////////////////////////
//  File version info:
// 
//  $Archive: /ALN Development/libaln/src/splitlfn.cpp $
//  $Workfile: splitlfn.cpp $
//  $Revision: 10 $
//  $Date: 8/18/07 3:04p $
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

int ALNAPI SplitLFN(ALN* pALN, ALNNODE* pNode)
{
	if(fabs(LFN_SPLIT_T(pNode))> 0) // if this is true, the choice of max or min is clear
	{
		if (LFN_SPLIT_T(pNode) <0) // This is TRUE if the ALN surface tends to be below the function values far from the centroid on the piece
	  {
			return ALNAddLFNs(pALN, pNode, GF_MAX, 2, NULL);  // A max tends to move the surface up far from the centroid.
	  }
		else
		{
			return ALNAddLFNs(pALN, pNode, GF_MIN, 2, NULL);
		}
	}
	else // the choice of max or min split is not clear, and we have set LFN_SPLIT_T to 0 in findsplitlfn.cpp
	{
		ALNNODE* pParent = NODE_PARENT(pNode);
		if(pParent == NULL || MINMAX_ISMAX(pParent)) //this affine piece is the top node just choose MAX
		{
			return ALNAddLFNs(pALN, pNode, GF_MAX, 2, NULL); 
		}
		else
		{
			return ALNAddLFNs(pALN, pNode, GF_MIN, 2, NULL);
		}
	}
}
