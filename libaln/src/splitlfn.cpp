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
	// We only split if the piece doesn't fit within the noise limit
	// and the direction of split will likely not be close
	if (LFN_SPLIT_T(pNode) < 0) // This is TRUE if the ALN surface tends to be below
		// the training values far from the centroid.
	{
		return ALNAddLFNs(pALN, pNode, GF_MAX, 2, NULL);  // A max is convex down   \_/
	}
	else
	{
		return ALNAddLFNs(pALN, pNode, GF_MIN, 2, NULL); // A min is convex up      ^
	}
}
