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

// countlfns.cpp

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
// count number of LFNs in an ALN

void ALNAPI CountLFNs(const ALNNODE* pNode, int& nTotal, int& nAdapted)
{
	ASSERT(pNode != NULL);
	
  int nTotalThis = 0;

	// are we an LFN?
	if (NODE_ISLFN(pNode))
  {
    nTotal++;
		nAdapted += !(NODE_RESPCOUNT(pNode) == 0);
  }
  else
  {
	  // we're a minmax... iterate over children
	  ASSERT(NODE_ISMINMAX(pNode));
  
    CountLFNs(MINMAX_LEFT(pNode), nTotal, nAdapted);
    CountLFNs(MINMAX_RIGHT(pNode), nTotal, nAdapted);
  }
}

