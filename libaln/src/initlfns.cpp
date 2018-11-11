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

// initlfns.cpp

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
// init any uninitialized LFN's

void ALNAPI InitLFNs(ALNNODE* pNode, ALN* pALN, const double* adblX)
{
	ASSERT(pNode != NULL);
	ASSERT(pALN != NULL);
	ASSERT(adblX != NULL);
	
	// are we an LFN?
	if (NODE_ISLFN(pNode))
	{
    if (!LFN_ISINIT(pNode)) // not initialized?
    {
  		int nOutput = pALN->nOutput;
  		int nDim = LFN_VDIM(pNode);
      ASSERT(nDim == pALN->nDim); 
  		double* adblW = LFN_W(pNode) + 1; // weight vector... skip bias
  		double* adblC = LFN_C(pNode);			// centroid vector
  		double* adblD = LFN_D(pNode);			// ave sq dist from centroid vector
			// vector initialization
  	  for (int i = 0; i < nDim; i++)
  	  { 
  			ALNCONSTRAINT* pConstr = GetVarConstraint(NODE_REGION(pNode), pALN, i);
  			ASSERT(pConstr != NULL);
			
        // init weights			
  			adblW[i] = max(min(pConstr->dblWMax, ALNRandFloat() * 0.0002 - 0.0001),
  			               pConstr->dblWMin);
  	    adblC[i] = adblX[i];
  	    adblD[i] = pConstr->dblSqEpsilon;
  	  }
    
  	  // the hyperplane is on the centroid to start so W[0] = 0
  	  LFN_W(pNode)[0] = 0;
    
  		ASSERT(adblW[nOutput] == -1.0);// make sure output var is -1

  	  // successfully initialized
      LFN_FLAGS(pNode) |= LF_INIT;
    }
	}
  else
  {
	  // we're a minmax... iterate over children
	  ASSERT(NODE_ISMINMAX(pNode));

    InitLFNs(MINMAX_LEFT(pNode), pALN, adblX);
    InitLFNs(MINMAX_RIGHT(pNode), pALN, adblX);
  }
}
