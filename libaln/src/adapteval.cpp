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

// adapteval.cpp


#ifdef ALNDLL
#define ALNIMP __declspec(dllexport)
#endif

#include <aln.h>
#include "alnpriv.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

double ALNAPI AdaptEval(ALNNODE* pNode, ALN* pALN, const double* adblX, 
                        CCutoffInfo* pCutoffInfo, ALNNODE** ppActiveLFN)
{
  ASSERT(pNode);
  ASSERT(pALN);
  ASSERT(adblX);
  ASSERT(ppActiveLFN);
  
  // do a cutoff eval to get active LFN and distance
  ALNNODE* pActiveLFN = NULL;
  double dbl;
  CEvalCutoff cutoff;
 
  // check for cutoff info
  if (pCutoffInfo != NULL)
  {
    // set up cutoff
    ALNNODE* pEval = pCutoffInfo->pLFN;
    if (pEval != NULL)
    {
      BuildCutoffRoute(pEval);
    }

    // evaluate using cutoff
    dbl = AdaptEval(pNode, pALN, adblX, cutoff, &pActiveLFN);

    // set new cutoff info
    pCutoffInfo->pLFN = pActiveLFN;
    pCutoffInfo->dblValue = dbl;
  }
  else
  {
    // eval with expanded cutoff
    dbl = AdaptEval(pNode, pALN, adblX, cutoff, &pActiveLFN);
  }

#ifdef _DEBUG
  ALNNODE* pLFNCheck = NULL;
  double dblCheck = DebugEval(pNode, pALN, adblX, &pLFNCheck);
	ASSERT(dbl == dblCheck);
	ASSERT(pLFNCheck == pActiveLFN);
#endif

  *ppActiveLFN = pActiveLFN;

  return dbl;
}
