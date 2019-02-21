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

// debugevalminmax.cpp


#ifdef ALNDLL
#define ALNIMP __declspec(dllexport)
#endif

#include <aln.h>
#include "alnpriv.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

#ifdef _DEBUG

double ALNAPI DebugEvalMinMax(const ALNNODE* pNode, const ALN* pALN, 
                            const double* adblX, ALNNODE** ppActiveLFN)
{
  ASSERT(NODE_ISMINMAX(pNode));

  // set first child
  ALNNODE* pChild0;
  if (MINMAX_EVAL(pNode))
    pChild0 = MINMAX_EVAL(pNode);
  else
    pChild0 = MINMAX_LEFT(pNode);

  // set next child
  ALNNODE* pChild1;
  if (pChild0 == MINMAX_LEFT(pNode))
    pChild1 = MINMAX_RIGHT(pNode);
  else
    pChild1 = MINMAX_LEFT(pNode);

  // eval first child
  ALNNODE* pActiveLFN0;
  double dbl0 = DebugEval(pChild0, pALN, adblX, &pActiveLFN0);
  
	// eval second child
  ALNNODE* pActiveLFN1;
	double dbl1 = DebugEval(pChild1, pALN, adblX, &pActiveLFN1);

  // get reference to region for this node
  const ALNREGION& region = pALN->aRegions[NODE_REGION(pNode)];

  // calc active child, active child response, and distance
  double dblDist, dblRespActive;
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
  
  return dblDist;
}


#endif