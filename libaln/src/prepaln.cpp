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
// prepaln.cpp

//  Archive: /ALN Development/libaln/src/prepaln.cpp
//  Workfile: prepaln.cpp 
//  Date: October 24, 2018
//  Changes by: W. W. Armstrong

#ifdef ALNDLL
#define ALNIMP __declspec(dllexport)
#endif

#include <aln.h>
#include "alnpriv.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

BOOL ALNAPI DoPrepRegions(ALN* pALN);
BOOL ALNAPI DoPrepNode(ALN* pALN, ALNNODE* pNode);

BOOL ALNAPI PrepALN(ALN* pALN)
{
  ASSERT(pALN);

  if (!DoPrepRegions(pALN))
    return FALSE;

  if (!DoPrepNode(pALN, pALN->pTree))
    return FALSE;

  return TRUE;
}

BOOL ALNAPI DoPrepRegions(ALN* pALN)
{
  ASSERT(pALN);

  // calc all region and var quantities

  for (int i = 0; i < pALN->nRegions; i++)
  {
    ALNREGION* pRegion = pALN->aRegions + i;
    
    // iterate over var constraints
    for (int j = 0; j < pRegion->nConstr; j++)
    {
      double dblEpsilon = pRegion->aConstr[j].dblEpsilon;

      if (j == pALN->nOutput)
      {
        // make sure wmin and wmax are -1
        pRegion->aConstr[j].dblWMin = -1.0;
        pRegion->aConstr[j].dblWMax = -1.0;
      }

      // calc sq epsilon
      pRegion->aConstr[j].dblSqEpsilon = dblEpsilon * dblEpsilon; // this may be useless now
    }

    // calc smoothing epsilon quantities
    SetSmoothingEpsilon(pRegion);
  }

  return TRUE;
}

BOOL ALNAPI DoPrepNode(ALN* pALN, ALNNODE* pNode)
{
  // currently no prepping required here
  return TRUE;


  // otherwise, recurse over entire tree
  /*
  if (NODE_ISLFN(pNode))
  {
  }
  else
  {
    ASSERT(NODE_ISMINMAX(pNode));
		DoPrepNode(MINMAX_LEFT(pNode));
    DoPrepNode(MINMAX_RIGHT(pNode));
  }

  return TRUE;
  */
}

