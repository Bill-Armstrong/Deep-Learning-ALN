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

// getvarconstraint.cpp

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
// constraint retrieval

ALNCONSTRAINT* ALNAPI GetVarConstraint(int nRegion, const ALN* pALN, 
                                       int nVar)
{
  ASSERT(pALN);
  ASSERT(nRegion >= 0 && nRegion < pALN->nRegions);
  ASSERT(nVar >= 0 && nVar < pALN->nDim);
  ALNREGION* pRegion = &(pALN->aRegions[nRegion]);

  while (TRUE)
  {
    // if region constrains all vars, then get constraint directly
    if (pRegion->nConstr == pALN->nDim)
    {
      ASSERT(pRegion->aConstr[nVar].nVarIndex == nVar);
      return &(pRegion->aConstr[nVar]);
    }
    else if (pRegion->afVarMap == NULL || TESTMAP(pRegion->afVarMap, nVar))
    {
      // no var map, or map indicates region contains this var
      for (int i = pRegion->nConstr - 1; i >= 0; i--)
      {
        ASSERT(pRegion->aConstr[i].nVarIndex >= 0 && 
               pRegion->aConstr[i].nVarIndex < pALN->nDim);
      
        if (pRegion->aConstr[i].nVarIndex == nVar)
          return &(pRegion->aConstr[i]);
      }
    }
    
    // variable not constrained in this region, check parent region
    if (pRegion->nParentRegion == -1)
      return NULL;  // var not found!?
    
    pRegion = &(pALN->aRegions[pRegion->nParentRegion]); 
  }
}

