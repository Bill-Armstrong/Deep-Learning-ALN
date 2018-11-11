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

#ifdef ALNDLL
#define ALNIMP __declspec(dllexport)
#endif

#include <aln.h>
#include "alnpriv.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

int ALNAPI ValidateALNDataInfo(const ALN* pALN,
                               const ALNDATAINFO* pDataInfo,
                               const ALNCALLBACKINFO* pCallbackInfo)
{
  // parameter variance
  if (pDataInfo == NULL)
  {
    return ALN_GENERIC;
  }

  if (pALN == NULL)
  {
    return ALN_GENERIC;
  }

  // must have at least one training point
  if (pDataInfo->nPoints <= 0)
  {
    return ALN_GENERIC;
  }

  // need proc if no data
  if(pDataInfo->adblData == NULL && 
     (pCallbackInfo == NULL || pCallbackInfo->pfnNotifyProc == NULL)   )
  {
    return ALN_GENERIC;
  }

  // need AN_VECTORINFO if no data
  if(pDataInfo->adblData == NULL && 
     (pCallbackInfo == NULL || !(pCallbackInfo->nNotifyMask & AN_VECTORINFO)))
  {
    return ALN_GENERIC;
  }

  // make sure columns valid
  if (pDataInfo->aVarInfo == NULL && 
      pDataInfo->adblData != NULL && 
      pDataInfo->nCols < pALN->nDim)
  {
    return ALN_GENERIC;
  }
  
  // check var info column, delta validity
  const VARINFO* aVarInfo = pDataInfo->aVarInfo;
  if (aVarInfo != NULL)
  {
    const double* adblData = pDataInfo->adblData;
    int nCols = pDataInfo->nCols;
    int nDim = pALN->nDim;
    int nEnd = 0;   // points from end
    int nStart = 0; // points from start
    for (int i = 0; i < nDim; i++)
    {
      if (adblData != NULL && 
          (aVarInfo[i].nCol < 0 || aVarInfo[i].nCol >= nCols))
        return ALN_GENERIC;

      if (aVarInfo[i].nDelta < 0 && abs(aVarInfo[i].nDelta) > nStart)
        nStart = abs(aVarInfo[i].nDelta);
      else if (aVarInfo[i].nDelta > 0 && aVarInfo[i].nDelta > nEnd)
        nEnd = aVarInfo[i].nDelta;
    }

    if (nEnd + nStart > pDataInfo->nPoints)  // deltas exceed number of points!
      return ALN_GENERIC;
  }

  return ALN_NOERROR;
}

#ifdef _DEBUG
void ALNAPI DebugValidateALNDataInfo(const ALN* pALN,
                                     const ALNDATAINFO* pDataInfo,
                                     const ALNCALLBACKINFO* pCallbackInfo)
{
  ASSERT(pDataInfo != NULL);

  // valid aln pointer
	ASSERT(pALN != NULL);

  // valid number of points
  ASSERT(pDataInfo->nPoints > 0);

  // valid data cols
  ASSERT((pDataInfo->adblData != NULL && pDataInfo->nCols > 0) || 
         pDataInfo->adblData == NULL);

  // valid notify proc
  ASSERT(pDataInfo->adblData != NULL || 
         (pCallbackInfo != NULL && 
          pCallbackInfo->pfnNotifyProc != NULL && 
          (pCallbackInfo->nNotifyMask & AN_VECTORINFO)));
  
  // valid varinfo
  ASSERT(pDataInfo->aVarInfo != NULL || pDataInfo->nCols >= pALN->nDim);
  if (pDataInfo->aVarInfo != NULL)
  {
    for (int i = 0; i < pALN->nDim; i++)
    {
      // column validity
      ASSERT(!pDataInfo->adblData || 
             (pDataInfo->aVarInfo[i].nCol >= 0 || 
              pDataInfo->aVarInfo[i].nCol < pDataInfo->nCols));
    }
  }
}
#endif
