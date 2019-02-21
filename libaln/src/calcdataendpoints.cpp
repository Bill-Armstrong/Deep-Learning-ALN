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

// calcdataendpoints.cpp

#ifdef ALNDLL
#define ALNIMP __declspec(dllexport)
#endif

#include <aln.h>
#include "alnpriv.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

void ALNAPI CalcDataEndPoints(long& nStart, long& nEnd, const ALN* pALN,
                              const ALNDATAINFO* pDataInfo)
{
  ASSERT(pALN);
  ASSERT(pDataInfo);

  int nDim = pALN->nDim;
  const VARINFO* aVarInfo = pDataInfo->aVarInfo;

  nStart = 0; 
  nEnd = 0;
  if (aVarInfo != NULL)
  {
    // scan var info structs and account for any time shifts
    for (int i = 0; i < nDim; i++)
    {
      if (aVarInfo[i].nDelta < 0 && abs(aVarInfo[i].nDelta) > nStart)
        nStart = abs(aVarInfo[i].nDelta);
      else if (aVarInfo[i].nDelta > 0 && aVarInfo[i].nDelta > nEnd)
        nEnd = aVarInfo[i].nDelta;
    }

    ASSERT(nStart + nEnd <= pDataInfo->nPoints);
    
    // adjust nEnd to reflect end point
    nEnd = pDataInfo->nPoints - nEnd - 1;
  }
  else
  {
    nStart = 0;
    nEnd = pDataInfo->nPoints - 1;
  }

  ASSERT(nStart >= 0 && nStart < pDataInfo->nPoints && nStart <= nEnd);
  ASSERT(nEnd >= 0 && nEnd < pDataInfo->nPoints);
}
