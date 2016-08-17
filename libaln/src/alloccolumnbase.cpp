// ALN Library
// Copyright (C) 1995 - 2010 William W. Armstrong.
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

// alloccolumnbase.cpp

///////////////////////////////////////////////////////////////////////////////
//  File version info:
// 
//  $Archive: /ALN Development/libaln/src/alloccolumnbase.cpp $
//  $Workfile: alloccolumnbase.cpp $
//  $Revision: 5 $
//  $Date: 7/17/07 5:08p $
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

const double** ALNAPI AllocColumnBase(int nStart,
                                      const ALN* pALN,
                                      const ALNDATAINFO* pDataInfo)
{
  ASSERT(pALN);
  ASSERT(pDataInfo);

  int nDim = pALN->nDim;
  int nCols = pDataInfo->nCols;
  const VARINFO* aVarInfo = pDataInfo->aVarInfo;
  const double* adblData = pDataInfo->adblData;

  const double** apdblBase = NULL;

  if (aVarInfo != NULL && adblData != NULL)
  {
    apdblBase = new const double*[nDim];
    if (!apdblBase) ThrowALNMemoryException();
    
    for (int i = 0; i < nDim; i++)
    {
      int nCol = (aVarInfo[i].nCol == -1) ? i : aVarInfo[i].nCol;
      apdblBase[i] = adblData + ((nStart + aVarInfo[i].nDelta) * nCols) + nCol;
    }
  }

  return apdblBase;
}