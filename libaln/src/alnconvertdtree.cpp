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

// alnconvertdtree.cpp

#ifdef ALNDLL
#define ALNIMP __declspec(dllexport)
#endif

#include <aln.h>
#include "alnpriv.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

// conversion to dtree
// pointer to constructed DTREE returned in ppDtree
// returns ALN_* error code, (ALN_NOERROR on success)
// if ALN_ERRDTREE is returned, check dtree_errno
ALNIMP int ALNAPI ALNConvertDtree(const ALN* pALN, int nMaxDepth,
                                  DTREE** ppDtree)
{
  // parameter variance
  if (pALN == NULL)
    return ALN_GENERIC;

  if (ppDtree == NULL)
    return ALN_GENERIC;

  if (nMaxDepth < DTREE_MINDEPTH || nMaxDepth > DTREE_MAXDEPTH)
    return ALN_GENERIC;

  int nResult = ALN_NOERROR;

  // build dtree
  *ppDtree = BuildDtree(pALN, nMaxDepth);
  if (*ppDtree == NULL && dtree_errno != DTR_NOERROR)
    nResult = ALN_GENERIC;  // dtree lib error
  else if (*ppDtree == NULL)
    nResult = ALN_OUTOFMEM;
  
  return nResult;
}