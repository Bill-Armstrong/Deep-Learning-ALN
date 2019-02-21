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

// alnquickeval.cpp

#ifdef ALNDLL
#define ALNIMP __declspec(dllexport)
#endif

#include <aln.h>
#include "alnpriv.h"
#include <errno.h>

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

// evaluation of ALN on single vector
// ALNQuickEval returns result of evaluating a single input vector, which
// must contain pALN->nDim elements
// the pointer to the responsible LFN that forms the surface is returned in
//   ppActiveLFN, which must be non-NULL
// the evaluation does not change any members of the ALN structure, which is
//   useful for applications like reinforcement learning, which may require
//   intermediate, non-destructive evaluations of the ALN during training
// returns ALN_* error code, (ALN_NOERROR on success)
// NOTE: currently only supports evaluation using default output variable
// NOTE: for efficiency reasons, there is _no_ parameter checking performed
//   and no error return value
ALNIMP double ALNAPI ALNQuickEval(const ALN* pALN, const double* adblX,
                                  ALNNODE** ppActiveLFN)
{
  ASSERT(pALN);
  ASSERT(adblX);

	// CUTOFF_EVAL returns distance from surface to point in the direction of
	// the default output variable, so we need to add that in to get the actual
	// surface value

  ALNNODE* pActiveLFN;
  double dbl =  adblX[pALN->nOutput] + CutoffEval(pALN->pTree, pALN, adblX, 
                                                  CEvalCutoff(), &pActiveLFN);
  if (ppActiveLFN)
    *ppActiveLFN = pActiveLFN;

  return dbl;
}
