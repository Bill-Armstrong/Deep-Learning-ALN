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

// alnconfidenceplimit.cpp

#ifdef ALNDLL
#define ALNIMP __declspec(dllexport)
#endif

#include <aln.h>
#include "alnpriv.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

static int ALNAPI ValidateALNConfidencePLimit(const ALNCONFIDENCE* pConfidence,
                                              double dblSignificance,
                                              double* pdblPLimit);

#ifdef _DEBUG
static void DebugValidateALNConfidencePLimit(const ALNCONFIDENCE* pConfidence,
                                             double dblSignificance,
                                             double* pdblPLimit);
#endif

// much of this is derived from theory documented in Master95 p302-323
// and Press et al p228-229

ALNIMP int ALNAPI ALNConfidencePLimit(const ALNCONFIDENCE* pConfidence,
                                      double dblSignificance,
                                      double* pdblPLimit)
{
	int nReturn = ValidateALNConfidencePLimit(pConfidence, dblSignificance, pdblPLimit);
  if (nReturn != ALN_NOERROR)
    return nReturn;
  
  #ifdef _DEBUG
    DebugValidateALNConfidencePLimit(pConfidence, dblSignificance, pdblPLimit);
  #endif

  // calc number of samples in tail
  int nTailSamples = (int)floor((float)pConfidence->nSamples * pConfidence->dblP - 1);
  
  // calculate limit at desired significance level
	*pdblPLimit = PLimit(pConfidence->nSamples, nTailSamples + 1, dblSignificance);

	return nReturn;
}


// validate params
static int ALNAPI ValidateALNConfidencePLimit(const ALNCONFIDENCE* pConfidence,
                                              double dblSignificance,
                                              double* pdblPLimit)
{
  if (pConfidence == NULL || pdblPLimit == NULL)
    return ALN_GENERIC;

  return ALN_NOERROR;
}

// debug version ASSERTS if bad params
#ifdef _DEBUG
static void DebugValidateALNConfidencePLimit(const ALNCONFIDENCE* pConfidence,
                                             double dblSignificance,
                                             double* pdblPLimit)
{
  ASSERT(dblSignificance >= 0.0 && dblSignificance <= 1.0);
}
#endif
