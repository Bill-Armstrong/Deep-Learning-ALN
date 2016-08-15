// ALN Library
// Copyright (C) 1995 - 2010 William W. Armstrong.
// 
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
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

// jitter.cpp

////////////////////////////////////////////////////////////////////////////
//  File version info:
// 
//  $Archive: /ALN Development/libaln/src/jitter.cpp $
//  $Workfile: jitter.cpp $
//  $Revision: 6 $
//  $Date: 7/17/07 6:46p $
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

inline
double ALNAPI Noise()
{
/*
   Generate point in triangular distribution from -1 to 1

   We could do this by generating two random uniformly distributed
   points x, y in [-1, 1] and rejecting any |x|+|y| > 1.
   On average this method will reject half the points we pick!

   Instead, we generate points x, y in [-1/(sqrt(2)), 1/(sqrt(2))]
  
   if r1 and r2 are uniform random numbers in [0, 1], then
  
   x = (2/sqrt(2))(r1 - 1/2) -> x in [-1/sqrt(2), 1/sqrt(2)]
   y = (2/sqrt(2))(r2 - 1/2) -> y in [-1/sqrt(2), 1/sqrt(2)]
  
   Now rotate 45 degrees...
   
   1/sqrt(2) |  1  -1  | | x | = | x' |
             |  1   1  | | y | = | y' |
  
   x' = (1/sqrt(2))(x - y)  
   y' = (1/sqrt(2))(x + y)  
  
   substituting for x and y
  
   x' = (1/sqrt(2))(((2/sqrt(2))(r1 - 1/2)) - ((2/sqrt(2))(r2 - 1/2)))
      = r1 - r2
  
   y' = (1/sqrt(2))(((2/sqrt(2))(r1 - 1/2)) + ((2/sqrt(2))(r2 - 1/2)))
      = r1 + r2 - 1
*/

#ifdef _DEBUG
  double dbl = ALNRandFloat() - ALNRandFloat();
  ASSERT(dbl >= -1.0 && dbl <= 1.0);
  return dbl;
#else
  return ALNRandFloat() - ALNRandFloat();
#endif
}

void ALNAPI Jitter(ALN* pALN, double* adblX)
{
  ASSERT(pALN);
  ASSERT(adblX);

  int nDim = pALN->nDim;
  int nOutput = pALN->nOutput;
  
  // save output value
  double dblOutput = adblX[nOutput];

  for (int i = 0; i < nDim; i++)
  {
#ifdef _DEBUG
    double dbl = adblX[i];
#endif

    adblX[i] += Noise() * pALN->aRegions[0].aConstr[i].dblEpsilon;

#ifdef _DEBUG
    double dblEps = pALN->aRegions[0].aConstr[i].dblEpsilon;
    ASSERT(adblX[i] >= dbl - dblEps && 
           adblX[i] <= dbl + dblEps);
#endif
  }

  // restore output value
  adblX[nOutput] = dblOutput;
}