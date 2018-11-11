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

// plimit.cpp


#ifdef ALNDLL
#define ALNIMP __declspec(dllexport)
#endif

#include <aln.h>
#include "alnpriv.h"
#include <boost\math\special_functions\beta.hpp>
using namespace boost::math;

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

// calculate probability p of an event occuring, such that
// the probablity of m or less such events occuring in n trials
// is x; currently limited to accuracy of 1.e-7

// implementation based on Ridders method documented in Press et al

// method: advance p until cumulative binomial dist 0 to m events in n 
//         trials drops to x

double ALNAPI PLimit(int n, int m, double dblX)
{
  static const double dblInc = 0.1;     // coarse increment
  static const double dblAcc = 1.0e-7;  // maximum accuracy

  if (dblX < 0.0 || dblX > 1.0 || n < 0)
  {
    return NAN;
  }
  
  // known cases
  if (m < 0)
    return 0.0;
  if (m >= n)
    return 1.0;

 
  // P is desired probability, Y is difference between desired area
  // under tail and the area under the tail given P

  // therfore, Y goes to 0 as P approaches desired value

  // lower bound
  double dblP1 = 0.0;
  double dblY1 = dblX - 1.0;
  ASSERT(dblY1 <= 0.0);

  // begin coarse approximation of P
  
  // scan upward until Y3 is +ve
  double dblP3, dblY3;
  for (dblP3 = dblInc; dblP3 < 1.0; dblP3 += dblInc)
  {
    dblY3 = dblX - (1.0 - ibeta(m + 1, n - m, dblP3));  //ibet??
                   // cumulative binomial dist 0 to m events in n trials, 
                   // see Press et al p229

    // convergence test (unlikely at this point)
    if (fabs(dblY3) < dblAcc)
      return dblP3;

    // check for sign change
    if (dblY3 > 0.0)
      break;  // we've bracketed desired value

    // else, new lower bound
    dblP1 = dblP3;
    dblY1 = dblY3;
  }

  // P1 and P3 bracket desired value... refine using ridders method
  const int nMaxIt = 100;
  
  for (int i = 0; i < nMaxIt; i++)
  {
    // get mid-values
    double dblP2 = 0.5 * (dblP1 + dblP3);
    if ((dblP3 - dblP1) < dblAcc)   // convergence test
      return dblP2;
    
    double dblY2 = dblX - (1.0 - ibeta(m + 1, n - m, dblP2)); //ibeta??

    // convergence test
    if (fabs(dblY2) < dblAcc)
      return dblP2;

    double dblDenom = sqrt(dblY2 * dblY2 - dblY1 * dblY3);  // y1, y3 opposite sign
    double dblTrial = dblP2 + (dblP1 - dblP2) * dblY2 / dblDenom;

    double dblY = dblX - (1.0 - ibeta(m + 1, n - m, dblTrial));  //ibeta

    // convergence test
    if (fabs(dblY) < dblAcc)
      return dblTrial;

    // between mid and test point?
    if ((dblY2 < 0.0) && (dblY > 0.0))
    {
      dblP1 = dblP2;    // new lower bound
      dblY1 = dblY2;
      dblP3 = dblTrial; // new upper bound
      dblY3 = dblY;
    }
    else if ((dblY < 0.0) && (dblY2 > 0.0))
    {
      dblP1 = dblTrial; // new lower bound
      dblY1 = dblY;
      dblP3 = dblP2;    // new upper bound
      dblY3 = dblY2;
    }
    else if (dblY < 0.0)  // both negative
    {
      dblP1 = dblTrial;
      dblY1 = dblY;
    }
    else  // both positive
    {
      dblP3 = dblTrial;
      dblY3 = dblY;
    }
  }

  // convergence failed... return best guess?
  return 0.5 * (dblP1 + dblP3); 
}