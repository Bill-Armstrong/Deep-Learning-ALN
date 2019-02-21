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

// shuffle.cpp

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
// shuffle

void ALNAPI Shuffle(int nStart, int nEnd, int* anShuffle)
{
  ASSERT(anShuffle);

  if ((nEnd - nStart) > 1)
  {
  	for (int nSwap = nStart; nSwap <= nEnd; nSwap++)
  	{
  		// calc swap indexes
  	  int a, b;
  	  a = ALNRand() % (nEnd - nStart + 1);
  	  do { b = ALNRand() % (nEnd - nStart + 1); } while (a == b);

  	  // swap indexes
  		#define _SWAP(a, b) (a) = (a)^(b); (b) = (a)^(b); (a) = (a)^(b)
  		_SWAP(anShuffle[a], anShuffle[b]);
  		#undef _SWAP
  	} // end shuffle    
  }
}
