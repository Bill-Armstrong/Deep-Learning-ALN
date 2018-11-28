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

// calcactivechild.cpp

#ifdef ALNDLL
#define ALNIMP __declspec(dllexport)
#endif

#include <aln.h>
#include "alnpriv.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif


int ALNAPI CalcActiveChild(double& dblRespActive, double& dblDistance, 
                           double dbl0, double dbl1, const ALNNODE* pNode, 
                           double dblSE, double dbl4SE, double dblOV16SE)

{
  int nActive = -1;

  // MAX node handling
	if (MINMAX_ISMAX(pNode)) 
	{
		if(dbl1 > dbl0 + dbl4SE)	//  this puts child 1 100% active
		{
      nActive = 1;
      dblRespActive = 1.0;
      dblDistance = dbl1;
		}
		else if	(dbl1 > dbl0 - dbl4SE)	// there is a fillet
		{
			if(dbl1 > dbl0)	 // child 1 is more active
			{
        nActive = 1;
			}
			else	// child 0 is more active
			{
				nActive = 0;
			}

      double dblDiff = dbl1 - dbl0;
      dblDistance = 0.5 * (dbl1 + dbl0) + 
                    dblOV16SE * dblDiff * dblDiff + 
                    dblSE;
			
      dblRespActive = 0.5 * (1.0 + fabs(dblDiff) / dbl4SE);
		}
		else	// child 0 is 100% active
		{
			nActive = 0;
      dblRespActive = 1.0;
      dblDistance = dbl0;
		}	
	}
	else // this is a MIN node
	{    
  	if(dbl1 < dbl0 - dbl4SE) //  this puts child 1 100% active
		{
			nActive = 1;
      dblRespActive = 1.0;
      dblDistance = dbl1;
  	}
		else if (dbl1 < dbl0 + dbl4SE)	 // there is a fillet
		{
			if(dbl1 < dbl0)
			{
			  nActive = 1;
			}
			else
			{
				nActive = 0;
			}

      double dblDiff = dbl1 - dbl0;
      dblDistance = 0.5 * (dbl1 + dbl0) - 
                    dblOV16SE * (dblDiff) * (dblDiff) - 
                    dblSE;
			
      dblRespActive = 0.5 * (1.0 + fabs(dblDiff) / dbl4SE);
		}
		else	 // child 0 is 100% active
		{
			nActive = 0;
      dblRespActive = 1.0;
      dblDistance = dbl0;
		}
	}

  ASSERT(nActive != -1);
  return nActive;
}
