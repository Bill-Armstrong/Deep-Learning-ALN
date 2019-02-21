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

// calcrmserror.cpp


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
// calc RMS error on data set

ALNIMP int ALNAPI ALNCalcRMSError(const ALN* pALN,
                                  const ALNDATAINFO* pDataInfo,
                                  const ALNCALLBACKINFO* pCallbackInfo,
                                  double* pdblRMSErr)
{
  int nReturn = ValidateALNDataInfo(pALN, pDataInfo, pCallbackInfo);

  try
  {
    *pdblRMSErr = DoCalcRMSError(pALN, pDataInfo, pCallbackInfo);
  }
	catch(CALNUserException* e)
  {
    nReturn = ALN_USERABORT;
    e->Delete();
    *pdblRMSErr = -1.0;
  }
  catch (CALNMemoryException* e)	// memory specific exceptions
	{
		nReturn = ALN_OUTOFMEM;
    *pdblRMSErr = -1.0;
    e->Delete();
	}
	catch (CALNException* e)	      // anything other exception we recognize
	{
		nReturn = ALN_GENERIC;
    *pdblRMSErr = -1.0;
    e->Delete();
  }
  catch(...)
  {
    nReturn = ALN_GENERIC;
    *pdblRMSErr = -1.0;
  }

  return nReturn;
}

double ALNAPI DoCalcRMSError(const ALN* pALN,
                             const ALNDATAINFO* pDataInfo,
                             const ALNCALLBACKINFO* pCallbackInfo)
{
#ifdef _DEBUG
  DebugValidateALNDataInfo(pALN, pDataInfo, pCallbackInfo);
#endif

  long nStart, nEnd;
  CalcDataEndPoints(nStart, nEnd, pALN, pDataInfo);
  
  double dblRMSError = -1.0;
	ALNNODE* pTree = pALN->pTree;
  int nDim = pALN->nDim;
  double* adblX = NULL;
  const double** apdblBase = NULL;
  CCutoffInfo* aCutoffInfo = NULL;  
  
  try
  {
    // allocate eval vector
    adblX = new double[nDim];
    if (!adblX) ThrowALNMemoryException();
    memset(adblX, 0, sizeof(double) * nDim);

    // allocate column base vector
    apdblBase = AllocColumnBase(nStart, pALN, pDataInfo);

    // allocate and init cutoff info array
    aCutoffInfo = new CCutoffInfo[nEnd - nStart + 1];
    if (!aCutoffInfo) ThrowALNMemoryException();
		
    for (int i = nStart; i <= nEnd; i++)
			aCutoffInfo[i - nStart].pLFN = NULL;

    // calc rms error
    double dblSqErrorSum = 0;
    for (int nPoint = nStart; nPoint <= nEnd; nPoint++)
	  {
      // get vector (cvt to zero based point index)
      FillInputVector(pALN, adblX, nPoint - nStart, nStart, apdblBase, 
                      pDataInfo, pCallbackInfo);
		
      // do an eval to get active LFN and distance
      ALNNODE* pActiveLFN = NULL;
      CCutoffInfo& cutoffinfo = aCutoffInfo[nPoint - nStart];
      double dbl = CutoffEval(pTree, pALN, adblX, &cutoffinfo, &pActiveLFN);
      
		  // now add square of distance from surface to error
		  dblSqErrorSum += dbl * dbl;
	  }	// end for each point

	  dblRMSError = sqrt(dblSqErrorSum / (nEnd - nStart + 1));
  }
  catch(...)
  {
    delete[] adblX;
    delete[] aCutoffInfo;
    FreeColumnBase(apdblBase);  

    throw;
  }
  
  delete[] adblX;
  delete[] aCutoffInfo;
  FreeColumnBase(apdblBase);  

  return dblRMSError;
}
