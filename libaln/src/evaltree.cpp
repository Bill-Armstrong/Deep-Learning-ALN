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

// evaltree.cpp
// eval support routines


#ifdef ALNDLL
#define ALNIMP __declspec(dllexport)
#endif

#include <aln.h>
#include "alnpriv.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

#ifdef _DEBUG
static void DebugValidateEvalTreeInfo(const ALN* pALN,
                                      const ALNDATAINFO* pDataInfo,
                                      const ALNCALLBACKINFO* pCallbackInfo,
                                      double* adblResult,
                                      int* pnStart, int* pnEnd,
                                      ALNNODE** apActiveLFNs,
                                      double* adblInput,
                                      double* adblOutput);
#endif

// evaluation of ALN on data

int ALNAPI EvalTree(const ALNNODE* pNode, 
                    const ALN* pALN,
                    const ALNDATAINFO* pDataInfo,
                    const ALNCALLBACKINFO* pCallbackInfo,
                    double* adblResult,
                    int* pnStart, int* pnEnd,
                    BOOL bErrorResults /*= FALSE*/,
                    ALNNODE** apActiveLFNs /*= NULL*/,
                    double* adblInput /*= NULL*/,
                    double* adblOutput /*= NULL*/)
{
  ASSERT(pNode);
#ifdef _DEBUG
  DebugValidateEvalTreeInfo(pALN, pDataInfo, pCallbackInfo,
                            adblResult, pnStart, pnEnd, 
                            apActiveLFNs, adblInput, adblOutput);
#endif
  
  int nDim = pALN->nDim;
  int nPoints = pDataInfo->nPoints;

  // calc start and end points
  long nStart, nEnd; 
  CalcDataEndPoints(nStart, nEnd, pALN, pDataInfo);
  
  if (pnStart != NULL)
    *pnStart = nStart;
  if (pnEnd != NULL)
    *pnEnd = nEnd;

  // evaluation loop
  int nReturn = ALN_NOERROR;        // assume OK
  double* adblX = NULL;             // eval vector
  const double** apdblBase = NULL;  // column base ptr
  ALNNODE* pTree = pALN->pTree;		  // on stack for quicker access
  CCutoffInfo* aCutoffInfo = NULL;  
 
	try
 	{
    // reset active lfn array
    if (apActiveLFNs != NULL)
    {
      memset(apActiveLFNs, 0, pDataInfo->nPoints * sizeof(ALNNODE*));
    }

  	// allocate input vector     
   	adblX = new double[nDim];   
    if (!adblX) ThrowALNMemoryException();
    memset(adblX, 0, sizeof(double) * nDim);

   	// allocate column base vector
    apdblBase = AllocColumnBase(nStart, pALN, pDataInfo);

    // main loop
    ALNNODE* pActiveLFN = NULL;
    for (int i = nStart; i <= nEnd; i++)
    {
      // fill input vector
      FillInputVector(pALN, adblX, i - nStart, nStart, apdblBase, 
                      pDataInfo, pCallbackInfo);

      // copy input vector?
      if (adblInput)
      {
        // get the input row
        double* adblRow = adblInput + (i * pALN->nDim);

        // copy values
        memcpy(adblRow, adblX, pALN->nDim * sizeof(double));
        
        // set the bias value in the output var spot
        adblRow[pALN->nOutput] = 1.0;
      }

      // copy desired output 
      // ... do this before setting output value in input vector to zero below
      if (adblOutput)
      {
        adblOutput[i] = adblX[pALN->nOutput];
      }

      // CutoffEval returns distance from surface to point in the direction of
	    // the output variable, so we need to add that to the existing output value
      // to get the actual surface value
      if (!bErrorResults)
      {
        adblX[pALN->nOutput] = 0; // set output value to zero...

        // ... since output value is zero, the distance CutoffEval returns
        // is the value of the function surface
      }

      // get the distance from the point to the surface defined by the ALN
      adblResult[i] = CutoffEval(pTree, pALN, adblX, CEvalCutoff(),
                                 &pActiveLFN);
      
      // save the active LFN
      if (apActiveLFNs != NULL)
      {
        apActiveLFNs[i] = pActiveLFN;
      }
    }
  }
  catch (CALNUserException* e)
  {
  	nReturn = ALN_USERABORT;
    e->Delete();
  }
  catch (CALNMemoryException* e)
  {
  	nReturn = ALN_OUTOFMEM;
    e->Delete();
  }
  catch (CALNException* e)
  {
  	nReturn = ALN_GENERIC;
    e->Delete();
  }
  catch (...)
  {
  	nReturn = ALN_GENERIC;
  }

  // clear memory	
	delete[] adblX;
  FreeColumnBase(apdblBase);
	
  return nReturn;
}

// debug version ASSERTS if bad params
#ifdef _DEBUG
static void DebugValidateEvalTreeInfo(const ALN* pALN,
                                      const ALNDATAINFO* pDataInfo,
                                      const ALNCALLBACKINFO* pCallbackInfo,
                                      double* adblResult,
                                      int* pnStart, int* pnEnd,
                                      ALNNODE** apActiveLFNs,
                                      double* adblInput,
                                      double* adblOutput)
{
  DebugValidateALNDataInfo(pALN, pDataInfo, pCallbackInfo);
  
  ASSERT(adblResult != NULL);
}
#endif
