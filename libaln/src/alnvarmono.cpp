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

// alnvarmono.cpp

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

// ALN variable monotonicity type
// monotonicity type returned in pnMono
// returns ALN_* error code, (ALN_NOERROR on success)
ALNIMP int ALNAPI ALNVarMono(const ALN* pALN, int nVar, int* pnMono)
{
  // parameter variance
  if (pALN == NULL)
    return ALN_GENERIC;

  if (nVar < 0 || nVar >= pALN->nDim)
    return ALN_GENERIC;

  if (pnMono == NULL)
    return ALN_GENERIC;

  *pnMono = CheckMonotonicity(pALN->pTree, pALN, nVar);
  ASSERT(*pnMono == MONO_CONSTANT || *pnMono == MONO_FREE ||
         *pnMono == MONO_STRONGINC || *pnMono == MONO_STRONGDEC ||
         *pnMono == MONO_WEAKINC || *pnMono == MONO_WEAKDEC);

  return ALN_NOERROR;
}

///////////////////////////////////////////////////////////////////////////////
// workhorse of ALNVarMono

int ALNAPI CheckMonotonicity(const ALNNODE* pNode, const ALN* pALN, int nVar)
{
	ASSERT(pNode);
	ASSERT(pALN);
  ASSERT(nVar >= 0 && nVar < pALN->nDim);

  // traverse tree and examine weights on variable
  if (NODE_MINMAXTYPE(pNode) & NF_LFN)
  {
    // examine actual weight on var
    double dblW = LFN_W(pNode)[nVar + 1];   // ...skip bias weight
    if (dblW < 0)
    {
      return MONO_STRONGDEC;
    }
    else if (dblW > 0)
    {
      return MONO_STRONGINC;
    }
    // else has zero weight so return MONO_NONE
  
    return MONO_CONSTANT; // return monotonicity
  }
  
  ASSERT(NODE_MINMAXTYPE(pNode) & NF_MINMAX);

  // monotonicity initially set to undefined
  int nMono = -1;

  int nChildren = MINMAX_NUMCHILDREN(pNode);
  ALNNODE*const* apChildren = MINMAX_CHILDREN(pNode);
  for (int i = 0; i < nChildren; i++)
  { 
    const ALNNODE* pChild = apChildren[i];
    ASSERT(pChild);

    int nChildMono = CheckMonotonicity(pChild, pALN, nVar);
                  
    if (nChildMono == MONO_FREE)
      return MONO_FREE;   // child is free, so we're free
    
    else if (i == 0)  
      nMono = nChildMono; // first time through, so take on child mono

    else if (nMono == nChildMono)
      continue;           // no change

    // move from const to weak if child is weak or strong                
    else if (nMono == MONO_CONSTANT && 
             (nChildMono == MONO_WEAKINC || nChildMono == MONO_STRONGINC))
      nMono = MONO_WEAKINC;
    
    else if (nMono == MONO_CONSTANT && 
             (nChildMono == MONO_WEAKDEC || nChildMono == MONO_STRONGDEC))
      nMono = MONO_WEAKDEC;

    // remain weak if child is strong or constant
    else if (nMono == MONO_WEAKINC &&
             (nChildMono == MONO_STRONGINC || nChildMono == MONO_CONSTANT))
      continue;
    
    else if (nMono == MONO_WEAKDEC &&
             (nChildMono == MONO_STRONGDEC || nChildMono == MONO_CONSTANT))
      continue;
          
    // move from strong to weak if we're strong and child is weak or constant
    else if (nMono == MONO_STRONGINC && 
             (nChildMono == MONO_WEAKINC || nChildMono == MONO_CONSTANT))
      nMono = MONO_WEAKINC;
             
    else if (nMono == MONO_STRONGDEC && 
             (nChildMono == MONO_WEAKDEC || nChildMono == MONO_CONSTANT))
      nMono = MONO_WEAKDEC;

    // opposite child monotonicities      
    else return MONO_FREE;
  }

  // no conflicting child monotonicities
  ASSERT(nMono != -1);
  return nMono;           
}

