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

// adaptminmax.cpp


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
// minmax specific adapt

static const double dblRespThresh = 1.0;//.9999999999; 
static const double dblRespMin = 0.0;

void ALNAPI AdaptMinMax(ALNNODE* pNode, ALN* pALN, const double* adblX, 
                      double dblResponse, BOOL bUsefulAdapt,
                      const TRAINDATA* ptdata)
{
  ASSERT(NODE_ISMINMAX(pNode));
 	ASSERT(NODE_ISEVAL(pNode));
	ASSERT(ptdata != NULL);

	if (bUsefulAdapt)
	{
    NODE_RESPCOUNT(pNode)++;
	}

  // get output var constraint
	ALNCONSTRAINT* pConstrOutput = GetVarConstraint(NODE_REGION(pNode), pALN, pALN->nOutput);
  
  // get children
  ALNNODE* pChild0 = MINMAX_LEFT(pNode);
	ALNNODE* pChild1 = MINMAX_RIGHT(pNode);
  BOOL bUsefulAdapt0 = FALSE;
  BOOL bUsefulAdapt1 = FALSE;

  // get resp counts
  int nResp0 = NODE_RESPCOUNT(pChild0) + NODE_RESPCOUNTLASTEPOCH(pChild0);
  int nResp1 = NODE_RESPCOUNT(pChild1) + NODE_RESPCOUNTLASTEPOCH(pChild1);

  // calculate the responsibilities of the children
  double dblResp0, dblResp1;
	double dblRespActive = MINMAX_RESPACTIVE(pNode);

  ASSERT(MINMAX_ACTIVE(pNode) != NULL);
  if (MINMAX_ACTIVE(pNode) == pChild0) // child 0 active
  {
    bUsefulAdapt0 = bUsefulAdapt;

    double dblR;

    if((fabs(ptdata->dblGlobalError) > pConstrOutput->dblEpsilon) && 
       (dblRespActive > dblRespThresh) && (nResp1 < nResp0))
		{
      // bring in useless piece 1
			dblR = dblRespThresh;
      
      // eval if necessary before adapting
			if(!NODE_ISEVAL(pChild1))
			{
        ALNNODE* pActiveLFN1;
        int nCutoffs = 0;
				AdaptEval(pChild1, pALN, adblX, CEvalCutoff(), &pActiveLFN1);
			}
		}
		else
		{
			dblR = dblRespActive;
		}

    // divide each quantity by 1 - 2r(1-r)

    double dblFactor = 1.0 / (1 - 2 * dblR * (1 - dblR));

    dblResp0 = dblR * dblFactor;
		dblResp1 = (1.0 - dblR) * dblFactor;
  }
	else // child 1 is active
	{
    bUsefulAdapt1 = bUsefulAdapt;

    double dblR;

    if((fabs(ptdata->dblGlobalError) > pConstrOutput->dblEpsilon) && 
       (dblRespActive > dblRespThresh) && (nResp0 < nResp1))
    { 
      // bring in useless piece 0
			dblR = dblRespThresh;
			
      // eval child 0 if necessary before adapting
			if(!NODE_ISEVAL(pChild0))
			{
        ALNNODE* pActiveLFN0;
        int nCutoffs = 0;
				AdaptEval(pChild0, pALN, adblX, CEvalCutoff(), &pActiveLFN0);
			}
		}
		else
		{
			dblR = dblRespActive;
		}
		    
    // divide each quantity by 1 - 2r(1-r)

    double dblFactor = 1.0 / (1 - 2 * dblR * (1 - dblR));

    dblResp1 = dblR * dblFactor;
		dblResp0 = (1.0 - dblR) * dblFactor;
	}
  
  // adapt child 0
  dblResp0 *= dblResponse;
  if(dblResp0 > dblRespMin)
	{
    Adapt(pChild0, pALN, adblX, dblResp0, bUsefulAdapt0, ptdata);
	}

  // adapt child 1
  dblResp1 *= dblResponse;
  if(dblResp1 > dblRespMin)
	{
    Adapt(pChild1, pALN, adblX, dblResp1, bUsefulAdapt1, ptdata);
  }
}


