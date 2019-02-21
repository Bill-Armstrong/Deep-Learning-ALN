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

// adaptlfn.cpp

#ifdef ALNDLL
#define ALNIMP __declspec(dllexport)
#endif

#include <aln.h>
#include "alnpriv.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

extern BOOL bTrainNV_ALN;

///////////////////////////////////////////////////////////////////////////////
// LFN specific adapt

void ALNAPI AdaptLFN(ALNNODE* pNode, ALN* pALN, const double* adblX, 
                     double dblResponse, BOOL bUsefulAdapt, const TRAINDATA* ptdata)
{
  ASSERT(NODE_ISLFN(pNode));
	ASSERT(LFN_ISINIT(pNode));
  ASSERT(LFN_VARMAP(pNode) == NULL);      // var map not yet supported
  ASSERT(LFN_VDIM(pNode) == pALN->nDim);  // no different sized vectors yet
  ASSERT(NODE_ISEVAL(pNode));
	// constraining region
	ASSERT(NODE_REGION(pNode) >= 0 && NODE_REGION(pNode) < pALN->nRegions);
	ALNREGION& region = pALN->aRegions[NODE_REGION(pNode)];
  // get output var constraint
  int nOutput = pALN->nOutput;
  // count useful adapts
  if (bUsefulAdapt)
  {
    NODE_RESPCOUNT(pNode)++; 
  }
  if (!LFN_CANSPLIT(pNode) || NODE_ISCONSTANT(pNode))
  {
    return;                     // no adaptation of this constant subtree
  }
	// This procedure adapts the nDim centroid values and the nDim - 1 weight values so that the changes
	// are likely to be individually small but together correct about fraction dblLearnrate of the error
	// of the surface w. r. t. the training point. The error is the distance in the output axis 
  // from the current datapoint to the ALN surface (negative for points above surface)
	// We have to think of the error as that of the *surface*, not the data point!
	double dblError = ptdata->dblGlobalError;
	// notify begining of LFN adapt
	if (CanCallback(AN_LFNADAPTSTART, ptdata->pfnNotifyProc, ptdata->nNotifyMask))
	{
		LFNADAPTINFO lai;
		lai.adblX = adblX;
		lai.pLFN = pNode;
		lai.dblError = dblError;
		lai.dblResponse = dblResponse;
		Callback(pALN, AN_LFNADAPTSTART, &lai, ptdata->pfnNotifyProc, ptdata->pvData);
	}
  // track split stats
  // below included LFN_CANSPLIT(pNode) &&
	if (bUsefulAdapt)
	{
		ASSERT(LFN_SPLIT(pNode) != NULL);
		LFN_SPLIT_COUNT(pNode)++;
		LFN_SPLIT_SQERR(pNode) += dblError * dblError * dblResponse;
		LFN_SPLIT_RESPTOTAL(pNode) += dblResponse;
		// there was a left brace here, moved to the end.

	// copy LFN vector pointers onto the stack for faster access
		int nDim = LFN_VDIM(pNode);
		ASSERT(nDim == pALN->nDim);
		double* adblW = LFN_W(pNode);	    // weight vector
		double* adblC = LFN_C(pNode);			// centroid vector
		double* adblD = LFN_D(pNode);			// average square dist from centroid vector
		// calculate adblA = how far the linear piece is above adblX[nOutput].
		// Note:  the sum below added to the bias weight would add up to zero for a point *on* the linear piece
		// If adblX[nOutput] is greater than the value of output on the piece, dblA is *negative*
		// dblA is also used later for computing convexity, which must be done w.r.t. the linear piece
		double dblA = *adblW++; // get weight adblW[0] into dblA then point to next weight
		//IMPORTANT: We talk about not losing numerical accuracy because we use centroids of linear pieces.
		// Here, we may lower accuracy by not using adblX[kk] - adblC[kk].  Another version should test this idea!
		for (int kk = 0; kk < nDim; kk++)
		{
			dblA += adblX[kk] * adblW[kk];
		}
		// We need to measure the error taking into account fillets.  The thickness of fillets is
		// dblError - dblA, positive when it is a MAX fillet.
		// If response is < 1, then the adjustment of the centroid and weights is lessened.
		// We need two learning rates.  The first is for the centroids and weights which divides by 2*nDim -1, thus
		// putting them on an equal footing with respect to correcting a share of the error.
		double dblLrnRate = ptdata->dblLearnRate;
		double dblLearnRespParam = dblLrnRate * dblResponse * region.dblLearnFactor / (2.0*nDim - 1.0);

		// ADAPT CENTROID FOR OUTPUT
		// Summing up:
		// L is the value of the linear piece at the input components of X
		// S is the level of the ALN function at the surface, including fillets
		// A = L - X
		// error = S - X
		// the thickness of the fillet (positive above the linear piece is S - L = error - A
		// the target level for the centroid C[output] is X - fillet thickness = X - error + A
		// We use a flywheel that averages that target level
		adblC[nOutput] += (adblX[nOutput] - dblError + dblA - adblC[nOutput]) * dblLearnRespParam;
		// Because the fillets may change, we aren't sure how much the dblError has changed.  We press ahead with
		// the old value of dblError while changing weights.  Similarly, we won't update the
		// dblError value as we change the weights.
		// Some helping variables
		double dblXmC = 0; // adblX[i] - adblC[i] "X minus C" for the current axis i
		double dblWeightFactor = 0;
		// ADAPT CENTROID AND WEIGHT FOR EACH INPUT VARIABLE  
		for (int i = 0; i < nDim; i++)
		{
			// get pointer to variable constraints
			ALNCONSTRAINT* pConstr = GetVarConstraint(NODE_REGION(pNode), pALN, i);
			// skip any variables of constant monotonicity; W is constant and X is irrelevant
			if (pConstr->dblWMax == pConstr->dblWMin) continue;
			// The following was changed on March 24, 2015, to allow for real-time
			// inputs where X[i] and the previous value might be correlated.
			// Compute the distance of X from the old centroid in axis i
			dblXmC = adblX[i] - adblC[i];
			// UPDATE VARIANCE BY EXPONENTIAL SMOOTHING
			// We adapt this first so the adaptation of the centroid to this input will not affect it
			ASSERT(adblD[i] >= 0);
			adblD[i] += (dblXmC * dblXmC - adblD[i]) * dblLrnRate; // Is this the right rate???
			// The exponentially smoothed estimate of variance
			// adblD[i] is not allowed to go below dblSqEpsilon of the current input variable to
			// slow rotations along some axes and prevent division by 0.
			if (adblD[i] < pConstr->dblSqEpsilon)adblD[i] = pConstr->dblSqEpsilon;
			// UPDATE THE CENTROID BY EXPONENTIAL SMOOTHING
			adblC[i] += dblXmC * dblLearnRespParam; // changed March 26 to correct mistake which led to terrible learning
			// Update the distance of X[i] from the centroid (the compiler will optimize these steps)
			dblXmC = adblX[i] - adblC[i];
			// ADAPT WEIGHTS
			dblWeightFactor = dblXmC / adblD[i];
			adblW[i] -= dblError * dblWeightFactor * dblLearnRespParam;
			// Bound the weight
			adblW[i] = max(min(pConstr->dblWMax, adblW[i]), pConstr->dblWMin);
			// COLLECT DATA FOR LATER SPLITTING THIS PIECE:
			// We exponentially smooth the error for points on the piece which are
			// further from the centroid than the variance stdev of the points on the piece along the current axis.
			// If the error is positive away from the centre, then we need a split of the LFN into a MIN node.
			if (LFN_CANSPLIT(pNode) && bUsefulAdapt && (dblXmC*dblXmC >= adblD[i]))
				LFN_SPLIT_T(pNode) += (dblError - LFN_SPLIT_T(pNode))* dblLrnRate;  //Is this the right rate???
		} // end loop over all nDim dimensions
		// compress the weighted centroid info into W[0]       
		double *pdblW0 = LFN_W(pNode);
		*pdblW0 = adblC[nDim - 1]; // there is no stored weight -1 for the output
		for (int i = 0; i < nDim - 1; i++)
		{
			*pdblW0 -= adblW[i] * adblC[i]; // here the W pointer is shifted up by one double
		}
		// notify end of LFN adapt
		if (CanCallback(AN_LFNADAPTEND, ptdata->pfnNotifyProc, ptdata->nNotifyMask))
		{
			LFNADAPTINFO lai;
			lai.adblX = adblX;
			lai.pLFN = pNode;
			lai.dblError = dblError;
			lai.dblResponse = dblResponse;
			Callback(pALN, AN_LFNADAPTEND, &lai, ptdata->pfnNotifyProc, ptdata->pvData);
		}
	}
}
