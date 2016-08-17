// ALN Library (libaln)
// file split_ops.cpp
// 
// Copyright (C) 1995 - 2010 William W. Armstrong.
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






#ifdef ALNDLL
#define ALNIMP __declspec(dllexport)
#endif

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

// include files
#include <stdafx.h>
#include <alnpp.h>
#include <datafile.h>


// include classes
#include ".\cmyaln.h" 

extern double dblMax;
extern int nDim;
extern double dblSetTolerance;
extern long nRowsTR;
extern long nRowsALNinputValFile;
extern BOOL bEstimateRMSError;
extern CDataFile TRfile;
extern CDataFile ALNinputValFile;

void dozerospliterror(CMyAln* pALN, ALNNODE* pNode);
void splitcontrol(CMyAln* pALN, double dblLimit);
void dozerospliterror(CMyAln* pALN, ALNNODE* pNode);
void dodivideTR(CMyAln* pALN, ALNNODE* pNode);
void dodivideVL(CMyAln* pALN, ALNNODE* pNode);
void dosplitcontrol(CMyAln* pALN, ALNNODE* pNode, double dblLimit);
void spliterrorsetTR(CMyAln * pALN);
void spliterrorsetVL(CMyAln * pALN);



void splitcontrol(CMyAln* pALN, double dblLimit)  // routine
{
  ASSERT(pALN);
  ASSERT(pALN->GetTree());
	// initialize the SPLIT components to zero
	dozerospliterror(pALN, pALN->GetTree());
	// get square errors of pieces on training set
	spliterrorsetTR(pALN);
	// divide the training errors of the pieces by the hit counts, set counts to zero
	dodivideTR(pALN,pALN->GetTree());
	// get square errors of the pieces on the validation set
	if(bEstimateRMSError) 
	{
		// there is a validation set
		spliterrorsetVL(pALN);
		// divide the validation errors of the pieces by the hit count
		dodivideVL(pALN,pALN->GetTree());
	}
  dosplitcontrol(pALN, pALN->GetTree(), dblLimit);
  // reset the SPLIT components to zero
	dozerospliterror(pALN, pALN->GetTree());
}

// We use the first three fields in ALNLFNSPLIT (declared in aln.h)
// in two different ways: in between intervals in which the hyperplanes
// in leaf nodes change weights, there is an analysis of the errors of
// each piece on the training and validation sets.  If the training
// error is less than a specified fraction of the validation error,
// the piece is not split.
#define DBLSQERRORVAL dblRespTotal
void dozerospliterror(CMyAln* pALN, ALNNODE* pNode) // routine
{
	// initializes counters
  ASSERT(pNode);
  if (NODE_ISMINMAX(pNode))
  {
    dozerospliterror(pALN, MINMAX_LEFT(pNode));
    dozerospliterror(pALN, MINMAX_RIGHT(pNode));
  }
  else
  {
		ASSERT(NODE_ISLFN(pNode));
		(pNode->DATA.LFN.pSplit)->nCount = 0;
    (pNode->DATA.LFN.pSplit)->dblSqError = 0;
    (pNode->DATA.LFN.pSplit)->DBLSQERRORVAL = 0;		
	}
}

void dodivideTR(CMyAln* pALN, ALNNODE* pNode) // routine
{
	// divides the total square errors of the pieces by their hit count
  ASSERT(pNode);
  if (NODE_ISMINMAX(pNode))
  {
    dodivideTR(pALN, MINMAX_LEFT(pNode));
    dodivideTR(pALN, MINMAX_RIGHT(pNode));
  }
  else
  {
    ASSERT(NODE_ISLFN(pNode));
		if((pNode->DATA.LFN.pSplit)->nCount > 2 * nDim)  // factor two inserted WWA 2009.10.06
		{
			(pNode->DATA.LFN.pSplit)->dblSqError /=
										 (pNode->DATA.LFN.pSplit)->nCount;
			(pNode->DATA.LFN.pSplit)->nCount =0;  // inserted WWA 2009.10.06 IMPORTANT ERROR WAS HERE
			// this is zeroed in order to use nCount for spliterrorsetVL

		}
		else
		{
      // if the count of training samples hitting the piece is 2* nDim or less, we don't want to split
			// so we make the training error of the piece zero
			(pNode->DATA.LFN.pSplit)->dblSqError = 0;
    }
	}
}


void dodivideVL(CMyAln* pALN, ALNNODE* pNode) // routine
{
	// divides the total square errors of the pieces by their hit count
  ASSERT(pNode);
  if (NODE_ISMINMAX(pNode))
  {
    dodivideVL(pALN, MINMAX_LEFT(pNode));
    dodivideVL(pALN, MINMAX_RIGHT(pNode));
  }
  else
  {
    ASSERT(NODE_ISLFN(pNode));
		if((pNode->DATA.LFN.pSplit)->nCount > nDim)
		{
      (pNode->DATA.LFN.pSplit)->DBLSQERRORVAL /=
										 (pNode->DATA.LFN.pSplit)->nCount;
		}
    else
    {
      // if there are too few or no validation samples, make sure we don't divide
      (pNode->DATA.LFN.pSplit)->DBLSQERRORVAL = dblMax;
    } 
	}
}

void dosplitcontrol(CMyAln* pALN, ALNNODE* pNode, double dblLimit) // routine
{
	// this routine visits all the leaf nodes and determines whether the
	// training error is below a certain limit related to validation error or tolerance
	double dblSqErrorPieceTrain;
	double dblSqErrorPieceVal;
  ASSERT(pNode);
  if (NODE_ISMINMAX(pNode))
  {
    dosplitcontrol(pALN, MINMAX_LEFT(pNode), dblLimit);
    dosplitcontrol(pALN, MINMAX_RIGHT(pNode), dblLimit);
  }
  else
  {
    ASSERT(NODE_ISLFN(pNode));
    dblSqErrorPieceTrain = (pNode->DATA.LFN.pSplit)->dblSqError; // average square error on TR set
    if(bEstimateRMSError)
    {
			dblSqErrorPieceVal = (pNode->DATA.LFN.pSplit)->DBLSQERRORVAL;
    }
    else
    {
			dblSqErrorPieceVal = dblSetTolerance * dblSetTolerance; // no validation set!
    }
		if(dblSqErrorPieceTrain * (dblLimit * dblLimit) < dblSqErrorPieceVal)
		// if the training average square error of the piece times dblLimit^2 (> 1.0) is less than the square validation error 
		// or, if we are skipping validation,
		// the training average square error times dblLimit is less than tolerance squared
		{
			// stop all future splitting of the piece
       LFN_FLAGS(pNode) &= ~LF_SPLIT;  // ~ is the unary one's complement operator
		}
	}
}


void spliterrorsetTR(CMyAln * pALN) // routine
{
	// assign the square errors on the validation set to the leaf nodes of the ALN
	double * adblX = (double *) malloc((nDim) * sizeof(double));
	double desired = 0;
	double predict = 0;
  double      se = 0; // square error accumulator
	ALNNODE* pActiveLFN;
	for(int j=0; j<nRowsTR; j++)
	{
		for(int i=0; i<nDim; i++)
		{
			adblX[i] = TRfile.GetAt(j,i,0);
		}
		desired = adblX[nDim - 1]; // get the desired result
		adblX[nDim-1] = 0; // not used in evaluation by QuickEval
		predict = pALN->QuickEval(adblX, &pActiveLFN);
    se = (predict - desired) * (predict - desired);
		(pActiveLFN->DATA.LFN.pSplit)->nCount++;
		(pActiveLFN->DATA.LFN.pSplit)->dblSqError += se;
  } // end loop over TRset
	free(adblX);
} // END of spliterrorsetTR


void spliterrorsetVL(CMyAln * pALN) // routine
{
	// assign the square errors on the validation set to the leaf nodes of the ALN
	double * adblX = (double *) malloc((nDim) * sizeof(double));
	double desired = 0;
	double predict = 0;
  double      se = 0; // square error added to LFN DBLSQERRORVAL
	ALNNODE* pActiveLFN;
	for(int j=0; j<nRowsALNinputValFile; j++)   // this is expensive using the whole validation set, but more accurate
  {
    for(int i=0; i<nDim; i++)
		{
			adblX[i] = ALNinputValFile.GetAt(j,i,0); // the value at nDim - 1 is used only for desired
		}
		desired = adblX[nDim - 1];
		adblX[nDim - 1] = 0; // not used by QuickEval WWA 2009.10.06
		predict = pALN->QuickEval(adblX, &pActiveLFN);
    se = (predict - desired) * (predict - desired);
		(pActiveLFN->DATA.LFN.pSplit)->nCount++;
		(pActiveLFN->DATA.LFN.pSplit)->DBLSQERRORVAL += se;
  } // end loop over VLfile
	free(adblX);
} // END of spliterrorsetVL
