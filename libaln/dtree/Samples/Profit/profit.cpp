/*
// profit.c
// sample Dendronic Decisions Dtree application
// ALN Library
// Copyright (C) Dendronic Decisions Limited 1995 - 2007.
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
// Dendronic Decisions Limited
// 3624 - 108 Street NW
// Edmonton, Alberta, Canada  T6J 1B4

*/


#include <stdio.h>
#include <stdlib.h>
#include <dtree.h>

#define NUMPOINTS  11 

void main()
{    
  long lVersion;          /* DTREE library version */
  DTREE* pDtree;          /* pointer to DTREE structure */
  int nErrCode;           /* library function return value */
  char szErrMsg[256];     /* error message */
  int i;                  /* loop counter */
  double dblMin;          /* variable min */
  double dblMax;          /* variable max */
  double dblStep;         /* step size between quantization levels */

  /* get dtree version */
  lVersion = GetDtreeVersion();
  printf("DTREE library v%d.%d\n", lVersion >> 16, lVersion & 0x0000FFFF);
                        
  /* load file "profit.dtr" from Dendronic Learning Engine profit example */
  printf("Opening \"profit.dtr\"... ");
  nErrCode = ReadDtree("profit.dtr", &pDtree);

  /* check error return */  
  if (nErrCode == DTR_NOERROR)
  {
    printf("succesfully parsed!\n\n");
  }
  else
  {
    GetDtreeError(nErrCode, szErrMsg, sizeof(szErrMsg));
    printf("\nError (%d): %s\n", dtree_lineno, szErrMsg);
    return;
  }    


  /* succesfully loaded... evaluate on some test points */                
  
  /* make sure that the output variable index is 1 */
  if (pDtree->nOutputIndex != 1)
  {
    printf("Fatal error: Expected output variable to have index 1!\n");
    DestroyDtree(pDtree);
    return;
  }

  /* get range of first variable (advertising expense) */
  dblMax = pDtree->aVarDefs[0].bound.dblMax;
  dblMin = pDtree->aVarDefs[0].bound.dblMin;

  /* set step size between each point */
  dblStep = (dblMax - dblMin) / (NUMPOINTS - 1);
  
  /* output some header info */
  printf("Ad. (x$1000)\t\tProfit ($)\n");
  
  /* iterate over each test point */
  for (i = 0; i < NUMPOINTS; i++)
  {
    double dblProfit;       /* result of Dtree evaluation */
    double adblX[2];        /* 2 element evaluation vector */
                            /* 2nd element never used in this example */

    /* set evaluation point (ad expense) */
    /* ...no need to provide input value for output variable */
    adblX[0] = dblMin + (i * dblStep);  
    adblX[1] = 0; /* output variable */
    
          
    /* get dtree evaluation (profit) */
    if ((nErrCode = EvalDtree(pDtree, adblX, &dblProfit, NULL)) == DTR_NOERROR)
    {
      printf("%f\t\t%f\n", adblX[0], dblProfit);
    }
    else  /* error! */
    {
      GetDtreeError(nErrCode, szErrMsg, sizeof(szErrMsg));
      printf("\nError (%d): %s\n", dtree_lineno, szErrMsg);
    }
  }
  
  /* destroy dtree */
  DestroyDtree(pDtree);
}
