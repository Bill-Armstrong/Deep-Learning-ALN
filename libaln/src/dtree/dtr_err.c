// dtr_err.c
// ALNfit Learning Engine for approximation of functions defined by samples.
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

#ifdef DTREEDLL
#define DTRIMP __declspec(dllexport)
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
                  
#include <dtree.h>
#include "dtr_priv.h"

typedef struct tagDTR_ERR
{ 
  int nErrCode;
  char* szErr;
} DTR_ERR;                              
                              
DTR_ERR _DTRErr[] =  
{                                                 
  DTR_GENERIC,                "unspecified error (NULL pointer???)",
  DTR_FILEERR,                "file access error",
  DTR_FILEWRITEERR,           "file write failure",
  DTR_FILEREADERR,            "file read failure",
  DTR_FILETYPEERR,            "not a valid DTREE binary file",
  DTR_ENDIANERR,              "binary file endian mismatch",
  DTR_MALLOCFAILED,           "memory allocation error", 
                                     
  DTR_BADVERSIONDEF,          "bad version defintion statement", 
  DTR_BADVERSIONINT,          "bad version number", 
  DTR_MISSINGVERSIONSEMI,     "missing ';' at end of version statement", 
  DTR_UNKNOWNVERSION,         "version of dtree file is newer then dtree library", 
                                                             
  DTR_BADDIMDEF,              "bad variable dimension statement", 
  DTR_BADDIMINT,              "bad variable dimension", 
  DTR_BADDIMRANGE,            "variable dimension must be greater than 1", 
  DTR_MISSINGDIMSEMI,         "missing ';' at end of variable dimension statement", 
                        
  DTR_BADVARDEF,              "bad variable definition statement", 
  DTR_BADVARDEFIDENT,         "bad variable identifier", 
  DTR_MISSINGVARDEFCOLON,     "missing colon after variable identifier", 
  DTR_MISSINGVARBOUNDSTART,   "missing '[' to start variable bound", 
  DTR_MISSINGVARBOUNDCOMMA,   "missing ',' after variable bound minimum", 
  DTR_BADVARBOUND,            "bad variable bound", 
  DTR_MISSINGVARBOUNDEND,     "missing ']' to end variable bound", 
  DTR_NEGATIVEVARBOUNDRANGE,  "variable bound maximum must be greater than minimum", 
  DTR_MISSINGVARDEFSEMI,      "missing ';' at end of variable definition statement", 
                                                             
  DTR_BADOUTPUTDEF,           "bad output variable definition statement", 
  DTR_BADOUTPUTIDENT,         "bad output variable identifier", 
  DTR_BADOUTPUTRANGE,         "output variable index out of range", 
  DTR_ZEROOUTPUTWEIGHT,       "output variable weight cannot be zero",
  DTR_MISSINGOUTPUTSEMI,      "missing ';' at end of output variable statement", 
                        
  DTR_UNKNOWNVARIDENT,        "undeclared variable identifier", 
  DTR_UNDEFINEDREGION,        "undefined region index", 
  DTR_UNDEFINEDLINEARFORM,    "undefined linear form index", 
                        
  DTR_BADBLOCKDEF,            "bad block definition statement", 
  DTR_BADBLOCKINT,            "bad block count", 
  DTR_MISSINGBLOCKINDEX,      "expected block index", 
  DTR_DUPBLOCKINDEX,          "duplicate block index", 
  DTR_MISSINGBLOCKCOLON,      "missing ':' after block index", 
  DTR_BADMINMAXNODE,          "bad min/max tree node definition",  
  DTR_MISSINGMINMAXLISTSTART, "missing '[' to start min/max child list", 
  DTR_EMPTYMINMAXLIST,        "cannot have an empty min/max child list", 
  DTR_BADMINMAXLISTELEM,      "bad min/max child list element",        
  DTR_BADBLOCKINDEXRANGE,     "block index out of range",
  DTR_MISSINGBLOCKSEMI,       "missing ';' at end of block statement", 
                        
  DTR_BADLINEARDEF,           "bad linear form definition statement", 
  DTR_BADLINEARINT,           "bad linear form count", 
  DTR_MISSINGLINEARINDEX,     "expected linear form index", 
  DTR_DUPLINEARINDEX,         "duplicate linear form index", 
  DTR_MISSINGLINEARCOLON,     "missing ':' after linear form index", 
  DTR_BADLINEARWEIGHT,        "bad linear form weight", 
  DTR_MISSINGCENTROIDSTART,   "missing '(' after weight",
  DTR_MISSINGCENTROIDEND,     "missing ')' after centroid",
  DTR_MISSINGCENTROIDIDENT,   "missing variable name after '('",
  DTR_MISSINGLINEARSIGN,      "missing sign in linear form",
  DTR_DUPVARCENTROID,         "duplicate variable weight/centroid definition",
  DTR_BADLINEARCENTROID,      "bad linear form centroid", 
  DTR_BADLINEARINDEXRANGE,    "linear index out of range",
  DTR_MISSINGLINEARSEMI,      "missing ';' at end of linear form statement", 
                                                             
  DTR_BADDTREEDEF,            "bad dtree definition statement", 
  DTR_BADDTREEINT,            "bad dtree node count", 
  DTR_MISSINGDTREEINDEX,      "expected dtree node index", 
  DTR_DUPDTREEINDEX,          "duplicate dtree node index", 
  DTR_MISSINGDTREECOLON,      "missing ':' after dtree node index identifier", 
  DTR_MISSINGDTREEVAR,        "missing variable name in dtree statement",
  DTR_MISSINGDTREELE,         "missing '<=' after variable name in dtree statement",
  DTR_MISSINGDTREETHRESHOLD,  "missing threshold after '<=' in dtree statement",
  DTR_BADTHRESHOLDRANGE,      "threshold out of range in dtree statement",
  DTR_MISSINGDTREEPAREN,      "missing ')' after threshold in dtree statement",
  DTR_MISSINGDTREEQUESTION,   "missing '?' after ')' in dtree statement",
  DTR_MISSINGDTREECHILDINDEX, "missing child dtree index in dtree statement",
  DTR_MISSINGDTREECHILDSEP,   "missing ':' after left child index in dtree statement",
  DTR_CIRCULARDTREECHILDREF,  "child dtree index must be greater than parent index",
  DTR_BADDTREEINDEXRANGE,     "dtree index out of range",
  DTR_MISSINGDTREEBLOCKINDEX, "missing block index after 'block' in dtree statement",
  DTR_MISSINGDTREESEMI,       "missing ';' at end of dtree statement", 

  0,                          "Unknown error code", 
};

void GetErrMsg(int nErrno, char* pBuf, int nMaxBufLen)
{
  int i;
  if (nErrno == 0)
  { 
    static char _szNoErr[] = "No error";
    strncpy(pBuf, _szNoErr, nMaxBufLen);
    pBuf[nMaxBufLen - 1] = 0;
    return;
  }    
  for (i = 0; _DTRErr[i].nErrCode != 0 && _DTRErr[i].nErrCode != nErrno; i++);
  strncpy(pBuf, _DTRErr[i].szErr, nMaxBufLen);
  pBuf[nMaxBufLen - 1] = 0;
}
