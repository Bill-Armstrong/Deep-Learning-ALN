/*
// dtr_priv.h
// DTREE implementation support routines
// Copyright (C) 1995 Dendronic Decisions Limited

// This source code is provided for education and evaluation purposes only.
// Use of this source code for any other purpose requires an Atree 3.0 
// Commercial License.  Please contact Dendronic Decisions Limited at
// +1 (403) 438 8285 for further details.
*/

#ifndef __DTR_PRIV_H__
#define __DTR_PRIV_H__
          
/*          
/////////////////////////////////////////////////////////////////////
// parser - returns DTE_NOERR on success                           
// returns new DTREE in *ppDtree - use DestroyDtree to delete
*/
int ParseDtreeFile(FILE* pFile, DTREE** ppDtree);
  
/*  
/////////////////////////////////////////////////////////////////////
// exporter - returns DTE_NOERR on success                           
// exports pDtree to a text file
*/
int ExportDtreeFile(FILE* pFile, const char* pszFileName, DTREE* pDtree);

/*          
/////////////////////////////////////////////////////////////////////
// binary file import - returns DTE_NOERR on success                           
// returns new DTREE in *ppDtree - use DestroyDtree to delete
*/
int ReadBinDtreeFile(FILE* pFile, DTREE** ppDtree);
  
/*  
/////////////////////////////////////////////////////////////////////
// binary file export - returns DTE_NOERR on success                           
// exports pDtree to a binary file
*/
int WriteBinDtreeFile(FILE* pFile, DTREE* pDtree);

  
/*                  
/////////////////////////////////////////////////////////////////////
// error code messages
// returns new DTREE in *ppDtree - use DestroyDtree to delete
*/
void GetErrMsg(int nErrno, char* pBuf, int nMaxBufLen);

                  
#endif  /* __DTR_PRIV_H__ */
