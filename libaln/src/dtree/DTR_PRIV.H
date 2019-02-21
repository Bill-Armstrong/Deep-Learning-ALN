// dtr_priv.h
// DTREE implementation support routines

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
