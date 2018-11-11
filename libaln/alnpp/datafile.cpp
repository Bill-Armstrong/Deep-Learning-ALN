// ALN Library C++ data file class
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

// datafile.cpp


#ifdef __GNUC__
#include <typeinfo>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <memory.h>
#include <ctype.h>

#include <datafile.h>

#ifdef sun
#include <floatingpoint.h>
#include <unistd.h>
#endif

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif


/////////////////////////////////////////////////////////////////////////////
// class CDataFile

CDataFile::CDataFile()
{
  m_pBuffer = NULL;
  m_lBufferLen = 0;
  m_lColumns = 0;
  m_lRows = 0;
}

CDataFile::CDataFile(const CDataFile& datafile)
{
  m_pBuffer = NULL;
  m_lBufferLen = 0;
  m_lColumns = 0;
  m_lRows = 0;

  *this = datafile;
}

// dtor
CDataFile::~CDataFile()
{               
  Destroy();
}

double CDataFile::GetColMax(long lCol) const
{
  ASSERT(lCol >=0 && lCol < m_lColumns);
  ASSERT(m_pBuffer);
               
  double* pdblBase = m_pBuffer + lCol;

  double dblMax = *pdblBase;               
  for(long i = 0; i < m_lRows; i++, pdblBase += m_lColumns)
  {
    double dbl = *pdblBase;               
    if (dbl > dblMax)
    {
      dblMax = dbl;
    }
  }
  
  return dblMax;
}

double CDataFile::GetColMin(long lCol) const
{
  ASSERT(lCol >=0 && lCol < m_lColumns);
  ASSERT(m_pBuffer);
  
  double* pdblBase = m_pBuffer + lCol;

  double dblMin = m_pBuffer[lCol];               
  for(long i = 0; i < m_lRows; i++, pdblBase += m_lColumns)
  {
    double dbl = *pdblBase;               
    if (dbl < dblMin)
    {
      dblMin = dbl;
    }
  }
  
  return dblMin;
}
         
BOOL CDataFile::Create(long lRows, long lColumns)
{                    
  Destroy();
  ASSERT(m_pBuffer == NULL);
  
  m_lColumns = lColumns;
  m_lRows = lRows;
   
  // allocate array mem         
  long lElements = m_lRows * m_lColumns * 2;
  if (lElements > 0)
  {
    if (!Grow(lElements * sizeof(double)))
    {
      m_lColumns = 0;
      m_lRows = 0;

      return FALSE;
    }

    memset(m_pBuffer, 0, m_lBufferLen);
  }

  return TRUE;
}

void CDataFile::Destroy()
{   
  if (m_pBuffer != NULL)
  {
    free(m_pBuffer);
    m_pBuffer = NULL;
  }
  
  m_lBufferLen = 0;
  m_lColumns = 0;
  m_lRows = 0;
}

CDataFile& CDataFile::operator = (const CDataFile& datafile)
{
  // don't copy ourself!
  if (&datafile == this)
    return *this;
  
  // try to grow ourself
  if (!Grow(datafile.m_lBufferLen))
  {
    Destroy();
    return *this; // failed to grow file... no error return!
  }
  
  m_lRows = datafile.m_lRows;
  m_lColumns = datafile.m_lColumns;
  memcpy(m_pBuffer, datafile.m_pBuffer, m_lBufferLen);

  return *this;
}

BOOL CDataFile::Grow(long lNewLen)
{
  static const long lGrowBytes = 1024;

	if (lNewLen > m_lBufferLen)
	{
		// grow the buffer
		long lNewBufferSize = (long)m_lBufferLen;

		// watch out for buffers which cannot be grown!
		ASSERT(lGrowBytes != 0);

		// determine new buffer size
		while (lNewBufferSize < lNewLen)
			lNewBufferSize += lGrowBytes;

		// allocate new buffer
		BYTE* pNew;
		if (m_pBuffer == NULL)
			pNew = (BYTE*)malloc(lNewBufferSize);
		else
			pNew = (BYTE*)realloc(m_pBuffer, lNewBufferSize);

		if (pNew == NULL)
			return FALSE;

		m_pBuffer = (double*)pNew;
		m_lBufferLen = lNewBufferSize;
	}

  return TRUE;
}
 
static BOOL SimpleFloatParse(char* pszText, double& d)
{
	ASSERT(pszText != NULL);
	while (*pszText == ' ' || *pszText == '\t')
		pszText++;

	char chFirst = pszText[0];
	d = strtod(pszText, &pszText);
	if (d == 0.0 && chFirst != '0')
		return FALSE;   // could not convert
	while (*pszText == ' ' || *pszText == '\t')
		pszText++;

	if (*pszText != '\0')
		return FALSE;   // not terminated properly

	return TRUE;
}

BOOL CDataFile::Append(const CDataFile& datafile)
{
  // grow memory to desired length
  long lLength = sizeof(double) * ((datafile.m_lRows + m_lRows) * m_lColumns);
  if (!Grow(lLength))
    return FALSE;

  // test if column count is same 
  if (m_lColumns == datafile.m_lColumns)
  {
    // just copy mem!
    long lDestPoints = m_lRows * m_lColumns;
    long lSrcPoints = datafile.m_lRows * m_lColumns;
    memcpy(m_pBuffer + lDestPoints, 
           datafile.m_pBuffer, 
           lSrcPoints * sizeof(double));
  }
  else
  {
    // test if src has less cols than dest
    if (m_lColumns > datafile.m_lColumns)
    {
      // zero out mem
      long lDestPoints = m_lRows * m_lColumns;
      long lSrcPoints = datafile.m_lRows * m_lColumns;
      memset(m_pBuffer + lDestPoints, 
             0, 
             lSrcPoints * sizeof(double));
    }

    // truncate number of copied columns if necessary
    long lMaxCol = (datafile.m_lColumns < m_lColumns) ?
                   datafile.m_lColumns : 
                   m_lColumns;

    long i = 0;        // row index of source
    long w = m_lRows;  // row index of dest
    for (; i < datafile.m_lRows; i++, w++)
    {
      register long lSrcOffset = i * datafile.m_lColumns;
      register long lDestOffset = w * m_lColumns;
      for (long j = 0; j < lMaxCol; j++)
      {
        m_pBuffer[lDestOffset + j] = datafile.m_pBuffer[lSrcOffset + j];
      }
    }
  }

  // set new row count
  m_lRows += datafile.m_lRows;
  
  return TRUE;
}
 
BOOL CDataFile::Read(const char* pszFileName)
{
  // clear existing data
  Destroy();

  // open file
  FILE* f = NULL;
  if (fopen_s(&f, pszFileName, "r") != 0)
    return FALSE;
  
  // alloc line buffer
  char* pBuf = new char[4096];
  if (pBuf == NULL)
    return FALSE;

  BOOL bFirstRow = TRUE;
  
  long lRow = 0;
  while (!feof(f))
  {   
    // read a line
    char* p = fgets(pBuf, 4096, f);
    if (p == NULL)
      break;

    long lCol = 0;
    
    // parse tokens
    char* pTok = NULL;
    while ((pTok = strtok(p, " ,\t\n")))
    {
      p = NULL;   // continue parsing from original string next 
                  //   time thru loop
     
      // detect comments (any punctuation except '-' and '.')
      if ((*pTok != '-') && (*pTok != '.') && ispunct(*pTok))
      {
        break;
      }

      // pTok is now null terminated, so we can parse out a double
      double dbl;
      if (!SimpleFloatParse(pTok, dbl))
      {
        // invalid input!
        delete[] pBuf;
        Destroy();
        return FALSE;
      }
     
      // calc current data block index
      long lIndex = (long)lRow * m_lColumns + lCol;
      
      // grow memory if necessary
      if (lIndex >= m_lBufferLen / (long)sizeof(double) &&
          !Grow((lIndex + 1) * sizeof(double)))
      {
        // failed allocation
        delete[] pBuf;
        Destroy();
        return FALSE;
      }

      // store data
      m_pBuffer[lIndex] = dbl;

      // increase column count
      lCol++;  
    }

    // increase row count
    lRow++;

    // finished parsing this row... check column count
    if (lCol == 0)
    {
      lRow--;  // no data on this row, fixup row count
    }
    else if (bFirstRow)
    {
      m_lColumns = lCol;
      bFirstRow = FALSE;
    }
    else if (lCol != m_lColumns)
    {
      // too many columns in this row!
      delete[] pBuf;
      Destroy();
      return FALSE;
    }
  }         

  m_lRows = lRow;

  fclose(f);  
  delete[] pBuf;

  return TRUE;
}

BOOL CDataFile::ReadBinary(const char* pszFileName)
{
  // open file
  FILE* f = NULL;
  if (fopen_s(&f, pszFileName, "rb") != 0)
    return FALSE;

  unsigned int nDataPoints;
  long lRows;
  long lColumns;

  // check signature
  char szSig[4] = "";
  static unsigned int nSigLen = sizeof(szSig);
  if (fread(szSig, 1, nSigLen, f) != nSigLen)
    goto error;
  else if (strcmp(szSig, "CDF") != 0)
    goto error;

  // read rows, cols
  if (fread(&lRows, sizeof(lRows), 1, f) != 1)
    goto error;

  if (fread(&lColumns, sizeof(lColumns), 1, f) != 1)
    goto error;
  
  // create space for data
  if (!Create(lRows, lColumns))
    goto error;

  // read data
  nDataPoints = m_lRows * m_lColumns;
  if (fread(m_pBuffer, sizeof(double), nDataPoints, f) != nDataPoints)
    goto error;

  // success
  fclose(f);
  return TRUE;

error:
  fclose(f);
  return FALSE;
}

BOOL CDataFile::ReadAppend(const char* pszFileName)
{
  CDataFile datafile;
  if (!datafile.Read(pszFileName))
    return FALSE;

  return Append(datafile);
}

BOOL CDataFile::ReadAppendBinary(const char* pszFileName)
{
  CDataFile datafile;
  if (!datafile.ReadBinary(pszFileName))
    return FALSE;

  return Append(datafile);
}

BOOL CDataFile::Write(const char* pszFileName, int nDelim /*= '\t'*/)
{
  // open file
  FILE* f = NULL;
  if (fopen_s(&f, pszFileName, "w") != 0)
    return FALSE;
  
  for (long i = 0; i < m_lRows; i++)
  {
    for(long j = 0; j < m_lColumns; j++)
    {
      long lIndex = i * m_lColumns + j;
      double dbl = m_pBuffer[lIndex];     

      // value
      fprintf(f, "%0.19g", dbl);
     
      // delim
      if (j < (m_lColumns - 1))
        fputc(nDelim, f);
    } 

    // new row
    fprintf(f, "\n");
  }
  
  fclose(f);
  return TRUE;
}

BOOL CDataFile::WriteBinary(const char* pszFileName)
{
  // open file
  FILE* f = NULL;
  if (fopen_s(&f, pszFileName, "wb") != 0)
    return FALSE;
 
  unsigned int nDataPoints;

  static char szSig[] = "CDF";
  static unsigned int nSigLen = sizeof(szSig);
  if (fwrite(szSig, 1, nSigLen, f) != nSigLen)
    goto error;

  // write rows, cols
  if (fwrite(&m_lRows, sizeof(m_lRows), 1, f) != 1)
    goto error;

  if (fwrite(&m_lColumns, sizeof(m_lColumns), 1, f) != 1)
    goto error;
  
  // write data
  nDataPoints = m_lRows * m_lColumns;
  if (fwrite(m_pBuffer, sizeof(double), nDataPoints, f) != nDataPoints)
    goto error;

  // success
  fclose(f);
  return TRUE;

error:
  fclose(f);
  remove(pszFileName);
  return FALSE;
}

