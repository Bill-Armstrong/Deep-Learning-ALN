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

// datafile.h


#ifndef __DATAFILE_H__
#define __DATAFILE_H__

#ifndef ASSERT
#ifndef ALNASSERT
#define ASSERT assert
#else
#define ASSERT ALNASSERT
#endif
#endif

#ifndef assert
#include <assert.h>
#endif

/////////////////////////////////////////////////////////////////////////////
// library files (Microsoft compilers only)

#ifndef ALNPP_NOFORCE_LIBS

#ifdef _MSC_VER

#ifdef _DEBUG
  #ifdef _MT
    #ifdef _DLL
    #pragma comment(lib, "libalnppdmtd.lib")   /* MT DLL */
    #else  /* !_DLL */
    #pragma comment(lib, "libalnppmtd.lib")    /* MT */
    #endif /* _DLL */
  #else
    #pragma comment(lib, "libalnppd.lib")
  #endif /* _MT */
#else
  #ifdef _MT
    #ifdef _DLL
    #pragma comment(lib, "libalnppdmt.lib")    /* MT DLL */
    #else  /* !_DLL */
    #pragma comment(lib, "libalnppmt.lib")     /* MT */
    #endif /* _DLL */
  #else
    #pragma comment(lib, "libalnpp.lib")
  #endif /* _MT */
#endif

#endif  // _MSC_VER

#endif  // ALNPP_NOFORCE_LIBS



/////////////////////////////////////////////////////////////////////////////
// data type definitions

typedef unsigned char BYTE;
typedef int BOOL;

#ifndef TRUE
#define TRUE 1
#endif

#ifndef FALSE
#define FALSE 0
#endif


///////////////////////////////////////////////////////////////////////////////
// class CDataFile

class CDataFile
{ 
// Constructors
public:
  CDataFile();    
  CDataFile(const CDataFile& datafile);
    
  BOOL Create(long lRows, long lColumns);

// Attributes
public:
  
  long ColumnCount() const
    { return m_lColumns; }
  long RowCount() const
    { return m_lRows; }
  
  double GetColMax(long lCol) const;
  double GetColMin(long lCol) const;
  
  double operator[](long lIndex) const
    {
      ASSERT(lIndex < (m_lBufferLen / (long)sizeof(double)));
      return m_pBuffer[lIndex];
    }
  double& operator[](long lIndex)
    {
      ASSERT(lIndex < (m_lBufferLen / (long)sizeof(double)));
      return m_pBuffer[lIndex];
    }

  long CalcDataIndex(long lRow, long lColumn, long lDelta = 0) const
    {
      ASSERT((lRow + lDelta) >= 0 && (lRow + lDelta) < m_lRows && 
             lColumn >=0 && lColumn < m_lColumns);
      ASSERT((lRow + lDelta) * m_lColumns + lColumn < (m_lBufferLen / (long)sizeof(double)));
      return (lRow + lDelta) * m_lColumns + lColumn;
    }

  const double* GetRowAt(long lRow) const
    {
      ASSERT(m_pBuffer);
      return m_pBuffer + CalcDataIndex(lRow, 0, 0);
    }

  double* GetRowAt(long lRow)
    {
      ASSERT(m_pBuffer);
      return m_pBuffer + CalcDataIndex(lRow, 0, 0);
    }

  double GetAt(long lRow, long lColumn, long lDelta = 0) const
    { 
      ASSERT(m_pBuffer);
      return m_pBuffer[CalcDataIndex(lRow, lColumn, lDelta)];
    }
  void SetAt(long lRow, long lColumn, double dbl, long lDelta = 0)
    { 
      ASSERT(m_pBuffer);
      m_pBuffer[CalcDataIndex(lRow, lColumn, lDelta)] = dbl;
    }

	const double* GetDataPtr() const
		{	return m_pBuffer; }

  double* GetDataPtr()
		{	return m_pBuffer; }

// Operations:
public:

  BOOL Append(const CDataFile& datafile);
    // appends datafile to end of this... truncates or adds columns
    // from datafile as necessary to match our columns

  BOOL Read(const char* pszFileName);
  BOOL ReadBinary(const char* pszFileName);
    // read data from a file, erase current contents
  
  BOOL ReadAppend(const char* pszFileName);
  BOOL ReadAppendBinary(const char* pszFileName);
    // read data from a file, appending to current contents
  
  BOOL Write(const char* pszFileName, int nDelim = '\t');
  BOOL WriteBinary(const char* pszFileName);
 
  void Destroy();
  
// implementation
public:
  virtual ~CDataFile();  
  CDataFile& operator = (const CDataFile& datafile);

protected:  

  // growing the data file
  BOOL Grow(long lNewLen);
  
  double* m_pBuffer;      // data block
  long m_lBufferLen;      // length of block
  long m_lColumns;        // number of columns
  long m_lRows;           // number of rows
};

///////////////////////////////////////////////////////////////////////////////

#endif  // __DATAFILE_H__
