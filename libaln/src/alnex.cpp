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

// alnex.cpp
// exception handling routines


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
// CALNException

CALNException::CALNException(BOOL bAutoDelete /*= TRUE*/, 
                             const char* pszReason /*= NULL*/)
{
  m_bAutoDelete = bAutoDelete;
	
  m_szReason[0] = '\0';	  // set NULL at beginning of string
	if (pszReason)
		strncpy(m_szReason, pszReason, sizeof(m_szReason) - 1);
  m_szReason[255] = '\0';	// set NULL at end of string
}
	
void CALNException::Delete()
{
  if (m_bAutoDelete)
    delete this;
}

void ALNAPI ThrowALNException()
{
  static CALNException gALNException(FALSE);
  throw &gALNException;
}
 
///////////////////////////////////////////////////////////////////////////////
// CALNMemoryException

CALNMemoryException::CALNMemoryException(BOOL bAutoDelete /*= TRUE*/)
  : CALNException(bAutoDelete)
{
}

void ALNAPI ThrowALNMemoryException()
{
  static CALNMemoryException gALNMemoryException(FALSE);
  throw &gALNMemoryException;
}

///////////////////////////////////////////////////////////////////////////////
// CALNUserException

CALNUserException::CALNUserException(BOOL bAutoDelete /*= TRUE*/)
  : CALNException(bAutoDelete)
{
}

void ALNAPI ThrowALNUserException()
{
  static CALNUserException gALNUserException(FALSE);
  throw &gALNUserException;
}
