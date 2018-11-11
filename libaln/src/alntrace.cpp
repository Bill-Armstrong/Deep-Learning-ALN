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

// alntrace.cpp


#ifdef ALNDLL
#define ALNIMP __declspec(dllexport)
#endif

#include <aln.h>
#include "alnpriv.h"
#include <stdarg.h>

#ifdef _WIN32
#define WIN32_EXTRA_LEAN
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#endif

#ifdef _DEBUG   // entire file

// NOTE: in separate module so it can replaced if needed

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

/////////////////////////////////////////////////////////////////////////////
// Diagnostic Trace

void __cdecl ALNTrace(const char* pszFormat, ...)
{
	va_list args;
	va_start(args, pszFormat);

	int nBuf;
	char szBuffer[512];

	nBuf = vsprintf(szBuffer, pszFormat, args);
	ASSERT(nBuf < sizeof(szBuffer));

	va_end(args);

  #ifdef _WIN32
  OutputDebugString(szBuffer);
  #else
  perror(szBuffer);
  #endif
}

#endif //_DEBUG


