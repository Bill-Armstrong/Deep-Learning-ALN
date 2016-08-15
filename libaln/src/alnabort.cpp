// ALN Library
// Copyright (C) 1995 - 2010 William W. Armstrong.
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

// 3624 - 108 Street NW
// Edmonton, Alberta, Canada  T6J 1B4
// alnabort.cpp

///////////////////////////////////////////////////////////////////////////////
//  File version info:
// 
//  $Archive: /ALN Development/libaln/src/alnabort.cpp $
//  $Workfile: alnabort.cpp $
//  $Revision: 5 $
//  $Date: 7/17/07 5:08p $
//  $Author: Arms $
//
///////////////////////////////////////////////////////////////////////////////

#ifdef ALNDLL
#define ALNIMP __declspec(dllexport)
#endif

#include <aln.h>
#include "alnpriv.h"

// NOTE: in separate module so it can replaced if needed

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

// abort proc address
static ALNABORTPROC _pfnAbort = NULL;

// set ALN abort procedure: if ALN needs to abort, it will call this
//   function; if not set, ALN will simply call abort()
// you can use this to chain an abort proc for another library to the ALN
//  library, and simply call ALNAbort
ALNIMP void ALNAPI ALNSetAbortProc(ALNABORTPROC pfnAbortProc)
{
	_pfnAbort = pfnAbortProc;
}

// ALN library abort... cleans up internal ALN library, then calls
// abort procedure defined by call to ALNSetAbortProc()
ALNIMP void __stdcall ALNAbort(void)
{
	// any ALN lib cleanup here

	if (!_pfnAbort || _pfnAbort == ALNAbort)
		abort();

	else
		(*_pfnAbort)();
}
