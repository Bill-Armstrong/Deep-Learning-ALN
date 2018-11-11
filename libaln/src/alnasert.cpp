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

// alnasert.cpp

#ifdef ALNDLL
#define ALNIMP __declspec(dllexport)
#endif

#include <aln.h>
#include "alnpriv.h"

#ifdef _DEBUG   // entire file

#ifdef _WIN32   // entire file
#define WIN32_EXTRA_LEAN
#define WIN32_LEAN_AND_MEAN
#include <windows.h>

// NOTE: in separate module so it can replaced if needed

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

LONG alnAssertBusy = -1;

int __stdcall ALNAssertFailed(const char* pszFileName, int nLine)
{
  char szMessage[512];
  
  // format message into buffer
  sprintf(szMessage, "File %hs, Line %d", pszFileName, nLine);
  
  // assume the debugger or auxiliary port
  char szT[512 + 20];
  sprintf(szT, "Assertion Failed: %s\n", szMessage);
  
  OutputDebugString(szT);
  
  if (InterlockedIncrement(&alnAssertBusy) > 0)
  {
    InterlockedDecrement(&alnAssertBusy);
    
    // assert within assert (examine call stack to determine first one)
    ALNDebugBreak();
    return FALSE;
  }
  
  // active popup window for the current thread  HOW CAN THIS BE IMPLEMENTED ?  WWA
  //HWND hWndParent = GetActiveWindow();
  //if (hWndParent != NULL)
  //hWndParent = GetLastActivePopup(hWndParent);
  
  // display the assert
  //int nCode = ::MessageBox(hWndParent, szMessage, "Assertion Failed!",
  //  MB_TASKMODAL|MB_ICONHAND|MB_ABORTRETRYIGNORE|MB_SETFOREGROUND);

  // cleanup
  InterlockedDecrement(&alnAssertBusy);

  //if (nCode == IDIGNORE)
    return FALSE;   // ignore
  
  //if (nCode == IDRETRY)
  //  return TRUE;    // will cause ALNDebugBreak

  ALNAbort();     // should not return (but otherwise ALNDebugBreak)
  return TRUE;
}

#endif  // _WIN32

#endif  // _DEBUG
