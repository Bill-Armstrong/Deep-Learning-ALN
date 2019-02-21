// OptionsDialog.cpp : implementation file
//
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


#include "stdafx.h"
#include "ALNfitDeep.h"
#include "ALNfitDeepDoc.h"
#include "ALNfitDeepView.h"
#include "OptionsDialog.h"
#include "DataFile.h"
#include "alnextern.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

/////////////////////////////////////////////////////////////////////////////
// COptionsDialog dialog


COptionsDialog::COptionsDialog(CWnd* pParent /*=NULL*/)
	: CDialog(COptionsDialog::IDD, pParent) , m_strSetTolerance(_T(""))
{
	//{{AFX_DATA_INIT(COptionsDialog)
	m_nPercentForTest = 10;
	m_nEstimateRMSError = 0;          // 0 means TRUE, compare half of remaining data with values from an overfitted ALN on the the other half
	m_strSetTolerance = _T("9999.9"); // this is not used now, but may be used eg for the probability in the F-test
	m_nZeroSmoothing = 0;             // 0 means TRUE, no smoothing is applied during approximation (it is never applied during least squares or noise level estimation)
	m_nNoJitter = 0;                  // 0 means TRUE, i.e we don't add jitter
	m_nDTREEDepth = 1;                // If greater than 1, parts of the input space are repeatedly split into two boxes to this depth, with ALNs on the blocks
	m_bTimePrefixes = TRUE;           // the time prefixed are prefixed to the output files, so all the files of a run stay together
	m_bReplaceUndefined = FALSE;      // set to true if we have to replace undefined values
	m_bDiagnostics = FALSE;           // if TRUE, diagnostics are printed
	//}}AFX_DATA_INIT
}


void COptionsDialog::DoDataExchange(CDataExchange* pDX)
{
  CDialog::DoDataExchange(pDX);
  //{{AFX_DATA_MAP(COptionsDialog)
  DDX_Text(pDX, IDC_EDITPERCENTFORTEST, m_nPercentForTest);
  DDX_Radio(pDX, IDC_RADIOVALUSEDATA, m_nEstimateRMSError);
  DDX_Text(pDX, IDC_EDIT_TOLERANCE, m_strSetTolerance);
  DDX_Radio(pDX, IDC_RADIO_ZEROSMOOTHING, m_nZeroSmoothing);
	DDX_Radio(pDX, IDC_RADIONOJITTER, m_nNoJitter);
	DDX_Text(pDX, IDC_EDITDEPTH, m_nDTREEDepth);
  DDX_Check(pDX, IDC_CHECK_TIME_PREFIXES, m_bTimePrefixes);
	DDX_Check(pDX, IDC_CHECK_REPLACE_DATA, m_bReplaceUndefined);
	DDX_Check(pDX, IDC_CHECK_DIAGNOSTICS, m_bDiagnostics);
  DDV_MinMaxInt(pDX, m_nDTREEDepth, 1, 25);
	DDV_MaxChars(pDX, m_strSetTolerance,20);
	//}}AFX_DATA_MAP
}

BEGIN_MESSAGE_MAP(COptionsDialog, CDialog)
	//{{AFX_MSG_MAP(COptionsDialog)
	//}}AFX_MSG_MAP
	ON_BN_CLICKED(IDC_RADIOVALUSEDATA, &COptionsDialog::OnBnClickedRadiovalusedata)
	ON_BN_CLICKED(IDC_RADIO_ZEROSMOOTHING, &COptionsDialog::OnBnClickedRadioZerosmoothing)
END_MESSAGE_MAP()


BOOL COptionsDialog::OnInitDialog() 
{
  return CDialog::OnInitDialog(); 
}



void COptionsDialog::OnBnClickedRadiovalusedata()
{
	// TODO: Add your control notification handler code here
}


void COptionsDialog::OnBnClickedRadioZerosmoothing()
{
	// TODO: Add your control notification handler code here
}
