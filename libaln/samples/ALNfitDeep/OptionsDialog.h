// OptionsDialog.h : header file
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

#if !defined(AFX_OPTIONSDIALOG_H__D4C60515_32E7_43E8_80F1_035BD2D8E8C1__INCLUDED_)
#define AFX_OPTIONSDIALOG_H__D4C60515_32E7_43E8_80F1_035BD2D8E8C1__INCLUDED_

#if _MSC_VER >= 1000
#pragma once
#endif // _MSC_VER >= 1000

// COptionsDialog dialog

class COptionsDialog : public CDialog
{
// Construction
public:
	COptionsDialog(CWnd* pParent = NULL);   // standard constructor

// Dialog Data
	//{{AFX_DATA(COptionsDialog)
	enum { IDD = IDD_OPTIONS };
	int		m_nPercentForTest;
	int		m_nEstimateRMSError;
  CString m_strSetTolerance;
	int		m_nZeroSmoothing;
	int		m_nNoJitter;
	int		m_nDTREEDepth;
	BOOL	m_bTimePrefixes;
	BOOL	m_bReplaceUndefined;
	BOOL	m_bDiagnostics;
	//}}AFX_DATA


// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(COptionsDialog)
	protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support
	//}}AFX_VIRTUAL

// Implementation
protected:
	// Generated message map functions
	//{{AFX_MSG(COptionsDialog)
	virtual BOOL OnInitDialog();
	//}}AFX_MSG

public:
  // This keeps the tolerance input able to be converted to a double
	void OnEnUpdateEditTolerance(void);
	DECLARE_MESSAGE_MAP()
	afx_msg void OnBnClickedRadiovalusedata();
	afx_msg void OnBnClickedRadioZerosmoothing();
};

//{{AFX_INSERT_LOCATION}}
// Microsoft Developer Studio will insert additional declarations immediately before the previous line.

#endif // !defined(AFX_OPTIONSDIALOG_H__D4C60515_32E7_43E8_80F1_035BD2D8E8C1__INCLUDED_)
