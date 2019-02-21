// ALNfitDeepView.h : interface of the CALNfitDeepView class
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



#if !defined(AFX_ALNFITDEEPVIEW_H__5CD205CA_22D6_4675_88DD_FB3789D96D18__INCLUDED_)
#define AFX_ALNFITDEEPVIEW_H__5CD205CA_22D6_4675_88DD_FB3789D96D18__INCLUDED_

#if _MSC_VER >= 1000
#pragma once
#endif // _MSC_VER >= 1000

//User-defined Messages where they are going to be sent in the worker thread
#define WM_UPDATESCREEN WM_USER + 100
#define WM_THREADFINISHED WM_USER + 101

// set up globals from the worker thread procedure
void PassBackStatus(int nMes, int nPer);

// worker thread procedure
UINT ActionsProc(LPVOID pParam);

class CALNfitDeepView : public CFormView
{
protected: // create from serialization only
	CALNfitDeepView();
	DECLARE_DYNCREATE(CALNfitDeepView)

public:
	//{{AFX_DATA(CALNfitDeepView)
	enum { IDD = IDD_ALNFITDEEP_FORM };
	CProgressCtrl	m_Progress;
	CString	m_strDataFileName;
	CString	m_strDTREEFileName;
	int		m_nFit;
	int		m_nTrain;
	CString	m_strReport;
	CString	m_strALNinputColName;
	int		m_nColsUniv;
	CString	m_strIOprop;
	int		m_nALNinputs;
	int		m_nSortCols;
	CString	m_strImportance;
	CString	m_strMinWeight;
	CString	m_strMaxWeight;
	CString	m_strVersionInfo;
	//}}AFX_DATA
  int m_nMessage;
  int m_nPercentProgress;
  CString m_strIOpropMessageArray[4];

// Attributes
public:
	CALNfitDeepDoc* GetDocument();
	void UpdateControlsFromDoc();

	// Operations
public:

// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CALNfitDeepView)
	public:
	virtual void OnInitialUpdate();
	virtual BOOL Create(LPCTSTR lpszClassName, LPCTSTR lpszWindowName, DWORD dwStyle, const RECT& rect, CWnd* pParentWnd, UINT nID, CCreateContext* pContext = NULL);
	protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support
	virtual BOOL OnPreparePrinting(CPrintInfo* pInfo);
	virtual void OnBeginPrinting(CDC* pDC, CPrintInfo* pInfo);
	virtual void OnEndPrinting(CDC* pDC, CPrintInfo* pInfo);
	virtual void OnUpdate(CView* pSender, LPARAM lHint, CObject* pHint);
	virtual BOOL PreCreateWindow(CREATESTRUCT& cs);
	virtual void CalcWindowRect(LPRECT lpClientRect, UINT nAdjustType = adjustBorder);
	virtual void OnDraw(CDC* pDC);
	virtual void OnPrint(CDC* pDC, CPrintInfo* pInfo);
	//}}AFX_VIRTUAL
  void DoActions();

// Implementation
public:
	void PrintPageHeader(CDC* pDC);
	int m_nPage;
	virtual ~CALNfitDeepView();
#ifdef _DEBUG
	virtual void AssertValid() const;
	virtual void Dump(CDumpContext& dc) const;
#endif

// Generated message map functions

protected:
	//{{AFX_MSG(CALNfitDeepView)
	afx_msg void OnButtonStart();
	afx_msg void OnButtonData();
	afx_msg void OnButtonHelp();
	afx_msg void OnShowWindow(BOOL bShow, UINT nStatus);
  afx_msg LRESULT OnUpdateScreen(WPARAM wParam,LPARAM lParam);
  afx_msg LRESULT OnThreadFinished(WPARAM wParam,LPARAM lParam);
	afx_msg void OnFileSave();
	afx_msg void OnFileSaveAs();
	afx_msg void OnButtonviewprev();
	afx_msg void OnButtonviewnext();
	afx_msg void OnSelchangeAlninputcolname();
	afx_msg void OnVScroll(UINT nSBCode, UINT nPos, CScrollBar* pScrollBar);
	afx_msg void OnButtonremovealninput();
	afx_msg void OnButtonaddalninput();
	afx_msg void OnRadioFit();
	afx_msg void OnRadioClass();
	afx_msg void OnRadioTrain();
	afx_msg void OnRadioEval();
	afx_msg void OnEditClearall();
	afx_msg void OnRadiosortimportance();
	afx_msg void OnRadiosortcols();
	afx_msg void OnUpdateEditmaxwt();
	afx_msg void OnUpdateEditminwt();
	afx_msg void OnKillfocusEditmaxwt();
	afx_msg void OnKillfocusEditminwt();
	afx_msg void OnButtonoptions();
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
private:
	CRect m_rectPrint;
public:
	afx_msg void OnEnChangeDatafilecolumns();
};

#ifndef _DEBUG  // debug version in ALNfitDeepView.cpp
inline CALNfitDeepDoc* CALNfitDeepView::GetDocument()
   { return (CALNfitDeepDoc*)m_pDocument; }
#endif

/////////////////////////////////////////////////////////////////////////////

//{{AFX_INSERT_LOCATION}}
// Microsoft Developer Studio will insert additional declarations immediately before the previous line.

#endif // !defined(AFX_ALNFITDEEPVIEW_H__5CD205CA_22D6_4675_88DD_FB3789D96D18__INCLUDED_)
