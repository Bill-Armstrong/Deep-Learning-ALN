// ALNfitDeepDoc.h : interface of the CALNfitDeepDoc class
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


#pragma once

#define INPUTDEC      0
#define INPUT         1
#define INPUTINC      2
#define OUTPUT        3

class CALNfitDeepDoc : public CDocument
{
protected: // create from serialization only
	CALNfitDeepDoc();
	DECLARE_DYNCREATE(CALNfitDeepDoc)

// Attributes
public:
	CALNfitDeepDoc* m_pDocument;   // pointer to the document object
  CString m_strTentativePathName;  // enables us to suggest a file name for saving doc
  CString m_strDocumentFileName; // name given to the document file
	CString	m_strDTREEFileName;   // the path of the DTREE file
	CString	m_strDataFileName;    // the path of the data file
	CString m_strProtocolFileName;// "Train"/"Eval" + "Protocol" + time +.txt
	CString m_strScatterFileName;  // "Train"/"Eval" + "Scatterplot"   + time +.txt
  CString m_strR_or_E_FileName; // contains replacements of missing values or an added column
	CString	m_strReport;          // string with current action message
	int		m_nFit;                 // 0 for regression, 1 for classification
	int		m_nTrain;               // 0 for training, 1 for evaluation
	int		m_nALNs;                // number of ALNs to be trained and averaged
  int   m_nMessage;             // the number of the message to be displayed
  int   m_nPercentProgress;     // integer percent progress
  // properties of the data file
  int m_nColsUniv;              // number of columns in the first line of the data file
  CString m_strVarname[101];    // variable names for the columns of the data file
  // specification of the ALN inputs and the output
  int m_nALNinputs;  // the number of inputs to the ALNs, including the output
  // properties of the input variables (column number,delay, importance)
  int m_nInputCol[101]; // which column of the data file the input or output comes from (0-based)
  int m_nLag[101];    // lags for the respective entries (default 0, positive means towards start of file)
  int m_nIOpropMessageNumber[101]; // gives input with monotonicity or output with where it goes
  double m_dblMinWeight[101];      // lower bound on weight of the ALN inputs
  double m_dblMaxWeight[101];      // upper bound on weight of the ALN inputs
  double m_dblImportance[101];     // importance values for the respective inputs (default 0)
  // the ALN output is always the highest numbered entry after the ALN inputs; its lag is 0
  int m_nDisplayedALNinput;        // this is the input or output being displayed in the view
  CStringArray m_stringArray;      // for printing parameters
  enum{nLinesPerPage = 56};
  BOOL m_bSortCols;                // TRUE when the view previous/next is in order by the column/lag
  BOOL m_bWarnEqual;               // true when two ALN inputs (col/lag) are the same,
                                   // or if the output at 0 lag is also an input
  BOOL m_bStarted;                  // TRUE if the start button has been pushed until actions end
  int m_nSamplesForTest;

  // document values set by OptionsDialog
	int m_nPercentForTest;
	int m_nEstimateRMSError;
  double m_dblSetTolerance;
	int m_nZeroSmoothing;
	int m_nNoJitter;
	int m_nDTREEDepth;
	BOOL m_bTimePrefixes;
  BOOL m_bReplaceUndefined;
  BOOL m_bDiagnostics;
  // Operations
public:
  void EditClearall();
  void GenerateReportString();     // called by OnUpdateScreen to generate an action report
  void DocUpdateScreen();          // handler of action messages calls this
  void DocThreadFinished();        // handler of action messages calls this
  void GetColumnNames(int ncols);  // gets the column names from an analysis of the data file
  void RemoveInput();              // removes the input currently displayed
  void AddInput();                 // adds a new input to the right of the one being displayed (if not the output)
  void GeneratePrint();
  void OnSort();                    // sorts according to column/lag if bSortCol is TRUE, else by importance
  // Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CALNfitDeepDoc)
	public:
	virtual BOOL OnNewDocument();
	virtual void Serialize(CArchive& ar);
	virtual void DeleteContents();
	//}}AFX_VIRTUAL

// Implementation
public:
	virtual ~CALNfitDeepDoc();
#ifdef _DEBUG
	virtual void AssertValid() const;
	virtual void Dump(CDumpContext& dc) const;
#endif

// Generated message map functions
// protected: removed so these are public
	//{{AFX_MSG(CALNfitDeepDoc)
	afx_msg void OnFileSave();
	afx_msg void OnFileSaveAs();
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
};


/////////////////////////////////////////////////////////////////////////////

//{{AFX_INSERT_LOCATION}}
// Microsoft Developer Studio will insert additional declarations immediately before the previous line.
