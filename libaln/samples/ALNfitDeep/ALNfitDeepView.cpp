// ALNfitDeepView.cpp : implementation of the CALNfitDeepView class
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


#ifdef ALNDLL
#define ALNIMP __declspec(dllexport)
#endif



#include "aln.h"
#include "alnpp.h"
#include "cmyaln.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

#include "stdafx.h"
#include "ALNfitDeep.h"
#include "ALNfitDeepDoc.h"
#include "ALNfitDeepView.h"
#include "OptionsDialog.h"
#include <fstream>
#include <math.h>
#include <datafile.h>
#include "alnextern.h"

extern CMyAln * pALN;

#define INPUTDEC      0
#define INPUT         1
#define INPUTINC      2
#define OUTPUT        3

/////////////////////////////////////////////////////////////////////////////
// CALNfitDeepView

IMPLEMENT_DYNCREATE(CALNfitDeepView, CFormView)

BEGIN_MESSAGE_MAP(CALNfitDeepView, CFormView)
	//{{AFX_MSG_MAP(CALNfitDeepView)
	ON_BN_CLICKED(IDC_BUTTON_START, OnButtonStart)
	ON_BN_CLICKED(IDC_BUTTON_DATA, OnButtonData)
	ON_BN_CLICKED(IDC_BUTTON_HELP, OnButtonHelp)
	ON_WM_SHOWWINDOW()
	ON_MESSAGE(WM_UPDATESCREEN, OnUpdateScreen)
	ON_MESSAGE(WM_THREADFINISHED, OnThreadFinished)
	ON_COMMAND(ID_FILE_SAVE, OnFileSave)
	ON_COMMAND(ID_FILE_SAVE_AS, OnFileSaveAs)
	ON_BN_CLICKED(IDC_BUTTONVIEWPREV, OnButtonviewprev)
	ON_BN_CLICKED(IDC_BUTTONVIEWNEXT, OnButtonviewnext)
	ON_CBN_SELCHANGE(IDC_ALNINPUTCOLNAME, OnSelchangeAlninputcolname)
	ON_WM_VSCROLL()
	ON_BN_CLICKED(IDC_BUTTONREMOVEALNINPUT, OnButtonremovealninput)
	ON_BN_CLICKED(IDC_BUTTONADDALNINPUT, OnButtonaddalninput)
	ON_BN_CLICKED(IDC_RADIO_FIT, OnRadioFit)
	ON_BN_CLICKED(IDC_RADIO_CLASS, OnRadioClass)
	ON_BN_CLICKED(IDC_RADIO_TRAIN, OnRadioTrain)
	ON_BN_CLICKED(IDC_RADIO_EVAL, OnRadioEval)
	ON_COMMAND(ID_EDIT_CLEARALL, OnEditClearall)
	ON_BN_CLICKED(IDC_RADIOSORTIMPORTANCE, OnRadiosortimportance)
	ON_BN_CLICKED(IDC_RADIOSORTCOLS, OnRadiosortcols)
	ON_EN_UPDATE(IDC_EDITMAXWT, OnUpdateEditmaxwt)
	ON_EN_UPDATE(IDC_EDITMINWT, OnUpdateEditminwt)
	ON_EN_KILLFOCUS(IDC_EDITMAXWT, OnKillfocusEditmaxwt)
	ON_EN_KILLFOCUS(IDC_EDITMINWT, OnKillfocusEditminwt)
	ON_BN_CLICKED(IDC_BUTTONOPTIONS, OnButtonoptions)
	//}}AFX_MSG_MAP
	// Standard printing commands
	ON_COMMAND(ID_FILE_PRINT, CFormView::OnFilePrint)
	ON_COMMAND(ID_FILE_PRINT_DIRECT, CFormView::OnFilePrint)
	ON_COMMAND(ID_FILE_PRINT_PREVIEW, CFormView::OnFilePrintPreview)
	ON_EN_CHANGE(IDC_DATAFILECOLUMNS, &CALNfitDeepView::OnEnChangeDatafilecolumns)
END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////////
// CALNfitDeepView construction/destruction

CALNfitDeepView::CALNfitDeepView()
	: CFormView(CALNfitDeepView::IDD),m_rectPrint(0,0,11520,-15120)
{
	//{{AFX_DATA_INIT(CALNfitDeepView)
	m_strDataFileName = "";
	m_strDTREEFileName = "";
	m_nFit = 0;  // 0 means TRUE, a .fit file will be used to evaluate a file 
	m_nTrain = 0;  // 0 means TRUE, the net will be trained
	m_strReport = "";
	m_strALNinputColName = _T("");
	m_nColsUniv = 0;
	m_strIOprop = _T("");
	m_nALNinputs = 0;
	m_nSortCols = -1;
	m_strImportance = _T("");
	m_strMinWeight = _T("");
	m_strMaxWeight = _T("");
	m_strVersionInfo = _T("");
	//}}AFX_DATA_INIT
  m_strIOpropMessageArray[0] ="INPUT  \\";
  m_strIOpropMessageArray[1] ="INPUT   ";
  m_strIOpropMessageArray[2] ="INPUT  /";
  m_strIOpropMessageArray[3] ="OUTPUT (lag is 0)";
}

CALNfitDeepView::~CALNfitDeepView()
{
}

void CALNfitDeepView::DoDataExchange(CDataExchange* pDX)
{
	CFormView::DoDataExchange(pDX);
	//{{AFX_DATA_MAP(CALNfitDeepView)
	DDX_Control(pDX, IDC_PROGRESS, m_Progress);
	DDX_Text(pDX, IDC_EDIT_DATA, m_strDataFileName);
	DDX_Radio(pDX, IDC_RADIO_FIT, m_nFit);
	DDX_Radio(pDX, IDC_RADIO_TRAIN, m_nTrain);
	DDX_Text(pDX, IDC_EDIT_REPORT, m_strReport);
	DDX_CBString(pDX, IDC_ALNINPUTCOLNAME, m_strALNinputColName);
	DDX_Text(pDX, IDC_DATAFILECOLUMNS, m_nColsUniv);
	DDX_Text(pDX, IDC_IOPROP, m_strIOprop);
	DDX_Text(pDX, IDC_NALNINPUTS, m_nALNinputs);
	DDX_Radio(pDX, IDC_RADIOSORTCOLS, m_nSortCols);
	DDX_Text(pDX, IDC_EDITIMPORTANCE, m_strImportance);
	DDX_Text(pDX, IDC_EDITMINWT, m_strMinWeight);
	DDV_MaxChars(pDX, m_strMinWeight, 32);
	DDX_Text(pDX, IDC_EDITMAXWT, m_strMaxWeight);
	DDV_MaxChars(pDX, m_strMaxWeight, 32);
	//}}AFX_DATA_MAP
}

/////////////////////////////////////////////////////////////////////////////
// CALNfitDeepView printing

BOOL CALNfitDeepView::OnPreparePrinting(CPrintInfo* pInfo)
{
  CALNfitDeepDoc* pDoc = (CALNfitDeepDoc*)GetDocument();
  if(pDoc->m_bStarted)
  {
    MessageBox("Print cancelled. Please wait until the run is finished for updated values.");
	  return 0;
  } 
  // the total printout is
  // 10 lines info
  // ceil(nColsUniv /5) lines for the column names
  // 1 line for the number of ALN inputs
  // 1 header line for the input descriptions
  // nALNinputs: two lines each except the last
  // We also allow two lines for the header on each page
  int nTotalLines = 12 + (int) ceil(((float)pDoc->m_nColsUniv)/5) + 2 * pDoc->m_nALNinputs;
  int nMaxPages = (int) ceil((double)nTotalLines/(double)CALNfitDeepDoc::nLinesPerPage);
  pInfo->SetMaxPage(nMaxPages);
	return DoPreparePrinting(pInfo);
}

void CALNfitDeepView::OnBeginPrinting(CDC* /*pDC*/, CPrintInfo* /*pInfo*/)
{

}

void CALNfitDeepView::OnEndPrinting(CDC* /*pDC*/, CPrintInfo* /*pInfo*/)
{
	// TODO: add cleanup after printing
}

void CALNfitDeepView::OnDraw(CDC* pDC)
{

}


/////////////////////////////////////////////////////////////////////////////
// CALNfitDeepView diagnostics

#ifdef _DEBUG
void CALNfitDeepView::AssertValid() const
{
	CFormView::AssertValid();
}

void CALNfitDeepView::Dump(CDumpContext& dc) const
{
	CFormView::Dump(dc);
}

CALNfitDeepDoc* CALNfitDeepView::GetDocument() // non-debug version is inline
{
	ASSERT(m_pDocument->IsKindOf(RUNTIME_CLASS(CALNfitDeepDoc)));
	return (CALNfitDeepDoc*)m_pDocument;
}
#endif //_DEBUG

/////////////////////////////////////////////////////////////////////////////
// CALNfitDeepView message handlers

void CALNfitDeepView::OnInitialUpdate() 
{
	// called on startup
	CALNfitDeepDoc* pDoc = (CALNfitDeepDoc*)GetDocument();
	CProgressCtrl* pProg = (CProgressCtrl*) GetDlgItem(IDC_PROGRESS);
	pProg->SetRange(0,100);
	CSpinButtonCtrl* pSpin = (CSpinButtonCtrl* ) GetDlgItem(IDC_SPIN_NALNS);
	pSpin->SetRange(1,10);
  pSpin->SetPos(pDoc->m_nALNs);
  CComboBox* pCLB = (CComboBox*) GetDlgItem(IDC_ALNINPUTCOLNAME);
  pCLB->SetDroppedWidth( 200);
  pCLB->ResetContent();
  for(int ii=0; ii< pDoc->m_nColsUniv; ii++)
  {
    pCLB->InsertString(ii,pDoc->m_strVarname[ii]);
  }
  // setting the initial value 0 and the range for the lag spin button
  pSpin = (CSpinButtonCtrl* ) GetDlgItem(IDC_SPINLAG);
  pSpin->SetRange(0,52);
  pSpin->SetPos(pDoc->m_nLag[pDoc->m_nDisplayedALNinput]);
  pSpin = (CSpinButtonCtrl* ) GetDlgItem(IDC_SPINALNIOPROP);
  pSpin->SetRange(0,4);
  pSpin->SetPos(pDoc->m_nIOpropMessageNumber[pDoc->m_nDisplayedALNinput]);
  CFormView::OnInitialUpdate();
}

void CALNfitDeepView::UpdateControlsFromDoc()
{
	// called from OnInitialUpdate and OnEditClearAll
  // and by OnUpdate which is called by UpdateAllViews
  // assumes the doc is updated completely
	CALNfitDeepDoc* pDoc = (CALNfitDeepDoc*)GetDocument();
	m_strDataFileName= pDoc->m_strDataFileName;
	m_strDTREEFileName= pDoc->m_strDTREEFileName;
	strcpy(szDTREEFileName, LPCTSTR( pDoc->m_strDTREEFileName)); // does this work?
  m_nMessage = pDoc->m_nMessage;
	m_strReport = pDoc->m_strReport;
	CEdit* pEdit = (CEdit* ) GetDlgItem(IDC_EDIT_REPORT);
	pEdit->SetWindowText(m_strReport);
	CSpinButtonCtrl* pSpin = (CSpinButtonCtrl* ) GetDlgItem(IDC_SPIN_NALNS);
	pSpin->SetPos(pDoc->m_nALNs);
	m_nFit = pDoc->m_nFit;
	m_nTrain = pDoc->m_nTrain;
	m_nPercentProgress = pDoc->m_nPercentProgress;
	CProgressCtrl* pProg = (CProgressCtrl*) GetDlgItem(IDC_PROGRESS);
	pProg->SetPos(m_nPercentProgress);
  // pro version updates
  m_nColsUniv = pDoc->m_nColsUniv;
  m_strALNinputColName = pDoc->m_strVarname[pDoc->m_nInputCol[pDoc->m_nDisplayedALNinput]];
  CComboBox* pCLB = (CComboBox*) GetDlgItem(IDC_ALNINPUTCOLNAME);
  pCLB->SetCurSel(pDoc->m_nInputCol[pDoc->m_nDisplayedALNinput]);
  m_nALNinputs = nALNinputs = pDoc->m_nALNinputs;
  pSpin = (CSpinButtonCtrl* ) GetDlgItem(IDC_SPINLAG);
	pSpin->SetPos(pDoc->m_nLag[pDoc->m_nDisplayedALNinput]);
  int nTemp = pDoc->m_nIOpropMessageNumber[pDoc->m_nDisplayedALNinput];
  m_strIOprop =  (CString) m_strIOpropMessageArray[nTemp];
	pSpin = (CSpinButtonCtrl* ) GetDlgItem(IDC_SPINALNIOPROP);
	pSpin->SetPos(nTemp);
  if(pDoc->m_nDisplayedALNinput < pDoc->m_nALNinputs -1) 
  {
    if((nTemp == INPUTDEC)&&(pDoc->m_dblMaxWeight[pDoc->m_nDisplayedALNinput]>0)) pDoc->m_dblMaxWeight[pDoc->m_nDisplayedALNinput] = 0;
    if((nTemp == INPUTINC)&&(pDoc->m_dblMinWeight[pDoc->m_nDisplayedALNinput]<0)) pDoc->m_dblMinWeight[pDoc->m_nDisplayedALNinput] = 0;
  }
  else
  {
    pDoc->m_dblMaxWeight[pDoc->m_nDisplayedALNinput] = -1;
    pDoc->m_dblMinWeight[pDoc->m_nDisplayedALNinput] = -1;
  }
  CString strVal;
  strVal.Format("%lf", (double)pDoc->m_dblMinWeight[pDoc->m_nDisplayedALNinput]);
  m_strMinWeight = strVal;
  strVal.Format("%lf", (double)pDoc->m_dblMaxWeight[pDoc->m_nDisplayedALNinput]);
  m_strMaxWeight = strVal;
  strVal.Format("%lf", (double)pDoc->m_dblImportance[pDoc->m_nDisplayedALNinput]);
  m_strImportance = strVal;
  if(pDoc->m_bSortCols == TRUE)
  {
    m_nSortCols = 0;
  }
  else
  {
    m_nSortCols = 1;
  }
	UpdateData(FALSE); // calls DDX
}

void CALNfitDeepView::OnButtonStart() 
{
	CALNfitDeepDoc* pDoc = (CALNfitDeepDoc*)GetDocument();
  if(pDoc->m_bStarted)
  {
    MessageBox("Already started");
    return;
  }
  m_strDataFileName = ""; 
  m_nFit = 0; // if a .fit file has been opened, the GUI will show this as 1
  m_nTrain = 0;
	bTrain = TRUE;
  UpdateData(TRUE); // get what is in the GUI into the view members and transfer to doc
  if(pDoc->m_strDataFileName != m_strDataFileName)
  {
    MessageBox("Please Browse to the data file before clicking Start.");
    return;
  }
  // set up all the other files
  CFile fData,fDTREE,fDocument;
  CFileException e;
  if( !fData.Open( m_strDataFileName, CFile::modeRead, &e ) )
  {
		MessageBox("Data file could not be opened. Please Browse to an existing data file *.txt.");
    return;
  }
  fData.Close();
  // just a test
  // start has been clicked, so we have the ALN variables defined
  // we can construct the files for training and testing that get input values
  // from various columns at various delays.
  
  nALNs = pDoc->m_nALNs;
  // set up the naming convention for the output files
  CTime theTime = CTime::GetCurrentTime();
	// We can get more fields: theTime.Format( "_%b_%d_%Y_%H_%M" );
	CString strPrefix, strInfix, strSuffix, strMessage;
  int nLen, nBSpos;
  TCHAR backslash = TCHAR(92);
  nLen = m_strDataFileName.GetLength(); 
  strSuffix = m_strDataFileName.Left(nLen-4);
  nBSpos = strSuffix.ReverseFind(backslash);
  strSuffix = strSuffix.Right(nLen - 4 - nBSpos -1);  // chop the path off
  if(bTimePrefixes)
  {
    strPrefix = theTime.Format( "%H%M" );
  }
  else
  {
    strPrefix = "";
  }
  strPrefix = strPrefix + strSuffix; // time and data file name -- good as prefix for all files
  if(m_nTrain == 0)
  {

    strSuffix = ".fit";
    pDoc->m_strTentativePathName = strPrefix + strSuffix;
    /*
    strMessage = "Save to a name like ";
    strSuffix = "???.fit";
    strMessage = strMessage + strPrefix + strSuffix;
    strSuffix = " where ??? is a comment";
    strMessage = strMessage + strSuffix;
    MessageBox(LPCTSTR(strMessage));
    */
    OnFileSaveAs();
  }
  // we name the other files accordingly
  // now check out the DTREE file name in the fit file
  if(m_nTrain == 1) // suppose we have a .fit file and this is evaluation
  {
    strSuffix = m_strDTREEFileName.Right(4);
    if((strSuffix != ".dtr") || !fDTREE.Open(pDoc->m_strDTREEFileName, CFile::modeRead))
    {
      // The file doesn't exist
		  MessageBox("The DTREE file to be used for evaluating the data file does not exist");
      return;
    }
  }
  else // we have not opened a .fit file, and so are training a new DTREE
  {
		strSuffix = "DTREE.dtr";
    pDoc->m_strDTREEFileName = m_strDTREEFileName = strPrefix + strSuffix;
    fDTREE.Open(m_strDTREEFileName, CFile::modeCreate|CFile::modeWrite);
    fDTREE.Close();
  }

  // the data, .fit and .dtr files exist, and we know if we are training or evaluating
  // we only know if we are classifying or not if we are evaluating, otherwise
  // we have to wait until analyzeTV
  UpdateData(FALSE);
  bTrain = (m_nTrain == 0)?TRUE:FALSE;
  bClassify = (m_nFit != 0)?TRUE:FALSE; // if training, we will set this later
  // create names for the Protocol, Output, and Replacement files based on the data file name
  if(m_nTrain == 0)
  {
    strSuffix = "TrainProtocol.txt";
  }
  else
  {
    strSuffix = "EvalProtocol.txt";
  }
	pDoc->m_strProtocolFileName = strPrefix + strSuffix;
  if(m_nTrain == 0)
  {
    strSuffix = "TrainScatterPlot.txt";
  }
  else
  {
    strSuffix = "EvalScatterPlot.txt";
  }
	pDoc->m_strScatterFileName = strPrefix + strSuffix;
  // now we replace some missing values in the input file
	if(bReplaceUndefined == TRUE)
	{
		strSuffix = "R.txt"; // this is kept short so there can be many replacements
	}
	else
	{
		strSuffix = "E.txt"; // this is for files with an extra output column at the right
	}
  pDoc->m_strR_or_E_FileName = strPrefix + strSuffix;
  // check that we have one input at the highest index
  int nBadInput;
  //CString strMessage;
  BOOL bIOpropErrorFlag = FALSE;
  for(int ii = 0; ii < pDoc->m_nALNinputs - 1; ii++)
  {
    if(pDoc->m_nIOpropMessageNumber[ii] > INPUTINC)
    {
      bIOpropErrorFlag = TRUE;
      nBadInput = ii;
    }
  }
  int ii = pDoc->m_nALNinputs - 1;
  if(pDoc->m_nIOpropMessageNumber[ii] <= INPUTINC)
  {
    bIOpropErrorFlag = TRUE;
    nBadInput = ii;
  }
  if(bIOpropErrorFlag)
  {
    strMessage.Format("Error: the input or output property number %d is not correct! Please correct and restart.",nBadInput + 1);
		MessageBox(strMessage);
    return;
  }
  BOOL bWeightMonotonicityConflict = FALSE;
  for(ii = 0; ii < pDoc->m_nALNinputs - 1; ii++)
  {
    if((pDoc->m_nIOpropMessageNumber[ii] == INPUTINC) && (pDoc->m_dblMaxWeight[ii] < 0))
    {
      bWeightMonotonicityConflict = TRUE;
      nBadInput = ii;
    }
    if((pDoc->m_nIOpropMessageNumber[ii] == INPUTDEC) && (pDoc->m_dblMinWeight[ii] > 0)) 
    {
      bWeightMonotonicityConflict = TRUE;
      nBadInput = ii;
    }
    if(pDoc->m_dblMaxWeight[ii] < pDoc->m_dblMinWeight[ii]) 
    {
      bWeightMonotonicityConflict = TRUE;
      nBadInput = ii;
    }
  }
  if(bWeightMonotonicityConflict)
  {
    strMessage.Format("Error: the monotonicity of ALN input number %d and the weight bounds for it are inconsistent! Please correct and restart.",nBadInput + 1);
    MessageBox(strMessage);
    return;
  }
  pDoc->m_bSortCols = TRUE;
  pDoc->OnSort();
  pDoc->m_bStarted = TRUE;
  // now get on with training or evaluation
  DoActions();
}

void CALNfitDeepView::DoActions()
{
  CALNfitDeepDoc* pDoc = (CALNfitDeepDoc*)GetDocument();
	// Convert the task type and the train/use input to booleans for use in ALN code
	// bClassify = (m_nFit != 0)?TRUE:FALSE; Now we don't get this from the interface
	bTrain = (m_nTrain == 0)?TRUE:FALSE; // set the number of ALNs in the program
	nALNs = pDoc->m_nALNs;
  // copy filenames into global strings
  strcpy(szDTREEFileName,pDoc->m_strDTREEFileName);
  strcpy(szDataFileName,pDoc->m_strDataFileName);
	strcpy(szProtocolFileName,pDoc->m_strProtocolFileName);
  strcpy(szR_or_E_FileName,pDoc->m_strR_or_E_FileName);
	strcpy(szScatterFileName,pDoc->m_strScatterFileName);
	// now start program actions
  ALNfitSetup();
  if(bTrain)
  {
    analyzeTV(); // moved from training thread, but only needed for training
  }
  m_nFit = pDoc->m_nFit = bClassify?1:0;
  UpdateControlsFromDoc();
  // start the training/evaluation thread
  CWinThread* pThread= 
  AfxBeginThread(ActionsProc,GetSafeHwnd(),THREAD_PRIORITY_NORMAL);
}


void CALNfitDeepView::OnButtonData() 
{
	CALNfitDeepDoc* pDoc = (CALNfitDeepDoc*)GetDocument();
  if(pDoc->m_bStarted == FALSE)
  {
	  CFileDialog dlg(TRUE,"txt","*.txt");
	  CFileException e;
	  CFile file;
    if(dlg.DoModal() ==IDOK)
	  {
		  BOOL bSuccess;
      bSuccess = file.Open(dlg.GetPathName(),CFile::modeRead, &e);
		  if(bSuccess)
      {
        // upon success with opening the data file, the DiagnoseFileSetupProtocol
        // can be opened in the same directory
         fpFileSetupProtocol = fopen("DiagnoseFileSetupProtocol.txt","w");
        if(fpFileSetupProtocol)
        {
          if(bPrint)fprintf(fpFileSetupProtocol,"The Data file was successfully opened.\n");
	        if(bPrint)fflush(fpFileSetupProtocol);
        }
        else
        {
          MessageBox("Problem opening DiagnoseFileSetupProtocol.txt");
					return;
        }
	      // since this handler will update in this case, we must transfer
        // other unsaved GUI data to the document and to the external variables
        UpdateData(TRUE);
	      pDoc->m_nFit = m_nFit;
        pDoc->m_nTrain = m_nTrain;
        // pDoc->m_nALNs is updated by OnVScroll, not by UpdateData(TRUE)
        nALNs = pDoc->m_nALNs;
        pDoc->m_nALNinputs = m_nALNinputs;
        bClassify = (m_nFit == 1)?TRUE:FALSE;
        bTrain = (m_nTrain == 0)?TRUE:FALSE;
			  pDoc->m_strDataFileName = m_strDataFileName;
        // now we get the new data file name from the dialog
			  m_strDataFileName = dlg.GetPathName();
			  pDoc->m_strDataFileName = m_strDataFileName;
        file.Close(); // just testing the first time
        // now we analyze the datafile and get the columns and default ALN variables
	      int nheaderlines = 0;
        long nrows = 0;
        int ncols = 0;
	      if(!analyzeInputFile((char *) LPCTSTR(m_strDataFileName), &nheaderlines, &nrows, &ncols,TRUE))
	      {
          MessageBox("Stopping: Input file analysis failed");
		      exit(0);
	      }
        if(bPrint)fprintf(fpFileSetupProtocol,"Data file analysis returned: header lines %d, rows %ld, columns %d \n", nheaderlines, nrows, ncols);
	      if(bPrint)fflush(fpFileSetupProtocol);
        BOOL bColumnsSame = TRUE;
        CString strColName;
        if(m_nTrain == 1)  // if we are evaluating, we have opened a .fit file
        {
          if(pDoc->m_nColsUniv == ncols)
          {
            pDoc->m_nDisplayedALNinput = pDoc->m_nALNinputs - 1; // display the output variable by default
          }
          else
          {
            MessageBox("Wrong number of columns in the data file. Please browse to a new data file.");
            return;
          }
        }
	      if(nheaderlines > 9)
	      {
          MessageBox("Stopping: You can't have more than nine headerlines in your data file! Please modify file.");
		      return;
	      }
				nheaderlinewithvars = -1;  //set the headerline with the variables to an impossible value
				if (nheaderlines > 0)
				{
					//does either of the last two headerlines have the correct number of strings to be variable names?
					if (nheaderlineitems[nheaderlines - 1] == ncols)
					{
						nheaderlinewithvars = nheaderlines - 1;
					}
					else
					{
						if ((nheaderlines > 1) && (nheaderlineitems[nheaderlines - 2] == ncols)) nheaderlinewithvars = nheaderlines - 2;
					}
				}
				else
				{
					// there is no acceptable list of column names and we have to synthesize them
					MessageBox("Can't find column names for the data file. Using defaults.");
					for (int nn = 0; nn < ncols; nn++)
					{
						char buff[32] = "_A";
						if (nn / 26 > 0)buff[0] = 'A' - 1 + nn / 26;
						buff[1] = 'A' + nn % 26;
						buff[2] = '\0';
						strColName = buff;
						if ((m_nTrain == 1) && (pDoc->m_strVarname[nn] != strColName)) bColumnsSame = FALSE;
						pDoc->m_strVarname[nn] = buff;
						for (int charno = 0; charno < 32; charno++)
						{
							varname[nn][charno] = buff[charno];
						}
					}
				}
				if(nheaderlinewithvars > 0) // nheaderlinewithvars is nheaderlines minus 1 or 2 and we have to extract the names of the variables
				{
					char item[101][32];
					// get the names of variables from the header line with vars
					int itemcount = sscanf(hdrline[nheaderlinewithvars], "%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s",
						item + 0, item + 1, item + 2, item + 3, item + 4, item + 5, item + 6, item + 7, item + 8, item + 9,
						item + 10, item + 11, item + 12, item + 13, item + 14, item + 15, item + 16, item + 17, item + 18, item + 19,
						item + 20, item + 21, item + 22, item + 23, item + 24, item + 25, item + 26, item + 27, item + 28, item + 29,
						item + 30, item + 31, item + 32, item + 33, item + 34, item + 35, item + 36, item + 37, item + 38, item + 39,
						item + 40, item + 41, item + 42, item + 43, item + 44, item + 45, item + 46, item + 47, item + 48, item + 49,
						item + 50, item + 51, item + 52, item + 53, item + 54, item + 55, item + 56, item + 57, item + 58, item + 59,
						item + 60, item + 61, item + 62, item + 63, item + 64, item + 65, item + 66, item + 67, item + 68, item + 69,
						item + 70, item + 71, item + 72, item + 73, item + 74, item + 75, item + 76, item + 77, item + 78, item + 79,
						item + 80, item + 81, item + 82, item + 83, item + 84, item + 85, item + 86, item + 87, item + 88, item + 89,
						item + 90, item + 91, item + 92, item + 93, item + 94, item + 95, item + 96, item + 97, item + 98, item + 99, item + 100);
					// set the variable names in the document and the external variable varname
					if (bPrint)fprintf(fpFileSetupProtocol, "The line read from the file with variable names is: \n %s \n", hdrline[nheaderlinewithvars]);
					if (bPrint)fprintf(fpFileSetupProtocol, "The line above is: \n %s \n", hdrline[nheaderlinewithvars-1]);
					if (bPrint)fprintf(fpFileSetupProtocol, "The line below is: \n %s \n", hdrline[nheaderlinewithvars+1]);    // THESE LINES SEEM EMPTY!!!!
					if (bPrint)fprintf(fpFileSetupProtocol, "First variable name %s \n", item[0]);
					for (int nn = 0; nn < ncols; nn++)
					{
						strColName = item[nn];
						if ((m_nTrain == 1) && (pDoc->m_strVarname[nn] != strColName)) bColumnsSame = FALSE;
						pDoc->m_strVarname[nn] = item[nn];
						// copy the names read to the external variable varname
						for (int charno = 0; charno < 32; charno++)
						{
							varname[nn][charno] = item[nn][charno];
						}
					}
				}
				if(bColumnsSame == FALSE)
        {
          MessageBox("CAUTION: the names of columns in the data file do not match those in the .fit file. Continuing...");
        }
        CComboBox* pCLB = (CComboBox*) GetDlgItem(IDC_ALNINPUTCOLNAME);
        pCLB->ResetContent();
        for(int ii=0; ii < ncols; ii++)
        {
          pCLB->InsertString(ii,LPCTSTR(pDoc->m_strVarname[ii]));
        }
        // one way or another we have an itemcount >0  and we can set the defaults
        //int nresponse = AfxMessageBox("Set up the default input variables, output, etc.?",MB_YESNO,0);
        //if(nresponse == IDNO)
        if(m_nTrain == 0)
        {
          // YES - use the default ALN inputs and output
          pDoc->m_nColsUniv = nColsUniv = ncols;
          m_nALNinputs = pDoc->m_nALNinputs = nALNinputs = nColsUniv;
          pDoc->m_nDisplayedALNinput = m_nALNinputs - 1; // display the output variable by default
          pCLB->SetCurSel(pDoc->m_nDisplayedALNinput);
          for(int nn = 0; nn < nColsUniv; nn++)
          {
            pDoc->m_nInputCol[nn] = nInputCol[nn] = nn;
            pDoc->m_nLag[nn] = nLag[nn] = 0;
            if(nn == nColsUniv -1)
            {
              // the output variable always has weight -1
              pDoc->m_nIOpropMessageNumber[nColsUniv - 1] = OUTPUT;
              pDoc->m_dblMinWeight[nn] = dblMinWeight[nn] = -1;
              pDoc->m_dblMaxWeight[nn] = dblMaxWeight[nn] = -1;
            }
            else
            {
              // otherwise set the defaults as generously as possible
              pDoc->m_nIOpropMessageNumber[nn] = INPUT;
              pDoc->m_dblMinWeight[nn] = dblMinWeight[nn] = -dblMax;
              pDoc->m_dblMaxWeight[nn] = dblMaxWeight[nn] =  dblMax;
            }
            pDoc->m_dblImportance[nn] = dblImportance[nn] = 0;
          }
        }
        pDoc->UpdateAllViews(NULL);
		  }
		  else
		  {
			  MessageBox("Can't open file.");
		  }
	  }
  }
  else
  {
    MessageBox("You can't change the input data file until the run is finished.");
  }
  fflush(fpFileSetupProtocol);
}

BOOL CALNfitDeepView::Create(LPCTSTR lpszClassName, LPCTSTR lpszWindowName, DWORD dwStyle, const RECT& rect, CWnd* pParentWnd, UINT nID, CCreateContext* pContext) 
{
	// TODO: Add your specialized code here and/or call the base class
	
	return CFormView::Create(lpszClassName,lpszWindowName, dwStyle, rect, pParentWnd, nID, pContext);
}

void CALNfitDeepView::OnUpdate(CView* pSender, LPARAM lHint, CObject* pHint) 
{
	UpdateControlsFromDoc();
}

void CALNfitDeepView::OnButtonHelp() 
{
	MessageBox("Getting started doing regression:\n\
1. Ignore the Task and Action boxes at the top of the screen.  Browse to a tab-separated data file of numbers. It can have column names in one of the last two lines of a file header of up to 9 lines.\n\
2. Initially, connections are 1-to-1 from the ALN inputs to the data file columns, with zero lags. The column for the desired output should be the rightmost ALN connection. If it isn't, change it using the pull-down list of column names. \n\
3. In the latter case, use the Previous and Next connection buttons to get to the original connection to the the desired output column and Remove that connection as an input. Also Remove columns which are constant or useless (e.g. sample numbers).\n\
4. Click the start button, save the .fit file for later use on new data files and wait for output.  The files belonging to the same run have the same time prefix. Examine the ... TrainProtocol.txt file, noting the estimated RMS noise level for future use.\n\
4. The ALN-estimated output is in the right-hand column in the ...E.txt output file, which extends the data file on the right.\n\
5. The Scatterplot output file is like the E file, but is for the test set. Use a spreadsheet to make a scatterplot of the two rightmost columns.\n\
\n\
   Additional things to try:\n\
A. The tolerance (options dialog) used to be the constant RMS error level below which a flat piece won't split into two. It is now not used since the noise variance need not be constant. It is determined using estimates of noise variance.\n\
B. You can maybe improve results on new data sets by averaging over more ALNs (bagging up to 10).\n\
C. You can use an earlier setup by opening a .fit file (menu Open). Then choose a data file from the same source. Click Start to evaluate the new file using the previous training.\n\
D. This program is not optimized for classification, which is done by regression with integer outputs. If you try it, it may work.\n\
E. If the data samples are a time series, you can select appropriate columns and lags for multiple time series analysis.\n\
F. If you don't know which connections are useful for computing a certain output, remove some connections with the lowest importance values after training. \n\
G. If theory says the function to be learned must be monotonic in some input variable(s), move to that connection and enter the constraint which is shown by a slash or backslash.\
More generally, enter a known bound on the rate of change with respect to any input variable (connection) as a weight (i.e. partial derivative) constraint.\n\
H. From the set of samples (rows of the data file), a percentage specified in Processing options will be removed for testing.  All of the rest will be used for determining the noise level and training.\n\
I. If there are missing values in some column of the the data, represent them by 99999 or -9999 and train on some columns without missing values to replace them by ALN-computed values in the R output file. (Options dialog check Replace) Browse to the R file as the data file to replace missing values another column.\n\
J. Smoothing with quadratic fillets allows fewer flat pieces to be used to fit a function to a given accuracy. This program uses the option of zero smoothing to greatly increase speed and allow DEEP LEARNING with an ALN of many layers.\n\
K. To check linear regression and stats using ALNfitDeep, look at the How to .. spreadsheet in the libaln\\samples\\realestate folder.");
}


BOOL CALNfitDeepView::PreCreateWindow(CREATESTRUCT& cs)
{
	// we got rid of the PreCreateWindow in the MainFrame!
	// THIS DOESN"T WORK!!!
	cs.style &= ~FWS_ADDTOTITLE;
	cs.lpszClass =
		AfxRegisterWndClass(CS_DBLCLKS | CS_HREDRAW | CS_VREDRAW, 0, ::CreateSolidBrush(RGB(255, 255, 255)), 0);
	if (cs.lpszClass == NULL || !CWnd::PreCreateWindow(cs))
	{
		return FALSE;
	}
	else
	{
		return TRUE;
	}
}

// The following functions are global, not part of the view class
// takes two numbers and assigns them to the globals for message number and percent progress
void PassBackStatus(int nMess, int nPer)
{
  //  This updates global variables
  nMessageNumber = nMess;
  nPercentProgress = nPer;
}

LRESULT CALNfitDeepView::OnUpdateScreen(WPARAM wParam,LPARAM lParam)
{
  CALNfitDeepDoc* pDoc = (CALNfitDeepDoc*)GetDocument();
  if(nMessageNumber == 4)
  {
    for(int input = 0; input < pDoc->m_nALNinputs; input++)
    {
      pDoc->m_dblImportance[input] = dblImportance[input];
    }
  }
  pDoc->DocUpdateScreen();
  return 0L;
}

LRESULT CALNfitDeepView::OnThreadFinished(WPARAM wParam,LPARAM lParam)
{	    
  CALNfitDeepDoc* pDoc = (CALNfitDeepDoc*)GetDocument();
  pDoc->DocThreadFinished();
  pDoc->m_bStarted = FALSE;
  if(bTrain)
  {
    OnFileSave(); // this allows us to save the importances
    MessageBox("The following text files have been generated: \n" +
    pDoc->m_strProtocolFileName +
    "\n" + pDoc->m_strDTREEFileName  +
    "\n" + pDoc->m_strScatterFileName  +
    "\n" + pDoc->m_strR_or_E_FileName);
  }
  else
  {
    MessageBox("The following text files have been generated: \n" +
    pDoc->m_strProtocolFileName +
    "\n" + pDoc->m_strScatterFileName  +
    "\n" + pDoc->m_strR_or_E_FileName);
  }
  return 0L;
}

UINT ActionsProc(LPVOID pParam)  // the actions thread
{
  if(bTrain) // procedure for training
  {
    // WM_UPDATESCREEN is user-defined message
    // sent to the view,whose thread created the window
    // and from where this worker thread was launched.
    // The view calls the real handler in the doc
    PassBackStatus(0,5);
    PostMessage((HWND) pParam, WM_UPDATESCREEN,0,0);
    if(bClassify)
    {
      fprintf(fpProtocol, "\n**************  The problem is to classify samples into two or a few classes  *****\n");
    }
    else
    {
      fprintf(fpProtocol, "\n**************  The problem is to fit samples with a smooth function  *****\n");
    }
    fflush(fpProtocol);
    if(bEstimateNoiseVariance)  // We are doing RMS Error estimation (works also with two-class classification)
    {
      // We do the following if we are doing regression and estimating noise variance.
			// Noise estimation may work for classification, but more study is required.
      PassBackStatus(1,10);  
      PostMessage((HWND) pParam, WM_UPDATESCREEN,0,0);
			doLinearRegression();
      PassBackStatus(2,15);  
      ::PostMessage((HWND) pParam, WM_UPDATESCREEN,0,0);
			createNoiseVarianceFile();
			//trainNoiseVarianceALN();
		}
    PassBackStatus(3,30);  
    ::PostMessage((HWND) pParam, WM_UPDATESCREEN,0,0);
		approximate();
    PassBackStatus(7,75);   
    ::PostMessage((HWND) pParam, WM_UPDATESCREEN,0,0);
    outputTrainingResults();
    PassBackStatus(4,80);  
    ::PostMessage((HWND) pParam, WM_UPDATESCREEN,0,0);
    trainAverage();
    PassBackStatus(5,85);  
    ::PostMessage((HWND) pParam, WM_UPDATESCREEN,0,0);
    constructDTREE(nDTREEDepth);
  } //end of actions for training

  // continue with actions for evaluation
  if(!bTrain || nPercentForTest >= 1)
  {
    PassBackStatus(6,90);  
    ::PostMessage((HWND) pParam, WM_UPDATESCREEN,0,0);
    evaluate();
  }
  PassBackStatus(9,95);  
  ::PostMessage((HWND) pParam, WM_UPDATESCREEN,0,0);

  // WM_THREADFINISHED is user-defined message
  // sent to the view which calls the real handler in the doc
  PassBackStatus(12,100);
  ::PostMessage((HWND) pParam, WM_THREADFINISHED,0,0);
  cleanup();
  return 0;
}

void CALNfitDeepView::OnShowWindow(BOOL bShow, UINT nStatus) 
{
	CFormView::OnShowWindow(bShow, nStatus);
}

void CALNfitDeepView::CalcWindowRect(LPRECT lpClientRect, UINT nAdjustType) 
{
	CFormView::CalcWindowRect(lpClientRect, nAdjustType);
}

void CALNfitDeepView::OnFileSave() 
{
	CALNfitDeepDoc* pDoc = (CALNfitDeepDoc*)GetDocument();
  if(pDoc->m_bStarted)
  {
    MessageBox("Please wait until the run is finished to do this.");
  }
  else
  {
    UpdateData(TRUE);
    pDoc->m_nFit = m_nFit;  
    pDoc->m_nTrain = m_nTrain;
    pDoc->m_strDTREEFileName = m_strDTREEFileName; // update the document
    pDoc->m_strDataFileName = m_strDataFileName; // update the document
    pDoc->OnFileSave();
  }
}

void CALNfitDeepView::OnFileSaveAs() 
{
	CALNfitDeepDoc* pDoc = (CALNfitDeepDoc*)GetDocument();
  if(pDoc->m_bStarted)
  {
    MessageBox("Please wait until the run is finished to do this.");
  }
  else
  {
    UpdateData(TRUE);
    pDoc->m_nFit = m_nFit;  
    pDoc->m_nTrain = m_nTrain;
    pDoc->m_strDTREEFileName = m_strDTREEFileName; // update the document
    pDoc->m_strDataFileName = m_strDataFileName; // update the document
    pDoc->OnFileSaveAs();
  }
}


void CALNfitDeepView::OnButtonviewprev() 
{
	CALNfitDeepDoc* pDoc = (CALNfitDeepDoc*)GetDocument();
  if(pDoc->m_nDisplayedALNinput > 0) pDoc->m_nDisplayedALNinput--;
  pDoc->UpdateAllViews(NULL);
}

void CALNfitDeepView::OnButtonviewnext() 
{
	CALNfitDeepDoc* pDoc = (CALNfitDeepDoc*)GetDocument();
  if(pDoc->m_nDisplayedALNinput < pDoc->m_nALNinputs - 1)
  {
    pDoc->m_nDisplayedALNinput++;
    pDoc->UpdateAllViews(NULL);
  }
}

void CALNfitDeepView::OnSelchangeAlninputcolname() 
{
	CALNfitDeepDoc* pDoc = (CALNfitDeepDoc*)GetDocument();
  if(!pDoc->m_bStarted)
  {
  	CComboBox* pCLB = (CComboBox*) GetDlgItem(IDC_ALNINPUTCOLNAME);
    pDoc->m_nInputCol[pDoc->m_nDisplayedALNinput] =nInputCol[pDoc->m_nDisplayedALNinput]= pCLB->GetCurSel();
    pDoc->UpdateAllViews(NULL);
  }
  else
  {
    MessageBox("Please wait until the run is finished to do this.");
    pDoc->UpdateAllViews(NULL);
  }
}

void CALNfitDeepView::OnVScroll(UINT nSBCode, UINT nPos, CScrollBar* pScrollBar) 
{
  if(nSBCode == SB_ENDSCROLL)
  {
    return; // reject spurious messages that could cause doubled message boxes
  }
  CALNfitDeepDoc* pDoc = (CALNfitDeepDoc*)GetDocument();
  CSpinButtonCtrl* pSpin;
  if(!pDoc->m_bStarted)
  {
    if(pScrollBar == NULL)
    {
      CScrollView::OnVScroll(nSBCode, nPos, NULL);
    }
    else
    {
      if((pScrollBar->GetDlgCtrlID() == IDC_SPINLAG)&&
         (pDoc->m_nDisplayedALNinput != m_nALNinputs - 1))
      {
	      pDoc->m_nLag[pDoc->m_nDisplayedALNinput] = nLag[pDoc->m_nDisplayedALNinput] = nPos;
      }
      else
      if(pScrollBar->GetDlgCtrlID() == IDC_SPINALNIOPROP)
      {
        if(pDoc->m_nDisplayedALNinput < pDoc->m_nALNinputs - 1)
        {
          if(nPos == OUTPUT)
          {
            nPos = INPUTINC;
            MessageBox("Only the rightmost variable can be the output");
            pSpin = (CSpinButtonCtrl* ) GetDlgItem(IDC_SPINALNIOPROP);
	          pSpin->SetPos(INPUTINC);
          }
          else
          {
	          pDoc->m_nIOpropMessageNumber[pDoc->m_nDisplayedALNinput] = nPos;
            if(nPos == INPUT)
            {
              if(pDoc->m_dblMaxWeight[pDoc->m_nDisplayedALNinput] <= 0)
                pDoc->m_dblMaxWeight[pDoc->m_nDisplayedALNinput]  = dblMax;
                dblMaxWeight[pDoc->m_nDisplayedALNinput]  = dblMax;
              if(pDoc->m_dblMinWeight[pDoc->m_nDisplayedALNinput] >=0)
                pDoc->m_dblMinWeight[pDoc->m_nDisplayedALNinput] = -dblMax;
                dblMinWeight[pDoc->m_nDisplayedALNinput] = -dblMax;
            }
            else if(nPos ==INPUTDEC)
            {
              if(pDoc->m_dblMaxWeight[pDoc->m_nDisplayedALNinput] > 0)
                pDoc->m_dblMaxWeight[pDoc->m_nDisplayedALNinput] = 0;
                dblMaxWeight[pDoc->m_nDisplayedALNinput] = 0;
              if(pDoc->m_dblMinWeight[pDoc->m_nDisplayedALNinput] > 0)
                pDoc->m_dblMinWeight[pDoc->m_nDisplayedALNinput]= -dblMax;
                dblMinWeight[pDoc->m_nDisplayedALNinput]= -dblMax;
            }
            else if(nPos == INPUTINC)
            {
              if(pDoc->m_dblMaxWeight[pDoc->m_nDisplayedALNinput] < 0)
                pDoc->m_dblMaxWeight[pDoc->m_nDisplayedALNinput]= dblMax;
                dblMaxWeight[pDoc->m_nDisplayedALNinput]= dblMax;
              if(pDoc->m_dblMinWeight[pDoc->m_nDisplayedALNinput] < 0)
                pDoc->m_dblMinWeight[pDoc->m_nDisplayedALNinput]= 0;
                dblMinWeight[pDoc->m_nDisplayedALNinput]= 0;
            }
          }
        }
        else // the displayed ALNinput is the output
        {
          if(nPos != OUTPUT)
          {
            nPos = OUTPUT;
            MessageBox("The rightmost variable must be the output");
            pSpin = (CSpinButtonCtrl* ) GetDlgItem(IDC_SPINALNIOPROP);
	          pSpin->SetPos(OUTPUT);
          }
          pDoc->m_dblMaxWeight[pDoc->m_nDisplayedALNinput] = -1;
          pDoc->m_dblMinWeight[pDoc->m_nDisplayedALNinput] = -1;
          dblMaxWeight[pDoc->m_nDisplayedALNinput] = -1;
          dblMinWeight[pDoc->m_nDisplayedALNinput] = -1;
        }
      }
      else
      if(pScrollBar->GetDlgCtrlID() == IDC_SPIN_NALNS)
      {
        pDoc->m_nALNs = nALNs = nPos;
      }
    }
  }
  else // not started yet
  {
    MessageBox("Please wait until the run is finished to do this.");
  }
  pDoc->UpdateAllViews(NULL);
}

void CALNfitDeepView::OnButtonremovealninput() 
{
	CALNfitDeepDoc* pDoc = (CALNfitDeepDoc*)GetDocument();
  if(!pDoc->m_bStarted)
  {
    if(pDoc->m_nDisplayedALNinput == pDoc->m_nALNinputs - 1)
    {
		  MessageBox("You can't remove the output. Instead you can change the column to a new one which doesn't already appear as an input with lag 0.");
      return;
    }
    if(pDoc->m_nALNinputs == 2)
    {
		  MessageBox("You can't eliminate the last input variable.  However you can edit the column and lag.");
    }
    else
    {
      pDoc->RemoveInput();
    }
  }
  else
  {
    MessageBox("Please wait until the run is finished to do this.");
  }
}

void CALNfitDeepView::OnButtonaddalninput() 
{
	CALNfitDeepDoc* pDoc = (CALNfitDeepDoc*)GetDocument();
  if(!pDoc->m_bStarted)
  {
    if(pDoc->m_nALNinputs == 101)
    {
      MessageBox("You can't add an input. You already have 101 ALN inputs (including one output), which is the upper limit.");
      return;
    }
    if(pDoc->m_nDisplayedALNinput == pDoc->m_nALNinputs - 1)
    {
		  MessageBox("You can add an input which has the same column as the output, but the lag can't be 0. Please continue...");
    }
    pDoc->AddInput();
  }
  else
  {
    MessageBox("Please wait until the run is finished to do this.");
  }
}

void CALNfitDeepView::OnPrint(CDC* pDC, CPrintInfo* pInfo) 
{
  CALNfitDeepDoc* pDoc = (CALNfitDeepDoc*)GetDocument();
  pDC->SetMapMode(MM_TWIPS);
  m_nPage = pInfo->m_nCurPage;
  CFont font;
  TEXTMETRIC tm;
  int nHeight;
  font.CreateFont(-200,0,0,0,400,FALSE,FALSE,0,ANSI_CHARSET,
    OUT_DEFAULT_PRECIS,CLIP_DEFAULT_PRECIS,
    DEFAULT_QUALITY,DEFAULT_PITCH | FF_ROMAN,"Times New Roman");
  CFont* pOldFont = (CFont*) pDC->SelectObject(&font);
  pDC->GetTextMetrics(&tm);
  nHeight = tm.tmHeight + tm.tmExternalLeading;
  pDoc->GeneratePrint();  // this generates the printed page no matter what
  int nStringArraySize = pDoc->m_stringArray.GetSize();
  int nStart, nEnd;
  nStart = (m_nPage - 1) * CALNfitDeepDoc::nLinesPerPage;
  nEnd = nStart + CALNfitDeepDoc::nLinesPerPage;
  PrintPageHeader(pDC);
  if(nStringArraySize < nEnd) nEnd = nStringArraySize;
  for(int i = nStart; i < nEnd; i++)
  {
    pDC->TextOut(720, -(i-nStart +2) * nHeight -720, pDoc->m_stringArray[i]);
  }
  pDC->SelectObject(pOldFont);
}


void CALNfitDeepView::PrintPageHeader(CDC* pDC)
{
  CALNfitDeepDoc* pDoc = (CALNfitDeepDoc*)GetDocument();
  CString Temp = "Parameters used for " + pDoc->m_strProtocolFileName;;
  CString Str;
  Str.Format("   Page %d",m_nPage);
  if(m_nPage > 1) Temp = Temp + Str;
  pDC->TextOut(720, -720, Temp);
}

void CALNfitDeepView::OnRadiosortimportance() 
{
  CALNfitDeepDoc* pDoc = (CALNfitDeepDoc*)GetDocument();
  // sort according to importance
  if(pDoc->m_bStarted)
  {
    MessageBox("Please wait until the run is finished to do this.");
    m_nSortCols = pDoc->m_bSortCols?0:1;
    UpdateControlsFromDoc();
  }
  else
  {
    m_nSortCols = 1;
    pDoc->m_bSortCols = FALSE;
    pDoc->OnSort();
  }	
}

void CALNfitDeepView::OnRadiosortcols() 
{
  CALNfitDeepDoc* pDoc = (CALNfitDeepDoc*)GetDocument();
  if(pDoc->m_bStarted)
  {
    MessageBox("Please wait until the run is finished to do this.");
    m_nSortCols = pDoc->m_bSortCols?0:1;
    UpdateControlsFromDoc();
  }
  else
  {
    // sort first according to the order of columns in the original data, then by increasing lag
    // set flag to warn of two column/lag combinations which are equal,
    // or an input which is the same as the output column with lag 0
    pDoc->m_bWarnEqual = FALSE;
    pDoc->m_bSortCols = TRUE;
    pDoc->OnSort();
    if(pDoc->m_bWarnEqual == TRUE) MessageBox("Warning: duplicate ALN inputs/output.  Please check!");
  }
}




void CALNfitDeepView::OnRadioFit() 
{
  CALNfitDeepDoc* pDoc = (CALNfitDeepDoc*)GetDocument();
  if(pDoc->m_bStarted)
  {
    MessageBox("Please wait until the run is finished to do this.");
    m_nFit = pDoc->m_nFit;
	  UpdateData(FALSE); // calls DDX
  }
  else
  {
    m_nFit = 0;
    pDoc->m_nFit = m_nFit;
  }
}

void CALNfitDeepView::OnRadioClass() 
{
  CALNfitDeepDoc* pDoc = (CALNfitDeepDoc*)GetDocument();
  if(pDoc->m_bStarted)
  {
    MessageBox("Please wait until the run is finished to do this.");
    m_nFit = pDoc->m_nFit;
	  UpdateData(FALSE); // calls DDX
  }
  else
  {
    m_nFit = 1;
    pDoc->m_nFit = m_nFit;
  }
}

void CALNfitDeepView::OnRadioTrain() 
{
  CALNfitDeepDoc* pDoc = (CALNfitDeepDoc*)GetDocument();
  if(pDoc->m_bStarted)
  {
    MessageBox("Please wait until the run is finished to do this.");
    m_nTrain = pDoc->m_nTrain;
	  UpdateData(FALSE); // calls DDX
  }
  else
  {
    m_nTrain = 0;
    pDoc->m_nTrain = m_nTrain;
  }
}

void CALNfitDeepView::OnRadioEval() 
{
  CALNfitDeepDoc* pDoc = (CALNfitDeepDoc*)GetDocument();
  if(pDoc->m_bStarted)
  {
    MessageBox("Please wait until the run is finished to do this.");
    m_nTrain = pDoc->m_nTrain;
	  UpdateData(FALSE); // calls DDX
  }
  else
  {
    m_nTrain = 1;
    pDoc->m_nTrain = m_nTrain;
  }
}

void CALNfitDeepView::OnEditClearall() 
{
  CALNfitDeepDoc* pDoc = (CALNfitDeepDoc*)GetDocument();
  if(pDoc->m_bStarted)
  {
    MessageBox("The run has started and you can't clear all the parameters!");
  }
  else
  {
    pDoc->EditClearall();
  }	
}

void CALNfitDeepView::OnUpdateEditmaxwt() 
{
  CString strWeight;
  double value;
  CALNfitDeepDoc* pDoc = (CALNfitDeepDoc*)GetDocument();
  if(pDoc->m_bStarted)
  {
    MessageBox("Please wait until the run is finished to do this.");
    value = pDoc->m_dblMaxWeight[pDoc->m_nDisplayedALNinput];
    m_strMaxWeight.Format("%fl",value);
	  UpdateData(FALSE); // calls DDX
  }
  else
  {
    UpdateData(TRUE);
    double value = atof(LPCTSTR(m_strMaxWeight));
    if(pDoc->m_nDisplayedALNinput != m_nALNinputs - 1)
    {
      pDoc->m_dblMaxWeight[pDoc->m_nDisplayedALNinput] = dblMaxWeight[pDoc->m_nDisplayedALNinput] =  value;
    }
    else
    {
      MessageBox("The weight on the output is always -1.");
      pDoc->m_dblMaxWeight[pDoc->m_nDisplayedALNinput] = dblMaxWeight[pDoc->m_nDisplayedALNinput] = -1;
      pDoc->m_dblMinWeight[pDoc->m_nDisplayedALNinput] = dblMinWeight[pDoc->m_nDisplayedALNinput] = -1;
    }
  }
}

void CALNfitDeepView::OnKillfocusEditmaxwt() 
{
  // here we set the monotonicity property to INPUTDEC if the upper bound on weight is < 0
  CString strWeight;
  double value;
  UpdateData(TRUE);
  value = atof(LPCTSTR(m_strMaxWeight));  // value is 0.0 if conversion is not possible
  CALNfitDeepDoc* pDoc = (CALNfitDeepDoc*)GetDocument();
  if(pDoc->m_nDisplayedALNinput != m_nALNinputs - 1)
  {
    if(value < dblMinWeight[pDoc->m_nDisplayedALNinput])
    {
      pDoc->m_dblMaxWeight[pDoc->m_nDisplayedALNinput] = dblMaxWeight[pDoc->m_nDisplayedALNinput] =  dblMax;
      MessageBox("You can't set the upper bound on a weight to be less than the lower bound! Upper bound has been reset to maximum.");
    }
    else
    {
      if(value < 0) pDoc->m_nIOpropMessageNumber[pDoc->m_nDisplayedALNinput] = INPUTDEC;
    }
  }
  UpdateControlsFromDoc();
}

void CALNfitDeepView::OnUpdateEditminwt() 
{
  CString strWeight;
  double value;
  CALNfitDeepDoc* pDoc = (CALNfitDeepDoc*)GetDocument();
  if(pDoc->m_bStarted)
  {
    MessageBox("Please wait until the run is finished to do this.");
    value = pDoc->m_dblMinWeight[pDoc->m_nDisplayedALNinput];
    m_strMinWeight.Format("%fl",value);
	  UpdateData(FALSE); // calls DDX
  }
  else
  {
    UpdateData(TRUE);
    double value = atof(LPCTSTR(m_strMinWeight));  // a bad string will be converted to 0.0
    if(pDoc->m_nDisplayedALNinput != m_nALNinputs - 1)
    {
      pDoc->m_dblMinWeight[pDoc->m_nDisplayedALNinput] = dblMinWeight[pDoc->m_nDisplayedALNinput] =  value;
    }
    else
    {
      MessageBox("The weight on the output is always -1.");
      pDoc->m_dblMaxWeight[pDoc->m_nDisplayedALNinput] = dblMaxWeight[pDoc->m_nDisplayedALNinput] = -1;
      pDoc->m_dblMinWeight[pDoc->m_nDisplayedALNinput] = dblMinWeight[pDoc->m_nDisplayedALNinput] = -1;
    }
  }
}

void CALNfitDeepView::OnKillfocusEditminwt() 
{
  // here we set the monotonicity property to INPUTDEC if the upper bound on weight is < 0
  CString strWeight;
  double value;
  UpdateData(TRUE);
  value = atof(LPCTSTR(m_strMinWeight));  // value is 0.0 if conversion is not possible
  CALNfitDeepDoc* pDoc = (CALNfitDeepDoc*)GetDocument();
  if(pDoc->m_nDisplayedALNinput != m_nALNinputs - 1)
  {
    if(value > dblMaxWeight[pDoc->m_nDisplayedALNinput])
    {
      pDoc->m_dblMinWeight[pDoc->m_nDisplayedALNinput] = dblMinWeight[pDoc->m_nDisplayedALNinput] =  -dblMax;
      MessageBox("You can't set the lower bound on a weight to be greater than the upper bound! The lower bound has been reset to minimum.");
    }
    else
    {
      if(value > 0) pDoc->m_nIOpropMessageNumber[pDoc->m_nDisplayedALNinput] = INPUTINC;
    }
  }
  UpdateControlsFromDoc();
}


void CALNfitDeepView::OnButtonoptions()
{
	CALNfitDeepDoc* pDoc = (CALNfitDeepDoc*)GetDocument();
	if (pDoc->m_bStarted || (fpFileSetupProtocol == NULL))
	{
		MessageBox("This must be done after setting up the Data file and before Start");
	}
	else
	{
		COptionsDialog dlg;
		// update all options from the doc
		dlg.m_nPercentForTest = pDoc->m_nPercentForTest;;
		dlg.m_nEstimateRMSError = pDoc->m_nEstimateRMSError;
		char szSetTolerance[128];
		sprintf(szSetTolerance, "%11.7lf", pDoc->m_dblSetTolerance);
		dlg.m_strSetTolerance = szSetTolerance;
		dlg.m_nZeroSmoothing = pDoc->m_nZeroSmoothing;
		dlg.m_nNoJitter = pDoc->m_nNoJitter;
		dlg.m_nDTREEDepth = pDoc->m_nDTREEDepth;
		dlg.m_bTimePrefixes = pDoc->m_bTimePrefixes;
		dlg.m_bReplaceUndefined = pDoc->m_bReplaceUndefined;
		dlg.m_bDiagnostics = pDoc->m_bDiagnostics;
		UpdateData(FALSE);
		if (dlg.DoModal() == IDOK)
		{
			// update the doc from the options, including some integers for options not selected
			UpdateData(TRUE);

			if ((dlg.m_nPercentForTest < 0) || (dlg.m_nPercentForTest > 50))
			{
				dlg.m_nPercentForTest = pDoc->m_nPercentForTest;
				MessageBox("The percentage used for testing must be between 0 and 50.  Value left unchanged.");
			}
			else
			{
				pDoc->m_nPercentForTest = nPercentForTest = dlg.m_nPercentForTest;
			}
			dblFracTest = (double)dlg.m_nPercentForTest / (double)100;
			pDoc->m_nEstimateRMSError = dlg.m_nEstimateRMSError;

			if (dlg.m_nEstimateRMSError == 1)
			{
				// we have highlighted the radiobutton to set the tolerance without
				// using a variance set, we convert the field to a double and
				// transmit that to the doc and its global variable
				pDoc->m_dblSetTolerance = dblSetTolerance = atof(LPCTSTR(dlg.m_strSetTolerance));
				bEstimateNoiseVariance = FALSE;
			}
			else
			{
				// the tolerance will be set by estimating RMS noise
				bEstimateNoiseVariance = TRUE;
			}
			
			pDoc->m_nZeroSmoothing = dlg.m_nZeroSmoothing;
			pDoc->m_nNoJitter = dlg.m_nNoJitter;
			bJitter = (dlg.m_nNoJitter == 0) ? FALSE : TRUE;
			pDoc->m_nDTREEDepth = nDTREEDepth = dlg.m_nDTREEDepth;
			pDoc->m_bReplaceUndefined =	bReplaceUndefined = dlg.m_bReplaceUndefined;
			pDoc->m_bTimePrefixes = bTimePrefixes = dlg.m_bTimePrefixes;
			pDoc->m_bDiagnostics = bDiagnostics = dlg.m_bDiagnostics;
			if (!bPrint) // where is bPrint defined??  let's leave it as TRUE ( it is not the same as bDiagnostics)
			{
				fclose(fpFileSetupProtocol);
			}
			UpdateData(FALSE);
		}
		else
		{
			// operation cancelled: update the dialog members, which may have changed, from the doc
			dlg.m_nPercentForTest = pDoc->m_nPercentForTest;
			dlg.m_nEstimateRMSError = pDoc->m_nEstimateRMSError;
			if (pDoc->m_nEstimateRMSError == 0) pDoc->m_dblSetTolerance = 0.1905255; // reset if we are estimating RMS error
			sprintf(szSetTolerance, "%11.7lf", pDoc->m_dblSetTolerance);
			dlg.m_strSetTolerance = szSetTolerance;
			dlg.m_nZeroSmoothing = pDoc->m_nZeroSmoothing;
			dlg.m_nNoJitter = pDoc->m_nNoJitter;
			dlg.m_nDTREEDepth = pDoc->m_nDTREEDepth;
			dlg.m_bTimePrefixes = pDoc->m_bTimePrefixes;
			dlg.m_bReplaceUndefined = pDoc->m_bReplaceUndefined;
			dlg.m_bDiagnostics = pDoc->m_bDiagnostics;
			UpdateData(FALSE);
		}
	}
}


void CALNfitDeepView::OnEnChangeDatafilecolumns()
{
	// TODO:  If this is a RICHEDIT control, the control will not
	// send this notification unless you override the CFormView::OnInitDialog()
	// function and call CRichEditCtrl().SetEventMask()
	// with the ENM_CHANGE flag ORed into the mask.

	// TODO:  Add your control notification handler code here
}
