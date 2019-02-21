// ALNfitDeepDoc.cpp : implementation of the CALNfitDeepDoc class
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
#include "datafile.h"
#include "alnextern.h"
//#include "alnintern.h"
#include "ALNfitDeep.h"
#include "ALNfitDeepDoc.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

#define		MESSAGE0  "Analyzing the data file"
#define		MESSAGE1  "Doing linear regression to find some useful starting values."
#define		MESSAGE2  "Creating noise variance samples; training an ALN on the noise variance."
#define		MESSAGE3  "Training ALNs using the noise variance function to prevent overtraining."
#define		MESSAGE4  "Bagging: training an ALN on samples of the average of several trained ALNs"
#define		MESSAGE5  "Constructing a DTREE from the average ALN and writing the .dtr file "
#define		MESSAGE6  "Loading and evaluating the DTREE on the test data file"
#define		MESSAGE7  "Determining importance of input variables using the trained ALNs"
#define		MESSAGE8  "Please enter parameters for a run above; then click on Start"
#define		MESSAGE9  "Finished. Please examine TrainProtocol and Output text files for results."
#define   MESSAGE10  "Opened.fit file."
#define   MESSAGE11  "Please save the current parameters to a .fit file, choose new ones, or exit."
#define   MESSAGE12  "Opened setup file. IMPORTANT: Now browse to a data file like the one above."
#define   MESSAGE13  "ERROR: the message number is out of range"

typedef struct
{
  int nIndex;
  int nInputCol;
  int nLag;
  int nIOpropMessageNumber;
  double dblMinWeight;
  double dblMaxWeight;
  double dblImportance;
} SRTND;  
SRTND SortArray[101];        // relates the index of the ALN input to its index in sorted order
int  SNCompareCols(void * p1, void * p2); // sort node comparison operation
int  SNCompareImportance(void * p1, void * p2); // sort node comparison operation
BOOL bWarnEqual = FALSE;

/////////////////////////////////////////////////////////////////////////////
// CALNfitDeepDoc

IMPLEMENT_DYNCREATE(CALNfitDeepDoc, CDocument)

BEGIN_MESSAGE_MAP(CALNfitDeepDoc, CDocument)
	//{{AFX_MSG_MAP(CALNfitDeepDoc)
	ON_COMMAND(ID_FILE_SAVE, OnFileSave)
	ON_COMMAND(ID_FILE_SAVE_AS, OnFileSaveAs)
	//}}AFX_MSG_MAP
END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////////
// CALNfitDeepDoc construction/destruction

CALNfitDeepDoc::CALNfitDeepDoc()
{
	m_pDocument = this;
}

CALNfitDeepDoc::~CALNfitDeepDoc()
{
}

void CALNfitDeepDoc::GetColumnNames(int ncols)
{
}

void CALNfitDeepDoc::GenerateReportString()
{
  // takes the message number in the doc's m_nMessage, gets the text,
  // and transfers it to the report string
	const int nMessageCount = 14;
	const char * apMessages[nMessageCount] = 
	{
		MESSAGE0,  //"Analyzing data files",
		MESSAGE1,  //"Using linear regression to support further training",
		MESSAGE2,  //"Creating noise variance samples and training an ALN on their logarithms",
		MESSAGE3,  //"Training one or more ALNs using the noise variance samples to guide ALN growth",
		MESSAGE4,  //"Training a new ALN with samples of the average of the trained ALNs",
		MESSAGE5,  //"Constructing a DTREE from the average ALN and writing the .dtr file ",
		MESSAGE6,  //"Loading and evaluating the DTREE on the test data file",
		MESSAGE7,  //"Determining importance of input variables using TV file and the trained ALNs",
		MESSAGE8,  //"Please enter parameters for a run above; then click on Start",
		MESSAGE9,  //"Finished. Please examine ..Protocol.. and ..Output.. text files for results.",
    MESSAGE10,  //"Opened .fit file",
    MESSAGE11,  //"Please save the current parameters to a .fit file, choose new ones, or exit.",
    MESSAGE12,  //"Opened setup file. IMPORTANT: Now browse to a data file like the one above."
    MESSAGE13   //"ERROR: the message number is out of range!"

  };
  if(m_nMessage < 0 || m_nMessage > 12)
  {
    m_nMessage = 13; // this is the error message
  }
  // get the string for the message
	m_strReport = apMessages[m_nMessage];
}

BOOL CALNfitDeepDoc::OnNewDocument()
{
	if (!CDocument::OnNewDocument())
		return FALSE;
 
  // (SDI documents will reuse this document)
  m_strDTREEFileName = "";
	m_strDataFileName = "";
	m_nFit = 0;   // 0 indicates regression, 1 classification
	m_nTrain = 0; // 0 indicates training, 1 evaluation
	m_nALNs = 1;  // default
	m_nMessage = 8; // asks to set up parameters then click on Start
  GenerateReportString();
	CString m_strProtocolFileName = "";   // set up by start button
	CString m_strScatterFileName = "";     // set up by start button
  m_bSortCols = TRUE;
	m_nPercentProgress = 0; // integer value for progress indicator
  m_nColsUniv = nColsUniv =101;
  m_nALNinputs = nALNinputs =101;
  m_nDisplayedALNinput = 100;
  for(int nn = 0; nn < m_nColsUniv; nn++) // for the training default m_nColsUniv equals m_nALNinputs
  {
    char buff[32] = "__A";
    //_itoa(nn + 1,buff+7,10);
    if(nn/26>0)buff[1] = 'A'- 1 + nn/26;
    buff[2] = 'A' + nn%26;
    buff[3] = '\0';
    m_strVarname[nn] = buff;
    for(int charno = 0; charno < 32; charno++)
    {
      varname[nn][charno]=buff[charno];
    }
    m_nInputCol[nn] = nInputCol[nn] = nn;
    m_nLag[nn] = nLag[nn] = 0;
    m_nIOpropMessageNumber[nn] = INPUT;
    if(nn < m_nColsUniv -1)
    {
      dblMinWeight[nn] = m_dblMinWeight[nn] = -dblMax;
      dblMaxWeight[nn] = m_dblMaxWeight[nn] = dblMax;
    }
    else
    {
      dblMinWeight[nn] = m_dblMinWeight[nn] = -1;
      dblMaxWeight[nn] = m_dblMaxWeight[nn] = -1;
    }
    m_dblImportance[nn] = dblImportance[nn] = 0;
  }
  if(m_nTrain == 0) // if this is a training run, we show the output column as such
  {
    m_nIOpropMessageNumber[m_nALNinputs-1] = OUTPUT;
  }
  m_bStarted = FALSE;
  // Set all option values
  m_nPercentForTest = 10;
  m_nEstimateRMSError = 0;   // 0 = TRUE; we overfit an ALN on some data and compare to a separate set of data points
  m_dblSetTolerance = 0.1905255;  // this value is good for classification when two classes are one unit apart
	m_nZeroSmoothing = 0;
	m_nNoJitter = 0;
  m_nDTREEDepth = 1;
	m_bTimePrefixes = TRUE;
  m_bReplaceUndefined = FALSE;
  m_bDiagnostics = FALSE;
  UpdateAllViews(NULL);
	return TRUE;
}

void CALNfitDeepDoc::RemoveInput()
{
  ASSERT(nALNinputs == m_nALNinputs);
  ASSERT(m_nALNinputs > 2);
  ASSERT(m_nDisplayedALNinput < m_nALNinputs -1);
  // shift entries above the displayed one down one step
  for(int i = m_nDisplayedALNinput; i < m_nALNinputs - 1; i++)
  {
    m_nInputCol[i] = nInputCol[i] =m_nInputCol[i+1]; 
    m_nLag[i] = nLag[i] = m_nLag[i+1];
    m_nIOpropMessageNumber[i] = m_nIOpropMessageNumber[i+1];
    m_dblImportance[i] = dblImportance[i] = m_dblImportance[i+1];
  }
  m_nALNinputs = --nALNinputs;
  if(m_nDisplayedALNinput == m_nALNinputs) m_nDisplayedALNinput--;
  UpdateAllViews(NULL);
}

void CALNfitDeepDoc::AddInput()
{
  //pull the entries above the displayed one up by one starting at the new output
  for(int i = m_nALNinputs; i > m_nDisplayedALNinput; i--)
  {
    m_nInputCol[i] = nInputCol[i] =m_nInputCol[i-1]; 
    m_nLag[i] = nLag[i] = m_nLag[i-1];
    m_nIOpropMessageNumber[i] = m_nIOpropMessageNumber[i-1];
    m_dblImportance[i] = dblImportance[i] = m_dblImportance[i-1];
    m_dblMinWeight[i] = dblMinWeight[i] = m_dblMinWeight[i-1];
    m_dblMaxWeight[i] = dblMaxWeight[i] = m_dblMaxWeight[i-1];
  }
  if(m_nDisplayedALNinput < m_nALNinputs -1)
  {
    m_nDisplayedALNinput++; // unless we were displaying the output,display the next higher one
  }
  // increment the number of inputs
  m_nALNinputs = ++nALNinputs;
  // display the input to be edited
  // m_nInputCol[m_nDisplayedALNinput] left unchanged, lag is increased by 1
  m_nLag[m_nDisplayedALNinput] = nLag[m_nDisplayedALNinput] = m_nLag[m_nDisplayedALNinput] +1;
  m_nIOpropMessageNumber[m_nDisplayedALNinput] = INPUT;
  m_dblImportance[m_nDisplayedALNinput] = dblImportance[m_nDisplayedALNinput] = 0;
  dblMinWeight[m_nDisplayedALNinput] = m_dblMinWeight[m_nDisplayedALNinput] = -dblMax;
  dblMaxWeight[m_nDisplayedALNinput] = m_dblMaxWeight[m_nDisplayedALNinput] = dblMax;
  UpdateAllViews(NULL);
}

void CALNfitDeepDoc::DeleteContents() 
{
  // in an SDI application, the document object is reused
	m_strDTREEFileName = "";
	m_strDataFileName = "";
	CString m_strProtocolFileName = "";   // set up by start button
	CString m_strScatterFileName = "";     // set up by start button
	m_strReport = "";
	m_nFit = 0;   // 0 indicates regression, 1 classification
	m_nTrain = 0; // 0 indicates training, 1 evaluation
	m_nALNs = 1;  // default
	m_nMessage = 0; // asks to set up parameters then click on Start
	m_nPercentProgress = 0; // integer value for progress indicator
  m_nColsUniv = 101;
  for(int nn = 0; nn < m_nColsUniv; nn++)
  {
    char buff[32] = "__A";
    //_itoa(nn + 1,buff+7,10);
    if(nn/26>0)buff[1] = 'A'- 1 + nn/26;
    buff[2] = 'A' + nn%26;
    buff[3] = '\0';
    m_strVarname[nn] = buff;
    m_nInputCol[nn] = nn;
    m_nLag[nn] = 0;
    m_dblImportance[nn] = 0;
  }
  m_stringArray.RemoveAll();
}

/////////////////////////////////////////////////////////////////////////////
// CALNfitDeepDoc serialization

void CALNfitDeepDoc::Serialize(CArchive& ar)
{
	if (ar.IsStoring())
	{
		ar << m_strDTREEFileName;
		ar <<	m_strDataFileName;
		ar <<	m_nFit;
		ar <<	m_nTrain;
		ar <<	m_nALNs;
    ar << m_nColsUniv;
    for(int ii = 0; ii < m_nColsUniv; ii++)
    {
      ar << m_strVarname[ii];
    }
    ar << m_nALNinputs;
    for(int jj = 0; jj < m_nALNinputs; jj++)
    {
      ar << m_nInputCol[jj];
      ar << m_nLag[jj];
      ar << m_nIOpropMessageNumber[jj];
      ar << m_dblMaxWeight[jj];
      ar << m_dblMinWeight[jj];
      ar << m_dblImportance[jj];
    }
    ar << m_bReplaceUndefined;
    ar << m_bTimePrefixes;
    ar << m_bDiagnostics;
    ar << m_nNoJitter;
    ar << m_nPercentForTest;
    ar << m_nSamplesForTest;
    ar << m_nZeroSmoothing;
    ar << m_nEstimateRMSError;
	}
	else
	{
		// this is executed after OnFileOpen from File menu etc.
    // it changes the doc and must be followed by UpdateControlsFromDoc
		ar >> m_strDTREEFileName;
    CString strTemp;
		ar >>	m_strDataFileName;
		ar >>	m_nFit;
		ar >>	m_nTrain;
    m_nTrain = 1; // always 1, i.e. no training, when a fit file is read in; we are evaluating
                  // this directs a lot of the processing so Browse Data and .fit work together
		ar >>	m_nALNs;
    ar >> m_nColsUniv;
		BOOL bColumnsSame = (nColsUniv == m_nColsUniv)?TRUE:FALSE;
    for(int ii = 0; ii < m_nColsUniv; ii++)
    {
      ar >> m_strVarname[ii];  // now the doc has the new variable name
      if(bColumnsSame)
      {
        strTemp = varname[ii]; // retains the old variable name from the document after contents deleted.
        if(strTemp != m_strVarname[ii]) bColumnsSame = FALSE;
      }
    }
    ar >> m_nALNinputs;
    for(int jj = 0; jj < m_nALNinputs; jj++)
    {
      ar >> m_nInputCol[jj];
      nInputCol[jj] = m_nInputCol[jj];
      ar >> m_nLag[jj];
      nLag[jj] = m_nLag[jj];
      ar >> m_nIOpropMessageNumber[jj];
      ar >> m_dblMaxWeight[jj];
      dblMaxWeight[jj] = m_dblMaxWeight[jj];
      ar >> m_dblMinWeight[jj];
      dblMinWeight[jj] = m_dblMinWeight[jj];
      ar >> m_dblImportance[jj];
      dblImportance[jj] = m_dblImportance[jj];
    }
    // m_nDisplayedALNinput not stored, set to rightmost entry (the output)
    m_nDisplayedALNinput = m_nALNinputs - 1;
    CString m_strProtocolFileName = "";   // set up by start button
	  CString m_strScatterFileName = "";     // set up by start button
		m_nMessage = 10;
    if(!bColumnsSame)m_nMessage = 12; // report a mismatch of names this way
    GenerateReportString();
    m_nPercentProgress = 0;
    m_bReplaceUndefined = FALSE;
    m_bTimePrefixes = TRUE;
    m_bDiagnostics = FALSE;
    m_nNoJitter = 0;
    m_nPercentForTest = 10;
    m_nZeroSmoothing = 0;
    m_nEstimateRMSError = 0;
    ar >> m_bReplaceUndefined;
    ar >> m_bTimePrefixes;
    ar >> m_bDiagnostics;
    ar >> m_nNoJitter;
    ar >> m_nPercentForTest;
    ar >> m_nSamplesForTest;
    ar >> m_nZeroSmoothing;
    ar >> m_nEstimateRMSError;
	}
  UpdateAllViews(NULL);
}

void CALNfitDeepDoc::GeneratePrint()
{
  m_stringArray.RemoveAll();
  m_stringArray.SetSize( 150, -1);
  CString Temp,Temp1,Temp2,Temp3,Temp4,Temp5,Temp6;
  m_stringArray[0] = "ALNfitDeep Parameters";
  CTime time;
  time = CTime::GetCurrentTime();
  CString timestr = time.Format( "%A, %B %d, %Y at %H:%M" );
  m_stringArray[1] ="Printed on " + timestr; 
  m_stringArray[2] ="DTREE file:    " + m_strDTREEFileName;
  m_stringArray[3] ="Data file:     " + m_strDataFileName;
  m_stringArray[4] ="Protocol file: " + m_strProtocolFileName;
  m_stringArray[5] ="Output file:   " + m_strScatterFileName;
  if(m_nFit == 0)
  {
    Temp = "The task is regression";
  }
  else
  {
    Temp = "The task is classification";
  }
  m_stringArray[6] = Temp;
  if(m_nTrain == 0)
  {
    Temp = "Training is to produce a new DTREE solution";
  }
  else
  {
    Temp = "An existing DTREE solution will be used on the data file";
  }
  m_stringArray[7] = Temp;
  char string[32];
  _itoa(m_nALNs,string,10);
  Temp = string;
  m_stringArray[8] = "The number of ALNs used to form the bagged average is " + Temp;
  _itoa(m_nColsUniv, string,10);
  Temp = string;
  m_stringArray[9] = "The number of columns in the data file is " + Temp;
  int k=10;
  int nlimit;
  for(int ii = 0; ii < m_nColsUniv; ii=ii+5)
  {
    m_stringArray[k] = "    " + m_strVarname[ii];
    nlimit = ii + 5;
    if(m_nColsUniv < nlimit) nlimit = m_nColsUniv;
    for(int jj = ii+1; jj < nlimit; jj++)
    {
       m_stringArray[k] = m_stringArray[k] + "  " + m_strVarname[jj];
    }
    k++;
  }
  _itoa(m_nALNinputs,string,10);
  Temp = string;
  m_stringArray[k] = "The number of ALN inputs is " + Temp;
  k++;
  m_stringArray[k] = "The ALN input column, lag, IO message, min weight, max weight, importance:";
  k++;
  CString strIOpropMessageArray[4] ={
    "input  --  monotonic decreasing",
    "input  --  no monotonicity constraint",
    "input  --  monotonic increasing",
    "output --  weight is always -1 and importance is always 0"};

  for(int jj = 0; jj < m_nALNinputs; jj++)
  {
    _itoa(jj+1,string,10);
    Temp = string;
    Temp1 = m_strVarname[m_nInputCol[jj]];
    _itoa(m_nLag[jj],string,10);
    Temp2 = string;
    Temp3 = strIOpropMessageArray[m_nIOpropMessageNumber[jj]];
    Temp4.Format((LPCTSTR)"%lf",m_dblMinWeight[jj]);
    Temp5.Format((LPCTSTR)"%lf",m_dblMaxWeight[jj]);
    Temp6.Format((LPCTSTR)"%lf",m_dblImportance[jj]);
    m_stringArray[k] = "   Input " + Temp + ": col= " + Temp1 + "   lag= "+ Temp2 + "   msg= "+  Temp3;
    k++;
    if(jj < m_nALNinputs -1)
    {
      m_stringArray[k] = "      min wt = " + Temp4 + "   max wt = " + Temp5 +"   imp = "+ Temp6;
      k++;
    }
  }
  k++;
  m_stringArray[k]= "OPTIONS";
	k++;
  _itoa(m_nPercentForTest,string,10);
  Temp = string;
  m_stringArray[k]= "Testing on " + Temp + " of the data set. To test on a new data file, restart using a .fit file ";
  k++;
  if(m_nEstimateRMSError == 0)
  {
    m_stringArray[k]= "We are estimating RMS Error in the data file";
  }
  k++;
  if(m_nNoJitter == 0)
  {
    m_stringArray[k]= "Inputs not jittered";
  }
  else
  {
    m_stringArray[k]= "Inputs jittered using triangular distribution";
  }
	k++;
	if (m_nZeroSmoothing == 0)
	{
		m_stringArray[k] = "Zero smoothing applied";
	}
	else
	{
		m_stringArray[k] = "Constant smoothing applied during approximation";
	}
  k++;
  if(m_nDTREEDepth == 1)
  {
		m_stringArray[k] = "Generating one-layer DTREE; i.e. without partitioning the input space";
  }
  else
  {
    _itoa(m_nDTREEDepth,string,10);
    Temp = string;
    m_stringArray[k]= "Generating DTREE with " + Temp + " layer(s)";
  }
  k++;
  if(m_bReplaceUndefined)
  {
     m_stringArray[k]= "Replacing undefined outputs in Data file with computed values";
	}
	else
	{
     m_stringArray[k]= "Adding a new column to right of Data file with computed values";
  }
  k++;
  if(m_bTimePrefixes)
  {
     m_stringArray[k]= "Output file names will have time prefix eg for 2:23PM prefix is 1423";
     k++;
  }
  if(m_bDiagnostics)
  {
     m_stringArray[k]= "Some internal files will be printed out to diagnose errors and bugs";
     k++;
  }
}

/////////////////////////////////////////////////////////////////////////////
// CALNfitDeepDoc diagnostics

#ifdef _DEBUG
void CALNfitDeepDoc::AssertValid() const
{
	CDocument::AssertValid();
}

void CALNfitDeepDoc::Dump(CDumpContext& dc) const
{
	CDocument::Dump(dc);
}
#endif //_DEBUG

/////////////////////////////////////////////////////////////////////////////
// CALNfitDeepDoc commands

void CALNfitDeepDoc::EditClearall() 
{ 
 	m_strDataFileName = ".txt";
	m_strDTREEFileName = ".dtr";
	m_nALNs = 1;
	m_nFit = 0;
	m_nTrain = 0;
	m_nMessage = 8;
	CString m_strProtocolFileName = "";   // set up by start button
	CString m_strScatterFileName = "";     // set up by start button
	GenerateReportString();
	m_nPercentProgress = 0;
  m_nColsUniv = nColsUniv =101;
  m_nALNinputs = nALNinputs =101;
  m_nDisplayedALNinput = 100;
  for(int nn = 0; nn < m_nColsUniv; nn++) // for the training default m_nColsUniv equals m_nALNinputs
  {
    char buff[32] = "__A";
    //_itoa(nn + 1,buff+7,10);
    if(nn/26>0)buff[1] = 'A'- 1 + nn/26;
    buff[2] = 'A' + nn%26;
    buff[3] = '\0';
    m_strVarname[nn] = buff;
    for(int charno = 0; charno < 32; charno++)
    {
      varname[nn][charno]=buff[charno];
    }
    m_nInputCol[nn] = nInputCol[nn] = nn;
    m_nLag[nn] = nLag[nn] = 0;
    m_nIOpropMessageNumber[nn] = INPUT;
    if(nn < m_nColsUniv -1)
    {
      dblMinWeight[nn] = m_dblMinWeight[nn] = -dblMax;
      dblMaxWeight[nn] = m_dblMaxWeight[nn] = dblMax;
    }
    else
    {
      dblMinWeight[nn] = m_dblMinWeight[nn] = -1;
      dblMaxWeight[nn] = m_dblMaxWeight[nn] = -1;
    }
    m_dblImportance[nn] = dblImportance[nn] = 0;
  }
  if(m_nTrain == 0) // if this is a training run, we show the output column as such
  {
    m_nIOpropMessageNumber[m_nALNinputs-1] = OUTPUT;
  }
  m_bStarted = FALSE;
  // Set all option values
  m_nPercentForTest = 10;
  m_nEstimateRMSError = 0;
	m_nZeroSmoothing = 0;
  m_nNoJitter = 0;
  m_nDTREEDepth = 1;
	m_bTimePrefixes = TRUE;
  m_bReplaceUndefined = FALSE;
  m_bDiagnostics = FALSE;
  UpdateAllViews(NULL);
}

void CALNfitDeepDoc::DocUpdateScreen()
{
  // takes global values for message number and percent progress,
  // updates the doc and the screen
  m_nMessage = nMessageNumber;
  GenerateReportString();
  m_nPercentProgress = nPercentProgress;
  UpdateAllViews(NULL);
}

void CALNfitDeepDoc::DocThreadFinished()
{
  // this time we don't pass back the message 
  // from the array; we get it from the result message message string
  m_nMessage = nMessageNumber;
  // Model m_strReport.Format("LFs: %d/%d  Train RMSE: %.3f",nActiveLFNs,nLFNs,dblRMSE);
  if(nPercentForTest == 0)
  {
    m_strReport.Format("Evaluation was not done, since no test set was provided.");
  }
  else
  {
    if(bClassify)
    {
      m_strReport.Format("Number of classification errors: %d, which is %.3f percent of the test set ",nEvalMisclassifications,dblEvalMisclassificationPercent);
    }
    else
    {
      m_strReport.Format("RMS difference from desired value (or from 0 if not given in test set) is %.3f",dblEvalRMSError);
    }
  }
  m_nPercentProgress = nPercentProgress;
  UpdateAllViews(NULL);
}

void CALNfitDeepDoc::OnFileSave() 
{
	// TODO: Add your command handler code her

  CDocument::OnFileSave();
}

void CALNfitDeepDoc::OnFileSaveAs() 
{
	// TODO: Add your command handler code here
  m_strPathName = m_strTentativePathName;
  CDocument::OnFileSaveAs();
}


typedef int (* CFT) (const void *, const void *);

void CALNfitDeepDoc::OnSort()
{
  bWarnEqual = FALSE;
  // load the sort array with the current ALN inputs
  for(int i = 0; i < m_nALNinputs - 1; i++) // don't sort the output
  {
    SortArray[i].nIndex = i;
    SortArray[i].nInputCol = m_nInputCol[i] ;
    SortArray[i].nLag = m_nLag[i];
    SortArray[i].nIOpropMessageNumber = m_nIOpropMessageNumber[i];
    SortArray[i].dblMinWeight = m_dblMinWeight[i];
    SortArray[i].dblMaxWeight = m_dblMaxWeight[i];
    SortArray[i].dblImportance = m_dblImportance[i];
  }
  if(m_bSortCols)
  {
    qsort((void **)SortArray, m_nALNinputs - 1, sizeof(SRTND), (CFT) SNCompareCols);
    // check to see if an input is the same as the output and with lag 0 too
    for(int i = 0; i < m_nALNinputs -1; i++)
    {
      if((m_nInputCol[i] == m_nInputCol[m_nALNinputs -1]) && (m_nLag[i] == 0))
      {
        bWarnEqual = TRUE;
      }
    }
  }
  else
  {
    qsort((void **) SortArray,m_nALNinputs - 1, sizeof(SRTND), (CFT) SNCompareImportance);
  }
  for(int i = 0; i < m_nALNinputs - 1; i++) // replace the sorted values
  {
    m_nInputCol[i] = nInputCol[i] = SortArray[i].nInputCol;
    m_nLag[i] = nLag[i] = SortArray[i].nLag;
    m_nIOpropMessageNumber[i] = SortArray[i].nIOpropMessageNumber;
    m_dblMinWeight[i] = SortArray[i].dblMinWeight;
    m_dblMaxWeight[i] = SortArray[i].dblMaxWeight;
    m_dblImportance[i] = dblImportance[i] = SortArray[i].dblImportance;
  }
  m_nDisplayedALNinput = 0;
  m_bWarnEqual = bWarnEqual;
  UpdateAllViews(NULL);
}

int SNCompareCols(void * p1, void * p2)
{
  if(((SRTND *) p1)->nInputCol < ((SRTND *) p2)->nInputCol)
  {
    return -1;
  }
  else if(((SRTND *) p1)->nInputCol > ((SRTND *) p2)->nInputCol)
  {
    return 1;
  }
  else // the input columns are equal, sort by lags
  {
    if(((SRTND *) p1)->nLag < ((SRTND *) p2)->nLag)
    {
      return -1;
    }
    else if(((SRTND *) p1)->nLag > ((SRTND *) p2)->nLag)
    {
      return 1;
    }
    else // the lags are also equal
    {
      bWarnEqual = TRUE;
      return 0;
    }
  }
}

int SNCompareImportance(void * p1, void * p2)
{
  if(((SRTND *) p1)->dblImportance < ((SRTND *) p2)->dblImportance)
  {
    return -1;
  }
  else if(((SRTND *) p1)->dblImportance > ((SRTND *) p2)->dblImportance)
  {
    return 1;
  }
  else // the importances are equal
  {
    return 0;
  }
}

