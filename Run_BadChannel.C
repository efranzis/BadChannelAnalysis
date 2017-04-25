/////////////////////////////////////////////////
//
// This macro has been developed to find bad cell candidates in EMCal and DCal based on cell amplitude distributions
// Input needed can be either outputs QA from AliAnaCalorimeterQA task (BadChannelAnalysis() function)
// Or from merged output of AliAnalysisTaskCaloCellsQA (use BCAnalysis() function)
//
// See https://twiki.cern.ch/twiki/bin/view/ALICE/BadChannelAnalysis
// for documentation
// ---------------------
//  Running the macro
// ---------------------
// use root -b to speed up (no canvas drawn)
// root [2] .L Run_BadChannel.C++
// root [2] Run_BadChannel(1,"LHC15o","Train_758","AnyINTnoBC",244411,"","runlist")
// root [2] SummarizeRunByRun("LHC15o","Train_758","AnyINTnoBC",".","runList45.txt")
//
//.x Run_BadChannel.C+(0,"LHC15o","Train_758","AnyINTnoBC",245554,"245554_Filtered.root","",".")
/////////////////////////////////////////////////

// --- ROOT system ---
#include <TString.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TH2D.h>
#include <TH1D.h>
#include <TList.h>
#include <TSystem.h>
#include <TStopwatch.h>
#include <TError.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TEnv.h>
#include <TGaxis.h>
#include <TLegend.h>

// --- ANALYSIS system ---
#include "AliAnaCaloChannelAnalysis.h" //include when compile
#include "AliEMCALGeometry.h"          //include when compile
#include "AliCalorimeterUtils.h"       //include when compile
#include "AliAODEvent.h"               //include when compile
#include "AliOADBContainer.h"          //include when compile

//colors
const Int_t colors[]       = {kBlack, kRed+1 , kBlue+1, kGreen+3, kMagenta+1, kOrange-1,kCyan+2,kYellow+2, kRed-9, kGreen-10, kBlue - 8};
const Int_t RainbowColors[]= {kRed, kRed-4, kRed-7, kRed-9, kRed-10, kYellow, kYellow-4, kYellow-7, kYellow-9, kYellow-10, kGreen, kGreen-4 , kGreen-7, kGreen-9, kGreen-10, kCyan, kCyan-4, kCyan-7, kCyan-9, kCyan-10, kBlue, kBlue-4, kBlue-7, kBlue-9, kBlue-10, kMagenta, kMagenta-4, kMagenta-7, kMagenta-9, kMagenta-10};

//definition of methods
TH1D* GetAbsIDHistogram(Int_t celln, TFile *fin);
TH1D** GetAllAbsIDHistogram(TFile *fin, Int_t& totbadCells);
TList* GetListAbsIDHistograms(TFile *fin);
TString RunSpecificHistogramFileName(Int_t runnumb, TString period, TString train = "Train_641", TString trigger= "AnyINTnoBC", TString workDir=".");
TH2F* CompressHistogram(TH2 *Histo,Int_t totalCells, Int_t badCells,Int_t nRuns);
void BuildMaxMinHisto(TH1D* inHisto, TH1D* minHist,TH1D* maxHist);
void PlotLowFractionCells(TString pdfName, std::vector<Int_t> cellVector,TH2F* badVsCell[],Int_t nRuns,TH2F* ampID[],TH1D* hCellGoodMean[]);
Bool_t IsItReallyBadRatio(TH1D* minHistoRatio,TH1D* maxHistoRatio,TH1D* meanHistoRatior,TString& crit);
Bool_t isHIRun(Int_t runID);
void SetHisto(TH2 *Histo,TString Xtitel,TString Ytitel,Bool_t longhisto);
void SetHisto(TH1 *Histo,TString Xtitel,TString Ytitel,Bool_t longhisto);
//void GetBestPeriodSplitting(TString period = "LHC15o", TString train = "Train_641",Int_t noOfSplits=4);

//________________________________________________________________________
void Run_BadChannel(Int_t nversion = 1, TString period = "LHC15n", TString train = "Train_603", TString trigger= "AnyINT", Int_t runNum= 245683, TString externalFile= "",TString listName="runList.txt",TString workDir=".")
{
	TStopwatch watch;
	watch.Start();
	//gErrorIgnoreLevel= 5000; //=1  =kWarning  //kInfo
	gROOT->ProcessLine("gErrorIgnoreLevel = kWarning;"); //Does not work -

	AliAnaCaloChannelAnalysis* Analysis;

	nversion=runNum; //..If you do the analysis run by run - this might be helpful
	Analysis=new AliAnaCaloChannelAnalysis(period,train,trigger,runNum,nversion,workDir,listName);

	//..Settings
	Analysis->SetExternalMergedFile(externalFile);
	//Analysis->SetQAChecks(1);  //1=Perform QA checks - takes a long time! Prints all good cells for cross check
	//Analysis->SetPrintOutput(1);  //1= prints more information about excluded cells

	//. . . . . . . . . . . . . . . . . . . . . . . .
	//. . Add different period analyses
	//. . . . . . . . . . . . . . . . . . . . . . . .
	//..the range of sigmas should be selected such
	//..that one does not cut into the natural fluctuation over the modules
/*This is version 1
    Analysis->AddPeriodAnalysis(2, 5.5,0.1,0.3);  // hits in cell in range Emin Emax
	Analysis->AddPeriodAnalysis(1, 4.0,0.1,0.3);  // energy/hit in range Emin Emax
	Analysis->AddPeriodAnalysis(2, 5.5,0.2,0.5);  // hits in cell range Emin Emax
	Analysis->AddPeriodAnalysis(1, 4.5,0.2,0.5);  // energy/hit in range Emin Emax
	Analysis->AddPeriodAnalysis(2, 4.0,0.5,1.0);  // hits in cell range Emin Emax
	Analysis->AddPeriodAnalysis(1, 4.5,0.5,1.0);  // energy/hit in range Emin Emax
	Analysis->AddPeriodAnalysis(2, 4.5,1.0,4.0);  // hits in cell range Emin Emax
	Analysis->AddPeriodAnalysis(1, 4.5,1.0,4.0);  // mean energy in range Emin Emax
	Analysis->AddPeriodAnalysis(2, 4.5,1.0,10.0); // hits in cell in range Emin Emax
	Analysis->AddPeriodAnalysis(1, 4.5,1.0,10.0); // energy/hit in range Emin Emax

	//..special test for extra high energy fluctuations
	Analysis->AddPeriodAnalysis(2, 5.5,3.0,40.0); // hits in cell in range Emin Emax
*/
/*This is version 3
 * 	Analysis->AddPeriodAnalysis(2, 5.5,0.1,0.3);  // hits in cell in range Emin Emax
	Analysis->AddPeriodAnalysis(1, 4.0,0.1,0.3);  // energy/hit in range Emin Emax
	Analysis->AddPeriodAnalysis(2, 5.5,0.2,0.5);  // hits in cell range Emin Emax
	Analysis->AddPeriodAnalysis(1, 4.5,0.2,0.5);  // energy/hit in range Emin Emax
	Analysis->AddPeriodAnalysis(2, 5.0,0.5,1.0);  // hits in cell range Emin Emax
	Analysis->AddPeriodAnalysis(1, 4.5,0.5,1.0);  // energy/hit in range Emin Emax
	Analysis->AddPeriodAnalysis(2, 5.0,1.0,4.0);  // hits in cell range Emin Emax
	Analysis->AddPeriodAnalysis(1, 4.5,1.0,4.0);  // mean energy in range Emin Emax
	Analysis->AddPeriodAnalysis(2, 5.0,1.0,10.0); // hits in cell in range Emin Emax
	Analysis->AddPeriodAnalysis(1, 4.5,1.0,10.0); // energy/hit in range Emin Emax

	//..special test for extra high energy fluctuations
	Analysis->AddPeriodAnalysis(2, 5.0,3.0,40.0); // hits in cell in range Emin Emax
*/
	Analysis->AddPeriodAnalysis(2, 5.5,0.1,0.3);  // hits in cell in range Emin Emax
	Analysis->AddPeriodAnalysis(1, 4.0,0.1,0.3);  // energy/hit in range Emin Emax
	Analysis->AddPeriodAnalysis(2, 5.5,0.2,0.5);  // hits in cell range Emin Emax
	Analysis->AddPeriodAnalysis(1, 4.5,0.2,0.5);  // energy/hit in range Emin Emax
	//Analysis->AddPeriodAnalysis(2, 5.5,0.3,0.6);  //neu* hits in cell range Emin Emax
	//Analysis->AddPeriodAnalysis(1, 4.5,0.3,0.6);  //neu* energy/hit in range Emin Emax
	Analysis->AddPeriodAnalysis(2, 5.5,0.5,1.0);  // hits in cell range Emin Emax
	Analysis->AddPeriodAnalysis(1, 4.5,0.5,1.0);  // energy/hit in range Emin Emax
	//Analysis->AddPeriodAnalysis(2, 5.5,1.0,2.0);  //neu* hits in cell range Emin Emax
	//Analysis->AddPeriodAnalysis(1, 4.5,1.0,2.0);  //neu* mean energy in range Emin Emax
	Analysis->AddPeriodAnalysis(2, 5.5,1.0,4.0);  // hits in cell range Emin Emax
	Analysis->AddPeriodAnalysis(1, 4.5,1.0,4.0);  // mean energy in range Emin Emax
	Analysis->AddPeriodAnalysis(2, 5.5,1.0,10.0); // hits in cell in range Emin Emax
	Analysis->AddPeriodAnalysis(1, 4.5,1.0,10.0); // energy/hit in range Emin Emax

	//..special test for extra high energy fluctuations
	Analysis->AddPeriodAnalysis(2, 4.5,3.0,40.0); // hits in cell in range Emin Emax


	///*test time stuff*/	Analysis->AddPeriodAnalysis(3, 6,-20,+20);// energy/hit in range Emin Emax

	//..Start the bad channel analysis
	Bool_t mergeOnly=0;//.. =1 do only merge and filter
	Analysis->Run(mergeOnly);

	watch.Stop();
	watch.Print();
}
//
// little helper function
// check where the bad cell is (row collumn), if you have only its ID
//________________________________________________________________________
void Get_RowCollumnID(Int_t runNum= 244411,Int_t inputCellID=-1,Int_t inputRow=-1,Int_t inputCollumn=-1,Int_t inputSM=-1)
{
	//......................................................
	//..Initialize EMCal/DCal geometry
	AliCalorimeterUtils* fCaloUtils = new AliCalorimeterUtils();
	//..Create a dummy event for the CaloUtils
	AliAODEvent* aod = new AliAODEvent();
	fCaloUtils->SetRunNumber(runNum);
	fCaloUtils->AccessGeometry(aod);
	AliEMCALGeometry * geom = fCaloUtils->GetEMCALGeometry();

	//..get row collumn from cell ID
	Int_t cellColumn=0,cellRow=0;
	Int_t cellColumnAbs=0,cellRowAbs=0;
	Int_t cellID=0;
	Int_t trash;

	cout<<"...................................................."<<endl;
	cout<<""<<endl;
	if(inputCellID!=-1)
	{
		fCaloUtils->GetModuleNumberCellIndexesAbsCaloMap(inputCellID,0,cellColumn,cellRow,trash,cellColumnAbs,cellRowAbs);
		cout<<"Cell Id provided: "<<inputCellID<<endl;
		cout<<"This corresponds to absolute row: "<<cellRowAbs<<" and absolute collumn: "<<cellColumnAbs<<endl;
	}
	if(inputRow!=-1 && inputCollumn!=-1)
	{
		cout<<"Supermodule provided: "<<inputSM<<endl;
		cout<<"Absolute row provided: "<< inputRow<<" and absolute collumn provided : "<<inputCollumn<<endl;
		cellID=geom->GetAbsCellIdFromCellIndexes(inputSM,inputRow,inputCollumn);
		cout<<"This corresponds to Cell Id: "<<cellID<<endl;
	}
	cout<<""<<endl;
	cout<<"...................................................."<<endl;
}
//
// little helper function
// Test if the file committed to OADB is correct
//________________________________________________________________________
void Test_OADB(TString period="LHC15n",Int_t trainNo=603,Int_t version=5,Int_t runNumber=244411)
{
	//......................................................
	// Test if the file committed to OADB is correct
	//......................................................
	Int_t nSM = 20;
	TH1C *h[20];
	Int_t cellColumnAbs,cellRowAbs,trash,cellID;

	//......................................................
	//..Get the OADB information
	TString fBasePath="/Users/Eliane/Software/alice/sw/osx_x86-64/AliPhysics/0-1/OADB/EMCAL";
	AliOADBContainer *cont=new AliOADBContainer("");
	cont->InitFromFile(Form("%s/EMCALBadChannels.root",fBasePath.Data()),"AliEMCALBadChannels");

	//......................................................
	//..Get the .root file with the original histogram to compare if they coincide
//	TString path        = Form("/Users/Eliane/Software/BadChannelAnalysis/AnalysisOutput/%s/Version%d",period.Data(),version);
	TString path        = Form("/Users/Eliane/Software/BadChannelAnalysis/AnalysisOutput/%s/Train_%i/VersionINT7Glob",period.Data(),trainNo);
//	TString rootFileName= Form("Train_%dAnyINTnoBC_Histograms_V%d.root",trainNo,version);
	TString rootFileName= Form("INT7_Histograms_V265630.root");
	TFile* outputRoot   = TFile::Open(Form("%s/%s",path.Data(),rootFileName.Data()));

	if(!outputRoot)cout<<"File "<<outputRoot->GetName()<<" does not exist"<<endl;

	TH2F* h2DChannelMap_FlagBad =(TH2F*)outputRoot->Get("2DChannelMap_Flag2");
	TH2F* h2DChannelMap_FlagDead=(TH2F*)outputRoot->Get("2DChannelMap_Flag1");

	/*cout<<"Open root file: "<<outputRoot->GetName()<<endl;
	if(h2DChannelMap_FlagBad)h2DChannelMap_FlagBad->SetName("BadMapFromFile");
	else cout<<"No valid file found"<<endl;*/
	//......................................................
	//..Initialize EMCal/DCal geometry
	AliCalorimeterUtils* fCaloUtils = new AliCalorimeterUtils();
	//..Create a dummy event for the CaloUtils
	AliAODEvent* aod = new AliAODEvent();
	fCaloUtils->SetRunNumber(runNumber);
	fCaloUtils->AccessGeometry(aod);
	AliEMCALGeometry * geom = fCaloUtils->GetEMCALGeometry();

	//.......................................................
	//..build two dimensional histogram with values row vs. column
	//..with info from OADB
	Int_t fNMaxCols    = 48;  //eta direction
	Int_t fNMaxRows    = 24;  //phi direction

	Int_t fNMaxColsAbs = 2*fNMaxCols;
	Int_t fNMaxRowsAbs = Int_t (nSM/2)*fNMaxRows; //multiply by number of supermodules
	TString histoName;
	histoName = Form("2DChannelMap_Flag_Bad");
	TH2F *plot2D_Bad_OADB = new TH2F(histoName,histoName,fNMaxColsAbs+1,-0.5,fNMaxColsAbs+0.5, fNMaxRowsAbs+1,-0.5,fNMaxRowsAbs+0.5);
	///*OADB  looks like Marcels figures*/	TH2F *plot2D_Bad_OADB = new TH2F(histoName,histoName,fNMaxColsAbs,0,fNMaxColsAbs, fNMaxRowsAbs,0,fNMaxRowsAbs);
	plot2D_Bad_OADB->GetXaxis()->SetTitle("cell column (#eta direction)");
	plot2D_Bad_OADB->GetYaxis()->SetTitle("cell row (#phi direction)");
	histoName = Form("2DChannelMap_Flag_Dead");
	TH2F *plot2D_Dead_OADB = new TH2F(histoName,histoName,fNMaxColsAbs+1,-0.5,fNMaxColsAbs+0.5, fNMaxRowsAbs+1,-0.5,fNMaxRowsAbs+0.5);
	///*OADB looks like Marcels figures*/	TH2F *plot2D_Dead_OADB = new TH2F(histoName,histoName,fNMaxColsAbs,0,fNMaxColsAbs, fNMaxRowsAbs,0,fNMaxRowsAbs);
	plot2D_Dead_OADB->GetXaxis()->SetTitle("cell column (#eta direction)");
	plot2D_Dead_OADB->GetYaxis()->SetTitle("cell row (#phi direction)");

	//.......................................................
	//.. Read the Bad Channel map from OADB
	TObjArray *recal=(TObjArray*)->GetObject(runNumber);
	//recal->ls();
	TCanvas* C1 = new TCanvas();
	C1->cd();
	//nSM=1;
	for(Int_t iSM = 0; iSM < nSM; iSM++)
	{
		h[iSM]=(TH1C *)recal->FindObject(Form("EMCALBadChannelMap_Mod%d",iSM));
		h[iSM]->GetXaxis()->SetTitle("Column"); // tmp debug
		h[iSM]->GetYaxis()->SetTitle("Row"); // tmp debug
		h[iSM]->DrawCopy("colz"); // tmp debug

		//..Loop though the SM to set which cells are bad
		for(Int_t column=0;column<48;column++)
		{
			for(Int_t row=0;row<24;row++)
			{
				Int_t inRow=row;
				Int_t inCol=column;
				cellID=geom->GetAbsCellIdFromCellIndexes(iSM,inRow,inCol);
				fCaloUtils->GetModuleNumberCellIndexesAbsCaloMap(cellID,0,inCol,inRow,trash,cellColumnAbs,cellRowAbs);
				if(h[iSM]->GetBinContent(column,row)==2)
					//if(h[iSM]->GetBinContent(column,row)>1)//..bad and warm
				{
					plot2D_Bad_OADB->SetBinContent(cellColumnAbs,cellRowAbs,1);
				}
				if(h[iSM]->GetBinContent(column,row)==1)
				{
					plot2D_Dead_OADB->SetBinContent(cellColumnAbs,cellRowAbs,1);
				}
			}
		}
	}

	//..................................................................
	TCanvas* C2 = new TCanvas("Info from OADB","Info from OADB",1);
	C2->Divide(2);
	C2->cd(1);
	plot2D_Bad_OADB->DrawCopy("colz");
	C2->cd(2);
	plot2D_Dead_OADB->DrawCopy("colz");
	//..................................................................
	TCanvas* C3 = new TCanvas("Orig Maps from Analyis","Orig Maps from Analyis",1);
	C3->Divide(2);
	C3->cd(1);
	h2DChannelMap_FlagBad->DrawCopy("colz");
	C3->cd(2);
	h2DChannelMap_FlagDead->DrawCopy("colz");
	//..................................................................
	TCanvas* C4 = new TCanvas("Subtraction of OADB-Orig.","Subtraction of OADB-Orig.",1);
	C4->Divide(2);
	C4->cd(1);
	plot2D_Bad_OADB->Add(h2DChannelMap_FlagBad,-1);
	plot2D_Bad_OADB->DrawCopy("colz");
	plot2D_Bad_OADB->Add(h2DChannelMap_FlagBad,+1);
	C4->cd(2);
	plot2D_Dead_OADB->Add(h2DChannelMap_FlagDead,-1);
	plot2D_Dead_OADB->DrawCopy("colz");
	plot2D_Dead_OADB->Add(h2DChannelMap_FlagDead,+1);

	//..................................................................
	TCanvas* C5 = new TCanvas("Division of OADB/Orig.","Division of OADB/Orig.",1);
	C5->Divide(2);
	C5->cd(1);
	plot2D_Bad_OADB->Divide(h2DChannelMap_FlagBad);
	plot2D_Bad_OADB->DrawCopy("colz");
	C5->cd(2);
	plot2D_Dead_OADB->Divide(h2DChannelMap_FlagDead);
	plot2D_Dead_OADB->DrawCopy("colz");
}
// TEST TESTE TEST ELI ELI ELI
/// Draw the good, bad, dead channel maps, and the amplitude distribution per each run and save them in pdf files in analysisOutput/train/RunByRunSummary
/// -- can improve: write different pages in one pdf
/// -- decide how to treat this info
//________________________________________________________________________
void SummarizeRunByRun(TString period = "LHC15o", TString train = "Train_641", TString trigger= "AnyINTnoBC", TString workDir=".", TString listName="runList")
{
	gStyle->SetOptTitle(0);
	gStyle->SetOptStat(0);
	//gStyle->SetPalette(53);  //standard is 1
	gStyle->SetCanvasColor(10);
	TGaxis::SetMaxDigits(4);
	gStyle->SetPadTopMargin(0.07);//0.05
	gStyle->SetPadBottomMargin(0.18);//0.15
	gStyle->SetPadRightMargin(0.10);
	gStyle->SetPadLeftMargin(0.15);
	gStyle->SetFrameFillColor(10);
	gStyle->SetLabelSize(0.05,"X");
	gStyle->SetLabelSize(0.05,"Y");
	gStyle->SetTitleSize(5.0,"X");
	gStyle->SetTitleSize(5.0,"Y");
	gEnv->SetValue("Canvas.ShowEventStatus",1);  //shows the status bar in the canvas

	//..............................................
	//..manually disable cells
	std::vector<Int_t> badcellsBlock1;
	std::vector<Int_t> badcellsBlock2;
	std::vector<Int_t> badcellsBlock3;
	badcellsBlock1.push_back(13483);
//	badcellsBlock2.push_back(12360);
	badcellsBlock2.push_back(13483);
	badcellsBlock3.push_back(13483);

	//..select runs after which a new bad map is built
	Int_t splitRuns1=44;
	Int_t splitRuns2=71;
	//..............................................

	TString analysisInput  = Form("AnalysisInput/%s",period.Data());
//	TString analysisOutput = Form("AnalysisOutput/%s/%s/CutRangesV1",period.Data(),train.Data());  // CutRangesV3
	TString analysisOutput = Form("AnalysisOutput/%s/%s",period.Data(),train.Data());  // CutRangesV3
	TString runList        = Form("%s/%s/%s/%s.txt",workDir.Data(), analysisInput.Data(), train.Data(),listName.Data());
	gSystem->mkdir(TString::Format("%s/RunByRunSummary/", analysisOutput.Data()));

	//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	//..open the text file and save the run IDs into the RunId[] array
	cout<<"o o o Open .txt file with run indices. Name = " << runList << endl;
	FILE *pFile = fopen(runList.Data(), "r");
	if(!pFile)
	{
		cout<<"couldn't open file "<<runList<<"!"<<endl;
		return;
	}
	Int_t q;
	Int_t ncols;
	Int_t nlines = 0 ;
	Int_t RunId[500] ;
	while (1)
	{
		ncols = fscanf(pFile,"  %d ",&q);
		if (ncols< 0) break;
		RunId[nlines]=q;
		nlines++;
	}
	fclose(pFile);

	//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	//..Create different canvases for run-by-runcomparision
	cout<<"o o o Found " << nlines <<" files in list"<< endl;
	const Int_t nRun = nlines;
	Int_t nRunsUsed = nlines;
	//ELI for MartinInt_t totalperCv = 4;
	Int_t totalperCv = 16;
	Int_t nPad = TMath::Sqrt(totalperCv);
	Int_t nCv = nRun/totalperCv+1;


	if(nCv<1)nCv=1;

	//..canvases per run
	TCanvas **cBad  = new TCanvas*[nCv];
	TCanvas **cGood = new TCanvas*[nCv];
	TCanvas **cDead = new TCanvas*[nCv];
	TCanvas **cAmp  = new TCanvas*[nCv];

	for(Int_t ic = 0; ic<nCv; ic++)
	{
		cBad [ic] = new TCanvas(TString::Format("badcells%d", ic), TString::Format("I) badcells  (%d/%d)", ic+1, nCv), 1000,750);
		cGood[ic] = new TCanvas(TString::Format("goodcells%d", ic),TString::Format("I) goodcells (%d/%d)", ic+1, nCv),1000,750);
		cDead[ic] = new TCanvas(TString::Format("deadcells%d", ic),TString::Format("I) deadcells (%d/%d)", ic+1, nCv),1000,750);
		cAmp [ic] = new TCanvas(TString::Format("Amplitide%d", ic),TString::Format("I) Amplitide (%d/%d)", ic+1, nCv),1000,750);

		cBad [ic] ->Divide(nPad,nPad,0.001,0.001);
		cGood[ic] ->Divide(nPad,nPad,0.001,0.001);
		cDead[ic] ->Divide(nPad,nPad,0.001,0.001);
		cAmp [ic] ->Divide(nPad,nPad,0.001,0.001);
	}

	//..summary figures for all runs
	Int_t nFlags = 3;
	TH2F** hFlagvsRun = new TH2F*[nFlags]; 
	TH2F* hFlagNew  = 0x0;
	TH2F* hFlagNewClean  = 0x0;
	TH2F** ampID         = new TH2F*[nRun];
	TH2F** ampIDCl       = new TH2F*[nRun];
	TH2F** ampIDCl3Block = new TH2F*[nRun];
	TH2F** ampIDCl1Block = new TH2F*[nRun];
	TH1D** hCellGoodMean = new TH1D*[nRun];
	TH1D** hNEvent       = new TH1D*[nRun];
	TH2F* hBadVsEvent    = new TH2F("hBadVsEvent","hBadVsEvent",100,100000,25000000,60,700,1800);
	TH2F* hDeadBadVsEvent= new TH2F("hDeadBadVsEvent","hDeadBadVsEvent",100,100000,25000000,60,700,1800);
	TH2F* hBadVsEventHI    = new TH2F("hBadVsEventHI","hBadVsEventHI",100,100000,25000000,60,700,1800);
	TH2F* hDeadBadVsEventHI= new TH2F("hDeadBadVsEventHI","hDeadBadVsEventHI",100,100000,25000000,60,700,1800);
	TH1D* deadbadCellsVsRun;
	TH1D* deadCellsVsRun;
	TH1D* badCellsVsRun;
	TH1D* deadCellsVsRunC;
	TH1D* badCellsVsRunC;
	TH1D* projSum;
	TH1D* projSumC;
	TH1D* projSumC3Blocks;
	TH1D* projSumC3BlocksA;
	TH1D* projSumC3BlocksB;
	TH1D* projSumC3BlocksC;
	TH1D* projSumC1Block;
	TH1D* nEventsVsRuns;
	TH2F* Sum2DSingleMask;
	TH2F* Sum2D3BlockMask;
	TH2F* Sum2DOrig;
	TH2F* Sum2DIdeal;
	TH1D* hgoodMean;

	for(Int_t i = 0; i<nFlags; i++)
	{
		hFlagvsRun[i] = 0x0;
	}
	TCanvas *cFlagDeadBad     = new TCanvas("cFlagDeadBad", "II) Flag dead or bad", 1600, 500);
	TCanvas *cFlagSum         = new TCanvas("cFlagSum", "II) Flag dead&bad", 1600, 500);
	TCanvas *cFlagSumCleaned  = new TCanvas("cFlagSumClean", "III) cleanded Flag dead&bad", 1600, 500);
	TCanvas *cFlagSumCompAll  = new TCanvas("cFlagSumCompAll", "III) compressed Flag dead&bad", 1600, 500);
	TCanvas *cFlagSumCompClean= new TCanvas("cFlagSumComp", "III) compressed&cleaned Flag dead&bad", 1600, 500);

	TCanvas *cFlagNew         = new TCanvas("cFlagNew", "IV) Frac dead&bad 2D", 1600, 800);
	TCanvas *cellSummaryCan   = new TCanvas("cellSummaryCan", "I) run overview", 1600, 800);
	TCanvas *cellSummaryCan2  = new TCanvas("cellSummaryCan2","I) run overview II",1600,800);
	TCanvas *cAmpSum          = new TCanvas("SumOfAmplitudes","I) Sum of Amplitides over runs",1500,750);
	TCanvas *cAmpSum2D        = new TCanvas("SumOf2DAmplitudes","I) Sum of 2D Amplitides over runs",1500,750);

	//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	//..defining variables and histograms
	TString rootFileName;
    TString badChannelOutput;
	/*TString badChannelOutput[4];
	badChannelOutput[0]="AnalysisOutput/LHC16h/Version2/Train_622AnyINTnoBC_Histograms_V2.root";
	badChannelOutput[1]="AnalysisOutput/LHC16i/Version2/Train_623AnyINTnoBC_Histograms_V2.root";
	badChannelOutput[2]="AnalysisOutput/LHC16k/Version0/Train_658AnyINTnoBC_Histograms_V0.root";
	badChannelOutput[3]="AnalysisOutput/LHC16o/Version3/Train_663AnyINTnoBC_Histograms_V3.root";
    */
	Int_t   noOfCells=0;
	Int_t usedRuns=0;
	Bool_t iFirstRun=0;

	AliCalorimeterUtils *fCaloUtils = new AliCalorimeterUtils();
	//..Create a dummy event for the CaloUtils
	AliAODEvent* aod = new AliAODEvent();
	fCaloUtils->SetRunNumber(RunId[0]);
	fCaloUtils->AccessGeometry(aod);
    noOfCells=fCaloUtils->GetEMCALGeometry()->GetNCells(); //..Very important number, never change after that point!

	hFlagvsRun[0] = new TH2F("hFlag1vsRun", "hFlag1vsRun (?); cell ID; Run number", noOfCells+10, 0, noOfCells+10, nRun, 0, nRun); // update this axis, need to have the run number
	hFlagvsRun[1] = new TH2F("hFlag2vsRun", "hFlag2vsRun (?); cell ID; Run number", noOfCells+10, 0, noOfCells+10, nRun, 0, nRun);
	hFlagvsRun[2] = new TH2F("hFlag3vsRun", "hFlag3vsRun (?); cell ID; Run number", noOfCells+10, 0, noOfCells+10, nRun, 0, nRun);
	nEventsVsRuns = new TH1D("nEventVsRun", "number of events in run", nRun, 0, nRun);

	//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	//..Open the different .root files with help of the run numbers from the text file
	//..loop over the amount of run numbers found in the previous text file.
	for(Int_t i = 0 ; i < nRun ; i++)  //Version%i" LHC16i_muon_caloLego_Histograms_V255539
	{
		rootFileName      = Form("%s_Histograms_V%i.root", trigger.Data(), RunId[i]);
		badChannelOutput  = Form("%s/Version%i/%s", analysisOutput.Data(), RunId[i],rootFileName.Data());

		//Martin cout<<"Open root file: "<<badChannelOutput[i]<<endl;
		cout<<"Open root file No: "<<i+1<<" - "<<badChannelOutput<<" - "<<flush;
		TFile *f = TFile::Open(badChannelOutput);
		if(!f)
		{
			cout<<"Couldn't open/find .root file: "<<badChannelOutput<<endl;
			cout<<endl;
			break;
		}

		//..you may introduce a cut here, selecting only
		//..runs with a certain number of events
		hNEvent[usedRuns]      = (TH1D*)f->Get("hNEvents");
		cout<<hNEvent[usedRuns]->Integral()<<" evt."<<endl;
//		if(hNEvent[usedRuns]->Integral()>1000000)continue;

		TH2F *goodCells = (TH2F*)f->Get("2DChannelMap_Flag0");
		TH2F *deadCells = (TH2F*)f->Get("2DChannelMap_Flag1");
		TH2F *badCells  = (TH2F*)f->Get("2DChannelMap_Flag2");
		TH1F *hCellFlag = (TH1F*)f->Get("fhCellFlag");
		ampID[usedRuns]        = (TH2F*)f->Get("hCellAmplitude");
		hCellGoodMean[usedRuns]= (TH1D*)f->Get("hgoodMean");

		if(!badCells || !goodCells || !deadCells || !ampID[usedRuns] || !hCellFlag)
		{
			if(!badCells) Printf("2DChannelMap_Flag2 not found");
			if(!goodCells)Printf("2DChannelMap_Flag0 not found");
			if(!deadCells)Printf("2DChannelMap_Flag1 not found");
			if(!ampID[usedRuns]) Printf("hCellAmplitude not found");
			if(!hCellFlag)Printf("fhCellFlag not found");
			cout<<endl;
			continue;
		}

		nEventsVsRuns->SetBinContent(usedRuns+1,hNEvent[usedRuns]->Integral());
		ampID[usedRuns]        ->SetName(Form("hCellAmplitudeRun%i",usedRuns));
		ampIDCl[usedRuns]      = (TH2F*)ampID[usedRuns]->Clone(Form("hCellAmplitudeClRun%i",usedRuns));
		ampIDCl3Block[usedRuns]= (TH2F*)ampID[usedRuns]->Clone(Form("hCellAmplitudeCl3RunBlock%i",usedRuns));
		ampIDCl1Block[usedRuns]= (TH2F*)ampID[usedRuns]->Clone(Form("hCellAmplitudeCl1RunBlock%i",usedRuns));

		hCellGoodMean[usedRuns]->SetLineColor(kGreen);
		if(iFirstRun==0)
		{
			Sum2DOrig      = (TH2F*)ampID[usedRuns]->Clone("Sum2DOrig");
			Sum2DSingleMask= (TH2F*)ampID[usedRuns]->Clone("Sum2DSingleMask");
			Sum2D3BlockMask= (TH2F*)ampID[usedRuns]->Clone("Sum2D3BlockMask");
			Sum2DIdeal     = (TH2F*)ampID[usedRuns]->Clone("Sum2DIdealDistr");
			Sum2DOrig      ->Reset();
			Sum2DSingleMask->Reset();
			Sum2D3BlockMask->Reset();
			Sum2DIdeal     ->Reset();
			hgoodMean      = (TH1D*)hCellGoodMean[usedRuns]->Clone("MeanSpectrumAllRuns");
			hgoodMean     ->Reset();
		}
		//if(i<30)hCellGoodMean[i]->SetLineColor(RainbowColors[i]);
		hgoodMean->Add(hCellGoodMean[usedRuns]); //..add all good distributions to build a mean of all runs
		hCellGoodMean[usedRuns]->SetLineWidth(3);

		//..fill the histo bad cell vs. run number
		Int_t percBad=0;
		for(Int_t icell = 0; icell < noOfCells; icell++)
		{
			Int_t flag =  hCellFlag->GetBinContent(icell+1);
			//..dead
			if(flag == 1) hFlagvsRun[0]->Fill(icell, usedRuns, 1);
			//..bad or warm
			if(flag>1)    hFlagvsRun[1]->Fill(icell, usedRuns, 1); //fill, use the x, y values
			//..dead+bad
			if(flag>0)    hFlagvsRun[2]->Fill(icell, usedRuns, 1);
			if(flag>0)    percBad++;
		}
		if(1.0*percBad/noOfCells>0.3)cout<<"Problem in this run detected. Large number of bad+dead cells (>30%) - please double check!"<<endl;

		if(!hFlagNew)
		{
			hFlagNew = (TH2F*)goodCells->Clone(TString::Format("h2DChannelMapNew_FlagAll"));
			hFlagNew->Reset();
			hFlagNew->SetTitle("Selected flag greater than 0; cell column (#eta direction); cell raw (#phi direction)");
			hFlagNewClean = (TH2F*)goodCells->Clone(TString::Format("h2DChannelMapNew_FlagAllClean"));
			hFlagNewClean->Reset();
			hFlagNewClean->SetTitle("Selected flag greater than 0; cell column (#eta direction); cell raw (#phi direction)");
		}

		// Drawing histograms for each run
		//....................................
		cBad[usedRuns/totalperCv]->cd(usedRuns%totalperCv+1);
		SetHisto(badCells,"","",0);
		badCells->Draw("colz");
		TLatex* text = new TLatex(0.2,0.85,Form("Bad Cells - Run %i",RunId[i]));
		text->SetTextSize(0.06);
		text->SetNDC();
		text->SetTextColor(1);
		text->Draw();
		//....................................
		cGood[usedRuns/totalperCv]->cd(usedRuns%totalperCv+1);
		SetHisto(goodCells,"","",0);
		goodCells->Draw("colz");
		TLatex* text1 = new TLatex(0.2,0.85,Form("Good Cells - Run %i",RunId[i]));
		text1->SetTextSize(0.06);
		text1->SetNDC();
		text1->SetTextColor(1);
		text1->Draw();
		//....................................
		cDead[usedRuns/totalperCv]->cd(usedRuns%totalperCv+1);
		SetHisto(deadCells,"","",0);
		deadCells->Draw("colz");
		TLatex* text2 = new TLatex(0.2,0.85,Form("Dead Cells - Run %i",RunId[i]));
		text2->SetTextSize(0.06);
		text2->SetNDC();
		text2->SetTextColor(1);
		text2->Draw();

		//....................................
		cAmp[usedRuns/totalperCv]->cd(usedRuns%totalperCv+1)->SetLogy();
		TH1D* proj = ampID[usedRuns]->ProjectionX(TString::Format("hampIDProj_Run%d",RunId[i]));
		SetHisto(proj,"","",0);
		proj->SetLineColor(kCyan+2);
		proj->GetYaxis()->SetRangeUser(0.0000001,100);
		proj->Draw("hist");
		TLatex* text3 = new TLatex(0.2,0.85,Form("Amplitudes - Run %i",RunId[i]));
		text3->SetTextSize(0.06);
		text3->SetNDC();
		text3->SetTextColor(1);
		text3->Draw();
		//..create a summ version
		if(iFirstRun==0)projSum = ampID[usedRuns]->ProjectionX("hampIDProj_Sum");
		if(iFirstRun>0) projSum->Add(proj);
		Sum2DOrig->Add(ampID[usedRuns]);

		iFirstRun=1;
		usedRuns++;
	}
	nRunsUsed=usedRuns;
	//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	// Count number of cells that are bad in at least one run
	//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	Int_t nBadCells=0;
	for(Int_t ic = 0; ic < noOfCells; ic++)
	{
		TH1D *htmpCell = hFlagvsRun[2]->ProjectionY(TString::Format("hIDProj_cell%d", ic), ic+1, ic+1);
		Double_t sumRun = 0;
		for(Int_t ir = 0; ir < nRunsUsed ; ir++)
		{
			sumRun += htmpCell->GetBinContent(ir+1);
			if(htmpCell->GetBinContent(ir+1)==1)
			{
				//..mask the bad cells for this run
				//..Direction of amplitude (Checks energies from 0-nBins GeV)
				for (Int_t amp = 1; amp <= ampIDCl[ir]->GetNbinsX(); amp++)
				{
					ampIDCl[ir]->SetBinContent(amp,ic+1,0);
				}
			}
		}
		if(sumRun!=0)nBadCells++; //only count for the dead+bad case
	}

	hgoodMean->Scale(1.0/nRunsUsed);

    //..create an ideal 2D energy distribution for a later division
	//..helps to identify where cells have been unmasked and whether
	//..this was a good or bad unmasking desicion (e.g. creating spikes)
	for(Int_t eBin=0;eBin<Sum2DIdeal->GetNbinsX();eBin++)
	{
		Double_t binVal=hgoodMean->GetBinContent(eBin+1);
		for(Int_t icell=0;icell<Sum2DIdeal->GetNbinsY();icell++)
		{
			Sum2DIdeal->SetBinContent(eBin+1,icell+1,binVal);
		}
	}
	//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	// Draw masked cell amplitude by masking cells that were identified bad or dead in this specific run
	//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	for(Int_t ir = 0; ir < nRunsUsed ; ir++)
	{
		cAmp[ir/totalperCv]->cd(ir%totalperCv+1)->SetLogy();
		TH1D* proj = ampIDCl[ir]->ProjectionX(TString::Format("hampIDMaskedProj_Run%d",RunId[ir]));
		SetHisto(proj,"","",0);
		proj->SetLineColor(kSpring-2);
		proj->Draw("hist same");
		//..create a sum version
		if(ir==0)projSumC = ampIDCl[ir]->ProjectionX("hampIDProj_SumMasked");
		if(ir>0) projSumC->Add(proj);
		Sum2DSingleMask->Add(ampIDCl[ir]);
	}
	//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	// Draw summary histogram with dead and bad cells vs run
	//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	deadCellsVsRun    =hFlagvsRun[0]->ProjectionY("deadCellsVsRun");
	badCellsVsRun     =hFlagvsRun[1]->ProjectionY("badCellsVsRun");
	deadbadCellsVsRun =hFlagvsRun[2]->ProjectionY("badDeadCellsVsRun");
	cellSummaryCan->Divide(2);
	cellSummaryCan->cd(1);
	SetHisto(badCellsVsRun,"Run","No. of cells",0);
	badCellsVsRun->GetYaxis()->SetTitleOffset(1.7);
	badCellsVsRun->GetYaxis()->SetRangeUser(0,1300);
	badCellsVsRun->SetLineColor(kCyan+2);
	badCellsVsRun->SetLineWidth(2);
	badCellsVsRun->DrawCopy("hist");
	deadCellsVsRun->SetLineColor(kMagenta-2);
	deadCellsVsRun->SetLineWidth(2);
  	deadCellsVsRun->DrawCopy("same");

  	cellSummaryCan->cd(2);
  	SetHisto(nEventsVsRuns,"Run","No. of Events",0);
  	nEventsVsRuns->DrawCopy("hist");

  	cellSummaryCan2->Divide(2,2);
  	for(Int_t iRun=0;iRun<nRunsUsed;iRun++)
  	{
  		hBadVsEvent    ->Fill(nEventsVsRuns->GetBinContent(iRun+1),badCellsVsRun->GetBinContent(iRun+1));
  		hDeadBadVsEvent->Fill(nEventsVsRuns->GetBinContent(iRun+1),deadbadCellsVsRun->GetBinContent(iRun+1),1);
  		if(isHIRun(RunId[iRun]))
  		{
  			hBadVsEventHI    ->Fill(nEventsVsRuns->GetBinContent(iRun+1),badCellsVsRun->GetBinContent(iRun+1));
  			hDeadBadVsEventHI->Fill(nEventsVsRuns->GetBinContent(iRun+1),deadbadCellsVsRun->GetBinContent(iRun+1),1);
  		}
  	}
  	cellSummaryCan2->cd(1);
  	SetHisto(hBadVsEvent,"events in run","bad cells in run",0);
	hBadVsEvent->DrawCopy("colz");
 	cellSummaryCan2->cd(2);
	SetHisto(hDeadBadVsEvent,"events in run","bad+dead cells in run",0);
	hDeadBadVsEvent->DrawCopy("colz");
  	cellSummaryCan2->cd(3);
  	SetHisto(hBadVsEventHI,"events in HI run","bad cells in HI run",0);
  	hBadVsEventHI->DrawCopy("colz");
 	cellSummaryCan2->cd(4);
	SetHisto(hDeadBadVsEventHI,"events in HI run","bad+dead cells in HI run",0);
	hDeadBadVsEventHI->DrawCopy("colz");
	//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	//Draw bad & dead cells vs. ID
	//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	Color_t histCol=kCyan-8;
    //..Draw summary for all runs
	cFlagDeadBad->Divide(0,2);
	cFlagDeadBad->cd(1)->SetLeftMargin(0.05);
	cFlagDeadBad->cd(1)->SetRightMargin(0.05);
	SetHisto(hFlagvsRun[0],"dead cell ID","Run in list",1);
	hFlagvsRun[0]->SetFillColor(histCol);
	hFlagvsRun[0]->Draw("BOX");
	cFlagDeadBad->cd(2)->SetLeftMargin(0.04);
	cFlagDeadBad->cd(2)->SetRightMargin(0.04);
	SetHisto(hFlagvsRun[1],"bad cell ID","Run in list",1);
	hFlagvsRun[1]->SetFillColor(histCol);
	hFlagvsRun[1]->Draw("BOX");

	cFlagSum->Divide(0,2);
	cFlagSum->cd(1)->SetLeftMargin(0.04);
	cFlagSum->cd(1)->SetRightMargin(0.04);
	SetHisto(hFlagvsRun[2],"dead+bad cell ID","Run in list",1);
	hFlagvsRun[2]->GetXaxis()->SetRangeUser(0,8837);
	hFlagvsRun[2]->SetFillColor(histCol);
	hFlagvsRun[2]->DrawCopy("BOX");
	cFlagSum->cd(2)->SetLeftMargin(0.04);
	cFlagSum->cd(2)->SetRightMargin(0.04);
	SetHisto(hFlagvsRun[2],"dead+bad cell ID","Run in list",1);
	hFlagvsRun[2]->GetXaxis()->SetRangeUser(8837,17674);
	hFlagvsRun[2]->DrawCopy("BOX");

	//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	//..compress the histogram for visibility since
	//..90% of the space is filled by empty good cells
	//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	TH2F* CompressedAll = CompressHistogram(hFlagvsRun[2],noOfCells , nBadCells,nRunsUsed);
	cFlagSumCompAll->Divide(0,2);
	cFlagSumCompAll->cd(1)->SetLeftMargin(0.04);
	cFlagSumCompAll->cd(1)->SetRightMargin(0.04);
	SetHisto(CompressedAll,"certain dead+bad cell","Run in list",1);
	CompressedAll->GetXaxis()->SetRangeUser(0,nBadCells/2);
	CompressedAll->SetFillColor(histCol);
	CompressedAll->DrawCopy("BOX");
	cFlagSumCompAll->cd(2)->SetLeftMargin(0.04);
	cFlagSumCompAll->cd(2)->SetRightMargin(0.04);
	SetHisto(CompressedAll,"certain dead+bad cell","Run in list",1);
	CompressedAll->GetXaxis()->SetRangeUser(nBadCells/2,nBadCells+2);
	CompressedAll->DrawCopy("BOX");

	//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	//..................................................................
	// Find cells that are bad in a low fraction of runs
	//..................................................................
	//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	cout<<"    o Summary: "<<nBadCells<<" bad cells of "<<noOfCells<<" total cells and "<<nRunsUsed<<" runs"<<endl;
	cout<<"    o 1 bad run out of "<<nRunsUsed<<" is "<<1.0/nRunsUsed<<endl;
	Double_t percbad = 0.20; //..if a cell is only bad at 10% of the runs, declare it as good
	Double_t setBadCompletley=0.8; //..If a cell is bad in >80% of the runs then set it bad completley
	cout<<"    o Cells with "<<percbad<<" are double checked. These are  "<<nRunsUsed*percbad<<" runs"<<endl;
	cout<<"    o cell id with low fraction of bad runs:"<<endl;

	std::vector<Int_t> cellVector; //Filled with cells that onl have a small fraction of bad runs
	Double_t fracRun = 0, sumRun = 0;
	for(Int_t ic = 0; ic < noOfCells; ic++)
	{
		TH1D *htmpCell = hFlagvsRun[2]->ProjectionY(TString::Format("hIDProj_cell%d", ic), ic+1, ic+1);
		fracRun = 0, sumRun = 0;
		for(Int_t ir = 0; ir < nRunsUsed ; ir++)
		{
			sumRun += htmpCell->GetBinContent(ir+1);
		}
		fracRun = sumRun/(Double_t)(nRunsUsed);

		//..loose selection criteria to remove the runs with zero entries
		if(fracRun>0)
		{
			Int_t cellColumn=0, cellRow=0;
			Int_t cellColumnAbs=0, cellRowAbs=0;
			Int_t trash = 0 ;
			fCaloUtils->GetModuleNumberCellIndexesAbsCaloMap(ic,0,cellColumn,cellRow,trash,cellColumnAbs,cellRowAbs);
			hFlagNew->Fill(cellColumnAbs, cellRowAbs, fracRun*100);
			//..If a cell is bad in >80% of the runs then set it bad completley
			if(fracRun>setBadCompletley)
			{
				for(Int_t j = 0 ; j < nRunsUsed; j++)
				{
					hFlagvsRun[2]->SetBinContent(ic+1,j+1,1);
				}
			}
			//..If a cell is bad in a low fraction of runs double check it
			if(fracRun<percbad)
			{
				cout<<ic<<", "<<flush;
				cellVector.push_back(ic);
			}
		}
	}
	cout<<endl;
	cout<<"    o In total "<<cellVector.size()<<" cells fall under this category"<<endl;
	//..................................................................
	// Plot cells with low fraction if bad runs
	// Double checks if they are really bad and re-includes them
	//..................................................................
	TString pdfName  = Form("%s/RunByRunSummary/%s_LowFractionCells",analysisOutput.Data(),listName.Data());
	PlotLowFractionCells(pdfName,cellVector,hFlagvsRun,nRunsUsed,ampID,hCellGoodMean);

	//..................................................................
	// Draw masked cell amplitude by masking cells that were identified bad or dead in a certain runblock
	//..................................................................
	for(Int_t ic = 0; ic < noOfCells; ic++)
	{
		TH1D *htmpCellAllRuns     =hFlagvsRun[2]->ProjectionY(TString::Format("hIDProj_cell%d", ic), ic+1, ic+1);
		Double_t integralBlock1   =htmpCellAllRuns->Integral(0,splitRuns1);
		Double_t integralBlock2   =htmpCellAllRuns->Integral(splitRuns1+1,splitRuns2);
		Double_t integralBlock3   =htmpCellAllRuns->Integral(splitRuns2+1,nRunsUsed);

		//..manually mask cells
		if(ic==badcellsBlock1.at(0))integralBlock1=1;
		if(ic==badcellsBlock2.at(0))integralBlock2=1;
		if(ic==badcellsBlock3.at(0))integralBlock3=1;

		Double_t integralBlock3Sum=integralBlock1+integralBlock2+integralBlock3;
		//..only if the cell is bad in 1 run we will start the
		//..masking proceedure
		if(integralBlock3Sum>0)
		{
			for(Int_t ir = 0; ir < nRunsUsed ; ir++)
			{
				if((integralBlock1>0 && ir<=splitRuns1) ||
				   (integralBlock2>0 && ir>splitRuns1 && ir<=splitRuns2) ||
				   (integralBlock3>0 && ir>splitRuns2))
				{
					for (Int_t amp = 1; amp <= ampIDCl3Block[ir]->GetNbinsX(); amp++)
					{
						//..mask the cell, if it is once bad in the specific block
						ampIDCl3Block[ir]->SetBinContent(amp,ic+1,0);
					}
				}
                //..mask the cell if it is bad in any run in the period
				for (Int_t amp = 1; amp <= ampIDCl1Block[ir]->GetNbinsX(); amp++)
				{
					ampIDCl1Block[ir]->SetBinContent(amp,ic+1,0);
				}
			}
		}
	}
	//..build projections of amplitudes
	for(Int_t ir = 0; ir < nRunsUsed ; ir++)
	{
		cAmp[ir/totalperCv]->cd(ir%totalperCv+1)->SetLogy();
		TH1D* proj      = ampIDCl3Block[ir]->ProjectionX(TString::Format("hampIDMaskedCleanedProj_Run%d",RunId[ir]));
		TH1D* projBlock = ampIDCl1Block[ir]->ProjectionX(TString::Format("hampIDMaskedCleaned1blockProj_Run%d",RunId[ir]));
		SetHisto(proj,"","",0);
		proj->SetLineColor(kBlue-1);
		proj->Draw("hist same");
		//..create a summ version
		if(ir==0)projSumC3Blocks = ampIDCl3Block[ir]->ProjectionX("hampIDProj_SumMaskedMultipleBlock");
		if(ir==0)projSumC3BlocksA= ampIDCl3Block[ir]->ProjectionX("hampIDProj_SumMaskedMultipleBlockA");
		if(ir==0)projSumC3BlocksB= ampIDCl3Block[ir]->ProjectionX("hampIDProj_SumMaskedMultipleBlockB");
		if(ir==0)projSumC3BlocksC= ampIDCl3Block[ir]->ProjectionX("hampIDProj_SumMaskedMultipleBlockC");
		if(ir>0)                                    projSumC3Blocks->Add(proj);
		if(ir>0 && ir<=splitRuns1)                  projSumC3BlocksA->Add(proj);
		if(ir>0 && ir>splitRuns1 && ir<=splitRuns2) projSumC3BlocksB->Add(proj);
		if(ir>0 && ir>splitRuns2)                   projSumC3BlocksC->Add(proj);
		if(ir==0)projSumC1Block = ampIDCl1Block[ir]->ProjectionX("hampIDProj_SumMaskedOneBlock");
		if(ir>0) projSumC1Block ->Add(projBlock);
		Sum2D3BlockMask->Add(ampIDCl3Block[ir]);
	}
	//..................................................................
	// Draw the cells masked in each respective run and summed amplitudes
	//..................................................................
	cAmpSum->Divide(2);
	cAmpSum->cd(1)->SetLogy();
	SetHisto(projSum,"","No. of hits/event",0);
	projSum->GetYaxis()->SetRangeUser(0.0000001,1000);
	projSum->SetLineColor(kCyan+3);
	projSum->SetLineWidth(2);
	projSum->Draw("hist");

	projSumC->SetLineColor(kSpring-2);
	projSumC->SetLineWidth(2);
	projSumC->DrawCopy("hist same");

	projSumC3Blocks->SetLineColor(kCyan-8);
	projSumC3Blocks->SetLineWidth(2);
	projSumC3Blocks->DrawCopy("hist same");
	projSumC1Block->SetLineColor(2);
	projSumC1Block->DrawCopy("hist same");
	projSumC3BlocksA->SetLineColor(4);
	projSumC3BlocksA->DrawCopy("hist same");
	projSumC3BlocksB->SetLineColor(5);
	projSumC3BlocksB->DrawCopy("hist same");
	projSumC3BlocksC->SetLineColor(6);
	projSumC3BlocksC->DrawCopy("hist same");
	TLegend *legSum = new TLegend(0.35,0.70,0.55,0.85);
	legSum->AddEntry(projSum,"Original engery distr.","l");
	legSum->AddEntry(projSumC,"Cells masked in each run","l");
	legSum->AddEntry(projSumC3Blocks,"Cells reincluded and masked in 3 blocks","l");
	legSum->SetBorderSize(0);
	legSum->SetTextSize(0.03);
	legSum->Draw("same");

	cAmpSum->cd(2);
	projSumC3Blocks->Divide(projSumC);
	SetHisto(projSumC3Blocks,"","block masked/single masked",0);
	projSumC3Blocks->SetLineColor(30);
	projSumC3Blocks->DrawCopy("hist");
	projSumC1Block->SetLineColor(4);
	projSumC1Block->Divide(projSumC);
	projSumC1Block->DrawCopy("same");

	projSumC1Block->SetLineColor(6);
	projSumC1Block->Divide(projSumC3Blocks);
	projSumC1Block->DrawCopy("same");

	TLegend *legRatio = new TLegend(0.35,0.70,0.55,0.85);
	legRatio->AddEntry(projSumC1Block,"re-incl & masked as one big block","l");
	legRatio->AddEntry(projSumC3Blocks,"re-incl & masked in three blocks","l");
	legRatio->AddEntry(projSumC1Block,"ratio 1Block/3Blocks","l");
	legRatio->SetBorderSize(0);
	legRatio->SetTextSize(0.03);
	legRatio->Draw("same");

	//..................................................................
	// Draw the Amp vs E histogram and the ratio to the masked 2D one
	//..................................................................
	cAmpSum2D->Divide(2);
	cAmpSum2D->cd(1)->SetLogz();
	SetHisto(Sum2D3BlockMask,"","",0);
	Sum2D3BlockMask->GetZaxis()->SetRangeUser(10e-8,10);
	Sum2D3BlockMask->DrawCopy("colz");

	cAmpSum2D->cd(2);
	Sum2D3BlockMask->Divide(Sum2DIdeal);
	SetHisto(Sum2D3BlockMask,"","",0);
	Sum2D3BlockMask->GetZaxis()->SetRangeUser(2,20);
	Sum2D3BlockMask->DrawCopy("colz");

	//..................................................................
	//..Plot the cleaned up histogram again
	//..................................................................
	cFlagSumCleaned->Divide(0,2);
	cFlagSumCleaned->cd(1)->SetLeftMargin(0.04);
	cFlagSumCleaned->cd(1)->SetRightMargin(0.04);
	SetHisto(hFlagvsRun[2],"dead+bad cell ID (cleaned)","Run in list",1);
	hFlagvsRun[2]->GetXaxis()->SetRangeUser(0,8837);
	hFlagvsRun[2]->SetFillColor(histCol);
	hFlagvsRun[2]->DrawCopy("BOX");
	cFlagSumCleaned->cd(2)->SetLeftMargin(0.04);
	cFlagSumCleaned->cd(2)->SetRightMargin(0.04);
	SetHisto(hFlagvsRun[2],"dead+bad cell ID (cleaned)","Run in list",1);
	hFlagvsRun[2]->GetXaxis()->SetRangeUser(8837,17674);
	hFlagvsRun[2]->DrawCopy("BOX");


	//..................................................................
	// Draw summary histogram with dead and bad cells vs run
	//..................................................................
	badCellsVsRunC  =hFlagvsRun[1]->ProjectionY("badCellsVsRunC");
	cellSummaryCan->cd(1);
	badCellsVsRunC->SetLineColor(kGreen-9);
	badCellsVsRunC->SetLineWidth(2);
	badCellsVsRunC->SetFillColorAlpha(10,0);
	badCellsVsRunC->DrawCopy("same");

	//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	//..compress the histogram for visibility since
	//..90% of the space is filled by empty good cells
	//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	cout<<endl;
	TH2F* CompressedClean = CompressHistogram(hFlagvsRun[2],noOfCells , nBadCells,nRunsUsed);

	cFlagSumCompClean->Divide(0,2);
	cFlagSumCompClean->cd(1)->SetLeftMargin(0.04);
	cFlagSumCompClean->cd(1)->SetRightMargin(0.04);
	SetHisto(CompressedClean,"certain dead+bad cell (cleaned)","Run in list",1);
	CompressedClean->GetXaxis()->SetRangeUser(0,nBadCells/2);
	CompressedClean->SetFillColor(histCol);
	CompressedClean->DrawCopy("BOX");
	cFlagSumCompClean->cd(2)->SetLeftMargin(0.04);
	cFlagSumCompClean->cd(2)->SetRightMargin(0.04);
	SetHisto(CompressedClean,"certain dead+bad cell (cleaned)","Run in list",1);
	CompressedClean->GetXaxis()->SetRangeUser(nBadCells/2,nBadCells+2);
	CompressedClean->DrawCopy("BOX");

	//..................................................................
	//..build and draw a new 2D histogram after all the cleaning/resetting was done
	//..................................................................
	for(Int_t ic = 0; ic < noOfCells; ic++)
	{
		TH1D *htmpCell = hFlagvsRun[2]->ProjectionY(TString::Format("hIDProj_cell%d", ic), ic+1, ic+1);
		fracRun = 0, sumRun = 0;
		for(Int_t ir = 0; ir < nRunsUsed ; ir++)
		{
			sumRun += htmpCell->GetBinContent(ir+1);
		}
		fracRun = sumRun/(Double_t)(nRunsUsed);
		if(fracRun>0)
		{
			Int_t cellColumn=0, cellRow=0;
			Int_t cellColumnAbs=0, cellRowAbs=0;
			Int_t trash = 0 ;
			fCaloUtils->GetModuleNumberCellIndexesAbsCaloMap(ic,0,cellColumn,cellRow,trash,cellColumnAbs,cellRowAbs);
			hFlagNewClean->Fill(cellColumnAbs, cellRowAbs, fracRun*100);
		}
	}
    //..2D fraction
    cFlagNew->Divide(2);
    cFlagNew->cd(1);
	SetHisto(hFlagNew,"","",0);
	hFlagNew->Draw("colz");
    cFlagNew->cd(2);
	SetHisto(hFlagNewClean,"","",0);
	hFlagNewClean->Draw("colz");

	//..................................................................
	// Print out an overview of this run-by-run analysis
	//..................................................................
	Int_t bcBolckSum=0;

	TH1D *htmpCellSum = hFlagvsRun[2]->ProjectionX("countSum");
	for(Int_t ic = 0; ic < noOfCells; ic++)
	{
		if(htmpCellSum->GetBinContent(ic+1)>0)bcBolckSum++;
	}

    cout<<"...................................."<<endl;
	cout<<"Summary: "<<nBadCells<<" bad cells of "<<noOfCells<<" total cells and "<<nRunsUsed<<" runs"<<endl;
    cout<<"o All bad cells: "<<nBadCells<<endl;
    cout<<"o low fraction bad cells: "<<cellVector.size()<<endl;
    cout<<"o bad cells after reinclusion: "<<bcBolckSum<<endl;
    cout<<"o rescued cells: "<<nBadCells-bcBolckSum<<endl;
    cout<<"o These are: "<<100*(nBadCells-bcBolckSum)/cellVector.size()<<"% of low frac cells"<<endl;
    cout<<"...................................."<<endl;

	//_______________________________________________________________________
	//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	//. . Save histograms to file
	//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	//_______________________________________________________________________
	TString fileName       = Form("%s/RunByRunSummary/%s_Results.root",analysisOutput.Data(),listName.Data());
	TFile* rootFile        = new TFile(fileName,"recreate");

	for(Int_t ic = 0; ic<nCv; ic++)
	{
		rootFile->WriteObject(cBad [ic],cBad [ic]->GetName());
		rootFile->WriteObject(cGood[ic],cGood[ic]->GetName());
		rootFile->WriteObject(cDead[ic],cDead[ic]->GetName());
		rootFile->WriteObject(cAmp [ic],cAmp [ic]->GetName());
	}
	rootFile->WriteObject(cellSummaryCan,cellSummaryCan->GetName());
	rootFile->WriteObject(cAmpSum,cAmpSum->GetName());
	rootFile->WriteObject(cFlagDeadBad,cFlagDeadBad->GetName());
	rootFile->WriteObject(cFlagSum,cFlagSum->GetName());
	rootFile->WriteObject(cFlagSumCleaned,cFlagSumCleaned->GetName());
	rootFile->WriteObject(cFlagSumCompAll,cFlagSumCompAll->GetName());
	rootFile->WriteObject(cFlagSumCompClean,cFlagSumCompClean->GetName());
	rootFile->WriteObject(cFlagNew,cFlagNew->GetName());
	rootFile->WriteObject(cAmpSum2D,cAmpSum2D->GetName());
	//..histograms
	rootFile->WriteObject(hFlagNew,hFlagNew->GetName());
	rootFile->WriteObject(CompressedAll,CompressedAll->GetName());
	rootFile->WriteObject(CompressedClean,CompressedClean->GetName());
	rootFile->WriteObject(projSum,projSum->GetName());
	rootFile->WriteObject(projSumC,projSumC->GetName());
	rootFile->WriteObject(projSumC3Blocks,projSumC3Blocks->GetName());
	rootFile->WriteObject(nEventsVsRuns,nEventsVsRuns->GetName());

	//..plot the canvases of cells into a .pdf file
	TString badCanvasName =Form("%s/RunByRunSummary/%s_BadCells.pdf", analysisOutput.Data(),listName.Data());
	TString goodCanvasName=Form("%s/RunByRunSummary/%s_GoodCells.pdf", analysisOutput.Data(),listName.Data());
	TString deadCanvasName=Form("%s/RunByRunSummary/%s_DeadCells.pdf", analysisOutput.Data(),listName.Data());
	TString ampCanvasName =Form("%s/RunByRunSummary/%s_Amplitudes.pdf", analysisOutput.Data(),listName.Data());
	TString summaryPDF    =Form("%s/RunByRunSummary/%s_FigureCollection.pdf", analysisOutput.Data(),listName.Data());

	for(Int_t ic = 0; ic<nCv; ic++)
	{
		if(ic==0 && nCv>1)
		{
			//..first pad
			cBad [ic] ->Print(Form("%s(",badCanvasName.Data()));
			cGood[ic] ->Print(Form("%s(",goodCanvasName.Data()));
			cDead[ic] ->Print(Form("%s(",deadCanvasName.Data()));
			cAmp [ic] ->Print(Form("%s(",ampCanvasName.Data()));
		}
		else if(ic==(nCv-1) && nCv>1)
		{
			//..last pad
			cBad [ic] ->Print(Form("%s)",badCanvasName.Data()));
			cGood[ic] ->Print(Form("%s)",goodCanvasName.Data()));
			cDead[ic] ->Print(Form("%s)",deadCanvasName.Data()));
			cAmp [ic] ->Print(Form("%s)",ampCanvasName.Data()));
		}
		else
		{
			//..all pads in between
			cBad [ic] ->Print(Form("%s",badCanvasName.Data()));
			cGood[ic] ->Print(Form("%s",goodCanvasName.Data()));
			cDead[ic] ->Print(Form("%s",deadCanvasName.Data()));
			cAmp [ic] ->Print(Form("%s",ampCanvasName.Data()));
		}
	}

	//..Add figures to the summary canvas
	cellSummaryCan   ->Print(Form("%s(",summaryPDF.Data()));
	cellSummaryCan2  ->Print(Form("%s(",summaryPDF.Data()));
	cAmpSum          ->Print(Form("%s",summaryPDF.Data()));
	cAmpSum2D        ->Print(Form("%s",summaryPDF.Data()));
	cFlagNew         ->Print(Form("%s",summaryPDF.Data()));
	cFlagDeadBad     ->Print(Form("%s",summaryPDF.Data()));
	cFlagSum         ->Print(Form("%s",summaryPDF.Data()));
	cFlagSumCleaned  ->Print(Form("%s",summaryPDF.Data()));
	cFlagSumCompAll  ->Print(Form("%s",summaryPDF.Data()));
	cFlagSumCompClean->Print(Form("%s)",summaryPDF.Data()));

}
//
//________________________________________________________________________
void BuildMaxMinHisto(TH1D* inHisto, TH1D* minHist,TH1D* maxHist)
{
	Double_t ref;
	Double_t min;
	Double_t max;
	for(Int_t bin=0;bin<50;bin++)
	{
		ref = inHisto->GetBinContent(bin);
	    max = maxHist->GetBinContent(bin);
	    min = minHist->GetBinContent(bin);
	    if((ref!=0 && ref<min) || min==0)minHist->SetBinContent(bin,ref);
	    if(ref>max)maxHist->SetBinContent(bin,ref);
	}
}
//
//________________________________________________________________________
Bool_t IsItReallyBadRatio(TH1D* minHistoRatio,TH1D* maxHistoRatio,TH1D* meanHistoRatio, TString& crit)
{
	//minHistoRatio  ->Smooth();
	//maxHistoRatio  ->Smooth();
	//meanHistoRatio ->Smooth();

	//..do the check only until 1 GeV because
	//..then we are running out of statistic
	/*Int_t OneGeVBin  =minHistoRatio->FindBin(1.47);  //=30bins
	Int_t HalfGeVBin =minHistoRatio->FindBin(0.73);  //=15bins
	Int_t ThirdGeVBinA =minHistoRatio->FindBin(0.35);
	Int_t ThirdGeVBinB =minHistoRatio->FindBin(0.7);
	Int_t ThirdGeVBinC =minHistoRatio->FindBin(1.05);
	*/
    Double_t mBlock1=0,mBlock2=0,mBlock3=0,mBlock4=0,mBlock5=0;
    Double_t zBlock1=0,zBlock2=0,zBlock3=0,zBlock4=0,zBlock5=0;

	Double_t minMean=0,    maxMean=0,    meanMean=0;
	Double_t zminMean=0,   zmaxMean=0,   zmeanMean=0;
/*	Double_t meanMean1=0,  meanMean2=0;
	Double_t meanMeanA3=0, meanMeanB3=0, meanMeanC3=0;
	Int_t zeroBinsHalf1=0, zeroBinsHalf2=0;
	Int_t zeroBinsThird1=0, zeroBinsThird2=0, zeroBinsThird3=0;
*/
	//for(Int_t bin=0;bin<30;bin++)
	for(Int_t bin=3;bin<53;bin++)
	{
		if(bin<33)
		{
			minMean  += minHistoRatio ->GetBinContent(bin);
			maxMean  += maxHistoRatio ->GetBinContent(bin);
			meanMean += meanHistoRatio->GetBinContent(bin);
			if(minHistoRatio->GetBinContent(bin)==0)zminMean++;
			if(maxHistoRatio->GetBinContent(bin)==0)zmaxMean++;
			if(meanHistoRatio->GetBinContent(bin)==0)zmeanMean++;
		}

/*		if(bin<18)      meanMean1+=meanHistoRatio->GetBinContent(bin);
		else if(bin<33) meanMean2+=meanHistoRatio->GetBinContent(bin);

*/
		/*		if(bin<10)     meanMeanA3+=meanHistoRatio->GetBinContent(bin+1);
		else if(bin<20)meanMeanB3+=meanHistoRatio->GetBinContent(bin+1);
		else           meanMeanC3+=meanHistoRatio->GetBinContent(bin+1);
		 */
		if(bin<13)      mBlock1+=meanHistoRatio->GetBinContent(bin);
		else if(bin<23) mBlock2+=meanHistoRatio->GetBinContent(bin);
		else if(bin<33) mBlock3+=meanHistoRatio->GetBinContent(bin);
		else if(bin<43) mBlock4+=meanHistoRatio->GetBinContent(bin);
		else if(bin<53) mBlock5+=meanHistoRatio->GetBinContent(bin);
		//..count zero bins
		if(meanHistoRatio->GetBinContent(bin)==0)
		{
			if(bin<13)      zBlock1++;
			else if(bin<23) zBlock2++;
			else if(bin<33) zBlock3++;
			else if(bin<43) zBlock4++;
			else if(bin<53) zBlock5++;
		}
	/*	//..correct for zero bins
		if(meanHistoRatio->GetBinContent(bin+1)==0 && bin<15) zeroBinsHalf1++;
		else if(meanHistoRatio->GetBinContent(bin+1)==0)      zeroBinsHalf2++;

		if(meanHistoRatio->GetBinContent(bin+1)==0 && bin<10)     zeroBinsThird1++;
		else if(meanHistoRatio->GetBinContent(bin+1)==0 && bin<20)zeroBinsThird2++;
		else if(meanHistoRatio->GetBinContent(bin+1)==0)          zeroBinsThird3++;
	*/}
	//.......................................
	//..correct for zero bins
	//.......................................
/*
	if(zeroBinsHalf1!=0)meanMean1=meanMean1/(1.0-1.0*zeroBinsHalf1/15);
	if(zeroBinsHalf2!=0)meanMean2=meanMean2/(1.0-1.0*zeroBinsHalf2/15);

	if(zeroBinsThird1!=0)meanMeanA3=meanMeanA3/(1.0-1.0*zeroBinsThird1/10);
	if(zeroBinsThird2!=0)meanMeanB3=meanMeanB3/(1.0-1.0*zeroBinsThird2/10);
	if(zeroBinsThird3!=0)meanMeanC3=meanMeanC3/(1.0-1.0*zeroBinsThird3/10);
*/
	if(zBlock1<10 && zBlock1!=0)mBlock1=mBlock1/(1.0-1.0*zBlock1/10);
	if(zBlock2<10 && zBlock2!=0)mBlock2=mBlock2/(1.0-1.0*zBlock2/10);
	if(zBlock3<10 && zBlock3!=0)mBlock3=mBlock3/(1.0-1.0*zBlock3/10);
	if(zBlock4<10 && zBlock4!=0)mBlock4=mBlock4/(1.0-1.0*zBlock4/10);
	if(zBlock5<10 && zBlock5!=0)mBlock5=mBlock5/(1.0-1.0*zBlock5/10);

	if(zminMean!=0)minMean  =minMean/(1.0-1.0*zminMean/30);
	if(zmaxMean!=0)maxMean  =maxMean/(1.0-1.0*zmaxMean/30);
	if(zmeanMean!=0)meanMean=meanMean/(1.0-1.0*zmeanMean/30);

	//..if more than half of the bins in the block were 0 exclude block
	if(zBlock1>5)mBlock1=0;
	if(zBlock2>5)mBlock2=0;
	if(zBlock3>5)mBlock3=0;
	if(zBlock4>5)mBlock4=0;
	if(zBlock5>5)mBlock5=0;
	//.......................................
	//..check criteria
	//.......................................

	//..bad channel is 5times lower than max distr.
	crit = "spectr. too low";
	if(meanMean/30.0>5) return 1;
	if(minMean/30.0 >5) return 1;   //5 times lower than the lowest run
	//..bad channel is 5times higher than max distr.
	crit = "spectr. too high";
	if(meanMean/30.0<0.2) return 1;
	if(maxMean/30.0 <0.2) return 1; //5 times higher than the highest run

	//..if there is a slope down (should be more than 10% decrease)
	crit = "slope down";
	Int_t down=0;
	if(mBlock1!=0 && mBlock2!=0 && mBlock1>mBlock2 && mBlock1>1.4*mBlock2)down++;
	if(mBlock2!=0 && mBlock3!=0 && mBlock2>mBlock3 && mBlock2>1.4*mBlock3)down++;
	if(mBlock3!=0 && mBlock4!=0 && mBlock3>mBlock4 && mBlock3>1.4*mBlock4)down++;
	if(mBlock4!=0 && mBlock5!=0 && mBlock4>mBlock5 && mBlock4>1.4*mBlock5)down++;
	if(down>=2)return 1; //..means at least three blocks have to be staggered

	crit = "slope up";
	Int_t up=0;
	if(mBlock1!=0 && mBlock2 !=0 && mBlock1<mBlock2 && 1.4*mBlock1<mBlock2)up++;
	if(mBlock2!=0 && mBlock3 !=0 && mBlock2<mBlock3 && 1.4*mBlock2<mBlock3)up++;
	if(mBlock3!=0 && mBlock4 !=0 && mBlock3<mBlock4 && 1.4*mBlock3<mBlock4)up++;
	if(mBlock4!=0 && mBlock5 !=0 && mBlock4<mBlock5 && 1.4*mBlock4<mBlock5)up++;
	if(up>=2)return 1; //..means at least three blocks have to be staggered

	//..if there is a step at 2.1 GeV
    crit = "step 2.1 GeV";
    if(mBlock4!=0 && mBlock5!=0 && mBlock4>mBlock5 && mBlock4>20*mBlock5) return 1;

    //..if there is a steep step at 1.1 GeV (can only be performed if this is not dominated by "zero" bins)
//    crit = "step 1.1";
//    if(zeroBinsThird3<5 && meanMeanB3>meanMeanC3 && meanMeanB3>3.5*meanMeanC3) return 1;
	//..is good
    crit = "";
	return 0;
}
//
//________________________________________________________________________
void PlotLowFractionCells(TString pdfName, std::vector<Int_t> cellVector,TH2F* badVsCell[],Int_t nRuns,TH2F* ampID[],TH1D* hCellGoodMean[])
{
	Int_t nRescuableChannels=cellVector.size();
	Int_t totalperCv = 16;
	Int_t nPad = TMath::Sqrt(totalperCv);
	Int_t nCv  = nRescuableChannels/totalperCv+1;
	Int_t lastGood=0;

	TLatex* text = new TLatex(0.45,0.6,"*Indeed bad*");//..
	text->SetTextSize(0.07);
	text->SetTextColor(kAzure-8);
	text->SetNDC();

	TH1D *maxHisto  = ampID[0]->ProjectionX("hMaxCells", 1, 1);
	TH1D *minHisto  = ampID[0]->ProjectionX("hMinCells", 1, 1);
	TH1D* hgoodMean = ampID[0]->ProjectionX("hMeanofRuns", 1, 1);
	if(nCv<1)nCv=1;

	cout<<"    o create: "<<nCv<<" Canvases with "<<nPad*nPad<<" pads"<<endl;
	//..to compare specific cells over the runs
	TCanvas **cComp     = new TCanvas*[nCv];
	TCanvas **cCompDiv  = new TCanvas*[nCv];
	for(Int_t i=0;i<nCv;i++)
	{
		cComp[i]    = new TCanvas(TString::Format("CompareGood%d", i), TString::Format("V) Candidates (%d/%d)", i+1, nCv), 1000,750);
		cComp[i]    ->Divide(nPad,nPad,0.001,0.001);
		cCompDiv[i] = new TCanvas(TString::Format("CompareGood Ratio%d", i), TString::Format("V) Candidates Ratio (%d/%d)", i+1, nCv), 1000,750);
		cCompDiv[i] ->Divide(nPad,nPad,0.001,0.001);
	}

	Int_t notBadCounter=0;
	//TH1F** hCellSpectr = new TH1F*[nRescuableChannels];
	for(Int_t cell=0;cell< (Int_t)cellVector.size();cell++)
	{
		if(cell%400==0)cout<<"cell No."<<cell<<endl;
		if(cell%20==0) cout<<"."<<flush;
		maxHisto->Reset();
		minHisto->Reset();
		hgoodMean->Reset();
		TH1D* declaredBad;
		std::vector<TH1D*> badHistVector;
		Int_t badRun=-1;
		for(Int_t i = 0 ; i < nRuns ; i++)
		{
			TH1D *htmpCell = ampID[i]->ProjectionX(TString::Format("hIDProj_cell%dRun%i", cellVector.at(cell),i), cellVector.at(cell)+1, cellVector.at(cell)+1);
			htmpCell->SetLineColor(kGray+1);
			if(badVsCell[2]->GetBinContent(cellVector.at(cell)+1,i+1)==1)
			{
				htmpCell->SetLineColor(2);
				htmpCell->SetFillColor(2);
				htmpCell->SetFillStyle(3002);
				declaredBad = (TH1D*)htmpCell->Clone("saveForLater");
				badHistVector.push_back(htmpCell);
				badRun=i;
				if(htmpCell->Integral()==0)cout<<"cell "<<cell<<" is dead for run: "<<i<<endl;
			}
			else
			{
				BuildMaxMinHisto(htmpCell,minHisto,maxHisto);
				hgoodMean->Add(htmpCell);
			}
			//..go to the last pad and draw the mean of all good cell distribution
			//..for all the tested runs
			if(i==0)
			{
				cComp[nCv-1]   ->cd(nPad*nPad)->SetLogy();
				SetHisto(hCellGoodMean[i],"","Number of Hits",0);
				hCellGoodMean[i]->GetXaxis()->SetRangeUser(0,3);
				hCellGoodMean[i]->Draw("hist");
			}
			else
			{
				cComp[nCv-1]   ->cd(nPad*nPad)->SetLogy();
				hCellGoodMean[i]->Draw("same hist");
			}
			//..go to pads and draw cell in all runs
			lastGood=(cell-notBadCounter)/totalperCv;//..you can overwrite good canvases to lessen the amount of canvases
			if(i==0)
			{
				cComp[(cell-notBadCounter)/totalperCv]->cd(((cell-notBadCounter)%totalperCv)+1)->SetLogy();
				SetHisto(htmpCell,Form("Energy of cell %i",cellVector.at(cell)),"Number of Hits",0);
				htmpCell->GetXaxis()->SetRangeUser(0,3);
				htmpCell->Draw("hist");
			}
			else
			{
				cComp[(cell-notBadCounter)/totalperCv]->cd(((cell-notBadCounter)%totalperCv)+1)->SetLogy();
				htmpCell->DrawCopy("same hist");
			}
		}//..end of loop over runs
		hgoodMean->Scale(1.0/nRuns);

		hCellGoodMean[badRun]->DrawCopy("same hist"); //..draw the mean of all good cells for the run where the cell was bad
		minHisto->SetLineColor(1);
		maxHisto->SetLineColor(1);
		minHisto->DrawCopy("same");     //..draw the combined minimum of that cell for all the runs
		maxHisto->DrawCopy("same");     //..draw the combined maximum of that cell for all the runs

		//..Draw bad again
		for(Int_t j=0;j< (Int_t)badHistVector.size();j++)
		{
			badHistVector[j]->SetLineColor(2);
			badHistVector[j]->DrawCopy("same");
		}

		TLegend *leg = new TLegend(0.35,0.65,0.75,0.85);
		leg->AddEntry(hCellGoodMean[badRun],"mean of good cells","l");
		leg->AddEntry(declaredBad,"Cell in the ''bad'' run","l");
		leg->AddEntry(maxHisto,"max and min values","l");
		leg->SetBorderSize(0);
		leg->SetTextSize(0.07);
		leg->Draw("same");

		//. . . . . . . . . . . . . . . . . . . . . . . . .
		//..Fill the ratio canvas
		cCompDiv[(cell-notBadCounter)/totalperCv]->cd(((cell-notBadCounter)%totalperCv)+1);

		//..Draw bad again
		Int_t runTry=0;
		Bool_t badRatio=0;
		TString failCrit;
		for(Int_t j=0;j< (Int_t)badHistVector.size();j++)
		{
			TH1D* minHistoCopy =(TH1D*)minHisto->Clone("minHistoCopy");
			TH1D* maxHistoCopy =(TH1D*)maxHisto->Clone("maxHistoCopy");
			TH1D* meanHistoCopy=(TH1D*)hgoodMean->Clone("meanHistoCopy");
			minHistoCopy ->Divide(badHistVector[j]);
			maxHistoCopy ->Divide(badHistVector[j]);
			meanHistoCopy->Divide(badHistVector[j]);

			SetHisto(maxHistoCopy,Form("Energy of cell %i",cellVector.at(cell)),"max,min,mean / bad cell",0);
			maxHistoCopy->GetXaxis()->SetRangeUser(0,3);
			maxHistoCopy->DrawCopy("hist");
			minHistoCopy->DrawCopy("same hist");
			meanHistoCopy->SetLineColor(kBlue-8);
			meanHistoCopy->DrawCopy("same hist");
			//cout<<"cell: "<<cell<<endl;
			badRatio = IsItReallyBadRatio(minHistoCopy,maxHistoCopy,meanHistoCopy,failCrit);
			runTry=j;
			if(badRatio==1)break; //if its bad for one of the runs its enough
		}
		if(badRatio==1)
		{
			gPad->SetFrameFillColor(kRed-10);
			text->DrawLatexNDC(0.45,0.8,Form("#splitline{(%i)Excluded: }{%s}",runTry,(const char*)failCrit));;
		}
		else
		{
			//de-mask cells not declared as bad. (re-inclusion)
			for(Int_t j = 0 ; j < nRuns ; j++)
			{
				badVsCell[0]->SetBinContent(cellVector.at(cell)+1,j+1,0);
				badVsCell[1]->SetBinContent(cellVector.at(cell)+1,j+1,0);
				badVsCell[2]->SetBinContent(cellVector.at(cell)+1,j+1,0);
			}
		}
		if(cell==(Int_t)cellVector.size()-1)text->SetText(0.45,0.8,"test");

	}
	//..plot the canvases of cells into a .pdf file
	for(Int_t can=0;can<nCv;can++)
	{
		TString internal_pdfName1=pdfName+"Low.pdf";
		TString internal_pdfName2=pdfName+"High.pdf";
		TString internal_pdfName3=pdfName+"Ratio.pdf";
		if(can==0)
		{
			//..first pad
			internal_pdfName1 +="(";
			internal_pdfName2 +="(";
			internal_pdfName3 +="(";
			cComp[can]    ->Print(internal_pdfName1.Data());
			cCompDiv[can] ->Print(internal_pdfName3.Data());
		}
		else if(can==(nCv-1))
		{
			//..last pad
			internal_pdfName1 +=")";
			internal_pdfName2 +=")";
			internal_pdfName3 +=")";
			cComp[can]    ->Print(internal_pdfName1.Data());
			cCompDiv[can] ->Print(internal_pdfName3.Data());
		}
		else
		{
			//..all pads in between
			cComp[can]    ->Print(internal_pdfName1.Data());
			cCompDiv[can] ->Print(internal_pdfName3.Data());
		}
	}
	for(Int_t i=lastGood+1;i<nCv-1;i++)
	{
		delete cComp[i];
	}
	cout<<endl;

}//..plot selected cells END

//________________________________________________________________________

TString RunSpecificHistogramFileName(Int_t runnumb, TString period, TString train, TString trigger, TString workDir){

	TString analysisOutput = Form("%s/AnalysisOutput/%s", workDir.Data(), period.Data());

	TString badChannelOutput  = TString::Format("%s/Version%i/%s%s_Histograms_V%i.root", analysisOutput.Data(), runnumb, train.Data(), trigger.Data(), runnumb);

	return badChannelOutput;
}

//________________________________________________________________________
TH1D* GetAbsIDHistogram(Int_t celln, TFile *fin){

	TString listname = "BadCell_Amplitudes";
	TList *listBad = (TList*)fin->Get(listname);
	if(!listBad){
		Printf("List %s not found, return 0x0", listname.Data());
		fin->ls();
		return 0x0;
	}

	TH1D* hIDAmpl = (TH1D*)listBad->FindObject(TString::Format("Cell %d", celln));
	if(!hIDAmpl){
		Printf("Cell %d is not among the bad cells", celln);
		return 0x0;
	}

	return hIDAmpl;

}

//________________________________________________________________________
TH1D** GetAllAbsIDHistogram(TFile *fin, Int_t& totbadCells){

	if(!fin->IsOpen()){
		Printf("File not found, return 0x0");
		return 0x0;

	}
	TString listname = "BadCell_Amplitudes";
	TList *listBad = (TList*)fin->Get(listname);
	if(!listBad){
		Printf("List %s not found, return 0x0", listname.Data());
		fin->ls();
		return 0x0;
	}
	//Printf("Debugging List found");
	totbadCells = listBad->GetEntries();
	TH1D** hIDAmpl = new TH1D*[totbadCells];

	for(Int_t i = 0; i < totbadCells; i++){
		hIDAmpl[i] = (TH1D*)listBad->At(i);
		//Printf("Pointer histo %d, %p", i, hIDAmpl[i]);
	}
	Printf("End, returning");
	return hIDAmpl;

}

//________________________________________________________________________

TList* GetListAbsIDHistograms(TFile *fin){
	if(!fin->IsOpen()){
		Printf("File not found, return 0x0");
		return 0x0;

	}
	TString listname = "BadCell_Amplitudes";
	TList *listBad = (TList*)fin->Get(listname);
	if(!listBad){
		Printf("List %s not found, return 0x0", listname.Data());
		fin->ls();
		return 0x0;
	}
	return listBad;
}

///
/// compresses the bad cell histogram
/// to delete entries of good cells (only 10-15% of cells are bad)
//________________________________________________________________________
TH2F* CompressHistogram(TH2 *Histo,Int_t totalCells, Int_t badCells,Int_t nRuns)
{
	TH2F* cpmpressed = new TH2F(Form("%s_Comp",Histo->GetName()),Form("%s_Comp",Histo->GetName()),badCells+2, 0,badCells+2,nRuns,0,nRuns);

	Histo->GetXaxis()->UnZoom();
	TH1D *htmpCell = Histo->ProjectionX(TString::Format("%s_proj",Histo->GetName()),0,nRuns);
	Int_t sumRun=0,newHistoBin=0;
	for(Int_t icell = 0; icell < totalCells ; icell++)
	{
		sumRun = htmpCell->GetBinContent(icell+1);
		//cout<<"enties cell("<<icell<<"): "<<sumRun<<endl;
		//..Fill non zero entries into the new histogram
		if(sumRun>0)
		{
			newHistoBin++;
			if(newHistoBin>badCells)cout<<"PROBLEM"<<endl;
			for(Int_t iRun = 0; iRun < nRuns ; iRun++)
			{
				cpmpressed->Fill(newHistoBin,iRun,Histo->GetBinContent(icell+1,iRun+1));
				//cout<<"fill bin "<<newHistoBin<<" with value: 2"<<endl;
			}
		}
	}
	return cpmpressed;
}
///
/// Function to set TH1 histograms to a similar style
///
//________________________________________________________________________
void SetHisto(TH2 *Histo,TString Xtitel,TString Ytitel,Bool_t longhisto)
{
	Histo->SetStats(0);
	Histo->SetTitle("");
	if(longhisto==0)	Histo->GetYaxis()->SetTitleOffset(1.4);
	if(longhisto==0)	Histo->GetXaxis()->SetTitleOffset(1.4);
	if(longhisto==1)	Histo->GetYaxis()->SetTitleOffset(0.2);
	if(longhisto==1)	Histo->GetXaxis()->SetTitleOffset(1.0);
	//if(big==1)	Histo->GetYaxis()->SetLabelOffset(0.015);
	//if(big==1)	Histo->GetXaxis()->SetLabelOffset(0.015);
	if(longhisto==0)	Histo->GetXaxis()->SetLabelSize(0.05);
	if(longhisto==0)Histo->GetYaxis()->SetLabelSize(0.05);
	if(longhisto==1)	Histo->GetXaxis()->SetLabelSize(0.07);
	if(longhisto==1)Histo->GetYaxis()->SetLabelSize(0.07);
	if(longhisto==0)	Histo->GetXaxis()->SetTitleSize(0.045);
	if(longhisto==0)	Histo->GetYaxis()->SetTitleSize(0.045);
	if(longhisto==1)Histo->GetXaxis()->SetTitleSize(0.08);
	if(longhisto==1)	Histo->GetYaxis()->SetTitleSize(0.08);
	//Histo->GetXaxis()->CenterTitle();
	//Histo->GetYaxis()->CenterTitle();
	if(longhisto==1)
	{
		Histo->GetXaxis()->SetNdivisions(520);
		Histo->GetYaxis()->SetNdivisions(10);
	}
	else
	{
		Histo->GetXaxis()->SetNdivisions(505);
		Histo->GetYaxis()->SetNdivisions(505);
	}
	//..make nice font
	Histo->GetXaxis()->SetLabelFont(42);
	Histo->GetYaxis()->SetLabelFont(42);
	Histo->GetXaxis()->SetTitleFont(42);
	Histo->GetYaxis()->SetTitleFont(42);
	if(Xtitel!="")Histo->GetXaxis()->SetTitle(Xtitel);
	if(Ytitel!="")Histo->GetYaxis()->SetTitle(Ytitel);

	Histo->SetLineColor(1);
	Histo->SetMarkerColor(1);
	Histo->SetMarkerStyle(20);
	Histo->SetMarkerSize(0.5);
}
///
/// Funtion to set TH1 histograms to a similar style
///
//________________________________________________________________________
void SetHisto(TH1 *Histo,TString Xtitel,TString Ytitel,Bool_t longhisto)
{
	Histo->SetStats(0);
	Histo->SetTitle("");
	if(longhisto==0)	Histo->GetYaxis()->SetTitleOffset(1.4);
	if(longhisto==0)	Histo->GetXaxis()->SetTitleOffset(1.4);
	if(longhisto==1)	Histo->GetYaxis()->SetTitleOffset(0.2);
	if(longhisto==1)	Histo->GetXaxis()->SetTitleOffset(1.0);
	//if(big==1)	Histo->GetYaxis()->SetLabelOffset(0.015);
	//if(big==1)	Histo->GetXaxis()->SetLabelOffset(0.015);
	if(longhisto==0)	Histo->GetXaxis()->SetLabelSize(0.05);
	if(longhisto==0)Histo->GetYaxis()->SetLabelSize(0.05);
	if(longhisto==1)	Histo->GetXaxis()->SetLabelSize(0.07);
	if(longhisto==1)Histo->GetYaxis()->SetLabelSize(0.07);
	if(longhisto==0)	Histo->GetXaxis()->SetTitleSize(0.045);
	if(longhisto==0)	Histo->GetYaxis()->SetTitleSize(0.045);
	if(longhisto==1)Histo->GetXaxis()->SetTitleSize(0.08);
	if(longhisto==1)	Histo->GetYaxis()->SetTitleSize(0.08);
	//Histo->GetXaxis()->CenterTitle();
	//Histo->GetYaxis()->CenterTitle();
	if(longhisto==1)
	{
		Histo->GetXaxis()->SetNdivisions(520);
		Histo->GetYaxis()->SetNdivisions(10);
	}
	else
	{
		Histo->GetXaxis()->SetNdivisions(505);
		Histo->GetYaxis()->SetNdivisions(505);
	}	//..make nice font
	Histo->GetXaxis()->SetLabelFont(42);
	Histo->GetYaxis()->SetLabelFont(42);
	Histo->GetXaxis()->SetTitleFont(42);
	Histo->GetYaxis()->SetTitleFont(42);
	if(Xtitel!="")Histo->GetXaxis()->SetTitle(Xtitel);
	if(Ytitel!="")Histo->GetYaxis()->SetTitle(Ytitel);

	Histo->SetLineColor(1);
	Histo->SetMarkerColor(1);
	Histo->SetMarkerStyle(20);
	Histo->SetMarkerSize(0.5);
}
///
/// Calculate how to best split period into run blocks
/// So that the number of masked cells is minimized
//________________________________________________________________________
void GetBestPeriodSplitting(TString period = "LHC15o", TString train = "Train_641",Int_t noOfSplits=4)
{
	if(noOfSplits>4)cout<<"Error: so far only implemented for 1-4 splits"<<endl;

	gStyle->SetOptTitle(0);
	gStyle->SetOptStat(0);
	//gStyle->SetPalette(53);  //standard is 1
	gStyle->SetCanvasColor(10);
	TGaxis::SetMaxDigits(4);
	gStyle->SetPadTopMargin(0.07);//0.05
	gStyle->SetPadBottomMargin(0.18);//0.15
	gStyle->SetPadRightMargin(0.04);
	gStyle->SetPadLeftMargin(0.06);
	gStyle->SetFrameFillColor(10);
	gStyle->SetLabelSize(0.05,"X");
	gStyle->SetLabelSize(0.05,"Y");
	gStyle->SetTitleSize(5.0,"X");
	gStyle->SetTitleSize(5.0,"Y");
	gEnv->SetValue("Canvas.ShowEventStatus",1);  //shows the status bar in the canvas

	//Int_t Nruns=99;
	Int_t Nruns=45;
	Int_t runNumber=245145;

	AliCalorimeterUtils *fCaloUtils = new AliCalorimeterUtils();
	//..Create a dummy event for the CaloUtils
	AliAODEvent* aod = new AliAODEvent();
	fCaloUtils->SetRunNumber(runNumber);
	fCaloUtils->AccessGeometry(aod);
    Int_t noOfCells=fCaloUtils->GetEMCALGeometry()->GetNCells(); //..Very important number, never change after that point!

	TString fileName = Form("AnalysisOutput/%s/%s/RunByRunSummary/runList%i_Results.root",period.Data(),train.Data(),Nruns);

	cout<<"Open root file: "<<fileName<<endl;
	TFile *f = TFile::Open(fileName);
	if(!f)
	{
		cout<<"Couldn't open/find .root file: "<<fileName<<endl;
	}
	//..Get the histogram
	TH2F *hbadAndDeadvsRun = (TH2F*)f->Get("hFlag3vsRun_Comp");
	TH2F *hbadAndDeadvsRunw = (TH2F*)hbadAndDeadvsRun->Clone("hbadAndDeadvsRunWeight");
	if(!hbadAndDeadvsRun)cout<<"couldn't find histogram - hFlag3vsRun_Comp"<<endl;
	hbadAndDeadvsRun->GetXaxis()->UnZoom();
	TH1D* hEventsPerRun = (TH1D*)f->Get("nEventVsRun");
	if(!hEventsPerRun)cout<<"couldn't find histogram - nEventVsRun"<<endl;

	//..weight bad cells by numbers of events
	for(Int_t icell = 0; icell < noOfCells ; icell++)
	{
		for(Int_t iRun = 0; iRun < Nruns ; iRun++)
		{
			if(hbadAndDeadvsRun->GetBinContent(icell+1,iRun+1)==1)
			{
				hbadAndDeadvsRunw->SetBinContent(icell+1,iRun+1,hEventsPerRun->GetBinContent(iRun+1));
			}
		}
	}
	//..Draw the histogram
	TCanvas *canEvt= new TCanvas("canEvt", "canEvt", 500, 500);
	canEvt->cd(1);
	SetHisto(hEventsPerRun,"No of events","Run",1);
	//hEventsPerRun->GetYaxis()->SetTitleOffset(0.35);
	hEventsPerRun->DrawCopy("hist"); //box


	TCanvas *can= new TCanvas("compressedIDs", "compressed cell ID's", 1600, 500);
	can->Divide(0,2);
	can->cd(1);
	SetHisto(hbadAndDeadvsRun,"","compressed ID",1);
	hbadAndDeadvsRun->GetYaxis()->SetTitleOffset(0.35);
	hbadAndDeadvsRun->DrawCopy("colz"); //box
	can->cd(2);
	SetHisto(hbadAndDeadvsRunw,"","compressed ID",1);
	hbadAndDeadvsRunw->GetYaxis()->SetTitleOffset(0.35);
	hbadAndDeadvsRunw->DrawCopy("colz");

	//..Analyze the histogram
	//..split into three blocks and see what performs better
	Int_t splitRun1=0;
	Int_t splitRun2=0;
	Int_t splitRun3=0;
	Int_t endBlock3;
	Double_t totalCellsBadRun=0;
	Double_t totalCellsBadEvt=0;
	Double_t nCellRunsBlock1=0;
	Double_t nCellRunsBlock2=0;
	Double_t nCellRunsBlock3=0;
	Double_t nCellRunsBlock4=0;
	Double_t nCellEvtBlock1=0;
	Double_t nCellEvtBlock2=0;
	Double_t nCellEvtBlock3=0;
	Double_t nCellEvtBlock4=0;

	for(Int_t iRun=1; iRun<=Nruns; iRun++)
	{
		cout<<"Round "<<iRun<<" of "<<Nruns<<endl;
		for(Int_t iRun2=iRun+1; iRun2<=Nruns; iRun2++)
		{
			if(noOfSplits<4)
			{
				endBlock3=Nruns-1;
			}
			else
			{
				endBlock3=iRun2+1;
			}

			for(Int_t iRun3=endBlock3; iRun3<=Nruns; iRun3++)
			{
				hbadAndDeadvsRun->GetXaxis()->UnZoom();
				TH1D *htmpCell1 = hbadAndDeadvsRun->ProjectionX(TString::Format("%s_proj1",hbadAndDeadvsRun->GetName()),0,iRun);
				TH1D *htmpCell2 = hbadAndDeadvsRun->ProjectionX(TString::Format("%s_proj2",hbadAndDeadvsRun->GetName()),iRun+1,iRun2);
				TH1D *htmpCell3 = hbadAndDeadvsRun->ProjectionX(TString::Format("%s_proj3",hbadAndDeadvsRun->GetName()),iRun2+1,iRun3);
				TH1D *htmpCell4 = hbadAndDeadvsRun->ProjectionX(TString::Format("%s_proj4",hbadAndDeadvsRun->GetName()),iRun3+1,Nruns);

				Double_t nEvtBlock1 = hEventsPerRun->Integral(0,iRun);
				Double_t nEvtBlock2 = hEventsPerRun->Integral(iRun+1,iRun2);
				Double_t nEvtBlock3 = hEventsPerRun->Integral(iRun2+1,iRun3);
				Double_t nEvtBlock4 = hEventsPerRun->Integral(iRun3+1,Nruns);

				Int_t sumRunBlock1=0;
				Int_t sumRunBlock2=0;
				Int_t sumRunBlock3=0;
				Int_t sumRunBlock4=0;
				Int_t nBlock1=0;
				Int_t nBlock2=0;
				Int_t nBlock3=0;
				Int_t nBlock4=0;

				for(Int_t icell = 0; icell < htmpCell1->GetNbinsX(); icell++)
				{
					sumRunBlock1 = htmpCell1->GetBinContent(icell+1);
					sumRunBlock2 = htmpCell2->GetBinContent(icell+1);
					sumRunBlock3 = htmpCell3->GetBinContent(icell+1);
					sumRunBlock4 = htmpCell4->GetBinContent(icell+1);
					if(sumRunBlock1>0)                nBlock1++;
					if(sumRunBlock2>0 && noOfSplits>1)nBlock2++;
					if(sumRunBlock3>0 && noOfSplits>2)nBlock3++;
					if(sumRunBlock4>0 && noOfSplits>3)nBlock4++;
				}
				nCellRunsBlock1=nBlock1*iRun;
				nCellRunsBlock2=nBlock2*(iRun2-iRun+1);
				nCellRunsBlock3=nBlock3*(iRun3-iRun2+1);
				nCellRunsBlock4=nBlock4*(Nruns-iRun3+1);

				/*cout<<". . . . . . . . . . . . . . . . . . ."<<endl;
				cout<<" events block1 : "<<nEvtBlock1<<endl;
				cout<<" events block2 : "<<nEvtBlock2<<endl;
				cout<<" events block3 : "<<nEvtBlock3<<endl;
				cout<<" events block4 : "<<nEvtBlock4<<endl;
				cout<<" cells block1 : "<<nBlock1<<endl;
				cout<<" cells block2 : "<<nBlock2<<endl;
				cout<<" cells block3 : "<<nBlock3<<endl;
				cout<<" cells block4 : "<<nBlock4<<endl;
*/
				nCellEvtBlock1 =nBlock1*nEvtBlock1;
				nCellEvtBlock2 =nBlock2*nEvtBlock2;
				nCellEvtBlock3 =nBlock3*nEvtBlock3;
				nCellEvtBlock4 =nBlock4*nEvtBlock4;
/*
				cout<<" cells*evt block1 : "<<nCellEvtBlock1<<endl;
				cout<<" cells*evt block2 : "<<nCellEvtBlock2<<endl;
				cout<<" cells*evt block3 : "<<nCellEvtBlock3<<endl;
				cout<<" cells*evt block4 : "<<nCellEvtBlock4<<endl;
*/
				//cout<<"split int: 0-"<<iRun<<", and "<<iRun+1<<"-"<<iRun2<<" and "<<iRun2+1<<"-"<<Nruns<<endl;
				//cout<<"bad cells total :"<<nBlock1+nBlock2+nBlock3+nBlock4<<endl;

				//..not weighted by nuber of events in run
				if(totalCellsBadRun==0 || (nCellRunsBlock1+nCellRunsBlock2+nCellRunsBlock3+nCellRunsBlock4)<totalCellsBadRun)
				{
					totalCellsBadRun=nCellRunsBlock1+nCellRunsBlock2+nCellRunsBlock3+nCellRunsBlock4;
					/*splitRun1=iRun;
					splitRun2=iRun2;
					splitRun3=iRun3;
					cout<<"update"<<endl;
					cout<<"totalCellsBadRun: "<<totalCellsBadRun<<", nCellRunsBlock1: "<<nCellRunsBlock1<<", nCellRunsBlock2: "<<nCellRunsBlock2<<", nCellRunsBlock3: "<<nCellRunsBlock3<<", nCellRunsBlock4: "<<nCellRunsBlock4<<endl;
					*/
				}
				//..weighted by nuber of events in run
				if(totalCellsBadEvt==0 || (nCellEvtBlock1+nCellEvtBlock2+nCellEvtBlock3+nCellEvtBlock4)<totalCellsBadEvt)
				{
					totalCellsBadEvt=nCellEvtBlock1+nCellEvtBlock2+nCellEvtBlock3+nCellEvtBlock4;
					splitRun1=iRun;
					splitRun2=iRun2;
					splitRun3=iRun3;
					//cout<<"update"<<endl;
					//cout<<"split int: 0-"<<iRun<<", and "<<iRun+1<<"-"<<iRun2<<" and "<<iRun2+1<<"-"<<iRun3<<" and "<<iRun3+1<<"-"<<Nruns<<endl;
					//cout<<"totalCellsBadRun: "<<totalCellsBadEvt<<", nCellEvtBlock1: "<<nCellEvtBlock1<<", nCellEvtBlock2: "<<nCellEvtBlock2<<", nCellEvtBlock3: "<<nCellEvtBlock3<<", nCellEvtBlock4: "<<nCellEvtBlock4<<endl;
				}
			}
		}
	}
	hbadAndDeadvsRun->GetXaxis()->UnZoom();
	TH1D *htmpCell1p = hbadAndDeadvsRun->ProjectionX(TString::Format("%s_proj1",hbadAndDeadvsRun->GetName()),0,splitRun1);
	TH1D *htmpCell2p = hbadAndDeadvsRun->ProjectionX(TString::Format("%s_proj2",hbadAndDeadvsRun->GetName()),splitRun1+1,splitRun2);
	TH1D *htmpCell3p = hbadAndDeadvsRun->ProjectionX(TString::Format("%s_proj3",hbadAndDeadvsRun->GetName()),splitRun2+1,splitRun3);
	TH1D *htmpCell4p = hbadAndDeadvsRun->ProjectionX(TString::Format("%s_proj4",hbadAndDeadvsRun->GetName()),splitRun3+1,Nruns);

	TCanvas *canSplit= new TCanvas("canSplit", "Split compressed cell ID's", 1600, 500);
	canSplit->Divide(2,2);
	canSplit->cd(1);
	SetHisto(htmpCell1p,"","nruns bad",1);
	htmpCell1p->GetYaxis()->SetTitleOffset(0.35);
	htmpCell1p->DrawCopy("hist");
	canSplit->cd(2);
	SetHisto(htmpCell2p,"","nruns bad",1);
	htmpCell2p->GetYaxis()->SetTitleOffset(0.35);
	htmpCell2p->DrawCopy("hist");
	canSplit->cd(3);
	SetHisto(htmpCell3p,"","nruns bad",1);
	htmpCell3p->GetYaxis()->SetTitleOffset(0.35);
	htmpCell3p->DrawCopy("hist");
	canSplit->cd(4);
	SetHisto(htmpCell4p,"","nruns bad",1);
	htmpCell4p->GetYaxis()->SetTitleOffset(0.35);
	htmpCell4p->DrawCopy("hist");

	cout<<"Best results are achieved by splitting into:"<<endl;
	cout<<"Run: 0-"<<splitRun1<<endl;
	cout<<"Run: "<<splitRun1+1<<"-"<<splitRun2<<endl;
	cout<<"Run: "<<splitRun2+1<<"-"<<splitRun3<<endl;
	cout<<"Run: "<<splitRun3+1<<"-"<<Nruns<<endl;
	cout<<"Number of bad cells*events ="<<totalCellsBadEvt<<endl;
	cout<<"Number of bad cells*runs   ="<<totalCellsBadRun<<endl;

	cout<<" events block1 : "<<hEventsPerRun->Integral(0,splitRun1)<<endl;
	cout<<" events block2 : "<<hEventsPerRun->Integral(splitRun1+1,splitRun2)<<endl;
	cout<<" events block3 : "<<hEventsPerRun->Integral(splitRun2+1,splitRun3)<<endl;
	cout<<" events block4 : "<<hEventsPerRun->Integral(splitRun3+1,Nruns)<<endl;
}
// is this a High intensity run?
//_____________________________________________________________________
Bool_t isHIRun(Int_t runID)
{
	Bool_t hiRun=0;
	Int_t hiRunList[57] = {245145,245146,245151,245152,245231,245232,245683,245700,245702,245705,
					  	245829,245831,245833,245949,245952,245954,245963,246001,246003,246037,
					  	246042,246052,246053,246087,246089,246113,246115,246217,246222,246225,
					  	246271,246272,246424,246434,246487,246488,246493,246495,246750,246751,
					  	246757,246758,246759,246760,246765,246766,246804,246805,246807,246808,
					  	246809,246810,246844,246845,246846,246928,246945};
	for(Int_t i=0;i<57;i++)
	{
		if(hiRunList[i]==runID)
		{
			hiRun=1;
			break;
		}
	}
	return hiRun;
}
