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
// root [2] Run_BadChannel(1,"LHC15n","Train_603","AnyINTnoBC",244411)
//  
/////////////////////////////////////////////////

#include <TString.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TH2D.h>
#include <TH1D.h>
#include <TList.h>
#include <TSystem.h>
#include <TStopwatch.h>

#include "AliAnaCaloChannelAnalysis.h" //include when compile
#include "AliEMCALGeometry.h"          //include when compile
#include "AliCalorimeterUtils.h"       //include when compile
#include "AliAODEvent.h"               //include when compile
#include "AliOADBContainer.h"          //include when compile

//colors
const Int_t colors[]     = {kBlack, kRed+1 , kBlue+1, kGreen+3, kMagenta+1, kOrange-1,kCyan+2,kYellow+2, kRed-9, kGreen-10, kBlue - 8};

//definition of methods
TH1D* GetAbsIDHistogram(Int_t celln, TFile *fin);
TH1D** GetAllAbsIDHistogram(TFile *fin, Int_t& totbadCells);
TList* GetListAbsIDHistograms(TFile *fin);
TString RunSpecificHistogramFileName(Int_t runnumb, TString period, TString train = "Train_641", TString trigger= "AnyINTnoBC", TString workDir=".");

//________________________________________________________________________
void Run_BadChannel(Int_t nversion = 1, TString period = "LHC15n", TString train = "Train_603", TString trigger= "AnyINTnoBC", Int_t runNum= 245683, TString externalFile= "",TString listName="runList.txt",TString workDir=".")
{
	TStopwatch watch;
	watch.Start();
	
	AliAnaCaloChannelAnalysis* Analysis;
	//..If you do the analysis run by run - this might be helpful
	Analysis=new AliAnaCaloChannelAnalysis(period,train,trigger,runNum,runNum,workDir,listName);

	//..If you do it with merged files use this
	//Analysis=new AliAnaCaloChannelAnalysis(period,train,trigger,runNum,nversion,workDir,listName);

	//..Settings
	Analysis->SetExternalMergedFile(externalFile);
	//Analysis->SetQAChecks(1);  //1=Perform QA checks - takes a long time! Prints all good cells for cross check

	//. . . . . . . . . . . . . . . . . . . . . . . .
	//. . Add different period analyses
	//. . . . . . . . . . . . . . . . . . . . . . . .
/*	Analysis->AddPeriodAnalysis(2, 4.,0.1,0.3); // hit/event in range Emin Emax
	Analysis->AddPeriodAnalysis(1, 6.,0.1,0.3); // energy/hit in range Emin Emax
	Analysis->AddPeriodAnalysis(2, 4.,0.2,0.5); // hit/event in range Emin Emax
	Analysis->AddPeriodAnalysis(1, 6.,0.2,0.5); // energy/hit in range Emin Emax
	Analysis->AddPeriodAnalysis(2, 4.,0.5,1.0); // hit/event in range Emin Emax
	Analysis->AddPeriodAnalysis(1, 6.,0.5,1.0); // energy/hit in range Emin Emax

	//..If there is enough statistic add also these:
	Analysis->AddPeriodAnalysis(2, 4.,1.0,4.0); // hit/event in range Emin Emax
	Analysis->AddPeriodAnalysis(1, 6.,1.0,4.0); // mean energy in range Emin Emax
	Analysis->AddPeriodAnalysis(2, 4.,1.0,10.0);// hit/event in range Emin Emax
	Analysis->AddPeriodAnalysis(1, 5.,1.0,10.0);// energy/hit in range Emin Emax
//	Analysis->AddPeriodAnalysis(2, 5.,3.0,10.0);// (PbPb extra range) hit/event in range Emin Emax
//	Analysis->AddPeriodAnalysis(1, 5.,3.0,10.0);// (PbPb extra range) energy/hit in range Emin Emax
*/

	///*test time stuff*/	Analysis->AddPeriodAnalysis(3, 6,-20,+20);// energy/hit in range Emin Emax

	Analysis->AddPeriodAnalysis(2, 5.5,0.1,0.3);  // hit/event in range Emin Emax
	Analysis->AddPeriodAnalysis(1, 4.5,0.1,0.3);  // energy/hit in range Emin Emax
	Analysis->AddPeriodAnalysis(2, 5.5,0.2,0.5);  // hit/event in range Emin Emax
	Analysis->AddPeriodAnalysis(1, 4.5,0.2,0.5);  // energy/hit in range Emin Emax
	Analysis->AddPeriodAnalysis(2, 4.0,0.5,1.0);  // hit/event in range Emin Emax
	Analysis->AddPeriodAnalysis(1, 4.5,0.5,1.0);  // energy/hit in range Emin Emax
	Analysis->AddPeriodAnalysis(2, 4.0,1.0,4.0);  // hit/event in range Emin Emax
	Analysis->AddPeriodAnalysis(1, 4.5,1.0,4.0);  // mean energy in range Emin Emax
	Analysis->AddPeriodAnalysis(2, 4.0,1.0,10.0); // hit/event in range Emin Emax
	Analysis->AddPeriodAnalysis(1, 4.5,1.0,10.0); // energy/hit in range Emin Emax
//	Analysis->AddPeriodAnalysis(2, 5.,3.0,10.0);// (PbPb extra range) hit/event in range Emin Emax
//	Analysis->AddPeriodAnalysis(1, 5.,3.0,10.0);// (PbPb extra range) energy/hit in range Emin Emax


	//..Start the bad channel analysis
	Analysis->Run();
	
	watch.Stop();
	watch.Print();
}
//
//check where the bad cell is (row collumn), if you have only its ID
//
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
	TString path        = Form("/Users/Eliane/Software/BadChannelAnalysis/AnalysisOutput/%s/Version%d",period.Data(),version);
	TString rootFileName= Form("Train_%dAnyINTnoBC_Histograms_V%d.root",trainNo,version);
	TFile* outputRoot   = TFile::Open(Form("%s/%s",path.Data(),rootFileName.Data()));

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
	TObjArray *recal=(TObjArray*)cont->GetObject(runNumber);
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
//________________________________________________________________________
/// Draw the good, bad, dead channel maps, and the amplitude distribution per each run and save them in pdf files in analysisOutput/train/RunByRunSummary
/// -- can improve: write different pages in one pdf
/// -- decide how to treat this info
void SummarizeRunByRun(TString period = "LHC15o", TString train = "Train_641", TString trigger= "AnyINTnoBC", TString workDir=".", TString listName="runList.txt")
{
	//..Open the text file with the run list numbers and run index
	cout<<"o o o Open .txt file with run indices. Name = " << listName << endl;
	TString analysisInput  = Form("AnalysisInput/%s",period.Data());
	TString analysisOutput = Form("AnalysisOutput/%s",period.Data());
	TString runList        = Form("%s/%s/%s/%s",workDir.Data(), analysisInput.Data(), train.Data(),listName.Data());

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


	//..Open the different .root files with help of the run numbers from the text file
	const Int_t nRun = nlines ;
	TString rootFileName;
	TString badChannelOutput;
	TString infoText;

	Int_t totalperCv = 25;
	Int_t nPad = TMath::Sqrt(totalperCv);
	Int_t nCv = nRun/totalperCv;
	
	TCanvas **cBad  = new TCanvas*[nCv];
	TCanvas **cGood = new TCanvas*[nCv];
	TCanvas **cDead = new TCanvas*[nCv];
	TCanvas **cAmp  = new TCanvas*[nCv];
	for(Int_t ic = 0; ic<nCv; ic++){
		cBad [ic] = new TCanvas(TString::Format("badcells%d", ic), TString::Format("badcells  (%d/%d)", ic+1, nCv), 1000,750);
		cGood[ic] = new TCanvas(TString::Format("goodcells%d", ic),TString::Format("goodcells (%d/%d)", ic+1, nCv),1000,750);
		cDead[ic] = new TCanvas(TString::Format("deadcells%d", ic),TString::Format("deadcells (%d/%d)", ic+1, nCv),1000,750);
		cAmp [ic] = new TCanvas(TString::Format("Amplitide%d", ic),TString::Format("Amplitide (%d/%d)", ic+1, nCv),1000,750);
		
		cBad [ic] ->Divide(nPad,nPad);
		cGood[ic] ->Divide(nPad,nPad);
		cDead[ic] ->Divide(nPad,nPad);
		cAmp [ic] ->Divide(nPad,nPad);
	}
	
	Int_t nFlags = 3;
	TH2F** hFlagvsRun = new TH2F*[nFlags]; 
	TH2F** hFlagNew = new TH2F*[nFlags];
	hFlagNew[0] = 0x0;
	hFlagNew[1] = 0x0;
	hFlagNew[2] = 0x0;
	TH1F** hFlagPerID = new TH1F*[nFlags];
	
	TCanvas **cFlag = new TCanvas*[nFlags];
	TCanvas **cFlagNewID = new TCanvas*[nFlags];
	// init
	for(Int_t i = 0; i<nFlags; i++){
		hFlagvsRun[i] = 0x0;
		hFlagNew  [i] = 0x0;
		hFlagPerID[i] = 0x0;
		
		cFlag[i] = new TCanvas(TString::Format("cFlag%d", i), TString::Format("cFlag %d ", i), 1200, 500);
		
		cFlagNewID[i] = new TCanvas(TString::Format("cFlagNewID%d", i), "cFlag new vs cell ID", 1200, 800);
	}
	
	TCanvas *cFlagNew = new TCanvas("cFlagNew", "cFlag new", 800, 800);
	cFlagNew->Divide(2, 2);
	
	
	AliCalorimeterUtils *fCaloUtils = new AliCalorimeterUtils();
	//..Create a dummy event for the CaloUtils
	AliAODEvent* aod = new AliAODEvent();
	fCaloUtils->SetRunNumber(RunId[0]);
	fCaloUtils->AccessGeometry(aod);
	
	Int_t    ncells  = 0, runNotFound = 0;
	//..loop over the amount of run numbers found in the previous text file.
	for(Int_t i = 0 ; i < nRun ; i++)  //Version%i" LHC16i_muon_caloLego_Histograms_V255539
	{
		rootFileName      = Form("%s%s_Histograms_V%i.root", train.Data(), trigger.Data(), RunId[i]);
		badChannelOutput  = Form("%s/Version%i/%s", analysisOutput.Data(), RunId[i],rootFileName.Data());
		
		cout<<"Open root file: "<<badChannelOutput<<endl;
		TFile *f = TFile::Open(badChannelOutput);
		if(!f)
		{
			cout<<"Couldn't open/find .root file: "<<badChannelOutput<<endl;
			RunId[i] = -1;
			runNotFound++;
			continue;
		}
		
		TH2F *badCells  = (TH2F*)f->Get("2DChannelMap_Flag2");
		if(!badCells) {
			Printf("2DChannelMap_Flag2 not found");
			continue;
		}
		
		
		TH2F *goodCells = (TH2F*)f->Get("2DChannelMap_Flag0");
		if(!goodCells) {
			Printf("2DChannelMap_Flag0 not found");
			continue;
		}
		
		
		TH2F *deadCells = (TH2F*)f->Get("2DChannelMap_Flag1");
		if(!deadCells) {
			Printf("2DChannelMap_Flag1 not found");
			continue;
		}
		
		TH2F *ampID     = (TH2F*)f->Get("hCellAmplitude");
		if(!ampID) {
			Printf("hCellAmplitude not found");
			continue;
		}
		
		TH1F *hCellFlag = (TH1F*)f->Get("fhCellFlag");
		if(!hCellFlag) {
			Printf("fhCellFlag not found");
			continue;
		}
		ncells  = hCellFlag->GetNbinsX();
		if(!hFlagvsRun[0]) {
			Printf("Defining hFlagvsRun");
			
			Double_t mincell = hCellFlag->GetXaxis()->GetBinLowEdge(1);
			Double_t maxcell = hCellFlag->GetXaxis()->GetBinLowEdge(ncells+1);
			
			hFlagvsRun[0] = new TH2F("hFlag1vsRun", "hFlag1vsRun (?); cell ID; Run number", ncells, mincell, maxcell, nRun-runNotFound, 0, nRun-runNotFound-1); // update this axis, need to have the run number
			hFlagvsRun[1] = new TH2F("hFlag2vsRun", "hFlag2vsRun (?); cell ID; Run number", ncells, mincell, maxcell, nRun-runNotFound, 0, nRun-runNotFound-1);
			hFlagvsRun[2] = new TH2F("hFlag3vsRun", "hFlag3vsRun (?); cell ID; Run number", ncells, mincell, maxcell, nRun-runNotFound, 0, nRun-runNotFound-1);

			
		}
		for(Int_t ic = 0; ic < ncells; ic++){
			Int_t flag = hCellFlag->GetBinContent(ic+1);
			
			if(flag == 1) hFlagvsRun[0]->Fill(ic, i, 1);
			if(flag>1){
				hFlagvsRun[1]->Fill(ic, i, 1); //fill, use the x, y values
			}
			if(flag>0){
				hFlagvsRun[2]->Fill(ic, i, 1);
			}
		}
		
		if(!hFlagNew[0]){
			hFlagNew[0] = (TH2F*)goodCells->Clone(TString::Format("2DChannelMapNew_Flag1"));
			hFlagNew[0]->Reset();
			hFlagNew[0]->SetTitle("Selected flag 1; cell column (#eta direction); cell raw (#phi direction)");
			hFlagNew[1] = (TH2F*)goodCells->Clone(TString::Format("2DChannelMapNew_Flagg1"));
			hFlagNew[1]->Reset();
			hFlagNew[1]->SetTitle("Selected flag greater than 1; cell column (#eta direction); cell raw (#phi direction)");
			hFlagNew[2] = (TH2F*)goodCells->Clone(TString::Format("2DChannelMapNew_FlagAll"));
			hFlagNew[2]->Reset();
			hFlagNew[2]->SetTitle("Selected flag greater than 0; cell column (#eta direction); cell raw (#phi direction)");
			
			hFlagPerID[0] = (TH1F*)hCellFlag->Clone(TString::Format("FlagCellID_Flag1"));
			hFlagPerID[0]->SetTitle("Selected flag 1; Abs. Cell ID; Fraction of Runs");
			hFlagPerID[0]->Reset();
			hFlagPerID[1] = (TH1F*)hCellFlag->Clone(TString::Format("FlagCellID_Flagg1"));
			hFlagPerID[1]->SetTitle("Selected flag greater than 1; Abs. Cell ID; Fraction of Runs");
			hFlagPerID[1]->Reset();
			hFlagPerID[2] = (TH1F*)hCellFlag->Clone(TString::Format("FlagCellID_FlagAll"));
			hFlagPerID[2]->SetTitle("Selected flag greater than 0; Abs. Cell ID; Fraction of Runs");
			hFlagPerID[2]->Reset();
		}
		
		// Drawing
		//....................................
		cBad[i/totalperCv]->cd(i/totalperCv + i%totalperCv+1);
		badCells->Draw("colz");
		infoText=Form("Bad Cells - Run %i",RunId[i]);
		TLatex* text = new TLatex(0.2,0.8,infoText);
		text->SetTextSize(0.06);
		text->SetNDC();
		text->SetTextColor(1);
		text->Draw();
		//....................................
		cGood[i/totalperCv]->cd(i/totalperCv + i%totalperCv+1);
		goodCells->Draw("colz");
		infoText=Form("Good Cells - Run %i",RunId[i]);
		TLatex* text1 = new TLatex(0.2,0.8,infoText);
		text1->SetTextSize(0.06);
		text1->SetNDC();
		text1->SetTextColor(1);
		text1->Draw();
		//....................................
		cDead[i/totalperCv]->cd(i/totalperCv + i%totalperCv+1);
		deadCells->Draw("colz");
		infoText=Form("Dead Cells - Run %i",RunId[i]);
		TLatex* text2 = new TLatex(0.2,0.8,infoText);
		text2->SetTextSize(0.06);
		text2->SetNDC();
		text2->SetTextColor(1);
		text2->Draw();
		//....................................
		cAmp[i/totalperCv]->cd(i/totalperCv + i%totalperCv+1)->SetLogz();
		//			ampID->GetYaxis()->SetRangeUser(7500,9500); //LHC16i
		//			ampID->GetYaxis()->SetRangeUser(1400,1600); //LHC16i
		ampID->GetYaxis()->SetRangeUser(4750,5050); //LHC16h
		ampID->Draw("colz");
		infoText=Form("Amplitudes - Run %i",RunId[i]);
		TLatex* text3 = new TLatex(0.2,0.8,infoText);
		text3->SetTextSize(0.06);
		text3->SetNDC();
		text3->SetTextColor(1);
		text3->Draw();
		
	}
	
	
	
	//cFlag->cd(2);
	//hFlagvsRun[1]->Draw("colz");
	//cFlag->cd(3);
	//hFlagvsRun[2]->Draw("colz");
	
	// define cut for bad channel declaration, not used for the moment
	Double_t percbad = 0.9;

	Printf("Nbins = %d, %d, total = %d", hFlagvsRun[0]->GetNbinsX(), hFlagvsRun[0]->GetNbinsY(), hFlagvsRun[0]->GetBin(hFlagvsRun[0]->GetNbinsX(), hFlagvsRun[0]->GetNbinsY()));
	for(Int_t iflag = 0; iflag < nFlags; iflag++){//
		Printf("%p", cFlag[iflag]);
		cFlag[iflag]->cd();
		hFlagvsRun[iflag]->Draw("colz");
		for(Int_t ic = 0; ic < ncells; ic++){//
			TH1D *htmpCell = hFlagvsRun[iflag]->ProjectionY(TString::Format("hIDProj_cell%d", ic), ic+1, ic+1);
			Double_t fracRun = 0, sumRun = 0;
			//if(htmpCell->Integral() > 0) Printf("Integral cell %d %f", ic,  htmpCell->Integral());
			for(Int_t ir = 0; ir < nRun ; ir++){
				sumRun += htmpCell->GetBinContent(ir);
			}
			fracRun = sumRun/(Double_t)(nRun-runNotFound);
			//loose selection criteria to remove the runs with zero entries
			if (fracRun>0) {
				//Printf("Frac = %f, %f", sumRun, fracRun);
			
			//if(fracRun > percbad) {
				Int_t cellColumn=0, cellRow=0;
				Int_t cellColumnAbs=0, cellRowAbs=0;
				Int_t trash = 0 ;
				fCaloUtils->GetModuleNumberCellIndexesAbsCaloMap(ic,0,cellColumn,cellRow,trash,cellColumnAbs,cellRowAbs);
				//Printf("Cell %d -> %d, %d     %d, %d", ic, cellColumn, cellRow, cellColumnAbs, cellRowAbs);
				hFlagNew[iflag]->Fill(cellColumnAbs, cellRowAbs, fracRun);
				
				hFlagPerID[iflag]->SetBinContent(ic, fracRun);
			}
		}
		cFlagNewID[iflag]->cd();
		hFlagPerID[iflag]->Draw();
	}
	cFlagNew->cd(1);
	hFlagNew[0]->Draw("colz");
	cFlagNew->cd(2);
	hFlagNew[1]->Draw("colz");
	cFlagNew->cd(3);
	hFlagNew[2]->Draw("colz");
	
	
	//cFlagNewID->cd(2);  
	//hFlagPerID[1]->Draw();
	//cFlagNewID->cd(3);  
	//hFlagPerID[2]->Draw();
	
	gSystem->mkdir(TString::Format("%s/%s/", analysisOutput.Data(), train.Data()));
	gSystem->mkdir(TString::Format("%s/%s/RunByRunSummary/", analysisOutput.Data(), train.Data()));
	for(Int_t ic = 0; ic<nCv; ic++){
		cBad [ic] ->SaveAs(TString::Format("%s/%s/RunByRunSummary/%s.pdf", analysisOutput.Data(), train.Data(), cBad [ic]->GetName()));
		cGood[ic] ->SaveAs(TString::Format("%s/%s/RunByRunSummary/%s.pdf", analysisOutput.Data(), train.Data(), cGood[ic]->GetName()));
		cDead[ic] ->SaveAs(TString::Format("%s/%s/RunByRunSummary/%s.pdf", analysisOutput.Data(), train.Data(), cDead[ic]->GetName()));
		cAmp [ic] ->SaveAs(TString::Format("%s/%s/RunByRunSummary/%s.pdf", analysisOutput.Data(), train.Data(), cAmp [ic]->GetName()));
	}
	for(Int_t iflag = 0; iflag < nFlags; iflag++){
		cFlag[iflag]->SaveAs(TString::Format("%s/%s/RunByRunSummary/%s.pdf", analysisOutput.Data(), train.Data(), cFlag[iflag]->GetName()));
		cFlagNewID[iflag]->SaveAs(TString::Format("%s/%s/RunByRunSummary/%s.pdf", analysisOutput.Data(), train.Data(), cFlagNewID[iflag]->GetName()));
	}
	cFlagNew->SaveAs(TString::Format("%s/%s/RunByRunSummary/%s.pdf", analysisOutput.Data(), train.Data(), cFlagNew->GetName()));
	
}

//________________________________________________________________________

void PlotSelectionAbsID(TString listName, TString period = "LHC15o", TString train = "Train_641", TString trigger= "AnyINTnoBC", TString workDir="."){
	
	
	cout<<"o o o Open .txt file with run indices. Name = " << listName << endl;
	TString analysisInput  = Form("AnalysisInput/%s",period.Data());
	TString analysisOutput = Form("AnalysisOutput/%s",period.Data());
	TString runList        = Form("%s/%s/%s/%s",workDir.Data(), analysisInput.Data(), train.Data(),listName.Data());

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
	
	//..Initialize EMCal/DCal geometry
	AliCalorimeterUtils* fCaloUtils = new AliCalorimeterUtils();
	//..Create a dummy event for the CaloUtils
	AliAODEvent* aod = new AliAODEvent();
	
	fCaloUtils->AccessGeometry(aod);
	
	
	//prepare Canvas
	Int_t bigNumb = 180, nPad = 9, ncvfilled = 0;
	TCanvas **cAmpl = new TCanvas*[bigNumb];
	for(Int_t i = 0; i< bigNumb; i++){
		cAmpl[i]= new TCanvas(TString::Format("cbadcells%d", i), TString::Format("badcells, cv %d", i),1000,750);
		cAmpl[i]->Divide(TMath::Sqrt(nPad),TMath::Sqrt(nPad));
	}
	std::map<Int_t, Int_t> mapChPadIdx;
	Int_t counter = 0, nc = 0, np = 0;
	for(Int_t ir = 0 ; ir < nlines; ir++){ //
		Printf("--------------- Run %d", RunId[ir]);
		TFile *fin = new TFile(RunSpecificHistogramFileName(RunId[ir], period, train, trigger, workDir));
		if(!fin->IsOpen())continue;
		// get the histograms
		Int_t nh = 0;
		TList *listAmp = GetListAbsIDHistograms(fin);
		
		fCaloUtils->SetRunNumber(RunId[ir]);
		Int_t ncells = fCaloUtils->GetEMCALGeometry()->GetNCells();
		
		for(Int_t ic = 0; ic < ncells; ic++){ //
			
			TH1D* hTmpBad = (TH1D*)listAmp->FindObject(TString::Format("Cell %d", ic));
			if(!hTmpBad) continue;
			
			if(mapChPadIdx.count(ic) == 0){ // if this map position was never initialized
				mapChPadIdx[ic] = counter; //map to the position in the canvas where this cell is drawn
				//Printf(" -> becomes %d", mapChPadIdx[ic]);
				nc = counter/nPad;
				np = counter%nPad;
				//Printf("Setting ic %d in cv %d pad %d", ic, nc, np);
				if(nc > bigNumb - 1) {
					Printf("too few pads foreseen");
					continue;
				}
						
				counter++;
				
			} else {
				// check the position where this cell is drawn
				nc = mapChPadIdx[ic]/nPad;
				np = mapChPadIdx[ic]%nPad;
				
				
			}
			//hTmpBad->SetLineColor(colors[ir]);
			hTmpBad->SetLineWidth(2);
			cAmpl[nc]->cd(np+1);
			gPad->SetLogy();
			hTmpBad->DrawClone("sames");
			
			//if(ic == 3778) Printf("Canvas %d, pad %d", mapChPadIdx[ic]/nPad, mapChPadIdx[ic]%nPad);
			if(nc > ncvfilled) ncvfilled = nc;
		} //loop on cells
	
	} //loop on runs
	
	//save to pdf file
	TString filenameout = TString::Format("%s/%s/cAmplitudeSeveralRuns.pdf", workDir.Data(), analysisOutput.Data());
		cAmpl[0]->Print(TString::Format("%s(", filenameout.Data()), "pdf");
		for(Int_t inc = 1; inc < ncvfilled+1; inc++){
			cAmpl[inc]->Print(TString::Format("%s", filenameout.Data()), "pdf");
		}
		//last page will be empty
		cAmpl[ncvfilled+1]->Print(TString::Format("%s)", filenameout.Data()), "pdf");
		
		for(Int_t inc = ncvfilled; inc < bigNumb; inc++){
			delete cAmpl[inc];
		}
	
}

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
