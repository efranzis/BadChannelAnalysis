/////////////////////////////////////////////////
//
// This macro has been developed to find badcell candidates in EMCal and DCal based on cell amplitude distributions
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
// root [2] Run_BadChannel("LHC15n","Train_603","AnyINTnoBC",244411)
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

//_______________________________________________________________________
// some helper methods
void CalculatePads(Int_t n, Int_t&nx, Int_t&ny, Int_t&dx, Int_t&dy, Int_t perrow = 5, Int_t stdd = 400);

//________________________________________________________________________
void Run_BadChannel(TString period = "LHC15n", TString train = "Train_603", TString trigger= "AnyINTnoBC", Int_t runNum= 245683, TString externalFile= "",TString listName="runList.txt",TString workDir=".", Int_t nversion = 1)
{
	
	TStopwatch watch;
	watch.Start();
	
	AliAnaCaloChannelAnalysis* Analysis;
	//..If you do the analysis run by run - this might be helpful
	//Analysis=new AliAnaCaloChannelAnalysis(period,train,trigger,runNum,runNum,workDir,listName);

	//..If you do it with merged files use this
	Analysis=new AliAnaCaloChannelAnalysis(period,train,trigger,runNum,nversion,workDir,listName);

	//..Settings
	Analysis->SetExternalMergedFile(externalFile);
	//Analysis->SetQAChecks(1);  //1=Perform QA checks - takes a long time! Prints all good cells for cross check

	//. . . . . . . . . . . . . . . . . . . . . . . .
	//. . Add different period analyses
	//. . . . . . . . . . . . . . . . . . . . . . . .
	Analysis->AddPeriodAnalysis(2, 4.,0.1,0.3); // hit/event in range Emin Emax
	Analysis->AddPeriodAnalysis(1, 6.,0.1,0.3); // energy/hit in range Emin Emax
	//Analysis->AddPeriodAnalysis(2, 4.,0.2,0.5); // hit/event in range Emin Emax
	//Analysis->AddPeriodAnalysis(1, 6.,0.2,0.5); // energy/hit in range Emin Emax
	//Analysis->AddPeriodAnalysis(2, 4.,0.5,1.0); // hit/event in range Emin Emax
	//Analysis->AddPeriodAnalysis(1, 6.,0.5,1.0); // energy/hit in range Emin Emax

	//..If there is enough statistic add also these:
	//Analysis->AddPeriodAnalysis(2, 4.,1.0,4.0); // hit/event in range Emin Emax
	//Analysis->AddPeriodAnalysis(1, 6.,1.0,4.0); // mean energy in range Emin Emax
//	Analysis->AddPeriodAnalysis(2, 4.,1.0,10.0);// hit/event in range Emin Emax
//	Analysis->AddPeriodAnalysis(1, 5.,1.0,10.0);// energy/hit in range Emin Emax
	//Analysis->AddPeriodAnalysis(2, 5.,3.0,10.0);// (PbPb extra range) hit/event in range Emin Emax
	//Analysis->AddPeriodAnalysis(1, 5.,3.0,10.0);// (PbPb extra range) energy/hit in range Emin Emax

///*test time stuff*/	Analysis->AddPeriodAnalysis(3, 6,-20,+20);// energy/hit in range Emin Emax


	//..Start the bad channel analysis
	Analysis->Run();
	
	watch.Stop();
	watch.Print();
}
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

	TH2F* h2DChannelMap_FlagBad=(TH2F*)outputRoot->Get("2DChannelMap_Flag2");
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
	TCanvas* C4 = new TCanvas("Division of OADB/Orig.","Division of OADB/Orig.",1);
	C4->Divide(2);
	C4->cd(1);
	plot2D_Bad_OADB->Add(h2DChannelMap_FlagBad,-1);
	plot2D_Bad_OADB->DrawCopy("colz");
	C4->cd(2);
	plot2D_Dead_OADB->Add(h2DChannelMap_FlagDead,-1);
	plot2D_Dead_OADB->DrawCopy("colz");
}
//________________________________________________________________________
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
	
	TCanvas *cBad [nCv];
	TCanvas *cGood[nCv];
	TCanvas *cDead[nCv];
	TCanvas *cAmp [nCv];
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
			continue;
		}

		//TH2F *badCells  = (TH2F*)f->Get("HitRowColumn_Flag2");
		//TH2F *goodCells = (TH2F*)f->Get("HitRowColumn_Flag0");
		//TH2F *deadCells = (TH2F*)f->Get("HitRowColumn_Flag1");
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
	gSystem->mkdir(TString::Format("%s/%s/", analysisOutput.Data(), train.Data()));
	gSystem->mkdir(TString::Format("%s/%s/RunByRunSummary/", analysisOutput.Data(), train.Data()));
	for(Int_t ic = 0; ic<nCv; ic++){
		cBad [ic] ->SaveAs(TString::Format("%s/%s/RunByRunSummary/%s.gif", analysisOutput.Data(), train.Data(), cBad [ic]->GetName()));
		cGood[ic] ->SaveAs(TString::Format("%s/%s/RunByRunSummary/%s.gif", analysisOutput.Data(), train.Data(), cGood[ic]->GetName()));
		cDead[ic] ->SaveAs(TString::Format("%s/%s/RunByRunSummary/%s.gif", analysisOutput.Data(), train.Data(), cDead[ic]->GetName()));
		cAmp [ic] ->SaveAs(TString::Format("%s/%s/RunByRunSummary/%s.gif", analysisOutput.Data(), train.Data(), cAmp [ic]->GetName()));
	}
}
//________________________________________________________________________
void CheckListDeadChannels(TString qaoutputPath = "/data/Work/EMCAL/BadChannels/LHC16h/510muon_caloLego/", TString runlist = "/data/Work/EMCAL/BadChannels/LHC16h/rulListLHC16h_EMCGood_Fieldp30.txt", TString listname = "CaloQA_AnyINT", TString hname = "EMCAL_hAmpId")
{

	// current list of dead cells from Martin
	const Int_t ndeadSM4 = 32;
	Int_t listDeadSM4[ndeadSM4] = {4833, 4832, 4881, 4880, 4835, 4834, 4883, 4882, 4837, 4836, 4885, 4884, 4839, 4838, 4887, 4886, 4841, 4840, 4889, 4888, 4843, 4842, 4891, 4890, 4845, 4844, 4893, 4892, 4847, 4846, 4895, 4894};

	const Int_t ndeadSM7 = 32;
	Int_t listDeadSM7[ndeadSM7] = {8481, 8480, 8529, 8528, 8483, 8482, 8531, 8530, 8485, 8484, 8533, 8532, 8487, 8486, 8535, 8534, 8489, 8488, 8537, 8536, 8491, 8490, 8539, 8538, 8493, 8492, 8541, 8540, 8495, 8494, 8543, 8542};

	ifstream read(runlist.Data());
	Int_t nruns = 0;
	TString currentRun = "";
	while(read){
		read>> currentRun;
		Printf("Run %s", currentRun.Data());

		TFile *fin = new TFile(Form("%s%s.root", qaoutputPath.Data(), currentRun.Data()));

		if(!fin->IsOpen()){
			Printf("%s%s.root not found", qaoutputPath.Data(), currentRun.Data());
			continue;
		}

		TDirectoryFile *dir = (TDirectoryFile*)fin->Get(listname);
		if(!dir){
			Printf("%s not found", listname.Data());
			fin->ls();
			continue;
		}

		TList *list = (TList*)dir->Get(listname);
		if(!list){
			Printf("%s not found", listname.Data());
			dir->ls();
			continue;
		}

		TH2D *hAmpCellId = (TH2D*)list->FindObject(hname);
		if(!hAmpCellId){
			Printf("hAmpCellId not found");
			continue;
		}

		for(Int_t ideadsm4 = 0; ideadsm4 < ndeadSM4; ideadsm4++){
			TH1D* hAmpForSpecificID = hAmpCellId->ProjectionX(Form("hAmpForSpecificID"), listDeadSM4[ideadsm4], listDeadSM4[ideadsm4]);
			Int_t entries = hAmpForSpecificID->GetEntries();
			if(entries > 0) {
				Printf("Error! %d entries found in Cell %d, run %s", entries, listDeadSM4[ideadsm4], currentRun.Data());
			}

			delete hAmpForSpecificID;
		}

		for(Int_t ideadsm7 = 0; ideadsm7 < ndeadSM7; ideadsm7++){
			TH1D* hAmpForSpecificID = hAmpCellId->ProjectionX(Form("hAmpForSpecificID"), listDeadSM7[ideadsm7], listDeadSM7[ideadsm7]);
			Int_t entries = hAmpForSpecificID->GetEntries();
			if(entries > 0) {
				Printf("Error! %d entries found in Cell %d, run %s", entries, listDeadSM7[ideadsm7], currentRun.Data());
			}
			delete hAmpForSpecificID;
		}

		nruns++;

	}
}



//---------------------------------------------------------------------------------------------------------------------------------------------------------
void CalculatePads(Int_t n, Int_t&nx, Int_t&ny, Int_t&dx, Int_t&dy, Int_t perrow, Int_t stdd){
   
   Int_t percolumn = n/perrow;
   if(n%perrow > 0) percolumn++;
   nx = percolumn;
   ny = perrow;
   dx = stdd*nx;
   dy= stdd*ny;
   
   return;
}
