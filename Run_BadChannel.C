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
/*#include "AliAnaCaloChannelAnalysis.h" //include when compile
#include "AliEMCALGeometry.h"          //include when compile
#include "AliCalorimeterUtils.h"       //include when compile
#include "AliAODEvent.h"               //include when compile
*/
//________________________________________________________________________
void Run_BadChannel(TString period = "LHC15n", TString train = "Train_603", TString trigger= "AnyINTnoBC", Int_t runNum= 244411, TString externalFile= "",TString workDir=".", TString listName="runList.txt")
{
	AliAnaCaloChannelAnalysis* Analysis;
	//..If you do the analysis run by run - this might be helpful
	//Analysis=new AliAnaCaloChannelAnalysis(period,train,trigger,runNum,runNum,workDir,listName);

	//..If you do it with merged files use this
	Analysis=new AliAnaCaloChannelAnalysis(period,train,trigger,runNum,2,workDir,listName);

	//..Settings
	Analysis->SetExternalMergedFile(externalFile);
	Analysis->SetQAChecks(0);  //1=Perform QA checks
	//Analysis->MaskEdgesSM(1);  //Mask the supermodule edges as bad - not yet implemented

	//. . . . . . . . . . . . . . . . . . . . . . . .
	//. . Add different period analyses
	//. . . . . . . . . . . . . . . . . . . . . . . .
	Analysis->AddPeriodAnalysis(2, 4.,0.1,0.3); // hit/event in range Emin Emax
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

	//..Start the bad channel analysis
	Analysis->Run();
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
//________________________________________________________________________
void SummarizeRunByRun(TString period = "LHC15f", TString pass = "pass2",TString workDir=".", TString listName="runList.txt")
{
	//..Open the text file with the run list numbers and run index
	cout<<"o o o Open .txt file with run indices. Name = " << listName << endl;
	TString analysisInput  = Form("AnalysisInput/%s",period.Data());
	TString analysisOutput = Form("AnalysisOutput/%s",period.Data());
	TString runList        = Form("%s/%s/%s/%s",workDir.Data(), analysisInput.Data(), pass.Data(), listName.Data());

	FILE *pFile = fopen(runList.Data(), "r");
	if(!pFile)
	{
		cout<<"couldn't open file!"<<endl;
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

	TCanvas *cBad  = new TCanvas("badcells","badcells",1000,750);
	TCanvas *cGood = new TCanvas("goodcells","goodcells",1000,750);
	TCanvas *cDead = new TCanvas("deadcells","deadcells",1000,750);
	TCanvas *cAmp  = new TCanvas("Amplitide","Amplitide",1000,750);
	if(nRun > 16)
	{
		cBad->Divide(5,5);
		cGood->Divide(5,5);
		cDead->Divide(5,5);
		cAmp->Divide(5,5);
	}
	else if(nRun > 9)
	{
		cBad->Divide(4,4);
		cGood->Divide(4,4);
		cDead->Divide(4,4);
		cAmp->Divide(4,4);
	}
	else if(nRun > 6)
	{
		cBad->Divide(3,3);
		cGood->Divide(3,3);
		cDead->Divide(3,3);
		cAmp->Divide(3,3);
	}

	//..loop over the amount of run numbers found in the previous text file.
	for(Int_t i = 0 ; i < nRun ; i++)  //Version%i" LHC16i_muon_caloLego_Histograms_V255539
	{
		rootFileName      = Form("%s_%s_Histograms_V%i.root",period.Data(),pass.Data(),RunId[i]);
		badChannelOutput  = Form("%s/Version%i/%s", analysisOutput.Data(), RunId[i],rootFileName.Data());

		cout<<"Open root file: "<<badChannelOutput<<endl;
		TFile *f = TFile::Open(badChannelOutput);
		if(!f)
		{
			cout<<"Couldn't open/find .root file: "<<badChannelOutput<<endl;
			continue;
		}

		TH2F *badCells  = (TH2F*)f->Get("HitRowColumn_Flag2");
		TH2F *goodCells = (TH2F*)f->Get("HitRowColumn_Flag0");
		TH2F *deadCells = (TH2F*)f->Get("HitRowColumn_Flag1");
		TH2F *ampID     = (TH2F*)f->Get("hCellAmplitude");
		if(i<25)
		{
			//....................................
			cBad->cd(i+1);
			badCells->Draw("colz");
			infoText=Form("Bad Cells - Run %i",RunId[i]);
			TLatex* text = new TLatex(0.2,0.8,infoText);
			text->SetTextSize(0.06);
			text->SetNDC();
			text->SetTextColor(1);
			text->Draw();
			//....................................
			cGood->cd(i+1);
			goodCells->Draw("colz");
			infoText=Form("Good Cells - Run %i",RunId[i]);
			TLatex* text1 = new TLatex(0.2,0.8,infoText);
			text1->SetTextSize(0.06);
			text1->SetNDC();
			text1->SetTextColor(1);
			text1->Draw();
			//....................................
			cDead->cd(i+1);
			deadCells->Draw("colz");
			infoText=Form("Dead Cells - Run %i",RunId[i]);
			TLatex* text2 = new TLatex(0.2,0.8,infoText);
			text2->SetTextSize(0.06);
			text2->SetNDC();
			text2->SetTextColor(1);
			text2->Draw();
			//....................................
			cAmp->cd(i+1)->SetLogz();
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
	}
}
