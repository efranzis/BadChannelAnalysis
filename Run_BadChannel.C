// This macro has been developed to find badcell candidates in EMCal and DCal based on cell amplitude distributions
// Input needed can be either outputs QA from AliAnaCalorimeterQA task (BadChannelAnalysis() function)
// Or from merged output of AliAnalysisTaskCaloCellsQA (use BCAnalysis() function)
//
// Author : Alexis Mas (SUBATECH) & M. Germain, based on getCellsRunQA.C from Olga Driga (SUBATECH)
//
// ---------------------
//  Running the macro
// ---------------------
//   root [2] .L Run_BadChannel.C++
//   root [2] Run_BadChannel("EMCAL","LHC15o","muon_caloLego","AnyINTnoBC",trial=0)
//  
//  !!! pay attention the trigger name depends on the caloQA_triggername you want to analyse check first in QAresults.root what is abvailable
//  !!! it is generally not good to run it on triggered data for the following reasons:
//    
// --------------------
//  outputs
// --------------------
// the output of this analysis povides you:
// intermediate steps files: (those will be recreated each time you rerun so pa attention to save the different files when changing period/listof runs....
//  - Criterum-xx_Emin-xx_Emax-xx.txt : list of identified bad for the different test/Emin/Emax
//  - <period><pass>.txt file with list of dead/bad cells identified 
//  - a pdf file with all energy distributions plots of all bad cells candidates (compared to a reference one (hard coded see Draw function to change)
//
/////////////////////////////////////////////////

#include <TString.h>
#include <TFile.h>
#include <TH2D.h>
#include <TH1D.h>
#include <TList.h>
#include "AliAnaCaloChannelAnalysis.h"


//________________________________________________________________________
void Run_BadChannel(TString period = "LHC15f", TString pass = "pass2", TString trigger= "default", Int_t runNum= 254381, TString externalFile= "",TString workDir="./", TString listName="runList.txt")
{
	AliAnaCaloChannelAnalysis* Analysis=new AliAnaCaloChannelAnalysis(period,pass,trigger,runNum,workDir,listName);
	Analysis->SetExternalMergedFile(externalFile);
	Analysis->SetQAChecks(0);  //1=Perform QA checks
	Analysis->SetNTrial(5);  //Perform QA checks
	//	Analysis->MaskEdgesSM(1);  //Mask the supermodule edges as bad
	if(trigger=="default"||trigger=="INT7"||trigger=="DMC7"||trigger=="AnyINTnoBC")
	{
		Analysis->AddPeriodAnalysis(2, 4.,0.2,0.5); // mean hit in range Emin Emax
		Analysis->AddPeriodAnalysis(2, 4.,0.5, 1.); // mean hit in range Emin Emax
		Analysis->AddPeriodAnalysis(1, 6.,0.5, 1.); // mean energy in range Emin Emax
		Analysis->AddPeriodAnalysis(2, 4., 1., 2.); // mean hit in range Emin Emax
		Analysis->AddPeriodAnalysis(1, 6., 1., 2.); // mean energy in range Emin Emax
		Analysis->AddPeriodAnalysis(2, 4., 1.,10.); // mean hit in range Emin Emax
		Analysis->AddPeriodAnalysis(1, 6., 1.,10.); // mean energy in range Emin Emax
	}
	else
	{
		Analysis->AddPeriodAnalysis(2, 6.,0.5, 2.); // mean hit in range Emin Emax
		Analysis->AddPeriodAnalysis(1, 6.,0.5, 2.); // mean energy in range Emin Emax
		Analysis->AddPeriodAnalysis(2, 6., 2., 5.); // mean hit in range Emin Emax
		Analysis->AddPeriodAnalysis(1, 6., 2., 5.); // mean energy in range Emin Emax
//		Analysis->AddPeriodAnalysis(2, 6., 5., 10.); // mean hit in range Emin Emax
//		Analysis->AddPeriodAnalysis(1, 6., 5., 10.); // mean energy in range Emin Emax

		//inverted order
		/*Analysis->AddPeriodAnalysis(1, 6.,0.5, 2.); // mean energy in range Emin Emax
		Analysis->AddPeriodAnalysis(2, 6.,0.5, 2.); // mean hit in range Emin Emax
		Analysis->AddPeriodAnalysis(1, 6., 2., 5.); // mean energy in range Emin Emax
		Analysis->AddPeriodAnalysis(2, 6., 2., 5.); // mean hit in range Emin Emax
*/
	}

	Analysis->Run();
}
//________________________________________________________________________

void CheckListDeadChannels(TString qaoutputPath = "/data/Work/EMCAL/BadChannels/LHC16h/510muon_caloLego/", TString runlist = "/data/Work/EMCAL/BadChannels/LHC16h/rulListLHC16h_EMCGood_Fieldp30.txt", TString listname = "CaloQA_AnyINT", TString hname = "EMCAL_hAmpId"){

	// current list of dead cells from Martin
	Int_t ndeadSM4 = 32;
	Int_t listDeadSM4[ndeadSM4] = {4833, 4832, 4881, 4880, 4835, 4834, 4883, 4882, 4837, 4836, 4885, 4884, 4839, 4838, 4887, 4886, 4841, 4840, 4889, 4888, 4843, 4842, 4891, 4890, 4845, 4844, 4893, 4892, 4847, 4846, 4895, 4894};
	
	Int_t ndeadSM7 = 32;
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
