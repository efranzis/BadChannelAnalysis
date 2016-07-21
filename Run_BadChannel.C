// This macro has been developed to find badcell candidates in EMCal and DCal based on cell amplitude distributions
// Input needed can be either outputs QA from AliAnaCalorimeterQA task (BadChannelAnalysis() function)
// Or from merged output of AliAnalysisTaskCaloCellsQA (use BCAnalysis() function)
//
// Author : Alexis Mas (SUBATECH) & M. Germain, based on getCellsRunQA.C from Olga Driga (SUBATECH)
//
// ---------------------
//  Running the macro
// ---------------------
//   root [2] .L BadChannelAnalysis.C++
//   root [2] BadChannelAnalysis("EMCAL","LHC15o","muon_caloLego","AnyINTnoBC",trial=0)
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
//    Further improvement: implement tests 1 and 2 on time distribution histogram
//
/////////////////////////////////////////////////

//#include <TString.h>
//#include <AliCaloChannelAnalysis.cxx>
//#include <AliCaloChannelAnalysis.h>

using namespace std;

//________________________________________________________________________
void Run_BadChannel(TString period = "LHC15f", TString pass = "pass2", TString trigger= "default", Int_t runNum= 254381, TString externalFile= "",TString workDir="./", TString listName="runList.txt")
{
	AliAnaCaloChannelAnalysis* Analysis=new AliAnaCaloChannelAnalysis(period,pass,trigger,runNum,workDir,listName);
	Analysis->SetExternalMergedFile(externalFile);
	Analysis->SetQAChecks(0);  //1=Perform QA checks
	Analysis->SetNTrial(4);  //Perform QA checks
	Analysis->MaskEdgesSM(1);  //Mask the supermodule edges as bad

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

		//inverted order
		/*Analysis->AddPeriodAnalysis(1, 6.,0.5, 2.); // mean energy in range Emin Emax
		Analysis->AddPeriodAnalysis(2, 6.,0.5, 2.); // mean hit in range Emin Emax
		Analysis->AddPeriodAnalysis(1, 6., 2., 5.); // mean energy in range Emin Emax
		Analysis->AddPeriodAnalysis(2, 6., 2., 5.); // mean hit in range Emin Emax
*/

		//Analysis->AddPeriodAnalysis(2, 6., 5.,10.);
		//Analysis->AddPeriodAnalysis(1, 6., 5.,10.);
	}

	Analysis->Run();
}
//________________________________________________________________________
