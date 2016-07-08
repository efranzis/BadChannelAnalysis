// This macro has been developed to find badcell candidates in EMCal and DCal based on cell amplitude distributions
// Input needed can be either outputs QA from AliAnaCalorimeterQA task (BadChannelAnalysis() function)
// Or from merged output of AliAnalysisTaskCaloCellsQA (use BCAnalysis() function)
//
// Author : Alexis Mas (SUBATECH) & M. Germain, based on getCellsRunQA.C from Olga Driga (SUBATECH)
//
//-----------------
// Main method:
//---------------
// BadChannelAnalysis:
//
//    step 1 : Convert()
//       read list of mergeable runs  in your working directory
//       (in example below the $workdir is "/scratch/alicehp2/germain/QANew2/"
//       The QAresults.root files should be aleady copied from alien and be in $workdir/<period>/<pass>/runnb.root
//       read/merge the histos"EMCAL_hAmpId" and "EMCAL_hTimeId"  from QAresults.root file and write them in 
//       $workdir/period/pass/<period><pass>Runlist0New.root"  !!! this is hardcoded !!!!!
//       step 1 has to be called only the first time runing on a new list
//
//    step 2 BCanalysis() main method to analyse previously created file (hardcoded)
// 
//       call of different Periodanalysis(criterium,..) functions according to the wanted tests (critreria) 
//           1 : average E for E>Emin
//           2 : entries for E>Emin
//           3 : ki²/ndf  (from fit of each cell Amplitude between Emin and Emax) 
//           4 : A parameter (from fit of each cell Amplitude between Emin and Emax) 
//           5 : B parameter (from fit of each cell Amplitude between Emin and Emax) 
//           6 : 
//           7 : give bad + dead list
//
//       Mainly used: 1 and 2 (with different settings (chi2, intervals of energy : this is quite dependent of the stat you may have )
//       Further improvement: implement tests 1 and 2 on time distribustion histogram
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
void Run_BadChannel(TString period = "LHC15f", TString pass = "pass2", TString trigger= "default", Int_t runNum= 254381, TString externalFile= "")
{
	AliCaloChannelAnalysis* Analysis=new AliCaloChannelAnalysis(period,pass,trigger,runNum);
	Analysis->SetExternalMergedFile(externalFile);

	if(fTrigger=="default"||fTrigger=="INT7"||fTrigger=="DMC7"||fTrigger=="AnyINTnoBC")
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
		Analysis->AddPeriodAnalysis(2, 6.,0.5, 2.);
		Analysis->AddPeriodAnalysis(1, 6.,0.5, 2.);
		Analysis->AddPeriodAnalysis(2, 6., 2., 5.);
		Analysis->AddPeriodAnalysis(1, 6., 2., 5.);
		//Analysis->AddPeriodAnalysis(2, 6., 5.,10.);
		//Analysis->AddPeriodAnalysis(1, 6., 5.,10.);
	}

	Analysis->Run();
}
//________________________________________________________________________
