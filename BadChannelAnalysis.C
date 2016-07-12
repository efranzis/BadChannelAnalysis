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


#if !defined(__CINT__) || defined(__MAKECINT__) 
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TF1.h>
#include <TROOT.h>
#include <TLine.h>
#include <Riostream.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TLegend.h>
#include <TString.h>
#include <AliCalorimeterUtils.h>
#include <AliAODEvent.h>
#include <TLatex.h>
#include <stdio.h>
#include <fstream>
#include <iostream>
#include <alloca.h>
#include <string>
#include <cstring>
#endif
using namespace std;

//___________________________________________________________________________________________________
/*//copied from AliCalorimeterUtils
Int_t GetModuleNumberCellIndexesAbsCaloMap(Int_t absId, Int_t calo,Int_t & icol   , Int_t & irow,
										  Int_t & iRCU,Int_t & icolAbs, Int_t & irowAbs)
{
	Int_t imod = GetModuleNumberCellIndexes(absId, calo, icol, irow,iRCU);

	icolAbs = icol;
	irowAbs = irow;
	//
	// EMCal/DCal
	//
	// Shift collumns in even SM
	Int_t shiftEta = 48;

	// Shift collumn even more due to smaller acceptance of DCal collumns
	if ( imod >  11 && imod < 18) shiftEta+=48/3;

	icolAbs = (imod % 2) ? icol + shiftEta : icol;

	//
	// Shift rows per sector
	irowAbs = irow + 24 * Int_t(imod / 2);

	// Shift row less due to smaller acceptance of SM 10 and 11 to count DCal rows
	if ( imod >  11 && imod < 20) irowAbs -= (2*24 / 3);

	return imod ;
}*/
void Draw2(Int_t cell, Int_t cellref=400)
{
    //ELI not used anywhere
	gROOT->SetStyle("Plain");
	gStyle->SetOptStat(0);
	gStyle->SetFillColor(kWhite);
	gStyle->SetTitleFillColor(kWhite);
	gStyle->SetPalette(1);
	char out[120]; char title[100]; char name[100];char name2[100];
	TString slide(Form("Cells %d-%d",cell,cell));

	sprintf(out,"%d.gif",cell);
	TH2 *hCellAmplitude = (TH2*) gFile->Get("hCellAmplitude");
	TH1 *hCellref = hCellAmplitude->ProjectionX("badcells",cellref+1,cellref+1);

	TCanvas *c1 = new TCanvas("badcells","badcells",600,600) ;
	c1->SetLogy();

	// hCellref->Rebin(3);
	TLegend *leg = new TLegend(0.7, 0.7, 0.9, 0.9);

	sprintf(name,"Cell %d",cell) ;
	TH1 *hCell = hCellAmplitude->ProjectionX(name,cell+1,cell+1);

	sprintf(title,"Cell %d      Entries : %d Enties ref: %d",cell, (Int_t)hCell->GetEntries(),(Int_t)hCellref->GetEntries()) ;
	hCell->SetLineColor(2)  ;
	// cout<<title<<endl ;
	hCell->SetMaximum(1e5);
	// hCell->Rebin(3);
	hCell->SetAxisRange(0.,10.);
	hCell->GetXaxis()->SetTitle("E (GeV)");
	hCell->GetYaxis()->SetTitle("N Entries");
	hCellref->SetAxisRange(0.,10.);
	hCell->SetLineWidth(1) ;
	hCellref->SetLineWidth(1) ;
	hCell->SetTitle(title);
	hCellref->SetLineColor(1)  ;
	leg->AddEntry(hCellref,"reference","l");
	leg->AddEntry(hCell,"current","l");
	hCell->Draw() ;
	hCellref->Draw("same") ;
	leg->Draw();
	sprintf(name2,"Cell%dLHC13MB.gif",cell) ;
	c1->SaveAs(name2);

}

void SaveBadCellsToPDF(Int_t cell[], Int_t iBC, Int_t nBC, TString PdfName, const Int_t cellref=2377)
{
	//Allow to produce a pdf file with badcells candidates (red) compared to a refence cell (black)
	gROOT->SetStyle("Plain");
	gStyle->SetOptStat(0);
	gStyle->SetFillColor(kWhite);
	gStyle->SetTitleFillColor(kWhite);
	gStyle->SetPalette(1);
	//char out[120];
	char title[100];
	char name[100];
	if(cell[iBC]==-1)cout<<"### strange shouldn't happen 1, cell id: "<<iBC<<endl;
	if(cell[iBC+nBC-1]==-1)cout<<"### strange shouldn't happen 2"<<endl;

	TString slide     = Form("Cells %d-%d",cell[iBC],cell[iBC+nBC-1]);
    TString reflegend = Form("reference Cell %i",cellref);
	//sprintf(out,"%d-%d.gif",cell[iBC],cell[iBC+nBC-1]);

	TH2 *hCellAmplitude = (TH2*) gFile->Get("hCellAmplitude");
	TH1 *hCellref = hCellAmplitude->ProjectionX("badcells",cellref+1,cellref+1);

	TCanvas *c1 = new TCanvas("badcells","badcells",1000,750);
	if(nBC > 6) c1->Divide(3,3);
	else if (nBC > 3)  c1->Divide(3,2);
	else  c1->Divide(3,1);

	TLegend *leg = new TLegend(0.7, 0.7, 0.9, 0.9);
	for(Int_t i=0; i<nBC ; i++)
	{
		sprintf(name, "Cell %d",cell[iBC+i]) ;
		TH1 *hCell = hCellAmplitude->ProjectionX(name,cell[iBC+i]+1,cell[iBC+i]+1);
		sprintf(title,"Cell %d      Entries : %d  Ref : %d",cell[iBC+i], (Int_t)hCell->GetEntries(), (Int_t)hCellref->GetEntries() ) ;

		c1->cd(i%9 + 1);
		c1->cd(i%9 + 1)->SetLogy();
		hCell->SetLineColor(2);
		hCell->SetMaximum(1e6);
		hCell->SetAxisRange(0.,10.);
		hCell->GetXaxis()->SetTitle("E (GeV)");
		hCell->GetYaxis()->SetTitle("N Entries");
		hCell->SetLineWidth(1) ;
		hCell->SetTitle(title);
		hCellref->SetAxisRange(0.,8.);
		hCellref->SetLineWidth(1);
		hCellref->SetLineColor(1);

		if(i==0)
		{
			leg->AddEntry(hCellref,reflegend,"l");
			leg->AddEntry(hCell,"current","l");
		}
		hCell->Draw() ;
		hCellref->Draw("same") ;
		leg->Draw();
	}

	//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	//..Store the created canvas in a .pdf file
	//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	if(nBC<9)
	{
		PdfName +=")";  //ELI this is strange
		c1->Print(PdfName.Data());
	}
	else if(iBC==0)
	{
		PdfName +="("; //ELI this is strange
		c1->Print(PdfName.Data());
	}
	else  c1->Print(PdfName.Data());

	delete hCellref;
	delete c1;
	delete leg;
}
//_________________________________________________________________________
//_________________________________________________________________________
TString Convert(TString period = "LHC11h", TString pass = "pass1_HLT", TString trigger= "default", TString runlistFileName = "runList.txt", TString inputPath = "AnalysisInput", TString outPath = "ConvertOutput")
{
	//parameters for folder sturcture: period, pass
	//parameter very important for file name: trigger

	//..Create one file for the analysis from several outputs QA files listed in runlist.txt
	//..You need :
	//..runlist.txt with runs listed
	//..outputsQA  e.g  period/pass/123456.root
    cout<<"o o o Start conversion process o o o"<<endl;
    cout<<"o o o period: " << period << ", pass: " << pass << ",  trigger: "<<trigger<< endl;

    //.. Create histograms needed for...
    TH1D *hNEventsProcessedPerRun = new TH1D("hNEventsProcessedPerRun","Number of processed events vs run number",200000,100000,300000);
    //ELI a little problematic to hard code properties of histograms??
    TH2F *hCellAmplitude          = new TH2F("hCellAmplitude","Cell Amplitude",200,0,10,23040,0,23040);
    TH2F *hCellTime               = new TH2F("hCellTime","Cell Time",250,-275,975,23040,0,23040);

    //..Open the text file with the run list numbers and run index
    /*ELI*/TString file = Form("%s/%s/%s/%s", inputPath.Data(), period.Data(), pass.Data(), runlistFileName.Data());
    cout<<"o o o Open .txt file with run indices. Name = " << file << endl;
    FILE *pFile = fopen(file.Data(), "r");
    if(!pFile)cout<<"count't open file!"<<endl;
    Int_t Nentr;
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
	TString base;
	
	gSystem->mkdir(outPath);
	
	TString BCfile= Form("%s/%s%sConverted.root", outPath.Data(), period.Data(),pass.Data());

	TString direct(Form("CaloQA_%s",trigger.Data()));
	TString infile;

    cout<<"o o o Start merging process of " << nRun <<" files"<< endl;
	//..loop over the amount of run numbers found in the previous text file.
	for(Int_t i = 0 ; i < nRun ; i++)
	{
		base  = Form("%s/%s/%s/%d", inputPath.Data(), period.Data(), pass.Data(), RunId[i]);

		// ICI on met le nom period/pass/runblabla.root
		if ((pass=="cpass1_pass2")||(pass=="cpass1-2"))
		{
			if (trigger=="default")
			{
				infile = Form("%s_barrel.root",base.Data());
			}
			else
			{
				infile = Form("%s_outer.root",base.Data());
			}
		}
		else
		{   //..This is our run2 case
			infile = Form("%s.root",base.Data()) ;
		}

		cout<<"    o Open .root file with name: "<<infile<<endl;
		TFile *f = TFile::Open(infile);

		//base=Form("%s/%s",base.Data(),trigger.Data()); //CHIARA: why changing the variable base and then not using it?
		//..Do some basic checks
		if(!f)
		{
			cout<<"Couldn't open/find .root file: "<<infile<<endl;
			continue;
		}
		TDirectoryFile *dir = (TDirectoryFile *)f->Get(direct);
		if(!dir)
		{
			cout<<"Couln't open directory file in .root file: "<<infile<<", directory: "<<direct<<endl;
			continue;
		}
		TList *outputList = (TList*)dir->Get(direct);
		if(!outputList)
		{
			cout << "Couln't get TList from directory file: "<<direct<<endl;
			continue;
		}

		TH2F *hAmpId;
		TH2F *hTimeId;
		TH2F *hNEvents;

		hAmpId =(TH2F *)outputList->FindObject("EMCAL_hAmpId");
		if(!hAmpId)
		{
			Printf("hAmpId not found");
			outputList->ls();
			continue;
		}
		hTimeId =(TH2F *)outputList->FindObject("EMCAL_hTimeId");
		if(!hTimeId)
		{
			Printf("hTimeId not found");
			outputList->ls();
			continue;
		}
		hNEvents =(TH2F *)outputList->FindObject("hNEvents");
		if(!hNEvents)
		{
			Printf("hNEvents not found");
			outputList->ls();
			continue;
		}

		Nentr =  (Int_t)hNEvents->GetEntries();

		//..does that mean do not merge small files?
		if (Nentr<100)
		{
			cout <<"    o File to small to be merged. Only N entries " << Nentr << endl;
			continue ;
		}
		cout <<"    o File with N entries " << Nentr<<" will be merged"<< endl;

		hNEventsProcessedPerRun->SetBinContent(RunId[i]-100000,(Double_t)Nentr);
		hCellAmplitude->Add(hAmpId);
		hCellTime->Add(hTimeId);

		//if(i==0){ cout<<"Merging/Converting procedure ..." ; cout.flush();}
		//else { cout<<"..." ; cout.flush();}
		outputList->Delete();
		dir->Delete();
		f->Close();
		delete f;
	}

    //.. Save the merged histograms
	cout<<"o o o Save the merged histogramms to .root file with name: "<<BCfile<<endl;
	TFile *BCF = TFile::Open(BCfile,"recreate");
	hNEventsProcessedPerRun->Write();
	hCellAmplitude->Write();
	hCellTime->Write();
	BCF->Close();
    cout<<"o o o End conversion process o o o"<<endl;
    return BCfile;
}

//_________________________________________________________________________
//_________________________________________________________________________

void Process(Int_t *pflag[23040][7], TH1* inhisto, Double_t Nsigma = 4., Int_t dnbins = 200, Double_t dmaxval = -1., TString dirOut = "BadChannelOutput")
{  
	//  1) create a distribution for the input histogram;
	//  2) fit the distribution with a gaussian
	//  3) define good area within +-Nsigma to identfy badcells
	//
	// inhisto -- input histogram;
	// dnbins  -- number of bins in distribution;
	// dmaxval -- maximum value on distribution histogram.

	gStyle->SetOptStat(1); // MG modif
	gStyle->SetOptFit(1);  // MG modif
	Int_t crit = *pflag[0][0] ; //identify the criterum processed
	if(crit==1)cout<<"    o Fit average energy per hit distribution"<<endl;
	if(crit==2)cout<<"    o Fit average hit per event distribution"<<endl;

	//..setings for the 2D histogram
	Int_t fNMaxCols = 48;  //eta direction
	Int_t fNMaxRows = 24;  //phi direction
	Int_t fNMaxColsAbs=2*fNMaxCols;
	Int_t fNMaxRowsAbs=Int_t (20/2)*fNMaxRows; //multiply by number of supermodules (20)
	Int_t CellColumn=0,CellRow=0;
	Int_t CellColumnAbs=0,CellRowAbs=0;
	Int_t Trash;
	//..load necessary libraries
    AliCalorimeterUtils* fCaloUtils = new AliCalorimeterUtils();
    //..AccessGeometry needs an input event to retrieve the run number, name, GetPHOSMatrix, GetEMCALMatrix
	//..
    AliAODEvent* aod = new AliAODEvent();
    aod->SetRunNumber(254381); //will not work
    cout<<"current run number: "<<aod->GetRunNumber()<<" , name: "<<aod->GetName()<<endl;
    fCaloUtils->SetRunNumber(254381);
    fCaloUtils->AccessGeometry(aod); // InputEvent()->GetRunNumber()
    //..Set the AODB calibration, bad channels etc. parameters at least once
    //fCaloUtils->AccessOADB(aod);
    //..apparently not initialized correctly like eg in AliEMCALGeometry!
    //fCaloUtils->SetNumberOfSuperModulesUsed(20);
    cout<<"get number of cells: "<<fCaloUtils->GetEMCALGeometry()->GetNCells()<<endl;
    cout<<"get number of supermod: "<<fCaloUtils->GetEMCALGeometry()->GetNumberOfSuperModules()<<endl;
    cout<<"get number of supermod utils: "<<fCaloUtils->GetNumberOfSuperModulesUsed()<<endl;

	TString HistoName=inhisto->GetName();
	Double_t goodmax= 0. ;
	Double_t goodmin= 0. ;
	*pflag[0][0] =1; //ELI why is this done??
	if (dmaxval < 0.)
	{
		dmaxval = inhisto->GetMaximum()*1.01;  // 1.01 - to see the last bin
		if(crit==2 && dmaxval > 1) dmaxval =1. ;
	}

	//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	//. . .build the distribution of average values
	//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	TH1 *distrib = new TH1F(Form("%sDistr",(const char*)HistoName), "", dnbins, inhisto->GetMinimum(), dmaxval);
	distrib->SetXTitle(inhisto->GetYaxis()->GetTitle());
	distrib->SetYTitle("Entries");
	//distrib->GetXaxis()->SetNdivisions(505);
	//..fill the distribution of avarge cell values
	for (Int_t c = 1; c <= inhisto->GetNbinsX(); c++)
	{
		distrib->Fill(inhisto->GetBinContent(c));
	}
	//..build two dimensional histogram with values row vs. column
	TH2F *Plot2D = new TH2F(Form("%s_HitRowColumn",(const char*)HistoName),Form("%s_HitRowColumn",(const char*)HistoName),fNMaxColsAbs+2,-1.5,fNMaxColsAbs+0.5, fNMaxRowsAbs+2,-1.5,fNMaxRowsAbs+0.5);
	Plot2D->GetXaxis()->SetTitle("cell column (#eta direction)");
	Plot2D->GetYaxis()->SetTitle("cell row (#phi direction)");


	for (Int_t c = 1; c <= inhisto->GetNbinsX(); c++)
	{
		//..Do that only for cell ids also accepted by the
		if(!fCaloUtils->GetEMCALGeometry()->CheckAbsCellId(c-1))continue;
		//..Get Row and Collumn for cell ID c
		fCaloUtils->GetModuleNumberCellIndexesAbsCaloMap(c-1,0,CellColumn,CellRow,Trash,CellColumnAbs,CellRowAbs);
		if(CellColumnAbs> fNMaxColsAbs || CellRowAbs>fNMaxRowsAbs)
		{
			cout<<"Problem! wrong calculated number of max col and max rows"<<endl;
			cout<<"current col: "<<CellColumnAbs<<", max col"<<fNMaxColsAbs<<endl;
			cout<<"current row: "<<CellRowAbs<<", max row"<<fNMaxRowsAbs<<endl;
		}
		Plot2D->SetBinContent(CellColumnAbs,CellRowAbs,inhisto->GetBinContent(c));
	}

	//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	//. . .draw histogram + distribution
	//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

	//Produ
	TCanvas *c1 = new TCanvas(HistoName,HistoName,900,900);
	c1->ToggleEventStatus();
	TPad*    upperPad    = new TPad("upperPad", "upperPad",.005, .5, .995, .995);
	TPad*    lowerPadLeft = new TPad("lowerPadL", "lowerPadL",.005, .005, .5, .5);
	TPad*    lowerPadRight = new TPad("lowerPadR", "lowerPadR",.5, .005, .995, .5);
	upperPad->Draw();
	lowerPadLeft->Draw();
	lowerPadRight->Draw();

	upperPad->cd();
	upperPad->SetLeftMargin(0.045);
	upperPad->SetRightMargin(0.03);
	upperPad->SetLogy();
	inhisto->SetTitleOffset(0.6,"Y");
	inhisto->GetXaxis()->SetRangeUser(0,17000);

	inhisto->SetLineColor(kBlue+1);
	inhisto->Draw();

	lowerPadRight->cd();
	lowerPadRight->SetLeftMargin(0.09);
	lowerPadRight->SetRightMargin(0.06);
	Plot2D->Draw("colz");

	lowerPadLeft->cd();
	lowerPadLeft->SetLeftMargin(0.09);
	lowerPadLeft->SetRightMargin(0.06);
	lowerPadLeft->SetLogy();
	distrib->SetLineColor(kBlue+1);
	distrib->Draw();

	//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	//. . .fit histogram
	//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	Int_t higherbin=0,i;
	for(i = 2; i <= dnbins; i++)
	{
		if(distrib->GetBinContent(higherbin) < distrib->GetBinContent(i))  higherbin = i ;
	}
	//..good range is around the max value as long as the
	//..bin content is larger than 2 entries
	for(i = higherbin ; i<=dnbins ; i++)
	{
		if(distrib->GetBinContent(i)<2) break ;
		goodmax = distrib->GetBinCenter(i);
	}
	for(i = higherbin ; i>1 ; i--)
	{
		if(distrib->GetBinContent(i)<2) break ;
		goodmin = distrib->GetBinLowEdge(i);
	}
    //cout<<"higherbin : "<<higherbin<<endl;
    //cout<<"good range : "<<goodmin<<" - "<<goodmax<<endl;

	TF1 *fit2 = new TF1("fit2", "gaus");
	//..start the fit with a mean of the highest value
	fit2->SetParameter(1,higherbin);

	distrib->Fit(fit2, "0LQEM", "", goodmin, goodmax);
	Double_t sig, mean, chi2ndf;
	// Marie midif to take into account very non gaussian distrig
	mean    = fit2->GetParameter(1);
	sig     = fit2->GetParameter(2);
	chi2ndf = fit2->GetChisquare()/fit2->GetNDF();

	if (mean <0.) mean=0.; //ELI is this not a highly problematic case??

	goodmin = mean - Nsigma*sig ;
	goodmax = mean + Nsigma*sig ;

	if (goodmin<0) goodmin=0.;

	cout<<"    o Result of fit: "<<endl;
	cout<<"    o  "<<endl;
	cout<<"    o Mean: "<<mean <<" sigma: "<<sig<<endl;
	cout<<"    o good range : "<<goodmin <<" - "<<goodmax<<endl;

	//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	//. . .Add info to histogram
	//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	TLine *lline = new TLine(goodmin, 0, goodmin, distrib->GetMaximum());
	lline->SetLineColor(kGreen+2);
	lline->SetLineStyle(7);
	lline->Draw();

	TLine *rline = new TLine(goodmax, 0, goodmax, distrib->GetMaximum());
	rline->SetLineColor(kGreen+2);
	rline->SetLineStyle(7);
	rline->Draw();

	TLegend *leg = new TLegend(0.60,0.82,0.9,0.88);
	leg->AddEntry(lline, "Good region boundary","l");
	leg->Draw("same");

	fit2->SetLineColor(kOrange-3);
	fit2->SetLineStyle(1);//7
	fit2->Draw("same");

	TLatex* text = 0x0;
	if(crit==1) text = new TLatex(0.2,0.8,Form("Good range: %.2f-%.2f",goodmin,goodmax));
	if(crit==2) text = new TLatex(0.2,0.8,Form("Good range: %.2f-%.2fx10^-5",goodmin*100000,goodmax*100000));
	text->SetTextSize(0.06);
	text->SetNDC();
	text->SetTextColor(1);
	//text->SetTextAngle(angle);
	text->Draw();
	//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	//. . .Save histogram
	//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	c1->Update();
	gSystem->mkdir(dirOut);
	TString name   =Form("%s/criteria-_%d.gif", dirOut.Data(), crit);
	if(crit==1)name=Form("%s/AverageEperHit_%s.gif", dirOut.Data(), (const char*)HistoName);
	if(crit==2)name=Form("%s/AverageHitperEvent_%s.gif", dirOut.Data(), (const char*)HistoName);
	c1->SaveAs(name);


	//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	//. . . Mark the bad cells in the pflag array
	//. . .(0= bad because cell average value lower than min allowed)
	//. . .(2= bad because cell average value higher than max allowed)
	//. . .(1 by default)
	//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	cout<<"    o Flag bad cells that are outside the good range "<<endl;
	Int_t cel;
	//..*pflag[1][0] stores the number of cells
	for(Int_t c = 1; c <= *pflag[1][0]; c++)
	{
		cel=(Int_t)(inhisto->GetBinLowEdge(c)+0.1);  //ELI what does that 0.1 stand for?
		//cel=0 and c=1, cel=1 and c=2
		if (inhisto->GetBinContent(c) <= goodmin)
		{
			*pflag[cel][crit]=0;
		}
		else if (inhisto->GetBinContent(c) > goodmax)
		{
			*pflag[cel][crit]=2;
		}
	}
	cout<<"    o "<<endl;

}
//_________________________________________________________________________
//_________________________________________________________________________

void TestCellEandN(Int_t *pflag[23040][7], Double_t Emin = 0.1, Double_t Emax=2., Double_t Nsigma = 4.)
{
	//..here the average hit per event and the average energy per hit is caluclated for each cell.
	Int_t dnbins = 200;
	TH2 *hCellAmplitude = (TH2*) gFile->Get("hCellAmplitude");

	//..binning parameters
	Int_t ncells  = hCellAmplitude->GetNbinsY();
	Double_t amin = hCellAmplitude->GetYaxis()->GetXmin();
	Double_t amax = hCellAmplitude->GetYaxis()->GetXmax();
	cout<<"    o Calculate average cell hit per event and average cell energy per hit "<<endl;

	TH1F *hCellNtotal = new TH1F(Form("hCellNtotal_E%.2f-%.2f",Emin,Emax),Form("Number of hits per events, %.2f < E < %.2f GeV",Emin,Emax), ncells,amin,amax);
	hCellNtotal->SetXTitle("AbsId");
	hCellNtotal->SetYTitle("Av. hits per events");
	hCellNtotal->GetXaxis()->SetNdivisions(505);

	TH1F *hCellEtoNtotal = new TH1F(Form("hCellEtoNtotal_E%.2f-%.2f",Emin,Emax),Form("Average energy per hit, %.2f < E < %.2f GeV",Emin,Emax), ncells,amin,amax);
	hCellEtoNtotal->SetXTitle("AbsId");
	hCellEtoNtotal->SetYTitle("Av. energy per hit, GeV");
	hCellEtoNtotal->GetXaxis()->SetNdivisions(505);

	TH1* hNEventsProcessedPerRun = (TH1*) gFile->Get("hNEventsProcessedPerRun");
	Double_t totalevents = hNEventsProcessedPerRun->Integral(1, hNEventsProcessedPerRun->GetNbinsX());

	//..here the average hit per event and the average energy per hit is caluclated for each cell.
	for (Int_t c = 1; c <= ncells; c++)
	{
		Double_t Esum = 0;
		Double_t Nsum = 0;

		for (Int_t j = 1; j <= hCellAmplitude->GetNbinsX(); j++)
		{
			Double_t E = hCellAmplitude->GetXaxis()->GetBinCenter(j);
			Double_t N = hCellAmplitude->GetBinContent(j, c);
			if (E < Emin || E > Emax) continue;
			Esum += E*N;
			Nsum += N;
		}
		if(totalevents> 0.)hCellNtotal   ->SetBinContent(c, Nsum/totalevents);  //..number of hits per event
		if(Nsum > 0.)      hCellEtoNtotal->SetBinContent(c, Esum/Nsum);         //..average energy per hit
		//ELI maybe plot 2-dimensional hit/event eta vs. phi??
	}
	delete hCellAmplitude;

	if(*pflag[0][0]==1) Process(pflag,hCellEtoNtotal,Nsigma,dnbins,-1);
	if(*pflag[0][0]==2 && Emin==0.5) Process(pflag,hCellNtotal,   Nsigma,dnbins*9000,-1);//ELI I did massivley increase the binning now but it helps a lot
	if(*pflag[0][0]==2 && Emin>0.5)  Process(pflag,hCellNtotal,   Nsigma,dnbins*17,-1);
}

//_________________________________________________________________________
//_________________________________________________________________________

void TestCellShapes(Int_t *pflag[23040][7], Double_t fitEmin, Double_t fitEmax, Double_t Nsigma =4.)
{
	//ELI this method is currently not used
	// Test cells shape using fit function f(x)=A*exp(-B*x)/x^2.
	// Produce values per cell + distributions for A,B and chi2/ndf parameters.

	TString hname= "hCellAmplitude";
	Int_t dnbins = 1000;
	TH2 *hCellAmplitude = (TH2*) gFile->Get(Form("%s",(const char*)hname));

	// binning parameters
	Int_t  ncells = hCellAmplitude->GetNbinsY();
	Double_t amin = hCellAmplitude->GetYaxis()->GetXmin();
	Double_t amax = hCellAmplitude->GetYaxis()->GetXmax();
	cout << "ncells " << ncells << " amin = " << amin << "amax = " << amax<< endl;

	// initialize histograms
	TH1 *hFitA = new TH1F(Form("hFitA_%s",(const char*)hname),"Fit A value", ncells,amin,amax);
	hFitA->SetXTitle("AbsId");
	hFitA->SetYTitle("A");

	TH1 *hFitB = new TH1F(Form("hFitB_%s",(const char*)hname),"Fit B value", ncells,amin,amax);
	hFitB->SetXTitle("AbsId");
	hFitB->SetYTitle("B");

	TH1 *hFitChi2Ndf = new TH1F(Form("hFitChi2Ndf_%s",(const char*)hname),"Fit #chi^{2}/ndf value", ncells,amin,amax);
	hFitChi2Ndf->SetXTitle("AbsId");
	hFitChi2Ndf->SetYTitle("#chi^{2}/ndf");

	Double_t maxval1=0., maxval2=0., maxval3=0.;
	Double_t prev=0., MSA=0., AvA = 0. ; //those param are used to automaticaly determined a reasonable maxval1
	Double_t prev2=0., MSB=0., AvB = 0.  ; //those param are used to automaticaly determined a reasonable maxval2
	Double_t prev3=0., MSki2=0., Avki2 = 0. ; //those param are used to automaticaly determined a reasonable maxval3
	Double_t ki2=0.0 ;
	for (Int_t k = 1; k <= ncells; k++)
	{
		TF1 *fit = new TF1("fit", "[0]*exp(-[1]*x)/x^2");
		TH1 *hCell = hCellAmplitude->ProjectionX("",k,k);
		if (hCell->GetEntries() == 0) continue;
		// hCell->Rebin(3);
		hCell->Fit(fit, "0QEM", "", fitEmin, fitEmax);
		delete hCell;

		if(fit->GetParameter(0) < 5000.)
		{
			hFitA->SetBinContent(k, fit->GetParameter(0));
			if(k<3000)
			{
				AvA +=  fit->GetParameter(0);
				if(k==2999)  maxval1  = AvA/3000. ;
				if (prev < fit->GetParameter(0)) MSA += fit->GetParameter(0) - prev;
				else MSA -= (fit->GetParameter(0) - prev) ;
				prev = fit->GetParameter(0);
			}
			else
			{
				if((fit->GetParameter(0) - maxval1) > 0. && (fit->GetParameter(0) - maxval1) < (MSA/1000.))
				{
					maxval1 = fit->GetParameter(0);
				}
			}
		}
		else hFitA->SetBinContent(k, 5000.);

		if(fit->GetParameter(1) < 5000.)
		{
			hFitB->SetBinContent(k, fit->GetParameter(1));
			if(k<3000)
			{
				AvB +=  fit->GetParameter(1);
				if(k==2999)  maxval2  = AvB/3000. ;
				if (prev2 < fit->GetParameter(1)) MSB += fit->GetParameter(1) - prev2;
				else MSB -= (fit->GetParameter(1) - prev2) ;
				prev2 = fit->GetParameter(1);
			}
			else
			{
				if((fit->GetParameter(1) - maxval2) > 0. && (fit->GetParameter(1) - maxval2) < (MSB/1000.))
				{
					maxval2 = fit->GetParameter(1);
				}
			}
		}
		else hFitB->SetBinContent(k, 5000.);


		if (fit->GetNDF() != 0 ) ki2 =  fit->GetChisquare()/fit->GetNDF();
		else ki2 = 1000.;

		if(ki2 < 1000.)
		{
			hFitChi2Ndf->SetBinContent(k, ki2);
			if(k<3000)
			{
				Avki2 +=  ki2;
				if(k==2999)  maxval3  = Avki2/3000. ;
				if (prev3 < ki2) MSki2 += ki2 - prev3;
				else MSki2 -= (ki2 - prev3) ;
				prev3 = ki2;
			}
			else
			{
				if((ki2 - maxval3) > 0. && (ki2 - maxval3) < (MSki2/1000.))
				{
					maxval3 = ki2;
				}
			}
		}
		else hFitChi2Ndf->SetBinContent(k, 1000.);

		delete fit ;
	}

	delete hCellAmplitude;

	// if you have problem with automatic parameter :
	//  maxval1 =
	//  maxval2 =
	//  maxval3 =
	if(*pflag[0][0]==3)
		Process(pflag, hFitChi2Ndf, Nsigma, dnbins, maxval3);
	if(*pflag[0][0]==4)
		Process(pflag, hFitA, Nsigma, dnbins,  maxval1);
	if(*pflag[0][0]==5)
		Process(pflag, hFitB, Nsigma, dnbins, maxval2);
}
//_________________________________________________________________________
//_________________________________________________________________________
void ExcludeCells(Int_t *pexclu[23040], Int_t NrCells)
{
	//..This function finds cells with zero entries
	//..to exclude them from the analysis
	//cout<<"In ExcludeCells: Name of current file: "<<gFile->GetName()<<endl;
	TH2 *hCellAmplitude = (TH2*) gFile->Get("hCellAmplitude");

	Int_t SumOfExcl=0;

	//..Direction of cell ID
    for (Int_t c = 1; c <= NrCells; c++)
	{
		Double_t Nsum = 0;
		//..Direction of amplitude
		for (Int_t amp = 1; amp <= hCellAmplitude->GetNbinsX(); amp++)
		{
			Double_t N = hCellAmplitude->GetBinContent(amp,c);
			Nsum += N;
		}
		//..If the amplitude in one cell is basically 0
		//..mark the cell as excluded
		//ELI I just wonder how you can have less than one but more than 0.5
		//shouldnt everything below 1 be excluded?
		if(Nsum >= 0.5 && Nsum < 1)cout<<"-----------------------small but non zero!!!!"<<endl;
		if(Nsum < 0.5 && Nsum != 0)cout<<"-----------------------non zero!!!!"<<endl;

		if(Nsum < 0.5 && *pexclu[c-1]!=5)
		{
			//..histogram bin=cellID+1
			*pexclu[c-1]=1;
			SumOfExcl++;
		}
		else *pexclu[c-1]=0;
	}
	delete hCellAmplitude;
	cout<<"    o Number of excluded cells: "<<SumOfExcl<<endl;
	cout<<"     ("<<SumOfExcl<<")"<<endl;
}

//_________________________________________________________________________
//_________________________________________________________________________
void KillCells(Int_t filter[], Int_t nbc, TString outPath = "ConvertOutput")
{
	cout<<"    o Kill cells -> set bin content of "<<nbc<<" bad cells to 0 "<<endl;
	// kill a cell : put its entry to 0
	TH2 *hCellAmplitude          = (TH2*) gFile->Get("hCellAmplitude");
	TH1* hNEventsProcessedPerRun = (TH1*) gFile->Get("hNEventsProcessedPerRun");

	//..loop over number of identified bad cells. ID is stored in filer[] array
	for(Int_t i =0; i<nbc; i++)
	{
		if(filter[i]==-1)cout<<"#### That is strange - shouln't happen"<<endl;
		//..set all amplitudes for a given cell to 0
		for(Int_t amp=0; amp<= hCellAmplitude->GetNbinsX() ;amp++)
		{
			                        //(amplitiude,cellID,new value)
			//..CellID=0 is stored in bin1 so we need to shift+1
			hCellAmplitude->SetBinContent(amp,filter[i]+1,0);
		}
	}
	gSystem->mkdir(outPath);
	TFile *tf = new TFile(Form("%s/filter.root", outPath.Data()),"recreate");
	hCellAmplitude->Write();
	hNEventsProcessedPerRun->Write();
	tf->Write();
	tf->Close();
	delete hCellAmplitude;
	delete hNEventsProcessedPerRun;
}
//_________________________________________________________________________
//_________________________________________________________________________
void PeriodAnalysis(Int_t criterum=7, Double_t Nsigma = 4.0, Double_t Emin=0.1, Double_t Emax=2.0, TString period = "LHC15f", TString pass = "pass2", Int_t trial=0, TString Infilefile ="none", TString outPath = "ConvertOutput")
{
	cout<<""<<endl;
	cout<<""<<endl;
	cout<<""<<endl;
	cout<<"o o o o o o o o o o o o o o o o o o o o o o  o o o"<<endl;
	cout<<"o o o PeriodAnalysis for flag "<<criterum<<" o o o"<<endl;
	cout<<"o o o Done in the energy range E "<<Emin<<"-"<<Emax<<endl;
	TH2 *hCellAmplitude = (TH2*) gFile->Get("hCellAmplitude");

    //..This function does perform different checks depending on the given criterium variable
	//..diffrent possibilities for criterium are:
	// 1 : average E for E>Emin
	// 2 : entries for E>Emin
	// 3 : ki²/ndf  (from fit of each cell Amplitude between Emin and Emax)
	// 4 : A parameter (from fit of each cell Amplitude between Emin and Emax)
	// 5 : B parameter (from fit of each cell Amplitude between Emin and Emax)
	// 6 :
	// 7 : give bad + dead list
	//ELI Number of cells - (24*48)  16 full modules and 4 1/3 modules == 19,968 check that number!!
    static const Int_t NrCells=17663;//19968; //17665;//23040;  //ELI this is in fact a very important number!!

    //ELI a comment about the array positions
    //..In the histogram: bin 1= cellID 0, bin 2= cellID 1 etc
    //..In the arrays: array[cellID]= some information
    Int_t newBC[NrCells];       // starts at newBC[0] stores cellIDs  (cellID = bin-1)
    Int_t newDC[NrCells];       // starts at newDC[0] stores cellIDs  (cellID = bin-1)
	Int_t *pexclu[NrCells];     // starts at 0 pexclu[CellID] stores 0 not excluded, 1 excluded
	Int_t exclu[NrCells];       // is the same as above
	Int_t *pflag[NrCells][7];   // pflag[cellID][crit] = 1(ok),2(bad),0(bad)     start at 0 (cellID 0 = histobin 1)
	Int_t flag[NrCells][7];     // is the same as above
	//..set all fields to -1
	memset(newBC,-1, NrCells *sizeof(int));
	memset(newDC,-1, NrCells *sizeof(int));

	Int_t CellID, nb1=0, nb2=0, nbcemc = 0, nbcdca = 0, ndcemc = 0, ndcdca = 0;
	//INIT
	TString output, bilan, DeadPdfName, BadPdfName;
	for(CellID=0;CellID<NrCells;CellID++)
	{
		exclu[CellID] =0;
		pexclu[CellID]=&exclu[CellID];
		for(Int_t j=0;j<7;j++)
		{
			flag[CellID][j] =1;
			pflag[CellID][j]=&flag[CellID][j];
		}
	}
	//..Side note [x][y=0] is never used as there is no criterium 0
	//..use the [0][0] to store the criterum that is tested
	flag[0][0]=criterum ;
	//..use the [0][1] to store the number of cells tested
	flag[1][0]=NrCells;

	//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	//.. CELLS EXCLUDED
	//.. exclude cells from analysis (will not appear in results)
	//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	cout<<"o o o Exclude Cells o o o"<<endl;
	ExcludeCells(pexclu,NrCells);


	//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	//.. ANALYSIS
	//.. Build average distributions and fit them
	//.. Three tests for bad cells:
	//.. 1) Average energy per hit;
	//.. 2) Average hit per event;
	//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	if(criterum < 6)cout<<"o o o Analyze average cell distributions o o o"<<endl;
	//..For case 1 or 2
	if (criterum < 3)      TestCellEandN(pflag, Emin, Emax,Nsigma);
	//..For case 3, 4 or 5
	else if (criterum < 6) TestCellShapes(pflag, Emin, Emax, Nsigma);


	//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	//.. RESULTS
	//.. 1) Print the bad cells
	//..    and write the results to a file
	//.. 2) Kill cells function...
	//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	if(criterum < 6)
	{
		//..Print the results on the screen and
		//..write the results in a file
		output.Form("%s/Criterion%d_Emin-%.2f_Emax-%.2f.txt", outPath.Data(), criterum,Emin,Emax);
		ofstream file(output, ios::out | ios::trunc);
		if(!file)
		{
			cout<<"#### Major Error. Check the textfile!"<<endl;
		}
		file<<"Criterion : "<<criterum<<", Emin = "<<Emin<<" GeV"<<", Emax = "<<Emax<<" GeV"<<endl;
		file<<"Bad by lower value : "<<endl;
		cout<<"    o bad cells by lower value (for cell E between "<<Emin<<"-"<<Emax<<")"<<endl;
		cout<<"      ";
		nb1=0;
		for(CellID=0;CellID<NrCells;CellID++)
		{
			if(flag[CellID][criterum]==0 && exclu[CellID]==0)
			{
				newBC[nb1]=CellID;
				nb1++;
				file<<CellID<<", ";
				cout<<CellID<<",";
			}
		}
		file<<"("<<nb1<<")"<<endl;
		cout<<"("<<nb1<<")"<<endl;
		file<<"Bad by higher value : "<<endl;
		cout<<"    o bad cells by higher value (for cell E between "<<Emin<<"-"<<Emax<<")"<<endl;
		cout<<"      ";
		nb2=0;
		for(CellID=0;CellID<NrCells;CellID++)
		{
			if(flag[CellID][criterum]==2 && exclu[CellID]==0)
			{
				newBC[nb1+nb2]=CellID;
				nb2++;
				file<<CellID<<", ";
				cout<<CellID<<",";
			}
		}
		file<<"("<<nb2<<")"<<endl;
		cout<<"("<<nb2<<")"<<endl;

		file<<"Total number of bad cells"<<endl;
		file<<"("<<nb1+nb2<<")"<<endl;
		file.close();
		cout<<"    o Total number of bad cells "<<endl;
		cout<<"      ("<<nb1+nb2<<")"<<endl;

		//..create a filtered file
		KillCells(newBC,nb1+nb2) ;
	}

	//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	//.. CRITERUM 7 : FINAL RESULT
	//.. 1) summarize all dead and bad cells in a text file
	//.. 2) plot all bad cell E distributions in a .pdf file
	//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	if(criterum ==7)
	{

		DeadPdfName = Form("%s/%s%sDC_SummaryResults_V%i.pdf", outPath.Data(), period.Data(), pass.Data(), trial);
	    BadPdfName  = Form("%s/%s%sBC_SummaryResults_V%i.pdf", outPath.Data(), period.Data(), pass.Data(), trial);
		bilan   = Form("%s/%s%sBC_SummaryResults_V%i.txt", outPath.Data(), period.Data(), pass.Data(), trial); ;
		cout<<"    o Final results o "<<endl;
		cout<<"    o write results into .txt file: "<<bilan<<endl;
		cout<<"    o write results into .pdf file: "<<BadPdfName<<endl;
		ofstream file(bilan, ios::out | ios::trunc);
		if(file)
		{
			file<<"Dead cells : "<<endl;
			cout<<"    o Dead cells : "<<endl;
			
			nb1 =0;
			for(CellID=0; CellID<NrCells; CellID++)
			{
				if(exclu[CellID]==1)
				{
					newDC[nb1]=CellID;
					file<<CellID<<"\n" ;
					cout<<CellID<<"," ;
					exclu[CellID]=5;
					nb1++;
					
					// count how many belong to emcal or dcal
					if(CellID < 12288) {
						nbcemc++;
						ndcemc++;
					}
					else {
						nbcdca++;
						ndcdca++;
					}
				}
			}
			file<<"("<<nb1<<")"<<endl;
			cout<<"("<<nb1<<")"<<endl;

			TFile::Open(Form("%s/filter.root", outPath.Data()));
			ExcludeCells(pexclu,NrCells);
			file<<"Bad cells (excluded): "<<endl;
			cout<<"    o Total number of bad cells (excluded): "<<endl;
			nb2=0;
			for(CellID=0;CellID<NrCells;CellID++)
			{
				if(exclu[CellID]==1)
				{
					newBC[nb2]=CellID;
					file<<CellID<<"\n" ;
					cout<<CellID<<"," ;
					nb2++;
					// count how many belong to emcal or dcal
					if(CellID < 12288) nbcemc++;
					else nbcdca++;
				}
			}
			file<<"("<<nb2<<")"<<endl;
			cout<<"("<<nb2<<")"<<endl;
			
			file<<"Tot in EMCal = "<<nbcemc<<" (dead = "<<ndcemc<<")"<<endl;
			file<<"Tot in DCal = "<< nbcdca<<" (dead = "<<ndcdca<<")"<<endl;
		}
		file.close();

		if(Infilefile!="none")
		{
			cout<<"    o Open original file: "<<Infilefile<<endl;
			TFile::Open(Infilefile);
			Int_t c;
			//..loop over the bad cells in packages of 9
			/*cout<<"    o Save the Dead channel spectra to a .pdf file"<<endl;
			for(Int_t w=0; (w*9)<nb1; w++)
			{
				if(9<=(nb1-w*9)) c = 9 ;
				else c = nb1-9*w ;
				SaveBadCellsToPDF(newDC, w*9, c,DeadPdfName);
			}*/
			cout<<"    o Save the bad channel spectra to a .pdf file"<<endl;
			/*for(Int_t w=0; (w*9)<nb2; w++)
			//for(Int_t w=0; (w*9)<10; w++)
			{
				if(9<=(nb2-w*9)) c = 9 ;
				else c = nb2-9*w ;
				SaveBadCellsToPDF(newBC, w*9, c,BadPdfName) ;
			}*/
		}
	}

}
//_________________________________________________________________________
//_________________________________________________________________________

void BCAnalysis(TString file, TString trigger = "default",TString period = "LHC15f", TString pass = "pass2",Int_t trial = 0){

	//..Configure a complete analysis with different criteria, it provides bad+dead cells lists
	//..You can manage criteria used and their order, the first criteria will use the original
	//..output file from AliAnalysisTaskCaloCellsQA task, then after each criteria it will use a
	//..filtered file without the badchannel previously identified
    cout<<"o o o Bad channel analysis o o o"<<endl;

	Int_t criter;
	Double_t  Emini, Emaxi, Nsig;
    //..Default Configuration:
	if(trigger=="default"||trigger=="INT7"||trigger=="DMC7"||trigger=="AnyINTnoBC")
	{
		TFile *fin = new TFile(file);
		if(!fin->IsOpen()){
			Printf("File %s not found", file.Data());
			return;
		}
		PeriodAnalysis(2, 4.,0.2,0.5,period,pass,trial); // nb ent emin emax
		TFile::Open("ConvertOutput/filter.root");
		PeriodAnalysis(2, 4.,0.5, 1.,period,pass,trial); // nb ent emin emax
		TFile::Open("ConvertOutput/filter.root");
		PeriodAnalysis(1, 6.,0.5, 1.,period,pass,trial); // energy mea emin emax
		TFile::Open("ConvertOutput/filter.root");
		PeriodAnalysis(2, 4., 1., 2.,period,pass,trial); // nb ent emin emax
		TFile::Open("ConvertOutput/filter.root");
		PeriodAnalysis(1, 6., 1., 2.,period,pass,trial); // energy mea emin emax
		TFile::Open("ConvertOutput/filter.root");
		PeriodAnalysis(2, 4., 1.,10.,period,pass,trial); //nb ent emin emax
		TFile::Open("ConvertOutput/filter.root");
		PeriodAnalysis(1, 6., 1.,10.,period,pass,trial); //energy mea emin emax
	}
	else
	{
		//..you have the possibility to change analysis configuration  in function of trigger type
		//..PeriodAnalysis(Int_t criterum=7, Double_t Nsigma = 4.0, Double_t Emin=0.1, Double_t Emax=2.0, Int_t compteur = 1, TString period = "LHC15f", TString pass = "pass2", Int_t trial=0, TString Infilefile ="none"){
		//..Criterium 1,2
		//..low energies 0.5-2
		TFile *fin = new TFile(file);
		if(!fin->IsOpen()){
			Printf("File %s not found", file.Data());
			return;
		}
		PeriodAnalysis(2, 6.,0.5, 2.,period,pass,trial); // nb ent emin emax
		TFile::Open("ConvertOutput/filter.root");
		PeriodAnalysis(1, 6.,0.5, 2.,period,pass,trial); // energy mea emin emax
		//..mid energies
        TFile::Open("ConvertOutput/filter.root");
		PeriodAnalysis(2, 6., 2., 5.,period,pass,trial); // nb ent emin emax
		TFile::Open("ConvertOutput/filter.root");
		PeriodAnalysis(1, 6., 2., 5.,period,pass,trial); // energy mea emin emax
//		TFile::Open("ConvertOutput/filter.root");
//		PeriodAnalysis(1, 6., 2., 5.01,period,pass,trial); // test test test

		//..high energies
		//ELI this is not working properly
/*		TFile::Open("ConvertOutput/filter.root");
		PeriodAnalysis(2, 6., 5.,10.,period,pass,trial); // nb ent emin emax
		TFile::Open("ConvertOutput/filter.root");
		PeriodAnalysis(1, 6., 5.,10.,period,pass,trial); // energy mea emin emax
*/
	}
	//provide dead cells list from original file and draw bad cells candidate from indicated file
	TFile::Open(file);
	PeriodAnalysis(7,0.,0.,0.,period,pass,trial,file);

	cout<<"o o o End of bad channel analysis o o o"<<endl;
}
//_________________________________________________________________________
//________________________________________________________________________
//void BadChannelAnalysis()
//{
//	cout<<"Error needs input arguments! Try:"<<endl;
//	cout<<".x BadChannelAnalysis.C(''LHC16h'',''muon_caloLego'',''AnyINT'')"<<endl;
//}
//_________________________________________________________________________
//________________________________________________________________________
void BadChannelAnalysis(TString period = "LHC15f", TString pass = "pass2", TString trigger= "default", TString workingDir = "./", TString externalFile= "", Int_t trial=0)
{
	///..Define externalFile if the Conversion has already been performed, namely if a file named <period><pass>Converted.root already exists, by default in the directory ConvertOutput. externalFile must contain the full path and filename.
	///..First use Convert() to merge historgrams from a runlist .txt file
    ///..The merged outputfile contains 3 different histograms
	///..In a second step analyse these merged histograms
	///..by calling BCAnalysis()
	TString inputfile = "";

	if(externalFile=="")
	{
		cout<<". . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ."<<endl;
		cout<<". . .Start process by converting files. . . . . . . . . . . ."<<endl;
		cout<<endl;
		inputfile = Convert(period, pass, trigger, workingDir);
		cout<<endl;
	}
	else
	{
		//inputfile="ConvertOutput/";
		inputfile+=externalFile;
	}

	cout<<". . .Load inputfile with name: "<<inputfile<<" . . . . . . . ."<<endl;
	cout<<". . .Continue process by . . . . . . . . . . . ."<<endl;
	//inputfile="LHC15omuon_caloLegoRunlist0New.root";
	cout<<endl;
	BCAnalysis(inputfile,trigger,period,pass,trial);
	cout<<endl;
	cout<<". . .End of process . . . . . . . . . . . . . . . . . . . . ."<<endl;
	cout<<". . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ."<<endl;
}

