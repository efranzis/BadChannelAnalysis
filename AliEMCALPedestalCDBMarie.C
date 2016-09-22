// Script to create dead channel map and store them into CDB
//  - 4 sets of maps parameters can be created, with now dead channels 
// and 5%, 10%, 20% and 30% of dead channels 
//  - it reads the stored map in a given file. 
// Author: Gustavo Conesa

//.x $ALICE_ROOT/EMCAL/macros/CalibrationDB/AliEMCALSetTowerStatusCDB.C

#if !defined(__CINT__)
#include <TControlBar.h>
#include <TString.h>
#include <TRandom.h>
#include <TStyle.h>
#include <TH2.h>
#include <TF1.h>
#include <TCanvas.h>

#include "AliRun.h"
#include "AliCaloCalibPedestal.h"
#include "AliEMCALGeoParams.h"
#include "AliCDBMetaData.h"
#include "AliCDBId.h"
#include "AliCDBEntry.h"
#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#endif


void AliEMCALPedestalCDBMarie()
{
	TControlBar *menu = new TControlBar("vertical","EMCAL CDB");
	menu->AddButton("Help to run EMCAL CDB","Help()","Explains how to use EMCAL CDS menus");
	menu->AddButton("Equal Tower Status Map, all Alive","SetTowerStatusMap(0)","Set all channels to alive");
	menu->AddButton("Create Random Status Map, 5% dead","SetTowerStatusMap(5)","Set randomly 5% of the channels dead");
	menu->AddButton("Create Random Status Map, 10% dead","SetTowerStatusMap(10)","Set randomly 10% of the channels dead");
	menu->AddButton("Create Random Status Map, 20% dead","SetTowerStatusMap(20)","Set randomly 20% of the channels dead");
	menu->AddButton("Create Random Status Map, 30% dead","SetTowerStatusMap(30)","Set randomly 30% of the channels dead");
	menu->AddButton("Set Map from txt file","SetTowerStatusMap(\"map.txt\")","Read bad channels from txt file and set them in root file");
	menu->AddButton("Read Tower Status Map","GetTowerStatusMap()","Read initial equal calibration coefficients");
	menu->Show();
}

//------------------------------------------------------------------------
void Help()
{
	char *string =
			"\nSet tower status map (dead, hot, alive) and write them into ALICE CDB. Press button \"Equal kAlive\" to set all channels alive. Press button \"Random, 5% dead\" to create random dead channel map at 5%\n";
	printf(string);
}

//------------------------------------------------------------------------
void SetTowerStatusMap(Int_t percent=0)
{
	// Writing status of all the channels in the OCDB with equal value
	// except a percent to be not alive. Right now only "alive" or "dead",
	// we need to implement the other cases like "hot"

	TString sDBFolder ="local://PedestalsDB/2015";
	//  TString sDBFolder ="local://PedestalsDB/2011";
	//  TString sDBFolder ="local:///scratch/alicehp2/germain/PedestalsDB";
	Int_t firstRun   =  0;
	Int_t lastRun    =  999999999;
	Int_t beamPeriod =  1;
	char* objFormat = Form("%d percent of bad channels", percent);

	AliCaloCalibPedestal *caloped=new AliCaloCalibPedestal(AliCaloCalibPedestal::kEmCal);
	caloped->Init();

	TObjArray map = caloped->GetDeadMap();
	printf("MAP entries %d\n",map.GetEntries());

	TRandom rn;
	//for(Int_t iSM = 0; iSM < AliEMCALGeoParams::fgkEMCALModules; iSM ++){
	for(Int_t iSM = 0; iSM < map.GetEntries(); iSM ++)
	{
		Int_t ndead = 0;
		printf(" >>> SM %d <<< Entries %d, NbinsX %d, NbinsY %d\n",iSM,((TH2D*)map[iSM])->GetEntries(),((TH2D*)map[iSM])->GetNbinsX(),((TH2D*)map[iSM])->GetNbinsY());
		for(Int_t i = 0; i < ((TH2D*)map[iSM])->GetNbinsX() ; i++)
		{
			for(Int_t j = 0; j < ((TH2D*)map[iSM])->GetNbinsY() ; j++)
			{
				//printf("Bin (%d-%d) Content, before: %d ",i,j,((TH2D*)map[iSM])->GetBinContent(i, j));	

				if(rn.Uniform(0,100) > percent)
					caloped->SetChannelStatus(iSM, i, j, AliCaloCalibPedestal::kAlive);
				else
				{
					caloped->SetChannelStatus(iSM, i, j, AliCaloCalibPedestal::kDead);
					ndead++;
				}
				//printf("; after: %d \n",((TH2D*)map[iSM])->GetBinContent(i, j));	
			}	
		}
		caloped->SetDeadTowerCount(caloped->GetDeadTowerCount()+ndead);
		printf("--- dead %d\n",ndead);
	}

	printf("--- total dead %d\n",caloped->GetDeadTowerCount());

	//Store map into database

	AliCDBMetaData md;
	md.SetComment(objFormat);
	md.SetBeamPeriod(beamPeriod);
	md.SetResponsible("Gustavo Conesa");

	AliCDBId id("EMCAL/Calib/Pedestals",firstRun,lastRun); // create in EMCAL/Calib/Pedestal sDBFolder

	AliCDBManager* man = AliCDBManager::Instance();
	AliCDBStorage* loc = man->GetStorage(sDBFolder.Data());
	loc->Put(caloped, id, &md);

}

//____________________________________________




//------------------------------------------------------------------------
void SetTowerStatusMap(char * file = "map.txt",Int_t firstRun = 0, Int_t lastRun =999999)
//void SetTowerStatusMap(char * file = "LHC11c153570-154733.txt")
{
	// Get the list of dead/hot channels from file and set them in OCDB

	TString sDBFolder ="local://PedestalsDB/2015";
	//TString sDBFolder ="local://PedestalsDB/2013";
	//  TString sDBFolder ="local:///scratch/alicehp2/germain/PedestalsDB";
	// Int_t firstRun   =  187534;
	// Int_t lastRun    =  187562;
	Int_t beamPeriod =  1;
	char* objFormat = Form("bad channels extracted from file %s", file);
	cout << "run range"<< firstRun << "-"<< lastRun << endl;
	cout << "file : " << file << endl;

	AliCaloCalibPedestal *caloped=new AliCaloCalibPedestal(AliCaloCalibPedestal::kEmCal);
	caloped->Init();

	// Read parameter file line-by-line
	ifstream f;
	f.open(file);

	Int_t iabsId =-1, iSM=-1, icol=-1, irow=-1, istatus=-1, ndead=0 ;
	TString string;
	if (f.good())
	{
		while(string.ReadLine(f, kFALSE) && !f.eof())
		{
			sscanf(string.Data(), "%d %d %d %d %d",&iabsId,&iSM,&icol,&irow,&istatus);
			cout<<"absId= "<<iabsId<<",SM= "<<iSM<<", col= "<<icol<<", row= "<<irow<<", status="<<istatus<<endl;
			if(iSM==-1) continue;
			caloped->SetChannelStatus(iSM, icol, irow, istatus);
			ndead++;
		}
	}
	//  caloped->SetDeadTowerCount(ndead-2);
	//   printf("--- dead %d\n",ndead-2);
	caloped->SetDeadTowerCount(ndead);
	printf("--- dead %d\n",ndead);

	printf("--- total dead %d\n",caloped->GetDeadTowerCount());

	//Store map into database

	AliCDBMetaData md;
	md.SetComment(objFormat);
	md.SetBeamPeriod(beamPeriod);
	md.SetResponsible("Gustavo Conesa");

	AliCDBId id("EMCAL/Calib/Pedestals",firstRun,lastRun); // create in EMCAL/Calib/Pedestal sDBFolder

	AliCDBManager* man = AliCDBManager::Instance();
	AliCDBStorage* loc = man->GetStorage(sDBFolder.Data());
	loc->Put(caloped, id, &md);

}



//------------------------------------------------------------------------
void GetTowerStatusMap(Int_t runNumber=0)
{
	// Read status map

	// const char* geoType="EMCAL_COMPLETE12SMV1";
	const char* geoType="EMCAL_COMPLETE12SMv1_DCAL_8SM";
	//const char* geoType="EMCAL_COMPLETEV1";
	AliEMCALGeometry* geom = new AliEMCALGeometry(geoType,"EMCAL");
	Int_t nSupMod, nModule, nIphi, nIeta;
	Int_t iphi, ieta,kk;

	//cout<<"AliEMCALGeoUtils('"<<geoType<< ",'EMCAL') !!"<<endl;
	//cout <<" how many SMs?:"<<(geom->GetEMCGeometry())->GetNumberOfSuperModules()<< endl;

	Int_t absId, status;
	// TString sDBFolder ="local://PedestalsDB/2011";
	TString sDBFolder ="local://PedestalsDB/2016";
	//  TString sDBFolder ="local:///scratch/alicehp2/germain/PedestalsDB";
	//Int_t runNumber   =  148600;

	AliCaloCalibPedestal* caloped  = (AliCaloCalibPedestal*)
  			(AliCDBManager::Instance()
	->GetStorage(sDBFolder.Data())
	->Get("EMCAL/Calib/Pedestals",
			runNumber)->GetObject());
	//  gSystem->Load("libOADB");
	//TGrid::Connect("alien://");

	AliCDBManager* man = AliCDBManager::Instance();
	// man->SetDefaultStorage("raw://");
	//  man->SetDefaultStorage("local://PedestalsDB/2013/");
	man->SetDefaultStorage("local://PedestalsDB/2016");

	//  man->SetDefaultStorage("local://PedestalsDB/2011");
	man->SetRun(runNumber);
	AliCDBStorage *storage = man->GetDefaultStorage();
	AliCaloCalibPedestal* caloped  = (AliCaloCalibPedestal*) (storage->Get("EMCAL/Calib/Pedestals", runNumber)->GetObject());


	cout<<" run number" << runNumber <<endl;
	TObjArray map = caloped->GetDeadMap();
	printf("MAP entries %d\n",map.GetEntries());
	Int_t Totdead = 0;
	Int_t Totwarm = 0;
	Int_t Tothot = 0;
	Int_t TotTOT = 0;
	for(Int_t iSM = 0; iSM < 20; iSM ++)
	{ // flag -2 marie to limit  the 2 last SM
		//   for(Int_t iSM = 0; iSM < map.GetEntries()-2; iSM ++){ // flag -2 marie to limit  the 2 last SM
		TCanvas *cMap   = new TCanvas(Form("cMap%d",iSM),Form("SM %d dead map",iSM), 12,12,400,400);
		cMap->Divide(1,1);
		Int_t ndead = 0;
		Int_t ndeadVrai = 0;
		Int_t nhot = 0;
		Int_t nwarm = 0;
		Int_t absId;
		//	  printf(" >>> SM %d <<< Entries %d, NbinsX %d, NbinsY %d\n",iSM,((TH2D*)map[iSM])->GetEntries(),((TH2D*)map[iSM])->GetNbinsX(),((TH2D*)map[iSM])->GetNbinsY());
		for(Int_t i = 0; i < ((TH2D*)map[iSM])->GetNbinsX() ; i++)
		{
			for(Int_t j = 0; j < ((TH2D*)map[iSM])->GetNbinsY() ; j++)
			{
				if(((TH2D*)map[iSM])->GetBinContent(i, j)!=AliCaloCalibPedestal::kAlive)
					// //     if(((TH2D*)map[iSM])->GetBinContent(i, j)==AliCaloCalibPedestal::kHot)
					// 	{
					// 		printf("Bin (%d-%d) Content: %d \n",i,j,((TH2D*)map[iSM])->GetBinContent(i, j));
					// 		ndead++;}

					//if(((TH2D*)map[iSM])->GetBinContent(i, j)==AliCaloCalibPedestal::kDead ||
					// ((TH2D*)map[iSM])->GetBinContent(i, j)==AliCaloCalibPedestal::kHot)
				{
					//	printf("Bin (%d-%d) Content: %d \n",i,j,((TH2D*)map[iSM])->GetBinContent(i, j));
					absId = geom->GetAbsCellIdFromCellIndexes(iSM, j, i) ;

					printf("%d \t %d \t %d \t %d \t %d \n",absId, iSM,i,j,((TH2D*)map[iSM])->GetBinContent(i, j));

					ndead++;
				}
				if(((TH2D*)map[iSM])->GetBinContent(i, j)==AliCaloCalibPedestal::kHot){nhot++;}
				if(((TH2D*)map[iSM])->GetBinContent(i, j)==AliCaloCalibPedestal::kWarning){nwarm++;}
				if(((TH2D*)map[iSM])->GetBinContent(i, j)==AliCaloCalibPedestal::kDead){ndeadVrai++;}
			}
		}
		//	  printf("--- dead :%d warm : %d hot: %d  Tot: %d\n",ndeadVrai,nwarm,nhot,ndead);
		Totdead+=ndeadVrai;
		Tothot+=nhot;
		Totwarm+=nwarm;
		TotTOT+=ndead;
		cMap->cd(iSM);
		(TH2D*)map[iSM])->Draw("lego2");
	}

	printf("Total DEAD %d\n", caloped->GetDeadTowerCount());
	printf("Total DEAD Vrai counted %d Hot %d  warm %d Total %d  \n", Totdead, Tothot,Totwarm, TotTOT);

}
//------------------------------------------------------------------------
void PlotTowerStatusMap(Int_t runNumber=0)
{
	// Read status map

	gStyle->SetOptStat(0);
	TString sDBFolder ="local://PedestalsDB/2015";

	//  TString sDBFolder ="local:///scratch/alicehp2/germain/PedestalsDB";
	//Int_t runNumber   =  148600;

	AliCaloCalibPedestal* caloped  = (AliCaloCalibPedestal*)
			(AliCDBManager::Instance()
	->GetStorage(sDBFolder.Data())
	->Get("EMCAL/Calib/Pedestals",
			runNumber)->GetObject());

	cout<<" run nmber" << runNumber <<endl;
	TObjArray map = caloped->GetDeadMap();
	printf("MAP entries %d\n",map.GetEntries());
	TCanvas *cMap   = new TCanvas("cMap","dead map", 12,12,1000,1000);
	cMap->Divide(2,5,0.001,0.001);
	//for(Int_t iSM = 0; iSM < map.GetEntries(); iSM ++){
	for(Int_t iSM = 0; iSM < 20; iSM ++)
	{
		//         Form("cMap_%d",iSM)->cd();

		Int_t ndead = 0;
		printf(" >>> SM %d <<< Entries %d, NbinsX %d, NbinsY %d\n",iSM,((TH2D*)map[iSM])->GetEntries(),((TH2D*)map[iSM])->GetNbinsX(),((TH2D*)map[iSM])->GetNbinsY());
		for(Int_t i = 0; i < ((TH2D*)map[iSM])->GetNbinsX() ; i++)
		{
			for(Int_t j = 0; j < ((TH2D*)map[iSM])->GetNbinsY() ; j++)
			{
				if(((TH2D*)map[iSM])->GetBinContent(i, j)!=AliCaloCalibPedestal::kAlive)
					printf("Bin (%d-%d) Content: %d \n",i,j,((TH2D*)map[iSM])->GetBinContent(i, j));
				//	ndead++;}

				if(((TH2D*)map[iSM])->GetBinContent(i, j)==AliCaloCalibPedestal::kDead ||
						((TH2D*)map[iSM])->GetBinContent(i, j)==AliCaloCalibPedestal::kHot)
					ndead++
					;
			}
		}
		printf("--- dead %d\n",ndead);
		//  for(Int_t iSM = 0; iSM < map.GetEntries(); iSM ++){
		if (iSM==0)cMap_1->cd();
		if (iSM==1)cMap_2->cd();
		if (iSM==2)cMap_3->cd();
		if (iSM==3)cMap_4->cd();
		if (iSM==4)cMap_5->cd();
		if (iSM==5)cMap_6->cd();
		if (iSM==6)cMap_7->cd();
		if (iSM==7)cMap_8->cd();
		if (iSM==8)cMap_9->cd();
		if (iSM==9)cMap_10->cd();
		//   if (iSM==10)cMap_11->cd();
		//if (iSM==11)cMap_12->cd();
		(TH2D*)map[iSM])->GetXaxis()->SetTitle("col(eta)");
		(TH2D*)map[iSM])->GetYaxis()->SetTitle("row(phi)");
		(TH2D*)map[iSM])->GetXaxis()->SetLabelSize(0.09);
		(TH2D*)map[iSM])->GetXaxis()->SetTitleSize(0.06);
		(TH2D*)map[iSM])->GetYaxis()->SetTitleSize(0.06);
		(TH2D*)map[iSM])->GetXaxis()->SetTitleOffset(0.4);
		(TH2D*)map[iSM])->GetYaxis()->SetTitleOffset(0.5);
		(TH2D*)map[iSM])->GetYaxis()->SetTitleSize(0.06);
		(TH2D*)map[iSM])->GetYaxis()->SetLabelSize(0.08);
		(TH2D*)map[iSM])->SetMaximum(3);
		(TH2D*)map[iSM])->Draw("colz");

	}
	TString base = "/scratch/alicehp2/germain/NewPedestalsDB/2015/PedestalOCDBrun";
	base + =runNumber;
	TString outfile; outfile= base + ".gif";

	cMap->SaveAs(outfile);

	printf("Total DEAD %d\n", caloped->GetDeadTowerCount());

}
