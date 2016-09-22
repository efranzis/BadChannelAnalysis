
#if !defined(__CINT__) || defined(__MAKECINT__)
//Root include files
#include <Riostream.h>
#include <TFile.h>
#include <TChain.h>
#include <TParticle.h>
#include <TNtuple.h>
#include <TCanvas.h>
#include <TObjArray.h>
#include <TSystem.h>
#include <TString.h> 
#include <TH1F.h>
#include <TH1D.h> 
#include <TH2F.h> 
#include <TH2D.h> 
#include <TVector.h> 
#include <TParticle.h> 
#include <TRefArray.h>
#include <TArrayS.h>
#include <Riostream.h>
#include <TFile.h>
#include <TChain.h>
#include <TParticle.h>
#include <TNtuple.h>
#include <TCanvas.h>
#include <TObjArray.h>
#include <TSystem.h>
#include <TString.h> 
#include <TH1I.h> 
#include <TH1F.h> 
#include <TVector.h> 
#include <TParticle.h> 
#include <TRefArray.h>
#include <TArrayS.h>
#include <TObject.h> //ajouté
#include <TMath.h>  //add
#include "AliRunLoader.h"
#include "AliStack.h"
#include "AliESDEvent.h"
#include "AliESDVertex.h"
#include "AliESDCaloCluster.h"
#include "AliESDCaloCells.h"
#include "AliPID.h"
#include "AliLog.h" 
#include "AliEMCALPID.h"
#include "AliEMCALRecParam.h"
#include "AliEMCALReconstructor.h"
#include "AliCDBManager.h"
#include "AliEMCALRawUtils.h"
#include "AliCDBEntry.h"
#include "AliVertex.h"//add
#include "AliESDRun.h"
#include <AliEMCALGeometry.h> 
#endif

void etaphi(){ 
	Double_t eta, phi;
	Int_t  absId=5430;
	AliEMCALGeometry *geo = new AliEMCALGeometry() ;
	geo->EtaPhiFromIndex(absId, eta, phi);
	cout<<eta<<phi<<endl;
	
}
void writeBadChannelStatus(TString period = "LHC15f", TString train = "", TString trigger= "default", Int_t runNum= 254381,Int_t trial,Int_t istat,TString workDir=".")
//void writeBadChannelStatus(Int_t istat, TString pathinput = "LHC15omuon_calopass1BC_SummaryResults_V0.txt", const char* geoType="EMCAL_COMPLETE12SMV1_DCAL_8SM")
{
	TString analysisOutput  = Form("AnalysisOutput/%s/Version%i",period.Data(),trial);
	TString cellSummaryFile = Form("%s/%s/%s%s_Bad_Amplitudes_V%i.txt",fWorkdir.Data(), fAnalysisOutput.Data(), train.Data(), trigger.Data(),trial); ;

	ifstream fdata(pathinput);


	//const char* geoType="EMCAL_COMPLETE12SMV1";
	AliEMCALGeometry* geom = new AliEMCALGeometry(geoType,"EMCAL");
	//cout<<"AliEMCALGeoUtils('"<<geoType<< ",'EMCAL') !!"<<endl;
	//cout <<" how many SMs?:"<<(geom->GetEMCGeometry())->GetNumberOfSuperModules()<< endl;

	Int_t absId, status;

	Int_t p = -1;//, q,r;
	Int_t ncols;
	//Int_t nlines = 0 ;
	string line;

	while (fdata)
	{
		getline(fdata, line);
		ncols = atoi(line.c_str());
		Printf("Line = %s, number = %d", line.data(), ncols);
		// = fscanf(fdata,"%d",&p);
		if (ncols< 0) break;
		absId=p;

		// run=q;
		status=istat;

		Int_t nSupMod, nModule, nIphi, nIeta;
		Int_t iphi, ieta,kk;

		// kk is the AbsId.
		geom->GetCellIndex(absId,  nSupMod, nModule, nIphi, nIeta);
		geom->GetCellPhiEtaIndexInSModule(nSupMod,nModule,nIphi,nIeta, iphi,ieta);
		//cout << absId <<  nSupMod  << ieta << iphi << status <<endl;
		printf("%d \t %d \t %d \t %d \t %d \n",absId,nSupMod,ieta,iphi,status);
		// printf("%d \t %d %d %d %d \n",absId,nSupMod,ieta,iphi,status);
	}
}
