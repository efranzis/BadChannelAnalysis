
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

void writeBadChannelStatus(Int_t istat){


  //   fdata = fopen("/scratch/alicehp2/germain/PedestalDB/LHC11aLE-period2.txt","r");

  //    fdata = fopen("/scratch/alicehp2/germain/QA/LHC11c/LHC11crun153353-153560And154763-154796.txtWithComments","r");
  //  fdata = fopen("/scratch/alicehp2/germain/QA/LHC11c/LHC11crun153570-154480.txtWithComments","r");
  // fdata = fopen("/scratch/alicehp2/germain/QA/LHC11b/LHC11bFull4pass2.txt","r");
  // fdata = fopen("/scratch/alicehp2/germain/PedestalDB/LHC11dBadChannels.txt","r");

  //  fdata = fopen("/scratch/alicehp2/germain/PedestalDB/LHC12a176326-177295.txt","r");

  //   fdata = fopen("/scratch/alicehp2/germain/PedestalDB/BadCellCandidateLHC12ccpass1FullEMCAL.txt","r");
  //    fdata = fopen("/scratch/alicehp2/germain/PedestalDB/BadCellsLHC12efromcpass1.txt","r");
  //fdata = fopen("/scratch/alicehp2/germain/PedestalDB/LHC12dvpass1run185456to185784.txt","r");
  //fdata = fopen("/scratch/alicehp2/germain/PedestalDB/LHC12hcpass1.txt","r");
  //  fdata = fopen("/scratch/alicehp2/germain/PedestalDB/LHC11hFinal.txt","r");
  //    fdata = fopen("/scratch/alicehp2/germain/PedestalDB/BCLHC11a2.txt","r");      
  // fdata = fopen("/scratch/alicehp2/germain/PedestalDB/BC2013-3.txt","r");
  // fdata = fopen("/scratch/alicehp2/germain/PedestalDB/additionnalBC2013New.txt","r");
  // fdata = fopen("/scratch/alicehp2/germain/PedestalDB/dead197348-197388.txt","r");
  // fdata = fopen("/scratch/alicehp2/germain/PedestalDB/addBCLHC13b-fFromMarta.txt","r");
  //    fdata = fopen("/scratch/alicehp2/germain/PedestalDB/2012FromDaniel/new/LHC12badcellsDaniel/LHC12i-SortedByRun.log","r");
  //      fdata = fopen("/scratch/alicehp2/germain/PedestalDB/2012FromDaniel/deadLHC12i193750-193766.txt","r");
  //  fdata = fopen("/scratch/alicehp2/germain/QANew2/BadChannelNew/ResultsLHC15omuon_calopass1allrunsDEAD.txt","r");
  //  fdata = fopen("/scratch/alicehp2/germain/QANew2/BadChannelNew/LHC15oBadINT7.txt","r");
  fdata = fopen("/scratch/alicehp2/germain/QANew2/BadChannelNew/LHC15jMBDEAD.txt","r");
  //fdata = fopen("/scratch/alicehp2/germain/QANew2/BadChannelNew/ResultsLHC15iallruns4sigTestallenergyWARM.txt","r");
  // fdata = fopen("/scratch/alicehp2/germain/PedestalDB/2012FromDaniel/LHC12cspecialRuns182730-182744.txt","r");

const char* geoType="EMCAL_COMPLETE12SMV1_DCAL_8SM";
//const char* geoType="EMCAL_COMPLETE12SMV1";
AliEMCALGeometry* geom = new AliEMCALGeometry(geoType,"EMCAL");
//cout<<"AliEMCALGeoUtils('"<<geoType<< ",'EMCAL') !!"<<endl;
//cout <<" how many SMs?:"<<(geom->GetEMCGeometry())->GetNumberOfSuperModules()<< endl;

 Int_t absId, status;

 Int_t p, q,r;
  Int_t ncols;
  Int_t nlines = 0 ;


  while (1){
    ncols = fscanf(fdata,"%d",&p);
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
