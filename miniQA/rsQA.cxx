/*************************************************************************
* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
*                                                                        *
* Author: The ALICE Off-line Project.                                    *
* Contributors are mentioned in the code where appropriate.              *
*                                                                        *
* Permission to use, copy, modify and distribute this software and its   *
* documentation strictly for non-commercial purposes is hereby granted   *
* without fee, provided that the above copyright notice appears in all   *
* copies and that both the copyright notice and this permission notice   *
* appear in the supporting documentation. The authors make no claims     *
* about the suitability of this software for any purpose. It is          *
* provided "as is" without express or implied warranty.                  * 
**************************************************************************/


#include <TFile.h>
#include <TChain.h>
#include <TTree.h>
#include <TSystem.h>
#include <TROOT.h>
#include "AliAnalysisManager.h"
#include "AliESDInputHandler.h"
#include "rsQA.h"
#include "AliLog.h"
#include "AliMagF.h"
#include "AliMultiplicity.h"
#include "TGeoGlobalMagField.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliMCParticle.h"
#include "AliStack.h"
#include "AliITSgeomTGeo.h"
#include "TParticle.h"
#include "TParticlePDG.h"
#include "TStopwatch.h"
#include "AliPID.h"

ClassImp(rsQA)

const Double_t rsQA::kMaxZVtx=8., rsQA::kMaxDCAY=6.0,rsQA::kMaxDCAZ=4.0, rsQA::kMaxQ2Pt=5, rsQA::kMaxEta=1.0;
const Double_t rsQA::kPtRange[rsQA::kNPtRanges] = {0.,1.,2.};

//________________________________________________________________________
rsQA::rsQA(const char *name) 
: AliAnalysisTaskSE(name)
,fOutput(0)
,fTreeOut(0)
,fFileOut(0)
//
,fHistManArr(0)
,fOutFName(name)
,fHManMtcEta(0)
  //
,fBz(0)
,fUseMC(kFALSE)
,fMCEvent(0)
,fStack(0)
{
  // Constructor
  fHistManArr.SetOwner(kTRUE);
  memset(fHManDCA,0,kNCharges*kNPtRanges*sizeof(HistoManager*));
  memset(fHManMtcSect,0,kNSides*sizeof(HistoManager*));
  if (fOutFName.IsNull()) fOutFName = "rsQA.root";
  if (!fOutFName.EndsWith(".root")) fOutFName += ".root";

  DefineOutput(1, TList::Class());
  //
}

//________________________________________________________________________
rsQA::~rsQA()
{
  // Destructor
  if (fOutput && !AliAnalysisManager::GetAnalysisManager()->IsProofMode()) {  //RRR
    printf("Deleteing output\n");
    delete fOutput;
    fOutput = 0;
  }
  //
}

//________________________________________________________________________
void rsQA::UserCreateOutputObjects() 
{
  //
  const char *kQName[kNCharges] = {"P","N"};
  const char *kQTitle[kNCharges] = {"+","-"};
  const char *kSideName[kNSides] = {"A","C"};

  TH2* h2 = 0;
  TH1* h1 = 0;
  TString hnm, htt;
  double dEta = 2*kMaxEta/kNEtaRanges;
  for (int iEta=0;iEta<kNEtaRanges;iEta++) {
    TString hnmEb = Form("eta%d",iEta);
    TString httEb = Form("%+.1f<#eta<%+.1f",-kMaxEta+iEta*dEta,-kMaxEta+(iEta+1)*dEta);
    //
    for (int iQ=0;iQ<kNCharges;iQ++) {
      TString hnmQ = Form("q%s",kQName[iQ]);
      TString httQ = Form("q%s",kQTitle[iQ]);
      //
      for (int iPt=0;iPt<kNPtRanges;iPt++) {
	//
	TString hnmpt = Form("ptbin%d",iPt);
	TString httpt = Form("|q/p_{T}|>%.1f",kPtRange[iPt]);      
	//
	if (!fHManDCA[iQ][iPt]) {
	  fHManDCA[iQ][iPt] = new HistoManager(Form("dca_%s_%s",hnmQ.Data(),hnmpt.Data()), fOutFName.Data());
	  fHistManArr.Add(fHManDCA[iQ][iPt]);
	}
	//
	// DCA Y vs phi
	hnm = Form("%s_%s_%s_DCAY",hnmEb.Data(),hnmQ.Data(),hnmpt.Data());
	htt = Form("#DeltaY %s %s %s",hnmEb.Data(),hnmQ.Data(),hnmpt.Data());
	h2 = new TH2F(hnm.Data(),htt.Data(), kNBinsPerSect*18, 0, 18., kNBinsDCAY, -kMaxDCAY,kMaxDCAY);
	h2->SetXTitle("sector");
	h2->SetYTitle("#DeltaY");
	fHManDCA[iQ][iPt]->AddHisto(h2, GetHistoID(iEta,kHDCAY)*10);
	//
	// DCA Z vs phi
	hnm = Form("%s_%s_%s_DCAZ",hnmEb.Data(),hnmQ.Data(),hnmpt.Data());
	htt = Form("#DeltaZ %s %s %s",hnmEb.Data(),hnmQ.Data(),hnmpt.Data());
	h2 = new TH2F(hnm.Data(),htt.Data(), kNBinsPerSect*18, 0, 18., kNBinsDCAZ, -kMaxDCAZ,kMaxDCAZ);
	h2->SetXTitle("sector");
	h2->SetYTitle("#DeltaZ");
	fHManDCA[iQ][iPt]->AddHisto(h2, GetHistoID(iEta,kHDCAZ)*10);
	//
      }
      //
    }
    //
  }
  //
  for (int iSide=0;iSide<kNSides;iSide++) {
    //
    TString hnmSd = Form("Side%s",kSideName[iSide]);
    TString httSd = Form("Side%s",kSideName[iSide]);
    //    
    if (!fHManMtcSect[iSide]) {
      fHManMtcSect[iSide] = new HistoManager(Form("mtc_%s",hnmSd.Data()),fOutFName.Data());
      fHistManArr.Add(fHManMtcSect[iSide]);
    }
    //
    // TPC reference tracks for matching vs phi
    hnm = Form("%s_TPCvsSect",hnmSd.Data());
    htt = Form("TPC tracks vs Sector %s",httSd.Data());
    h2 = new TH2F(hnm.Data(),htt.Data(), kNBinsPerSect*18, 0, 18., kNBinsQ2Pt, -kMaxQ2Pt,kMaxQ2Pt);
    h2->SetXTitle("sector");
    h2->SetYTitle("q/p_{T}");
    fHManMtcSect[iSide]->AddHisto(h2, GetHistoID(iSide,kHTrcTPCSect)*10);
    //
    //
    // TPC-ITS tracks for matching vs phi
    hnm = Form("%s_TPCITSvsSect",hnmSd.Data());
    htt = Form("TPC-ITS tracks vs Sector %s",httSd.Data());
    h2 = new TH2F(hnm.Data(),htt.Data(), kNBinsPerSect*18, 0, 18., kNBinsQ2Pt, -kMaxQ2Pt,kMaxQ2Pt);
    h2->SetXTitle("sector");
    h2->SetYTitle("q/p_{T}");
    fHManMtcSect[iSide]->AddHisto(h2, GetHistoID(iSide,kHTrcTPCITSSect)*10);
    //
    // same with SPD request
    hnm = Form("%s_TPCITSSPDvsSect",hnmSd.Data());
    htt = Form("TPC-ITS(spd) tracks vs Sector %s",httSd.Data());
    h2 = new TH2F(hnm.Data(),htt.Data(), kNBinsPerSect*18, 0, 18., kNBinsQ2Pt, -kMaxQ2Pt,kMaxQ2Pt);
    h2->SetXTitle("sector");
    h2->SetYTitle("q/p_{T}");
    fHManMtcSect[iSide]->AddHisto(h2, GetHistoID(iSide,kHTrcTPCITSSPDSect)*10);
    //

  }
  //
  if (!fHManMtcEta) {
    fHManMtcEta = new HistoManager(Form("mtceta"),fOutFName.Data());
    fHistManArr.Add(fHManMtcEta);
  }

  // TPC reference tracks for matching vs eta
  hnm = Form("TPCvsEta");
  htt = Form("TPC tracks vs #eta");
  h2 = new TH2F(hnm.Data(),htt.Data(), kNBinsEta, -kMaxEta, kMaxEta, kNBinsQ2Pt, -kMaxQ2Pt,kMaxQ2Pt);
  h2->SetXTitle("#eta");
  h2->SetYTitle("q/p_{T}");
  fHManMtcEta->AddHisto(h2, kHTrcTPCEta*10);
  //
  //
  // TPC-ITS tracks for matching vs eta
  hnm = Form("TPCITSvsEta");
  htt = Form("TPC-ITS tracks");
  h2 = new TH2F(hnm.Data(),htt.Data(), kNBinsEta, -kMaxEta, kMaxEta, kNBinsQ2Pt, -kMaxQ2Pt,kMaxQ2Pt);
  h2->SetXTitle("#eta");
  h2->SetYTitle("q/p_{T}");
  fHManMtcEta->AddHisto(h2, kHTrcTPCITSEta*10);
  //
  // TPC-ITS tracks for matching vs eta
  hnm = Form("TPCITSspdEta");
  htt = Form("TPC-ITS(spd) tracks");
  h2 = new TH2F(hnm.Data(),htt.Data(), kNBinsEta, -kMaxEta, kMaxEta, kNBinsQ2Pt, -kMaxQ2Pt,kMaxQ2Pt);
  h2->SetXTitle("#eta");
  h2->SetYTitle("q/p_{T}");
  fHManMtcEta->AddHisto(h2, kHTrcTPCITSSPDEta*10);
  //
  //  PostData(1, fOutput);
  //  
  AliInfoF("Runing in %s mode",fUseMC ? "MC":"Data");
}


//________________________________________________________________________
void rsQA::UserExec(Option_t *) 
{
  // Main loop
  //
  AliAnalysisManager* anMan = AliAnalysisManager::GetAnalysisManager();
  AliESDInputHandler *handler = (AliESDInputHandler*)anMan->GetInputEventHandler();
  AliESDEvent* esdEv = (AliESDEvent*)handler->GetEvent();
  if (fUseMC) {
    AliMCEventHandler* eventHandler = (AliMCEventHandler*)anMan->GetMCtruthEventHandler();
    if (!eventHandler) { printf("ERROR: Could not retrieve MC event handler\n"); return; }
    fMCEvent = eventHandler->MCEvent();
    fStack = fMCEvent->Stack();
  }
  //
  if(!esdEv) {AliInfo("no ESD"); return;} 
  //
  ProcessEvent(esdEv);
  //
}      

//________________________________________________________________________
void rsQA::Terminate(Option_t *) 
{
  Printf("Terminating...");
  //
  if (fFileOut) {
    fFileOut->cd();
    if (fTreeOut) fTreeOut->Write();
    delete fTreeOut;
    fFileOut->Close();
    delete fFileOut;
  }
  //
  //
  for (int ih=0;ih<fHistManArr.GetEntries();ih++) {
    ((HistoManager*)fHistManArr[ih])->Write();
  }

}

//________________________________________________________________________
void rsQA::ProcessEvent(AliESDEvent* ev)
{

  AliMagF* fld = (AliMagF*) TGeoGlobalMagField::Instance()->GetField();
  if (!fld) ev->InitMagneticField();

  const AliESDVertex* vtx = ev->GetPrimaryVertex();
  if (!vtx->GetStatus()) return;
  if (TMath::Abs(vtx->GetZ())>kMaxZVtx) return;
  //
  const AliESDVertex* vtxTPC = ev->GetPrimaryVertexTPC();
  if (!vtxTPC->GetStatus()) return;
  if (TMath::Abs(vtxTPC->GetZ()-vtx->GetZ())>1) return;

  fBz = ev->GetMagneticField();
  int ntr = ev->GetNumberOfTracks();
  for (int itr=0;itr<ntr;itr++) {
    AliESDtrack* trc = ev->GetTrack(itr);
    //    if (!trc->IsOn(AliESDtrack::kTPCrefit)) continue;
    Bool_t res = ProcessTrack(trc,vtx);
  }
}

//________________________________________________________________________
Bool_t rsQA::ProcessTrack(const AliESDtrack* trc, const AliESDVertex* vtx)
{
  // 
  static int trCount=0;
  trCount++;

  // reference TPC track
  if (!trc->IsOn(AliESDtrack::kTPCrefit)) return kFALSE;
  if (trc->GetKinkIndex(0)>0) return kFALSE;

  const AliExternalTrackParam* tpcPar = trc->GetTPCInnerParam();
  if (!tpcPar || TMath::Abs(tpcPar->GetX())>3) return kFALSE;
  const AliExternalTrackParam* tpcIP = trc->GetInnerParam();
  if (!tpcIP) return kFALSE;
  //
  double eta = tpcPar->Eta();
  if (TMath::Abs(eta)>kMaxEta) return kFALSE;
  int ietaBin = (eta + kMaxEta)/(2.*kMaxEta/kNEtaRanges);
  if (ietaBin<0) return kFALSE;
  if (ietaBin>=kNEtaRanges) return kFALSE;
  //
  if (trc->GetTPCncls()<70) return kFALSE;
  //  if (trc->GetTPCchi2()/trc->GetTPCncls() > 5) return kFALSE;
  //
  Double_t dz[2],cov[3];
  AliExternalTrackParam tpcTr(*tpcPar);
  if (!tpcTr.PropagateToDCA(vtx, fBz, 99, dz, cov)) return kFALSE;
  double dcaY = dz[0], dcaZ = dz[1];
  if (TMath::Abs(dcaY)>kMaxDCAY) return kFALSE;
  if (TMath::Abs(dcaZ)>kMaxDCAZ) return kFALSE;
  // reference track cut

  int side = tpcIP->GetTgl()>0 ? kSideA : kSideC;
  int sideZ = tpcIP->GetZ()>0  ? kSideA : kSideC;
  if (side!=sideZ) return kFALSE;
  //

  double sector = tpcIP->PhiPos();
  if (sector<0) sector+= TMath::Pi()*2;
  sector *= 9./TMath::Pi();
  //
  double pt  = tpcPar->Pt();
  double q2pt = tpcPar->GetSigned1Pt();
  int sgq = q2pt>0 ? kQPlus:kQMinus;  
  //
  // fill DCA histos
  for (int iPt=0;iPt<kNPtRanges;iPt++) {
    if (pt < kPtRange[iPt]) break;
    //    if (iPt>1) printf("Fill %d to %d side:%d sideZ:%d iEta:%d q2pt: %+.1e\n",trCount,GetHistoID(ietaBin,kHDCAY)*10,side,sideZ,ietaBin,q2pt);
    fHManDCA[sgq][iPt]->GetHisto2F( GetHistoID(ietaBin,kHDCAY)*10 )->Fill(sector,dcaY);
    fHManDCA[sgq][iPt]->GetHisto2F( GetHistoID(ietaBin,kHDCAZ)*10 )->Fill(sector,dcaZ);
  }

  // TPC reference for matching
  fHManMtcSect[side]->GetHisto2F( GetHistoID(side,kHTrcTPCSect)*10 )->Fill(sector, q2pt);
  fHManMtcEta->GetHisto2F( kHTrcTPCEta*10 )->Fill(eta, q2pt);
  //
  //
  // ITS request
  if (!trc->IsOn(AliESDtrack::kITSrefit) || trc->GetITSNcls()<3) return kFALSE;
  fHManMtcSect[side]->GetHisto2F( GetHistoID(side,kHTrcTPCITSSect)*10 )->Fill(sector, q2pt);
  fHManMtcEta->GetHisto2F( kHTrcTPCITSSect*10 )->Fill(eta, q2pt);
  //
  // SPD request
  if (!(trc->HasPointOnITSLayer(0)||trc->HasPointOnITSLayer(1))) return kFALSE;
  fHManMtcSect[side]->GetHisto2F( GetHistoID(side,kHTrcTPCITSSPDSect)*10 )->Fill(sector, q2pt);
  fHManMtcEta->GetHisto2F( kHTrcTPCITSSPDSect*10 )->Fill(eta, q2pt);
  //
  return kTRUE;
}

Int_t rsQA::PdgToPid(int pdg) const
{
  // convert PDG to Alice PID
  pdg = TMath::Abs(pdg);
  for (int i=0;i<AliPID::kSPECIESC;i++) if (AliPID::ParticleCode(i)==pdg) return i;
  return -1;
}
