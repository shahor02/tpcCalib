#include "HistoManager.h"
#include "TF1.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TString.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TStyle.h"
#include "SaveCanvas.h"

HistoManager* ProcDCA(HistoManager* hman, int ietaMin,int ietaMax);

enum {kQPlus,kQMinus,kNCharges};
enum {kSideA,kSideC,kNSides};
enum {kHDCAY, kHDCAZ};
enum {kHTrcTPCSect,kHTrcTPCITSSect,kHTrcTPCITSSPDSect};
enum {kHTrcTPCEta,kHTrcTPCITSEta,kHTrcTPCITSSPDEta};  
enum {kNBinsPerSect=10, kNBinsEta=10, kNBinsQ2Pt=40, kNBinsDCA=160, kNPtRanges=3};
const char *kQName[kNCharges] = {"P","N"};
const char *kQTitle[kNCharges] = {"+","-"};
const char *kSideName[kNSides] = {"A","C"};

// special selections
Bool_t modeY = kTRUE;
Int_t iEtaSideMin = 0;  // eta bin starting from the center
Int_t iEtaSideMax = kNBinsEta/2 - 1; // and towards positive eta's

void SetModeY() {modeY = kTRUE;}
void SetModeZ() {modeY = kFALSE;}

Bool_t LoadHM(const char* fname, const char* pref);
int GetHistoID(int side, int kH) {return side*10+kH;}

HistoManager 
*fHManDCA[kNCharges+1][kNPtRanges], 
  *fHManMtcSect[kNSides], 
  *fHManMtcEta;
 
void ProcMiniQA(const char* fname, const char* pref="", int color=kBlack, int mStyle=20, float mSize=0.6)
{
  if (!LoadHM(fname,pref)) return;
  for (int iPt=0;iPt<kNPtRanges;iPt++) {
    for (int iQ=0;iQ<kNCharges+1;iQ++) {
      fHManDCA[iQ][iPt]->SetColor(color);
      fHManDCA[iQ][iPt]->SetMarkerStyle(mStyle);
      fHManDCA[iQ][iPt]->SetMarkerSize(mSize);
    }
  }
  //
}

TCanvas* DrawDCAPhi(int iq, int ipt, TCanvas* cnv=0, 
		    float rangeMin=1, float rangeMax=-1,
		    float rangeMinS=1, float rangeMaxS=-1)
{
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  Bool_t optSame = cnv;
  if (!cnv) {
    TString cname = Form("cnvq%dp%d",iq,ipt);
    cnv = new TCanvas(cname, cname, 1000,700);
    cnv->Divide(2,2);
  }
  HistoManager* hman = fHManDCA[iq][ipt];
  for (int iside=0;iside<kNSides;iside++) {

    int bEtaMin,bEtaMax;
    if (iside==kSideA) {
      bEtaMin = kNBinsEta/2 + iEtaSideMin;
      bEtaMax = kNBinsEta/2 + iEtaSideMax;
    }
    else {
      bEtaMin = kNBinsEta/2 - iEtaSideMax-1;
      bEtaMax = kNBinsEta/2 - iEtaSideMin-1;
    }
    
    HistoManager *sdArr = ProcDCA(hman, bEtaMin, bEtaMax);
    cnv->cd(iside+1);
    if (!optSame) gPad->SetLeftMargin(0.15);
    int hid = modeY ? 1 : 11;
    TH1* h1 = sdArr->GetHisto1F(hid);
    if (!optSame && rangeMin<rangeMax) {
      h1->SetMinimum(rangeMin);
      h1->SetMaximum(rangeMax);
    }
    h1->Draw(optSame ? "same":"");
    gPad->SetGrid();
    gPad->Modified();
    h1->GetYaxis()->SetTitleSize(0.05);
    h1->GetYaxis()->SetTitleOffset(1.4);
    if (!optSame) AddLabel(Form("Side %s",kSideName[iside]),0.2, 0.87,kBlack,0.05);

    gPad->Update();
    //
    //
    cnv->cd(iside+1+2);
    gPad->SetLeftMargin(0.15);
    hid = modeY ? 2 : 12;
    h1 = sdArr->GetHisto1F(hid);
    if (!optSame && rangeMinS<rangeMaxS) {
      h1->SetMinimum(rangeMinS);
      h1->SetMaximum(rangeMaxS);
    }
    h1->Draw(optSame ? "same":"");
    gPad->SetGrid();
    gPad->Modified();
    h1->GetYaxis()->SetTitleSize(0.05);
    h1->GetYaxis()->SetTitleOffset(1.4);
    gPad->Update();
    //   
    //    
  }
  if (!optSame) {
    cnv->cd(1);
    TString lbl;
    if (iq!=2) lbl += Form("q%c ",iq==kQPlus ? '+':'-');
    else       lbl += "q+,q- ";
    if (ipt) lbl += Form("p_{T}>%d",ipt);
    else     lbl += "all p_{T}";
    AddLabel(lbl.Data(),0.8, 0.96,kBlack,0.05);
  }

  if (rangeMin>rangeMax) { // auto rescaling
    float mxA=0,mnA=0,mxC=0,mnC=0;
    cnv->cd(1);
    GetHistosMinMaxRange(gPad,mnA,mxA);
    cnv->cd(2);
    GetHistosMinMaxRange(gPad,mnC,mxC);
    float mn = TMath::Min(mnA,mnC);
    float mx = TMath::Max(mxA,mxC);
    if (mn<mx) {
      cnv->cd(1);
      SetHistosMinMaxRange(gPad,mn,mx,0.1,0.1);
      cnv->cd(2);
      SetHistosMinMaxRange(gPad,mn,mx,0.1,0.1);
    }
    else printf("Failed to extract min/max of residual histos");
  }

  if (rangeMin>rangeMax) { // auto rescaling
    float mxA=0,mnA=0,mxC=0,mnC=0;
    cnv->cd(3);
    GetHistosMinMaxRange(gPad,mnA,mxA);
    cnv->cd(4);
    GetHistosMinMaxRange(gPad,mnC,mxC);
    float mn = TMath::Min(mnA,mnC);
    float mx = TMath::Max(mxA,mxC);
    if (mn<mx) {
      cnv->cd(3);
      SetHistosMinMaxRange(gPad,0,mx,0.1,0.0);
      cnv->cd(4);
      SetHistosMinMaxRange(gPad,0,mx,0.1,0.0);
    }
    else printf("Failed to extract min/max of residual histos");
  }


  cnv->cd();
  gPad->Modified();
  gPad->Update();
  return cnv;
}


Bool_t LoadHM(const char* fname, const char* pref)
{
  
  for (int iPt=0;iPt<kNPtRanges;iPt++) {
    TString hnmpt = Form("ptbin%d",iPt);
    //
    for (int iQ=0;iQ<kNCharges;iQ++) {
      TString hnmQ = Form("q%s",kQName[iQ]);
      fHManDCA[iQ][iPt] = new HistoManager(Form("dca_%s_%s",hnmQ.Data(),hnmpt.Data()), fname, 1, pref);
      if (fHManDCA[iQ][iPt]->GetEntriesFast()<1) return kFALSE;
    }
    //
    // charges added
    fHManDCA[kNCharges][iPt] = fHManDCA[kQPlus][iPt]->CreateClone("PM");
    fHManDCA[kNCharges][iPt]->AddHistos(fHManDCA[kQMinus][iPt]);
  }
  

  //
  for (int iSide=0;iSide<kNSides;iSide++) {
    TString hnmSd = Form("Side%s",kSideName[iSide]);
    fHManMtcSect[iSide] = new HistoManager(Form("mtc_%s",hnmSd.Data()),fname, 1, pref);
    if (fHManMtcSect[iSide]->GetEntriesFast()<1) return kFALSE;
    //
  }
  //
  fHManMtcEta = new HistoManager(Form("mtceta"),fname,1,pref);
  if (fHManMtcEta->GetEntriesFast()<1) return kFALSE;
  //
  return kTRUE;
}


HistoManager* ProcDCA(HistoManager* hman, int ietaMin,int ietaMax)
{
  printf("Adding eta bin %d - %d\n",ietaMin,ietaMax);
  HistoManager* harr = new HistoManager();
  TObjArray arrTmp;
  TH2* h2y = hman->GetHisto2F(GetHistoID(ietaMin,kHDCAY)*10);
  TString nms = Form("%s_eta%d_%d",h2y->GetName(),ietaMin,ietaMax);
  h2y = (TH2*)h2y->Clone(nms.Data());
  TH2* h2z = hman->GetHisto2F(GetHistoID(ietaMin,kHDCAZ)*10);
  nms = Form("%s_eta%d_%d",h2z->GetName(),ietaMin,ietaMax);
  h2z = (TH2*)h2z->Clone(nms.Data());
  //
  double maxDCA = h2y->GetYaxis()->GetXmax();
  static TF1* gs = new TF1("gs","gaus", -maxDCA, maxDCA);
  for (int ieta=ietaMin+1;ieta<=ietaMax;ieta++) {
    int hid = 0;
    hid = GetHistoID(ieta,kHDCAY)*10;
    h2y->Add(hman->GetHisto2F(hid));
    hid = GetHistoID(ieta,kHDCAZ)*10;
    h2z->Add(hman->GetHisto2F(hid));
    //
  }  
  harr->AddHisto(h2y,0);
  h2y->FitSlicesY(gs,0,-1,0,"qnrs0", &arrTmp);
  ((TH1*)arrTmp[1])->SetYTitle("<#DeltaY>, cm"); //Form("<#DeltaY> Side %s, cm", iSide==0 ? "A":"C"));
  ((TH1*)arrTmp[2])->SetYTitle("<#sigmaY>, cm"); ///Form("<#sigmaY> Side %s, cm", iSide==0 ? "A":"C"));    
  harr->AddHisto((TH1*)arrTmp[1],1);
  harr->AddHisto((TH1*)arrTmp[2],2);
  arrTmp.SetOwner(kFALSE);arrTmp.Clear();
  //
  harr->AddHisto(h2z,10);
  h2z->FitSlicesY(gs,0,-1,0,"qnrs0", &arrTmp);
  ((TH1*)arrTmp[1])->SetYTitle("<#DeltaZ>, cm"); //Form("<#DeltaY> Side %s, cm", iSide==0 ? "A":"C"));
  ((TH1*)arrTmp[2])->SetYTitle("<#sigmaZ>, cm"); ///Form("<#sigmaY> Side %s, cm", iSide==0 ? "A":"C"));    
  harr->AddHisto((TH1*)arrTmp[1],11);
  harr->AddHisto((TH1*)arrTmp[2],12);
  arrTmp.SetOwner(kFALSE);arrTmp.Clear();
  //
  harr->SetColor(hman->GetHisto2F(0)->GetMarkerColor());
  harr->SetMarkerStyle(hman->GetHisto2F(0)->GetMarkerStyle());  
  harr->SetMarkerSize(hman->GetHisto2F(0)->GetMarkerSize());    
  //
  return harr;
}
