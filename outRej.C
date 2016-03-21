#include <TSystem.h>
#include <TTree.h>
#include <TBranch.h>
#include <TFile.h>
#include <TVectorF.h>
#include <TVectorD.h>
#include <TString.h>
#include <TMath.h>
#include <TGeoMatrix.h>
#include <TGeoGlobalMagField.h>
#include <TGrid.h>
#include <TNDArray.h>
#include <TPad.h>
#include <THn.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TEnv.h>
#include <TObjArray.h>
#include <TLine.h>
#include <TLatex.h>
#include <TStopwatch.h>
#include <TGraphErrors.h>
#include "AliExternalTrackParam.h"
#include "AliTPCcalibAlignInterpolation.h"
#include "AliGeomManager.h"
#include "AliCDBManager.h"
#include "AliGRPManager.h"
#include "AliMagF.h"
#include "AliSysInfo.h"
#include "TStatToolkit.h"
#include "AliSymMatrix.h"
#include "AliTPCChebCorr.h"

enum {kNSect=18,kNSect2=2*kNSect,kNROC=4*kNSect,kNPadRows=159, kNRowIROC=63, kNRowOROC1=64, kNRowOROC2=32};
const float kSecDPhi=20.f*TMath::DegToRad(), kSecDPhiH = kSecDPhi*0.5f;
const float kMaxY2X = TMath::Tan(kSecDPhiH);  // max Y/X in sector coordinates (w/o excluding dead zone)
const char* kDriftFileName= "fitDrift";
const float kMaxZ2X = 1.0f, kZLim = 250.0f; // kMaxZ2X is Z/X range, while kZLim is endcap position

TFile *chunkFile = 0;
TTree* tree = 0;

TVectorF vecLocalDelta(kNPadRows);
TVectorF *vecDY=0,*vecDZ=0,*vecZ=0,*vecR=0,*vecSec=0,*vecPhi=0, *vecDYITS=0, *vecDZITS=0;
UShort_t npValid = 0;
Int_t timeStamp = 0;
AliExternalTrackParam* param = 0;
int fRun = 0;
float fBz=0;
float fMaxQ2Pt = 3;

TVectorD     *fVDriftParam=0;
TGraphErrors *fVDriftGraph=0;  
Float_t fCorrTime = 0;


TCanvas* cnv = 0, *cnvRS=0;
TGraph *grTrY=0,*grTrZ=0;
TGraph *grYX=0, *grZX=0;
TGraph *grYXRejMI=0, *grZXRejMI=0;
TGraph *grYXRejRS=0, *grZXRejRS=0;

TGraph *grYdfMA=0,*grZdfMA=0;
TGraph *grYdfLL=0,*grZdfLL=0;
TGraph *grYdfML=0,*grZdfML=0;

TH2F *hBaseYX=0,*hBaseZX=0;
TObjArray auxObj;

TH2F *hOutTrMaxDY=0, *hOutTrMaxDZ=0;
TH2F *hOutClDYZMA=0, *hOutClDYZLL=0;

float sTrc[kNPadRows];
float arrXtrack[kNPadRows],arrYtrack[kNPadRows],arrZtrack[kNPadRows];
float arrPhi[kNPadRows],arrX[kNPadRows],arrY[kNPadRows],arrZ[kNPadRows],
  arrDY[kNPadRows],arrDZ[kNPadRows],arrDYITS[kNPadRows],arrDZITS[kNPadRows];

float arrOKE[16*kNPadRows];


float yDiffMA[kNPadRows],zDiffMA[kNPadRows];
float yDiffLL[kNPadRows],zDiffLL[kNPadRows];
float yDiffML[kNPadRows],zDiffML[kNPadRows];

Bool_t rejMI[kNPadRows];
Bool_t rejRS[kNPadRows];
int arrSectID[kNPadRows], nCl=0;


TLatex* AddTxtLabel(const char*txt,float x,float y,int color=kRed,float size=0.04);

void CheckTrack(float q2pt, int np, const float *x, const float* y, const float *z, 
		const float* resy, const float *resz, const int* sect36,
		float *xTrc,float *yTrc, float *zTrc);

int CheckResiduals(int np, const float *x, const float *y, const float *z, const int *sec36, Bool_t* kill,
		   int nVois=3,float cut=16.);

void FitCircle(int np, const float* x, const float* y, 
	       float &xc, float &yc, float &r2, float* dy=0);

void LoadVDrift();
float GetDriftCorrection(float z, float x, float phi, int rocID);

//void InitForBugFix(const char* ocdb="local:///cvmfs/alice.cern.ch/calibration/data/2015/OCDB");
void InitForBugFix(const char* ocdb="raw://");
void Process(int ev);
void DrawSectorEdges(int nc,const int* sec, const float* x, const float* y, float mn, float mx, Bool_t verbose=kTRUE);
Bool_t FitPoly2(const float* x,const float* y, const float* w, int np, float *res, float *err);
Bool_t FitPoly1(const float* x,const float* y, const float* w, int np, float *res, float *err);

void MaxLeftRightDev(int np, const float* x, const float *y, const int nVoisin, 
		     float *diffY, float *diffLR, float *kink);

void DiffToMA(int np, const float* x, const float *y, const int winLR, float* diffMA);

int DiffToLocLine(int np, const float* x, const float *y, const int nVoisin, float *diffY);
int DiffToMedLine(int np, const float* x, const float *y, const int nVoisin, float *diffY);

void FindOutliersKinks(int np, const float* x, const float *y, const int nVoisin, 
		       float *diffY, float *diffYE, float *kink, float *kinkE);

void FixAlignmentBug(int sect, float q2pt, float bz,
		     float& alp, float& xCl, float &zTr, float &deltaY, float &deltaZ);
float tgpXY(float x, float y, float q2p, float bz);

//------------------------------------
float RoFunc(int np, const float* x, const float* y, float b, float &aa);
Float_t SelKthMin(int k, int np, float* arr);
void medFit(int np, const float* x, const float* y, float &a, float &b, float delI=0.f);
//------------------------------------

void Init(
	  //const char* fname="/data/testTPC/res245231/ResidualTrees.root",
	  //const char* fname="/data/testTPC/res244918/ResidualTrees.root",
	  //const char* fname="resTree244918/alice/data/2015/LHC15o/000244918/cpass0_pass1/ResidualMerge/002/ResidualTrees.root",
	  const char* fname="/data/testTPC/resTree000245231/alice/data/2015/LHC15o/000245231/cpass0_pass1/ResidualMerge/001/ResidualTrees.root",
	  //const char* fname="/data/testTPC/res245231/ResidualTrees.root",
	  //const char* fname="/data/testTPC/resTree244418/ResidualTrees.root",
	  int run=245231
	  //int run=244418
	  //int run=244918
)
{

  fRun = run;

  chunkFile = TFile::Open(fname);
  tree = (TTree*)chunkFile->Get("delta");
  tree->SetBranchAddress("timeStamp",&timeStamp);
  tree->SetBranchAddress("vecR.",&vecR);
  tree->SetBranchAddress("vecSec.",&vecSec);
  tree->SetBranchAddress("vecPhi.",&vecPhi);
  tree->SetBranchAddress("vecZ.",&vecZ);
  tree->SetBranchAddress("track.",&param);
  tree->SetBranchAddress("npValid",&npValid);
  tree->SetBranchAddress("trd0.",&vecDY);
  tree->SetBranchAddress("trd1.",&vecDZ);
  tree->SetBranchAddress("its0.",&vecDYITS);
  tree->SetBranchAddress("its1.",&vecDZITS);
  //
  InitForBugFix();
  //
  hOutTrMaxDY = new TH2F("hOutTrMaxDY","maxDY to helix",10,0,fMaxQ2Pt,100,0,2.);
  hOutTrMaxDZ = new TH2F("hOutTrMaxDZ","maxDY to helix",10,0,fMaxQ2Pt,100,0,2.);

  hOutClDYZMA = new TH2F("hOutClDYZMA","DyDz cluster MovingAverage",100,-10.,10.,100,-10.,10.);  
  hOutClDYZLL = new TH2F("hOutClDYZLL","DyDz cluster MovingLine"   ,100,-10.,10.,100,-10.,10.);  
  LoadVDrift(); //!!!

}

// RS settings
float kMaxDevHelixY = 0.3;
float kMaxDevHelixZ = 0.3;
float kMaxStdDev    = 16.0;
float kMaxRejFrac   = 0.15;
int   fMinNCl       = 50;
int   kNVoisin      = 3;


void Process(int ev)
{
  const Float_t *vLocalDelta=vecLocalDelta.GetMatrixArray();
  if (!tree) Init();
  tree->GetEntry(ev);
  //
  float q2pt = param->GetParameter()[4];
  float tgLam = param->GetParameter()[3];
  if (TMath::Abs(q2pt)>fMaxQ2Pt) {printf("bad pt\n"); return;}
  //
  //-----------------------------------------------------

  //-----------------------------------------------------
  // MI settings
  const Int_t   kMaxSkippedCluster=10;  // 10 cluster
  const Float_t kMaxRMSTrackCut=2.0;    // maximal RMS (cm) between the tracks 
  const Float_t kMaxRMSClusterCut=0.3;    // maximal RMS (cm) between the cluster and local mean
  const Float_t kMaxDeltaClusterCut=0.5;    // maximal delta(cm) between the cluster and local mean

  Float_t rmsTrack=3, rmsCluster=1;
  Int_t nSkippedCluster=AliTPCcalibAlignInterpolation::CalculateDistance(*vecDY,*vecDYITS, 
									 *vecSec, vecLocalDelta, 
									 npValid, rmsTrack, rmsCluster,1.5);
  printf("SkippedCl: %d RMStr: %6.2f RMScl: %6.2f\n",nSkippedCluster,rmsTrack,rmsCluster);
  Bool_t trRejMI = kFALSE;
  if (nSkippedCluster>kMaxSkippedCluster) {printf("Kill SkipCl\n"); trRejMI = kTRUE;}
  if (rmsTrack>kMaxRMSTrackCut)           {printf("Kill RMStr\n");  trRejMI = kTRUE;}
  if (rmsCluster>kMaxRMSClusterCut)       {printf("Kill RMScl\n");  trRejMI = kTRUE;}
  //-----------------------------------------------------


  const Float_t *vSec= vecSec->GetMatrixArray();
  const Float_t *vPhi= vecPhi->GetMatrixArray();
  const Float_t *vR  = vecR->GetMatrixArray();
  const Float_t *vZ  = vecZ->GetMatrixArray();
  const Float_t *vDY = vecDY->GetMatrixArray();
  const Float_t *vDZ = vecDZ->GetMatrixArray();
  const Float_t *vDYITS = vecDYITS->GetMatrixArray();
  const Float_t *vDZITS = vecDZITS->GetMatrixArray();
  //

  nCl = 0;

  float fCorrTime = (fVDriftGraph!=NULL) ? fVDriftGraph->Eval(timeStamp):0; // for VDrift correction


  int nRejMI = 0;
  for (int ip=0;ip<npValid;ip++) { // 1st fill selected track data to buffer for eventual outlier rejection
    if (vR[ip]<1 || vDY[ip]<-900 || vDYITS[ip]<-900) continue;
    //
    arrX[nCl] = vR[ip];
    arrZ[nCl] = vZ[ip];
    arrDY[nCl] = vDY[ip];
    arrDZ[nCl] = vDZ[ip];
    arrDYITS[nCl] = vDYITS[ip];
    arrDZITS[nCl] = vDZITS[ip];
    arrPhi[nCl] = vPhi[ip];

    rejMI[nCl] = TMath::Abs(vLocalDelta[ip])> kMaxDeltaClusterCut;

    if (rejMI[nCl]) nRejMI++;

    int rocID = TMath::Nint(vSec[ip]);
    //
    if (1) {
      float phiITS = arrPhi[nCl];
      float xITS   = arrX[nCl];
      float zITS   = arrZ[nCl];  // ITS track Z was stored!!!
      arrZ[nCl] += arrDZ[nCl] - arrDZITS[nCl]; // recover ITS-TRD track position from ITS and deltas
      FixAlignmentBug(rocID, q2pt, fBz, arrPhi[nCl], arrX[nCl], arrZ[nCl], arrDY[nCl],arrDZ[nCl]);
      FixAlignmentBug(rocID, q2pt, fBz, phiITS     ,      xITS, zITS     , arrDYITS[nCl],arrDZITS[nCl]);
    }
    else {
      arrZ[nCl] += arrDZ[nCl] - arrDZITS[nCl]; // recover ITS-TRD track position from ITS and deltas
    }
    if (arrPhi[nCl]<0) arrPhi[nCl] += 2.*TMath::Pi();
    
    //
    // apply drift velocity calibration if available
    arrDZ[nCl] += GetDriftCorrection(arrZ[nCl],arrX[nCl],arrPhi[nCl],rocID);
    //
    arrSectID[nCl] = rocID%kNSect2; // 0-36 for sectors from A0 to C17
    float sna = TMath::Sin(arrPhi[nCl]-(0.5f +rocID%kNSect)*kSecDPhi);
    float csa = TMath::Sqrt((1.f-sna)*(1.f+sna));
    //
    // by using propagation in cluster frame in AliTPCcalibAlignInterpolation::Process,
    // the X of the track is evaluated not at the pad-row x=r*csa but at x=r*sca-dy*sna
    double xrow = arrX[nCl]*csa;
    double dx   = arrDY[nCl]*sna;
    double xtr = xrow - dx;
    double ycl = arrX[nCl]*sna;           // the cluster Y in the sector frame is 
    double ytr = ycl + arrDY[nCl]*csa;    // and the track Y in the sector frame at x=xtr is 
    //
    double ztr = arrZ[nCl];               // Z of the track at x=xtr
    double zcl = ztr - arrDZ[nCl];        // and the Z of the cluster is Ztr-deltaZ
    //
    // Now we need to take the track to real pad-row X
    // 1) get track angle in the sector frame at xtr
    double tgXtr = tgpXY(xtr,ytr,q2pt,fBz);
    double csXtrInv = TMath::Sqrt(1.+tgXtr*tgXtr); // (inverse cosine of track angle)
    // 2) use linear extrapolation:
    ytr += dx*tgXtr;
    ztr += dx*tgLam*csXtrInv;
    //
    // assign to arrays and recalculate residuals
    arrX[nCl] = xrow;
    arrY[nCl] = ycl;
    arrZ[nCl] = zcl;//ztr;
    arrDY[nCl] = ytr - ycl;
    arrDZ[nCl] = ztr - zcl;
    //
    nCl++;
  }  
  //
  printf("Nclusters: %d\n",nCl);
  if (nCl<3) return;

  Bool_t trRejRS = kFALSE;
  int nRejRS = CheckResiduals(nCl,arrX,arrDY,arrDZ,arrSectID,rejRS, kNVoisin, kMaxStdDev);
  CheckTrack(q2pt, nCl, arrX,arrY,arrZ, arrDY,arrDZ, arrSectID, arrXtrack,arrYtrack,arrZtrack);
  float hmnY=1e9,hmxY=-1e9,hmnZ=1e9,hmxZ=-1e9;
  for (int i=nCl;i--;) {
    float val = arrYtrack[i];
    if (val<hmnY) hmnY = val;
    if (val>hmxY) hmxY = val;
    val = arrZtrack[i];
    if (val<hmnZ) hmnZ = val;
    if (val>hmxZ) hmxZ = val;
  }
  float maxDevY = TMath::Abs(hmxY-hmnY);
  float maxDevZ = TMath::Abs(hmxZ-hmnZ);

  if (nCl<fMinNCl) trRejRS = kTRUE;
  if (maxDevY>kMaxDevHelixY || maxDevZ>kMaxDevHelixZ) trRejRS = kTRUE;
  if (float(nRejRS)/nCl>kMaxRejFrac) trRejRS = kTRUE;

  // track deviation from parabola
  hOutTrMaxDY->Fill(TMath::Abs(q2pt),maxDevY);
  hOutTrMaxDZ->Fill(TMath::Abs(q2pt),maxDevZ);
  //
  // cluster deviations
  for (int icl=0;icl<nCl;icl++) {
    hOutClDYZMA->Fill(yDiffMA[icl],zDiffMA[icl]);
    hOutClDYZLL->Fill(yDiffLL[icl],zDiffLL[icl]);
  }
  //
  delete grTrZ;
  delete grTrY;
  grTrY = new TGraph(nCl,arrXtrack,arrYtrack);
  //  grTrZ = new TGraph(nCl,arrXtrack,arrZtrack);
  grTrZ = new TGraph(nCl,sTrc,arrZtrack);
  grTrY->SetMarkerStyle(20);
  grTrZ->SetMarkerStyle(20);

  delete grYX;
  delete grZX;

  delete grYdfMA;
  delete grZdfMA;
  grYdfMA = new TGraph(nCl,arrX,yDiffMA);
  grZdfMA = new TGraph(nCl,arrX,zDiffMA);
  grYdfMA->SetMarkerStyle(7);
  grZdfMA->SetMarkerStyle(7);

  delete grYdfLL;
  delete grZdfLL;
  grYdfLL = new TGraph(nCl,arrX,yDiffLL);
  grZdfLL = new TGraph(nCl,arrX,zDiffLL);
  grYdfLL->SetMarkerStyle(7);
  grZdfLL->SetMarkerStyle(7);
 
  delete grYXRejMI; grYXRejMI = 0;
  delete grZXRejMI; grZXRejMI = 0;
  
  delete grYXRejRS; grYXRejRS = 0;
  delete grZXRejRS; grZXRejRS = 0;

  auxObj.Delete();

  if (nRejRS) {
    grYXRejRS = new TGraph(nRejRS);
    grZXRejRS = new TGraph(nRejRS);
    grYXRejRS->SetMarkerStyle(25);
    grZXRejRS->SetMarkerStyle(25);
    grYXRejRS->SetMarkerColor(kBlue);
    grZXRejRS->SetMarkerColor(kBlue);
    grYXRejRS->SetMarkerSize(0.9);
    grZXRejRS->SetMarkerSize(0.9);
    int kr = 0;
    for (int i=0;i<nCl;i++) {
      if (!rejRS[i]) continue;
      grYXRejRS->SetPoint(kr, arrX[i],arrDY[i]);
      grZXRejRS->SetPoint(kr, arrX[i],arrDZ[i]);
      kr++;
    }
  }

  if (nRejMI) {
    grYXRejMI = new TGraph(nRejMI);
    grZXRejMI = new TGraph(nRejMI);
    grYXRejMI->SetMarkerStyle(24);
    grZXRejMI->SetMarkerStyle(24);
    grYXRejMI->SetMarkerColor(kGreen+2);
    grZXRejMI->SetMarkerColor(kGreen+2);
    grYXRejMI->SetMarkerSize(1.2);
    grZXRejMI->SetMarkerSize(1.2);
    int kr = 0;
    for (int i=0;i<nCl;i++) {
      if (!rejMI[i]) continue;
      grYXRejMI->SetPoint(kr, arrX[i],arrDY[i]);
      grZXRejMI->SetPoint(kr, arrX[i],arrDZ[i]);
      kr++;
    }
  }

  grYX = new TGraph(nCl,arrX,arrDY);
  grZX = new TGraph(nCl,arrX,arrDZ);
  grYX->SetMarkerStyle(7);
  grZX->SetMarkerStyle(7);
  grYX->SetMarkerColor(kRed);
  grZX->SetMarkerColor(kRed);
  //
  double mnX,mxX, mnY,mxY, mnZ,mxZ, mnYITS,mxYITS, mnZITS,mxZITS;
  TStatToolkit::GetMinMax(arrX,nCl, mnX,mxX);
  TStatToolkit::GetMinMax(arrDY,nCl, mnY,mxY);
  TStatToolkit::GetMinMax(arrDZ,nCl, mnZ,mxZ);
  //
  double d;
  d = 0.1*(mxX-mnX);
  mnX -= d;
  mxX += d; 
  //
  d = 0.1*(mxY-mnY);
  mnY -= d;
  mxY += d; 
  //
  d = 0.1*(mxZ-mnZ);
  mnZ -= d;
  mxZ += d; 
  //
  delete hBaseYX; 
  delete hBaseZX; 
  hBaseYX = new TH2F("hBaseYX","YX",100,mnX,mxX,100,mnY,mxY);
  hBaseYX->SetXTitle("X");   
  hBaseYX->SetYTitle("DY");
  hBaseZX = new TH2F("hBaseZX","ZX",100,mnX,mxX,100,mnZ,mxZ);
  hBaseZX->SetXTitle("X");   
  hBaseZX->SetYTitle("DZ");

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  if (!cnv) cnv = new TCanvas("cnv","cnv",1400,700);
  if (!cnvRS) cnvRS = new TCanvas("cnvRS","cnvRS",1400,700);
  cnv->Clear();
  cnv->Divide(2,2);
  cnv->cd(1);
  hBaseYX->Draw();
  grYX->Draw("p");
  if (grYXRejMI) grYXRejMI->Draw("p");
  if (grYXRejRS) grYXRejRS->Draw("p");
  DrawSectorEdges(nCl,arrSectID,arrX,arrDY,mnY,mxY,kTRUE);
  gPad->SetGrid();

  TString txtMI = Form("MI: TagCl: %d SkpCl: %d(%c) RMStr: %5.2f(%c) RMScl: %5.2f(%c) -> %s\n",nRejMI, 
		       nSkippedCluster,nSkippedCluster>kMaxSkippedCluster ? 'x':' ',
		       rmsTrack,rmsTrack>kMaxRMSTrackCut ? 'x':' ',  
		       rmsCluster, rmsCluster>kMaxRMSClusterCut ? 'x':' ',
		       trRejMI ? "REJ":"ACC");
  TLatex* lab = AddTxtLabel(txtMI.Data(),0.05,0.94,kBlue,0.03);
  auxObj.Add(lab);

  //
  cnv->cd(2);
  hBaseZX->Draw();
  grZX->Draw("p");
  if (grZXRejMI) grZXRejMI->Draw("p");
  if (grZXRejRS) grZXRejRS->Draw("p");
  DrawSectorEdges(nCl,arrSectID,arrX,arrDZ,mnZ,mxZ,kFALSE);
  gPad->SetGrid();
  //

  TString txtRS = Form("RS: TagCl: %d of %d (%.3f) DiffToHelix: Y:%.3f Z:%.3f -> %s\n",
		       nRejRS,nCl,float(nRejRS)/nCl,
		       maxDevY,maxDevZ,trRejRS ? "REJ":"ACC");
  lab = AddTxtLabel(txtRS.Data(),0.05,0.94,kRed,0.03);
  auxObj.Add(lab);


  cnv->cd(3);
  grTrY->Draw("ap");
  gPad->SetGrid();
  //
  cnv->cd(4);
  grTrZ->Draw("ap");
  gPad->SetGrid();
  //
  //
  cnvRS->Clear();
  cnvRS->Divide(2,2);

  cnvRS->cd(1);
  grYdfMA->Draw("ap"); gPad->SetGrid();

  cnvRS->cd(2);
  grZdfMA->Draw("ap"); gPad->SetGrid();

  cnvRS->cd(3);
  grYdfLL->Draw("ap"); gPad->SetGrid();

  cnvRS->cd(4);
  grZdfLL->Draw("ap"); gPad->SetGrid();

}

//_________________________________________________
void DrawSectorEdges(int nc,const int* sec, const float* x, const float* y, float mn, float mx, Bool_t verbose)
{
  const char *sname[2] = {"A","C"};
  const char *ioname[2] = {"IROC","OROC"};
  for (int i=1;i<nc;i++) {
    if (sec[i]==sec[i-1]) continue;
    TLine* ln = new TLine(0.5*(x[i]+x[i-1]),mn, 0.5*(x[i]+x[i-1]),mx);
    
    int sc0 = sec[i-1];
    int sc1 = sec[i];
    Bool_t side0 = (sc0/18)&0x1;
    Bool_t side1 = (sc1/18)&0x1;
    Bool_t Oroc0 = sc0/36;
    Bool_t Oroc1 = sc1/36;

    if (side0!=side1) ln->SetLineColor(kGreen+2);
    if (Oroc0!=Oroc1) ln->SetLineStyle(2);
    if (sc0%18 != sc1%18) ln->SetLineWidth(2);
    if (verbose) {
      printf("Sector change %d (%s%02d %s) -> %d (%s%02d %s) at X=%.2f\n",
	     sc0, sname[side0],sc0%18,ioname[Oroc0],
	     sc1, sname[side1],sc1%18,ioname[Oroc1],
	     x[i]);
    }
    ln->Draw();
    auxObj.Add(ln);
  }
}

//_________________________________________________
void InitForBugFix(const char* ocdb)
{
  ::Info(" AliTPCcalibAlignInterpolation::InitForBugFix","Alignment bug fix is requested\n");
  //
  // this requires the field and the geometry ...
  if (fRun<1) ::Fatal("tstw","InitForBugFix: Run number is not provided");
  Bool_t geomOK = AliGeomManager::GetGeometry() != 0;
  AliMagF* fld = (AliMagF*)TGeoGlobalMagField::Instance()->GetField();
  if (!geomOK || !fld) { // need to setup ocdb?
    AliCDBManager* man = AliCDBManager::Instance();
    if (!man->IsDefaultStorageSet()) man->SetDefaultStorage(ocdb);
    if (man->GetRun()!=fRun) man->SetRun(fRun);
  }
  if (!geomOK) {
    AliGeomManager::LoadGeometry();
    AliGeomManager::ApplyAlignObjsFromCDB("TPC");
  }
  if (!fld) {
    AliGRPManager grpMan;
    grpMan.ReadGRPEntry();
    grpMan.SetMagField();
    fld = (AliMagF*)TGeoGlobalMagField::Instance()->GetField();
  }
  fBz = fld->SolenoidField();
}


TLatex* AddTxtLabel(const char*txt,float x,float y,int color,float size)
{
  TLatex* lt = new TLatex(x,y,txt); 
  lt->SetNDC(); 
  lt->SetTextColor(color);
  lt->SetTextSize(size);
  lt->Draw();
  return lt;
}

void CheckTrack(float q2pt, int np, const float *x, const float* y, const float *z, 
		const float* resy, const float *resz, const int* sect36,
		float *xTrc,float *yTrc, float *zTrc)
{
  float ssum2=0;
  int sectCPrev=-1,sect0 = sect36[0]%kNSect; // align to 1st point
  double sn=0,cs=0;
  //  float sTrc[kNPadRows];
  sTrc[0] = 0.f;
  float crv = TMath::Abs(q2pt*fBz*(-0.299792458e-3f));
  for (int ip=0;ip<np;ip++) {
    double xp = x[ip];
    double yp = y[ip] + resy[ip];
    double zp = z[ip] + resz[ip]; // track Z coordinate is provided
    int sect = sect36[ip]%kNSect;
    if (sect!=sect0) { // rotate to reference sector
      if (sect!=sectCPrev) {
	double dalp = (sect - sect0)*20*TMath::DegToRad();
	sn = TMath::Sin(dalp);
	cs = TMath::Cos(dalp);
	sectCPrev = sect; // to not recalculate sin,cos every time
      }
      float xtmp = xp;
      xp = xtmp*cs-yp*sn;
      yp = yp*cs+xtmp*sn;
    }
    xTrc[ip] = xp;
    yTrc[ip] = yp;
    zTrc[ip] = zp;
    if (ip) {
      float dx = xp-xTrc[ip-1];
      float dy = yp-yTrc[ip-1];
      float ds2 = dx*dx+dy*dy;
      float ds  = TMath::Sqrt(ds2); // circular path
      if (ds*crv>0.01) { 
	// account for the arc-chord difference as 1st 2 terms of asin expansion	
	ds *= (1.f+ds2*crv*crv/24.f);
      }
      sTrc[ip] = sTrc[ip-1]+ds;
    }
  }
  //
  float xc=0,yc=0,r2=0;
  FitCircle(np,xTrc,yTrc,xc,yc,r2,yTrc);
  float pol1z[2],pol1zE[4] ;
  Bool_t resfZ = FitPoly1(sTrc, zTrc, 0, np, pol1z, pol1zE);
  //
  for (int ip=0;ip<np;ip++) zTrc[ip] -= pol1z[0]+sTrc[ip]*pol1z[1];
  // 
}

//_____________________________________________________
Bool_t FitPoly2(const float* x,const float* y, const float* w, int np, float *res, float *err)
{
  // poly2 fitter
  if (np<3) return kFALSE; // no enough points
  double sumW[5]={0},sumY[3]={0};
  for (int ip=np;ip--;) {
    double ww = w ? w[ip] : 1.0;
    sumW[0] += ww;
    sumY[0] += ww*y[ip];
    sumW[1] += (ww*=x[ip]); 
    sumY[1] += ww*y[ip];
    sumW[2] += (ww*=x[ip]); 
    sumY[2] += ww*y[ip];
    sumW[3] += (ww*=x[ip]); 
    sumW[4] += (ww*=x[ip]); 
  }
  double min00 = (sumW[2]*sumW[4]-sumW[3]*sumW[3]);
  double min01 = (sumW[1]*sumW[4]-sumW[3]*sumW[2]);
  double min02 = (sumW[1]*sumW[3]-sumW[2]*sumW[2]);
  double min11 = (sumW[0]*sumW[4]-sumW[2]*sumW[2]);
  double min12 = (sumW[0]*sumW[3]-sumW[2]*sumW[1]);
  double min22 = (sumW[0]*sumW[2]-sumW[1]*sumW[1]);
  double det  = sumW[0]*min00
    -           sumW[1]*min01
    +           sumW[2]*min02;
  if (TMath::Abs(det)<1e-12) return kFALSE;
  double detI = 1./det;
  double det0 = sumY[0]*min00
    -           sumW[1]*(sumY[1]*sumW[4]-sumW[3]*sumY[2]) 
    +           sumW[2]*(sumY[1]*sumW[3]-sumW[2]*sumY[2]);
  double det1 = sumW[0]*(sumY[1]*sumW[4]-sumW[3]*sumY[2]) 
    -           sumY[0]*min01
    +           sumW[2]*(sumW[1]*sumY[2]-sumY[1]*sumW[2]);
  double det2 = sumW[0]*(sumW[2]*sumY[2]-sumY[1]*sumW[3]) 
    -           sumW[1]*(sumW[1]*sumY[2]-sumY[1]*sumW[2]) 
    +           sumY[0]*min02;
  res[0] = det0*detI;
  res[1] = det1*detI;
  res[2] = det2*detI;
  //
  err[0] = min00*detI; // e00
  err[1] =-min01*detI; // e10
  err[2] = min11*detI; // e11
  err[3] = min02*detI; // e20
  err[4] =-min12*detI; // e21
  err[5] = min22*detI; // e21
  //
  return kTRUE;
}

//_____________________________________________________
Bool_t FitPoly1(const float* x,const float* y, const float* w, int np, float *res, float *err)
{
  // poly1 fitter
  if (np<2) return kFALSE; // no enough points
  double sumW[3]={0},sumY[2]={0};
  for (int ip=np;ip--;) {
    double ww = w ? w[ip]:1.0;
    sumW[0] += ww;
    sumY[0] += ww*y[ip];
    sumW[1] += (ww*=x[ip]); 
    sumY[1] += ww*y[ip];
    sumW[2] += (ww*=x[ip]);
  }
  double det  = sumW[0]*sumW[2] - sumW[1]*sumW[1];
  if (TMath::Abs(det)<1e-12) return kFALSE;
  double detI = 1./det;
  double det0 = sumY[0]*sumW[2] - sumW[1]*sumY[1];
  double det1 = sumW[0]*sumY[1] - sumY[0]*sumW[1];
  res[0] = det0*detI;
  res[1] = det1*detI;
  //
  err[0] = sumW[2]*detI; // e00
  err[1] =-sumW[1]*detI; // e10
  err[2] = sumW[0]*detI; // e11
  //
  return kTRUE;
}


//_______________________________________________
void FindOutliersKinks(int np, const float* x, const float *y, const int nVoisin, 
		       float *diffY, float *diffYE, float *kink, float *kinkE)
{
  // calculate the difference between the point and linear extrapolation from neigbourhood
  const float kEps = 1e-9;
  double sumX1b[kNPadRows+1],sumX2b[kNPadRows+1],sumYb[kNPadRows+1],sumYXb[kNPadRows+1];
  double *sumX1 = sumX1b+1, *sumX2 = sumX2b+1, *sumY0 = sumYb+1, *sumYX = sumYXb+1;
  //
  sumX1[-1]=0.f;
  sumX2[-1]=0.f;
  sumY0[-1]=0.f;
  sumYX[-1]=0.f;
  // accumulate matrix elements for whole array, element -1 is at 0
  double sX1=0., sX2=0., sY0=0., sYX=0.;

  for (int ip=0;ip<np;ip++) {
    sumX1[ip] = sX1 += x[ip];
    sumX2[ip] = sX2 += x[ip]*x[ip];
    sumY0[ip] = sY0 += y[ip];
    sumYX[ip] = sYX += y[ip]*x[ip];
    diffY[ip] = diffYE[ip] = 0.f;
    kink[ip] = kinkE[ip] = 0.f;
  }
  // 

  for (int ip=0;ip<np;ip++) {

    float slopL=0,slopLE=0,slopR=0,slopRE=0;
    float yLEst=0.f,yLEstE2I=0.f, yREst=0.f,yREstE2I=0.f;
    float &yEst = diffY[ip], &yEstE = diffYE[ip];
    float &tgKink = kink[ip], &tgKinkE = kinkE[ip];
    int ip0,ip1;
    //
    // estimate from the left
    ip0 = ip-nVoisin;
    if (ip0>-1) { 
      ip1 = ip-1;    // extract sum from ip0 to ip1 from cumulant, S00=nVoisin
      sX1 = sumX1[ip1] - sumX1[ip0-1]; // S01
      sX2 = sumX2[ip1] - sumX2[ip0-1]; // S11
      sY0 = sumY0[ip1] - sumY0[ip0-1]; // RHS0
      sYX = sumYX[ip1] - sumYX[ip0-1]; // RHS1
      double det = nVoisin*sX2 - sX1*sX1;
      if (det>kEps) {
	double detI = 1./det;
	// yLEst = offs + slop*x[ip] 
	// with offs=(sY0*sX2-sYX*sX1)/det and slop=(nVoisin*sYX-sY0*sX1)/det
	// inverse err^2 = 1/(errOffs+x^2*errSlop+2*x*errSlpOffs) with
	// errOffs = S11/det, errSlp=S00/det, errSlpOffs=-S01/det
	yLEst    = ((sY0*sX2-sYX*sX1) + (nVoisin*sYX-sY0*sX1)*x[ip])*detI;
	yLEstE2I = det/(sX2+x[ip]*x[ip]*nVoisin-2.f*x[ip]*sX1); // 1/err2
	slopL  = (nVoisin*sYX-sY0*sX1)*detI;
	slopLE = nVoisin*detI;
      }
    }
    // estimate from the right
    ip1 = ip+nVoisin;
    if (ip1<np) {   // extract sum from ip+1 to ip1 from cumulant, S00=nVoisin
      sX1 = sumX1[ip1] - sumX1[ip]; // S01
      sX2 = sumX2[ip1] - sumX2[ip]; // S11
      sY0 = sumY0[ip1] - sumY0[ip]; // RHS0
      sYX = sumYX[ip1] - sumYX[ip]; // RHS1
      double det = nVoisin*sX2 - sX1*sX1;
      if (det>kEps) {
	double detI = 1./det;
	yREst    = ((sY0*sX2-sYX*sX1) + (nVoisin*sYX-sY0*sX1)*x[ip])*detI;
	yREstE2I = det/(sX2+x[ip]*x[ip]*nVoisin-2.f*x[ip]*sX1); // 1/err2
	slopR  = (nVoisin*sYX-sY0*sX1)*detI;
	slopRE = nVoisin*detI;
      }
    }
    // joint estimate of expected value
    yEstE = yREstE2I + yLEstE2I;
    yEst  = yLEst*yLEstE2I + yREst*yREstE2I;
    if (yEstE>kEps) {
      yEstE = 1./yEstE; // err^2
      yEst  *= yEstE;   // weighted estimate
      yEstE = TMath::Sqrt(yEstE);
    }
    yEst -= y[ip]; // take difference
    //
    // tangent of the kink from 2 slopes
    if (slopRE && slopLE) {
      float tLRone  = 1.f+slopR*slopL;
      if (TMath::Abs(tLRone)>kEps) {
	float tLRoneI = 1.f/tLRone;
	float tL2one = 1.f+slopL*slopL;
	float tR2one = 1.f+slopR*slopR;
	tgKink  = (slopL-slopR)*tLRoneI;
	tgKinkE = slopLE*tR2one*tR2one + slopRE*tL2one*tL2one;
	tgKinkE = TMath::Sqrt(tgKinkE)*tLRoneI*tLRoneI;
      }
      else {
	tgKink = 1.0f;
	tgKinkE = 1.0f;
      }
    }

  }
  //
}

//____________________________________________________________________
void FitCircle(int np, const float* x, const float* y, 
	       float &xc, float &yc, float &r2, float* dy)
{
  // fit points to circle, if dy!=0, fill residuals
  double x0=0.,y0=0.;
  for (int i=np;i--;) {
    x0 += x[i];
    y0 += y[i];
  } 
  x0 /= np;
  y0 /= np;
  double su2=0,sv2=0,suv=0,su3=0,sv3=0,su2v=0,suv2=0;
  for (int i=np;i--;) {
    double ui = x[i]-x0, ui2 = ui*ui;
    double vi = y[i]-y0, vi2 = vi*vi;
    suv += ui*vi;
    su2 += ui2;
    sv2 += vi2;
    su3 += ui2*ui;
    sv3 += vi2*vi;
    su2v += ui2*vi;
    suv2 += ui*vi2;
  } 
  double rhsU = 0.5*(su3+suv2), rhsV = 0.5*(sv3+su2v);
  double det = su2*sv2-suv*suv;
  double uc  = (rhsU*sv2 - rhsV*suv)/det;
  double vc  = (su2*rhsV - suv*rhsU)/det;
  r2  = uc*uc + vc*vc + (su2+sv2)/np;
  xc = uc + x0;
  yc = vc + y0;
  //
  if (dy) {
    for (int i=np;i--;) {
      double dx = x[i]-xc;
      double dxr = r2 - dx*dx;
      double ys = dxr>0 ? TMath::Sqrt(dxr) : 0;
      double dy0 = y[i]-yc;
      double dysp = dy0-ys;
      double dysm = dy0+ys;
      dy[i] = TMath::Abs(dysp)<TMath::Abs(dysm) ? dysp : dysm;
    }
  }
}


//______________________________________________
void FixAlignmentBug(int sect, float q2pt, float bz,
		     float& alp, float& xCl, float &zTr, float &deltaY, float &deltaZ)
{
  // fix alignment bug: https://alice.its.cern.ch/jira/browse/ATO-339?focusedCommentId=170850&page=com.atlassian.jira.plugin.system.issuetabpanels:comment-tabpanel#comment-170850
  //
  static TGeoHMatrix *mCache[72] = {0};
  if (sect<0||sect>=72) {
    ::Error("FixAlignmentBug","Invalid sector %d",sect);
    return;
  }
  int lr = sect/36 ? (AliGeomManager::kTPC2) : (AliGeomManager::kTPC1);
  TGeoHMatrix* mgt = mCache[sect];
  if (!mgt) {
    int volID = AliGeomManager::LayerToVolUIDSafe(lr,sect%36);
    mgt = new TGeoHMatrix(*AliGeomManager::GetTracking2LocalMatrix(volID));
    mgt->MultiplyLeft(AliGeomManager::GetMatrix(volID));
    mCache[sect] = mgt;
    printf("Caching matrix for sector %d\n",sect);
  }  
  double alpSect = ((sect%18)+0.5)*20.*TMath::DegToRad();

  // cluster in its proper alpha frame with alignment bug. Note that the Track Z is supplied!!!
  double xyzClUse[3] = {xCl,0,zTr}; // this is what we read from the residual tree
  double xyzTrUse[3] = {xCl, deltaY, zTr}; // track in bad cluster frame
  //
  // recover cluster Z position by adding deltaZ
  xyzClUse[2] -= deltaZ;
  static AliExternalTrackParam trDummy;
  trDummy.Local2GlobalPosition(xyzClUse,alp); // misaligned cluster in global frame
  double xyz0[3]={xyzClUse[0],xyzClUse[1],xyzClUse[2]};
  mgt->MasterToLocal(xyz0,xyzClUse);
  // we got ideal cluster in the sector tracking frame
  //
  // go to ideal cluster frame
  trDummy.Local2GlobalPosition(xyzClUse,alpSect); // ideal global
  double alpFix = TMath::ATan2(xyzClUse[1],xyzClUse[0]);    // fixed cluster phi
  trDummy.Global2LocalPosition(xyzClUse,alpFix);     // fixed cluster in its frame
  //
  trDummy.Local2GlobalPosition(xyzTrUse,alp); // track in global frame
  trDummy.Global2LocalPosition(xyzTrUse,alpFix); // track in cluster frame
  alp = alpFix;
  //
  double dx = xyzTrUse[0] - xyzClUse[0]; // x might not be the same after alignment fix
  // deduce track slopes assuming it comes from the vertex
  double tgphi = tgpXY(xyzClUse[0],xyzTrUse[1],q2pt,bz);
  xyzTrUse[1] -= dx*tgphi;
  xyzTrUse[2] -= dx*xyzClUse[2]/xyzClUse[0]; // z2x
  //
  xCl = xyzClUse[0];
  zTr = xyzTrUse[2]; // we still use track Z as a reference ...
  deltaY = xyzTrUse[1]-xyzClUse[1];
  deltaZ = xyzTrUse[2]-xyzClUse[2];
  //
}

//_____________________________________________________
float tgpXY(float x, float y, float q2p, float bz)
{
  // get the tg of primary track inclination wrt x axis given
  // that it was registered at X,Y coordinates
  float c = q2p*bz*(-0.299792458e-3f);
  if (TMath::Abs(c)<1e-9) return y/x;
  float r2 = x*x+y*y;
  float det = 4./r2 - c*c;
  float snp  = 0;
  if (det<0) {
    snp = TMath::Sign(-0.8f,c);
    printf("track of q2p=%f cannot reach x:%f y:%f\n",q2p,x,y);
  }
  else {
    snp = 0.5f*(y*TMath::Sqrt(det)-c*x); // snp at vertex
    snp += x*c;  // snp at x,y
  }
  return snp/TMath::Sqrt((1.f-snp)*(1.f+snp));
}



//_______________________________________________
void MaxLeftRightDev(int np, const float* x, const float *y, const int nVoisin, 
		     float *diffY, float *diffLR, float *kink)
{
  // calculate the difference between the point and linear extrapolation from neigbourhood
  const float kEps = 1e-9;
  double sumX1b[kNPadRows+1],sumX2b[kNPadRows+1],sumYb[kNPadRows+1],sumYXb[kNPadRows+1];
  double *sumX1 = sumX1b+1, *sumX2 = sumX2b+1, *sumY0 = sumYb+1, *sumYX = sumYXb+1;
  //
  sumX1[-1]=0.f;
  sumX2[-1]=0.f;
  sumY0[-1]=0.f;
  sumYX[-1]=0.f;
  // accumulate matrix elements for whole array, element -1 is at 0
  double sX1=0., sX2=0., sY0=0., sYX=0.;

  for (int ip=0;ip<np;ip++) {
    sumX1[ip] = sX1 += x[ip];
    sumX2[ip] = sX2 += x[ip]*x[ip];
    sumY0[ip] = sY0 += y[ip];
    sumYX[ip] = sYX += y[ip]*x[ip];
    diffY[ip]  = 0.0f;
    diffLR[ip] = 0.0f;
    kink[ip]   = 0.0f;
  }
  // 

  for (int ip=0;ip<np;ip++) {

    float slopL=0,slopR=0;
    float yLEst=0.f,yREst=0.f;
    float &yEst = diffY[ip];
    float &estDiff = diffLR[ip];
    float &tgKink = kink[ip];
    int ip0,ip1;
    Bool_t okL=kFALSE,okR=kFALSE;
    //
    // estimate from the left
    ip0 = ip-nVoisin;
    if (ip0>-1) { 
      ip1 = ip-1;    // extract sum from ip0 to ip1 from cumulant, S00=nVoisin
      sX1 = sumX1[ip1] - sumX1[ip0-1]; // S01
      sX2 = sumX2[ip1] - sumX2[ip0-1]; // S11
      sY0 = sumY0[ip1] - sumY0[ip0-1]; // RHS0
      sYX = sumYX[ip1] - sumYX[ip0-1]; // RHS1
      double det = nVoisin*sX2 - sX1*sX1;
      if (det>kEps) {
	double detI = 1./det;
	// yLEst = offs + slop*x[ip] 
	// with offs=(sY0*sX2-sYX*sX1)/det and slop=(nVoisin*sYX-sY0*sX1)/det
	// inverse err^2 = 1/(errOffs+x^2*errSlop+2*x*errSlpOffs) with
	// errOffs = S11/det, errSlp=S00/det, errSlpOffs=-S01/det
	yLEst = ((sY0*sX2-sYX*sX1) + (nVoisin*sYX-sY0*sX1)*x[ip])*detI;
	slopL = (nVoisin*sYX-sY0*sX1)*detI;
	okL   = kTRUE;
      }
    }
    // estimate from the right
    ip1 = ip+nVoisin;
    if (ip1<np) {   // extract sum from ip+1 to ip1 from cumulant, S00=nVoisin
      sX1 = sumX1[ip1] - sumX1[ip]; // S01
      sX2 = sumX2[ip1] - sumX2[ip]; // S11
      sY0 = sumY0[ip1] - sumY0[ip]; // RHS0
      sYX = sumYX[ip1] - sumYX[ip]; // RHS1
      double det = nVoisin*sX2 - sX1*sX1;
      if (det>kEps) {
	double detI = 1./det;
	yREst    = ((sY0*sX2-sYX*sX1) + (nVoisin*sYX-sY0*sX1)*x[ip])*detI;
	slopR  = (nVoisin*sYX-sY0*sX1)*detI;
	okR   = kTRUE;
      }
    }
    //
    float dfL = okL ? y[ip] - yLEst : 0;
    float dfR = okR ? y[ip] - yREst : 0;
    //
    yEst = TMath::Abs(dfL)>TMath::Abs(dfR) ? dfL : dfR;
    if (okL && okR) {
      float tLRone  = 1.f+slopR*slopL;
      tgKink = TMath::Abs(tLRone)>kEps ? (slopL-slopR)/tLRone : 0;
      estDiff = yLEst - yREst;
    }
    //
    
  
  }
  //
}

int CheckResiduals(int np, const float *x, const float *y, const float *z, const int *sec36, Bool_t* kill, 
		   int nVois,float cut)
{

  int ip0=0,ip1;
  int sec0 = sec36[ip0];
  int npLast = np-1;
  //
  const int nMinAcc = 30;

  memset(yDiffMA,0,np*sizeof(float));
  memset(zDiffMA,0,np*sizeof(float));
  memset(yDiffLL,0,np*sizeof(float));
  memset(zDiffLL,0,np*sizeof(float));
  memset(kill,0,np*sizeof(Bool_t));
  float absDevY[kNPadRows],absDevZ[kNPadRows];
  for (int i=0;i<np;i++) {
    if (sec36[i]==sec0 && i<npLast) continue;
    //
    // sector change or end of input reached
    // run estimators for the points in the same sector
    int npSec = i-ip0;
    if (i==npLast) npSec++;
    //
    DiffToLocLine(npSec, x+ip0, y+ip0, nVois, yDiffLL+ip0);
    DiffToLocLine(npSec, x+ip0, z+ip0, nVois, zDiffLL+ip0);
    DiffToMA(npSec, x+ip0, y+ip0, nVois, yDiffMA+ip0);
    DiffToMA(npSec, x+ip0, z+ip0, nVois, zDiffMA+ip0);
    //
    ip0 = i;
    sec0 = sec36[ip0];
  }
  // store abs deviations
  int naccY=0,naccZ=0;
  for (int i=np;i--;) {
    //    if (yDiffMA[i]) absDevY[naccY++] = TMath::Abs(yDiffMA[i]);
    //    if (zDiffMA[i]) absDevZ[naccZ++] = TMath::Abs(zDiffMA[i]);
    if (yDiffLL[i]) absDevY[naccY++] = TMath::Abs(yDiffLL[i]);
    if (zDiffLL[i]) absDevZ[naccZ++] = TMath::Abs(zDiffLL[i]);
  }
  //
  // estimate rms on 90% smallest deviations
  int kmnY = 0.9*naccY,kmnZ = 0.9*naccZ;
  if (naccY<nMinAcc || naccZ<nMinAcc) { // kill all
    for (int i=np;i--;) kill[i] = kTRUE;
    return np;
  }

  SelKthMin(kmnY, naccY, absDevY);
  SelKthMin(kmnZ, naccZ, absDevZ);
  float rmsKY=0,rmsKZ=0;
  for (int i=kmnY;i--;) rmsKY += absDevY[i]*absDevY[i];
  for (int i=kmnZ;i--;) rmsKZ += absDevZ[i]*absDevZ[i];
  rmsKY = TMath::Sqrt(rmsKY/kmnY);
  rmsKZ = TMath::Sqrt(rmsKZ/kmnZ);
  printf("RMSY %d min of %d: %f | RMSZ %d min of %d: %f\n",kmnY,naccY,rmsKY, kmnZ,naccZ,rmsKZ);
  //
  //
  int nKill = 0;
  for (int ip=0;ip<np;ip++) {

    yDiffLL[ip] /= rmsKY;
    zDiffLL[ip] /= rmsKZ;
    yDiffMA[ip] /= rmsKY;
    zDiffMA[ip] /= rmsKZ;
    float dy = yDiffLL[ip], dz = zDiffLL[ip];
    if (dy*dy+dz*dz>cut) {
      kill[ip] = kTRUE;
      nKill++;
    }
  }
  return nKill;
  //
}

/*
void CheckResiduals(int np, const float *x, const float *y, const float *z, const int *sec36)
{

  int segmStart[kNSect],segmNClus[kNSect],nAccSegmY[kNSect],nAccSegmZ[kNSect],nSegments=0;
  int ip0=0,ip1;
  int sec0 = sec36[ip0];
  int npLast = np-1;
  //
  float absDevY[kNPadRows],absDevZ[kNPadRows];
  int naccADY=0,naccADZ=0;

  
  for (int i=0;i<np;i++) {
    if (sec36[i]==sec0 && i<npLast) continue;
    //
    // sector change or end of input reached
    // run estimators for the points in the same sector
    int npSec = i-ip0;
    if (i==npLast) npSec++;
    //
    // register track segment
    segmStart[nSegments] = ip0;
    segmNClus[nSegments] = npSec;
    //
    DiffToMA(npSec, x+ip0, y+ip0, 5, yDiffMA+ip0);
    DiffToMA(npSec, x+ip0, z+ip0, 5, zDiffMA+ip0);
    //
    float* yDiffLLp = yDiffLL+ip0;
    float* zDiffLLp = zDiffLL+ip0;
    nAccSegmY[nSegments] = DiffToLocLine(npSec, x+ip0, y+ip0, 4, yDiffLLp);
    nAccSegmZ[nSegments] = DiffToLocLine(npSec, x+ip0, z+ip0, 4, zDiffLLp);
    nSegments++;
    //
    // store abs deviations
    for (int ia=npSec;ia--;) {
      float valy = TMath::Abs(yDiffLLp[ia]);
      float valz = TMath::Abs(zDiffLLp[ia]);
      if (valy) absDevY[naccADY++] = valy;
      if (valz) absDevZ[naccADZ++] = valz;
    }
    //
    int naccyM = DiffToMedLine(npSec, x+ip0, y+ip0, 5, yDiffML+ip0);
    int nacczM = DiffToMedLine(npSec, x+ip0, z+ip0, 5, zDiffML+ip0);

    ip0 = i;
    sec0 = sec36[ip0];
  }
  //
  float meanY=0.f,meanZ=0.f,rmsY=0.f,rmsZ=0.f;
  //
  int naccY=0,naccZ=0;
  for (int ip=np;ip--;) {
    float devY = yDiffMA[ip];
    float devZ = zDiffMA[ip];   
    if (devY) {
      meanY += devY;
      rmsY  += devY*devY;
      naccY++;
    }
    if (devZ) {
      meanZ += devZ;
      rmsZ  += devZ*devZ;
      naccZ++;
    }
    //
  }
  if (naccY) {
    meanY /= naccY;
    rmsY = TMath::Sqrt(rmsY/naccY - meanY*meanY);
  }
  if (naccZ) {
    meanZ /= naccZ;
    rmsZ = TMath::Sqrt(rmsZ/naccZ - meanZ*meanZ);
  }
  printf("Y: %d %.3f  Z: %d %.3f\n",naccY,rmsY, naccZ,rmsZ);
  //
  int kmnY = 0.9*naccADY;
  int kmnZ = 0.9*naccADZ;
  SelKthMin(kmnY, naccADY, absDevY);
  float rmsKY=0,rmsKZ=0;
  for (int i=kmnY;i--;) rmsKY += absDevY[i]*absDevY[i];
  for (int i=kmnZ;i--;) rmsKZ += absDevZ[i]*absDevZ[i];
  rmsKY = TMath::Sqrt(rmsKY/kmnY);
  rmsKZ = TMath::Sqrt(rmsKZ/kmnZ);
  printf("RMSY %d min of %d: %f | RMSZ %d min of %d: %f\n",kmnY,naccADY,rmsKY, kmnZ,naccADZ,rmsKZ);
  for (int i=np;i--;) {
    yDiffLL[i]/=rmsKY;
    zDiffLL[i]/=rmsKZ;
  }
  //
  Bool_t killed[kNPadRows] = {0};
  float xmed[kNPadRows],ymed[kNPadRows];
  int nVois = 5;
  for (int ip=0;ip<np;ip++) {
    if (abs(yDiffLL[ip])<3) continue;
    //
    // check large value with median filter
    int nacc=0,naccR=0,naccL=0,il=ip-1;
    while (naccL<nVois && il>-1 && !killed[il]) {
      xmed[naccL] = x[il];
      ymed[naccL] = y[il];
      il--;
      naccL++;
    }
    il = ip+1;
    nacc = naccL;
    while (naccR<nVois && il<np) {
      xmed[nacc] = x[il];
      ymed[nacc] = y[il];
      il++;
      naccR++;
      nacc++;
    }
    nacc = naccL + naccR;
    float offs=0,slp=0;
    if (nacc<nVois) {killed[ip] = kTRUE; continue;}
    medFit(nacc, xmed, ymed, offs,slp);
    float dev = (y[ip]-(offs+slp*x[ip]))/rmsKY;
    if (TMath::Abs(dev)>4.) killed[ip] = kTRUE;
    yDiffML[ip] = dev;
  }
  //
}
*/

//_________________________________________
void DiffToMA(int np, const float* x, const float *y, const int winLR, float* diffMA)
{
  // difference to moving average, excluding central element
  //
  double arrSumO[kNPadRows+1], *arrSum=arrSumO+1;
  arrSum[-1] = 0.;
  for (int ip=0;ip<np;ip++) arrSum[ip] = arrSum[ip-1]+y[ip];
  for (int ip=0;ip<np;ip++) {
    diffMA[ip] = 0;
    int ipmn = ip-winLR;
    int ipmx = ip+winLR;
    if (ipmn<0)   ipmn=0;
    if (ipmx>=np) ipmx=np-1;
    int nrm = (ipmx-ipmn);
    if (nrm<winLR) continue;
    double ma = (arrSum[ipmx]-arrSum[ipmn-1] - (arrSum[ip]-arrSum[ip-1]))/nrm;
    diffMA[ip] = y[ip] - ma;
  }

}


//_______________________________________________
int DiffToLocLine(int np, const float* x, const float *y, const int nVoisin, float *diffY)
{
  // calculate the difference between the point and linear extrapolation from neigbourhood
  const float kEps = 1e-9;
  double sumX1b[kNPadRows+1],sumX2b[kNPadRows+1],sumYb[kNPadRows+1],sumYXb[kNPadRows+1];
  double *sumX1 = sumX1b+1, *sumX2 = sumX2b+1, *sumY0 = sumYb+1, *sumYX = sumYXb+1;
  //
  sumX1[-1]=0.f;
  sumX2[-1]=0.f;
  sumY0[-1]=0.f;
  sumYX[-1]=0.f;
  // accumulate matrix elements for whole array, element -1 is at 0
  double sX1=0., sX2=0., sY0=0., sYX=0.;

  for (int ip=0;ip<np;ip++) {
    sumX1[ip] = sX1 += x[ip];
    sumX2[ip] = sX2 += x[ip]*x[ip];
    sumY0[ip] = sY0 += y[ip];
    sumYX[ip] = sYX += y[ip]*x[ip];
    diffY[ip]  = 0.0f;
  }
  // 
  int nAcc = 0;
  for (int ip=0;ip<np;ip++) {
    float &yEst = diffY[ip];
    // estimate from the left
    int ip0 = ip-nVoisin;
    int ip1 = ip+nVoisin;
    if (ip0<0)   ip0=0;
    if (ip1>=np) ip1=np-1;
    int nrm = (ip1-ip0);
    if (nrm<nVoisin) continue;
    int ip0m = ip0-1, ipm = ip-1;
    // extract sum from ip0 to ip1 from cumulant, S00=nrm, excluding tested point
    sX1 = sumX1[ip1] - sumX1[ip0m] - (sumX1[ip]-sumX1[ipm]); // S01
    sX2 = sumX2[ip1] - sumX2[ip0m] - (sumX2[ip]-sumX2[ipm]); // S11
    sY0 = sumY0[ip1] - sumY0[ip0m] - (sumY0[ip]-sumY0[ipm]); // RHS0
    sYX = sumYX[ip1] - sumYX[ip0m] - (sumYX[ip]-sumYX[ipm]); // RHS1
    double det = nrm*sX2 - sX1*sX1;
    if (det<kEps) continue;
    double detI = 1./det;
    // yLEst = offs + slop*x[ip] 
    // with offs=(sY0*sX2-sYX*sX1)/det and slop=(nVoisin*sYX-sY0*sX1)/det
    // inverse err^2 = 1/(errOffs+x^2*errSlop+2*x*errSlpOffs) with
    // errOffs = S11/det, errSlp=S00/det, errSlpOffs=-S01/det
    yEst = y[ip]-((sY0*sX2-sYX*sX1) + (nrm*sYX-sY0*sX1)*x[ip])*detI;
    nAcc++;
    //    
  }
  return nAcc;
  //
}

//====================================================================
int DiffToMedLine(int np, const float* x, const float *y, const int nVoisin, float *diffY)
{
  int nAcc = 0;
  float offs=0,slp=0;
  float buff[kNPadRows];
  for (int ip=0;ip<np;ip++) {
    float &yEst = diffY[ip];
    yEst = 0;
    // estimate from the left
    int ip0 = ip-nVoisin;
    int ip1 = ip+nVoisin;
    if (ip0<0)   ip0=0;
    if (ip1>=np) ip1=np-1;
    int nrm = (ip1-ip0+1);
    if (nrm<nVoisin) continue;
    int ip0m = ip0-1, ipm = ip-1;
    const float *arrx = x+ip0, *arry = y+ip0;
    medFit(nrm, arrx, arry, offs,slp);
    /*
    float asum = 0;
    for (int i=nrm;i--;) {
      buff[i] = arry[i] - (offs + slp*arrx[i]);
      if (i+ip0!=ip) asum += TMath::Abs(buff[i]);
    }
    asum /= nrm-1;
    yEst = buff[ip-ip0]/asum;
    */
    yEst = y[ip] - (offs + slp*x[ip]);
    nAcc++;
  }
  return nAcc;
}

//___________________________________________________________________
void medFit(int np, const float* x, const float* y, float &a, float &b, float delI)
{
  // Median linear fit: minimizes abs residuals instead of squared ones
  // Adapted from "Numerical Recipes in C"
  float aa,bb,b1,b2,f,f1,f2,sigb,chisq=0.0f;
  if (!delI) {
    float sx=0.0f,sxx=0.0f,sy=0.0f,sxy=0.0f,del;
    //
    for (int j=np;j--;) { sx += x[j]; sxx += x[j]*x[j];}
    del = np*sxx-sx*sx;
    //
    for (int j=np;j--;) { sy += y[j]; sxy += x[j]*y[j];}
    //
    float delI = 1./del;
    aa = (sxx*sy-sx*sxy)*delI;
    bb = (np*sxy-sx*sy)*delI;
  }
  else { // initial values provided
    aa = a;
    bb = b;
  }
  //
  for (int j=np;j--;) {
    float temp = y[j]-(aa+bb*x[j]);
    chisq += temp*temp;
  }
  //
  sigb = TMath::Sqrt(chisq*delI);
  b1=bb;
  f1 = RoFunc(np,x,y,b1,aa);
  if (sigb>0) {
    b2 = bb+TMath::Sign(float(3.0f*sigb),f1);
    f2 = RoFunc(np,x,y,b2,aa);
    if (f1==f2) {
      a = aa;
      b = bb;
      return;
    }
    while (f1*f2 > 0.0f) { // bracketing
      bb = b2 + 1.6f*(b2-b1);
      b1 = b2;
      f1 = f2;
      b2 = bb;
      f2 = RoFunc(np,x,y,b2,aa);
    }
    sigb = 0.01*sigb;
    while (fabs(b2-b1)>sigb) {
      bb = b1 + 0.5f*(b2-b1);
      if (bb==b1 || bb==b2) break;
      f = RoFunc(np,x,y,bb,aa);
      if (f*f1 >= 0.0f) {
	f1=f;
	b1=bb;
      } 
      else {
	f2 = f;
	b2 = bb;
      }
    }
  }
  a = aa;
  b = bb;
  //
}

//___________________________________________________________________
float RoFunc(int np, const float* x, const float* y, float b, float &aa)
{
  const float kEPS = 1.0e-7f; 
  static float* arrTmp = 0;
  static int nBook = 0;
  if (np>nBook) { // make sure the buffer is ok
    nBook = np;
    delete[] arrTmp;
    arrTmp = new float[nBook];
  }
  float d,sum=0.0f;
  for (int j=np;j--;) arrTmp[j] = y[j]-b*x[j];
  //
  int nph = np>>1;
  if (np<20) {  // it is faster to do insertion sort 
    for (int i=1;i<np;i++) {
      float v = arrTmp[i];
      int j;
      for (j=i;j--;) if (arrTmp[j]>v) arrTmp[j+1]=arrTmp[j]; else break;
      arrTmp[j+1] = v;
    }
    aa = (np&0x1) ? arrTmp[nph] : 0.5f*(arrTmp[nph-1]+arrTmp[nph]);
  }
  else {
    aa = (np&0x1) ? SelKthMin(nph,np,arrTmp) : 
      0.5f*(SelKthMin(nph-1,np,arrTmp)+SelKthMin(nph,np,arrTmp));
  }
  for (int j=np;j--;) {
    d = y[j] - (b*x[j] + aa);
    if (y[j] != 0.0f) d /= TMath::Abs(y[j]);
    if (TMath::Abs(d) > kEPS) sum += (d >= 0.0f ? x[j] : -x[j]);
  }
  return sum;
}

//___________________________________________________________________
Float_t SelKthMin(int k, int np, float* arr)
{
  // Returns the k th smallest value in the array. The input array will be rearranged
  // to have this value in location arr[k] , with all smaller elements moved before it
  // (in arbitrary order) and all larger elements after (also in arbitrary order).
  // From Numerical Recipes in C++

  int i,ir,j,l,mid;
  float a;
  l=0;ir=np-1;
  for (;;) {
    if (ir<=l+1) {
      if (ir==l+1 && arr[ir]<arr[l]) swap(arr[l],arr[ir]);
      return arr[k];
    } 
    else {
      int mid = (l+ir)>>1, i=l+1;
      swap(arr[mid],arr[i]);
      if (arr[i]>arr[ir]) swap(arr[i],arr[ir]);
      if (arr[l]>arr[ir]) swap(arr[l]  ,arr[ir]);
      if (arr[i]>arr[l])  swap(arr[i],arr[l]);
      j=ir;
      a=arr[l];
      for (;;) {
	do i++; while (arr[i]<a);
	do j--; while (arr[j]>a);
	if (j<i) break;
	swap(arr[i],arr[j]);
      }
      arr[l]=arr[j];
      arr[j]=a;
      if (j>=k) ir=j-1;
      if (j<=k) l=i;
    }
  }
}
//=========================================================

//_________________________________________
void LoadVDrift()
{
  // load vdrift params
  fVDriftGraph = 0;
  fVDriftParam = 0;
  TFile *fdrift = TFile::Open(Form("%s.root",kDriftFileName));
  if (fdrift) {
    TTree * tree = (TTree*)fdrift->Get("fitTimeStat");
    if (tree==NULL) {
      ::Error("LoadDriftCalibration FAILED", "tree fitTimeStat not avaliable in file %s.root",kDriftFileName);
    }
    else {      
      tree->SetBranchAddress("grTRDReg.",&fVDriftGraph);
      tree->SetBranchAddress("paramRobust.",&fVDriftParam);
      tree->GetEntry(0);
      if (fVDriftGraph==NULL || fVDriftGraph->GetN()<=0) {
	::Info("LoadDriftCalibration FAILED", "ITS/TRD drift calibration not availalble. Trying ITS/TOF");
	tree->SetBranchAddress("grTOFReg.",&fVDriftGraph);
	tree->GetEntry(0);
      }
      /*
      else {
	::Info("LoadDriftCalibration", "tree fitTimeStat not avaliable in file %s.root",kDriftFileName);
      }
      */
    }
    delete tree;
  }
  else {
    ::Error("LoadDriftCalibration FAILED", "fitDrift.root not present");
  }
  if (fdrift) fdrift->Close();
  delete fdrift;
}

//_________________________________________
float GetDriftCorrection(float z, float x, float phi, int rocID)
{
  // apply vdrift correction
  int side = ((rocID/kNSect)&0x1) ? -1:1; // C:A
  float drift = side>0 ? kZLim-z : z+kZLim;
  float gy    = TMath::Sin(phi)*x;
  Double_t pvecFit[3];
  pvecFit[0]= side;             // z shift (cm)
  pvecFit[1]= drift*gy/kZLim;   // global y gradient
  pvecFit[2]= drift;            // drift length
  float expected = (fVDriftParam==NULL) ? 0:
    (*fVDriftParam)[0]+
    (*fVDriftParam)[1]*pvecFit[0]+
    (*fVDriftParam)[2]*pvecFit[1]+
    (*fVDriftParam)[3]*pvecFit[2];
  return -side*(expected+fCorrTime*drift);
  //
}
