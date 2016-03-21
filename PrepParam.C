#include <TString.h>
#include <TFile.h>
#include <TSystem.h>
#include <THn.h>
#include <TH2.h>
#include <TROOT.h>
#include <TStopwatch.h>
#include <TKey.h>
#include <TStyle.h>
#include <TCanvas.h>

#include "AliCDBManager.h"
#include "AliCDBEntry.h"
#include "AliGRPManager.h"
#include "AliGRPObject.h"
#include "AliNDLocalRegression.h"
#include "AliTPCcalibDB.h"
#include "AliTPCRecoParam.h"
#include "AliTPCTransform.h"
#include "AliTPCParam.h"
#include "AliMagF.h"
#include "AliTPCChebCorr.h"
#include "AliTPCParamSR.h"
#include "AliTPCCorrection.h"
#include "AliTPCExB.h"
#include <TGeoGlobalMagField.h>

//TString dataPath = "/data/testTPC/rawparam/000245231/Time1448614500";
//TString dataPath = "/data/testTPC/rawparam/000245231/Time1448613598";
//TString dataPath = "/data/testTPC/ref244918/mapsNDL/TimeAll";
//TString dataPath = "/data/testTPC/ref245231/mapsNDL/TimeAll";   // !!! where to look for hisX directories !!!
TString dataPath = "/data/testTPC/ref244918/mapsNDL/MI_Marh11/Time1448447533";   // !!! where to look for hisX directories !!!

//TString dataPath = "/data/testTPC/rawparam/%09d/fit";
TString mapTemplate = "delta%s%sDistFit%s_sec%d_%d_theta%d_%d_SmoothConst1;1";  // !!! parameterization to use !!!
//TString mapTemplate = "delta%s%sDistFit%s_sec%d_%d_theta%d_%d;1";
TString fitVersion = "0";                                                       // !!! fit version to use !!!!

//TString mapTemplate = "delta%s%sDistFit%s_sec%d_%d_theta%d_%d_SmoothConst3;1";
//TString fitVersion = "1";

//TString mapTemplate = "delta%s%sDistFit%s_sec%d_%d_theta%d_%d;1";
TString fileTemplate = "delta%s%sDistFit%s_sec%d_%d_theta%d_%d";               /// !!! template of file in hisX dir. !!!

enum {kDTRD,kDTOF};
enum {kNInnerSectors=36, kInnerNRow=63};

//TString dimS[2] = {"RPhi","RPhi"};
TString dimS[2] = {"RPhi","Z"};
enum {kndResRPhi,kndResZ2R};
Int_t useHistID[2] = {1,4};
TString detPairS[2] = {"TPCITSTRD","TPCITSTOF"};
int scCov = 18; //9; // number of sectors per file !!! grouping of sectors !!!
int curRun = 233678, curDetComb = 0;


const int kNQ2PtQuery = 3;
const float kQ2PtQuery[kNQ2PtQuery] = {-2.,0.,2.};
const float kq2ptCent = 0; // q/pt for the central value of residual
//
const Bool_t corrX = kFALSE;//kTRUE;

const Bool_t hackZ2RtoZ2X = kFALSE; // set to true if NDLocal used R instead of X
//
Bool_t useIniTransform = kFALSE;
Bool_t usedZTrack = kTRUE; // were residuals evaluated at Ztrack instead of Zcluster?
const double inversionEps = 20e-4; // when inverting, stop Newton-Raphson iterations at this eps
const int    inversionMaxIt = 3; // when inverting, stop Newton-Raphson after some numbers of iterations
const double kSectWindow = TMath::Pi()/9.;


double bZ=0;
UInt_t timeMin=0,timeMax=0,fCurrentTimeStamp=0;
const int kRecoParamID = 0;
AliTPCRecoParam* fCurrentRecoParam = 0;

AliNDLocalRegression *regMap[2] = {0,0};
Bool_t statusMap[18*2] = {kFALSE};
enum {kDimR,kDimSect,kDimZ2R,kDimQ2Pt,kNDimReg};
Double_t bounds[kNDimReg][2];
enum {kndR, kndSec, kndZ2R, kndQ2PT};
TFile* treeOutFile = 0;
TTree* treeOut = 0;
TCanvas* cnv=0;

Bool_t LoadDistMap(int run, int sect, Bool_t zPos, Int_t detComb);
Bool_t GetDistortion(int sector, int row, float y2x, float z2x, float xyz[3]);
void SetRun(int run);
void SetupTPCCalibDB(int run);
double GetTgPhi(double x, double y2x, double q2p, double b);
void Transform(int sector,int row, double x[3]);
void EvalNDLocal(double *vecInp,double &drp,double &dz);

void TrainCorr(int row, float* tzLoc, float* corrLoc);
AliTPCChebCorr* MakeCheb(const char* name, int nz=1, int nsc=2, Bool_t useS=kTRUE);
void extTreeTPC(AliTPCChebCorr* par, int nby=10,int nbz=200,const char* treeFile="corrTree.root", Bool_t update=kFALSE);
void BookTree(const char* treeFile, Bool_t update=kFALSE, const char* treeName="dtr", const char* treeTitle="SummaryTree");
void CloseTree();
void BenchChTPC(AliTPCChebCorr* chp, int n=100000);
//void BenchTPCCheb(AliTPCChebCorr* chp,int nby=50,int nbz=125,Bool_t benchTree = kFALSE); /// !!! if benchTree is true, large QA tree is produced !!!
void BenchTPCCheb(AliTPCChebCorr* chp,int nby=50,int nbz=125,Bool_t benchTree = kTRUE); /// !!! if benchTree is true, large QA tree is produced !!!

// Precision of Cheb param for each dimension
float precD[3]={100e-4,100e-4,100e-4};
//float precD[3]={50e-4,50e-4,50e-4};
//
// number of nodes to parameterize
//  int np[3][2] = {{10,10},{10,10},{10,10}};
int np[3][2] = {{15,15},{15,15},{15,15}};
//int np[3][2] = {{25,25},{25,25},{25,25}};
////int np[3][2] = {{20,15},{20,15},{20,15}};
//int np[3][2] = {{60,60},{60,60},{60,60}};
typedef struct {
  Float_t rL[3]; // sector coordinates
  Float_t dL[3]; // sector frame residuals, original
  Float_t dC[3]; // sector frame residuals from Cheb[Short] param
  Float_t tgphi; // tangent of track az.angle at padrow
  UChar_t sect;   // sector
  UChar_t row;    // row
  Short_t iz;     // Z bin
  Short_t iy2x;   // sector Y/X bin
  //
} Dist_t;

Dist_t dst;
double tgAngCent=0;

//________________________________________________________________
AliTPCChebCorr* PrepParam(Int_t run=245231, Bool_t useTransform=kFALSE)
{
  useIniTransform = useTransform;

  SetRun(run);
  // 
  // modify recoparam according to what we want to parameterize
  AliTPCChebCorr* cc = MakeCheb(Form("cheb%d",run));
  BenchTPCCheb(cc);
  return cc;
}

//________________________________________________________________
AliTPCChebCorr* MakeCheb(const char* name, int nz, int nsc, Bool_t useS)
{
  //
  AliTPCcalibDB* calib = AliTPCcalibDB::Instance();
  AliTPCParam* param   = calib->GetParameters();
  //
  float zMax =  1.0;
  //
  AliTPCChebCorr* chp = new AliTPCChebCorr(name,name,nsc,nz,zMax);
  chp->SetUseFloatPrec(!useS);
  chp->Parameterize(TrainCorr,3,np,precD);
  chp->SetTimeStampStart(timeMin);
  chp->SetTimeStampEnd(timeMax);
  chp->SetTimeDependent(kFALSE);
  chp->SetUseZ2R(kTRUE);
  return chp;
}

//_____________________________________________________________________
void TrainCorr(int row, float* tzLoc, float* corrLoc)
{
  // compute correction for the point
  //
  // xtzLoc: y2x z or z2x in sector frame
  // corrLoc: x, y, z corrections in sector frame
  static int sector=0;
  static int side = 0;
  if (!tzLoc || !corrLoc) {
    sector = row%18; 
    side   = (row/18)&0x1;
    printf("Set to Sector %d Side: %d\n",sector,side);
    return;
  }
  //
  Bool_t res = GetDistortion(sector, row, tzLoc[0], tzLoc[1], corrLoc);
  if (!res) {
    printf("Failed to evaluate original distortion\n");
    exit(1);
  }
  //
}
//_______________________________________________________________________
void extTreeTPC(AliTPCChebCorr* par, int nby,int nbz,const char* treeFile, Bool_t update)
{
  printf("Extracting corrections for full tpc %s\n",treeFile);
  //
  BookTree(treeFile,update);
  //
  AliTPCcalibDB* calib = AliTPCcalibDB::Instance();
  AliTPCParam* tpcParam   = calib->GetParameters();
  //
  float zmin = par->GetZMin();//
  float zmax = par->GetZMax();//
  //
  float spanY2X = par->GetMaxY2X()*2;
  int nrowL = tpcParam->GetNRowLow();
  int nrowU = tpcParam->GetNRowUp();
  int nrow  = nrowL + nrowU;
  for (int i=3;i--;) dst.dC[i] = dst.dL[i] = 0;
  //
  Float_t distOr[3],distPar[3];
  TStopwatch sw;
  sw.Start();
  for (int iz=0;iz<=nbz;iz++) {
    dst.iz = iz;
    for (int isc=0;isc<18;isc++) {
      dst.sect = isc;    
      for (int ir=0;ir<nrow;ir++) {
	double xL = tpcParam->GetPadRowRadii(ir<nrowL ? 0 : 36 ,ir<nrowL ? ir : ir-nrowL);
	dst.row = ir;
	//
	for (int iy=0;iy<=nby;iy++) {
	  double yL =  (iy*spanY2X/nby - spanY2X/2.)*xL;
	  dst.iy2x = iy;
	  // go to sector/row in 0-71 convention
	  int sect72 = isc;
	  int row72  = ir;
	  if (ir>=nrowL) {
	    sect72 += 36; 
	    row72 -= nrowL;
	  }
	  //
	  double zL = zmin + iz*(zmax-zmin)/nbz;  // sector local measured
	  if (par->GetUseZ2R()) zL *= xL;
	  if (zL<0) sect72 += 18;
	  //
	  dst.rL[0] = xL; dst.rL[1] = yL; dst.rL[2] = zL;
	  //
	  Float_t tz[2] = {yL/xL, zL/xL};
	  //printf("%+f %+f %+f | %+e %+e\n",xL,yL,zL,tz[0],tz[1]);
	  if (!GetDistortion(isc, ir, tz[0], tz[1], distOr)) {
	    printf("Failed to evaluate original distortion\n");
	    exit(1);
	  }
	  dst.tgphi = tgAngCent;
	  //
	  par->Eval(sect72,ir>=kInnerNRow ? ir-kInnerNRow : ir, tz,distPar);
	  for (int i=3;i--;) {
	    dst.dC[i] = distPar[i];
	    dst.dL[i] = distOr[i];
	  }
	  treeOut->Fill();
	}
      }
    }
  }
  sw.Stop();
  sw.Print();
  //
  CloseTree();
  //
}

//_______________________________________________
void BookTree(const char* treeFile, Bool_t update, const char* treeName, const char* treeTitle)
{
  treeOutFile = TFile::Open(treeFile,update ? "update" : "recreate");
  treeOut = new TTree(treeName,treeTitle);
  //
  treeOut->Branch("rL", &dst.rL ,"rL[3]/F");
  treeOut->Branch("dL", &dst.dL ,"dL[3]/F");
  treeOut->Branch("dC", &dst.dC ,"dC[3]/F");
  treeOut->Branch("tgphi", &dst.tgphi ,"tgphi/F");
  treeOut->Branch("sect", &dst.sect ,"sect/b");  
  treeOut->Branch("row", &dst.row ,"row/b");  
  treeOut->Branch("iz", &dst.iz ,"iz/b");  
  treeOut->Branch("iy2x", &dst.iy2x ,"iy2x/b");  

  treeOut->SetAlias("sector","sect+iy2x/51.");
  treeOut->SetAlias("dX","dL[0]");
  treeOut->SetAlias("dY","dL[1]");
  treeOut->SetAlias("dZ","dL[2]");
  treeOut->SetAlias("dXch","dC[0]");
  treeOut->SetAlias("dYch","dC[1]");
  treeOut->SetAlias("dZch","dC[2]");
  treeOut->SetAlias("z2x","rL[2]/rL[0]");
  treeOut->SetAlias("y2x","rL[1]/rL[0]");
  //
}

//_______________________________________________
void CloseTree()
{
  treeOutFile->cd();
  treeOut->Write();
  delete treeOut;
  treeOutFile->Close();
  delete treeOutFile;
}

//_______________________________________________
void BenchChTPC(AliTPCChebCorr* chp, int n)
{
  float zmin = chp->GetZMin();
  float zmax = chp->GetZMax();
  float x[2],res[3];
  float y2xmin = -chp->GetMaxY2X(), y2xspan = 2*chp->GetMaxY2X();
  TStopwatch sw;
  sw.Start();
  for (int i=n;i--;) {
    int sect  = gRandom->Integer(72);
    int row   = gRandom->Integer(sect<kNInnerSectors ? kInnerNRow : 159-kInnerNRow);
    x[0] = y2xmin + y2xspan*gRandom->Rndm();
    x[1] = zmin + (zmax-zmin)*gRandom->Rndm();
    chp->Eval(sect,row,x,res);
  }
  sw.Stop();
  printf("TimePerCall: %e\n",sw.CpuTime()/n);
  sw.Print();
  //
}

//__________________________________________________________________________
void BenchTPCCheb(AliTPCChebCorr* chTPC, int nby,int nbz, Bool_t benchTree)
{
  //
  TFile* dumpf = 0;
  dumpf = TFile::Open("dumpChebTPC.root","recreate");
  chTPC->Write();
  dumpf->Close();
  delete dumpf;
  dumpf = TFile::Open("dumpChebTPC.root");
  TKey* kch = dumpf->GetKey(chTPC->GetName());
  //
  printf("\nBenchmark for TPC Cheb param\n");
  BenchChTPC(chTPC);
  printf("Size: Mem: %.3e Disk: %.3e kB\n",kch->GetObjlen()/1024.,kch->GetNbytes()/1024.);
  //
  dumpf->Close();
  delete dumpf;
  //
  printf("\n");
  if (!benchTree) return;
  extTreeTPC(chTPC,nby,nbz,"dumpChebTPC.root",kTRUE);
  //
  dumpf = TFile::Open("dumpChebTPC.root");
  TTree* tr = (TTree*)dumpf->Get("dtr");
  //
  gStyle->SetTitleW(0.3);
  cnv = new TCanvas("cnv","AliTPCChebCorr Benchmark",1400,1000);
  cnv->Divide(3,2);
  TH2* htemp;

  printf("plotting residuals between original and Cheb2D parameterizations\n");

  const char* kAxName[3]={"X","Y","Z"};
  for (int i=0;i<3;i++) {
    cnv->cd(i+1);
    gPad->SetLeftMargin(0.17);
    gPad->SetBottomMargin(0.15);
    tr->SetMarkerColor(kBlue);
    tr->Draw(Form("dL[%i]-dC[%i] : dL[%i]",i,i,i));
    htemp = (TH2F*)gPad->GetPrimitive("htemp");
    htemp->SetXTitle(Form("#Delta%s_{orig}, cm",kAxName[i]));
    htemp->SetYTitle(Form("#Delta%s_{orig}-#Delta%s_{param}, cm",kAxName[i],kAxName[i]));
    htemp->SetTitle(Form("#Delta%s_{orig}-#Delta%s_{param}, cm",kAxName[i],kAxName[i]));
    htemp->GetYaxis()->SetTitleSize(0.05);
    htemp->GetYaxis()->SetTitleOffset(1.4);
    htemp->GetXaxis()->SetTitleSize(0.05);
    htemp->GetXaxis()->SetTitleOffset(1.2);
    gPad->SetGrid();
    gPad->Modified();
    gPad->Update();
    //
  }
  //
  for (int i=0;i<3;i++) {
    cnv->cd(i+1+3);
    gPad->SetLeftMargin(0.17);
    gPad->SetBottomMargin(0.15);
    tr->SetLineColor(kRed);
    tr->Draw(Form("dL[%i]-dC[%i]",i,i));
    htemp = (TH2F*)gPad->GetPrimitive("htemp");
    htemp->SetXTitle(Form("#Delta%s_{orig}-#Delta%s_{param}, cm",kAxName[i],kAxName[i]));
    htemp->SetTitle(Form("#Delta%s_{orig}-#Delta%s_{param}, cm",kAxName[i],kAxName[i]));
    htemp->GetYaxis()->SetTitleSize(0.05);
    htemp->GetYaxis()->SetTitleOffset(1.4);
    htemp->GetXaxis()->SetTitleSize(0.05);
    htemp->GetXaxis()->SetTitleOffset(1.2);
    gPad->SetLogy();
    gPad->SetGrid();
    gPad->Modified();
    gPad->Update();
    //
  }
  //
  cnv->Print("benchTPCCheb.gif");
}

//________________________________________________________________
Bool_t LoadDistMap(int run, int sect, Bool_t zPos, Int_t detComb)
{
  static int lastRun = -1;
  static int lastDet = -1;
  static Bool_t firstRep[2*18] = {0};
  int sec0 = (sect/scCov)*scCov;
  int sec1 = sec0+scCov;
  if (run==lastRun && statusMap[sect+(zPos?0:18)] && detComb==lastDet) return kTRUE;  
  int th0 = zPos ? 0:-1;
  int th1 = zPos ? 1:0;
  const char* det = detPairS[detComb].Data();
  //
  lastRun = run;
  lastDet = detComb;
  memset(statusMap,0,2*18*sizeof(Bool_t));
  //
  TString path = Form(dataPath.Data(),run);
  //
  for (int id=0;id<2;id++) {delete regMap[id]; regMap[id] = 0;}
  for (int id=0;id<2;id++) {
    TString fileName = Form(Form("%%s/his%%d/%s.root",fileTemplate.Data()),path.Data(),useHistID[id],dimS[id].Data(),det,"",sec0,sec1,th0,th1);
    if (gSystem->AccessPathName(fileName.Data())) {
      printf("File %s not found\n",fileName.Data());
      return kFALSE;
    }
    TString objName = Form(mapTemplate.Data(),dimS[id].Data(),det,fitVersion.Data(),sec0,sec1,th0,th1);
    //printf("Need to open new file %s\n",fileName.Data());
    TFile* flinp = TFile::Open(fileName.Data());
    regMap[id] = (AliNDLocalRegression*)flinp->Get(objName.Data());
    if (regMap[id]) {
      if (!firstRep[sect+(zPos?0:18)]) printf("Loaded %s\n",objName.Data());
    }
    else {
      printf("Failed to load %s\n",objName.Data());
      return kFALSE;
    }
    //
    flinp->Close();
    delete flinp;
  }
  // registor loaded sectors
  for (int is=sec0;is<sec1;is++) {
    statusMap[is+(zPos?0:18)] = kTRUE;
    firstRep[is+(zPos?0:18)] = kTRUE;
  }
  //
  // define boundaries
  const THn* hn = regMap[0]->GetHistogram();
  for (int i=0;i<kNDimReg;i++) {
    TAxis* ax = hn->GetAxis(i);
    bounds[i][0] = ax->GetXmin();
    bounds[i][1] = ax->GetXmax();
  }
  //
  return kTRUE;
}

void SetRun(int run)
{
  curRun = run;
  SetupTPCCalibDB(run);
}


//__________________________________________________________________________
Bool_t GetDistortion(int sector, int row, float y2x, float z2x, float xyz[3])
{
  // computes NDLocal parameterization, i.e. gets distortion at measured cluster position
  // e.g. true = cluster + xyz
  // The input is in sector coordinates, sector in 0:17, row in 0:158 format

  AliTPCcalibDB* calib = AliTPCcalibDB::Instance();
  AliTPCParam* tpcParam = calib->GetParameters();
  if (!LoadDistMap(curRun,sector,z2x>0,curDetComb)) return kFALSE;
  double vecInp[kNDimReg]={0};
  int roc72 = sector;
  int row72 = row;
  if (z2x<0) roc72 += 18;
  if (row72>=63) {
    roc72 += 36;
    row72 -= 63;
  }
  double xRow = tpcParam->GetPadRowRadii(roc72,row72);
  double xyz0[3] = {xRow, y2x*xRow, z2x*xRow};

  vecInp[kndR] = xRow;
  //
  // convert y to "sector" coordinate
  vecInp[kndSec]  = sector + TMath::ATan2(xyz0[1],xyz0[0])/kSectWindow + 0.5;
  vecInp[kndZ2R]  = z2x;

  if (hackZ2RtoZ2X) { // in the old maps Z/R is used
    vecInp[kndZ2R] = xyz0[2]/TMath::Sqrt(xyz0[0]*xyz0[0]+xyz0[1]*xyz0[1]);
  }
  //
  double dRPhi,dZ;
  Bool_t xok = kFALSE;
  while (corrX && TMath::Abs(bZ)>1e-4) {
    // assume deltaY = DeltaY - DeltaX * tg(phi)
    // where deltaY is measured track-cluster residual at the pad-row
    // DeltaY and DeltaX are real distortions in Y and X: track (ionization point) - cluster

    double deltSum=0,tgsum=0,tg2sum=0,deltgsum=0, dzsum=0;
    for (int iqpt=kNQ2PtQuery;iqpt--;) {
      vecInp[kndQ2PT] = kQ2PtQuery[iqpt];
      EvalNDLocal(vecInp,dRPhi,dZ);
      double tg = GetTgPhi(vecInp[kndR],y2x, vecInp[kndQ2PT],bZ);
      deltSum += dRPhi;
      dzsum   += dZ;
      tgsum   += tg;
      tg2sum  += tg*tg;
      deltgsum += dRPhi*tg;
    }
    deltSum /= kNQ2PtQuery;
    tgAngCent = tgsum /= kNQ2PtQuery;
    tg2sum /= kNQ2PtQuery;
    deltgsum /= kNQ2PtQuery;
    dzsum /= kNQ2PtQuery;
    double det = tgsum*tgsum - tg2sum;
    if (TMath::Abs(det)<1e-6) {
      printf("q/pt evaluation at Sector: %.1f X: %.1f Z/X: %+.1f gives %e determinant\n",
	     vecInp[kndSec],vecInp[kndR],vecInp[kndZ2R],det);
      printf("sums are: sDY: %+e sTg: %+e sTg2: %+e sTgDY: %+e sDZ: %+e\n\n", deltSum,tgsum,tg2sum,deltgsum, dzsum);
      break;
    }
    xyz[0] = (deltgsum - deltSum*tgsum)/det;
    xyz[1] = (deltgsum*tgsum - deltSum*tg2sum)/det;
    xyz[2] = (dzsum + xyz[0]*z2x);
    xok = kTRUE;
    break;
  }
  //
  if (!xok) {
    vecInp[kndQ2PT] = kq2ptCent;
    EvalNDLocal(vecInp,dRPhi,dZ);
    //
    tgAngCent = GetTgPhi(vecInp[kndR],y2x, kq2ptCent ,bZ);
    //
    xyz[0] = 0;
    xyz[1] = dRPhi;
    xyz[2] = dZ;
  }
  //
  if (useIniTransform) {
    double xyz0Tr[3] = {xyz0[0],xyz0[1],xyz0[2]};
    //xyz[0]=xyz[1]=xyz[2]=0;
    Transform(roc72,row72, xyz0Tr);
    for (int i=3;i--;) xyz[i] += xyz0Tr[i]-xyz0[i]; // default correction, if distortions are not done wrt raw data  
  }
  //printf("--> %+8.4f %+8.4f %+8.4f\n",xyz[0],xyz[1],xyz[2]);
  return kTRUE;
}

/*

//__________________________________________________________________________
Bool_t GetDistortion(int sector, int row, float y2x, float z2x, float xyz[3])
{
  // computes NDLocal parameterization, i.e. gets distortion at measured cluster position
  // e.g. true = cluster + xyz
  // The input is in sector coordinates, sector in 0:17, row in 0:158 format

  AliTPCcalibDB* calib = AliTPCcalibDB::Instance();
  AliTPCParam* tpcParam = calib->GetParameters();
  if (!LoadDistMap(curRun,sector,z2x>0,curDetComb)) return kFALSE;
  double vecInp[kNDimReg]={0};
  int roc72 = sector;
  int row72 = row;
  if (z2x<0) roc72 += 18;
  if (row72>=63) {
    roc72 += 36;
    row72 -= 63;
  }
  double xRow = tpcParam->GetPadRowRadii(roc72,row72);
  double xyz0[3] = {xRow, y2x*xRow, z2x*xRow};

  vecInp[kndR] = xRow;
  //
  // convert y to "sector" coordinate
  vecInp[kndSec]  = sector + TMath::ATan2(xyz0[1],xyz0[0])/kSectWindow + 0.5;
  vecInp[kndZ2R]  = z2x;

  if (hackZ2RtoZ2X) { // in the old maps Z/R is used
    vecInp[kndZ2R] = xyz0[2]/TMath::Sqrt(xyz0[0]*xyz0[0]+xyz0[1]*xyz0[1]);
  }
  //
  vecInp[kndQ2PT] = kq2ptCent;
  double dRPhi0,dZ0;
  EvalNDLocal(vecInp,dRPhi0,dZ0);
  //
  double tg0 = GetTgPhi(vecInp[kndR],y2x, kq2ptCent ,bZ);
  // tmp
  tgAngCent = tg0;
  //
  xyz[0] = 0;
  xyz[1] = dRPhi0;
  xyz[2] = dZ0;
  if (corrX && TMath::Abs(bZ)>1e-4) {
    // assume deltaY = DeltaY - DeltaX * tg(phi)
    // where deltaY is measured track-cluster residual at the pad-row
    // DeltaY and DeltaX are real distortions in Y and X: track (ionization point) - cluster

    vecInp[kndQ2PT] = kq2ptBin;
    double dRPhiP,dZP,dRPhiM,dZM;
    EvalNDLocal(vecInp,dRPhiP,dZP);
    vecInp[kndQ2PT] = -kq2ptBin;
    EvalNDLocal(vecInp,dRPhiM,dZM);
    //
    //    double tg0 = GetTgPhi(vecInp[kndR],y2x, 0,bZ);
    double tgP = GetTgPhi(vecInp[kndR],y2x, kq2ptBin,bZ);
    double tgM = GetTgPhi(vecInp[kndR],y2x,-kq2ptBin,bZ);
    //
    double dtgP = tgP - tg0, dtgM = tgM - tg0, dx=0,dy=0;
    int nmeas = 0;
    if (TMath::Abs(dtgP)>1e-4) {
      dx+= (dRPhi0-dRPhiP)/dtgP;          // dx from q2pt=0 and q2pt>0 pair
      dy+= (dRPhi0*tgP-dRPhiP*tg0)/dtgP;  // dy from q2pt=0 and q2pt>0 pair
      nmeas++;
    }
    //
    if (TMath::Abs(dtgM)>1e-4) {
      dx+= (dRPhi0-dRPhiM)/dtgM;          // dx from q2pt=0 and q2pt<0 pair
      dy+= (dRPhi0*tgM-dRPhiM*tg0)/dtgM;  // dy from q2pt=0 and q2pt<0 pair
      nmeas++;
    }
    //
    if (nmeas==2) {
      dx *= 0.5;
      dy *= 0.5;
    }

    ////    dx *=-1; // TMP Inversion is not needed
    if (nmeas) {
      xyz[0] = dx;
      xyz[1] = dy;
      xyz[2] = dZ0 + dx*z2x;
    }
  }
  //
  if (useIniTransform) {
    double xyz0Tr[3] = {xyz0[0],xyz0[1],xyz0[2]};
    //xyz[0]=xyz[1]=xyz[2]=0;
    Transform(roc72,row72, xyz0Tr);
    for (int i=3;i--;) xyz[i] += xyz0Tr[i]-xyz0[i]; // default correction, if distortions are not done wrt raw data  
  }
  //printf("--> %+8.4f %+8.4f %+8.4f\n",xyz[0],xyz[1],xyz[2]);
  return kTRUE;
}

 */

//__________________________________________________________________________
Bool_t GetDistortionInv(int sector, float xSect, float ySect, float zSect, float xyz[3])
{
  // computes inverse of NDLocal parameterization, i.e. gets distortion for true position of the point
  // e.g. cluster = true + xyz
  // The input is in sector coordinates, sector in 0:17, row in 0:158 format
  //
  return kTRUE;
}

//_______________________________________________
double GetTgPhi(double x, double y2x, double q2p, double b)
{
  // calculate tangent of primary track at sector frame at given x,y
  double y = y2x*x;
  double c = q2p*b*(-0.299792458e-3);
  if (TMath::Abs(c)<1e-9) return y2x;
  double r2 = y*y+x*x;
  double det = 4./r2 - c*c;
  if (det<0) printf("track of q2p=%f cannot reach x:%f y:%f\n",q2p,x,y);
  double snp = 0.5*(y*TMath::Sqrt(det)-c*x); // snp at vertex
  snp += x*c;  // snp at x,y
  return snp/TMath::Sqrt((1-snp)*(1+snp));
}

//_______________________________________________
void SetupTPCCalibDB(int run)
{
  AliCDBManager* man = AliCDBManager::Instance();
  //
  man->SetDefaultStorage("local:///cvmfs/alice.cern.ch/calibration/data/2015/OCDB"); //local://$ALICE_ROOT/../src/OCDB");
  if (gSystem->AccessPathName("OCDB.root", kFileExists)==0) {        
    man->SetSnapshotMode("OCDB.root");
  }
  if (gSystem->AccessPathName("localOCDBaccessConfig.C", kFileExists)==0) {        
    gROOT->ProcessLine(".x localOCDBaccessConfig.C");
  }
  //
  //  man->SetRaw(1);
  man->SetRun(run);
  AliGRPManager gm;
  gm.ReadGRPEntry();
  gm.SetMagField();
  const AliGRPObject* grp = gm.GetGRPData();
  timeMin = grp->GetTimeStart();
  timeMax = grp->GetTimeEnd();
  //
  AliTPCcalibDB* calib = AliTPCcalibDB::Instance();
  //
  TObjArray* rpa = 0;
  AliCDBEntry* ent = man->Get("TPC/Calib/RecoParam");
  rpa = (TObjArray*)ent->GetObject();
  fCurrentRecoParam = (AliTPCRecoParam*)(rpa->At(kRecoParamID))->Clone();
  fCurrentRecoParam->SetUseCorrectionMap(kFALSE);
  fCurrentRecoParam->SetUseComposedCorrection(kTRUE);
  fCurrentRecoParam->SetUseAlignmentTime(kTRUE); // !!! not for future object
  fCurrentRecoParam->SetUseTOFCorrection(kTRUE);
  //
  AliTPCTransform* trans = calib->GetTransform();
  trans->SetCurrentRecoParam(fCurrentRecoParam);
  trans->SetCurrentTimeStamp(fCurrentTimeStamp);
  trans->UpdateTimeDependentCache();
  //
  AliMagF* magF= (AliMagF*)TGeoGlobalMagField::Instance()->GetField();
  bZ = magF->SolenoidField(); //field in kGaus
  //
}


//_______________________________________________
void Transform(int sector,int row, double x[3])
{
  // emulate transform behaviour for what we need, starting from
  // input in global frame, the output is in sector frame
  //
  // sector, row is in 0-71 ROC conventions
  //
  //
  Bool_t isInRotated = kTRUE;
  //
  AliTPCcalibDB* calib = AliTPCcalibDB::Instance();
  AliTPCTransform* trans = calib->GetTransform();
  AliTPCParam  * param    = calib->GetParameters();
  //
  AliMagF* magF= (AliMagF*)TGeoGlobalMagField::Instance()->GetField();
  Double_t bzField = magF->SolenoidField(); //field in kGaus
  //  
  // //
  // if (fCurrentRecoParam->GetUseCorrectionMap()) {
  //   const AliTPCChebCorr* fCorrMapCache0 = trans->GetCorrMapCache0();
  //   const AliTPCChebCorr* fCorrMapCache1 = trans->GetCorrMapCache1();
  //   //
  //   float delta0[3], y2x=x[1]/x[0], z2x = fCorrMapCache0->GetUseZ2R() ? x[2]/x[0] : x[2];
  //   fCorrMapCache0->Eval(row,sector,y2x,z2x,delta0);
  //   // 
  //   // for time dependent correction need to evaluate 2 maps, assuming linear dependence
  //   if (fCorrMapCache1) {
  //     float delta1[3];
  //     fCorrMapCache1->Eval(row,sector,y2x,z2x,delta1);   
  //     UInt_t t0 = fCorrMapCache0->GetTimeStampCenter();
  //     UInt_t t1 = fCorrMapCache1->GetTimeStampCenter();
  //     // possible division by 0 is checked at upload of maps
  //     double dtScale = (t1-fCurrentTimeStamp)/double(t1-t0);
  //     for (int i=3;i--;) delta0[i] += (delta1[i]-delta0[i])*dtScale;
  //   }
  //   for (int i=3;i--;) x[i] += delta0[i];
  // }
  //
  if(fCurrentRecoParam->GetUseExBCorrection()) {
    if (!isInRotated) {
      trans->RotatedGlobal2Global(sector,x);
      isInRotated = kFALSE;    
    }
    Double_t xx[3];
    calib->GetExB()->Correct(x,xx);   // old ExB correction
    for (int i=3;i--;) x[i] = xx[i];
  }
  //
   //
  if(fCurrentRecoParam->GetUseComposedCorrection()) {
    //
    if (isInRotated) {
      trans->RotatedGlobal2Global(sector,x);
      isInRotated = kFALSE;          
    }
    // HACK: composed correction has jumps exactly on the boundary
    // recalculate sector
    double hackDX=0,hackDY=0,hackDZ=0;
    //
    /*
    double phi = TMath::ATan2(x[1],x[0]);
    if (phi<0) phi+=TMath::Pi()*2;
    int sct = phi/(TMath::Pi()/9);
    double phiSect = (sector*20+10)*TMath::DegToRad();
    double hspan = 10*TMath::DegToRad();

    const double kTolAng = 2e-4, kTolZ = 2e-4;
    double rad = TMath::Sqrt(x[0]*x[0]+x[1]*x[1]);
    if (phi>phiSect) { // is it on the upper boundary
      if (TMath::Abs(phiSect+hspan - phi)<kTolAng) { // move point inside the sector
	hackDX = TMath::Cos(phi-kTolAng)*rad - x[0];
	hackDY = TMath::Sin(phi-kTolAng)*rad - x[1];
      }
    }
    else {
      if (TMath::Abs(phiSect-hspan - phi)<kTolAng) { // move point inside the sector
	hackDX = TMath::Cos(phi+kTolAng)*rad - x[0];
	hackDY = TMath::Sin(phi+kTolAng)*rad - x[1];
      }
    }
    // check if side is consistent with Z
    int side = (sector/18)&0x1;
    if (TMath::Abs(x[2])<kTolZ) {
      if (side==0) hackDZ = kTolZ;
      else         hackDZ = -kTolZ;
    }
    //
    x[0] += hackDX;
    x[1] += hackDY;
    x[2] += hackDZ;
    //
    */
    AliTPCCorrection * correction = calib->GetTPCComposedCorrection();   // first user defined correction  // if does not exist  try to get it from calibDB array
    if (!correction) correction = calib->GetTPCComposedCorrection(bzField);
    AliTPCCorrection * correctionDelta = calib->GetTPCComposedCorrectionDelta();
    if (correction) {
      Float_t distPoint[3]={static_cast<Float_t>(x[0]),static_cast<Float_t>(x[1]),static_cast<Float_t>(x[2])};
      correction->CorrectPoint(distPoint, sector);
      for (int i=3;i--;) x[i]=distPoint[i];
      if (correctionDelta&&fCurrentRecoParam->GetUseAlignmentTime()){  // appply time dependent correction if available and enabled
	Float_t distPointDelta[3]={static_cast<Float_t>(x[0]),static_cast<Float_t>(x[1]),static_cast<Float_t>(x[2])};
	correctionDelta->CorrectPoint(distPointDelta, sector);
	for (int i=3;i--;) x[i]=distPointDelta[i];
      }
    }
    // undo the hack
    x[0] -= hackDX;
    x[1] -= hackDY;
    x[2] -= hackDZ;
  }
  //
  // Time of flight correction
  //
  if (fCurrentRecoParam->GetUseTOFCorrection()){
    const Int_t kNIS=param->GetNInnerSector(), kNOS=param->GetNOuterSector();
    Float_t sign=1;
    if (sector < kNIS) {
      sign = (sector < kNIS/2) ? 1 : -1;
    } else {
      sign = ((sector-kNIS) < kNOS/2) ? 1 : -1;
    }
    Float_t deltaDr =0;
    Float_t dist=0;
    dist+=x[0]*x[0];  // (fPrimVtx[0]-x[0])*(fPrimVtx[0]-x[0]); // RS the vertex is anyway close to 0
    dist+=x[1]*x[1];  // (fPrimVtx[1]-x[1])*(fPrimVtx[1]-x[1]);
    dist+=x[2]*x[2]; //(fPrimVtx[2]-x[2])*(fPrimVtx[2]-x[2]);
    dist = TMath::Sqrt(dist);
    // drift length correction because of TOF
    // the drift velocity is in cm/s therefore multiplication by 0.01
    deltaDr = (dist*(0.01*param->GetDriftV()))/TMath::C();
    x[2]+=sign*deltaDr;
  }
  //
  //
  if (!isInRotated) {
    trans->Global2RotatedGlobal(sector,x);
    isInRotated = kTRUE;
  }
  //
  /*
  // Apply non linear distortion correction
  //
  int useFieldCorr = fCurrentRecoParam->GetUseFieldCorrection();
  if (useFieldCorr&(0x2|0x4|0x8)) {
    AliTPCCalPad * distortionMapY=0,*distortionMapZ=0,*distortionMapR=0;
    //
    // wt - to get it form the OCDB
    // ignore T1 and T2
    Double_t vdrift = param->GetDriftV()/1000000.; // [cm/us]   // From dataBase: to be updated: per second (ideally)
    Double_t ezField = 400; // [V/cm]   // to be updated: never (hopefully)
    if (sector%36<18) ezField*=-1;
    //    Double_t wt = -10.0 * (bzField*10) * vdrift / ezField ;
    Double_t wt = -10.0 * (bzField) * vdrift / ezField ; // RS: field is already in kGauss
    Double_t c0=1./(1.+wt*wt);
    Double_t c1=wt/c0;
    
    //can be switch on for each dimension separatelly
    if (useFieldCorr&0x2 && (distortionMapY=calib->GetDistortionMap(0))) {
      x[1]-= c0*distortionMapY->GetCalROC(sector)->GetValue(row,pad);
      x[0]-= c1*distortionMapY->GetCalROC(sector)->GetValue(row,pad);
    }
    if (useFieldCorr&0x4 && (distortionMapZ=calib->GetDistortionMap(1))) {
      x[2]-=distortionMapZ->GetCalROC(sector)->GetValue(row,pad);
    }
    if (useFieldCorr&0x8 && (distortionMapR=calib->GetDistortionMap(2))) {
      x[0]-= c0*distortionMapR->GetCalROC(sector)->GetValue(row,pad);
      x[1]-=-c1*distortionMapR->GetCalROC(sector)->GetValue(row,pad)*wt;
    }
    //
  }
  */
  //
  /*
  if (isInRotated) {
    trans->RotatedGlobal2Global(sector,x);
    isInRotated = kFALSE;
  }
  */
}


void EvalNDLocal(double *vec,double &drp,double &dz)
{
  drp = regMap[kndResRPhi]->Eval(vec);
  dz  = regMap[kndResZ2R]->Eval(vec);
  if (usedZTrack) { // the resudals were defined at track position, we need at cluster position
    double zc = vec[kndZ2R]*vec[kndR];
    double zt = zc + dz; // 1st guess on true zt at measured zc
    // use Newton-Raphson method for NDLocal inversion to get zt = F(zc) from 
    // dz == zt - zc = NDLoc(zt)   ->  zc - [zt - NDLoc(zt)]=0 
    // ->  zt_{i+1} = zt_i - (zc - [zt - NDLoct(zt_i)])/(-d[zt_i-NDLoc(zt_i)]/dz)
    double deriv[kNDimReg], dztmp; 
    int it=0;
    double sav = vec[kndZ2R], change = 0, rInv = 1./vec[kndR];
    do {
      vec[kndZ2R] = zt*rInv;
      regMap[kndResZ2R]->EvalAndDerivative(vec,dztmp,deriv);
      double bot = 1. - deriv[kndZ2R]*rInv;  // dF(zt_i)/dz
      if (TMath::Abs(bot)<1e-6) break;
      double top = zc - (zt - dztmp); // zc - F(zt_i) 
      double change = top/bot;
      //  printf("It %d Eps:%+e, Zc:%+e, Zt:%+e Ztn:%+e DZ:%+e | dztmp: :%+e DD:%+e\n",
      //     it,change,zc,zt,zt+change,dz,dztmp,deriv[kndZ2R]*rInv);
      zt += change;
      it++;
    } while(it<inversionMaxIt && TMath::Abs(change)>inversionEps);
    vec[kndZ2R] = sav;
    dz = zt - zc;
    //    printf("--> %+e\n",dz);
  }
  //
}
