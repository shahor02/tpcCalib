#ifndef ALITPCDISTORTIONEXTRACTOR_H
#define ALITPCDISTORTIONEXTRACTOR_H
#include <TSystem.h>
#include <TTree.h>
#include <TBranch.h>
#include <TFile.h>
#include <TVectorF.h>
#include <TVectorD.h>
#include <TString.h>
#include <TMath.h>
#include <TGeoGlobalMagField.h>
#include <TGrid.h>
#include <TNDArray.h>
#include <THn.h>
#include <TH1F.h>
#include <TF1.h>
#include <TEnv.h>
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

class AliTPCDistortionExtraction: public TObject
{
 public:
  enum {kEpanechnikovKernel, kGausianKernel};  // defined kernels
  enum {kNSect=18,kNSect2=2*kNSect,kNROC=4*kNSect,kNPadRows=159, kNRowIROC=63, kNRowOROC1=64, kNRowOROC2=32};

  //
  // the voxels are defined in following space
  enum {kVoxQ,   // tg of track inclination wrt pad row, each voxel has its own range according to F,X,Y
	kVoxF,   // y/x in sector coordinates
	kVoxX,   // sector X coordinate
	kVoxZ,   // Z/X sector coordinates
	kVoxV,   // variable within the voxel (delta, stat, etc): last dimension of all THn histos
	kVoxHDim, kVoxDim=kVoxHDim-1};

  enum {kResX,kResY,kResZ,kResDim}; // output dimensions

  // content of processed voxel
  enum {kEstNorm,kEstMean,kEstSig,kEstMax,  // statistics
	kEstMeanL,kEstMeanEL, kEstSigL,kEstSigEL, // LTM
	kEstNormG,kEstMeanG,kEstMeanEG,kEstSigG,kEstSigEG,kEstChi2G, // gaussian fit
	kNEstPar};

  struct dts_t {  // struct for basic residual
    UChar_t dy;   // Y residual
    UChar_t dz;   // Z residual
    UChar_t bvox[kVoxDim]; // voxel bin info: kVoxQ,kVoxF,kVoxX,kVoxZ
  };

  struct bres_t  { // final result for the voxel
    Float_t D[kResDim];      // values of extracted distortions
    Float_t E[kResDim];      // their errors
    Float_t stat[kVoxHDim];  // statistics: averages (weigted over Q bins) of each voxel dimension + entries
    Float_t DS[kResDim];     // smoothed value
    UChar_t bvox[kVoxDim];   // voxel identifier, here the bvox[0] shows number of Q bins used for Y
    UChar_t bsec;            // sector ID (0-35)
    UChar_t smooth;          // smoother flag
  };

  struct bstat_t {           // stat info on the voxel
    Float_t stat[kVoxHDim];  // statistics: averages of each voxel dimension + entries
    Float_t distY[kNEstPar]; // distortion estimators for Y
    Float_t distZ[kNEstPar]; // distortion estimators for Z
    UChar_t bvox[kVoxDim];   // voxel identifier
    UChar_t bsec;            // sector ID (0-35)
  };

  struct voxDef_t {          // to dumpe the voxel definition (within sector)
    UChar_t bvox[kVoxDim];   // voxel identifier
    Float_t vmin[kVoxDim];   // min boundary
    Float_t vmax[kVoxDim];   // max boundary
  };


public:

  AliTPCDistortionExtraction();
  virtual ~AliTPCDistortionExtraction();
  
  void Init(int run,const char * residualList,Long64_t tmin,Long64_t tmax,float maxDY,float maxDZ
	    ,float maxQ2Pt,int nY2XBins,int nZ2XBins,int nXBins,int nDeltaBinsY,int nDeltaBinsZ
	    ,int maxTracks,
	    Bool_t fixAlignmentBug,int cacheInp,int learnSize,Bool_t switchCache);
  void CollectData();
  void WriteVoxelDefinitions();
  void ProcessSectorResiduals(int is, bstat_t &voxStat);
  void ExtractVoxelData(bstat_t &stat,const TNDArrayT<short>* harrY,
			const TNDArrayT<short>* harrZ,const TNDArrayT<float>* harrStat);
  void ExtractDistortionsData(TH1F* histo, float est[kNEstPar],float minNorm, 
			      float fracLTM, float fitNSig);

  void InitForBugFix(const char* ocdb);
  THnF* CreateVoxelStatHisto(int sect);
  THn* CreateSectorResidualsHisto(int sect, int nbDelta,float range, const char* pref);

  void    LoadVDrift();
  Float_t GetDriftCorrection(float z, float x, float phi, int rocID);
  Float_t tgpXY(float x, float y, float q2p, float bz);
  
  void WriteStatHistos();
  void LoadStatHistos();
  void WriteResTree();

  Float_t ExtractResidualHisto(const TNDArrayT<short>* harr, const Long64_t bprod[kVoxHDim], 
			       const UChar_t vox[kVoxDim], TH1F* dest);

  Float_t ExtractResidualHisto(const TNDArrayT<short>* harr, const Long64_t bprod[kVoxHDim], 
			       const UChar_t voxMin[kVoxDim], const UChar_t voxMax[kVoxDim], TH1F* dest);

  TH1F* ExtractResidualHisto(Bool_t y, int sect, const UChar_t vox[kVoxDim]);
  TH1F* ExtractResidualHisto(Bool_t y, int sect, const UChar_t vox[kVoxDim], const UChar_t vox1[kVoxDim]);
  void  ExtractXYZDistortions();
  Bool_t ExtractVoxelXYZDistortions(const bstat_t voxIQ[kNQBins],bres_t &res, int minStat, 
				    float maxGChi2, int minYBinsOK);

  void FillHoles(int isect, bres_t *sectData, const int fNBProdSectG[2], int minGoodPoints);
  
  Bool_t FitPoly2(const float* x,const float* y, const float* w, int np, float *res, float *err);
  Bool_t FitPoly1(const float* x,const float* y, const float* w, int np, float *res, float *err);

  Int_t Smooth0(int isect);
  Bool_t GetSmoothEstimate(int isect, float x, float p, float z, float *res, float *deriv);
  void SetKernelType(int tp, float bwX, float bwP, float bwZ, float scX,float scP,float scZ);
  
  void CreateCorrectionObject();
  void InitBinning(int nbx, int nby, int nbz);
  Int_t GetXBin(float x);
  Int_t GetRowID(float x);
  
  Bool_t FindVoxelBin(int sectID,float q2pt,float x,float y,float z,UChar_t bin[kVoxHDim],float voxVars[kVoxHDim]);
  
  Int_t   GetXBinExact(float x);
  Float_t GetY2X(int ix, int iy);
  Float_t GetY2XLow(int ix, int iy);
  Float_t GetDY2X(int ix);
  Float_t GetDY2XI(int ix);
  Float_t GetX(int i);
  Float_t GetXLow(int i);
  Float_t GetDX(int i);
  Float_t GetDXI(int i);
  Int_t   GetY2XBinExact(float y2x, int ix);
  Int_t   GetY2XBin(float y2x, int ix);
  Int_t   GetZ2XBinExact(float z2x);
  Int_t   GetZ2XBin(float z2x);
  Float_t GetZ2XLow(int iz);
  Float_t GetDZ2X();
  Float_t GetDZ2XI();
  void    FindVoxel(float x, float y2x, float z2x, int &ix,int &ip, int &iz);
  void    FindVoxel(float x, float y2x, float z2x, UChar_t &ix,UChar_t &ip, UChar_t &iz);
  void    GetVoxelCoordinates(int isec, int ix, int ip, int iz,float &x, float &p, float &z);
  Double_t GetKernelWeight(double u2);
  //  Int_t  GetQBin(float tgp, int binX, int binY);
  Int_t  GetQBin(float tgp);
  Long64_t GetBin2Fill(const Long64_t bprod[kVoxHDim],const UChar_t binVox[kVoxDim], UShort_t bVal);
  Int_t  GetVoxGBin(int ix, int ip, int iz);
  Int_t  GetVoxGBin(UChar_t bvox[kVoxDim]);

protected:
  //
  Bool_t   fInitDone = kFALSE;
  Bool_t   fUseErrInSmoothing;                      // weight kernel by point error
  Bool_t   fSwitchCache;                            // reset the cache when the reading mode is changing
  Bool_t   fFixAlignmentBug;                        // flag to apply the fix
  Bool_t   fApplyZt2Zc;                             // Apply fix for using Z_track instead of Z_cluster in the data


  // --------------------------------Chebyshev object creation 
  Int_t    fChebZSlicePerSide1;                     // z partitions per side
  Int_t    fChebPhiSlicePerSector;                  // azimuthal partitions per sector
  Int_t    fNPCheb[3][2];                           // cheb. nodes per slice

  Float_t  fChebPrecD[3];                           // nominal precision per output dimension
  AliTPCChebCorr* fChebCorr;                        // final Chebyshev object

  // -------------------------------Task defintion
  Int_t    fRun;     // run numbet 
  Long64_t fTMin;    // time start
  Long64_t fTMax;    // time stop
  Int_t    fMaxTracks;  // max tracks to accept
  Int_t    fCacheInp;      // input trees cache in MB
  Int_t    fLearnSize;     // event to learn for the cache
  Float_t  fBz;            // B field
  TString  fResidualList;   // list of residuals tree

  // -------------------------------Binning
  Float_t  fMaxDY;   // max residual in Y
  Float_t  fMaxDZ;   // max residual in Z
  Float_t  fMaxQ2Pt; // max |q/pt|
  Float_t  fMidQ2Pt; // middle |q/pt| for slopes binning 
  Int_t    fNY2XBins=-1;  // y/x bins per sector
  Int_t    fNZ2XBins;    // z/x bins per sector
  Int_t    fNXBins=-1;    // n bins in radial dim.
  Int_t    fNDeltaYBins; // n bins in Y residual space
  Int_t    fNDeltaZBins; // n bins in Z residual space
  Bool_t   fUniformBins[kVoxDim]; // uniform binning? Currently only X may be non-uniform (per pad-row)


  Float_t  fDZ2X;            // Z2X bin size
  Float_t  fDX;            // X bin size
  Float_t  fDZ2XI;           // inverse Z2X bin size 
  Float_t  fDXI;           // inverse X bin size 
  Float_t  fDeltaYbinI;    // inverse deltaY bin size
  Float_t  fDeltaZbinI;    // inverse deltaZ bin size

  Int_t    fNGVoxPerSector; // total number of geometrical voxels per sector (excluding Q binning)

  Float_t  *fMaxY2X;        // max Y/X at each X bin, account for dead zones
  Float_t  *fDY2X;          // Y/X bin size at given X bin
  Float_t  *fDY2XI;         // inverse of Y/X bin size at given X bin
  Float_t  *fBinMinQ;       // min value of tg(inclination) at given X,Y bin
  Float_t  *fBinDQ;         // tg(inclination) bin size at given X,Y bin
  Float_t  *fBinDQI;        // inverse of tg(inclination) bin size at given X,Y bin

  Long64_t fNBProdSt[kVoxHDim]; // aux arrays for fast bin calculation
  Long64_t fNBProdDY[kVoxHDim];
  Long64_t fNBProdDZ[kVoxHDim];
  Int_t    fNBProdSectG[2];   // aux info for fast bin index calculation in geom voxel space


  // ------------------------------Smoothing
  Int_t    fNMaxNeighb;     // max neighbours to loop for smoothing
  Int_t    fKernelType;     // kernel type
  Int_t    fStepKern[kVoxDim] = {0,2,2,4};
  Float_t  fKernelScaleEdge[kVoxDim] = {1, 1.,1.,1.}; // scaling factor for edge points

  Double_t fLastSmoothingRes[kResDim*4];  // result of last kernel minimization

  // ------------------------------Selection Stats
  Int_t    fNTrSelTot;      // selected tracks
  Int_t    fNTrSelTotWO;    // would be selected w/o outliers rejection
  Int_t    fNReadCallTot;   // read calls from input trees
  Int_t    fNBytesReadTot;  // total bytes read


  // ------------------------------VDrift correction
  TVectorD     *fVDriftParam;
  TGraphErrors *fVDriftGraph;  
  Float_t      fCorrTime;

  // -----------------------------Results of processing
  bres_t *fSectGVoxRes[kNSect2];         //[fNGVoxPerSector] sectors results for geometric voxel
  TTree* fStatTree;                      //! tree with voxels statistics
  TTree* fTmpTree[kNSect2];              //! IO tree per sector
  TFile* fTmpFile[kNSect2];              //! file for fTmpTree
  THnF*  fStatHist[kNSect2];             //! histos for statistics bins
  TNDArrayT<float> *farrNDstat[kNSect2]; //! alias arrays for fast access to fStatHist

  TH1F* fHDelY;                          //! work histo for delta Y fits
  TH1F* fHDelZ;                          //! work histo for delta Z fits

  // ---------------------------------- outliers rejection
  Bool_t  fFilterOutliers;               // reject outliers
  Int_t   fMaxSkippedCluster=10;  // 10 cluster
  Float_t fMaxRMSTrackCut=2.0;    // maximal RMS (cm) between the tracks 
  Float_t fMaxRMSClusterCut=0.3;    // maximal RMS (cm) between the cluster and local mean
  Float_t fMaxDeltaClusterCut=0.5;    // maximal delta(cm) between the cluster and local mean
  
  //
  static const float kSecDPhi;
  static const float kSecDPhiH;
  static const float kMaxY2X; // max Y/X in sector coordinates (w/o excluding dead zone)
  static const float kMinX;   // min X to cover
  static const float kMaxX;   // max X to cover
  static const float kMaxZ2X;   // max z/x
  static const float kZLim;   // endcap position
  static const char* kTmpFileName;
  static const char* kStatOut;
  static const char* kResOut;
  static const char* kDriftFileName;
  static const float kDeadZone;  // dead zone on sector edges in cm
  static const int   kNQBins;     // number of bins in voxQ variable
  static const ULong64_t kMByte;
  static const Float_t kZeroK; // zero kernel weight

  static const char* kVoxName[];
  static const char* kResName[];
  static const char* kEstName[];
  
  static const Float_t kTPCRowX[]; // X of the pad-row
  static const Float_t kTPCRowDX[]; // pitch in X

  ClassDef(AliTPCDistortionExtraction,1);
};

//________________________________________________________________
inline Int_t AliTPCDistortionExtraction::GetXBinExact(float x) 
{
  // convert X to bin ID, following pad row widths
  if (fUniformBins[kVoxX]) {
    int ix = (x-kMinX)*fDXI;
    return (ix<0 || ix>=fNXBins) ? -2 : ix;
  }
  else return GetRowID(x);
}

//________________________________________________________________
inline Float_t AliTPCDistortionExtraction::GetY2X(int ix, int iy)
{
  // get Y2X bin center for ix,iy bin
  return (0.5f+iy)*fDY2X[ix] - fMaxY2X[ix];
}

//________________________________________________________________
inline Float_t AliTPCDistortionExtraction::GetY2XLow(int ix, int iy)
{
  // get Y2X bin low edge for ix,iy bin
  return iy*fDY2X[ix] - fMaxY2X[ix];
}

//________________________________________________________________
inline Float_t AliTPCDistortionExtraction::GetDY2X(int ix)
{
  // get Y2X bin size value for ix bin
  return fDY2X[ix];
}

//________________________________________________________________
inline Float_t AliTPCDistortionExtraction::GetDY2XI(int ix)
{
  // get Y2X inverse bin size  for ix bin
  return fDY2XI[ix];
}

//________________________________________________________________
inline Float_t AliTPCDistortionExtraction::GetX(int i)
{
  // low edge of i-th X bin
  return (fUniformBins[kVoxX]) ? kMinX+(0.5+i)*fDX : kTPCRowX[i];
}

//________________________________________________________________
inline Float_t AliTPCDistortionExtraction::GetXLow(int i)
{
  // low edge of i-th X bin
  return fUniformBins[kVoxX] ? kMinX+i*fDX : kTPCRowX[i] - 0.5*kTPCRowDX[i];
}

//________________________________________________________________
inline Float_t AliTPCDistortionExtraction::GetDX(int i)
{
  // width of i-th X bin
  return fUniformBins[kVoxX] ? fDX : kTPCRowDX[i];
}

//________________________________________________________________
inline Float_t AliTPCDistortionExtraction::GetDXI(int i)
{
  // inverse width of i-th X bin
  return (fUniformBins[kVoxX]) ? fDXI : 1.f/kTPCRowDX[i];
}

//________________________________________________________________
inline Int_t AliTPCDistortionExtraction::GetY2XBinExact(float y2x, int ix) 
{
  // get exact y2x bin at given x range
  float bf = ( y2x + fMaxY2X[ix] ) * GetDY2XI(ix);
  if (bf<0) return -1;
  else if (bf>=fNY2XBins) return fNY2XBins;
  return int(bf);
}

//________________________________________________________________
inline Int_t AliTPCDistortionExtraction::GetY2XBin(float y2x, int ix) 
{
  // get closest y2x bin at given x range
  int bf = ( y2x + fMaxY2X[ix] ) * GetDY2XI(ix);
  if (bf<0) bf = 0;
  else if (bf>=fNY2XBins) fNY2XBins-1;
  return bf;
}

//________________________________________________________________
inline Int_t AliTPCDistortionExtraction::GetZ2XBinExact(float z2x)
{
  // get exact z2x bin at given x range
  float bz = TMath::Abs(z2x)*GetDZ2XI();
  if (bz>=fNZ2XBins) return -1;
  return int(bz);
}

//________________________________________________________________
inline Int_t AliTPCDistortionExtraction::GetZ2XBin(float z2x) 
{
  // get closest z2x bin
  int bz = TMath::Abs(z2x)*GetDZ2XI();
  return bz<fNZ2XBins ? bz : fNZ2XBins-1;
}

//________________________________________________________________
inline Float_t GetZ2X(int iz)
{
  // get Z2X bin center for iz, !! always positive
  return (0.5f+iz)*GetDZ2X();
}

//________________________________________________________________
inline Float_t AliTPCDistortionExtraction::GetZ2XLow(int iz)
{
  // get Z2X bin low edge for iz !! bin positive
  return iz*GetDZ2X();
}

//________________________________________________________________
inline Float_t AliTPCDistortionExtraction::GetDZ2X()
{
  // get Z2X bin size value
  return fDZ2X;
}

//________________________________________________________________
inline Float_t AliTPCDistortionExtraction::GetDZ2XI()
{
  // get Z2X inverse bin size
  return fDZ2XI;
}

//_____________________________________
inline void AliTPCDistortionExtraction::FindVoxel(float x, float y2x, float z2x, int &ix, int &ip, int &iz)
{
  // calculate voxel center sector coordinates (wrt sector)
  ix = GetXBin(x);
  ip = GetY2XBin(y2x,ix);
  iz = GetZ2XBin(z2x);
  //
}

//_____________________________________
inline void AliTPCDistortionExtraction::FindVoxel(float x, float y2x, float z2x, UChar_t &ix, UChar_t &ip, UChar_t &iz)
{
  // calculate voxel center sector coordinates (wrt sector)
  ix = GetXBin(x);
  ip = GetY2XBin(y2x,ix);
  iz = GetZ2XBin(z2x);
  //
}

//_____________________________________
inline void AliTPCDistortionExtraction::GetVoxelCoordinates(int isec, int ix, int ip, int iz, float &x, float &p, float &z)
{
  // calculate voxel center sector coordinates (wrt sector)
  x = GetX(ix);
  p = GetY2X(ix,ip);
  z = GetZ2X(iz);
  if (isec>=kNSect) z = -z;
}

//_____________________________________
inline Double_t AliTPCDistortionExtraction::GetKernelWeight(double u2)
{
  if (fKernelType == kEpanechnikovKernel) {
    if (u2>1) return 0.;
    return 3./4.*(1.-u2);
  }
  else if (fKernelType == kGaussianKernel) {
    return u2<5 ? TMath::Exp(-u2)/TMath::Sqrt(2.*TMath::Pi()) : 0;
  }
  else {
    ::Fatal("GetKernelWeight","Kernel type %d is not defined",fKernelType);
  }
}

//_____________________________________________________
inline Int_t AliTPCDistortionExtraction::GetQBin(float q2pt)
{
  // get binID in track Q variable (tg of inclination) for given X,Y bin
  //
  float q2ptA = TMath::Abs(q2pt);
  if (q2pTA>=fMaxQ2Pt) return -1;
  int bin = (q2pTA<fMidQ2Pt) ?  0 : 1;
  if (q2pt>0) bin = kNQBins-1 - bin;
}
//_____________________________________________________

/*
inline Int_t AliTPCDistortionExtraction::GetQBin(float tgp, int binX, int binY)
{
  // get binID in track Q variable (tg of inclination) for given X,Y bin
  //
  int id = binX*fNY2XBins+binY;
  float bf = (tgp - fBinMinQ[id])*fBinDQI[id];
  if (bf<0 || bf>kNQBins) return -1;
  return int(bf);
}
*/

//_____________________________________________________
inline Long64_t AliTPCDistortionExtraction::GetBin2Fill(const Long64_t bprod[kVoxHDim],
							const UChar_t binVox[kVoxDim], UShort_t bVal) 
{
  // TH5 bin calculation, bval is the last dimention binID
  ULong64_t binToFill = bVal+1; // 0 bin is undeflow
  for (int id=kVoxDim;id--;) binToFill += bprod[id]*(1+binVox[id]);
  return binToFill;
}

//_____________________________________________________
inline Int_t AliTPCDistortionExtraction::GetVoxGBin(int ix, int ip, int iz) 
{
  // index of geometrix voxel (no Q info)
  return iz+fNBProdSectG[1]*ip+fNBProdSectG[0]*ix;
}

//_____________________________________________________
inline Int_t AliTPCDistortionExtraction::GetVoxGBin(UChar_t bvox[kVoxDim]) 
{
  // index of geometrix voxel (no Q info)
  return bvox[kVoxZ]+fNBProdSectG[1]*bvox[kVoxF]+fNBProdSectG[0]*bvox[kVoxX];
}


#endif
