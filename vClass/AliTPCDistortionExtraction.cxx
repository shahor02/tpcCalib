#include "AliTPCDistortionExtraction.h"


const char* AliTPCDistortionExtraction::kVoxName[kVoxDim] = {"tgSlp","y2x","x","z2x"};
const char* AliTPCDistortionExtraction::kResName[kResDim] = {"dX","dY","dZ"};
const char* AliTPCDistortionExtraction::kEstName[kNEstPar] = {
  "Nrm","Mean","Sig","Max",
  "MeanL","MeanEL","SigL","SigEL",
  "NormG","MeanG","MeanEG","SigG","SigEG","Chi2G"};

const float AliTPCDistortionExtraction::kSecDPhi=20.f*TMath::DegToRad();
const float AliTPCDistortionExtraction::kSecDPhiH = AliTPCDistortionExtraction::kSecDPhi*0.5f;
const float AliTPCDistortionExtraction::kMaxY2X = TMath::Tan(kSecDPhiH);
const float AliTPCDistortionExtraction::kMinX = 85.0f;
const float AliTPCDistortionExtraction::kMaxX = 246.0f;
const float AliTPCDistortionExtraction::kMaxZ2X = 1.0f;
const float AliTPCDistortionExtraction::kZLim = 250.0f;
const char* AliTPCDistortionExtraction::kTmpFileName  = "tmpDeltaSect";
const char* AliTPCDistortionExtraction::kStatOut      = "voxelStat";
const char* AliTPCDistortionExtraction::kResOut       = "voxelRes";
const char* AliTPCDistortionExtraction::kDriftFileName= "fitDrift";
const float AliTPCDistortionExtraction::kDeadZone = 1.5;
const int   AliTPCDistortionExtraction::kNQBins = 4;
const ULong64_t AliTPCDistortionExtraction::kMByte = 1024LL*1024LL;
const Float_t AliTPCDistortionExtraction::kZeroK = 1e-6;


const Float_t AliTPCDistortionExtraction::kTPCRowX[AliTPCDistortionExtraction::kNPadRows] = { // pad-row center X
  85.225, 85.975, 86.725, 87.475, 88.225, 88.975, 89.725, 90.475, 91.225, 91.975, 92.725, 93.475, 94.225, 94.975, 95.725,
  96.475, 97.225, 97.975, 98.725, 99.475,100.225,100.975,101.725,102.475,103.225,103.975,104.725,105.475,106.225,106.975,
  107.725,108.475,109.225,109.975,110.725,111.475,112.225,112.975,113.725,114.475,115.225,115.975,116.725,117.475,118.225,
  118.975,119.725,120.475,121.225,121.975,122.725,123.475,124.225,124.975,125.725,126.475,127.225,127.975,128.725,129.475,
  130.225,130.975,131.725,135.100,136.100,137.100,138.100,139.100,140.100,141.100,142.100,143.100,144.100,145.100,146.100,
  147.100,148.100,149.100,150.100,151.100,152.100,153.100,154.100,155.100,156.100,157.100,158.100,159.100,160.100,161.100,
  162.100,163.100,164.100,165.100,166.100,167.100,168.100,169.100,170.100,171.100,172.100,173.100,174.100,175.100,176.100,
  177.100,178.100,179.100,180.100,181.100,182.100,183.100,184.100,185.100,186.100,187.100,188.100,189.100,190.100,191.100,
  192.100,193.100,194.100,195.100,196.100,197.100,198.100,199.350,200.850,202.350,203.850,205.350,206.850,208.350,209.850,
  211.350,212.850,214.350,215.850,217.350,218.850,220.350,221.850,223.350,224.850,226.350,227.850,229.350,230.850,232.350,
  233.850,235.350,236.850,238.350,239.850,241.350,242.850,244.350,245.850
};
const Float_t AliTPCDistortionExtraction::kTPCRowDX[AliTPCDistortionExtraction::kNPadRows] = { // pad-row pitch in X
  0.750,0.750,0.750,0.750,0.750,0.750,0.750,0.750,0.750,0.750,0.750,0.750,0.750,0.750,0.750,0.750,0.750,0.750,0.750,0.750,0.750,
  0.750,0.750,0.750,0.750,0.750,0.750,0.750,0.750,0.750,0.750,0.750,0.750,0.750,0.750,0.750,0.750,0.750,0.750,0.750,0.750,0.750,
  0.750,0.750,0.750,0.750,0.750,0.750,0.750,0.750,0.750,0.750,0.750,0.750,0.750,0.750,0.750,0.750,0.750,0.750,0.750,0.750,0.750,
  1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,
  1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,
  1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,
  1.500,1.500,1.500,1.500,1.500,1.500,1.500,1.500,1.500,1.500,1.500,1.500,1.500,1.500,1.500,1.500,1.500,1.500,1.500,1.500,1.500,
  1.500,1.500,1.500,1.500,1.500,1.500,1.500,1.500,1.500,1.500,1.500
};


AliTPCDistortionExtraction::AliTPCDistortionExtraction( ) : 
  
  fInitDone(kFALSE)
  ,fUseErrInSmoothing(kTRUE)
  ,fSwitchCache(kFALSE)
  ,fFixAlignmentBug(kTRUE)
  ,fApplyZt2Zc(kTRUE)

  ,fChebZSlicePerSide(1)
  ,fChebPhiSlicePerSector(2)
  ,fChebCorr(0)

  ,fRun(0)
  ,fTMin(0)
  ,fTMax(0xffffffffffffffffLL)
  ,fMaxTracks(9999999)
  ,fCacheInp(300)
  ,fLearnSize(1)
  ,fBz(0)
  ,fResidualList()


  ,fMaxDY(6.f)
  ,fMaxDZ(6.f)
  ,fMaxQ2Pt(3.f)
  ,fMidQ2Pt(1.22f)
  ,fNY2XBins(15)
  ,fNZ2XBins(10)
  ,fNXBins(-1)
  ,fNDeltaYBins(120)
  ,fNDeltaZBins(120)  
  ,fDZ2X(0)
  ,fDX(0)
  ,fDZ2XI(0)
  ,fDXI(0)
  ,fDeltaYbinI(0)
  ,fDeltaZbinI(0)
  ,fNGVoxPerSector(0)
  //
  ,fMaxY2X(0)
  ,fDY2X(0)
  ,fDY2XI(0)
  ,fBinMinQ(0)
  ,fBinDQ(0)
  ,fBinDQI(0)

  ,fNMaxNeighb(0)
  ,fKernelType(kGausianKernel)
  

  ,fNTrSelTot(0)
  ,fNTrSelTotWO(0)
  ,fNReadCallTot(0)
  ,fNBytesReadTot(0)

  ,fVDriftParam(0)
  ,fVDriftGraph(0)
  ,fCorrTime(0)

  ,fStatTree(0)
  ,fHDelY(0)
  ,fHDelZ(0)

  ,fFilterOutliers(kTRUE) 
  ,fMaxSkippedCluster(10)  // 10 cluster
  ,fMaxRMSTrackCut(2.0)    // maximal RMS (cm) between the tracks 
  ,fMaxRMSClusterCut(0.3)    // maximal RMS (cm) between the cluster and local mean
  ,fMaxDeltaClusterCut(0.5)    // maximal delta(cm) between the cluster and local mean


{
  for (int i=0;i<kResDim;i++) {
    for (int j=0;i<2;j++) fNPCheb[i][j] = 15;
    fChebPrecD[i] = 100e-4;
  }

  memset(fLastSmoothingRes,0,kResDim*4*sizeof(double));
  fStepKern[kVoxQ] = 1;
  fStepKern[kVoxF] = 2;
  fStepKern[kVoxX] = 2;
  fStepKern[kVoxZ] = 4;
  Float_t  fKernelScaleEdge[kVoxDim] = {1, 1.,1.,1.}; // scaling factor for edge points

  for (int i=kVoxDim;i--;) fUniformBins[i] = kTRUE;

  for (int i=kVoxHDim;i--;) {
    fNBProdSt[i] = 0;
    fNBProdDY[i] = 0;
    fNBProdDZ[0] = 0;
  }
  fNBProdSectG[0] = fNBProdSectG[1] = 0;
  //
  for (int i=0;i<kNSect2;i++) {
    fSectGVoxRes[i] = 0;
    fTmpTree[i] = 0;
    fStatHist[i] = 0;
    farrNDstat[i] = 0;
    fTmpFile[i] = 0;
  }
}

//==================================================================================


void AliTPCDistortionExtraction::Init(int run
				      ,const char * residualList
				      ,Long64_t tmin
				      ,Long64_t tmax
				      ,float maxDY
				      ,float maxDZ
				      ,float maxQ2Pt
				      ,int nY2XBins
				      ,int nZ2XBins
				      ,int nXBins
				      ,int nDeltaBinsY
				      ,int nDeltaBinsZ
				      ,int maxTracks
				      ,Bool_t fixAlignmentBug
				      ,int cacheInp
				      ,int learnSize
				      ,Bool_t switchCache
				      )
{          
  AliSysInfo::AddStamp("ProjStart",0,0,0,0);
  fRun = run;
  if (fRun<1) {
    fRun = TString(gSystem->Getenv("runNumber")).Atoi();
    gSystem->Setenv("runNumber","245231"); // for test only
  }
  InitForBugFix();
  //
  //
  fResidualList = residualList;
  fTMin = tmin;
  fTMax = tmax;
  fMaxDY = maxDY;
  fMaxDZ = maxDZ;
  fMaxQ2Pt = maxQ2Pt;
  if (fMidQ2Pt<0) fMidQ2Pt = fMaxQ2Pt/2.f;
  //  
  //
  fNDeltaYBins = nDeltaBinsY;
  fNDeltaZBins = nDeltaBinsZ;
  if (fNDeltaYBins>0xff) {
    fNDeltaYBins = 0xff;
    AliErrorF("N DeltaY bins %d exceeds max allowed, setting to %d",nDeltaBinsY,fNDeltaYBins);
  }
  if (fNDeltaZBins>0xff) {
    fNDeltaZBins = 0xff;
    AliErrorF("N DeltaZ bins %d exceeds max allowed, setting to %d",nDeltaBinsZ,fNDeltaZBins);
  }
  fMaxTracks = maxTracks;
  fFixAlignmentBug = fixAlignmentBug;
  fCacheInp = cacheInp;
  fLearnSize = learnSize;
  fSwitchCache = switchCache;
  //
  fApplyZt2Zc = kFALSE;//kTRUE; //RRR
  //                                   
  // define boundaries
  InitBinning(nXBins,nY2XBins,nZ2XBins);
  
  //
  //
  LoadVDrift(); //!!!
  //
  // prepare aux info for stat and residuals histo bin calculation, see doc of TNDArray bin calculation
  fNBProdDY[kVoxHDim-1] = 1;
  fNBProdDZ[kVoxHDim-1] = 1;
  fNBProdSt[kVoxHDim-1] = 1;
  int nbh[kVoxDim];
  nbh[kVoxQ] = kNQBins;
  nbh[kVoxF] = fNY2XBins;
  nbh[kVoxX] = fNXBins;
  nbh[kVoxZ] = fNZ2XBins;
  for (int i=kVoxDim;i--;) {   // +2 to account for under/over-flows
    fNBProdSt[i] = fNBProdSt[i+1]*(2 + ((i==kVoxDim-1) ? kVoxHDim     : nbh[i+1])); 
    fNBProdDY[i] = fNBProdDY[i+1]*(2 + ((i==kVoxDim-1) ? fNDeltaYBins : nbh[i+1]));
    fNBProdDZ[i] = fNBProdDZ[i+1]*(2 + ((i==kVoxDim-1) ? fNDeltaZBins : nbh[i+1]));
  }
  //
  AliSysInfo::AddStamp("Init",0,0,0,0);
  //
  // kernel smoothing parameters
  SetKernelType();
  //
  fInitDone = kTRUE;
}

//_____________________________________________________
void AliTPCDistortionExtraction::CollectData() 
{
  // process residual trees and write local binned data trees
  const float kEps = 1e-6;
  if (!fInitDone) {printf("Init not done\n"); return;}
  //  gEnv->SetValue("TFile.AsyncPrefetching", 1);
  TVectorF *vecDY=0,*vecDZ=0,*vecZ=0,*vecR=0,*vecSec=0,*vecPhi=0, *vecDYITS=0;
  UShort_t npValid = 0;
  Int_t timeStamp = 0;
  AliExternalTrackParam* param = 0;
  //
  TVectorF *vecLocalDelta = new TVectorF(kNPadRows);
  TStopwatch swTot;
  swTot.Start();
  
  // temporary trees for local delta's storage
  dts_t dts, *dtsP = &dts; 
  //
  for (int is=0;is<kNSect2;is++) {
    tmpFile[is] = TFile::Open(Form("%s%d.root",kTmpFileName,is),"recreate");
    tmpTree[is] = new TTree(Form("ts%d",is),"");
    tmpTree[is]->Branch("dts",&dtsP);
    //tmpTree[is]->SetAutoFlush(150000);
    //
    statHist[is] = CreateVoxelStatHisto(is);
    arrNDstat[is] = (TNDArrayT<float>*)&statHist[is]->GetArray();
  }
  //
  // prepare input tree
  TString  chunkList = gSystem->GetFromPipe(TString::Format("cat %s",fResidualList.Data()).Data());
  TObjArray *chunkArray= chunkList.Tokenize("\n");  
  Int_t nChunks = chunkArray->GetEntriesFast();  
  //
  AliSysInfo::AddStamp("ProjInit",0,0,0,0);

  for (int ichunk=0;ichunk<nChunks;ichunk++) {
    //
    int ntrSelChunkWO=0, ntrSelChunk=0,nReadCallsChunk=0,nBytesReadChunk=0;
    //
    TStopwatch swc;
    swc.Start();
    TString fileNameString(chunkArray->At(ichunk)->GetName());
    if (fileNameString.Contains("alien://") && (!gGrid || (gGrid && !gGrid->IsConnected()))) TGrid::Connect("alien://");
    TFile *chunkFile = TFile::Open(fileNameString.Data());
    if (!chunkFile) continue;
    TTree *tree = (TTree*)chunkFile->Get("delta");
    if (!tree) {printf("No delta tree in %s\n",fileNameString.Data());continue;}
    tree->SetCacheLearnEntries(fLearnSize);
    tree->SetCacheSize(0);
    tree->SetCacheSize(fCacheInp*kMByte);
    //
    tree->SetBranchStatus("*",kFALSE);
    tree->SetBranchStatus("timeStamp",kTRUE);
    tree->SetBranchStatus("vecR.",kTRUE);
    tree->SetBranchStatus("vecSec.",kTRUE);
    tree->SetBranchStatus("vecPhi.",kTRUE);
    tree->SetBranchStatus("vecZ.",kTRUE);
    tree->SetBranchStatus("track.*",kTRUE);      
    tree->SetBranchStatus("npValid",kTRUE);
    tree->SetBranchStatus("trd0.",kTRUE);
    tree->SetBranchStatus("trd1.",kTRUE);
    //
    if (fFilterOutliers) {
      tree->SetBranchStatus("its0.",kTRUE);
      tree->SetBranchAddress("its0.",&vecDYITS);
    }
    //
    tree->SetBranchAddress("timeStamp",&timeStamp);
    tree->SetBranchAddress("vecR.",&vecR);
    tree->SetBranchAddress("vecSec.",&vecSec);
    tree->SetBranchAddress("vecPhi.",&vecPhi);
    tree->SetBranchAddress("vecZ.",&vecZ);
    tree->SetBranchAddress("track.",&param);
    tree->SetBranchAddress("npValid",&npValid);
    tree->SetBranchAddress("trd0.",&vecDY);
    tree->SetBranchAddress("trd1.",&vecDZ);
    //
    TBranch* brTime = tree->GetBranch("timeStamp");
    //
    int nTracks = tree->GetEntries();
    printf("Processing %d tracks of %s\n",nTracks,fileNameString.Data());

    float arrPhi[kNPadRows],arrX[kNPadRows],arrY[kNPadRows],arrZ[kNPadRows],arrDY[kNPadRows],arrDZ[kNPadRows];
    int arrSectID[kNPadRows], nCl=0;
    
    Bool_t lastReadMatched = kFALSE; // reset the cache when swithching between the timeStamp and Event read modes
    for (int itr=0;itr<nTracks;itr++) {
      nBytesReadChunk += brTime->GetEntry(itr);
      if (timeStamp<fTMin  || timeStamp>fTMax) {
	if (lastReadMatched && fSwitchCache) { // reset the cache
	  tree->SetCacheSize(0);
	  tree->SetCacheSize(fCacheInp*kMByte);
	  lastReadMatched = kFALSE;
	}
	continue;	
      }
      //
      if (!lastReadMatched && fSwitchCache) { // reset the cache before switching to event reading mode
	tree->SetCacheSize(0);
	tree->SetCacheSize(fCacheInp*kMByte);
      }
      lastReadMatched = kTRUE;
      nBytesReadChunk += tree->GetEntry(itr);
      float q2pt = param->GetParameter()[4];
      if (TMath::Abs(q2pt)>fMaxQ2Pt) continue;
      //
      const Float_t *vSec= vecSec->GetMatrixArray();
      const Float_t *vPhi= vecPhi->GetMatrixArray();
      const Float_t *vR  = vecR->GetMatrixArray();
      const Float_t *vZ  = vecZ->GetMatrixArray();
      const Float_t *vDY = vecDY->GetMatrixArray();
      const Float_t *vDZ = vecDZ->GetMatrixArray();
      //
      ntrSelChunkWO++;
      if (fFilterOutliers) { 
	// at the moment use Marian's implementation >>>>>
	Float_t rmsTrack=3, rmsCluster=1;
	Int_t nSkippedCluster=AliTPCcalibAlignInterpolation::CalculateDistance(*vecDY,*vecDYITS, *vecSec, *vecLocalDelta, npValid, rmsTrack, rmsCluster,1.5);
	if (nSkippedCluster>kMaxSkippedCluster) continue;
	if (rmsTrack>kMaxRMSTrackCut) continue;
	if (rmsCluster>kMaxRMSClusterCut) continue;	
      }
      //
      // at the moment use Marian's implementation <<<<<
      //
      fCorrTime = (fVDriftGraph!=NULL) ? fVDriftGraph->Eval(timeStamp):0; // for VDrift correction
      //
      ntrSelChunk++;
      nCl = 0;
      for (int ip=npValid;ip--;) { // 1st fill selected track data to buffer for eventual outlier rejection
	if (vR[ip]<1 || vDY[ip]<-900 ) continue;
	//
	arrX[nCl] = vR[ip];
	arrZ[nCl] = vZ[ip];
	arrDY[nCl] = vDY[ip];
	arrDZ[nCl] = vDZ[ip];
	arrPhi[nCl] = vPhi[ip];
	int rocID = TMath::Nint(vSec[ip]);
	//
	if (fFixAlignmentBug) {
	  AliTPCcalibAlignInterpolation::FixAlignmentBug(rocID, q2pt, fBz, arrPhi[nCl], arrX[nCl], arrZ[nCl], arrDY[nCl],arrDZ[nCl]);
	}
	if (arrPhi[nCl]<0) arrPhi[nCl] += 2.*TMath::Pi();

	if (TMath::Abs(arrDY[nCl])>fMaxDY-kEps) continue;
	//
	// apply drift velocity calibration if available
	arrDZ[nCl] += GetDriftCorrection(arrZ[nCl],arrX[nCl],arrPhi[nCl],rocID);

	if (TMath::Abs(arrDZ[nCl])>fMaxDZ-kEps) continue;
	//
	arrSectID[nCl] = rocID%kNSect2; // 0-36 for sectors from A0 to C17
	arrY[nCl] = arrX[nCl]*TMath::Sin(arrPhi[nCl]-(0.5f +rocID%kNSect)*kSecDPhi); // Y in sector frame
	arrX[nCl] *= TMath::Cos(arrPhi[nCl]-(0.5f +rocID%kNSect)*kSecDPhi); // R -> X in sector frame
	if (arrX[nCl]<kMinX || arrX[nCl]>kMaxX) continue;
	if (TMath::Abs(arrZ[nCl])>kZLim) continue;;
	nCl++;
      }
      // now fill the local trees and statistics
      float voxVars[kVoxHDim]={0}; // voxel variables (unbinned)
      for (int icl=nCl;icl--;) {
	int sectID = arrSectID[icl]; // 0-35 numbering
	// 
	// calculate voxel variables and bins
	// 
	if (!FindVoxelBin(sectID, q2pt, arrX[icl], arrY[icl], arrZ[icl], dts.bvox, voxVars)) continue;
	dts.dy   = (arrDY[icl]+fMaxDY)*fDeltaYbinI;
	dts.dz   = (arrDZ[icl]+fMaxDZ)*fDeltaZbinI;
	tmpTree[sectID]->Fill();
	//
	// fill statistics on distribution within the voxel, last dimension, kVoxV is for Nentries
	ULong64_t binToFill = GetBin2Fill(fNBProdSt,dts.bvox,kVoxV); // bin of sector stat histo
	float &binEntries = arrNDstat[sectID]->At(binToFill); // entries in the voxel
	float oldEntries  = binEntries++;
	float norm        = 1.f/binEntries;
	for (int iv=kVoxDim;iv--;) {
	  float &mean = arrNDstat[sectID]->At(binToFill+iv-kVoxV);
	  mean = ( mean*oldEntries + voxVars[iv]) * norm; // account new bin entry in averages calculation
	}
	//
      } // loop over clusters
    } // loop over tracks
    //
    swc.Stop();
    nReadCallsChunk =  chunkFile->GetReadCalls();
    printf("Selected %d tracks (%d with outliers) from chunk %d | %.1f MB read in %d read calls\n",
	   ntrSelChunk,ntrSelChunkWO, ichunk,float(nBytesReadChunk)/kMByte,nReadCallsChunk); swc.Print();
    fNTrSelTot += ntrSelChunk;
    fNTrSelTotWO += ntrSelChunkWO;
    fNReadCallTot += nReadCallsChunk;
    fNBytesReadTot += nBytesReadChunk;
    //
    delete tree;
    chunkFile->Close();
    delete chunkFile;
    //
    AliSysInfo::AddStamp("ProjTreeLoc", ichunk ,fNTrSelTot,fNTrSelTot,fNReadCallTot );
    //
    if (fNTrSelTot > fMaxTracks) {
      printf("Max number of tracks exceeded\n");
      break;
    }
    //
  } // loop over chunks
  //
  // write/close local trees
  for (int is=0;is<kNSect2;is++) {
    tmpFile[is]->cd();
    tmpTree[is]->Write("", TObject::kOverwrite);
    delete tmpTree[is];
    tmpFile[is]->Close();
    delete tmpFile[is];
  }
  //
  printf("Selected %d tracks (with outliers: %d) | %.1f MB read in %d read calls\n",
	 fNTrSelTot,fNTrSelTotWO,float(fNBytesReadTot)/kMByte,fNReadCallTot); 
  swTot.Print();

  AliSysInfo::AddStamp("ProjTreeLocSave");

  WriteStatHistos();
  //
  delete vecLocalDelta;
}

//________________________________________________
void AliTPCDistortionExtraction::ProcessResiduals()
{
  // project local trees, extract distortions
  if (!fInitDone) {printf("Init not done\n"); return;}

  bstat_t voxStat, *statP = &voxStat;
  AliSysInfo::AddStamp("ProcResid",0,0,0,0);
  TFile* flOut = new TFile(Form("%sTree.root",kStatOut),"recreate");
  fStatTree = new TTree("voxStat","");
  fStatTree->Branch("bins",&statP);
  //
  for (int i=0;i<kNEstPar;i++) {
    fStatTree->SetAlias(Form("Y%s",kEstName[i]),Form("distY[%d]",i));
    fStatTree->SetAlias(Form("Z%s",kEstName[i]),Form("distZ[%d]",i));
  }
  for (int i=0;i<kVoxDim;i++) {
    fStatTree->SetAlias(kVoxName[i],Form("bvox[%d]",i));
    fStatTree->SetAlias(Form("%sAV",kVoxName[i]),Form("stat[%d]",i));
  }
  fStatTree->SetAlias("N",Form("stat[%d]",kVoxV)); // entries
  //
  for (int is=0;is<kNSect2;is++) {
    ProcessSectorResiduals(is, voxStat);
    AliSysInfo::AddStamp("ProjResid",is);
    //
    TString sectFileName = Form("%s%d.root",kTmpFileName,is);
    ::Info(" AliTPCcalibAlignInterpolation::ProcessResidualsInTimeBin","Deleting %s\n",sectFileName.Data());
    //    unlink(sectFileName.Data());
  }
  //
  flOut->cd();
  fStatTree->Write("", TObject::kOverwrite);
  delete fStatTree;
  fStatTree = 0;
  flOut->Close();
  delete flOut;
  //
  AliSysInfo::AddStamp("ProcResid",1,0,0,0);
  //
}

//_________________________________
void  AliTPCDistortionExtraction::WriteVoxelDefinitions()
{
  // Store voxel boundaries
  if (!fInitDone) {printf("Init not done\n"); return;}

  voxDef_t vdef, *vdefP = &vdef;
  //
  TFile* flOut = new TFile(Form("%sTree.root",kStatOut),"update");
  //
  TTree* trDef = new TTree("voxDef","Voxel Boundaries definition");
  trDef->Branch("vDef",&vdef);
  for (int ix=0;ix<fNXBins;ix++) {
    vdef.bvox[kVoxX] = ix;
    vdef.vmin[kVoxX] = GetXLow(ix);
    vdef.vmax[kVoxX] = vdef.vmin[kVoxX] + GetDX(ix);
    //
    for (int ip=0;ip<fNY2XBins;ip++) {
      vdef.bvox[kVoxF] = ip;
      vdef.vmin[kVoxF] = GetY2XLow(ix,ip);
      vdef.vmax[kVoxF] = vdef.vmin[kVoxF] + GetDY2X(ix);
      //
      for (int iq=0;iq<kNQBins;iq++) {
	vdef.bvox[kVoxQ] = iq;
	int idxy = ix*fNY2XBins + ip;
	if (fBinDQI[idxy]>0) {
	  vdef.vmin[kVoxQ] = fBinMinQ[idxy] + iq*fBinDQ[idxy];
	  vdef.vmax[kVoxQ] = vdef.vmin[kVoxQ] + fBinDQ[idxy];
	}
	else {
	  vdef.vmin[kVoxQ] = -1;
	  vdef.vmax[kVoxQ] =  1;	  
	}
	//
	for (int iz=0;iz<fNZ2XBins;iz++) {
	  vdef.bvox[kVoxZ] = iz;
	  vdef.vmin[kVoxZ] = GetZ2XLow(iz);
	  vdef.vmax[kVoxZ] = vdef.vmin[kVoxZ] + GetDZ2X();
	  //
	  trDef->Fill();
	}
      }
    }
  }
  //
  for (int i=0;i<kVoxDim;i++) {
    trDef->SetAlias(kVoxName[i],Form("bvox[%d]",i));
    trDef->SetAlias(Form("%sMN",kVoxName[i]),Form("vmin[%d]",i));
    trDef->SetAlias(Form("%sMX",kVoxName[i]),Form("vmax[%d]",i));
  }
  //
  trDef->Write("", TObject::kOverwrite);
  delete trDef;
  //
  flOut->Close();
  delete flOut;
  //
  AliSysInfo::AddStamp("WriteVoxDef",0,0,0,0);
  //
}


//_________________________________________________
void AliTPCDistortionExtraction::ProcessSectorResiduals(int is, bstat_t &voxStat)
{
  // process residuals for single sector and store in the tree
  //
  TStopwatch sw;  sw.Start();
  printf("ProcessingSectoResiduals %d\n",is);
  AliSysInfo::AddStamp("ProcSectRes",is,0,0,0);
  //
  TString sectFileName = Form("%s%d.root",kTmpFileName,is);
  TFile* sectFile = TFile::Open(sectFileName.Data());
  if (!sectFile) {::Fatal("ProcessSectorResiduals",
			  "file %s not found",sectFileName.Data());}
  TString treeName = Form("ts%d",is);
  TTree *tmpTree = (TTree*) sectFile->Get(treeName.Data());
  if (!tmpTree) {::Fatal("ProcessSectorResiduals",
			 "tree %s is not found in file %s",treeName.Data(),sectFileName.Data());}
  //
  dts_t dts, *dtsP = &dts; 
  tmpTree->SetBranchAddress("dts",&dtsP);
  int npoints = tmpTree->GetEntries();
  //
  THnS *hisY = (THnS*)CreateSectorResidualsHisto(is, fNDeltaYBins, fMaxDY, "DY");
  THnS *hisZ = (THnS*)CreateSectorResidualsHisto(is, fNDeltaZBins, fMaxDZ, "DZ");  
  TNDArrayT<short>* ndArrY = (TNDArrayT<short>*)&hisY->GetArray();
  TNDArrayT<short>* ndArrZ = (TNDArrayT<short>*)&hisZ->GetArray();
  //
  for (int ip=0;ip<npoints;ip++) {
    tmpTree->GetEntry(ip);
    UShort_t  binDY = dts.dy;
    UShort_t  binDZ = dts.dz;
    ULong64_t binToFillY = GetBin2Fill(fNBProdDY,dts.bvox,binDY);
    ULong64_t binToFillZ = GetBin2Fill(fNBProdDZ,dts.bvox,binDZ);
    //
    ndArrY->At(binToFillY)++;
    ndArrZ->At(binToFillZ)++;
    //
  }
  //
  hisY->SetEntries(npoints);
  hisZ->SetEntries(npoints);
  //
  delete tmpTree;
  sectFile->Close(); // to reconsider: reuse the file
  delete sectFile;
  //
  AliSysInfo::AddStamp("ProjSectRes",is,0,0,0);
  //
  // extract full info about the voxel and write in the tree
  voxStat.bsec = is;
  ExtractVoxelData(voxStat, ndArrY, ndArrZ, arrNDstat[is]);

  // at the moment save histos ...
  TFile* flOut = TFile::Open(Form("residualSect%d.root",is),"update"); // RS shall we dump all sector histos to 1 file?
  hisY->Write("", TObject::kOverwrite);
  hisZ->Write("", TObject::kOverwrite);
  flOut->Close();
  //
  delete hisY;
  delete hisZ;
  //
  sw.Stop(); 
  sw.Print();
  AliSysInfo::AddStamp("ProcSectRes",is,1,0,0);
  //
}

//_________________________________________________
void AliTPCDistortionExtraction::ExtractVoxelData(bstat_t &stat, 
						  const TNDArrayT<short>* harrY, 
						  const TNDArrayT<short>* harrZ,
						  const TNDArrayT<float>* harrStat)
{
  // Extract distortion estimators from each voxel histo
  if (!fHDelY) fHDelY = new TH1F("dy","dy",fNDeltaYBins,-fMaxDY,fMaxDY);
  if (!fHDelZ) fHDelZ = new TH1F("dz","dz",fNDeltaZBins,-fMaxDZ,fMaxDZ);
  //
  UChar_t bvox1[kVoxDim];
  for (stat.bvox[kVoxZ]=0;stat.bvox[kVoxZ]<fNZ2XBins;stat.bvox[kVoxZ]++) {
    for (stat.bvox[kVoxX]=0;stat.bvox[kVoxX]<fNXBins;stat.bvox[kVoxX]++) { 
      for (stat.bvox[kVoxF]=0;stat.bvox[kVoxF]<fNY2XBins;stat.bvox[kVoxF]++) {
	//
	// in Z we integrate over all Q bins, since no dependence is expected
	stat.bvox[kVoxQ]=0;
	for (int i=kVoxDim;i--;) bvox1[i] = stat.bvox[i];
	bvox1[kVoxQ] = kNQBins-1;
	ExtractResidualHisto(harrZ,fNBProdDZ,stat.bvox,bvox1,fHDelZ); // integrate over Q bins
	ExtractDistortionsData(fHDelZ, stat.distZ);
	//	
	for (stat.bvox[kVoxQ]=0;stat.bvox[kVoxQ]<kNQBins;stat.bvox[kVoxQ]++) {
	  //
	  ExtractResidualHisto(harrY,fNBProdDY,stat.bvox,fHDelY);
	  ExtractDistortionsData(fHDelY, stat.distY);
	  //
	  //ExtractResidualHisto(harrZ,fNBProdDZ,stat.bvox,fHDelZ); // integrate over Q bins
	  //ExtractDistortionsData(fHDelZ, stat.distZ);
	  //
	  // extract voxel statistics: COG for each dimension
	  Long64_t bglo = GetBin2Fill(fNBProdSt,stat.bvox,kVoxV);
	  for (int i=0;i<kVoxHDim;i++) stat.stat[i] = harrStat->At(bglo+i-kVoxV);
	  fStatTree->Fill();
	}
      }
    }
  }
  //
}

//______________________________________________________________________________
void AliTPCDistortionExtraction::ExtractDistortionsData(TH1F* histo, float est[kNEstPar], 
							float minNorm, float fracLTM, float fitNSig)
{
  const float kMinEntries=30, kUseLLFrom=20;
  static TF1 fgaus("fgaus","gaus",-10,10);
  float nrm=0,mean=0,mom2=0,rms=0,maxVal=0;
  //
  memset(est,0,kNEstPar*sizeof(float));
  //
  float *w = (float*) histo->GetArray();
  w++; // skip underflows
  int nb = histo->GetNbinsX();
  float x = histo->GetXaxis()->GetXmax();
  float dx = (x+x)/nb;
  x -= dx*0.5f;
  for (int ip=nb;ip--;) {
    nrm  += w[ip];
    mean += x*w[ip];
    mom2 += x*x*w[ip];
    x -= dx;
    if (maxVal<w[ip]) maxVal = w[ip];
  }
  if (nrm>0) {
    mean /= nrm;
    mom2 /= nrm;
    rms = mom2 - mean*mean;
    rms = rms>0 ? TMath::Sqrt(rms):0;
  }
  est[kEstNorm] = nrm;
  est[kEstMean] = mean;
  est[kEstSig]  = rms;
  est[kEstMax]  = maxVal;
  if (nrm<minNorm) return;
  //
  TVectorF vecLTM(10);
  TStatToolkit::LTMHisto(histo, vecLTM, fracLTM); 
  est[kEstMeanL]  = vecLTM[1];
  est[kEstSigL]   = vecLTM[2];
  est[kEstMeanEL] = vecLTM[3];
  est[kEstSigEL]  = vecLTM[4];
  //
  if (nrm>=kMinEntries && est[kEstSigL]>0) {
    fgaus.SetParameters(nrm/(est[kEstMeanL]/dx),est[kEstMeanL],est[kEstSigL]);
    if (fitNSig>0) fgaus.SetRange(est[kEstMeanL]-fitNSig*est[kEstSigL], // limit fit range
				  est[kEstMeanL]+fitNSig*est[kEstSigL]);
    //
    TFitResultPtr fitPtr= histo->Fit(&fgaus,maxVal<kUseLLFrom ? "qnrlS":"qnrS");
    TFitResult * result = fitPtr.Get();
    if (result!=NULL) {
      est[kEstNormG] = fgaus.GetParameter(0);
      est[kEstMeanG] = fgaus.GetParameter(1);
      est[kEstSigG]  = fgaus.GetParameter(2);
      est[kEstMeanEG] = fgaus.GetParError(1);
      est[kEstSigEG]  = fgaus.GetParError(2);
      est[kEstChi2G] = fgaus.GetChisquare()/fgaus.GetNumberFreeParameters();
    }
  }
  //
}


//_________________________________________________
void AliTPCDistortionExtraction::InitForBugFix(const char* ocdb)
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


//___________________________________________________________________________
THnF* AliTPCDistortionExtraction::CreateVoxelStatHisto(int sect)
{
  // prepare histogram to store the means of distributions within the voxel

  // create binning for voxels statistics histograms
  Int_t    voxNBins[kVoxHDim];
  Double_t voxBinMin[kVoxHDim],voxBinMax[kVoxHDim];
  TString  voxAxisName[kVoxHDim];

  voxAxisName[kVoxQ] = "tgphi_Bin";
  voxNBins[kVoxQ]    = kNQBins;
  voxBinMin[kVoxQ]   = 0;
  voxBinMax[kVoxQ]   = kNQBins;
  //
  voxAxisName[kVoxF] = "Y2X_Bin";
  voxNBins[kVoxF]    = fNY2XBins;
  voxBinMin[kVoxF]   = 0;
  voxBinMax[kVoxF]   = fNY2XBins;
  //
  voxAxisName[kVoxX]   = "X_Bin";
  voxNBins[kVoxX]      = fNXBins;
  voxBinMin[kVoxX]     = 0;
  voxBinMax[kVoxX]     = fNXBins;
  //
  voxAxisName[kVoxZ]   = "Z2X_Bin";
  voxNBins[kVoxZ]      = fNZ2XBins;
  voxBinMin[kVoxZ]     = 0;
  voxBinMax[kVoxZ]     = fNZ2XBins;
  //
  voxAxisName[kVoxV] = "Stat_Bin";
  voxNBins[kVoxV]    = kVoxHDim;
  voxBinMin[kVoxV]   = 0;
  voxBinMax[kVoxV]   = kVoxHDim;
  //
  THnF* h = new THnF(Form("hs%d",sect),"",kVoxHDim,voxNBins,voxBinMin,voxBinMax);
  for (int i=0;i<kVoxHDim;i++) h->GetAxis(i)->SetName(voxAxisName[i].Data());
  h->SetEntries(1); // otherwise drawing does not work well
  return h;
}

//___________________________________________________________________________
THn* AliTPCDistortionExtraction::CreateSectorResidualsHisto(int sect, int nbDelta,float range, const char* pref)
{
  // prepare histogram to store the residuals withing the sector

  // create binning for voxels and residuals
  Int_t voxNBins[kVoxHDim];
  Double_t voxBinMin[kVoxHDim],voxBinMax[kVoxHDim];
  TString  voxAxisName[kVoxHDim];

  voxAxisName[kVoxQ] = "tgphi_Bin";
  voxNBins[kVoxQ]    = kNQBins;
  voxBinMin[kVoxQ]   = 0;
  voxBinMax[kVoxQ]   = kNQBins;
  //
  voxAxisName[kVoxF] = "Y2X_Bin";
  voxNBins[kVoxF]    = fNY2XBins;
  voxBinMin[kVoxF]   = 0;
  voxBinMax[kVoxF]   = fNY2XBins;
  //
  voxAxisName[kVoxX]   = "X_Bin";
  voxNBins[kVoxX]      = fNXBins;
  voxBinMin[kVoxX]     = 0;
  voxBinMax[kVoxX]     = fNXBins;
  //
  voxAxisName[kVoxZ]   = "Z2X_Bin";
  voxNBins[kVoxZ]      = fNZ2XBins;
  voxBinMin[kVoxZ]     = 0;
  voxBinMax[kVoxZ]     = fNZ2XBins;
  //
  //
  voxAxisName[kVoxV]   = pref;
  voxNBins[kVoxV]      = nbDelta;
  voxBinMin[kVoxV]     =-range;
  voxBinMax[kVoxV]     = range;
  //
  THnS* h = new THnS(Form("delta%d_%s",sect,pref),"",kVoxHDim,voxNBins,voxBinMin,voxBinMax);
  for (int i=0;i<kVoxHDim;i++) {
    h->GetAxis(i)->SetName(voxAxisName[i].Data());
    h->GetAxis(i)->SetTitle(voxAxisName[i].Data());
  }
  return h;
}


//_________________________________________
void AliTPCDistortionExtraction::LoadVDrift()
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
float AliTPCDistortionExtraction::GetDriftCorrection(float z, float x, float phi, int rocID)
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

//_____________________________________________________
float AliTPCDistortionExtraction::tgpXY(float x, float y, float q2p, float bz)
{
  // get the tg of primary track inclination wrt padrow given
  // that it was registered at X,Y sector coordinates
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


//=========================================================================
//
//   IO related methods
//
//=========================================================================
//
//___________________________________________________________________
void AliTPCDistortionExtraction::WriteStatHistos()
{
  // write stat histos
  TString statOutName = Form("%s.root",kStatOut);
  TFile* statOutFile = TFile::Open(statOutName.Data(),"recreate");
  for (int is=0;is<kNSect2;is++) statHist[is]->Write("", TObject::kOverwrite);
  statOutFile->Close();
  delete statOutFile;
  //
}

//___________________________________________________________________
void AliTPCDistortionExtraction::LoadStatHistos()
{
  // load bin stat histos
  TString statOutName = Form("%s.root",kStatOut);
  TFile* statOutFile = TFile::Open(statOutName.Data());
  if (!statOutFile) ::Fatal("tstw","LoadStatHistos: failed to read file %s",statOutName.Data()); 
  for (int is=0;is<kNSect2;is++) {
    statHist[is] = (THnF*) statOutFile->Get(Form("hs%d",is));
    if (!statHist[is]) ::Fatal("tstw","LoadStatHistos: failed to read secto %d histo from %s",is,statOutName.Data());
    arrNDstat[is] = (TNDArrayT<float>*)&statHist[is]->GetArray();
  }
  statOutFile->Close();
  delete statOutFile;
  //
}

//___________________________________________________________________
void AliTPCDistortionExtraction::WriteResTree()
{
  // output file for results tree
  TStopwatch sw;
  sw.Start();
  bres_t voxRes, *voxResP=&voxRes;

  AliSysInfo::AddStamp("ResTree",0,0,0,0);

  TFile* flOut = new TFile(Form("%sTree.root",kResOut),"recreate");
  TTree* resTree = new TTree("voxRes","final distortions, see GetListOfAliases");
  resTree->Branch("res",&voxRes);
  for (int i=0;i<kVoxDim;i++) {
    resTree->SetAlias(kVoxName[i],Form("bvox[%d]",i));
    resTree->SetAlias(Form("%sAV",kVoxName[i]),Form("stat[%d]",i));
  }
  for (int i=0;i<kResDim;i++) {
    resTree->SetAlias(kResName[i],Form("D[%d]",i));
    resTree->SetAlias(Form("%sE",kResName[i]),Form("E[%d]",i));
  }
  //
  for (int is=0;is<kNSect2;is++) { 

    bres_t* sectData = fSectGVoxRes[is];
    if (!sectData) {
      printf("No processed data for sector %d\n",is);
      continue;
    }

    // now store the data for the sector
    for (int iz=0;iz<fNZ2XBins;iz++) {
      for (int ix=0;ix<fNXBins;ix++) { 
	for (int ip=0;ip<fNY2XBins;ip++) {
	  int binGlo = GetVoxGBin(ix,ip,iz);
	  bres_t *voxel = &sectData[binGlo];
	  if (fChebCorr) {
	    int row = GetRowID(voxel->stat[kVoxX]);
	    int roc = is;
	    if (row>=kNRowIROC) {roc += kNSect2; row -= kNRowIROC;}
	    fChebCorr->Eval(roc,row, voxel->stat[kVoxF],voxel->stat[kVoxZ],voxel->DC);
	  }
	  memcpy(&voxRes,voxel,sizeof(bres_t)); // store in the sector data array
	  resTree->Fill();
	}
      }
    }    
    //
  } // end of sector loop
  //
  flOut->cd();
  resTree->Write("", TObject::kOverwrite);
  delete resTree;
  //
  if (fChebCorr) {
    fChebCorr->Write();
  }
  //
  flOut->Close();
  delete flOut;
  //
  sw.Stop();
  sw.Print();
  AliSysInfo::AddStamp("ResTree",1,0,0,0);

}



//=========================================================================
//
//   voxels processing related methods
//
//=========================================================================

//___________________________________________________________________________________
float AliTPCDistortionExtraction::ExtractResidualHisto(const TNDArrayT<short>* harr, const Long64_t bprod[kVoxHDim], 
						       const UChar_t vox[kVoxDim], TH1F* dest)
{
  // extract residuals for the voxel from THn to 1D histo
  // Here we use shortcut to fill the histo, don't modify neither this part, not the histo def.
  float* arr = (float*)dest->GetArray();
  int nb = dest->GetNbinsX();
  arr++; // skip the underflows bin
  ULong64_t binGlo = GetBin2Fill(bprod,vox,0);
  float sum = 0;
  for (int i=0;i<nb;i++) sum += arr[i] = harr->At(binGlo++);
  dest->SetEntries(sum);
  return sum;
}

//___________________________________________________________________________________
float AliTPCDistortionExtraction::ExtractResidualHisto(const TNDArrayT<short>* harr, const Long64_t bprod[kVoxHDim], 
						       const UChar_t voxMin[kVoxDim], const UChar_t voxMax[kVoxDim], TH1F* dest)
{
  // extract residuals for the voxel from THn to 1D histo
  // Here we use shortcut to fill the histo, don't modify neither this part, not the histo def.
  float* arr = (float*)dest->GetArray();
  int nb = dest->GetNbinsX();
  memset(arr,0,sizeof(float)*(nb+2));
  arr++; // skip the underflows bin
  float sum = 0;
  UChar_t vox[kVoxDim];
  
  for (vox[kVoxQ]=voxMin[kVoxQ];vox[kVoxQ]<=voxMax[kVoxQ];vox[kVoxQ]++) 
    for (vox[kVoxF]=voxMin[kVoxF];vox[kVoxF]<=voxMax[kVoxF];vox[kVoxF]++) 
      for (vox[kVoxX]=voxMin[kVoxX];vox[kVoxX]<=voxMax[kVoxX];vox[kVoxX]++) 
	for (vox[kVoxZ]=voxMin[kVoxZ];vox[kVoxZ]<=voxMax[kVoxZ];vox[kVoxZ]++) {
	  ULong64_t binGlo = GetBin2Fill(bprod,vox,0);
	  for (int i=0;i<nb;i++) arr[i] += harr->At(binGlo++);
	}
  //
  for (int i=0;i<nb;i++) sum += arr[i];
  dest->SetEntries(sum);
  return sum;
}

//___________________________________________________________________________________
TH1F* AliTPCDistortionExtraction::ExtractResidualHisto(Bool_t y, int sect, const UChar_t vox[kVoxDim])
{
  // create Y or Z residuals for the voxel in sector sect,  from THn to 1D histo
  TString fln = Form("residualSect%d.root",sect);
  TFile* fl = TFile::Open(fln.Data());
  if (!fl) {printf("Failed to open %s\n",fln.Data()); return 0;}
  THnS* hn = (THnS*)fl->Get(Form("delta%d_D%s",sect,y ? "Y":"Z"));
  if (!hn) {printf("Did not find histo %s in %s\n",Form("delta%d_D%s",sect,y ? "Y":"Z"),fln.Data());return 0;}
  TString ttl = Form("d%s_S%d_vox%d_%d_%d_%d",y ? "Y":"Z",sect,vox[kVoxQ],vox[kVoxF],vox[kVoxX],vox[kVoxZ]);
  TH1F* h1 = new TH1F(ttl.Data(),ttl.Data(),
		      y ? fNDeltaYBins : fNDeltaZBins ,
		      y ? -fMaxDY : -fMaxDZ,
		      y ?  fMaxDY :  fMaxDZ);
  h1->SetDirectory(0);
  const TNDArrayT<short>* harr = (TNDArrayT<short>*)&hn->GetArray();
  ExtractResidualHisto(harr, y ? fNBProdDY : fNBProdDZ, vox, h1);
  delete hn;
  fl->Close();
  delete fl;
  printf("Don't forget to delete (TH1F*)%p %s\n",h1,h1->GetName());
  return h1;
}

//___________________________________________________________________________________
TH1F* AliTPCDistortionExtraction::ExtractResidualHisto(Bool_t y, int sect, const UChar_t vox[kVoxDim], const UChar_t vox1[kVoxDim])
{
  // create Y or Z residuals for the voxel in sector sect,  from THn to 1D histo
  TString fln = Form("residualSect%d.root",sect);
  TFile* fl = TFile::Open(fln.Data());
  if (!fl) {printf("Failed to open %s\n",fln.Data()); return 0;}
  THnS* hn = (THnS*)fl->Get(Form("delta%d_D%s",sect,y ? "Y":"Z"));
  if (!hn) {printf("Did not find histo %s in %s\n",Form("delta%d_D%s",sect,y ? "Y":"Z"),fln.Data());return 0;}
  TString ttl = Form("d%s_S%d_vox%d_%d_%d_%d__%d_%d_%d_%d",y ? "Y":"Z",sect,
		     vox[kVoxQ],vox[kVoxF],vox[kVoxX],vox[kVoxZ],
		     vox1[kVoxQ],vox1[kVoxF],vox1[kVoxX],vox1[kVoxZ]);
  TH1F* h1 = new TH1F(ttl.Data(),ttl.Data(),
		      y ? fNDeltaYBins : fNDeltaZBins ,
		      y ? -fMaxDY : -fMaxDZ,
		      y ?  fMaxDY :  fMaxDZ);
  h1->SetDirectory(0);
  const TNDArrayT<short>* harr = (TNDArrayT<short>*)&hn->GetArray();
  ExtractResidualHisto(harr, y ? fNBProdDY : fNBProdDZ, vox, vox1, h1);
  delete hn;
  fl->Close();
  delete fl;
  printf("Don't forget to delete (TH1F*)%p %s\n",h1,h1->GetName());
  return h1;
}

//__________________________________________________________________
void AliTPCDistortionExtraction::ExtractXYZDistortions()
{
  if (!fInitDone) {printf("Init not done\n"); return;}

  TStopwatch sw;
  // extract XYZ distortions from fitted Y,Z residuals vs Q variable
  bstat_t voxStat, *statP = &voxStat;
  bstat_t voxIQ[kNQBins];
  bres_t voxRes, *voxResP=&voxRes;
  //
  AliSysInfo::AddStamp("ExtXYZ",0,0,0,0);
  printf("ExtractXYZDistortions\n");
  TFile* flStat = 0;
  if (!fStatTree) {
    TString fname = Form("%sTree.root",kStatOut);
    flStat = new TFile(fname.Data());
    if (!flStat) {::Fatal("ExtractXYZDistortions","file %s not found",fname.Data());}
    fStatTree = (TTree*)flStat->Get("voxStat");
    if (!fStatTree) {::Fatal("ExtractXYZDistortions","voxStat tree not found in %s",fname.Data());}
  }
  fStatTree->SetBranchAddress("bins",&statP);
  
  //
  int ent = 0;

  // 1st loop over good voxels, extract X distortions
  for (int is=0;is<kNSect2;is++) { 
    
    bres_t* sectData = fSectGVoxRes[is] = new bres_t[fNGVoxPerSector]; // here we keep main result
    voxRes.bsec = is;
    float cntGood = 0;
    for (voxRes.bvox[kVoxZ]=0;voxRes.bvox[kVoxZ]<fNZ2XBins;voxRes.bvox[kVoxZ]++) {
      for (voxRes.bvox[kVoxX]=0;voxRes.bvox[kVoxX]<fNXBins;voxRes.bvox[kVoxX]++) { 
	for (voxRes.bvox[kVoxF]=0;voxRes.bvox[kVoxF]<fNY2XBins;voxRes.bvox[kVoxF]++) {
	  for (int iq=0;iq<kNQBins;iq++) {
	    fStatTree->GetEntry(ent++);
	    //
	    // check
	    if (voxStat.bvox[kVoxZ]!=voxRes.bvox[kVoxZ] || voxStat.bvox[kVoxX]!=voxRes.bvox[kVoxX] ||
		voxStat.bvox[kVoxF]!=voxRes.bvox[kVoxF] || voxStat.bvox[kVoxQ]!=iq) {
	      printf("voxel QFXZ : Expected %d %2d %3d %2d | Read %d %2d %3d %2d",iq,
		     voxRes.bvox[kVoxF],voxRes.bvox[kVoxX],voxRes.bvox[kVoxZ],
		     voxStat.bvox[kVoxQ],voxStat.bvox[kVoxF],voxStat.bvox[kVoxX], voxStat.bvox[kVoxZ]);
	      ::Fatal("ExtractXYZDistortions","Mismatch between expected and obtained voxel %d",ent-1);
	    }
	    memcpy(&voxIQ[iq],&voxStat,sizeof(bstat_t)); // save for the analysis vs Q
	  }
	  ExtractVoxelXYZDistortions(voxIQ,voxRes);      // extract residuals deconvoluted for the slopes etc
	  cntGood += voxRes.bvox[kVoxQ]>0;
	  int binGlo = GetVoxGBin(voxRes.bvox);
	  memcpy(&sectData[binGlo],&voxRes,sizeof(bres_t)); // store in the sector data array
	  //
	}
      }
    } 
    //
    int cntSmooth = Smooth0(is); // smooth sector data
    //FillHoles(is, sectData, fNBProdSectG);

    printf("Sector%2d: voxels with data %6d (%4.1f%%) smoothed %6d (%4.1f%%) of %d\n",is,int(cntGood),
	   cntGood/fNGVoxPerSector*100.,cntSmooth,float(cntSmooth)/fNGVoxPerSector*100.,fNGVoxPerSector);
  }

  delete fStatTree;
  fStatTree = 0;
  flStat->Close();
  delete flStat;
  //

}

//_____________________________________________
Bool_t AliTPCDistortionExtraction::ExtractVoxelXYZDistortions(const bstat_t voxIQ[kNQBins], 
							      bres_t &res, int minStat, 
							      float maxGChi2, int minYBinsOK)
{
  // extract XYZ distortions from voxel fitted Y,Z residuals vs Q variable
  //
  const float kMaxSigG2L = 5.0f; // max ratio between Gaussian sigma and RMS LTM
  const float kMinNormG2M = 0.2f; // min ratio between Gaussian amplitude and max Value
  const float kZeroSigma = 1e-4; 
  Bool_t okG=kFALSE,okL;
  //
  int nyOK = 0;
  float av[kNQBins],meas[kNQBins],wgh[kNQBins],resFit[2],errFit[3];
  //
  for (int i=kResDim;i--;) res.D[i] = res.DS[i] = res.DC[i] = res.E[i] = 0;
  for (int i=kVoxHDim;i--;) res.stat[i] = 0;
  //
  for (int iq=kNQBins;iq--;) {
    const bstat_t &vox = voxIQ[iq];
    float ent = vox.stat[kVoxV];
    if (ent<minStat) continue;
    //
    okG = okL = vox.distY[kEstSigL]>kZeroSigma;
    //
    if (okG && 
	vox.distY[kEstNormG]<kMinNormG2M*vox.distY[kEstMax] || // gaussian norm should not be negligible
	vox.distY[kEstChi2G]>maxGChi2 || vox.distY[kEstSigG]<kZeroSigma ||
	vox.distY[kEstSigG]>kMaxSigG2L) okG = kFALSE;
    //
    // assume that measured Y resodual dy is related to real residuals DY and DX as
    // dy = DY - DX*tg(slope) 
    // where the slope is average track inclination angle at the pad-row
    // For DZ calculate simple weighted mean
    //
    if (okL) { // collect data for linear fit
      av[nyOK] = vox.stat[kVoxQ]; // mean value of tg(slope) for this Q bin
      meas[nyOK] = okG ? vox.distY[kEstMeanG]  : vox.distY[kEstMeanL];
      float wy   = okG ? vox.distY[kEstMeanEG] : vox.distY[kEstMeanEL];
      wgh[nyOK]  = 1./(wy*wy);
      //
      float st = vox.stat[kVoxV]; // statistics of the voxel
      res.stat[kVoxV] += st;
      for (int i=kVoxDim;i--;) res.stat[i] += vox.stat[i]*st;
      nyOK++;
    }
    //
  }
  //
  if (res.stat[kVoxV]>0) {
    float stI = 1.0f/res.stat[kVoxV];
    for (int i=kVoxDim;i--;) res.stat[i] *= stI; // average of each voxel dimension for selected bins    
  }
  if (nyOK>=minYBinsOK && FitPoly1(av,meas,wgh,nyOK,resFit,errFit)) {
    //
    res.D[kResY] = resFit[0];
    res.D[kResX] =-resFit[1];
    res.E[kResY] = errFit[0];
    res.E[kResX] = errFit[2];
    res.bvox[kVoxQ] = nyOK; // number of points used
  }
  else { // estimation impossible
    res.bvox[kVoxQ] = 0;
    // set coordinates to bin center
    GetVoxelCoordinates(res.bsec,res.bvox[kVoxX],res.bvox[kVoxF],res.bvox[kVoxZ],
			res.stat[kVoxX],res.stat[kVoxF],res.stat[kVoxZ]);
  }
  //
  // Z fits were integrated over Q, use just 1st bin (they are all the same), correcting for X shift
  // as measured DZ -> dZ + DX*<Z/X>
  const bstat_t &vox0 = voxIQ[0];
  okG = okL = vox0.distZ[kEstSigL]>kZeroSigma;  
  if (okG && 
      vox0.distZ[kEstNormG]<kMinNormG2M*vox0.distZ[kEstMax] || // gaussian norm should not be negligible
      vox0.distZ[kEstChi2G]>maxGChi2 || vox0.distZ[kEstSigG]<kZeroSigma ||
      vox0.distZ[kEstSigG]>kMaxSigG2L) okG = kFALSE;
    //
  if (okL) {
    res.D[kResZ] = (okG ? vox0.distZ[kEstMeanG]  : vox0.distZ[kEstMeanL]) + res.stat[kVoxZ]*res.D[kResX];
    res.E[kResZ] = okG ? vox0.distZ[kEstMeanEG] : vox0.distZ[kEstMeanEL];
  }
  //
  for (int i=0;i<kResDim;i++) res.E[i] = res.E[i]>0 ? TMath::Sqrt(res.E[i]) : 0;
  //
}


//=========================================================================
//
//   fitting/smoothing related methods
//
//=========================================================================
//________________________________
void AliTPCDistortionExtraction::FillHoles(int isect, bres_t *sectData, const int fNBProdSectG[2], int minGoodPoints)
{
  /// RRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR DELETE?
  // fill holes within 1 sector
  // scan transverse plane and extra/interpolate in Z/X 
  const float kZeroError=1e-6, kDummyError = 1.0;
  float resFit0[3],errFit0[6],resFit1[3],errFit1[6];
  float dz = GetDZ2X();
  int maxDim = TMath::Max(fNXBins,TMath::Max(fNY2XBins,fNZ2XBins));
  bres_t **currLine = new bres_t*[maxDim];
  Bool_t missY[maxDim],missZ[maxDim];
  Float_t val0[maxDim],pos0[maxDim],wgh0[maxDim];
  Float_t val1[maxDim],pos1[maxDim],wgh1[maxDim];
  printf("FillHoles Sector%d\n",isect);
  for (int ix=0;ix<fNXBins;ix++) {
    for (int ip=0;ip<fNY2XBins;ip++) {
      int nmissY=0,nmissZ=0;

      for (int iz=0;iz<fNZ2XBins;iz++) {  // extract line in z
	int binGlo = iz + fNBProdSectG[1]*ip + fNBProdSectG[0]*ix; // global bin
	currLine[iz] = &sectData[binGlo];
	missY[iz] = missZ[iz] = kFALSE;
	if (currLine[iz]->E[kResY]<kZeroError) {
	  missY[iz] = kTRUE;
	  nmissY++;
	}
	if (currLine[iz]->E[kResZ]<kZeroError) {
	  missZ[iz] = kTRUE;
	  nmissZ++;
	}
      }
      // 
      if (nmissY && (fNZ2XBins-nmissY)>=minGoodPoints) { // recover points in Y and X
	int npGood = 0;
	for (int iz=0;iz<fNZ2XBins;iz++) {
	  if (missY[iz]) continue;
	  val0[npGood] = currLine[iz]->D[kResY];
	  pos0[npGood] = currLine[iz]->stat[kVoxZ]; // average Z position
	  wgh0[npGood] = 1./(currLine[iz]->E[kResY]*currLine[iz]->E[kResY]);
	  // X distortion goes with Y
	  val1[npGood] = currLine[iz]->D[kResX];
	  pos1[npGood] = currLine[iz]->stat[kVoxZ]; // average Z position
	  wgh1[npGood] = 1./(currLine[iz]->E[kResX]*currLine[iz]->E[kResX]);
	  npGood++;
	}
	Bool_t res = FitPoly2(pos0,val0,wgh0,npGood, resFit0, errFit0)
	  &&         FitPoly2(pos1,val1,wgh1,npGood, resFit1, errFit1);
	if (res) {
	  for (int iz=0;iz<fNZ2XBins;iz++) {
	    if (!missY[iz]) continue;
	    // evaluate in the center of the bin
	    double z = (isect>=kNSect ? -1.0f:1.0f)*(iz+0.5)*dz, z2=z*z, z3=z2*z, z4=z3*z;
	    currLine[iz]->D[kResY] = resFit0[0]+z*resFit0[1]+z2*resFit0[2];
	    currLine[iz]->D[kResX] = resFit1[0]+z*resFit1[1]+z2*resFit1[2];
	    double evErr0= errFit0[0] + errFit0[2]*z2 + errFit0[5]*z4
	      +            2.*(errFit0[1]*z + errFit0[3]*z2 + errFit0[4]*z3);
	    double evErr1= errFit1[0] + errFit1[2]*z2 + errFit1[5]*z4
	      +            2.*(errFit1[1]*z + errFit1[3]*z2 + errFit1[4]*z3);
	    //
	    currLine[iz]->E[kResY] = evErr0>0 ? TMath::Sqrt(evErr0) : kDummyError;
	    currLine[iz]->E[kResX] = evErr1>0 ? TMath::Sqrt(evErr1) : kDummyError;
	    //		
	  }
	  printf("Sect%2d bX=%3d bF=%3d DY vs Z: filled %d holes using %d values\n",isect,ix,ip, nmissY,npGood);
	}
	else printf("Sect%2d bX=%3d bF=%3d DY vs Z: FAILED to fill %d holes using %d values\n",isect,ix,ip, nmissY,npGood);	
      }
      //
      //
      if (nmissZ && (fNZ2XBins-nmissZ)>=minGoodPoints) { // recover points in Z
	int npGood = 0;
	for (int iz=0;iz<fNZ2XBins;iz++) {
	  if (missY[iz]) continue;
	  val0[npGood] = currLine[iz]->D[kResZ];
	  pos0[npGood] = currLine[iz]->stat[kVoxZ]; // average Z position
	  wgh0[npGood] = 1./(currLine[iz]->E[kResZ]*currLine[iz]->E[kResZ]);
	  npGood++;
	}
	Bool_t res = FitPoly2(pos0,val0,wgh0,npGood, resFit0, errFit0);
	if (res) {
	  for (int iz=0;iz<fNZ2XBins;iz++) {
	    if (!missZ[iz]) continue;
	    // evaluate in the center of the bin
	    double z = (isect>=kNSect ? -1.0f:1.0f)*(iz+0.5)*dz, z2=z*z, z3=z2*z, z4=z3*z;
	    currLine[iz]->D[kResZ] = resFit0[0]+z*resFit0[1]+z2*resFit0[2];
	    double evErr = errFit0[0] + errFit0[2]*z2 + errFit0[5]*z4
	      +            2.*(errFit0[1]*z + errFit0[3]*z2 + errFit0[4]*z3);
	    currLine[iz]->E[kResZ] = evErr>0 ? TMath::Sqrt(evErr) : kDummyError;
	  }
	  printf("Sect%2d bX=%3d bF=%3d DZ vs Z: filled %d holes using %d values\n",isect,ix,ip, nmissZ,npGood);
	}
	else printf("Sect%2d bX=%3d bF=%3d DZ vs Z: FAILED to fill %d holes using %d values\n",isect,ix,ip, nmissY,npGood);	
      }
      //
    } // loop in phi bins
  } // loop in x bins
  delete[] currLine;
  //
}

//_____________________________________________________
Bool_t AliTPCDistortionExtraction::FitPoly2(const float* x,const float* y, const float* w, int np, float *res, float *err)
{
  // poly2 fitter
  if (np<3) return kFALSE; // no enough points
  double sumW[5]={0},sumY[3]={0};
  for (int ip=np;ip--;) {
    double ww = w[ip];
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
Bool_t AliTPCDistortionExtraction::FitPoly1(const float* x,const float* y, const float* w, int np, float *res, float *err)
{
  // poly1 fitter
  if (np<2) return kFALSE; // no enough points
  double sumW[3]={0},sumY[2]={0};
  for (int ip=np;ip--;) {
    double ww = w[ip];
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

//________________________________
Int_t AliTPCDistortionExtraction::Smooth0(int isect)
{
  // apply linear regression kernel smoother 
  int cnt = 0;
  bres_t* sectData = fSectGVoxRes[isect];
  for (int ix=0;ix<fNXBins;ix++) {
    for (int ip=0;ip<fNY2XBins;ip++) {
      for (int iz=0;iz<fNZ2XBins;iz++) {  // extract line in z
	int binGlo = GetVoxGBin(ix,ip,iz);
	bres_t *vox = &sectData[binGlo];
	Bool_t res = GetSmoothEstimate(vox->bsec,vox->stat[kVoxX],vox->stat[kVoxF],vox->stat[kVoxZ],vox->DS);
	vox->smooth = res;
	if (res) cnt++;
      }
    }
  }
  //
  return cnt;
}

//________________________________________________________________
Bool_t AliTPCDistortionExtraction::GetSmoothEstimate(int isect, float x, float p, float z, float *res, float *deriv)
{
  // get smooth estimate for point in sector coordinates (x,y/x,z/x)
  // smoothing results also saved in the fLastSmoothingRes (allow derivative calculation)
  //
  const int kMinPointsTot = 4; // we fit 12 paremeters, each point provides 3 values
  const int kMaxTrials = 3; // max allowed iterations if neighbours are missing

  res[kResX]=res[kResY]=res[kResZ] = 0;
  //
  double cmatX[10],cmatY[10],cmatZ[10];
  double &m00x=cmatX[0], 
    &m10x=cmatX[1], &m11x=cmatX[2], 
    &m20x=cmatX[3], &m21x=cmatX[4], &m22x=cmatX[5], 
    &m30x=cmatX[6], &m31x=cmatX[7], &m32x=cmatX[8], &m33x=cmatX[9];
  double &m00y=cmatY[0], 
    &m10y=cmatY[1], &m11y=cmatY[2], 
    &m20y=cmatY[3], &m21y=cmatY[4], &m22y=cmatY[5], 
    &m30y=cmatY[6], &m31y=cmatY[7], &m32y=cmatY[8], &m33y=cmatY[9];
  double &m00z=cmatZ[0], 
    &m10z=cmatZ[1], &m11z=cmatZ[2], 
    &m20z=cmatZ[3], &m21z=cmatZ[4], &m22z=cmatZ[5], 
    &m30z=cmatZ[6], &m31z=cmatZ[7], &m32z=cmatZ[8], &m33z=cmatZ[9];

  //loop over neighbours which can contribute
  //
  bres_t *currClus[fNMaxNeighb];
  double *rhsX = &fLastSmoothingRes[0], *rhsY = &fLastSmoothingRes[4] , *rhsZ = &fLastSmoothingRes[8];
  //
  int ix0,ip0,iz0;
  FindVoxel(x,p,z, ix0,ip0,iz0); // find nearest voxel
  bres_t* sectData = fSectGVoxRes[isect];
  int binCen = GetVoxGBin(ix0,ip0,iz0);  // global bin of nearest voxel
  bres_t* voxCen = &sectData[binCen]; // nearest voxel
  //
  int trial = 0, nbOK = 0;
  while(1)  {
    //
    memset(fLastSmoothingRes,0,kResDim*4*sizeof(double));
    if (trial>kMaxTrials) {printf("Trials limit reached\n"); return kFALSE;}

    memset(cmatX,0,10*sizeof(double));
    if (fUseErrInSmoothing) { // in this case we cannot use the same matrix for x,y,z
      memset(cmatY,0,10*sizeof(double));
      memset(cmatZ,0,10*sizeof(double));
    }
    //
    nbOK=0; // accounted neighbours
    //
    int stepX = fStepKern[kVoxX] + trial;
    int stepF = fStepKern[kVoxF] + trial;
    int stepZ = fStepKern[kVoxZ] + trial;
    //
    if (!voxCen->bvox[kVoxQ]) { // closest voxel has no data, increase smoothing step
      stepX+=1;
      stepF+=1;
      stepZ+=1;
    }
    //
    // effective kernel widths accounting for the increased bandwidth at the edges
    float kWXI = GetDXI(ix0)/stepX;
    float kWFI = GetDY2XI(ix0)/stepF;
    float kWZI = GetDZ2XI()/stepZ;

    // for edge bins increase kernel size and neighbours search
    int ixMn=ix0-stepX,ixMx=ix0+stepX;
    if (ixMn<0) {
      ixMn = 0;
      ixMx = TMath::Min(TMath::Nint(ix0+stepX*fKernelScaleEdge[kVoxX]),fNXBins-1);
      kWXI /= fKernelScaleEdge[kVoxX];
    }
    if (ixMx>=fNXBins) {
      ixMx = fNXBins-1;
      ixMn = TMath::Max(TMath::Nint(ix0-stepX*fKernelScaleEdge[kVoxX]),0);
      kWXI /= fKernelScaleEdge[kVoxX];
    }
    //
    int ipMn=ip0-stepF,ipMx=ip0+stepF;
    if (ipMn<0) {
      ipMn = 0;
      ipMx = TMath::Min(TMath::Nint(ip0+stepF*fKernelScaleEdge[kVoxF]),fNY2XBins-1);
      kWFI /= fKernelScaleEdge[kVoxF];
    }
    if (ipMx>=fNY2XBins) {
      ipMx = fNY2XBins-1; 
      ipMn = TMath::Max(TMath::Nint(ip0-stepF*fKernelScaleEdge[kVoxF]),0);
      kWFI /= fKernelScaleEdge[kVoxF];
    }
    //
    int izMn=iz0-stepZ,izMx=iz0+stepZ;
    if (izMn<0) {
      izMn = 0;
      izMx = TMath::Min(TMath::Nint(iz0+stepZ*fKernelScaleEdge[kVoxZ]),fNZ2XBins-1);
      kWZI /= fKernelScaleEdge[kVoxZ];
    }
    if (izMx>=fNZ2XBins) {
      izMx = fNZ2XBins-1;
      izMn = TMath::Max(TMath::Nint(iz0-stepZ*fKernelScaleEdge[kVoxZ]),0);
      kWZI /= fKernelScaleEdge[kVoxZ];
    }
    //
    int nOccX[ixMx-ixMn+1],nOccF[ipMx-ipMn+1],nOccZ[izMx-izMn+1];
    //
    for (int i=ixMx-ixMn+1;i--;) nOccX[i]=0;
    for (int i=ipMx-ipMn+1;i--;) nOccF[i]=0;
    for (int i=izMx-izMn+1;i--;) nOccZ[i]=0;

    for (int ix=ixMn;ix<=ixMx;ix++) {
      for (int ip=ipMn;ip<=ipMx;ip++) {
	for (int iz=izMn;iz<=izMx;iz++) {
	  //
	  int binNb = GetVoxGBin(ix,ip,iz);  // global bin
	  bres_t* voxNb = &sectData[binNb];
	  if (!voxNb->bvox[kVoxQ]) continue; // skip voxels w/o data
	  // estimate weighted distance
	  float dx = voxNb->stat[kVoxX]-x;
	  float df = voxNb->stat[kVoxF]-p;
	  float dz = voxNb->stat[kVoxZ]-z;
	  float dxw = dx*kWXI, dfw = df*kWFI, dzw = dz*kWZI;
	  double u2 = dxw*dxw+dfw*dfw+dzw*dzw;
	  double kernW = GetKernelWeight(u2);
	  if (kernW<kZeroK) continue;
	  nOccX[ix-ixMn]++;
	  nOccF[ip-ipMn]++;
	  nOccZ[iz-izMn]++;
	  //
	  if (fUseErrInSmoothing) { // apart from the kernel value, account for the point error
	    double kernWX = kernW/(voxNb->E[kResX]*voxNb->E[kResX]);
	    double kernWY = kernW/(voxNb->E[kResY]*voxNb->E[kResY]);
	    double kernWZ = kernW/(voxNb->E[kResZ]*voxNb->E[kResZ]);
	    //
	    m00x += kernWX;
	    m10x += kernWX*dx;   m11x += kernWX*dx*dx;
	    m20x += kernWX*df;   m21x += kernWX*dx*df;  m22x += kernWX*df*df;
	    m30x += kernWX*dz;   m31x += kernWX*dx*dz;  m32x += kernWX*df*dz;   m33x += kernWX*dz*dz;
	    //
	    m00y += kernWY;
	    m10y += kernWY*dx;   m11y += kernWY*dx*dx;
	    m20y += kernWY*df;   m21y += kernWY*dx*df;  m22y += kernWY*df*df;
	    m30y += kernWY*dz;   m31y += kernWY*dx*dz;  m32y += kernWY*df*dz;   m33y += kernWY*dz*dz;
	    //
	    m00z += kernWZ;
	    m10z += kernWZ*dx;   m11z += kernWZ*dx*dx;
	    m20z += kernWZ*df;   m21z += kernWZ*dx*df;  m22z += kernWZ*df*df;
	    m30z += kernWZ*dz;   m31z += kernWZ*dx*dz;  m32z += kernWZ*df*dz;   m33z += kernWZ*dz*dz;
	    //
	    rhsX[0] += kernWX*voxNb->D[kResX];
	    rhsY[0] += kernWY*voxNb->D[kResY];
	    rhsZ[0] += kernWZ*voxNb->D[kResZ];
	    //
	    rhsX[1] += kernWX*voxNb->D[kResX]*dx;
	    rhsY[1] += kernWY*voxNb->D[kResY]*dx;
	    rhsZ[1] += kernWZ*voxNb->D[kResZ]*dx;
	    //
	    rhsX[2] += kernWX*voxNb->D[kResX]*df;
	    rhsY[2] += kernWY*voxNb->D[kResY]*df;
	    rhsZ[2] += kernWZ*voxNb->D[kResZ]*df;
	    //
	    //
	    rhsX[3] += kernWX*voxNb->D[kResX]*dz;
	    rhsY[3] += kernWY*voxNb->D[kResY]*dz;
	    rhsZ[3] += kernWZ*voxNb->D[kResZ]*dz;	    
	  }
	  else { // single matrix is used
	    //
	    m00x += kernW;
	    m10x += kernW*dx;   m11x += kernW*dx*dx;
	    m20x += kernW*df;   m21x += kernW*dx*df;  m22x += kernW*df*df;
	    m30x += kernW*dz;   m31x += kernW*dx*dz;  m32x += kernW*df*dz;   m33x += kernW*dz*dz;
	    //
	    rhsX[0] += kernW*voxNb->D[kResX];
	    rhsY[0] += kernW*voxNb->D[kResY];
	    rhsZ[0] += kernW*voxNb->D[kResZ];
	    //
	    rhsX[1] += kernW*voxNb->D[kResX]*dx;
	    rhsY[1] += kernW*voxNb->D[kResY]*dx;
	    rhsZ[1] += kernW*voxNb->D[kResZ]*dx;
	    //
	    rhsX[2] += kernW*voxNb->D[kResX]*df;
	    rhsY[2] += kernW*voxNb->D[kResY]*df;
	    rhsZ[2] += kernW*voxNb->D[kResZ]*df;
	    //
	    rhsX[3] += kernW*voxNb->D[kResX]*dz;
	    rhsY[3] += kernW*voxNb->D[kResY]*dz;
	    rhsZ[3] += kernW*voxNb->D[kResZ]*dz;
	  }
	  //
	  currClus[nbOK] = voxNb;
	  nbOK++;
	}
      }
    }
  
    // check if we have enough points in every dimension
    int npx=0,npp=0,npz=0;
    for (int i=ixMx-ixMn+1;i--;) if (nOccX[i]) npx++; 
    for (int i=ipMx-ipMn+1;i--;) if (nOccF[i]) npp++;
    for (int i=izMx-izMn+1;i--;) if (nOccZ[i]) npz++;
    if (npx<2 || npp<2 || npz<2 || nbOK<kMinPointsTot) {
      trial++;
      printf("Sector:%2d x=%.3f y/x=%.3f z/x=%.3f (iX:%d iY2X:%d iZ2X:%d)\n"
	     "not enough neighbours (need min %d) %d %d %d (tot: %d) | Steps: %d %d %d\n"
	     "trying to increase filter bandwidth (trial%d)\n",
	     isect,x,p,z,ix0,ip0,iz0,2,npx,npp,npz,nbOK,stepX,stepF,stepZ,trial);
      continue;
    }
    //
    Bool_t fitRes = kTRUE;
    //
    // solve system of linear equations
    AliSymMatrix mat(4);
    if (fUseErrInSmoothing) {
      mat(0,0) = m00x;
      mat(1,0) = m10x;   mat(1,1) = m11x;
      mat(2,0) = m20x;   mat(2,1) = m21x;  mat(2,2) = m22x; 
      mat(3,0) = m30x;   mat(3,1) = m31x;  mat(3,2) = m32x;  mat(3,3) = m33x;
      fitRes &= mat.SolveChol(rhsX);
      //
      mat.Reset();
      mat(0,0) = m00y;
      mat(1,0) = m10y;   mat(1,1) = m11y;
      mat(2,0) = m20y;   mat(2,1) = m21y;  mat(2,2) = m22y; 
      mat(3,0) = m30y;   mat(3,1) = m31y;  mat(3,2) = m32y;  mat(3,3) = m33y;
      fitRes &= mat.SolveChol(rhsY);
      //
      mat.Reset();
      mat(0,0) = m00z;
      mat(1,0) = m10z;   mat(1,1) = m11z;
      mat(2,0) = m20z;   mat(2,1) = m21z;  mat(2,2) = m22z; 
      mat(3,0) = m30z;   mat(3,1) = m31z;  mat(3,2) = m32z;  mat(3,3) = m33z;
      fitRes &= mat.SolveChol(rhsZ);
    }
    else {
      mat(0,0) = m00x;
      mat(1,0) = m10x;   mat(1,1) = m11x;
      mat(2,0) = m20x;   mat(2,1) = m21x;  mat(2,2) = m22x; 
      mat(3,0) = m30x;   mat(3,1) = m31x;  mat(3,2) = m32x;  mat(3,3) = m33x;
      fitRes &= mat.SolveCholN(fLastSmoothingRes,kResDim);
    //
    }
    if (!fitRes) {
      trial++;
      printf("Sector:%2d x=%.3f y/x=%.3f z/x=%.3f (iX:%d iY2X:%d iZ2X:%d)\n"
	     "neighbours range used %d %d %d (tot: %d) | Steps: %d %d %d\n"
	     "Solution for smoothing Failed, trying to increase filter bandwidth (trial%d)\n",
	     isect,x,p,z,ix0,ip0,iz0,npx,npp,npz,nbOK,stepX,stepF,stepZ,trial);
      continue;
    }
    //
    break; // success
  } // end of loop over allowed trials
  res[0] = rhsX[0];
  res[1] = rhsY[0];
  res[2] = rhsZ[0];
  //
  if (deriv) { // derivatives are requested
    for (int i=0;i<kResDim;i++) {
      for (int j=0;j<3;j++) deriv[i*3+j] = fLastSmoothingRes[i*4+j];
    }
  }
  return kTRUE;
}

//_____________________________________
void SetKernelType(int tp, float bwX, float bwP, float bwZ, float scX,float scP,float scZ)
{
  // set kernel type and widths in terms of binning in X,Y/X and Z/X, define aux variables
  fKernelType = tp;
  fStepKern[kVoxX] = TMath::Nint(bwX);
  fStepKern[kVoxF] = TMath::Nint(bwP);
  fStepKern[kVoxZ] = TMath::Nint(bwZ);
  //
  fKernelScaleEdge[kVoxX] = scX;
  fKernelScaleEdge[kVoxF] = scP;
  fKernelScaleEdge[kVoxZ] = scZ;

  if (fKernelType == kEpanechnikovKernel) { // bandwidth 1
  }
  else if (kGaussianKernel) {
    for (int i=0;i<kVoxDim;i++) fStepKern[i] *= 5.; // look in ~5 sigma
  }
  else {
    ::Fatal("GetKernelWeight","Kernel type %d is not defined",fKernelType);
  }
  fNMaxNeighb = 2*(2*fStepKern[kVoxX]+1)*(2*fStepKern[kVoxF]+1)*(2*fStepKern[kVoxZ]+1);
}

//=========================================================================

//_____________________________________________
void AliTPCDistortionExtraction::CreateCorrectionObject()
{
  // create correction object for given time slice

  AliSysInfo::AddStamp("MakeCheb",0,0,0,0);


  TString name = Form("run%d_%lld_%lld",fRun,fTMin,fTMax);
  fChebCorr = new AliTPCChebCorr(name.Data(),name.Data(),
				 fChebPhiSlicePerSector,fChebZSlicePerSide,1.0f);
  fChebCorr->SetUseFloatPrec(kFALSE);
  fChebCorr->SetTimeStampStart(fTMin);
  fChebCorr->SetTimeStampEnd(fTMax);
  fChebCorr->SetTimeDependent(kFALSE);
  fChebCorr->SetUseZ2R(kTRUE);
  //
  fChebCorr->Parameterize(trainCorr,3,fNPCheb,fChebPrecD);
  //
  AliSysInfo::AddStamp("MakeCheb",1,0,0,0);
}

//________________________________________________________________
void AliTPCDistortionExtraction::InitBinning(int nbx, int nby, int nbz)
{
  // initialize binning structures
  //
  // X binning
  if (nbx>0 && nbx<kNPadRows) {
    fNXBins = nbx;
    printf("X-binning: uniform %d bins from %.2f to %.2f\n",fNXBins,kMinX,kMaxX);
    fDXI         = fNXBins/(kMaxX-kMinX);
    fDX          = 1.0f/fDXI;
    fUniformBins[kVoxX] = kTRUE;
  }
  else {
    fNXBins = kNPadRows;
    printf("X-binning: bin per pad-row\n");
    fUniformBins[kVoxX] = kFALSE;
    fDX = kTPCRowDX[0];
    fDXI = 1.f/fDX; // should not be used
  }
  //
  // Y binning
  if (fNXBins<1) ::Fatal("InitYBins","X bins must be initialized first");
  fMaxY2X = new Float_t[fNXBins];        // max Y/X at each X bin, account for dead zones
  fDY2XI  = new Float_t[fNXBins];        // inverse of Y/X bin size at given X bin
  fDY2X   = new Float_t[fNXBins];        // Y/X bin size at given X bin
  //
  fNY2XBins = nby;
  int nxy = fNXBins*fNY2XBins;
  fBinMinQ = new Float_t[nxy];
  fBinDQI  = new Float_t[nxy];
  fBinDQ   = new Float_t[nxy];
  //
  for (int ix=0;ix<fNXBins;ix++) {
    float x = GetX(ix);
    fMaxY2X[ix] = kMaxY2X - kDeadZone/x;
    fDY2XI[ix] = fNY2XBins / (2.f*fMaxY2X[ix]);
    fDY2X[ix] = 1.f/fDY2XI[ix];
    for (int iy=0;iy<fNY2XBins;iy++) {
      float y = GetY2X(ix,iy)*x;
      float tgMn = tgpXY(x,y,-fMaxQ2Pt,fBz);
      float tgMx = tgpXY(x,y, fMaxQ2Pt,fBz);
      if (tgMn>tgMx) swap(tgMn,tgMx);
      int ixy = ix*fNY2XBins + iy;
      fBinMinQ[ixy] = TMath::Abs(fBz)>0.01 ? tgMn : -0.5;
      fBinDQ[ixy]   = TMath::Abs(fBz)>0.01 ? (tgMx-tgMn)/kNQBins : 1.;
      fBinDQI[ixy]  = 1./fBinDQ[ixy];
    }
  }
  //
  fNZ2XBins = nbz;
  fDZ2XI = fNZ2XBins/kMaxZ2X;
  fDZ2X  = 1.0f/fDZ2XI;
  //
  // inverse bin sizes for residuals
  fDeltaYbinI  = fNDeltaYBins/(2.0f*fMaxDY);
  fDeltaZbinI  = fNDeltaZBins/(2.0f*fMaxDZ);
  //
  fNGVoxPerSector = fNY2XBins*fNZ2XBins*fNXBins;
  fNBProdSectG[0] = fNY2XBins*fNZ2XBins;
  fNBProdSectG[1] = fNZ2XBins;

}

//________________________________________________________________
Int_t AliTPCDistortionExtraction::GetXBin(float x) 
{
  // convert X to bin ID, following pad row widths
  if (fUniformBins[kVoxX]) {
    int ix = (x-kMinX)*fDXI;
    if (ix<0) return 0;
    else if (ix>=fNXBins) return fNXBins-1;
    return ix;
  }
  else {
    int ix;
    if (x<kTPCRowX[kNRowIROC-1]+0.5*kTPCRowDX[kNRowIROC-1]) {     // uniform pad size in IROC
      ix = (x-(kTPCRowX[0]-kTPCRowDX[0]*0.5))/kTPCRowDX[0];
      if (ix<0) ix = 0;
    }
    // uniform pad size in OROC2
    else if ( x>= kTPCRowX[kNRowIROC+kNRowOROC1]-0.5*kTPCRowDX[kNRowIROC+kNRowOROC1] ) {
      ix = (x-(kTPCRowX[kNRowIROC+kNRowOROC1]-0.5*kTPCRowDX[kNRowIROC+kNRowOROC1]))/kTPCRowDX[kNPadRows-1]
	+ kNRowIROC + kNRowOROC1;
      if (ix>=kNPadRows) ix = kNPadRows-1;
    }
    else { // uniform pad size in OROC1
      ix = (x-(kTPCRowX[kNRowIROC]-0.5*kTPCRowDX[kNRowIROC]))/kTPCRowDX[kNRowIROC] + kNRowIROC;
      if (ix<kNRowIROC) { // we go between IROC and OROC?
	if (x> 0.5*(kTPCRowX[kNRowIROC-1]+kTPCRowX[kNRowIROC]))  ix = kNRowIROC; // 1st OROC1 row
	else ix = kNRowIROC-1;
      }
      
    }
    return ix;
    // 
  }
}

//________________________________________________________________
Int_t AliTPCDistortionExtraction::GetRowID(float x)
{
  // return row ID
  int ix;
  if (x<kTPCRowX[kNRowIROC-1]+0.5*kTPCRowDX[kNRowIROC-1]) {     // uniform pad size in IROC
    ix = (x-(kTPCRowX[0]-kTPCRowDX[0]*0.5))/kTPCRowDX[0];
    if (ix<0) ix = -1;
  }
  // uniform pad size in OROC2
  else if ( x>= kTPCRowX[kNRowIROC+kNRowOROC1]-0.5*kTPCRowDX[kNRowIROC+kNRowOROC1] ) {
    ix = (x-(kTPCRowX[kNRowIROC+kNRowOROC1]-0.5*kTPCRowDX[kNRowIROC+kNRowOROC1]))/kTPCRowDX[kNPadRows-1]
      + kNRowIROC + kNRowOROC1;
    if (ix>=kNPadRows) ix = -2;
  }
  else { // uniform pad size in OROC1
    ix = (x-(kTPCRowX[kNRowIROC]-0.5*kTPCRowDX[kNRowIROC]))/kTPCRowDX[kNRowIROC] + kNRowIROC;
    if (ix<kNRowIROC) { // we go between IROC and OROC?
      ix = -3;
    }
    
  }
  return ix;
}

//_____________________________________________________
Bool_t AliTPCDistortionExtraction::FindVoxelBin(int sectID, float q2pt, float x, float y, float z, 
						UChar_t bin[kVoxHDim],float voxVars[kVoxHDim])
{
  // define voxel variables and bin
  //  
  //
  // track Z/X bin
  voxVars[kVoxZ] = z/x;
  if (TMath::Abs(voxVars[kVoxZ])>kMaxZ2X) return kFALSE;
  int binZ       = GetZ2XBinExact( sectID<kNSect ? voxVars[kVoxZ] : -voxVars[kVoxZ] );
  if (binZ<0) return kFALSE;
  bin[kVoxZ] = binZ;
  //
  // track X bin
  voxVars[kVoxX] = x;
  int binX = GetXBinExact(x); // exact matching to pad-rows
  if (binX<0 || binX>=fNXBins) return kFALSE;
  bin[kVoxX] = binX;
  //
  // track Y/X bin accounting for sector edges dead zones
  voxVars[kVoxF] = y/x;
  int binF = GetY2XBinExact(voxVars[kVoxF],binX);
  if (binF<0||binF>=fNY2XBins) return kFALSE;
  bin[kVoxF] = binF;
  //
  // track inclination at pad row
  voxVars[kVoxQ] = tgpXY(x,y,q2pt,fBz);
  int binQ = GetQBin(q2pt); //GetQBin(voxVars[kVoxQ],binX,binF);
  if (binQ<0) return kFALSE;
  bin[kVoxQ] = binQ;
  //
  return kTRUE;
}

//=====================================================================
///////////////////////////////////////////////////
//
// This is a special non-meber f-n for cheb. learning
//
void trainCorr(int row, float* tzLoc, float* corrLoc)
{
  // Cheb. object training f-n: compute correction for the point
  //
  // xtzLoc: y2x and z2x in sector frame
  // corrLoc: x, y, z corrections in sector frame
  // when called with pointer at 0, row will set the sector/side 
  // (row should be sector in 0-35 or 0-71 format)
  static int sector=0;
  if (!tzLoc || !corrLoc) {
    sector = row%36; 
    printf("training Sector %d\n",sector);
    return;
  }
  //
  float x = kTPCRowX[row];
  float dist[3], deriv[9];

  float y2x = tzLoc[0];
  float z2x = tzLoc[1];
  //
  Bool_t res = GetSmoothEstimate(sector, x, y2x, z2x, dist);
  if (!res) { printf("Failed to evaluate smooth distortion\n"); exit(1); }

  // Marian stored Z track coordinate instead of cluster one, need to correct for this
  if (fApplyZt2Zc) {
    const double inversionEps = 20e-4; // when inverting, stop Newton-Raphson iterations at this eps
    const int    inversionMaxIt = 3; // when inverting, stop Newton-Raphson after some numbers of iterations
    double change = 0, xInv = 1./x;
    float zc = z2x*x;
    float zt = zc + dist[kResZ]; // 1st guess on true zt at measured zc
    //
    // use Newton-Raphson method for NDLocal inversion to get zt = F(zc) from 
    // dz == zt - zc = NDLoc(zt)   ->  zc - [zt - NDLoc(zt)]=0 
    // ->  zt_{i+1} = zt_i - (zc - [zt - NDLoct(zt_i)])/(-d[zt_i-NDLoc(zt_i)]/dz)
    //
    int it = 0;
    do {
      z2x = zt*xInv;
      res = GetSmoothEstimate(sector, x, y2x, z2x, dist, deriv);
      if (!res) {printf("Failed to evaluate smooth distortion\n");exit(1);}
      double bot = 1. - deriv[kResZ*3+2]*xInv;  // dF(zt_i)/dz
      if (TMath::Abs(bot)<1e-6) break;
      double top = zc - (zt - dist[kResZ]); // zc - F(zt_i) 
      double change = top/bot;
      //  printf("It %d Eps:%+e, Zc:%+e, Zt:%+e Ztn:%+e DZ:%+e | dztmp: :%+e DD:%+e\n",
      //     it,change,zc,zt,zt+change,dz,dztmp,deriv[kndZ2R]*xInv);
      zt += change;
      it++;
    } while(it<inversionMaxIt && TMath::Abs(change)>inversionEps);
    //
    // now query at fixed Z2X
    z2x = zt*xInv;
    res = GetSmoothEstimate(sector, x, y2x, z2x, dist);
    if (!res) {printf("Failed to evaluate smooth distortion\n");exit(1);}
  }
  corrLoc[kResX] = dist[kResX];
  corrLoc[kResY] = dist[kResY];
  corrLoc[kResZ] = dist[kResZ];
  //
}
//======================================================================================
