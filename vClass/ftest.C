#include "TMath.h"
#include "TVectorF.h"
#include "TStatToolkit.h"
#include "TGraphErrors.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TPad.h"

TGraphErrors* grMLTM=0;
TGraphErrors* grSLTM=0;
TGraphErrors* grMG=0;
TGraphErrors* grSG=0;
TGraph* grLogL=0;
TGraph* grLogL0=0;
TGraph* grLogLRat=0;
TH1F *hbaseM=0,*hbaseS=0,*hbaseL=0;
TCanvas* cnv=0;
TH1F* histoClone = 0;

Double_t GetLogL(TH1F* histo, int bin0, int bin1, double &mu, double &sig, double &logL0);
void TruncNormMod(double a, double b, double mu0, double sig0, double &muCf, double &sigCf);
Bool_t GetTruncNormMuSig(double a, double b, double &mean, double &sig);


//______________________________________________________________________________
void ftest(TH1F* histo, int kMinN=10, float maxLTM=1.0, float minLTM=0.5, float dLTM=0.05, float sclSLTM=1.5)
{
  static TF1 fgaus("fgaus","gaus",-10,10);
  static float est[20];
  //
  memset(est,0,20*sizeof(float));
  //
  if (!cnv) {
    cnv = new TCanvas("cnv","cnv",1200,900);
    cnv->Draw();
  }
  cnv->Clear();
  cnv->Divide(3,2);
  cnv->cd(1);
  delete histoClone;
  histoClone = (TH1F*) histo->Clone("histoClone");
  histoClone->Draw();
  gPad->SetGrid();
  gPad->cd(2);
  //
  int nb = histo->GetNbinsX();

  int nLTM = (maxLTM-minLTM)/dLTM;
  if (nLTM<1) nLTM=1;
  dLTM = (maxLTM-minLTM)/nLTM;

  delete grMLTM; grMLTM = new TGraphErrors(nLTM);
  delete grSLTM; grSLTM = new TGraphErrors(nLTM);
  delete grMG; grMG = new TGraphErrors(nLTM);
  delete grSG; grSG = new TGraphErrors(nLTM);
  delete grLogL;  grLogL  = new TGraph(nLTM); 
  delete grLogL0; grLogL0 = new TGraph(nLTM); 
  delete grLogLRat; grLogLRat = new TGraph(nLTM); 

  TVectorF vecLTM(10);
  int nfit = 0;
  float minM=1e9,maxM=-1e9, minS=1e9,maxS=-1e9, minL=1e99,maxL=-1e99;
  Bool_t drawDone = kFALSE;

  int bminPrev = -1, bmaxPrev = -1;
  float dxh = histo->GetBinWidth(1)*0.5;
  //
  double prevLRat = 1e6;
  double maxVal = histo->GetMaximum();
  for (int iltm=0;iltm<nLTM;iltm++) {
    float cut = maxLTM - iltm*dLTM;
    TStatToolkit::LTMHisto(histo, vecLTM, cut);
    int bmin = int(vecLTM[5]), bmax = int(vecLTM[6]);
    printf("Test at %f %d %d\n",cut, bmin,bmax);
    //
    if (bmin==bminPrev && bmax==bmaxPrev) {
      printf("Same bin sel at %f %d %d\n",cut, bmin,bmax);
      //      vecLTM.Print();
      continue;  // don't check twice same truncation
    }
    bminPrev = bmin;
    bmaxPrev = bmax;    
    // skip empty bins from edges
    while (!histo->GetBinContent(bmin)) bmin++;
    while (!histo->GetBinContent(bmax)) bmax--;
    double muEst=0,sigEst=1;
    //
    if (bmax==bmin) { // do not allow single bins ?
      if (bmin>1) bmin--;
      if (bmax<nb) bmax++; 
    }
    // get estimate of real mean and sigma from truncated mean and corresponding
    // sample and reference log-likelihood
    double logL0=0, logL = GetLogL(histo,bmin,bmax,muEst,sigEst,logL0);
    //
    double logLRat = (logL-logL0)/logL0;
    grLogL->SetPoint(nfit, 1.f-cut, logL);
    grLogL0->SetPoint(nfit, 1.f-cut, logL0);
    grLogLRat->SetPoint(nfit,1.f-cut,logLRat); //(logL+logL0));
    //
    fgaus.SetParameters(maxVal ,muEst, sigEst);
    
    float rngMin = histo->GetBinCenter(bmin) - dxh;
    float rngMax = histo->GetBinCenter(bmax) + dxh;
    //
    TFitResultPtr fitPtr= histo->Fit(&fgaus,"qnrLS","",rngMin,rngMax); //maxVal<kUseLLFrom ? "qnrlS":"qnrS");
    TFitResult * result = fitPtr.Get();
    float estMG=0,estSG=0,estMGE=0,estSGE=0,chi2=0;
    if (result!=NULL) {
      estMG = fgaus.GetParameter(1);
      estSG  = fgaus.GetParameter(2);
      estMGE = fgaus.GetParError(1);
      estSGE  = fgaus.GetParError(2);
      chi2 = fgaus.GetChisquare()/fgaus.GetNumberFreeParameters();
      printf("#%2d          Mean: %+6.3f RMS: %.3f MErr: %+.3f SErr: %.3f | Chi2: %.2f |in %+6.3f %+6.3f\n", 
	     nfit,estMG,estSG,estMGE,estSGE,chi2,rngMin,rngMax);
      
      grMG->SetPoint(nfit, 1.f-cut, estMG);
      grSG->SetPoint(nfit, 1.f-cut, estSG);
      grMG->SetPointError(nfit, dLTM/2., TMath::Max(1e-6f,estMGE));
      grSG->SetPointError(nfit, dLTM/2., TMath::Max(1e-6f,estSGE));

      if (minM>estMG) minM = estMG;
      if (minS>estSG) minS = estSG;      
      if (maxM<estMG) maxM = estMG;
      if (maxS<estSG) maxS = estSG;
      //
      if (!drawDone && abs(logLRat)<0.5 && abs(logLRat-prevLRat)<0.2 && abs(logL-logL0)<5.) { //RRRR
	cnv->cd(2);
	fgaus.SetParameters(maxVal ,muEst, sigEst);
	printf("Fitting in %e %e IniParams: %e %e %e\n",rngMin,rngMax,maxVal ,muEst, sigEst);
	histo->Fit(&fgaus,"L","",rngMin,rngMax);
	histo->GetXaxis()->SetRangeUser(estMG-3*estSG,estMG+3*estSG);
	gPad->SetGrid();
	gPad->Modified();
	drawDone = kTRUE;
      }
      prevLRat = logLRat;
      //
    }

    float estMLTM = vecLTM[1];
    float estSLTM = vecLTM[2];
    float estMLTME = vecLTM[3];
    float estSLTME = vecLTM[4];
    int nb0 = vecLTM[5];
    int nb1 = vecLTM[6];
    int nbd = nb1-nb0+1;
    
    grMLTM->SetPoint(nfit, 1.f-cut, estMLTM);
    grSLTM->SetPoint(nfit, 1.f-cut, estSLTM);
    grMLTM->SetPointError(nfit, dLTM/2., TMath::Max(1e-6f,estMLTME));
    grSLTM->SetPointError(nfit, dLTM/2., TMath::Max(1e-6f,estSLTME));
    
    if (minM>estMLTM) minM = estMLTM;
    if (minS>estSLTM) minS = estSLTM;

    if (maxM<estMLTM) maxM = estMLTM;
    if (maxS<estSLTM) maxS = estSLTM;

    if (minL>logL0) minL=logL0;
    if (minL>logL)  minL=logL;
    if (maxL<logL0) maxL=logL0;
    if (maxL<logL)  maxL=logL;
    //
    printf("#%2d Cut:%4.2f Mean: %+6.3f RMS: %.3f MErr: %+.3f SErr: %.3f | %3d %3d -> %3d |LL %+e vs %+e\n",
	   nfit,cut,estMLTM,estSLTM,estMLTME,estSLTME,nb0,nb1,nb1-nb0,logL,logL0);
    nfit++;
  }
  //
  printf("%f-%f  %f-%f\n",minM,maxM,minS,maxS);
  delete hbaseM;
  delete hbaseS;
  //
  hbaseM = new TH1F("hbaseM","M",2*nLTM,1.f-maxLTM-dLTM,1.f-minLTM+dLTM);
  hbaseS = new TH1F("hbaseS","S",2*nLTM,1.f-maxLTM-dLTM,1.f-minLTM+dLTM);
  hbaseL = new TH1F("hbaseL","L",2*nLTM,1.f-maxLTM-dLTM,1.f-minLTM+dLTM);
  //
  hbaseM->SetMinimum( minM - 0.15*(maxM-minM) );
  hbaseM->SetMaximum( maxM + 0.15*(maxM-minM) );
  hbaseS->SetMinimum( minS - 0.15*(maxS-minS) );
  hbaseS->SetMaximum( maxS + 0.15*(maxS-minS) );
  hbaseL->SetMinimum( minL - 0.15*(maxL-minL) );
  hbaseL->SetMaximum( maxL + 0.15*(maxL-minL) );
  //
  cnv->cd(3);
  hbaseM->Draw();
  grMLTM->SetMarkerStyle(24);
  grMLTM->SetMarkerColor(kBlue);
  grMLTM->SetLineColor(kBlue);
  grMLTM->Draw("p");
  grMG->SetMarkerStyle(20);
  grMG->SetMarkerColor(kRed);
  grMG->SetLineColor(kRed);
  grMG->Draw("p");
  gPad->SetGrid();
  //
  cnv->cd(4);
  hbaseS->Draw();
  grSLTM->SetMarkerStyle(24);
  grSLTM->SetMarkerColor(kBlue);
  grSLTM->SetLineColor(kBlue);
  grSLTM->Draw("p");
  grSG->SetMarkerStyle(20);
  grSG->SetMarkerColor(kRed);
  grSG->SetLineColor(kRed);
  grSG->Draw("p");
  gPad->SetGrid();
  //
  cnv->cd(5);
  hbaseL->Draw();
  grLogL0->SetMarkerStyle(24);
  grLogL0->SetMarkerColor(kBlue);
  grLogL0->SetLineColor(kBlue);
  grLogL0->Draw("p");
  grLogL->SetMarkerStyle(20);
  grLogL->SetMarkerColor(kRed);
  grLogL->SetLineColor(kRed);
  grLogL->Draw("p");
  gPad->SetGrid();
  //
  cnv->cd(6);
  grLogLRat->SetMarkerStyle(20);
  grLogLRat->SetMarkerColor(kRed);
  grLogLRat->SetLineColor(kRed);
  grLogLRat->Draw("ap");
  gPad->SetGrid();
}


Double_t GetLogL(TH1F* histo, int bin0, int bin1, double &mu, double &sig, double &logL0)
{
  // Calculate log likelihood of normal distribution for the histo between boundaries 
  // bin0 and bin for given mu and sigma assumption. Exact Poisson statistics is assumed
  // Also the approximate "reference" log-likelihood logLO is calculated in the following way:
  // if the Poisson prob. for given bin is "m", then the logL0 gets contribution for this bin
  // log("reference probability"), which I define as a geometric mean of probabilities for
  // observing [M] and [M]+1 entries if m is large: 
  // P_ref = exp(-m)/m! M^m sqrt( m/(M+1) )
  // -> ln(P_ref) = -(1/2+M)*ln(M/m) + M-m + 1/2 ln(M+1) - 1/2 ln(2 pi) -> -1/2 ln(2 pi m) for m>>1 (take m>5)
  //                (precise up to 1/2 ln(2 pi m) term of Stirling formula
  // or           = -m + m*log(m) - log( Gamma(1.+m) )                                     for m<~1
  // 
  // integral
  const double kNuLarge = 5.0, kMinSig2BinH = 0.01;
  double dxh = 0.5*histo->GetBinWidth(1);
  if ((sig/dxh)<kMinSig2BinH) {
    printf("Too small sigma %.4e is provided for bin width %.4e\n",sig,dxh);
    logL0 = -1;
    return -1e9;
  }
  double sum=0, sum1=0, sum2=0;
  
  for (int ib=bin0;ib<=bin1;ib++) {
    double w = histo->GetBinContent(ib);
    double x = histo->GetBinCenter(ib);
    sum += w;
    sum1 += w*x;
    sum2 += w*x*x;
  }  
  //
  double xb0 = histo->GetBinCenter(bin0)-dxh;
  double xb1 = histo->GetBinCenter(bin1)+dxh;
  //
  if (sum<1e-6) {logL0 = -1e6; return -1e9;}
  mu = sum1/sum;
  sig = TMath::Sqrt(sum2/sum - mu*mu);
  //printf("Sample mu : %e sig: %e in %e %e\n",mu,sig,xb0,xb1);

  // estimated sig, mu are from the truncated sample, try to recover the truth
  GetTruncNormMuSig(xb0,xb1, mu, sig);
  //
  xb0 -= mu;
  xb1 -= mu;
  double sqri2 = 1./(TMath::Sqrt(2.)*sig);
  // normalization constant
  double norm = 2.*sum / (TMath::Erf(xb1*sqri2) - TMath::Erf(xb0*sqri2));
  //
  //  printf("Norm: %e\n",norm);
  // likelihood
  double logL = 0;
  logL0 = 0;
  const double kMinExp = 1e-100;
  for (int i=bin0;i<=bin1;i++) {
    double x = histo->GetBinCenter(i)-mu;
    double w = histo->GetBinContent(i);
    xb0 = x-dxh;
    xb1 = x+dxh;
    // bin expectation: normal integral within the bin
    double nu = 0.5*norm*(TMath::Erf(xb1*sqri2) - TMath::Erf(xb0*sqri2));  
    if (nu<kMinExp) nu = kMinExp;
    double logNFac = w<100 ? TMath::Log(TMath::Factorial(w)) : w*TMath::Log(w)-w + TMath::Log( sqrt(2*TMath::Pi()*w));
    double logNu = TMath::Log(nu);
    double logc = -nu + w*logNu - logNFac;  // contribution of this bin to log-likelihood
    logL += logc;
    // now get the reference contribution
    double logc0 = 0;
    if (nu>kNuLarge) logc0 = -0.5*TMath::Log(2.*TMath::Pi()*nu);
    else {
      logc0 = -nu + nu*logNu - TMath::Log( TMath::Gamma(1.+nu) );
    }
    logL0 += logc0;  // reference LL update
    //printf("b: %d x:%+.2e nstd:%+.2e Exp:%e Obs:%e logc: %e logc0: %e\n",i,x,(x-mu)/sig, nu,w,logc, logc0);

  }
  //  printf("LogL: %e LogL0: %e\n",logL,logL0);
  //
  return logL;
}

//___________________________________________________________________________
void TruncNormMod(double a, double b, double mu0, double sig0, double &muCf, double &sigCf)
{
  // calculate truncated mean and sigma of normal distribution as 
  // mu_tr  = mu0 + sig0*muCf
  // sig_tr = sig0 * sigCf
  //
  const double sqrt2PiI = 1./TMath::Sqrt(TMath::Pi()*2.);
  double sigI = 1./(sig0*TMath::Sqrt(2.));
  double ra = (a-mu0)*sigI, rb = (b-mu0)*sigI;
  double ra2 = ra*ra, rb2 = rb*rb;
  double af = ra2<100 ? sqrt2PiI*TMath::Exp(-ra2) : 0;
  double bf = rb2<100 ? sqrt2PiI*TMath::Exp(-rb2) : 0;
  //  double aF = 0.5*(1.+TMath::Erf(ra)), bF = 0.5*(1.+TMath::Erf(rb)), deltaF = bF-aF
  double deltaF = 0.5*( TMath::Erf(rb) - TMath::Erf(ra) );
  double deltaf = af - bf;
  muCf = deltaf / deltaF;
  sigCf = 1./TMath::Sqrt(1. + TMath::Sqrt(2)*(ra*af-rb*bf)/deltaF - muCf*muCf); 
  //
}

//_____________________________________________________
Bool_t GetTruncNormMuSig(double a, double b, double &mean, double &sig)
{
  // get estimate of real mu and sigma of normal distribution provided
  // the mean and rms of sample truncated between a and b
  const double kMinWindow=1e-2,kEpsRMS = 1e-4, kEpsMu = 1e-4;
  const int kMaxIter = 200;
  //
  if (sig<1e-12) {
    printf("Input sigma %e is too small\n",sig);
    return kFALSE;
  }
  if ( (b-a)/sig<kMinWindow ) {
    printf("Truncation window %e-%e is %e sigma only\n",a,b,(b-a)/sig);
    return kFALSE;
  }
  //
  double sig0=sig, mean0=mean; // initial values
  // for protection, don't allow the sigma to grow above twice of the flat distribution
  double sigMax = 2*(b-a)/TMath::Sqrt(12.);
  //
  double m = mean, s = sig;
  for (int i=0;i<kMaxIter;i++) {
    double sclRMS,sclMU;
    TruncNormMod(a,b,m,s, sclMU,sclRMS);
    //
    s = sig * sclRMS;
    double mPrev = m, sPrev = s;
    m = mean - sclMU*s;
    //printf("%d -> M: %e S: %e\n",i,m, s);
    if ( s>sigMax) {
      printf("Iteration took sigma to twice of the flat distribution for "
	     "mu0=%+.3e sig0=%.3e in %+.3e:%+.3e interval\n",mean0,sig0, a,b);
      printf("Abandoning and returning input sigma and mean\n");
      m = mean0; s = sig0; break;
    }
    if (TMath::Abs(1.-sPrev/s)<kEpsRMS && TMath::Abs(m-mPrev)<kEpsMu ) break;
  }
  //
  mean = m;
  sig  = s;

  return kTRUE;
}
