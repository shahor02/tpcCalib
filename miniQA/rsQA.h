#ifndef RSQA_H
#define RSQA_H

#include "AliAnalysisTaskSE.h"
#include <TString.h>
#include <TStopwatch.h>
#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "HistoManager.h"
class AliAlgSteer;
class TTree;
class TFile;
class AliStack;
class AliMCEvent;


class rsQA : public AliAnalysisTaskSE {
 public:
  //
  enum {kQPlus,kQMinus,kNCharges};
  enum {kSideA,kSideC,kNSides};
  enum {kHDCAY, kHDCAZ};
  enum {kHTrcTPCSect,kHTrcTPCITSSect,kHTrcTPCITSSPDSect};
  enum {kHTrcTPCEta,kHTrcTPCITSEta,kHTrcTPCITSSPDEta};  
  enum {kNBinsPerSect=10, kNBinsEta=20, kNBinsQ2Pt=40, kNBinsDCAY=120, kNBinsDCAZ=80, kNPtRanges=3, kNEtaRanges=10};

  rsQA(const char *name = "rsQA");
  virtual ~rsQA(); 

  int GetHistoID(int side, int kH) {return side*10+kH;}

  virtual void  UserCreateOutputObjects();
  virtual void  UserExec(Option_t *option);
  virtual void  Terminate(Option_t *);
  void ProcessEvent(AliESDEvent* ev);
  Bool_t ProcessTrack(const AliESDtrack* trc, const AliESDVertex* vtx);
  //
  void SetUseMC(Bool_t v) {fUseMC=v;}
  Int_t PdgToPid(int pdg) const;
  //
 protected:
  //
  TList*    fOutput;                   // output list send on output slot 1
  TTree* fTreeOut;
  TFile* fFileOut;
  //
  TObjArray fHistManArr;
  TString fOutFName;
  //
  HistoManager 
    *fHManDCA[kNCharges][kNPtRanges], 
    *fHManMtcSect[kNSides], 
    *fHManMtcEta;
  //
  static const Double_t kMaxZVtx,kMaxDCAY,kMaxDCAZ,kMaxQ2Pt,kMaxEta,kPtRange[kNPtRanges];
  Double_t fBz;
  //
  // unused stuff
  Bool_t fUseMC;
  AliMCEvent*  fMCEvent; //!
  AliStack*    fStack;   //!

 private:    
  rsQA(const rsQA&); // not implemented
  rsQA& operator=(const rsQA&); // not implemented 
  //  
  ClassDef(rsQA, 1);  
};


#endif
