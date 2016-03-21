
{
  //  Bool_t ymode = kFALSE;//
  Bool_t ymode = kTRUE;

  if (ymode) SetModeY();
  else       SetModeZ();


  //  TString destName = "figProp/lhc11h_168511_vs_245231_fixF22_M16_EP";
  //TString destName = "figProp/lhc11h_168511_vs_244918_fixF22_M16_EP";
  TString destName = "figProp/lhc15o_244918_M14_M14Ref";

  TCanvas* cnv[3][3], *ctmp;

  //ProcMiniQA("out/coll_LHC11h_168511Std.root","r168511Std",kGreen,20,0.6);
  //  ProcMiniQA("out/coll_245231_Feb22Std.root","r245231std",kGreen,20,0.6);
  //  ProcMiniQA("out/coll_LHC15o_244918Std.root","r244918std",kGreen,20,0.6);
  ProcMiniQA("out/coll_245231_Feb22bfix.root","r245231corr",kGreen,24,0.6);

  for (int iq=0;iq<=2;iq++) {
    for (int ip=0;ip<=2;ip++) {
      ctmp = DrawDCAPhi(iq,ip,0);
      cnv[iq][ip] = ctmp;
      cnv[iq][ip]->cd(2);
      //AddLabel("LHC11h_168511", 0.15, 0.96, kBlack, 0.05);
      //      AddLabel("LHC15o 245231", 0.15, 0.96, kBlack, 0.05);
      AddLabel("LHC15o 245231", 0.15, 0.96, kBlack, 0.05);
      AddLabel("Corr", 0.38, 0.96, kBlack, 0.05);
      AddLabel("Feb22", 0.45, 0.96, kGreen, 0.05);
      //AddLabel("NoRef", 0.4, 0.96, kBlue, 0.05);
      cnv[iq][ip]->cd();
    }
  }
  //
  
  //  ProcMiniQA("out/coll_245231_Feb22bfix.root","r245231corr",kBlue,24,0.6);
  //  ProcMiniQA("out/coll_244918_Feb22bfix.root","r244918corrF22",kBlue,24,0.6);
  ProcMiniQA("/data/testTPC/miniQA/out/coll_245231_March14fix.root","r244918corrM14",kBlue,20,0.6);

  for (int iq=0;iq<=2;iq++) {
    for (int ip=0;ip<=2;ip++) {
      cnv[iq][ip]->cd();
      DrawDCAPhi(iq,ip,cnv[iq][ip]);
      cnv[iq][ip]->cd(2);
      //      AddLabel("LHC15o 244918", 0.15, 0.96, kBlack, 0.05);
      //AddLabel("Corr", 0.46, 0.96, kBlack, 0.05);      
      AddLabel("14/03/16", 0.55, 0.96, kBlack, 0.05);      
      AddLabel("w/o Ref", 0.7, 0.96, kBlue, 0.05);      

      cnv[iq][ip]->cd();
      //      SaveCanvas(cnv[iq][ip],Form("%s_%s_q%d_pt%d",destName.Data(),ymode ? "Y":"Z",iq,ip),"cg");
    }
  }

  //  ProcMiniQA("/data/testTPC/miniQA/out/coll_245231_March14fix.root","r244918corrM14",kBlue,20,0.6);

  //  ProcMiniQA("out/coll_244918_EP_NOref_March16.root","r244918corrM16",kRed,20,0.6);

  ProcMiniQA("out/coll_March14_WithRef_000245231.root","r244918corrM14Ref",kRed,20,0.6);

  //  ProcMiniQA("/data/testTPC/miniQA/out/coll_245231_March14fix.root","r244918corrM14",kRed,20,0.6);
  for (int iq=0;iq<=2;iq++) {
    for (int ip=0;ip<=2;ip++) {
      cnv[iq][ip]->cd();
      DrawDCAPhi(iq,ip,cnv[iq][ip]);
      cnv[iq][ip]->cd(2);
      //      AddLabel("Corr", 0.46, 0.96, kBlack, 0.05);      
      // AddLabel("16/03/16", 0.7, 0.96, kRed, 0.05);           
      AddLabel("w/Ref03", 0.85, 0.96, kRed, 0.05);      
      //      AddLabel("245231", 0.85, 0.96, kRed, 0.05);           
      cnv[iq][ip]->cd();
      SaveCanvas(cnv[iq][ip],Form("%s_%s_q%d_pt%d",destName.Data(),ymode ? "Y":"Z",iq,ip),"cg");
    }
  }

}
