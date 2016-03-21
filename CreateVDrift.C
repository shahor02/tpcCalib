
void CreateVDrift(int run
		  ,const char* ocdbDest="local://"
		  ,const char* fitDrift="fitDrift.root"
		  ,const char* ocdbSrc="local:///cvmfs/alice.cern.ch/calibration/data/2015/OCDB")
{
  AliCDBManager * man = AliCDBManager::Instance();
  man->SetDefaultStorage(ocdbSrc ? ocdbSrc : "raw://");
  man->SetRun(run);
  AliGRPManager gm;
  gm.ReadGRPEntry();
  gm.SetMagField();
  AliTPCcalibAlignInterpolation::MakeVDriftOCDB(fitDrift,run,ocdbDest,0);
}
