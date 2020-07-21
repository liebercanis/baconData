{
  printf("\n this is rootlogon for bacon \n");
  TString arch=gSystem->GetBuildArch();
  cout << " arch is " << arch << endl; 
  gSystem->AddIncludePath(" -I. -I./bobj/");
  gSystem->AddDynamicPath("./bobj/");
  printf(" include path %s \n\n",gSystem->GetIncludePath());
  cout << "DYNAMIC PATH "  << gSystem->GetDynamicPath() << endl;
  cout << "LINKED LIBS "  << gSystem->GetLinkedLibs() << endl;
  gROOT->LoadMacro("util.C");
  int iload = gSystem->Load("bobj/libBaconData.so");
  printf(" loaded libBaconData = %i zero is success! \n",iload);
}

