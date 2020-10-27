////////////////////////////////////////////////////////
enum {NDET=4};
TString tag;
TFile *fin;
TDirectory *evdir;
TCanvas *can;

void openFile() {
  // open ouput file and make some histograms
  TString fileName; fileName.Form("ana-%s.root",tag.Data());
  printf(" looking for file %s\n",fileName.Data());
  fin = new TFile(fileName,"readonly");
  if(fin->IsZombie()) {
    printf(" couldnt open file %s\n",fileName.Data());
    return;
  }
  fin->GetObject("EvDir",evdir);
  evdir->ls();
}


void base(Int_t ievent=0, Int_t idet=0) {

  TH1D* hist=NULL;
  TH1D* bhist=NULL;

  if (!fin->IsOpen()) return;

  TString search;
  search.Form("EvWave%1i_DET%1i",ievent,idet);
  TString bsearch;
  bsearch.Form("EvBase%1i_DET%1i",ievent,idet);
  printf(" looking for %s %s \n",search.Data(),bsearch.Data());

  
  TList* list = evdir->GetListOfKeys() ;
  if (!list) { printf("<E> No keys found in file\n") ; exit(1) ; }
  TIter next(list);
  TKey* key;
  TObject* obj;

  while ( (key = (TKey*)next()) ) {
    obj = key->ReadObj() ;
    if (!obj->InheritsFrom("TH1D")) continue;
    TString hname(obj->GetName());
    if( hname.Contains(search)  ) hist = (TH1D*) obj;
    if( hname.Contains(bsearch)) bhist = (TH1D*) obj;
  }

  if(!(hist&&bhist)) { 
    printf(" cannot find hist for event %i \n",ievent);
    return;
  }

  printf("found %s  %s \n",hist->GetName(),  bhist->GetName() );

  hist->SetLineColor(kBlack);
  hist->SetFillStyle(0);
  

  bhist->SetLineColor(kGreen);
  bhist->SetLineWidth(2);



  TString canName;
  canName.Form("%s-Base-Event%iDet%i",tag.Data(),ievent,idet);

  TCanvas * can = new TCanvas(canName,canName);
  hist->Draw("");
  bhist->Draw("same");
  can->Print(".png");
}


void plote(int ievent=0) {

  TH1D* hist[NDET];
  TH1D* phist[NDET];
  TH1D* dhist[NDET];

  if (!fin->IsOpen()) return;

  TString search[NDET];
  TString psearch[NDET];
  TString dsearch[NDET];
  for(int id=0; id<NDET; ++id) {
    hist[id]=NULL;
    search[id].Form("EvSignal%1i_DET%1i",ievent,id);
    phist[id]=NULL;
    psearch[id].Form("EvPeaks%1i_DET%1i",ievent,id);
    dhist[id]=NULL;
    dsearch[id].Form("EvDWave%1i_DET%1i",ievent,id);
  }

  
  TList* list = evdir->GetListOfKeys() ;
  if (!list) { printf("<E> No keys found in file\n") ; exit(1) ; }
  TIter next(list);
  TKey* key;
  TObject* obj;

  while ( (key = (TKey*)next()) ) {
    obj = key->ReadObj() ;
    if (!obj->InheritsFrom("TH1D")) continue;
    TString hname(obj->GetName());
     for(int id=0; id<NDET; ++id) {
       if( hname.Contains(search[id]) ) hist[id] = (TH1D*) obj;;
       if( hname.Contains(psearch[id]) ) phist[id] = (TH1D*) obj;;
       if( hname.Contains(dsearch[id]) ) dhist[id] = (TH1D*) obj;;
     }
  }

  for(int id=0; id<NDET; ++id) {
    if(!(hist[id])) { 
      printf(" cannot find hist %s for event %i \n",search[id].Data(),ievent);
      return;
    }
    if(!(phist[id])) { 
      printf(" cannot find hist %s for event %i \n",psearch[id].Data(),ievent);
      return;
    }
    if(!(dhist[id])) { 
      printf(" cannot find hist %s for event %i \n",dsearch[id].Data(),ievent);
      return;
    }
  }

  TString canName;
  canName.Form("%s-Raw-All-Event%i",tag.Data(),ievent);
  TCanvas * can = new TCanvas(canName,canName);
  can->Divide(1,4);
  for(int id=0; id<NDET; ++id) {
    can->cd(id+1);
    hist[id]->Draw();
    phist[id]->SetLineColor(kRed);
    phist[id]->SetFillColor(kRed);
    phist[id]->SetFillStyle(3001);
    phist[id]->Draw("sames");
  }
  can->Print(".pdf");

  return;

  canName.Form("%s-DWave-All-Event%i",tag.Data(),ievent);
  TCanvas * dcan = new TCanvas(canName,canName);
  dcan->Divide(1,4);
  for(int id=0; id<NDET; ++id) {
    dcan->cd(id+1);
    dhist[id]->Draw();
  }


  canName.Form("%s-Peaks-All-Event%i",tag.Data(),ievent);
  TCanvas * pcan = new TCanvas(canName,canName);
  pcan->Divide(1,4);
  for(int id=0; id<NDET; ++id) {
    pcan->cd(id+1);
    phist[id]->Draw();
  }
  //can->Print(".pdf");

}
/** look at fft ***/
void plotfft(Int_t ievent=0, Int_t idet=0) {

  TH1D* hist=NULL;
  TH1D* fhist=NULL;

  if (!fin->IsOpen()) return;

  TString search;
  search.Form("EvRawWave%1i_DET%1i",ievent,idet);
  TString fsearch;
  fsearch.Form("EvInvFFT%1i_DET%1i",ievent,idet);
  printf(" looking for %s %s \n",search.Data(),fsearch.Data());


  
  TList* list = evdir->GetListOfKeys() ;
  if (!list) { printf("<E> No keys found in file\n") ; exit(1) ; }
  TIter next(list);
  TKey* key;
  TObject* obj;

  while ( (key = (TKey*)next()) ) {
    obj = key->ReadObj() ;
    if (!obj->InheritsFrom("TH1D")) continue;
    TString hname(obj->GetName());
    if( hname.Contains(search)  ) hist = (TH1D*) obj;
    if( hname.Contains(fsearch)) fhist = (TH1D*) obj;
  }

  if(!(hist&&fhist)) { 
    printf(" cannot find hist for event %i \n",ievent);
    return;
  }

  printf("found %s  %s \n",hist->GetName(), fhist->GetName() );

  TString histName;
  histName.Form("run-%s-det%i-Event%i",tag.Data(),idet,ievent);
  hist->SetTitle(histName);

  hist->SetLineColor(kBlack);
  fhist->SetLineColor(kRed);


  TString canName;
  canName.Form("FFTrun-%s-det%i-Event%i",tag.Data(),idet,ievent);

  TCanvas * can = new TCanvas(canName,canName);
  hist->Draw("");
  fhist->Draw("same");
  can->Print(".png");
}


/** looked at fft ***/

void plotb(Int_t ievent=0, Int_t idet=0) {

  TH1D* hist=NULL;
  TH1D* phist=NULL;
  TH1D* bhist=NULL;

  if (!fin->IsOpen()) return;

  TString search;
  search.Form("EvSignal%1i_DET%1i",ievent,idet);
  TString psearch;
  psearch.Form("EvPeaks%1i_DET%1i",ievent,idet);
  printf(" looking for %s %s \n",search.Data(),psearch.Data());
  TString bsearch;
  bsearch.Form("EvBase%1i_DET%1i",ievent,idet);
  printf(" looking for %s %s \n",search.Data(),bsearch.Data());


  
  TList* list = evdir->GetListOfKeys() ;
  if (!list) { printf("<E> No keys found in file\n") ; exit(1) ; }
  TIter next(list);
  TKey* key;
  TObject* obj;

  while ( (key = (TKey*)next()) ) {
    obj = key->ReadObj() ;
    if (!obj->InheritsFrom("TH1D")) continue;
    TString hname(obj->GetName());
    if( hname.Contains(search)  ) hist = (TH1D*) obj;
    if( hname.Contains(psearch)) phist = (TH1D*) obj;
    if( hname.Contains(bsearch)) bhist = (TH1D*) obj;
  }

  if(!(hist&&phist&&bhist)) { 
    printf(" cannot find hist for event %i \n",ievent);
    return;
  }

  printf("found %s  %s %s \n",hist->GetName(), phist->GetName(), bhist->GetName() );

  TString histName;
  histName.Form("run-%s-det%i-Event%i",tag.Data(),idet,ievent);
  phist->SetTitle(histName);

  hist->SetLineColor(kBlack);
  hist->SetFillColor(kBlack);
  hist->SetFillStyle(0);
  
  phist->SetLineColor(kRed);
  phist->SetFillColor(kRed);
  phist->SetFillStyle(3001);

  bhist->SetLineColor(kGreen);



  TString canName;
  canName.Form("run-%s-det%i-Event%i",tag.Data(),idet,ievent);

  TCanvas * can = new TCanvas(canName,canName);
  phist->Draw("");
  hist->Draw("same");
  bhist->Draw("same");
  can->Print(".png");
}

//
void plotp(Int_t ievent=0, Int_t idet=0) {

  TH1D* hist=NULL;
  TH1D* phist=NULL;
  TH1D* dhist=NULL;
  TH1D* dphist=NULL;
  TH1D* bhist=NULL;

  if (!fin->IsOpen()) return;

  TString search;
  search.Form("EvSignal%1i_DET%1i",ievent,idet);
  TString psearch;
  psearch.Form("EvPeaks%1i_DET%1i",ievent,idet);
  TString dsearch;
  dsearch.Form("EvDWave%1i_DET%1i",ievent,idet);
  TString bsearch;
  bsearch.Form("EvBase%1i_DET%1i",ievent,idet);
  TString dpsearch;
  dpsearch.Form("EvDPeaks%1i_DET%1i",ievent,idet);


  printf(" looking for %s %s %s %s %s \n",search.Data(),dsearch.Data(),psearch.Data(),bsearch.Data(),dpsearch.Data());

  
  TList* list = evdir->GetListOfKeys() ;
  if (!list) { printf("<E> No keys found in file\n") ; exit(1) ; }
  TIter next(list);
  TKey* key;
  TObject* obj;

  while ( (key = (TKey*)next()) ) {
    obj = key->ReadObj() ;
    if (!obj->InheritsFrom("TH1D")) continue;
    TString hname(obj->GetName());
    //cout << hname << " " << endl; 
    if( hname.Contains(search))  hist = (TH1D*) obj;
    if( hname.Contains(dsearch)) dhist = (TH1D*) obj;
    if( hname.Contains(psearch)) phist = (TH1D*) obj;
    if( hname.Contains(bsearch)) bhist = (TH1D*) obj;
    if( hname.Contains(dpsearch)) dphist = (TH1D*) obj;
  }

  if(!hist) printf(" cannot find %s for event %i \n",search.Data(),ievent);
  if(!dhist) printf(" cannot find %s for event %i \n",dsearch.Data(),ievent);
  if(!phist) printf(" cannot find %s for event %i \n",psearch.Data(),ievent);
  if(!bhist) printf(" cannot find %s for event %i \n",bsearch.Data(),ievent);
  if(!dphist) printf(" cannot find %s for event %i \n",dpsearch.Data(),ievent);

  if(!hist||!dhist||!phist||!bhist||!dphist) {
    printf(" missing histo \n");
    return;
  }

  printf("found %s  %s \n",hist->GetName(), dhist->GetName() );

  dhist->SetLineColor(kRed);
  //dhist->SetFillColor(kRed);
  bhist->SetLineColor(kGreen);
  bhist->SetLineWidth(2);

  //dhist->SetFillStyle(3004);

  hist->SetLineColor(kBlue);
  hist->SetFillColor(kBlue);
  hist->SetFillStyle(0);


  TString canName;
  canName.Form("%s-Event%i",tag.Data(),ievent);

  TCanvas * can = new TCanvas(canName,canName);
  can->Divide(1,2);
  can->cd(1);
  hist->Draw();
  //bhist->Draw("same");
  phist->SetFillColor(kRed);
  phist->SetFillStyle(3003);
  phist->Draw("same");
  //if(shist) shist->Print("all");
  can->cd(2);
  dhist->Draw("");
  dphist->SetFillColor(kRed);
  dphist->SetFillStyle(3003);
  dphist->Draw("same");
  can->Print(".png");
}




void next(int ievent=0) {
  plote(ievent);
}

void plots(int nplot=10) {

  for(int i=0; i<nplot; ++i ){
    next(i);
    can->Print(".png");
  }

}

void plotEv(TString rtag = "9_25_2020") {
  tag = rtag;
  openFile();
  TFile *fout = new TFile( Form("plotEv%s",tag.Data()),"recreate");
  plote(0);
}
