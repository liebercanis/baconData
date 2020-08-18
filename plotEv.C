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



void plote(int ievent=0) {

  TH1D* hist[NDET];
  TH1D* phist[NDET];

  if (!fin->IsOpen()) return;

  TString search[NDET];
  TString psearch[NDET];
  for(int id=0; id<NDET; ++id) {
    hist[id]=NULL;
    search[id].Form("EvWave%1i_DET_%1i",ievent,id);
    phist[id]=NULL;
    psearch[id].Form("EvPeaks%1i_DET_%1i",ievent,id);

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

  canName.Form("%s-Peaks-All-Event%i",tag.Data(),ievent);
  TCanvas * pcan = new TCanvas(canName,canName);

  pcan->Divide(1,4);

  for(int id=0; id<NDET; ++id) {
    pcan->cd(id+1);
    phist[id]->SetLineColor(kRed);
    phist[id]->SetFillColor(kRed);
    phist[id]->SetFillStyle(3001);
    phist[id]->Draw();
  }
  //can->Print(".pdf");

}


void plotb(Int_t ievent=0, Int_t idet=0) {

  TH1D* hist=NULL;
  TH1D* phist=NULL;

  if (!fin->IsOpen()) return;

  TString search;
  search.Form("EvSignal%1i_DET_%1i",ievent,idet);
  TString psearch;
  psearch.Form("EvPeaks%1i_DET_%1i",ievent,idet);
  printf(" looking for %s %s \n",search.Data(),psearch.Data());

  
  TList* list = evdir->GetListOfKeys() ;
  if (!list) { printf("<E> No keys found in file\n") ; exit(1) ; }
  TIter next(list);
  TKey* key;
  TObject* obj;

  while ( (key = (TKey*)next()) ) {
    obj = key->ReadObj() ;
    if (!obj->InheritsFrom("TH1D")) continue;
    TString hname(obj->GetName());
    if( hname.Contains(search) && !hname.Contains(psearch) ) hist = (TH1D*) obj;;
    if( hname.Contains(psearch)) phist = (TH1D*) obj;;
  }

  if(!(hist&&phist)) { 
    printf(" cannot find hist for event %i \n",ievent);
    return;
  }

  printf("found %s  %s \n",hist->GetName(), phist->GetName() );

  hist->SetLineColor(kBlack);
  hist->SetFillColor(kBlack);
  hist->SetFillStyle(0);
  
  phist->SetLineColor(kRed);
  phist->SetFillColor(kRed);
  phist->SetFillStyle(3001);


  TString canName;
  canName.Form("%s-Raw-Event%i",tag.Data(),ievent);

  TCanvas * can = new TCanvas(canName,canName);
  phist->Draw("");
  hist->Draw("same");
  //hist->Draw("same");
  can->Print(".png");
}

//
void plotp(Int_t ievent=0, Int_t idet=0) {

  TH1D* hist=NULL;
  TH1D* phist=NULL;
  TH1D* dhist=NULL;

  if (!fin->IsOpen()) return;

  TString search;
  search.Form("EvSignal%1i_DET_%1i",ievent,idet);
  TString psearch;
  psearch.Form("EvPeaks%1i_DET_%1i",ievent,idet);
  TString dsearch;
  dsearch.Form("EvDWave%1i_DET_%1i",ievent,idet);

  printf(" looking for %s %s %s \n",search.Data(),dsearch.Data(),psearch.Data());

  
  TList* list = evdir->GetListOfKeys() ;
  if (!list) { printf("<E> No keys found in file\n") ; exit(1) ; }
  TIter next(list);
  TKey* key;
  TObject* obj;

  while ( (key = (TKey*)next()) ) {
    obj = key->ReadObj() ;
    if (!obj->InheritsFrom("TH1D")) continue;
    TString hname(obj->GetName());
    cout << hname << " " << endl; 
    if( hname.Contains(search))  hist = (TH1D*) obj;
    if( hname.Contains(dsearch)) dhist = (TH1D*) obj;
    if( hname.Contains(psearch)) phist = (TH1D*) obj;
  }

  if(!hist) printf(" cannot find %s for event %i \n",search.Data(),ievent);
  if(!dhist) printf(" cannot find %s for event %i \n",dsearch.Data(),ievent);
  if(!phist) printf(" cannot find %s for event %i \n",psearch.Data(),ievent);

  if(!hist||!dhist||!phist) return;


  printf("found %s  %s \n",hist->GetName(), dhist->GetName() );

  dhist->SetLineColor(kRed);
  dhist->SetFillColor(kRed);
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
  phist->SetFillColor(kRed);
  phist->SetFillStyle(3003);
  phist->Draw("same");
  //if(shist) shist->Print("all");
  can->cd(2);
  dhist->Draw("");
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

void plotEv(TString rtag = "MuTrigger") {
  tag = rtag;
  openFile();
  TFile *fout = new TFile( Form("plotEv%s",tag.Data()),"recreate");
  plote(0);
}
