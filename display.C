////////////////////////////////////////////////////////
enum {NDET=4};
TTree *btree;
TBRawEvent *detList[NDET];
TH1D *hevent[NDET];
TString runTag;

void openFile() {
  // open ouput file and make some histograms
  TString fileName; fileName.Form("rootData/%s.root",runTag.Data());
  printf(" looking for file %s\n",fileName.Data());
  TFile *fin = new TFile(fileName,"readonly");
  if(fin->IsZombie()) {
    printf(" couldnt open file %s\n",fileName.Data());
    return;
  }

  fin->GetObject("BTree",btree);
  TObjArray *brList = btree->GetListOfBranches();

  TIter next(brList);
  TBranch *aBranch=NULL;
  int idet=0;

  while( ( aBranch = (TBranch *) next() ) ) {
    detList[idet] = new TBRawEvent(aBranch->GetName());
    btree->SetBranchAddress(aBranch->GetName() , &detList[idet]);
    ++idet;
  }

  for(unsigned id=0; id<NDET; ++id) hevent[id]=NULL;
  cout << " number of entries in BTree is " <<  btree->GetEntries() << endl;
  btree->GetListOfBranches()->ls();

}

void eventDisplay(Int_t  ievent=0)
{
  btree->GetEntry(ievent);
  int nsamples = (int) detList[0]->digi.size();
  cout << detList[0]->GetName() << " nsamples " << nsamples << endl;
  if(!hevent[0]) for(unsigned id=0; id<NDET; ++id)  hevent[id] = 
    new TH1D(Form("Event%s",detList[id]->GetName()),Form("Event%s",detList[id]->GetName()), nsamples,0, nsamples);
  
  // fill histograms
  for(unsigned id=0; id<NDET; ++id) for(unsigned i=0; i< nsamples ; ++i) hevent[id]->SetBinContent(i+1,detList[id]->digi[i]);
  
  TCanvas *can = new TCanvas(Form("%sEv%i",runTag.Data(),ievent),Form(" run %s Event %i",runTag.Data(),ievent));
  can->Divide(1,NDET);
  for(unsigned id=0; id<NDET; ++id) {
    can->cd(id+1);
    hevent[id]->Draw();
  }


}

void next(int ievent=0) {
  eventDisplay(ievent);
}

void display(TString tag = "7_27_2020") {
  runTag = tag;
  openFile();
  TFile *fout = new TFile( Form("display%s",tag.Data()),"recreate");
  next(0);
}
