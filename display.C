////////////////////////////////////////////////////////
enum {NDET=3};
TTree *btree;
TBRawEvent *detList[NDET];
TH1D *hevent[NDET];
TString runTag;

void openFile() {
  // open ouput file and make some histograms
  TString fileName; fileName.Form("%s.root",runTag.Data());
  printf(" looking for file %s\n",fileName.Data());
  TFile *fin = new TFile(fileName,"readonly");
  if(fin->IsZombie()) {
    printf(" couldnt open file %s\n",fileName.Data());
    return;
  }

  fin->GetObject("BTree",btree);
  detList[0] = new TBRawEvent("det6");
  detList[1] = new TBRawEvent("det4");
  detList[2]= new TBRawEvent("det2");
  btree->SetBranchAddress("det6", &detList[0]);
  btree->SetBranchAddress("det4", &detList[1]);
  btree->SetBranchAddress("det2", &detList[2]);
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

void display(TString tag = "run1") {
  runTag = tag;
  openFile();
  TFile *fout = new TFile( Form("display%s",tag.Data()),"recreate");
  next(0);
}
