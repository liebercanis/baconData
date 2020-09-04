////////////////////////////////////////////////////////
#include "TBRun.cxx"
enum {NDET=4};
TTree *btree;
TBRun *brun;
TBRawEvent *detList[NDET];
TH1D *hevent[NDET];
TString runTag;
TCanvas *can;
TBEvent *bevent;
TH1D* hHitTime[NDET];
TH1D* hHitQ[NDET];
TNtuple *nt;


void openFile() {
  // open ouput file and make some histograms
  TString fileName; fileName.Form("rootData/ana-%s.root",runTag.Data());
  printf(" looking for file %s\n",fileName.Data());
  TFile *fin = new TFile(fileName,"readonly");
  if(fin->IsZombie()) {
    printf(" couldnt open file %s\n",fileName.Data());
    return;
  }

  fin->GetObject("BRunTree",btree);

    
  for(unsigned i=0; i<NDET; ++i ) {
    cout << " btree adding branch " << brun->detList[i]->GetName() << endl;
    btree->SetBranchAddress(brun->detList[i]->GetName(),&(brun->detList[i]));
  }
  cout << " btree adding branch " << brun->bevent->GetName() << endl;
  btree->SetBranchAddress("bev",&(brun->bevent));

}

void anaPost(Long64_t nev) {
  double timeUnit = 2.0 ; // ns per count
  double microSec=1.0E-3;  
  Long64_t nevents = btree->GetEntries();
  if(nev!=0) nevents = nev;
  printf(" loop over  %lld out of %lld total \n",nevents, btree->GetEntries());
  for(Long64_t entry=0; entry<nevents; ++entry) {
    btree->GetEntry(entry);
    if(entry%100==0) printf(" ... %lld \n",entry);
    //int sipmHits=0;
    for(unsigned id=0; id<NDET; ++id) {
      //if( id>0) sipmHits += brun->detList[id]->hits.size();
      for(unsigned ihit=0; ihit<brun->detList[id]->hits.size(); ++ihit) {
        TDetHit dhit = brun->detList[id]->hits[ihit];
        hHitTime[id]->Fill(dhit.startTime*timeUnit*microSec);
        hHitQ[id]->Fill(dhit.qsum*1.0E-6);
        nt->Fill(id,dhit.qsum*1.0E-6, dhit.startTime*timeUnit*microSec);
      }
    }
    //if(sipmHits>0) printf(" ... %lld sipmhits %d \n",entry,sipmHits);

  }
}



void post(Long64_t nev=0, TString tag = "3.0VLED_51.0VSiPM") {
  
  brun = new TBRun(tag);
  brun->det0.SetName(Form("wave%i",1));
  brun->det1.SetName(Form("wave%i",4));
  brun->det2.SetName(Form("wave%i",2));
  brun->det3.SetName(Form("wave%i",6));


  runTag = tag;
  openFile();
  TFile *fout = new TFile( Form("post-%s.root",tag.Data()),"recreate");
  int icolor[NDET] = {kBlack,kBlue,kGreen,kRed};
  for(unsigned id=0; id<NDET; ++id) {
    TString htit;
    htit.Form("hitTime-%s",brun->detList[id]->GetName());
    hHitTime[id]= new TH1D(htit,htit,80,0,16);
    hHitTime[id]->GetXaxis()->SetTitle(" micro-seconds ");
    hHitTime[id]->SetLineColor( icolor[id]) ;
    htit.Form("hitQ-%s",brun->detList[id]->GetName());
    hHitQ[id]= new TH1D(htit,htit,200,0,1);
    hHitQ[id]->GetXaxis()->SetTitle(" hit Q  ");
    hHitQ[id]->SetLineColor( icolor[id]) ;
  }
  nt = new TNtuple("nt","hits","id:q:t");

  anaPost(nev);

  TCanvas *can = new TCanvas(Form("hitTime%s",tag.Data()),Form("hitTime%s",tag.Data()));
  
  hHitTime[1]->Draw(""); 
  for(unsigned id=0; id<NDET; ++id) { 
    hHitTime[id]->Draw("sames");
    gPad->Update();
    TPaveStats *st = (TPaveStats*) hHitTime[id]->FindObject("stats");
    st->SetTextColorAlpha(icolor[id],1);
  }

  TCanvas *canq = new TCanvas(Form("hitQ%s",tag.Data()),Form("hitQ%s",tag.Data()));
  hHitQ[1]->Draw(""); 
  for(unsigned id=0; id<NDET; ++id) { 
    hHitQ[id]->Draw("sames");
    gPad->Update();
    TPaveStats *stq = (TPaveStats*) hHitQ[id]->FindObject("stats");
    stq->SetTextColorAlpha(icolor[id],1);
  }


}
