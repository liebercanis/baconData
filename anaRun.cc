//////////////////////////////////////////////////////////
//  M.Gold August 2020
//////////////////////////////////////////////////////////
#include <sstream>
#include <unistd.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <complex> //includes std::pair, std::make_pair
#include <valarray>
//
#include <TROOT.h>
#include <TVirtualFFT.h>
#include <TChain.h>
#include <TMath.h>
#include <TNtuple.h>
#include <TFile.h>
#include <Rtypes.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>
#include <TFormula.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <algorithm> // std::sort
#include "TSpectrum.h"
#include "TRandom3.h"
//
#include "TBRawEvent.hxx"

class anaRun
{
  public:
    enum {NDET=4};
    TTree *btree;
    TBRawEvent *detList[NDET];
    TFile *fin;
    TFile *fout;
    TString runTag;

    anaRun(TString tag = "7_27_2020", Int_t maxEvents = 0);
    virtual ~anaRun() { ; }
    bool openFile();
    void analyze(Long64_t maxEvents);
    
    TH1D *hEvent[NDET];
    TH1D *hNoise[NDET];
    TH1D *hBase[NDET];
    TNtuple *nt;
};

////////////////////////////////////////////////////////
anaRun::anaRun(TString tag, Int_t maxEvents)
{
  runTag = tag;
  fin=NULL;
  if(!openFile()) { cout << " file with tag not found " << endl; return; }
  cout <<  " opened file with tag " << endl;

  fout = new TFile( Form("ana-%s.root",tag.Data()),"recreate");
  btree->GetEntry(0);
  int nsamples = (int) detList[0]->digi.size();
  cout << tag << " samples " << nsamples << endl;
  for(unsigned id=0; id<NDET; ++id) { 
    hEvent[id] = new TH1D(Form("Event%s",detList[id]->GetName()),Form("Sum%s",detList[id]->GetName()), nsamples,0, nsamples);
    hBase[id] = new TH1D(Form("Base%s",detList[id]->GetName()),Form("Base det %s",detList[id]->GetName()), 3000,0, 3000);
    hNoise[id] = new TH1D(Form("Noise%s",detList[id]->GetName()),Form("Noiese det %s",detList[id]->GetName()), 500,0, 200);
    hNoise[id]->GetXaxis()->SetTitle(" ADC counts ");
    hNoise[id]->GetYaxis()->SetTitle(" frequency ");
  }
  nt = new TNtuple("nt",Form("nt%s",tag.Data()),"b0:n0:b1:n1:b2:n2");
 

  analyze(maxEvents);

  fin->Close();
  fout->Write();
  //fout->Close();

}

void anaRun::analyze(Long64_t nmax)
{
  Long64_t nentries = btree->GetEntries();
  if(nmax>0) nentries=nmax;
  cout << " analyze " << nentries << endl;
  for(Long64_t ievent=0; ievent< nentries; ++ievent) {
    btree->GetEntry(ievent);
    // fill histograms
    double baseline[NDET];
    double noise[NDET];
    for(unsigned id=0; id<NDET; ++id) {
      for(unsigned long i=0; i< detList[id]->digi.size() ; ++i) hEvent[id]->SetBinContent(i+1, hEvent[id]->GetBinContent(1+i)+detList[id]->digi[i]);
      vector<double> digi = detList[id]->digi;
      std::sort(digi.begin(), digi.end());
      double nsamples = (double) detList[id]->digi.size();
      baseline[id] = digi[0.5*nsamples];
      noise[id] = std::abs( digi[0.68*nsamples] - baseline[id]);
      hBase[id]->Fill(baseline[id]);
      hNoise[id]->Fill(noise[id]);
    }
    nt->Fill( baseline[0],noise[0],baseline[1],noise[1],baseline[2],noise[2]);
  }

}

bool  anaRun::openFile() 
{
  // open ouput file and make some histograms
  TString fileName; fileName.Form("rootData/%s.root",runTag.Data());
  printf(" looking for file %s\n",fileName.Data());
  fin = new TFile(fileName,"readonly");
  if(fin->IsZombie()) {
    printf(" couldnt open file %s\n",fileName.Data());
    return false;
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

  for(unsigned id=0; id<NDET; ++id) hEvent[id]=NULL;
  cout << " number of entries in BTree is " <<  btree->GetEntries() << endl;
  btree->GetListOfBranches()->ls();
  return true;
}

