/////////////////////////////////////////////////////////
//  make data sum versus time plots
//  M.Gold Sept 2020
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
#include "TBRun.hxx"

class anaSum
{
  public:
    enum {NDET=4};
    enum {UPCROSS,DOWNCROSS};
    double vsign[NDET];
    TTree *btree;
    TBRawEvent *detList[NDET];
    TFile *fin;
    TFile *fout;
    TBRun *brun;
    TString runTag;
    enum { minLength=7};
    enum { baseLineHalfWindow = 200 }; // even integer

    anaSum(Int_t maxEvents = 0, TString tag = "MuTrigger");
    virtual ~anaSum() { ; }
    bool openFile();
    void analyze(Long64_t maxEvents);

    // summed wave histograms
    TH1D *hSumWave[NDET];
    TH1D *hSumEvent[NDET];

};

////////////////////////////////////////////////////////
anaSum::anaSum(Int_t maxEvents, TString tag)
{
  runTag = tag;
  fin=NULL;
  if(!openFile()) { cout << " file with tag not found " << endl; return; }
  cout <<  " opened file with tag " << endl;

  fout = new TFile( Form("anaSum-%s.root",tag.Data()),"recreate");
  fout->cd();

  // get sample size 
  btree->GetEntry(0);
  int nsamples = (int) detList[0]->digi.size();
  cout << tag << " samples " << nsamples << endl;

  for(unsigned id=0; id<NDET; ++id) hSumWave[id] = new TH1D(Form("SumWave%s",detList[id]->GetName()),Form("Sum Wave%s",detList[id]->GetName()), nsamples,0, nsamples);

  analyze(maxEvents);

  fin->Close();
  fout->Write();

  TString cname;
  cname.Form("EventWave%s",tag.Data());
  TCanvas *can2 = new TCanvas(cname,cname);
  can2->Divide(1,4);
  for(unsigned id=0; id<NDET; ++id) {
    can2->cd(id+1);
    hSumWave[id]->Draw(); 
  }
  can2->Print(".pdf");
  //fout->Close()


  cname.Form("EventSum%s",tag.Data());
  TCanvas *can = new TCanvas(cname,cname);
  can->Divide(2,2);
  gStyle->SetOptLogy();
  for(unsigned id=0; id<NDET; ++id) {
    can->cd(id+1);
    hSumEvent[id]->Draw(); 
  }
  can->Print(".pdf");

 
  //fout->Close()

}

void anaSum::analyze(Long64_t nmax)
{
  double vsign[NDET]={ 1., 1., 1., 1. };
  //
  Long64_t nentries = btree->GetEntries();
  if(nmax>0) nentries=nmax;
  cout << " analyze " << nentries << endl;
  for(Long64_t ievent=0; ievent< nentries; ++ievent) {
    btree->GetEntry(ievent);
    if (ievent%100==0) printf(" ... %lld \n", ievent);
    // loop over detectors
    //
    if(ievent==0) {
      printf(" det names \n");
      for(unsigned id=0; id<NDET; ++id) { 
        if(TString(detList[id]->GetName()).Contains(TString("1"))) {
          vsign[id]=-1.0;
          hSumEvent[id] = new TH1D(Form("SumEvent%s",detList[id]->GetName()),Form("SumEv Wave%s",detList[id]->GetName()),600,-10,50);
        } else {
          hSumEvent[id] = new TH1D(Form("SumEvent%s",detList[id]->GetName()),Form("SumEv Wave%s",detList[id]->GetName()),300,-1,2);
        }
        for(unsigned id=0; id<NDET; ++id) printf("%u %s sign %0.1f \n",id, detList[id]->GetName(), vsign[id] );
      }
    }

    for(unsigned id=0; id<NDET; ++id) {
      // get baseline copy vector to sort
      vector<double> sdigi = detList[id]->digi;
      std::sort(sdigi.begin(), sdigi.end());
      double nsamples = (double) detList[id]->digi.size();
      double baseline = sdigi[0.5*nsamples];
      double evSum=0;
      //cout << " idet " << id << " baseline " << baseline << endl;
      for (int ibin=0; ibin < nsamples ; ++ibin) { 
        double qbin = vsign[id]*(detList[id]->digi[ibin] -  baseline);
        evSum += qbin;
      }
      evSum /= double(nsamples);
      hSumEvent[id]->Fill(evSum);
      //cout << " idet " << id << " baseline " << baseline << endl;
      if(evSum>0.5) {
        for (int ibin=0; ibin < nsamples ; ++ibin) { 
          double qbin = vsign[id]*(detList[id]->digi[ibin] -  baseline);
          hSumWave[id]->SetBinContent(ibin, hSumWave[id]->GetBinContent(ibin)+qbin);
          hSumWave[id]->SetBinError(ibin, sqrt(pow(hSumWave[id]->GetBinError(ibin), 2)+pow(qbin, 2)));
        }
      }

      //cout << ievent << "  " << id << " "   << detList[id]->GetName() << " " << evSum/double(nsamples) << endl;
    } // loop over dets
  } // sum over events
}

bool  anaSum::openFile() 
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
    detList[idet]->SetName(aBranch->GetName());
    ++idet;
  }

  cout << " number of entries in BTree is " <<  btree->GetEntries() << endl;
  btree->GetListOfBranches()->ls();
  return true;
}

