////////////////////////////////////////////////////////
// analyze anaRun output data structure TTree BRunTree
// data structures are all in directory bobj "bacon object" 
#include "TBRun.cxx"
enum {NDET=4};
TTree *btree;  // pointer to TTree BRunTree on input file
TBRun *brun;   // TTRee to hold data internally
TH1D *hevent[NDET];
TString runTag;
TH1D* hHitTime[NDET];
TH1D* hHitQ[NDET];
TH1D* hHitMult;
TNtuple *nt;
TFile *fin;

using namespace TMath;  // lots of useful functions defined here


// open file and attach to local data structure brun 
void openFile() {
  // open ouput file and make some histograms
  TString fileName; fileName.Form("ana-%s.root",runTag.Data());
  printf(" looking for file %s\n",fileName.Data());
  fin = new TFile(fileName,"readonly");
  if(fin->IsZombie()) {
    printf(" couldnt open file %s\n",fileName.Data());
    return;
  }

  fin->GetObject("BRunTree",btree);
  brun = new TBRun(btree,runTag);
    
  for(unsigned i=0; i<NDET; ++i ) {
    cout << " btree adding branch det " << brun->detList[i]->GetName() << endl;
    btree->SetBranchAddress(brun->detList[i]->GetName(),&(brun->detList[i]));
  }
  cout << " btree adding branch event " << brun->bevent->GetName() << endl;
  btree->SetBranchAddress("bev",&(brun->bevent));

}

void anaPost(Long64_t nev=0, TString tag = "9_25_2020") {

  runTag = tag;
  openFile();

  // make output file
  TFile *fout = new TFile( Form("anaPost-%s.root",tag.Data()),"recreate");
  // define ntuple
  nt = new TNtuple("nt","hits","id:ibin:wid:time:qhit");
  hHitMult= new TH1D("hitMult"," hit multiplicity by det ",4,0,4);
  hHitMult->GetXaxis()->SetTitle(" det ");
  hHitMult->GetYaxis()->SetTitle(" hit multiplicity  ");

  //histograms by detector
  int icolor[NDET] = {kBlack,kBlue,kGreen,kRed};
  for(unsigned id=0; id<NDET; ++id) {
    TString htit;
    htit.Form("hitTime-%s",brun->detList[id]->GetName());
    hHitTime[id]= new TH1D(htit,htit,4000,0,16);
    hHitTime[id]->Sumw2();
    hHitTime[id]->GetXaxis()->SetTitle(" micro-seconds ");
    hHitTime[id]->GetYaxis()->SetTitle(" hit charge ");
    hHitTime[id]->SetLineColor( icolor[id]) ;
    htit.Form("hitQ-%s",brun->detList[id]->GetName());
    // event1 is the PMT
    if(TString(brun->detList[id]->GetName()).Contains("1")  ) 
      hHitQ[id]= new TH1D(htit,htit,2000,-.2,1.8);
    else 
      hHitQ[id]= new TH1D(htit,htit,2000,-.005,0.015);
    hHitQ[id]->GetXaxis()->SetTitle(" hit Q  ");
    hHitQ[id]->SetLineColor( icolor[id]) ;
  }

  double qsigmaDet1 = 5.67949e-04;
  double timeUnit = 2.0 ; // ns per count
  double microSec=1.0E-3;  
  Long64_t nevents = btree->GetEntries();
  if(nev!=0) nevents = nev;
  printf(" loop over  %lld out of %lld total \n",nevents, btree->GetEntries());

  // loop over events
  for(Long64_t entry=0; entry<nevents; ++entry) {
    btree->GetEntry(entry); // read event into local data structure brun
    if(entry%1000==0) printf(" ... %lld \n",entry);
    // loop over detectors
    for(unsigned id=0; id<NDET; ++id) {
      bool isPMT=false;
      if(TString(brun->detList[id]->GetName()).Contains("1")  ) isPMT=true;
      hHitMult->Fill(id,brun->detList[id]->hits.size());
      // loop over hits
      for(unsigned ihit=0; ihit<brun->detList[id]->hits.size(); ++ihit) {
        TDetHit dhit = brun->detList[id]->hits[ihit]; // hit structure 
        double hitCharge = dhit.qsum*1.0E-6;
        double hitChargeError = dhit.qerr*1.0E-6;
        hHitQ[id]->Fill(hitCharge);
        double  hitqerr = sqrt(abs(hitCharge));
        double  hitTime = dhit.startTime*timeUnit*microSec;
        int hitBin =  hHitTime[id]->FindBin(hitTime); 
        if(isPMT&&hitCharge>3*qsigmaDet1) {
          hHitTime[id]->SetBinContent(hitBin, hHitTime[id]->GetBinContent(hitBin)+hitCharge);
          hHitTime[id]->SetBinError(hitBin, sqrt( pow(hHitTime[id]->GetBinError(hitBin),2)+pow(hitChargeError,2)));
        } else if(!isPMT&&hitCharge>0) {
          //printf("det %i hit %i bin %i q %f \n",id,ihit,hitBin,hitCharge);
          hHitTime[id]->SetBinContent(hitBin, hHitTime[id]->GetBinContent(hitBin)+hitCharge);
          hHitTime[id]->SetBinError(hitBin, sqrt( pow(hHitTime[id]->GetBinError(hitBin),2)+pow(hitChargeError,2)));
        }
        // fill ntuple
        nt->Fill(id,dhit.firstBin,dhit.lastBin-dhit.firstBin,dhit.startTime*timeUnit*microSec,dhit.qsum*1.0E-6);
      }
    }
    //if(sipmHits>0) printf(" ... %lld sipmhits %d \n",entry,sipmHits);

  }
   fout->Write();
}


