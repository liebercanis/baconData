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
TH1D* hHitTimeUn[NDET];
TH1D* hHitQ[NDET];
TH1D* hHitMult;
TNtuple *nt;
TFile *fin;

using namespace TMath;

static double fexp(double *xx, double *par)
{
  double x = xx[0];
  double binwidth = par[2];
  double tau = par[0];
  //double cnorm = Exp(-tailStart/tau) - Exp(-tailStop/tau);
  double cnorm = 1.0;
  double f = par[1]/tau*binwidth/cnorm*Exp(-x/tau);
  return f; 
}

void fitIt(Int_t idet) 
{
  double binwidth = hHitTime[idet]->GetBinWidth(0);
  double xlow = 1.2;
  double xhigh = 8.0;
  int lowBin = hHitTime[idet]->FindBin(xlow);
  int highBin = hHitTime[idet]->FindBin(xhigh);
  double xmax = hHitTime[idet]->GetBinLowEdge(hHitTime[idet]->GetNbinsX());
  double integral = hHitTime[idet]->Integral(lowBin,highBin);
  printf(" fit for det %i from %f to %f integral %f \n",idet,xlow,xhigh,integral);
  TF1 * fe = new TF1(Form("fexpDet-%i",idet),fexp,xlow,xhigh,3);
  fe->SetNpx(1000); // numb points for function
  fe->SetParNames("lifetime","integral","binwidth");
  fe->SetParameters(1.4,integral,binwidth);
  fe->FixParameter(2,binwidth);

  fe->Print();
  for(int ii=0; ii<3; ++ii) {
    printf("\t  param %i %s %.4f +/- %.4f \n",ii,fe->GetParName(ii),fe->GetParameter(ii),fe->GetParError(ii));
  }

  TCanvas *can = new TCanvas(Form("tripletFit%s-det%i",runTag.Data(),idet),Form("tripletFit%s-det%i",runTag.Data(),idet));
  gStyle->SetOptStat(11);
  gStyle->SetOptDate();
  gPad->SetLogy();
  hHitTime[idet]->SetLineColor(kBlack);
  hHitTime[idet]->SetMarkerStyle(21);
  hHitTime[idet]->SetMarkerSize(0.2);
  hHitTime[idet]->Fit(fe,"RLE+","",xlow,xhigh);
  hHitTime[idet]->Draw("samesE1");
  //fe->Draw("sames");
  TPaveStats *ps = (TPaveStats*)can->GetPrimitive("stats");
  ps->SetTextColor(kBlack);
  gStyle->SetOptFit(111);
  can->Update();

  can->Print(".png");
}


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

void anaPost(Long64_t nev) {

  // for wave 1 
  // 2  Mean         6.10091e-04   1.12858e-06   9.50574e-08   7.16740e+01
  // 3  Sigma        5.67949e-04   8.68690e-07   2.26459e-05  -8.97484e-01
  double qsigmaDet1 = 5.67949e-04;
  double timeUnit = 2.0 ; // ns per count
  double microSec=1.0E-3;  
  Long64_t nevents = btree->GetEntries();
  if(nev!=0) nevents = nev;
  printf(" loop over  %lld out of %lld total \n",nevents, btree->GetEntries());
  for(Long64_t entry=0; entry<nevents; ++entry) {
    btree->GetEntry(entry);
    if(entry%1000==0) printf(" ... %lld \n",entry);
    double timeOffset = 1;
    // loop to get PMT time offset 
    for(unsigned id=0; id<NDET; ++id) {
      bool isPMT=false;
      if(TString(brun->detList[id]->GetName()).Contains("1")  ) isPMT=true;
      if(isPMT&&brun->detList[id]->hits.size()>0) timeOffset= brun->detList[id]->hits[0].startTime*timeUnit*microSec; 
    }
    //int sipmHits=0;
    for(unsigned id=0; id<NDET; ++id) {
      bool isPMT=false;
      if(TString(brun->detList[id]->GetName()).Contains("1")  ) isPMT=true;
      hHitMult->Fill(id,brun->detList[id]->hits.size());
      for(unsigned ihit=0; ihit<brun->detList[id]->hits.size(); ++ihit) {
        TDetHit dhit = brun->detList[id]->hits[ihit];
        double hitCharge = dhit.qsum*1.0E-6;
        double hitChargeError = dhit.qerr*1.0E-6;
        hHitQ[id]->Fill(hitCharge);
        double  hitqerr = sqrt(abs(hitCharge));
        // hit time relative to first hit 
        double  hitTime = dhit.startTime*timeUnit*microSec - timeOffset + 1.;
        int hitBin =  hHitTime[id]->FindBin(hitTime); 
        hHitTimeUn[id]->SetBinContent(hitBin, hHitTime[id]->GetBinContent(hitBin)+1);
        if(isPMT&&hitCharge>0) {
          hHitTime[id]->SetBinContent(hitBin, hHitTime[id]->GetBinContent(hitBin)+hitCharge);
          hHitTime[id]->SetBinError(hitBin, sqrt( pow(hHitTime[id]->GetBinError(hitBin),2)+pow(hitChargeError,2)));
        } else if(!isPMT&&hitCharge>0) {
          //printf("det %i hit %i bin %i q %f \n",id,ihit,hitBin,hitCharge);
          hHitTime[id]->SetBinContent(hitBin, hHitTime[id]->GetBinContent(hitBin)+hitCharge);
          hHitTime[id]->SetBinError(hitBin, sqrt( pow(hHitTime[id]->GetBinError(hitBin),2)+pow(hitChargeError,2)));
        }
        nt->Fill(id,dhit.firstBin,dhit.lastBin-dhit.firstBin,dhit.startTime*timeUnit*microSec,dhit.qsum*1.0E-6);
      }
    }
    //if(sipmHits>0) printf(" ... %lld sipmhits %d \n",entry,sipmHits);

  }
}

void plotHist() {

  TH1D* hbase[NDET]={0,0,0,0};
  TH1D* hnoise[NDET]={0,0,0,0};
  TList* list = fin->GetListOfKeys() ;
  if (!list) { printf("<E> No keys found in file\n") ; exit(1) ; }
  TIter next(list);
  TKey* key;
  TObject* obj;
  TH1D* hist;

  int igot=0;
  while ( (key = (TKey*)next()) && igot<NDET) {
    obj = key->ReadObj() ;
    if (!obj->InheritsFrom("TH1D")) continue;
    TString hname(obj->GetName());
    if( hname.Contains(TString("RawBase"))  ) {
      hbase[igot] = (TH1D*) obj;
    }
    if( hname.Contains(TString("RawNoise"))  ) {
      hnoise[igot] = (TH1D*) obj;
      ++igot;
    }
  }

  if(igot<NDET) {
    cout << " igot " << igot << " file  " << fin->GetName() << endl; return;
  }

  TCanvas *can = new TCanvas(Form("baseLine%s",runTag.Data()),Form("baseLine%s",runTag.Data()));
  can->Divide(2,2);
  for(unsigned i=0; i<NDET; ++i) { 
    can->cd(i+1);
    gPad->SetLogy();
    hbase[i]->Draw("");
  }

  
  TCanvas *canN = new TCanvas(Form("noise%s",runTag.Data()),Form("noise%s",runTag.Data()));
  canN->Divide(2,2);
  for(unsigned i=0; i<NDET; ++i) { 
    canN->cd(i+1); 
    gPad->SetLogy();
    hnoise[i]->Draw("");
  }
}


// main routine
void post(Long64_t nev=0, TString tag = "9_25_2020") {

  runTag = tag;
  openFile();

  plotHist();


  TFile *fout = new TFile( Form("post-%s.root",tag.Data()),"recreate");
  nt = new TNtuple("nt","hits","id:ibin:wid:time:qhit");
  hHitMult= new TH1D("hitMult"," hit multiplicity by det ",4,0,4);
  hHitMult->GetXaxis()->SetTitle(" det ");
  hHitMult->GetYaxis()->SetTitle(" hit multiplicity  ");

  int icolor[NDET] = {kBlack,kBlue,kGreen,kRed};
  for(unsigned id=0; id<NDET; ++id) {
    TString htit;
    htit.Form("hitTime-%s",brun->detList[id]->GetName());
    hHitTime[id]= new TH1D(htit,htit,4000,0,16);
    hHitTime[id]->Sumw2();
    hHitTime[id]->GetXaxis()->SetTitle(" micro-seconds ");
    hHitTime[id]->GetYaxis()->SetTitle(" hit charge ");
    hHitTime[id]->SetLineColor( icolor[id]) ;

    htit.Form("hitTimeUn-%s",brun->detList[id]->GetName());
    hHitTimeUn[id]= new TH1D(htit,htit,4000,0,16);
    hHitTimeUn[id]->Sumw2();
    hHitTimeUn[id]->GetXaxis()->SetTitle(" micro-seconds ");
    hHitTimeUn[id]->GetYaxis()->SetTitle(" hit charge ");
    hHitTimeUn[id]->SetLineColor( icolor[id]) ;

    htit.Form("hitQ-%s",brun->detList[id]->GetName());
    if(TString(brun->detList[id]->GetName()).Contains("1")  ) 
      hHitQ[id]= new TH1D(htit,htit,2000,-.2,1.8);
    else 
      hHitQ[id]= new TH1D(htit,htit,2000,-.005,0.015);
    hHitQ[id]->GetXaxis()->SetTitle(" hit Q  ");
    hHitQ[id]->SetLineColor( icolor[id]) ;
  }

  anaPost(nev);

  TCanvas *cant = new TCanvas(Form("hitTime%s",tag.Data()),Form("hitTime%s",tag.Data()));
  cant->Divide(2,2);
  for(unsigned id=0; id<NDET; ++id) { 
    cant->cd(id+1);
    hHitTime[id]->Draw("");
    gPad->Update();
    TPaveStats *st = (TPaveStats*) hHitTime[id]->FindObject("stats");
    st->SetTextColorAlpha(icolor[id],1);
  }

  TCanvas *canq = new TCanvas(Form("hitQ-Run%s",tag.Data()),Form("hitQRun%s",tag.Data()));
  gStyle->SetOptLogy();
  canq->Divide(2,2);
  for(unsigned id=0; id<NDET; ++id) { 
    canq->cd(id+1);
    hHitQ[id]->Draw("");
    gPad->Update();
    TPaveStats *stq = (TPaveStats*) hHitQ[id]->FindObject("stats");
    stq->SetTextColorAlpha(icolor[id],1);
  }

  
  //for(unsigned id=0; id<NDET; ++id) fitIt(id);
  fitIt(1);
  fout->Write();

}
