////////////////////////////////////////////////////////
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
#include "TBRun.hxx"

typedef std::map<Double_t,TDetHit,std::less<Double_t> >  hitMap;
typedef std::map<Double_t,TDetHit,std::less<Double_t> >::iterator  hitMapIter;

typedef std::vector<std::pair<unsigned,unsigned> >  peakType;
typedef std::vector<std::pair<unsigned,unsigned> >::iterator  peakTypeIter;

const Double_t qnorm=1.0;
class anaRun
{
  public:
    enum {NDET=4};
    enum {UPCROSS,DOWNCROSS};
    double vsign[NDET];
    TTree *btree;
    TBRawEvent *detList[NDET];
    ofstream dumpFile;
    TFile *fin;
    TFile *fout;
    TBRun *brun;
    TString runTag;
    enum { minLength=7};
    enum { baseLineHalfWindow = 200 }; // even integer

    anaRun(Int_t maxEvents = 0, TString tag = "9_25_2020");
    virtual ~anaRun() { ; }
    bool openFile();
    void analyze(Long64_t maxEvents);
    std::vector<Double_t> differentiate(std::vector<Double_t> v, unsigned diffStep); 
    void getAverage(std::vector<Double_t> digi, Double_t& ave, Double_t& sigma,unsigned max=0);
    void getAverages(std::vector<Double_t> digi, Double_t& ave, Double_t& sigma,unsigned max=0);
    void plotWave(Long64_t ievent,unsigned idet);
    std::vector<Double_t> getBaselineWeights(unsigned arraySize, peakType peakList, Int_t& maxWeightWidth);
    std::vector<Double_t> getBaselineWMARecursive(Double_t ave, std::vector<Double_t> sig, std::vector<Double_t> weight, Int_t NWindow);
    peakType  derivativePeaks(std::vector<Double_t> v, Int_t idet , Int_t nwindow, Double_t rms, std::vector<Int_t>& peakKind);
    peakType  simplePeaks(std::vector<Double_t> v, Int_t idet , unsigned minWidth, unsigned maxWidth, Double_t rms, std::vector<Int_t>& peakKind);
    void extendPeaks(std::vector<Double_t> v, peakType &peakList );


    hitMap makeHits(peakType peakList, std::vector<Int_t> peakKind, std::vector<Double_t> digi, Double_t sigma, Double_t& firstTime, Double_t& firstCharge);

    double timeUnit;
    double derivativeSigma[NDET];
    double microSec;
    double maxLife;
    Int_t nSigma;
    Int_t aveWidth;
    unsigned diffStep;


    TH1D *hEvent[NDET];
    TH1D *hNoise[NDET];
    TH1D *hBase[NDET];
    TH1D *hDerivative[NDET];
    TH1D *hEvWave[NDET];
    TH1D *hEvDWave[NDET];
    TH1D *hEvSignal[NDET];
    TH1D *hEvPeaks[NDET];
    TH1D *hBaselineWMA[NDET];

    TH1D *hWave[NDET];
    TH1D *hEvWeight[NDET];
    TH1D *hDWave[NDET];
    TH1D *hDNoise[NDET];
    TH1D *hSNoise[NDET];
    TH1D *hLife[NDET];
    TH1I* hPeakNWidth;
    TH1I* hHitLength;
    TH1D *hPeakCount;
    TNtuple *nt;
    TNtuple *ntDer;
    TNtuple *ntSimplePeaks;
    TNtuple *ntWave;
    TNtuple *ntCal;
    TNtuple *ntHit;

    TDirectory *evDir;

};

////////////////////////////////////////////////////////
anaRun::anaRun(Int_t maxEvents, TString tag)
{
  runTag = tag;
  fin=NULL;
  timeUnit = 2.0 ; // ns per count
  microSec=1.0E-3;
  diffStep=7;
  nSigma=5;
  aveWidth=20;

  if(!openFile()) { cout << " file with tag not found " << endl; return; }
  cout <<  " opened file with tag " << endl;

  btree->GetEntry(0);
  for(unsigned id=0; id< NDET; ++id) {
    derivativeSigma[id]=3.0;
    vsign[id]=1.0;
    TString dname = detList[id]->GetName();
    if(dname.Contains( TString("1")))  {
      derivativeSigma[id]=5.0;
      vsign[id]=-1.0;
    }
    printf(" setting  derivative Sigma for %s to %.3f \n",dname.Data(),derivativeSigma[id]);
  }
  double baseOffset = 8100;
  int nsamples = (int) detList[0]->digi.size();
  maxLife = double(nsamples)*timeUnit*microSec; // microseconds
  cout << tag << " samples " << nsamples << endl;
  for(unsigned id=0; id<NDET; ++id) { 
    hEvWave[id] = new TH1D(Form("EvWave%s",detList[id]->GetName()),Form("Wave%s",detList[id]->GetName()), nsamples,0, nsamples);
    hEvDWave[id] = new TH1D(Form("EvDWave%s",detList[id]->GetName()),Form("DWave%s",detList[id]->GetName()), nsamples,0, nsamples);
    hEvSignal[id] = new TH1D(Form("EvSignal%s", detList[id]->GetName()), Form("DetSignal%s", detList[id]->GetName()), nsamples, 0, nsamples);
    hEvPeaks[id] = new TH1D(Form("EvPeaks%s", detList[id]->GetName()), Form("DetPeaks%s", detList[id]->GetName()), nsamples, 0, nsamples);
    hEvWeight[id] = new TH1D(Form("EvWeight%s",detList[id]->GetName()),Form("baseline weight det %s",detList[id]->GetName()), nsamples,0, nsamples);
    hBaselineWMA[id]=new TH1D(Form("BaselineWMA%s", detList[id]->GetName()), Form("BaselineWMA%s", detList[id]->GetName()), nsamples, 0, nsamples);
  }


  fout = new TFile( Form("ana-%s.root",tag.Data()),"recreate");
  evDir = fout->mkdir("EvDir");
  fout->cd();

  brun = new TBRun(btree,tag);

  printf(" ***** assigned det names *****  %lu \n", brun->detList.size() );
  for(unsigned id=0; id<brun->detList.size() ; ++id) printf(" %i %s %s \n",id,brun->detList[id]->GetName() , detList[id]->GetName() );

  // open dump file
  dumpFile.open (Form("dump-%s.txt",brun->bevent->GetName()),ios::out);
  dumpFile << "C****** open  dump for run " <<  brun->bevent->GetName() << endl;
  cout << " ****** open  dump for run " <<  brun->bevent->GetName() << endl;


  for(unsigned id=0; id<NDET; ++id) { 
    hEvent[id] = new TH1D(Form("Event%s",detList[id]->GetName()),Form("Sum%s",detList[id]->GetName()), nsamples,0, nsamples);
    hBase[id] = new TH1D(Form("RawBase%s",detList[id]->GetName()),Form("Base det %s",detList[id]->GetName()), 400,baseOffset-200,baseOffset+200);
    hNoise[id] = new TH1D(Form("RawNoise%s",detList[id]->GetName()),Form("Noise det %s",detList[id]->GetName()), 50,0, 50);
    hNoise[id]->GetXaxis()->SetTitle(" ADC counts ");
    hNoise[id]->GetYaxis()->SetTitle(" frequency ");
    hDNoise[id] = new TH1D(Form("DNoise%s",detList[id]->GetName()),Form("DNoise%s",detList[id]->GetName()), 1000,-200,200);
    hSNoise[id] = new TH1D(Form("SNoise%s", detList[id]->GetName()), Form("SNoise det %s", detList[id]->GetName()), 50, -50, 50);
    hLife[id] = new TH1D(Form("Life%i", id), Form(" lifetime DET %i ", id), 500, 0, maxLife);
    hLife[id]->Sumw2();
  }
  hPeakCount = new TH1D("PeakCount"," peaks by det ",NDET,0,NDET);
  hHitLength = new TH1I("HitLength", " hit length", 100, 0, 100);
  hPeakNWidth = new TH1I("PeakNWidth", "PeakNWidth", 100, 0, 100);
  nt = new TNtuple("ntBase", Form("ntBase%s", tag.Data()), "id:rawBase:rawNoise:rawAve:rawSigma");
  ntDer =  new TNtuple("ntDer"," deriviative ","t:sigma:d0:kover:type");
  ntSimplePeaks =  new TNtuple("ntSimplePeaks"," simple peaks ","idet:rms:v:t:width");
  ntWave = new TNtuple("ntWave"," wave ","event:v:d");
  ntCal =  new TNtuple("ntCal","ntuple Cal","iev:ipmt:base:sigma:dbase:dsigma:sbase:ssigma");
  ntHit = new TNtuple("ntHit", "ntuple Hit", "det:hits:order:istart:time:q:width:peak:good:kind");



  analyze(maxEvents);

  fin->Close();
  fout->Write();
  if (dumpFile.is_open()) { 
    cout << " ****** close dump for run " <<  brun->bevent->GetName() << endl;
    dumpFile << "C****** close dump for run " <<  brun->bevent->GetName() << endl;
    dumpFile.close();
  }

//fout->Close()

}

void anaRun::analyze(Long64_t nmax)
{
  // reset wave histos
  for (unsigned id=0; id<NDET; ++id) {
    hEvWave[id]->Reset();
    hEvDWave[id]->Reset();
    hEvSignal[id]->Reset();
    hEvPeaks[id]->Reset();
    hEvWeight[id]->Reset();
    hBaselineWMA[id]->Reset();
  }
  //
  Long64_t nentries = btree->GetEntries();
  if(nmax>0) nentries=nmax;
  cout << " analyze " << nentries << endl;
  for(Long64_t ievent=0; ievent< nentries; ++ievent) {
    if (ievent%100==0) printf(" ... %lld \n", ievent);
    brun->clear();
    btree->GetEntry(ievent);
    // fill raw ntuple... these are not used
    double rawBaseline[NDET];
    double noise[NDET];
    for(unsigned id=0; id<NDET; ++id) {
      for(unsigned long i=0; i< detList[id]->digi.size() ; ++i) hEvent[id]->SetBinContent(i+1, hEvent[id]->GetBinContent(1+i)+detList[id]->digi[i]);
      vector<double> adigi = detList[id]->digi;
      std::sort(adigi.begin(), adigi.end());
      double nsamples = (double) detList[id]->digi.size();
      rawBaseline[id] = adigi[0.5*nsamples];
      noise[id] = std::abs( adigi[0.68*nsamples] - rawBaseline[id]);
    }

    fout->cd();
    double trigTime=0;

    for(int id = 0; id < NDET; id++){
      Double_t derAve, derSigma;
      Double_t rawAve, rawSigma;
      getAverages(detList[id]->digi, rawAve, rawSigma, 600); // first 600 bins
      hBase[id]->Fill(rawAve);
      hNoise[id]->Fill(rawSigma);

      // first subtract rawBaseline
      nt->Fill(id,rawBaseline[id],noise[id],rawBaseline[id],noise[id],rawAve,rawSigma);
      //printf(" det %i raw ave %f base %f noise %f %f \n",id,rawAve,rawBaseline[id] ,rawSigma, noise[id] );
      std::vector<double> ddigi;
      for (unsigned i = 0; i < detList[id]->digi.size(); i++) ddigi.push_back(detList[id]->digi[i]-rawAve); 
 
      // make the derivative waveform 
      std::vector<double> vderiv;
      vderiv = differentiate(ddigi,diffStep);
      getAverages(vderiv,derAve,derSigma);
      for(unsigned isample = 0; isample < vderiv.size(); isample++){
        hEvWave[id]->SetBinContent(isample+1,ddigi[isample]);
        hEvDWave[id]->SetBinContent(isample+1,vderiv[isample]);
        ntWave->Fill(id,ddigi[isample],vderiv[isample]);
      }

      // derAVe, derSigma  best to just fit Gaussian
      hDNoise[id]->Reset();
      for(unsigned isample = 0; isample < vderiv.size(); isample++) hDNoise[id]->Fill(vderiv[isample]);
      hDNoise[id]->Fit("gaus","0Q");
      TF1 *fit1 = (TF1*)hDNoise[id]->GetFunction("gaus");
      derAve = fit1->GetParameter(1);
      derSigma = fit1->GetParameter(2);

      // find  peaks 
      peakType peakList;
      std::vector<Int_t> peakKind;
      TString dname = detList[id]->GetName();
      // for derivativePeaks, window in time is timeUnit*windowSize (ns) . timeUnit = 2
      Int_t windowSize=10;              
      // min, max width in time bins for simple peaks
      unsigned maxWidth = 100000; 
      unsigned minWidth = 10;
      if(dname.Contains( TString("1"))) peakList= derivativePeaks(vderiv, id, windowSize, derSigma, peakKind);
      else peakList= simplePeaks(ddigi, id, minWidth, maxWidth, rawSigma, peakKind);
      hPeakCount->Fill(id,peakList.size());

       // get WMARecursive baseline weigthts 
      int maxWeightWidth=0;
      std::vector<Double_t> weight = getBaselineWeights(ddigi.size(), peakList,maxWeightWidth);

      // set early weights to zero
      int wbin = 1000;
      for (unsigned i = 0; i < wbin; i++) weight[i]=0; 
      for (unsigned i = 0; i < weight.size(); i++) hEvWeight[id]->SetBinContent(i+1, weight[i]); 


      // get WMARecursive baseline 
      int NWindow  = (int)detList[id]->digi.size()/8200*400;  // size is 8200
      std::vector<Double_t> baselineDigi = getBaselineWMARecursive(0,ddigi, weight, NWindow);

      // zero basseline until after firt pulse
      if(peakList.size()>0) for(unsigned i=0; i<= std::get<1>(peakList[0]) ; ++i ) baselineDigi[i]=0;

      // subract base 
      std::vector<double> sdigi;
      for (unsigned i = 0; i < detList[id]->digi.size(); i++) {
          hBaselineWMA[id]->SetBinContent(i, baselineDigi[i]); 
          double sval = ddigi[i] - baselineDigi[i];
          sdigi.push_back(sval);
          hEvSignal[id]->SetBinContent(i+1,sval);
      }

      // extend peaks to baseline subtracted wave 
      extendPeaks(sdigi,peakList);
      /*** peak finding */
      for(unsigned isample = 0; isample < sdigi.size(); isample++) hSNoise[id]->Fill(sdigi[isample]);
      hSNoise[id]->Fit("gaus","0Q");
      TF1 *fit2 = (TF1*)hSNoise[id]->GetFunction("gaus");
      double sAve = 0;
      double sSigma = 0;
      if(fit2) {
        sAve=fit2->GetParameter(1);
        sSigma=fit2->GetParameter(2);
      } else 
        printf(" no SNoise fit event %lld det %i ",ievent,id);

      ntCal->Fill(ievent,id,rawAve,rawSigma,derAve,derSigma,sAve,sSigma);

      /** make hits **/
     double triggerTime=0;
     double firstCharge=0;
     int hitCount=0;
     hitMap  detHits = makeHits(peakList,peakKind,vderiv,sSigma,triggerTime,firstCharge);
     if (id==0) trigTime =  triggerTime;
     for (hitMapIter hitIter=detHits.begin(); hitIter!=detHits.end(); ++hitIter) {
        TDetHit  hiti = hitIter->second;
        Double_t hitQ = hiti.qsum;
        brun->detList[id]->hits.push_back(hiti);
        //Double_t phitQErr = phiti.qerr*timeUnit*1E9;
        Int_t width = hiti.lastBin - hiti.firstBin +1;
        for (int ibin=hiti.firstBin; ibin <= hiti.lastBin; ++ibin) 
          hEvPeaks[id]->SetBinContent(ibin+1, sdigi[ibin]);
        Double_t hitTime =  hiti.startTime*timeUnit*microSec-trigTime;
        Int_t istartBin =  hLife[id]->FindBin(hitTime);
        ntHit->Fill(id,detHits.size(), ++hitCount, istartBin, hiti.startTime*timeUnit*microSec, hitQ, width, hiti.qpeak, hiti.good, hiti.kind);
        hLife[id]->SetBinContent(istartBin, hLife[id]->GetBinContent(istartBin)+hitQ);
        hLife[id]->SetBinError(istartBin, sqrt(pow(hLife[id]->GetBinError(istartBin), 2)+pow(hiti.qerr, 2)));
     }
     if (ievent<10) plotWave(ievent, id);
    } // end of det loop
    // fill BRun
    brun->bevent->event =ievent;
    brun->bevent->trigTime =trigTime;
    brun->fill();

    if(dumpFile.is_open() ) brun->dumpEvent(dumpFile);
  } // loop over events
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
    detList[idet]->SetName(aBranch->GetName());
    ++idet;
  }

  for(unsigned id=0; id<NDET; ++id) hEvent[id]=NULL;
  cout << " number of entries in BTree is " <<  btree->GetEntries() << endl;
  btree->GetListOfBranches()->ls();
  return true;
}

// WMARecursive weights
std::vector<Double_t> anaRun::getBaselineWeights(unsigned arraySize, peakType peakList, Int_t& maxWeightWidth)
{
    std::vector<Double_t> weight(arraySize, 0);
    if (peakList.size() == 0) return weight;

    // construct the weights
    int alpha;
    int beta;
    maxWeightWidth=0;
    for (int ip=0; ip<int(peakList.size()); ++ip) {
        if (ip==0) beta = -1;
        else beta = std::get<1>(peakList[ip-1]);
        if (ip==peakList.size()) alpha = weight.size();
        else alpha = std::get<0>(peakList[ip]);
        int width = std::get<1>(peakList[ip])-std::get<0>(peakList[ip]);
        //printf(" \t peak %i at %i width %i\n",ip,std::get<0>(peakList[ip]),width);
        if (width>maxWeightWidth) maxWeightWidth=width;
        for (int jp=beta+1; jp<alpha; ++jp)  weight[jp]=alpha-beta-1;
    }
    return weight;
}
//Weighted Moving Average baseline from Zugec et al https://arxiv.org/pdf/1601.04512.pdf 
std::vector<Double_t> anaRun::getBaselineWMARecursive(Double_t ave, std::vector<Double_t> sig, std::vector<Double_t> weight,Int_t NWindow)
{
  // default baseline is zero
  std::vector<Double_t> baseline(sig.size(),ave);

  // starting values 
  Double_t KW=0;
  Double_t CW=0;
  Double_t SW=0;
  Double_t KS=0;
  Double_t CS=0;
  Double_t SS=0;

  // recursive sums
  Double_t kcos = TMath::Cos( TMath::Pi()/Double_t(NWindow));
  Double_t ksin = TMath::Sin( TMath::Pi()/Double_t(NWindow));
  //printf(" kcos %f ksin %f \n",kcos,ksin);
  Int_t resetBin = 100;

  for(int iw=0; iw<int(sig.size()); ++iw) {
    // reset
    if(iw%resetBin ==  0) {
      int low =  TMath::Max(0,iw-NWindow);
      int high = TMath::Min(iw+NWindow,int(sig.size())-1);
      KW=0;CW=0;SW=0;KS=0; CS=0;SS=0;
      for(int jp =low; jp < high; ++jp) {
        Double_t cosj = TMath::Cos(  Double_t(jp-iw)*TMath::Pi()/Double_t(NWindow) );
        Double_t sinj = TMath::Sin(  Double_t(jp-iw)*TMath::Pi()/Double_t(NWindow) );
        KW+= weight[jp];
        KS+= sig[jp]*weight[jp];
        CW+= weight[jp]*cosj;
        CS+= sig[jp]*weight[jp]*cosj;
        SW+= weight[jp]*sinj;
        SS+= sig[jp]*weight[jp]*sinj;
      }
    } else {
      //save previous iteraton
      int i1 = iw+NWindow;
      int i2 = iw-1-NWindow;
      Double_t CWlast=CW;
      Double_t CSlast=CS;

      // case sums 
      CW = kcos*CWlast+ksin*SW;
      CS = kcos*CSlast+ksin*SS;
      SW = kcos*SW-ksin*CWlast;
      SS = kcos*SS-ksin*CSlast;

      if(iw  <= NWindow && iw + NWindow <= sig.size() - 1){
        KW += weight[i1];
        KS += sig[i1]*weight[i1];
        CW += -1.*weight[i1];
        CS += -1.*sig[i1]*weight[i1];
      } else if(iw > NWindow && iw + NWindow <= sig.size() - 1){
        KW += -1.*weight[i2] + weight[i1];
        KS += -1.* sig[i2]*weight[i2] +  sig[i1]*weight[i1];
        CW +=  kcos*weight[i2] - weight[i1];
        CS +=  kcos*sig[i2]*weight[i2] - sig[i1]*weight[i1];
        SW += -1.*ksin*weight[i2] ;
        SS += -1.*ksin*sig[i2]*weight[i2] ;
      } else if(iw  > NWindow && iw + NWindow > sig.size() - 1){
        KW += -1.*weight[i2];
        KS +=  -1.*sig[i2]*weight[i2];
        CW +=  kcos*weight[i2];
        CS +=  kcos*sig[i2]*weight[i2];
        SW += -1.*ksin*weight[i2];
        SS += -1.*ksin*sig[i2]*weight[i2];
      }  // sum baseline 
    }
    if( TMath::Abs(KW+CW)>1E-12) baseline[iw]=(KS+CS)/(KW+CW);
  }
  return baseline;
}

peakType anaRun::simplePeaks(std::vector<Double_t> v, Int_t idet, unsigned minWidth, unsigned maxWidth, Double_t rms,std::vector<Int_t>& peakKind ) 
{
  peakType peakList;
  peakKind.clear();
  std::vector<unsigned> crossings;
  std::vector<unsigned> crossingBin;
  std::vector<double> crossingTime;
  unsigned vsize = v.size();
  Double_t cut = derivativeSigma[idet]*rms;
  //cout << " >>>> rms " << rms << " cut " << cut << endl;

  /* loop over bins finding v>cut.  continue until v<cut */
  for( unsigned ibin=1; ibin< vsize-1; ++ibin ) {
    Double_t vi = vsign[idet]*v[ibin];
    Double_t vmax=0;
    Double_t vj=0;
    unsigned jbin=ibin+1;
    if(vi>cut) {
      //if(idet==2) printf(" at  %i ibin %u vi %f\n",idet,ibin,vi);
      for(jbin=ibin+1; jbin< vsize; ++jbin ) {
        vj = vsign[idet]*v[jbin];
        if(vj>vmax) vmax=vj;
        //if(idet==2) printf(" \t   jbin %u vj %f \n",jbin,vj);
        if(vj<cut) break;
      }
      if(jbin-ibin>minWidth&&jbin-ibin<maxWidth) {
        peakList.push_back( std::make_pair(ibin,jbin-1) );
        peakKind.push_back(1);
        ntSimplePeaks->Fill(idet,rms,vmax,float(ibin),float(jbin-ibin));
        //printf(" found peak det %i vi %f vj %f start bin %u end  %u \n",idet,vi,vsign[idet]*v[jbin-1],ibin,jbin-ibin);
        ibin = jbin+1;
      }
    }
  }
  // return list
  return peakList;
}

peakType anaRun::derivativePeaks(std::vector<Double_t> v,  Int_t idet, Int_t nsum, Double_t rms,std::vector<Int_t>& peakKind ) 
{
  peakType peakList;
  peakKind.clear();
  std::vector<unsigned> crossings;
  std::vector<unsigned> crossingBin;
  std::vector<double> crossingTime;
  unsigned vsize = v.size();
  Double_t cut = derivativeSigma[idet]*rms;
  //cout << " >>>> rms " << rms << " cut " << cut << endl;
  Double_t ncut = -derivativeSigma[idet]*rms;
  // find all crossings
  for( unsigned ibin=1; ibin< vsize; ++ibin ) {
    Double_t u = double(ibin)*timeUnit*microSec;
    Double_t vi = vsign[idet]*v[ibin];
    Double_t vj = vsign[idet]*v[ibin-1];
    unsigned ctype=10;
    if( vi>cut && vj <cut ) {
      crossings.push_back(UPCROSS);
      ctype=UPCROSS;
      crossingBin.push_back(ibin);
      crossingTime.push_back(u);
    } else if( vi<cut && vj > cut ) {
      crossings.push_back(UPCROSS);
      ctype=UPCROSS;
      crossingBin.push_back(ibin);
      crossingTime.push_back(u);
    } else if( vi<ncut && vj > ncut ) {
      crossings.push_back(DOWNCROSS);
      ctype=DOWNCROSS;
      crossingBin.push_back(ibin);
      crossingTime.push_back(u);
    } else if( vi >ncut && vj < ncut ) {
      crossings.push_back(DOWNCROSS);
      ctype=DOWNCROSS;
      crossingBin.push_back(ibin);
      crossingTime.push_back(u);
    }
    //if (idet==0&&ctype<10)  printf("\t %u vj %f vi %f cut %f cross type %u \n", ibin, vj, vi, cut, ctype );
    //if (idet==0)  printf("\t %u vj %f vi %f cut %f  not cross  \n", ibin, vj, vi, cut );
  }

  if(crossings.size()<4) return peakList;

  // label found crossings, intially all false
  std::vector<bool> crossingFound;
  crossingFound.resize(crossings.size());
  for(unsigned jc=0; jc<crossings.size(); ++jc)  {
    //printf(" det %i crossing %i bin %i time %f type %i \n",idet,jc,crossingBin[jc],crossingTime[jc],crossings[jc]);
    crossingFound[jc]=false;
  }


  // parse crossings to make pairs 
  unsigned ip =0; 
  //printf("crossings %zu %E %E \n",crossings.size(),cut,ncut);
  while ( ip<= crossings.size() -4 ) {
    // UP UP DOWN DOWN 
    if(crossings[ip]==UPCROSS&&crossings[ip+1]==UPCROSS&&crossings[ip+2]==DOWNCROSS&&crossings[ip+3]==DOWNCROSS) {
      // check zero crossings
      unsigned nzero=0;
      for(unsigned ibin = crossingBin[ip+1]; ibin<crossingBin[ip+2]; ++ibin) if(v[ibin]>0&&v[ibin+1]<0) ++nzero;
      if(nzero<=3) {
        //printf(" peak %i time %f (%i %i %i %i ) nzero %u \n",ip,crossingTime[ip],crossings[ip],crossings[ip+1],crossings[ip+2],crossings[ip+3],nzero);
        peakList.push_back( std::make_pair(crossingBin[ip],crossingBin[ip+3]) );
        peakKind.push_back(0);
        ntDer->Fill(rms,v[crossingBin[ip]],double(crossingBin[ip+3]-crossingBin[ip]),double(0));//sigma:d0:step:dstep
        crossingFound[ip]=true;
        crossingFound[ip+1]=true;
        crossingFound[ip+2]=true;
        crossingFound[ip+3]=true;
        ip=ip+4;
      } else 
        ++ip;
    } else {
      ++ip;
    }
  }
  //printf("peaks with fours %zu  \n",peakList.size());

  // pick up UP UP, DOWN DOWN cases
  ip =0; 
  while ( ip<= crossings.size() -2 ) {  
    if(!crossingFound[ip]&&!crossingFound[ip+1]) {
      // UP UP
      if(crossings[ip]==UPCROSS&&crossings[ip+1]==UPCROSS) {
        //printf(" up up peak %i time %f (%i %i )\n",ip,crossingTime[ip],crossings[ip],crossings[ip+1]); 
        peakList.push_back( std::make_pair(crossingBin[ip],crossingBin[ip+1]) );
        peakKind.push_back(1);
        ntDer->Fill(rms,v[crossingBin[ip]],double(crossingBin[ip+1]-crossingBin[ip]),double(1));//sigma:d0:step:dstep
        // DOWN DOWN 
      } else if(crossings[ip]==DOWNCROSS&&crossings[ip+1]==DOWNCROSS) {
        //printf(" down down peak %i time %f (%i %i )\n",ip,crossingTime[ip],crossings[ip],crossings[ip+1]); 
        peakList.push_back( std::make_pair(crossingBin[ip],crossingBin[ip+1]) );
        peakKind.push_back(2);
        ntDer->Fill(rms,v[crossingBin[ip]],double(crossingBin[ip+1]-crossingBin[ip]),double(2));//sigma:d0:step:dstep
      }
    }
    ip=ip+2;
  }
 
  // return list
  return peakList;
}
void anaRun::getAverages(std::vector<Double_t> digi, Double_t& ave, Double_t& sigma, unsigned max) 
{
  if(max==0) max = digi.size();
  std::vector<double> vsort;
  for (unsigned is=0; is< max; ++is) vsort.push_back(digi[is]); 
  std::sort(vsort.begin(), vsort.end());
  ave = vsort[unsigned(0.5*double(vsort.size()))];
  sigma = TMath::Abs(vsort[0.659*vsort.size()]-ave);
}

void anaRun::getAverage(std::vector<Double_t> digi, Double_t& ave, Double_t& sigma, unsigned max) 
{
  
  if(max==0) max = digi.size();
  double sum=0;
  double sum2=0;
  for (unsigned is=0; is< max; ++is) {
    sum += digi[is];
    sum2 += pow(digi[is],2.);
  }
  ave = sum/double(max);
  sigma = sqrt( sum2/double(max) - pow(ave,2.));
}
//implimented as in zugec et al, arXiv:1601.04512v1
std::vector<Double_t> anaRun::differentiate(std::vector<Double_t> v, unsigned diffStep)
{
  std::vector<Double_t> d;
  unsigned nsamples = v.size();
  Double_t sump=0;
  Double_t summ=0;
  d.push_back(0); // first entry is zero
  for(unsigned i=1; i<nsamples; ++i) {
    unsigned i2 = 2*i;
    unsigned max = TMath::Min( diffStep, i);
    max = TMath::Min(max, nsamples - 1 -i);
    // beginning, middle, end cases
    if(i<=diffStep && i2 <= nsamples -1 ) {
      sump = sump - v[i] + v[i2-1] + v[i2];
      summ = summ + v[i-1];
    }
    else if(i>diffStep && i+diffStep <= nsamples -1){
      sump = sump - v[i] + v[i+diffStep];
      summ = summ + v[i-1]-v[i-1-diffStep];
    }
    else if(i+diffStep> nsamples-1 && i2 > nsamples-1) {
      sump = sump - v[i];
      summ = summ + v[i-1] - v[i2-nsamples-1] - v[i2-nsamples];
    }
    d.push_back(sump-summ);
  }
  return d;
}

hitMap anaRun::makeHits(peakType peakList, std::vector<Int_t> peakKind, std::vector<Double_t> ddigi, Double_t sigma, Double_t& triggerTime, Double_t& firstCharge)
{
    triggerTime=1E9;
    firstCharge=0;
    hitMap detHits;
    if (peakList.size()<1) return detHits;
    Double_t qmax=0;

    unsigned minLength=5;
    for (unsigned ip=0; ip<peakList.size(); ++ip) {
        unsigned klow  = std::get<0>(peakList[ip]);
        unsigned khigh = std::get<1>(peakList[ip]);
        //printf(" hit  %u (%u,%u) kind %i length %u \n",ip,klow,khigh,peakKind[ip],khigh-klow+1);
        hHitLength->Fill(khigh-klow+1);
        if (khigh-klow+1<minLength) {
            continue;
        }
        Double_t qhit=0;
        UInt_t peakt=0;
        Double_t qpeak=0;
        Double_t qsum = 0;
        for (unsigned k=klow; k<khigh; ++k) {
            double qdigik = -1.*ddigi[k];
            qsum += qdigik;
            if (qdigik>qpeak) {
                peakt=k;
                qpeak = qdigik;
            }
        }

        TDetHit dhit;
        dhit.peakBin=Int_t(peakt);
        dhit.qsum=qsum;
        dhit.qpeak=qpeak;
        dhit.firstBin = klow;
        dhit.lastBin = khigh;
        dhit.peakMaxTime=peakt;
        dhit.peakt=peakt;
        dhit.startTime=klow;
        dhit.peakWidth=khigh-klow;
        // this is N= q/qnorm and delta q = root(n)*qnorm;
        dhit.qerr = sqrt(pow(sigma*Double_t(dhit.peakWidth), 2)+ qnorm*qsum);
        dhit.kind = peakKind[ip];

        // just use the biggest pulse 
        if (qsum>qmax) {
            qmax=qsum;
            triggerTime=dhit.startTime*timeUnit*microSec;
            firstCharge = qsum;
        }
        Double_t hitTime = dhit.startTime*timeUnit*microSec;
        //printf(" insert  %lu (%u,%u)  \n",detHits.size(),dhit.firstBin,dhit.lastBin);
        //for (unsigned k=klow; k<khigh; ++k) printf(" \t %u %f ; ", k, ddigi[k]);
        //cout << endl;
        detHits.insert(std::pair<Double_t, TDetHit>(hitTime, dhit));
        hPeakNWidth->Fill(dhit.lastBin-dhit.firstBin+1);
    }

    // first time, charge from map
    /*
    hitMapIter hitIter;
    hitIter=detHits.begin();
    TDetHit dhit0 = hitIter->second;
    triggerTime = dhit0.startTime*microSec;
    firstCharge = dhit0.qsum;
    */
    return  detHits;
}

// extend peaks to zero of waveform  
void anaRun::extendPeaks(std::vector<Double_t> v, peakType &peakList ) {
  for(unsigned ip=0; ip<peakList.size(); ++ip)  {
    // high direction
    unsigned high = std::get<1>(peakList[ip]);
    unsigned next = v.size();
    if(ip<peakList.size()-1) next = std::get<1>(peakList[ip+1]);
    for(unsigned kp= high; kp < next ; ++kp ) {
      if( v[kp] <0) break;
      std::get<1>(peakList[ip])=kp;
    }

    // low direction
    unsigned low = std::get<0>(peakList[ip]);
    unsigned prev = 1;
    if(ip>0) prev = TMath::Max(prev,std::get<0>(peakList[ip-1]));
    for(unsigned kp= low; kp > prev ; --kp ) {
      if( v[kp] <0 ) break;
      std::get<0>(peakList[ip])=kp;
    }
  }
}



void anaRun::plotWave(Long64_t ievent, unsigned idet) {
    evDir->cd();
    TString histName;
    TString detName = detList[idet]->GetName();

    histName.Form("EvWave%lli_DET%1i_%s", ievent,idet, detName.Data());
    TH1D* hwave = (TH1D*)hEvWave[idet]->Clone(histName);

    histName.Form("EvBase%lli_DET%1i_%s", ievent,idet,detName.Data());
    TH1D* hbase = (TH1D*)hBaselineWMA[idet]->Clone(histName);

    histName.Form("EvDWave%lli_DET%1i_%s", ievent,idet,detName.Data());
    TH1D* hdwave = (TH1D*)hEvDWave[idet]->Clone(histName);

    histName.Form("EvSignal%lli_DET%1i_%s", ievent,idet,detName.Data());
    TH1D* hsignal = (TH1D*)hEvSignal[idet]->Clone(histName);

    histName.Form("EvPeaks%lli_DET%1i_%s", ievent,idet,detName.Data());
    TH1D* hpeaks = (TH1D*)hEvPeaks[idet]->Clone(histName);

    histName.Form("EvWeight%lli_DET%1i_%s", ievent,idet,detName.Data());
    TH1D* hweight = (TH1D*)hEvWeight[idet]->Clone(histName);

    fout->cd();
}
