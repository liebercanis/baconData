#include "TBRun.hxx"
ClassImp(TBRun)

TBRun::TBRun(TString runName ): TNamed(runName,runName)
{
  btree = new TTree("BRunTree"," bacon data " );
  bevent = new TBEvent(runName);
  btree->Branch("bev",&bevent);
  det0.SetName(Form("wave%i",1));
  det1.SetName(Form("wave%i",4));
  det2.SetName(Form("wave%i",2));
  det3.SetName(Form("wave%i",6));
  detList.push_back(&det0);
  detList.push_back(&det1);
  detList.push_back(&det2);
  detList.push_back(&det3);

  for(unsigned i=0; i<NDET; ++i ) {
    cout << " btree adding branch " << detList[i]->GetName() << endl;
    btree->Branch(detList[i]->GetName(),detList[i]);
  }
  cout << " TBRun tree " << btree->GetName() << endl;
  btree->GetListOfBranches()->ls();
  clear();
}


//TBRun::~TBRun(){}

void TBRun::clear()
{
  bevent->clear();
  detListClear();
}

void TBRun::detListClear() 
{
  for(unsigned i=0; i<NDET; ++i ) detList[i]->clear(); 
}


