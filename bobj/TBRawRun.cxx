#include "TBRawRun.hxx"
ClassImp(TBRawRun)

TBRawRun::TBRawRun(TString runName ): TNamed(runName,runName)
{
  btree = new TTree("BTree"," bacon data " );
  det0.SetName(Form("det%i",6));
  det1.SetName(Form("det%i",4));
  det2.SetName(Form("det%i",1));
  det3.SetName(Form("det%i",2));
  detList.push_back(&det0);
  detList.push_back(&det1);
  detList.push_back(&det2);
  detList.push_back(&det3);

  for(unsigned i=0; i<NDET; ++i ) {
    btree->Branch(detList[i]->GetName(),detList[i]);
  }
  cout << " TBRawRun tree " << btree->GetName() << endl;
  btree->GetListOfBranches()->ls();
  clear();
}


//TBRawRun::~TBRawRun(){}

void TBRawRun::clear()
{
  run=0;
  timeUnit=2; // ns, default.
  detListClear();
}

void TBRawRun::detListClear() 
{
  for(unsigned i=0; i<3; ++i ) detList[i]->clear(); 
}


