/**
** MG, July 2020 
**/
#ifndef TBRAWRUN_DEFINED
#define TBRAWRUN_DEFINED
#include <iostream>
#include <string>
#include <map>
#include <TNamed.h>
#include <TTree.h>
#include "TBRawEvent.hxx"

using namespace std;

// class to store info for the run

class TBRawRun: public TNamed {
  public:
    TBRawRun(TString runName = "run0");
    //~TBRawRun();

    void clear();
    enum {NDET=4};
    Int_t run;
    TTree *btree;
    int timeUnit; //ns
    TBRawEvent det3;
    TBRawEvent det2;
    TBRawEvent det1;
    TBRawEvent det0;

    vector<TBRawEvent*> detList;
    void detListClear(); 
    void fill() {
      btree->Fill();
    }

    void print() {
      printf(" %s  entries %lld \n",this->GetName(), btree->GetEntries());
      for(unsigned i=0; i<detList.size(); ++i) printf(" %s ev %i %lu \n", detList[i]->description.Data(), detList[i]->event, detList[i]->digi.size() );
      btree->GetListOfBranches()->ls();
    }

    ClassDef(TBRawRun,1)
};
#endif

