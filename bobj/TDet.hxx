/**
** MG, March 2020 
**/
#ifndef TDET_DEFINED
#define TDET_DEFINED
#include <iostream>
#include <string>
#include <TNamed.h>
#include <vector>
#include "TDetHit.hxx"

using namespace std;

// class to store info for the event

class TDet: public TNamed {
	public:
		TDet();
    //TDet::~TDet(){}

		// data elements
    Long64_t event;
    Int_t    flags;
    Int_t    nspe;
    Double_t energy;
    Double_t qsum;
    std::vector<TDetHit> hits; 

    unsigned nhits() { return hits.size(); }

    void clear()
    {
      event=0;
      flags=0;
      nspe=0;
      energy=0;
      qsum=0;
      hits.clear();
    }

    ClassDef(TDet,1)
};
#endif

