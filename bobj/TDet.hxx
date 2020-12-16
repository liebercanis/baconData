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
    Double_t qPrompt;
    Double_t qSum;
    Double_t hitPrompt;
    Double_t hitSum;

    std::vector<TDetHit> hits; 

    unsigned nhits() { return hits.size(); }

    void clear()
    {
      event=0;
      flags=0;
      nspe=0;
      qPrompt=0;
      qSum=0;
      hitPrompt=0;
      hitSum=0;
      hits.clear();
    }

    ClassDef(TDet,2)
};
#endif

