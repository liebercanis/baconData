/**
** MG, July 2020 
**/
#ifndef TBRAWEVENT_DEFINED
#define TBRAWEVENT_DEFINED
#include <iostream>
#include <string>
#include <TNamed.h>
#include <vector>

using namespace std;

// class to store info for the event

class TBRawEvent: public TNamed {
	public:
    TBRawEvent(TString detName = "det0");
    //		~TBRawEvent();

		void clear();
		// data elements
    TString description;
    Int_t   event;
		std::vector<Double_t> digi;		 
		ClassDef(TBRawEvent,1)
};
#endif

