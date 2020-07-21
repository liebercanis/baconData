#include "TBRawEvent.hxx"
ClassImp(TBRawEvent)

TBRawEvent::TBRawEvent(TString detName): TNamed(detName,detName)
{
  description=TString("none");
  cout << " new TBRawEvent "<< detName << "  "  << this->GetName() << " " << description << endl;
  clear();
}

//TBRawEvent::~TBRawEvent(){}

void TBRawEvent::clear()
{
  event=0;
  digi.clear();	 
}

