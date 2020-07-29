#include <fstream>
#include "TString.h"
#include "TObjString.h"
#include "TSystemDirectory.h"
#include "TFile.h"
#include "TBRawRun.hxx"

TBRawRun *brun;
vector<ifstream*> streams;

using namespace std;

void readEvent() 
{
  char line[256];

  for(unsigned is =0; is<streams.size(); ++is ) {
    TBRawEvent* rawEv = brun->detList[is];
    //cout << is << " for "  << rawEv->description  << endl;
    while( streams[is]->good() ) {
      streams[is]->getline(line,256);
      TString tline(line);
      TObjArray* tokenArray = tline.Tokenize(' ');
      if(tokenArray->GetEntries()>1) {
        TObjString *sv = (TObjString*) tokenArray->At(0);
        if(sv->GetString()==TString("DC") && rawEv->digi.size()>0 ) {
          //cout << " end of event " << rawEv->event << "  " << rawEv->GetName() << endl;
          break;
        }
      } else {
        TObjString *sv = (TObjString*) tokenArray->At(0);
        if(sv==NULL) continue;
        rawEv->event =  brun->btree->GetEntries() ;
        rawEv->digi.push_back(sv->GetString().Atof());
      }
    }
    //printf (" readEvent for %s digis  %lu \n",rawEv->GetName(),rawEv->digi.size());
  }
  brun->fill();
  brun->detListClear();
}


void openFiles() 
{
  TString dirName(Form("data/%s",brun->GetName()));
  cout << dirName << endl;
  TSystemDirectory dir("dataFiles",dirName);
  TList *files = dir.GetListOfFiles();
  files->ls();

  TIter next(files);
  TSystemFile *file;
  while( (file = (TSystemFile*) next()) ) {
    string name = string(file->GetName());
    string exten  = name.substr( name.find_last_of(".")+1 );
    if(exten!=string("txt")) continue;
    string tag = name.substr( 0, name.find_last_of(".") );
    cout << " file is " << name << " tag  " << tag << endl;
    string fullName =  string( dirName.Data())  + string("/")+name;
    ifstream* in = new ifstream(fullName,std::ios::in);
    if(in->is_open()) streams.push_back(in);
    brun->detList[streams.size()-1]->description=TString(tag.c_str());

  }
  printf("opened %ld files \n defined detectors: \n",streams.size());
  for(unsigned  i=0; i< brun->detList.size(); ++i)  brun->detList[i]->print();

}


void closeFiles() 
{
  for(unsigned i=0; i<streams.size(); ++i ) streams[i]->close();
}


void readRaw(TString runName = "pmtCh1")
{
  TFile* fout = new TFile(Form("%s.root",runName.Data()),"recreate");
  brun = new TBRawRun(runName);
  brun->timeUnit=2; // ns
  cout << " run name " << brun->GetName() << endl;
  openFiles();
  while( streams[0]->good()) { 
    readEvent(); 
    if(int(brun->btree->GetEntries())%100==0) printf("... %lld \n",brun->btree->GetEntries() );
  }
  closeFiles();
  fout->ls();
  fout->Write();
}
