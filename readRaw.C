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
  brun->detListClear();
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
    //printf (" readEvent %i for %s digi.size  %lu \n",is,rawEv->GetName(),rawEv->digi.size());
  }
  brun->fill();
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
    TBranch *branch = brun->btree->GetBranch(brun->detList[streams.size()-1]->GetName());
    branch->SetName(tag.c_str());
    branch->SetTitle(tag.c_str());
    brun->detList[streams.size()-1]->SetName(tag.c_str());
    brun->detList[streams.size()-1]->description=TString(tag.c_str());

  }
  printf("opened %ld files \n defined detectors: \n",streams.size());

}


void closeFiles() 
{
  for(unsigned i=0; i<streams.size(); ++i ) streams[i]->close();
}


void readRaw( Long64_t maxRead=0, TString runName = "7_27_2020" )
{
  TFile* fout = new TFile(Form("rootData/%s.root",runName.Data()),"recreate");
  brun = new TBRawRun(runName);
  brun->timeUnit=2; // ns
  cout << " run name " << brun->GetName() << endl;
  openFiles();
  while( streams[0]->good()) { 
    readEvent(); 
    if(int(brun->btree->GetEntries())%100==0) printf("... %lld \n",brun->btree->GetEntries() );
    if(maxRead>0 && brun->btree->GetEntries()> maxRead ) break;
  }
  printf(" finished  %lld \n",brun->btree->GetEntries() );
  closeFiles();
  brun->btree->GetListOfBranches()->ls();
  fout->ls();
  fout->Write();
}
