//Generates a TTree with particles of interest in each event, and stores it in a .root file

// Stdlib header file for input and output.
#include <iostream>
#include <stdlib.h>
#include <vector>
#include <string>
#include <sstream>

// Header file to access Pythia 8 program elements.
#include "Pythia8/Pythia.h"
#include "TF1.h"

#include <cmath>
#include <ctime>
// FastJet interface
#include "Pythia8/FastJet3.h"

// ROOT, for histogramming.
#include "TROOT.h"
#include "TH1.h"
#include "TProfile.h"
#include "THStack.h"
#include "TStyle.h"

// ROOT, for interactive graphics.
#include "TVirtualPad.h"
#include "TApplication.h"

// ROOT, for saving file.
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"

// tdrStyle
#include "tdrstyle_mod1.C"
// scripts
#include "jetsorter_auxiliary.h"

#include "TTree.h"

using namespace Pythia8;

int main(int argc, char* argv[]) 
{

  int  nEvent = 100;
  if (argc > 1)
  {
    nEvent = atoi(argv[1]);
  }

  // Pythia setup
  Pythia pythia;
  Event& event = pythia.event;
  Info& info = pythia.info;
  pythia.readFile("pythiaSettings.cmnd");
  pythia.init();
  pythia.settings.listChanged();

  std::stringstream outputFilename("");
  outputFilename << nEvent <<"eventsZjet.root";

  //ROOT TTree setup  
  TFile outFile(outputFilename.str().c_str(), "RECREATE"); //output file. change the name in physicsDef.cc too
  
  UShort_t size = 2000;     //CHECK: expected maximum number of particles to be stored for each event. Will lead to SEGFAULT if small.
  
  UShort_t n;               //exact number of particles stored for each event 
  Int_t id[size];
  Float_t weight, pT[size], eta[size], phi[size], m[size];
  UChar_t status[size];


  TTree *tree = new TTree("Events","TTree of pythia events");
  tree->Branch("n", &n, "n/s");
  tree->Branch("weight", &weight, "weight/F");
  tree->Branch("id", id, "id[n]/I");
  tree->Branch("status", status, "status[n]/b");
  tree->Branch("pT", pT, "pT[n]/F");
  tree->Branch("eta", eta, "eta[n]/F");
  tree->Branch("phi", phi, "phi[n]/F");
  tree->Branch("m", m, "m[n]/F");

  TH1D* type = new TH1D("type","type",26,0,26);
  
  // Float_t rad;
  // tree->Branch("rad", &rad, "rad/F");
  // TH1D * fsr = new TH1D("fsr","fsr",300,0,1.5);
 
  int state, count, temp,  index, identity;
  // Float_t num, den;
  vector<int> leptonListFinal;

  /****************************************END OF SET-UP**************************************************/
  for (int iEvent = 0; iEvent != nEvent; ++iEvent) 
  {
    // makeTree(event);
    if (!pythia.next()) continue;
    // cout<<"----------\n";

    leptonListFinal.resize(0);
    weight = info.weight();
    count = -1;

    for(int t=0; t != event.size(); ++t)
    {
      assert(count<size);
      state = abs(event[t].status());
      identity = abs(event[t].id());

      if(event[t].isFinal())
      {
        if(identity == 13) 
        {
          leptonListFinal.push_back(t);
          continue;
        }
        else
        {
          ++count;
          status[count] = 1;
        }
      }

      else if(state == 23 && identity != 13)
      {
        ++count;
        status[count] = 3;
        type->Fill(identity);
      }

      else continue;

      id[count] = event[t].id();
      pT[count] = event[t].pT();
      eta[count] = event[t].eta();
      phi[count] = event[t].phi();
      m[count] = event[t].m();
    }

    temp = 0;

    for(unsigned int x = 0; x != leptonListFinal.size(); ++x)
    {
      index = leptonListFinal[x];
      while(abs(event[index].id())==13)
      {
        index = event[index].mother1();
      }
      
      if(event[index].id()==23)
      {
        
        ++count;
        status[count] = 2;
        ++temp;
        if(temp==1) index = event[index].daughter1();//, cout<<index<<" ";
        if(temp==2) index = event[index].daughter2();//, cout<<index<<" ";
        // assert(temp<2);
        // cout<<event[index].phi()<<" "<<event[index].eta()<<endl;
        // if(temp == 1) daughter1 = event[index].daughter1(), a = leptonListFinal[x];
        // if(temp == 2) daughter2 = event[index].daughter2(), b = leptonListFinal[x];

        // while(abs(event[index].id())==23)
        // {
          // index = event[index].mother1();
          // cout<<index<<" ";
      }
      
      else
      {
        ++count;
        status[count] = 1;
      }
      // cout<<endl;
      

      id[count] = event[index].id();
      pT[count] = event[index].pT();
      eta[count] = event[index].eta();
      phi[count] = event[index].phi();
      m[count] = event[index].m();
      // else cout<<"Dead end.";
      // cout<<endl;
    }
    // cout<<endl;
    assert(temp==2);

    // den = (event[daughter1].p()+event[daughter2].p()).pT();
    // num = (event[a].p()+event[b].p()).pT();
    
    // cout<<num<<" "<<den<<endl;
    // if(weight>1) cout<<weight<<" ";
    // weight = num/den;
    n = UShort_t(count+1);
    tree->Fill();
    // fsr->Fill(num/den);
  }

  // fsr->Write();

  // tree->Print();
  type->Write();
  tree->AutoSave("Overwrite");
  outFile.Close();
  // makeTree();
  // cout<<"Done.\n";
  // cout<<"TTree is saved in "<<outputFilename<<endl;
  return 0;
}
