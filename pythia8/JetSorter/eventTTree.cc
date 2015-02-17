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
  outputFilename << nEvent <<"events.root";

  //ROOT TTree setup  
  TFile outFile(outputFilename.str().c_str(), "recreate"); //output file. change the name in physicsDef.cc too
  
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
  
  Float_t rad;
  tree->Branch("rad", &rad, "rad/F");
  TH1D * fsr = new TH1D("fsr","fsr",150,0,1.5);
 
  int state, count, daughter1, daughter2, temp,  index;
  Float_t num, den;
  vector<int> leptonListFinal;

  /****************************************END OF SET-UP**************************************************/
  for (int iEvent = 0; iEvent != nEvent; ++iEvent) 
  {
    if (!pythia.next()) continue;
    
    leptonListFinal.resize(0);
    // weight = info.weight();
    // count = -1;
    num = 0; den = 0;
    // cout<<endl;
    temp = 0;
    for(int t=0; t != event.size(); ++t)
    {
      if(abs(event[t].id()) == 13)
      {
        if(event[t].isFinal()) 
        {
          // cout<<t<<" ";
          leptonListFinal.push_back(t);
          // temp++;
        }
      }
    }

    for(unsigned int x = 0; x != leptonListFinal.size(); ++x)
    {
      index = leptonListFinal[x];
      // cout<<index<<": ";
      while(abs(event[index].id())==13)
      {
        index = event[index].mother1();
        // cout<<index<<" ";
      }
      if(event[index].id()==23)
      {
        // cout<<"(Z0) ";
        num += event[leptonListFinal[x]].pT();
        temp++;
        if(temp == 1) daughter1 = event[index].daughter1();
        if(temp == 2) daughter2 = event[index].daughter2();

        while(abs(event[index].id())==23)
        {
          index = event[index].mother1();
          // cout<<index<<" ";
        }
      }
      // else cout<<"Dead end.";
      // cout<<endl;
    }
    assert(temp==2);
    den = event[daughter1].pT()+event[daughter2].pT();
    weight = num/den;
    if(weight>1) cout<<iEvent<<" ";
    tree->Fill();
    fsr->Fill(num/den);
  }

  fsr->Write();

  // tree->Print();
  tree->AutoSave("Overwrite"); 
  outFile.Close();
  // cout<<"Done.\n";
  // cout<<"TTree is saved in "<<outputFilename<<endl;
  return 0;
}
