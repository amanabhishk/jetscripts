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
  pythia.readFile("pythia_Zjet.cmnd");
  pythia.init();

  std::stringstream outputFilename("");
  outputFilename << nEvent <<"events_Zjet.root";

  TFile outFile(outputFilename.str().c_str(), "NEW");
  
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

  int state, count, temp,  index, identity;
  vector<int> leptonListFinal;

  std::clock_t start = std::clock();

  /****************************************END OF SET-UP**************************************************/
  for (int iEvent = 0; iEvent != nEvent; ++iEvent) 
  {
    if (!pythia.next()) continue;
    leptonListFinal.resize(0);
    weight = pythia.info.weight();
    count = -1;

    for(int t=0; t != event.size(); ++t)
    {
      assert(count<size);
      Particle& p = pythia.event.at(t);
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
          fillTree(p, count, id, pT, eta, phi, m);
        }
      }

      if(state == 23 && identity != 13)
      {
        ++count;
        status[count] = 3;
        fillTree(p, count, id, pT, eta, phi, m);
      }

      hadronic_definition_status_codes(event, t, count, status, id, pT, eta, phi, m);
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
        if(temp==1) index = event[index].daughter1();
        if(temp==2) index = event[index].daughter2();
      }
      
      else
      {
        ++count;
        status[count] = 1;
      }

      fillTree(event.at(index), count, id, pT, eta, phi, m);
    }
    //assert(temp==2);

    n = UShort_t(count+1);
    tree->Fill();
  }

  tree->AutoSave("Overwrite");
  outFile.Close();
  printTime((std::clock()-start),nEvent);
  cout<<"TTree is saved in "<<outputFilename.str().c_str()<<endl;
  return 0;
}
