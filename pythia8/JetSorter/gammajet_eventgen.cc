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
  pythia.readFile("pythia_gammajet.cmnd");
  pythia.init();

  std::stringstream outputFilename("");
  outputFilename << nEvent <<"events_gammaJet.root";

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

  int count, gammaIndex = -1;
  bool skip;
  std::clock_t start = std::clock();

  /****************************************END OF SET-UP**************************************************/
  for (int iEvent = 0; iEvent != nEvent; ++iEvent) 
  {
    skip = false;
    if (!pythia.next()) continue;
    weight = pythia.info.weight();
    count = -1;

    for(int t = 0; t != event.size(); ++t) //find gamma from the hard process
    {
      if(event[t].status() == -23 && event[t].id() == 22)
      {
        gammaIndex = t;
        break;
      }
    }

    //assert(gammaIndex != -1);

    while(!event[gammaIndex].isFinal()) //find exclude all the complicated cases with pair production
    {                                   
      if(event[gammaIndex].daughter1() != event[gammaIndex].daughter2()) 
      {
        skip = true;
        break;
      }
      
      else gammaIndex = event[gammaIndex].daughter1();
    }

    if(skip) continue;

    ++count;
    status[count] = 2;
    fillTree(event.at(gammaIndex), count, id, pT, eta, phi, m);
    // id[count] = event[gammaIndex].id();
    // pT[count] = event[gammaIndex].pT();
    // eta[count] = event[gammaIndex].eta();
    // phi[count] = event[gammaIndex].phi();
    // m[count] = event[gammaIndex].m();


    for(int t=0; t != event.size(); ++t)
    {
      assert(count<size);
      Particle& p = pythia.event.at(t);

      if(event[t].isFinal() && t != gammaIndex)
      {
        ++count;
        status[count] = 1;
        fillTree(p, count, id, pT, eta, phi, m);
      }

      if(event[t].status() == -23 && abs(event[t].id()) != 22)
      {
        ++count;
        status[count] = 3;
        fillTree(p, count, id, pT, eta, phi, m); 
      }

      hadronic_definition_status_codes(event, t, count, status, id, pT, eta, phi, m);
    }

    n = UShort_t(count+1);
    tree->Fill();
  }

  tree->AutoSave("Overwrite");
  outFile.Close();
  printTime((std::clock()-start),nEvent);
  cout<<"TTree is saved in "<<outputFilename.str().c_str()<<endl;
  return 0;
}
