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
  //Info& info = pythia.info;
  pythia.readFile("pythia_dijet.cmnd");
  pythia.init();
  pythia.settings.listChanged();

  std::stringstream outputFilename("");
  outputFilename << nEvent <<"events_dijet.root";

  //ROOT TTree setup  
  TFile outFile(outputFilename.str().c_str(), "NEW"); //output file. change the name in physicsDef.cc too
  
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
 
  int state, count;
  
  /****************************************END OF SET-UP**************************************************/
  for (int iEvent = 0; iEvent != nEvent; ++iEvent) 
  {
    if (!pythia.next()) continue;
    
    weight = pythia.info.weight();
    count = -1;

    for(int t=0; t != event.size(); ++t)
    {
      assert(count<size);
      //state = abs(event[t].status());      
      
      if(event[t].status() == -23) 
      {
        ++count;
        status[count] = 3;
      }
      
      else if(event[t].isFinal())
      {
        ++count;
        status[count] = 1;
      }

      else continue;

      id[count] = event[t].id();
      pT[count] = event[t].pT();
      eta[count] = event[t].eta();
      phi[count] = event[t].phi();
      m[count] = event[t].m();
    }

    n = UShort_t(count+1);
    tree->Fill();
  
  }

  tree->Print();
  tree->AutoSave("Overwrite"); 
  outFile.Close();
  cout<<"Done.\n";
  cout<<"TTree is saved in "<<outputFilename.str().c_str()<<endl;
  return 0;
}
