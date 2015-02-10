// This class sorts pythia8 jets with the fastjet algorithm. See READMEi_ScriptInfo for further details.

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
  if (argc > 1){
    nEvent = atoi(argv[1]);
  }


  // Pythia setup
  Pythia pythia;
  Event& event = pythia.event;
  Info& info = pythia.info;

  pythia.readString("HardQCD:all = on");
  pythia.readString("PhaseSpace:pTHatMin = 30.");
  pythia.readString("Next:numberShowInfo = 0");
  pythia.readString("Next:numberShowProcess = 0");
  pythia.readString("Next:numberShowEvent = 0");
  pythia.readString("PhaseSpace:bias2Selection = on");
  pythia.readString("PhaseSpace:bias2SelectionPow = 4.5");
  pythia.readString("Beams:eCM = 8000.");

  pythia.init();
  pythia.settings.listChanged();

  //ROOT TTree setup  
  UShort_t size = 5000;
  
  UShort_t n;
  Int_t id[size];
  Float_t weight, pT[size], eta[size], phi[size], m[size];

  TFile outFile("output1.root", "RECREATE");

  TTree *tree = new TTree("Events","TTree of pythia events");
  tree->Branch("n", &n, "n/s");
  tree->Branch("weight", &weight, "weight/F");
  tree->Branch("id", id, "id[n]/I");
  tree->Branch("pT", pT, "pT[n]/F");
  tree->Branch("eta", eta, "eta[n]/F");
  tree->Branch("phi", phi, "phi[n]/F");
  tree->Branch("m", m, "m[n]/F");
 
  /****************************************END OF SET-UP**************************************************/
  for (int iEvent = 0; iEvent != nEvent; ++iEvent) 
  {
    if (!pythia.next()) continue;
    
    n = event.size();
    weight = info.weight();

    for(int t=0; t != event.size(); ++t)
    {
      assert(t<size);

      id[t] = event[t].id();
      pT[t] = event[t].pT();
      eta[t] = event[t].eta();
      phi[t] = event[t].phi();
      m[t] = event[t].m();
    }
    tree->Fill();
  }

  tree->AutoSave("Overwrite"); 
  outFile.Close();
  return 0;
}
