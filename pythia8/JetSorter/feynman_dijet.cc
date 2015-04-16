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



int vertexCheck(int x, int y)
{
  if( x == 21 && y == 21) return 1;
  else if( x == 21 && abs(y) < 7 && y > 0) return 2;
  else if( x == 21 && abs(y) < 7 && y < 0) return 3;
  else if( abs(x) < 7 && x > 0 && abs(y) < 7 && y < 0) return 4;
  else if( abs(x) < 7 && x > 0 && abs(y) < 7 && y > 0) return 5;
  else if( abs(x) < 7 && x < 0 && abs(y) < 7 && y < 0) return 6;
  else vertexCheck(y,x);
}


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
  pythia.readFile("pythia_dijet.cmnd");
  pythia.init();

  std::stringstream outputFilename("");
  outputFilename << nEvent <<"dijet_feynman.root";

  TFile outFile(outputFilename.str().c_str(), "RECREATE");
  
  UShort_t size = 2000;     //CHECK: expected maximum number of particles to be stored for each event. Will lead to SEGFAULT if small.
  
  int inType, outType;               //exact number of particles stored for each event 
  float pT, weight;

  TTree *tree = new TTree("Events","TTree of Feynman diagrams");
  tree->Branch("inType", &inType, "inType/I");
  tree->Branch("outType", &outType, "outType/I");
  tree->Branch("weight", &weight, "weight/F");
  tree->Branch("pT1", &pT, "pT1/F");
  
  int lim, count;
  vector<int> out;
  vector<int> in;
  bool back2back;

  std::clock_t start = std::clock();
  
  /****************************************END OF SET-UP**************************************************/
  for (int iEvent = 0; iEvent != nEvent; ++iEvent) 
  {
    //cout<<iEvent<<endl;
    if (!pythia.next()) continue;

    out.resize(0);
    in.resize(0);
    
    for(int t=0; t != event.size(); ++t)
    {
      assert(count<size);
      Particle& p = pythia.event.at(t);

      if(event[t].status() == -23) 
      {
        out.push_back(t);
        //fillTree(p, count, id, pT, eta, phi, m);
      }
      
      //hadronic_definition_status_codes(event, t, count, status, id, pT, eta, phi, m);
      else if(event[t].status() == -21)
      {
        in.push_back(t);
        //++count;
        //status[count] = 2;
        //fillTree(p, count, id, pT, eta, phi, m);
      }
    }

    assert(out.size()==2);
    assert(in.size()==2);
    
    weight = pythia.info.weight();
    inType = vertexCheck(event[in[0]].id(),event[in[1]].id());
    outType = vertexCheck(event[out[0]].id(),event[out[1]].id());
    
    back2back = deltaPhi(event[out[0]].phi(),event[out[1]].phi())>2.8;
    
    if(back2back && fabs(event[out[0]].eta()<etaMax));
    {
      
      pT = event[out[0]].pT();
      tree->Fill();
    }

    //back2back = deltaPhi(event[out[0]].phi(),event[out[1]].phi())>2.8;
    if(back2back  && fabs(event[out[1]].eta()<etaMax));
    {
      pT = event[out[1]].pT();
      tree->Fill();
    }
  }

  tree->AutoSave("Overwrite"); 
  outFile.Close();
  cout<<"TTree is saved in "<<outputFilename.str().c_str()<<endl;
  return 0;
}
