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
  TH1D * fsr = new TH1D("fsr","fsr",10,0,2);
 
  int state, count, daughter1, daughter2, temp;
  Float_t num, den;
  
  /****************************************END OF SET-UP**************************************************/
  for (int iEvent = 0; iEvent != nEvent; ++iEvent) 
  {
    if (!pythia.next()) continue;
    
    // weight = info.weight();
    // count = -1;
    num = 0; den = 0;
    cout<<endl;
    for(int t=0; t != event.size(); ++t)
    {
      if(event[t].id() == 23)
      {
        daughter1 = event[t].daughter1();
        daughter2 = event[t].daughter2();

        while(daughter1 == daughter2 && daughter1 != 0)
        {
          daughter1 = event[daughter1].daughter1();
          daughter2 = event[daughter2].daughter2();
          cout<<daughter1<<" "<<daughter2<<endl;
        }
        break;
      }
    }
    
    den = event[daughter1].pT()+event[daughter2].pT();

    temp = daughter2;

    cout<<"daughter1:"<<daughter1<<"\n";
    // cout<<daughter1<<" "<<daughter2<<endl;


    daughter2 = event[daughter1].daughter2(); //dont change the
    daughter1 = event[daughter1].daughter1(); //order of these ;)
    
    cout<<daughter1<<" "<<daughter2<<endl;


    while(daughter1 == daughter2 && daughter1 != 0)
    {
      daughter1 = event[daughter1].daughter1();
      daughter2 = event[daughter2].daughter2();
      cout<<daughter1<<" "<<daughter2<<endl;
    }

    if(!event[daughter1].isFinal() && daughter1 != 0)cout<<daughter1<<" was not final.\n";
    while(!event[daughter1].isFinal() && daughter1 != 0) 
    {
      cout<< daughter1<< " ";
      daughter1 = event[daughter1].daughter1();
    }
    cout<<"----"<<endl;

    if(!event[daughter2].isFinal() && daughter2 != 0)cout<<daughter2<<" was not final.\n";
    while(!event[daughter2].isFinal() && daughter2 != 0) 
    {
      cout<< daughter2<<" ";
      daughter2 = event[daughter2].daughter1();
    }
    cout<<"----"<<endl;


    cout<<"daughter2:\n";
    daughter1 = event[temp].daughter1();
    daughter2 = event[temp].daughter2();
    cout<<daughter1<<" "<<daughter2<<endl;
    
    
    while(daughter1 == daughter2 && daughter1 != 0)
    {
      // cout<<daughter1<<" "<<daughter2<<endl;
      daughter1 = event[daughter1].daughter1();
      daughter2 = event[daughter2].daughter2();
      cout<<daughter1<<" "<<daughter2<<endl;
      
    }
    
    if(!event[daughter1].isFinal() && daughter1 != 0)cout<<daughter1<<" was not final.\n";
    while(!event[daughter1].isFinal() && daughter1 != 0) 
    {
      cout<< daughter1<< " ";
      daughter1 = event[daughter1].daughter1();
    }
    cout<<"----"<<endl;

    if(!event[daughter2].isFinal() && daughter2 != 0) cout<<daughter2<<" was not final.\n";
    while(!event[daughter2].isFinal() && daughter2 != 0) 
    {
      cout<< daughter2<<" ";
      daughter2 = event[daughter2].daughter1();
    }
    cout<<"----"<<endl;






      // state = abs(event[t].status());      
      
      // if(event[t].id() < 30 && event[t].id() > 20) 
      // {
      //   ++count;
      //   status[count] = 3;
      // }
      
      // else if(state == 71 || state == 72 || state == 61 || state == 62 || state == 63) 
      // { 
      //   ++count;
      //   status[count] = 2;
      // }

    //   if(event[t].isVisible() && abs(event[t].id()) == 11)
    //   {
    //     ++count;
    //     if(abs(event[event[t].mother1()].id() == 11) || abs(event[event[t].mother2()].id() == 11))
    //     {
    //       num + = 
    //       cout<<"found: "<<t<<endl;

    //     }
    //     // else status[count] = 1;
    //   }

    //   else continue;

    //   id[count] = event[t].id();
    //   pT[count] = event[t].pT();
    //   eta[count] = event[t].eta();
    //   phi[count] = event[t].phi();
    //   m[count] = event[t].m();
    // }

    // n = UShort_t(count+1);
    // tree->Fill();
  
  }

  tree->Print();
  tree->AutoSave("Overwrite"); 
  outFile.Close();
  cout<<"Done.\n";
  cout<<"TTree is saved in "<<outputFilename<<endl;
  return 0;
}
