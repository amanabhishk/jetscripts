//Takes .root file as input and applies the physics definition for flavor tagging

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
#include "fastjet/tools/GridMedianBackgroundEstimator.hh"

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
#include "TLorentzVector.h"


using namespace Pythia8;

int main(int argc, char* argv[]) 
{
  std::clock_t start = std::clock();

  char* options[argc];
  if(argc != 4)
  {
    cout<<"Incorrect number of arguments. Correct input is:\n" ;
    cout<<"$ ./physics_def.exe [inputfile.root] [outputfile.root] [1(dijet) 2(Zjet) 3(gamma)]\n";
    exit(0);
  }
  for(unsigned short int t = 0; t < argc ; ++t) options[t] = argv[t];
  unsigned short int sample = atoi(options[3]);

  TApplication theApp("event_generation", &argc, argv);

  // fastjet setup
  fastjet::JetDefinition jetDef(fastjet::antikt_algorithm, R, fastjet::E_scheme, fastjet::Best);
  std::vector <fastjet::PseudoJet> fjInputs; //particles that will be clustered into jets

  TFile *f = new TFile(options[1]);              //input file with events
  assert(f && !f->IsZombie());
  TTree *Events = (TTree*)f->Get("Events");
  assert(Events && !Events->IsZombie());
  
  TFile outFile(options[2], "RECREATE");    //output file 

  //Extracting information from the TTree
  unsigned int nEvent = (unsigned int)Events->GetEntries(); 
  cout<<nEvent<<" events found in "<<options[1]<<"."<<endl;

  unsigned short int eventParticleCount;
  int id[size];
  float weight, pT[size], eta[size], phi[size], m[size];
  unsigned char status[size];
  
  Events->SetBranchAddress("n", &eventParticleCount);
  Events->SetBranchAddress("weight", &weight);
  Events->SetBranchAddress("id", id);
  Events->SetBranchAddress("status", status);
  Events->SetBranchAddress("pT", pT);
  Events->SetBranchAddress("eta", eta);
  Events->SetBranchAddress("phi", phi);
  Events->SetBranchAddress("m", m);


  jet_data phy_def_jets; //class for storing jet data

  TLorentzVector v;       
  vector<int> leptonList;   //Z+jet
  vector<int> partonList;   
  int gamma;                //gamma+jet

  /**************************************END OF SET-UP**************************************/

    cout<<"Using physics definition.\n";
    cout<<"Progress:\n";
    for (unsigned int iEvent = 0; iEvent < nEvent; ++iEvent) 
    {
      
      if(100*iEvent == 25*nEvent) cout<< "25%\n";
      if(100*iEvent == 50*nEvent) cout<< "50%\n";
      if(100*iEvent == 75*nEvent) cout<< "75%\n";
      if(iEvent == nEvent-1) cout<< "100%.......Done.\n";

      Events->GetEntry(iEvent);

      leptonList.resize(0);
      fjInputs.resize(0);
      partonList.resize(0);

      //select relevant events and make partonList and fjInputs vectors
      for (unsigned int i = 0; i != eventParticleCount; ++i) 
      {
        if ( status[i] == 1 ) 
        {   
          v.SetPtEtaPhiM(pT[i],eta[i],phi[i],m[i]);
          fastjet::PseudoJet particleTemp = v;
          particleTemp.set_user_index(-1);
          fjInputs.push_back( particleTemp );
        }
        else if(status[i] == 2)
        {
          leptonList.push_back(i);
          gamma = i;
        }
        else if(status[i] == 3)
        {
          partonList.push_back(i);
        }
      }//Event selector loop
    
      if (fjInputs.size() == 0) continue;
      
      //clustering using fastjet
      vector <fastjet::PseudoJet> unsortedJets, sortedJets;
      fastjet::ClusterSequence jetCluster(fjInputs, jetDef);

      unsortedJets = jetCluster.inclusive_jets( pTMin );
      sortedJets = sorted_by_pt(unsortedJets);

      if(sortedJets.size()==0)  continue;

      struct clustered_info clusteredData;
      clusteredData.leptonList = leptonList;
      clusteredData.gamma = gamma;
      clusteredData.sortedJets = sortedJets;
      clusteredData.weight = weight;

      if(!is_good_event(pT,eta,phi,m,clusteredData,sample)) continue;

      vector <int> jetFlavor(sortedJets.size(),0);
      
      //flavor tagging begins
      for (unsigned int i = 0; i != sortedJets.size(); ++i) 
      {      
        vector<fastjet::PseudoJet> jetParts = sortedJets[i].constituents();
        if ( jetParts.size() == 1 ) continue;
        if (i > 1) break;
      
        for(unsigned int k = 0; k != partonList.size(); ++k)  
        {

          double dR = deltaR( phi[partonList[k]], sortedJets[i].phi(), eta[partonList[k]],sortedJets[i].eta());
          if ( dR < R ) 
          {
            if(jetFlavor[i]!=0)
            {
              jetFlavor[i] = 0;
              break;
            }    
            else jetFlavor[i] = abs(id[partonList[k]]);
          }
        }//tagging loop
      }//Loop over leading jets
      
      //store jet data
      clusteredData.jetFlavor = jetFlavor;
      phy_def_jets.fill(clusteredData,sample);
    }//Event loop
    
  phy_def_jets.save();
  outFile.Close();
  printTime((std::clock()-start));
  cout<<"Analysis stored in "<<options[2]<<endl;
  
  return 0;
}
