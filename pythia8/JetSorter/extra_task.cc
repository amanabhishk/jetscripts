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
    cout<<"$ ./physicsDef.exe [inputfile.root] [outputfile.root] [1(dijet) 2(Zjet) 3(gamma)]\n";
    exit(0);
  }
  for(unsigned short int t = 0; t < argc ; ++t) options[t] = argv[t];
  unsigned short int sample = atoi(options[3]);

  TApplication theApp("event_generation", &argc, argv);
  unsigned int size = 2000;
  
  int ptBins = 48.;
  const double ptRange[]=
    {18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84,
    97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468,
    507, 548, 592, 638, 686, 737, 790, 846, 905, 967,
    1032, 1101, 1172, 1248, 1327, 1410, 1497, 1588, 1684, 1784, 1890, 2000};

  double R      = 0.5;    // Jet size.
  double pTMin  = 10.0;   // Min jet pT
  double etaMax = 1.3;    // Pseudorapidity range
  int count;              //keeping track of tagged jets;
  bool dijetCriteria;     //selection of good dijet events

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


  float Jw, JpT, JpTD, Jsigma2[]={0,0};
  unsigned int Jmul[]={0,0}, hadronCount;
  unsigned char Jflavor;

  TTree* tree = new TTree("tree","tree");
  tree->Branch("jet_weight", &Jw , "jet_weight/F" );
  tree->Branch("jet_pT", &JpT , "jet_pT/F" );
  tree->Branch("jet_pTD", &JpTD , "jet_pTD/F" );
  tree->Branch("jet_sigma2", Jsigma2 , "jet_sigma2[2]/F" );
  tree->Branch("jet_multiplicity", Jmul , "jet_multiplicity[2]/i" );
  tree->Branch("jet_flavor", &Jflavor , "jet_flavor/b" );
  tree->Branch("hadron_count",&hadronCount,"hadron_count/i");

  TLorentzVector v, v1, v2;       //for converting to cartesian coordinates
  vector<int> leptonList;
  int gamma = -1;

  /**************************************END OF SET-UP**************************************/

    cout<<"Catching hadrons.\n";
    cout<<"Progress:\n";
    for (unsigned int iEvent = 0; iEvent < nEvent; ++iEvent) 
    {
      
      if(100*iEvent == 25*nEvent) cout<< "25%\n";
      if(100*iEvent == 50*nEvent) cout<< "50%\n";
      if(100*iEvent == 75*nEvent) cout<< "75%\n";
      if(iEvent == nEvent-1) cout<< "100%.......Done.\n";

      Events->GetEntry(iEvent);

      //cout<<iEvent<<endl;
      fjInputs.resize(0);
      //select relevant events and make partonList and fjInputs vectors
      for (unsigned int i = 0; i != eventParticleCount; ++i) 
      {
        if ( status[i] == 1 ) 
        {   
          v.SetPtEtaPhiM(pT[i],eta[i],phi[i],m[i]);
          fastjet::PseudoJet particleTemp = v;
          particleTemp.set_user_index(id[i]);
          fjInputs.push_back( particleTemp );
        }
      }//Event selector loop
    
      if (fjInputs.size() == 0) continue;
      //clustering using fastjet
      vector <fastjet::PseudoJet> unsortedJets, sortedJets;
      fastjet::ClusterSequence jetCluster(fjInputs, jetDef);

      unsortedJets = jetCluster.inclusive_jets( pTMin );
      sortedJets = sorted_by_pt(unsortedJets);

      if(sortedJets.size()==0)  continue;

      JpT = sortedJets[0].pt();
      hadronCount = hadron_count(sortedJets);
      tree->Fill();

    }//Event loop
    
  tree->AutoSave("Overwrite");
  outFile.Close();
  printTime((std::clock()-start));
  cout<<"Analysis stored in "<<options[2]<<endl;
  
  return 0;
}
