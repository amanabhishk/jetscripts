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
  if(argc != 3)
  {
    cout<<"Incorrect number of arguments. Correct input is:\n $ ./physicsDef.exe [inputfile.root] [outputfile.root]\n";
    exit(0);
  }
  for(unsigned short int t = 0; t < argc ; ++t) options[t] = argv[t];

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
  bool ZjetCriteria;     //selection of Z+jet events

  // fastjet setup
  fastjet::JetDefinition jetDef(fastjet::antikt_algorithm, R, fastjet::E_scheme, fastjet::Best);
  std::vector <fastjet::PseudoJet> fjInputs; //particles that will be clustered into jets

  TFile *f = new TFile(options[1]);              //input file with events
  TTree *Events = (TTree*)f->Get("Events");
  
  TFile outFile(options[2], "RECREATE");    //output file 

  //Extracting information from the TTree
  unsigned int nEvent = (unsigned int)Events->GetEntries(); 
  cout<<nEvent<<" events found in the ROOT file."<<endl;

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
  unsigned int Jmul[]={0,0};
  unsigned char Jflavor;

  TTree* tree = new TTree("tree","tree");
  tree->Branch("jet_weight", &Jw , "jet_weight/F" );
  tree->Branch("jet_pT", &JpT , "jet_pT/F" );
  tree->Branch("jet_pTD", &JpTD , "jet_pTD/F" );
  tree->Branch("jet_sigma2", Jsigma2 , "jet_sigma2[2]/F" );
  tree->Branch("jet_multiplicity", Jmul , "jet_multiplicity[2]/i" );
  tree->Branch("jet_flavor", &Jflavor , "jet_flavor/b" );

  TLorentzVector v;//,v1,v2;       //for converting to cartesian coordinates
  int parton = -1, gamma = -1;
  bool check1, check2;
  /**************************************END OF SET-UP**************************************/
  cout<<"Progress:\n";
  for (unsigned int iEvent = 0; iEvent < nEvent; ++iEvent) 
  {
    
    if(100*iEvent == 25*nEvent) cout<< "25%\n";
    if(100*iEvent == 50*nEvent) cout<< "50%\n";
    if(100*iEvent == 75*nEvent) cout<< "75%\n";
    if(iEvent == nEvent-1) cout<< "100%.......Done.\n";

    cout << std::setprecision(10);
    Events->GetEntry(iEvent);

    fjInputs.resize(0);
    //vector<int> leptonList;
    check1 = true, check2 = true;

    //select relevant events and make partonList and fjInputs vectors
    for (unsigned int i = 0; i != eventParticleCount; ++i) 
    {
      if( status[i] == 3 ) 
      {
        parton = i;
        assert(check1);
        check1 = false;
      }
      else if ( status[i] == 1 ) 
      {   
        v.SetPtEtaPhiM(pT[i],eta[i],phi[i],m[i]);
        fastjet::PseudoJet particleTemp = v;
        fjInputs.push_back( particleTemp );
      }
      else if(status[i] == 2) 
      {
        gamma = i;
        assert(check2);
        check2 = false;
      }
    }//Event selector loop
  
    if (fjInputs.size() == 0) assert(false);
    
    vector <fastjet::PseudoJet> unsortedJets, sortedJets;
    fastjet::ClusterSequence jetCluster(fjInputs, jetDef);

    unsortedJets = jetCluster.inclusive_jets( pTMin );
    sortedJets = sorted_by_pt(unsortedJets);

    if(sortedJets.size()==0)  continue;
    if(sortedJets[1].pt()>0.3*pT[gamma]) continue;
    if(deltaR(phi[gamma],sortedJets[0].phi(),eta[gamma],sortedJets[0].eta()) < R) continue;

    unsigned short int jetFlavor = 0;

    if ( sortedJets[0].constituents().size() == 1 ) continue;
    //if(parton == -1) continue;
    assert(parton != -1);
    //match with status 23 particles and assign flavor to leading jet only
    double dR = deltaR( phi[parton], sortedJets[0].phi(), eta[parton],sortedJets[0].eta());

    if ( dR < R ) 
    {  
      jetFlavor = abs(id[parton]);
    }//tag tagging loop
    
    if(fabs(sortedJets[0].eta()) > etaMax) continue;

    Jw = weight;
    JpT = sortedJets[0].pt();
    Jmul[0] = multiplicity(sortedJets[0],0);
    Jmul[1] = multiplicity(sortedJets[0],2);
    Jflavor = jetFlavor;
    JpTD = pTD(sortedJets[0]);
    sigma2(sortedJets[0],Jsigma2);
    tree->Fill();
  }//Event loop
  
  tree->AutoSave("Overwrite");
  outFile.Close();
  printTime((std::clock()-start));
  cout<<"Analysis stored in "<<options[2]<<endl;
  
  return 0;
}
