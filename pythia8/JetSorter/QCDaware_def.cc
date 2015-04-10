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
#include "/home/aabhishe/fastjet-3.1.0/plugins/QCDaware/QCDAware.hh"

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
using namespace fastjet::contrib;

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

  // fastjet setup
  //QCD aware clustering of status 70 particles
  AntiKtMeasure *antikt_measure = new AntiKtMeasure(R);
  QCDAware *qcd_aware_clustering_definition = new QCDAware(antikt_measure);
  std::vector <fastjet::PseudoJet> qcd_aware_clustering_input; //particles that will be clustered into jets

  //Normal clustering of final state particles
  fastjet::JetDefinition normal_jet_clustering_definition(fastjet::antikt_algorithm, R, fastjet::E_scheme, fastjet::Best);
  std::vector <fastjet::PseudoJet> clustering_input;

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


  jet_data QCDaware_def_jets; // class for storing jet data

  TH1D* efficiency = new TH1D("tagged jets", "tagged jets", 100,0,10);
  //TH1D* jet_count = new TH1D("jet_count","jet_count",100,0,50);

  TLorentzVector v;       //for converting to cartesian coordinates
  vector<int> leptonList;
  int gamma;
  fastjet::PseudoJet particleTemp;

  int index;
  bool gluonPresence;
  double hardest_pT = 0;
  /**************************************END OF SET-UP**************************************/

  cout<<"Using QCD-aware algorithm.\n";
  cout<<"Progress:\n";
  for (unsigned int iEvent = 0; iEvent < nEvent; ++iEvent) 
  {

    Events->GetEntry(iEvent);

    cout<<iEvent<<endl;
    leptonList.resize(0);
    qcd_aware_clustering_input.resize(0);
    clustering_input.resize(0);

    //select relevant events and make jets_qcd and qcd_aware_clustering_input vectors
    for (unsigned int i = 0; i != eventParticleCount; ++i) 
    {
      //if(isCharm(id[i]) && abs(id[i]) != 4)cout<<" "<<id[i]<<endl;
      if (status[i] == 70) 
      {   
        v.SetPtEtaPhiM(pT[i],eta[i],phi[i],m[i]);
        particleTemp = v;
        particleTemp.set_user_index(id[i]);
        qcd_aware_clustering_input.push_back( particleTemp );
      }

      else if (status[i] == 1) 
      {   
        v.SetPtEtaPhiM(pT[i],eta[i],phi[i],m[i]);
        particleTemp = v;
        particleTemp.set_user_index(0);
        clustering_input.push_back( particleTemp );
      }

      else if(status[i] == 2)
      {
        leptonList.push_back(i);
        gamma = i;
      }
    }//Event selector loop
  
    if (qcd_aware_clustering_input.size() == 0) continue;
    //if (clustering_input.size() == 0) continue; //obsolete since jets_qcd.size()==0 is checked
    
    //QCD aware clustering of status 70 particles using fastjet
    vector <fastjet::PseudoJet> jets_qcd, unsortedJets, sortedJets;
    fastjet::ClusterSequence qcd_aware_clustering(qcd_aware_clustering_input, qcd_aware_clustering_definition);
    jets_qcd = qcd_aware_clustering.inclusive_jets( pTMin );

    //Adding QCD clustered jets as ghost particles to normal clustering which will follow.
    if(jets_qcd.size()==0) continue;

    for(unsigned int k = 0; k != jets_qcd.size(); ++k)
    {
      v.SetPtEtaPhiM(jets_qcd[k].pt(),jets_qcd[k].eta(),jets_qcd[k].phi(),jets_qcd[k].m());
      particleTemp = v*pow(10,-18);
      particleTemp.set_user_index( jets_qcd[k].user_index() );
      clustering_input.push_back( particleTemp );
    }

    //clustering final state particles with ghosts obtained from QCD aware
    fastjet::ClusterSequence jetCluster(clustering_input, normal_jet_clustering_definition);
    unsortedJets = jetCluster.inclusive_jets( pTMin );
    sortedJets = sorted_by_pt(unsortedJets);

    if(sortedJets.size()==0)  continue;

    struct clustered_info clusteredData;
    clusteredData.leptonList = leptonList;
    clusteredData.gamma = gamma;
    clusteredData.sortedJets = sortedJets;

    if(!is_good_event(pT,eta,phi,m,clusteredData,sample)) continue;
    
    //flavor tagging begins
    vector <int> jetFlavor(sortedJets.size(),0);

    for (unsigned int i = 0; i != sortedJets.size(); ++i) 
    {      
      vector<fastjet::PseudoJet> jetParts = sortedJets[i].constituents();
      if ( jetParts.size() == 1 ) continue;
      if( i > 1 ) break;
      gluonPresence = false;
      hardest_pT = 0;

      for(unsigned int k = 0; k != jetParts.size(); ++k)  
      {
        
        index = abs(jetParts[k].user_index());
        if(index != 0 && hardest_pT < pow(10,18)*jetParts[k].pt())
        {
          jetFlavor[i] = index;
          hardest_pT = pow(10,18)*jetParts[k].pt();
        }
      }
    }//tagging loop


    //store jet data
    clusteredData.weight = weight;
    clusteredData.jetFlavor = jetFlavor;
    QCDaware_def_jets.fill(clusteredData,sample);
  }//Event loop
    
  //efficiency->Write();
  //jet_count->Write();
  //cout<<jetcount<<endl;
  //tree->AutoSave("Overwrite");
  QCDaware_def_jets.save();
  outFile.Close();
  printTime((std::clock()-start));
  cout<<"Analysis stored in "<<options[2]<<endl;
  
  return 0;
}
