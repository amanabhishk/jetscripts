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
  //QCD aware clustering of status 70 particles
  AntiKtMeasure *antikt_measure = new AntiKtMeasure(0.4);
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


  jet_data QCDaware_def_jets;

  //TH1D* efficiency = new TH1D("tagged jets", "tagged jets", 100,0,10);
  //TH1D* jet_count = new TH1D("jet_count","jet_count",100,0,50);

  TLorentzVector v, v1, v2;       //for converting to cartesian coordinates
  vector<int> leptonList;
  int gamma = -1;
  fastjet::PseudoJet particleTemp;

  /**************************************END OF SET-UP**************************************/

  cout<<"Using QCD-aware algorithm.\n";
  cout<<"Progress:\n";
  for (unsigned int iEvent = 0; iEvent < nEvent; ++iEvent) 
  {
    
    if(100*iEvent == 25*nEvent) cout<< "25%\n";
    if(100*iEvent == 50*nEvent) cout<< "50%\n";
    if(100*iEvent == 75*nEvent) cout<< "75%\n";
    if(iEvent == nEvent-1) cout<< "100%.......Done.\n";

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

      if (status[i] == 1) 
      {   
        v.SetPtEtaPhiM(pT[i],eta[i],phi[i],m[i]);
        particleTemp = v;
        particleTemp.set_user_index(0);
        clustering_input.push_back( particleTemp );
      }

      if(status[i] == 2)
      {
        leptonList.push_back(i);
        gamma = i;
      }
    }//Event selector loop
  
    if (qcd_aware_clustering_input.size() == 0) continue;
    if (clustering_input.size() == 0) continue;
    
    //QCD aware clustering of status 70 particles using fastjet
    vector <fastjet::PseudoJet> jets_qcd, unsortedJets, sortedJets;
    fastjet::ClusterSequence qcd_aware_clustering(qcd_aware_clustering_input, qcd_aware_clustering_definition);
    jets_qcd = qcd_aware_clustering.inclusive_jets( pTMin );

    //Adding QCD clustered jets as ghost particles to normal clustering which will follow.
    //jet_count->Fill(jets_qcd.size());
    if(jets_qcd.size()==0)
    {
      continue;
    }

    // for(unsigned int k = 0; k != jets_qcd.size(); ++k)
    // {
    //   v.SetPtEtaPhiM(jets_qcd[k].pt(),jets_qcd[k].eta(),jets_qcd[k].phi(),jets_qcd[k].m());
    //   v *= pow(10,-18);
    //   particleTemp.set_user_index( jets_qcd[k].user_index() );
    //   clustering_input.push_back( particleTemp );
    // }

    //clustering final state particles with ghosts obtained from QCD aware
    fastjet::ClusterSequence jetCluster(clustering_input, normal_jet_clustering_definition);
    unsortedJets = jetCluster.inclusive_jets( pTMin );
    sortedJets = sorted_by_pt(unsortedJets);

    if(sortedJets.size()==0)  continue;

    struct clustered_info clusteredData;
    clusteredData.leptonList = leptonList;
    clusteredData.gamma = gamma;
    clusteredData.sortedJets = sortedJets;

    if(!is_good_event(pT,eta,phi,m,clusteredData,sample));
    
    //flavor tagging begins
    vector <int> jetFlavor(sortedJets.size(),0);

    for (unsigned int i = 0; i != sortedJets.size(); ++i) 
    {      
      //if(pTD(sortedJets) != sortedJets.size()) cout<< "NO.\n";
      vector<fastjet::PseudoJet> jetParts = sortedJets[i].constituents();
      if ( jetParts.size() == 1 ) continue;
      
      //match with status 23 particles and assign flavor to (2) leading jets
      for(unsigned int k = 0; k != jets_qcd.size(); ++k)  
      {

        double dR = deltaR( jets_qcd[k].phi(), sortedJets[i].phi(), jets_qcd[k].eta(),sortedJets[i].eta());
        if ( dR < R/2 ) 
        {
          // taggedJets-break1);
          if(jetFlavor[i]!=0) 
          {
            jetFlavor[i] = 0;
            break;    
          }
          jetFlavor[i] = abs(jets_qcd[k].user_index());
        }
      }//tag tagging loop
    }//Loop over leading jets


    //store jet data
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
