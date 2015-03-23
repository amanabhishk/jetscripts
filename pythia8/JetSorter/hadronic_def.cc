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

// tdrStyle
#include "tdrstyle_mod1.C"
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
    cout<<"Incorrect number of arguments. Correct input is:\n" ;
    cout<<"$ ./physicsDef.exe [inputfile.root] [outputfile.root]\n";
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
  unsigned int Jmul[]={0,0};
  unsigned char Jflavor;

  TTree* tree = new TTree("tree","tree");
  tree->Branch("jet_weight", &Jw , "jet_weight/F" );
  tree->Branch("jet_pT", &JpT , "jet_pT/F" );
  tree->Branch("jet_pTD", &JpTD , "jet_pTD/F" );
  tree->Branch("jet_sigma2", Jsigma2 , "jet_sigma2[2]/F" );
  tree->Branch("jet_multiplicity", Jmul , "jet_multiplicity[2]/i" );
  tree->Branch("jet_flavor", &Jflavor , "jet_flavor/b" );

  TLorentzVector v;       //for converting to cartesian coordinates

  /**************************************END OF SET-UP**************************************/

    cout<<"Using hadron-based definition.\n";
    cout<<"Progress:\n";
    for (unsigned int iEvent = 0; iEvent < nEvent; ++iEvent) 
    {
      
      if(100*iEvent == 25*nEvent) cout<< "25%\n";
      if(100*iEvent == 50*nEvent) cout<< "50%\n";
      if(100*iEvent == 75*nEvent) cout<< "75%\n";
      if(iEvent == nEvent-1) cout<< "100%.......Done.\n";

      Events->GetEntry(iEvent);

      fjInputs.resize(0);

      //select relevant events and make partonList and fjInputs vectors
      for (unsigned int i = 0; i != eventParticleCount; ++i) 
      {
        //if(isCharm(id[i]) && abs(id[i]) != 4)cout<<" "<<id[i]<<endl;
        if (status[i] == 70) 
        {   
          //if(abs(id[i])>99)cout<<abs(id[i])<<" ";
          v.SetPtEtaPhiM(pT[i],eta[i],phi[i],m[i]);
          v *= pow(10,-18);
          fastjet::PseudoJet particleTemp = v;
          particleTemp.set_user_index(id[i]*100);
          fjInputs.push_back( particleTemp );
        }
        else if (status[i] == 4 || status[i] == 5 || status[i] == 6) 
        {   
          //if(abs(id[i])>99)cout<<abs(id[i])<<" ";
          v.SetPtEtaPhiM(pT[i],eta[i],phi[i],m[i]);
          v *= pow(10,-18);
          fastjet::PseudoJet particleTemp = v;
          particleTemp.set_user_index( (status[i]==6)? 3:status[i] );
          fjInputs.push_back( particleTemp );
        }
        else if ( status[i] == 1 ) 
        {   
          v.SetPtEtaPhiM(pT[i],eta[i],phi[i],m[i]);
          fastjet::PseudoJet particleTemp = v;
          particleTemp.set_user_index(-1);
          fjInputs.push_back( particleTemp );
        }
      }//Event selector loop
    
      if (fjInputs.size() == 0) continue;
      //clustering using fastjet
      vector <fastjet::PseudoJet> unsortedJets, sortedJets;
      fastjet::ClusterSequence jetCluster(fjInputs, jetDef);

      unsortedJets = jetCluster.inclusive_jets( pTMin );
      sortedJets = sorted_by_pt(unsortedJets);

      //selecting good dijet events
      if(sortedJets.size()<2) dijetCriteria = false;
      else if(sortedJets.size()>2)
      {
       dijetCriteria = deltaPhi(sortedJets[0].phi(),sortedJets[1].phi())>2.8 && 0.1*abs(sortedJets[0].pt()+sortedJets[1].pt())>sortedJets[2].pt();
      }
      else
      {
       dijetCriteria = deltaPhi(sortedJets[0].phi(),sortedJets[1].phi())>2.8;
      }
      
      if(!dijetCriteria) continue;

      vector <int> jetFlavor(sortedJets.size(),0);

      cout << std::setprecision(10);
      
      //flavor tagging begins
      int flavor_from_hadron, flavor_from_parton, index;
      for (unsigned int i = 0; i != sortedJets.size(); ++i) 
      {      
        vector<fastjet::PseudoJet> jetParts = sortedJets[i].constituents();
        jetParts = sorted_by_pt(jetParts);

        if ( jetParts.size() == 1 ) continue;
        if(i>1) break;
        
        flavor_from_parton = 0;
        flavor_from_hadron = 0;
        
        for(int p = 0; p != jetParts.size(); ++p)
        {
          index = jetParts[p].user_index();
          if(index == 4 || index == 5)
          {
            if(flavor_from_hadron < index) flavor_from_hadron = index;  
          }
          
          else if(index == 400 || index == 500)
          {
            if(flavor_from_parton < index) flavor_from_parton = index/100;  
          }
        }

        if(flavor_from_parton == 0)
        {
          for(int p = 0; p != jetParts.size(); ++p)
          {
            index = jetParts[p].user_index();
            if(index > 99) 
            {
              flavor_from_parton = index/100;
              break;
            }
          }
        }

        jetFlavor[i]=(flavor_from_hadron == 0)? flavor_from_parton:flavor_from_hadron;
      }//Loop over leading jets
      
      //store jet data
      for(int k = 0; k != jetFlavor.size(); ++k)
      {
        if(k > 1) break;
        if(sortedJets[k].eta() > etaMax) continue;
        
        Jw = weight;
        JpT = sortedJets[0].pt();
        Jmul[0] = multiplicity(sortedJets[0],0);
        Jmul[1] = multiplicity(sortedJets[0],2);
        Jflavor = jetFlavor[0];
        JpTD = pTD(sortedJets[0]);
        sigma2(sortedJets[0],Jsigma2);
        tree->Fill();
      }
    }//Event loop
    
  tree->AutoSave("Overwrite");
  outFile.Close();
  printTime((std::clock()-start));
  cout<<"Analysis stored in "<<options[2]<<endl;
  
  return 0;
}