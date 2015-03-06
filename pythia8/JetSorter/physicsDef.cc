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
  bool dijetCriteria;     //selection of good dijet events

  // fastjet setup
  fastjet::JetDefinition jetDef(fastjet::antikt_algorithm, R, fastjet::E_scheme, fastjet::Best);
  std::vector <fastjet::PseudoJet> fjInputs; //particles that will be clustered into jets

  // Book histograms
  TProfile gluonFrac("g","g",ptBins,ptRange);
  TProfile lightquarkFrac("lq","lq",ptBins,ptRange);
  TProfile charmFrac("c","c",ptBins,ptRange);
  TProfile bottomFrac("b","b",ptBins,ptRange);
  TProfile unmatchedFrac("unmatched","unmatched",ptBins,ptRange);
  TH1D* taggedJets =  new TH1D("taggedJets","taggedJets",10, 0.5, 10.5);


  //std::stringstream inputFilename("");
  //inputFilename << argv[1];
  //cout<<options[1]<<endl;
  TFile *f = new TFile(options[1]);              //input file with events
  TTree *Events = (TTree*)f->Get("Events");
  
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


  double num, den;
  TLorentzVector v;       //for converting to cartesian coordinates

  /**************************************END OF SET-UP**************************************/
  cout<<"Progress:\n";
  for (unsigned int iEvent = 0; iEvent < nEvent; ++iEvent) 
  {
    
    if(100*iEvent == 25*nEvent) cout<< "25%\n";
    if(100*iEvent == 50*nEvent) cout<< "50%\n";
    if(100*iEvent == 75*nEvent) cout<< "75%\n";
    if(iEvent == nEvent-1) cout<< "100%.......Done.\n";

    Events->GetEntry(iEvent);

    fjInputs.resize(0);
    vector<int> partonList; //pick out status 23 particles

    //select relevant events and make partonList and fjInputs vectors
    for (unsigned int i = 0; i != eventParticleCount; ++i) 
    {
      if( status[i] == 3 ) partonList.push_back(i);
      if ( status[i] == 1 ) 
      {   
        v.SetPtEtaPhiM(pT[i],eta[i],phi[i],m[i]);
        fastjet::PseudoJet particleTemp = v;
        fjInputs.push_back( particleTemp );
      }
    }//Event selector loop
  
    if (fjInputs.size() == 0) continue;
    assert(partonList.size() == 2);
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

    count = 0;
    
    //flavor tagging begins
    for (unsigned int i = 0; i != sortedJets.size(); ++i) 
    {      
      //if(pTD(sortedJets) != sortedJets.size()) cout<< "NO.\n";
      vector<fastjet::PseudoJet> jetParts = sortedJets[i].constituents();
      if ( jetParts.size() == 1 ) continue;
      
      //match with status 23 particles and assign flavor to (2) leading jets
      for(unsigned int k = 0; k != partonList.size(); ++k)  
      {

        double dR = deltaR( phi[partonList[k]], sortedJets[i].phi(), eta[partonList[k]],sortedJets[i].eta());
        if ( dR < R ) 
        {
          count += 1;
          // taggedJets-break1);
          assert(jetFlavor[i]==0);     
          jetFlavor[i] = abs(id[partonList[k]]);
        }
      }//tag tagging loop
    }//Loop over leading jets
    
    taggedJets->Fill(count);

    //fill histograms
    for(int k = 0; k != jetFlavor.size(); ++k)
    {
      if(k == partonList.size()) break;
      if(sortedJets[k].eta() > etaMax) continue;
      
      Jw = weight;
      JpT = sortedJets[0].pt();
      Jmul[0] = multiplicity(sortedJets[0],0);
      Jmul[1] = multiplicity(sortedJets[0],2);
      Jflavor = jetFlavor[0];
      JpTD = pTD(sortedJets[0]);
      sigma2(sortedJets[0],Jsigma2);
      tree->Fill();

      gluonFrac.Fill(sortedJets[k].pt(), (jetFlavor[k] == 21)? 1:0, weight);
      lightquarkFrac.Fill(sortedJets[k].pt(), (jetFlavor[k] == 1 || jetFlavor[k] == 2 || jetFlavor[k] == 3)? 1:0, weight);
      charmFrac.Fill(sortedJets[k].pt(), (jetFlavor[k] == 4)? 1:0, weight);
      bottomFrac.Fill(sortedJets[k].pt(), (jetFlavor[k] == 5)? 1:0, weight);
      unmatchedFrac.Fill(sortedJets[k].pt(), (jetFlavor[k] == 0)? 1:0, weight);
    }
  }//Event loop
  
  TH1D *lqF = lightquarkFrac.ProjectionX("light quarks","");
  lqF->Write();

  TH1D *gF = gluonFrac.ProjectionX("gluons","");
  gF->Write();

  TH1D *cF = charmFrac.ProjectionX("charm","");
  cF->Write();

  TH1D *bF = bottomFrac.ProjectionX("bottom","");
  bF->Write();

  TH1D *uF = unmatchedFrac.ProjectionX("unmatched","");
  uF->Write();

  tree->AutoSave("Overwrite");
  taggedJets->Write();
  outFile.Close();
  cout<<"Done in "<<(std::clock()-start)/CLOCKS_PER_SEC<<" seconds. "<<endl;
  cout<<"Analysis stored in "<<options[2]<<endl;
  
  return 0;
}
