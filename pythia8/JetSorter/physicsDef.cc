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

  // Book histograms
  TProfile gluonFrac("g","g",ptBins,ptRange);
  TProfile lightquarkFrac("lq","lq",ptBins,ptRange);
  TProfile charmFrac("c","c",ptBins,ptRange);
  TProfile bottomFrac("b","b",ptBins,ptRange);
  TProfile unmatchedFrac("unmatched","unmatched",ptBins,ptRange);
  TH1D* taggedJets =  new TH1D("taggedJets","taggedJets",10, 0.5, 10.5);
  TH1D* jetMultiplicity =  new TH1D("jetMultiplicity","jetMultiplicity",100, 0, 20);
  TH1D* visSize = new TH1D("visSize","visSize",500,0,1000);

  TFile *f = new TFile("50000events_gammaJet.root");              //input file with events
  TTree *Events = (TTree*)f->Get("Events");
  
  TFile outFile("jets_output.root", "RECREATE");    //output file 

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

  // unsigned int limit = 20;

  float Jphi[20], Jeta[20], JpT[20], Jm[20], Je[20];
  float Gphi, Geta, GpT, Gm, Ge;
  float Zphi, Zeta, ZpT, Zm, Ze;
  unsigned short int jet_multiplicity;
  
  TTree* tree = new TTree("tree","tree");
  // TTree* leptons = new TTree("Leptons","Leptons");
  // TTree* Z = new TTree("Z","Z");

  tree->Branch("jet_multiplicity", &jet_multiplicity, "jet_multiplicity/s");
  tree->Branch("jet_pT", JpT, "jet_pT[jet_multiplicity]/F");
  tree->Branch("jet_eta", Jeta, "jet_eta[jet_multiplicity]/F");
  tree->Branch("jet_phi", Jphi, "jet_phi[jet_multiplicity]/F");
  tree->Branch("jet_m", Jm, "jet_m[jet_multiplicity]/F");
  tree->Branch("jet_e", Je, "jet_e[jet_multiplicity]/F");
  
  tree->Branch("weight", &weight, "weight/F");

  tree->Branch("gamma_pT", &GpT, "gamma_pT/F");
  tree->Branch("gamma_eta", &Geta, "gamma_eta/F");
  tree->Branch("gamma_phi", &Gphi, "gamma_phi/F");
  //tree->Branch("gamma_m", &Gm, "gamma_m/F");
  tree->Branch("gamma_e", &Ge, "gamma_e/F");
  
  //tree->Branch("Z_pT", &ZpT, "Z_pT/F");
  //tree->Branch("Z_eta", &Zeta, "Z_eta/F");
  //tree->Branch("Z_phi", &Zphi, "Z_phi/F");
  //tree->Branch("Z_m", &Zm, "Z_m/F");
  //tree->Branch("Z_e", &Ze, "Z_e/F");
  
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
      else assert(false);

    }//Event selector loop
  
    assert(gamma != -1);

    if (fjInputs.size() == 0) assert(false);
    
    vector <fastjet::PseudoJet> unsortedJets, sortedJets;
    fastjet::ClusterSequence jetCluster(fjInputs, jetDef);

    unsortedJets = jetCluster.inclusive_jets( pTMin );
    sortedJets = sorted_by_pt(unsortedJets);

    jetMultiplicity->Fill(sortedJets.size());
    jet_multiplicity = sortedJets.size();
    if(sortedJets.size()==0)  continue;

    

    // cout<<"Checkpoint1\n";

    for(unsigned int t = 0; t != sortedJets.size(); ++t)
    {
      Jphi[t] = sortedJets[t].phi();
      Jeta[t] = sortedJets[t].eta();
      JpT[t] = sortedJets[t].pt();
      Jm[t] = sortedJets[t].m();
      Je[t] = sortedJets[t].e();
      // jets->Fill();
    }

    
    Gphi = phi[gamma];
    Geta = eta[gamma];
    GpT = pT[gamma];
    Gm = m[gamma];
      
    v.SetPtEtaPhiM(pT[gamma],eta[gamma],phi[gamma],m[gamma]);
    Ge = v.E();
      // leptons->Fill();
    

    // for(unsigned int t = 0; t != le)
    //v1.SetPtEtaPhiM(pT[leptonList[0]],eta[leptonList[0]],phi[leptonList[0]],m[leptonList[0]]);
    //v2.SetPtEtaPhiM(pT[leptonList[1]],eta[leptonList[1]],phi[leptonList[1]],m[leptonList[1]]);
    //v = v1+v2;
    //Zphi = v.Phi();
    //Zeta = v.Eta();
    //ZpT = v.Pt();
    //Zm = v.M();
    //Ze = v.E();
    // Z->Fill();
    tree->Fill();

    // cout<<"Checkpoint2\n";
    if(deltaR(phi[gamma],sortedJets[0].phi(),eta[gamma],sortedJets[0].eta()) < R) continue;
    //if(deltaR(phi[leptonList[1]],sortedJets[0].phi(),eta[leptonList[1]],sortedJets[0].eta()) < R) continue;

    unsigned short int jetFlavor = 0;

    //vector<fastjet::PseudoJet> jetParts = sortedJets[0].constituents();
    if ( sortedJets[0].constituents().size() == 1 ) continue;
    //if(parton == -1) continue;
    assert(parton != -1);
    //match with status 23 particles and assign flavor to leading jet only
    double dR = deltaR( phi[parton], sortedJets[0].phi(), eta[parton],sortedJets[0].eta());
    // cout<<eta[parton]<<" "<<phi[parton]<<endl;
    // cout<<sortedJets[0].eta()<<" "<<sortedJets[0].phi() <<endl;
    // cout<<dR<<endl<<endl;

    // cout<<"Checkpoint3\n";
    if ( dR < R ) 
    {  
      jetFlavor = abs(id[parton]);
    }//tag tagging loop
    
    if(sortedJets[0].eta() > etaMax) continue;
    


    gluonFrac.Fill(sortedJets[0].pt(), (jetFlavor == 21)? 1:0, weight);
    lightquarkFrac.Fill(sortedJets[0].pt(), (jetFlavor == 1 || jetFlavor == 2 || jetFlavor == 3)? 1:0, weight);
    charmFrac.Fill(sortedJets[0].pt(), (jetFlavor == 4)? 1:0, weight);
    bottomFrac.Fill(sortedJets[0].pt(), (jetFlavor == 5)? 1:0, weight);
    unmatchedFrac.Fill(sortedJets[0].pt(), (jetFlavor == 0)? 1:0, weight);
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

  jetMultiplicity->Write();
  taggedJets->Write();
  visSize->Write();
  
  tree->AutoSave("Overwrite");
  // leptons->AutoSave("Overwrite");
  // Z->AutoSave("Overwrite");
  outFile.Close();

  cout<<"Done in "<<(std::clock()-start)/CLOCKS_PER_SEC<<" seconds. "<<endl;
  
  return 0;
}
