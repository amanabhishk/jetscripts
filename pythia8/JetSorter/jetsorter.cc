// This class sorts pythia8 jets with the fastjet algorithm. See READMEi_ScriptInfo for further details.

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


using namespace Pythia8;

int main(int argc, char* argv[]) 
{

  bool verbose = false;
  // Create the ROOT application environment.
  TApplication theApp("event_generation", &argc, argv);

  // Settings
  int  nEvent = 100;
  // 0 for gluon jet, 1 for all quarks, 2 for light quarks, 3 for heavy quarks
  if (argc > 1){
    nEvent = atoi(argv[1]);
  }
  int ptBins = 48.;
  const double ptRange[]=
    {18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84,
    97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468,
    507, 548, 592, 638, 686, 737, 790, 846, 905, 967,
    1032, 1101, 1172, 1248, 1327, 1410, 1497, 1588, 1684, 1784, 1890, 2000};

  // int power     = -1;     // -1 = ant-kT; 0 = C/A; 1 = kT
  double R      = 0.5;    // Jet size.
  double pTMin  = 10.0;   // Min jet pT
  double etaMax = 1.3;    // Pseudorapidity range
  // int weightedPt = 1;

  // Create Pythia instance and set it up to generate hard QCD processes
  // above pTHat = 20 GeV for pp collisions at 14 TeV.
  Pythia pythia;
  Event& event = pythia.event;
  Info& info = pythia.info;

  pythia.readString("HardQCD:all = on");
  pythia.readString("PhaseSpace:pTHatMin = 30.");

  // pythia.readString("particleDecays:limitTau0=on");
  // pythia.readString("particleDecays:tauMax=10.");
  pythia.readString("Next:numberShowInfo = 0");
  pythia.readString("Next:numberShowProcess = 0");
  pythia.readString("Next:numberShowEvent = 0");
  pythia.readString("PhaseSpace:bias2Selection = on");
  pythia.readString("PhaseSpace:bias2SelectionPow = 4.5");
  // pythia.readString("PartonLevel:ISR=on");
  //pythia.particleData.listAll();

  // PtHatReweightUserHook ptGenReweight;

  // if (weightedPt){
  //   pythia.setUserHooksPtr( &ptGenReweight );
  // }


  pythia.readString("Beams:eCM = 8000.");
  pythia.init();
  pythia.settings.listChanged();

  // Create file on which histogram(s) can be saved.
  TFile outFile("output.root", "RECREATE");
  TProfile gluonFrac("g","g",ptBins,ptRange);
  // TProfile gluonFrac2("g + U","g + U",ptBins,ptRange);
  TProfile lightquarkFrac("lq","lq",ptBins,ptRange);
  TProfile charmFrac("c","c",ptBins,ptRange);
  TProfile bottomFrac("b","b",ptBins,ptRange);
  TProfile unmatchedFrac("unmatched","unmatched",ptBins,ptRange);



  // fastjet::JetDefinition jetDef(fastjet::genkt_algorithm, R, power);
  fastjet::JetDefinition jetDef(fastjet::antikt_algorithm, R, fastjet::E_scheme, fastjet::Best);

  // Fastjet input
  std::vector <fastjet::PseudoJet> fjInputs;

  // std::vector<std::vector> frac,

  std::clock_t start = std::clock();
  double time_processor = 0; int hours; int minutes; int seconds; 
  
  TH1D* taggedJets =  new TH1D("taggedJets","taggedJets",10, 0.5, 10.5);
  TH1D* distribution =  new TH1D("distribution","distribution",ptBins,ptRange);
  // TH1D *pt_spectrum_0 = new TH1D("pt without weight","pt without weight",ptBins,ptRange);
  // TH1D *pt_spectrum = new TH1D("pt spectrum","pt spectrum",ptBins,ptRange);

  int count; //etaBreak;
  bool dijetCriteria;
  double weight;// weightSum = 0;
  //////////////////////////END OF SET-UP////////////////////////////
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) 
  {
    // cout<<
    
    if(verbose)
    {
      cout<<"-------------------\n";
      cout<<"Event#"<<iEvent+1<<endl;
    }
    if (!pythia.next()) 
    {
      if(verbose)cout<<"pythia.next():skip event.";
      continue;
    }

    if (iEvent!=0&&iEvent%100==0 && verbose)
    {
      time_processor = (std::clock() - start)/(( (double) CLOCKS_PER_SEC ) );
      time_processor = time_processor*( ((double) nEvent)/iEvent-1); 
      minutes =  time_processor/60; hours = minutes/60;
      seconds = time_processor-60*minutes;
      minutes = minutes - hours*60;
      cout << iEvent << " events analyzed. Eta : " << hours << "h" <<
        minutes << "m" << seconds << "s." << endl;
    }  
    
    fjInputs.resize(0);
    vector<int> partonList;
  
    for (int i = 0; i != event.size(); ++i) 
    {
      double status = abs( event[i].status() ); 
      if( status == 23 ) partonList.push_back(i);
      if ( event[i].isFinal() && event[i].isVisible() ) 
      {   
        fastjet::PseudoJet particleTemp = event[i];
        fjInputs.push_back( particleTemp );
        distribution->Fill(event[i].pT());
      }
    }//Event selector loop
  
    if (fjInputs.size() == 0) 
    {
      if(verbose) cout << "Error: event with no final state particles" << endl;
      continue;
    }
    
    vector <fastjet::PseudoJet> unsortedJets, sortedJets;
    fastjet::ClusterSequence jetCluster(fjInputs, jetDef);

    unsortedJets = jetCluster.inclusive_jets( pTMin );
    sortedJets = sorted_by_pt(unsortedJets);

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
    if(verbose)
    {  
      cout<<"parton1: "<<event[partonList[0]].eta()<<" "<<event[partonList[0]].phi()<</*" "<<event[partonList[0]].eT()<<*/endl;
      cout<<"parton2: "<<event[partonList[1]].eta()<<" "<<event[partonList[1]].phi()<</*" "<<event[partonList[1]].eT()<<*/endl;
    }

    count = 0;
    for (unsigned int i = 0; i != sortedJets.size(); ++i) 
    {      
      
      

      weight = info.weight();
      // cout<< sortedJets[i].pt()<<"  "<<weight<<endl;
      // pt_spectrum->Fill(abs(sortedJets[i].pt()));//,weight);

      if(abs(sortedJets[i].eta())>etaMax)
      {
        if(verbose) cout<< "etaMax condition violated.\n";
        continue;
      }

      vector<fastjet::PseudoJet> jetParts = sortedJets[i].constituents();
      if ( jetParts.size() == 1 ) 
      {
        if(verbose) cout<<"Trivial jet found.\n";
        continue;
      }
 
      if(verbose) cout<<"jet#"<<i+1<<":  "<<sortedJets[i].eta()<<" "<<sortedJets[i].phi()<</*" "<<sortedJets[i].Et()<<*/endl;
      
      
      for(unsigned int k = 0; k != partonList.size(); ++k)  
      {

        double dR = deltaR( event[partonList[k]].phi(), sortedJets[i].phi(),event[partonList[k]].eta(),sortedJets[i].eta());
        // cout<<dR<<" ";
        // count += 1;
        if ( dR < R ) 
        {
          count += 1;
          assert(jetFlavor[i]==0);     
          // cout<<k<<":"<<event[partonList[k]].id()<<endl;
          jetFlavor[i] = abs(event[partonList[k]].id());
        }
      }//tag tagging loop
      if(verbose) cout<<endl;
      
    }//Loop over each jet
    
    taggedJets->Fill(count);

    if(verbose)
    {
      cout<<"jetFlavor: ";for(int i=0; i!= jetFlavor.size(); ++i) cout<<jetFlavor[i]<<" ";
      cout<<endl;
    }

    //fill histograms
    
    // weightSum += weight;
    for(int k = 0; k != 2; ++k)
    {
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

  taggedJets->Write();
  // distribution->Write();
  // pt_spectrum->Write();

  return 0;
}//Done
