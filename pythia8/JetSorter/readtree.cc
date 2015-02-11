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

#include "TTree.h"
#include "TLorentzVector.h"
// #include "TSystem.h"

int main(int argc, char* argv[])
{
	// gSystem->Load("libTree");
	TApplication theApp("event_generation", &argc, argv);
	TFile *f = new TFile("output.root");
	TTree *Events = (TTree*)f->Get("Events");
	TFile outFile("canvas.root", "RECREATE");

	Int_t event_count = (Int_t)Events->GetEntries(); 
	cout<<event_count<<endl;
	Float_t phi[5000]; //with nmax greater or equal to the max value for nGenjet5.
	
	Events->SetBranchAddress("phi",phi);
	UShort_t size;
	Events->SetBranchAddress("n",&size);
	TH1D *phi_hist = new TH1D("phi","phi",1000,0,5000);
	for (Int_t i=0; i != event_count; ++i) 
	{
	// 	UShort_t size;
 //  		
  		Events->GetEntry(i);
  		// cout<<size<<endl;
 //  		// UShort_t id[size];
		Float_t phi[size], phi[size], phi[size], m[size];
  		// cout<< size<< " ";
  		// phi_hist->Fill(size);
  		for(UShort_t k = 0; k != size; ++k)
  		{
			phi_hist->Fill(phi[k]);
	// // // 		// cout<<k+1<<" "<<id[k]<<" "<<phi[k]<<" "<<phi[k]<<" "<<m[k]<<endl;


  		}
  		
  		// cout<<"*************************************\n";

	}
	phi_hist->Draw();
	phi_hist->Write();

	return 0;
	
}

