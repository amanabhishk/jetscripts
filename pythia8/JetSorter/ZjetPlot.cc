#include "tdrstyle_mod14.C"
#include "TFile.h"
#include "THStack.h"
#include <iostream>
#include "TTree.h"
// #include "TH1.h"


using namespace std;

void ZjetPlot()
{
	TDirectory *curdir = gDirectory;
	TFile *f21 = new TFile("jets_output21.root","READ");
	assert(f21 && !f21->IsZombie());
	TTree *tree21 = (TTree*)f21->Get("tree");
	unsigned int N21 = (unsigned int)tree21->GetEntries(); 
	cout<<N21<<" events of 2->1 process."<<endl;
	float Z_pt_21;
	tree21->SetBranchAddress("Z_pT",&Z_pt_21);


	TFile *f22 = new TFile("jets_output22.root","READ");
	assert(f22 && !f22->IsZombie());
	TTree *tree22 = (TTree*)f22->Get("tree");
	unsigned int N22 = (unsigned int)tree22->GetEntries(); 
	cout<<N22<<" events of 2->2 process."<<endl;
	float Z_pt_22, weight;
	tree22->SetBranchAddress("Z_pT",&Z_pt_22);
	tree22->SetBranchAddress("weight",&weight);
	
	TH1D* Z21 = new TH1D("21sample","21sample",100,0,2200);
	Z21->Sumw2();
	TH1D* Z22 = new TH1D("22sample","22sample",100,0,2200);
	Z22->Sumw2();

	for(unsigned int x=0; x != N21; ++x)
	{
		tree21->GetEntry(x);
		Z21->Fill(Z_pt_21);
	}


	for(unsigned int x=0; x != N22; ++x)
	{
		tree22->GetEntry(x);
		Z22->Fill(Z_pt_22,weight);
	}

	Z22->Draw("SAMEE");
	
	Z21->SetLineColor(kRed);
	Z21->Scale(0.001);
	Z21->Draw("SAMEE");
	
	gPad->SetLogy();


}