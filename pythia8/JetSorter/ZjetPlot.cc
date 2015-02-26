#include "tdrstyle_mod14.C"
#include "TFile.h"
#include "THStack.h"
#include <iostream>
#include "TTree.h"
// #include "TH1.h"


using namespace std;

void ZjetPlot()
{
	int ptBins = 48.;
	const double ptRange[]=
    {18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84,
    97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468,
    507, 548, 592, 638, 686, 737, 790, 846, 905, 967,
    1032, 1101, 1172, 1248, 1327, 1410, 1497, 1588, 1684, 1784, 1890, 2000};

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
	
	TH1D* Z21 = new TH1D("21sample","21sample",ptBins,ptRange);
	Z21->Sumw2();
	TH1D* Z22 = new TH1D("22sample","22sample",ptBins,ptRange);
	Z22->Sumw2();
	TH1D* ratio = new TH1D("ratio","ratio",ptBins,ptRange);
	ratio->Sumw2();

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
	ratio->Divide(Z22,Z21,1,0.001);


	// TCanvas *c1 = new TCanvas("c1", "c1",900, 700);
	// c1->Divide(1,2);
	
	// c1->cd(1);
	// gPad->SetLogy();
	// Z22->Draw("SAMEE");

	// Z21->Draw("SAMEE");
	// Z21->SetLineColor(kRed);
	Z21->Scale(0.0009);
	
	// c1->cd(2);
	// gPad->SetLogy();
	
	// ratio->Draw();
	// ratio->GetXaxis()->SetRange(0,300);
	// ratio->SetLineColor(kGreen+3);

	// 
	setTDRStyle();
	TH1D *h1 = new TH1D("h1",";p_{T} (GeV);Cross section",100,0,1000);
	h1->SetMinimum(2*1e-4);
	h1->SetMaximum(2*1e2);
	h1->GetYaxis()->SetNoExponent();
	h1->GetXaxis()->SetNoExponent();
	h1->GetXaxis()->SetMoreLogLabels(kTRUE);
	h1->GetXaxis()->SetRangeUser(20,1000);
	// h1->GetYaxis()->SetMoreLogLabels(kTRUE);
	TH1D *h2 = new TH1D("h2",";p_{T} (GeV);(241+242)/221",100,0,1000);
	h2->SetMaximum(1.6);
	h2->SetMinimum(0.3);
	h2->GetXaxis()->SetNoExponent();
	h2->GetYaxis()->SetNoExponent();
	h2->GetXaxis()->SetMoreLogLabels(kTRUE);
	h2->GetXaxis()->SetRangeUser(20,1000);
	// h2->GetYaxis()->SetMoreLogLabels(kTRUE);
	TCanvas *c2 = tdrDiCanvas("c2",h1,h2,0,33);
	c2->cd(1);
	gPad->SetLogy();
	gPad->SetLogx();
	// gStyle->SetOptStat(kFALSE);
	tdrDraw(Z22,"HISTE",kDot,kBlue,kSolid,-1,1001,kBlue-6);
	tdrDraw(Z21,"PE",kMultiply, kBlue+4);
	
	c2->cd(2);
	gPad->SetLogx();
	tdrDraw(ratio,"P",kStar,kGreen+4);

}