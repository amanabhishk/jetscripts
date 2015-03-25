#include "tdrstyle_mod14.C"
#include "TFile.h"
#include "THStack.h"
#include <iostream>
#include "TTree.h"
#include "TLorentzVector.h"
// #include "TH1.h"


using namespace std;

void plot_Zjet_channels()
{
	unsigned int size = 2000;
	int ptBins = 48+6;
	const double ptRange[]=
    {0,5,8,10,12,15,18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84,
    97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468,
    507, 548, 592, 638, 686, 737, 790, 846, 905, 967,
    1032, 1101, 1172, 1248, 1327, 1410, 1497, 1588, 1684, 1784, 1890, 2000};

	TDirectory *curdir = gDirectory;
		
	TH1D* Z21 = new TH1D("21sample","21sample",ptBins,ptRange);
	Z21->Sumw2();
	TH1D* Z22 = new TH1D("22sample","22sample",ptBins,ptRange);
	Z22->Sumw2();
	TH1D* ratio = new TH1D("ratio","ratio",ptBins,ptRange);
	ratio->Sumw2();

	TFile *f21 = new TFile("40001events_Zjet.root","READ");
	assert(f21 && !f21->IsZombie());
	TTree *tree21 = (TTree*)f21->Get("Events");
	unsigned int N21 = (unsigned int)tree21->GetEntries(); 
	cout<<N21<<" events of 2->1 process."<<endl;

  	unsigned short int eventParticleCount;
  	int id[size];
  	float weight, pT[size], eta[size], phi[size], m[size];
  	unsigned char status[size];
  
  	tree21->SetBranchAddress("n", &eventParticleCount);
  	tree21->SetBranchAddress("weight", &weight);
  	tree21->SetBranchAddress("id", id);
  	tree21->SetBranchAddress("status", status);
  	tree21->SetBranchAddress("pT", pT);
  	tree21->SetBranchAddress("eta", eta);
  	tree21->SetBranchAddress("phi", phi);
  	tree21->SetBranchAddress("m", m);

	TLorentzVector v1;    
	for(unsigned int x=0; x != N21; ++x)
	{
		tree21->GetEntry(x);
		TLorentzVector v;
		for(unsigned int y = 0; y< eventParticleCount; ++y)
		{
			if(status[y]==2)
			{
				v1.SetPtEtaPhiM(pT[y],eta[y],phi[y],m[y]);
				v += v1;
			}
		}
		Z21->Fill(v.Pt(),weight);
		//cout<<v.Pt()*weight<<endl;
	}

	TFile *f22 = new TFile("40000events_Zjet.root","READ");
	assert(f22 && !f22->IsZombie());
	TTree *tree22 = (TTree*)f22->Get("Events");
	unsigned int N22 = (unsigned int)tree22->GetEntries(); 
	cout<<N22<<" events of 2->1 process."<<endl;
  
  	tree22->SetBranchAddress("n", &eventParticleCount);
  	tree22->SetBranchAddress("weight", &weight);
  	tree22->SetBranchAddress("id", id);
  	tree22->SetBranchAddress("status", status);
  	tree22->SetBranchAddress("pT", pT);
  	tree22->SetBranchAddress("eta", eta);
  	tree22->SetBranchAddress("phi", phi);
  	tree22->SetBranchAddress("m", m);

 
	for(unsigned int x=0; x != N22; ++x)
	{
		tree22->GetEntry(x);
		TLorentzVector v2;
		for(unsigned int y = 0; y< eventParticleCount; ++y)
		{
			if(status[y]==2)
			{
				v1.SetPtEtaPhiM(pT[y],eta[y],phi[y],m[y]);
				v2 += v1;
			}
		}
		Z22->Fill(v2.Pt(),weight);
		
	}


	// TCanvas *c1 = new TCanvas("c1", "c1",900, 700);
	// c1->Divide(1,2);
	
	// c1->cd(1);
	// gPad->SetLogy();
	// Z22->Draw("SAMEE");

	// Z21->Draw("SAMEE");
	// Z21->SetLineColor(kRed);
	ratio->Divide(Z22,Z21,1,0.001);
	Z21->Scale(0.0009);
	//Z22->Draw();
	
	setTDRStyle();
	TH1D *h1 = new TH1D("h1",";p_{T} (GeV);Cross section",100,0,1000);
	h1->SetMinimum(0.5*1e-7);
	h1->SetMaximum(1*1e1);
	//h1->GetYaxis()->SetNoExponent();
	//h1->GetXaxis()->SetNoExponent();
	//h1->GetXaxis()->SetMoreLogLabels(kTRUE);
	h1->GetXaxis()->SetRangeUser(10,1000);
	//h1->GetXaxis()->LabelsOption("v");
	// h1->GetYaxis()->SetMoreLogLabels(kTRUE);
	TH1D *h2 = new TH1D("h2",";p_{T} (GeV);(241+242)/221",100,0,1000);
	h2->SetMaximum(1.6);
	h2->SetMinimum(0.3);
	h2->GetXaxis()->SetNoExponent();
	h2->GetYaxis()->SetNoExponent();
	h2->GetXaxis()->SetMoreLogLabels(kTRUE);
	h2->GetXaxis()->SetRangeUser(10,1000);
	// h2->GetYaxis()->SetMoreLogLabels(kTRUE);
	TCanvas *c2 = tdrDiCanvas("c2",h1,h2,0,33);
	c2->cd(1);
	gPad->SetLogy();
	gPad->SetLogx();
	// gStyle->SetOptStat(kFALSE);
	tdrDraw(Z22,"PE",kFullCircle,kBlue);
	tdrDraw(Z21,"PE",kFullTriangleUp, kRed-6);
	
	c2->cd(2);
	gPad->SetLogx();
	tdrDraw(ratio,"P",kStar,kGreen+4);

}