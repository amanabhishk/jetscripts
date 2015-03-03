#include "tdrstyle_mod14.C"
#include "TFile.h"
#include "THStack.h"
#include <iostream>
#include "TTree.h"
// #include "TH1.h"


using namespace std;

void dijetPlot()
{
	int ptBins = 48.;
	const double ptRange[]=
    {18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84,
    97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468,
    507, 548, 592, 638, 686, 737, 790, 846, 905, 967,
    1032, 1101, 1172, 1248, 1327, 1410, 1497, 1588, 1684, 1784, 1890, 2000};

	TDirectory *curdir = gDirectory;
	TFile *f = new TFile("jets_outputL.root","READ");
	assert(f && !f->IsZombie());
	TTree *tree = (TTree*)f->Get("tree");
	unsigned int N = (unsigned int)tree->GetEntries(); 
	cout<<N<<" jets in total."<<endl;
	
	float pT, weight, pTD;
	unsigned int constituents = 99;
	unsigned char flavor;
	
	tree->SetBranchAddress("jet_pT",&pT);
	tree->SetBranchAddress("jet_pTD",&pTD);
	tree->SetBranchAddress("jet_flavor",&flavor);
	tree->SetBranchAddress("jet_weight",&weight);
	tree->SetBranchAddress("jet_multiplicity",&constituents);
	
	TH1D* multiplicity_g = new TH1D("multiplicity_g","multiplicity_g",60,0,60);
	TH1D* multiplicity_q = new TH1D("multiplicity_q","multiplicity_q",60,0,60);
	TH1D* multiplicity_u = new TH1D("multiplicity_u","multiplicity_u",60,0,60);
	TH1D* pTD_g = new TH1D("pTD_g","pTD_g",210,0,70);
	TH1D* pTD_q = new TH1D("pTD_q","pTD_q",210,0,70);
	TH1D* pTD_u = new TH1D("pTD_u","pTD_u",210,0,70);
	multiplicity_q->Sumw2();
	multiplicity_g->SetLineColor(kRed);
	multiplicity_u->SetLineColor(kGreen);
	multiplicity_g->Sumw2();
	multiplicity_u->Sumw2();
	pTD_u->Sumw2();
	pTD_g->Sumw2();
	pTD_q->Sumw2();
	pTD_g->SetLineColor(kRed);
	pTD_u->SetLineColor(kGreen);
	//TH1D* ratio = new TH1D("ratio","ratio",ptBins,ptRange);
	//ratio->Sumw2();

	for(unsigned int x=0; x != N; ++x)
	{
		tree->GetEntry(x);
		if(pT>100 || pT<80) continue;
		//cout<<constituents<<endl;
		//multiplicity_g->Fill(constituents);
		if(flavor == 21) multiplicity_g->Fill(constituents,weight), pTD_g->Fill(pTD,weight);
		else if(flavor == 0) multiplicity_u->Fill(constituents,weight),pTD_u->Fill(pTD,weight);
		else multiplicity_q->Fill(constituents,weight),pTD_q->Fill(pTD,weight);
	}

	//multiplicity_g->Scale(1/multiplicity_g->Integral());
	//multiplicity_q->Scale(1/multiplicity_q->Integral());
	
	multiplicity_g->Draw("HISTSAME");
	multiplicity_q->Draw("HISTSAME");
	multiplicity_u->Draw("HISTSAME");

	//pTD_g->Draw("HIST");
	//pTD_q->Draw("HISTSAME");
	//pTD_u->Draw("HISTSAME");
	
	// TCanvas *c1 = new TCanvas("c1", "c1",900, 700);
	// c1->Divide(1,2);
	
	// c1->cd(1);
	// gPad->SetLogy();
	// Z22->Draw("SAMEE");

	// Z21->Draw("SAMEE");
	// Z21->SetLineColor(kRed);
	
	// c1->cd(2);
	// gPad->SetLogy();
	
	// ratio->Draw();
	// ratio->GetXaxis()->SetRange(0,300);
	// ratio->SetLineColor(kGreen+3);

	// 
	// setTDRStyle();
	// TH1D *h1 = new TH1D("h1",";p_{T} (GeV);Cross section",100,0,1000);
	// h1->SetMinimum(2*1e-4);
	// h1->SetMaximum(2*1e2);
	// h1->GetYaxis()->SetNoExponent();
	// h1->GetXaxis()->SetNoExponent();
	// h1->GetXaxis()->SetMoreLogLabels(kTRUE);
	// h1->GetXaxis()->SetRangeUser(20,1000);
	// // h1->GetYaxis()->SetMoreLogLabels(kTRUE);
	// TH1D *h2 = new TH1D("h2",";p_{T} (GeV);(241+242)/221",100,0,1000);
	// h2->SetMaximum(1.6);
	// h2->SetMinimum(0.3);
	// h2->GetXaxis()->SetNoExponent();
	// h2->GetYaxis()->SetNoExponent();
	// h2->GetXaxis()->SetMoreLogLabels(kTRUE);
	// h2->GetXaxis()->SetRangeUser(20,1000);
	// // h2->GetYaxis()->SetMoreLogLabels(kTRUE);
	// TCanvas *c2 = tdrDiCanvas("c2",h1,h2,0,33);
	// c2->cd(1);
	// gPad->SetLogy();
	// gPad->SetLogx();
	// // gStyle->SetOptStat(kFALSE);
	// tdrDraw(Z22,"HISTE",kDot,kBlue,kSolid,-1,1001,kBlue-6);
	// tdrDraw(Z21,"PE",kMultiply, kBlue+4);
	
	// c2->cd(2);
	// gPad->SetLogx();
	// tdrDraw(ratio,"P",kStar,kGreen+4);

}