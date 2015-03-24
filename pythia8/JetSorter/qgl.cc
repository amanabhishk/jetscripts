#include "tdrstyle_mod14.C"
#include "TFile.h"
#include "THStack.h"
#include <iostream>
#include "TTree.h"
#include "TProfile.h"


using namespace std;

void sigma()
{

	int ptBins = 48.;
	const double ptRange[]=
    {18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84,
    97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468,
    507, 548, 592, 638, 686, 737, 790, 846, 905, 967,
    1032, 1101, 1172, 1248, 1327, 1410, 1497, 1588, 1684, 1784, 1890, 2000};

	TDirectory *curdir = gDirectory;
	TFile *f = new TFile("d.root","READ");
	assert(f && !f->IsZombie());
	TTree *tree = (TTree*)f->Get("tree");
	unsigned int N = (unsigned int)tree->GetEntries(); 
	cout<<N<<" jets in total."<<endl;

	TFile *f1 = new TFile("d1.root","READ");
	assert(f1 && !f1->IsZombie());
	TTree *tree1 = (TTree*)f1->Get("tree");
	unsigned int N1 = (unsigned int)tree1->GetEntries(); 
	cout<<N<<" jets in total."<<endl;
	
	float pT, weight, pTD, sigma2[2];
	unsigned int constituents[2];
	unsigned char flavor;
	
	tree->SetBranchAddress("jet_pT",&pT);
	tree->SetBranchAddress("jet_sigma2",sigma2);
	tree->SetBranchAddress("jet_pTD",&pTD);
	tree->SetBranchAddress("jet_flavor",&flavor);
	tree->SetBranchAddress("jet_weight",&weight);
	tree->SetBranchAddress("jet_multiplicity",constituents);

	TH1D* sigma2_g = new TH1D("sigma2_g","sigma2_g",100,0,0.2);
	TH1D* sigma2_q = new TH1D("sigma2_q","sigma2_q",100,0,0.2);
	TH1D* sigma2_u = new TH1D("sigma2_u","sigma2_u",100,0,0.2);
	sigma2_u->Sumw2();
	sigma2_g->Sumw2();
	sigma2_q->Sumw2();
	sigma2_g->SetLineColor(kRed);
	sigma2_u->SetLineColor(kGreen);

	TH1D* pTD_g = new TH1D("pTD_g","pTD_g",100,0,1);
	TH1D* pTD_q = new TH1D("pTD_q","pTD_q",100,0,1);
	TH1D* pTD_u = new TH1D("pTD_u","pTD_u",100,0,1);
	pTD_u->Sumw2();
	pTD_g->Sumw2();
	pTD_q->Sumw2();
	pTD_g->SetLineColor(kRed);
	pTD_u->SetLineColor(kGreen);

	TH1D* multiplicity_g = new TH1D("multiplicity_g","multiplicity_g",60,0,60);
	TH1D* multiplicity_q = new TH1D("multiplicity_q","multiplicity_q",60,0,60);
	TH1D* multiplicity_u = new TH1D("multiplicity_u","multiplicity_u",60,0,60);
	multiplicity_q->Sumw2();
	multiplicity_g->Sumw2();
	multiplicity_u->Sumw2();
	multiplicity_g->SetLineColor(kRed);
	multiplicity_u->SetLineColor(kGreen);


	for(unsigned int x=0; x != N; ++x)
	{
		tree->GetEntry(x);

		if(pT>100 || pT<80) continue;
		if(flavor == 21) 
		{
			sigma2_g->Fill(sigma2[1],weight);
			multiplicity_g->Fill(constituents[1],weight);
			pTD_g->Fill(pTD,weight);
		}
		else if(flavor == 0) 
		{
			sigma2_u->Fill(sigma2[1],weight);
			multiplicity_u->Fill(constituents[1],weight);
			sigma2_u->Fill(sigma2[1],weight);	
		}
		else 
		{
			sigma2_q->Fill(sigma2[1],weight);
			multiplicity_q->Fill(constituents[1],weight);
			sigma2_q->Fill(sigma2[1],weight);	
		}
	}


	
	float pT1, weight1, pT1D, sigma21[2];
	unsigned int constituents1[2];
	unsigned char flavor1;
	
	tree1->SetBranchAddress("jet_pT",&pT1);
	tree1->SetBranchAddress("jet_sigma2",sigma21);
	tree1->SetBranchAddress("jet_pTD",&pT1D);
	tree1->SetBranchAddress("jet_flavor",&flavor1);
	tree1->SetBranchAddress("jet_weight",&weight1);
	tree1->SetBranchAddress("jet_multiplicity",constituents1);

	TH1D* sigma21_g = new TH1D("sigma21_g","sigma21_g",100,0,0.2);
	TH1D* sigma21_q = new TH1D("sigma21_q","sigma21_q",100,0,0.2);
	TH1D* sigma21_u = new TH1D("sigma21_u","sigma21_u",100,0,0.2);
	sigma21_u->Sumw2();
	sigma21_g->Sumw2();
	sigma21_q->Sumw2();
	sigma21_g->SetLineColor(kRed);
	sigma21_u->SetLineColor(kGreen);

	for(unsigned int x=0; x != N1; ++x)
	{
		tree1->GetEntry(x);

		if(pT1>100 || pT1<80) continue;
		if(flavor1 == 21) 
		{
			sigma21_g->Fill(sigma21[1],weight1);
		}
		else if(flavor1 == 0) 
		{
			sigma21_u->Fill(sigma21[1],weight1);	
		}
		else 
		{
			sigma21_q->Fill(sigma21[1],weight1);	
		}
	}
	/****************sigma2 with CUT****************/

	//SEPARATE
	setTDRStyle();
	TH1F *us = (TH1F*)sigma2_u->Clone("us");
	TH1F *gs = (TH1F*)sigma2_g->Clone("gs");
	TH1F *qs = (TH1F*)sigma2_q->Clone("qs");

	TH1F *us1 = (TH1F*)sigma21_u->Clone("us1");
	TH1F *gs1 = (TH1F*)sigma21_g->Clone("gs1");
	TH1F *qs1 = (TH1F*)sigma21_q->Clone("qs1");

	TH1D *h9 = new TH1D("h9",";#sigma_{2}(with cut);Events",100,0,0.2);
	h9->SetMinimum(0);
	h9->SetMaximum(0.08);
	h9->GetYaxis()->SetNoExponent();
	h9->GetXaxis()->SetNoExponent();
	h9->GetXaxis()->SetRangeUser(0,0.2);

	TCanvas *c9 = tdrCanvas("c9",h9,0,33);
	tdrDraw(gs,"HIST",kDot,kRed-3,kSolid,-1,3003,kRed-3);
	tdrDraw(qs,"HIST",kDot,kBlue,kSolid,-1,3005,kBlue);
	qs->Scale(1/qs->Integral());
	gs->Scale(1/gs->Integral());
	
	tdrDraw(gs1,"HIST",kDot,kRed-3,kSolid,-1,3003,kRed-3);
	tdrDraw(qs1,"HIST",kDot,kBlue,kSolid,-1,3005,kBlue);
	qs1->Scale(1/qs1->Integral());
	gs1->Scale(1/gs1->Integral());
}

