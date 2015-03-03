#include "tdrstyle_mod14.C"
#include "TFile.h"
#include "THStack.h"
#include <iostream>
#include "TTree.h"
#include "TProfile.h"


using namespace std;

void plot()
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
	
	float pT, weight, pTD, sigma2;
	unsigned int constituents;
	unsigned char flavor;
	
	tree->SetBranchAddress("jet_pT",&pT);
	tree->SetBranchAddress("jet_sigma2",&sigma2);
	tree->SetBranchAddress("jet_pTD",&pTD);
	tree->SetBranchAddress("jet_flavor",&flavor);
	tree->SetBranchAddress("jet_weight",&weight);
	tree->SetBranchAddress("jet_multiplicity",&constituents);
	
	TProfile gluonFrac("g","g",ptBins,ptRange);
  	TProfile lightquarkFrac("lq","lq",ptBins,ptRange);
  	TProfile charmFrac("c","c",ptBins,ptRange);
  	TProfile bottomFrac("b","b",ptBins,ptRange);
  	TProfile unmatchedFrac("unmatched","unmatched",ptBins,ptRange);

	TH1D* multiplicity_g = new TH1D("multiplicity_g","multiplicity_g",60,0,60);
	TH1D* multiplicity_q = new TH1D("multiplicity_q","multiplicity_q",60,0,60);
	TH1D* multiplicity_u = new TH1D("multiplicity_u","multiplicity_u",60,0,60);
	multiplicity_q->Sumw2();
	multiplicity_g->Sumw2();
	multiplicity_u->Sumw2();
	multiplicity_g->SetLineColor(kRed);
	multiplicity_u->SetLineColor(kGreen);
	

	TH1D* pTD_g = new TH1D("pTD_g","pTD_g",210,0,70);
	TH1D* pTD_q = new TH1D("pTD_q","pTD_q",210,0,70);
	TH1D* pTD_u = new TH1D("pTD_u","pTD_u",210,0,70);
	pTD_u->Sumw2();
	pTD_g->Sumw2();
	pTD_q->Sumw2();
	pTD_g->SetLineColor(kRed);
	pTD_u->SetLineColor(kGreen);

	TH1D* sigma2_g = new TH1D("sigma2_g","sigma2_g",210,0,70);
	TH1D* sigma2_q = new TH1D("sigma2_q","sigma2_q",210,0,70);
	TH1D* sigma2_u = new TH1D("sigma2_u","sigma2_u",210,0,70);
	sigma2_u->Sumw2();
	sigma2_g->Sumw2();
	sigma2_q->Sumw2();
	sigma2_g->SetLineColor(kRed);
	sigma2_u->SetLineColor(kGreen);


	for(unsigned int x=0; x != N; ++x)
	{
		tree->GetEntry(x);
		
		gluonFrac.Fill(pT, (flavor == 21)? 1:0, weight);
    	lightquarkFrac.Fill(pT, (flavor == 1 || flavor == 2 || flavor == 3)? 1:0, weight);
    	charmFrac.Fill(pT, (flavor == 4)? 1:0, weight);
    	bottomFrac.Fill(pT, (flavor == 5)? 1:0, weight);
    	unmatchedFrac.Fill(pT, (flavor == 0)? 1:0, weight);

		if(pT>100 || pT<80) continue;
		if(flavor == 21) 
		{
			multiplicity_g->Fill(constituents,weight);
			pTD_g->Fill(pTD,weight);
			sigma2_g->Fill(sigma2,weight);
		}
		else if(flavor == 0) 
		{
			multiplicity_u->Fill(constituents,weight);
			pTD_u->Fill(pTD,weight);
			sigma2_u->Fill(sigma2,weight);	
		}
		else 
		{
			multiplicity_q->Fill(constituents,weight);
			pTD_q->Fill(pTD,weight);
			sigma2_q->Fill(sigma2,weight);	
		}
	}

	TH1D *light_quarks = lightquarkFrac.ProjectionX("light quarks","");
  	TH1D *gluons = gluonFrac.ProjectionX("gluons","");
  	TH1D *charm = charmFrac.ProjectionX("charm","");
  	TH1D *bottom = bottomFrac.ProjectionX("bottom","");
  	TH1D *unmatched = unmatchedFrac.ProjectionX("unmatched","");
	
  	TH1D *h = new TH1D("h",";p_{T} (GeV);Fraction",1000,20,2000);

	tdrDraw(unmatched,"",kOpenCircle,kGray+2,kSolid,-1,1001,kGray);
	tdrDraw(gluons,"",kPlus,kBlue+2,kSolid,-1,1001,kBlue-9);
	tdrDraw(light_quarks,"",kFullCircle,kGreen-1,kSolid,-1,1001,kYellow-9);
	tdrDraw(charm,"",kFullCircle,kGreen-1,kSolid,-1,1001,kGreen-9);
	tdrDraw(bottom,"",kFullCircle,kRed-2,kSolid,-1,1001,kRed-9);

	THStack *hs  = new THStack("hs","test stacked histograms");

	TCanvas *c1 = tdrCanvas("c1",h,2,0,kSquare);

	hs->Add(bottom);
	hs->Add(charm);
	hs->Add(light_quarks);
	hs->Add(gluons);
	hs->Add(unmatched);

	hs->Draw();

	hs->GetXaxis()->SetNoExponent();
	// hs->GetXaxis()->SetNoExponent(kTRUE);
	hs->GetXaxis()->SetRange(9,47);
	hs->GetXaxis()->SetNdivisions(5,kTRUE);
	hs->GetXaxis()->SetTitle("p_{T} (GeV)");
	hs->GetYaxis()->SetTitle("Flavor fraction");
	// hs->SetLogx();
	TLegend *leg = new TLegend(0.175,0.50,0.5,0.78);


	leg->SetFillStyle(kNone);
	leg->SetBorderSize(0);
 	leg->SetTextSize(0.045);

	
	leg->AddEntry(bottom,"Bottom","f");
	leg->AddEntry(charm,"Charm","f");
	leg->AddEntry(light_quarks,"Light","f");
	leg->AddEntry(gluons,"Gluon","f");
	leg->AddEntry(unmatched,"None","f");
	
	leg->Draw();
	
	gPad->SetLogx();
	
	/////////////////////////////////////

	setTDRStyle();
	TH1D *h1 = new TH1D("h1",";Number of constituents;Events",100,0,1000);
	h1->SetMinimum(0);
	//h1->SetMaximum(1);
	h1->GetYaxis()->SetNoExponent();
	h1->GetXaxis()->SetNoExponent();
	//h1->GetXaxis()->SetMoreLogLabels(kTRUE);
	h1->GetXaxis()->SetRangeUser(0,60);
	// h1->GetYaxis()->SetMoreLogLabels(kTRUE);
	// gStyle->SetOptStat(kFALSE);
	TCanvas *c2 = tdrCanvas("c2",h1,0,33);
	tdrDraw(multiplicity_g,"P",kDot,kRed-9,kSolid,-1,1001,kBlue-6);
	tdrDraw(multiplicity_u,"P",kDot,kGray,kSolid,-1,1001,kBlue-6);
	tdrDraw(multiplicity_q,"P",kDot,kBlue,kSolid,-1,1001,kBlue-6);
	multiplicity_g->Scale(1/multiplicity_g->Integral());
	multiplicity_u->Scale(1/multiplicity_u->Integral());
	multiplicity_q->Scale(1/multiplicity_q->Integral());	
	////////////////////////////////////

	setTDRStyle();
	TH1D *h2 = new TH1D("h2",";p_{T}D;Events",100,0,1000);
	//h2->SetMinimum(0);
	//h2->SetMaximum(1);
	h2->GetYaxis()->SetNoExponent();
	h2->GetXaxis()->SetNoExponent();
	//h2->GetXaxis()->SetMoreLogLabels(kTRUE);
	h2->GetXaxis()->SetRangeUser(0,1);
	// h2->GetYaxis()->SetMoreLogLabels(kTRUE);
	// gStyle->SetOptStat(kFALSE);
	TCanvas *c3 = tdrCanvas("c3",h2,0,33);
	tdrDraw(pTD_g,"P",kDot,kRed-9,kSolid,-1,1001,kBlue-6);
	tdrDraw(pTD_u,"P",kDot,kGray,kSolid,-1,1001,kBlue-6);
	tdrDraw(pTD_q,"P",kDot,kBlue,kSolid,-1,1001,kBlue-6);
	pTD_g->Scale(1/pTD_g->Integral());
	pTD_u->Scale(1/pTD_u->Integral());
	pTD_q->Scale(1/pTD_q->Integral());

	////////////////////////////////////
	setTDRStyle();
	TH1D *h3 = new TH1D("h3",";#sigma_{2};Events",100,0,1000);
	//h3->SetMinimum(0);
	//h3->SetMaximum(1);
	h3->GetYaxis()->SetNoExponent();
	h3->GetXaxis()->SetNoExponent();
	//h3->GetXaxis()->SetMoreLogLabels(kTRUE);
	h3->GetXaxis()->SetRangeUser(0,0.2);
	// h3->GetYaxis()->SetMoreLogLabels(kTRUE);
	// gStyle->SetOptStat(kFALSE);
	TCanvas *c4 = tdrCanvas("c4",h3,0,33);
	tdrDraw(sigma2_g,"P",kDot,kRed-9,kSolid,-1,1001,kBlue-6);
	tdrDraw(sigma2_u,"P",kDot,kGray,kSolid,-1,1001,kBlue-6);
	tdrDraw(sigma2_q,"P",kDot,kBlue,kSolid,-1,1001,kBlue-6);
	sigma2_q->Scale(1/sigma2_q->Integral());
	sigma2_u->Scale(1/sigma2_u->Integral());
	sigma2_g->Scale(1/sigma2_g->Integral());
}

