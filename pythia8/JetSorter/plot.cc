#include "tdrstyle_mod14.C"
#include "TFile.h"
#include "THStack.h"
#include <iostream>


using namespace std;

void plot()
{
	// if(merge) cout<< "works!\n";
	TDirectory *curdir = gDirectory;
	
	TFile *fin = new TFile("jets_output.root","READ");
	assert(fin && !fin->IsZombie());
	curdir->cd();
	// 
	TH1D *unmatched = (TH1D*)fin->Get("unmatched");
	assert(unmatched);
	unmatched->UseCurrentStyle();
	// unmatched->

	setTDRStyle();
	TH1D *gluons = (TH1D*)fin->Get("gluons");
	assert(gluons);
	gluons->UseCurrentStyle();

	TH1D *light_quarks = (TH1D*)fin->Get("light quarks");
	assert(light_quarks);
	light_quarks->UseCurrentStyle();

	TH1D *charm = (TH1D*)fin->Get("charm");
	assert(charm);
	charm->UseCurrentStyle();
	
	TH1D *bottom = (TH1D*)fin->Get("bottom");
	assert(bottom);
	bottom->UseCurrentStyle();
	
	
	TH1D *h = new TH1D("h",";p_{T} (GeV);Fraction",1000,20,2000);

	tdrDraw(unmatched,"",kOpenCircle,kGray+2,kSolid,-1,1001,kGray);
	tdrDraw(gluons,"",kPlus,kBlue+2,kSolid,-1,1001,kBlue-9);
	tdrDraw(light_quarks,"",kFullCircle,kGreen-1,kSolid,-1,1001,kYellow-9);
	tdrDraw(charm,"",kFullCircle,kGreen-1,kSolid,-1,1001,kGreen-9);
	tdrDraw(bottom,"",kFullCircle,kRed-2,kSolid,-1,1001,kRed-9);
	
	// gluons->ResetAttMarker();
	// gluons->SetMarkerStyle(20);

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
	// tdrDraw(gluons,"");

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
	// TPad->GetXAxis()->SetNoExponent();
	// gPad->RedrawAxis();

	// c2->cd(2);
	// unmatched->Draw();
	// gStyle->SetPalette(1);
    // gStyle->SetOptStat(0);
	// gPad->Update();
}

