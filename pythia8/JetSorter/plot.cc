#include "tdrstyle_mod14.C"
#include "TFile.h"
#include "THStack.h"
#include <iostream>


using namespace std;

void plot(bool merge = false)
{
	// if(merge) cout<< "works!\n";
	TDirectory *curdir = gDirectory;
	
	TFile *fin = new TFile("output.root","READ");
	assert(fin && !fin->IsZombie());
	curdir->cd();
	// 
	TH1D *unmatched = (TH1D*)fin->Get("unmatched");
	assert(unmatched);
	unmatched->UseCurrentStyle();

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
	

	// TH1D *h = new TH1D("h","");
	if(merge) gluons->Add(unmatched);
	
	TH1D *h = new TH1D("h",";p_{T} (GeV);Fraction",1000,20,2000);
	// h->SetMinimum(0);
 //    h->SetMaximum(1.2);
    //unmatched->GetXaxis()->SetMoreLogLabels();

	// unmatched->Draw();
	
	// // curdir->cd();
	// // tdrDraw(h,"");

	tdrDraw(unmatched,"",kOpenCircle,kOrange+7,kSolid,-1,1001,kOrange-3);
	tdrDraw(gluons,"",kFullCircle,kBlue+2,kSolid,-1,1001,kBlue-9);
	tdrDraw(light_quarks,"",kFullCircle,kGreen-1,kSolid,-1,1001,kGreen-8);
	tdrDraw(charm,"",kFullCircle,kYellow-2,kSolid,-1,1001,kYellow-9);
	tdrDraw(bottom,"",kFullCircle,kMagenta-1,kSolid,-1,1001,kMagenta-9);
	// // TLegend *leg = new TLegend(0.25,0.25,0.55,0.30);
	// leg->AddEntry(bottom,"test","w");
	// leg->SetHeader("The Legend Title");

	THStack *hs  = new THStack("hs","test stacked histograms");
    
    // hs->GetXaxis()->SetNoExponent();
	

	TCanvas *c1 = tdrCanvas("c1",h,2,0,kSquare);
	 // TCanvas *c1 = new TCanvas("c1","multipads",900,700);
	// c1->Divide(1,0,0.2,0.2);
	
	// c1->cd(1);
	// gluons->Draw();
	hs->Add(bottom);
	hs->Add(charm);
	hs->Add(light_quarks);
	hs->Add(gluons);
	if(!merge) hs->Add(unmatched);

	// leg->Draw();
    // unmatched->GetXaxis()->SetLogx(1);
    //unmatched->GetYaxis()->SetMoreLogLabels();
    // hs->GetYaxis()->SetNoExponent();

	// hs->Add(unmatched);
	// leg->Draw();
	
	
	// TCanvas c2("c2","stacked hists",10,10,700,900);

	hs->Draw();
	// tdrDraw(gluons,"");

	gPad->SetLogx();
	gPad->RedrawAxis();

	// c2->cd(2);
	// unmatched->Draw();
	// gStyle->SetPalette(1);
    // gStyle->SetOptStat(0);
	// gPad->Update();
}

