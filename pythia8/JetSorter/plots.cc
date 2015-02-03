#include "tdrstyle_mod14.C"
#include "TFile.h"
#include "THStack.h"
#include <iostream>


using namespace std;

void plots()
{
	TDirectory *curdir = gDirectory;
	
	TFile *fin = new TFile("output.root","READ");
	assert(fin && !fin->IsZombie());
	curdir->cd();
	// 
	setTDRStyle();
	TH1D *h2 = (TH1D*)fin->Get("gluons");
	assert(h2);
	h2->UseCurrentStyle();
	TH1D *h3 = (TH1D*)fin->Get("light quarks");
	assert(h3);
	h3->UseCurrentStyle();
	TH1D *h4 = (TH1D*)fin->Get("charm");
	assert(h4);
	h4->UseCurrentStyle();
	
	TH1D *h5 = (TH1D*)fin->Get("bottom");
	assert(h5);
	h5->UseCurrentStyle();
	
	TH1D *h1 = (TH1D*)fin->Get("unmatched");
	assert(h1);
	h1->UseCurrentStyle();
	// TH1D *h = new TH1D("h","");
	h2->Add(h1);
	
	TH1D *h = new TH1D("h",";p_{T} (GeV);Fraction",1000,20,2000);
	// h->SetMinimum(0);
 //    h->SetMaximum(1.2);
    //h1->GetXaxis()->SetMoreLogLabels();

	// h1->Draw();
	
	// // curdir->cd();
	// // tdrDraw(h,"");

	tdrDraw(h1,"",kOpenCircle,kOrange+7,kSolid,-1,1001,kOrange-3);
	tdrDraw(h2,"",kFullCircle,kBlue+2,kSolid,-1,1001,kBlue-9);
	tdrDraw(h3,"",kFullCircle,kGreen-1,kSolid,-1,1001,kGreen-8);
	tdrDraw(h4,"",kFullCircle,kYellow-2,kSolid,-1,1001,kYellow-9);
	tdrDraw(h5,"",kFullCircle,kMagenta-1,kSolid,-1,1001,kMagenta-9);
	// // TLegend *leg = new TLegend(0.25,0.25,0.55,0.30);
	// leg->AddEntry(h5,"test","w");
	// leg->SetHeader("The Legend Title");

	THStack *hs  = new THStack("hs","test stacked histograms");
    
    // hs->GetXaxis()->SetNoExponent();
	

	TCanvas *c1 = tdrCanvas("c1",h,2,0,kSquare);
	

	// h2->Draw();
	hs->Add(h5);
	hs->Add(h4);
	hs->Add(h3);
	hs->Add(h2);
	// hs->Add(h1);

	// leg->Draw();
    // h1->GetXaxis()->SetLogx(1);
    //h1->GetYaxis()->SetMoreLogLabels();
    // hs->GetYaxis()->SetNoExponent();

	// hs->Add(h1);
	// leg->Draw();
	
	
	// TCanvas c2("c2","stacked hists",10,10,700,900);

	hs->Draw();
	// tdrDraw(h2,"");

	gPad->SetLogx();
	gPad->RedrawAxis();
	// gStyle->SetPalette(1);
    // gStyle->SetOptStat(0);
	// gPad->Update();
}

