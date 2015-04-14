#include "tdrstyle_mod14.C"
#include "TFile.h"
#include "THStack.h"
#include <iostream>
#include "TTree.h"
#include "TProfile.h"
#include <string>
#include <sstream>

using namespace std;

void cont()
{

	bool merge = false;
	int ptBins = 48.;
	const double ptRange[]=
    {18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84,
    97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468,
    507, 548, 592, 638, 686, 737, 790, 846, 905, 967,
    1032, 1101, 1172, 1248, 1327, 1410, 1497, 1588, 1684, 1784, 1890, 2000};

	TDirectory *curdir = gDirectory;
	TFile *f = new TFile("5000dijet_feynman.root","READ");
	assert(f && !f->IsZombie());
	TTree *tree = (TTree*)f->Get("Events");
	unsigned int N = (unsigned int)tree->GetEntries(); 
	cout<<N<<" events."<<endl;
	
	int in, out;               //exact number of particles stored for each event 
  	float pt1, pt2, wt;
	
	tree->SetBranchAddress("pT1",&pt1);
	tree->SetBranchAddress("pT2",&pt2);
	tree->SetBranchAddress("weight",&wt);
	tree->SetBranchAddress("outType",&out);
	tree->SetBranchAddress("inType",&in);

	vector<TProfile*> diag(0);
	
	for(unsigned int x = 0; x < 7; x++)
	{
		for(unsigned int y = 0; y < 7; y++)
		{
			std::stringstream temp("");
			temp<<x<<y;
			diag.push_back(new TProfile(temp.str().c_str(),temp.str().c_str(),ptBins,ptRange));
		}
	}

	for(unsigned int x=0; x != N; ++x)
	{
		tree->GetEntry(x);
		
		for(unsigned int k = 0; k < diag.size(); k++)
		{
			if(abs(10*in+out) != k) 
			{
				diag[k]->Fill(pt1,0,wt);
				diag[k]->Fill(pt2,0,wt);
			}
			else 
			{
				diag[k]->Fill(pt1,1,wt);
				diag[k]->Fill(pt2,1,wt);
			}
		}
	}
	cout<<"check";
	vector<TH1D*> diag_proj(0);
	
	for(unsigned int x = 0; x < 7; x++)
	{
		for(unsigned int y = 0; y < 7; y++)
		{
			std::stringstream t("");
			t<<x<<y;
			TProfile* god = diag[x*y];
			diag_proj.push_back(god->ProjectionX(t.str().c_str(),""));
		}
	}

	// TH1D *light_quarks = lightquarkFrac.ProjectionX("light quarks","");
 //  	TH1D *gluons = gluonFrac.ProjectionX("gluons","");
 //  	TH1D *strange = strangeFrac.ProjectionX("strange","");
 //  	if(merge) light_quarks->Add(strange);
 //  	TH1D *charm = charmFrac.ProjectionX("charm","");
 //  	TH1D *bottom = bottomFrac.ProjectionX("bottom","");
 //  	TH1D *unmatched = unmatchedFrac.ProjectionX("unmatched","");
	
   	TH1D *h = new TH1D("h",";p_{T} (GeV);Fraction",1000,20,2000);

	for(unsigned int a = 0; a < diag_proj.size(); ++a) 
	{
		tdrDraw(diag_proj[a],"",kOpenCircle,kGray+2,kSolid,-1,1001,a%5);
	}		
	// gStyle->SetOptStat(kFALSE); //removes old legend
	// tdrDraw(gluons,"",kPlus,kBlue+2,kSolid,-1,1001,kBlue-9);
	// tdrDraw(light_quarks,"",kFullCircle,kGreen-1,kSolid,-1,1001,kYellow-9);
	// tdrDraw(strange,"",kFullCircle,kAzure-6,kSolid,-1,1001,kAzure-8);
	// tdrDraw(charm,"",kFullCircle,kGreen-1,kSolid,-1,1001,kGreen-9);
	// tdrDraw(bottom,"",kFullCircle,kRed-2,kSolid,-1,1001,kRed-9);

	THStack *hs  = new THStack("hs","test stacked histograms");

	TCanvas *c1 = tdrCanvas("c1",h,2,0,kSquare);

	for(unsigned int a = 0; a < diag_proj.size(); ++a) 
	{
		hs->Add(diag_proj[a]);
	}
	// //light_quarks->Add(strange);
	// hs->Add(bottom);
	// hs->Add(charm);
	// if(!merge)hs->Add(strange);
	// hs->Add(light_quarks);
	// hs->Add(gluons);
	// hs->Add(unmatched);

	hs->Draw();

	hs->GetXaxis()->SetNoExponent();
	// hs->GetXaxis()->SetNoExponent(kTRUE);
	hs->GetXaxis()->SetRange(9,47);
	hs->GetXaxis()->SetMoreLogLabels(kTRUE);
	//hs->GetXaxis()->SetNdivisions(5,kTRUE);
	hs->GetXaxis()->SetTitle("p_{T} (GeV)");
	hs->GetYaxis()->SetTitle("Flavor fraction");
	// hs->SetLogx();

	// double x0, y0;
	// x0 = 0.45;
	// y0 = 0.05;
	// //TLegend *leg = tdrLeg(0.5,0.82-0.1,0.175,0.50-0.1); 				//physics def
	// //TLegend *leg = tdrLeg(0.5,0.82+0.07,0.175,0.50+0.07);				//hadronic def
	// TLegend *leg = tdrLeg(0.5+0.5,0.82-0.2,0.175+0.5,0.50-0.2);		//QCDaware def
	

	
	// TLegend *sample = tdrLeg(0.675-0.05-x0,0.50-y0,0.775-0.05-x0,0.505-y0);				//QCDaware
	// TLegend *alphacut = tdrLeg(0.77-x0,0.50-0.05-y0,0.87-x0,0.505-0.05-y0);				//goes
	// TLegend *etacut = tdrLeg(0.61-x0,0.50-0.05-y0,0.71-x0,0.505-0.05-y0);				//here

	// TLegend *heading = tdrLeg(0.675-0.4,0.50+0.44,0.775-0.4,0.505+0.44); 	
	// // TLegend *sample = tdrLeg(0.675-0.05,0.50,0.775-0.05,0.505);			//everything except
	// // TLegend *alphacut = tdrLeg(0.77,0.50-0.05,0.87,0.505-0.05);			//QCD
	// // TLegend *etacut = tdrLeg(0.61,0.50-0.05,0.71,0.505-0.05);			//aware

	// sample->SetHeader("Z+jet sample");
	// heading->SetHeader("Hadronic Definition, #sqrt{s} = 8 TeV");
	// alphacut->SetHeader("#alpha<0.3");
	// etacut->SetHeader("#left|#eta#right|< 1.3,");

	
	// leg->AddEntry(bottom,"Bottom","f");
	// leg->AddEntry(charm,"Charm","f");
	// leg->AddEntry(strange,"Strange","f");
	// leg->AddEntry(light_quarks,"Light","f");
	// leg->AddEntry(gluons,"Gluon","f");
	// leg->AddEntry(unmatched,"None","f");
	
	// leg->Draw();
	// heading->Draw();
	// sample->Draw();
	// alphacut->Draw();
	// etacut->Draw();
	
	
	gPad->SetLogx();       
	// //c1->SaveAs("output.pdf");
}

