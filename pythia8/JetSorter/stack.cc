#include "tdrstyle_mod14.C"
#include "TFile.h"
#include "THStack.h"
#include <iostream>
#include "TTree.h"
#include "TProfile.h"


using namespace std;

void stack()
{

	bool merge = false;
	int ptBins = 48.;
	const double ptRange[]=
    {18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84,
    97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468,
    507, 548, 592, 638, 686, 737, 790, 846, 905, 967,
    1032, 1101, 1172, 1248, 1327, 1410, 1497, 1588, 1684, 1784, 1890, 2000};

	TDirectory *curdir = gDirectory;
	TFile *f = new TFile("qcd_g.root","READ");
	assert(f && !f->IsZombie());
	TTree *tree = (TTree*)f->Get("tree");
	unsigned int N = (unsigned int)tree->GetEntries(); 
	cout<<N<<" jets in total."<<endl;
	
	float pT, weight;
	unsigned char flavor;
	
	tree->SetBranchAddress("jet_pT",&pT);
	tree->SetBranchAddress("jet_flavor",&flavor);
	tree->SetBranchAddress("jet_weight",&weight);
	
	TProfile gluonFrac("g","g",ptBins,ptRange);
  	TProfile lightquarkFrac("lq","lq",ptBins,ptRange);
  	TProfile strangeFrac("s","s",ptBins,ptRange);
  	TProfile charmFrac("c","c",ptBins,ptRange);
  	TProfile bottomFrac("b","b",ptBins,ptRange);
  	TProfile unmatchedFrac("unmatched","unmatched",ptBins,ptRange);

	for(unsigned int x=0; x != N; ++x)
	{
		tree->GetEntry(x);
		
		gluonFrac.Fill(pT, (flavor == 21)? 1:0, weight);
    	lightquarkFrac.Fill(pT, (flavor == 1 || flavor == 2)? 1:0, weight);
    	strangeFrac.Fill(pT, (flavor == 3)? 1:0, weight);
    	charmFrac.Fill(pT, (flavor == 4)? 1:0, weight);
    	bottomFrac.Fill(pT, (flavor == 5)? 1:0, weight);
    	unmatchedFrac.Fill(pT, (flavor == 0)? 1:0, weight);
    	//if(weight != 0) cout<< weight << endl;
	}

	TH1D *light_quarks = lightquarkFrac.ProjectionX("light quarks","");
  	TH1D *gluons = gluonFrac.ProjectionX("gluons","");
  	TH1D *strange = strangeFrac.ProjectionX("strange","");
  	if(merge) light_quarks->Add(strange);
  	TH1D *charm = charmFrac.ProjectionX("charm","");
  	TH1D *bottom = bottomFrac.ProjectionX("bottom","");
  	TH1D *unmatched = unmatchedFrac.ProjectionX("unmatched","");
	
  	TH1D *h = new TH1D("h",";p_{T} (GeV);Fraction",1000,20,2000);

	tdrDraw(unmatched,"",kOpenCircle,kGray+2,kSolid,-1,1001,kGray);
	gStyle->SetOptStat(kFALSE); //removes old legend
	tdrDraw(gluons,"",kPlus,kBlue+2,kSolid,-1,1001,kBlue-9);
	tdrDraw(light_quarks,"",kFullCircle,kGreen-1,kSolid,-1,1001,kYellow-9);
	tdrDraw(strange,"",kFullCircle,kAzure-6,kSolid,-1,1001,kAzure-8);
	tdrDraw(charm,"",kFullCircle,kGreen-1,kSolid,-1,1001,kGreen-9);
	tdrDraw(bottom,"",kFullCircle,kRed-2,kSolid,-1,1001,kRed-9);

	THStack *hs  = new THStack("hs","test stacked histograms");

	TCanvas *c1 = tdrCanvas("c1",h,2,0,kSquare);

	//light_quarks->Add(strange);
	hs->Add(bottom);
	hs->Add(charm);
	if(!merge)hs->Add(strange);
	hs->Add(light_quarks);
	hs->Add(gluons);
	hs->Add(unmatched);

	hs->Draw();

	hs->GetXaxis()->SetNoExponent();
	// hs->GetXaxis()->SetNoExponent(kTRUE);
	hs->GetXaxis()->SetRange(9,47);
	hs->GetXaxis()->SetMoreLogLabels(kTRUE);
	hs->GetXaxis()->SetTitle("p_{T} (GeV)");
	hs->GetYaxis()->SetTitle("Flavor fraction");
	hs->GetYaxis()->SetTitleSize(0.05);
	hs->GetYaxis()->SetTitleOffset(1.2);
	hs->GetXaxis()->SetTitleSize(0.05);
	hs->GetXaxis()->SetTitleOffset(1);
	hs->SetMaximum(0.95);
	// hs->SetLogx();

	double x0, y0;
	x0 = 0.4;
	y0 = 0.05;

	// //hadronic def
	// TLegend *leg = tdrLeg(0.642617,0.304878,0.968121,0.623693);				
	// TLegend *sample = tdrLeg(0.173658,0.3223-0.01,0.729027,0.54878-0.01);			//
	// TLegend *alphacut = tdrLeg(0.162752+0.16,0.336237,0.708054+0.16,0.409408);			//hadronic
	// TLegend *etacut = tdrLeg(0.166107,0.334495 ,0.709732,0.407666);			

	//QCDaware def
	TLegend *leg = tdrLeg(0.5+0.5,0.82-0.2,0.175+0.5,0.50-0.2);			
	TLegend *sample = tdrLeg(0.675-0.05-x0,0.50-y0,0.775-0.05-x0,0.505-y0);				//QCDaware
	TLegend *alphacut = tdrLeg(0.77-x0,0.50-0.05-y0,0.87-x0,0.505-0.05-y0);				//goes
	TLegend *etacut = tdrLeg(0.61-x0,0.50-0.05-y0,0.71-x0,0.505-0.05-y0);				//here

	////physics def	
	// TLegend *leg = tdrLeg(0.5,0.82-0.1,0.175,0.50-0.1); 				
	// TLegend *sample = tdrLeg(0.675-0.05,0.50,0.775-0.05,0.505);			//
	// TLegend *alphacut = tdrLeg(0.77,0.50-0.05,0.87,0.505-0.05);			//physics
	// TLegend *etacut = tdrLeg(0.61,0.50-0.05,0.71,0.505-0.05);			//



	sample->SetHeader("#gamma+jet sample");
	//TLegend *heading = tdrLeg(0.675-0.4,0.50+0.5,0.775-0.4,0.505+0.5); 	
	//heading->SetHeader("Hadronic Definition, #sqrt{s} = 8 TeV");
	alphacut->SetHeader("#alpha<0.3");
	etacut->SetHeader("#left|#eta#right|< 1.3,");

	leg->AddEntry(unmatched,"None","f");	
	leg->AddEntry(gluons,"Gluon","f");
	leg->AddEntry(light_quarks,"Light","f");
	leg->AddEntry(strange,"Strange","f");
	leg->AddEntry(charm,"Charm","f");
	leg->AddEntry(bottom,"Bottom","f");
	
	leg->Draw();
	// heading->Draw();
	sample->Draw();
	alphacut->Draw();
	etacut->Draw();
	
	
	gPad->SetLogx();       
	c1->SaveAs("/home/aman/Documents/Thesis/Figures/stack_final/qcd_g.pdf");
}

