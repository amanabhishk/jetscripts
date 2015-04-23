#include "tdrstyle_mod14.C"
#include "TFile.h"
#include "THStack.h"
#include <iostream>
#include "TTree.h"
#include "TProfile.h"
#include <string>
#include <sstream>

using namespace std;

void draw_Zgamma_fd()
{

	bool merge = false;
	int ptBins = 48.;
	const double ptRange[]=
    {18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84,
    97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468,
    507, 548, 592, 638, 686, 737, 790, 846, 905, 967,
    1032, 1101, 1172, 1248, 1327, 1410, 1497, 1588, 1684, 1784, 1890, 2000};

	TDirectory *curdir = gDirectory;
	TFile *f = new TFile("1000000Zjet_massive_feynman.root","READ");
	assert(f && !f->IsZombie());
	TTree *tree = (TTree*)f->Get("Events");
	unsigned int N = (unsigned int)tree->GetEntries(); 
	cout<<N<<" events."<<endl;
	
	int in,out;               //exact number of particles stored for each event 
  	float pt, wt;
	
	tree->SetBranchAddress("pT",&pt);
	tree->SetBranchAddress("weight",&wt);
	tree->SetBranchAddress("outType",&out);
	tree->SetBranchAddress("inType",&in);

	vector<TProfile*> diag(0);
	
	int i = 0;
	for(unsigned int x = 11; x < 70; x++)
	{
		if(x%10 > 6)
		{
			std::stringstream temp("");
			temp<<x;
			diag.push_back(new TProfile(temp.str().c_str(),temp.str().c_str(),ptBins,ptRange));
			//cout<<"TH1D* d"<< x << "= diag["<<i<<"]->ProjectionX(d"<<x<<","");"<<endl;
			i++;
		}
		

	}
	// 	cout<<"check1";
	for(unsigned int x=0; x != N; ++x)
	{
		tree->GetEntry(x);
		for(unsigned int k = 0; k < diag.size(); k++)
		{
			diag[k]->Fill(pt,(abs(3*in - 10 + out) == k)? 1:0,wt) ;
		}
	}
	
	
	
	TH1D* d17= diag[0]->ProjectionX("d17","");
	TH1D* d18= diag[1]->ProjectionX("d18","");
	TH1D* d19= diag[2]->ProjectionX("d19","");

	TH1D* d27= diag[3]->ProjectionX("d27","");
	TH1D* d28= diag[4]->ProjectionX("d28","");
	TH1D* d29= diag[5]->ProjectionX("d29","");
	
	TH1D* d37= diag[6]->ProjectionX("d37","");
	TH1D* d38= diag[7]->ProjectionX("d38","");
	TH1D* d39= diag[8]->ProjectionX("d39","");
	
	TH1D* d47= diag[9]->ProjectionX("d47","");
	TH1D* d48= diag[10]->ProjectionX("d48","");
	TH1D* d49= diag[11]->ProjectionX("d49","");
	
	TH1D* d57= diag[12]->ProjectionX("d57","");
	TH1D* d58= diag[13]->ProjectionX("d58","");
	TH1D* d59= diag[14]->ProjectionX("d59","");
	
	TH1D* d67= diag[15]->ProjectionX("d67","");
	TH1D* d68= diag[16]->ProjectionX("d68","");
	TH1D* d69= diag[17]->ProjectionX("d69","");	



	vector <TH1D*> all(0);
	all.push_back(d17);
	all.push_back(d18);
	all.push_back(d19);
	all.push_back(d27);
	all.push_back(d28);
	all.push_back(d29);
	all.push_back(d37);
	all.push_back(d38);
	all.push_back(d39);
	all.push_back(d47);
	all.push_back(d48);
	all.push_back(d49);
	all.push_back(d57);
	all.push_back(d58);
	all.push_back(d59);
	all.push_back(d67);
	all.push_back(d68);
	all.push_back(d69);


 	TFile* outFile = new TFile("out.root", "RECREATE");
 	TH1D *h = new TH1D("h",";p_{T} (GeV);Fraction",1000,20,2000);
 	gStyle->SetOptStat(kFALSE);
		
	// // gStyle->SetOptStat(kFALSE); //removes old legend
	// // tdrDraw(gluons,"",kPlus,kBlue+2,kSolid,-1,1001,kBlue-9);
	// // tdrDraw(light_quarks,"",kFullCircle,kGreen-1,kSolid,-1,1001,kYellow-9);
	// // tdrDraw(strange,"",kFullCircle,kAzure-6,kSolid,-1,1001,kAzure-8);
	// // tdrDraw(charm,"",kFullCircle,kGreen-1,kSolid,-1,1001,kGreen-9);
	// // tdrDraw(bottom,"",kFullCircle,kRed-2,kSolid,-1,1001,kRed-9);

	curdir->cd();
	THStack *hs  = new THStack("hs","Feynman diagrams for Z+jet event");
	TCanvas *c1 = tdrCanvas("c1",h,2,0,kSquare);

	for(unsigned int a = 0; a < all.size(); ++a) 
	{
		if(all[a]->GetMaximum() != 0) 
		{
			cout<<a<<endl;
		}
	}

	tdrDraw(d28,"",kPlus,kBlue+2,kSolid,-1,1001,kYellow-10);
	hs->Add(d28);
	tdrDraw(d39,"",kPlus,kBlue+2,kSolid,-1,1001,kYellow-9);
	hs->Add(d39);
	tdrDraw(d47,"",kPlus,kBlue+2,kSolid,-1,1001,kBlue-7);
	hs->Add(d47);
	
	hs->Draw();

	hs->GetXaxis()->SetNoExponent();
	hs->GetXaxis()->SetRange(9,47);
	hs->GetXaxis()->SetMoreLogLabels(kTRUE);
	hs->GetXaxis()->SetTitle("p_{T} (GeV)");
	hs->GetYaxis()->SetTitle("Diagram contribution");
	hs->SetMaximum(0.95);
	hs->GetYaxis()->SetTitleSize(0.05);
	hs->GetYaxis()->SetTitleOffset(1.2);
	hs->GetXaxis()->SetTitleSize(0.05);
	hs->GetXaxis()->SetTitleOffset(1);

	TLegend *sample = tdrLeg(0.675-0.05,0.50,0.775-0.05,0.505); 
	sample->SetHeader("#gamma+jet sample");
	sample->Draw();
	TLegend *leg = tdrLeg(0.5,0.2,0.275,0.40);		
	
	// leg->AddEntry(d47,"q#bar{q}#rightarrowZg","f");
	// leg->AddEntry(d39,"g#bar{q}#rightarrowZ#bar{q}","f");
	// leg->AddEntry(d28,"gq#rightarrowZq","f");

	leg->AddEntry(d47,"q#bar{q}#rightarrow#gammag","f");
	leg->AddEntry(d39,"g#bar{q}#rightarrow#gamma#bar{q}","f");
	leg->AddEntry(d28,"gq#rightarrow#gammaq","f");

	leg->Draw();
	
	gPad->SetLogx();    
	outFile->cd();

	curdir->cd();
}

