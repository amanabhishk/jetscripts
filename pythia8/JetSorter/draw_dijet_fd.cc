#include "tdrstyle_mod14.C"
#include "TFile.h"
#include "THStack.h"
#include <iostream>
#include "TTree.h"
#include "TProfile.h"
#include <string>
#include <sstream>

using namespace std;

void draw_dijet_fd()
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
	
	int in,out;               //exact number of particles stored for each event 
  	float pt, wt;
	
	tree->SetBranchAddress("pT",&pt);
	tree->SetBranchAddress("weight",&wt);
	tree->SetBranchAddress("outType",&out);
	tree->SetBranchAddress("inType",&in);

	vector<TProfile*> diag(0);
	
	for(unsigned int x = 11; x < 67; x++)
	{
		if(x%10 < 7 && x%10 != 0)
		{
			std::stringstream temp("");
			temp<<x;
			diag.push_back(new TProfile(temp.str().c_str(),temp.str().c_str(),ptBins,ptRange));
		}
	}
	// 	cout<<"check1";
	for(unsigned int x=0; x != N; ++x)
	{
		tree->GetEntry(x);
		for(unsigned int k = 0; k < diag.size(); k++)
		{
			diag[k]->Fill(pt,(abs(6*in - 7 + out) == k)? 1:0,wt) ;
		}
	}
	
	
	
	TH1D* d11 = diag[0]->ProjectionX("d11","");
	TH1D* d12 = diag[1]->ProjectionX("d12","");
	TH1D* d13 = diag[2]->ProjectionX("d13","");
	TH1D* d14 = diag[3]->ProjectionX("d14","");
	TH1D* d15 = diag[4]->ProjectionX("d15","");
	TH1D* d16 = diag[5]->ProjectionX("d16","");

	TH1D* d21 = diag[6]->ProjectionX("d21","");
	TH1D* d22 = diag[7]->ProjectionX("d22","");
	TH1D* d23 = diag[8]->ProjectionX("d23","");
	TH1D* d24 = diag[9]->ProjectionX("d24","");
	TH1D* d25 = diag[10]->ProjectionX("d25","");
	TH1D* d26 = diag[11]->ProjectionX("d26","");

	TH1D* d31 = diag[12]->ProjectionX("d31","");
	TH1D* d32 = diag[13]->ProjectionX("d32","");
	TH1D* d33 = diag[14]->ProjectionX("d33","");
	TH1D* d34 = diag[15]->ProjectionX("d34","");
	TH1D* d35 = diag[16]->ProjectionX("d35","");
	TH1D* d36 = diag[17]->ProjectionX("d36","");

	TH1D* d41 = diag[18]->ProjectionX("d41","");
	TH1D* d42 = diag[19]->ProjectionX("d42","");
	TH1D* d43 = diag[20]->ProjectionX("d43","");
	TH1D* d44 = diag[21]->ProjectionX("d44","");
	TH1D* d45 = diag[22]->ProjectionX("d45","");
	TH1D* d46 = diag[23]->ProjectionX("d46","");	

	TH1D* d51 = diag[24]->ProjectionX("d51","");	
	TH1D* d52 = diag[25]->ProjectionX("d52","");	
	TH1D* d53 = diag[26]->ProjectionX("d53","");	
	TH1D* d54 = diag[27]->ProjectionX("d54","");	
	TH1D* d55 = diag[28]->ProjectionX("d55","");	
	TH1D* d56 = diag[29]->ProjectionX("d56","");

	TH1D* d61 = diag[30]->ProjectionX("d61","");	
	TH1D* d62 = diag[31]->ProjectionX("d62","");	
	TH1D* d63 = diag[32]->ProjectionX("d63","");	
	TH1D* d64 = diag[33]->ProjectionX("d64","");	
	TH1D* d65 = diag[34]->ProjectionX("d65","");	
	TH1D* d66 = diag[35]->ProjectionX("d66","");	



	vector <TH1D*> all(0);

	all.push_back(d11);
	all.push_back(d12);
	all.push_back(d13);
	all.push_back(d14);
	all.push_back(d15);
	all.push_back(d16);

	all.push_back(d21);
	all.push_back(d22);
	all.push_back(d23);
	all.push_back(d24);
	all.push_back(d25);
	all.push_back(d26);

	all.push_back(d31);
	all.push_back(d32);
	all.push_back(d33);
	all.push_back(d34);
	all.push_back(d35);
	all.push_back(d36);

	all.push_back(d41);
	all.push_back(d42);
	all.push_back(d43);
	all.push_back(d44);
	all.push_back(d45);
	all.push_back(d46);

	all.push_back(d51);
	all.push_back(d52);
	all.push_back(d53);
	all.push_back(d54);
	all.push_back(d55);
	all.push_back(d56);

	all.push_back(d61);
	all.push_back(d62);
	all.push_back(d63);
	all.push_back(d64);
	all.push_back(d65);
	all.push_back(d66);

	// // TH1D *light_quarks = lightquarkFrac.ProjectionX("light quarks","");
 // //  	TH1D *gluons = gluonFrac.ProjectionX("gluons","");
 // //  	TH1D *strange = strangeFrac.ProjectionX("strange","");
 // //  	if(merge) light_quarks->Add(strange);
 // //  	TH1D *charm = charmFrac.ProjectionX("charm","");
 // //  	TH1D *bottom = bottomFrac.ProjectionX("bottom","");
 // //  	TH1D *unmatched = unmatchedFrac.ProjectionX("unmatched","");
	
 	assert(all.size()==diag.size());
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
	THStack *hs  = new THStack("hs","stack");

	// TH1D *h = new TH1D("h",";p_{T} (GeV);Fraction",1000,20,2000);
	TCanvas *c1 = tdrCanvas("c1",h,2,0,kSquare);

	int colors[]={kBlue-9,kYellow-9,kAzure-8,kGreen-9,kRed-9,kYellow-5,kCyan-6};
	for(unsigned int a = 0; a < all.size(); ++a) 
	{
		if(all[a]->GetMaximum() != 0) 
		{
			cout<<a<<endl;
		}
	}

	tdrDraw(d55,"",kPlus,kBlue+2,kSolid,-1,1001,kYellow-10);
	hs->Add(d55);
	tdrDraw(d44,"",kPlus,kBlue+2,kSolid,-1,1001,kYellow-9);
	hs->Add(d44);
	tdrDraw(d14,"",kPlus,kBlue+2,kSolid,-1,1001,kYellow-6);
	hs->Add(d14);
	tdrDraw(d66,"",kPlus,kBlue+2,kSolid,-1,1001,kYellow-2);
	hs->Add(d66);
	
	tdrDraw(d33,"",kPlus,kBlue+2,kSolid,-1,1001,kRed-10);
	hs->Add(d33);
	tdrDraw(d22,"",kPlus,kBlue+2,kSolid,-1,1001,kRed-9);
	hs->Add(d22);
	
	tdrDraw(d41,"",kPlus,kBlue+2,kSolid,-1,1001,kBlue-7);
	hs->Add(d41);
	tdrDraw(d11,"",kPlus,kBlue+2,kSolid,-1,1001,kBlue-9);
	hs->Add(d11);


	// // //light_quarks->Add(strange);
	// // hs->Add(bottom);
	// // hs->Add(charm);
	// // if(!merge)hs->Add(strange);
	// // hs->Add(light_quarks);
	// // hs->Add(gluons);
	// // hs->Add(unmatched);

	
	hs->Draw();

	hs->GetXaxis()->SetNoExponent();
	// hs->GetXaxis()->SetNoExponent(kTRUE);
	hs->GetXaxis()->SetRange(9,47);
	hs->GetXaxis()->SetMoreLogLabels(kTRUE);
	//hs->GetXaxis()->SetNdivisions(5,kTRUE);
	hs->GetXaxis()->SetTitle("p_{T} (GeV)");
	hs->GetYaxis()->SetTitle("Diagram contribution");
	hs->GetYaxis()->SetTitleSize(0.05);
	hs->GetYaxis()->SetTitleOffset(1.2);
	hs->GetXaxis()->SetTitleSize(0.05);
	hs->GetXaxis()->SetTitleOffset(1);

	TH1D* mid = new TH1D("line","line",ptBins,ptRange);
	mid->Add(d55);
	mid->Add(d14);
	mid->Add(d66);
	mid->Add(d44);
	mid->Add(d22,0.5);
	mid->Add(d33,0.5);
	tdrDraw(mid,"HIST",kNone,kBlack,kDashed,-1,kNone,0);



	
	//mid->Draw();
	// hs->SetLogx();

	// // double x0, y0;
	// // x0 = 0.45;
	// // y0 = 0.05;
	// // //TLegend *leg = tdrLeg(0.5,0.82-0.1,0.175,0.50-0.1); 				//physics def
	// // //TLegend *leg = tdrLeg(0.5,0.82+0.07,0.175,0.50+0.07);				//hadronic def
	TLegend *leg = tdrLeg(0.199664,0.437282,0.426174,0.87108,0.03);		//QCDaware def
	

	
	// // TLegend *sample = tdrLeg(0.675-0.05-x0,0.50-y0,0.775-0.05-x0,0.505-y0);				//QCDaware
	// // TLegend *alphacut = tdrLeg(0.77-x0,0.50-0.05-y0,0.87-x0,0.505-0.05-y0);				//goes
	// // TLegend *etacut = tdrLeg(0.61-x0,0.50-0.05-y0,0.71-x0,0.505-0.05-y0);				//here

	// // TLegend *heading = tdrLeg(0.675-0.4,0.50+0.44,0.775-0.4,0.505+0.44); 	
	// // // TLegend *sample = tdrLeg(0.675-0.05,0.50,0.775-0.05,0.505);			//everything except
	// // // TLegend *alphacut = tdrLeg(0.77,0.50-0.05,0.87,0.505-0.05);			//QCD
	// // // TLegend *etacut = tdrLeg(0.61,0.50-0.05,0.71,0.505-0.05);			//aware

	// // sample->SetHeader("Z+jet sample");
	// // heading->SetHeader("Hadronic Definition, #sqrt{s} = 8 TeV");
	// // alphacut->SetHeader("#alpha<0.3");
	// // etacut->SetHeader("#left|#eta#right|< 1.3,");

	
	leg->AddEntry(d11,"gg#rightarrowgg","f");
	leg->AddEntry(d41,"q#bar{q}#rightarrowgg","f");
	leg->AddEntry(d22,"gq#rightarrowgq","f");
	leg->AddEntry(d33,"g#bar{q}#rightarrowg#bar{q}","f");
	leg->AddEntry(d66,"#bar{q}#bar{q}#rightarrow#bar{q}#bar{q}","f");
	leg->AddEntry(d14,"gg#rightarrowq#bar{q}","f");
	leg->AddEntry(d44,"q#bar{q}#rightarrowq#bar{q}","f");
	leg->AddEntry(d55,"qq#rightarrowqq","f");
	// // leg->AddEntry(charm,"Charm","f");
	// // leg->AddEntry(strange,"Strange","f");
	// // leg->AddEntry(light_quarks,"Light","f");
	// // leg->AddEntry(gluons,"Gluon","f");
	// // leg->AddEntry(unmatched,"None","f");
	
	leg->Draw();
	// // heading->Draw();
	// // sample->Draw();
	// // alphacut->Draw();
	// // etacut->Draw();
	
	
	
	gPad->SetLogx();    
	outFile->cd();
	//hs->Write(); 

	//outFile->Close();
	curdir->cd();

	// // //c1->SaveAs("output.pdf");
}

