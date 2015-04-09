#include "tdrstyle_mod14.C"
#include "TFile.h"
#include "THStack.h"
#include <iostream>
#include "TTree.h"
#include "TProfile.h"
#include <string>
#include <sstream>


using namespace std;

void compare()
{
	TDirectory *curdir = gDirectory;

	//vector<vector<string>> files;
	vector<string> files  = {"had_z.root","had_g.root","had_d.root","phy_z.root","phy_g.root","phy_d.root"};
	//vector<string> qcd_def  = {"z1.root","g1.root","d1.root"};
	//file.push_back(physics_def);
	//file.push_back(qcd_def);
	vector<string> algo  = {"qcd","qcd","qcd","physics","physics","physics"};
    vector<int> markers = {kFullSquare,kFullCircle,kFullTriangleUp,kOpenSquare,kOpenCircle,kOpenTriangleUp};
    vector<int> color = {kRed-3,kBlue,kGreen-6};
    vector<string> variable = {"sigma2","pTD","multiplicity"};
    vector<string> partons = {"gluon","quark","unmatched"};
    vector<float> x_min = {0, 0, 0};
    vector<float> x_max = {0.2, 1, 60};
    vector<float> y_min = {0, 0, 0};
    vector<float> y_max = {0.08, 0.08, 0.1};
    vector<int> bins = {100, 100, 60};

    vector<vector<vector<TH1D*>>> plots;
    //for(unsigned int t=0; t != variable.size(); ++t) //loop over variables
	//{
	for(unsigned int q=0; q != files.size(); ++q) //loop over files
	{
		std::stringstream tmpString("");
		tmpString << files[q];
		TFile *f = new TFile(tmpString.str().c_str(),"READ");
		assert(f && !f->IsZombie());
		
		TTree *tree = (TTree*)f->Get("tree");
		assert(tree && !tree->IsZombie());
		unsigned int N = (unsigned int)tree->GetEntries(); 

		float pT, weight, pTD, sigma2;
		unsigned int constituents;
		unsigned char flavor;
		
		tree->SetBranchAddress("jet_pT",&pT);
		tree->SetBranchAddress("jet_sigma2",&sigma2);
		tree->SetBranchAddress("jet_pTD",&pTD);
		tree->SetBranchAddress("jet_flavor",&flavor);
		tree->SetBranchAddress("jet_weight",&weight);
		tree->SetBranchAddress("jet_multiplicity",&constituents);

		
		vector<vector<TH1D*>> one_sample;
		for(unsigned int p = 0; p < partons.size(); ++p) //book histograms for gluons, quarks, and unmatched
		{
			vector<TH1D*> s;
			for(unsigned int k=0; k < variable.size(); ++k) //QGL variables
			{
				std::stringstream temp("");
				temp << variable[k] << "_" << partons[p];
				s.push_back(new TH1D(temp.str().c_str(),temp.str().c_str(),bins[k],x_min[k],x_max[k]));
			}
			one_sample.push_back(s);
		}

		plots.push_back(one_sample);
		for(unsigned int x=0; x != N; ++x)
		{
			tree->GetEntry(x);

			if(pT>100 || pT<80) continue;
			if(flavor == 21)
			{
				plots[q][0][0]->Fill(sigma2,weight);
				plots[q][0][1]->Fill(pTD,weight);
				plots[q][0][2]->Fill(constituents,weight);
			}
			else if(flavor == 0) 
			{
				plots[q][2][0]->Fill(sigma2,weight);
				plots[q][2][1]->Fill(pTD,weight);
				plots[q][2][2]->Fill(constituents,weight);
			}
			else 
			{
				plots[q][1][0]->Fill(sigma2,weight);
				plots[q][1][1]->Fill(pTD,weight);
				plots[q][1][2]->Fill(constituents,weight);
			}
		}
		
		
	}
	
	setTDRStyle();
	
	// vector<vector<TH1D*>> final;
	// for(int )
	vector<TH1D*> xkcd;
	TH1D *h1 = new TH1D("h1",";#sigma_{2};Events",bins[0],x_min[0],x_max[0]);
	TH1D *h2 = new TH1D("h2",";pTD;Events",bins[1],x_min[1],x_max[1]);
	TH1D *h3 = new TH1D("h3",";multiplicity;Events",bins[2],x_min[2],x_max[2]);
	xkcd.push_back(h1);
	xkcd.push_back(h2);
	xkcd.push_back(h3);
		
	for(unsigned int i = 0; i < xkcd.size(); ++i)
	{
		xkcd[i] ->SetMinimum(y_min[i]);
		xkcd[i] ->SetMaximum(y_max[i]);
		xkcd[i] ->GetYaxis()->SetNoExponent();
		xkcd[i] ->GetXaxis()->SetNoExponent();
		xkcd[i] ->GetXaxis()->SetRangeUser(x_min[i],x_max[i]);
	}

	vector<TCanvas*> c(6);

	for(int l = 0; l<3 ; ++l)
	{
		std::stringstream c_head("");
		//c_head<<"c";

		for(int t = 0; t < 6; ++t)
		{
			//c_head = "";
			//c_head << t+1;
			if(t%3==0) 
			{
				c_head << variable[l] <<"_"<< algo[t];
				c[t] = tdrCanvas(c_head.str().c_str(),xkcd[l],0,33);
				c_head.str("");
				//c_head << "c";
				
			}
			tdrDraw(plots[t][0][l] ,"P",markers[t%3] ,color[t%3]);
			tdrDraw(plots[t][1][l],"P",markers[3+t%3],color[t%3]);
			plots[t][0][l]->Scale(1/plots[t][0][l]->Integral());
			plots[t][1][l]->Scale(1/plots[t][1][l]->Integral());

			//TLegend *leg = tdrLeg(0.5,0.82,0.175,0.50);
			
			if(t%3==2) 
			{
				TLegend* leg = tdrLeg(0.88+0.05,0.72,0.675+0.05,0.30+0.045*4);
				leg->SetHeader("gluon");
				leg->AddEntry(plots[t-2][0][l],"Z+jet","P");
				leg->AddEntry(plots[t-1][0][l],"#gamma+jet","P");
				leg->AddEntry(plots[t][0][l],"dijet","P");
				leg->Draw();

				TLegend* gel = tdrLeg(0.88+0.05,0.48,0.675+0.05,0.24);
				gel->SetHeader("quark");
				gel->AddEntry(plots[t-2][1][l],"Z+jet","P");
				gel->AddEntry(plots[t-1][1][l],"#gamma+jet","P");
				gel->AddEntry(plots[t][1][l],"dijet","P");
				gel->Draw();				
			}

			//c[t]->SaveAs("test.jpg");
	
		}
	}
}

