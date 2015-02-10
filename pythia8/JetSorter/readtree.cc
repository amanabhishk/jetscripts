{
	TFile *f = new TFile("output.root");
	TTree *Events = (TTree*)f->Get("Events");

	Int_t event_count = (Int_t)Events->GetEntries(); 
	cout<<event_count<<endl;
	Float_t eta[5000]; //with nmax greater or equal to the max value for nGenjet5.
	
	Events->SetBranchAddress("eta",eta);
	UShort_t size;
	Events->SetBranchAddress("n",&size);
	TH1D *eta_hist = new TH1D("eta","eta",100,-60,60);
	for (Int_t i=0; i != event_count; ++i) 
	{
	// 	UShort_t size;
 //  		
  		Events->GetEntry(i);
  		// cout<<size<<endl;
 //  		// UShort_t id[size];
	// 	// Float_t eta[size], eta[size], eta[size], m[size];
  		// eta_hist->Fill(size);
  		for(UShort_t k = 0; k != size; ++k)
  		{
			eta_hist->Fill(eta[k]);
	// // 		// cout<<k+1<<" "<<id[k]<<" "<<eta[k]<<" "<<eta[k]<<" "<<m[k]<<endl;


  		}
  		
  		// cout<<"*************************************\n";

	}
	eta_hist->Draw();
	
}

