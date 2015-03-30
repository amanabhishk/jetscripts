// This holds all the auxiliary functions and classes of jetsorter

// Stdlib header file for input and output.
#include <iostream>
#include <stdlib.h>
#include <vector>
#include <string>
#include <sstream>

// Header file to access Pythia 8 program elements.
#include "Pythia8/Pythia.h"
#include "TF1.h"

#include <math.h>
// FastJet interface
#include "Pythia8/FastJet3.h"

// ROOT, for histogramming.
#include "TROOT.h"
#include "TH1.h"
#include "TProfile.h"
#include "THStack.h"
#include "TStyle.h"

// ROOT, for interactive grdPhics.
#include "TVirtualPad.h"
#include "TApplication.h"

// ROOT, for saving file.
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"

#include "TMatrix.h"
#include "TMatrixDSymEigen.h"

#include "constants.h"
#include "TTree.h"
#include "TLorentzVector.h"

using namespace Pythia8;




double deltaPhi(double phi1, double phi2){
  double pi = 3.141592;
  
  while ( phi1 < 0 ) phi1 += 2.*pi;
  while ( phi1 > 2*pi ) phi1 -= 2.*pi;
  while ( phi2 < 0 ) phi2 += 2.*pi;
  while ( phi2 > 2*pi ) phi2 -= 2.*pi;

  double dPhi = abs(phi1 - phi2);
  if(dPhi>pi) dPhi = 2*pi - dPhi;
  return dPhi;
}

double deltaEta(double eta1, double eta2){
  return eta1-eta2;
}

double deltaR( double phi1, double phi2, double eta1, double eta2 ){
  
  double dEta = deltaEta(eta1,eta2);
  
  double dPhi = deltaPhi(phi1,phi2);

  double dR = pow( pow( dPhi, 2 ) + pow( dEta, 2 ) , 0.5 );

  // cout<<"dEta:"<<dEta<<" dPhi:"<<dPhi<<" dR:"<<dR<<endl;
  return dR;
}

bool isCharged(const int& c)
{
  int pdgid = abs(c), digit, charge = 0;
   //photon and neutrinos
  if(pdgid == 22 || pdgid == 12 ||pdgid == 14 ||pdgid == 16 ) return false;
  //charged leptons
  if(pdgid == 11 ||pdgid == 13 ||pdgid == 15 ) return true; 

  pdgid = (pdgid/10)%1000;
  if(pdgid < 100) //Meson
  {
    if((pdgid%10)%2 == 0) charge += 2;
    else charge += -1;
    
    if((pdgid/10)%2 == 0) charge += -2;
    else charge += 1;
    
    if(charge == 0) return false;
    else return true;
  }
  else //Baryon
  {
    while(pdgid != 0)
    {
      digit = pdgid%10;
      pdgid = pdgid/10;
      if(digit%2 == 0) charge += 2;
      else charge += -1;  
    }
    if(charge == 0) return false;
    else return true; 
  } 
}

bool isHadron(const int& c)
{
  if(abs(c)>99) return true;
  else return false;
}

bool isMeson(const int& c)
{
  if(!isHadron(c))
  {
    cout<<"A non-hadron was tested. Abort.\n";
    exit(0);
  }
  int pdgid = abs(c);
  pdgid = (pdgid/10)%1000;
  if(pdgid<100) return true;
  else return false;
}

double pTD(const fastjet::PseudoJet& jet){
  
  double num = 0, den = 0;
  vector<fastjet::PseudoJet> jetParts = jet.constituents();

  for(unsigned int q = 0; q != jetParts.size(); ++q) 
  {
    num += pow(jetParts[q].pt(),2);
    den += jetParts[q].pt();
  }
  //for(unsigned int q = 0; q != jetParts.size(); ++q) 

  num = pow(num, 0.5);
  return num/den;
}

void implicit_cuts(const vector<fastjet::PseudoJet>& jet_ref, vector<fastjet::PseudoJet>& jet)
{
  int id;
  for(unsigned int q = 0; q != jet_ref.size(); ++q) 
  {
    id = jet_ref[q].user_index();
    if(isHadron(id))
    {
      if(isCharged(id) && jet_ref[q].pt() > 0.3) jet.push_back(jet_ref[q]);
      else if(!isCharged(id) && jet_ref[q].pt() > 3) jet.push_back(jet_ref[q]);
      else continue;
    }
    else jet.push_back(jet_ref[q]);
  }  
}

void explicit_cuts(const vector<fastjet::PseudoJet>& jet_ref, vector<fastjet::PseudoJet>& jet)
{
  int id;
  for(unsigned int q = 0; q != jet_ref.size(); ++q) 
  {
    id = jet_ref[q].user_index();
    if(id == 22 && jet_ref[q].pt() < 1) continue;
    else if(isHadron(id) && !isCharged(id) && jet_ref[q].pt() < 1) continue;
    else jet.push_back(jet_ref[q]);
  }  
}

unsigned int multiplicity(const fastjet::PseudoJet& jet, const unsigned char& cut)
{
  if(cut == 0) return jet.constituents().size();
  else if(cut == 1)
  {
    vector <fastjet::PseudoJet> jetParts(0);
    explicit_cuts(jet.constituents(),jetParts);
    return jetParts.size(); 
  }
  else if(cut == 2)
  {
    vector <fastjet::PseudoJet> jetParts(0);
    vector <fastjet::PseudoJet> jetParts2(0);
    explicit_cuts(jet.constituents(), jetParts);
    implicit_cuts(jetParts,jetParts2);
    return jetParts2.size();
  }
  
}

void sigma2(const fastjet::PseudoJet& jet, float* output){

  vector <fastjet::PseudoJet> jetParts = jet.constituents(); 
  double M11 = 0, M22 = 0, M12 = 0, M11_cut = 0, M22_cut = 0, M12_cut = 0;
  double phi(0), eta(0), pT2(0);
  
  for(unsigned int q = 0; q != jetParts.size(); ++q) 
  {
    pT2 += pow(jetParts[q].pt(),2);
    eta += pow(jetParts[q].pt(),2)*jetParts[q].eta();
    phi += pow(jetParts[q].pt(),2)*jetParts[q].phi();
  }

  eta = eta/pT2;
  phi = phi/pT2;

  int id;

  for(unsigned int q = 0; q != jetParts.size(); ++q) 
  {
    id = jetParts[q].user_index();
    if((jetParts[q].pt()>1 && !isCharged(id) && isHadron(id)) || !isHadron(id) || (isHadron(id) && isCharged(id)))
    {
      M11_cut += pow(jetParts[q].pt()*deltaEta(jetParts[q].eta(),eta),2);
      M22_cut += pow(jetParts[q].pt()*deltaPhi(jetParts[q].phi(),phi),2);
      M12_cut += -pow(jetParts[q].pt(),2)*deltaEta(jetParts[q].eta(),eta)*deltaPhi(jetParts[q].phi(),phi);          
    }
    //else cout<<id<<endl;

    M11 += pow(jetParts[q].pt()*deltaEta(jetParts[q].eta(),eta),2);
    M22 += pow(jetParts[q].pt()*deltaPhi(jetParts[q].phi(),phi),2);
    M12 += -pow(jetParts[q].pt(),2)*deltaEta(jetParts[q].eta(),eta)*deltaPhi(jetParts[q].phi(),phi);    
  }
  
  double e[4] = {
    M11, M12,
    M12, M22
  };

  double e_cut[4] = {
    M11_cut, M12_cut,
    M12_cut, M22_cut
  };

  TMatrixDSym m(2, e);
  TMatrixDSymEigen me(m);

  TMatrixDSym m_cut(2, e_cut);
  TMatrixDSymEigen me_cut(m_cut);


  TVectorD eigenval = me.GetEigenValues();
  TVectorD eigenval_cut = me_cut.GetEigenValues();

  output[0] = pow(eigenval[1]/pT2,0.5);
  output[1] = pow(eigenval_cut[1]/pT2,0.5);
  //return pow(eigenval[1]/pT2,0.5);
}

void fillTree(const Particle& p, const int& count, Int_t* id, Float_t* pT, Float_t* eta, Float_t* phi, Float_t* m){
  id[count] = p.id();
  pT[count] = p.pT();
  eta[count] = p.eta();
  phi[count] = p.phi();
  m[count] = p.m();
}

bool isStrange( const int& c ) {
  int id = abs(c);
  int digit;
  
  if(id < 100)
  {
    if(id==3) return true;
    else return false;
  }

  id = id/10;
  bool result = false;
  while(id != 0)
  {
    digit = id%10;
    id = id/10;
    if(id==3) result = true;
    if(id>3) return false;  
  }
  return result;
}

bool isCharm( const int& c ) {
  int id = abs(c);
  int digit;
  
  if(id < 100)
  {
    if(id==4) return true;
    else return false;
  }

  id = id/10;
  bool result = false;
  while(id != 0)
  {
    digit = id%10;
    id = id/10;
    if(id==4) result = true;
    if(id>4) return false;  
  }
  return result;
}

bool isBottom( const int& c ) {
  int id = abs(c);
  int digit;
  
  if(id < 100)
  {
    if(id==5) return true;
    else return false;
  }

  id = id/10;
  bool result = false;
  while(id != 0)
  {
    digit = id%10;
    id = id/10;
    if(id==5) result = true;
    if(id>5) return false;  
  }
  return result;
}


void hadronic_definition_status_codes( Event& event, const int& index, int& count, UChar_t* status, Int_t* id, Float_t* pT, Float_t* eta, Float_t* phi, Float_t* m){
  if((abs(event[index].status())==71 || abs(event[index].status())==72) && abs(event[index].id()) < 100)
  {
    ++count;
    status[count] = 70;
    fillTree(event.at(index), count, id, pT, eta, phi, m);
  }

  if(isBottom(event[index].id()) && event[index].id() > 99)
  {
    if(!isBottom(event[event[index].daughter1()].id()) && !isBottom(event[event[index].daughter2()].id()))
    {
      ++count;
      status[count] = 5;
      fillTree(event.at(index), count, id, pT, eta, phi, m);
    }
  }

  if(isCharm(event[index].id()) && event[index].id() > 99)
  {
    if(!isCharm(event[event[index].daughter1()].id()) && !isCharm(event[event[index].daughter2()].id()))
    {
      ++count;
      status[count] = 4;
      fillTree(event.at(index), count, id, pT, eta, phi, m);
    }
  }
  
  if(isStrange(event[index].id()) && event[index].id() > 99)
  {
    if(!isStrange(event[event[index].daughter1()].id()) && !isStrange(event[event[index].daughter2()].id()))
    {
      ++count;
      status[count] = 6;
      fillTree(event.at(index), count, id, pT, eta, phi, m);
    }
  } 
}

void printTime(int total){
  total = total/CLOCKS_PER_SEC;

  int h = total/3600;
  int min = (total%3600)/60;
  int sec = total - 3600*h - 60*min;

  cout<<"Done in "<<h<<" hours "<<min<<" minutes "<<sec<<" seconds."<<endl;
}

void printTime(int total, const int& events){
  total = total/CLOCKS_PER_SEC;

  int h = total/3600;
  int min = (total%3600)/60;
  int sec = total - 3600*h - 60*min;

  cout<<"Done in "<<h<<" hours "<<min<<" minutes "<<sec<<" seconds."<<endl;
  cout<<total*100/double(events)<<"seconds/1000 events.\n";
}

void qcd_aware_clustering_flavor(const fastjet::PseudoJet& jet, TH1D*& h){
  vector <fastjet::PseudoJet> componenets;
  int temp = 0;
  componenets = jet.constituents();
  for(unsigned int m = 0; m != componenets.size(); ++m)
  {
    if(componenets[m].user_index()==0) continue;
    else
    {
      temp++;
    }
  }
  h->Fill(temp);
}

unsigned int hadron_count(const vector<fastjet::PseudoJet>& jets){
  unsigned int output = 0, check = 0;
  
  for(unsigned int x = 0; x != jets.size(); ++x)
  {
    if(jets[x].pt()<2000) break;
    if(jets[x].pt()>2500) continue;
    ++check;
    vector<fastjet::PseudoJet> parts = jets[x].constituents();

    for(unsigned int y = 0; y != parts.size(); ++y)
    {
      if(isHadron(parts[y].user_index()) && parts[y].pt()>100) ++output;
    }
  }
  //cout<<output<<endl;
  return output;
}

struct clustered_info
{
    vector<int> leptonList;
    int gamma;
    vector <fastjet::PseudoJet> sortedJets;
    vector <int> jetFlavor;
    float weight;
};


bool is_good_event(const float* pT, const float* eta, const float* phi, const float* m, const struct clustered_info& data, const int& sample)
{
  if(sample == 1)//selecting good dijet events
  {
    bool dijetCriteria = false;
    if(data.sortedJets.size()<2) dijetCriteria = false;
    else if(data.sortedJets.size()>2)
    {
     dijetCriteria = deltaPhi(data.sortedJets[0].phi(),data.sortedJets[1].phi())>2.8 && 0.15*fabs(data.sortedJets[0].pt()+data.sortedJets[1].pt())>data.sortedJets[2].pt();
    }
    else
    {
     dijetCriteria = deltaPhi(data.sortedJets[0].phi(),data.sortedJets[1].phi())>2.8;
    }
    return dijetCriteria;
  }

  else if(sample == 2)//selecting good Z-jet events
  {
    // Checking sufficient resolution
    TLorentzVector v1, v2; 
    if(deltaR(phi[data.leptonList[0]],data.sortedJets[0].phi(),eta[data.leptonList[0]],data.sortedJets[0].eta()) < R) return false;
    if(deltaR(phi[data.leptonList[1]],data.sortedJets[0].phi(),eta[data.leptonList[1]],data.sortedJets[0].eta()) < R) return false;

    //the pT of the muons are required to be greater than 20 and 10 GeV, respectively
    if(!((pT[data.leptonList[0]]>20 && pT[data.leptonList[1]]>10) || (pT[data.leptonList[1]]>20 && pT[data.leptonList[0]]>10))) return false;

    //the subleading jet in the event is required to have a pT smaller than 30% of that of the dimuon system.
    v1.SetPtEtaPhiM(pT[data.leptonList[0]],eta[data.leptonList[0]],phi[data.leptonList[0]],m[data.leptonList[0]]);
    v2.SetPtEtaPhiM(pT[data.leptonList[1]],eta[data.leptonList[1]],phi[data.leptonList[1]],m[data.leptonList[1]]);
    if(data.sortedJets[1].pt()>0.3*(v1+v2).Pt()) return false;

    //the dimuon invariant mass is required to fall in the 70-110 GeV range
    if(fabs((v1+v2).M())<70 || fabs((v1+v2).M())>110) return false;
    return true;
  }

  else if(sample == 3)//selecting good data.gamma-jet events
  {
    if(data.sortedJets[1].pt()>0.3*pT[data.gamma]) return false;
    if(deltaR(phi[data.gamma],data.sortedJets[0].phi(),eta[data.gamma],data.sortedJets[0].eta()) < R) return false;
    return true;
  }
  
  else return false;
}



class jet_data
{
  public:

    void fill(const struct clustered_info& data, const int& sample)
    {
      for(int k = 0; k != data.sortedJets.size(); ++k)
      {
        if(sample == 1 && k == 2) break;
        if((sample == 2 || sample == 3) && k == 1) break; 
        
        if(fabs(data.sortedJets[k].eta()) > etaMax) continue;
        
        Jw = data.weight;
        JpT = data.sortedJets[k].pt();
        Jmul[0] = multiplicity(data.sortedJets[k],0);
        Jmul[1] = multiplicity(data.sortedJets[k],2);
        Jflavor = data.jetFlavor[k];
        JpTD = pTD(data.sortedJets[k]);
        sigma2(data.sortedJets[k],Jsigma2);
        tree->Fill();
      }
    }

    jet_data()
    {
      tree->Branch("jet_weight", &Jw , "jet_weight/F" );
      tree->Branch("jet_pT", &JpT , "jet_pT/F" );
      tree->Branch("jet_pTD", &JpTD , "jet_pTD/F" );
      tree->Branch("jet_sigma2", Jsigma2 , "jet_sigma2[2]/F" );
      tree->Branch("jet_multiplicity", Jmul , "jet_multiplicity[2]/i" );
      tree->Branch("jet_flavor", &Jflavor , "jet_flavor/b" ); 
    }

    void save()
    {
      tree->AutoSave("Overwrite");
    }
  
  protected:
    float Jw, JpT, JpTD, Jsigma2[2];
    unsigned int Jmul[2];
    unsigned char Jflavor;
    TTree* tree = new TTree("tree","tree");

};