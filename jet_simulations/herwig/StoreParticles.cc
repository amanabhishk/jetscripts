#include "StoreParticles.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/PDT/StandardMatchers.h"

using namespace jetanalysis;
using namespace ThePEG;

#ifndef LWH_AIAnalysisFactory_H
#ifndef LWH 
#define LWH ThePEGLWH
#endif
#include "ThePEG/Analysis/LWH/AnalysisFactory.h"
#endif

#include <iostream>
//#include "../generic/help_functions.h"
//#include "herwig_functions.h"

using std::cout;
using std::endl;

int StoreParticles::getStatusCode(tPPtr part) const
{
   int status = 1;
   size_t nChildren = part->children().size();
   if ( nChildren > 0 || part->next() ) {
      tStepPtr step = part->birthStep();
      if ((!step || (step && (!step->handler() || step->handler() == eh))) && part->id() != 82)
         status = 3;
      else
         status = 2;
   }
   /* This is the official ThePEG method for obtaining status, but it does not
    * work with the CMSSW interpretation of jet flavours. */
   // if ( nChildren > 0 || part->next() ) status = 11; 
   // if ( nChildren > 1 ) { 
   //    long id = part->data().id();
   //    if ( BaryonMatcher::Check(id) || MesonMatcher::Check(id) ||
   //       id == ParticleID::muminus || id == ParticleID::muplus ||
   //       id == ParticleID::tauminus || id == ParticleID::tauplus )
   //          if ( part->mass() <= part->data().massMax() &&
   //             part->mass() >= part->data().massMin() ) status = 2;
   // }
   return status;
}
  
void StoreParticles::analyze(tEventPtr event, long ieve, int loop, int state) 
{
   /* Rotate to CMS, extract final state particles and call analyze(particles). */
   AnalysisHandler::analyze(event, ieve, loop, state);
   if ( loop > 0 || !event ) return;
   // if ( loop > 0 || status != 0 || !event ) return;

   //tPVector parts=event->getFinalState();
   tPVector parts;
   event->select(std::back_inserter(parts),SelectAll());
   eh = event->primaryCollision()->handler();
   
   int index, count = -1; 
   
   /* Loop over all particles. */ 
   for (tPVector::const_iterator pit = parts.begin(); pit != parts.end(); ++pit) {
      int absId = abs( (*pit)->id() );
      if ( absId==2101 || absId ==2203 || absId==82 ) continue; // uu, ud, p+rem
      int pStatus = getStatusCode( *pit );
      //if ( pStatus==3 ) continue; // Beam particles and partons, not of interest
      index = (*pit)->number();
      if(index==6 || index==7) pStatus = 3;
      /* Normal end state particles */
      if (pStatus == 1 || pStatus == 3) 
      {
         //pEvent->AddPrtcl((*pit)->momentum().x(),(*pit)->momentum().y(),(*pit)->momentum().z(),
         //   (*pit)->momentum().t(), (*pit)->id(), (*pit)->data().charge(), pStatus);
        count++;
        id[count] = (*pit)->id();
        pT[count] = (*pit)->momentum().t();
        eta[count] = (*pit)->momentum().x();
        phi[count] = (*pit)->momentum().y();
        m[count] = (*pit)->momentum().z();
        status[count] = pStatus;
        N = UShort_t(count+1);
      }
      
      /* Ghost particles */
      // int ghostStatus = 0;
      
      // if ( (pStatus==2) && ( (absId<=6 && absId>0) || absId==21 ) ) {
      //       ghostStatus = 11; /* Partons */
      // } else if (absId >= 100) { /* Status codes below this are not conventional hadrons */
      //    if (hasBottom(absId) && !isExcitedHadronState(*pit,5) ) {
      //       ghostStatus = 12; /* b Hadrons */
      //    }
      //    if (hasCharm(absId) && !isExcitedHadronState(*pit,4) ) {
      //       ghostStatus = 13; /* c Hadrons */
      //    }
      // } /* Add Leptons here if needed */
      
      // if (ghostStatus) {
      //    pEvent->AddPrtcl((*pit)->momentum().x(),(*pit)->momentum().y(),(*pit)->momentum().z(),
      //       (*pit)->momentum().t(), (*pit)->id(), (*pit)->data().charge(),ghostStatus);
      // }
   }
   // TODO ghost particles/partons
   
   herwigTree->Fill();
   //pEvent->Clear();
}

void StoreParticles::dofinish() 
{
  AnalysisHandler::dofinish();
  
  herwigTree->GetCurrentFile();
  herwigTree->AutoSave("Overwrite");
  herwigFile->Close();
  cout << "StoreParticles: root tree has been written to a file" << endl;  
}

void StoreParticles::doinitrun() 
{
  cout<<"Check3.\n";
  AnalysisHandler::doinitrun();

  cout<<"Check2.\n";
  // create ROOT File
  herwigFile = new TFile ("herwig_particles.root","RECREATE");
  herwigFile->SetCompressionLevel(1); // by default file is compressed 

  if (!herwigFile) {
    cout << "StoreParticles: root file has not been created..." << endl;
    return;
  }
  
  // create ROOT Tree
  herwigTree = new TTree ("HerwigTree","Tree filled with herwig data.");
  if (!herwigTree) {
    cout << "StoreParticles: root tree has not been created..." << endl;
    return;
  }
  //herwigTree->SetAutoSave(1000000000); /* autosave when 1 Gbyte written */
  //herwigTree->SetCacheSize(10000000);  /* set a 10 MBytes cache (useless when writing local files) */
  cout<<"Check1.\n";
  //TTree::SetBranchStyle(1); /* new style by default */
  //pEvent = new PrtclEvent;
  herwigTree->Branch("n", &N, "n/s");
  herwigTree->Branch("weight", &weight, "weight/F");
  herwigTree->Branch("id", id, "id[n]/I");
  herwigTree->Branch("status", status, "status[n]/b");
  herwigTree->Branch("pT", pT, "pT[n]/F");
  herwigTree->Branch("eta", eta, "eta[n]/F");
  herwigTree->Branch("phi", phi, "phi[n]/F");
  herwigTree->Branch("m", m, "m[n]/F");
  //TBranch *branch = herwigTree->Branch("event", &pEvent, 32000,4);
  //branch->SetAutoDelete(kFALSE);
  //herwigTree->BranchRef();
}

/* *** Attention *** The following static variable is needed for the type
 * description system in ThePEG. Please check that the template arguments
 * are correct (the class and its base class), and that the constructor
 * arguments are correct (the class name and the name of the dynamically
 * loadable library where the class implementation can be found). */
DescribeNoPIOClass<StoreParticles,AnalysisHandler>
  describejetanalysisStoreParticles("jetanalysis::StoreParticles", "../lib/libStoreParticles.so");

