#ifndef PRTCLEVENT_H
#define PRTCLEVENT_H

///////////////////////////////////////////////////////////////////////
// A generic event class for storing particle data from simulations. //
// Hannu Siikonen 11.3.2015                                          //
// (special thanks to Rene Brun's ROOT examples)                     //
///////////////////////////////////////////////////////////////////////

#include <iostream>
#include <cmath>
#include <cassert>

#include "TROOT.h"
#include "TObject.h"
#include "TClonesArray.h"
#include "TProcessID.h"
#include "TMath.h"
#include "TClass.h"
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"


// The information concerning one particle, stored to a ROOT tree.
class PrtclData : public TObject 
{
public:
   PrtclData() { Class()->IgnoreTObjectStreamer(); }
   virtual ~PrtclData() { }

   void SetFourVector(double,double,double,double);
   void SetParams(int,double,int);
   

private:
   /* Use a pure ROOT LorentzVector so that for instance Pt can be found out
    * even without the sources of this event class. This is a slightly better 
    * format than TLorentzVector and is in use for instance in the KKousouris 
    * scripts (indirectly, through CMSSW). */ 
   //ROOT::Math::LorentzVector< ROOT::Math::PxPyPzE4D<double> > fP4;
  
   int fPDGCode;
   int fChargeTimes3;
   float pT,eta,phi,m;
   
   /* Indicates why the particle was saved and/or other properties it holds.
    * 1: generic particle (stable)
    * 10: pi0 photon (stable)
    * 11: parton/gluon for jet flavour studies
    * 12: b hadron for jet flavour studies
    * 13: c hadron for jet flavour studies
    * 14: leptons for jet flavour studies*/
   int fAnalysisStatus;
   
public:
   ClassDef(PrtclData,1)
};


// An event class for the particle data that is stored
class PrtclEvent : public TObject 
{
public:
   PrtclEvent(size_t = 10000);
   virtual ~PrtclEvent() { Reset(); };

   void AddPrtcl(double,double,double,double,int,double,int);
   PrtclData *InitPrtcl();

   void Clear(Option_t *option ="");
   void Reset(Option_t *option ="");

   TClonesArray *GetParts() const { return fPrtcls; }

private:
   size_t fN_Prtcl; //! Present amount of particles in the tree 
   size_t fSizeLim; //! Maximal amount of particles within an event
  
   TClonesArray *fPrtcls;
   static TClonesArray *fgPrtcls;

public:
   ClassDef(PrtclEvent, 1)
};


#endif // PRTCLEVENT_H