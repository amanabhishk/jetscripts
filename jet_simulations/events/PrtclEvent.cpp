#include "PrtclEvent.h"

ClassImp(PrtclData)
ClassImp(PrtclEvent)


TClonesArray *PrtclEvent::fgPrtcls = 0;

using std::cout;
using std::endl;


void PrtclData::SetFourVector(double a, double b, double c, double d)
{
   pT = a;
   eta = b;
   phi = c;
   m = d;
}


void PrtclData::SetParams(int id, double charge, int status)
{
   fPDGCode = id;
   fChargeTimes3 = TMath::Nint(3*charge);
   fAnalysisStatus = status;
}


PrtclEvent::PrtclEvent(size_t tmpStore)
{
   Class()->IgnoreTObjectStreamer();
   fN_Prtcl = 0;
   if (!fgPrtcls) fgPrtcls = new TClonesArray("PrtclData",tmpStore);
   fSizeLim = tmpStore;
   fPrtcls = fgPrtcls;
}


void PrtclEvent::AddPrtcl(double pT, double eta, double phi, double m, int id, 
   double charge, int status)
{
   int ObjectNumber = TProcessID::GetObjectCount();
  
   PrtclData *part;
   part = InitPrtcl();
   part->SetFourVector(pT,eta,phi,m);
   part->SetParams(id,charge,status);
   
   TProcessID::SetObjectCount(ObjectNumber);
}


PrtclData* PrtclEvent::InitPrtcl()
{
   assert(fSizeLim>fN_Prtcl);
   PrtclData *part = (PrtclData*) fPrtcls->ConstructedAt(fN_Prtcl++);
   return part;
}


void PrtclEvent::Clear(Option_t* /* option */)
{
   fPrtcls->Clear("C");
   fN_Prtcl = 0;
}


void PrtclEvent::Reset(Option_t* option)
{
   delete fgPrtcls; fgPrtcls = 0;
   fN_Prtcl = 0;
}
