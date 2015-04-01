
#include "pythia6_functions.h"


int main(int argc, char** argv)
{
   size_t nEvent = 400;
   if (argc > 1) nEvent = atoi(argv[1]); assert( nEvent > 0 );
   bool addLeptons = false;
   
   if (nEvent > 0) {
      return makeEventSample(nEvent);
   } else {
      return 1;
   }
}
