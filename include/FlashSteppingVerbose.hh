
class FlashSteppingVerbose;

#ifndef FlashSteppingVerbose_h
#define FlashSteppingVerbose_h 1

#include "G4SteppingVerbose.hh"



class FlashSteppingVerbose : public G4SteppingVerbose
{
 public:

   FlashSteppingVerbose();
   virtual ~FlashSteppingVerbose();

   virtual void StepInfo();
   virtual void TrackingStarted();

};



#endif