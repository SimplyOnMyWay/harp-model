/******************************************/ 
/* Simplified string model.               */ 
/* Compatible with STK version 4.2.1.     */
/* Julius Smith, 2000-2005                */
/******************************************/

#if !defined(__SimpString_h)
#define __SimpString_h

#include "Instrmnt.h"
#include "Delay.h"
#include "DelayA.h"
#include "Filter.h"

class SimpString : public Instrmnt
{
 protected:  
  DelayA *delayLine;  // main string
  Delay *combDelay;   // pick-position simulation
  Filter *loopFilter; // string losses
  long length;
  StkFloat loopGain;
  StkFloat lastFreq;
  StkFloat lastLength;
  StkFloat pluckAmp;
  StkFloat pluckPos;
  StkFloat computeSample(StkFloat stringInput);
  StkFloat computeSample(void);
 public:
  SimpString(StkFloat lowestFreq);
  ~SimpString();
  void clear();
  void setPluckPos(StkFloat position);
  void pluck(StkFloat amplitude);
  void setFreq(StkFloat frequency);
  void pluck(StkFloat amplitude,StkFloat position);
  void noteOn(StkFloat freq, StkFloat amp);
  void noteOn(StkFloat freq, StkFloat amp, StkFloat position);
  StkFloat tick(StkFloat stringInput);
  StkFloat tick(void);
  void noteOff(StkFloat amp);
  void setBearing(StkFloat angle); /* -45 to +45 degrees */
};

#endif

