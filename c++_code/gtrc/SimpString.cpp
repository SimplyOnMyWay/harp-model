/******************************************/  
/* Simplified string model                */
/* Compatible with STK version 4.2.1.     */
/* Julius Smith, 2000-2005                */
/******************************************/

#include "SimpString.h"
#include <iostream>

using namespace std;

SimpString :: SimpString(StkFloat lowestFreq)
{
  length = (long) (SRATE / lowestFreq + 1);
  lastLength = length;
  lastFreq = lowestFreq;
  loopGain = (StkFloat) 0.999;
  delayLine = new DelayA(length,length);
  combDelay = new Delay(length,length);
  pluckAmp = (StkFloat) 0.2;
  pluckPos = (StkFloat) 0.25;
  combDelay->setDelay((int) (pluckPos * lastLength + 0.5));

  loopFilter = new Filter;

  // MOC 01.08.2022
  // a and b coeffs copied from outputs of Octave/Matlab file estimateGdb.m 
  StkFloat aCoeffs[26] = {1,1.7168,1.5387,0.85268,0.18961,0.93938,1.7727,1.6546,0.64077,-0.26411,-0.21655,0.23863,0.44181,0.0062434,-0.40779,-0.43513,-0.31064,-0.13673,-0.23406,-0.33322,-0.21455,-0.15972,-0.17757,-0.19113,-0.2026,-0.091798};
  
  std::vector<StkFloat> a(aCoeffs, aCoeffs+26);
  loopFilter->setDenominator(a);

  StkFloat bCoeffs[26] = {-0.51064,0.21721,0.040176,0.064167,0.15054,-0.36574,0.29526,-0.23235,0.28405,0.062146,-0.20363,0.24276,-0.36267,0.30918,-0.056316,-0.0016206,0.10411,-0.22981,0.17852,-0.12055,0.077652,-0.095642,-0.0056801,0.16772,-0.077353,0.067852};
  
  std::vector<StkFloat> b(bCoeffs, bCoeffs+26);
  loopFilter->setNumerator(b);

 /*
  //original loop filter in SimpString.cpp: 
  StkFloat bCoeffs[5] = {0.4,0.2,0.2,0.1,0.1}; 
  std::vector<StkFloat> b(bCoeffs, bCoeffs+5); 
  loopFilter->setNumerator(b);
  */
  
}

SimpString :: ~SimpString()
{
  delete delayLine;
  delete combDelay;
  delete loopFilter;
}

void SimpString :: clear()
{
  delayLine->clear();
  combDelay->clear();
  loopFilter->clear();
}

void SimpString :: setFreq(StkFloat frequency)
{
  StkFloat loopFilterDelay = 0.5; /* just a guess */
  lastFreq = frequency;
  lastLength = ((StkFloat) SRATE / lastFreq); /* length - delays */
  delayLine->setDelay(lastLength - loopFilterDelay);
}

void SimpString :: setPluckPos(StkFloat position)
{ /* Set Pick Position  */
  pluckPos = position;
  combDelay->setDelay((int) (pluckPos * lastLength + 0.5));
}

void SimpString :: noteOff(StkFloat amp)
{
  loopGain = (StkFloat) (1.0 - amp) * 0.5; // Accelerate decay
#if defined(_debug_)        
  printf("SimpString : NoteOff: Amp=%lf\n",amp);
#endif    
}

void SimpString :: pluck(StkFloat amplitude)
{
  pluckAmp = amplitude; /* actual excitation is by tick(input) */
}

void SimpString :: pluck(StkFloat amplitude, StkFloat position)
{ 
  pluckPos = ((position<0) ? 0.0 : ((position>1) ? 1.0 : position));
  combDelay->setDelay((int) (pluckPos * lastLength + 0.5)); 
  this->pluck(amplitude);
}

void SimpString :: noteOn(StkFloat freq, StkFloat amp)
{
  this->setFreq(freq);
  this->pluck(amp);
#if defined(_debug_)        
  printf("SimpString : NoteOn: Freq=%lf Amp=%lf\n",freq,amp);
#endif    
}

void SimpString :: noteOn(StkFloat freq, StkFloat amp, StkFloat position)
{
  this->setFreq(freq);
  this->pluck(amp,position);
#if defined(_debug_)        
  printf("SimpString : NoteOn: Freq=%lf Amp=%lf Pos=%lf\n",
         freq,amp,position);
#endif    
}

StkFloat SimpString :: computeSample(StkFloat stringInput)
{
  /* pluck-position comb filtering: */
  stringInput -= combDelay->tick(stringInput);
  lastOutput_ = delayLine->tick(loopFilter->tick( stringInput + (delayLine->lastOut() * loopGain)));

//  lastOutput_ = delayLine->tick(stringInput + (delayLine->lastOut()) * loopGain);    
  return lastOutput_;
}

StkFloat SimpString :: tick(StkFloat stringInput)
{
  return computeSample(stringInput);
}

StkFloat SimpString :: computeSample(void)
{
  return this->tick(0);
}
