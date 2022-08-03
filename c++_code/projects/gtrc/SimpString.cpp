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
  StkFloat aCoeffs[26] = {1,1.0121,-0.0033335,0.020044,0.065425,0.010869,0.27115,0.36892,0.099736,-0.31453,-0.09778,0.55836,0.35303,0.1377,0.11322,0.10453,0.28679,0.26397,0.12531,0.035282,-0.12458,0.069552,0.2519,0.20875,0.23293,0.07192};
  
  std::vector<StkFloat> a(aCoeffs, aCoeffs+26);
  loopFilter->setDenominator(a);

  StkFloat bCoeffs[26] = {-0.41846,0.2904,0.25022,-0.37259,0.25113,0.15289,-0.27758,-0.00741,0.23392,-0.023239,-0.25541,0.014386,0.18875,-0.17571,-0.0048375,0.17101,-0.019006,-0.11101,0.0031186,0.24299,-0.069115,-0.2156,0.13133,-6.1139e-05,-0.1062,0.10544};
  
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
