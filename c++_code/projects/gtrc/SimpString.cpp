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
  loopGain = (StkFloat) 1.0;
  delayLine = new DelayA(length,length);
  combDelay = new Delay(length,length);
  pluckAmp = (StkFloat) 0.2;
  pluckPos = (StkFloat) 0.25;
  combDelay->setDelay((int) (pluckPos * lastLength + 0.5));

  loopFilter = new Filter;

  // MOC 01.08.2022
  // a and b coeffs copied from outputs of Octave/Matlab file estimateGdb.m 
  StkFloat aCoeffs[26] = {1,0.84957,-0.23412,-0.086182,0.15575,0.14016,0.29239,0.30141,-0.083032,-0.40584,-0.021339,0.58943,0.24425,-0.045438,-0.0074434,0.079901,0.23753,0.24937,0.1673,-0.021739,-0.19302,0.04546,0.20634,0.10795,0.11098,0.0023504};
  
  std::vector<StkFloat> a(aCoeffs, aCoeffs+26);
  loopFilter->setDenominator(a);

  StkFloat bCoeffs[26] = {0.9113,0.86472,-0.13591,-0.11715,0.11981,0.16158,0.29417,0.28849,-0.048681,-0.37634,-0.067495,0.54502,0.28484,-0.021092,-0.028164,0.060289,0.2464,0.26189,0.14685,-0.0020269,-0.16457,0.016445,0.19475,0.12125,0.10997,0.0093465};
  
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
  StkFloat loopFilterDelay = -10; /* just a guess */
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
