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

 StkFloat aCoeffs[31] = {1,0.17558,-0.33058,0.40997,0.41358,0.073151,-0.14008,0.11177,-0.097517,-0.027241,0.15791,-0.04225,-0.07367,0.19572,0.025484,0.051416,0.19739,0.10054,0.092412,0.043431,0.10815,0.12631,0.030902,0.028377,0.024654,0.023142,0.031067,0.028542,0.026187,-0.02292,0.016288};
  // StkFloat aCoeffs[31] = {1,0.52874,-0.064736,0.48365,0.65922,0.40633,0.10369,0.25905,0.11832,0.10181,0.23729,0.11339,0.026769,0.24361,0.11236,0.13613,0.29448,0.22862,0.20701,0.17915,0.23339,0.2578,0.18605,0.15629,0.1266,0.10381,0.10485,0.084,0.0772,0.019029,0.016185};
  
  std::vector<StkFloat> a(aCoeffs, aCoeffs+31);
  loopFilter->setDenominator(a);

 StkFloat bCoeffs[31] = {0.90452,0.24799,-0.28885,0.32865,0.4194,0.1261,-0.12793,0.062536,-0.047502,-0.025819,0.11429,-0.012592,-0.073636,0.18,0.018402,0.056438,0.20416,0.095699,0.078081,0.06785,0.10485,0.10778,0.043839,0.035201,0.017689,0.022559,0.047436,0.015413,0.024893,-0.0059326,0.0035131};
  // StkFloat bCoeffs[31] = {0.99002,0.53323,-0.059946,0.47646,0.6579,0.41096,0.10609,0.25464,0.1224,0.1032,0.23355,0.1154,0.027333,0.24254,0.11144,0.13616,0.29518,0.22837,0.20541,0.1811,0.23351,0.25601,0.18682,0.1572,0.12634,0.1038,0.10661,0.083271,0.077115,0.020829,0.01552};
  
  std::vector<StkFloat> b(bCoeffs, bCoeffs+31);
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
  StkFloat loopFilterDelay = -10.0; /* just a guess */
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
