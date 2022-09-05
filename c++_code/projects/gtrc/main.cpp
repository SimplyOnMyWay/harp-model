/********************************************/  
/* main.cpp                                 */
/*  Test main program for SimpString.cpp    */
/*  Compatible with STK version 4.2.1.      */
/*  Julius Smith, 2000-2005.                */
/********************************************/  

#include "FileWvIn.h"
#include "FileWvOut.h"
#include "SimpString.h"
#include "ADSR.h"
#include "Stk.h"
#include <stdlib.h>
#include <ctype.h>
#include <iostream>

int main(int argc,char *argv[])
{

  Stk::setSampleRate(48000);

  long i;
  FileWvOut output(argv[0]); /* creates output soundfile */
  SimpString *simpString =  new SimpString(/* lowest pitch */ 50);

  //StkFloat srate = Stk::sampleRate();

  simpString->noteOn(370,
                     1.0 /* amplitude */, 
                     0.75 /* pluck position (0:1) */
                     );

  // generate excitation as noise-burst with IR and FR shaped according to body recording - see body3.m
  ADSR env;

  env.keyOn();
  env.setAllTimes(0.001,0.02,0.03,0.0);


  FileWvIn imp_temp("../../../octave_code/e_sig.wav");
  StkFloat fs = (long) imp_temp.getFileRate();
  Stk::setSampleRate(fs); // set sampling rate to that of input
  std::cout << "srate = " << fs << std::endl;
  FileWvIn imp("../../../octave_code/e_sig.wav");
  StkFloat computeSeconds = 2;
  long nSamps = (long) (fs*computeSeconds);
  //  StkFloat forOutput[nSamps]; // output buffer


  StkFloat noiseburst[nSamps];
  for (i=0;i<nSamps;i++){
    noiseburst[i] = (rand()-0.5)*2; //*env.tick();
}
  env.keyOff();


  output.tick(simpString->tick(1));
  
  for (i=0;i<nSamps;i++)   {
       // output.tick(simpString->tick(imp.tick()));
    //output.tick((rand()-0.5)*2);
     output.tick(simpString->tick(0));
  }

  delete simpString;
  return(0);
}
