/********************************************/  
/* main.cpp                                 */
/*  Test main program for SimpString.cpp    */
/*  Compatible with STK version 4.2.1.      */
/*  Julius Smith, 2000-2005.                */
/********************************************/  

#include "FileWvOut.h"
#include "SimpString.h"
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

  StkFloat srate = Stk::sampleRate();

  std::cout << "srate = " << srate << std::endl;
  
  simpString->noteOn(370,
                     1.0 /* amplitude */, 
                     0.75 /* pluck position (0:1) */
                     );

  output.tick(simpString->tick(2.0)); /* impulse */


  for (i=1;i<srate*5;i++)   {
    output.tick(simpString->tick(0.0));
  }

  delete simpString;
  return(0);
}
