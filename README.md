# harp-model

sound file cliag_G4_3.wav contains pluck of my early Gaelic harp, note G4 370Hz

octave file estimateGdb.m estimates the loop filter gains gdb, and fits a order 25 IIR filter to it using invfreq.
<br/>
all other .m files are dependencies to estimateGdb.m

NOTE: currently the results of converting gdb to minimum phase looks strange!

the outputs of estimateGdb.m are the A and B coefficients for the IIR filter
<br/>
these are copied / pasted into the c++ file SimpString.cpp (line 26 - 36); the file remains unchanged otherwise

NOTE: currently running "make test" produces a WAV file which sounds completely wrong! I suspect this is either due to<br/>
a) the minimum phase conversion of the loop filter freq response above being incorrect<br/>
b) something I am overlooing in how I have adjusted the c++ DWG model of a string with order 3 FIR loop filter, to an order 25 IIR loop filter.

<br/>

**Instructions:**
## octave code
open octave<br/>
on octave command line, enter "run estimateGdb.m" <br/>
to suppress line fit plots set variable doplot = 0 on line 122<br/>

## c code
navigate to folder projects/gtrc <br/>
on command line enter "make test" <br/>

