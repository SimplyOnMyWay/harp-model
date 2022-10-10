// Based on <faust-0.9.8.6>/faust-examples/faust/karplus.dsp
import("stdfaust.lib");

// MIDI-driven parameters:
freq = nentry("freq", 370, 20, 20000, 1); // Hz
gain = nentry("gain", 1, 0, 10, 0.01);    // 0 to 1
gate = button("gate");                    // 0 or 1

// Excitation window (convert gate to a one-period pulse):
diffgtz(x) = (x-x') > 0;
decay(n,x) = x - (x>0)/n;
release(n) = + ~ decay(n);
trigger(n) = diffgtz : release(n) : > (0.0);

// Resonator:
average(x) = (x+x')/2;
bv = 0.99002,0.53323,-0.059946,0.47646,0.6579,0.41096,0.10609,0.25464,0.1224,0.1032,0.23355,0.1154,0.027333,0.24254,0.11144,0.13616,0.29518,0.22837,0.20541,0.1811,0.23351,0.25601,0.18682,0.1572,0.12634,0.1038,0.10661,0.083271,0.077115,0.020829,0.01552;
// note av[0] = 1 is assumed by Faust!
av = 0.52874,-0.064736,0.48365,0.65922,0.40633,0.10369,0.25905,0.11832,0.10181,0.23729,0.11339,0.026769,0.24361,0.11236,0.13613,0.29448,0.22862,0.20701,0.17915,0.23339,0.2578,0.18605,0.15629,0.1266,0.10381,0.10485,0.084,0.0772,0.019029,0.016185;

//bv = 0.5,0.5;
//av = 0;

loopFilter = fi.iir(bv,av);
P = (ma.SR/freq)-0.5;
//resonator = (+ : de.delay(4096, P)) ~ (average);
dwg = (+ : de.delay(4096, P)) ~ (loopFilter);

imp = 1-1';

process = no.noise : *(gain) : *(gate : trigger(P)) : dwg <: _,_;
//process = imp : dwg <: _,_;