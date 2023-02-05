## test algorithm for adsr for use in my artblocks project
## plotting shapes to ensure correct reasoning in the code


close all; clear all;

## usage: [env,t]= ad(a_ms,d_ms,fs)
## inputs:
## a_ms - attack duration in ms
## d_ms - decay duration in ms
## fs - sampling freq
## outpus:
## env - vector of adsr envelope values
## t - vector of time values in ms

function [env,t] = ad (ams,dms,fs)
  a = ams/1000*fs;             # attack duration in samples
  d = dms/1000*fs;             # decay duration in samples

  aslope = 1/a;
  dslope = -1/d;

  for i = 1:a
    env(i) = aslope*i;
    t(i) = i/fs*1000;
  end

  for i = 1:d
    env(a+i) = 1 + dslope*i;
    t(a+i) = (a+i)/fs*1000;
  end

  
endfunction


[env,t] = ad(20,50,48000);

plot(t,env);
xlabel("time (ms)");




