## NOT QUITE THERE!  FUNCTIONS ARE AT FIRST REV, NEEDS A THROUGH REVIEW AS OUTPUTS NOT CORRECT, YET!



clear all;




## asrfe envelope

global SR = 44100;
global fRec0 = [0,0];

## usage: pole = tau2pole (tau)
##
##
function pole = tau2pole (tau)
  SR = 44100;
  pole = exp(-1 / (tau * SR));
endfunction


## usage: y = smooth (s)
##
##
function y = smooth (s,x)
  fRec0 = [0,0];
  for n = 1:length(x)
    fRec0(1) = (1 - s)*x(n) + s*fRec0(2);
    y(n) = fRec0(1);
    fRec0(2) = fRec0(1);
  endfor
endfunction


## usage: envelope = asrfe (attT60,susLvl,relT60,finLvl,gate)
##
##
function envelope = asrfe (attT60,susLvl,relT60,finLvl,gate)
  ugate = gt(gate,0);
  if (ugate == 0)
    target = finLvl;
    t60 = attT60;
  elseif (ugate == 1)
    target = susLvl;
    t60 = relT60;
  endif
  pole = tau2pole(t60/6.91);
  envelope = smooth(pole, target);
endfunction


gateSignal1 = ones(1, 10);
gateSignal0 = zeros(1, 10);
gateSignal = horzcat(gateSignal1,gateSignal0);

for i = 1:length(gateSignal)
  env(i) = asrfe(0.01,0.5,0.02,0.25,gateSignal(i));
endfor


plot(1:length(gateSignal),env)
