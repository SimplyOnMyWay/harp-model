## NOT QUITE THERE!  FUNCTIONS ARE AT FIRST REV, NEEDS A THROUGH REVIEW AS OUTPUTS NOT CORRECT, YET!



clear all;
## close all;                      

load "faustout.mat" faustout;


printf("running...!\n");
## asrfe envelope

global SR = 48000;
global fRec0 = [0,0];


global zeroCounter = 0;
global oneCounter = 0;  


## output_precision(16);



## usage: pole = tau2pole (tau)
##
##
function pole = tau2pole (tau)
  global SR;
  pole = exp(-1.0 / (tau * SR));
endfunction


## usage: y = smooth (s)
##
##
function y = smooth (s,x)
  global fRec0
  for n = 1:length(x)
    fRec0(1) = (1 - s)*x(n) + s*fRec0(2);
    ## fRec0(1) = s*(fRec0(2) - x(n)) + x(n);
    y(n) = fRec0(1);
    fRec0(2) = fRec0(1);
  endfor
endfunction


## usage: envelope = asrfe (attT60,susLvl,relT60,finLvl,gate)
##
##
function envASRFE = asrfe (attT60,susLvl,relT60,finLvl,gate)
  global zeroCounter;
  global oneCounter;  
  ugate = gt(gate,0);
  if (ugate == 0)
    zeroCounter++;
    target = finLvl;
    t60 = relT60;
  elseif (ugate == 1)
    target = susLvl;
    t60 = attT60;
    oneCounter++;
  endif
  pole = tau2pole(t60/6.91);
  envASRFE = smooth(pole, target);
endfunction


## usage: envARFE = arfe (attT60,relT60,fv,gate)
##
##
function envARFE = arfe (attT60,relT60,fv,gate)
  envARFE = asrfe(attT60,1.0,relT60,fv,gate);
endfunction

## usage: envARE = are (attT60,relT60,gate)
##
##
function envARE = are (attT60,relT60,gate)
  envARE = arfe(attT60,relT60,0.0,gate);
endfunction


len = 8487;
a_ = 0.005;
r_ = (len - (a_*SR))/SR; #0.1718;
s_ = 0.1;
f_ = 0.05;
gateSignal1 = ones(1, round(a_*SR));
gateSignal0 = zeros(1, round(r_*SR));
gateSignal = horzcat(gateSignal1,gateSignal0);

for i = 1:length(gateSignal)
#  envASRFE(i) = asrfe(a_,s_,r_,f_,gateSignal(i));
#  envARFE(i) = arfe(a_,r_,f_,gateSignal(i));
  envARE(i) = are(a_,r_,gateSignal(i));  
endfor

tvec = [1:length(gateSignal)]';
plot(tvec,[gateSignal',envARE', faustout]);
ylim([-1.0,1.0]);
xlim([0,len]);
grid on; grid minor on;


printf("zeroCounter = %d\n",zeroCounter);
printf("oneCounter = %d\n",oneCounter);
