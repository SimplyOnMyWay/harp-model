# asrfe envelope

global SR = 44100

## usage: pole = tau2pole (tau)
##
##
function pole = tau2pole (tau)
  SR = 44100;
  pole = exp(-1 / (tau * SR));
endfunction
p = tau2pole(0.68);
disp(["p = " p]);
printf("p  =  %d\n",p)

## usage: y = smooth (s)
##
##
function y = smooth (s)
  
endfunction


## usage: asrfe (attT60,susLvl,relT60,finLvl,gate)
##
##
function asrfe (attT60,susLvl,relT60,finLvl,gate)
  
endfunction
