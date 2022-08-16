## usage: [B,A] = modes_invfreq (S,f,fc,bw,bwfactor)
## S: Complex FR, non-negative freqs (0:fs/2)
## f: frequency vector up to fs/2
## fc: peak center frequency
## bw: estimated bandwidth
function [B,A] = modes_invfreq (S,f,fc,bw,bwfactor)
  if nargin < 5
  bwfactor = 2;  #seems to be optimal, 1 and 3 tested, not so good
  endif
  if nargin < 4
    bw = 10;
    endif
  df = f(2)-f(1);
  ibw = find(abs(f - (fc - bwfactor*bw)) < df/2) : ...
         find(abs(f - (fc + bwfactor*bw)) < df/2);
  S_ = S(ibw);
  wkk_ = f(ibw)./f(end).*pi;
  [B,A] = invfreqz(S_,wkk_,2,2);              
endfunction
