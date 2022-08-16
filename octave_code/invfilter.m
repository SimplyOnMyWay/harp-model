## usage: [A,B] = invfilter (f,bw,r,fs)
## f = center frequency of peak
## bw = bandwidth of peak
## r = zero/pole factor (notch isolation)
## fs = sampling rate (Hz)
function [A,B] = invfilter (f,bw,r,fs)
  R = exp( - pi * bw / fs);            % pole radius
  z = R * exp(j * 2 * pi * f / fs); % pole itself
  B = [1, -(z + conj(z)), z * conj(z)]; % inverse filter numerator
  ## r = 0.95;     % zero/pole factor (notch isolation)
  A = B .* (r .^ [0 : length(B)-1]);   % inverse filter denominator
endfunction
