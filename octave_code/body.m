## model of harp body

clear all; close all;

## usage: plotIR (t,signal,fs);
## set plotlims to 0 to show full extents of plot axes
function plotIR (signal,fs,plotlims)
  T = 1/fs;
  t = T.*[1:length(signal)].*1000;
  ## plot impulse response of body
  set(gcf, 'Position', get(0, 'Screensize'));
  plot(t,signal);
  if plotlims
    xlim([0,50]);
  endif
  xlabel("time (ms)");
  ylabel("body impulse response (amplitude)");
  grid minor on;

endfunction

## usage: plotFR (f,s,plotlims)
## set plotlims to 0 to show full extents of plot axes
function plotFR (signal,fs,plotlims)
  ## frequency magnitude response of body
  Nfft = 2^16;
  fk = linspace(0,fs/2,Nfft/2+1);
  sfull = fft(signal, Nfft);
  s = sfull(1:Nfft/2+1);
  sdB = 20*log(abs(s));
  sdBnorm = sdB - max(sdB);

  set(gcf, 'Position', get(0, 'Screensize'));
  plot(fk,sdBnorm,".-k");
  if plotlims
  xlim([0 400]);
  ylim([-60,0]);
  end
  xlabel("frequency (Hz)");
  ylabel("");
  grid minor on;
endfunction

## usage: [residual,A,B] = peakremove (f,bw,signal,fs)
## localised peak removal
## https://ccrma.stanford.edu/realsimple/BodyFactoring/Matlab_Localized_Peak_Removal.html
function [residual,A,B] = peakremove (f,bw,signal,fs)
  R = exp( - pi * bw / fs);            % pole radius
  z = R * exp(j * 2 * pi * f / fs); % pole itself
  B = [1, -(z + conj(z)), z * conj(z)]; % inverse filter numerator
  r = 0.95;  % zero/pole factor (notch isolation)
  A = B .* (r .^ [0 : length(B)-1]);   % inverse filter denominator
  residual = filter(B,A,signal);       % apply inverse filter
endfunction


## ######################
## main part of the code!
## ######################

filename = "gtrbody.wav";
[signal,fs] = audioread(filename);
sound(signal);

## find portion of response where most modes have decayed considerably, ie after 125ms
T = 1/fs;
t = 0:T:(length(signal)-1)*T;
i125ms = find(abs(t-0.025) < T/2); #0.25
i200ms = find(abs(t-0.3) < T/2); #0.3

## figure(9);
## set(gcf, 'Position', get(0, 'Screensize'));
## plot(t,signal);
## hold on; grid minor on;
## plot(t(i125ms:end),signal(i125ms:end));

## figure(10);
## plotlims = 1;
## Nfft = 2^17;
## fk = linspace(0,fs/2,Nfft/2+1);

## ## bark warping
## W = fk./fk(end);
## [WB,G] = lin2bark(W,fs);
## fkb = WB.*fk(end);

## sfull = fft(signal(i125ms:i200ms), Nfft);
## s = sfull(1:Nfft/2+1);
## sdB = 20*log(abs(s));
## sdBnorm = sdB - max(sdB);

## plot(fkb,sdBnorm,".-k");
## set(gcf, 'Position', get(0, 'Screensize'));
## if plotlims
##   xlim([0 1400]);
##   ylim([-60,0]);
## end
## xlabel("Bark warped frequency (Hz)");
## ylabel("");
## grid on;



idx = i125ms:i200ms;

figure(1);plotIR(signal(idx),fs,0);
figure(2);plotFR(signal(idx),fs,1);


## remove 1st mode

freq1 = 104.6;  % estimated peak frequency in Hz
bw1 = 10;        % peak bandwidth estimate in Hz

[res1,A1,B1] = peakremove(freq1,bw1,signal(idx),fs);

figure(3);plotIR(res1,fs,1);title("mode1 removed");
figure(4);plotFR(res1,fs,1);title("mode1 removed");

[H1,w1] = freqz(B1,A1);
f50 = figure(50); plot(w1./pi*fs/2,abs(H1));

## invfreq on 1st peak
# redo the setup
Nfft = 2^16;
fk = linspace(0,fs/2,Nfft/2+1);
sfull = fft(signal(idx), Nfft);
s = sfull(1:Nfft/2+1);
sdB = 20*log(abs(s));
sdBnorm = sdB - max(sdB);

wt = zeros(length(fk));

p


             


## ## ## remove 2nd mode
## freq2 = 205.4;  % estimated peak frequency in Hz
## bw2 = 10;        % peak bandwidth estimate in Hz

## [res2,A2,B2] = peakremove(freq2,bw2,res1,fs);

## [H2,w2] = freqz(B2,A2);
## figure(50); hold on; plot(w2./pi*fs/2,abs(H2));

## figure(5);plotIR(res2,fs,1);title("mode2 removed");
## figure(6);plotFR(res2,fs,1);title("mode2 removed");

## ## remove 4th mode
## freq4 = 394.1;  % estimated peak frequency in Hz
## bw4 = 10;        % peak bandwidth estimate in Hz

## [res4,A4,B4] = peakremove(freq4,bw4,res2,fs);

## figure(7);plotIR(res4,fs,1);title("mode4 removed");
## figure(8);plotFR(res4,fs,1);title("mode4 removed");
