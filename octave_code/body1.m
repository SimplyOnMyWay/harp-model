clear all; close all;

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
  wkk_ = f(ibw)'./f(end).*pi;
  S_ = S(ibw);
  [B,A] = invfreqz(S_,wkk_,2,2);              
endfunction




filename = "gtrbody.wav";
[signal,fs] = audioread(filename);
sound(signal);
# time axes...
T = 1/fs;
t = T:T:T*length(signal);
istart = find(abs(t-0.001)<T/2);
iend = find(abs(t-0.300)<T/2);

Nfft = 2^16;  #I think this is oversampling by factor of 2^3 = 8, but it seems to need it as FR looks jaggedy below 2^16
divider = Nfft/length(signal); #to remove zero pads when plotting ifft recovered signals
Sfull = fft(signal(istart:iend),Nfft);
iposFreq = 1:Nfft/2+1;
S = Sfull(iposFreq);
SdB = 20*log10(abs(S));
offset = max(SdB);
SdBNorm = SdB - offset;
fk = linspace(0,fs,Nfft);
fkk = fk(iposFreq);

## ## bark warping
WKK = fkk./fkk(end);
[WKKB,GG] = lin2bark(WKK,fs);
fkkB = WKKB.*fkk(end);

df = fkk(2)-fkk(1);
SdBNormi = interp1(fkkB,SdBNorm,fkk)';

SdBi = SdBNormi+offset;
Simag = 10.^(SdBi./20);
pS = unwrap(angle(S));
pSi = interp1(fkkB,pS,fkk)';
Si = Simag.*(cos(pSi)+i*sin(pSi));
Sfulli = [Si' flip(Si(2:end-1)')]';

## recover original spectrum from SdB, for referance
Smag = 10.^(SdB./20);
Sr = Smag.*(cos(pS)+i*sin(pS));
Sfullr = [Sr' flip(Sr(2:end-1)')]';

## ## plot the original against the Bark FR
## figure;
## set(gcf, 'Position', get(0, 'Screensize'));
## plot(fk,[Sfull,Sfulli]);
## legend("linear freq","Bark freq");

# inverse fft the full (up to fs) "bark" warped FR
signali = ifft(Sfulli);
signalr = ifft(Sfull);

## ## plot both the recovered and the bark warped against the original
## figure;
## set(gcf, 'Position', get(0, 'Screensize'));

## plot(t,[signal,signalr(1:length(signalr)/divider),signali(1:length(signalr)/divider)]);
## grid minor on;
## legend("original","recovered","Bark warped!");


## ## plot both linear and Bark-warped FR for comparison
## xlimval = 400;
## df = fkk(2) - fkk(1);
## ixlim = find(abs(fkk - xlimval) < df/2);
## figure;
## set(gcf, 'Position', get(0, 'Screensize'));
## subplot(211);
## plot(fkk,SdBNorm);
## ## xlim([0,fkk(ixlim)]);
## ylim([-60,0]);
## grid minor on;
## subplot(212);
## plot(fkkB,SdBNorm);
## ## xlim([0,fkkB(ixlim)]);
## ylim([-60,0]);
## grid minor on;

## ## zoom in on LF portion of Bark-warped FR, this time using fkk as x-axis
## figure;
## set(gcf, 'Position', get(0, 'Screensize'));
## plot(fkk./1000,SdBNormi);
## axis([0 1.4 -60 0]);
## grid on;


## ########################
## inverse filtering
## ########################

## inverse filtering mode 1

freq1 = 469;  % estimated peak frequency in Hz
bw1 = 10;        % peak bandwidth estimate in Hz
r1 = 0.97;
[A1,B1] = invfilter(freq1,bw1,r1,fs);
res1 = filter(B1,A1,signali);       % apply inverse filter

Sfull1 = fft(res1,Nfft);
S1 = Sfull1(iposFreq);
S1dB = 20*log10(abs(S1));
S1dBNorm = S1dB - offset;

figure;
set(gcf, 'Position', get(0, 'Screensize'));
plot(fkk./1000,[SdBNormi,S1dBNorm]);
axis([0 1.4 -60 0]);
grid on;


## inverse filtering mode 2

freq2 = 915;  % estimated peak frequency in Hz
bw2 = 20;        % peak bandwidth estimate in Hz
r2 = 0.9;

[A2,B2] = invfilter(freq2,bw2,r2,fs);

res2 = filter(B2,A2,res1);       % apply inverse filter

Sfull2 = fft(res2,Nfft);
S2 = Sfull2(iposFreq);
S2dB = 20*log10(abs(S2));
## offset = max(S1dB);
S2dBNorm = S2dB - offset;

figure;
set(gcf, 'Position', get(0, 'Screensize'));
plot(fkk./1000,[SdBNormi,S1dBNorm,S2dBNorm]);
axis([0 1.4 -60 0]);
grid on;


## inverse filtering mode 3

freq3 = 1102;  % estimated peak frequency in Hz
bw3 = 20;        % peak bandwidth estimate in Hz
r3 = 0.9;

[A3,B3] = invfilter(freq3,bw3,r3,fs);

res3 = filter(B3,A3,res2);       % apply inverse filter

Sfull3 = fft(res3,Nfft);
S3 = Sfull3(iposFreq);
S3dB = 20*log10(abs(S3));
## offset = max(S1dB);
S3dBNorm = S3dB - offset;

figure;
set(gcf, 'Position', get(0, 'Screensize'));
plot(fkk./1000,[SdBNormi,S1dBNorm,S2dBNorm,S3dBNorm]);
axis([0 1.4 -60 0]);
grid on;
title("this one!");






##  overlay filters and their inverse
[H1,w1] = freqz(A1,B1,Nfft/2+1);
H1dB = 20.*log10(abs(H1));
H1dBNorm = H1dB - offset; #max(H1dB);
[H1inv,w1inv] = freqz(B1,A1,Nfft/2+1);
H1invdB = 20.*log10(abs(H1inv));
H1invdBNorm = H1invdB - offset; #max(H1dB);

[H2,w2] = freqz(A2,B2,Nfft/2+1);
H2dB = 20.*log10(abs(H2));
H2dBNorm = H2dB - offset; #max(H2dB);
[H2inv,w2inv] = freqz(B2,A2,Nfft/2+1);
H2invdB = 20.*log10(abs(H2inv));
H2invdBNorm = H2invdB - offset; #max(H2dB);
## hold on;
figure;
plot(fkk/1000,[H1dBNorm,H1invdBNorm,H2dBNorm,H2invdBNorm]);



## plot filtered signals in time domain
figure;
set(gcf, 'Position', get(0, 'Screensize'));
itime = 1:length(res1)/divider;
## plot(t,[signali(itime), res1(itime), res2(itime)]);
plot(t,[signali(itime), res1(itime), res2(itime), res3(itime)]);
legend("warped original","mode1 removed","mode2 removed","mode3 removed");






## ## test inverse fft on SdB and S
## signalinv1 = ifft(Sfull);

## T = 1/fs;
## t = T.*[1:length(signal)];
## Tinv = t(end)/length(signalinv1);
## tinv = 1/(fs*8):1/(fs*8):t(end);
## figure;
## plot(t,signal,tinv,signalinv1);

## figure;
## set(gcf, 'Position', get(0, 'Screensize'));
## hax1 = subplot(211);
## plot(t,signal);
## grid minor on;
## hax2 = subplot(212);
## ii = 1:round(length(tinv)/8);
## plot(tinv(ii),signalinv1(ii));
## grid minor on;
## linkaxes([hax1, hax2]);
## xlim([0,0.1])


## sound(signalinv1,fs)

## sound(signal,fs)



## ###################################
## alternative approach using invfreqz
## ###################################

[B1_,A1_] = modes_invfreq(Si,fkk,freq1,bw1,2);
H1_ = freqz(B1_,A1_,fkk,fs)';
H1_dB = 20*log10(abs(H1_));
H1_dBNorm = H1_dB - offset; 
H1_inv = freqz(A1_,B1_,fkk,fs)';
H1_invdB = 20*log10(abs(H1_inv));
H1_invdBNorm = H1_invdB - offset; 
figure;
set(gcf, 'Position', get(0, 'Screensize'));
plot(fkk,[SdBNormi, H1_dBNorm,H1invdBNorm]);
grid minor on;

## culprit!!!! not so straight forward to invert the biquad filter it seems, hopefully a workaround is possible!
res1_ = filter(A1_,B1_,signali(1:length(signal)));       % apply inverse filter

Sfull1_ = fft(res1_,Nfft);
S1_ = Sfull1_(iposFreq);
S1_dB = 20*log10(abs(S1_));
S1_dBNorm = S1_dB - offset;

figure;
set(gcf, 'Position', get(0, 'Screensize'));
plot(fkk./1000,[SdBNormi,S1_dBNorm]);
axis([0 1.4 -60 0]);
grid on;
title("Inverse filtering using invfreq!");



