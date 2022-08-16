clear all; close all;

doplot_if1 = 0; # whether (1) or not (0) to plot invfilter method 1 plots

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
  S_ = S(ibw);
  wkk_ = f(ibw)./f(end).*pi;
  [B,A] = invfreqz(S_,wkk_,2,2);              
endfunction

## usage: A_ = stabilise (A)
##
##
function A_ = stabilise (A)
  A_roots = roots(A);
  P1 = A_roots(1); 
  P1refl = P1; #only changes if unstable
  P2 = A_roots(2); 
  P2refl = P2; #only changes if unstable
  if abs(P1) > 1.0
    P1refl = 1/conj(P1);  
    disp("P1 reflected!");
  endif
  if abs(P2) > 1.0
    P2refl = 1/conj(P2);
    disp("P2 reflected!");
  endif
  A_stabilised_roots = [P1refl,P2refl];
  A_ = poly(A_stabilised_roots); ; #omitting any scaling (see stackexchange.com/questions/26114/how-to-stabilize-a-filter              
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
SfullidB = 20*log10(Sfulli);

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

if (doplot_if1)
  figure;
  set(gcf, 'Position', get(0, 'Screensize'));
  plot(fkk./1000,[SdBNormi,S1dBNorm]);
  axis([0 1.4 -60 0]);
  grid on;
endif


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

if (doplot_if1)
  figure;
  set(gcf, 'Position', get(0, 'Screensize'));
  plot(fkk./1000,[SdBNormi,S1dBNorm,S2dBNorm]);
  axis([0 1.4 -60 0]);
  grid on;
endif


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

if (doplot_if1)
  figure;
  set(gcf, 'Position', get(0, 'Screensize'));
  plot(fkk./1000,[SdBNormi,S1dBNorm,S2dBNorm,S3dBNorm]);
  axis([0 1.4 -60 0]);
  grid on;
  title("this one!");
endif





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
if (doplot_if1)
  figure;
  plot(fkk/1000,[H1dBNorm,H1invdBNorm,H2dBNorm,H2invdBNorm]);
endif


if (doplot_if1)
  ## plot filtered signals in time domain
  figure;
  set(gcf, 'Position', get(0, 'Screensize'));
  itime = 1:length(res1)/divider;
  plot(t.*1000,[signali(itime), res1(itime), res2(itime), res3(itime)]);
  legend("warped original","mode1 removed","mode2 removed","mode3 removed");
  axis([0 30 -1 1]);
  grid minor on;
endif




## ###################################
## alternative approach using invfreqz
## ###################################


# mode 1

## avoid fitting the -ve values of Si up to ca 300Hz
## i300Hz = find(abs(fkk-300)<df/2); 

## convert Bark warped FR to minimum phase. Seems to really help fit, and also appears to make TF invertible without adjusting poles for stability!
Ns = Nfft/2+1;
## Si_mpp = minphase(fk,Nfft,Ns,SfullidB);
## Si_mp = mps(Sfulli);

Si_mpp = Si'; #mp(1:Nfft/2+1)';

## figure;
## set(gcf, 'Position', get(0, 'Screensize'));
## Si_mppdB = 20*log10(abs(Si_mpp));
## subplot(211);
## plot(fkk,Si_mppdB);
## hold on;
## plot(fkk,SdBi);
## ylabel("Magnitude [dB]");
## legend("Si_{mpp}","SdBi");
## subplot(212);
## plot(fkk,rad2deg(unwrap(angle(Si_mpp))));
## hold on;
## plot(fkk,rad2deg(unwrap(angle(Si))));
## ylabel("Phase [degrees]");
## legend("Si_{mpp}","SdBi");


[B1_,A1_] = modes_invfreq(Si_mpp,fkk,freq1,bw1,2);
H1_ = freqz(B1_,A1_,fkk,fs)';
H1_dB = 20*log10(abs(H1_));
H1_dBNorm = H1_dB - offset; 
H1_inv = freqz(A1_,B1_,fkk,fs)';
H1_invdB = 20*log10(abs(H1_inv));
H1_invdBNorm = H1_invdB - offset; 
## figure;
## set(gcf, 'Position', get(0, 'Screensize'));
## plot(fkk,[SdBNormi, H1_dBNorm,H1_invdBNorm]);
## grid minor on;

## B1_ is a culprit as not stable - zero outside unit circle!!!! inverting not good :) Reflect inside to stabilise...
## may not be needed due to converting FR to min phase!
B1_s = stabilise(B1_);

## H1_sinv = freqz(A1_,B1_s,fkk,fs)'; #originally used pre mps
H1_sinv = freqz(A1_,B1_,fkk,fs)';
H1_sinvdB = 20*log10(abs(H1_sinv));
H1_sinvdBNorm = H1_sinvdB - offset;

H1_s = freqz(B1_,A1_,fkk,fs)';
H1_sdB = 20*log10(abs(H1_s));
H1_sdBNorm = H1_sdB - offset;

res1_ = filter(A1_,B1_,signali(1:length(signal)));       % apply inverse filter
Sfull1_ = fft(res1_,Nfft);
S1_ = Sfull1_(iposFreq);
S1_dB = 20*log10(abs(S1_));
S1_dBNorm = S1_dB - max(S1_dB); #offset;

figure;
set(gcf, 'Position', get(0, 'Screensize'));
plot(fkk./1000,[SdBNormi,S1_dBNorm,H1_dBNorm,H1_sdBNorm,H1_invdBNorm,H1_sinvdBNorm]);
legend("SdBNormi","S1_{dBNorm}","H1_{dBNorm}","H1_{sdBNorm}","H1_{invdBNorm}","H1_{sinvdBNorm}");
axis([0 1.4 -60 0]);
grid minor on;
title("Inverse filtering using invfreq!");


## mode 2


[B2_,A2_] = modes_invfreq(Si_mpp,fkk,freq2,bw2,0.5);
H2_ = freqz(B2_,A2_,fkk,fs)';
H2_dB = 20*log10(abs(H2_));
H2_dBNorm = H2_dB - offset; 
H2_inv = freqz(A2_,B2_,fkk,fs)';
H2_invdB = 20*log10(abs(H2_inv));
H2_invdBNorm = H2_invdB - offset; 
## figure;
## set(gcf, 'Position', get(0, 'Screensize'));
## plot(fkk,[SdBNormi, H2_dBNorm,H2_invdBNorm]);
## grid minor on;

## B2_ is a culprit as not stable - zero outside unit circle!!!! inverting not good :) Reflect inside to stabilise...
B2_s = stabilise(B2_);

H2_sinv = freqz(A2_,B2_,fkk,fs)';
H2_sinvdB = 20*log10(abs(H2_sinv));
H2_sinvdBNorm = H2_sinvdB - offset;

H2_s = freqz(B2_,A2_,fkk,fs)';
H2_sdB = 20*log10(abs(H2_s));
H2_sdBNorm = H2_sdB - offset;

res2_ = filter(A2_,B2_,res1_(1:length(signal)));       % apply inverse filter
Sfull2_ = fft(res2_,Nfft);
S2_ = Sfull2_(iposFreq);
S2_dB = 20*log10(abs(S2_));
S2_dBNorm = S2_dB - offset;

figure;
set(gcf, 'Position', get(0, 'Screensize'));
plot(fkk./1000,[SdBNormi,S2_dBNorm,H2_dBNorm,H2_sdBNorm,H2_invdBNorm,H2_sinvdBNorm]);
axis([0 1.4 -60 0]);
grid minor on;
legend("SdBNormi","S2_{dBNorm}","H2_{dBNorm}","H2_{sdBNorm}","H2_{invdBNorm}","H2_{sinvdBNorm}");
title("Inverse filtering using invfreq!");


## ## mode 3


## [B3_,A3_] = modes_invfreq(Si_mpp,fkk,freq3,bw3,0.8);
## H3_ = freqz(B3_,A3_,fkk,fs)';
## H3_dB = 20*log10(abs(H3_));
## H3_dBNorm = H3_dB - offset; 
## H3_inv = freqz(A3_,B3_,fkk,fs)';
## H3_invdB = 20*log10(abs(H3_inv));
## H3_invdBNorm = H3_invdB - offset; 
## ## figure;
## ## set(gcf, 'Position', get(0, 'Screensize'));
## ## plot(fkk,[SdBNormi, H3_dBNorm,H3_invdBNorm]);
## ## grid minor on;

## ## B3_ is a culprit as not stable - zero outside unit circle!!!! inverting not good :) Reflect inside to stabilise...
## B3_s = stabilise(B3_);

## H3_sinv = freqz(A3_,B3_,fkk,fs)';
## H3_sinvdB = 20*log10(abs(H3_sinv));
## H3_sinvdBNorm = H3_sinvdB - offset;

## H3_s = freqz(B3_,A3_,fkk,fs)';
## H3_sdB = 20*log10(abs(H3_s));
## H3_sdBNorm = H3_sdB - offset;

## res3_ = filter(A3_,B3_,res2_(1:length(signal)));       % apply inverse filter
## Sfull3_ = fft(res3_,Nfft);
## S3_ = Sfull3_(iposFreq);
## S3_dB = 20*log10(abs(S3_));
## S3_dBNorm = S3_dB - offset;

## figure;
## set(gcf, 'Position', get(0, 'Screensize'));
## plot(fkk./1000,[SdBNormi,S3_dBNorm,H3_dBNorm,H3_sdBNorm,H3_invdBNorm,H3_sinvdBNorm]);
## axis([0 1.4 -60 0]);
## grid minor on;
## title("Inverse filtering using invfreq!");



## ## plot filtered signals in time domain
## figure;
## set(gcf, 'Position', get(0, 'Screensize'));
## itime = 1:length(res1)/divider;
## plot(t.*1000,[signali(itime), res1_(itime), res2_(itime), res3_(itime)]);
## legend("warped original","mode1 removed","mode2 removed","mode3 removed");
## title("Alternate method using invfreq");
## axis([0 30 -1 1]);
## grid minor on;




## ######
## unwarp
## ######

## FRs

## [S1lin_, W1_, G1_] = bark2lin(S1_, WKKB, fs);
## S1lin_dB = 20*log10(S1lin_);
## S1lin_dBNorm = S1lin_dB - offset;

## [S2lin_, W2_, G2_] = bark2lin(S2_, WKKB, fs);
## pS2lin_dB = 20*log10(S2lin_);
## S2lin_dBNorm = S2lin_dB - offset;

## [S3lin_, W3_, G3_] = bark2lin(S3_, WKKB, fs);
## f3kk = W3_.*fkk(end);
## S3lin_dB = 20*log10(S3lin_);
## S3lin_dBNorm = S3lin_dB - offset;

## [Z1lin_, P1lin_, K1lin_, G1lin_] = bark2lin_(roots(B1), roots(A1_), 1, fs, WKKB);
## ## H1lin_ = freqz(poly(Z1lin_),poly(P1lin_),fkk,fs)';
## B1lin_ = poly(Z1lin_);
## A1lin_ = poly(P1lin_);
## H1lin_ = freqz(K1lin_.*B1lin_,A1lin_,fkk,fs)';
## H1lin_dB = 20*log10(abs(H1lin_));
## H1lin_dBNorm = H1lin_dB - offset;

## [Z2lin_, P2lin_, K2lin_, G2lin_] = bark2lin_(roots(B2_), roots(A2_), 1, fs, WKKB);
## B2lin_ = poly(Z2lin_);
## A2lin_ = poly(P2lin_);
## H2lin_ = freqz(K2lin_.*B2lin_,A2lin_,fkk,fs)';
## H2lin_dB = 20*log10(abs(H2lin_));
## H2lin_dBNorm = H2lin_dB - offset;

## [Z3lin_, P3lin_, K3lin_, G3lin_] = bark2lin_(roots(B3_), roots(A3_), 1, fs, WKKB);
## B3lin_ = poly(Z3lin_);
## A3lin_ = poly(P3lin_);
## H3lin_ = freqz(K3lin_.*B3lin_, A3lin_,fkk,fs)';

## H3lin_dB = 20*log10(abs(H3lin_));
## H3lin_dBNorm = H3lin_dB - offset;

## figure;
## set(gcf, 'Position', get(0, 'Screensize'));
## ## plot(f3kk,S3lin_dBNorm)
## ## hold on;
## plot(fkk./1000,[SdBNorm,S3lin_dBNorm,H1lin_dBNorm,H2lin_dBNorm,H3lin_dBNorm]);
## axis([0 1.4 -60 0]);
## grid minor on;






## IRs
