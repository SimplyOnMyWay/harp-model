clear all; close all;

## usage: phase = ph (S)
##
##
function phase = ph (S)
  phase = rad2deg(unwrap(angle(S)));
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
SdBN = SdB - offset;
fk = linspace(0,fs,Nfft);
fkk = fk(iposFreq);

mps_method = 1;
if mps_method == 1
  Smp = mps(Sfull);
  Smpp = Smp(iposFreq);
elseif mps_method == 2
  Ns = length(iposFreq)
  Smpp = minphase(fk,Nfft,Ns,Sfull);
endif

SmppdB = 20*log10(abs(Smpp));
SmppdBN = SmppdB - offset;

figure;
set(gcf, 'Position', get(0, 'Screensize'));
subplot(211);
plot(fkk,[SdBN,SmppdBN]);
grid minor on;
subplot(212);
plot(fkk,[ph(S),ph(Smpp)]);
grid minor on;


fc1 = 104.7;
bw1 = 10;
fc2 = 205.4;
bw2 = 10;
fc3 = 247.5;
bw3 = 0.4;


[B1_,A1_] = modes_invfreq(S',fkk,fc1,bw1,2);
H1 = freqz(B1_,A1_,fkk,fs);
H1 = mps(H1);
[B1,A1] = invfreqz(H1,fkk./fkk(end).*pi,2,2);
H1dB = 20*log10(abs(H1));
H1dBN = H1dB - offset;

[B2_,A2_] = modes_invfreq(S',fkk,fc2,bw2,2);
H2 = freqz(B2_,A2_,fkk,fs);
H2 = mps(H2);
H2dB = 20*log10(abs(H2));
H2dBN = H2dB - offset;

[B3_,A3_] = modes_invfreq(S',fkk,fc3,bw3,1);
H3 = freqz(B3_,A3_,fkk,fs);
H3 = mps(H3);
H3dB = 20*log10(abs(H3));
H3dBN = H3dB - offset;


figure;
set(gcf, 'Position', get(0, 'Screensize'));
plot(fkk,[SmppdBN,H1dBN',H2dBN',H3dBN']);



