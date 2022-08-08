clear all; close all;

filename = "gtrbody.wav";
[signal,fs] = audioread(filename);
sound(signal);

frameSizeMS = 100;
minFrameLen = fs*frameSizeMS/1000; 
frameLenPow = nextpow2(minFrameLen);
frameLen = 2^frameLenPow; % frame length = fft size
[B,F,T] = specgram(signal,frameLen,fs);

BdB = 10*log10(abs(B));

figure;
set(gcf, 'Position', get(0, 'Screensize'));
mesh(T,F/1000,BdB);
view(130,30);
title('Spectrogram');
xlabel('Time (s)');ylabel('Frequency (kHz)');zlabel('Magnitude (dB)');
axis tight;zoom on;



## ## bark warping
W = F./F(end);
[WB,G] = lin2bark(W,fs);
FB = WB.*F(end);

iT = 2;

figure;
set(gcf, 'Position', get(0, 'Screensize'));
subplot(211);
plot(F,BdB(:,iT));
grid minor on;
subplot(212);
plot(FB,BdB(:,iT));
grid minor on;





Nfft = 2^10;
Sfull = fft(signal,Nfft);
iposFreq = 1:Nfft/2+1;
S = Sfull(iposFreq);
SdB = 20*log10(abs(S));
SdB10 = 10*log10(abs(S));
offset = max(SdB);
SdBNorm = SdB - offset;
SdBNorm10 = SdB10 - max(SdB10);
fk = linspace(0,fs,Nfft);
fkk = fk(iposFreq);

## ## bark warping
WKK = fkk./fkk(end);
[WKKB,GG] = lin2bark(WKK,fs);
fkkB = WKKB.*fkk(end);

figure;
set(gcf, 'Position', get(0, 'Screensize'));
subplot(211);
plot(fkk,[SdBNorm SdBNorm10]);
grid minor on;
subplot(212);
plot(fkkB,SdBNorm);
grid minor on;


