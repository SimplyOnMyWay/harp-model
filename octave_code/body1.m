clear all; close all;

filename = "gtrbody.wav";
[signal,fs] = audioread(filename);
sound(signal);

Nfft = 2^16;
Sfull = fft(signal,Nfft);
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
## Smag = 10.^(SdB./20);
## Sr = Smag.*(cos(pS)+i*sin(pS));
## Sfullr = [Sr' flip(Sr(2:end-1)')]';

figure;
set(gcf, 'Position', get(0, 'Screensize'));
plot(fk,[Sfull,Sfulli]);
legend("linear freq","Bark freq");

# inverse fft the full (up to fs) "bark" warped FR
signali = ifft(Sfulli);
signalr = ifft(Sfull);

# time axes...
T = 1/fs;
t = T:T:T*length(signal);
figure;
set(gcf, 'Position', get(0, 'Screensize'));
plot(t,[signal,signalr(1:length(signalr)/8),signali(1:length(signalr)/8)]);
grid minor on;
legend("original","recovered","Bark warped!");

xlimval = 400;
df = fkk(2) - fkk(1);
ixlim = find(abs(fkk - xlimval) < df/2);
figure;
set(gcf, 'Position', get(0, 'Screensize'));
subplot(211);
plot(fkk,SdBNorm);
## xlim([0,fkk(ixlim)]);
ylim([-60,0]);
grid minor on;
subplot(212);
plot(fkkB,SdBNorm);
## xlim([0,fkkB(ixlim)]);
ylim([-60,0]);
grid minor on;








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



