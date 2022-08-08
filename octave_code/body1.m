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

## recover original spectrum from SdB, for referance
Smag = 10.^(SdB./20);
Sr = Smag.*(cos(pS)+i*sin(pS));
Sfullr = [Sr' flip(Sr(2:end-1)')]';

## plot the original against the Bark FR
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

## plot both the recovered and the bark warped against the original
figure;
set(gcf, 'Position', get(0, 'Screensize'));
plot(t,[signal,signalr(1:length(signalr)/8),signali(1:length(signalr)/8)]);
grid minor on;
legend("original","recovered","Bark warped!");


## plot both linear and Bark-warped FR for comparison
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

## zoom in on LF portion of Bark-warped FR, this time using fkk as x-axis
figure;
set(gcf, 'Position', get(0, 'Screensize'));
plot(fkk./1000,SdBNormi);
axis([0 1.4 -60 0]);
grid on;

## inverse filtering mode 1

freq = 469;  % estimated peak frequency in Hz
bw = 10;        % peak bandwidth estimate in Hz

R = exp( - pi * bw / fs);            % pole radius
z = R * exp(j * 2 * pi * freq / fs); % pole itself
B = [1, -(z + conj(z)), z * conj(z)] % inverse filter numerator
r = 0.95;     % zero/pole factor (notch isolation)
A = B .* (r .^ [0 : length(B)-1]);   % inverse filter denominator

residual = filter(B,A,signali);       % apply inverse filter

Sfull1 = fft(residual,Nfft);
S1 = Sfull1(iposFreq);
S1dB = 20*log10(abs(S1));
## offset = max(S1dB);
S1dBNorm = S1dB - offset;

figure;
set(gcf, 'Position', get(0, 'Screensize'));
plot(fkk./1000,[SdBNormi,S1dBNorm]);
axis([0 1.4 -60 0]);
grid on;


## inverse filtering mode 2

freq = 915;  % estimated peak frequency in Hz
bw = 40;        % peak bandwidth estimate in Hz

R = exp( - pi * bw / fs);            % pole radius
z = R * exp(j * 2 * pi * freq / fs); % pole itself
B2 = [1, -(z + conj(z)), z * conj(z)] % inverse filter numerator
r = 0.9;     % zero/pole factor (notch isolation)
A2 = B2 .* (r .^ [0 : length(B2)-1]);   % inverse filter denominator

residual2 = filter(B2,A2,residual);       % apply inverse filter

Sfull2 = fft(residual2,Nfft);
S2 = Sfull2(iposFreq);
S2dB = 20*log10(abs(S2));
## offset = max(S1dB);
S2dBNorm = S2dB - offset;

figure;
set(gcf, 'Position', get(0, 'Screensize'));
plot(fkk./1000,[SdBNormi,S1dBNorm,S2dBNorm]);
axis([0 1.4 -60 0]);
grid on;


## plot filtered signals in time domain
figure;
set(gcf, 'Position', get(0, 'Screensize'));
itime = 1:length(residual)/8;
plot(t,[signali(itime), residual(itime), residual2(itime)]);
legend("warped original","mode1 removed","mode2 removed");











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



