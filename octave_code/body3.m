clear all; close all;

##################
## local functions
## ###############

## usage: y = lpf (n,wc,x)
##
##
function y = lpf (n,wc,x)
  [b,a] = butter(n, wc); % low pass Butterworth filter with cutoff pi*Wc radians - choose the order of the filter n and cut-off frequency Wc to suit
  y = filter(b,a,x);
endfunction


## usage: phase = ph (S)
##
##
function phase = ph (S)
  phase = rad2deg(unwrap(angle(S)));
endfunction

## usage: a = adsr (target,gain,duration,fs)
## inputs:
## target - vector of attack, sustain, release target values
## gain - vector of attack, sustain, release gain values
## duration - vector of attack, sustain, release durations in ms
## outpus:
## a - vector of adsr envelope values
## from https://web.nmsu.edu/~pdeleon/Research/Publications/ASEE_GSW_2000.pdf
function a = adsr (target,gain,duration,fs,envdur)
  a = zeros(round(fs*envdur),1); 
  duration = round(duration./1000.*fs); % envelope duration in samp
  ## Attack phase
  start = 2;
  stop = duration(1);
  for n = [start:stop]
    a(n) = target(1)*gain(1) + (1.0 - gain(1))*a(n-1);
  end

  ## Sustain phase
  start = stop + 1;
  stop = start + duration(2);
  for n = [start:stop]
    a(n) = target(2)*gain(2) + (1.0 - gain(2))*a(n-1);
  end
  ## Release phase
  start = stop + 1;
  stop = start + duration(3);
  for n = [start:stop]
    a(n) = target(3)*gain(3) + (1.0 - gain(3))*a(n-1);
  end;
endfunction


## ###############################################
## read in body impulse response and set up its FR
## ###############################################

filename = "gtrbody.wav";
[signal,fs] = audioread(filename);
sound(signal,fs);
# time axes...
T = 1/fs;
t = T:T:T*length(signal);
istart = find(abs(t-0.001)<T/2);
iend = find(abs(t-0.300)<T/2);

Nfft = 2^16;  #I think this is oversampling by factor of 2^3 = 8, but it seems to need it as FR looks jaggedy below 2^16
Sfull = fft(signal(istart:iend),Nfft);
iposFreq = 1:Nfft/2+1;
S = Sfull(iposFreq);
SdB = 20*log10(abs(S));
offset = max(SdB);
SdBN = SdB - offset;
fk = linspace(0,fs,Nfft);
fkk = fk(iposFreq);


## identify central freq and bandwidth from observing the full FR
fc1 = 104.7;
bw1 = 3;
fc2 = 205.4;
bw2 = 10;
fc3 = 247.5;
bw3 = 5;


## ############################
## apply the inverse filters...
## ############################

## mode 1
r = 0.98;
[B1,A1] = invfilter(fc1,bw1,r,fs);
H1 = freqz(B1,A1,fkk,fs);
H1dB = 20*log10(H1);
H1dBN = H1dB - offset;
H1inv = freqz(A1,B1,fkk,fs);
H1invdB = 20*log10(H1inv);
H1invdBN = H1invdB - offset;

## mode 2
r = 0.9;
[B2,A2] = invfilter(fc2,bw2,r,fs);
H2 = freqz(B2,A2,fkk,fs);
H2dB = 20*log10(H2);
H2dBN = H2dB - offset;
H2inv = freqz(A2,B2,fkk,fs);
H2invdB = 20*log10(H2inv);
H2invdBN = H2invdB - offset;

## mode 3
r = 0.98;
[B3,A3] = invfilter(fc3,bw3,r,fs);
H3 = freqz(B3,A3,fkk,fs);
H3dB = 20*log10(H3);
H3dBN = H3dB - offset;
H3inv = freqz(A3,B3,fkk,fs);
H3invdB = 20*log10(H3inv);
H3invdBN = H3invdB - offset;

## plot biquad fits and their inverses, against original SdB
figure;
set(gcf, 'Position', get(0, 'Screensize'));
plot(fkk,[SdBN,H1dBN',H1invdBN',H2dBN',H2invdBN',H3dBN',H3invdBN']);
grid minor on;
axis([0 500 -80 0]);


## #####################################
## apply inverse filters to get residues
## #####################################

res1 = filter(A1,B1,signal);
res2 = filter(A2,B2,res1);
res3 = filter(A3,B3,res2);

figure;
set(gcf, 'Position', get(0, 'Screensize'));
plot(t.*1000,[signal, res1, res2, res3]);
axis([0, 80, -1, 1]);
grid minor on;
legend("signal","res1","res2","res3");


## ######################################################
## check FR after applying each biquad as inverse of mode
## ######################################################

## FR after applying biquad 1
S1full = fft(res1,Nfft);
S1 = S1full(iposFreq);
S1dB = 20*log10(abs(S1));
S1dBN = S1dB - offset;

## FR after applying biquad 2
S2full = fft(res2,Nfft);
S2 = S2full(iposFreq);
S2dB = 20*log10(abs(S2));
S2dBN = S2dB - offset;

## FR after applying biquad 3
S3full = fft(res3,Nfft);
S3 = S3full(iposFreq);
S3dB = 20*log10(abs(S3));
S3dBN = S3dB - offset;

## plot FR of original signal, and after applying each inverse filter
figure;
set(gcf, 'Position', get(0, 'Screensize'));
plot(fkk,[SdBN,S1dBN,S2dBN,S3dBN]);
axis([0 500 -60 0]);
grid minor on;
legend("SdBN","S1dBN","S2dBN","S3dBN");


## #######################################
## fit adsr to res3 in time and freqdomain
## #######################################

target = [0.999;0.05;0.00];
gain = [0.05;0.006;0.003];
duration = [1;20;50];
env = adsr(target,gain,duration,fs,t(end));

figure;
set(gcf, 'Position', get(0, 'Screensize'));
plot(t.*1000,env);
axis([0 80 -1 1]);
grid minor on;

S3mp = mps(S3full);
S3mpp = S3mp(iposFreq);
wk = iposFreq.*pi./iposFreq(end);
wkk = wk(1:2^7:end);
wt = 1./(wkk+1);
S3mppi = interp1(wk,S3mpp,wkk);
S3mppidBN = 20*log10(S3mppi) - offset;
iforder = 10;
[B3if,A3if] = invfreqz(S3mppi,wkk,iforder,iforder,wt);
S3if = freqz(B3if,A3if,fkk,fs);
S3ifdB = 20*log10(S3if);
S3ifdBN = S3ifdB - offset;

figure;
set(gcf, 'Position', get(0, 'Screensize'));
plot(wkk./pi.*fkk(end),S3mppidBN,"o-");hold on;
plot(fkk,[S3dBN, S3ifdBN']);
grid minor on;
legend("S3mppidBN","S3dBN","S3ifdBN");

## generate synthesised noise burst, with IR shaped by adsr
rand("seed",1);
x = (rand(round(fs*t(end)),1) - 0.5).*2;
# xlp = lpf(5,0.9,x);
y =  env .* x; % Modulate

## shape the FR of the noise burst using the LPF designed using invfreqz,ie [B3if,A3if]
## res3synth = lpf(2,0.9,y);
res3synth_ = filter(B3if,A3if,y);
res3synth = res3synth_./max(res3synth_);
## y = y./max(y);
sound(res3synth,fs);
audiowrite("adsrimpulse.wav",res3synth,fs);

figure;
set(gcf, 'Position', get(0, 'Screensize'));
plot(t.*1000,[res3, res3synth]);
axis([0, 80, -1, 1]);
grid minor on;
legend("res3","res3synth");

## FR after applying biquad 3
S3synthfull = fft(res3synth,Nfft);
S3synth = S3synthfull(iposFreq);
S3synthdB = 20*log10(abs(S3synth));
S3synthdBN = S3synthdB - offset;

figure;
set(gcf, 'Position', get(0, 'Screensize'));
plot(fkk,[S3dBN, S3synthdBN]);
#axis([0, 80, -1, 1]);
grid minor on;
legend("S3dBN","S3synthdBN");


## #####################################
## reconstitute res2, res1 and signal...
## #####################################

res2synth = filter(B3,A3,res3synth);
res1synth = filter(B2,A2,res2synth);
signalsynth = filter(B1,A1,res1synth);

figure;
set(gcf, 'Position', get(0, 'Screensize'));
## plot(t.*1000,[res2,res2synth]);
plot(t.*1000,[signalsynth,signal, res1synth, res1, res2synth, res2]);
axis([0 80 -1 1]);
grid minor on;
legend("signalsynth","signal", "res1synth", "res1", "res2synth", "res2");



## #########
## last bits
## #########

## enjoy the filtered sounds!
sound(signal,fs);
sound(res1,fs);
sound(res2,fs);
sound(res3,fs);

## save to WAV files
audiowrite("res1.wav",res1,fs);
audiowrite("res2.wav",res2,fs);
audiowrite("res3.wav",res3,fs);




