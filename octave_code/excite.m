## ################################################################
## excite.m
## synthesise excitation of my harp string model
## based on https://ccrma.stanford.edu/realsimple/vguitar/Computing_Excitation_Signal.html
## thanks and respect Julius!!!
## ################################################################


## close all; clear all;
## run estimateGdb.m;

## find the excitation signal. 
## FUND_F corresponds to the fundamental freq
N = round(fs/f0);
y = signal;

del_y = [zeros(N,1);y(1:end-N)];
filt_y = filter(Ai,Bi,del_y);

e_sig = y-filt_y;
audiowrite("e_sig.wav",e_sig,fs);


## plots
T = 1/fs;
t = T:T:T*length(y);
ydb = 20*log10(abs(y(2:end)));
offset = max(ydb);
ydbNorm = ydb - offset; #max(ydb);
e_sigdb = 20*log10(abs(e_sig(2:end)));
e_sigdbNorm = e_sigdb - offset; #max(e_sigdb);


figure;
set(gcf, 'Position', get(0, 'Screensize'));
subplot(311);
plot(t,[y,e_sig]);
grid minor on;
xlabel("t [s]");
ylabel("tone amplitude");
legend("recorded","extracted excitation");
subplot(312);
plot(t(2:end),[ydbNorm,e_sigdbNorm]);
grid minor on;
xlabel("t [s]");
ylabel("tone amplitude [dB]");
legend("recorded","extracted excitation");
subplot(313);
plot(t,[y,e_sig]);
grid minor on;
xlabel("t [s]");
ylabel("tone amplitude");
legend("recorded","extracted excitation");
axis([0 1 -1 1]);


## play my new excitation!
iend = find(abs(t-2.0)<T/2);
sound(e_sig(1:iend),fs);
sound(y(1:iend),fs);


## fit adsr with LP noise to e_sig
## per https://www.dsprelated.com/freebooks/pasp/Approximating_Shortened_Excitations_Noise.html

