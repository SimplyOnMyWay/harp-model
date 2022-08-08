% edrAndPlot.m  
% calculate and plot the Energy Decay Relief of a signal
%
% Music 421, Spring Quarter 2003-2004

% function B_EDRdb = edrAndPlot(s,mode)
% 

fileName = "mcliag_G4_3.wav";
minPlotDB = -200; % minimum value in DB to be plotted
frameSizeMS = 10; % minimum frame length, in ms
overlap = 0.5; % fraction of frame overlapping
windowType = 'hann'; % type of windowing used for each frame

mode = 1;

[signal,fs] = audioread(fileName);
if (mode == 1)
  signal = signal(:,1);
end

%% calculate STFT frames
minFrameLen = fs*frameSizeMS/1000; 
frameLenPow = nextpow2(minFrameLen);
frameLen = 2^frameLenPow; % frame length = fft size
% eval(['frameWindow = ' windowType '(frameLen);']);
[B,F,T] = specgram(signal,frameLen,fs);

Bdb = 20*log10(abs(B));





%%% WRITE CODE HERE %%%
     % Complete the rest of the code in order to compute the signal
     % spectrogram. Use the Matlab function "specgram", and
     % save the result into vectors [B,F,T] (this will be used
     % later for plotting the time and frequency axes).

%B_EDRdb = zeros(size(B));

% calculate the EDR
     %%% WRITE CODE HERE %%%
     % The variable "B_EDRdb" is the EDR in units of decibels.
     % Keep the matrix orientation of the EDR the same as the spectrogram. 
     % Implementation hint: It is faster to calculate the EDR going
     % backward in time rather than forward.

[nBins, nFrames] = size(B);

B_energy = B.*conj(B);
B_EDR = zeros(nBins, nFrames);
for i = 1:nBins
          B_EDR(i,:) = fliplr(cumsum(fliplr(B_energy(i,:))));
end

## convert to dB
B_EDRdb = 20*log10(abs(B_EDR));

% normalize EDR to 0 dB and truncate the plot below a given dB threshold
offset = max(max(B_EDRdb));
B_EDRdbN = B_EDRdb - offset;
[nBins,nFrames] = size(B_EDRdbN);
for i=1:nFrames
  I = find(B_EDRdbN(:,i) < minPlotDB);
  if (I)
    B_EDRdbN(I,i) = minPlotDB;
  end
end

% plot the energy decay relief
figure;
set(gcf, 'Position', get(0, 'Screensize'));
mesh(T,F/1000,B_EDRdbN);
view(130,30);
title('Normalized Energy Decay Relief (EDR)');
xlabel('Time (s)');ylabel('Frequency (kHz)');zlabel('Magnitude (dB)');
axis tight;zoom on;

#{
% plot the spectrum 
figure(gcf);clf;
mesh(T,F/1000,Bdb);
view(130,30);
title('Spectrum');
xlabel('Time (s)');ylabel('Frequency (kHz)');zlabel('Magnitude (dB)');
axis tight;zoom on;
#}


% %%%%%%%%%%%%%%%%
% find harmonic peaks
% %%%%%%%%%%%%%%%%

% calculate the min number of frequency samples between peaks - assumes peaks are not much less than f0 apart
minPeakDistance = floor(f0/(F(2)-F(1)));
% add 100 to B_EDRdbN so no neg values - not exact, so increase if error thrown "Data contains negative values..."
[pks,loc,extra]=findpeaks(2*offset+B_EDRdbN(:,1),"MinPeakDistance",minPeakDistance);

figure;
set(gcf, 'Position', get(0, 'Screensize'));
plot(F,B_EDRdbN(:,1)); hold on;
plot(F(loc),B_EDRdbN(loc,1),"r*"); grid minor on;
xlabel("Frequency (Hz)");
ylabel("B_{EDRdbN}(:,1) (dB)"); 
title("harmonics identified using function findpeaks (ignore first airmode harmonic at 50Hz)");



