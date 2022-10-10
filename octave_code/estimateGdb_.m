% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% estimateGdb.m
% estimate T60s and use to identify loop filter gains
% for harp string model
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;
clear all;

## functions

## usage: dbval = mydb (val)
##
##
function dbval = mydb (val)
  dbval = 20*log10(val);
endfunction

## ##############################################################
## f0 is the fundamental frequency of the string being identified
## ##############################################################
f0 = 370; % McLiag G4 is flat!!  Just G3 = 196Hz, G4 = 392Hz


## #################
## call edrAndPlot.m
## #################
run edrAndPlot;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% first estimate T60 - taken from estimateT60.m REAL SIMPLE LAB on Resonators
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimate the decay time (time it takes to decay 60 dB) at different 
% frequency bins by fitting a line to log magnitude data across time frames
%
% Variables are re-used from edrAndPlot.m, so this program should
% be run after calculating the EDR. You may need to need to modify
% the "dBfit" and "startTimeMS" values in order to obtain a good 
% region over which to fit the line.
%
% Music 421, Spring Quarter 2003-2004

Sdb = B_EDRdbN; % either the spectrogram or the EDR, in dB

slopes = zeros(nBins,1);
offsets = zeros(nBins,1);
RT = zeros(nBins,1); % array of t60s


%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% find reverberation times - for harmonics 1 to 10
%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for ii = 1:length(loc); % use harmonic indices in variables "loc", previously determined using findpeaks fn within edrAndPlot.m.

  %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% line fit using polyfit, adjusting dBfit for various freq bands...
  band1 = floor(length(loc)*0.2); %first 25%
  band2 = floor(length(loc)*0.4); %first 20%
  band3 = floor(length(loc)*0.66); %first 50% (above ca 10Hz, where the gain should start to drop dramatically)                
  
  if (ii <= band1) 
    dBfit = 30; % fit line to this amount of decay in dB
    startTimeMS = 100; % start line fitting after this offset in ms      
  elseif (ii > band1 && ii <= band2)   
    dBfit = 20; % inital decay steeper in this bandb
    startTimeMS = 50; % start line fitting after this offset in ms      
  elseif ((ii > band2) && (ii <= band3))
    dBfit = 10; % fit to  short steep decay in this band
    startTimeMS = 10; % start line fitting after this offset in ms      
  elseif (ii > band3)
    dBfit = 1; % fit to very short steep decay in this band
    startTimeMS = 10; % start line fitting after this offset in ms
  ## elseif (ii == 1) # || ii == 29)
  ##   dBfit = 50; % fit to very short steep decay in this band
  ##   startTimeMS = 6000; % start line fitting after this offset in ms  
  endif;


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% slopes / gdb values are super sensitive to startTimeMS
%% tuned by observing shape of EDR decays close to zero time 
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[minDiffMS minIdx] = min(abs(T*1000-startTimeMS));
startFrame(loc(ii)) = minIdx(1); % starting frame corresponds to the time offset

      
i = loc(ii);
idx = find(Sdb(i,:) < (Sdb(i,startFrame(i))-dBfit));
if (length(idx) == 0)
  endFrame(i) = nFrames; % last frame
else 
  endFrame(i) = idx(1); % first frame idx after dBfit (30 or so) dB of decay
end

segmentFrames = (startFrame(i):endFrame(i));
overlapSamp = endFrame(i) - startFrame; %my additional line!
segmentEDRdb = Sdb(i,segmentFrames); % EDR segment to fit line on

P = polyfit(T(segmentFrames), segmentEDRdb, 1);
slopes(i) = P(1);
offsets(i) = P(2);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% reveration time RT60, time to decay 60dB
  % https://ccrma.stanford.edu/realsimple/phys_mod_overview/Loop_Filter_Estimation.html
  % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  RT(i) = -60 / slopes(i);


  %% %%%%%%%%%
  % gdb
  % %%%%%%%%%

  gdb(i) = -60/(RT(i)*f0);

% above equiv to 
% g(i) = 10.^(-3./(RT(i)*f0));
% gdb(i) = 20*log10(g(i));

end


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot fitted slopes 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
doplot = 0;  % set to zero to suppress RT line fit plots
if (doplot)
  batchsize = 20;
  nobatches = floor(length(loc)/batchsize);
  res = mod(length(loc),batchsize*nobatches);
  batch = [0:batchsize:batchsize*nobatches (batchsize*nobatches + res)]; % group plots of slope fits into batches of ca. 20

  for i = 2:length(batch)
    figure;
    set(gcf, 'Position', get(0, 'Screensize'));
    batchvec = (batch(i-1)+1):batch(i);
    disp(["batchno = " num2str(i)]);
    for j = batchvec(1):batchvec(end)
      disp(["harmonic no. = " num2str(j)]);
      j_ = j - batch(i-1);
      disp(["plot no. = " num2str(j_)]);
      sp1_ix = loc(j);
      segmentFrames_1 = (startFrame(sp1_ix):endFrame(sp1_ix));
      subplot(4,5,j_);
      plot(T,Sdb(sp1_ix,:),'.-b');
      hold on;
      plot(T,(T*slopes(sp1_ix)+offsets(sp1_ix,:)),'r--');
      grid;
      plot(T(segmentFrames_1),Sdb(sp1_ix,segmentFrames_1),'g','linewidth',3);
            % legend(['EDC at ' num2str(F(sp1_ix)) ' Hz'],'line fit');
      xlabel('Time(s)');ylabel('magnitude (dB)');
% title(['harm=' num2str(ii(i)) '; freq=' num2str(F(loc(ii(i)))) '  t60 = ' num2str(RT(sp1_ix))]);
      title(['harm=' num2str(j) '; freq=' num2str(F(loc(j))) '  t60 = ' num2str(RT(sp1_ix))]);
      ylim([-200,0]);
    end
  end
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% prepare gdb to be in full "fft buffer format", ahead of converting to minimum phase
% The gain at half the sampling-rate was set arbitrarily to -4.08 dB per
% https://ccrma.stanford.edu/realsimple/phys_mod_overview/Loop_Filter_Estimation.html
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gdb_peaks = [gdb(loc(1)) gdb(loc) -4.0];
F_peaks = [F(1) F(loc) F(end)];

% figure;
% plot(F_peaks,gdb_peaks,'o-','color','r','linewidth',2);
% grid minor on;


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% smooth ahead of invfreq
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lambda = 0.000001;
## power of 2 plus 1 for the dc value, when neg freqs added later (ie not incl fs/2 or dc) total freq samples will be a power of 2
Nsmooth = 2^9 + 1; 
F_smooth = linspace(0,F(end),Nsmooth);
gdb_smooth = regdatasmooth(F_peaks, gdb_peaks, "lambda", lambda, "xhat", F_smooth');

## figure;
## set(gcf, 'Position', get(0, 'Screensize'));
## plot(F_peaks, gdb_peaks); hold on;
## plot(F_smooth, gdb_smooth); grid minor on;
## legend("gdb_{peaks}","gdb_{smooth}");


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% convert gdb into whole spectrum "fft buffer format",
% i.e., dc followed by positive freq vals, followed by neg freq vals, and power of 2 in length
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## gdb_whole = [gdb_smooth' flip(gdb_smooth')];
gdb_s = gdb_smooth';
negfreqi = Nsmooth-1:-1:2;
gdb_whole = [gdb_s gdb_s(negfreqi)];
## add on neg frequencies, excl dc and fs/2
F_whole = [F_smooth F_smooth(end) + F_smooth(2:end-1)];



% %%%%%%%%%%%%%%%%%%%%%%%%%%
% CONVERT TO MINIMUM PHASE
% see https://ccrma.stanford.edu/~jos/filters/Conversion_Minimum_Phase.html
% calls mps.m, downloaded from:
% https://ccrma.stanford.edu/~jos/filters/Matlab_listing_mps_m_test.html#sec:tmps
% note alternative code from mp also here:
% https://ccrma.stanford.edu/~jos/pasp/Converting_Desired_Amplitude_Response.html
% contribution of minimum-phase zeros to the complex cepstrum was described in:
% https://ccrma.stanford.edu/~jos/filters/Poles_Zeros_Cepstrum.html#sec:cepstrum
% %%%%%%%%%%%%%%%%%%%%%%%%%%

fk = F_smooth;
Nfft = 2^10;                     # 512
Ns = Nsmooth; 
gdb_s = gdb_smooth';
Sdb = [gdb_s,gdb_s(Ns-1:-1:2)]; % install negative-frequencies

Smpp = minphase(fk,Nfft,Ns,Sdb);

wt = 1 ./ (fk+1); % typical weight fn for audio
wk = 2*pi*fk/fs;
NZ = 30;
NP = 30;
[Bi,Ai] = invfreqz(Smpp,wk,NZ,NP,wt);
Hh = freqz(Bi,Ai,Ns);


figure;
set(gcf, 'Position', get(0, 'Screensize'));
plot(fk,mydb([Smpp(:),Hh(:)]));
hold on;
plot(F_peaks, gdb_peaks, "*m");
grid minor on;
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
title('Magnitude Frequency Response');
legend('Desired','Filter','Measured');

[H,w] = freqz(Bi,Ai);
figure;
set(gcf, 'Position', get(0, 'Screensize'));
subplot(211);
plot(fk,abs(Smpp)); hold on;
plot(w./pi.*fk(end),abs(H));
ylabel("Magnitude response (-)");
grid minor on;
title("Comparison of desired (minimum phase) frequency response, with filter frequency response");
legend("desired","filter","location","SouthWest");
subplot(212);
plot(fk,rad2deg(angle(Smpp))); hold on;
plot(w./pi.*fk(end),rad2deg(angle(H)));
ylabel("phase response (Hz)");
xlabel("frequency (Hz)");
grid minor on;
legend("desired","filter","location","SouthWest");

## ###############################################
## Ai and Bi in rows to copy over to c++ DWG model
## ###############################################
disp(length(Ai));       
for i = 1:length(Ai); printf([',' num2str(Ai(i))]);end;
disp('\n');
for i = 1:length(Bi); printf([',' num2str(Bi(i))]);end;


