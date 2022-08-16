## usage: Smp = minphase (fk,Nfft,Ns,Sdb)
## alternative mp code!
## see
## fk: frequency values up to fs/2
## Nfft: number of freq blocks in fft
## Ns: number of samples in "measured" and smoothed FR to be minimum phase'd
## Sdb: the FR to be min-phase'd, in decibels
function Smpp = minphase (fk,Nfft,Ns,Sdb)
  #ensure Sdb is a row vector
  [mm,nn] = size(Sdb);
  if mm>nn
  Sdb = Sdb'; # ensure row vector
  endif
  S = 10 .^ (Sdb/20); % convert to linear magnitude
  if Ns~=Nfft/2+1, error("confusion"); end
  s = ifft(S); % desired impulse response
  s = real(s); % any imaginary part is quantization noise
  tlerr = 100*norm(s(round(0.9*Ns:1.1*Ns)))/norm(s);
  disp(sprintf(['Time-limitedness check: Outer 20%% of impulse ' ...
                  'response is %0.2f %% of total rms'],tlerr));
                       % = 0.02 percent
  if tlerr>1.0 % arbitrarily set 1% as the upper limit allowed
    error('Increase Nfft and/or smooth Sdb');
  end

  c = ifft(Sdb); % compute real cepstrum from log magnitude spectrum
        % Check aliasing of cepstrum (in theory there is always some):
  caliaserr = 100*norm(c(round(Ns*0.9:Ns*1.1)))/norm(c);
  disp(sprintf(['Cepstral time-aliasing check: Outer 20%% of ' ...
                  'cepstrum holds %0.2f %% of total rms'],caliaserr));
                                % = 0.09 percent
  if caliaserr>1.0 % arbitrary limit
    error('Increase Nfft and/or smooth Sdb to shorten cepstrum');
  end
    % Fold cepstrum to reflect non-min-phase zeros inside unit circle:
    % If complex:
% cf = [c(1), c(2:Ns-1)+conj(c(Nfft:-1:Ns+1)), c(Ns), zeros(1,Nfft-Ns)];
  cf = [c(1), c(2:Ns-1)+c(Nfft:-1:Ns+1), c(Ns), zeros(1,Nfft-Ns)];
  Cf = fft(cf); % = dB_magnitude + j * minimum_phase
  Smp = 10 .^ (Cf/20); % minimum-phase spectrum

  Smpp = Smp(1:Ns)'; % nonnegative-frequency portion
  
endfunction
