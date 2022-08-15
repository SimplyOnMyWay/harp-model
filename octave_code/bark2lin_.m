function [out1, out2, out3, out4] = bark2lin_(in1, in2, in3, in4, in5)
% BARK2LIN - inverse Bark frequency warping via first-order conformal map.
%
% G = bark2lin(Fs) returns the parameter G such that the conformal
% map z --> (z - G) / (1 - z*G) closely aproximates an inverse Bark
% frequency warping for the input sample rate Fs.  Note that the
% sample rate Fs must lie in the interval (0Hz, 50kHz]. Also, 0 < G < 1.
%
% [W, G] = bark2lin(WB, Fs) warps input frequencies WB on Bark
% according to the conformal map generated by Fs to form linear
% frequencies W.  The input WB must lie in the range [-1.0, 1.0].
% Note that the corresponding output linear frequency W will always be
% less than the Bark input frequency WB for the inverse Bark mapping.
%
% [H, W, G] = bark2lin(HB, WB, Fs) warps the transfer function HB(WB)
% from Bark so that the values of HB at the input frequencies WB get mapped 
% to the new function H at linear frequencies W.
%
% [Z, P, K, G] = bark2lin(ZB, PB, KB, Fs) maps an input transfer
% function in zero-pole-gain -- ZB, PB, KB -- format from Bark to linear.
%
% See also LIN2BARK.

% Notes:
% (1) The map z --> (z - G)/(1 - z*G) provides the best fit to Bark
% at lower sample rates, though good fits are obtained through 50kHz.
%
% (2) Note that the output warped transfer function H is computed
% via H = spline(W, HB, WB).  If the input frequency scale WB is not
% sufficiently dense, H may contain spurious sample values.

% References:
% [1] Julius O. Smith III, Techniques for Digital Filter Design and System
% Identification with Application to the Violin, Stanford University Ph.D.
% Thesis, June 1983, pp. 133--137.
%
% [2] E. Zwicker, Psychoacoustics, Springer-Verlag: New York, 1990, page 142.

% Copyright (c) 1994 Abel Innovations.  All rights reserved.
%
% Jonathan S. Abel
% Created: August 28, 1994
% Revised: September 9, 1994    - pole-zero-gain and H(w) formats added
% Revised: October 23, 1998 by Yoon Kim
%          - Corrected formula for the weight v(\omega_k) in calculating rho
%          - Revised the warping of the frequency from WB to W so that 
%            WB is the Bark frequency that is always larger than the linear
%            frequency W in the inverse Bark mapping. 
%          - Corrected the formula for calculating zeros and poles ZB, PB. 
%          - Revised the outputs so that rho is not negated.
%          - Revised comments regarding the definition of the warped transfer
%            function HB so that it is mathematically feasible.
% Version: 1.0
%


%% assign and verify inputs
%%
if (nargin < 1),
    disp(['lin2bark: ERROR - not enough input arguments.']);
    return;
elseif (nargin == 1),
    fs = in1(1,1);
elseif (nargin == 2),
    w = in1(:);
    fs = in2(1,1);
elseif (nargin == 3),
    H = in1(:);
    w = in2(:);
    fs = in3(1,1);
elseif (nargin == 5),
    z = in1(:);
    p = in2(:);
    k = in3(1,1);
    fs = in4(1,1);
    w = in5(:);
else,
    disp(['lin2bark: ERROR - too many input arguments.']);
    return;
end;

if (fs ~= max([eps min([50000 fs])])),
    disp(['lin2bark: ERROR - input sample rate out of range.']),
    return;
end;
if (w ~= max(-1.0, min(1.0, w))),
    disp(['lin2bark: ERROR - some input frequencies fall outside [-1.0, 1.0]']),
    return;
end;

%% determine conformal map parameter RHO
%%
linear = [0:24]'/24;
lambda = pi * linear;

bark = [0 100 200 300 400 510 630 770 920 1080 1270 1480 1720 2000 2320, ...
        2700 3150 3700 4400 5300 6400 7700 9500 12000 15500 20500 27000]';
bs = spline(bark, [0:length(bark)-1], fs/2);
beta = pi * spline([0:length(bark)-1], bark, bs * linear) / (fs/2);

rhoU = sum(cos(beta) - cos(lambda)) / sum(1 - cos(beta+lambda));
weight = ones(size(lambda)) ./ (1 + rhoU^2 - 2*rhoU*cos(lambda));
rho =  (cos(beta) - cos(lambda))'*weight / ((1 - cos(beta+lambda))'*weight);

%% compute HB and/or WB if needed
%%
if (nargin == 2 | nargin == 3),
    cmap = (exp(sqrt(-1)*lambda)+ rho) ./ (1 + rho*exp(sqrt(-1)*lambda));
    wb = spline(linear, angle(cmap)/pi, abs(w)) .* sign(w);
    if (nargin == 3),
        Hb = spline(wb, H, w);
    end;
end;

%% compute ZB, PB, and KB if needed
%%
if (nargin == 5),
    zb = (z + rho) ./ (1 + rho*z);
    pb = (p + rho) ./ (1 + rho*p);
    npz = length(p)-length(z);
    if (npz > 0),
        zb = [zb; ones(npz,1)*rho];
    elseif (npz < 0),
        pb = [pb; ones(-npz,1)*rho];
    end;
    kb = k * prod(1 + z*rho) / prod(1 + p*rho);
end;


%% return results
%%
if (nargin == 1),
    out1 = rho;
elseif (nargin == 2),
    out1 = wb;
    out2 = rho;
elseif (nargin == 3),
    out1 = Hb;
    out2 = wb;
    out3 = rho;
elseif (nargin == 5),
    out1 = zb;
    out2 = pb;
    out3 = kb;
    out4 = rho;
end;
