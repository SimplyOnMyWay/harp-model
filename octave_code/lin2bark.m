function [out1, out2, out3, out4] = lin2bark(in1, in2, in3, in4)
% LIN2BARK - Bark frequency warping via first-order conformal map.
%
% G = lin2bark(Fs) returns the parameter G such that the conformal
% map z --> (z + G) / (z*G + 1) closely aproximates a Bark frequency
% warping for the input sample rate Fs.  Note that the sample rate
% Fs must lie in the interval (0Hz, 50kHz]. Also, 0 < G < 1.
%
% [WB, G] = lin2bark(W, Fs) warps input frequencies W according to
% the conformal map generated by Fs to form Bark frequencies WB.  The
% input W must lie in the range [-1.0, 1.0]. Note that the corresponding
% Bark frequency will always be greater than the input frequency W for
% our forward Bark mapping (0 < G < 1)
%
% [HB, WB, G] = lin2bark(H, W, Fs) warps the transfer function H(W)
% onto Bark so that the values of H at the input frequencies W get mapped to 
% the new function HB at Bark frequencies W'=WB.
%
% [ZB, PB, KB, G] = lin2bark(Z, P, K, Fs) maps an input transfer
% function in zero-pole-gain -- Z, P, K -- format to Bark.
%
% See also BARK2LIN.

% Notes:
% (1) The map z --> (z + G)/(z*G + 1) provides the best fit to Bark
% at lower sample rates, though good fits are obtained through 50kHz.
%
% (2) Note that the output warped transfer function HB is computed
% via HB = spline(WB, H, W).  If the input frequency scale W is not
% sufficiently dense, HB may contain spurious sample values.

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
% Revised: September 9, 1994            - pole-zero-gain and H(w) formats added
% Revised: October 23, 1998 by Yoon Kim
%          - Corrected formula for the weight v(\omega_k) in calculating rho
%          - Reversed the warping of the frequency from W to WB so that 
%            WB is the Bark frequency that is always larger than the linear
%            frequency W in the forward Bark mapping. 
%          - Corrected the formula for calculating zeros and poles ZB, PB. 
%          - Revised comments regarding the definition of the warped transfer
%            function HB so that it is mathematically feasible.

% Version: 1.0




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
elseif (nargin == 4),
    z = in1(:);
    p = in2(:);
    k = in3(1,1);
    fs = in4(1,1);
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
rho = (cos(beta) - cos(lambda))'*weight / ((1 - cos(beta+lambda))'*weight);

%% compute HB and/or WB if needed
%%
if (nargin == 2 || nargin == 3),
    cmap = (exp(sqrt(-1)*lambda) + rho) ./ (1 + rho*exp(sqrt(-1)*lambda));
    wb = spline(angle(cmap)/pi, linear, abs(w)) .* sign(w);
    if (nargin == 3),
        Hb = spline(wb, H, w);
    end;
end;

%% compute ZB, PB, and KB if needed
%%
if (nargin == 4),
    zb = (z - rho) ./ (1 - rho*z);
    pb = (p - rho) ./ (1 - rho*p);
    npz = length(p)-length(z);
    if (npz > 0),
        zb = [zb; ones(npz,1)*(-rho)];
    elseif (npz < 0),
        pb = [pb; ones(-npz,1)*(-rho)];
    end;
    kb = k * prod(1 - z*rho) / prod(1 - p*rho);
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
elseif (nargin == 4),
    out1 = zb;
    out2 = pb;
    out3 = kb;
    out4 = rho;
end;
