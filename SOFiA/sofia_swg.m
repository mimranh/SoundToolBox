% /// ASAR Research Group
%
% Cologne University of Applied Sciences
% Technical University of Berlin
% Deutsche Telekom Laboratories
% WDR Westdeutscher Rundfunk
%
% SOFiA sound field analysis
%
% S/W/G Sampled Wave Generator Wrapper R11-1220
%
% Copyright (C)2011 by bBrn - benjamin Bernschütz
%                             rockzentrale 'AT' me.com
%                             +49 171 4176069 Germany
%
% This file is part of the SOFiA toolbox under GNU General Public License
%
%
% [fftdata, kr] = sofia_swg(r, gridData, ac, FS, NFFT, ...
%                                       AZ, EL, Nlim, t, c, wavetype, ds)
% ------------------------------------------------------------------------
% fftData  Complex sound pressures                   [(N+1)^2 x NFFT]
% kr       kr-Vector
%          Can also be a matrix [krm; krs] for rigid sphere configurations:
%           [1,:] => krm referring to the Microphone Radius
%           [2,:] => krs referring to the Sphere/Microphone2 Radius
% ------------------------------------------------------------------------
% r        Microphone Radius
%          Can also be a vector for rigid sphere configurations:
%           [1,1] => rm  Microphone Radius
%           [2,1] => rs  Sphere Radius (Scatterer)
% gridData Quadrature grid                           [default LEB110]
%           Columns : Position Number 1...M
%           Rows    : [AZ EL Weight]
%           Angles AZ, EL in [RAD]
% ac       Array Configuration
%           0  Open Sphere with p Transducers
%           1  Open Sphere with pGrad Transducers
%           2  Rigid Sphere with p Transducers
%           3  Rigid Sphere with pGrad Transducers (Thx to Nils Peters!)
%           4  Dual Open Sphere with p Transducers (Thx to Nils Peters!)
% FS       Sampling Frequency
% NFFT     FFT Order (Number of bins) should be 2^x, x=1,2,3,...
% AZ       Azimuth   angle in [DEG] 0-2pi
% EL       Elevation angle in [DEG] 0-pi
% Nlim     Internal generator transform order limit
% c        Speed of sound in [m/s] (Default: 343m/s)
% t        Time Delay in s
% wavetype Type of the Wave. 0: Plane Wave (default) 1: Spherical Wave
% ds       Distance of the source in [m] (For wavetype = 1 only)
%
% This file is a wrapper generating the complex pressures at the
% positions given in 'gridData' for a full spectrum 0-FS/2 Hz (NFFT Bins)
% wave impinging to an array. The wrapper involves the W/G/C wave
% generator core and the I/T/C spatial transform core.
%
% S/W/G emulates discrete sampling. You can observe alias artifacts.
%


% CONTACT AND LICENSE INFORMATION:
%
% /// ASAR Research Group
%
%     [1] Cologne University of Applied Sciences
%     [2] Technical University of Berlin
%     [3] Deutsche Telekom Laboratories
%     [4] WDR Westdeutscher Rundfunk
%
% SOFiA sound field analysis
%
% Copyright (C)2011 bBrn - benjamin Bernschütz [1,2] et al.(§)
%
% Contact ------------------------------------
% Cologne University of Applied Sciences
% Institute of Communication Systems
% Betzdorfer Street 2
% D-50679 Germany (Europe)
%
% phone       +49 221 8275 -2496
% cell phone  +49 171 4176069
% mail        rockzentrale 'at' me.com
% --------------------------------------------
%
% This file is part of the SOFiA sound field analysis toolbox
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program. If not, see <http://www.gnu.org/licenses/>.
%
%
% (§) Christoph Pörschmann [1]   christoph.poerschmann 'at' fh-koeln.de
%     Sascha Spors         [2,3] sascha.spors 'at' telekom.de
%     Stefan Weinzierl     [2]   stefan.weinzierl 'at' tu-berlin.de


function [fftData, kr] = sofia_swg(r, gridData, ac, FS, NFFT, AZ, EL, Nlim, t, c, wavetype, ds)

disp('SOFiA S/W/G Sampled Wave Generator Wrapper R11-1220');

if nargin < 12
    ds = 1;       %Default 1m
end

if nargin < 11
    wavetype = 0; %Default 0: Plane Wave
end

if nargin < 10
    c=343;        %Default 343m/s
end

if nargin < 9
    t=0;          %Default no timeshift
end

if nargin < 8
    Nlim = 120;   %Internal modal resolution limit of the wavefield generator
end

if nargin < 7
    EL = pi/2;    %Default 90°
end

if nargin < 6
    AZ = 0;       %Default 0°
end

if nargin < 5
    NFFT = 512;   %Default 512 Blocks
end

if nargin < 4
    FS = 48000;   %Default 48KHz
end

if nargin < 3
    ac = 0;       %Default Open Sphere
end

if nargin < 2
    gridData = sofia_lebedev(110, 0);  %EXAMPLE GRID LEB110, No Plot
end

if nargin == 0
    r = 0.1;
end

pause(1e-12);

minOrderLim                     = 70; %Default minimum
kr                              = linspace(0,r*pi*FS/c,NFFT/2+1);
rqOrders                        = ceil(kr*2);
maxReqOrder                     = max(rqOrders);
rqOrders(rqOrders<=minOrderLim) = minOrderLim;
rqOrders(rqOrders>Nlim)         = Nlim;
Ng                              = max(rqOrders);
Pnm                             = zeros((Ng+1)^2,(NFFT/2+1));


index  = 1;
ctr    = -1;

if (maxReqOrder > Ng)
    disp(' ');
    disp(['   -> Warning: The requested wave needs a generator order of ', num2str(maxReqOrder),' to']);
    disp(['      deliver a fully valid result for the current configuration.']);
    disp(['      The generator is actually limited to a maximum order of ',num2str(Ng),'.']);
    disp(' ');
else
    if Ng == minOrderLim
        disp(['      Full spectrum generator order: ', num2str(Ng)]);
    else
        disp(['      Segmented generator orders: ', num2str(minOrderLim),' to ',num2str(Ng)]);
    end
end

while(1)  %Segmentation
    
    if mod(ctr,36) == 0 && ctr>0
        fprintf('\n      ');
    end
    fOrders = find(rqOrders == rqOrders(index));
    try
        Pnm = Pnm+sofia_wgc(Ng, r, ac, FS, NFFT, AZ, EL, t, c, wavetype, ds, min(fOrders), max(fOrders), rqOrders(index));
    catch ME
        if strcmp(ME.identifier,'MATLAB:unexpectedCPPexception')
            error('SOFiA:WGC:WaveUnstable','SOFiA W/G/C: Result is unstable. Change simulation parameters (r, FS, Nlim).');
        else
            error(ME.identifier,ME.message);
        end
        fftData = zeros(size(gridData,1),(NFFT/2+1));
        return
    end
    index = index+size(fOrders,2);
    if index > size(rqOrders,2)
        break;
    end
    ctr=ctr+1;
    pause(1e-12)
end

pause(1e-12);

% Inverse spatial transform at the discrete sampling positions:

try
    fftData = sofia_itc(Pnm, gridData);
catch ME
    if strcmp(ME.identifier,'MATLAB:unexpectedCPPexception')
        error('SOFiA:WGC:WaveUnstable','SOFiA W/G/C: Result is unstable. Change simulation parameters (r, FS, Nlim).');
    else
        error(ME.identifier,ME.message);
    end
    fftData = zeros(size(gridData,1),(NFFT/2+1));
    return
end

