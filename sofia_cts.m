% /// ASAR-MARA Research Group
%
% Cologne University of Applied Sciences
% Berlin University of Technology
% University of Rostock
% Deutsche Telekom Laboratories
% WDR Westdeutscher Rundfunk
% IOSONO GmbH
%
% SOFiA sound field analysis
%
% C/T/S Center Transducer Simulation R13-0306
%
% Produces the corresponding center transducer signals for the wave
% generators. The center transducer signal is e.g. required for the
% SOFiA B/S/A BEMA Spatial Anti Aliasing module.
%
% Copyright (C)2012-2013 Benjamin Bernschütz
%                        rockzentrale 'AT' me.com
%                        +49 171 4176069 Germany
%
% This file is part of the SOFiA toolbox under GNU General Public License
%
%
% center = sofia_cts(FS, NFFT, AZ, EL, t, c)
% ------------------------------------------------------------------------
% center  Complex sound pressures (Center Transducer)  [1 x NFFT]
% ------------------------------------------------------------------------
% FS      Sampling rate [1/s]                       [default = 48000]
% NFFT    Number of FFT bins                        [default = 512]
% AZ      Azimuth angle in [RAD] 0-2pi              [default = 0]
% EL      Elevation angle in [RAD] 0-pi             [default = pi/2]
% t       Timeshift, for (FS/t) Samples             [default = 0]
% c       Speed of Sound [m/s]                      [default = 343]
%


% CONTACT AND LICENSE INFORMATION:
%
% /// ASAR-MARA Research Group
%
%     [1] Cologne University of Applied Sciences
%     [2] Berlin University of Technology
%     [3] Deutsche Telekom Laboratories
%     [4] WDR Westdeutscher Rundfunk
%     [5] University of Rostock
%     [6] IOSONO GmbH
%
% SOFiA sound field analysis
%
% Copyright (C)2011-2013 Benjamin Bernschütz [1,2] et al.(§)
%
% Contact -------------------------------------
% Cologne University of Applied Sciences
% Institute of Communication Systems
% Betzdorfer Street 2
% D-50679 Germany (Europe)
%
% phone +49 221 8275 -2496
% mail  benjamin.bernschuetz@fh-koeln.de
% ---------------------------------------------
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
% (§) Christoph Pörschmann [1]   christoph.poerschmann 'at' fh-koeln.de
%     Stefan Weinzierl     [2]   stefan.weinzierl 'at' tu-berlin.de
%     Sascha Spors         [5]   sascha.spors 'at' uni-rostock.de


function [center] = sofia_cts(FS, NFFT, AZ, EL, t, c)

if nargin < 6
    c=343;
end

if nargin < 5
    t=0;
end

if nargin < 4
    EL = pi/2;
end

if nargin < 3
    AZ = 0;
end

if nargin < 2
    NFFT = 512;
end

if nargin < 1
    FS = 48000;
end


if (FS*t >= NFFT)
    warning('Timeshift t > NFFT ! (Cyclic convolution)');
    disp(' ');
elseif (FS*t >= NFFT/2)
    warning('Timeshift t > NFFT/2 - Remember to guard headroom for filter responses.');
    disp(' ');
end

disp('SOFiA C/T/S - Center Transducer Simulation R13-0306');

EL = pi-EL;

f  = linspace(0,FS/2,NFFT/2+1);
k  = 2*pi*f/c;
w  = 2*pi*f;

fftBins = size(k,2);

center = zeros(1,fftBins);

[Xo,Yo,Zo]=sph2cart(0,0,0);
for bin=1:fftBins
    center(bin)=exp(1i*(k(bin)*Xo*sin(EL)*cos(AZ)+k(bin)...
        *Yo*sin(AZ)*sin(EL)+k(bin)*Zo*cos(EL)))*exp(-1i*w(bin)*t);
end


