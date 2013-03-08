% /// ASAR Research Group
% 
%     Cologne University of Applied Sciences
%     Technical University of Berlin
%     Deutsche Telekom Laboratories
%     WDR Westdeutscher Rundfunk
% 
% SOFiA sound field analysis 
% 
% F/D/T frequency domain transform R11-1220
% 
% Copyright (C)2011 by bBrn - benjamin Bernsch�tz   
% 
% This file is part of the SOFiA toolbox under GNU General Public License
% 
%  
% [fftData, kr, f] = sofia_fdt(timeData, FFToversize, 
%                                   firstSample, lastSample)
% ------------------------------------------------------------------------
% fftData           Frequency domain data ready for the Spatial 
%                   Fourier Transform (SOFiA S/T/C Core)
%                   ! To save processing power the SOFiA chain always
%                     works on the half-sided FFT spectrum only (NFFT/2+1)
% kr                kr-Values of the delivered data (required e.g. for the 
%                   modal radial filters SOFiA M/F)
% f                 Absolute frequency scale. Not really required but good 
%                   for control purposes or scaling a spectrum plot.
% ------------------------------------------------------------------------ 
% timeData          Struct with minimum fields: 
% 
%                   * .impulseResponses     [Channels X Samples]          
%                   * .FS
%                   * .radius               Array radius
%                   * .averageAirTemp       Temperature in [�C] 
% 
% 
% FFToversize       FFToversize rises the FFT Blocksize.   [default = 2] 
%                   A FFT of the blocksize (FFToversize*NFFT) is applied 
%                   to the time domain data,  where  NFFT is determinded  
%                   as the next power of two of the signalSize  which is
%                   signalSize = (lastSample-firstSample).                   
% 
%                   The function will pick a window of 
%                   (lastSample-firstSample) for the FFT:
% 
% firstSample       First time domain sample to be included 
% lastSample        Last time domain sample to be included
%                   If firstSample and lastSample are not defined 
%                   the full IR will be transformed:
%                   [ default: firstSample = 1;
%                     lastSample = size(timeData.impulseResponses,2);]
% 
% Call this function with a running window (firstSample+td->lastSample+td) 
% iteration increasing td to obtain time slices. This way you resolve the
% temporal information within the captured sound field.                  



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
% Copyright (C)2011 bBrn - benjamin Bernsch�tz [1,2] et al.(�)   
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
% (�) Christoph P�rschmann [1]   christoph.poerschmann 'at' fh-koeln.de
%     Sascha Spors         [2,3] sascha.spors 'at' telekom.de  
%     Stefan Weinzierl     [2]   stefan.weinzierl 'at' tu-berlin.de
%     Nils Peters                nils 'at' icsi.berkeley.edu


function [fftData, kr, f] = sofia_fdt(timeData, FFToversize, firstSample, lastSample)

disp('SOFiA F/D/T - Frequency Domain Transform R11-1220');

try
    impulseResponses = timeData.impulseResponses;
    FS               = timeData.FS;
    temperature      = timeData.averageAirTemp;
    radius           = timeData.radius;
catch
    fftData=[]; 
    kr=[]; 
    f=[];
    error('Input data not valid');
    return
end

if nargin < 2
   FFToversize = 1;
end

FFToversize = round (FFToversize);
if FFToversize < 1 
   error('FFToversize must be greater or equal 1')   
end

if nargin < 4
   firstSample = 1;
   lastSample  = size(impulseResponses,2);
end

if lastSample < firstSample 
   error('lastSample must be greater than firstSample')   
end

if firstSample < 1
   error('firstSample must be greater than zero')   
end

if lastSample > size(impulseResponses,2)
   error('lastSample out of range')    
end


totalSamples     = lastSample - firstSample + 1;
impulseResponses = impulseResponses(:, firstSample:lastSample);
NFFT             = 2^(nextpow2(totalSamples));
fftData          = fft(impulseResponses, NFFT * FFToversize, 2);
fftData          = fftData(:,1:end/2+1);
fftData          = cast(fftData ,'double');

f  = linspace(0, FS/2, (NFFT*FFToversize)/2+1);
c  = 331.5+0.6*temperature;
k  = 2*pi*f/c;
kr = k*radius;


