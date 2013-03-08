% /// ASAR Research Group
% 
%     Cologne University of Applied Sciences
%     Technical University of Berlin
%     Deutsche Telekom Laboratories
%     WDR Westdeutscher Rundfunk
% 
% SOFiA sound field analysis
% 
% Impulse response reconstruction R11-1220
% 
% Copyright (C)2011 by bBrn - benjamin Bernschütz   
% 
% This file is part of the SOFiA toolbox under GNU General Public License
% 
%  
% impulseResponses = sofia_makeIR(Y, [win], [resampleFactor])
% ------------------------------------------------------------------------     
% impulseResponses   Reconstructed impulse response
%                    Columns : Index / Channel: IR1, IR2, ... IRn
%                    Rows    : Impulse responses (time domain)
% ------------------------------------------------------------------------              
% Y                  Frequency domain FFT data for multiple channels
%                    Columns : Index / Channel
%                    Rows    : FFT data (frequency domain)
% 
% [win]              Window IR tail [0...1] with a HANN window
%                    0    off
%                    0-1  window coverage (1 full, 0 off)
%                    [default 1/8: 1/8 of the IR length is windowed]
%                    ! Signal Processing Toolbox required
%                   
% [resampleFactor]   Optional resampling: Resampling factor 
%                    e.g. FS_target/FS_source
%                    Resampling is done using MATLAB RESAMPLE 
%                    (See MATLAB documentation for more details)
%                    ! Signal Processing Toolbox required  
% 
%  This function recombines impulse responses for multiple channels from
% frequency domain data. It is made to work with half-sided spectrum FFT 
% data.  The impulse responses can be windowed.  The IFFT blocklength is 
% determined by the Y data itself:
% 
% Y should have a size [NumberOfChannels x ((2^n)/2)+1] with n=[1,2,3,...] 
% and the function returns [NumberOfChannels x resampleFactor*2^n] samples.
% 
% Dependencies: MATLAB Signal Processing Toolbox required for 
%               windowing and resampling.    
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
% 



function impulseResponses = sofia_makeIR(Y, win, resampleFactor)

disp('SOFiA Impulse Response Reconstruction R11-1220');

if nargin == 0
   error('Arguments missing: impulseResponses = sofia_makeIR(Y, [win], [resampleFactor])');
end

if nargin < 2
   win = 1/8;
end

if win > 1
   error('Argument win must be in the range 0 to 1.')
end

signalProcessingExists = 0;
toolboxes = ver;
for i = 1: size(toolboxes,2)
    if ~isempty(strfind(toolboxes(i).Name,'Signal Processing Toolbox'))
       signalProcessingExists = 1;   
       break
    end
end

if ~signalProcessingExists 
   disp('WARNING: No signal processing toolbox found. Windowing and resampling disabled.') 
end

Y(:,1) = Y(:,2); 
y = [Y(:,:), conj(fliplr(Y(:,2:end-1)))];
y = real(ifft(y,[],2)); 

if signalProcessingExists && win ~= 0
   winfkt = window(@hann,round(size(Y,2)*2*win))';
   winfkt = winfkt(round(end/2+1):end);
   winfkt = [ones(1,2*size(Y,2)-size(winfkt,2)-2), winfkt];
   winfkt = repmat(winfkt, size(Y,1), 1);     
   y = y.*winfkt;
end

if nargin == 3 && signalProcessingExists 
    impulseResponses = resample(y',resampleFactor,1)';
else
    impulseResponses = y;   
end

