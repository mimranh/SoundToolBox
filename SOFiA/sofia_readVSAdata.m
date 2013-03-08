% /// ASAR Research Group
% 
%     Cologne University of Applied Sciences
%     Technical University of Berlin
%     Deutsche Telekom Laboratories
%     WDR Westdeutscher Rundfunk
% 
% SOFiA sound field analysis
% 
% VariSphear -> SOFiA data import R11-1220
% 
% For more information on the VariSphear Array 
% system visit: http://varisphear.fh-koeln.de/
%
% Copyright (C)2011 by bBrn - benjamin Bernschütz   
% 
% This file is part of the SOFiA toolbox under GNU General Public License
% 
%  
% [timeDataCH1, timeDataCH2] = readVSA_data([downSample], ...
%                                    [normalize], [directory])
% ------------------------------------------------------------------------     
% timeDataCH1/CH2     Structs with fields:
%
%                     .impulseResponses     [Channels x Samples]          
%                     .FS
%                     .quadratureGrid       [AZ1 EL1 W1; ...; AZn ELn Wn]
%                     .metaData             Cell Array, VSA Metadata
%                     .downSample
%                     .averageAirTemp       Temperature in DEG 
%                     .irOverlay            Plot this for a good total 
%                                           overview of the dataset.
% ------------------------------------------------------------------------              
% downSample         Downsampling factor  [default = 1: No downsampling]
%                    Downsampling is done using DECIMATE and a FIR low 
%                    pass filter of order 30. See MATLAB documentation 
%                    for more information. 
%                    !!! MATLAB Signal Processing Library required
%
% normalize          Normalize flag 1:on, 0:off         [default = 1: on]       
%                    Normalizes the impulse responses with respect to the 
%                    absolute maximum value within the complete dataset.
%
% directory          VSA dataset directory. If not defined a user dialog
%                    opens to pick a directory.
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


function [timeDataCH1, timeDataCH2] = sofia_readVSAdata(downSample, normalize, directory)

disp('SOFiA VariSphear -> SOFiA data import R11-1220');

air_temperature = 20; %DEFAULT if no temperatures found in vsd-dataset
    
if nargin == 0
   downSample = 1; 
end

if nargin < 2
   normalize = 1; 
end

if nargin < 3
    directory = uigetdir();
    if directory == 0
       error('No directory picked');
    end
end

files=dir(fullfile(directory, '*.mat'));        
arraycounterCH1=0;
arraycounterCH2=0;

Temp1=[];
Temp2=[];

    for filecounter=1:size(files,1)       

        imported = files(filecounter).name;
        vsd=[];

        if ~isempty(imported)             
           if mod(filecounter,5) == 0
              fprintf('|');        
           end           
           if mod(filecounter, 200) == 0
              fprintf('\n'); 
           end
           
           load(fullfile(directory, files(filecounter).name), 'vsd')
           
           if ~isempty(vsd)
               
               if vsd.CH == 1 || vsd.HalfSphere == 1% -------------------------------------------- CH1
                   
                   if strcmp(vsd.StorageType,'MAT')                     
                        if downSample ~= 1
                           timeDataCH1.impulseResponses(arraycounterCH1+1,:)=decimate(cast(vsd.ImpulseResponse,'double'),downSample ,'FIR');
                           timeDataCH1.FS=vsd.FS/downSample;
                        else
                           timeDataCH1.impulseResponses(arraycounterCH1+1,:)=vsd.ImpulseResponse;
                           timeDataCH1.FS=vsd.FS;                                
                        end                   
                   else 
                        try
                            irdata = wavread(fullfile(directory,vsd.ImpulseResponse(5:end)))';
                        catch
                            error(['ERROR - Data file missing: ', vsd.ImpulseResponse(5:end)]);
                        end
                        
                        if downSample ~= 1
                           timeDataCH1.impulseResponses(arraycounterCH1+1,:) = decimate(cast(irdata,'double'),downSample ,'FIR');
                           timeDataCH1.FS=vsd.FS/downSample;
                        else
                           timeDataCH1.impulseResponses(arraycounterCH1+1,:) = cast(irdata,'double');
                           timeDataCH1.FS=vsd.FS;                                
                        end
                   end
                   
                   timeDataCH1.metaData{arraycounterCH1+1}=rmfield(vsd,'ImpulseResponse'); 
                   timeDataCH1.radius = vsd.Radius;
                   timeDataCH1.quadratureGrid(arraycounterCH1+1,1)=vsd.Azimuth*pi/180;
                   timeDataCH1.quadratureGrid(arraycounterCH1+1,2)=vsd.Elevation*pi/180;
                   timeDataCH1.quadratureGrid(arraycounterCH1+1,3)=vsd.GridWeight;
                   
                   if isfield(vsd,'Temperature')
                      Temp1(arraycounterCH1+1)=vsd.Temperature;                           
                   else
                      Temp1(arraycounterCH1+1)=air_temperature;
                   end   
                   arraycounterCH1=arraycounterCH1+1;
               
               else               % -------------------------------------------- CH2
                   
                   if strcmp(vsd.StorageType,'MAT')                     
                        if downSample ~= 1
                           timeDataCH2.impulseResponses(arraycounterCH2+1,:)=decimate(cast(vsd.ImpulseResponse,'double'),downSample ,'FIR');
                           timeDataCH2.FS=vsd.FS/downSample;
                        else
                           timeDataCH2.impulseResponses(arraycounterCH2+1,:)=vsd.ImpulseResponse;
                           timeDataCH2.FS=vsd.FS;                                
                        end
                   else 
                        try
                            if (ismac)
                               irdata = wavread([directory,'/',vsd.ImpulseResponse(5:end)])';
                            else
                               irdata = wavread([directory,'\',vsd.ImpulseResponse(5:end)])';
                            end
                        catch
                            error(['ERROR - Data file missing: ', vsd.ImpulseResponse(5:end)]);
                        end
                        
                        if downSample ~= 1
                           timeDataCH2.impulseResponses(arraycounterCH2+1,:) = decimate(cast(irdata,'double'),downSample ,'FIR');
                           timeDataCH2.FS=vsd.FS/downSample;
                        else
                           timeDataCH2.impulseResponses(arraycounterCH2+1,:) = cast(irdata,'double');
                           timeDataCH2.FS=vsd.FS;                                
                        end
                   end
                
                   timeDataCH2.metaData{arraycounterCH2+1}=rmfield(vsd,'ImpulseResponse'); 
                   timeDataCH2.radius = vsd.Radius;
                   timeDataCH2.quadratureGrid(arraycounterCH2+1,1)=vsd.Azimuth*pi/180;
                   timeDataCH2.quadratureGrid(arraycounterCH2+1,2)=vsd.Elevation*pi/180;
                   timeDataCH2.quadratureGrid(arraycounterCH2+1,3)=vsd.GridWeight;
                   
                   if isfield(vsd,'Temperature')
                      Temp2(arraycounterCH2+1)=vsd.Temperature;                           
                   else
                      Temp2(arraycounterCH2+1)=air_temperature;
                   end     
                   
                   arraycounterCH2=arraycounterCH2+1; 
               end                      
          end
        end   
    end
    
  

if arraycounterCH1>0 
   
   fprintf('\n\n');
   disp([num2str(arraycounterCH1),' spatial sampling points imported.']);   
   
   if normalize == 1
      timeDataCH1.impulseResponses = timeDataCH1.impulseResponses./max(max(abs(timeDataCH1.impulseResponses)));     
   end
   
   timeDataCH1.downSample=downSample;
   timeDataCH1.averageAirTemp = mean(Temp1);    
   timeDataCH1.irOverlay=sum(timeDataCH1.impulseResponses,1);
   timeDataCH1.irOverlay=abs(timeDataCH1.irOverlay/max(abs(timeDataCH1.irOverlay)));
   
  
   if arraycounterCH2>0       
      if normalize == 1
          timeDataCH2.impulseResponses=timeDataCH2.impulseResponses./max(max(abs(timeDataCH2.impulseResponses))); 
      end
      timeDataCH2.downSample=downSample;
      timeDataCH2.averageAirTemp = mean(Temp2);    
      timeDataCH2.irOverlay=sum(timeDataCH2.impulseResponses,1);
      timeDataCH2.irOverlay=abs(timeDataCH2.irOverlay/max(abs(timeDataCH2.irOverlay)));

   else
      timeDataCH2=[]; 
   end
      
else
   disp(['Nothing found to import.']); 
   timeDataCH1=[];
   timeDataCH2=[];
end    

fprintf('\n');
