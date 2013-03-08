% /// ASAR Research Group
% 
% Cologne University of Applied Sciences
% Technical University of Berlin
% Deutsche Telekom Laboratories
% WDR Westdeutscher Rundfunk
% 
% SOFiA sound field analysis
% 
% Make Matrix 3D Data (Pressures) R11-1220
%
% Preparing pressure data for sofia_visual3d() 
% 
% Copyright (C)2011 by bBrn - benjamin Bernschütz   
%                             rockzentrale 'AT' me.com
%                             +49 171 4176069 Germany  
% 
% This file is part of the SOFiA toolbox under GNU General Public License
% 
%  
% mtxData = sofia_makeMtx(N, Pnm, dn, krIndex) 
% ------------------------------------------------------------------------     
% mtxData   SOFiA 3D-matrix-data in 1° steps
% ------------------------------------------------------------------------              
% N         Order of the spatial fourier transform     [default = 3]
% Pnm       Spatial Fourier Coefficients (from S/T/C)
% dn        Modal Radial Filters (from M/F)
% krindex   Index of kr Vector                         [default = 1]
%
% Dependencies: SOFiA P/D/C
% 
% The file generates a SOFiA mtxData Matrix of 181x360 pixels for the
% visualisation with sofia_visual3d() in 1° Steps (65160 plane waves).
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


function mtxData = sofia_makeMTX(N, Pnm, dn, krIndex)  

disp('SOFiA Visual Matrix Generator Wrapper R11-1220');

angle = zeros(65160,2);
cnt=1;
for Elev=0:180
    for Azim=0:359
        angle(cnt,:) = [Azim Elev];
        cnt = cnt+1;
    end
end

angle=angle*pi/180;
Y = sofia_pdc(N, angle, Pnm(:,krIndex), dn(:,krIndex));

cnt=1;
mtxData=zeros(181,360);
for Elev=1:181
    for Azim=1:360
        mtxData(Elev, Azim)=Y(cnt);
        cnt=cnt+1;
    end
end

